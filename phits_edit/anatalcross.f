!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_crsreg(m,ntf,
     &                             np,ne,nt,na,nr,nm,
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
      integer           nt
      integer           na
      integer           nr
      integer           nm

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,na,nt,nr,nm,2)
      integer           nfile
      double precision  tranatal(np,ne,na,nt,nr,nm,2*nfile)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ip
      integer           ie
      integer           it
      integer           ia
      integer           ir
      integer           im

      double precision  Xa, sigx
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do im = 1, nm
      do ir = 1, nr
       do ia = 1, na
        do it = 1, nt
         do ie = 1, ne
          do ip = 1, np

             if( trRES(ip,ie,ia,it,ir,im,1) .gt. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(ip,ie,ia,it,ir,im,1),
     &                        trRES(ip,ie,ia,it,ir,im,2),
     &                        1.0d+0)
              trRES(ip,ie,ia,it,ir,im,1) = Xa
              trRES(ip,ie,ia,it,ir,im,2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(ip,ie,ia,it,ir,im,iat(1,ntf)) =
     &           trRES(ip,ie,ia,it,ir,im,1)

              ! (sig_xj)**2
              tranatal(ip,ie,ia,it,ir,im,iat(2,ntf)) =
     &        (  trRES(ip,ie,ia,it,ir,im,2)
     &         * trRES(ip,ie,ia,it,ir,im,1) )**2

              ! sig_x
              tranatal(ip,ie,ia,it,ir,im,iat(2,ntf)) = dsqrt(
     &        tranatal(ip,ie,ia,it,ir,im,iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(ip,ie,ia,it,ir,im,iat(2,ntf)) =
     &       (tranatal(ip,ie,ia,it,ir,im,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(ip,ie,ia,it,ir,im,iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(ip,ie,ia,it,ir,im,iat(1,ntf)) =
     &        tranatal(ip,ie,ia,it,ir,im,iat(1,ntf)) * resc2(m)

            if( tranatal(ip,ie,ia,it,ir,im,iat(1,ntf)) > cmax )
     &            cmax = tranatal(ip,ie,ia,it,ir,im,iat(1,ntf))
            if( tranatal(ip,ie,ia,it,ir,im,iat(1,ntf)) < cmin )
     &            cmin = tranatal(ip,ie,ia,it,ir,im,iat(1,ntf))

          end do     ! ip loop end
         end do      ! ie loop end
        end do       ! it loop end
       end do        ! ia loop end
      end do         ! ir loop end
      end do         ! im loop end

! sumover
      call anatal_calc_anova_sub(m,ntf,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return
      end subroutine anatal_calc_anova_crsreg



!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_crsrz(m,ntf,
     &                             np,ne,na,nt,nr,nz,nm,
     &                             sumfactor,weightRate,
     &                             trRES,tranatal,nfile,
     &                             tzRES,tzanatal,
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
      integer           na
      integer           nt
      integer           nr
      integer           nz
      integer           nm

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,na,nt,(nr+1)*nz,nm,2)
      integer           nfile
      double precision  tranatal(np,ne,na,nt,(nr+1)*nz,nm,2*nfile)
      double precision  tzRES   (np,ne,na,nt,nr*(nz+1),nm,2)
      double precision  tzanatal(np,ne,na,nt,nr*(nz+1),nm,2*nfile)


      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           iz
      integer           ir
      integer           it
      integer           ia
      integer           ie
      integer           ip
      integer           im

      double precision  Xa, sigx
*-----------------------------------------------------------------------
      integer           irf
      integer           izf
      integer           iat, iad, iaf

      irf(ir,iz) = ir + ( iz - 1 ) * ( nr + 1 )
      izf(ir,iz) = iz + ( ir - 1 ) * ( nz + 1 )

      iat(iad,iaf) = iad + (iaf-1) * 2

!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
*        z-crossing
*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do im = 1, nm
      do iz = 1, nz + 1
       do ir = 1, nr
        do it = 1, nt
         do ia = 1, na
          do ie = 1, ne
           do ip = 1, np

             if( tzRES(ip,ie,ia,it,izf(ir,iz),im,1) .gt. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        tzRES(ip,ie,ia,it,izf(ir,iz),im,1),
     &                        tzRES(ip,ie,ia,it,izf(ir,iz),im,2),
     &                        1.0d+0)
              tzRES(ip,ie,ia,it,izf(ir,iz),im,1) = Xa
              tzRES(ip,ie,ia,it,izf(ir,iz),im,2) = sigx
             end if

              ! X_bar = Xj_bar
              tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(1,ntf)) =
     &           tzRES(ip,ie,ia,it,izf(ir,iz),im,1)

              ! (sig_xj)**2
              tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(2,ntf)) =
     &        (  tzRES(ip,ie,ia,it,izf(ir,iz),im,2)
     &         * tzRES(ip,ie,ia,it,izf(ir,iz),im,1) )**2

              ! sig_x
              tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(2,ntf)) = dsqrt(
     &        tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(2,ntf)) =
     &       (tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(1,ntf)) =
     &        tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(1,ntf)) * resc2(m)

            if( tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(1,ntf)) > cmax )
     &            cmax = tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(1,ntf))
            if( tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(1,ntf)) < cmin )
     &            cmin = tzanatal(ip,ie,ia,it,izf(ir,iz),im,iat(1,ntf))

           end do    ! ip loop end
          end do     ! ie loop end
         end do      ! ia loop end
        end do       ! it loop end
       end do        ! ir loop end
      end do         ! iz loop end
      end do         ! im loop end

! sumover
      call anatal_calc_anova_tz_sub(m,ntf,
     &                               resc2(m),resc3(m) )

*-----------------------------------------------------------------------
*        r-crossing
*-----------------------------------------------------------------------

      do im = 1, nm
      do iz = 1, nz
       do ir = 1, nr + 1
        do it = 1, nt
         do ia = 1, na
          do ie = 1, ne
           do ip = 1, np

             if( trRES(ip,ie,ia,it,irf(ir,iz),im,1) .gt. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(ip,ie,ia,it,irf(ir,iz),im,1),
     &                        trRES(ip,ie,ia,it,irf(ir,iz),im,2),
     &                        1.0d+0)
              trRES(ip,ie,ia,it,irf(ir,iz),im,1) = Xa
              trRES(ip,ie,ia,it,irf(ir,iz),im,2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(ip,ie,ia,it,irf(ir,iz),im,iat(1,ntf)) =
     &           trRES(ip,ie,ia,it,irf(ir,iz),im,1)

              ! (sig_xj)**2
              tranatal(ip,ie,ia,it,irf(ir,iz),im,iat(2,ntf)) =
     &        (  trRES(ip,ie,ia,it,irf(ir,iz),im,2)
     &         * trRES(ip,ie,ia,it,irf(ir,iz),im,1) )**2

              ! sig_x
              tranatal(ip,ie,ia,it,irf(ir,iz),im,iat(2,ntf)) = dsqrt(
     &        tranatal(ip,ie,ia,it,irf(ir,iz),im,iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(ip,ie,ia,it,irf(ir,iz),im,iat(2,ntf)) =
     &       (tranatal(ip,ie,ia,it,irf(ir,iz),im,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(ip,ie,ia,it,irf(ir,iz),im,iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(ip,ie,ia,it,irf(ir,iz),im,iat(1,ntf)) =
     &        tranatal(ip,ie,ia,it,irf(ir,iz),im,iat(1,ntf)) * resc2(m)


           end do    ! ip loop end
          end do     ! ie loop end
         end do      ! ia loop end
        end do       ! it loop end
       end do        ! ir loop end
      end do         ! iz loop end
      end do         ! im loop end

! sumover
      call anatal_calc_anova_tr_sub(m,ntf,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return
      end subroutine anatal_calc_anova_crsrz



!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_crsxyz(m,ntf,
     &                             nt,np,na,ne,ny,nx,nz,nm,
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

      integer           nt
      integer           np
      integer           na
      integer           ne
      integer           ny
      integer           nx
      integer           nz
      integer           nm

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,na,nt,nx*ny*(nz+1),nm,2)
      integer           nfile
      double precision  tranatal(np,ne,na,nt,nx*ny*(nz+1),nm,2*nfile)

      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           iz
      integer           ix
      integer           iy

      integer           ie
      integer           ia
      integer           ip
      integer           it
      integer           im
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

      do im = 1, nm
      do iz = 1, nz + 1
       do ix = 1, nx

        do iy = 1, ny
         do ie = 1, ne
          do ia = 1, na
           do ip = 1, np
            do it = 1, nt

             if( trRES(ip,ie,ia,it,icf(ix,iy,iz),im,1) .gt. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &                        trRES(ip,ie,ia,it,icf(ix,iy,iz),im,2),
     &                        1.0d+0)
              trRES(ip,ie,ia,it,icf(ix,iy,iz),im,1) = Xa
              trRES(ip,ie,ia,it,icf(ix,iy,iz),im,2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf)) =
     &           trRES(ip,ie,ia,it,icf(ix,iy,iz),im,1)

              ! (sig_xj)**2
              tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(2,ntf)) =
     &        (  trRES(ip,ie,ia,it,icf(ix,iy,iz),im,2)
     &         * trRES(ip,ie,ia,it,icf(ix,iy,iz),im,1) )**2

              ! sig_x
              tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(2,ntf)) = dsqrt(
     &        tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(2,ntf)) =
     &       (tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf)) =
     &        tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf))*resc2(m)

            if(tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf)) > cmax)
     &          cmax = tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf))
            if(tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf)) < cmin)
     &          cmin = tranatal(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf))

            end do   ! it loop end
           end do    ! ip loop end
          end do     ! ia loop end
         end do      ! ie loop end
        end do       ! iy loop end
       end do        ! ix loop end
      end do         ! iz loop end
      end do         ! im loop end

! sumover
      call anatal_calc_anova_sub(m,ntf,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return
      end subroutine anatal_calc_anova_crsxyz



************************************************************************
*                                                                      *
      subroutine anatal_psufreg(m,np,nr,ne,na,nt,nm,mr,kr,ar,eb,ab,tb,
     &                   tr,nfile,weightRate,idasa,manatally)
*                                                                      *
*       output of region crossing surface tally of                     *
*       current spectrum and flux.                                     *
*       created by T.Miura on 2021/09/14                               *
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

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      double precision resc2   ! W
      double precision resc3   ! N
      common /restart/ resc2(itlmax), resc3(itlmax)

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)

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
      common /tall35/ iteps(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /istcut/ ist_cut, ist_bat
*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      dimension   kr(mr+1)
      dimension   eb(ne+1)
      dimension   ew(ne)
      dimension   ab(na+1)
      dimension   tb(nt+1)
!      dimension   aw(na)
!      dimension   tw(nt)
      dimension   ar(nr)
!      dimension   vl(nr,1,1)
      real(8),allocatable :: aw(:),tw(:),vl(:,:,:)

      dimension   tr(np,ne,na,nt,nr,nm,2*nfile)

      integer     nfile, manatally
!      dimension   weightRate(nfile), rdata(2,nfile)
      dimension   weightRate(nfile)
      real(8),allocatable :: rdata(:,:)

      integer irst,nrst, itmpdata
!!      dimension   anatalrst(np,ne+1,na+1,nt+1,nr*1*1,nm+1,2*nfile+3)
! sumover
!      dimension   anatalrst(np,ne+1,na+1,nt+1,(nr+1)*1*1,nm+1,2*nfile+3)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)

      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9  ! S.H. 2021.8.9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*21,cijaxs*21,cijaxs2*21
      integer iDaxis
!      dimension   dlr(nr+1) !  integer ar(nr) ---> double dlr(nr)
      real(8),allocatable :: dlr(:)

*-----------------------------------------------------------------------

      dimension tott(6,2)

      character hsunit(16)*32

      data hsunit( 1) / '[1/cm^2/source]                 '/
      data hsunit( 2) / '[1/cm^2/MeV/source]             '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]        '/
      data hsunit( 4) / '[1/cm^2/sr/source]              '/
      data hsunit( 5) / '[1/cm^2/MeV/sr/source]          '/
      data hsunit( 6) / '[1/cm^2/Lethargy/sr/source]     '/
      data hsunit(11) / '[1/cm^2/nsec/source]            '/
      data hsunit(12) / '[1/cm^2/MeV/nsec/source]        '/
      data hsunit(13) / '[1/cm^2/Lethargy/nsec/source]   '/
      data hsunit(14) / '[1/cm^2/sr/nsec/source]         '/
      data hsunit(15) / '[1/cm^2/MeV/sr/nsec/source]     '/
      data hsunit(16) / '[1/cm^2/Lethargy/sr/nsec/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

      character rpa*1
      data rpa /'}'/

      character dum1*10000
      character dum2*10000
      character dum3*10000

      character cname*7
      character dname*7

      character aname*3

      character yen*1

      dimension dt_one(1)
      data dt_one/1.0d0/
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
      yen  = char(92)

      allocate (aw(na),tw(nt),vl(nr,1,1))
      allocate (rdata(2,nfile))
      allocate (dlr(nr+1))

      nrst = 2*nfile+3 ! S.H. 2021.8.15

      allocate (anatalrst(np,ne+1,na+1,nt+1,(nr+1)*1*1,nm+1,nrst))

*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 2) = '[1/cm^2/(MeV/n)/source]         '
         hsunit( 5) = '[1/cm^2/(MeV/n)/sr/source]      '
         hsunit(12) = '[1/cm^2/(MeV/n)/nsec/source]    '
         hsunit(15) = '[1/cm^2/(MeV/n)/sr/nsec/source] '
      end if

*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------
*        current or flux
*-----------------------------------------------------------------------

            if( itout(m) .eq. 1 ) then

               cname = '  Flux '
               dname = '  flux '

            else if( itout(m) .eq. 2 ) then

               cname = 'Current'
               dname = 'current'

            else if( itout(m) .eq. 3 ) then

               cname = ' F-Curr'
               dname = ' f-curr'

            else if( itout(m) .eq. 4 ) then

               cname = ' B-Curr'
               dname = ' b-curr'

            else if( itout(m) .eq. 5 ) then

               cname = ' O-Curr'
               dname = ' o-curr'

            else if( itout(m) .eq. 6 ) then

               cname = 'OF-Curr'
               dname = 'of-curr'

            else if( itout(m) .eq. 7 ) then

               cname = 'OB-Curr'
               dname = 'ob-curr'

            else if( itout(m) .eq. 8 ) then

               cname = ' A-Curr'
               dname = ' a-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 9 ) then

               cname = 'OA-Curr'
               dname = 'oa-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

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
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq.  2 .or. itunt(m) .eq.  5 .or.
     &               itunt(m) .eq. 12 .or. itunt(m) .eq. 15 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 13 .or. itunt(m) .eq. 16 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            end if

*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

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
*           itunt(m) = 11, 12, 13, 14 : /.../nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .gt. 10 ) then

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

ccse 2021.08 add (use anatal_rearrange sub.)
               do ir = 1, nr
                  dlr(ir) = dble(ir)
               end do
               dlr(nr+1) = dble(nr)

               do ir = 1, nr
                  vl(ir,1,1) = ar(ir)
               end do
               ar_sum = sum(ar(:))

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*-----------------------------------------------------------------------
                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

            do 100 im = 1, nm
            do 100 ir = 1, nr
            do 100 ia = 1, na
            do 100 it = 1, nt
            do 100 ie = 1, ne
            do 100 ip = 1, np
            do ntf = 1, nfile

               if( tr(ip,ie,ia,it,ir,im,iat(1,ntf)) .gt. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,ia,it,ir,im,iat(1,ntf)),
     &                            tr(ip,ie,ia,it,ir,im,iat(2,ntf)),
     &                            rtfac(m)/ar(ir)/ew(ie)/aw(ia)/tw(it))

                  tr(ip,ie,ia,it,ir,im,iat(1,ntf)) = Xa
                  tr(ip,ie,ia,it,ir,im,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

               else

                  isdz = 1
                  tr(ip,ie,ia,it,ir,im,iat(2,ntf)) = 0.0

               end if

! sumover
               call psufreg_sumover_stdev_ntf(0,m,ntf,
     &               ip,ie,ia,it,ir,im,
     &               rtfac(m),ew(ie),aw(ia),tw(it),ar(ir),
     &               ew_sum,aw_sum,tw_sum,ar_sum)

            end do
  100       continue

!$OMP PARALLEL
!$OMP& private(ipomp,ir,ia,im,it,ie,ip,ntf,isdz,rdata,answer,rerr)
!$OMP& private(iax,ir_a,ia_a,it_a,ie_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
            do 101 im = 1, nm
            do 101 ir = 1, nr
            do 101 ia = 1, na
            do 101 it = 1, nt
            do 101 ie = 1, ne
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,ia,it,ir,im,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,ia,it,ir,im,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(ip,ie,ia,it,ir,im,1) = answer
              anatalrst(ip,ie,ia,it,ir,im,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,ia,it,ir,im,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,ia,it,ir,im,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(ip,ie,ia,it,ir,im,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,ia,it,ir,im,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,ia,it,ir,im,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,ia,it,ir,im,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,ia,it,ir,im,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,ia,it,ir,im,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,ir,im,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,ir,im,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,ir,im,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(ip,ie,ia,it,ir,im,iat(ioe,ntf))
     &              = tr(ip,ie,ia,it,ir,im,iat(ioe,ntf))
                end do
               end do

            end if


!sumover
            do iax = 1,itaxn(m)

              ie_a = ie
              ia_a = ia
              it_a = it
              ir_a = ir
              ido_ana = 0
              if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 14)
     &             .and. ie == 1) then
                ie_a = ne + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 2 .and. ir == 1) then
                ir_a = nr + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 9 .and. it == 1) then
                it_a = nt + 1
                ido_ana = 1
              elseif((itaxs(m,iax) == 8 .or. itaxs(m,iax) == 10)
     &             .and. ia == 1) then
                ia_a = na + 1
                ido_ana = 1
              endif

              if(ido_ana == 1) then

                 call get_crsreg_tr_sum_data(m,iax,nfile,
     &                 ip,ie,ia,it,ir,im,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   anatalrst(ip,ie_a,ia_a,it_a,ir_a,im,1) = answer
                   anatalrst(ip,ie_a,ia_a,it_a,ir_a,im,2) = rerr

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  anatalrst(ip,ie_a,ia_a,it_a,ir_a,im,1) ! mean
     &              = fmval
                 if ( fmval .gt. 0.0d0 ) then
                  anatalrst(ip,ie_a,ia_a,it_a,ir_a,im,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                 else
                   anatalrst(ip,ie_a,ia_a,it_a,ir_a,im,2) ! total uncertainty
     &                 = 0.0d0
                 end if

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      anatalrst(ip,ie_a,ia_a,it_a,ir_a,im,iat(ioe,ntf))
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

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do iax = 1, itfln(m)

         if( itmdp(m,0) .eq. 0 ) then

            if( ( itall .eq. 2 .and. nobch .lt. maxbch )
     &           .or. itall .eq. 4 ) then

               call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &              nobch,maxbch,npe)

            else

               fname = ctfln(m,iax)

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
*     input echo
*-----------------------------------------------------------------------

               call tcrsech(iot,m,iax,1)

*-----------------------------------------------------------------------
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        energy, reg, time, angle axes
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1
     &        .or. itaxs(m,iax) .eq. 2
     &        .or. itaxs(m,iax) .eq. 9
     &        .or. itaxs(m,iax) .eq. 8
     &        .or. itaxs(m,iax) .eq. 10
     &        .or. itaxs(m,iax) .eq. 14 ) then ! T.Sato 2023/01/09

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 1   ! energy
           case ( 2 )
            iDaxis = 2   ! reg axis
           case ( 9 )
            iDaxis = 9   ! time
           case ( 8 )
            iDaxis = 12  ! cos of angle(p) axis
           case ( 10 )
            iDaxis = 13  ! the of angle(p) axis
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.15

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(itaty(m),iMeVperu,1,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      na,    nt,     nr,  0,  0,  nm, nrst,
     &         eb,ew,  ab,aw,  tb,tw, dlr,  dt_one,  dt_one,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

               inum = 0

                  j = 0
                  igm = ( mmmax - 1 ) * 2 + 1

          do ianataldata = 1, nanataldata

                  j = j + 1
                  ntrn = kr(j)
                  j = j + 1
                  mtrn = kr(j)
                  call echrg2(mtrn,kr(j+1),dum1,lng1,icmb,igm)
                  j = j + mtrn
                  j = j + 1
                  ntrn = kr(j)
                  j = j + 1
                  mtrn = kr(j)
                  call echrg2(mtrn,kr(j+1),dum2,lng2,icmb,igm)
                  j = j + mtrn


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


                  l = 0
                  dum3(l+1:l+9) = '#   no. ='
                  l = l + 9
                  write(dum3(l+1:l+5),'(i5)') inum
                  l = l + 5
                  dum3(l+1:l+9) = '   reg = '
                  l = l + 9
                  dum3(l+1:l+lng1) = dum1(1:lng1)
                  l = l + lng1
                  dum3(l+1:l+3) = ' - '
                  l = l + 3
                  dum3(l+1:l+lng2) = dum2(1:lng2)
                  l = l + lng2

                  write(iot,'(600a1)') (dum3(i:i),i=1,l)

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
               write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

        if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 14 ) then ! T.Sato 2023/01/09
               if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 16 .or.
     &             itety(m) .eq.  3 .or. itety(m) .eq.  5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 9 ) then
               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then
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
c S.H. added IF statement below for epsout=2 (2016.7.25)
            if ( iteps(m) .ne. 2 ) then
             if( itaxs(m,iax) .eq. 2 ) then ! reg axis
               write(iot,'( "h: x",7x,"n",5x,"n",11x,
     &                    1000(a10,"),",a4," n   "))')
     &                     ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

             else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if
            else
             if( itaxs(m,iax) .eq. 2 ) then ! reg axis
               write(iot,'( "h: x",7x,"n",5x,"n",11x,
     &         1000(a1,i0,a9,"),",a4
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &         ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             else
               write(iot,'( "h: n",12x,"x",12x,
     &         1000(a1,i0,a9,"),hh0",a3
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &         ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

            if( itaxs(m,iax) .eq. 2 ) then ! reg axis
               write(iot,'( "#  num    reg     area    ",
     &                    1000(a1,2x,a8,4x,"r.err "))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            else
               write(iot,'( "#  lower        upper  ",3x,
     &                    1000(a1,2x,a8,4x,"r.err "))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            end if

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
             if( itaxs(m,iax) .eq. 2 ) then ! reg axis
               write(iot,'( "h: x",7x,"n",5x,"n",11x,
     &                    1000(a10,"),",a4," n n n "))')
     &                     ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

             else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n n n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if
            else
             if( itaxs(m,iax) .eq. 2 ) then ! reg axis
               write(iot,'("h: x",7x,"n",5x,"n",11x,
     &         1000(a1,i0,a9,"),",a4
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &         ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             else
               write(iot,'( "h: n",12x,"x",12x,
     &         1000(a1,i0,a9,"),hh0",a3
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &         ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

            if( itaxs(m,iax) .eq. 2 ) then ! reg axis
               write(iot,'( "#  num    reg     area    ",
     &                    1000(a1,2x,a8,4x,"r.err(tot, syst, stat)"))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            else
               write(iot,'( "#  lower        upper  ",3x,
     &                    1000(a1,2x,a8,4x,"r.err(tot, syst, stat)"))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            end if

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

           voll = 0.0d0

           njaxs = nijaxs + 1

           do ijaxs = 1, nijaxs

            if ( manatally .eq. 0 ) then ! user defined analysis
             if( itaxs(m,iax) .eq. 2 ) then
                write(iot,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0pf8.4))')
     &               ijaxs, int(dlr(ijaxs)), fgaxs(ijaxs),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05
             else
                write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05
             end if

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
             if( itaxs(m,iax) .eq. 2 ) then
                write(iot,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0p3f8.4))')
     &               ijaxs, int(dlr(ijaxs)), fgaxs(ijaxs),
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
            if( itaxs(m,iax) .eq.  1 ) then ! eng axis
                     if( itunt(m) .eq.  2 .or. itunt(m) .eq.  3 .or.
     &                   itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &                   itunt(m) .eq. 12 .or. itunt(m) .eq. 13 .or.
     &                   itunt(m) .eq. 15 .or. itunt(m) .eq. 16 ) then
                 tott(ip,1) = tott(ip,1) / voll
                     end if

            else if( itaxs(m,iax) .eq. 2 ) then ! reg axis
                     tott(ip,1) = tott(ip,1) / voll

            else if( itaxs(m,iax) .eq.  9 ) then ! time axis
                  if( itunt(m) .gt.  10 ) then
                     tott(ip,1) = tott(ip,1) / voll
                  end if

            else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then ! cos or the axis
                     if( itunt(m) .eq.  4 .or. itunt(m) .eq.  5 .or.
     &                   itunt(m) .eq.  6 .or.
     &                   itunt(m) .eq. 14 .or. itunt(m) .eq. 15 .or.
     &                   itunt(m) .eq. 16 ) then
                        tott(ip,1) = tott(ip,1) / voll
                     end if
            end if
               end do

! sumover
            do ip = 1, np
            do irst=1,2
              tott(ip,irst) = anataldata(ip,ianataldata,njaxs,irst)
            enddo
            enddo

           if ( manatally .ne. 2 ) then
               write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4))')
     &                     (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
           end if

                  l = 0
                  write(dum3(l+1:l+1),'(a1)') cha
                  l = l + 1
                  dum3(l+1:l+5) = 'no. ='
                  l = l + 5
                  write(dum3(l+1:l+5),'(i5)') inum
                  l = l + 5
                  dum3(l+1:l+9) = '   reg = '
                  l = l + 9
                  dum3(l+1:l+lng1) = dum1(1:lng1)
                  l = l + lng1
                  dum3(l+1:l+3) = ' - '
                  l = l + 3
                  dum3(l+1:l+lng2) = dum2(1:lng2)
                  l = l + lng2

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

            call prestart(m,iot) !OBINATA(2012.6.13)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

      end do

*-----------------------------------------------------------------------

      deallocate (anatalrst)
      deallocate (aw,tw,vl)
      deallocate (rdata)
      deallocate (dlr)

      return
      end subroutine anatal_psufreg



************************************************************************
*                                                                      *
      subroutine anatal_psufrz(m,np,nr,nz,ne,na,nt,nm,rm,zm,eb,ab,tb,
     &                  tr,tz,nfile,weightRate,idasa,manatally)
*                                                                      *
*       output r-z scoring mesh crossing current & flux                *
*       created by T.Miura on 2021/09/14                               *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
      common /tall09/ itpan(itlmax), itpat(itlmax,6,2),
     &                jtpat(itlmax,6,6,2)
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
      common /tall35/ iteps(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /istcut/ ist_cut, ist_bat
*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      dimension   rm(nr+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   ab(na+1)
      dimension   tb(nt+1)
!      dimension   ew(ne)
!      dimension   aw(na)
!      dimension   tw(nt)
      real(8),allocatable :: ew(:),aw(:),tw(:)

      dimension   tr(np,ne,na,nt,(nr+1)*nz,nm,2)
      dimension   tz(np,ne,na,nt,nr*(nz+1),nm,2)


      integer     nfile, manatally
      dimension   weightRate(nfile)

*-----------------------------------------------------------------------

      dimension tott(6,2)
      real(8),allocatable :: ar_r(:),ar_z(:)

      character hsunit(16)*32

      data hsunit( 1) / '[1/cm^2/source]                 '/
      data hsunit( 2) / '[1/cm^2/MeV/source]             '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]        '/
      data hsunit( 4) / '[1/cm^2/sr/source]              '/
      data hsunit( 5) / '[1/cm^2/MeV/sr/source]          '/
      data hsunit( 6) / '[1/cm^2/Lethargy/sr/source]     '/
      data hsunit(11) / '[1/cm^2/nsec/source]            '/
      data hsunit(12) / '[1/cm^2/MeV/nsec/source]        '/
      data hsunit(13) / '[1/cm^2/Lethargy/nsec/source]   '/
      data hsunit(14) / '[1/cm^2/sr/nsec/source]         '/
      data hsunit(15) / '[1/cm^2/MeV/sr/nsec/source]     '/
      data hsunit(16) / '[1/cm^2/Lethargy/sr/nsec/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character chl(6)*3
      data chl /'l  ','dr ','ub ','mg ','qrr','vbb'/
      character chm(6)*4
      data chm /'l3  ','d4r ','u5b ','m6g ','q7rr','v8bb'/

      character chp(6)*11  ! kitamura22/03/31
      character chq(6)*9   ! kitamura22/03/31

      character chb(6)*1
      data chb /6*' '/

*-----------------------------------------------------------------------

      character rpa*1
      data rpa /'}'/

      character cname*7
      character dname*7
      character aname*3

      character yen*1

      include 'samepage_include/samepage000.inc'

*-----------------------------------------------------------------------

         irf(ir,iz) = ir + ( iz - 1 ) * ( nr + 1 )
         izf(ir,iz) = iz + ( ir - 1 ) * ( nz + 1 )

*-----------------------------------------------------------------------
*        set mesh area
*-----------------------------------------------------------------------

               az(ir) = pi * ( rm(ir+1)**2 - rm(ir)**2 )

               ar(ir,iz) = 2.0 * pi * rm(ir)
     &                     * ( zm(iz+1) - zm(iz) )

               

*-----------------------------------------------------------------------

      include 'samepage_include/samepage001.inc'

      allocate (ew(ne),aw(na),tw(nt))

      az_sum = 0.0d0
      do ir=1,nr
        az_sum = az_sum + az(ir)
      enddo

      allocate (ar_r(nz),ar_z(nr))

      do iz=1,nz
         do ir=1,nr
           ar_r(iz)= ar_r(iz) + ar(ir,iz)
           ar_z(ir)= ar_z(ir) + ar(ir,iz)
         enddo
      enddo

      yen  = char(92)
      igsh = 0

*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 2) = '[1/cm^2/(MeV/n)/source]         '
         hsunit( 5) = '[1/cm^2/(MeV/n)/sr/source]      '
         hsunit(12) = '[1/cm^2/(MeV/n)/nsec/source]    '
         hsunit(15) = '[1/cm^2/(MeV/n)/sr/nsec/source] '
      end if

*-----------------------------------------------------------------------
*        current or flux
*-----------------------------------------------------------------------

            if( itout(m) .eq. 1 ) then

               cname = '  Flux '
               dname = '  flux '

            else if( itout(m) .eq. 2 ) then

               cname = 'Current'
               dname = 'current'

            else if( itout(m) .eq. 3 ) then

               cname = ' F-Curr'
               dname = ' f-curr'

            else if( itout(m) .eq. 4 ) then

               cname = ' B-Curr'
               dname = ' b-curr'

            else if( itout(m) .eq. 5 ) then

               cname = ' O-Curr'
               dname = ' o-curr'

            else if( itout(m) .eq. 6 ) then

               cname = 'OF-Curr'
               dname = 'of-curr'

            else if( itout(m) .eq. 7 ) then

               cname = 'OB-Curr'
               dname = 'ob-curr'

            else if( itout(m) .eq. 8 ) then

               cname = ' A-Curr'
               dname = ' a-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 9 ) then

               cname = 'OA-Curr'
               dname = 'oa-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

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
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq.  2 .or. itunt(m) .eq.  5 .or.
     &               itunt(m) .eq. 12 .or. itunt(m) .eq. 15 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 13 .or. itunt(m) .eq. 16 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            end if

*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.d+0

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
*           itunt(m) = 11, 12, 13, 14 : /.../nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .gt. 10 ) then

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
*-----------------------------------------------------------------------
*        z-crossing
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1


            do 100 im = 1, nm
            do 100 iz = 1, nz + 1
            do 100 ir = 1, nr
            do 100 it = 1, nt
            do 100 ia = 1, na
            do 100 ie = 1, ne
            do 100 ip = 1, np

               if( tz(ip,ie,ia,it,izf(ir,iz),im,1) .gt. 0.0d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                         tz(ip,ie,ia,it,izf(ir,iz),im,1),
     &                         tz(ip,ie,ia,it,izf(ir,iz),im,2),
     &                         rtfac(m)/az(ir)/ew(ie)/aw(ia)/tw(it))

                  tz(ip,ie,ia,it,izf(ir,iz),im,1) = Xa
                  tz(ip,ie,ia,it,izf(ir,iz),im,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tz(ip,ie,ia,it,izf(ir,iz),im,1) .gt. cmax )
     &                          cmax = tz(ip,ie,ia,it,izf(ir,iz),im,1)

                  if( tz(ip,ie,ia,it,izf(ir,iz),im,1) .lt. cmin )
     &                          cmin = tz(ip,ie,ia,it,izf(ir,iz),im,1)

               else

                  isdz = 1
                  tz(ip,ie,ia,it,izf(ir,iz),im,2) = 0.0

               end if

! sumover
               call psufrz_sumover_tz_stdev_ntf(0,m,ntf,
     &               ip,ie,ia,it,ir,iz,im,
     &               rtfac(m),ew(ie),aw(ia),tw(it),az(ir),
     &               ew_sum,aw_sum,tw_sum,az_sum)


  100       continue

*-----------------------------------------------------------------------
*        r-crossing
*-----------------------------------------------------------------------


            do 200 im = 1, nm
            do 200 iz = 1, nz
            do 200 ir = 1, nr + 1
            do 200 it = 1, nt
            do 200 ia = 1, na
            do 200 ie = 1, ne
            do 200 ip = 1, np

               if( tr(ip,ie,ia,it,irf(ir,iz),im,1) .gt. 0.d0 ) then

                call calc_stdev(m,Xa,sigx,
     &                       tr(ip,ie,ia,it,irf(ir,iz),im,1),
     &                       tr(ip,ie,ia,it,irf(ir,iz),im,2),
     &                       rtfac(m)/ar(ir,iz)/ew(ie)/aw(ia)/tw(it))

                  tr(ip,ie,ia,it,irf(ir,iz),im,1) = Xa
                  tr(ip,ie,ia,it,irf(ir,iz),im,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

               else

                  isdz = 1
                  tr(ip,ie,ia,it,irf(ir,iz),im,2) = 0.0

               end if
! sumover
               call psufrz_sumover_stdev_ntf(0,m,ntf,
     &               ip,ie,ia,it,ir,iz,im,
     &               rtfac(m),ew(ie),aw(ia),tw(it),ar(ir,iz),
     &               ew_sum,aw_sum,tw_sum,ar_r(iz),ar_z(ir))


  200       continue

            if( nobch .gt. ist_bat ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do iax = 1, itfln(m)

!OBINATA(2012.7.12): output *.err
         noe = 1
         if ( any( itaxs(m,iax) .eq. (/ 7 /) )
     &       .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &     noe = 2

         do ioe = 1, noe

         if( itall .eq. 2 .and. nobch .lt. maxbch ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.7.12): output *.err
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

!OBINATA(2012.7.12): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tcrsech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        z-crossing
*        energy axis
*-----------------------------------------------------------------------

      include 'samepage_include/samepage002_pmeatr1z1.inc'
      include 'samepage_include/samepagechp_pmeatrz_npmax.inc'
      include 'samepage_include/samepageseti.inc'

            nrstepi_0 = nrstepi
            nzstepi_0 = nzstepi

               inum = 0

        if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 14 ) then ! T.Sato 2023/01/09


            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz + 1, nzstepi
            do iri = 1, nr, nrstepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ir = iri
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3,3x,
     &            "ir  =",i3/
     &            "#   zmesh = ",1p1e13.4,/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, ir,
     &                        zm(iz),
     &                        rm(ir), rm(ir+nrstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia ="i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it ="i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'(/"x: Energy [MeV/n]")')
                else
                  write(iot,'(/"x: Energy [MeV]")')
                end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

           end if

               if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 16 .or.
     &             itety(m) .eq.  3 .or. itety(m) .eq.  5 ) then

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

      include 'samepage_include/peatrzm_e_ana.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,a1)')
     &                        cha, inum, iz, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,",  it =",i3,a1)')
     &                        cha, inum, iz, ir, it, cha
                  end if

               else

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ia =",i3,a1)')
     &                        cha, inum, iz, ir, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ia =",i3,",  it =",i3,a1)')
     &                        cha, inum, iz, ir, ia, it, cha
                  end if

               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,",  it =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, it, itmnt(m,im), cha
                  end if

               else

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ia =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ia =",i3,",  it =",i3,
     &                        ",  mset =",i3,a1)')
     &                     cha, inum, iz, ir, ia, it, itmnt(m,im), cha
                  end if

               end if

            end if

          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz surf =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(7) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------


                  areasum = 0.0d0
                  do i=iri,iri+nrstepi-1
                    areasum = areasum + az(i)
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  z surface crossing"/
     &                        "  z     &=&",1pe13.4," [cm]"/
     &                        "  area  &=&",1pe13.4," [cm^2]"/
     &                        "  rmin  &=&",1pe13.4," [cm]"/
     &                        "  rmax  &=&",1pe13.4," [cm]")')
     &                        yen,
     &                        zm(iz), areasum, rm(ir), rm(ir+nrstepi)
               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itout(m) .gt. 7 ) then
                  write(iot,'(
     &                        "  amin  &=&",1pe13.4,/
     &                        "  amax  &=&",1pe13.4)')
     &                        ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

                  write(iot,'("e:")')
            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        r axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 6 ) then

            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz + 1, nzstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ie = iei
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                if(nzstepi .eq. 1) then
                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3)')
     &                        inum, iz
                else
                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3," - ",i3)')
     &                        inum, iz,iz+nzstepi-1
                end if

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie ="i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia ="i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it ="i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: r [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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

      include 'samepage_include/peatrzm_r_ana.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",   ie =",i3,a1)')
     &                     cha, inum, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",   ie =",i3,",   it =",i3,a1)')
     &                     cha, inum, iz, ie, it, cha
                  end if

               else if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,a1)')
     &                     cha, inum, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,",   it =",i3,a1)')
     &                     cha, inum, iz, it, cha
                  end if

               else if( itout(m) .eq. 8 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",   it =",i3,a1)')
     &                     cha, inum, iz, ie, ia, it, cha
                  end if

               else if( itout(m) .eq. 9 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ia =",i3,",   it =",i3,a1)')
     &                     cha, inum, iz, ia, it, cha
                  end if

               end if

*-----------------------------------------------------------------------
            else

               if( itout(m) .le. 4 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",   ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",   ie =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .eq. 8 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, ia, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .eq. 9 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ia =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ia, it, itmnt(m,im), cha
                  end if

               end if

            end if

          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz surf =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(4) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  z     &=&",1pe13.4," [cm]")')
     &                     yen, zm(iz)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')
            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iei = 1, ne, nestepi
            do iri = 1, nr, nrstepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ie = iei
               ir = iri
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,"ir  =",i3/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir,
     &                        rm(ir), rm(ir+nrstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie ="i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia ="i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it ="i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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

      include 'samepage_include/peatrzm_z_ana.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ir =",i3,a1)')
     &                     cha, inum, ie, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ir =",i3,",   it =",i3,a1)')
     &                     cha, inum, ie, ir, it, cha
                  end if

               else if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,a1)')
     &                     cha, inum, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, it, cha
                  end if

               else if( itout(m) .eq. 8 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie ="i3,
     &                     ",  ia ="i3,
     &                     ",  ir =",i3,a1)')
     &                     cha, inum, ie, ia, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie ="i3,
     &                     ",  ia ="i3,
     &                     ",  ir =",i3,",   it =",i3,a1)')
     &                     cha, inum, ie, ia, ir, it, cha
                  end if

               else if( itout(m) .eq. 9 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia ="i3,
     &                     ",   ir =",i3,a1)')
     &                     cha, inum, ia, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia ="i3,
     &                     ",   ir =",i3,",   it =",i3,a1)')
     &                     cha, inum, ia, ir, it, cha
                  end if

               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ir =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ir =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ir, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .eq. 8 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie ="i3,
     &                     ",  ia ="i3,
     &                     ",  ir =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie ="i3,
     &                     ",  ia ="i3,
     &                     ",  ir =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, ir, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .eq. 9 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia ="i3,
     &                     ",   ir =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia ="i3,
     &                     ",   ir =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, ir, it, itmnt(m,im), cha
                  end if

               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(5) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  areasum = 0.0d0
                  do i=iri,iri+nrstepi-1
                    areasum = areasum + az(i)
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]")')
     &                     yen, areasum,
     &                     rm(ir), rm(ir+nrstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
               end if
               if( itout(m) .ge. 8 ) then

                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
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

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then

            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz + 1, nzstepi
            do iri = 1, nr, nrstepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ir = iri
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

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3,3x,
     &            "ir  =",i3/
     &            "#   zmesh = ",1p1e13.4,/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, ir,
     &                        zm(iz),
     &                        rm(ir), rm(ir+nrstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 8 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

           end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_a_ana.inc'

*-----------------------------------------------------------------------

           if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,a1)')
     &                        cha, inum, iz, ir, ie, cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,",   it =",i3,a1)')
     &                        cha, inum, iz, ir, ie, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,a1)')
     &                        cha, inum, iz, ir, cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,",   it =",i3,a1)')
     &                        cha, inum, iz, ir, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, ie, itmnt(m,im), cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,",   it =",i3,
     &                        ",  mset =",i3,a1)')
     &                     cha, inum, iz, ir, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, itmnt(m,im), cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,",   it =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz surf =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(6) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  areasum = 0.0d0
                  do i=iri,iri+nrstepi-1
                    areasum = areasum + az(i)
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     zm(iz), areasum, rm(ir), rm(ir+nrstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .eq. 8 ) then
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
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

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz + 1, nzstepi
            do iri = 1, nr, nrstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ir = iri
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

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3,3x,
     &            "ir  =",i3/
     &            "#   zmesh = ",1p1e13.4,/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, ir,
     &                        zm(iz),
     &                        rm(ir), rm(ir+nrstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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

      include 'samepage_include/peatrzm_t_ana.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,a1)')
     &                        cha, inum, iz, ir, ie, cha
               else if( itout(m) .le. 7 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,a1)')
     &                        cha, inum, iz, ir, cha
               else if( itout(m) .eq. 8 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,
     &                        ",  ia =",i3,a1)')
     &                        cha, inum, iz, ir, ie, ia, cha
               else if( itout(m) .eq. 9 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ia =",i3,a1)')
     &                        cha, inum, iz, ir, ia, cha
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, ie, itmnt(m,im), cha
               else if( itout(m) .le. 7 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, itmnt(m,im), cha
               else if( itout(m) .eq. 8 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,
     &                        ",  ia =",i3,
     &                        ",  mset =",i3,a1)')
     &                     cha, inum, iz, ir, ie, ia, itmnt(m,im), cha
               else if( itout(m) .eq. 9 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ia =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, ia, itmnt(m,im), cha
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz surf =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(8) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                  areasum = 0.0d0
                  do i=iri,iri+nrstepi-1
                     areasum = areasum + az(i)
                  end do
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     zm(iz), areasum, rm(ir), rm(ir+nrstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4/
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

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

         end if

*-----------------------------------------------------------------------
*        r-crossing
*        energy axis
*-----------------------------------------------------------------------

        if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 14 ) then ! T.Sato 2023/01/09

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
          nr_1 = 1
          if(rm(1) .le. 0.0) then
            nr_1 =2
            if(iloopmode .eq. 6) then
              nrstepi = nrstepi - 1
            end if
          end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iri = nr_1, nr + 1, nrstepi
            do izi = 1, nz, nzstepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               iz = izi
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ir surf =",i3,3x,
     &            "iz  =",i3/
     &            "#   rmesh = ",1p1e13.4,/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir, iz,
     &                        rm(ir),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#    ia =",i3/
     &            "#   ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#    it =",i3/
     &            "#     t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Energy [MeV]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

           end if

               if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 16 .or.
     &             itety(m) .eq.  3 .or. itety(m) .eq.  5 ) then

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

      include 'samepage_include/peatrzm_e_r_cross_ana.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ir, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,a1)')
     &                     cha, inum, ir, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, it, itmnt(m,im), cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir surf =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(7) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                 areasum = 0.0d0
                 do i=izi,izi+nzstepi-1
                   areasum = areasum + ar(iri,i)
                 end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  r surface crossing"/
     &                     "  r     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     rm(ir), areasum, zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
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

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        r-crossing
*        r axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 6 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iei = 1, ne, nestepi
            do izi = 1, nz, nzstepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ie = iei
               iz = izi
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iz  =",i3/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz,
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: r [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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


      include 'samepage_include/peatrzm_r_r_ana.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then
            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ie, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ie, iz, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, iz, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ie, ia, iz, cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ie, ia, iz, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ia, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ia, iz, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, iz, itmnt(m,im), cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, iz, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(4) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  r surface crossing"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        r-crossing
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
          nr_1 = 1
          if(rm(1) .le. 0.0) then
            nr_1 =2
            if(iloopmode .eq. 6) then
              nrstepi = nrstepi - 1
            end if
          end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iri = nr_1, nr + 1, nrstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               ie = iei
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ir surf =",i3/
     &            "#   rmesh = ",1p1e13.4)')
     &                        inum, ir, rm(ir)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                          ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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


      include 'samepage_include/peatrzm_z_r_cross_ana.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,a1)')
     &                     cha, inum, ir, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir surf =",i3,a1)')
     &                     cha, inum, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir surf =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,a1)')
     &                     cha, inum, ir, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ia =",i3,a1)')
     &                     cha, inum, ir, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ia =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir surf =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir surf =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ie, ia, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ia =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir surf =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(5) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  r surface crossing"/
     &                     "  r     &=&",1pe13.4," [cm]")')
     &                     yen, rm(ir)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        r-crossing
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
          nr_1 = 1
          if(rm(1) .le. 0.0) then
            nr_1 =2
            if(iloopmode .eq. 6) then
              nrstepi = nrstepi - 1
            end if
          end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iri = nr_1, nr + 1, nrstepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
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

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ir surf =",i3,3x,
     &            "iz  =",i3/
     &            "#   rmesh = ",1p1e13.4,/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir, iz,
     &                        rm(ir),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 8 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

           end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_a_r_cross_ana.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,a1)')
     &                     cha, inum, ir, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, ie, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ir, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir surf =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(6) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                 areasum = 0.0
                 do i=izi,izi+nzstepi-1
                   areasum = areasum + ar(iri,i)
                 end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  r surface crossing"/
     &                     "  r     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     rm(ir), areasum, zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .eq. 8 ) then
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
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

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        r-crossing
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
          nr_1 = 1
          if(rm(1) .le. 0.0) then
            nr_1 =2
            if(iloopmode .eq. 6) then
              nrstepi = nrstepi - 1
            end if
          end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iri = nr_1, nr + 1, nrstepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               im = imi
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

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ir surf =",i3,3x,
     &            "iz  =",i3/
     &            "#   rmesh = ",1p1e13.4,/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir, iz,
     &                        rm(ir),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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

      include 'samepage_include/peatrzm_t_r_cross_ana.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,a1)')
     &                     cha, inum, ir, iz, ie, cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ir, iz, cha
               else if( itout(m) .le. 8 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  ia = ",i3,a1)')
     &                     cha, inum, ir, iz, ie, ia, cha
               else if( itout(m) .eq. 9 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,a1)')
     &                     cha, inum, ir, iz, ia, cha
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, itmnt(m,im), cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
               else if( itout(m) .le. 8 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, ia, itmnt(m,im), cha
               else if( itout(m) .eq. 9 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, itmnt(m,im), cha
               end if

            end if

          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir surf =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(8) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  areasum = 0.0d0
                  do i=izi,izi+nzstepi-1
                    areasum = areasum + ar(iri,i)
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  r surface crossing"/
     &                     "  r     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     rm(ir), areasum, zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4/
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

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

            call prestart(m,iot) !OBINATA(2012.7.12)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

      end do
      end do


*-----------------------------------------------------------------------
      include 'samepage_include/samepage999.inc'

      deallocate (ar_r,ar_z)
      deallocate (ew,aw,tw)

      return
      end subroutine anatal_psufrz



************************************************************************
*                                                                      *
      subroutine anatal_psufxyz(m,np,nx,ny,nz,ne,na,nt,nm,xm,ym,zm,
     &                   eb,ab,tb,
     &                   tr,nfile,weightRate,igsh,idasa,manatally)
*                                                                      *
*       output xyz scoring mesh z surfacecrossing current & flux       *
*       created by T.Miura on 2021/09/14                               *
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

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      double precision resc2   ! W
      double precision resc3   ! N
      common /restart/ resc2(itlmax), resc3(itlmax)

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
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

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

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
      dimension   ab(na+1)
      dimension   tb(nt+1)
!      dimension   ew(ne)
!      dimension   aw(na)
!      dimension   tw(nt)
      real(8),allocatable :: ew(:),aw(:),tw(:)

      integer,allocatable :: ixyz(:)
      dimension   tr(np,ne,na,nt,nx*ny*(nz+1),nm,2*nfile)

      integer     nfile, manatally
!      dimension   weightRate(nfile), rdata(2,nfile)
      dimension   weightRate(nfile)
      real(8),allocatable :: rdata(:,:)

      integer irst,nrst, itmpdata
!!      dimension anatalrst(np,ne+1,na+1,nt+1,nx*ny*(nz+1),nm+1,2*nfile+3)
! sumover
!      dimension anatalrst(np,ne+1,na+1,nt+1,(nx+1)*(ny+1)*(nz+2),
!     &                    nm+1,2*nfile+3)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)
      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9  ! S.H. 2021.8.9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*21,cijaxs*21,cijaxs2*21
      integer iDaxis
!      real(8) vl(nx,ny,nz+1)
      real(8),allocatable ::  vl(:,:,:)

      real(8),allocatable :: ax_x(:),ax_y(:)

*-----------------------------------------------------------------------

      dimension tott(6,2)

      character hsunit(16)*32

      data hsunit( 1) / '[1/cm^2/source]                 '/
      data hsunit( 2) / '[1/cm^2/MeV/source]             '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]        '/
      data hsunit( 4) / '[1/cm^2/sr/source]              '/
      data hsunit( 5) / '[1/cm^2/MeV/sr/source]          '/
      data hsunit( 6) / '[1/cm^2/Lethargy/sr/source]     '/
      data hsunit(11) / '[1/cm^2/nsec/source]            '/
      data hsunit(12) / '[1/cm^2/MeV/nsec/source]        '/
      data hsunit(13) / '[1/cm^2/Lethargy/nsec/source]   '/
      data hsunit(14) / '[1/cm^2/sr/nsec/source]         '/
      data hsunit(15) / '[1/cm^2/MeV/sr/nsec/source]     '/
      data hsunit(16) / '[1/cm^2/Lethargy/sr/nsec/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

*-----------------------------------------------------------------------

      character dc2*4

      character rpa*1
      data rpa /'}'/

      character cname*7
      character dname*7
      character aname*3

      character yen*1
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2

*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny
      icf2(ix,iy,iz) = ix + ( iy - 1 ) * (nx + 1)
     &               + ( iz - 1 ) * (nx + 1) * (ny + 1)

*-----------------------------------------------------------------------
*        set mesh area
*-----------------------------------------------------------------------

               ax(ix,iy) = ( xm(ix+1) - xm(ix) )
     &                   * ( ym(iy+1) - ym(iy) )

*-----------------------------------------------------------------------

      nrst = 2*nfile+3 ! S.H. 2021.8.15
      allocate (anatalrst(np,ne+1,na+1,nt+1,(nx+1)*(ny+1)*(nz+2),
     &                    nm+1,nrst))

      allocate (ew(ne),aw(na),tw(nt))
      allocate (rdata(2,nfile))

      allocate (ax_x(ny),ax_y(nx))
      ax_x(:) = 0.0d0
      ax_y(:) = 0.0d0
      do iy=1,ny
        do ix=1,nx
          ax_x(iy) = ax_x(iy) + ax(ix,iy)
          ax_y(ix) = ax_y(ix) + ax(ix,iy)
        enddo
      enddo

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

      allocate (vl(nx,ny,nz+1))
      do ix = 1, nx
      do iy = 1, ny
      do iz = 1, nz+1
         vl(ix,iy,iz) = ax(ix,iy)
      end do
      end do
      end do

*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 2) = '[1/cm^2/(MeV/n)/source]         '
         hsunit( 5) = '[1/cm^2/(MeV/n)/sr/source]      '
         hsunit(12) = '[1/cm^2/(MeV/n)/nsec/source]    '
         hsunit(15) = '[1/cm^2/(MeV/n)/sr/nsec/source] '
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
            nag = na
            ntg = nt

         else

            npg = 1
            neg = 1
            nag = 1
            ntg = 1

         end if

*-----------------------------------------------------------------------
*        current or flux
*-----------------------------------------------------------------------

            if( itout(m) .eq. 1 ) then

               cname = '  Flux '
               dname = '  flux '

            else if( itout(m) .eq. 2 ) then

               cname = 'Current'
               dname = 'current'

            else if( itout(m) .eq. 3 ) then

               cname = ' F-Curr'
               dname = ' f-curr'

            else if( itout(m) .eq. 4 ) then

               cname = ' B-Curr'
               dname = ' b-curr'

            else if( itout(m) .eq. 5 ) then

               cname = ' O-Curr'
               dname = ' o-curr'

            else if( itout(m) .eq. 6 ) then

               cname = 'OF-Curr'
               dname = 'of-curr'

            else if( itout(m) .eq. 7 ) then

               cname = 'OB-Curr'
               dname = 'ob-curr'

            else if( itout(m) .eq. 8 ) then

               cname = ' A-Curr'
               dname = ' a-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 9 ) then

               cname = 'OA-Curr'
               dname = 'oa-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq.  2 .or. itunt(m) .eq.  5 .or.
     &               itunt(m) .eq. 12 .or. itunt(m) .eq. 15 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 13 .or. itunt(m) .eq. 16 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            end if

*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

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
*           itunt(m) = 11, 12, 13, 14 : /.../nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .gt. 10 ) then

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
*        c1 : nomalization for source
*-----------------------------------------------------------------------
*        z-crossing
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

            do 100 im = 1, nm
            do 100 iz = 1, nz + 1
            do 100 ix = 1, nx
            do 100 iy = 1, ny
            do 100 ie = 1, ne
            do 100 ia = 1, na
            do 100 ip = 1, np
            do 100 it = 1, nt
            do ntf = 1, nfile

               if( tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf))
     &              .gt. 0.d0 ) then

                call calc_stdev(m,Xa,sigx,
     &                 tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf)),
     &                 tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(2,ntf)),
     &                 rtfac(m)/ax(ix,iy)/ew(ie)/aw(ia)/tw(it))

                  tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf)) = Xa
                  tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf))
     &                 .gt. cmax )
     &                 cmax =tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf))

                  if( tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf))
     &                 .lt. cmin )
     &                 cmin =tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(2,ntf)) = 0.0

               end if

! sumover
               call psufxyz_sumover_stdev_ntf(0,m,ntf,
     &               ip,ie,ia,it,ix,iy,iz,im,
     &               rtfac(m),ew(ie),aw(ia),tw(it),ax(ix,iy),
     &               ew_sum,aw_sum,tw_sum,ax_x(iy),ax_y(ix))

            end do
  100       continue

!$OMP PARALLEL
!$OMP& private(ipomp,ix,iy,iz,ia,im,it,ie,ip,ntf,isdz,rdata,answer,rerr)
!$OMP& private(iax,ix_a,iy_a,iz_a,ia_a,it_a,ie_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
            do 101 im = 1, nm
            do 101 iz = 1, nz + 1
            do 101 ix = 1, nx
            do 101 iy = 1, ny
            do 101 ie = 1, ne
            do 101 ia = 1, na
            do 101 ip = 1, np
            do 101 it = 1, nt

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

               do ntf = 1, nfile
                  rdata(1,ntf)
     &                 = tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf))
                  rdata(2,ntf)
     &                 = tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),im,1) = answer
              anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),im,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf)
     &                 = tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(1,ntf))
                  rdata(2,ntf)
     &                 = tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),im,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),im,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),im,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),im,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),im,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),im,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),im,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),im,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),im,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),im,iat(ioe,ntf))
     &              = tr(ip,ie,ia,it,icf(ix,iy,iz),im,iat(ioe,ntf))
                end do
               end do

            end if

!sumover
            do iax = 1,itaxn(m)

              ie_a = ie
              ia_a = ia
              it_a = it
              ix_a = ix
              iy_a = iy
              iz_a = iz
              ido_ana = 0
              if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 14)
     &             .and. ie == 1) then
                ie_a = ne + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 3 .and. ix == 1) then
                ix_a = nx + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 4 .and. iy == 1) then
                iy_a = ny + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 5 .and. iz == 1) then
                iz_a = nz + 2
                ido_ana = 1
              elseif(itaxs(m,iax) == 9 .and. it == 1) then
                it_a = nt + 1
                ido_ana = 1
              elseif((itaxs(m,iax) == 8 .or. itaxs(m,iax) == 10)
     &             .and. ia == 1) then
                ia_a = na + 1
                ido_ana = 1
              endif

              if(ido_ana == 1) then

               call get_crsxyz_tr_sum_data(m,iax,nfile,
     &                 ip,ie,ia,it,ix,iy,iz,im,rdata)

              if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                call usranatal(nfile,rdata,answer,rerr)

                anatalrst(ip,ie_a,ia_a,it_a,icf2(ix_a,iy_a,iz_a),im,1)
     &             = answer
                anatalrst(ip,ie_a,ia_a,it_a,icf2(ix_a,iy_a,iz_a),im,2)
     &             = rerr

              else if( manatally .eq. 1 ) then

                call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                anatalrst(ip,ie_a,ia_a,it_a,icf2(ix_a,iy_a,iz_a),im,1) ! mean
     &              = fmval
                if ( fmval .gt. 0.0d0 ) then
                 anatalrst(ip,ie_a,ia_a,it_a,icf2(ix_a,iy_a,iz_a),im,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                else
                 anatalrst(ip,ie_a,ia_a,it_a,icf2(ix_a,iy_a,iz_a),im,2) ! total uncertainty
     &                 = 0.0d0
                 end if

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      anatalrst(ip,ie_a,ia_a,it_a,icf2(ix_a,iy_a,iz_a),
     &                          im,iat(ioe,ntf))
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

!OBINATA(2012.7.12): output *.err
         noe = 1
         if ( any( itaxs(m,iax) .eq. (/ 7 /) )
     &       .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &     noe = 2

         do 900 ioe = 1, noe

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or.
     &       ittwo(m) .eq. 4 .or. ittwo(m) .eq. 5 ) ) goto 900

         if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &       igsh .eq. 0 ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.7.2): output *.err
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

!OBINATA(2012.7.2): output *.err
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

               call tcrsech(iot,m,iax,1)

*-----------------------------------------------------------------------
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        z-crossing
*        energy, x, y, z, time, angle axes
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1
     &        .or. itaxs(m,iax) .eq. 3
     &        .or. itaxs(m,iax) .eq. 4
     &        .or. itaxs(m,iax) .eq. 5
     &        .or. itaxs(m,iax) .eq. 9
     &        .or. itaxs(m,iax) .eq. 8
     &        .or. itaxs(m,iax) .eq. 10
     &        .or. itaxs(m,iax) .eq. 14 ) then ! T.Sato 2023/01/09

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 1   ! energy
           case ( 3 )
            iDaxis = 3   ! x axis
           case ( 4 )
            iDaxis = 4   ! y axis
           case ( 5 )
            iDaxis = 5   ! z axis
           case ( 9 )
            iDaxis = 9   ! time
           case ( 8 )
            iDaxis = 12  ! cos of angle(p) axis
           case ( 10 )
            iDaxis = 13  ! the of angle(p) axis
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.15

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(itaty(m),iMeVperu,3,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      na,    nt,     nx, ny, nz+1,  nm, nrst,
     &         eb,ew,  ab,aw,  tb,tw, xm, ym, zm,1,
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


*-----------------------------------------------------------------------

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
               write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

        if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 14 ) then ! T.Sato 2023/01/09
               if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 16 .or.
     &             itety(m) .eq.  3 .or. itety(m) .eq.  5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 3 ) then
               if( itxty(m) .eq. 3 .or. itxty(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 4 ) then
               if( ityty(m) .eq. 3 .or. ityty(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 5 ) then
               if( itzty(m) .eq. 3 .or. itzty(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 9 ) then
               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then
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
c S.H. added IF statement below for epsout=2 (2016.7.25)
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
            else
               write(iot,'( "h: n",12x,"x",12x,
     &         1000(a1,i0,a9,"),hh0",a3
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &         ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if

               write(iot,'( "#  lower        upper  ",3x,
     &                    1000(a1,2x,a8,4x,"r.err "))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n n n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
            else
               write(iot,'( "h: n",12x,"x",12x,
     &         1000(a1,i0,a9,"),hh0",a3
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &         ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if

               write(iot,'( "#  lower        upper  ",3x,
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

           voll = 0.0d0

           njaxs = nijaxs + 1

           do ijaxs = 1, nijaxs

            if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
                write(iot,'(1p2e13.4,1000(1pe13.4,0p3f8.4))')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,4)
     &               ,ip=1,np) ! frtati 2021/10/05

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
            if( itaxs(m,iax) .eq.  1 ) then ! eng axis
                     if( itunt(m) .eq.  2 .or. itunt(m) .eq.  3 .or.
     &                   itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &                   itunt(m) .eq. 12 .or. itunt(m) .eq. 13 .or.
     &                   itunt(m) .eq. 15 .or. itunt(m) .eq. 16 ) then
                 tott(ip,1) = tott(ip,1) / voll
                     end if

            else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then ! cos or the axis
                     if( itunt(m) .eq.  4 .or. itunt(m) .eq.  5 .or.
     &                   itunt(m) .eq.  6 .or.
     &                   itunt(m) .eq. 14 .or. itunt(m) .eq. 15 .or.
     &                   itunt(m) .eq. 16 ) then
                 tott(ip,1) = tott(ip,1) / voll
                     end if

            else
                     tott(ip,1) = tott(ip,1) / voll
            end if
               end do

! sumover
            do ip = 1, np
             do irst=1,2
              tott(ip,irst) = anataldata(ip,ianataldata,njaxs,irst)
            enddo
            enddo

           if ( manatally .ne. 2 ) then
               write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4))')
     &                     (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
           end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

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
*        z-crossing
*        xy axis ( matrix)
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

          select case( itaxs(m,iax) )
           case ( 7 )
            iDaxis = 31   ! xy
          end select

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(itaty(m),iMeVperu,3,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      na,    nt,     nx, ny, nz+1,  nm, nrst,
     &         eb,ew,  ab,aw,  tb,tw, xm, ym, zm,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

               inum = 0

          do ianataldata = 1, nanataldata

               inum = inum + 1

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

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

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
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                 write(iot,'("y: ",a15)') cijaxs2

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

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
     &                    "  Flux       r.err")')

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

                  write(iot,'( "#   x = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7,/
     &                         "#   y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            axs2min + axs2del/2d0, axs2max - axs2del/2d0, axs2del,
     &            axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &          'ax2/ax1',( fgaxs(ijaxs) + axs1del/2d0, ijaxs=1,nijaxs )

               do iy = ny, 1, -1

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

               zval = fgaxs3(ij(5))
               none = 1
               iaxs = 1
               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                nijaxs+1,nijaxs2+1,none,fgaxs,fgaxs2,zval,ixyz(1),
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

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else

        if( itazl(m) .eq. 0 ) then

         write(iot,'("y: ",a7,1x,a32)') cname, hsunit(itunt(m))

        else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

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

         end if

*-----------------------------------------------------------------------


            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )

      deallocate (ew,aw,tw)
      deallocate (ax_x,ax_y)
      deallocate (vl,rdata)

      deallocate (anatalrst)
      
      return
      end subroutine anatal_psufxyz



************************************************************************
*                                                                      *
      subroutine anatal_psufxyz_rpp(m,np,nx,ny,nz,ne,na,nt,nm,
     &                                   xm,ym,zm,eb,ab,
     &                   tb,tr,nfile,weightRate,igsh,idasa,manatally)
*                                                                      *
*       output xyz scoring mesh current & flux for enclosure mode      *
*       created by T.Miura on 2021/09/14                               *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
      common /tall09/ itpan(itlmax), itpat(itlmax,6,2),
     &                jtpat(itlmax,6,6,2)
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
      common /tall28/ itgsh(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall73/ itenclo(itlmax), itangform(itlmax)

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

*-----------------------------------------------------------------------

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   ab(na+1)
      dimension   tb(nt+1)
!      dimension   ew(ne)
!      dimension   aw(na)
!      dimension   tw(nt)
      real(8),allocatable :: ew(:),aw(:),tw(:)

      integer,allocatable :: ixyz(:)
      dimension   tr(np,ne,na,nt,nx*ny*(nz+1),nm,2)

      integer     nfile, manatally
      dimension   weightRate(nfile)

*-----------------------------------------------------------------------

      dimension tott(6,2)

      character hsunit(16)*32

      data hsunit( 1) / '[1/cm^2/source]                 '/
      data hsunit( 2) / '[1/cm^2/MeV/source]             '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]        '/
      data hsunit( 4) / '[1/cm^2/sr/source]              '/
      data hsunit( 5) / '[1/cm^2/MeV/sr/source]          '/
      data hsunit( 6) / '[1/cm^2/Lethargy/sr/source]     '/
      data hsunit(11) / '[1/cm^2/nsec/source]            '/
      data hsunit(12) / '[1/cm^2/MeV/nsec/source]        '/
      data hsunit(13) / '[1/cm^2/Lethargy/nsec/source]   '/
      data hsunit(14) / '[1/cm^2/sr/nsec/source]         '/
      data hsunit(15) / '[1/cm^2/MeV/sr/nsec/source]     '/
      data hsunit(16) / '[1/cm^2/Lethargy/sr/nsec/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character chl(6)*3
      data chl /'l  ','dr ','ub ','mg ','qrr','vbb'/
      character chm(6)*4
      data chm /'l3  ','d4r ','u5b ','m6g ','q7rr','v8bb'/

      character chp(6)*11  ! kitamura22/03/31
      character chq(6)*9   ! kitamura22/03/31

      character chb(6)*1
      data chb /6*' '/

*-----------------------------------------------------------------------

      character dc2*4

      character rpa*1
      data rpa /'}'/

      character cname*7
      character dname*7
      character aname*3

      character yen*1
      
      real(8),allocatable :: ax_x(:,:),ax_y(:,:),ax_z(:,:)

      include 'samepage_include/samepage000.inc'

*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

*-----------------------------------------------------------------------
*        set mesh area for rpp
*-----------------------------------------------------------------------

      ax(ix,iy,iz) = 2.d0 * ( xm(ix+1)-xm(ix) ) * ( ym(iy+1)-ym(iy) )
     &             + 2.d0 * ( ym(iy+1)-ym(iy) ) * ( zm(iz+1)-zm(iz) )
     &             + 2.d0 * ( zm(iz+1)-zm(iz) ) * ( xm(ix+1)-xm(ix) )

      include 'samepage_include/samepage001.inc'

*-----------------------------------------------------------------------

      allocate  (ew(ne),aw(na),tw(nt))

! sumover

      allocate (ax_x(ny,nz),ax_y(nx,nz),ax_z(nx,ny))
      ax_x(:,:) = 0.0d0
      ax_y(:,:) = 0.0d0
      ax_z(:,:) = 0.0d0
      do iz = 1,nz
        do iy = 1,ny
          do ix = 1,nx
            ax_x(iy,iz) = ax_x(iy,iz) + ax(ix,iy,iz)
            ax_y(ix,iz) = ax_y(ix,iz) + ax(ix,iy,iz)
            ax_z(ix,iy) = ax_z(ix,iy) + ax(ix,iy,iz)
          enddo
        enddo
      enddo

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 2) = '[1/cm^2/(MeV/n)/source]         '
         hsunit( 5) = '[1/cm^2/(MeV/n)/sr/source]      '
         hsunit(12) = '[1/cm^2/(MeV/n)/nsec/source]    '
         hsunit(15) = '[1/cm^2/(MeV/n)/sr/nsec/source] '
      end if

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            nag = na
            ntg = nt

         else

            npg = 1
            neg = 1
            nag = 1
            ntg = 1

         end if

*-----------------------------------------------------------------------
*        current or flux
*-----------------------------------------------------------------------

            if( itout(m) .eq. 1 ) then

               cname = '  Flux '
               dname = '  flux '

            else if( itout(m) .eq. 2 ) then

               cname = 'Current'
               dname = 'current'

            else if( itout(m) .eq. 3 ) then

               cname = ' F-Curr'
               dname = ' f-curr'

            else if( itout(m) .eq. 4 ) then

               cname = ' B-Curr'
               dname = ' b-curr'

            else if( itout(m) .eq. 5 ) then

               cname = ' O-Curr'
               dname = ' o-curr'

            else if( itout(m) .eq. 6 ) then

               cname = 'OF-Curr'
               dname = 'of-curr'

            else if( itout(m) .eq. 7 ) then

               cname = 'OB-Curr'
               dname = 'ob-curr'

            else if( itout(m) .eq. 8 ) then

               cname = ' A-Curr'
               dname = ' a-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 9 ) then

               cname = 'OA-Curr'
               dname = 'oa-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq.  2 .or. itunt(m) .eq.  5 .or.
     &               itunt(m) .eq. 12 .or. itunt(m) .eq. 15 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 13 .or. itunt(m) .eq. 16 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            end if

*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

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
*           itunt(m) = 11, 12, 13, 14 : /.../nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .gt. 10 ) then

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
*        c1 : nomalization for source
*-----------------------------------------------------------------------
*        rpp-crossing
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

            do 100 im = 1, nm
            do 100 iz = 1, nz
            do 100 ix = 1, nx
            do 100 iy = 1, ny
            do 100 ie = 1, ne
            do 100 ia = 1, na
            do 100 ip = 1, np
            do 100 it = 1, nt

               if( tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) .gt. 0.d0 ) then

                call calc_stdev(m,Xa,sigx,
     &                          tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &                          tr(ip,ie,ia,it,icf(ix,iy,iz),im,2),
     &                       rtfac(m)/ax(ix,iy,iz)/ew(ie)/aw(ia)/tw(it))

                  tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) = Xa
                  tr(ip,ie,ia,it,icf(ix,iy,iz),im,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) .gt. cmax )
     &                     cmax = tr(ip,ie,ia,it,icf(ix,iy,iz),im,1)

                  if( tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) .lt. cmin )
     &                     cmin = tr(ip,ie,ia,it,icf(ix,iy,iz),im,1)

               else

                  isdz = 1
                  tr(ip,ie,ia,it,icf(ix,iy,iz),im,2) = 0.0

               end if

! sumover
               call psufxyz_sumover_rpp_stdev_ntf(0,m,ntf,
     &               ip,ie,ia,it,ia,ix,iy,iz,im,
     &               rtfac(m),ew(ie),aw(ia),tw(it),ax(ix,iy,iz),
     &               ew_sum,aw_sum,tw_sum,
     &               ax_x(iy,iz),ax_y(ix,iz),ax_z(ix,iy))

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

!OBINATA(2012.7.12): output *.err
         noe = 1
         if ( any( itaxs(m,iax) .eq. (/ 7 /) )
     &       .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &     noe = 2

         do 900 ioe = 1, noe

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or.
     &       ittwo(m) .eq. 4 .or. ittwo(m) .eq. 5 ) ) goto 900

         if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &       igsh .eq. 0 ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.7.2): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else

!OBINATA(2012.7.2): output *.err
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

               call tcrsech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        rpp-crossing
*        energy axis
*-----------------------------------------------------------------------
      include 'samepage_include/samepage002_pmeatxyz1.inc'
      include 'samepage_include/samepagechp_pmeatxyz_npmax.inc'
      include 'samepage_include/samepageseti.inc'

        if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 14 ) then ! T.Sato 2023/01/09

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if

            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz, nzstepi
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ix = ixi
               iy = iyi
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------

                  write(iot,'("#   no. =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iy  =",i3/
     &            "iz  =",i3,3x,
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iy, iz,
     &                        xm(ix), xm(ix+nxstepi),
     &                        ym(iy), ym(iy+nystepi),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'(/"x: Energy [MeV/n]")')
                else
                  write(iot,'(/"x: Energy [MeV]")')
                end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

           end if

               if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 16 .or.
     &             itety(m) .eq.  3 .or. itety(m) .eq.  5 ) then

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

      include 'samepage_include/peatxyzm_e_non_ana.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ix, iy, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iy, iz, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, it, itmnt(m,im), cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ia, it,
     &                     itmnt(m,im), cha
                  end if
               end if

            end if

*-----------------------------------------------------------------------
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(8) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  areasum = 0.0d0
                  do iii=1,izi,izi+nzstepi-1
                   do ii=1,iyi,iyi+nystepi-1
                    do i=1,ixi,ixi+nxstepi-1
                        areasum = areasum + ax(i,ii,iii)
                    end do
                   end do
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, areasum, xm(ix), xm(ix+nxstepi),
     &                 ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
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
            end do

*-----------------------------------------------------------------------
*        rpp-crossing
*        x axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 3 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if

            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz, nzstepi
            do iyi = 1, ny, nystepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               iy = iyi
               ie = iei
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iy  =",i3,3x
     &            "iz  =",i3/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iy, iz,
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: x [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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

      include 'samepage_include/peatxyzm_x_xyz_ana.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, iy, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, iz, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3, a1)')
     &                     cha, inum, iy, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, iz, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iy, iz, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, iz, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iy, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, ie, ia, it,
     &                     itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(4) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &            yen, ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')
            end do

            end do
            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        rpp-crossing
*        y axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 4 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz, nzstepi
            do ixi = 1, nx, nxstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ix = ixi
               ie = iei
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iz  =",i3/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iz,
     &                        xm(ix), xm(ix+nxstepi),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: y [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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


      include 'samepage_include/peatxyzm_y_xyz_ana.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, ix, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iz, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3, a1)')
     &                     cha, inum, ix, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iz, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ix, iz, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iz, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ix, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, ie, ia, it,
     &                     itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(5) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &               yen, xm(ix), xm(ix+nxstepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')
            end do

            end do
            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        rpp-crossing
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ie = iei
               ia = iai
               ix = ixi
               iy = iyi
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ix  =",i3,3x
     &            "iy  =",i3,/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iy,
     &                  xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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

      include 'samepage_include/peatxyzm_z_xyz_ana.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ie, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ix, iy, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iy, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ia, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ia, ix, iy, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ix, iy, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, it,
     &                     itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, ix, iy, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(6) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                  xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')
            end do

            end do
            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        rpp-crossing
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz, nzstepi
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ix = ixi
               iy = iyi
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

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
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

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 8 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

           end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatxyzm_a_ana.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ix, iy, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iy, iz, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, it,
     &                     itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(7) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, ax(ix,iy,iz), xm(ix), xm(ix+nxstepi),
     &                 ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .eq. 8 ) then
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
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
            end do

*-----------------------------------------------------------------------
*        rpp-crossing
*        t axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz, nzstepi
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ix = ixi
               iy = iyi
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

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ix =",i3,3x,
     &            "iy  =",i3,3x,
     &            "iz  =",i3/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iy, iz,
     &                        xm(ix), xm(ix+nxstepi),
     &                        ym(iy), ym(iy+nystepi),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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

      include 'samepage_include/peatxyzm_t_ana.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ix, iy, iz, cha
               else if( itout(m) .eq. 8 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  ia =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, ia, cha
               else if( itout(m) .eq. 9 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  ia =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ia, cha
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, itmnt(m,im), cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, itmnt(m,im), cha
               else if( itout(m) .eq. 8 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, ia,
     &                     itmnt(m,im), cha
               else if( itout(m) .eq. 9 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  ia =",i3,
     &                     ",  mset =",i3, a1)')
     &                     cha, inum, ix, iy, iz, ia, itmnt(m,im), cha
               end if

            end if

          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(9) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  areasum = 0.0d0
                  do iii = izi, izi+nzstepi-1
                    do ii = iyi, iyi+nystepi-1
                      do i = ixi, ixi+nxstepi-1
                        areasum = areasum + ax(i,ii,iii)
                      end do
                    end do
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, areasum, xm(ix), xm(ix+nxstepi),
     &                 ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
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
            end do

*-----------------------------------------------------------------------
*        rpp-crossing
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

            do im = 1, nm
            do iz = 1, nz
            do ip = 1, npg
            do ie = 1, neg
            do ia = 1, nag
            do it = 1, ntg

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------

                  write(iot,'("#   no. =",i3,3x,
     &            "iz =",i3/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, zm(iz), zm(iz+1)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+1)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, ia, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: x [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: y [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

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

               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), ix = 1, nx ),
     &             iy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# x          y        ",
     &                      "  flux       r.err")')

               do iy = 1, ny
               do ix = 1, nx

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               xm(ix)  + rtxdl(m)/2.0,
     &               ym(iy)  + rtydl(m)/2.0,
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'( "#   x = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7,/
     &                         "#   y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m),
     &            ym(1) + rtydl(m)/2.0, ym(ny) + rtydl(m)/2.0, rtydl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'y/x',( xm(ix) + rtxdl(m)/2.0, ix = 1, nx )

               do iy = ny, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            ym(iy) + rtydl(m)/2.0,
     &            ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), ix = 1, nx )

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

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else

        if( itazl(m) .eq. 0 ) then

         write(iot,'("y: ",a7,1x,a32)') cname, hsunit(itunt(m))

        else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

        end if

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  part  &=&  ",a8)')
     &                     yen, zm(iz), zm(iz+1), chq(ip)
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+1)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+1)
                  end if
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,it)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
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
            end do

*-----------------------------------------------------------------------
*        rpp-crossing
*        yz axis (matrix)
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

            do im = 1, nm
            do ix = 1, nx
            do ip = 1, npg
            do ie = 1, neg
            do ia = 1, nag
            do it = 1, ntg

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------

                  write(iot,'("#   no. =",i3,3x,
     &            "ix =",i3/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, xm(ix), xm(ix+1)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+1)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, ix, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,a1)')
     &                     cha, inum, ix, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ix, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ix, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, ie, ia, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: y [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

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

               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), iz = 1, nz ),
     &             iy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# y          z        ",
     &                      "  flux       r.err")')

               do iz = 1, nz
               do iy = 1, ny

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               ym(iy)  + rtydl(m)/2.0,
     &               zm(iz)  + rtzdl(m)/2.0,
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'( "#   y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7,/
     &                         "#   z = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            ym(1) + rtydl(m)/2.0, ym(ny) + rtydl(m)/2.0, rtydl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'y/z',( zm(iz) + rtzdl(m)/2.0, iz = 1, nz )

               do iy = ny, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            ym(iy) + rtydl(m)/2.0,
     &            ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), iz = 1, nz )

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

      else

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else

        if( itazl(m) .eq. 0 ) then

         write(iot,'("y: ",a7,1x,a32)') cname, hsunit(itunt(m))

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

                  write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  part  &=&  ",a8)')
     &                     yen, xm(ix), xm(ix+1), chq(ip)
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+1)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+1)
                  end if
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,it)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
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
            end do

*-----------------------------------------------------------------------
*        rpp-crossing
*        xz axis ( matrix)
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

            do im = 1, nm
            do iy = 1, ny
            do ip = 1, npg
            do ie = 1, neg
            do ia = 1, nag
            do it = 1, ntg

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------

                  write(iot,'("#   no. =",i3,3x,
     &            "iy =",i3/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iy, ym(iy), ym(iy+1)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+1)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, iy, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iy, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iy, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, ie, ia, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: x [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

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

               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), iz = 1, nz ),
     &             ix = nx, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# x          z        ",
     &                      "  flux       r.err")')

               do iz = 1, nz
               do ix = 1, nx

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               xm(ix)  + rtxdl(m)/2.0,
     &               zm(iz)  + rtzdl(m)/2.0,
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'( "#   y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7,/
     &                         "#   z = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'x/z',( zm(iz) + rtzdl(m)/2.0, iz = 1, nz )

               do ix = nx, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            xm(ix) + rtxdl(m)/2.0,
     &            ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), iz = 1, nz )

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
               iaxs = 2
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

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else

        if( itazl(m) .eq. 0 ) then

         write(iot,'("y: ",a7,1x,a32)') cname, hsunit(itunt(m))

        else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

        end if

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  part  &=&  ",a8)')
     &                     yen, ym(iy), ym(iy+1), chq(ip)
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+1)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+1)
                  end if
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,it)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
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
            end do

*-----------------------------------------------------------------------

         end if          ! itaxis

*-----------------------------------------------------------------------

            call prestart(m,iot) !OBINATA(2012.7.12)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      include 'samepage_include/samepage999.inc'

      deallocate (ax_x,ax_y,ax_z)
      deallocate (ew,aw,tw)


      return
      end subroutine anatal_psufxyz_rpp



************************************************************************
*                                                                      *
      subroutine anatal_psufrz_rcc(m,np,nr,nz,ne,na,nt,nm,
     &                                  rm,zm,eb,ab,tb,
     &                      tr,tz,nfile,weightRate,idasa,manatally)
*                                                                      *
*       output r-z scoring mesh current & flux for enclosure mode      *
*       created by T.Miura on 2021/09/14                               *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
      common /tall09/ itpan(itlmax), itpat(itlmax,6,2),
     &                jtpat(itlmax,6,6,2)
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
      common /tall35/ iteps(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall73/ itenclo(itlmax), itangform(itlmax)

      common /istcut/ ist_cut, ist_bat
*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

      common /tall67/ itstd(itlmax), rtstd(itlmax)

*-----------------------------------------------------------------------

      dimension   rm(nr+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   ab(na+1)
      dimension   tb(nt+1)
!      dimension   ew(ne)
!      dimension   aw(na)
!      dimension   tw(nt)
      real(8),allocatable :: ew(:),aw(:),tw(:)

      dimension   tr(np,ne,na,nt,(nr+1)*nz,nm,2)
      dimension   tz(np,ne,na,nt,nr*(nz+1),nm,2)

      real(8),allocatable :: ar_r(:),ar_z(:)

      integer     nfile, manatally
      dimension   weightRate(nfile)

*-----------------------------------------------------------------------

      dimension tott(6,2)

      character hsunit(16)*32

      data hsunit( 1) / '[1/cm^2/source]                 '/
      data hsunit( 2) / '[1/cm^2/MeV/source]             '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]        '/
      data hsunit( 4) / '[1/cm^2/sr/source]              '/
      data hsunit( 5) / '[1/cm^2/MeV/sr/source]          '/
      data hsunit( 6) / '[1/cm^2/Lethargy/sr/source]     '/
      data hsunit(11) / '[1/cm^2/nsec/source]            '/
      data hsunit(12) / '[1/cm^2/MeV/nsec/source]        '/
      data hsunit(13) / '[1/cm^2/Lethargy/nsec/source]   '/
      data hsunit(14) / '[1/cm^2/sr/nsec/source]         '/
      data hsunit(15) / '[1/cm^2/MeV/sr/nsec/source]     '/
      data hsunit(16) / '[1/cm^2/Lethargy/sr/nsec/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character chl(6)*3
      data chl /'l  ','dr ','ub ','mg ','qrr','vbb'/
      character chm(6)*4
      data chm /'l3  ','d4r ','u5b ','m6g ','q7rr','v8bb'/

      character chp(6)*11  ! kitamura22/03/31
      character chq(6)*9   ! kitamura22/03/31

      character chb(6)*1
      data chb /6*' '/

*-----------------------------------------------------------------------

      character rpa*1
      data rpa /'}'/

      character cname*7
      character dname*7
      character aname*3

      character yen*1

      character dc2*4

      include 'samepage_include/samepage000.inc'

*-----------------------------------------------------------------------

         irf(ir,iz) = ir + ( iz - 1 ) * ( nr + 1 )

*-----------------------------------------------------------------------
*        set mesh area ... two circules and two side walls
*-----------------------------------------------------------------------

         ar(ir,iz) = 2.d0 * pi * ( rm(ir+1)**2 - rm(ir)**2 )
     &             + 2.d0 * pi * ( rm(ir+1) + rm(ir) )
     &               * ( zm(iz+1) - zm(iz) )

*-----------------------------------------------------------------------
      include 'samepage_include/samepage001.inc'

      allocate (ew(ne),aw(na),tw(nt))
! sumover
       allocate (ar_r(nz),ar_z(nr))
       ar_r(:) = 0.0d0
       ar_z(:) = 0.0d0
       do iz=1,nz
          do ir=1,nr
            ar_r(iz) = ar_r(iz) + ar(ir,iz)
            ar_z(ir) = ar_z(ir) + ar(ir,iz)
          enddo
       enddo

      yen  = char(92)
      igsh = 0

*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 2) = '[1/cm^2/(MeV/n)/source]         '
         hsunit( 5) = '[1/cm^2/(MeV/n)/sr/source]      '
         hsunit(12) = '[1/cm^2/(MeV/n)/nsec/source]    '
         hsunit(15) = '[1/cm^2/(MeV/n)/sr/nsec/source] '
      end if

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            nag = na
            ntg = nt

         else

            npg = 1
            neg = 1
            nag = 1
            ntg = 1

         end if

*-----------------------------------------------------------------------
*        current or flux
*-----------------------------------------------------------------------

            if( itout(m) .eq. 1 ) then

               cname = '  Flux '
               dname = '  flux '

            else if( itout(m) .eq. 2 ) then

               cname = 'Current'
               dname = 'current'

            else if( itout(m) .eq. 3 ) then

               cname = ' F-Curr'
               dname = ' f-curr'

            else if( itout(m) .eq. 4 ) then

               cname = ' B-Curr'
               dname = ' b-curr'

            else if( itout(m) .eq. 5 ) then

               cname = ' O-Curr'
               dname = ' o-curr'

            else if( itout(m) .eq. 6 ) then

               cname = 'OF-Curr'
               dname = 'of-curr'

            else if( itout(m) .eq. 7 ) then

               cname = 'OB-Curr'
               dname = 'ob-curr'

            else if( itout(m) .eq. 8 ) then

               cname = ' A-Curr'
               dname = ' a-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 9 ) then

               cname = 'OA-Curr'
               dname = 'oa-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

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
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq.  2 .or. itunt(m) .eq.  5 .or.
     &               itunt(m) .eq. 12 .or. itunt(m) .eq. 15 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 13 .or. itunt(m) .eq. 16 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            end if

*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

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
*           itunt(m) = 11, 12, 13, 14 : /.../nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .gt. 10 ) then

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
*-----------------------------------------------------------------------
*        rcc-crossing
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

            do 200 im = 1, nm
            do 200 iz = 1, nz
            do 200 ir = 1, nr
            do 200 it = 1, nt
            do 200 ia = 1, na
            do 200 ie = 1, ne
            do 200 ip = 1, np

               if( tr(ip,ie,ia,it,irf(ir,iz),im,1) .gt. 0.d0 ) then

                call calc_stdev(m,Xa,sigx,
     &                       tr(ip,ie,ia,it,irf(ir,iz),im,1),
     &                       tr(ip,ie,ia,it,irf(ir,iz),im,2),
     &                       rtfac(m)/ar(ir,iz)/ew(ie)/aw(ia)/tw(it))

                  tr(ip,ie,ia,it,irf(ir,iz),im,1) = Xa
                  tr(ip,ie,ia,it,irf(ir,iz),im,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,ia,it,irf(ir,iz),im,1) .gt. cmax )
     &                     cmax = tr(ip,ie,ia,it,irf(ir,iz),im,1)

                  if( tr(ip,ie,ia,it,irf(ir,iz),im,1) .lt. cmin )
     &                     cmin = tr(ip,ie,ia,it,irf(ir,iz),im,1)

               else

                  isdz = 1
                  tr(ip,ie,ia,it,irf(ir,iz),im,2) = 0.0

               end if

! sumover
               call psufrz_sumover_stdev_ntf(0,m,ntf,
     &               ip,ie,ia,it,ir,iz,im,
     &               rtfac(m),ew(ie),aw(ia),tw(it),ar(ir,iz),
     &               ew_sum,aw_sum,tw_sum,
     &               ar_r(iz),ar_z(ir))

  200       continue

            if( nobch .gt. ist_bat ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do iax = 1, itfln(m)

!OBINATA(2012.7.12): output *.err
         noe = 1
         if ( any( itaxs(m,iax) .eq. (/ 7 /) )
     &       .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &     noe = 2

         do ioe = 1, noe

         if( itall .eq. 2 .and. nobch .lt. maxbch ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.7.12): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else

!OBINATA(2012.7.12): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tcrsech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        rcc-crossing
*        energy axis
*-----------------------------------------------------------------------
      include 'samepage_include/samepage002_pmeatrz.inc'
      include 'samepage_include/samepagechp_pmeatrz_npmax.inc'
      include 'samepage_include/samepageseti.inc'

               inum = 0

        if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 14 ) then ! T.Sato 2023/01/09

            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iri = 1, nr, nrstepi
            do izi = 1, nz, nzstepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               iz = izi
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------

                  write(iot,'("#   no. =",i3,3x,
     &            "ir  =",i3,3x,
     &            "iz  =",i3/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir, iz,
     &                        rm(ir), rm(ir+nrstepi),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#    ia =",i3/
     &            "#   ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#    it =",i3/
     &            "#     t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Energy [MeV]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

           end if

               if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 16 .or.
     &             itety(m) .eq.  3 .or. itety(m) .eq.  5 ) then

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

      include 'samepage_include/peatrzm_e_rcc_ana.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ir, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,a1)')
     &                     cha, inum, ir, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, it, itmnt(m,im), cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(7) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


*-----------------------------------------------------------------------
                  areasum = 0.0d0
                  do ii=izi,izi+nzstepi-1
                    do i=iri,iri+nrstepi-1
                      areasum = areasum + ar(i,ii)
                    end do
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rcc surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &         areasum, rm(ir), rm(ir+nrstepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
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
*        rcc-crossing
*        r axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 6 ) then

            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iei = 1, ne, nestepi
            do izi = 1, nz, nzstepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ie = iei
               iz = izi
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iz  =",i3/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz,
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: r [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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

      include 'samepage_include/peatrzm_r_rcc_ana.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ie, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ie, iz, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, iz, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ie, ia, iz, cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ie, ia, iz, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ia, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ia, iz, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, iz, itmnt(m,im), cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, iz, it, itmnt(m,im), cha
                  end if
               end if

            end if

          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(4) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rcc surface crossing"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')
            end do

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        rcc-crossing
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iri = 1, nr, nrstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               ie = iei
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ir  =",i3/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir,
     &                        rm(ir), rm(ir+nrstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                          ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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

      include 'samepage_include/peatrzm_z_rcc_ana.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ie =",i3,a1)')
     &                     cha, inum, ir, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ie =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir =",i3,a1)')
     &                     cha, inum, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,a1)')
     &                     cha, inum, ir, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ia =",i3,a1)')
     &                     cha, inum, ir, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ia =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, ia, it, cha
                  end if
               end if
            else

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, itmnt(m,im), cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
               else if( itout(m) .le. 8 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, ia, itmnt(m,im), cha
               else if( itout(m) .eq. 9 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, itmnt(m,im), cha
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(5) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rcc surface crossing"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]")')
     &                     yen, rm(ir), rm(ir+nrstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')
            end do

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        rcc-crossing
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then

            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iri = 1, nr, nrstepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
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

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
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

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 8 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

           end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_a_rcc_ana.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,a1)')
     &                     cha, inum, ir, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, ie, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ir, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(6) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  areasum = 0.0d0
                  do ii=izi,izi+nzstepi-1
                    do i=iri,iri+nrstepi-1
                      areasum = areasum + ar(i,ii)
                    end do
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rcc surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &         areasum, rm(ir), rm(ir+nrstepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .eq. 8 ) then
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
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
*        r-crossing
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iri = 1, nr, nrstepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               im = imi
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

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
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

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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

      include 'samepage_include/peatrzm_t_rcc_ana.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,a1)')
     &                     cha, inum, ir, iz, ie, cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ir, iz, cha
               else if( itout(m) .le. 8 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  ia = ",i3,a1)')
     &                     cha, inum, ir, iz, ie, ia, cha
               else if( itout(m) .eq. 9 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,a1)')
     &                     cha, inum, ir, iz, ia, cha
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, itmnt(m,im), cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
               else if( itout(m) .le. 8 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, ia, itmnt(m,im), cha
               else if( itout(m) .eq. 9 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, itmnt(m,im), cha
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(8) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if


*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  areasum = 0.0d0
                  do ii=izi,izi+nzstepi-1
                    do i=iri,iri+nrstepi-1
                      areasum = areasum + ar(i,ii)
                    end do
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rcc surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &           areasum, rm(ir), rm(ir+nrstepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4/
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
*        rcc-crossing
*        rz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 13 ) then

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

            do im = 1, nm
            do ip = 1, npg
            do ie = 1, neg
            do ia = 1, nag
            do it = 1, ntg

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+1)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,a1)')
     &                     cha, inum, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ", it =",i3,a1)')
     &                     cha, inum, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if

*-----------------------------------------------------------------------

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: r [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

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

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

                  write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            rm(nr) + rtrdl(m)/2.0, rm(1) + rtrdl(m)/2.0, rtrdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,irf(ir,iz),im,ioe),
     &              iz = 1, nz ), ir = nr, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# r          z        ",
     &                      "  flux       r.err")')

               do iz = 1, nz
               do ir = 1, nr

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               rm(ir)  + rtrdl(m)/2.0,
     &               zm(iz)  + rtzdl(m)/2.0,
     &               tr(ip,ie,ia,it,irf(ir,iz),im,1),
     &               tr(ip,ie,ia,it,irf(ir,iz),im,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'( "#   r = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7,/
     &                         "#   z = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            rm(1) + rtrdl(m)/2.0, rm(nr) + rtrdl(m)/2.0, rtrdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'r/z',( zm(iz) + rtzdl(m)/2.0, iz = 1, nz )

               do ir = nr, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            rm(ir) + rtrdl(m)/2.0,
     &            ( tr(ip,ie,ia,it,irf(ir,iz),im,ioe), iz = 1, nz )

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

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else

        if( itazl(m) .eq. 0 ) then

         write(iot,'("y: ",a7,1x,a32)') cname, hsunit(itunt(m))

        else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

        end if

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                  write(iot,'(/"wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  part  &=&  ",a8)')
     &                     yen, chq(ip)
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+1)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+1)
                  end if
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,it)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
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
            end do

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

            call prestart(m,iot) !OBINATA(2012.7.12)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

      end do
      end do

*-----------------------------------------------------------------------
      include 'samepage_include/samepage999.inc'

      deallocate (ar_r,ar_z)
      deallocate (ew,aw,tw)

      return
      end subroutine anatal_psufrz_rcc


!***********************************************************************
!                                                                      *
! sumover subroutine                                                   *
!                                                                      *
!***********************************************************************

!***********************************************************************
!                                                                      *
      subroutine get_crsreg_tr_sum_data(m,iax,nfile,
     &           ip,ie,ia,it,ir,im,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,ia,it,ir,im
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

        call get_crsreg_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &     ittnm_sum(m,iax),
     &     itrgn_sum(m,iax),itmst_sum(m,iax),
     &     ip,ie,ia,it,ir,im,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_crsreg_tr_sum_data_sub(tr_sum,
     &          np,ne,na,nt,nr,nm,ip,ie,ia,it,ir,im,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,na,nt,nr,nm,ip,ie,ia,it,ir,im
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,ne,na,nt,nr,nm,2)
 
      do i=1,2
        rdata(i)= tr_sum(ip,ie,ia,it,ir,im,i)
      enddo

      end

!***********************************************************************
!                                                                      *
      subroutine get_crsrz_tr_sum_data(m,iax,nfile,
     &           ip,ie,ia,it,ir,iz,im,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,ia,it,ir,iz,im
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

        call get_crsrz_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &     ittnm_sum(m,iax),
     &     itrnm_sum(m,iax),itznm_sum(m,iax),
     &     itmst_sum(m,iax),
     &     ip,ie,ia,it,ir,iz,im,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_crsrz_tr_sum_data_sub(tr_sum,
     &          np,ne,na,nt,nr,nz,nm,ip,ie,ia,it,ir,iz,im,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,na,nt,nr,nz,nm,ip,ie,ia,it,ir,iz,im
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,ne,na,nt,(nr+1)*nz,nm,2)

      icf(ir,iz) = ir + ( iz - 1 ) * (nr +1)

      do i=1,2
        rdata(i)= tr_sum(ip,ie,ia,it,icf(ir,iz),im,i)
      enddo

      end

!***********************************************************************
!                                                                      *
      subroutine get_crsrz_tz_sum_data(m,iax,nfile,
     &           ip,ie,ia,it,ir,iz,im,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,ia,it,ir,iz,im
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
!$        iskip = italsize0_2_sum(m,iax)
!$         p_sum => tr00_sum(ianatalm_sum(m,iax,ntf)+iskip:)
!$       else
           iskip = italsize_2_sum(m,iax)
           p_sum => tr0_sum(ianatalm_sum(m,iax,ntf)+iskip:)
C for nonshared_tally option
!$       end if

        call get_crsrz_tz_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &     ittnm_sum(m,iax),
     &     itrnm_sum(m,iax),itznm_sum(m,iax),
     &     itmst_sum(m,iax),
     &     ip,ie,ia,it,ir,iz,im,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_crsrz_tz_sum_data_sub(tz_sum,
     &          np,ne,na,nt,nr,nz,nm,ip,ie,ia,it,ir,iz,im,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,na,nt,nr,nz,nm,ip,ie,ia,it,ir,iz,im
      real(8) :: rdata(2)

      real(8) :: tz_sum(np,ne,na,nt,nr*(nz+1),nm,2)

      izf(ir,iz) = iz + ( ir - 1 ) * (nz +1)

      do i=1,2
        rdata(i)= tz_sum(ip,ie,ia,it,izf(ir,iz),im,i)
      enddo

      end

!***********************************************************************
!                                                                      *
      subroutine get_crsxyz_tr_sum_data(m,iax,nfile,
     &           ip,ie,ia,it,ix,iy,iz,im,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,ia,it,ix,iy,iz,im
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

        call get_crsxyz_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &     ittnm_sum(m,iax),
     &     itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &     itmst_sum(m,iax),
     &     ip,ie,ia,it,ix,iy,iz,im,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_crsxyz_tr_sum_data_sub(tr_sum,
     &          np,ne,na,nt,nx,ny,nz,nm,ip,ie,ia,it,ix,iy,iz,im,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,na,nt,nx,ny,nz,nm,ip,ie,ia,it,ix,iy,iz,im
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,ne,na,nt,nx*ny*(nz+1),nm,2)

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

      do i=1,2
        rdata(i)= tr_sum(ip,ie,ia,it,icf(ix,iy,iz),im,i)
      enddo

      end


!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_tr_sub(m, ntf,
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

      integer :: m, ntf

      double precision :: resc2,resc3

      real(8),pointer :: p_sum(:)
      real(8),pointer :: tran_sum(:)

      integer :: ln_sum,iax

      do iax =1,6

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(p_sum,m,iax)
!$          ln_sum = italsize0_2_sum(m,iax)/2
!$       else
            call GET_TR_HEAD_POINTER_SUM(p_sum,m,iax)
            ln_sum = italsize_2_sum(m,iax)/2
C for nonshared_tally option
!$       end if

        tran_sum => trANATAL_SUM(ianatalm_sum(m,iax,ntf):)

        if(ln_sum > 0 ) then 
          call anatal_calc_anova_2trAN(m,p_sum,tran_sum,
     &       ln_sum, resc2,resc3)

        endif
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_tz_sub(m, ntf,
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

      integer :: m, ntf

      double precision :: resc2,resc3

      real(8),pointer :: p_sum(:)
      real(8),pointer :: tran_sum(:)

      integer :: ln_sum,iax,iskip

      do iax =1,6

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TZ_HEAD_POINTER0_SUM(p_sum,m,iax)
!$          ln_sum = (italsize0_sum(m,iax) - italsize0_2_sum(m,iax))/2
!$          iskip = italsize0_2_sum(m,iax)
!$       else
            call GET_TZ_HEAD_POINTER_SUM(p_sum,m,iax)
            ln_sum = (italsize_sum(m,iax) - italsize_2_sum(m,iax))/2
            iskip = italsize_2_sum(m,iax)
C for nonshared_tally option
!$       end if

        tran_sum => trANATAL_SUM(ianatalm_sum(m,iax,ntf)+iskip:)

        if(ln_sum > 0 ) then 
          call anatal_calc_anova_2trAN(m,p_sum,tran_sum,
     &       ln_sum, resc2,resc3)

        endif
      enddo

      return
      end

