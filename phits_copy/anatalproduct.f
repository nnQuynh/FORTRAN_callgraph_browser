!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_productreg(m,ntf,
     &                             np,na,nt,ne,nr,
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
      integer           na
      integer           nt
      integer           ne
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,nt,na,nr,2)
      integer           nfile
      double precision  tranatal(np,ne,nt,na,nr,2*nfile)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir
      integer           ie
      integer           it
      integer           ia
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
        do it = 1, nt
         do ia = 1, na
          do ip = 1, np

             if( trRES(ip,ie,it,ia,ir,1) .gt. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(ip,ie,it,ia,ir,1),
     &                        trRES(ip,ie,it,ia,ir,2),
     &                        1.0d+0)
              trRES(ip,ie,it,ia,ir,1) = Xa
              trRES(ip,ie,it,ia,ir,2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(ip,ie,it,ia,ir,iat(1,ntf)) =
     &           trRES(ip,ie,it,ia,ir,1)

              ! (sig_xj)**2
              tranatal(ip,ie,it,ia,ir,iat(2,ntf)) =
     &        (  trRES(ip,ie,it,ia,ir,2)
     &         * trRES(ip,ie,it,ia,ir,1) )**2

              ! sig_x
              tranatal(ip,ie,it,ia,ir,iat(2,ntf)) = dsqrt(
     &        tranatal(ip,ie,it,ia,ir,iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(ip,ie,it,ia,ir,iat(2,ntf)) =
     &       (tranatal(ip,ie,it,ia,ir,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(ip,ie,it,ia,ir,iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(ip,ie,it,ia,ir,iat(1,ntf)) =
     &        tranatal(ip,ie,it,ia,ir,iat(1,ntf)) * resc2(m)

            if( tranatal(ip,ie,it,ia,ir,iat(1,ntf)) > cmax )
     &            cmax = tranatal(ip,ie,it,ia,ir,iat(1,ntf))
            if( tranatal(ip,ie,it,ia,ir,iat(1,ntf)) < cmin )
     &            cmin = tranatal(ip,ie,it,ia,ir,iat(1,ntf))

          end do     ! ip loop end
         end do      ! ie loop end
        end do       ! it loop end
       end do        ! ir loop end
      end do         ! im loop end

! sumover
      call anatal_calc_anova_sub(m,ntf,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return
      end subroutine anatal_calc_anova_productreg



!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_productrz(m,ntf,
     &                             np,na,nt,ne,nr,nz,
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
      integer           na
      integer           nt
      integer           ne
      integer           nr
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,nt,na,nr,nz,2)
      integer           nfile
      double precision  tranatal(np,ne,nt,na,nr,nz,2*nfile)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           iz
      integer           ir
      integer           ie
      integer           it
      integer           ia
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
         do it = 1, nt
          do ia = 1, na
           do ip = 1, np

             if( trRES(ip,ie,it,ia,ir,iz,1) .gt. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(ip,ie,it,ia,ir,iz,1),
     &                        trRES(ip,ie,it,ia,ir,iz,2),
     &                        1.0d+0)
              trRES(ip,ie,it,ia,ir,iz,1) = Xa
              trRES(ip,ie,it,ia,ir,iz,2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(ip,ie,it,ia,ir,iz,iat(1,ntf)) =
     &           trRES(ip,ie,it,ia,ir,iz,1)

              ! (sig_xj)**2
              tranatal(ip,ie,it,ia,ir,iz,iat(2,ntf)) =
     &        (  trRES(ip,ie,it,ia,ir,iz,2)
     &         * trRES(ip,ie,it,ia,ir,iz,1) )**2

              ! sig_x
              tranatal(ip,ie,it,ia,ir,iz,iat(2,ntf)) = dsqrt(
     &        tranatal(ip,ie,it,ia,ir,iz,iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(ip,ie,it,ia,ir,iz,iat(2,ntf)) =
     &       (tranatal(ip,ie,it,ia,ir,iz,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(ip,ie,it,ia,ir,iz,iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(ip,ie,it,ia,ir,iz,iat(1,ntf)) =
     &        tranatal(ip,ie,it,ia,ir,iz,iat(1,ntf)) * resc2(m)

            if( tranatal(ip,ie,it,ia,ir,iz,iat(1,ntf)) > cmax )
     &            cmax = tranatal(ip,ie,it,ia,ir,iz,iat(1,ntf))
            if( tranatal(ip,ie,it,ia,ir,iz,iat(1,ntf)) < cmin )
     &            cmin = tranatal(ip,ie,it,ia,ir,iz,iat(1,ntf))

           end do    ! ip loop end
          end do     ! ia loop end
         end do      ! it loop end
        end do       ! ie loop end
       end do        ! ir loop end
      end do         ! iz loop end

! sumover
      call anatal_calc_anova_sub(m,ntf,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return
      end subroutine anatal_calc_anova_productrz



!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_productxyz(m,ntf,
     &                             np,na,nt,ne,nx,ny,nz,
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
      integer           na
      integer           nt
      integer           ne
      integer           nx
      integer           ny
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,nt,na,nx*ny*nz,2)
      integer           nfile
      double precision  tranatal(np,ne,nt,na,nx*ny*nz,2*nfile)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           iz
      integer           iy
      integer           ix
      integer           ie
      integer           it
      integer           ia
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
          do it = 1, nt
           do ia = 1, na
            do ip = 1, np

             if( trRES(ip,ie,it,ia,icf(ix,iy,iz),1) .gt. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(ip,ie,it,ia,icf(ix,iy,iz),1),
     &                        trRES(ip,ie,it,ia,icf(ix,iy,iz),2),
     &                        1.0d+0)
              trRES(ip,ie,it,ia,icf(ix,iy,iz),1) = Xa
              trRES(ip,ie,it,ia,icf(ix,iy,iz),2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf)) =
     &           trRES(ip,ie,it,ia,icf(ix,iy,iz),1)

              ! (sig_xj)**2
              tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(2,ntf)) =
     &        (  trRES(ip,ie,it,ia,icf(ix,iy,iz),2)
     &         * trRES(ip,ie,it,ia,icf(ix,iy,iz),1) )**2

              ! sig_x
              tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(2,ntf)) = dsqrt(
     &        tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(2,ntf)) =
     &       (tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf)) =
     &        tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf)) * resc2(m)

            if( tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf)) > cmax )
     &            cmax = tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf))
            if( tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf)) < cmin )
     &            cmin = tranatal(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf))

            end do   ! ip loop end
           end do    ! ia loop end
          end do     ! it loop end
         end do      ! ie loop end
        end do       ! ix loop end
       end do        ! iy loop end
      end do         ! iz loop end

! sumover
      call anatal_calc_anova_sub(m,ntf,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return
      end subroutine anatal_calc_anova_productxyz



************************************************************************
*                                                                      *
      subroutine anatal_ppdctreg(m,np,nr,mr,ne,nt,na,kr,eb,tb,ab,tr,
     &                    nfile,weightRate,
     &                    nvl,ivl,rvl,
     &                    nx,ny,nz,xm,ym,zm,igsh,idasa,manatally)
*                                                                      *
*       output the product tally in region mesh                        *
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

      parameter ( pi  = 3.141592653589793d0 )

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

*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   ab(na+1)
!      dimension   vl(nr)
!      dimension   lr(nr)
!      dimension   ew(ne)
!      dimension   tw(nt)
!      dimension   aw(na)
      dimension   tr(np,ne,nt,na,nr,2*nfile)
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
!      dimension   anatalrst(np,ne+1,na+1,nt+1,nr*1*1,1,2*nfile+3)
! sumover
!      dimension   anatalrst(np,ne+1,na+1,nt+1,(nr+1)*1*1,1,2*nfile+3)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)

      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9  ! S.H. 2021.8.9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*21,cijaxs*21,cijaxs2*21
      integer iDaxis
!      dimension   dlr(nr+1) !  integer lr(nr) ---> double dlr(nr)
      real(8),allocatable :: vl(:),ew(:),tw(:),aw(:),rdata(:,:),
     &                       val(:),dlr(:)
      integer,allocatable :: lr(:)

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
!-----------------------------------------------------------------------
      character yen*1

      dimension dt_one(1)
      data dt_one /1.0d0/

      integer           iat, iad, iaf
      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------
      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )
      allocate (vl(nr),ew(ne),tw(nt),aw(na),rdata(2,nfile),
     &          val(nr),dlr(nr+1))
      allocate (lr(nr))

      nrst = 2*nfile+3 ! S.H. 2021.8.15
      allocate (anatalrst(np,ne+1,na+1,nt+1,(nr+1)*1*1,1,nrst))


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
*        set mesh volume and name
*-----------------------------------------------------------------------

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)
               vl_sum = sum(vl(:))

ccse 2021.08 add (use anatal_rearrange sub.)
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
            do 100 it = 1, nt
            do 100 ia = 1, na
            do 100 ip = 1, np
            do ntf = 1, nfile

               if( tr(ip,ie,it,ia,ir,iat(1,ntf)) .gt. 0.d0 ) then

                  cc = rtfac(m) / ew(ie) / tw(it) / aw(ia)

                  if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                itunt(m) .eq. 34 .or. itunt(m) .eq. 36 )
     &                cc = cc / vl(ir)

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,ia,ir,iat(1,ntf)),
     &                            tr(ip,ie,it,ia,ir,iat(2,ntf)),
     &                            cc)

                  tr(ip,ie,it,ia,ir,iat(1,ntf)) = Xa
                  tr(ip,ie,it,ia,ir,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,it,ia,ir,iat(1,ntf)) .gt. cmax )
     &                             cmax = tr(ip,ie,it,ia,ir,iat(1,ntf))

                  if( tr(ip,ie,it,ia,ir,iat(1,ntf)) .lt. cmin )
     &                             cmin = tr(ip,ie,it,ia,ir,iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,it,ia,ir,iat(2,ntf)) = 0.0

               end if

! sumover
               if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then
                 call ppdctreg_sumover_stdev_ntf(0,m,ntf,
     &               ip,ie,it,ia,ir,
     &               rtfac(m),ew(ie),tw(it),aw(ia),vl(ir),
     &               ew_sum,tw_sum,aw_sum,vl_sum)
               else
                 call ppdctreg_sumover_stdev_ntf(0,m,ntf,
     &               ip,ie,it,ia,ir,
     &               rtfac(m),ew(ie),tw(it),aw(ia),1.0d0,
     &               ew_sum,tw_sum,aw_sum,1.0d0)
               endif

            end do
  100       continue

!$OMP PARALLEL
!$OMP& private(ipomp,ir,ie,it,ia,ip,ntf,isdz,rdata,answer,rerr)
!$OMP& private(iax,ir_a,ie_a,it_a,ia_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
            do ir = 1, nr
            do 101 ie = 1, ne
            do 101 it = 1, nt
            do 101 ia = 1, na
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,it,ia,ir,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,it,ia,ir,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(ip,ie,ia,it,ir,1,1) = answer
              anatalrst(ip,ie,ia,it,ir,1,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,it,ia,ir,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,it,ia,ir,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(ip,ie,ia,it,ir,1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,ia,it,ir,1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,ia,it,ir,1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,ia,it,ir,1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,ia,it,ir,1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,ia,it,ir,1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,ir,1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,ir,1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,ir,1,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(ip,ie,ia,it,ir,1,iat(ioe,ntf))
     &              = tr(ip,ie,it,ia,ir,iat(ioe,ntf))
                end do
               end do

            end if

!sumover
            do iax = 1,itaxn(m)

              ie_a = ie
              it_a = it
              ir_a = ir
              ia_a = ia
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
              elseif((itaxs(m,iax) == 12 .or. itaxs(m,iax) == 13)
     &                .and. ia == 1) then
                ia_a = na + 1
                ido_ana = 1
              endif

              if(ido_ana == 1) then

                 call get_pdctreg_tr_sum_data(m,iax,nfile,
     &                 ip,ie,it,ia,ir,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   anatalrst(ip,ie_a,ia_a,it_a,ir_a,1,1) = answer
                   anatalrst(ip,ie_a,ia_a,it_a,ir_a,1,2) = rerr

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  anatalrst(ip,ie_a,ia_a,it_a,ir_a,1,1) ! mean
     &              = fmval
                 if ( fmval .gt. 0.0d0 ) then
                  anatalrst(ip,ie_a,ia_a,it_a,ir_a,1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                 else
                   anatalrst(ip,ie_a,ia_a,it_a,ir_a,1,2) ! total uncertainty
     &                 = 0.0d0
                 end if

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      anatalrst(ip,ie_a,ia_a,it_a,ir_a,1,iat(ioe,ntf))
     &                = rdata(ioe,ntf)
                    end do
                  end do

                end if

              endif

            enddo
  101       continue
      enddo  ! ir loop T.Sato 2021/04/15
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
     &        .or. itaxs(m,iax) .eq. 11
     &        .or. itaxs(m,iax) .eq. 12
     &        .or. itaxs(m,iax) .eq. 13
     &        .or. itaxs(m,iax) .eq. 15 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 1   ! energy
           case ( 2 )
            iDaxis = 2   ! reg axis
           case ( 11 )
            iDaxis = 9   ! time
           case ( 12 )
            iDaxis = 12  ! cos of angle(p) axis
           case ( 13 )
            iDaxis = 13  ! the of angle(p) axis
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.15

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(itaty(m),iMeVperu,1,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      na,    nt,     nr,  0,  0,   0, nrst,
     &         eb,ew,  ab,aw,  tb,tw, dlr,  dt_one,  dt_one,1,
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
               write(iot,'(/"x: ",a21)') cijaxs
            else
               write(iot,'(/"x: ",200a1)') (itaxt(m)(i:i),i=1,itaxl(m))
            end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
            if( itayl(m) .eq. 0 ) then
               write(iot,'( "y: Number ",a29)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

        if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 15 ) then ! T.Sato 2023/01/09
              if( itety(m) .eq.  3 .or. itety(m) .eq.  5 .or.
     &             itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &             itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &             itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then
                 write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
              else
                 write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
              end if

           else if( itaxs(m,iax) .eq. 11 ) then
              if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then
                 write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
              else
                 write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
              end if

           else if( itaxs(m,iax) .eq. 12 .or.
     &             itaxs(m,iax) .eq. 13  .or.
     &             itaxs(m,iax) .eq. 2 ) then
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
     &                     1000(a10,"),",a4," n   "))')
     &                     ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             else if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: x",12x,"n",12x,
     &                    6(a10,"),h0",a3," n "))')
     &                     ( chp(i), chl(i) , i = 1, np )
             else
                write(iot,'( "h: n",12x,"x",12x,
     &                     1000(a10,"),hh0",a3," n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if
            else
             if( itaxs(m,iax) .eq. 2 ) then ! reg axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                     1000(a1,i0,a9,"),",a4
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             else if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: n",12x,"x",12x,
     &           6(a1,i1,a9,"),h0",a3
     &          ," ny2",i1," dy",i1,"=[y",i1,"*y2",i1,"] "))')
     &           ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &           ,i,i,i,i, i = 1, np )
             else
                write(iot,'( "h: n",12x,"x",12x,
     &                     1000(a1,i0,a9,"),hh0",a3
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

            if( itaxs(m,iax) .eq. 2 ) then ! reg axis
               write(iot,'( "#  num    reg     volume  ",
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
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                     1000(a10,"),",a4," n n n "))')
     &                     ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             else if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: x",12x,"n",12x,
     &                    6(a10,"),h0",a3," n n n "))')
     &                     ( chp(i), chl(i) , i = 1, np )
             else
                write(iot,'( "h: n",12x,"x",12x,
     &                     1000(a10,"),hh0",a3," n n n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if

            else
             if( itaxs(m,iax) .eq. 2 ) then ! reg axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                     1000(a1,i0,a9,"),",a4
     &                ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             else if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: n",12x,"x",12x,
     &           6(a1,i1,a9,"),h0",a3
     &          ," ny2",i1," dy",i1,"=[y",i1,"*y2",i1,"] n n "))')
     &           ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &           ,i,i,i,i, i = 1, np )
             else
                write(iot,'( "h: n",12x,"x",12x,
     &                     1000(a1,i0,a9,"),hh0",a3
     &                ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

            if( itaxs(m,iax) .eq. 2 ) then ! reg axis
               write(iot,'( "#  num    reg     volume  ",
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
                write(iot,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0p3f8.4))')
     &               ijaxs, lr(ijaxs), fgaxs(ijaxs),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05
             else
                write(iot,'(1p2e13.4,1000(1pe13.4,0p3f8.4))')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05
             end if

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
             if( itaxs(m,iax) .eq. 2 ) then
                write(iot,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0p3f8.4))')
     &               ijaxs, lr(ijaxs), fgaxs(ijaxs),
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

           end do                 ! ijaxs = 1, nijaxs

           do ip = 1, np
            if( itaxs(m,iax) .eq.  1 ) then ! eng axis
              if( itunt(m) .eq.  3 .or. itunt(m) .eq.  4 .or.
     &             itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 14 .or.
     &             itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &             itunt(m) .eq. 23 .or. itunt(m) .eq. 24 .or.
     &             itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &             itunt(m) .eq. 33 .or. itunt(m) .eq. 34 .or.
     &             itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then
                 tott(ip,1) = tott(ip,1) / voll
              end if

            else if( itaxs(m,iax) .eq.  11 ) then ! time axis
               if( ( itunt(m) .gt. 10 .and. itunt(m) .lt. 20 ) .or.
     &              ( itunt(m) .gt. 30 .and. itunt(m) .lt. 40 ) ) then
                  tott(ip,1) = tott(ip,1) / voll
               end if

            else if( itaxs(m,iax) .eq. 12 .or.
     &            itaxs(m,iax) .eq. 13 ) then ! cos or the axis
              if( itunt(m) .gt.  20 ) then
                 tott(ip,1) = tott(ip,1) / voll
              end if

            else if( itaxs(m,iax) .eq. 2 ) then ! reg axis
               if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &              itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &              itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &              itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &              itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &              itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then
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
             if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 15 ) then ! T.Sato 2023/01/09
               write(iot,'(/"#   sum over ",   13x ,
     &               1000(1pe13.4,0pf8.4,16x))')
     &              (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
            else
               write(iot,'(/"#   sum over ",1pe13.4,
     6               1000(1pe13.4,0pf8.4,16x))')
     &              voll,(tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
            end if
           end if

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

            call prestart(m,iot) !OBINATA(2012.9.11)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      deallocate (anatalrst)
      deallocate (vl,ew,tw,aw,rdata,val,dlr)
      deallocate (lr)

      return
      end subroutine anatal_ppdctreg



************************************************************************
*                                                                      *
      subroutine anatal_ppdcttet(m,np,nr,mr,ne,nt,na,eb,tb,ab,tr,
     &                    nfile,weightRate,
     &                    nx,ny,nz,kr,xm,ym,zm,igsh,idasa,manatally)
*                                                                      *
*       output the product tally in tetra mesh                         *
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

      parameter ( pi  = 3.141592653589793d0 )

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

*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   ab(na+1)
!      dimension   ew(ne)
!      dimension   tw(nt)
!      dimension   aw(na)
      dimension   tr(np,ne,nt,na,nr,2*nfile)

      integer     nfile, manatally
!      dimension   weightRate(nfile), rdata(2,nfile)
      dimension   weightRate(nfile)

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)

      integer,allocatable :: ixyz(:)
      integer,allocatable :: lr(:)        !FURUTA20190204
      real(8),allocatable :: vl(:),val(:) !FURUTA20190204
c Dont know why but necessary to avoid segmentation fault

      integer irst,nrst, itmpdata
!      dimension   anatalrst(np,ne+1,na+1,nt+1,nr*1*1,1,2*nfile+3)
! sumover
!      dimension   anatalrst(np,ne+1,na+1,nt+1,(nr+1)*1*1,1,2*nfile+3)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)
      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9  ! S.H. 2021.8.9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*21,cijaxs*21,cijaxs2*21
      integer iDaxis
!      dimension   dlr(nr+1) !  integer lr(nr) ---> double dlr(nr)
      real(8),allocatable :: ew(:),tw(:),aw(:),rdata(:,:),dlr(:)

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

      dimension dt_one(1)
      data dt_one /1.0d0/

      integer           iat, iad, iaf
      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------
      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )
      allocate (ew(ne),tw(nt),aw(na),rdata(2,nfile),dlr(nr+1))

      nrst = 2*nfile+3 ! S.H. 2021.8.15
      allocate (anatalrst(np,ne+1,na+1,nt+1,(nr+1)*1*1,1,nrst))


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
               ew_sum =  eb(ne+1) - eb(1)

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
*        set mesh volume and name
*-----------------------------------------------------------------------

               allocate ( lr(nr),vl(nr),val(nr) ) !FURUTA20190204
               call ttetvl(mr,kr,nr,vl,lr)
               vl_sum = sum(vl(:))

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
            do 100 it = 1, nt
            do 100 ia = 1, na
            do 100 ip = 1, np
            do ntf = 1, nfile

               if( tr(ip,ie,it,ia,ir,iat(1,ntf)) .gt. 0.d0 ) then

                  cc = rtfac(m) / ew(ie) / tw(it) / aw(ia)

                  if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                itunt(m) .eq. 34 .or. itunt(m) .eq. 36 )
     &                cc = cc / vl(ir)

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,ia,ir,iat(1,ntf)),
     &                            tr(ip,ie,it,ia,ir,iat(2,ntf)),
     &                            cc)

                  tr(ip,ie,it,ia,ir,iat(1,ntf)) = Xa
                  tr(ip,ie,it,ia,ir,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,it,ia,ir,iat(1,ntf)) .gt. cmax )
     &                             cmax = tr(ip,ie,it,ia,ir,iat(1,ntf))

                  if( tr(ip,ie,it,ia,ir,iat(1,ntf)) .lt. cmin )
     &                             cmin = tr(ip,ie,it,ia,ir,iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,it,ia,ir,iat(2,ntf)) = 0.0

               end if

! sumover
               if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then
                 call ppdctreg_sumover_stdev_ntf(0,m,ntf,
     &               ip,ie,it,ia,ir,
     &               rtfac(m),ew(ie),tw(it),aw(ia),vl(ir),
     &               ew_sum,tw_sum,aw_sum,vl_sum)
               else
                 call ppdctreg_sumover_stdev_ntf(0,m,ntf,
     &               ip,ie,it,ia,ir,
     &               rtfac(m),ew(ie),tw(it),aw(ia),1.0d0,
     &               ew_sum,tw_sum,aw_sum,1.0d0)
               endif

            end do
  100       continue

!$OMP PARALLEL
!$OMP& private(ipomp,ir,ie,it,ia,ip,ntf,isdz,rdata,answer,rerr)
!$OMP& private(iax,ir_a,ie_a,it_a,ia_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
            do ir = 1, nr
            do 101 ie = 1, ne
            do 101 it = 1, nt
            do 101 ia = 1, na
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,it,ia,ir,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,it,ia,ir,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(ip,ie,ia,it,ir,1,1) = answer
              anatalrst(ip,ie,ia,it,ir,1,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,it,ia,ir,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,it,ia,ir,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(ip,ie,ia,it,ir,1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,ia,it,ir,1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,ia,it,ir,1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,ia,it,ir,1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,ia,it,ir,1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,ia,it,ir,1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,ir,1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,ir,1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,ir,1,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(ip,ie,ia,it,ir,1,iat(ioe,ntf))
     &              = tr(ip,ie,it,ia,ir,iat(ioe,ntf))
                end do
               end do

            end if

!sumover
            do iax = 1,itaxn(m)

              ie_a = ie
              it_a = it
              ir_a = ir
              ia_a = ia
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
              elseif((itaxs(m,iax) == 12 .or. itaxs(m,iax) == 13)
     &                .and. ia == 1) then
                ia_a = na + 1
                ido_ana = 1
              endif

              if(ido_ana == 1) then

                 call get_pdctreg_tr_sum_data(m,iax,nfile,
     &                 ip,ie,it,ia,ir,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   anatalrst(ip,ie_a,ia_a,it_a,ir_a,1,1) = answer
                   anatalrst(ip,ie_a,ia_a,it_a,ir_a,1,2) = rerr

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  anatalrst(ip,ie_a,ia_a,it_a,ir_a,1,1) ! mean
     &              = fmval
                 if ( fmval .gt. 0.0d0 ) then
                  anatalrst(ip,ie_a,ia_a,it_a,ir_a,1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                 else
                   anatalrst(ip,ie_a,ia_a,it_a,ir_a,1,2) ! total uncertainty
     &                 = 0.0d0
                 end if

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      anatalrst(ip,ie_a,ia_a,it_a,ir_a,1,iat(ioe,ntf))
     &                = rdata(ioe,ntf)
                    end do
                  end do

                end if

              endif

            enddo

  101       continue
      enddo  ! ir loop T.Sato 2021/04/15
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
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        energy, tet, time, angle axes
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1
     &        .or. itaxs(m,iax) .eq. 14
     &        .or. itaxs(m,iax) .eq. 11
     &        .or. itaxs(m,iax) .eq. 12
     &        .or. itaxs(m,iax) .eq. 13
     &        .or. itaxs(m,iax) .eq. 15 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 1   ! energy
           case ( 14 )
            iDaxis = 14  ! tet axis
           case ( 11 )
            iDaxis = 9   ! time
           case ( 12 )
            iDaxis = 12  ! cos of angle(p) axis
           case ( 13 )
            iDaxis = 13  ! the of angle(p) axis
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.15

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(itaty(m),iMeVperu,4,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      na,    nt,     nr,  0,  0,   0, nrst,
     &         eb,ew,  ab,aw,  tb,tw, dlr,  dt_one,  dt_one,1,
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
               write(iot,'(/"x: ",a21)') cijaxs
            else
               write(iot,'(/"x: ",200a1)') (itaxt(m)(i:i),i=1,itaxl(m))
            end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
            if( itayl(m) .eq. 0 ) then
               write(iot,'( "y: Number ",a29)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

        if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 15 ) then ! T.Sato 2023/01/09
              if( itety(m) .eq.  3 .or. itety(m) .eq.  5 .or.
     &             itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &             itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &             itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then
                 write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
              else
                 write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
              end if

           else if( itaxs(m,iax) .eq. 11 ) then
              if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then
                 write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
              else
                 write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
              end if

           else if( itaxs(m,iax) .eq. 12 .or.
     &             itaxs(m,iax) .eq. 13  .or.
     &             itaxs(m,iax) .eq. 14 ) then
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
     &                     1000(a10,"),",a4," n   "))')
     &                     ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             else if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: x",12x,"n",12x,
     &                    6(a10,"),h0",a3," n "))')
     &                     ( chp(i), chl(i) , i = 1, np )
             else
                write(iot,'( "h: n",12x,"x",12x,
     &                     1000(a10,"),hh0",a3," n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if
            else
             if( itaxs(m,iax) .eq. 14 ) then ! tet axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                     1000(a1,i0,a9,"),",a4
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             else if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: n",12x,"x",12x,
     &           6(a1,i1,a9,"),h0",a3
     &          ," ny2",i1," dy",i1,"=[y",i1,"*y2",i1,"] "))')
     &           ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &           ,i,i,i,i, i = 1, np )
             else
                write(iot,'( "h: n",12x,"x",12x,
     &                     1000(a1,i0,a9,"),hh0",a3
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

            if( itaxs(m,iax) .eq. 14 ) then ! tet axis
               write(iot,'( "#  num    reg     volume  ",
     &                    1000(a1,2x,a8,4x,"r.err "))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            else
               write(iot,'( "#  lower        upper  ",3x,
     &                    1000(a1,2x,a8,4x,"r.err "))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            end if

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
             if( itaxs(m,iax) .eq. 14 ) then ! tet axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                     1000(a10,"),",a4," n n n "))')
     &                     ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             else if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: x",12x,"n",12x,
     &                    6(a10,"),h0",a3," n n n "))')
     &                     ( chp(i), chl(i) , i = 1, np )
             else
                write(iot,'( "h: n",12x,"x",12x,
     &                     1000(a10,"),hh0",a3," n n n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if

            else
             if( itaxs(m,iax) .eq. 14 ) then ! tet axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                     1000(a1,i0,a9,"),",a4
     &                ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             else if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: n",12x,"x",12x,
     &           6(a1,i1,a9,"),h0",a3
     &          ," ny2",i1," dy",i1,"=[y",i1,"*y2",i1,"] n n "))')
     &           ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &           ,i,i,i,i, i = 1, np )
             else
                write(iot,'( "h: n",12x,"x",12x,
     &                     1000(a1,i0,a9,"),hh0",a3
     &                ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

            if( itaxs(m,iax) .eq. 14 ) then ! tet axis
               write(iot,'( "#  num    reg     volume  ",
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
             if( itaxs(m,iax) .eq. 14 ) then
                write(iot,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0p3f8.4))')
     &               ijaxs, lr(ijaxs), fgaxs(ijaxs),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05
             else
                write(iot,'(1p2e13.4,1000(1pe13.4,0p3f8.4))')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05
             end if

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
             if( itaxs(m,iax) .eq. 14 ) then
                write(iot,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0p3f8.4))')
     &               ijaxs, lr(ijaxs), fgaxs(ijaxs),
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

           end do                 ! ijaxs = 1, nijaxs

           do ip = 1, np
            if( itaxs(m,iax) .eq.  1 ) then ! eng axis
              if( itunt(m) .eq.  3 .or. itunt(m) .eq.  4 .or.
     &             itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 14 .or.
     &             itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &             itunt(m) .eq. 23 .or. itunt(m) .eq. 24 .or.
     &             itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &             itunt(m) .eq. 33 .or. itunt(m) .eq. 34 .or.
     &             itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then
                 tott(ip,1) = tott(ip,1) / voll
              end if

            else if( itaxs(m,iax) .eq.  11 ) then ! time axis
               if( ( itunt(m) .gt. 10 .and. itunt(m) .lt. 20 ) .or.
     &              ( itunt(m) .gt. 30 .and. itunt(m) .lt. 40 ) ) then
                  tott(ip,1) = tott(ip,1) / voll
               end if

            else if( itaxs(m,iax) .eq. 12 .or.
     &            itaxs(m,iax) .eq. 13 ) then ! cos or the axis
              if( itunt(m) .gt.  20 ) then
                 tott(ip,1) = tott(ip,1) / voll
              end if

            else if( itaxs(m,iax) .eq. 14 ) then ! tet axis
               if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &              itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &              itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &              itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &              itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &              itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then
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
            if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 15 ) then ! T.Sato 2023/01/09
               write(iot,'(/"#   sum over ",   13x ,
     6               1000(1pe13.4,0pf8.4,16x))')
     &              (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
            else
               write(iot,'(/"#   sum over ",1pe13.4,
     6               1000(1pe13.4,0pf8.4,16x))')
     &              voll,(tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
            end if
           end if

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

            call prestart(m,iot) !OBINATA(2012.9.11)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      deallocate( lr,vl,val ) !FURUTA20190204
      if(itfoam(m).eq.2) deallocate(xcm) !FURUTA20191028
      deallocate (anatalrst)
      deallocate (ew,tw,aw,rdata,dlr)

      return
      end subroutine anatal_ppdcttet



************************************************************************
*                                                                      *
      subroutine anatal_ppdctrz(m,np,nr,nz,ne,nt,na,rm,zm,eb,tb,ab,tr,
     &                   nfile,weightRate,
     &                   idasa,manatally)
*                                                                      *
*       output r-z scoring mesh product tally                          *
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

      parameter ( pi  = 3.141592653589793d0 )

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
      dimension   tb(nt+1)
      dimension   ab(na+1)
!      dimension   ew(ne)
!      dimension   tw(nt)
!      dimension   aw(na)
      dimension   tr(np,ne,nt,na,nr,nz,2*nfile)
      integer irst,nrst, itmpdata
!      dimension   anatalrst(np,ne+1,na+1,nt+1,nr*nz*1,1,2*nfile+3)
! sumover
!      dimension   anatalrst(np,ne+1,na+1,nt+1,(nr+1)*(nz+1)*1,
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
      real(8),allocatable :: ew(:),tw(:),aw(:),rdata(:,:)

*-----------------------------------------------------------------------

      dimension tott(6,2)

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


      real(8) vl(nr,nz,1)
      real(8),allocatable :: vl_r(:,:),vl_z(:,:)
      dimension dt_one(1)
      data dt_one /1.0d0/


*-----------------------------------------------------------------------

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

      yen  = char(92)
      igsh = 0

      allocate (ew(ne),tw(nt),aw(na),rdata(2,nfile))

      nrst = 2*nfile+3 ! S.H. 2021.8.15
      allocate (anatalrst(np,ne+1,na+1,nt+1,(nr+1)*(nz+1)*1,1,nrst))


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

*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------
! sumover
      allocate (vl_r(nz,1),vl_z(nr,1))
      if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &    itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &    itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &    itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &    itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &    itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then
        vl_r(:,:) = 0.0d0
        vl_z(:,:) = 0.0d0
        do ir = 1, nr
        do iz = 1, nz
                  vl(ir,iz,1) = pi * ( rm(ir+1)**2 - rm(ir)**2 )
     &                           * ( zm(iz+1) - zm(iz) )
          vl_r(iz,1) = vl_r(iz,1) + vl(ir,iz,1)
          vl_z(ir,1) = vl_z(ir,1) + vl(ir,iz,1)
        end do
        end do
      else
        vl(:,:,:) = 1.0d0
        vl_r(:,:) = 1.0d0
        vl_z(:,:) = 1.0d0
      endif


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
                  aw_sum = aw_sum +aw(i)

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

        do itmprz=1,nr*nz
            ir=(itmprz-1)/nz+1
            iz=itmprz-(ir-1)*nz
            do 100 ie = 1, ne
            do 100 it = 1, nt
            do 100 ia = 1, na
            do 100 ip = 1, np
            do ntf = 1, nfile

               if( tr(ip,ie,it,ia,ir,iz,iat(1,ntf)) .gt. 0.d0 ) then

                  if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then

                     vm= vl(ir,iz,1)

                  else

                     vm = 1.d+0

                  end if

                  cc = rtfac(m) / vm / ew(ie) / tw(it) / aw(ia)

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,ia,ir,iz,iat(1,ntf)),
     &                            tr(ip,ie,it,ia,ir,iz,iat(2,ntf)),
     &                            cc)

                  tr(ip,ie,it,ia,ir,iz,iat(1,ntf)) = Xa
                  tr(ip,ie,it,ia,ir,iz,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,it,ia,ir,iz,iat(1,ntf)) .gt. cmax )
     &                 cmax = tr(ip,ie,it,ia,ir,iz,iat(1,ntf))

                  if( tr(ip,ie,it,ia,ir,iz,iat(1,ntf)) .lt. cmin )
     &                 cmin = tr(ip,ie,it,ia,ir,iz,iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,it,ia,ir,iz,iat(2,ntf)) = 0.0

               end if

! sumover
               call ppdctrz_sumover_stdev_ntf(0,m,ntf,
     &               ip,ie,it,ia,ir,iz,
     &               rtfac(m),ew(ie),tw(it),aw(ia),vm,
     &               ew_sum,tw_sum,aw_sum,vl_r(iz,1),vl_z(ir,1))

            end do
  100       continue
      enddo  ! rz loop T.Sato 2021/04/15

!$OMP PARALLEL
!$OMP& private(ipomp,ir,iz,ia,it,ie,ip,ntf,itmprz,rdata,answer,rerr)
!$OMP& private(iax,ir_a,iz_a,ia_a,it_a,ie_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
        do itmprz=1,nr*nz
            ir=(itmprz-1)/nz+1
            iz=itmprz-(ir-1)*nz
            do 101 ie = 1, ne
            do 101 it = 1, nt
            do 101 ia = 1, na
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

               do ntf = 1, nfile
                  rdata(1,ntf)=tr(ip,ie,it,ia,ir,iz,iat(1,ntf))
                  rdata(2,ntf)=tr(ip,ie,it,ia,ir,iz,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(ip,ie,ia,it,icf(ir,iz),1,1) = answer
              anatalrst(ip,ie,ia,it,icf(ir,iz),1,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf)=tr(ip,ie,it,ia,ir,iz,iat(1,ntf))
                  rdata(2,ntf)=tr(ip,ie,it,ia,ir,iz,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(ip,ie,ia,it,icf(ir,iz),1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,ia,it,icf(ir,iz),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,ia,it,icf(ir,iz),1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,ia,it,icf(ir,iz),1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,ia,it,icf(ir,iz),1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,ia,it,icf(ir,iz),1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,icf(ir,iz),1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,icf(ir,iz),1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,icf(ir,iz),1,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(ip,ie,ia,it,icf(ir,iz),1,iat(ioe,ntf))
     &              = tr(ip,ie,it,ia,ir,iz,iat(ioe,ntf))
                end do
               end do

            end if

! sumover
            do iax = 1,itaxn(m)

              ie_a = ie
              it_a = it
              ia_a = ia
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
              elseif((itaxs(m,iax) == 12 .or. itaxs(m,iax) == 13)
     &               .and. ia == 1) then
                ia_a = na + 1
                ido_ana = 1
              endif

              if(ido_ana == 1) then

                 call get_pdctrz_tr_sum_data(m,iax,nfile,
     &                 ip,ie,it,ia,ir,iz,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   anatalrst(ip,ie_a,ia_a,it_a,icf(ir_a,iz_a),1,1)
     &             = answer
                   anatalrst(ip,ie_a,ia_a,it_a,icf(ir_a,iz_a),1,2)
     &             = rerr

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  anatalrst(ip,ie_a,ia_a,it_a,icf(ir_a,iz_a),1,1) ! mean
     &              = fmval
                 if ( fmval .gt. 0.0d0 ) then
                   anatalrst(ip,ie_a,ia_a,it_a,icf(ir_a,iz_a),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                 else
                   anatalrst(ip,ie_a,ia_a,it_a,icf(ir_a,iz_a),1,2) ! total uncertainty
     &                 = 0.0d0
                 end if

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      anatalrst(ip,ie_a,ia_a,it_a,1,icf(ir_a,iz_a),
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
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        angle(p), energy, r, z, time axes
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 12 .or. itaxs(m,iax) .eq. 13 .or.
     &       itaxs(m,iax) .eq.  1 .or.
     &       itaxs(m,iax) .eq.  6 .or.
     &       itaxs(m,iax) .eq.  5 .or.
     &       itaxs(m,iax) .eq. 11 ) then

          select case( itaxs(m,iax) )
           case ( 12 )
            iDaxis = 12  ! cos of angle(p) axis
           case ( 13 )
            iDaxis = 13  ! the of angle(p) axis
           case ( 1 )
            iDaxis = 1   ! energy
           case ( 6 )
            iDaxis = 6   ! r(-z) axis
           case ( 5 )
            iDaxis = 7   ! (r-)z axis
           case ( 11 )
            iDaxis = 9   ! time
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.15

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(itaty(m),iMeVperu,2,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      na,    nt,     nr, nz,  0,    0, nrst,
     &         eb,ew,    ab,aw,  tb,tw,  rm, zm, dt_one,1,
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
               write(iot,'( "y: Number ",a29)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

           if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 15 ) then ! T.Sato 2023/01/09
               if( itety(m) .eq.  3 .or. itety(m) .eq.  5 .or.
     &             itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &             itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &             itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 11 ) then ! time axis
               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 12
     &              .or. itaxs(m,iax) .eq. 13 ) then      ! angle axis
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

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
             if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: x",12x,"n",12x,
     &                    6(a10,"),h0",a3," n "))')
     &                     ( chp(i), chl(i) , i = 1, np )
             else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

             end if
            else
             if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: n",12x,"x",12x,
     &           6(a1,i1,a9,"),h0",a3
     &          ," ny2",i1," dy",i1,"=[y",i1,"*y2",i1,"] "))')
     &           ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &           ,i,i,i,i, i = 1, np )
             else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a1,i0,a9,"),hh0",a3
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                    ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                    ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

            write(iot,'( "#  lower        upper  ",3x,
     &                 1000(a1,2x,a8,4x,"r.err "))')
     &                 ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
             if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: x",12x,"n",12x,
     &                    6(a10,"),h0",a3," n n n "))')
     &                     ( chp(i), chl(i) , i = 1, np )
             else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n n n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if

            else
             if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: n",12x,"x",12x,
     &           6(a1,i1,a9,"),h0",a3
     &          ," ny2",i1," dy",i1,"=[y",i1,"*y2",i1,"] n n "))')
     &           ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &           ,i,i,i,i, i = 1, np )
             else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a1,i0,a9,"),hh0",a3
     &                ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &                    ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                    ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
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
                if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 15 ) then ! T.Sato 2023/01/09
                     if( itunt(m) .eq.  3 .or. itunt(m) .eq.  4 .or.
     &                   itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &                   itunt(m) .eq. 13 .or. itunt(m) .eq. 14 .or.
     &                   itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &                   itunt(m) .eq. 23 .or. itunt(m) .eq. 24 .or.
     &                   itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &                   itunt(m) .eq. 33 .or. itunt(m) .eq. 34 .or.
     &                   itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then
                        tott(ip,1) = tott(ip,1) / voll
                     end if

                  else if( itaxs(m,iax) .eq. 11 ) then ! time axis
                    if( ( itunt(m) .gt. 10 .and. itunt(m) .lt. 20 ) .or.
     &                    ( itunt(m) .gt. 30 .and. itunt(m) .lt. 40 ) )
     &                    then
                       tott(ip,1) = tott(ip,1) / voll
                    end if

                  else if( itaxs(m,iax) .eq. 12
     &                   .or. itaxs(m,iax) .eq. 13 ) then ! angle axis
                      if( itunt(m) .gt.  20 ) then
                       tott(ip,1) = tott(ip,1) / voll
                     end if

                  else if( itaxs(m,iax) .eq. 6
     &                    .or. itaxs(m,iax) .eq. 5 ) then ! r,z axis
                  if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then
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

           if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 11
     &     .or.itaxs(m,iax) .eq. 15 ) then ! T.Sato 2023/01/09
              write(iot,'(/"#   sum over ",   13x ,
     6              1000(1pe13.4,0pf8.4,16x))')
     &             (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
           else if( itaxs(m,iax) .eq. 6 .or. itaxs(m,iax) .eq. 5 ) then ! r,z axis
              if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &             itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &             itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &             itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &             itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &             itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then
               write(iot,'(/"#   sum over ",1pe13.4,
     6               1000(1pe13.4,0pf8.4,16x))')
     &                voll,(tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
            else
               write(iot,'(/"#   sum over ",   13x ,
     6               1000(1pe13.4,0pf8.4,16x))')
     &              (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
            end if
           else
              write(iot,'(/"#   sum over ",1pe13.4,
     6              1000(1pe13.4,0pf8.4,16x))')
     &             voll,(tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
           end if

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
          call anatal_rearrange_sum(itaty(m),iMeVperu,2,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      na,    nt,     nr, nz,  0,    0, nrst,
     &         eb,ew,    ab,aw,  tb,tw,  rm, zm, dt_one,1,
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

!OBINATA(2012.9.5): output *.err
               write(iot,'(1p10e11.3)')
     &       ( ( anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                ijaxs = 1, nijaxs ), ijaxs2 = nijaxs2, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

            write(iot,'(/"# axis2      axis1    ",
     &                      "  Number     r.err")')

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
*        gshow
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
         write(iot,'("y: Number ",a32)') hsunit(itunt(m))
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

            call prestart(m,iot) !OBINATA(2012.9.11)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

        end do

      end do

*-----------------------------------------------------------------------
      deallocate( vl_r,vl_z )
      deallocate (anatalrst)
      deallocate (ew,tw,aw,rdata)

      return
      end subroutine anatal_ppdctrz



************************************************************************
*                                                                      *
      subroutine anatal_ppdctxyz(m,np,nl,lt,
     &                    nx,ny,nz,ne,nt,na,xm,ym,zm,eb,tb,ab,tr,
     &                    nfile,weightRate,
     &                    igsh,idasa,manatally)
*                                                                      *
*       output xyz scoring mesh product tally                          *
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

      parameter ( pi  = 3.141592653589793d0 )

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
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   ab(na+1)
!      dimension   ew(ne)
!      dimension   tw(nt)
!      dimension   aw(na)
      real(8),allocatable :: ew(:),tw(:),aw(:)

      dimension   tr(np,ne,nt,na,nx*ny*nz,2*nfile)
      integer irst,nrst, itmpdata
!      dimension   anatalrst(np,ne+1,na+1,nt+1,nx*ny*nz,1,2*nfile+3)
! sumover
!      dimension   anatalrst(np,ne+1,na+1,nt+1,(nx+1)*(ny+1)*(nz+1),
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
      real(8),allocatable :: rdata(:,:)

      integer,allocatable :: ixyz(:)

*-----------------------------------------------------------------------

      dimension tott(6,2)

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


      real(8),allocatable :: vl(:,:,:)
! sumover
      real(8),allocatable :: vl_x(:,:),vl_y(:,:),vl_z(:,:)

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
      integer           iat, iad, iaf
      iat(iad,iaf) = iad + (iaf-1) * 2

*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx
     &              + ( iz - 1 ) * nx * ny
      icf2(ix,iy,iz) = ix + ( iy - 1 ) * (nx+1)
     &              + ( iz - 1 ) * (nx+1) * (ny+1)

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------
      allocate (ew(ne),tw(nt),aw(na))
      allocate (rdata(2,nfile))

      allocate( vl(nx,ny,nz) )
      allocate( vl_x(ny,nz),vl_y(nx,nz),vl_z(nx,ny) )
      if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &    itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &    itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &    itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &    itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &    itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then

         vl_x(:,:) = 0.0d0
         vl_y(:,:) = 0.0d0
         vl_z(:,:) = 0.0d0

         do ix = 1, nx
         do iy = 1, ny
         do iz = 1, nz
           vl(ix,iy,iz) = vls(nl,lt,itmcn(m),itvm(m),itmtr(m,4),
     &                      xm(ix),xm(ix+1),
     &                      ym(iy),ym(iy+1),
     &                      zm(iz),zm(iz+1))
           vl_x(iy,iz) = vl_x(iy,iz) + vl(ix,iy,iz)
           vl_y(ix,iz) = vl_y(ix,iz) + vl(ix,iy,iz)
           vl_z(ix,iy) = vl_z(ix,iy) + vl(ix,iy,iz)

        end do
        end do
        end do

      else
         vl(:,:,:) = 1.0d0
         vl_x(:,:) = 1.0d0
         vl_y(:,:) = 1.0d0
         vl_z(:,:) = 1.0d0

      endif

*-----------------------------------------------------------------------

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

      nrst = 2*nfile+3 ! S.H. 2021.8.15
      allocate(anatalrst(np,ne+1,na+1,nt+1,(nx+1)*(ny+1)*(nz+1),1,nrst))

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
            do 100 it = 1, nt
            do 100 ia = 1, na
            do 100 ip = 1, np
            do ntf = 1, nfile

           if( tr(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf)) .gt. 0.d0 ) then

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

                  cc = rtfac(m) / vm / ew(ie) / aw(ia)

                  call calc_stdev(m,Xa,sigx,
     &                 tr(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf)),
     &                 tr(ip,ie,it,ia,icf(ix,iy,iz),iat(2,ntf)),
     &                 cc)

                  tr(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf)) = Xa
                  tr(ip,ie,it,ia,icf(ix,iy,iz),iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                if( tr(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf)) .gt. cmax )
     &                   cmax = tr(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf))

                if( tr(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf)) .lt. cmin )
     &                   cmin = tr(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,it,ia,icf(ix,iy,iz),iat(2,ntf)) = 0.d+0

               end if

! sumover
               call ppdctxyz_sumover_stdev_ntf(0,m,ntf,
     %               ip,ie,it,ia,ix,iy,iz,
     &               rtfac(m),ew(ie),1.0d0,aw(ia),vm,
     &               ew_sum,1.0d0,aw_sum,
     &               vl_x(iy,iz),vl_y(ix,iz),vl_z(ix,iy))


            end do
  100       continue
      enddo  ! xyz loop T.Sato 2021/04/15

!$OMP PARALLEL
!$OMP& private(ipomp,ix,iy,iz,it,ie,ip,ia,ntf,itmpxyz,rdata,answer,rerr)
!$OMP& private(iax,ix_a,iy_a,iz_a,it_a,ie_a,ia_a,ido_ana,ioe)
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
            do 101 it = 1, nt
            do 101 ia = 1, na
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

               do ntf = 1, nfile
                  rdata(1,ntf)=tr(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf))
                  rdata(2,ntf)=tr(ip,ie,it,ia,icf(ix,iy,iz),iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),1,1) = answer
              anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),1,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf)=tr(ip,ie,it,ia,icf(ix,iy,iz),iat(1,ntf))
                  rdata(2,ntf)=tr(ip,ie,it,ia,icf(ix,iy,iz),iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),1,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(ip,ie,ia,it,icf2(ix,iy,iz),1,iat(ioe,ntf))
     &              = tr(ip,ie,it,ia,icf(ix,iy,iz),iat(ioe,ntf))
                end do
               end do

            end if

!sumover
            do iax = 1,itaxn(m)
              ie_a = ie
              it_a = it
              ia_a = ia
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
              elseif((itaxs(m,iax) == 12 .or. itaxs(m,iax) == 13)
     &               .and. ia == 1) then
                ia_a = na + 1
                ido_ana = 1
              endif

              if(ido_ana == 1) then

                 call get_pdctxyz_tr_sum_data(m,iax,nfile,
     &                 ip,ie,it,ia,ix,iy,iz,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   anatalrst(ip,ie_a,ia_a,it_a,icf2(ix_a,iy_a,iz_a),1,1)
     &             = answer
                   anatalrst(ip,ie_a,ia_a,it_a,icf2(ix_a,iy_a,iz_a),1,2)
     &             = rerr

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  anatalrst(ip,ie_a,ia_a,it_a,icf2(ix_a,iy_a,iz_a),1,1) ! mean
     &              = fmval
                 if ( fmval .gt. 0.0d0 ) then
                   anatalrst(ip,ie_a,ia_a,it_a,icf2(ix_a,iy_a,iz_a),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                 else
                   anatalrst(ip,ie_a,ia_a,it_a,icf2(iz_a,iy_a,iz_a),1,2) ! total uncertainty
     &                 = 0.0d0
                 end if

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      anatalrst(ip,ie_a,ia_a,it_a,icf2(ix_a,iy_a,iz_a),
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

            if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &          igsh .eq. 0 ) then

               write(fnume,'(i3.3)') nobch
C S.H. revised for Dump-Restart on 2014/5/7
               fname = ctfln(m,iax)(1:itfll(m,iax))//fnume

            else if ( itall.eq.4 .and. igsh.eq.0 ) then
               write(fnume,'(i3.3)') nobch
               if ( npe.gt.1 ) write (fnume,'(i3.3)') nobch/(npe-1)
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
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        angle(p), energy, x, y, z, time axes
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 12 .or. itaxs(m,iax) .eq. 13 .or.
     &       itaxs(m,iax) .eq.  1 .or.
     &       itaxs(m,iax) .eq.  3 .or.
     &       itaxs(m,iax) .eq.  4 .or.
     &       itaxs(m,iax) .eq.  5 .or.
     &       itaxs(m,iax) .eq. 11 ) then

          select case( itaxs(m,iax) )
           case ( 12 )
            iDaxis = 12  ! cos of angle(p) axis
           case ( 13 )
            iDaxis = 13  ! the of angle(p) axis
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
          call anatal_rearrange_sum(itaty(m),iMeVperu,3,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      na,    nt,     nx, ny, nz,    0, nrst,
     &         eb,ew,    ab,aw,  tb,tw,  xm, ym, zm,1,
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
               write(iot,'( "y: Number ",a29)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

        if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 15 ) then ! T.Sato 2023/01/09
               if( itety(m) .eq.  3 .or. itety(m) .eq.  5 .or.
     &             itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &             itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &             itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 11 ) then ! time axis
               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax).eq.12 .or. itaxs(m,iax).eq.13 ) then ! angle axis
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

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
             if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: x",12x,"n",12x,
     &                    6(a10,"),h0",a3," n "))')
     &                     ( chp(i), chl(i) , i = 1, np )
             else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

             end if
            else
             if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: n",12x,"x",12x,
     &           6(a1,i1,a9,"),h0",a3
     &          ," ny2",i1," dy",i1,"=[y",i1,"*y2",i1,"] "))')
     &           ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &           ,i,i,i,i, i = 1, np )
             else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a1,i0,a9,"),hh0",a3
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                    ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                    ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

            write(iot,'( "#  lower        upper  ",3x,
     &                 1000(a1,2x,a8,4x,"r.err "))')
     &                 ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
             if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: x",12x,"n",12x,
     &                    6(a10,"),h0",a3," n n n "))')
     &                     ( chp(i), chl(i) , i = 1, np )
             else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n n n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if

            else
             if( itaxs(m,iax) .eq. 13 ) then ! theta axis
                write(iot,'( "h: n",12x,"x",12x,
     &           6(a1,i1,a9,"),h0",a3
     &          ," ny2",i1," dy",i1,"=[y",i1,"*y2",i1,"] n n "))')
     &           ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &           ,i,i,i,i, i = 1, np )
             else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a1,i0,a9,"),hh0",a3
     &                ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &                    ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                    ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
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
                if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 15 ) then ! T.Sato 2023/01/09
                     if( itunt(m) .eq.  3 .or. itunt(m) .eq.  4 .or.
     &                   itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &                   itunt(m) .eq. 13 .or. itunt(m) .eq. 14 .or.
     &                   itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &                   itunt(m) .eq. 23 .or. itunt(m) .eq. 24 .or.
     &                   itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &                   itunt(m) .eq. 33 .or. itunt(m) .eq. 34 .or.
     &                   itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then
                        tott(ip,1) = tott(ip,1) / voll
                     end if

                  else if( itaxs(m,iax) .eq. 11 ) then ! time axis
                  if( ( itunt(m) .gt. 10 .and. itunt(m) .lt. 20 ) .or.
     &                ( itunt(m) .gt. 30 .and. itunt(m) .lt. 40 ) )
     &            then
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

               if( itaxs(m,iax) .eq. 1 .or. itaxs(m,iax) .eq. 11
     &              .or. itaxs(m,iax) .eq. 12
     &              .or. itaxs(m,iax) .eq. 13
     &              .or. itaxs(m,iax) .eq. 15 ) then  ! T.Sato 2023/01/09
                  write(iot,'(/"#   sum over",14x,6(1pe13.4,0pf8.4))')
     &                 (tott(ip,1),tott(ip,2),ip=1,np)

               else if( itaxs(m,iax).eq.3 .or. itaxs(m,iax).eq.4
     &                 .or. itaxs(m,iax).eq.5 ) then
                  if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then
              write(iot,'(/"#   sum over ",1pe13.4,
     &              1000(1pe13.4,0pf8.4,16x))')
     &             voll,(tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
                  else
               write(iot,'(/"#   sum over ",   13x ,
     &               1000(1pe13.4,0pf8.4,16x))')
     &              (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
                  end if

               end if

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
          call anatal_rearrange_sum(itaty(m),iMeVperu,3,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      na,    nt,     nx, ny, nz,    0, nrst,
     &         eb,ew,    ab,aw,  tb,tw,  xm, ym, zm,1,
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

!OBINATA(2012.9.11): output *.err
               write(iot,'(1p10e11.3)')
     &       ( ( anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                ijaxs = 1, nijaxs ), ijaxs2 = nijaxs2, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

              write(iot,'(/"# axis2      axis1    ",
     &                    "  Number       r.err")')

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

!OBINATA(2012.9.11): output *.err
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
               iaxs = 3
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
         write(iot,'("y: Number ",a32)') hsunit(itunt(m))
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

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 5
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               if ( itaxs(m,iax) .eq. 7 ) then ! when axis = xy (x-y)
                  bmpfIType = (/ 'a', 't', 'e', 'p', 'z' /)
               else if ( itaxs(m,iax) .eq. 8 ) then ! when axis = yz (z-y)
                  bmpfIType = (/ 'a', 't', 'e', 'p', 'x' /)
               else if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xz (z-x)
                  bmpfIType = (/ 'a', 't', 'e', 'p', 'y' /)
               end if
                  bmpfIndex = (/ nag, ntg, neg, npg, nijaxs3 /)
                  bmpWidth  = nijaxs
                  bmpHeight = nijaxs2

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
                  do ie = 1, neg
                  do ia = 1, nag

                     write(isunit_vtk_meta) ip, ie, ia

                   if ( itaxs(m,iax) .eq. 7 ) then ! when axis = xy (x-y)
                     write(isunit_vtk)
     &                    ( ( ( anataldata(ip,
     &                    im+((ijaxs3-1)+(it-1+(ie-1)*ntg)*nijaxs3)*nm,
     &                    ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                             ijaxs = 1, nijaxs ),
     &                             ijaxs2 = 1, nijaxs2 ),
     &                             ijaxs3 = 1, nijaxs3 )
                   else if ( itaxs(m,iax) .eq. 8 ) then ! when axis = yz (z-y)
                     write(isunit_vtk)
     &                    ( ( ( anataldata(ip,
     &                    im+((ijaxs3-1)+(it-1+(ie-1)*ntg)*nijaxs3)*nm,
     &                    ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                             ijaxs3 = 1, nijaxs3 ),
     &                             ijaxs2 = 1, nijaxs2 ),
     &                             ijaxs = 1, nijaxs )
                   else if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xz (z-x)
                     write(isunit_vtk)
     &                    ( ( ( anataldata(ip,
     &                    im+((ijaxs3-1)+(it-1+(ie-1)*ntg)*nijaxs3)*nm,
     &                    ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                             ijaxs2 = 1, nijaxs2 ),
     &                             ijaxs3 = 1, nijaxs3 ),
     &                             ijaxs = 1, nijaxs )
                   end if

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
             end if

            end if

            deallocate( fg )
            deallocate( fgaxs )
            deallocate( fgaxs2 )
            deallocate( fgaxs3 )
            deallocate( anataldata )
            deallocate( delvol )

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
      deallocate( vl )

      deallocate (vl_x,vl_y,vl_z)
      deallocate (anatalrst)
      deallocate (ew,tw,aw)
      deallocate (rdata)


      return
      end subroutine anatal_ppdctxyz

!***********************************************************************
!                                                                      *
! sumover subroutine                                                   *
!                                                                      *
!***********************************************************************

!***********************************************************************
!                                                                      *
      subroutine get_pdctreg_tr_sum_data(m,iax,nfile,
     &           ip,ie,it,ia,ir,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,it,ia,ir
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

        call get_pdctreg_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &     itanm_sum(m,iax),itrgn_sum(m,iax),
     &     ip,ie,it,ia,ir,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_pdctreg_tr_sum_data_sub(tr_sum,
     &          np,ne,nt,na,nr,ip,ie,it,ia,ir,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,nt,na,nr,ip,ie,it,ia,ir
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,ne,nt,na,nr,2)
 
      do i=1,2
        rdata(i)= tr_sum(ip,ie,it,ia,ir,i)
      enddo

      end

!***********************************************************************
!                                                                      *
      subroutine get_pdctrz_tr_sum_data(m,iax,nfile,
     &           ip,ie,it,ia,ir,iz,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,it,ia,ir,iz
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

        call get_pdctrz_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &     itanm_sum(m,iax),
     &     itrnm_sum(m,iax),itznm_sum(m,iax),
     &     ip,ie,it,ia,ir,iz,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_pdctrz_tr_sum_data_sub(tr_sum,
     &          np,ne,nt,na,nr,nz,ip,ie,it,ia,ir,iz,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,nt,na,nr,nz,ip,ie,it,ia,ir,iz
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,ne,nt,na,nr,nz,2)

      do i=1,2
        rdata(i)= tr_sum(ip,ie,it,ia,ir,iz,i)
      enddo

      end

!***********************************************************************
!                                                                      *
      subroutine get_pdctxyz_tr_sum_data(m,iax,nfile,
     &           ip,ie,it,ia,ix,iy,iz,rdata)
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

        call get_pdctxyz_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &     itanm_sum(m,iax),
     &     itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &     ip,ie,it,ia,ix,iy,iz,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_pdctxyz_tr_sum_data_sub(tr_sum,
     &          np,ne,nt,na,nx,ny,nz,ip,ie,it,ia,ix,iy,iz,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,nt,na,nx,ny,nz,ip,ie,it,ia,ix,iy,iz
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,ne,nt,na,nx*ny*nz,2)

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

      do i=1,2
        rdata(i)= tr_sum(ip,ie,it,ia,icf(ix,iy,iz),i)
      enddo

      end
