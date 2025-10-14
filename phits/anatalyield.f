!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_yieldreg(m,ntf,
     &                             mn,mz,mm,nr, ! S.H. added mm (2022.3.11)
     &                             sumfactor,weightRate,
     &                             trRES,tranatal,nfile,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     setting of tranatal of t-yield with reg mesh                     *
!     uploaded by S.Hashimoto on 2022/1/21                             *
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

      integer           mn
      integer           mz
      integer           nr
      integer           mm ! S.H. added mm (2022.3.11)

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nr,mz,mn,0:mm,2) ! S.H. 2-> mm (2022.3.11)
      integer           nfile
      double precision  tranatal(nr,mz,mn,0:mm,2*nfile) ! S.H. 2-> mm (2022.3.11)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           im

      integer           iz
      integer           in
      integer           il
      integer           ir
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
       do iz = 1, mz
        do in = 1, mn
         do il = 0, mm ! S.H. 2-> mm (2022.3.11)

             if( trRES(ir,iz,in,il,1) .ne. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(ir,iz,in,il,1),
     &                        trRES(ir,iz,in,il,2),
     &                        1.0d+0)
              trRES(ir,iz,in,il,1) = Xa
              trRES(ir,iz,in,il,2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(ir,iz,in,il,iat(1,ntf)) =
     &           trRES(ir,iz,in,il,1)

              ! (sig_xj)**2
              tranatal(ir,iz,in,il,iat(2,ntf)) =
     &        (  trRES(ir,iz,in,il,2)
     &         * trRES(ir,iz,in,il,1) )**2

              ! sig_x
              tranatal(ir,iz,in,il,iat(2,ntf)) = dsqrt(
     &        tranatal(ir,iz,in,il,iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(ir,iz,in,il,iat(2,ntf)) =
     &       (tranatal(ir,iz,in,il,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(ir,iz,in,il,iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(ir,iz,in,il,iat(1,ntf)) =
     &        tranatal(ir,iz,in,il,iat(1,ntf)) * resc2(m)

             if( tranatal(ir,iz,in,il,iat(1,ntf)) > cmax )
     &             cmax = tranatal(ir,iz,in,il,iat(1,ntf))
             if( tranatal(ir,iz,in,il,iat(1,ntf)) < cmin )
     &             cmin = tranatal(ir,iz,in,il,iat(1,ntf))

         end do     ! il loop end
        end do      ! in loop end
       end do       ! iz loop end
      end do        ! ir loop end

!-----------------------------------------------------------------------
      return

      end subroutine anatal_calc_anova_yieldreg


!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_yieldrz(m,ntf,
     &                             mn,mz,mm,nz,nr, ! S.H. added mm (2022.3.11)
     &                             sumfactor,weightRate,
     &                             trRES,tranatal,nfile,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     setting of tranatal of t-yield with r-z mesh                     *
!     uploaded by S.Hashimoto on 2022/1/21                             *
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

      integer           mn
      integer           mz
      integer           nz
      integer           nr
      integer           mm ! S.H. added mm (2022.3.11)

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nr,nz,mz,mn,0:mm,2) ! S.H. 2-> mm (2022.3.11)
      integer           nfile
      double precision  tranatal(nr,nz,mz,mn,0:mm,2*nfile) ! S.H. 2-> mm (2022.3.11)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           il
      integer           in
      integer           iz
      integer           jz
      integer           jr
      double precision  Xa, sigx
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do jr = 1, nr
       do jz = 1, nz
        do iz = 1, mz
         do in = 1, mn
          do il = 0, mm ! S.H. 2-> mm (2022.3.11)

             if( trRES(jr,jz,iz,in,il,1) .ne. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(jr,jz,iz,in,il,1),
     &                        trRES(jr,jz,iz,in,il,2),
     &                        1.0d+0)
              trRES(jr,jz,iz,in,il,1) = Xa
              trRES(jr,jz,iz,in,il,2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(jr,jz,iz,in,il,iat(1,ntf)) =
     &           trRES(jr,jz,iz,in,il,1)

              ! (sig_xj)**2
              tranatal(jr,jz,iz,in,il,iat(2,ntf)) =
     &        (  trRES(jr,jz,iz,in,il,2)
     &         * trRES(jr,jz,iz,in,il,1) )**2

              ! sig_x
              tranatal(jr,jz,iz,in,il,iat(2,ntf)) = dsqrt(
     &        tranatal(jr,jz,iz,in,il,iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(jr,jz,iz,in,il,iat(2,ntf)) =
     &       (tranatal(jr,jz,iz,in,il,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(jr,jz,iz,in,il,iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(jr,jz,iz,in,il,iat(1,ntf)) =
     &        tranatal(jr,jz,iz,in,il,iat(1,ntf)) * resc2(m)

           if( tranatal(jr,jz,iz,in,il,iat(1,ntf)) > cmax )
     &          cmax = tranatal(jr,jz,iz,in,il,iat(1,ntf))
           if( tranatal(jr,jz,iz,in,il,iat(1,ntf)) < cmin )
     &          cmin = tranatal(jr,jz,iz,in,il,iat(1,ntf))

          end do     ! il loop end
         end do      ! in loop end
        end do       ! iz loop end
       end do        ! jz loop end
      end do         ! jr loop end

!-----------------------------------------------------------------------
      return

      end subroutine anatal_calc_anova_yieldrz


!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_yieldxyz(m,ntf,
     &                             mn,mz,mm,nz,ny,nx, ! S.H. added mm (2022.3.11)
     &                             sumfactor,weightRate,
     &                             trRES,tranatal,nfile,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     setting of tranatal of t-yield with xyz mesh                     *
!     uploaded by S.Hashimoto on 2022/1/21                             *
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

      integer           mn
      integer           mz
      integer           nz
      integer           ny
      integer           nx
      integer           mm ! S.H. added mm (2022.3.11)

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nx,ny,nz,mz,mn,0:mm,2) ! S.H. 2-> mm (2022.3.11)
      integer           nfile
      double precision  tranatal(nx,ny,nz,mz,mn,0:mm,2*nfile) ! S.H. 2-> mm (2022.3.11)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           jx
      integer           jy
      integer           jz
      integer           iz
      integer           in
      integer           il
      double precision  Xa, sigx
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do jx = 1, nx
       do jy = 1, ny
        do jz = 1, nz
         do iz = 1, mz
          do in = 1, mn
           do il = 0, mm ! S.H. 2-> mm (2022.3.11)

             if( trRES(jx,jy,jz,iz,in,il,1) .ne. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(jx,jy,jz,iz,in,il,1),
     &                        trRES(jx,jy,jz,iz,in,il,2),
     &                        1.0d+0)
              trRES(jx,jy,jz,iz,in,il,1) = Xa
              trRES(jx,jy,jz,iz,in,il,2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(jx,jy,jz,iz,in,il,iat(1,ntf)) =
     &           trRES(jx,jy,jz,iz,in,il,1)

              ! (sig_xj)**2
              tranatal(jx,jy,jz,iz,in,il,iat(2,ntf)) =
     &        (  trRES(jx,jy,jz,iz,in,il,2)
     &         * trRES(jx,jy,jz,iz,in,il,1) )**2

              ! sig_x
              tranatal(jx,jy,jz,iz,in,il,iat(2,ntf)) = dsqrt(
     &        tranatal(jx,jy,jz,iz,in,il,iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(jx,jy,jz,iz,in,il,iat(2,ntf)) =
     &       (tranatal(jx,jy,jz,iz,in,il,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(jx,jy,jz,iz,in,il,iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(jx,jy,jz,iz,in,il,iat(1,ntf)) =
     &        tranatal(jx,jy,jz,iz,in,il,iat(1,ntf)) * resc2(m)

            if( tranatal(jx,jy,jz,iz,in,il,iat(1,ntf)) > cmax )
     &            cmax = tranatal(jx,jy,jz,iz,in,il,iat(1,ntf))
            if( tranatal(jx,jy,jz,iz,in,il,iat(1,ntf)) < cmin )
     &            cmin = tranatal(jx,jy,jz,iz,in,il,iat(1,ntf))

           end do    ! il loop end
          end do     ! in loop end
         end do      ! iz loop end
        end do       ! jz loop end
       end do        ! jy loop end
      end do         ! jx loop end

!-----------------------------------------------------------------------
      return

      end subroutine anatal_calc_anova_yieldxyz
************************************************************************


************************************************************************
*                                                                      *
      subroutine anatal_pyildreg(m,mz,mn,mm,nr,mr,nn,kr,nt,ikzz,iknn, ! S.H. added mm (2022.3.11)
     &                    tr,nfile,weightRate,
     &                    nvl,ivl,rvl,
     &                    nx,ny,nz,xm,ym,zm,igsh,idasa,manatally)
*                                                                      *
*       anatally calculation and output of t-yield with reg mesh       *
*                                                                      *
************************************************************************

      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

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
      real(8) bzz(mz+1), bnn(mn+1)
!      dimension   tm(maxnt+maxpt,2*nfile+3)
!      dimension   vl3d(nr,1,1)
!      dimension   vl(nr)
!      dimension   lr(nr)
      dimension   tr(nr,mz,mn,0:mm,2*nfile) ! S.H. 2-> mm (2022.3.11)
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
!      dimension   anatalrst(mm+1,mz+1,mn+1,1,nr*1*1,1,2*nfile+3) ! S.H. 3-> mm+1 (2022.3.11)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)
      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*21,cijaxs*21,cijaxs2*21
      integer iDaxis
!      dimension   dlr(nr+1) !  integer lr(nr) ---> double dlr(nr)
      real(8),allocatable :: tm(:,:),vl3d(:,:,:),vl(:),rdata(:,:),
     &                       val(:),dlr(:)
      integer,allocatable :: lr(:)

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

      dimension dummy(1)
      data dummy/1.0d0/

*-----------------------------------------------------------------------
      character(28) cfmt                !FURUTA20200615
      character(12) cir                 !FURUTA20200615
      common /redufmt/ iredufmt(itlmax) !FURUTA20200615
*-----------------------------------------------------------------------
      character rpa*1
      data rpa /'}'/
      character yen*1
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------
      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )
      allocate (vl3d(nr,1,1),vl(nr),rdata(2,nfile),val(nr),dlr(nr+1))
      allocate (lr(nr))

      nrst = 2*nfile+3
      allocate (anatalrst(mm+1,mz+1,mn+1,1,nr*1*1,1,nrst))


*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)

ccse 2021.08 add (use anatal_rearrange sub.)
               do ir = 1, nr
                  dlr(ir) = dble(lr(ir))
               end do
               dlr(nr+1) = dble(lr(nr))

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

            do 100 ir = 1, nr

               if( itunt(m) .eq. 1 ) then

                  cc = rtfac(m)

               else if( itunt(m) .eq. 2 ) then

                  cc = rtfac(m) / vl(ir)

               end if

            do 100 iz = 1, mz
            do 100 in = 1, mn
            do 100 il = 0, mm ! S.H. 2-> mm (2022.3.11)
            do ntf = 1, nfile

               if( tr(ir,iz,in,il,iat(1,ntf)) .ne. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ir,iz,in,il,iat(1,ntf)),
     &                            tr(ir,iz,in,il,iat(2,ntf)),
     &                            cc)

                  tr(ir,iz,in,il,iat(1,ntf)) = Xa
                  tr(ir,iz,in,il,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ir,iz,in,il,iat(1,ntf)) .gt. cmax )
     &                 cmax = tr(ir,iz,in,il,iat(1,ntf))

                  if( tr(ir,iz,in,il,iat(1,ntf)) .lt. cmin )
     &                 cmin = tr(ir,iz,in,il,iat(1,ntf))

               else

                  if( il .eq. 0 ) isdz = 1
                  tr(ir,iz,in,il,iat(2,ntf)) = 0.0

               end if

            end do
  100       continue

!$OMP PARALLEL
!$OMP& private(ipomp,ir,iz,in,il,ntf,isdz,rdata,answer,rerr)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
            do ir = 1, nr
            do 101 iz = 1, mz
            do 101 in = 1, mn
            do 101 il = 0, mm ! S.H. 2-> mm (2022.3.11)

            if( manatally .eq. 0 ) then ! user defined analysis

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ir,iz,in,il,iat(1,ntf))
                  rdata(2,ntf) = tr(ir,iz,in,il,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(il+1,iz,in,1,ir,1,1) = answer
              anatalrst(il+1,iz,in,1,ir,1,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then ! systematic uncertainty

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ir,iz,in,il,iat(1,ntf))
                  rdata(2,ntf) = tr(ir,iz,in,il,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(il+1,iz,in,1,ir,1,1) ! mean
     &              = fmval
               if ( fmval .ne. 0d0 ) then
                  anatalrst(il+1,iz,in,1,ir,1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(il+1,iz,in,1,ir,1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(il+1,iz,in,1,ir,1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(il+1,iz,in,1,ir,1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(il+1,iz,in,1,ir,1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(il+1,iz,in,1,ir,1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(il+1,iz,in,1,ir,1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(il+1,iz,in,1,ir,1,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then ! c-value dependence

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(il+1,iz,in,1,ir,1,iat(ioe,ntf))
     &              = tr(ir,iz,in,il,iat(ioe,ntf))
                end do
               end do

            end if

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
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        reg, mass, charge axes
*-----------------------------------------------------------------------

         do iz = 1, mz
            bzz(iz) = iz
         end do
         do in = 1, mn
            bnn(in) = in
         end do
         do ir = 1, nr
            vl3d(ir,1,1) = vl(ir)
         end do

*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1
     &        .or. itaxs(m,iax) .eq. 2
     &        .or. itaxs(m,iax) .eq. 7 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 15   ! mass axis
           case ( 2 )
            iDaxis = 16   ! reg axis
           case ( 7 )
            iDaxis = 17   ! charge axis
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.19

C rearrange data from anatalrst to anataldata
          call anatal_rearrange(0,0,1,
     &         iDaxis,ntaxis,anatalrst,
     &         mm+1,     mz,      mn,    0,      nr,  1,  1,    0, nrst, ! S.H. 3-> mm+1 (2022.3.11)
     &         bzz, 1,  bnn,1,  1, 1,  dlr,  dummy,  dummy,1,
     &         vl3d,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          allocate( tm(nijaxs,2*nfile+3) )

          inum = 0

          if( nn .eq. 0 ) then
             nc = 1
          else
             nc = nn
          end if

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

           ic = ij(2)

*-----------------------------------------------------------------------

           iz = 0
           ia = 0
           in = 0
           if( nn .gt. 0 ) then

              iz = nt(ic) / 1000
              ia = nt(ic) - iz * 1000

              if( ia .gt. 0 ) then
                 in = ia - iz
                 call chname(idum,ia,iz,chau)
              end if

           end if

*-----------------------------------------------------------------------

           do itmp2 = 1, nrst
           do itmp1 = 1, nijaxs
              tm(itmp1,itmp2) = 0.0d+0
           end do
           end do

           seka = 0.0
           sera = 0.0
           seva = 0.0
           voll = 0.0

           do ijaxs = 1, nijaxs ! (1)

            if( itunt(m) .eq. 1 ) then

               vm = 1.0

            else if( itunt(m) .eq. 2 ) then

               vm = delvol(ianataldata,ijaxs)

            end if

            voll = voll + delvol(ianataldata,ijaxs)

*-----------------------------------------------------------------------

            if( nn .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

              vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
              tm(ijaxs,1) = tm(ijaxs,1) + vn
              tm(ijaxs,2) = tm(ijaxs,2)
     &             + ( vn * anataldata(il+1,ianataldata,ijaxs,2) )**2
              if ( manatally .eq. 1 ) then ! systematic uncertainty
                 tm(ijaxs,3) = tm(ijaxs,3)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,3) )**2
                 tm(ijaxs,4) = tm(ijaxs,4)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,4) )**2
              end if

             end do

            else if( ia .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

              vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
              tm(ijaxs,1) = tm(ijaxs,1) + vn
              tm(ijaxs,2) = tm(ijaxs,2)
     &             + ( vn * anataldata(il+1,ianataldata,ijaxs,2) )**2
              if ( manatally .eq. 1 ) then ! systematic uncertainty
                 tm(ijaxs,3) = tm(ijaxs,3)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,3) )**2
                 tm(ijaxs,4) = tm(ijaxs,4)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,4) )**2
              end if

             end do

            else

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

              vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
              tm(ijaxs,1) = tm(ijaxs,1) + vn
              tm(ijaxs,2) = tm(ijaxs,2)
     &             + ( vn * anataldata(il+1,ianataldata,ijaxs,2) )**2
              if ( manatally .eq. 1 ) then ! systematic uncertainty
                 tm(ijaxs,3) = tm(ijaxs,3)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,3) )**2
                 tm(ijaxs,4) = tm(ijaxs,4)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,4) )**2
              end if
             end do

            end if

            seka = seka + tm(ijaxs,1)
            sera = sera + tm(ijaxs,2)
            seva = seva + vm

            if( tm(ijaxs,1) .gt. 0.0 ) then

             tm(ijaxs,2) = sqrt( tm(ijaxs,2) ) / tm(ijaxs,1)
             tm(ijaxs,1) = tm(ijaxs,1) / vm
             if ( manatally .eq. 1 ) then ! systematic uncertainty
                tm(ijaxs,3) = sqrt( tm(ijaxs,3) ) / tm(ijaxs,1)
                tm(ijaxs,4) = sqrt( tm(ijaxs,4) ) / tm(ijaxs,1)
             end if

            end if

           end do                 ! ijaxs = 1, nijaxs (1)

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
               write(iot,'( "y: Number ",a15)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

           write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

           if( itanl(m) .gt. 0 ) then
              write(iot,'( "p: ",200a1)')
     &             ( itang(m)(i:i),i = 1, itanl(m) )
           end if

           if( itsans(m) .gt. 0 ) then
              call write_sangel(iot,m,0)
           end if

           if ( manatally .eq. 0 ) then ! user defined analysis
            if( itaxs(m,iax) .eq. 1 ) then ! mass axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5    y,hl0       n")')

             else
                write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

             end if

             write(iot,'( "#  mass     number      r.err")')

            else if( itaxs(m,iax) .eq. 2 ) then ! reg axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                    "y,l3        n")')

             else
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                    "y1,l3       ny2 dy1=[y1*y2]")')

             end if

             write(iot,'( "#  num    reg     volume  ",
     &                    "   number      r.err")')

            else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5    y,hl0       n")')

             else
                write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

             end if

             write(iot,'( "#  charge   number      r.err")')

            end if

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if( itaxs(m,iax) .eq. 1 ) then ! mass axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5   y,hl0       n n n")')

             else
                write(iot,'( "h: x-0.5   y1,hl0      ",
     &                      "ny2 dy1=[y1*y2] n n")')

             end if

             write(iot,'( "#  mass    number      "
     &                   ,"r.err(tot, syst, stat)")')

           else if( itaxs(m,iax) .eq. 2 ) then ! reg axis
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h:",3x,"x",7x,"n",3x,"n",12x,
     &                    "y1,l3",7x,"n")')

            else
               write(iot,'( "h:",3x,"x",7x,"n",3x,"n",12x,
     &                    "y1,l3",7x,"ny2 dy1=[y1*y2] n n")')

            end if

            write(iot,'( "#  num     reg   volume  ",
     &                    "     number      r.err(tot, syst, stat)")')

           else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: x-0.5   y,hl0       n n n")')

            else
               write(iot,'( "h: x-0.5   y1,hl0      ",
     &                      "ny2 dy1=[y1*y2] n n")')

            end if

            write(iot,'( "#charge    number      "
     &                   ,"r.err(tot, syst, stat)")')

           end if

          else if ( manatally .eq. 2 ) then ! c-value dependence
           if ( iteps(m) .ne. 2 ) then
              write(iot,'( "h:   x",10x,
     &                    "y,l3        n")')

           else
              write(iot,'( "h:   x",10x,
     &                    "y1,l3       ny2 dy1=[y1*y2]")')

           end if

           write(iot,'( "#  c-value  ",
     &                    "   number      r.err")')

          end if

*-----------------------------------------------------------------------

           do ijaxs = 1, nijaxs
              if( tm(ijaxs,1) .gt. 0.0d0 ) goto 201
           end do

 201       im = ijaxs

           do ijaxs = nijaxs, im + 1, -1
              if( tm(ijaxs,1) .gt. 0.0d0 ) goto 211
           end do

           jm = im
 211       jm = ijaxs

           im = max( 1, im - 2 )
           jm = jm + 2

*-----------------------------------------------------------------------

           if( itaxs(m,iax).eq.1 .or. itaxs(m,iax).eq.7 ) then ! mass, charge axis
            if( im .eq. 1 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(3x,f4.1,1x,1pe13.4,0p1f8.4)')
     &               0.5, 0.0, 0.0

             else if ( manatally .eq. 1 ) then ! systematic uncertainty
                write(iot,'(3x,f4.1,1x,1pe13.4,0p3f8.4)')
     &               0.5, 0.0, 0.0, 0.0, 0.0

             else if ( manatally .eq. 2 ) then ! c-value dependence

             end if

            end if
           end if

           do ijaxs = 1, nijaxs ! (2)

            if ( manatally .eq. 0 ) then ! user defined analysis
             if( itaxs(m,iax) .eq. 1 ) then ! mass axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p1f8.4)')
     &                ijaxs+iz, (tm(ijaxs,irst),irst=1,2)
              end if

             else if( itaxs(m,iax) .eq. 2 ) then ! reg axis
                write(iot,'(1x,i5,1x,i7,1pe13.4,1pe13.4,0p1f8.4)')
     &               ijaxs, lr(ijaxs), delvol(ianataldata,ijaxs)
     &               ,(tm(ijaxs,irst),irst=1,2)

             else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p1f8.4)')
     &                ijaxs, (tm(ijaxs,irst),irst=1,2)

              end if

             end if

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
             if( itaxs(m,iax) .eq. 1 ) then ! mass axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                write(iot,'(3x,i4,1x,1pe13.4,0p3f8.4)')
     &               ijaxs+iz, (tm(ijaxs,irst),irst=1,4)
              end if

             else if( itaxs(m,iax) .eq. 2 ) then ! reg axis
                write(iot,'(1x,i5,1x,i7,1pe13.4,1pe13.4,0p3f8.4)')
     &               ijaxs, lr(ijaxs), delvol(ianataldata,ijaxs)
     &               ,(tm(ijaxs,irst),irst=1,4)

             else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p3f8.4)')
     &                ijaxs, (tm(ijaxs,irst),irst=1,4)
              end if

             end if

            else if ( manatally .eq. 2 ) then ! c-value dependence
               write(iot,'(3x,i4,1x,1pe13.4,0p1f8.4)')
     &              ijaxs, (tm(ijaxs,irst),irst=1,2)

            end if

           end do                 ! ijaxs = 1, nijaxs (2)

           if( seka .gt. 0.0 ) then
              sera = sqrt( sera ) / seka
              if( itunt(m) .eq. 2 ) seka = seka / seva
           end if

           if( itaxs(m,iax) .eq. 1 ) then ! mass axis
            write(iot,'(/"#   sum ",1pe13.4,0pf8.4)') seka,sera

            if( itunt(m) .eq. 2 ) then
               write(iot,'("#   vol sum = ",1p1e13.4,"  [cm^3]")') vm
            end if

           else if( itaxs(m,iax) .eq. 2 ) then ! reg axis
            write(iot,'(/"#   sum over  ",1pe13.4,1pe13.4,0pf8.4)')
     &             voll, seka, sera

            if( itunt(m) .eq. 2 ) then
               write(iot,'("#   vol sum = ",1p1e13.4,"  [cm^3]")') seva

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
              if ( ijtmp .eq. 2 ) then
                 if ( iz.gt.0)
     &                write(iot,'(a15,"&=&",i5,3x,"(",a4,")")')
     &                cij(ijtmp), iz, elmnt(iz)
              else if ( ijtmp .eq. 3 ) then
                 if ( ia .gt. 0 )
     &                write(iot,'(a15,"&=&",i5)') cij(ijtmp), ia
              else if ( ijtmp .eq. 5 ) then
                 write(iot,'(a15,"&=&",i5)')
     &                cij(ijtmp), idnint(fg(itmp))
                 write(iot,'(3x,"volume",6x,"&=&",1pe13.4," [cm^3]")')
     &                delvol(ianataldata,1)
              else
                 write(iot,'(a15,"&=&",i5)')
     &                cij(ijtmp), idnint(fg(itmp))
              end if
              itmp = itmp - ij(ijtmp) + nij(ijtmp)
             else
                write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                           /1pe13.4,2x,"$--$",1pe13.4)')
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
          deallocate( tm )

*-----------------------------------------------------------------------
*        chart axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 ) then

          iDaxis = 18         ! chart axis

C rearrange data from anatalrst to anataldata
          call anatal_rearrange(0,0,1,
     &         iDaxis,ntaxis,anatalrst,
     &         mm+1,     mz,      mn,    0,      nr,  1,  1,    0, nrst, ! S.H. 3-> mm+1 (2022.3.11)
     &         bzz, 1,  bnn,1,  1, 1,  dlr,  dummy,  dummy,1,
     &         vl3d,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          inum = 0

*-----------------------------------------------------------------------

          do ianataldata = 1, nanataldata

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

           ir = ij(5)

           do il = 0, mm         ! loop for isomeric level ! S.H. 2-> mm (2022.3.11)

*-----------------------------------------------------------------------

            itmax  = 0
            icmax  = 0
            inmax  = 0

            do iz = 1, maxpt
             do in = 1, maxnt

              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
                 ijaxs = in + (iz-1)*maxnt
              else
                 ijaxs = igetiznmp(iz,in,il,m)
              end if
              if( anataldata(il+1,ianataldata,ijaxs,1).gt.0d0 ) then
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

            if( itmax .eq. 0 ) exit

*-----------------------------------------------------------------------

            ireg = idnint(fg(ij(5)))

            write(iot,'(/"#",78("-"))')

            inum = inum + 1

            if( inum .eq. 1 .and. il .eq. 0 ) then
               write(iot,'( "#newpage:")')
            else
               write(iot,'( " newpage:")')
            end if

*-----------------------------------------------------------------------

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
                  write(iot,'(/"x: ",a21)') cijaxs
               else
                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))
               end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
               if( itayl(m) .eq. 0 ) then
                  write(iot,'(/"y: ",a21)') cijaxs2
               else
                  write(iot,'( "y: ",200a1)')
     &                 (itayt(m)(i:i),i=1,itayl(m))
               end if
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
     &              icmax+2, inmax+2

             do ijaxs2 = icmax+2, 1, -1
              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
               write(iot,'(1p10e11.3)')
     &       ( anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &               ijaxs = 1, inmax+2 )

              else
               write(iot,'(1p10e11.3)')
     &                ( anataldata(il+1,ianataldata,
     &                igetiznmp(ijaxs2,ijaxs,il,m),ioe),
     &               ijaxs = 1, inmax+2 )
              end if
             end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

             write(iot,'(/"#   Z    N  Mass    number    r.err")')

             do ijaxs2 = 1, icmax
              do ijaxs = 1, inmax

              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
               if( anataldata(il+1,ianataldata
     &                ,ijaxs+(ijaxs2-1)*nijaxs,1) .gt. 0.d0 ) then

                  write(iot,'(3i5,1pe13.4,0pf8.4)')
     &                 ijaxs2, ijaxs, ijaxs2+ijaxs,
     &           anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &           anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,2)

               end if

              else
               if( anataldata(il+1,ianataldata
     &                ,igetiznmp(ijaxs2,ijaxs,il,m),1) .gt. 0.d0 ) then

                  write(iot,'(3i5,1pe13.4,0pf8.4)')
     &                 ijaxs2, ijaxs, ijaxs2+ijaxs,
     &      anataldata(il+1,ianataldata,igetiznmp(ijaxs2,ijaxs,il,m),1),
     &      anataldata(il+1,ianataldata,igetiznmp(ijaxs2,ijaxs,il,m),2)
               end if
              end if

              end do
             end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

             write(iot,'(/"#   Z = 1 to ",i3/
     &                      "#   N = 1 to ",i3)')
     &              icmax+2, inmax+2

             write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'Z/N', ( dble( in ), in = 1, inmax+2 )

             do ijaxs2 = icmax+2, 1, -1
              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
                write(iot,'(1p1000e11.3)')
     &               dble( ijaxs2 ),
     &         ( anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &               ijaxs = 1, inmax+2 )

              else
                write(iot,'(1p1000e11.3)')
     &               dble( ijaxs2 ),
     &    ( anataldata(il+1,ianataldata,igetiznmp(ijaxs2,ijaxs,il,m),1),
     &               ijaxs = 1, inmax+2 )
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

      else

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else
        if( itazl(m) .eq. 0 ) then
           write(iot,'("y: Number ",a15)') hsunit(itunt(m))
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
       write(iot,'("wt: s(0.7)",/a1,"vspace{-3}")') yen
       itmp = 0
       do ijtmp = 2,8
        if (cij(ijtmp) .ne. 'F' ) then
         itmp = itmp + ij(ijtmp)
         if ( ibin(ijtmp) .eq. 0 ) then
            if ( ijtmp .eq. 5 ) then
               write(iot,'(a15,"&=&",i5)')
     &              cij(ijtmp), idnint(fg(itmp))
               write(iot,'(3x,"volume")')
               write(iot,'(5x,"&=&",1pe13.4," [cm^3]")')
     &              delvol(ianataldata,1)
            else
               write(iot,'(a15,"&=&",i5)')
     &              cij(ijtmp), idnint(fg(itmp))
            end if
            itmp = itmp - ij(ijtmp) + nij(ijtmp)
           else
              write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                           /1pe13.4,2x,"$--$",1pe13.4)')
     &             cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
              itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
           end if
        end if
       end do
       if(il .eq. 0) then       ! specify isomeric level at the bottom of figure
          write(iot,'("il = 0 (Ground state)}")')
       elseif(il .eq. 1) then
          write(iot,'("il = 1 (1st isomer)}")')
       else
          write(iot,'("il = 2 (2nd isomer)}")')
       endif

       write(iot,'("e:")')

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

           end do                    ! il = 0, 2 ! loop for isomeric level
          end do                    ! ianataldata = 1, nanataldata

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( anataldata )
          deallocate( delvol )

*-----------------------------------------------------------------------
*        dchain axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 13 ) then

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
                           t0_tr = anatalrst(il+1,iz,in,1,ir,1,1)
                        else
                           t0_tr =
     &                   anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,ir,1,1)
                        end if
                        if( t0_tr .ne. 0.0 ) goto 120

                     end do
                     end do

                     goto 170

  120             im = in

                     do in = maxnt, im + 1, -1
                     do ir = 1, nr

                        if( itnzn(m).eq.0 ) then
                           t0_tr = anatalrst(il+1,iz,in,1,ir,1,1)
                        else
                           t0_tr =
     &                   anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,ir,1,1)
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
     &           lr(ir), ( anatalrst(il+1,iz,i,1,ir,1,ioe), i = n3, n4 )
                   else
                     write(iot,'(i7,1p12e11.3)') lr(ir),
     &               ( anatalrst(0+1,igetiznmp(iz,i,il,m),1,1,ir,1,ioe),
     &                i = n3, n4 )
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
     &       lr(ir), ( anatalrst(il+1,iz,i,1,ir,1,2)*100.0, i = n3, n4 )
                      else
                        write(iou,'(i7,1p12e11.3)') lr(ir),
     &         ( anatalrst(0+1,igetiznmp(iz,i,il,m),1,1,ir,1,2)*100.0,
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
                  if(anatalrst(il+1,iz,in,1,ir,1,ioe).ne.0d0)then
                   write(iot,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                  (anatalrst(il+1,iz,in,1,ir,1,ioee),ioee=1,2)
                   if(itout(m).ne.0)then
                    write(iou,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                   anatalrst(il+1,iz,in,1,ir,1,2)
                   endif
                  endif
               else
                  if(anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,ir,1,ioe)
     &                .ne.0d0)then
                   write(iot,cfmt)ir,iz*10000+(in+iz)*10+il,
     &    (anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,ir,1,ioee),ioee=1,2)
                   if(itout(m).ne.0)then
                    write(iou,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                  anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,ir,1,2)
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

         end if   ! itaxs(m,iax)

*-----------------------------------------------------------------------

         close(iot)

         if( iteps(m) .ne. 0 .and. itaxs(m,iax) .ne. 13 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      deallocate (anatalrst)
      deallocate (vl3d,vl,rdata,val,dlr)
      deallocate (lr)

      return
      end


************************************************************************
*                                                                      *
      subroutine anatal_pyildrz(m,mz,mn,mm,nfr,nfz, ! S.H. added mm (2022.3.11)
     &                   nr,nz,nn,rm,zm,nt,ikzz,iknn,
     &                   tr,nfile,weightRate,
     &                   idasa,manatally)
*                                                                      *
*       anatally calculation and output of t-yield with r-z mesh       *
*                                                                      *
************************************************************************

      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

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
      real(8) bzz(mz+1), bnn(mn+1)
!      dimension   tm(maxnt+maxpt,2*nfile+3)
!      dimension   vl3d(nr,nz,1)
!      dimension   fm(nfz,nfr)
      dimension   tr(nr,nz,mz,mn,0:mm,2*nfile) ! S.H. 2-> mm (2022.3.11)
!      dimension   weightRate(nfile), rdata(2,nfile)
      dimension   weightRate(nfile)

      integer irst,nrst, itmpdata
!      dimension   anatalrst(mm+1,mz+1,mn+1,1,nr*nz*1,1,
!     &                     2*nfile+3) ! S.H. 3-> mm+1 (2022.3.11)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)
      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*21,cijaxs*21,cijaxs2*21
      integer iDaxis
      real(8),allocatable :: tm(:,:),vl3d(:,:,:),fm(:,:),rdata(:,:)

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

      data ipstep / 12 /
      character erfnm*100

      character chau*8
      character cha*1
      data cha /"'"/

      character rpa*1
      data rpa /'}'/
      character yen*1

      dimension dummy(1)
      data dummy/1.0d0/
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------

      vl(jr,jz) = pi * ( rm(jr+1)**2 - rm(jr)**2 )
     &     * ( zm(jz+1) - zm(jz) )

*-----------------------------------------------------------------------

      icf(jr,jz) = jr + ( jz - 1 ) * nr

*-----------------------------------------------------------------------

      yen  = char(92)
      igsh = 0
      allocate (vl3d(nr,nz,1),fm(nfz,nfr),rdata(2,nfile))

      nrst = 2*nfile+3
      allocate (anatalrst(mm+1,mz+1,mn+1,1,nr*nz*1,1,
     &                     nrst))


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

            do 100 jr = 1, nr
            do 100 jz = 1, nz

               if( itunt(m) .eq. 1 ) then

                  cc = rtfac(m)

               else if( itunt(m) .eq. 2 ) then

                  cc = rtfac(m) / vl(jr,jz)

               end if

            do 100 iz = 1, mz
            do 100 in = 1, mn
            do 100 il = 0, mm ! S.H. 2-> mm (2022.3.11)
            do ntf = 1, nfile

               if( tr(jr,jz,iz,in,il,iat(1,ntf)) .ne. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(jr,jz,iz,in,il,iat(1,ntf)),
     &                            tr(jr,jz,iz,in,il,iat(2,ntf)),
     &                            cc)

                  tr(jr,jz,iz,in,il,iat(1,ntf)) = Xa
                  tr(jr,jz,iz,in,il,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(jr,jz,iz,in,il,iat(1,ntf)) .gt. cmax )
     &                 cmax = tr(jr,jz,iz,in,il,iat(1,ntf))

                  if( tr(jr,jz,iz,in,il,iat(1,ntf)) .lt. cmin )
     &                 cmin = tr(jr,jz,iz,in,il,iat(1,ntf))

               else

                  if( il .eq. 0 ) isdz = 1
                  tr(jr,jz,iz,in,il,iat(2,ntf)) = 0.0

               end if

            end do
  100       continue

!$OMP PARALLEL
!$OMP& private(ipomp,jr,jz,iz,in,il,ntf,isdz,rdata,answer,rerr)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
            do 101 jr = 1, nr
            do 101 jz = 1, nz
            do 101 iz = 1, mz
            do 101 in = 1, mn
            do 101 il = 0, mm ! S.H. 2-> mm (2022.3.11)

            if( manatally .eq. 0 ) then ! user defined analysis

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(jr,jz,iz,in,il,iat(1,ntf))
                  rdata(2,ntf) = tr(jr,jz,iz,in,il,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(il+1,iz,in,1,icf(jr,jz),1,1) = answer
              anatalrst(il+1,iz,in,1,icf(jr,jz),1,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then ! systematic uncertainty

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(jr,jz,iz,in,il,iat(1,ntf))
                  rdata(2,ntf) = tr(jr,jz,iz,in,il,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(il+1,iz,in,1,icf(jr,jz),1,1) ! mean
     &              = fmval
               if ( fmval .ne. 0d0 ) then
                  anatalrst(il+1,iz,in,1,icf(jr,jz),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(il+1,iz,in,1,icf(jr,jz),1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(il+1,iz,in,1,icf(jr,jz),1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(il+1,iz,in,1,icf(jr,jz),1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(il+1,iz,in,1,icf(jr,jz),1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(il+1,iz,in,1,icf(jr,jz),1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(il+1,iz,in,1,icf(jr,jz),1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(il+1,iz,in,1,icf(jr,jz),1,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then ! c-value dependence

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(il+1,iz,in,1,icf(jr,jz),1,iat(ioe,ntf))
     &              = tr(jr,jz,iz,in,il,iat(ioe,ntf))
                end do
               end do

            end if

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

*-----------------------------------------------------------------------
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        r, z, mass, charge axes
*-----------------------------------------------------------------------

         do iz = 1, mz
            bzz(iz) = iz
         end do
         do in = 1, mn
            bnn(in) = in
         end do
         do jz = 1, nz
         do jr = 1, nr
            vl3d(jr,jz,1) = vl(jr,jz)
         end do
         end do

*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1
     &        .or. itaxs(m,iax) .eq. 5
     &        .or. itaxs(m,iax) .eq. 6
     &        .or. itaxs(m,iax) .eq. 7 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 15   ! mass axis
           case ( 5 )
            iDaxis = 20   ! z axis
           case ( 6 )
            iDaxis = 19   ! r axis
           case ( 7 )
            iDaxis = 17   ! charge axis
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.19

C rearrange data from anatalrst to anataldata
          call anatal_rearrange(0,0,2,
     &         iDaxis,ntaxis,anatalrst,
     &         mm+1,     mz,      mn,    0,      nr, nz,  1,    0, nrst, ! S.H. 3-> mm+1 (2022.3.11)
     &         bzz, 1,  bnn,1,  1, 1,   rm, zm,  dummy,1,
     &         vl3d,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          allocate( tm(nijaxs,2*nfile+3) )

          inum = 0

          if( nn .eq. 0 ) then
             nc = 1
          else
             nc = nn
          end if

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

           ic = ij(2)

*-----------------------------------------------------------------------

           iz = 0
           ia = 0
           in = 0
           if( nn .gt. 0 ) then

              iz = nt(ic) / 1000
              ia = nt(ic) - iz * 1000

              if( ia .gt. 0 ) then
                 in = ia - iz
                 call chname(idum,ia,iz,chau)
              end if

           end if

*-----------------------------------------------------------------------

           do itmp2 = 1, nrst
           do itmp1 = 1, nijaxs
              tm(itmp1,itmp2) = 0.0d+0
           end do
           end do

           seka = 0.0
           sera = 0.0
           seva = 0.0

           do ijaxs = 1, nijaxs ! (1)

            if( itunt(m) .eq. 1 ) then

               vm = 1.0

            else if( itunt(m) .eq. 2 ) then

               vm = delvol(ianataldata,ijaxs)

            end if

*-----------------------------------------------------------------------

            if( nn .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

              vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
              tm(ijaxs,1) = tm(ijaxs,1) + vn
              tm(ijaxs,2) = tm(ijaxs,2)
     &             + ( vn * anataldata(il+1,ianataldata,ijaxs,2) )**2
              if ( manatally .eq. 1 ) then ! systematic uncertainty
                 tm(ijaxs,3) = tm(ijaxs,3)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,3) )**2
                 tm(ijaxs,4) = tm(ijaxs,4)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,4) )**2
              end if

             end do

            else if( ia .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

              vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
              tm(ijaxs,1) = tm(ijaxs,1) + vn
              tm(ijaxs,2) = tm(ijaxs,2)
     &             + ( vn * anataldata(il+1,ianataldata,ijaxs,2) )**2
              if ( manatally .eq. 1 ) then ! systematic uncertainty
                 tm(ijaxs,3) = tm(ijaxs,3)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,3) )**2
                 tm(ijaxs,4) = tm(ijaxs,4)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,4) )**2
              end if

             end do

            else

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

              vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
              tm(ijaxs,1) = tm(ijaxs,1) + vn
              tm(ijaxs,2) = tm(ijaxs,2)
     &             + ( vn * anataldata(il+1,ianataldata,ijaxs,2) )**2
              if ( manatally .eq. 1 ) then ! systematic uncertainty
                 tm(ijaxs,3) = tm(ijaxs,3)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,3) )**2
                 tm(ijaxs,4) = tm(ijaxs,4)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,4) )**2
              end if
             end do

            end if

            seka = seka + tm(ijaxs,1)
            sera = sera + tm(ijaxs,2)
            seva = seva + vm

            if( tm(ijaxs,1) .gt. 0.0 ) then

             tm(ijaxs,2) = sqrt( tm(ijaxs,2) ) / tm(ijaxs,1)
             tm(ijaxs,1) = tm(ijaxs,1) / vm
             if ( manatally .eq. 1 ) then ! systematic uncertainty
                tm(ijaxs,3) = sqrt( tm(ijaxs,3) ) / tm(ijaxs,1)
                tm(ijaxs,4) = sqrt( tm(ijaxs,4) ) / tm(ijaxs,1)
             end if

            end if

           end do                 ! ijaxs = 1, nijaxs (1)

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
               write(iot,'( "y: Number ",a15)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

           write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

           if( itanl(m) .gt. 0 ) then
              write(iot,'( "p: ",200a1)')
     &             ( itang(m)(i:i),i = 1, itanl(m) )
           end if

           if( itsans(m) .gt. 0 ) then
              call write_sangel(iot,m,0)
           end if

           if ( manatally .eq. 0 ) then ! user defined analysis
            if( itaxs(m,iax) .eq. 1 ) then ! mass axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5    y,hl0       n")')

             else
                write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

             end if

             write(iot,'( "#  mass     number      r.err")')

            else if( itaxs(m,iax).eq.5 .or. itaxs(m,iax).eq.6 ) then ! z,r axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: n",12x,"x",12x,
     &                    "y,hhl0      n")')

             else
                write(iot,'( "h: n",12x,"x",12x,
     &                    "y1,hhl0     ny2 dy1=[y1*y2]")')

             end if

             write(iot,'( "#  lower        upper  ",3x,
     &                    "   number      r.err")')

            else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5    y,hl0       n")')

             else
                write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

             end if

             write(iot,'( "#  charge   number      r.err")')

            end if

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if( itaxs(m,iax) .eq. 1 ) then ! mass axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5   y,hl0       n n n")')

             else
                write(iot,'( "h: x-0.5   y1,hl0      ",
     &                      "ny2 dy1=[y1*y2] n n")')

             end if

             write(iot,'( "#  mass    number      "
     &                   ,"r.err(tot, syst, stat)")')

            else if( itaxs(m,iax).eq.5 .or. itaxs(m,iax).eq.6 ) then ! r,z axis
             if ( iteps(m) .ne. 2 ) then
              write(iot,'( "h: n",12x,"x",12x,
     &                    "y,hhl0      n n n")')

             else
              write(iot,'( "h: n",12x,"x",12x,
     &                    "y1,hhl0     ny2 dy1=[y1*y2] n n")')

             end if

             write(iot,'( "#  lower        upper  ",3x,
     &                    "   number      r.err(tot, syst, stat)")')

            else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5   y,hl0       n n n")')

             else
               write(iot,'( "h: x-0.5   y1,hl0      ",
     &                      "ny2 dy1=[y1*y2] n n")')

             end if

             write(iot,'( "#charge    number      "
     &                   ,"r.err(tot, syst, stat)")')

            end if

           else if ( manatally .eq. 2 ) then ! c-value dependence
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h:   x",10x,
     &                    "y,l3        n")')

            else
               write(iot,'( "h:   x",10x,
     &                    "y1,l3       ny2 dy1=[y1*y2]")')

            end if

            write(iot,'( "#  c-value  ",
     &                    "   number      r.err")')

           end if

*-----------------------------------------------------------------------

           do ijaxs = 1, nijaxs
              if( tm(ijaxs,1) .gt. 0.0d0 ) goto 201
           end do

 201       im = ijaxs

           do ijaxs = nijaxs, im + 1, -1
              if( tm(ijaxs,1) .gt. 0.0d0 ) goto 211
           end do

           jm = im
 211       jm = ijaxs

           im = max( 1, im - 2 )
           jm = jm + 2

*-----------------------------------------------------------------------

           if( itaxs(m,iax).eq.1 .or. itaxs(m,iax).eq.7 ) then ! mass, charge axis
            if( im .eq. 1 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(3x,f4.1,1x,1pe13.4,0p1f8.4)')
     &                                   0.5, 0.0, 0.0

             else if ( manatally .eq. 1 ) then ! systematic uncertainty
                write(iot,'(3x,f4.1,1x,1pe13.4,0p3f8.4)')
     &                                   0.5, 0.0, 0.0, 0.0, 0.0

             else if ( manatally .eq. 2 ) then ! c-value dependence

             end if

            end if
           end if

           do ijaxs = 1, nijaxs ! (2)

            if ( manatally .eq. 0 ) then ! user defined analysis
             if( itaxs(m,iax) .eq. 1 ) then ! mass axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p1f8.4)')
     &                ijaxs+iz, (tm(ijaxs,irst),irst=1,2)

              end if

             else if( itaxs(m,iax).eq.5 .or. itaxs(m,iax).eq.6 ) then ! z,r axis
                write(iot,'(1p2e13.4,1pe13.4,0p1f8.4)')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               (tm(ijaxs,irst),irst=1,2)

             else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p1f8.4)')
     &                ijaxs, (tm(ijaxs,irst),irst=1,2)

              end if

             end if

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
             if( itaxs(m,iax) .eq. 1 ) then ! mass axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p3f8.4)')
     &                ijaxs+iz, (tm(ijaxs,irst),irst=1,4)

              end if

             else if( itaxs(m,iax).eq.5 .or. itaxs(m,iax).eq.6 ) then ! z,r axis
                write(iot,'(1p2e13.4,1pe13.4,0p3f8.4)')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               (tm(ijaxs,irst),irst=1,4)

             else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                write(iot,'(3x,i4,1x,1pe13.4,0p3f8.4)')
     &               ijaxs, (tm(ijaxs,irst),irst=1,4)
              end if

             end if

            else if ( manatally .eq. 2 ) then ! c-value dependence
               write(iot,'(3x,i4,1x,1pe13.4,0p1f8.4)')
     &              ijaxs, (tm(ijaxs,irst),irst=1,2)


            end if

           end do                 ! ijaxs = 1, nijaxs (2)

           if( seka .gt. 0.0 ) then
              sera = sqrt( sera ) / seka
              if( itunt(m) .eq. 2 ) seka = seka / seva
           end if

           if ( manatally .eq. 0 ) then ! user defined analysis

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if( itaxs(m,iax) .eq. 1 ) then ! mass axis
             write(iot,'(/"#   sum  ",1pe13.4,0pf8.4)') seka,sera

             if( itunt(m) .eq. 2 ) then
                write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm
             end if

            else if( itaxs(m,iax).eq.5 .or. itaxs(m,iax).eq.6 ) then ! z,r axis
               write(iot,'(/"#   sum over",14x,1pe13.4,0pf8.4)')
     &               seka,sera

               if( itunt(m) .eq. 2 ) then
                  write(iot,'(
     &           "#   vol sum = ",1p1e13.4,"  [cm^3]")') seva
               end if

            else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
               write(iot,'(/"#   sum  ",1pe13.4,0pf8.4)') seka,sera

               if( itunt(m) .eq. 2 ) then
                  write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm
               end if

            end if

           else if ( manatally .eq. 2 ) then ! c-value dependence

           end if

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
              if ( ijtmp .eq. 2 ) then
                 if ( iz.gt.0)
     &                write(iot,'(a15,"&=&",i5,3x,"(",a4,")")')
     &                cij(ijtmp), iz, elmnt(iz)
              else if ( ijtmp .eq. 3 ) then
                 if ( ia .gt. 0 )
     &                write(iot,'(a15,"&=&",i5)') cij(ijtmp), ia
              else if ( ijtmp .eq. 5 ) then
                 write(iot,'(a15,"&=&",i5)')
     &                cij(ijtmp), idnint(fg(itmp))
                 write(iot,'(3x,"volume",6x,"&=&",1pe13.4," [cm^3]")')
     &                delvol(ianataldata,1)
              else
                 write(iot,'(a15,"&=&",i5)')
     &                cij(ijtmp), idnint(fg(itmp))
              end if
              itmp = itmp - ij(ijtmp) + nij(ijtmp)
             else
                write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                           /1pe13.4,2x,"$--$",1pe13.4)')
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
          deallocate( tm )

*-----------------------------------------------------------------------
*        chart axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 ) then

          iDaxis = 18         ! chart axis

C rearrange data from anatalrst to anataldata
          call anatal_rearrange(0,0,2,
     &         iDaxis,ntaxis,anatalrst,
     &         mm+1,     mz,      mn,    0,      nr, nz,  1,    0, nrst, ! S.H. 3-> mm+1 (2022.3.11)
     &         bzz, 1,  bnn,1,  1, 1,   rm, zm,  dummy,1,
     &         vl3d,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          inum = 0

*-----------------------------------------------------------------------

          do ianataldata = 1, nanataldata

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

           do il = 0, mm         ! loop for isomeric level ! S.H. 2-> mm (2022.3.11)

*-----------------------------------------------------------------------

            itmax  = 0
            icmax  = 0
            inmax  = 0

            do iz = 1, maxpt
             do in = 1, maxnt

              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
                 ijaxs = in + (iz-1)*maxnt
              else
                 ijaxs = igetiznmp(iz,in,il,m)
              end if
              if( anataldata(il+1,ianataldata,ijaxs,1).gt.0d0 ) then
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

            if( itmax .eq. 0 ) exit

*-----------------------------------------------------------------------

            write(iot,'(/"#",78("-"))')

            inum = inum + 1

            if( inum .eq. 1 .and. il .eq. 0 ) then
               write(iot,'( "#newpage:")')
            else
               write(iot,'( " newpage:")')
            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &              yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else
               if( itaxl(m) .eq. 0 ) then
                  write(iot,'(/"x: ",a21)') cijaxs
               else
                  write(iot,'(/"x: ",200a1)')
     &                 (itaxt(m)(i:i),i=1,itaxl(m))
               end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
               if( itayl(m) .eq. 0 ) then
                  write(iot,'(/"y: ",a21)') cijaxs2
               else
                  write(iot,'( "y: ",200a1)')
     &                 (itayt(m)(i:i),i=1,itayl(m))
               end if
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
     &              icmax+2, inmax+2

             do ijaxs2 = icmax+2, 1, -1
              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
              write(iot,'(1p10e11.3)')
     &       ( anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &               ijaxs = 1, inmax+2 )

              else
               write(iot,'(1p10e11.3)')
     &                ( anataldata(il+1,ianataldata,
     &                igetiznmp(ijaxs2,ijaxs,il,m),ioe),
     &               ijaxs = 1, inmax+2 )
              end if
             end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

             write(iot,'(/"#   Z    N  Mass  ",
     &                      "  number    r.err")')

             do ijaxs2 = 1, icmax
              do ijaxs = 1, inmax

              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
               if( anataldata(il+1,ianataldata
     &                ,ijaxs+(ijaxs2-1)*nijaxs,1) .gt. 0.d0 ) then

                  write(iot,'(3i5,1pe13.4,0pf8.4)')
     &                 ijaxs2, ijaxs, ijaxs2+ijaxs,
     &           anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &           anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,2)
               end if

              else
               if( anataldata(il+1,ianataldata
     &                ,igetiznmp(ijaxs2,ijaxs,il,m),1) .gt. 0.d0 ) then

                  write(iot,'(3i5,1pe13.4,0pf8.4)')
     &                 ijaxs2, ijaxs, ijaxs2+ijaxs,
     &      anataldata(il+1,ianataldata,igetiznmp(ijaxs2,ijaxs,il,m),1),
     &      anataldata(il+1,ianataldata,igetiznmp(ijaxs2,ijaxs,il,m),2)

               end if
              end if

              end do
             end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

             write(iot,'(/"#   Z = 1 to ",i3/
     &                      "#   N = 1 to ",i3)')
     &              icmax+2, inmax+2

             write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'Z/N', ( dble( in ), in = 1, inmax+2 )

             do ijaxs2 = icmax+2, 1, -1
              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
                write(iot,'(1p1000e11.3)')
     &               dble( ijaxs2 ),
     &         ( anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &               ijaxs = 1, inmax+2 )

              else
                write(iot,'(1p1000e11.3)')
     &               dble( ijaxs2 ),
     &    ( anataldata(il+1,ianataldata,igetiznmp(ijaxs2,ijaxs,il,m),1),
     &               ijaxs = 1, inmax+2 )
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

      else

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else
        if( itazl(m) .eq. 0 ) then
           write(iot,'("y: Number ",a15)') hsunit(itunt(m))
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
       write(iot,'("wt: s(0.7)",/a1,"vspace{-3}")') yen
       itmp = 0
       do ijtmp = 2,8
        if (cij(ijtmp) .ne. 'F' ) then
         itmp = itmp + ij(ijtmp)
         if ( ibin(ijtmp) .eq. 0 ) then
            if ( ijtmp .eq. 5 ) then
               write(iot,'(a15,"&=&",i5)')
     &              cij(ijtmp), idnint(fg(itmp))
               write(iot,'(3x,"volume")')
               write(iot,'(5x,"&=&",1pe13.4," [cm^3]")')
     &              delvol(ianataldata,1)
            else
               write(iot,'(a15,"&=&",i5)')
     &              cij(ijtmp), idnint(fg(itmp))
            end if
            itmp = itmp - ij(ijtmp) + nij(ijtmp)
           else
              write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                           /1pe13.4,2x,"$--$",1pe13.4)')
     &             cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
              itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
           end if
          end if
         end do
         if(il .eq. 0) then     ! specify isomeric level at the bottom of figure
            write(iot,'("il = 0 (Ground state)}")')
         elseif(il .eq. 1) then
            write(iot,'("il = 1 (1st isomer)}")')
         else
            write(iot,'("il = 2 (2nd isomer)}")')
         endif

         write(iot,'("e:")')

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

           end do                    ! il = 0, 2 ! loop for isomeric level
          end do                    ! ianataldata = 1, nanataldata

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( anataldata )
          deallocate( delvol )

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

          select case( itaxs(m,iax) )
           case ( 12 )
            iDaxis = 35   ! rz
          end select

C rearrange data from anatalrst to anataldata
          call anatal_rearrange(0,0,2,
     &         iDaxis,ntaxis,anatalrst,
     &         mm+1,     mz,      mn,    0,      nr, nz,  1,    0, nrst, ! S.H. 3-> mm+1 (2022.3.11)
     &         bzz, 1,  bnn,1,  1, 1,   rm, zm,  dummy,1,
     &         vl3d,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          allocate( tm(nijaxs*nijaxs2,2*nfile+3) )

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

          do ianataldata = 1, nanataldata
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

           ic = ij(2)

           iz = 0
           ia = 0
           in = 0
           if( nn .gt. 0 ) then
              iz = nt(ic) / 1000
              ia = nt(ic) - iz * 1000
              if( ia .gt. 0 ) in = ia - iz
           end if

           do ijaxs = 1, nijaxs*nijaxs2 ! (1)

            sek = 0.0

            if( itunt(m) .eq. 1 ) then

               vm = 1.0

            else if( itunt(m) .eq. 2 ) then

               vm = delvol(ianataldata,ijaxs)

            end if

*-----------------------------------------------------------------------

            if( nn .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

                vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
                sek = sek + vn

             end do

            else if( ia .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

                vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
                sek = sek + vn

             end do

            else

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

                vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
                sek = sek + vn
             end do

            end if

            fmsv = sek / vm

            if( fmsv .gt. fmax ) fmax = fmsv
            if( fmsv .gt. 0.0 .and.
     &           fmsv .lt. fmin ) fmin = fmsv

           end do                 ! ijaxs = 1, nijaxs*nijaxs2 (1)

          end do

*-----------------------------------------------------------------------

          do ianataldata = 1, nanataldata
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

           ic = ij(2)

           inum = inum + 1

*-----------------------------------------------------------------------

           iz = 0
           ia = 0
           in = 0
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

*-----------------------------------------------------------------------

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

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
                 write(iot,'(/"x: ",a21)') cijaxs
              else
                 write(iot,'(/"x: ",200a1)')
     &                (itaxt(m)(i:i),i=1,itaxl(m))
              end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
              if( itayl(m) .eq. 0 ) then
                 write(iot,'(/"y: ",a21)') cijaxs2
              else
                 write(iot,'( "y: ",200a1)')
     &                (itayt(m)(i:i),i=1,itayl(m))
              end if
           end if

*-----------------------------------------------------------------------

              if( inum .eq. 1 ) then

                  form  = ( fgaxs2(nijaxs2+1) - fgaxs2(1) )
     &                  / ( fgaxs(nijaxs+1) - fgaxs(1) )
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

           do itmp2 = 1, nrst
           do itmp1 = 1, nijaxs*nijaxs2
              tm(itmp1,itmp2) = 0.0d+0
           end do
           end do

           do ijaxs = 1, nijaxs*nijaxs2 ! (2)

            sek = 0.0

            if( itunt(m) .eq. 1 ) then

               vm = 1.0

            else if( itunt(m) .eq. 2 ) then

               vm = delvol(ianataldata,ijaxs)

            end if

*-----------------------------------------------------------------------

            if( nn .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

                vn  = vm * anataldata(il+1,ianataldata,ijaxs,ioe)
                tm(ijaxs,ioe) = tm(ijaxs,ioe) + vn

             end do

            else if( ia .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

                vn  = vm * anataldata(il+1,ianataldata,ijaxs,ioe)
                tm(ijaxs,ioe) = tm(ijaxs,ioe) + vn

             end do

            else

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

                vn  = vm * anataldata(il+1,ianataldata,ijaxs,ioe)
                tm(ijaxs,ioe) = tm(ijaxs,ioe) + vn

             end do

            end if

            tm(ijaxs,ioe) = tm(ijaxs,ioe) / vm

           end do                 ! ijaxs = 1, nijaxs*nijaxs2 (2)

*-----------------------------------------------------------------------

           write(iot,'("#  n2 = ",i3,"   n1 = ",i3)')
     &          nijaxs2, nijaxs

           if( ittwo(m) .ne. 4 ) then

              write(iot,'( "# ( ( data(i1,i2), i1 = 1, n1 ),",
     &                      " i2 = n2, 1, -1 )")')

           end if

*-----------------------------------------------------------------------

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              write(iot,'(/a4," y = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7," ; x = ",
     &         1p1g14.7," to ",1p1g14.7," by ",
     &         1p1g14.7," ;")') dc2,
     &         axs2max - axs2del/2d0, axs2min + axs2del/2d0, axs2del,
     &         axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

              write(iot,'(1p10e11.3)')
     &             ( ( tm(ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &             ijaxs = 1, nijaxs ), ijaxs2 = nijaxs2, 1, -1 )

*-----------------------------------------------------------------------

           else if( ittwo(m) .eq. 4 ) then

              write(iot,'(/"# axis2      axis1    ",
     &                      "  number")')

              do ijaxs = 1, nijaxs
                 do ijaxs2 = 1, nijaxs2

                    write(iot,'(1p10e11.3)')
     &                   fgaxs2(ijaxs2)  + axs2del/2d0,
     &                   fgaxs(ijaxs)  + axs1del/2d0,
     &                   tm(ijaxs+(ijaxs2-1)*nijaxs,ioe)

                 end do
              end do

*-----------------------------------------------------------------------

           else if( ittwo(m) .eq. 5 ) then

              write(iot,'("#  ax2= ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7/
     &                     "#  ax1= ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7)')
     &         axs2min + axs2del/2d0, axs2max - axs2del/2d0, axs2del,
     &         axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

              write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &          'ax2/ax1',( fgaxs(ijaxs) + axs1del/2d0, ijaxs=1,nijaxs )

              do ijaxs2 = nijaxs2, 1, -1

                 write(iot,'(1p1000e11.3)')
     &                fgaxs2(ijaxs2)  + axs2del/2d0,
     &                ( tm(ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                ijaxs = 1, nijaxs )

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

      else

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else
        if( itazl(m) .eq. 0 ) then
           write(iot,'("y: Number ",a15)') hsunit(itunt(m))
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
       write(iot,'("wt: s(0.7)",/a1,"vspace{-3}")') yen
       itmp = 0
       do ijtmp = 2,8
        if (cij(ijtmp) .ne. 'F' ) then
         itmp = itmp + ij(ijtmp)
         if ( ibin(ijtmp) .eq. 0 ) then
          if ( ijtmp .eq. 2 ) then
             if ( iz.gt.0)
     &            write(iot,'(a15,"&=&",i5,3x,"(",a4,")")')
     &            cij(ijtmp), iz, elmnt(iz)
          else if ( ijtmp .eq. 3 ) then
             if ( ia .gt. 0 )
     &            write(iot,'(a15,"&=&",i5)') cij(ijtmp), ia
          else if ( ijtmp .eq. 5 ) then
             write(iot,'(a15,"&=&",i5)')
     &            cij(ijtmp), idnint(fg(itmp))
             write(iot,'(3x,"volume",6x,"&=&",1pe13.4," [cm^3]")')
     &            delvol(ianataldata,1)
          else
             write(iot,'(a15,"&=&",i5)')
     &            cij(ijtmp), idnint(fg(itmp))
          end if
          itmp = itmp - ij(ijtmp) + nij(ijtmp)
         else
            write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                           /1pe13.4,2x,"$--$",1pe13.4)')
     &           cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
            itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
         end if
        end if
       end do

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

          end do                    ! ianataldata = 1, nanataldata

*-----------------------------------------------------------------------

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( fgaxs2 )
          deallocate( fgaxs3 )
          deallocate( anataldata )
          deallocate( delvol )
          deallocate( tm )

*-----------------------------------------------------------------------
*        dchain axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 13 ) then

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
                          t0_tr = anatalrst(il+1,iz,in,1,icf(jr,jz),1,1)
                        else
                           t0_tr =
     &          anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,icf(jr,jz),1,1)
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
                          t0_tr = anatalrst(il+1,iz,in,1,icf(jr,jz),1,1)
                        else
                          t0_tr =
     &          anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,icf(jr,jz),1,1)
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
     &                        jr,jz,
     &            (anatalrst(il+1,iz,i,1,icf(jr,jz),1,ioe), i = n3, n4 )
                      else
                        write(iot,'(2i7,1p12e11.3)')
     &        (anatalrst(0+1,igetiznmp(iz,i,il,m),1,1,icf(jr,jz),1,ioe),
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
     &                           ir,iz,
     &        (anatalrst(il+1,iz,i,1,icf(jr,jz),1,2)*100.0, i = n3, n4 )
                         else
                           write(iou,'(2i7,1p12e11.3)') ir,iz,
     &    (anatalrst(0+1,igetiznmp(iz,i,il,m),1,1,icf(jr,jz),1,2)*100.0,
     &                     i = n3, n4 )
                         end if
                        end do
                        end do

                     end if

                  end do

  170       continue

            if( itout(m) .ne. 0 ) close(iou)

*-----------------------------------------------------------------------

         end if                   ! itaxs(m,iax)

*-----------------------------------------------------------------------

         close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

 900  continue

      end do

*-----------------------------------------------------------------------
      deallocate (anatalrst)
      deallocate (vl3d,fm,rdata )

      return
      end


************************************************************************
*                                                                      *
      subroutine anatal_pyildxyz(m,mz,mn,mm,nf,nl,lt, ! S.H. added mm (2022.3.11)
     &                    nx,ny,nz,nn,xm,ym,zm,nt,ikzz,iknn,
     &                    tr,nfile,weightRate,
     &                    igsh,idasa,manatally)
*                                                                      *
*       anatally calculation and output of t-yield with xyz mesh       *
*                                                                      *
************************************************************************

      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

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
      real(8) bzz(mz+1), bnn(mn+1)
!      dimension   tm(maxnt+maxpt,2*nfile+3)
!      dimension   fm(nf,nf)
      dimension   tr(nx,ny,nz,mz,mn,0:mm,2*nfile) ! S.H. 2-> mm (2022.3.11)

      integer,allocatable :: ixyz(:)
!      dimension   weightRate(nfile), rdata(2,nfile)
      dimension   weightRate(nfile)

      integer irst,nrst, itmpdata
!      dimension   anatalrst(mm+1,mz+1,mn+1,1,nx*ny*nz,
!     &                      1,2*nfile+3) ! S.H. 3-> mm+1 (2022.3.11)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)
      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*21,cijaxs*21,cijaxs2*21
      integer iDaxis
      real(8),allocatable :: tm(:,:),fm(:,:),rdata(:,:)

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

*-----------------------------------------------------------------------
      character(28) cfmt                !FURUTA20200615
      character(12) cir                 !FURUTA20200615
      common /redufmt/ iredufmt(itlmax) !FURUTA20200615
*-----------------------------------------------------------------------
      integer           icf
      integer           iat, iad, iaf

      icf(jx,jy,jz) = jx + ( jy - 1 ) * nx + ( jz - 1 ) * nx * ny
      iat(iad,iaf) = iad + (iaf-1) * 2

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------

      real(8),allocatable :: vl(:,:,:)
      allocate( vl(nx,ny,nz) )
      do jx = 1, nx
      do jy = 1, ny
      do jz = 1, nz
         vl(jx,jy,jz) = vls(nl,lt,itmcn(m),itvm(m),itmtr(m,4),
     &                      xm(jx),xm(jx+1),
     &                      ym(jy),ym(jy+1),
     &                      zm(jz),zm(jz+1))
      end do
      end do
      end do

*-----------------------------------------------------------------------

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

      allocate (fm(nf,nf),rdata(2,nfile))

      nrst = 2*nfile+3
      allocate (anatalrst(mm+1,mz+1,mn+1,1,nx*ny*nz,
     &                      1,nrst))


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

            do 100 jx = 1, nx
            do 100 jy = 1, ny
            do 100 jz = 1, nz

               if( itunt(m) .eq. 1 ) then

                  cc = rtfac(m)

               else if( itunt(m) .eq. 2 ) then

                  cc = rtfac(m) / vl(jx,jy,jz)

               end if

            do 100 iz = 1, mz
            do 100 in = 1, mn
            do 100 il = 0, mm ! S.H. 2-> mm (2022.3.11)
            do ntf = 1, nfile

               if( tr(jx,jy,jz,iz,in,il,iat(1,ntf)) .ne. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(jx,jy,jz,iz,in,il,iat(1,ntf)),
     &                            tr(jx,jy,jz,iz,in,il,iat(2,ntf)),
     &                            cc)

                  tr(jx,jy,jz,iz,in,il,iat(1,ntf)) = Xa
                  tr(jx,jy,jz,iz,in,il,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(jx,jy,jz,iz,in,il,iat(1,ntf)) .gt. cmax )
     &                 cmax = tr(jx,jy,jz,iz,in,il,iat(1,ntf))

                  if( tr(jx,jy,jz,iz,in,il,iat(1,ntf)) .lt. cmin )
     &                 cmin = tr(jx,jy,jz,iz,in,il,iat(1,ntf))

               else

                  if( il .eq. 0 ) isdz = 1
                  tr(jx,jy,jz,iz,in,il,iat(2,ntf)) = 0.0

               end if

            end do
  100       continue

!$OMP PARALLEL
!$OMP& private(ipomp,jx,jy,jz,iz,in,il,ntf,isdz,rdata,answer,rerr)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
           do itmpxyz=1,nx*ny*nz
            jz=(itmpxyz-1)/(ny*nx)+1
            jy=(itmpxyz-1-(jz-1)*ny*nx)/nx+1
            jx=itmpxyz-(jy-1)*nx-(jz-1)*ny*nx
            do 101 iz = 1, mz
            do 101 in = 1, mn
            do 101 il = 0, mm ! S.H. 2-> mm (2022.3.11)

            if( manatally .eq. 0 ) then ! user defined analysis

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(jx,jy,jz,iz,in,il,iat(1,ntf))
                  rdata(2,ntf) = tr(jx,jy,jz,iz,in,il,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,1) = answer
              anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then ! systematic uncertainty

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(jx,jy,jz,iz,in,il,iat(1,ntf))
                  rdata(2,ntf) = tr(jx,jy,jz,iz,in,il,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,1) ! mean
     &              = fmval
               if ( fmval .ne. 0d0 ) then
                  anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then ! c-value dependence

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,iat(ioe,ntf))
     &              = tr(jx,jy,jz,iz,in,il,iat(ioe,ntf))
                end do
               end do

            end if

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

*-----------------------------------------------------------------------
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        x, y, z, mass, charge axes
*-----------------------------------------------------------------------

         do iz = 1, mz
            bzz(iz) = iz
         end do
         do in = 1, mn
            bnn(in) = in
         end do

*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 3
     &        .or. itaxs(m,iax) .eq. 4
     &        .or. itaxs(m,iax) .eq. 5
     &        .or. itaxs(m,iax) .eq. 1
     &        .or. itaxs(m,iax) .eq. 7 ) then

          select case( itaxs(m,iax) )
           case ( 3 )
            iDaxis = 21   ! x axis
           case ( 4 )
            iDaxis = 22   ! y axis
           case ( 5 )
            iDaxis = 23   ! z axis
           case ( 1 )
            iDaxis = 15   ! mass axis
           case ( 7 )
            iDaxis = 17   ! charge axis
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.19

C rearrange data from anatalrst to anataldata
          call anatal_rearrange(0,0,3,
     &         iDaxis,ntaxis,anatalrst,
     &         mm+1,     mz,      mn,    0,      nx, ny, nz,   0, nrst, ! S.H. 3-> mm+1 (2022.3.11)
     &         bzz, 1,  bnn,1,  1, 1,   xm, ym, zm,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          allocate( tm(nijaxs,2*nfile+3) )

          inum = 0

          if( nn .eq. 0 ) then
             nc = 1
          else
             nc = nn
          end if

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

           ic = ij(2)

*-----------------------------------------------------------------------

           iz = 0
           ia = 0
           in = 0
           if( nn .gt. 0 ) then

              iz = nt(ic) / 1000
              ia = nt(ic) - iz * 1000

              if( ia .gt. 0 ) then
                 in = ia - iz
                 call chname(idum,ia,iz,chau)
              end if

           end if

*-----------------------------------------------------------------------

           do itmp2 = 1, nrst
           do itmp1 = 1, nijaxs
              tm(itmp1,itmp2) = 0.0d+0
           end do
           end do

           seka = 0.0
           sera = 0.0
           seva = 0.0

           do ijaxs = 1, nijaxs ! (1)

            if( itunt(m) .eq. 1 ) then

               vm = 1.0

            else if( itunt(m) .eq. 2 ) then

               vm = delvol(ianataldata,ijaxs)

            end if

*-----------------------------------------------------------------------

            if( nn .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

              vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
              tm(ijaxs,1) = tm(ijaxs,1) + vn
              tm(ijaxs,2) = tm(ijaxs,2)
     &             + ( vn * anataldata(il+1,ianataldata,ijaxs,2) )**2
              if ( manatally .eq. 1 ) then ! systematic uncertainty
                 tm(ijaxs,3) = tm(ijaxs,3)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,3) )**2
                 tm(ijaxs,4) = tm(ijaxs,4)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,4) )**2
              end if

             end do

            else if( ia .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

              vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
              tm(ijaxs,1) = tm(ijaxs,1) + vn
              tm(ijaxs,2) = tm(ijaxs,2)
     &             + ( vn * anataldata(il+1,ianataldata,ijaxs,2) )**2
              if ( manatally .eq. 1 ) then ! systematic uncertainty
                 tm(ijaxs,3) = tm(ijaxs,3)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,3) )**2
                 tm(ijaxs,4) = tm(ijaxs,4)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,4) )**2
              end if

             end do

            else

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

              vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
              tm(ijaxs,1) = tm(ijaxs,1) + vn
              tm(ijaxs,2) = tm(ijaxs,2)
     &             + ( vn * anataldata(il+1,ianataldata,ijaxs,2) )**2
              if ( manatally .eq. 1 ) then ! systematic uncertainty
                 tm(ijaxs,3) = tm(ijaxs,3)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,3) )**2
                 tm(ijaxs,4) = tm(ijaxs,4)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,4) )**2
              end if
             end do

            end if

            seka = seka + tm(ijaxs,1)
            sera = sera + tm(ijaxs,2)
            seva = seva + vm

            if( tm(ijaxs,1) .gt. 0.0 ) then

             tm(ijaxs,2) = sqrt( tm(ijaxs,2) ) / tm(ijaxs,1)
             tm(ijaxs,1) = tm(ijaxs,1) / vm
             if ( manatally .eq. 1 ) then ! systematic uncertainty
                tm(ijaxs,3) = sqrt( tm(ijaxs,3) ) / tm(ijaxs,1)
                tm(ijaxs,4) = sqrt( tm(ijaxs,4) ) / tm(ijaxs,1)
             end if

            end if

           end do                 ! ijaxs = 1, nijaxs (1)

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
               write(iot,'( "y: Number ",a15)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

           write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

           if( itanl(m) .gt. 0 ) then
              write(iot,'( "p: ",200a1)')
     &             ( itang(m)(i:i),i = 1, itanl(m) )
           end if

           if( itsans(m) .gt. 0 ) then
              call write_sangel(iot,m,0)
           end if

           if ( manatally .eq. 0 ) then ! user defined analysis
            if( itaxs(m,iax) .eq. 1 ) then ! mass axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5    y,hl0       n")')

             else
                write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

             end if

             write(iot,'( "#  mass     number      r.err")')

            else if( itaxs(m,iax).ge.3 .and. itaxs(m,iax).le.5 ) then ! x,y,z axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: n",12x,"x",12x,
     &                    "y,hhl0      n")')

             else

                write(iot,'( "h: n",12x,"x",12x,
     &                    "y1,hhl0     ny2 dy1=[y1*y2]")')

             end if

             write(iot,'( "#  lower        upper  ",3x,
     &                    "   number      r.err")')

            else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5    y,hl0       n")')

             else
                write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

             end if

             write(iot,'( "#  charge   number      r.err")')

            end if

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if( itaxs(m,iax) .eq. 1 ) then ! mass axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5   y,hl0       n n n")')

             else
                write(iot,'( "h: x-0.5   y1,hl0      ",
     &                      "ny2 dy1=[y1*y2] n n")')

             end if

             write(iot,'( "#  mass    number      "
     &                   ,"r.err(tot, syst, stat)")')

            else if( itaxs(m,iax).ge.3 .and. itaxs(m,iax).le.5 ) then ! x,y,z axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: n",12x,"x",12x,
     &                    "y,hhl0      n n n")')

             else
                write(iot,'( "h: n",12x,"x",12x,
     &                    "y1,hhl0     ny2 dy1=[y1*y2] n n")')

             end if

             write(iot,'( "#  lower        upper  ",3x,
     &                    "   number      r.err(tot, syst, stat)")')

            else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5   y,hl0       n n n")')

             else
                write(iot,'( "h: x-0.5   y1,hl0      ",
     &                      "ny2 dy1=[y1*y2] n n")')

             end if

             write(iot,'( "#charge    number      "
     &                   ,"r.err(tot, syst, stat)")')

            end if

           else if ( manatally .eq. 2 ) then ! c-value dependence
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h:   x",10x,
     &                    "y,l3        n")')

            else
               write(iot,'( "h:   x",10x,
     &                    "y1,l3       ny2 dy1=[y1*y2]")')

            end if

            write(iot,'( "#  c-value  ",
     &                    "   number      r.err")')

           end if

*-----------------------------------------------------------------------

           do ijaxs = 1, nijaxs
              if( tm(ijaxs,1) .gt. 0.0d0 ) goto 201
           end do

 201       im = ijaxs

           do ijaxs = nijaxs, im + 1, -1
              if( tm(ijaxs,1) .gt. 0.0d0 ) goto 211
           end do

           jm = im
 211       jm = ijaxs

           im = max( 1, im - 2 )
           jm = jm + 2

*-----------------------------------------------------------------------

           if( itaxs(m,iax).eq.1 .or. itaxs(m,iax).eq.7 ) then ! mass, charge axis
            if( im .eq. 1 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(3x,f4.1,1x,1pe13.4,0p1f8.4)')
     &                                   0.5, 0.0, 0.0

             else if ( manatally .eq. 1 ) then ! systematic uncertainty
                write(iot,'(3x,f4.1,1x,1pe13.4,0p3f8.4)')
     &                                   0.5, 0.0, 0.0, 0.0, 0.0

             else if ( manatally .eq. 2 ) then ! c-value dependence
                write(iot,'(3x,f4.1,1x,1pe13.4,0p1f8.4)')
     &                                   0.5, 0.0, 0.0

             end if

            end if
           end if

           do ijaxs = 1, nijaxs ! (2)
            if ( manatally .eq. 0 ) then ! user defined analysis
             if( itaxs(m,iax) .eq. 1 ) then ! mass axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p1f8.4)')
     &                ijaxs+iz, (tm(ijaxs,irst),irst=1,2)
              end if

             else if( itaxs(m,iax).ge.3 .and. itaxs(m,iax).le.5 ) then ! x,y,z axis
                write(iot,'(1p2e13.4,1pe13.4,0p1f8.4)')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               (tm(ijaxs,irst),irst=1,2)

             else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p1f8.4)')
     &                ijaxs, (tm(ijaxs,irst),irst=1,2)
              end if

             end if

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
             if( itaxs(m,iax) .eq. 1 ) then ! mass axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                write(iot,'(3x,i4,1x,1pe13.4,0p3f8.4)')
     &               ijaxs+iz, (tm(ijaxs,irst),irst=1,4)
              end if

             else if( itaxs(m,iax).ge.3 .and. itaxs(m,iax).le.5 ) then ! x,y,z axis
                write(iot,'(1p2e13.4,1pe13.4,0p3f8.4)')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               (tm(ijaxs,irst),irst=1,4)

             else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p3f8.4)')
     &                ijaxs, (tm(ijaxs,irst),irst=1,4)
              end if

             end if

            else if ( manatally .eq. 2 ) then ! c-value dependence
               write(iot,'(3x,i4,1x,1pe13.4,0p1f8.4)')
     &              ijaxs, (tm(ijaxs,irst),irst=1,2)

            end if

           end do                 ! ijaxs = 1, nijaxs (2)

           if( seka .gt. 0.0 ) then
              sera = sqrt( sera ) / seka
              if( itunt(m) .eq. 2 ) seka = seka / seva
           end if


           if ( manatally .eq. 0 ) then ! user defined analysis

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if( itaxs(m,iax) .eq. 1 ) then ! mass axis
               write(iot,'(/"#   sum  ",1pe13.4,0pf8.4)') seka,sera

               if( itunt(m) .eq. 2 ) then
                  write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm
               end if

            else if( itaxs(m,iax).ge.3 .and. itaxs(m,iax).le.5 ) then ! x,y,z axis
               write(iot,'(/"#   sum over",14x,1pe13.4,0pf8.4)') 
     &               seka,sera

               if( itunt(m) .eq. 2 ) then
                  write(iot,'(
     &           "#   vol sum = ",1p1e13.4,"  [cm^3]")') seva
               end if

            else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
               write(iot,'(/"#   sum  ",1pe13.4,0pf8.4)') seka,sera

               if( itunt(m) .eq. 2 ) then
                  write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm
               end if

            end if

           else if ( manatally .eq. 2 ) then ! c-value dependence

           end if

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
              if ( ijtmp .eq. 2 ) then
                 if ( iz.gt.0)
     &                write(iot,'(a15,"&=&",i5,3x,"(",a4,")")')
     &                cij(ijtmp), iz, elmnt(iz)
              else if ( ijtmp .eq. 3 ) then
                 if ( ia .gt. 0 )
     &                write(iot,'(a15,"&=&",i5)') cij(ijtmp), ia
              else if ( ijtmp .eq. 5 ) then
                 write(iot,'(a15,"&=&",i5)')
     &                cij(ijtmp), idnint(fg(itmp))
                 write(iot,'(3x,"volume",6x,"&=&",1pe13.4," [cm^3]")')
     &                delvol(ianataldata,1)
              else
                 write(iot,'(a15,"&=&",i5)')
     &                cij(ijtmp), idnint(fg(itmp))
              end if
              itmp = itmp - ij(ijtmp) + nij(ijtmp)
             else
                write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                           /1pe13.4,2x,"$--$",1pe13.4)')
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
          deallocate( tm )

*-----------------------------------------------------------------------
*        chart axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 ) then

          iDaxis = 18         ! chart axis

C rearrange data from anatalrst to anataldata
          call anatal_rearrange(0,0,3,
     &         iDaxis,ntaxis,anatalrst,
     &         mm+1,     mz,      mn,    0,      nx, ny, nz,   0, nrst, ! S.H. 3-> mm+1 (2022.3.11)
     &         bzz, 1,  bnn,1,  1, 1,   xm, ym, zm,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          inum = 0

          do ianataldata = 1, nanataldata

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

           do il = 0, mm         ! loop for isomeric level ! S.H. 2-> mm (2022.3.11)

*-----------------------------------------------------------------------

            itmax  = 0
            icmax  = 0
            inmax  = 0

            do iz = 1, maxpt
             do in = 1, maxnt

              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
                 ijaxs = in + (iz-1)*maxnt
              else
                 ijaxs = igetiznmp(iz,in,il,m)
              end if
              if( anataldata(il+1,ianataldata,ijaxs,1).gt.0d0 ) then
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

            if( itmax .eq. 0 ) exit

*-----------------------------------------------------------------------

            write(iot,'(/"#",78("-"))')

            inum = inum + 1

            if( inum .eq. 1 .and. il .eq. 0 ) then
               write(iot,'( "#newpage:")')
            else
               write(iot,'( " newpage:")')
            end if

*-----------------------------------------------------------------------

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
                  write(iot,'(/"x: ",a21)') cijaxs
               else
                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))
               end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
               if( itayl(m) .eq. 0 ) then
                  write(iot,'(/"y: ",a21)') cijaxs2
               else
                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))
               end if
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

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

             write(iot,'(/"hc: y = ",i3," to 1 by -1 ;",
     &                         " x = 1 to ",i3," by 1 ;")')
     &              icmax+2, inmax+2

             do ijaxs2 = icmax+2, 1, -1
              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
              write(iot,'(1p10e11.3)')
     &       ( anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &               ijaxs = 1, inmax+2 )

              else
               write(iot,'(1p10e11.3)')
     &                ( anataldata(il+1,ianataldata,
     &                igetiznmp(ijaxs2,ijaxs,il,m),ioe),
     &               ijaxs = 1, inmax+2 )
              end if
             end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

             write(iot,'(/"#   Z    N  Mass  ",
     &                      "  number    r.err")')

             do ijaxs2 = 1, icmax
              do ijaxs = 1, inmax

              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
               if( anataldata(il+1,ianataldata
     &                ,ijaxs+(ijaxs2-1)*nijaxs,1) .gt. 0.d0 ) then

                  write(iot,'(3i5,1pe13.4,0pf8.4)')
     &                 ijaxs2, ijaxs, ijaxs2+ijaxs,
     &           anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &           anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,2)

               end if

              else
               if( anataldata(il+1,ianataldata
     &                ,igetiznmp(ijaxs2,ijaxs,il,m),1) .gt. 0.d0 ) then

                  write(iot,'(3i5,1pe13.4,0pf8.4)')
     &                 ijaxs2, ijaxs, ijaxs2+ijaxs,
     &      anataldata(il+1,ianataldata,igetiznmp(ijaxs2,ijaxs,il,m),1),
     &      anataldata(il+1,ianataldata,igetiznmp(ijaxs2,ijaxs,il,m),2)

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

             do ijaxs2 = icmax+2, 1, -1

              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
                write(iot,'(1p1000e11.3)')
     &               dble( ijaxs2 ),
     &         ( anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &            ijaxs = 1, inmax+2 )

              else
                write(iot,'(1p1000e11.3)')
     &               dble( ijaxs2 ),
     &    ( anataldata(il+1,ianataldata,igetiznmp(ijaxs2,ijaxs,il,m),1),
     &               ijaxs = 1, inmax+2 )
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

      else

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else
        if( itazl(m) .eq. 0 ) then
           write(iot,'("y: Number ",a15)') hsunit(itunt(m))
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
       write(iot,'("wt: s(0.7)",/a1,"vspace{-3}")') yen
       itmp = 0
       do ijtmp = 2,8
        if (cij(ijtmp) .ne. 'F' ) then
         itmp = itmp + ij(ijtmp)
         if ( ibin(ijtmp) .eq. 0 ) then
            if ( ijtmp .eq. 5 ) then
               write(iot,'(a15,"&=&",i5)')
     &              cij(ijtmp), idnint(fg(itmp))
               write(iot,'(3x,"volume")')
               write(iot,'(5x,"&=&",1pe13.4," [cm^3]")')
     &              delvol(ianataldata,1)
            else
               write(iot,'(a15,"&=&",i5)')
     &              cij(ijtmp), idnint(fg(itmp))
            end if
            itmp = itmp - ij(ijtmp) + nij(ijtmp)
         else
            write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                           /1pe13.4,2x,"$--$",1pe13.4)')
     &           cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
            itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
         end if
        end if
       end do
       if(il .eq. 0) then       ! specify isomeric level at the bottom of figure
          write(iot,'("il = 0 (Ground state)}")')
       elseif(il .eq. 1) then
          write(iot,'("il = 1 (1st isomer)}")')
       else
          write(iot,'("il = 2 (2nd isomer)}")')
       endif

       write(iot,'("e:")')

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

           end do                    ! il = 0, 2 ! loop for isomeric level
          end do                    ! ianataldata = 1, nanataldata

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( anataldata )
          deallocate( delvol )

*-----------------------------------------------------------------------
*        xy, yz, xz axes (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9
     &         .or. itaxs(m,iax) .eq. 10
     &         .or. itaxs(m,iax) .eq. 11 ) then


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
           case ( 9 )
            iDaxis = 36   ! xy
           case ( 10 )
            iDaxis = 37   ! yz
           case ( 11 )
            iDaxis = 38   ! xz axes
          end select

C rearrange data from anatalrst to anataldata
          call anatal_rearrange(0,0,3,
     &         iDaxis,ntaxis,anatalrst,
     &         mm+1,     mz,      mn,    0,      nx, ny, nz,   0, nrst, ! S.H. 3-> mm+1 (2022.3.11)
     &         bzz, 1,  bnn,1,  1, 1,   xm, ym, zm,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          allocate( tm(nijaxs*nijaxs2,2*nfile+3) )

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

          do ianataldata = 1, nanataldata
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

           ic = ij(2)

           iz = 0
           ia = 0
           in = 0
           if( nn .gt. 0 ) then
              iz = nt(ic) / 1000
              ia = nt(ic) - iz * 1000
              if( ia .gt. 0 ) in = ia - iz
           end if

           do ijaxs = 1, nijaxs*nijaxs2 ! (1)

            sek = 0.0

            if( itunt(m) .eq. 1 ) then

               vm = 1.0

            else if( itunt(m) .eq. 2 ) then

               vm = delvol(ianataldata,ijaxs)

            end if

*-----------------------------------------------------------------------

            if( nn .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

                vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
                sek = sek + vn

             end do

            else if( ia .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

                vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
                sek = sek + vn

             end do

            else

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)
                vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
                sek = sek + vn
             end do

            end if

            fmsv = sek / vm

            if( fmsv .gt. fmax ) fmax = fmsv
            if( fmsv .gt. 0.0 .and.
     &           fmsv .lt. fmin ) fmin = fmsv

           end do                 ! ijaxs = 1, nijaxs*nijaxs2 (1)

          end do

*-----------------------------------------------------------------------

          do ianataldata = 1, nanataldata
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

           ic = ij(2)

           inum = inum + 1

*-----------------------------------------------------------------------

           iz = 0
           ia = 0
           in = 0
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

*-----------------------------------------------------------------------

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

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
                 write(iot,'(/"x: ",a21)') cijaxs
              else
                 write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))
              end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
              if( itayl(m) .eq. 0 ) then
                 write(iot,'(/"y: ",a21)') cijaxs2
              else
                 write(iot,'( "y: ",200a1)')
     &                (itayt(m)(i:i),i=1,itayl(m))
              end if
           end if

*-----------------------------------------------------------------------

              if( inum .eq. 1 ) then

                  form  = ( fgaxs2(nijaxs2+1) - fgaxs2(1) )
     &                  / ( fgaxs(nijaxs+1) - fgaxs(1) )
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

            do itmp2 = 1, nrst
            do itmp1 = 1, nijaxs*nijaxs2
               tm(itmp1,itmp2) = 0.0d+0
            end do
            end do

            do ijaxs = 1, nijaxs*nijaxs2 ! (2)

             sek = 0.0

             if( itunt(m) .eq. 1 ) then

                vm = 1.0

             else if( itunt(m) .eq. 2 ) then

                vm = delvol(ianataldata,ijaxs)

             end if

*-----------------------------------------------------------------------

             if( nn .eq. 0 ) then

              do il = 0, mm ! S.H. 2-> mm (2022.3.11)

                 vn  = vm * anataldata(il+1,ianataldata,ijaxs,ioe)
                 tm(ijaxs,ioe) = tm(ijaxs,ioe) + vn

              end do

             else if( ia .eq. 0 ) then

              do il = 0, mm ! S.H. 2-> mm (2022.3.11)

                 vn  = vm * anataldata(il+1,ianataldata,ijaxs,ioe)
                 tm(ijaxs,ioe) = tm(ijaxs,ioe) + vn

              end do

             else

              do il = 0, mm ! S.H. 2-> mm (2022.3.11)

                 vn  = vm * anataldata(il+1,ianataldata,ijaxs,ioe)
                 tm(ijaxs,ioe) = tm(ijaxs,ioe) + vn

              end do

             end if

             tm(ijaxs,ioe) = tm(ijaxs,ioe) / vm

            end do                ! ijaxs = 1, nijaxs*nijaxs2 (2)

*-----------------------------------------------------------------------

            write(iot,'("#  n2 = ",i3,"   n1 = ",i3)')
     &           nijaxs2, nijaxs

            if( ittwo(m) .ne. 4 ) then

               write(iot,'( "# ( ( data(i1,i2), i1 = 1, n1 ),",
     &                      " i2 = n2, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

               write(iot,'(/a4," y = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7," ; x = ",
     &         1p1g14.7," to ",1p1g14.7," by ",
     &         1p1g14.7," ;")') dc2,
     &         axs2max - axs2del/2d0, axs2min + axs2del/2d0, axs2del,
     &         axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

               write(iot,'(1p10e11.3)')
     &       ( ( tm(ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                ijaxs = 1, nijaxs ), ijaxs2 = nijaxs2, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# axis2      axis1    ",
     &                      "  number")')

               do ijaxs = 1, nijaxs
                  do ijaxs2 = 1, nijaxs2

                     write(iot,'(1p10e11.3)')
     &                    fgaxs2(ijaxs2)  + axs2del/2d0,
     &                    fgaxs(ijaxs)  + axs1del/2d0,
     &                    tm(ijaxs+(ijaxs2-1)*nijaxs,ioe)

                  end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

               write(iot,'("#  ax2= ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7/
     &                     "#  ax1= ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7)')
     &         axs2min + axs2del/2d0, axs2max - axs2del/2d0, axs2del,
     &         axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

               write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &          'ax2/ax1',( fgaxs(ijaxs) + axs1del/2d0, ijaxs=1,nijaxs )

               do ijaxs2 = nijaxs2, 1, -1

                  write(iot,'(1p1000e11.3)')
     &                 fgaxs2(ijaxs2)  + axs2del/2d0,
     &                 ( tm(ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                 ijaxs = 1, nijaxs )

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

      else

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else
        if( itazl(m) .eq. 0 ) then
           write(iot,'("y: Number ",a15)') hsunit(itunt(m))
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
       write(iot,'("wt: s(0.7)",/a1,"vspace{-3}")') yen
       itmp = 0
       do ijtmp = 2,8
        if (cij(ijtmp) .ne. 'F' ) then
         itmp = itmp + ij(ijtmp)
         if ( ibin(ijtmp) .eq. 0 ) then
          if ( ijtmp .eq. 2 ) then
             if ( iz.gt.0)
     &            write(iot,'(a15,"&=&",i5,3x,"(",a4,")")')
     &            cij(ijtmp), iz, elmnt(iz)
          else if ( ijtmp .eq. 3 ) then
             if ( ia .gt. 0 )
     &            write(iot,'(a15,"&=&",i5)') cij(ijtmp), ia
          else if ( ijtmp .eq. 5 ) then
             write(iot,'(a15,"&=&",i5)')
     &            cij(ijtmp), idnint(fg(itmp))
             write(iot,'(3x,"volume",6x,"&=&",1pe13.4," [cm^3]")')
     &            delvol(ianataldata,1)
          else
             write(iot,'(a15,"&=&",i5)')
     &            cij(ijtmp), idnint(fg(itmp))
          end if
          itmp = itmp - ij(ijtmp) + nij(ijtmp)
         else
            write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                           /1pe13.4,2x,"$--$",1pe13.4)')
     &           cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
            itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
         end if
        end if
       end do

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

          end do                    ! ianataldata = 1, nanataldata

*-----------------------------------------------------------------------

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 2
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xy (x-y)
                  bmpfIType = (/ 'n', 'z' /)
               else if ( itaxs(m,iax) .eq. 10 ) then ! when axis = yz (z-y)
                  bmpfIType = (/ 'n', 'x' /)
               else if ( itaxs(m,iax) .eq. 11 ) then ! when axis = xz (z-x)
                  bmpfIType = (/ 'n', 'y' /)
               end if
               bmpfIndex = (/ nc, nijaxs3 /)
               bmpWidth  = nijaxs
               bmpHeight = nijaxs2

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

               if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xy (x-y)
                  write(isunit_vtk_meta) nijaxs,nijaxs2,nijaxs3
                  write(isunit_vtk_meta)
     &                 ( fgaxs(ijaxs), ijaxs=1,nijaxs+1 )
                  write(isunit_vtk_meta)
     &                 ( fgaxs2(ijaxs2), ijaxs2=1,nijaxs2+1 )
                  write(isunit_vtk_meta)
     &                 ( fgaxs3(ijaxs3), ijaxs3=1,nijaxs3+1 )
               else if ( itaxs(m,iax) .eq. 10 ) then ! when axis = yz (z-y)
                  write(isunit_vtk_meta) nijaxs3,nijaxs2,nijaxs
                  write(isunit_vtk_meta)
     &                 ( fgaxs3(ijaxs3), ijaxs3=1,nijaxs3+1 )
                  write(isunit_vtk_meta)
     &                 ( fgaxs2(ijaxs2), ijaxs2=1,nijaxs2+1 )
                  write(isunit_vtk_meta)
     &                 ( fgaxs(ijaxs), ijaxs=1,nijaxs+1 )
               else if ( itaxs(m,iax) .eq. 11 ) then ! when axis = xz (z-x)
                  write(isunit_vtk_meta) nijaxs2,nijaxs3,nijaxs
                  write(isunit_vtk_meta)
     &                 ( fgaxs2(ijaxs2), ijaxs2=1,nijaxs2+1 )
                  write(isunit_vtk_meta)
     &                 ( fgaxs3(ijaxs3), ijaxs3=1,nijaxs3+1 )
                  write(isunit_vtk_meta)
     &                 ( fgaxs(ijaxs), ijaxs=1,nijaxs+1 )
               end if

               do ianataldata = 1, nanataldata
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

                  ic = ij(2)

                  write(isunit_vtk_meta) ic

                   if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xy (x-y)
                     write(isunit_vtk)
     &                    ( ( ( anataldata(1,
     &                    ic+(ijaxs3-1)*nc,
     &                    ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                             ijaxs = 1, nijaxs ),
     &                             ijaxs2 = 1, nijaxs2 ),
     &                             ijaxs3 = 1, nijaxs3 )
                   else if ( itaxs(m,iax) .eq. 10 ) then ! when axis = yz (z-y)
                     write(isunit_vtk)
     &                    ( ( ( anataldata(1,
     &                    ic+(ijaxs3-1)*nc,
     &                    ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                             ijaxs3 = 1, nijaxs3 ),
     &                             ijaxs2 = 1, nijaxs2 ),
     &                             ijaxs = 1, nijaxs )
                   else if ( itaxs(m,iax) .eq. 11 ) then ! when axis = xz (z-x)
                     write(isunit_vtk)
     &                    ( ( ( anataldata(1,
     &                    ic+(ijaxs3-1)*nc,
     &                    ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                             ijaxs2 = 1, nijaxs2 ),
     &                             ijaxs3 = 1, nijaxs3 ),
     &                             ijaxs = 1, nijaxs )
                   end if
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

             if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xy (x-y)
               call vtk_set_cell_regmat(
     &                 isunit_vtk_rm,
     &                 isunit_vtk_geom, isunit_vtk_geom_meta,
     &                 nijaxs, nijaxs2, nijaxs3, fgaxs, fgaxs2, fgaxs3,
     &                 iaxs,iuni,ires, igser,
     &                 1,1,krr,vll,itmtr(m,4))
             else if ( itaxs(m,iax) .eq. 10 ) then ! when axis = yz (z-y)
               call vtk_set_cell_regmat(
     &                 isunit_vtk_rm,
     &                 isunit_vtk_geom, isunit_vtk_geom_meta,
     &                 nijaxs3, nijaxs2, nijaxs, fgaxs3, fgaxs2, fgaxs,
     &                 iaxs,iuni,ires, igser,
     &                 1,1,krr,vll,itmtr(m,4))
             else if ( itaxs(m,iax) .eq. 11 ) then ! when axis = xz (z-x)
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
            deallocate( tm )

*-----------------------------------------------------------------------
*        dchain axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 13 ) then

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
                       t0_tr = anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,1)
                        else
                          t0_tr =
     &       anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,icf(jx,jy,jz),1,1)
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
                       t0_tr = anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,1)
                        else
                          t0_tr =
     &       anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,icf(jx,jy,jz),1,1)
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
     &                        jx,jy,jz,
     &         (anatalrst(il+1,iz,i,1,icf(jx,jy,jz),1,ioe), i = n3, n4 )
                      else
                        write(iot,'(3i7,1p12e11.3)') jx,jy,jz,
     &  (anatalrst(0+1,igetiznmp(iz,i,il,m),1,1,icf(jx,jy,jz),1,ioe),
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
     &                           jx,jy,jz,
     &     (anatalrst(il+1,iz,i,1,icf(jx,jy,jz),1,2)*100.0, i = n3, n4 )
                         else
                           write(iou,'(3i7,1p12e11.3)') jx,jy,jz,
     &       (anatalrst(0+1,igetiznmp(iz,i,il,m),1,1,icf(jx,jy,jz),1,2)
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
                  if(anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,ioe)
     &                    .ne.0d0)then
                   write(iot,cfmt)ir,iz*10000+(in+iz)*10+il,
     &          (anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,ioee),ioee=1,2)
                   if(itout(m).ne.0)then
                    write(iou,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                   anatalrst(il+1,iz,in,1,icf(jx,jy,jz),1,2)
                   endif
                  endif
                 else
                  if(anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,
     &             icf(jx,jy,jz),1,ioe)
     &               .ne.0d0)then
                   write(iot,cfmt)ir,iz*10000+(in+iz)*10+il,
     &              (anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,
     &              icf(jx,jy,jz),1,ioee),
     &              ioee=1,2)
                   if(itout(m).ne.0)then
                    write(iou,cfmt)ir,iz*10000+(in+iz)*10+il,
     &              anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,
     &              icf(jx,jy,jz),1,2)
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

*-----------------------------------------------------------------------

         end if                   ! itaxs(m,iax)

*-----------------------------------------------------------------------

         close(iot)

*-----------------------------------------------------------------------

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      deallocate( vl )
      deallocate( fm,rdata )
      
      return
      end


************************************************************************
*                                                                      *
      subroutine anatal_pyildtet(m,mz,mn,mm, ! S.H. added mm (2022.3.11)
     &                    nr,mr,nn,nt,ikzz,iknn,
     &                    tr,nfile,weightRate,
     &                    nx,ny,nz,kr,xm,ym,zm,igsh,idasa)
*                                                                      *
*       anatally calculation and output of t-yield with tetra mesh     *
*                                                                      *
************************************************************************

      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

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
      real(8) bzz(mz+1), bnn(mn+1)
!      dimension   tm(maxnt+maxpt,2*nfile+3)
!      dimension   vl3d(nr,1,1)
      dimension   tr(nr,mz,mn,0:mm,2*nfile) ! S.H. 2-> mm (2022.3.11)
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
!      dimension   anatalrst(mm+1,mz+1,mn+1,1,(nr+1)*1*1,
!     &                     1,2*nfile+3) ! S.H. 3-> mm+1 (2022.3.11)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)
      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*21,cijaxs*21,cijaxs2*21
      integer iDaxis
!      dimension   dlr(nr+1) !  integer lr(nr) ---> double dlr(nr)
      real(8),allocatable :: tm(:,:),vl3d(:,:,:),rdata(:,:),dlr(:)

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

      dimension dummy(1)
      data dummy/1.0d0/
*-----------------------------------------------------------------------
      character(28) cfmt                !FURUTA20200615
      character(12) cir                 !FURUTA20200615
      common /redufmt/ iredufmt(itlmax) !FURUTA20200615
*-----------------------------------------------------------------------
      character yen*1
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------
      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )
      allocate (vl3d(nr,1,1),rdata(2,nfile),
     &          dlr(nr+1))

      nrst = 2*nfile+3
      allocate (anatalrst(mm+1,mz+1,mn+1,1,nr*1*1,
     &                     1,nrst))

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               allocate ( lr(nr),vl(nr),val(nr) ) !FURUTA20190204
               call ttetvl(mr,kr,nr,vl,lr)

ccse 2021.08 add (use anatal_rearrange sub.)
               do ir = 1, nr
                  dlr(ir) = dble(lr(ir))
               end do
               dlr(nr+1) = dble(lr(nr))

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

            do 100 ir = 1, nr

               if( itunt(m) .eq. 1 ) then

                  cc = rtfac(m)

               else if( itunt(m) .eq. 2 ) then

                  cc = rtfac(m) / vl(ir)

               end if

            do 100 iz = 1, mz
            do 100 in = 1, mn
            do 100 il = 0, mm ! S.H. 2-> mm (2022.3.11)
            do ntf = 1, nfile

               if( tr(ir,iz,in,il,iat(1,ntf)) .ne. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ir,iz,in,il,iat(1,ntf)),
     &                            tr(ir,iz,in,il,iat(2,ntf)),
     &                            cc)

                  tr(ir,iz,in,il,iat(1,ntf)) = Xa
                  tr(ir,iz,in,il,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ir,iz,in,il,iat(1,ntf)) .gt. cmax )
     &                 cmax = tr(ir,iz,in,il,iat(1,ntf))

                  if( tr(ir,iz,in,il,iat(1,ntf)) .lt. cmin )
     &                 cmin = tr(ir,iz,in,il,iat(1,ntf))

               else

                  if( il .eq. 0 ) isdz = 1
                  tr(ir,iz,in,il,iat(2,ntf)) = 0.0

               end if
            end do
  100       continue

!$OMP PARALLEL
!$OMP& private(ipomp,ir,iz,in,il,ntf,isdz,rdata,answer,rerr)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
            do ir = 1, nr
            do 101 iz = 1, mz
            do 101 in = 1, mn
            do 101 il = 0, mm ! S.H. 2-> mm (2022.3.11)

            if( manatally .eq. 0 ) then ! user defined analysis

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ir,iz,in,il,iat(1,ntf))
                  rdata(2,ntf) = tr(ir,iz,in,il,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(il+1,iz,in,1,ir,1,1) = answer
              anatalrst(il+1,iz,in,1,ir,1,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then ! systematic uncertainty

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ir,iz,in,il,iat(1,ntf))
                  rdata(2,ntf) = tr(ir,iz,in,il,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(il+1,iz,in,1,ir,1,1) ! mean
     &              = fmval
               if ( fmval .ne. 0d0 ) then
                  anatalrst(il+1,iz,in,1,ir,1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(il+1,iz,in,1,ir,1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(il+1,iz,in,1,ir,1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(il+1,iz,in,1,ir,1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(il+1,iz,in,1,ir,1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(il+1,iz,in,1,ir,1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(il+1,iz,in,1,ir,1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(il+1,iz,in,1,ir,1,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then ! c-value dependence

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(il+1,iz,in,1,ir,1,iat(ioe,ntf))
     &              = tr(ir,iz,in,il,iat(ioe,ntf))
                end do
               end do

            end if

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
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        tet, mass, charge axes
*-----------------------------------------------------------------------

         do iz = 1, mz
            bzz(iz) = iz
         end do
         do in = 1, mn
            bnn(in) = in
         end do
         do ir = 1, nr
            vl3d(ir,1,1) = vl(ir)
         end do

*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1
     &          .or. itaxs(m,iax) .eq. 7
     &          .or. itaxs(m,iax) .eq. 14 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 15   ! mass axis
           case ( 7 )
            iDaxis = 17   ! charge axis
           case ( 14 )
            iDaxis = 8    ! tet axis
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.19

C rearrange data from anatalrst to anataldata
          call anatal_rearrange(0,0,4,
     &         iDaxis,ntaxis,anatalrst,
     &         mm+1,     mz,      mn,    0,      nr,  1,  1,    0, nrst, ! S.H. 3-> mm+1 (2022.3.11)
     &         bzz, 1,  bnn,1,  1, 1,  dlr,  dummy,  dummy,1,
     &         vl3d,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          allocate( tm(nijaxs,2*nfile+3) )

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
             call ttetcm(itet0,nr,xcm)
            endif

*-----------------------------------------------------------------------

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

           ic = ij(2)

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

           iz = 0
           ia = 0
           in = 0
           if( nn .gt. 0 ) then

              iz = nt(ic) / 1000
              ia = nt(ic) - iz * 1000

              if( ia .gt. 0 ) then
                 in = ia - iz
                 call chname(idum,ia,iz,chau)
              end if

           end if

*-----------------------------------------------------------------------

           do itmp2 = 1, nrst
           do itmp1 = 1, nijaxs
              tm(itmp1,itmp2) = 0.0d+0
           end do
           end do

           seka = 0.0
           sera = 0.0
           seva = 0.0
           voll = 0.0

           do ijaxs = 1, nijaxs ! (1)

            if( itunt(m) .eq. 1 ) then

               vm = 1.0

            else if( itunt(m) .eq. 2 ) then

               vm = delvol(ianataldata,ijaxs)

            end if

            voll = voll + delvol(ianataldata,ijaxs)

*-----------------------------------------------------------------------

            if( nn .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

              vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
              tm(ijaxs,1) = tm(ijaxs,1) + vn
              tm(ijaxs,2) = tm(ijaxs,2)
     &             + ( vn * anataldata(il+1,ianataldata,ijaxs,2) )**2
              if ( manatally .eq. 1 ) then ! systematic uncertainty
                 tm(ijaxs,3) = tm(ijaxs,3)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,3) )**2
                 tm(ijaxs,4) = tm(ijaxs,4)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,4) )**2
              end if

             end do

            else if( ia .eq. 0 ) then

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

              vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
              tm(ijaxs,1) = tm(ijaxs,1) + vn
              tm(ijaxs,2) = tm(ijaxs,2)
     &             + ( vn * anataldata(il+1,ianataldata,ijaxs,2) )**2
              if ( manatally .eq. 1 ) then ! systematic uncertainty
                 tm(ijaxs,3) = tm(ijaxs,3)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,3) )**2
                 tm(ijaxs,4) = tm(ijaxs,4)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,4) )**2
              end if

             end do

            else

             do il = 0, mm ! S.H. 2-> mm (2022.3.11)

              vn  = vm * anataldata(il+1,ianataldata,ijaxs,1)
              tm(ijaxs,1) = tm(ijaxs,1) + vn
              tm(ijaxs,2) = tm(ijaxs,2)
     &             + ( vn * anataldata(il+1,ianataldata,ijaxs,2) )**2
              if ( manatally .eq. 1 ) then ! systematic uncertainty
                 tm(ijaxs,3) = tm(ijaxs,3)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,3) )**2
                 tm(ijaxs,4) = tm(ijaxs,4)
     &                + ( vn * anataldata(il+1,ianataldata,ijaxs,4) )**2
              end if
             end do

            end if

            seka = seka + tm(ijaxs,1)
            sera = sera + tm(ijaxs,2)
            seva = seva + vm

            if( tm(ijaxs,1) .gt. 0.0 ) then

             tm(ijaxs,2) = sqrt( tm(ijaxs,2) ) / tm(ijaxs,1)
             tm(ijaxs,1) = tm(ijaxs,1) / vm
             if ( manatally .eq. 1 ) then ! systematic uncertainty
                tm(ijaxs,3) = sqrt( tm(ijaxs,3) ) / tm(ijaxs,1)
                tm(ijaxs,4) = sqrt( tm(ijaxs,4) ) / tm(ijaxs,1)
             end if

            end if

           end do                 ! ijaxs = 1, nijaxs (1)

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
               write(iot,'( "y: Number ",a15)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

           write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

           if( itanl(m) .gt. 0 ) then
              write(iot,'( "p: ",200a1)')
     &             ( itang(m)(i:i),i = 1, itanl(m) )
           end if

           if( itsans(m) .gt. 0 ) then
              call write_sangel(iot,m,0)
           end if

           if ( manatally .eq. 0 ) then ! user defined analysis
            if( itaxs(m,iax) .eq. 1 ) then ! mass axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5    y,hl0       n")')

             else
                write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

             end if

             write(iot,'( "#  mass     number      r.err")')

            else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5    y,hl0       n")')

             else
                write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

             end if

             write(iot,'( "#  charge   number      r.err")')

            else if( itaxs(m,iax) .eq. 14 ) then ! tet axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                    "y,l3        n")')

             else
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                    "y1,l3       ny2 dy1=[y1*y2]")')

             end if

             write(iot,'( "#  num    tetra   volume  ",
     &                    "   number      r.err")')

            end if

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if( itaxs(m,iax) .eq. 1 ) then ! mass axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5   y,hl0       n n n")')

             else
                write(iot,'( "h: x-0.5   y1,hl0      ",
     &                      "ny2 dy1=[y1*y2] n n")')

             end if

             write(iot,'( "#  mass    number      "
     &                   ,"r.err(tot, syst, stat)")')

            else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: x-0.5   y,hl0       n n n")')

             else
                write(iot,'( "h: x-0.5   y1,hl0      ",
     &                      "ny2 dy1=[y1*y2] n n")')

             end if

             write(iot,'( "#charge    number      "
     &                   ,"r.err(tot, syst, stat)")')

            else if( itaxs(m,iax) .eq. 14 ) then ! tet axis
             if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h:",3x,"x",7x,"n",3x,"n",12x,
     &                    "y1,l3",7x,"n")')

             else
                write(iot,'( "h:",3x,"x",7x,"n",3x,"n",12x,
     &                    "y1,l3",7x,"ny2 dy1=[y1*y2] n n")')

             end if

             write(iot,'( "#  num     reg   volume  ",
     &                    "     number      r.err(tot, syst, stat)")')

            end if

           else if ( manatally .eq. 2 ) then ! c-value dependence
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h:   x",10x,
     &                    "y,l3        n")')

            else
               write(iot,'( "h:   x",10x,
     &                    "y1,l3       ny2 dy1=[y1*y2]")')

            end if

            write(iot,'( "#  c-value  ",
     &                    "   number      r.err")')

           end if

*-----------------------------------------------------------------------

           do ijaxs = 1, nijaxs
              if( tm(ijaxs,1) .gt. 0.0d0 ) goto 201
           end do

 201       im = ijaxs

           do ijaxs = nijaxs, im + 1, -1
              if( tm(ijaxs,1) .gt. 0.0d0 ) goto 211
           end do

           jm = im
 211       jm = ijaxs

           im = max( 1, im - 2 )
           jm = jm + 2

*-----------------------------------------------------------------------

           if( itaxs(m,iax).eq.1 .or. itaxs(m,iax).eq.7 ) then ! mass, charge axis
            if( im .eq. 1 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(3x,f4.1,1x,1pe13.4,0p1f8.4)')
     &                                   0.5, 0.0, 0.0

             else if ( manatally .eq. 1 ) then ! systematic uncertainty
                write(iot,'(3x,f4.1,1x,1pe13.4,0p3f8.4)')
     &                                   0.5, 0.0, 0.0, 0.0, 0.0

             else if ( manatally .eq. 2 ) then ! c-value dependence

             end if

            end if
           end if

           do ijaxs = 1, nijaxs ! (2)

            if ( manatally .eq. 0 ) then ! user defined analysis
             if( itaxs(m,iax) .eq. 1 ) then ! mass axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p1f8.4)')
     &                ijaxs+iz, (tm(ijaxs,irst),irst=1,2)
              end if

             else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p1f8.4)')
     &                ijaxs, (tm(ijaxs,irst),irst=1,2)
              end if

             else if( itaxs(m,iax) .eq. 14 ) then ! tet axis
                write(iot,'(i8,1x,i8,1pe13.4,1pe13.4,0pf8.4)')
     &               ir, ireg, vl(ir), sek, ser
                write(iot,'(1x,i5,1x,i7,1pe13.4,1pe13.4,0p1f8.4)')
     &               ijaxs, lr(ijaxs), delvol(ianataldata,ijaxs)
     &               ,(tm(ijaxs,irst),irst=1,2)
cFURUTA20190208 OpenFOAM output
                if ( itfoam(m) .eq. 1) then
                   write(iot2,'(1pe13.4)') tm(ijaxs,1)
cFURUTA20191028 CSV output
                elseif(itfoam(m) .eq. 2) then
                   write(buf,'(i8,5(",",1pe13.4),",",0pf8.4)')
     &                  lr(ijaxs), xcm(1:3,ijaxs)
     &                  ,delvol(ianataldata,ijaxs)
     &                  ,(tm(ijaxs,irst),irst=1,2)
                   call remove_spaces(buf)
                   write(iot2,'(a)')trim(buf)
                endif

             end if

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
             if( itaxs(m,iax) .eq. 1 ) then ! mass axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p3f8.4)')
     &                ijaxs+iz, (tm(ijaxs,irst),irst=1,4)
              end if

             else if( itaxs(m,iax) .eq. 7 ) then ! charge axis
              if( ijaxs.ge.im .and. ijaxs.le.jm ) then ! from im to jm
                 write(iot,'(3x,i4,1x,1pe13.4,0p3f8.4)')
     &                ijaxs, (tm(ijaxs,irst),irst=1,4)
              end if

             else if( itaxs(m,iax) .eq. 14 ) then ! tet axis
                write(iot,'(1x,i5,1x,i7,1pe13.4,1pe13.4,0p3f8.4)')
     &               ijaxs, lr(ijaxs), delvol(ianataldata,ijaxs)
     &               ,(tm(ijaxs,irst),irst=1,4)

             end if

            else if ( manatally .eq. 2 ) then ! c-value dependence
               write(iot,'(3x,i4,1x,1pe13.4,0p1f8.4)')
     &              ijaxs, (tm(ijaxs,irst),irst=1,2)

            end if

           end do                 ! ijaxs = 1, nijaxs (2)

           if( seka .gt. 0.0 ) then
              sera = sqrt( sera ) / seka
              if( itunt(m) .eq. 2 ) seka = seka / seva
           end if

           if( itaxs(m,iax) .eq. 1 ) then ! mass axis
            write(iot,'(/"#   sum ",1pe13.4,0pf8.4)') seka,sera

            if( itunt(m) .eq. 2 ) then
               write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm
            end if

           else if( itaxs(m,iax) .eq. 14 ) then ! tet axis
            write(iot,'(/"#   sum over  ",1pe13.4,1pe13.4,0pf8.4)')
     &             voll, seka, sera

            if( itunt(m) .eq. 2 ) then
               write(iot,'(
     &         "#   vol sum = ",1p1e13.4,"  [cm^3]")') seva
            end if

           end if

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
              if ( ijtmp .eq. 2 ) then
                 if ( iz.gt.0)
     &                write(iot,'(a15,"&=&",i5,3x,"(",a4,")")')
     &                cij(ijtmp), iz, elmnt(iz)
              else if ( ijtmp .eq. 3 ) then
                 if ( ia .gt. 0 )
     &                write(iot,'(a15,"&=&",i5)') cij(ijtmp), ia
              else if ( ijtmp .eq. 5 ) then
                 write(iot,'(a15,"&=&",i5)')
     &                cij(ijtmp), idnint(fg(itmp))
                 write(iot,'(3x,"volume",6x,"&=&",1pe13.4," [cm^3]")')
     &                delvol(ianataldata,1)
              else
                 write(iot,'(a15,"&=&",i5)')
     &                cij(ijtmp), idnint(fg(itmp))
              end if
              itmp = itmp - ij(ijtmp) + nij(ijtmp)
             else
                write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                           /1pe13.4,2x,"$--$",1pe13.4)')
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
          deallocate( tm )

*-----------------------------------------------------------------------
cFURUTA20190208 OpenFOAM output
             if ( itfoam(m) .eq. 1) then !FURUTA20191028
               write(iot2,'(a)')')'
               close(iot2)
             endif
*-----------------------------------------------------------------------

cFURUTA20190208 OpenFOAM output
            if ( itfoam(m) .ne. 0) then
              deallocate( foamfIType,foamfIndex )
            endif

*-----------------------------------------------------------------------
*        chart axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 ) then

          iDaxis = 18         ! chart axis

C rearrange data from anatalrst to anataldata
          call anatal_rearrange(0,0,1,
     &         iDaxis,ntaxis,anatalrst,
     &         mm+1,     mz,      mn,    0,      nr,  1,  1,    0, nrst, ! S.H. 3-> mm+1 (2022.3.11)
     &         bzz, 1,  bnn,1,  1, 1,  dlr,  dummy,  dummy,1,
     &         vl3d,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          inum = 0

*-----------------------------------------------------------------------

          do ianataldata = 1, nanataldata

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

           ir = ij(5)

           do il = 0, mm         ! loop for isomeric level ! S.H. 2-> mm (2022.3.11)

*-----------------------------------------------------------------------

            itmax  = 0
            icmax  = 0
            inmax  = 0

            do iz = 1, maxpt
             do in = 1, maxnt

              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
                 ijaxs = in + (iz-1)*maxnt
              else
                 ijaxs = igetiznmp(iz,in,il,m)
              end if
              if( anataldata(il+1,ianataldata,ijaxs,1).gt.0d0 ) then
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

            if( itmax .eq. 0 ) exit

*-----------------------------------------------------------------------

            ireg = idnint(fg(ij(5)))

            write(iot,'(/"#",78("-"))')

            inum = inum + 1

            if( inum .eq. 1 .and. il .eq. 0 ) then
               write(iot,'( "#newpage:")')
            else
               write(iot,'( " newpage:")')
            end if

*-----------------------------------------------------------------------

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
                  write(iot,'(/"x: ",a21)') cijaxs
               else
                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))
               end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
               if( itayl(m) .eq. 0 ) then
                  write(iot,'(/"y: ",a21)') cijaxs2
               else
                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))
               end if
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

             do ijaxs2 = icmax+2, 1, -1
              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
              write(iot,'(1p10e11.3)')
     &       ( anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                ijaxs = 1, inmax+2 )

              else
               write(iot,'(1p10e11.3)')
     &                ( anataldata(il+1,ianataldata,
     &                igetiznmp(ijaxs2,ijaxs,il,m),ioe),
     &               ijaxs = 1, inmax+2 )
              end if
             end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

             write(iot,'(/"#   Z    N  Mass  ",
     &                      "  number    r.err")')

             do ijaxs2 = 1, icmax
              do ijaxs = 1, inmax

              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
               if( anataldata(il+1,ianataldata
     &                ,ijaxs+(ijaxs2-1)*nijaxs,1) .gt. 0.d0 ) then

                  write(iot,'(3i5,1pe13.4,0pf8.4)')
     &                 ijaxs2, ijaxs, ijaxs2+ijaxs,
     &           anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &           anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,2)

               end if

              else
               if( anataldata(il+1,ianataldata
     &                ,igetiznmp(ijaxs2,ijaxs,il,m),1) .gt. 0.d0 ) then

                  write(iot,'(3i5,1pe13.4,0pf8.4)')
     &                 ijaxs2, ijaxs, ijaxs2+ijaxs,
     &      anataldata(il+1,ianataldata,igetiznmp(ijaxs2,ijaxs,il,m),1),
     &      anataldata(il+1,ianataldata,igetiznmp(ijaxs2,ijaxs,il,m),2)

               end if
              end if

              end do
             end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

             write(iot,'(/"#   Z = 1 to ",i3/
     &                      "#   N = 1 to ",i3)')
     &              icmax+2, inmax+2

             write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'Z/N', ( dble( in ), in = 1, inmax+2 )

             do ijaxs2 = icmax+2, 1, -1
              if( itnzn(m).eq.0 ) then ! S.H. added (2022.3.11)
                write(iot,'(1p1000e11.3)')
     &               dble( ijaxs2 ),
     &         ( anataldata(il+1,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &            ijaxs = 1, inmax+2 )

              else
                write(iot,'(1p1000e11.3)')
     &               dble( ijaxs2 ),
     &    ( anataldata(il+1,ianataldata,igetiznmp(ijaxs2,ijaxs,il,m),1),
     &               ijaxs = 1, inmax+2 )
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

      else

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else
        if( itazl(m) .eq. 0 ) then
           write(iot,'("y: Number ",a15)') hsunit(itunt(m))
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
       write(iot,'("wt: s(0.7)",/a1,"vspace{-3}")') yen
       itmp = 0
       do ijtmp = 2,8
        if (cij(ijtmp) .ne. 'F' ) then
         itmp = itmp + ij(ijtmp)
         if ( ibin(ijtmp) .eq. 0 ) then
            if ( ijtmp .eq. 5 ) then
               write(iot,'(a15,"&=&",i5)')
     &              cij(ijtmp), idnint(fg(itmp))
               write(iot,'(3x,"volume")')
               write(iot,'(5x,"&=&",1pe13.4," [cm^3]")')
     &              delvol(ianataldata,1)
            else
               write(iot,'(a15,"&=&",i5)')
     &              cij(ijtmp), idnint(fg(itmp))
            end if
            itmp = itmp - ij(ijtmp) + nij(ijtmp)
         else
            write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                           /1pe13.4,2x,"$--$",1pe13.4)')
     &           cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
            itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
         end if
        end if
       end do
       if(il .eq. 0) then       ! specify isomeric level at the bottom of figure
          write(iot,'("il = 0 (Ground state)}")')
       elseif(il .eq. 1) then
          write(iot,'("il = 1 (1st isomer)}")')
       else
          write(iot,'("il = 2 (2nd isomer)}")')
       endif

       write(iot,'("e:")')

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

           end do                    ! il = 0, 2 ! loop for isomeric level
          end do                    ! ianataldata = 1, nanataldata

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( anataldata )
          deallocate( delvol )

*-----------------------------------------------------------------------
*        dchain axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 13 ) then

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
                          t0_tr = anatalrst(il+1,iz,in,1,ir,1,1)
                        else
                          t0_tr =
     &                   anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,ir,1,1)
                        end if
                        if( t0_tr .ne. 0.0 ) goto 120

                     end do
                     end do

                     goto 170

  120             im = in

                     do in = maxnt, im + 1, -1
                     do ir = 1, nr

                        if( itnzn(m).eq.0 ) then
                          t0_tr = anatalrst(il+1,iz,in,1,ir,1,1)
                        else
                          t0_tr =
     &                   anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,ir,1,1)
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
     &           lr(ir), ( anatalrst(il+1,iz,i,1,ir,1,ioe), i = n3, n4 )
                   else
                     write(iot,'(i7,1p12e11.3)') lr(ir),
     &               ( anatalrst(0+1,igetiznmp(iz,i,il,m),1,1,ir,1,ioe),
     &                i = n3, n4 )
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
     &       lr(ir), ( anatalrst(il+1,iz,i,1,ir,1,2)*100.0, i = n3, n4 )
                      else
                        write(iou,'(i7,1p12e11.3)') lr(ir),
     &         ( anatalrst(0+1,igetiznmp(iz,i,il,m),1,1,ir,1,2)*100.0,
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
                  if(anatalrst(il+1,iz,in,1,ir,1,ioe).ne.0d0)then
                   write(iot,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                  (anatalrst(il+1,iz,in,1,ir,1,ioee),ioee=1,2)
                   if(itout(m).ne.0)then
                    write(iou,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                   anatalrst(il+1,iz,in,1,ir,1,2)
                   endif
                  endif
                 else
                  if(anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,ir,1,ioe)
     &                .ne.0d0)then
                   write(iot,cfmt)ir,iz*10000+(in+iz)*10+il,
     &    (anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,ir,1,ioee),ioee=1,2)
                   if(itout(m).ne.0)then
                    write(iou,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                  anatalrst(0+1,igetiznmp(iz,in,il,m),1,1,ir,1,2)
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

         end if   ! itaxs(m,iax)

*-----------------------------------------------------------------------

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

      deallocate( vl3d,rdata,dlr )

      return
      end
