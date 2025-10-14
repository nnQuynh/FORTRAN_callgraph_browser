!***********************************************************************
!                                                                      *
      subroutine sumtal_read_theat(m,iax,ntf,ierr)
!                                                                      *
!     original routine is read_theat in restheat.f                     *
!                                                                      *
!   m: the tally number, index of ital.                                *
!***********************************************************************

        use RESTALMOD
        use sumtallymod
        use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'
        include 'param01.inc'

*-----------------------------------------------------------------------

        common /mpi00/ npe, me

*-----------------------------------------------------------------------

        common /talout/ itall
        common /talmm/  nmmax, lmmax, itlmx

*-----------------------------------------------------------------------

        common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
        common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                  rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
        common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                  rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
        common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                  rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
        common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                  rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
        common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                  rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall17/ itndy(itlmax)
        common /tall18/ ithet(itlmax)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

        common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
        character crfln*100

        dimension     idas(mdas*2)
        equivalence ( das, idas )

*-----------------------------------------------------------------------

        character dsin(0:9,2)*200
        dimension idsi(0:9,2)
        dimension ill(0:9,2), ilf(0:9,2)
        dimension jsn(2), jsi(2)
        dimension ierrs(2)

        integer   ntf   ! Number of tallyfname (tally file name)
*-----------------------------------------------------------------------

        ierr = 0
        noe  = 1

        if ( any( itaxs(m,iax) .eq. (/ 7, 8, 9, 10 /) )
     &      .and. ittwo(m) .ne. 4 ) noe = 2

*-----------------------------------------------------------------------
*   open restart file
*-----------------------------------------------------------------------
        call sumtal_open_resfile
     &      (m,noe,ntf,jsn,jsi,dsin,idsi,ill,ilf,newtall,ierr)

        if ( newtall .ne. 0 ) goto 900  !! it's new tally
        if ( ierr    .ne. 0 ) goto 900

*-----------------------------------------------------------------------
*   check tally
*-----------------------------------------------------------------------
        do 800 ioe = 1, noe

        if(ireschk.eq.0) then  ! T.Sato 2013/10/19
        idsi(0,ioe) = ltallyfname(ntf,m)
        dsin(0,ioe)(1:idsi(0,ioe))=
     &                         tallyfname(ntf,m)(1:ltallyfname(ntf,m))
        call check_theat(m,iax,jsn(ioe),jsi(ioe),
     &                   dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                   ierr)
        endif

        if (ierr .ne. 0 ) goto 900

*-----------------------------------------------------------------------

        rewind(jsi(ioe))
        ill(jsn(ioe),ioe) = 1

*-----------------------------------------------------------------------

        idas0 = nmmax
        idas1 = lmmax
        idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
        idas3 = idas1 + itrgn(m) + ( itrgn(m) + mod(itrgn(m),2) ) / 2
        idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_hetreg(m,iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),
     &      itpan(m),itrgn(m),itrgm(m),
     &      idas_itreg(itreg(m)),itenm(m),das_iterg(iterg(m)),
     &      trRES(irestalm(m)),das(ithet(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call read_hetrz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),
     &      itpan(m),itrnm(m),itznm(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      itenm(m),das_iterg(iterg(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call read_hetxyz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),
     &      itpan(m),itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      itenm(m),das_iterg(iterg(m)),
     &      trRES(irestalm(m)),das(ithet(m)))

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_hetreg(m,
     &      itndy(m),itout(m),
     &      itpan(m),itrgn(m),itrgm(m),
     &      idas_itreg(itreg(m)),itenm(m),das_iterg(iterg(m)),
     &      trRES(irestalm(m)),das(ithet(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call restore_hetrz(m,
     &      itndy(m),itout(m),
     &      itpan(m),itrnm(m),itznm(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      itenm(m),das_iterg(iterg(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call restore_hetxyz(m,
     &      itndy(m),itout(m),
     &      itpan(m),itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      itenm(m),das_iterg(iterg(m)),
     &      trRES(irestalm(m)),das(ithet(m)))

        end if

  900   continue
*-----------------------------------------------------------------------
*   close restart file
*-----------------------------------------------------------------------

        do ioe = 1, noe

          close(jsi(ioe))

        end do

*-----------------------------------------------------------------------
  999   continue

      end subroutine sumtal_read_theat


!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_hetreg(m,ntf,
     &                             ne,nd,nr,
     &                             sumfactor,weightRate,
     &                             trRES,trSUMTAL,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by T.Miura on 2015/04/27                                 *
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
      integer           ne
      integer           nd
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nd,nr,0:ne,2)
      double precision  trSUMTAL(nd,nr,0:ne,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  calcfact1,calcfact2

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ie
      integer           ik
      integer           ir

      integer,save   :: maxcasd
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr
       do ik = 1, nd
        do ie = 0, ne

         if( ntf > 1 ) then

          ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
          !                 + F r2 sig_N2 xi wi
          trSUMTAL(ik,ir,ie,1) =
     &    trSUMTAL(ik,ir,ie,1) +
     &            ( sumfactor * weightRate(ntf) *
     &       trRES(ik,ir,ie,1) )

          ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
          !                    + F**2 r2**2 sig_N2(xiwi)**2
          trSUMTAL(ik,ir,ie,2) =
     &    trSUMTAL(ik,ir,ie,2) +
     &            ( sumfactor**2 * weightRate(ntf)**2 *
     &       trRES(ik,ir,ie,2) )

         else
          trSUMTAL(ik,ir,ie,1) =
     &            ( sumfactor * weightRate(ntf) *
     &       trRES(ik,ir,ie,1) )

          trSUMTAL(ik,ir,ie,2) =
     &            ( sumfactor**2 * weightRate(ntf)**2 *
     &       trRES(ik,ir,ie,2) )

         end if

         if( trSUMTAL(ik,ir,ie,1) > cmax )
     &     cmax = trSUMTAL(ik,ir,ie,1)
         if( trSUMTAL(ik,ir,ie,1) < cmin )
     &     cmin = trSUMTAL(ik,ir,ie,1)

        end do      ! ie loop end
       end do       ! ik loop end
      end do        ! ir loop end

! sumover
      calcfact1 =  sumfactor * weightRate(ntf)
      calcfact2 = sumfactor**2 * weightRate(ntf)**2
      call sumtal_calc_stdev_tracreg_sub(m, ntf,
     &                                 calcfact1,calcfact2)

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = resc2SUMTAL
     &            + weightRate(ntf) * resc3(m) * resc2(m)/resc3(m)

      ! N = N1 + N2
      resc3SUMTAL = resc3SUMTAL + resc3(m)

      resc2(m) = resc2SUMTAL
      resc3(m) = resc3SUMTAL

!-----------------------------------------------------------------------
      if( ntf > 1 ) then
       if( maxcasd /= maxcasres ) then
        ierr = 1
        m_err = 'maxcas is inconsistent.'
        ErrCha = ''
        ErrID = 'L:338/R:sumtal_calc_stdev_hetreg/F:sumtalheat.f' !E56_001_001


        write(*,'(/'' ***** Error in sumtal calc stdev hetreg *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_hetreg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_hetrz(m,ntf,
     &                             ne,nd,nz,nr,
     &                             sumfactor,weightRate,
     &                             trRES,trSUMTAL,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by T.Miura on 2015/04/27                                 *
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
      integer           ne
      integer           nd
      integer           nz
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nd,nr,nz,0:ne,2)
      double precision  trSUMTAL(nd,nr,nz,0:ne,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           im

      integer           ir
      integer           iz

      integer           ik
      integer           ie

      integer,save   :: maxcasd
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr
       do iz = 1, nz
        do ik = 1, nd
         do ie = 0, ne

          if( ntf > 1 ) then

           ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
           !                 + F r2 sig_N2 xi wi
           trSUMTAL(ik,ir,iz,ie,1) =
     &            trSUMTAL(ik,ir,iz,ie,1) +
     &            ( sumfactor * weightRate(ntf) *
     &            trRES(ik,ir,iz,ie,1) )

           ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
           !                    + F**2 r2**2 sig_N2(xiwi)**2
           trSUMTAL(ik,ir,iz,ie,2) =
     &          trSUMTAL(ik,ir,iz,ie,2) +
     &          ( sumfactor**2 * weightRate(ntf)**2 *
     &          trRES(ik,ir,iz,ie,2) )

          else
           trSUMTAL(ik,ir,iz,ie,1) =
     &            ( sumfactor * weightRate(ntf) *
     &            trRES(ik,ir,iz,ie,1) )

           trSUMTAL(ik,ir,iz,ie,2) =
     &          ( sumfactor**2 * weightRate(ntf)**2 *
     &          trRES(ik,ir,iz,ie,2) )

          end if

          if( trSUMTAL(ik,ir,iz,ie,1) > cmax )
     &         cmax = trSUMTAL(ik,ir,iz,ie,1)
          if( trSUMTAL(ik,ir,iz,ie,1) < cmin )
     &         cmin = trSUMTAL(ik,ir,iz,ie,1)

         end do      ! ie loop end
        end do       ! ik loop end
       end do        ! iz loop end
      end do         ! ir loop end

! sumover
      calcfact1 =  sumfactor * weightRate(ntf)
      calcfact2 = sumfactor**2 * weightRate(ntf)**2
      call sumtal_calc_stdev_tracreg_sub(m, ntf,
     &                                 calcfact1,calcfact2)

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = resc2SUMTAL
     &            + weightRate(ntf) * resc3(m) * resc2(m)/resc3(m)

      ! N = N1 + N2
      resc3SUMTAL = resc3SUMTAL + resc3(m)

      resc2(m) = resc2SUMTAL
      resc3(m) = resc3SUMTAL

!-----------------------------------------------------------------------
      if( ntf > 1 ) then
       if( maxcasd /= maxcasres ) then
        ierr = 1
        m_err = 'maxcas is inconsistent.'
        ErrCha = ''
        ErrID = 'L:487/R:sumtal_calc_stdev_hetrz/F:sumtalheat.f' !E56_001_002


        write(*,'(/'' ***** Error in sumtal calc stdev hetrz *****''/)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_hetrz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_hetxyz(m,ntf,
     &                             ne,nd,nz,ny,nx,
     &                             sumfactor,weightRate,
     &                             trRES,trSUMTAL,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by T.Miura on 2015/04/27                                 *
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
      integer           ne
      integer           nd
      integer           nz
      integer           ny
      integer           nx

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nd,nx,ny,nz,0:ne,2)
      double precision  trSUMTAL(nd,nx,ny,nz,0:ne,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           im

      integer           ix
      integer           iy
      integer           iz

      integer           ik
      integer           ie

      integer,save   :: maxcasd
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ix = 1, nx
       do iy = 1, ny
        do iz = 1, nz
         do ik = 1, nd
          do ie = 0, ne

           if( ntf > 1 ) then

            ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
            !                 + F r2 sig_N2 xi wi
            trSUMTAL(ik,ix,iy,iz,ie,1) =
     &      trSUMTAL(ik,ix,iy,iz,ie,1) +
     &        ( sumfactor * weightRate(ntf) *
     &         trRES(ik,ix,iy,iz,ie,1) )

            ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
            !                    + F**2 r2**2 sig_N2(xiwi)**2
            trSUMTAL(ik,ix,iy,iz,ie,2) =
     &      trSUMTAL(ik,ix,iy,iz,ie,2) +
     &        ( sumfactor**2 * weightRate(ntf)**2 *
     &         trRES(ik,ix,iy,iz,ie,2) )

           else
            trSUMTAL(ik,ix,iy,iz,ie,1) =
     &        ( sumfactor * weightRate(ntf) *
     &         trRES(ik,ix,iy,iz,ie,1) )

            trSUMTAL(ik,ix,iy,iz,ie,2) =
     &        ( sumfactor**2 * weightRate(ntf)**2 *
     &         trRES(ik,ix,iy,iz,ie,2) )

           end if

           if( trSUMTAL(ik,ix,iy,iz,ie,1) > cmax )
     &           cmax = trSUMTAL(ik,ix,iy,iz,ie,1)
           if( trSUMTAL(ik,ix,iy,iz,ie,1) < cmin )
     &           cmin = trSUMTAL(ik,ix,iy,iz,ie,1)

          end do    ! ie loop end
         end do     ! it loop end
        end do      ! ix loop end
       end do       ! iy loop end
      end do        ! iz loop end

! sumover
      calcfact1 =  sumfactor * weightRate(ntf)
      calcfact2 = sumfactor**2 * weightRate(ntf)**2
      call sumtal_calc_stdev_tracreg_sub(m, ntf,
     &                                 calcfact1,calcfact2)

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = resc2SUMTAL
     &            + weightRate(ntf) * resc3(m) * resc2(m)/resc3(m)

      ! N = N1 + N2
      resc3SUMTAL = resc3SUMTAL + resc3(m)

      resc2(m) = resc2SUMTAL
      resc3(m) = resc3SUMTAL

!-----------------------------------------------------------------------
      if( ntf > 1 ) then
       if( maxcasd /= maxcasres ) then
        ierr = 1
        m_err = 'maxcas is inconsistent.'
        ErrCha = ''
        ErrID = 'L:640/R:sumtal_calc_stdev_hetxyz/F:sumtalheat.f' !E56_001_003


        write(*,'(/'' ***** Error in sumtal calc stdev hetxyz *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_hetxyz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_hetreg(m,ntf,nfile,
     &                               ne,nd,nr,
     &                               sumfactor,weightRate,sumWR,
     &                               trRES,trSUMTAL,
     &                               resc2SUMTAL,resc3SUMTAL,
     &                               ierr)
!                                                                      *
!     created by T.Miura on 2015/04/27                                 *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------------------------------
      include 'param.inc'
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
      integer           m
      integer           ntf
      integer           nfile
      integer           ne
      integer           nd
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (nd,nr,0:ne,2)
      double precision  trSUMTAL(nd,nr,0:ne,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir
      integer           ik
      integer           ie

      double precision  Xa, sigx
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr
       do ik = 1, nd
        do ie = 0, ne

         if( trRES(ik,ir,ie,1) .gt. 0.d0 ) then
          call calc_stdev(m,Xa,sigx,
     &                    trRES(ik,ir,ie,1),
     &                    trRES(ik,ir,ie,2),
     &                    1.0d+0)
          trRES(ik,ir,ie,1) = Xa
          trRES(ik,ir,ie,2) = sigx
         end if

         if( ntf > 1 ) then

          ! X_bar = F Sigma rj/r Xj_bar
          trSUMTAL(ik,ir,ie,1) =
     &    trSUMTAL(ik,ir,ie,1) +
     &       weightRate(ntf)
     &     * trRES(ik,ir,ie,1)

          ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
          trSUMTAL(ik,ir,ie,2) =
     &    trSUMTAL(ik,ir,ie,2) +
     &    (  trRES(ik,ir,ie,2)
     &     * trRES(ik,ir,ie,1) )**2
     &     * weightRate(ntf)**2

         ! ntf = 1
         else
          ! X_bar = F Sigma rj/r Xj_bar
          trSUMTAL(ik,ir,ie,1) =
     &       weightRate(ntf)
     &     * trRES(ik,ir,ie,1)

          ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
          trSUMTAL(ik,ir,ie,2) =
     &    (  trRES(ik,ir,ie,2)
     &     * trRES(ik,ir,ie,1) )**2
     &     * weightRate(ntf)**2
         end if

        end do     ! ie loop end
       end do      ! ik loop end
      end do       ! ir loop end

!sumover  tr_sum -> trSUMTAL_sum
      call calc_stdev_tr_sum(m, 1.0d0)
      call sumtal_calc_average_sub(m, ntf,weightRate(ntf))

!-----------------------------------------------------------------------
  100 continue
      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = resc2SUMTAL
     &            + weightRate(ntf) * resc3(m) * resc2(m)/resc3(m)

      ! N = N1 + N2
      resc3SUMTAL = resc3SUMTAL + resc3(m)
!-----------------------------------------------------------------------
      if( ntf == nfile ) then
       do ir = 1, nr
        do ik = 1, nd
         do ie = 0, ne

          ! X_bar
          trSUMTAL(ik,ir,ie,1) =
     &    trSUMTAL(ik,ir,ie,1) * sumfactor/sumWR

          ! sig_x
          trSUMTAL(ik,ir,ie,2) = sqrt(
     &    trSUMTAL(ik,ir,ie,2)) * sumfactor/sumWR


          ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
          trSUMTAL(ik,ir,ie,2) =
     &   (trSUMTAL(ik,ir,ie,2)**2 *
     &    resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &    resc3SUMTAL *
     &    trSUMTAL(ik,ir,ie,1)**2) *
     &    (resc2SUMTAL/resc3SUMTAL)**2

          ! Sigma xi wi = X_bar W
          trSUMTAL(ik,ir,ie,1) =
     &    trSUMTAL(ik,ir,ie,1) * resc2SUMTAL

          if( trSUMTAL(ik,ir,ie,1) > cmax )
     &          cmax = trSUMTAL(ik,ir,ie,1)
          if( trSUMTAL(ik,ir,ie,1) < cmin )
     &          cmin = trSUMTAL(ik,ir,ie,1)

         end do     ! ie loop end
        end do      ! ik loop end
       end do       ! ir loop end

!sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_hetreg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_hetrz(m,ntf,nfile,
     &                               ne,nd,nz,nr,
     &                               sumfactor,weightRate,sumWR,
     &                               trRES,trSUMTAL,
     &                               resc2SUMTAL,resc3SUMTAL,
     &                               ierr)
!                                                                      *
!     created by T.Miura on 2015/04/27                                 *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------------------------------
      include 'param.inc'
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
      integer           m
      integer           ntf
      integer           nfile

      integer           ne
      integer           nd
      integer           nz
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (nd,nr,nz,0:ne,2)
      double precision  trSUMTAL(nd,nr,nz,0:ne,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir
      integer           iz
      integer           ik
      integer           ie

      double precision  Xa, sigx
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr
       do iz = 1, nz
        do ik = 1, nd
         do ie = 0, ne

          if( trRES(ik,ir,iz,ie,1) .gt. 0.d0 ) then
           call calc_stdev(m,Xa,sigx,
     &                     trRES(ik,ir,iz,ie,1),
     &                     trRES(ik,ir,iz,ie,2),
     &                     1.0d+0)
           trRES(ik,ir,iz,ie,1) = Xa
           trRES(ik,ir,iz,ie,2) = sigx
          end if

          if( ntf > 1 ) then

           ! X_bar = F Sigma rj/r Xj_bar
           trSUMTAL(ik,ir,iz,ie,1) =
     &     trSUMTAL(ik,ir,iz,ie,1) +
     &     weightRate(ntf)
     &      * trRES(ik,ir,iz,ie,1)

           ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
           trSUMTAL(ik,ir,iz,ie,2) =
     &     trSUMTAL(ik,ir,iz,ie,2) +
     &     (  trRES(ik,ir,iz,ie,2)
     &      * trRES(ik,ir,iz,ie,1) )**2
     &      * weightRate(ntf)**2

          ! ntf = 1
          else
           ! X_bar = F Sigma rj/r Xj_bar
           trSUMTAL(ik,ir,iz,ie,1) =
     &     weightRate(ntf)
     &      * trRES(ik,ir,iz,ie,1)

           ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
           trSUMTAL(ik,ir,iz,ie,2) =
     &     (  trRES(ik,ir,iz,ie,2)
     &      * trRES(ik,ir,iz,ie,1) )**2
     &      * weightRate(ntf)**2
          end if

         end do      ! ie loop end
        end do       ! ik loop end
       end do        ! iz loop end
      end do         ! ir loop end

!sumover  tr_sum -> trSUMTAL_sum
      call calc_stdev_tr_sum(m, 1.0d0)
      call sumtal_calc_average_sub(m, ntf,weightRate(ntf))

!-----------------------------------------------------------------------
  100 continue
      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = resc2SUMTAL
     &            + weightRate(ntf) * resc3(m) * resc2(m)/resc3(m)

      ! N = N1 + N2
      resc3SUMTAL = resc3SUMTAL + resc3(m)
!-----------------------------------------------------------------------
      if( ntf == nfile ) then
       do ir = 1, nr
        do iz = 1, nz
         do ik = 1, nd
          do ie = 0, ne

           ! X_bar
           trSUMTAL(ik,ir,iz,ie,1) =
     &     trSUMTAL(ik,ir,iz,ie,1) * sumfactor/sumWR

           ! sig_x
           trSUMTAL(ik,ir,iz,ie,2) = sqrt(
     &     trSUMTAL(ik,ir,iz,ie,2)) * sumfactor/sumWR


           ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
           trSUMTAL(ik,ir,iz,ie,2) =
     &    (trSUMTAL(ik,ir,iz,ie,2)**2 *
     &     resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &     resc3SUMTAL *
     &     trSUMTAL(ik,ir,iz,ie,1)**2) *
     &    (resc2SUMTAL/resc3SUMTAL)**2

           ! Sigma xi wi = X_bar W
           trSUMTAL(ik,ir,iz,ie,1) =
     &     trSUMTAL(ik,ir,iz,ie,1) * resc2SUMTAL

           if( trSUMTAL(ik,ir,iz,ie,1) > cmax )
     &         cmax = trSUMTAL(ik,ir,iz,ie,1)
           if( trSUMTAL(ik,ir,iz,ie,1) < cmin )
     &         cmin = trSUMTAL(ik,ir,iz,ie,1)

          end do      ! ie loop end
         end do       ! ik loop end
        end do        ! iz loop end
       end do         ! ir loop end

!sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_hetrz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_hetxyz(m,ntf,nfile,
     &                               ne,nd,nz,ny,nx,
     &                               sumfactor,weightRate,sumWR,
     &                               trRES,trSUMTAL,
     &                               resc2SUMTAL,resc3SUMTAL,
     &                               ierr)
!                                                                      *
!     created by T.Miura on 2015/04/27                                 *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------------------------------
      include 'param.inc'
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
      integer           m
      integer           ntf
      integer           nfile

      integer           ne
      integer           nd
      integer           nz
      integer           ny
      integer           nx

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (nd,nx,ny,nz,0:ne,2)
      double precision  trSUMTAL(nd,nx,ny,nz,0:ne,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ix
      integer           iy
      integer           iz

      integer           ik
      integer           ie

      double precision  Xa, sigx
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ix = 1, nx
       do iy = 1, ny
        do iz = 1, nz
         do ik = 1, nd
          do ie = 0, ne

           if( trRES(ik,ix,iy,iz,ie,1) .gt. 0.d0 ) then
            call calc_stdev(m,Xa,sigx,
     &                      trRES(ik,ix,iy,iz,ie,1),
     &                      trRES(ik,ix,iy,iz,ie,2),
     &                      1.0d+0)
            trRES(ik,ix,iy,iz,ie,1) = Xa
            trRES(ik,ix,iy,iz,ie,2) = sigx
           end if

           if( ntf > 1 ) then

            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(ik,ix,iy,iz,ie,1) =
     &      trSUMTAL(ik,ix,iy,iz,ie,1) +
     &         weightRate(ntf)
     &       * trRES(ik,ix,iy,iz,ie,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(ik,ix,iy,iz,ie,2) =
     &      trSUMTAL(ik,ix,iy,iz,ie,2) +
     &      (  trRES(ik,ix,iy,iz,ie,2)
     &       * trRES(ik,ix,iy,iz,ie,1) )**2
     &       * weightRate(ntf)**2

           ! ntf = 1
           else
            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(ik,ix,iy,iz,ie,1) =
     &         weightRate(ntf)
     &       * trRES(ik,ix,iy,iz,ie,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(ik,ix,iy,iz,ie,2) =
     &      (  trRES(ik,ix,iy,iz,ie,2)
     &       * trRES(ik,ix,iy,iz,ie,1) )**2
     &       * weightRate(ntf)**2
           end if

          end do   ! ie loop end
         end do    ! ik loop end
        end do     ! iz loop end
       end do      ! iy loop end
      end do       ! ix loop end

!sumover  tr_sum -> trSUMTAL_sum
      call calc_stdev_tr_sum(m, 1.0d0)
      call sumtal_calc_average_sub(m, ntf,weightRate(ntf))

!-----------------------------------------------------------------------
  100 continue
      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = resc2SUMTAL
     &            + weightRate(ntf) * resc3(m) * resc2(m)/resc3(m)

      ! N = N1 + N2
      resc3SUMTAL = resc3SUMTAL + resc3(m)
!-----------------------------------------------------------------------
      if( ntf == nfile ) then
       do ix = 1, nx
        do iy = 1, ny
         do iz = 1, nz
          do ik = 1, nd
           do ie = 0, ne

            ! X_bar
            trSUMTAL(ik,ix,iy,iz,ie,1) =
     &      trSUMTAL(ik,ix,iy,iz,ie,1) * sumfactor/sumWR

            ! sig_x
            trSUMTAL(ik,ix,iy,iz,ie,2) = sqrt(
     &      trSUMTAL(ik,ix,iy,iz,ie,2)) * sumfactor/sumWR


            ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
            trSUMTAL(ik,ix,iy,iz,ie,2) =
     &     (trSUMTAL(ik,ix,iy,iz,ie,2)**2 *
     &      resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &      resc3SUMTAL *
     &      trSUMTAL(ik,ix,iy,iz,ie,1)**2) *
     &      (resc2SUMTAL/resc3SUMTAL)**2

            ! Sigma xi wi = X_bar W
            trSUMTAL(ik,ix,iy,iz,ie,1) =
     &      trSUMTAL(ik,ix,iy,iz,ie,1) * resc2SUMTAL

            if( trSUMTAL(ik,ix,iy,iz,ie,1) > cmax )
     &            cmax = trSUMTAL(ik,ix,iy,iz,ie,1)
            if( trSUMTAL(ik,ix,iy,iz,ie,1) < cmin )
     &            cmin = trSUMTAL(ik,ix,iy,iz,ie,1)

           end do     ! ie loop end
          end do      ! ik loop end
         end do       ! iz loop end
        end do        ! iy loop end
       end do         ! ix loop end

!sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_hetxyz
