!***********************************************************************
!                                                                      *
      subroutine sumtal_read_tdpa(m,iax,ntf,ierr)
!                                                                      *
!     original routine is read_tdpa in restdpa.f                       *
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
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

        dimension     idas(mdas*2)
        equivalence ( das, idas )

        integer   ntf   ! Number of tallyfname (tally file name)
*-----------------------------------------------------------------------

        character dsin(0:9,2)*200
        dimension idsi(0:9,2)
        dimension ill(0:9,2), ilf(0:9,2)
        dimension jsn(2), jsi(2)
        dimension ierrs(2)

*-----------------------------------------------------------------------

        ierr = 0
        noe  = 1

        if ( any( itaxs(m,iax) .eq. (/ 6, 7, 8, 9 /) )
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
        call check_tdpa(m,iax,jsn(ioe),jsi(ioe),
     &                 dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                 ierr)
        endif

        if (ierr .ne. 0 ) goto 900

*-----------------------------------------------------------------------

        rewind(jsi(ioe))
        ill(jsn(ioe),ioe) = 1

*-----------------------------------------------------------------------

        if (itmsh(m) .eq. 1 ) then
          idas0 = nmmax
          idas1 = lmmax
          idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
          idas3 = idas1 + itrgn(m)
     &          + ( itrgn(m) + mod(itrgn(m),2) ) / 2
          idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
          idasa = idas3
        else if (itmsh(m) .eq. 2 ) then
          idas0 = nmmax
          idas1 = ( lmmax - 1 ) * 2 + 1
          idasa = lmmax
        else if (itmsh(m) .eq. 3 ) then
          idas0 = nmmax
          idasa = lmmax
cFURUTA20190121 ! copy from reg
        else if (itmsh(m) .eq. 4 ) then
          idas0 = nmmax
          idas1 = lmmax
          idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
          idas3 = idas1 + itrgn(m)
     &          + ( itrgn(m) + mod(itrgn(m),2) ) / 2
          idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
          idasa = idas3
        end if

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_dpareg(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),
     &      itpan(m),itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call read_dparz(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),
     &      itpan(m),itrnm(m),itznm(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call read_dpaxyz(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),
     &      itpan(m),itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then
          call read_dpatet(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),itpan(m),itrgn(m),
     &      trRES(irestalm(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_dpareg(m,
     &      itndy(m),itout(m),
     &      itpan(m),itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call restore_dparz(m,
     &      itndy(m),itout(m),
     &      itpan(m),itrnm(m),itznm(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call restore_dpaxyz(m,
     &      itndy(m),itout(m),
     &      itpan(m),itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then

          call restore_dpatet(m,
     &      itndy(m),itout(m),itpan(m),itrgn(m),itrgm(m),
     &      trRES(irestalm(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      idas_itreg(itreg(m)),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

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

      end subroutine sumtal_read_tdpa


!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_dpareg(m,ntf,
     &                             nd,np,nr,
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

      integer           nd
      integer           np
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nd,np,nr,2)
      double precision  trSUMTAL(nd,np,nr,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir
      integer           ip
      integer           ik

      integer,save   :: maxcasd
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr
       do ip = 1, np
        do ik = 1, nd

         if( ntf > 1 ) then

          ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
          !                 + F r2 sig_N2 xi wi
          trSUMTAL(ik,ip,ir,1) =
     &    trSUMTAL(ik,ip,ir,1) +
     &            ( sumfactor * weightRate(ntf) *
     &       trRES(ik,ip,ir,1) )

          ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
          !                    + F**2 r2**2 sig_N2(xiwi)**2
          trSUMTAL(ik,ip,ir,2) =
     &    trSUMTAL(ik,ip,ir,2) +
     &            ( sumfactor**2 * weightRate(ntf)**2 *
     &       trRES(ik,ip,ir,2) )

         else
          trSUMTAL(ik,ip,ir,1) =
     &            ( sumfactor * weightRate(ntf) *
     &       trRES(ik,ip,ir,1) )

          trSUMTAL(ik,ip,ir,2) =
     &            ( sumfactor**2 * weightRate(ntf)**2 *
     &       trRES(ik,ip,ir,2) )

         end if

         if( trSUMTAL(ik,ip,ir,1) > cmax )
     &        cmax = trSUMTAL(ik,ip,ir,1)
         if( trSUMTAL(ik,ip,ir,1) < cmin )
     &        cmin = trSUMTAL(ik,ip,ir,1)

        end do       ! ik loop end
       end do        ! ip loop end
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
        ErrID = 'L:371/R:sumtal_calc_stdev_dpareg/F:sumtaldpa.f' !E59_001_001


        write(*,'(/'' ***** Error in sumtal calc stdev dpareg *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_dpareg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_dparz(m,ntf,
     &                             nd,np,nr,nz,
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

      integer           nd
      integer           np
      integer           nr
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nd,np,nr,nz,2)
      double precision  trSUMTAL(nd,np,nr,nz,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           iz
      integer           ir
      integer           ip
      integer           ik

      integer,save   :: maxcasd
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do iz = 1, nz
       do ir = 1, nr
        do ip = 1, np
         do ik = 1, nd

          if( ntf > 1 ) then

           ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
           !                 + F r2 sig_N2 xi wi
           trSUMTAL(ik,ip,ir,iz,1) =
     &     trSUMTAL(ik,ip,ir,iz,1) +
     &             ( sumfactor * weightRate(ntf) *
     &        trRES(ik,ip,ir,iz,1) )

           ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
           !                    + F**2 r2**2 sig_N2(xiwi)**2
           trSUMTAL(ik,ip,ir,iz,2) =
     &     trSUMTAL(ik,ip,ir,iz,2) +
     &             ( sumfactor**2 * weightRate(ntf)**2 *
     &        trRES(ik,ip,ir,iz,2) )

          else
           trSUMTAL(ik,ip,ir,iz,1) =
     &             ( sumfactor * weightRate(ntf) *
     &        trRES(ik,ip,ir,iz,1) )

           trSUMTAL(ik,ip,ir,iz,2) =
     &             ( sumfactor**2 * weightRate(ntf)**2 *
     &        trRES(ik,ip,ir,iz,2) )

          end if

          if( trSUMTAL(ik,ip,ir,iz,1) > cmax )
     &         cmax = trSUMTAL(ik,ip,ir,iz,1)
          if( trSUMTAL(ik,ip,ir,iz,1) < cmin )
     &         cmin = trSUMTAL(ik,ip,ir,iz,1)

         end do     ! ik loop end
        end do      ! ip loop end
       end do       ! ir loop end
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
        ErrID = 'L:519/R:sumtal_calc_stdev_dparz/F:sumtaldpa.f' !E59_001_002


        write(*,'(/'' ***** Error in sumtal calc stdev dparz *****''/)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_dparz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_dpaxyz(m,ntf,
     &                             nd,np,nx,ny,nz,
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

      integer           nd
      integer           np
      integer           nx
      integer           ny
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nd,np,nx,ny,nz,2)
      double precision  trSUMTAL(nd,np,nx,ny,nz,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           iz
      integer           iy
      integer           ix
      integer           ip
      integer           ik

      integer,save   :: maxcasd
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do iz = 1, nz
       do iy = 1, ny
        do ix = 1, nx
         do ip = 1, np
          do ik = 1, nd

           if( ntf > 1 ) then

            ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
            !                 + F r2 sig_N2 xi wi
            trSUMTAL(ik,ip,ix,iy,iz,1) =
     &      trSUMTAL(ik,ip,ix,iy,iz,1) +
     &              ( sumfactor * weightRate(ntf) *
     &         trRES(ik,ip,ix,iy,iz,1) )

            ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
            !                    + F**2 r2**2 sig_N2(xiwi)**2
            trSUMTAL(ik,ip,ix,iy,iz,2) =
     &      trSUMTAL(ik,ip,ix,iy,iz,2) +
     &              ( sumfactor**2 * weightRate(ntf)**2 *
     &         trRES(ik,ip,ix,iy,iz,2) )

           else
            trSUMTAL(ik,ip,ix,iy,iz,1) =
     &              ( sumfactor * weightRate(ntf) *
     &         trRES(ik,ip,ix,iy,iz,1) )

            trSUMTAL(ik,ip,ix,iy,iz,2) =
     &              ( sumfactor**2 * weightRate(ntf)**2 *
     &         trRES(ik,ip,ix,iy,iz,2) )

           end if

           if( trSUMTAL(ik,ip,ix,iy,iz,1) > cmax )
     &           cmax = trSUMTAL(ik,ip,ix,iy,iz,1)
           if( trSUMTAL(ik,ip,ix,iy,iz,1) < cmin )
     &           cmin = trSUMTAL(ik,ip,ix,iy,iz,1)

          end do    ! ik loop end
         end do     ! ip loop end
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
        ErrID = 'L:671/R:sumtal_calc_stdev_dpaxyz/F:sumtaldpa.f' !E59_001_003


        write(*,'(/'' ***** Error in sumtal calc stdev dpaxyz *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_dpaxyz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_dpareg(m,ntf,nfile,
     &                               nd,np,nr,
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

      integer           nd
      integer           np
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (nd,np,nr,2)
      double precision  trSUMTAL(nd,np,nr,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir
      integer           ip
      integer           ik

      double precision  Xa, sigx
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr
       do ip = 1, np
        do ik = 1, nd

         if( trRES(ik,ip,ir,1) .gt. 0.d0 ) then
          call calc_stdev(m,Xa,sigx,
     &                    trRES(ik,ip,ir,1),
     &                    trRES(ik,ip,ir,2),
     &                    1.0d+0)
          trRES(ik,ip,ir,1) = Xa
          trRES(ik,ip,ir,2) = sigx
         end if

         if( ntf > 1 ) then

          ! X_bar = F Sigma rj/r Xj_bar
          trSUMTAL(ik,ip,ir,1) =
     &    trSUMTAL(ik,ip,ir,1) +
     &             weightRate(ntf)
     &     * trRES(ik,ip,ir,1)

          ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
          trSUMTAL(ik,ip,ir,2) =
     &    trSUMTAL(ik,ip,ir,2) +
     &    (  trRES(ik,ip,ir,2)
     &     * trRES(ik,ip,ir,1) )**2
     &           * weightRate(ntf)**2

         ! ntf = 1
         else
          ! X_bar = F Sigma rj/r Xj_bar
          trSUMTAL(ik,ip,ir,1) =
     &             weightRate(ntf)
     &     * trRES(ik,ip,ir,1)

          ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
          trSUMTAL(ik,ip,ir,2) =
     &    (  trRES(ik,ip,ir,2)
     &     * trRES(ik,ip,ir,1) )**2
     &           * weightRate(ntf)**2
         end if

        end do       ! ik loop end
       end do        ! ip loop end
      end do         ! ir loop end

! sumover  tr_sum -> trSUMTAL_sum
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
        do ip = 1, np
         do ik = 1, nd

          ! X_bar
          trSUMTAL(ik,ip,ir,1) =
     &    trSUMTAL(ik,ip,ir,1) * sumfactor/sumWR

          ! sig_x
          trSUMTAL(ik,ip,ir,2) = sqrt(
     &    trSUMTAL(ik,ip,ir,2)) * sumfactor/sumWR


          ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
          trSUMTAL(ik,ip,ir,2) =
     &   (trSUMTAL(ik,ip,ir,2)**2 *
     &             resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &             resc3SUMTAL *
     &    trSUMTAL(ik,ip,ir,1)**2) *
     &            (resc2SUMTAL/resc3SUMTAL)**2

          ! Sigma xi wi = X_bar W
          trSUMTAL(ik,ip,ir,1) =
     &    trSUMTAL(ik,ip,ir,1) * resc2SUMTAL

          if( trSUMTAL(ik,ip,ir,1) > cmax )
     &         cmax = trSUMTAL(ik,ip,ir,1)
          if( trSUMTAL(ik,ip,ir,1) < cmin )
     &         cmin = trSUMTAL(ik,ip,ir,1)

         end do       ! ik loop end
        end do        ! ip loop end
       end do         ! ir loop end

! sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_dpareg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_dparz(m,ntf,nfile,
     &                               nd,np,nr,nz,
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

      integer           nd
      integer           np
      integer           nr
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (nd,np,nr,nz,2)
      double precision  trSUMTAL(nd,np,nr,nz,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           iz
      integer           ir
      integer           ip
      integer           ik

      double precision  Xa, sigx
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do iz = 1, nz
       do ir = 1, nr
        do ip = 1, np
         do ik = 1, nd

          if( trRES(ik,ip,ir,iz,1) .gt. 0.d0 ) then
           call calc_stdev(m,Xa,sigx,
     &                     trRES(ik,ip,ir,iz,1),
     &                     trRES(ik,ip,ir,iz,2),
     &                     1.0d+0)
           trRES(ik,ip,ir,iz,1) = Xa
           trRES(ik,ip,ir,iz,2) = sigx
          end if

          if( ntf > 1 ) then

           ! X_bar = F Sigma rj/r Xj_bar
           trSUMTAL(ik,ip,ir,iz,1) =
     &     trSUMTAL(ik,ip,ir,iz,1) +
     &              weightRate(ntf)
     &      * trRES(ik,ip,ir,iz,1)

           ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
           trSUMTAL(ik,ip,ir,iz,2) =
     &     trSUMTAL(ik,ip,ir,iz,2) +
     &     (  trRES(ik,ip,ir,iz,2)
     &      * trRES(ik,ip,ir,iz,1) )**2
     &            * weightRate(ntf)**2

          ! ntf = 1
          else
           ! X_bar = F Sigma rj/r Xj_bar
           trSUMTAL(ik,ip,ir,iz,1) =
     &              weightRate(ntf)
     &      * trRES(ik,ip,ir,iz,1)

           ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
           trSUMTAL(ik,ip,ir,iz,2) =
     &     (  trRES(ik,ip,ir,iz,2)
     &      * trRES(ik,ip,ir,iz,1) )**2
     &            * weightRate(ntf)**2
          end if

         end do      ! ik loop end
        end do       ! ip loop end
       end do        ! ir loop end
      end do         ! iz loop end

! sumover  tr_sum -> trSUMTAL_sum
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
       do iz = 1, nz
        do ir = 1, nr
         do ip = 1, np
          do ik = 1, nd

           ! X_bar
           trSUMTAL(ik,ip,ir,iz,1) =
     &     trSUMTAL(ik,ip,ir,iz,1) * sumfactor/sumWR

           ! sig_x
           trSUMTAL(ik,ip,ir,iz,2) = sqrt(
     &     trSUMTAL(ik,ip,ir,iz,2)) * sumfactor/sumWR


           ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
           trSUMTAL(ik,ip,ir,iz,2) =
     &    (trSUMTAL(ik,ip,ir,iz,2)**2 *
     &              resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &              resc3SUMTAL *
     &     trSUMTAL(ik,ip,ir,iz,1)**2) *
     &             (resc2SUMTAL/resc3SUMTAL)**2

           ! Sigma xi wi = X_bar W
           trSUMTAL(ik,ip,ir,iz,1) =
     &     trSUMTAL(ik,ip,ir,iz,1) * resc2SUMTAL

           if( trSUMTAL(ik,ip,ir,iz,1) > cmax )
     &          cmax = trSUMTAL(ik,ip,ir,iz,1)
           if( trSUMTAL(ik,ip,ir,iz,1) < cmin )
     &          cmin = trSUMTAL(ik,ip,ir,iz,1)

          end do      ! ik loop end
         end do       ! ip loop end
        end do        ! ir loop end
       end do         ! iz loop end

! sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_dparz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_dpaxyz(m,ntf,nfile,
     &                               nd,np,nx,ny,nz,
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

      integer           nd
      integer           np
      integer           nx
      integer           ny
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (nd,np,nx,ny,nz,2)
      double precision  trSUMTAL(nd,np,nx,ny,nz,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           iz
      integer           iy
      integer           ix
      integer           ip
      integer           ik

      double precision  Xa, sigx
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do iz = 1, nz
       do iy = 1, ny
        do ix = 1, nx
         do ip = 1, np
          do ik = 1, nd

           if( trRES(ik,ip,ix,iy,iz,1) .gt. 0.d0 ) then
            call calc_stdev(m,Xa,sigx,
     &                      trRES(ik,ip,ix,iy,iz,1),
     &                      trRES(ik,ip,ix,iy,iz,2),
     &                      1.0d+0)
            trRES(ik,ip,ix,iy,iz,1) = Xa
            trRES(ik,ip,ix,iy,iz,2) = sigx
           end if

           if( ntf > 1 ) then

            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(ik,ip,ix,iy,iz,1) =
     &      trSUMTAL(ik,ip,ix,iy,iz,1) +
     &               weightRate(ntf)
     &       * trRES(ik,ip,ix,iy,iz,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(ik,ip,ix,iy,iz,2) =
     &      trSUMTAL(ik,ip,ix,iy,iz,2) +
     &      (  trRES(ik,ip,ix,iy,iz,2)
     &       * trRES(ik,ip,ix,iy,iz,1) )**2
     &             * weightRate(ntf)**2

           ! ntf = 1
           else
            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(ik,ip,ix,iy,iz,1) =
     &               weightRate(ntf)
     &       * trRES(ik,ip,ix,iy,iz,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(ik,ip,ix,iy,iz,2) =
     &      (  trRES(ik,ip,ix,iy,iz,2)
     &       * trRES(ik,ip,ix,iy,iz,1) )**2
     &             * weightRate(ntf)**2
           end if

          end do     ! ik loop end
         end do      ! ip loop end
        end do       ! ix loop end
       end do        ! iy loop end
      end do         ! iz loop end

! sumover  tr_sum -> trSUMTAL_sum
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
       do iz = 1, nz
        do iy = 1, ny
         do ix = 1, nx
          do ip = 1, np
           do ik = 1, nd

            ! X_bar
            trSUMTAL(ik,ip,ix,iy,iz,1) =
     &      trSUMTAL(ik,ip,ix,iy,iz,1) * sumfactor/sumWR

            ! sig_x
            trSUMTAL(ik,ip,ix,iy,iz,2) = sqrt(
     &      trSUMTAL(ik,ip,ix,iy,iz,2)) * sumfactor/sumWR


            ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
            trSUMTAL(ik,ip,ix,iy,iz,2) =
     &     (trSUMTAL(ik,ip,ix,iy,iz,2)**2 *
     &               resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &               resc3SUMTAL *
     &      trSUMTAL(ik,ip,ix,iy,iz,1)**2) *
     &              (resc2SUMTAL/resc3SUMTAL)**2

            ! Sigma xi wi = X_bar W
            trSUMTAL(ik,ip,ix,iy,iz,1) =
     &      trSUMTAL(ik,ip,ix,iy,iz,1) * resc2SUMTAL

            if( trSUMTAL(ik,ip,ix,iy,iz,1) > cmax )
     &            cmax = trSUMTAL(ik,ip,ix,iy,iz,1)
            if( trSUMTAL(ik,ip,ix,iy,iz,1) < cmin )
     &            cmin = trSUMTAL(ik,ip,ix,iy,iz,1)

           end do     ! ik loop end
          end do      ! ip loop end
         end do       ! ix loop end
        end do        ! iy loop end
       end do         ! iz loop end

! sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)


       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_dpaxyz
