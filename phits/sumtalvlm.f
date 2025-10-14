!***********************************************************************
!                                                                      *
      subroutine sumtal_read_tvlm(m,iax,ntf,ierr)
!                                                                      *
!                                                                      *
!   m: the tally number, index of ital.                                *
!***********************************************************************

        use RESTALMOD
        use sumtallymod
      use moddas_region

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /mpi00/ npe, me

*-----------------------------------------------------------------------

        common /talmm/  nmmax, lmmax, itlmx

*-----------------------------------------------------------------------

        common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)

        common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
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
        call check_tvlm(m,iax,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    ierr)
        endif

        if (ierr .ne. 0 ) goto 900

*-----------------------------------------------------------------------

        rewind(jsi(ioe))
        ill(jsn(ioe),ioe) = 1

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_vlmreg(m,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                    trRES(irestalm(m)),
     &                    itrnv(m),idas(itriv(m)),das(itrrv(m)) )

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_vlmreg(m,itrgn(m),itrgm(m),
     &                    idas_itreg(itreg(m)),trRES(irestalm(m)),
     &                    itrnv(m),idas(itriv(m)),das(itrrv(m)) )

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

      end subroutine sumtal_read_tvlm



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_vlmreg(m,ntf,nr,
     &                             sumfactor,weightRate,
     &                             trRES,trSUMTAL,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
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
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nr,2)
      double precision  trSUMTAL(nr,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------

      integer           ir

      integer,save   :: maxcasd
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

       do ir = 1, nr

           if( ntf > 1 ) then

            ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
            !                 + F r2 sig_N2 xi wi
            trSUMTAL(ir,1) =
     &             trSUMTAL(ir,1) +
     &             ( sumfactor * weightRate(ntf) *
     &             trRES(ir,1) )

            ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
            !                    + F**2 r2**2 sig_N2(xiwi)**2
            trSUMTAL(ir,2) =
     &           trSUMTAL(ir,2) +
     &           ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ir,2) )

           else
            trSUMTAL(ir,1) =
     &             ( sumfactor * weightRate(ntf) *
     &             trRES(ir,1) )

            trSUMTAL(ir,2) =
     &           ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ir,2) )

           end if

           if( trSUMTAL(ir,1) > cmax )
     &          cmax = trSUMTAL(ir,1)
           if( trSUMTAL(ir,1) < cmin )
     &          cmin = trSUMTAL(ir,1)

       end do        ! ir loop end

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
        ErrID = 'L:240/R:sumtal_calc_stdev_vlmreg/F:sumtalvlm.f' !E67_002_001


        write(*,'(/'' ***** Error in sumtal calc stdev vlmreg *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_vlmreg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_vlmreg(m,ntf,nfile,nr,
     &                               sumfactor,weightRate,sumWR,
     &                               trRES,trSUMTAL,
     &                               resc2SUMTAL,resc3SUMTAL,
     &                               ierr)
!                                                                      *
!     created by S.Hashimoto on 2015/1/20                              *
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
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (nr,2)
      double precision  trSUMTAL(nr,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------

      integer           ir

      double precision  Xa, sigx
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

       do ir = 1, nr

           if( trRES(ir,1) .gt. 0.d0 ) then
            call calc_stdev(m,Xa,sigx,
     &             trRES(ir,1),
     &             trRES(ir,2),
     &             1.0d+0)
            trRES(ir,1) = Xa
            trRES(ir,2) = sigx
           end if

           if( ntf > 1 ) then

            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(ir,1) =
     &             trSUMTAL(ir,1) +
     &             weightRate(ntf)
     &             * trRES(ir,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(ir,2) =
     &           trSUMTAL(ir,2) +
     &           (  trRES(ir,2)
     &           * trRES(ir,1) )**2
     &           * weightRate(ntf)**2

           ! ntf = 1
           else
            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(ir,1) =
     &             weightRate(ntf)
     &             * trRES(ir,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(ir,2) =
     &           (  trRES(ir,2)
     &           * trRES(ir,1) )**2
     &           * weightRate(ntf)**2
           end if

       end do        ! ir loop end

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

            ! X_bar
            trSUMTAL(ir,1) =
     &             trSUMTAL(ir,1) * sumfactor/sumWR

            ! sig_x
            trSUMTAL(ir,2) = sqrt(
     &           trSUMTAL(ir,2)) * sumfactor/sumWR


            ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
            trSUMTAL(ir,2) =
     &           (trSUMTAL(ir,2)**2 *
     &           resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &           resc3SUMTAL *
     &           trSUMTAL(ir,1)**2) *
     &           (resc2SUMTAL/resc3SUMTAL)**2

            ! Sigma xi wi = X_bar W
            trSUMTAL(ir,1) =
     &           trSUMTAL(ir,1) * resc2SUMTAL

            if( trSUMTAL(ir,1) > cmax )
     &           cmax = trSUMTAL(ir,1)
            if( trSUMTAL(ir,1) < cmin )
     &           cmin = trSUMTAL(ir,1)

        end do        ! ir loop end

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_vlmreg



