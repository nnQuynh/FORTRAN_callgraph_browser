!***********************************************************************
!                                                                      *
      subroutine sumtal_read_tyield(m,iax,ntf,ierr)
!                                                                      *
!     original routine is read_tyield in restyield.f                   *
!                                                                      *
!   m: the tally number, index of ital.                                *
!***********************************************************************

        use RESTALMOD
        use sumtallymod
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
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                  itnun(itlmax), itnuc(itlmax),
     &                  itndz(itlmax), itndn(itlmax),
     &                  itnkz(itlmax), itnkn(itlmax)
        common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
        common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/02/18

        dimension     idas(mdas*2) ! S.H. 2023.6.13 temporary for fbounds-check
        equivalence ( das, idas )

        integer   ntf   ! Number of tallyfname (tally file name)
*-----------------------------------------------------------------------

        character dsin(0:9,2)*200
        dimension idsi(0:9,2)
        dimension ill(0:9,2), ilf(0:9,2)
        dimension jsn(2), jsi(2)
        dimension ierrs(2)

*-----------------------------------------------------------------------
        common /redufmt/ iredufmt(itlmax) !FURUTA20200615
*-----------------------------------------------------------------------

        ierr = 0
        noe  = 1

        if ( any( itaxs(m,iax) .eq. (/ 8, 9, 10, 11, 12, 13 /) )
     &      .and. ittwo(m) .ne. 4 ) noe = 2

        if(itaxs(m,iax).eq.13.and.iredufmt(m).eq.1) noe=1 !FURUTA20200615

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
        call check_tyield(m,iax,jsn(ioe),jsi(ioe),
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
        idas2 = idas1 + ( maxnt + maxpt ) * 2
        if( itmsh(m) .eq. 1 ) then
          idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
          idas4 = idas2 + itrgn(m)
     &          + ( itrgn(m) + mod(itrgn(m),2) ) / 2
          idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
        else if ( itmsh(m) .eq. 3 ) then
          idas3 = ( idas2 + itnfn(m) * itnfn(m) - 1 ) * 2 + 1
cFURUTA20190121 ! copy from reg
        else if( itmsh(m) .eq. 4 ) then
          idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
          idas4 = idas2 + itrgn(m)
     &          + ( itrgn(m) + mod(itrgn(m),2) ) / 2
          idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
        end if

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_yieldreg(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndz(m),itndn(m),itndm(m),
     &      itrgn(m),itrgm(m),itnun(m),
     &      idas_itreg(itreg(m)),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call read_yieldrz(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndz(m),itndn(m),itndm(m),itnfr(m),itnfz(m),
     &      itrnm(m),itznm(m),itnun(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call read_yieldxyz(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndz(m),itndn(m),itndm(m),itnfn(m),
     &      itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),itnun(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then
          call read_yieldtet(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndz(m),itndn(m),itndm(m),itrgn(m),itnun(m),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_yieldreg(m,
     &      itndz(m),itndn(m),itndm(m),
     &      itrgn(m),itrgm(m),itnun(m),
     &      idas_itreg(itreg(m)),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call restore_yieldrz(m,
     &      itndz(m),itndn(m),itndm(m),itnfr(m),itnfz(m),
     &      itrnm(m),itznm(m),itnun(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call restore_yieldxyz(m,
     &      itndz(m),itndn(m),itndm(m),itnfn(m),
     &      itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),itnun(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then

          call restore_yieldtet(m,
     &      itndz(m),itndn(m),itndm(m),
     &      itrgn(m),itrgm(m),itnun(m),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
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

      end subroutine sumtal_read_tyield


!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_yieldreg(m,ntf,
     &                             mn,mz,mm,nr, ! frtati 2022/02/18 added mm
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

      integer           mn
      integer           mz
      integer           mm
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nr,mz,mn,0:mm,2)
      double precision  trSUMTAL(nr,mz,mn,0:mm,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           im

      integer           iz
      integer           in
      integer           il
      integer           ir

      integer,save   :: maxcasd
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr
       do iz = 1, mz
        do in = 1, mn
         do il = 0, mm ! frtati 2022/02/18

          if( ntf > 1 ) then

           ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
           !                 + F r2 sig_N2 xi wi
           trSUMTAL(ir,iz,in,il,1) =
     &     trSUMTAL(ir,iz,in,il,1) +
     &             ( sumfactor * weightRate(ntf) *
     &        trRES(ir,iz,in,il,1) )

           ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
           !                    + F**2 r2**2 sig_N2(xiwi)**2
           trSUMTAL(ir,iz,in,il,2) =
     &     trSUMTAL(ir,iz,in,il,2) +
     &             ( sumfactor**2 * weightRate(ntf)**2 *
     &        trRES(ir,iz,in,il,2) )

          else
           trSUMTAL(ir,iz,in,il,1) =
     &             ( sumfactor * weightRate(ntf) *
     &        trRES(ir,iz,in,il,1) )

           trSUMTAL(ir,iz,in,il,2) =
     &             ( sumfactor**2 * weightRate(ntf)**2 *
     &        trRES(ir,iz,in,il,2) )

          end if

          if( trSUMTAL(ir,iz,in,il,1) > cmax )
     &         cmax = trSUMTAL(ir,iz,in,il,1)
          if( trSUMTAL(ir,iz,in,il,1) < cmin )
     &         cmin = trSUMTAL(ir,iz,in,il,1)

         end do     ! il loop end
        end do      ! in loop end
       end do       ! iz loop end
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
        ErrID = 'L:386/R:sumtal_calc_stdev_yieldreg/F:sumtalyield.f' !E57_001_001


        write(*,'(/'' ***** Error in sumtal calc stdev yieldreg *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_yieldreg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_yieldrz(m,ntf,
     &                             mn,mz,mm,nz,nr, ! frtati 2022/02/18 added mm
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

      integer           mn
      integer           mz
      integer           mm
      integer           nz
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nr,nz,mz,mn,0:mm,2)
      double precision  trSUMTAL(nr,nz,mz,mn,0:mm,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           il
      integer           in
      integer           iz
      integer           jz
      integer           jr

      integer,save   :: maxcasd
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do jr = 1, nr
       do jz = 1, nz
        do iz = 1, mz
         do in = 1, mn
          do il = 0, mm ! frtati 2022/02/18

           if( ntf > 1 ) then

            ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
            !                 + F r2 sig_N2 xi wi
            trSUMTAL(jr,jz,iz,in,il,1) =
     &      trSUMTAL(jr,jz,iz,in,il,1) +
     &              ( sumfactor * weightRate(ntf) *
     &         trRES(jr,jz,iz,in,il,1) )

            ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
            !                    + F**2 r2**2 sig_N2(xiwi)**2
            trSUMTAL(jr,jz,iz,in,il,2) =
     &      trSUMTAL(jr,jz,iz,in,il,2) +
     &              ( sumfactor**2 * weightRate(ntf)**2 *
     &        trRES(jr,jz,iz,in,il,2) )

           else
            trSUMTAL(jr,jz,iz,in,il,1) =
     &              ( sumfactor * weightRate(ntf) *
     &         trRES(jr,jz,iz,in,il,1) )

            trSUMTAL(jr,jz,iz,in,il,2) =
     &              ( sumfactor**2 * weightRate(ntf)**2 *
     &         trRES(jr,jz,iz,in,il,2) )

           end if

           if( trSUMTAL(jr,jz,iz,in,il,1) > cmax )
     &          cmax = trSUMTAL(jr,jz,iz,in,il,1)
           if( trSUMTAL(jr,jz,iz,in,il,1) < cmin )
     &          cmin = trSUMTAL(jr,jz,iz,in,il,1)

          end do     ! il loop end
         end do      ! in loop end
        end do       ! iz loop end
       end do        ! jz loop end
      end do         ! jr loop end

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
        ErrID = 'L:538/R:sumtal_calc_stdev_yieldrz/F:sumtalyield.f' !E57_001_002


        write(*,'(/'' ***** Error in sumtal calc stdev yieldrz *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_yieldrz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_yieldxyz(m,ntf,
     &                             mn,mz,mm,nz,ny,nx, ! frtati 2022/02/18 added mm
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

      integer           mn
      integer           mz
      integer           mm
      integer           nz
      integer           ny
      integer           nx

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (nx,ny,nz,mz,mn,0:mm,2)
      double precision  trSUMTAL(nx,ny,nz,mz,mn,0:mm,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           jx
      integer           jy
      integer           jz
      integer           iz
      integer           in
      integer           il

      integer,save   :: maxcasd
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
           do il = 0, mm ! frtati 2022/02/18

            if( ntf > 1 ) then

             ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
             !                 + F r2 sig_N2 xi wi
             trSUMTAL(jx,jy,jz,iz,in,il,1) =
     &       trSUMTAL(jx,jy,jz,iz,in,il,1) +
     &               ( sumfactor * weightRate(ntf) *
     &          trRES(jx,jy,jz,iz,in,il,1) )

             ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
             !                    + F**2 r2**2 sig_N2(xiwi)**2
             trSUMTAL(jx,jy,jz,iz,in,il,2) =
     &       trSUMTAL(jx,jy,jz,iz,in,il,2) +
     &               ( sumfactor**2 * weightRate(ntf)**2 *
     &          trRES(jx,jy,jz,iz,in,il,2) )

            else
             trSUMTAL(jx,jy,jz,iz,in,il,1) =
     &               ( sumfactor * weightRate(ntf) *
     &          trRES(jx,jy,jz,iz,in,il,1) )

             trSUMTAL(jx,jy,jz,iz,in,il,2) =
     &               ( sumfactor**2 * weightRate(ntf)**2 *
     &          trRES(jx,jy,jz,iz,in,il,2) )

            end if

            if( trSUMTAL(jx,jy,jz,iz,in,il,1) > cmax )
     &            cmax = trSUMTAL(jx,jy,jz,iz,in,il,1)
            if( trSUMTAL(jx,jy,jz,iz,in,il,1) < cmin )
     &            cmin = trSUMTAL(jx,jy,jz,iz,in,il,1)

           end do    ! il loop end
          end do     ! in loop end
         end do      ! iz loop end
        end do       ! jz loop end
       end do        ! jy loop end
      end do         ! jx loop end

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
        ErrID = 'L:694/R:sumtal_calc_stdev_yieldxyz/F:sumtalyield.f' !E57_001_003


        write(*,'(/'' ***** Error in sumtal calc stdev yieldxyz *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_yieldxyz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_yieldreg(m,ntf,nfile,
     &                               mn,mz,mm,nr, ! frtati 2022/02/18 added mm
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

      integer           mn
      integer           mz
      integer           mm
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (nr,mz,mn,0:mm,2)
      double precision  trSUMTAL(nr,mz,mn,0:mm,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           il
      integer           in
      integer           iz
      integer           ir

      double precision  Xa, sigx
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr
       do iz = 1, mz
        do in = 1, mn
         do il = 0, mm ! frtati 2022/02/18

          if( trRES(ir,iz,in,il,1) .gt. 0.d0 ) then
           call calc_stdev(m,Xa,sigx,
     &                     trRES(ir,iz,in,il,1),
     &                     trRES(ir,iz,in,il,2),
     &                     1.0d+0)
           trRES(ir,iz,in,il,1) = Xa
           trRES(ir,iz,in,il,2) = sigx
          end if

          if( ntf > 1 ) then

           ! X_bar = F Sigma rj/r Xj_bar
           trSUMTAL(ir,iz,in,il,1) =
     &     trSUMTAL(ir,iz,in,il,1) +
     &              weightRate(ntf)
     &      * trRES(ir,iz,in,il,1)

           ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
           trSUMTAL(ir,iz,in,il,2) =
     &     trSUMTAL(ir,iz,in,il,2) +
     &     (  trRES(ir,iz,in,il,2)
     &      * trRES(ir,iz,in,il,1) )**2
     &             * weightRate(ntf)**2

          ! ntf = 1
          else
           ! X_bar = F Sigma rj/r Xj_bar
           trSUMTAL(ir,iz,in,il,1) =
     &              weightRate(ntf)
     &      * trRES(ir,iz,in,il,1)

           ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
           trSUMTAL(ir,iz,in,il,2) =
     &     (  trRES(ir,iz,in,il,2)
     &      * trRES(ir,iz,in,il,1) )**2
     &             * weightRate(ntf)**2
          end if

         end do      ! il loop end
        end do       ! in loop end
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
        do iz = 1, mz
         do in = 1, mn
          do il = 0, mm ! frtati 2022/02/18

           ! X_bar
           trSUMTAL(ir,iz,in,il,1) =
     &     trSUMTAL(ir,iz,in,il,1) * sumfactor/sumWR

           ! sig_x
           trSUMTAL(ir,iz,in,il,2) = sqrt(
     &     trSUMTAL(ir,iz,in,il,2)) * sumfactor/sumWR


           ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
           trSUMTAL(ir,iz,in,il,2) =
     &    (trSUMTAL(ir,iz,in,il,2)**2 *
     &              resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &              resc3SUMTAL *
     &     trSUMTAL(ir,iz,in,il,1)**2) *
     &             (resc2SUMTAL/resc3SUMTAL)**2

           ! Sigma xi wi = X_bar W
           trSUMTAL(ir,iz,in,il,1) =
     &     trSUMTAL(ir,iz,in,il,1) * resc2SUMTAL

           if( trSUMTAL(ir,iz,in,il,1) > cmax )
     &          cmax = trSUMTAL(ir,iz,in,il,1)
           if( trSUMTAL(ir,iz,in,il,1) < cmin )
     &          cmin = trSUMTAL(ir,iz,in,il,1)

          end do      ! il loop end
         end do       ! in loop end
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

      end subroutine sumtal_calc_average_yieldreg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_yieldrz(m,ntf,nfile,
     &                               mn,mz,mm,nz,nr, ! frtati 2022/02/18
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

      integer           mn
      integer           mz
      integer           mm
      integer           nz
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (nr,nz,mz,mn,0:mm,2)
      double precision  trSUMTAL(nr,nz,mz,mn,0:mm,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           jr
      integer           jz
      integer           iz
      integer           in
      integer           il

      double precision  Xa, sigx
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do jr = 1, nr
       do jz = 1, nz
        do iz = 1, mz
         do in = 1, mn
          do il = 0, mm ! frtati 2022/02/18

           if( trRES(jr,jz,iz,in,il,1) .gt. 0.d0 ) then
            call calc_stdev(m,Xa,sigx,
     &                      trRES(jr,jz,iz,in,il,1),
     &                      trRES(jr,jz,iz,in,il,2),
     &                      1.0d+0)
            trRES(jr,jz,iz,in,il,1) = Xa
            trRES(jr,jz,iz,in,il,2) = sigx
           end if

           if( ntf > 1 ) then

            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(jr,jz,iz,in,il,1) =
     &      trSUMTAL(jr,jz,iz,in,il,1) +
     &               weightRate(ntf)
     &       * trRES(jr,jz,iz,in,il,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(jr,jz,iz,in,il,2) =
     &      trSUMTAL(jr,jz,iz,in,il,2) +
     &      (  trRES(jr,jz,iz,in,il,2)
     &       * trRES(jr,jz,iz,in,il,1) )**2
     &              * weightRate(ntf)**2

           ! ntf = 1
           else
            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(jr,jz,iz,in,il,1) =
     &               weightRate(ntf)
     &       * trRES(jr,jz,iz,in,il,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(jr,jz,iz,in,il,2) =
     &      (  trRES(jr,jz,iz,in,il,2)
     &       * trRES(jr,jz,iz,in,il,1) )**2
     &              * weightRate(ntf)**2
           end if

          end do     ! il loop end
         end do      ! in loop end
        end do       ! iz loop end
       end do        ! jz loop end
      end do         ! jr loop end

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
       do jr = 1, nr
        do jz = 1, nz
         do iz = 1, mz
          do in = 1, mn
           do il = 0, mm ! frtati 2022/02/18

            ! X_bar
            trSUMTAL(jr,jz,iz,in,il,1) =
     &      trSUMTAL(jr,jz,iz,in,il,1) * sumfactor/sumWR

            ! sig_x
            trSUMTAL(jr,jz,iz,in,il,2) = sqrt(
     &      trSUMTAL(jr,jz,iz,in,il,2)) * sumfactor/sumWR


            ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
            trSUMTAL(jr,jz,iz,in,il,2) =
     &     (trSUMTAL(jr,jz,iz,in,il,2)**2 *
     &               resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &               resc3SUMTAL *
     &      trSUMTAL(jr,jz,iz,in,il,1)**2) *
     &              (resc2SUMTAL/resc3SUMTAL)**2

            ! Sigma xi wi = X_bar W
            trSUMTAL(jr,jz,iz,in,il,1) =
     &      trSUMTAL(jr,jz,iz,in,il,1) * resc2SUMTAL

            if( trSUMTAL(jr,jz,iz,in,il,1) > cmax )
     &           cmax = trSUMTAL(jr,jz,iz,in,il,1)
            if( trSUMTAL(jr,jz,iz,in,il,1) < cmin )
     &           cmin = trSUMTAL(jr,jz,iz,in,il,1)

           end do     ! il loop end
          end do      ! in loop end
         end do       ! iz loop end
        end do        ! jz loop end
       end do         ! jr loop end

!sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_yieldrz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_yieldxyz(m,ntf,nfile,
     &                               mn,mz,mm,nz,ny,nx, ! frtati 2022/02/18 added mm
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

      integer           mn
      integer           mz
      integer           mm
      integer           nz
      integer           ny
      integer           nx

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (nx,ny,nz,mz,mn,0:mm,2)
      double precision  trSUMTAL(nx,ny,nz,mz,mn,0:mm,2)
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

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do jx = 1, nx
       do jy = 1, ny
        do jz = 1, nz
         do iz = 1, mz
          do in = 1, mn
           do il = 0, mm ! frtati 2022/02/18

            if( trRES(jx,jy,jz,iz,in,il,1) .gt. 0.d0 ) then
             call calc_stdev(m,Xa,sigx,
     &                       trRES(jx,jy,jz,iz,in,il,1),
     &                       trRES(jx,jy,jz,iz,in,il,2),
     &                       1.0d+0)
             trRES(jx,jy,jz,iz,in,il,1) = Xa
             trRES(jx,jy,jz,iz,in,il,2) = sigx
            end if

            if( ntf > 1 ) then

             ! X_bar = F Sigma rj/r Xj_bar
             trSUMTAL(jx,jy,jz,iz,in,il,1) =
     &       trSUMTAL(jx,jy,jz,iz,in,il,1) +
     &                weightRate(ntf)
     &        * trRES(jx,jy,jz,iz,in,il,1)

             ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
             trSUMTAL(jx,jy,jz,iz,in,il,2) =
     &       trSUMTAL(jx,jy,jz,iz,in,il,2) +
     &       (  trRES(jx,jy,jz,iz,in,il,2)
     &        * trRES(jx,jy,jz,iz,in,il,1) )**2
     &               * weightRate(ntf)**2

            ! ntf = 1
            else
             ! X_bar = F Sigma rj/r Xj_bar
             trSUMTAL(jx,jy,jz,iz,in,il,1) =
     &                weightRate(ntf)
     &        * trRES(jx,jy,jz,iz,in,il,1)

             ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
             trSUMTAL(jx,jy,jz,iz,in,il,2) =
     &       (  trRES(jx,jy,jz,iz,in,il,2)
     &        * trRES(jx,jy,jz,iz,in,il,1) )**2
     &               * weightRate(ntf)**2
            end if

           end do    ! il loop end
          end do     ! in loop end
         end do      ! iz loop end
        end do       ! jz loop end
       end do        ! jy loop end
      end do         ! jx loop end

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
       do jx = 1, nx
        do jy = 1, ny
         do jz = 1, nz
          do iz = 1, mz
           do in = 1, mn
            do il = 0, mm ! frtati 2022/02/18

             ! X_bar
             trSUMTAL(jx,jy,jz,iz,in,il,1) =
     &       trSUMTAL(jx,jy,jz,iz,in,il,1) * sumfactor/sumWR

             ! sig_x
             trSUMTAL(jx,jy,jz,iz,in,il,2) = sqrt(
     &       trSUMTAL(jx,jy,jz,iz,in,il,2)) * sumfactor/sumWR


             ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
             trSUMTAL(jx,jy,jz,iz,in,il,2) =
     &      (trSUMTAL(jx,jy,jz,iz,in,il,2)**2 *
     &                resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &                resc3SUMTAL *
     &       trSUMTAL(jx,jy,jz,iz,in,il,1)**2) *
     &               (resc2SUMTAL/resc3SUMTAL)**2

             ! Sigma xi wi = X_bar W
             trSUMTAL(jx,jy,jz,iz,in,il,1) =
     &       trSUMTAL(jx,jy,jz,iz,in,il,1) * resc2SUMTAL

             if( trSUMTAL(jx,jy,jz,iz,in,il,1) > cmax )
     &             cmax = trSUMTAL(jx,jy,jz,iz,in,il,1)
             if( trSUMTAL(jx,jy,jz,iz,in,il,1) < cmin )
     &             cmin = trSUMTAL(jx,jy,jz,iz,in,il,1)

            end do    ! il loop end
           end do     ! in loop end
          end do      ! iz loop end
         end do       ! jz loop end
        end do        ! jy loop end
       end do         ! jx loop end

!sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_yieldxyz
