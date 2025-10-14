!***********************************************************************
!                                                                      *
      subroutine sumtal_read_twwg(m,iax,ntf,ierr)
!                                                                      *
!     original routine is read_twwg in resttrack.f                     *
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
        include 'err.inc'

*-----------------------------------------------------------------------

        common /mpi00/ npe, me

*-----------------------------------------------------------------------

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
        common /tall08/ rtrx0(itlmax), rtry0(itlmax)
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                  rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
        common /tall21/ rtfac(itlmax)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
        common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

        common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
        character crfln*100

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

      character m_err*200

*-----------------------------------------------------------------------
cKN 2018/01/08

        ierr = 0
        noe  = 1

        if(itaxs(m,iax).eq.7)then !FURUTA20240111
         goto 998
        endif

        if ( any( itaxs(m,iax) .eq. (/ 3, 4, 5 /) )
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
        call check_twwg(m,iax,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    ierr)
        endif

        if (ierr .ne. 0 ) goto 900

*-----------------------------------------------------------------------

        rewind(jsi(ioe))
        ill(jsn(ioe),ioe) = 1

*-----------------------------------------------------------------------

        idas1 = lmmax
        idas2 = idas1 + itenm(m)
        idas3 = idas2 + ittnm(m)
        idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_wwgreg(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &                    ittnm(m),
     &                    idas_itreg(itreg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)),
     &                    itrnv(m),idas(itriv(m)),das(itrrv(m)) )

        else if ( itmsh(m) .eq. 3 ) then

          call read_wwgxyz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itmtn(m),ismte(itmtt(m)),
     &                    itxnm(m),itynm(m),itznm(m),itenm(m),itmst(m),
     &                    ittnm(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then
cFURUTA20240110
!<-20211217murofushi update
!--          call read_tractet(m,iax,ioe,jsn(ioe),jsi(ioe),
!--
!--     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
!--     &                    itpan(m),itrgn(m),itenm(m),itmst(m),ittnm(m),
!--     &                    das(iterg(m)),das(ittrg(m)),
!--     &                    trRES(irestalm(m)))
          call read_wwgtet(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itrgn(m),itenm(m),itmst(m),ittnm(m),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))
!--->

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_wwgreg(m,
     &                    itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &                    ittnm(m),
     &                    idas_itreg(itreg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)),
     &                    itrnv(m),idas(itriv(m)),das(itrrv(m)) )

        else if ( itmsh(m) .eq. 3 ) then

          call restore_wwgxyz(m,
     &                    itpan(m),itmtn(m),ismte(itmtt(m)),
     &                    itxnm(m),itynm(m),itznm(m),itenm(m),itmst(m),
     &                    ittnm(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then

          call restore_wwgtet(m,
     &                    itpan(m),itrgn(m),itrgm(m),
     &                    itenm(m),itmst(m),ittnm(m),
     &                    idas_itreg(itreg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

        end if

  900   continue
*-----------------------------------------------------------------------
*   close restart file
*-----------------------------------------------------------------------

        do ioe = 1, noe

          close(jsi(ioe))

        end do

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  998    m_err = '1st axis of T-WWG shold not be wwg for sumtally'
         ErrCha = ''
         ErrID = 'L:236/R:sumtal_read_twwg/F:sumtalwwg.f'
         goto 999

*-----------------------------------------------------------------------
  999   continue
         ierr  = 1

         call ErrWrite(ErrID, ErrCha)
         write(*,*) 'Error: ' // m_err

*-----------------------------------------------------------------------

      return

      end subroutine sumtal_read_twwg


!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_wwgreg(m,ntf,
     &                             np,ne,nt,nr,nm,
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
      integer           np
      integer           ne
      integer           nt
      integer           nr
      integer           nm

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,nt,nr,nm,2)
      double precision  trSUMTAL(np,ne,nt,nr,nm,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           im

      integer           ir

      integer           ip
      integer           ie
      integer           it

      integer,save   :: maxcasd
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do im = 1, nm
       do ir = 1, nr
        do it = 1, nt
         do ie = 1, ne
          do ip = 1, np

           if( ntf > 1 ) then

            ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
            !                 + F r2 sig_N2 xi wi
            trSUMTAL(ip,ie,it,ir,im,1) =
     &             trSUMTAL(ip,ie,it,ir,im,1) +
     &             ( sumfactor * weightRate(ntf) *
     &             trRES(ip,ie,it,ir,im,1) )

            ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
            !                    + F**2 r2**2 sig_N2(xiwi)**2
            trSUMTAL(ip,ie,it,ir,im,2) =
     &           trSUMTAL(ip,ie,it,ir,im,2) +
     &           ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie,it,ir,im,2) )

           else
            trSUMTAL(ip,ie,it,ir,im,1) =
     &             ( sumfactor * weightRate(ntf) *
     &             trRES(ip,ie,it,ir,im,1) )

            trSUMTAL(ip,ie,it,ir,im,2) =
     &           ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie,it,ir,im,2) )

           end if

           if( trSUMTAL(ip,ie,it,ir,im,1) > cmax )
     &          cmax = trSUMTAL(ip,ie,it,ir,im,1)
           if( trSUMTAL(ip,ie,it,ir,im,1) < cmin )
     &          cmin = trSUMTAL(ip,ie,it,ir,im,1)

          end do     ! ip loop end
         end do      ! ie loop end
        end do       ! it loop end
       end do        ! ir loop end
      end do         ! im loop end

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
        ErrID = 'L:389/R:sumtal_calc_stdev_wwgreg/F:sumtalwwg.f' !E65_001_001


        write(*,'(/'' ***** Error in sumtal calc stdev wwgreg *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_wwgreg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_wwgxyz(m,ntf,
     &                             np,ne,nt,nx,ny,nz,nm,
     &                             sumfactor,weightRate,
     &                             trRES,trSUMTAL,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     create by T.Miura on 2014/11/30                                  *
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
      integer           nm

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,nt,nx*ny*nz,nm,2)
      double precision  trSUMTAL(np,ne,nt,nx*ny*nz,nm,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           im

      integer           ix
      integer           iy
      integer           iz

      integer           ip
      integer           ie
      integer           it

      integer,save   :: maxcasd
!-----------------------------------------------------------------------
      integer           icf

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do im = 1, nm
       do iz = 1, nz
        do iy = 1, ny
         do ix = 1, nx
          do it = 1, nt
           do ie = 1, ne
            do ip = 1, np

             if( ntf > 1 ) then

              ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
              !                 + F r2 sig_N2 xi wi
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) =
     &        trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) +
     &          ( sumfactor * weightRate(ntf) *
     &           trRES(ip,ie,it,icf(ix,iy,iz),im,1) )

              ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
              !                    + F**2 r2**2 sig_N2(xiwi)**2
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,2) =
     &        trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,2) +
     &          ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie,it,icf(ix,iy,iz),im,2) )

             else
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) =
     &          ( sumfactor * weightRate(ntf) *
     &           trRES(ip,ie,it,icf(ix,iy,iz),im,1) )

              trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,2) =
     &          ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie,it,icf(ix,iy,iz),im,2) )

             end if

             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) > cmax )
     &             cmax = trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1)
             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) < cmin )
     &             cmin = trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1)

            end do   ! ip loop end
           end do    ! ie loop end
          end do     ! it loop end
         end do      ! ix loop end
        end do       ! iy loop end
       end do        ! iz loop end
      end do         ! im loop end

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
        ErrID = 'L:554/R:sumtal_calc_stdev_wwgxyz/F:sumtalwwg.f' !E65_001_002



        write(*,'(/'' ***** Error in sumtal calc stdev tracxyz *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)

       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_wwgxyz


!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_wwgreg(m,ntf,nfile,
     &                               np,ne,nt,nr,nm,
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
      integer           np
      integer           ne
      integer           nt
      integer           nr
      integer           nm

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (np,ne,nt,nr,nm,2)
      double precision  trSUMTAL(np,ne,nt,nr,nm,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           im

      integer           ir

      integer           ip
      integer           ie
      integer           it
      double precision  Xa, sigx
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do im = 1, nm
       do ir = 1, nr
        do it = 1, nt
         do ie = 1, ne
          do ip = 1, np

           if( trRES(ip,ie,it,ir,im,1) .gt. 0.d0 ) then
            call calc_stdev(m,Xa,sigx,
     &             trRES(ip,ie,it,ir,im,1),
     &             trRES(ip,ie,it,ir,im,2),
     &             1.0d+0)
            trRES(ip,ie,it,ir,im,1) = Xa
            trRES(ip,ie,it,ir,im,2) = sigx
           end if

           if( ntf > 1 ) then

            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(ip,ie,it,ir,im,1) =
     &             trSUMTAL(ip,ie,it,ir,im,1) +
     &             weightRate(ntf)
     &             * trRES(ip,ie,it,ir,im,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(ip,ie,it,ir,im,2) =
     &           trSUMTAL(ip,ie,it,ir,im,2) +
     &           (  trRES(ip,ie,it,ir,im,2)
     &           * trRES(ip,ie,it,ir,im,1) )**2
     &           * weightRate(ntf)**2

           ! ntf = 1
           else
            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(ip,ie,it,ir,im,1) =
     &             weightRate(ntf)
     &             * trRES(ip,ie,it,ir,im,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(ip,ie,it,ir,im,2) =
     &           (  trRES(ip,ie,it,ir,im,2)
     &           * trRES(ip,ie,it,ir,im,1) )**2
     &           * weightRate(ntf)**2
           end if

          end do     ! ip loop end
         end do      ! ie loop end
        end do       ! it loop end
       end do        ! ir loop end
      end do         ! im loop end

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
       do im = 1, nm
        do ir = 1, nr
         do it = 1, nt
          do ie = 1, ne
           do ip = 1, np

            ! X_bar
            trSUMTAL(ip,ie,it,ir,im,1) =
     &             trSUMTAL(ip,ie,it,ir,im,1) * sumfactor/sumWR

            ! sig_x
            trSUMTAL(ip,ie,it,ir,im,2) = sqrt(
     &           trSUMTAL(ip,ie,it,ir,im,2)) * sumfactor/sumWR


            ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
            trSUMTAL(ip,ie,it,ir,im,2) =
     &           (trSUMTAL(ip,ie,it,ir,im,2)**2 *
     &           resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &           resc3SUMTAL *
     &           trSUMTAL(ip,ie,it,ir,im,1)**2) *
     &           (resc2SUMTAL/resc3SUMTAL)**2

            ! Sigma xi wi = X_bar W
            trSUMTAL(ip,ie,it,ir,im,1) =
     &           trSUMTAL(ip,ie,it,ir,im,1) * resc2SUMTAL

            if( trSUMTAL(ip,ie,it,ir,im,1) > cmax )
     &           cmax = trSUMTAL(ip,ie,it,ir,im,1)
            if( trSUMTAL(ip,ie,it,ir,im,1) < cmin )
     &           cmin = trSUMTAL(ip,ie,it,ir,im,1)

           end do     ! ip loop end
          end do      ! ie loop end
         end do       ! it loop end
        end do        ! ir loop end
       end do         ! im loop end

!sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_wwgreg


!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_wwgxyz(m,ntf,nfile,
     &                               np,ne,nt,nx,ny,nz,nm,
     &                               sumfactor,weightRate,sumWR,
     &                               trRES,trSUMTAL,
     &                               resc2SUMTAL,resc3SUMTAL,
     &                               ierr)
!                                                                      *
!     create by T.Miura on 2014/11/30                                  *
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
      integer           np
      integer           ne
      integer           nt
      integer           nx
      integer           ny
      integer           nz
      integer           nm

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (np,ne,nt,nx*ny*nz,nm,2)
      double precision  trSUMTAL(np,ne,nt,nx*ny*nz,nm,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           im

      integer           ix
      integer           iy
      integer           iz

      integer           ip
      integer           ie
      integer           it
      double precision  Xa, sigx
!-----------------------------------------------------------------------
      integer           icf

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do im = 1, nm
       do iz = 1, nz
        do iy = 1, ny
         do ix = 1, nx
          do it = 1, nt
           do ie = 1, ne
            do ip = 1, np

             if( trRES(ip,ie,it,icf(ix,iy,iz),im,1) .gt. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(ip,ie,it,icf(ix,iy,iz),im,1),
     &                        trRES(ip,ie,it,icf(ix,iy,iz),im,2),
     &                        1.0d+0)
              trRES(ip,ie,it,icf(ix,iy,iz),im,1) = Xa
              trRES(ip,ie,it,icf(ix,iy,iz),im,2) = sigx
             end if

             if( ntf > 1 ) then

              ! X_bar = F Sigma rj/r Xj_bar
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) =
     &        trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) +
     &           weightRate(ntf)
     &         * trRES(ip,ie,it,icf(ix,iy,iz),im,1)

              ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,2) =
     &        trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,2) +
     &        (  trRES(ip,ie,it,icf(ix,iy,iz),im,2)
     &         * trRES(ip,ie,it,icf(ix,iy,iz),im,1) )**2
     &         * weightRate(ntf)**2

             ! ntf = 1
             else
              ! X_bar = F Sigma rj/r Xj_bar
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) =
     &           weightRate(ntf)
     &         * trRES(ip,ie,it,icf(ix,iy,iz),im,1)

              ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,2) =
     &        (  trRES(ip,ie,it,icf(ix,iy,iz),im,2)
     &         * trRES(ip,ie,it,icf(ix,iy,iz),im,1) )**2
     &         * weightRate(ntf)**2
             end if

            end do   ! ip loop end
           end do    ! ie loop end
          end do     ! it loop end
         end do      ! ix loop end
        end do       ! iy loop end
       end do        ! iz loop end
      end do         ! im loop end

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
       do im = 1, nm
        do iz = 1, nz
         do iy = 1, ny
          do ix = 1, nx
           do it = 1, nt
            do ie = 1, ne
             do ip = 1, np

              ! X_bar
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) =
     &        trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) * sumfactor/sumWR

              ! sig_x
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,2) = sqrt(
     &        trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,2)) * sumfactor/sumWR


              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,2) =
     &       (trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,2)**2 *
     &        resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &        resc3SUMTAL *
     &        trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1)**2) *
     &        (resc2SUMTAL/resc3SUMTAL)**2

              ! Sigma xi wi = X_bar W
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) =
     &        trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) * resc2SUMTAL

              if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) > cmax )
     &              cmax = trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1)
              if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1) < cmin )
     &              cmin = trSUMTAL(ip,ie,it,icf(ix,iy,iz),im,1)

             end do   ! ip loop end
            end do    ! ie loop end
           end do     ! it loop end
          end do      ! ix loop end
         end do       ! iy loop end
        end do        ! iz loop end
       end do         ! im loop end

!sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_wwgxyz

