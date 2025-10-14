!***********************************************************************
!                                                                      *
      subroutine sumtal_read_tcross(m,iax,ntf,ierr)
!                                                                      *
!     original routine is read_tcross in restcross.f                   *
!                                                                      *
!   m: the tally number, index of ital.                                *
!***********************************************************************

        use RESTALMOD
        use sumtallymod
        use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas_mesh
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

        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
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
        common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)
        common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                  rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)
        common /tall26/ itrcs(itlmax), itrcc(itlmax), itrss(itlmax)

*-----------------------------------------------------------------------
      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall73/ itenclo(itlmax), itangform(itlmax)

*-----------------------------------------------------------------------

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

        ierr = 0
        noe  = 1

        if ( any( itaxs(m,iax) .eq. (/ 7 /) )
     &      .and. ittwo(m) .ne. 4 ) noe = 2

*-----------------------------------------------------------------------
*   open restart file
*-----------------------------------------------------------------------
        call sumtal_open_resfile
     &               (m,noe,ntf,jsn,jsi,dsin,idsi,ill,ilf,newtall,ierr)

        if ( newtall .ne. 0 ) goto 900  !! it's new tally
        if ( ierr    .ne. 0 ) goto 900

*-----------------------------------------------------------------------
*   check tally
*-----------------------------------------------------------------------
        do 800 ioe = 1, noe

        if(ireschk.eq.0) then ! T.Sato 2013/10/19
        idsi(0,ioe) = ltallyfname(ntf,m)
        dsin(0,ioe)(1:idsi(0,ioe))=
     &                         tallyfname(ntf,m)(1:ltallyfname(ntf,m))
        call check_tcross(m,iax,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    ierr)
        endif

        if (ierr .ne. 0 ) goto 900

*-----------------------------------------------------------------------

        rewind(jsi(ioe))
        ill(jsn(ioe),ioe) = 1

*-----------------------------------------------------------------------

        idas0 = nmmax
        idas1 = lmmax
        idas2 = idas1 + itenm(m)
        idas3 = idas2 + itanm(m)
        idas4 = ( idas3 + ittnm(m) - 1 ) * 2 + 1
        idasa = idas3 + ittnm(m)

        idasz = itenm(m) * itpan(m) * 2
     &        * ( itrnm(m) + 1 ) * itznm(m)
     &        * itanm(m)
     &        * ittnm(m)
     &        * itmst(m)

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_crsreg(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itrcn(m),itenm(m),itanm(m),ittnm(m),
     &                    itmst(m),
     &                    itrss(m),
     &                    idas_itrcc(itrcc(m)),das_itrca(itrca(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 2 ) then

         if( itenclo(m) .eq. 1 ) then

          call read_crsrz_rcc(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itrnm(m),itznm(m),itenm(m),itanm(m),
     &                    ittnm(m),itmst(m),
     &                    das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))
!     &                    trRES(irestalm(m)+idasz))

         else

          call read_crsrz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itrnm(m),itznm(m),itenm(m),itanm(m),
     &                    ittnm(m),itmst(m),
     &                    das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)),
     &                    trRES(irestalm(m)+idasz))

         endif

        else if ( itmsh(m) .eq. 3 ) then

         if( itenclo(m) .eq. 1 ) then

          call read_crsxyz_rpp(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itxnm(m),itynm(m),itznm(m),
     &                    itenm(m),itanm(m),ittnm(m),itmst(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

         else

          call read_crsxyz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itxnm(m),itynm(m),itznm(m),
     &                    itenm(m),itanm(m),ittnm(m),itmst(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

         endif

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_crsreg(m,
     &                    itpan(m),itrcn(m),itenm(m),itanm(m),ittnm(m),
     &                    itmst(m),
     &                    itrss(m),
     &                    idas_itrcc(itrcc(m)),das_itrca(itrca(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 2 ) then

         if( itenclo(m) .eq. 1 ) then

          call restore_crsrz_rcc(m,
     &                    itpan(m),itrnm(m),itznm(m),itenm(m),itanm(m),
     &                    ittnm(m),itmst(m),
     &                    das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)),
     &                    trRES(irestalm(m)+idasz))

         else

          call restore_crsrz(m,
     &                    itpan(m),itrnm(m),itznm(m),itenm(m),itanm(m),
     &                    ittnm(m),itmst(m),
     &                    das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)),
     &                    trRES(irestalm(m)+idasz))

         endif

        else if ( itmsh(m) .eq. 3 ) then

         if( itenclo(m) .eq. 1 ) then

          call restore_crsxyz_rpp(m,
     &                    itpan(m),itxnm(m),itynm(m),itznm(m),
     &                    itenm(m),itanm(m),ittnm(m),itmst(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

         else

          call restore_crsxyz(m,
     &                    itpan(m),itxnm(m),itynm(m),itznm(m),
     &                    itenm(m),itanm(m),ittnm(m),itmst(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

         endif

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

      end subroutine sumtal_read_tcross


!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_crsreg(m,ntf,
     &                             np,ne,nt,na,nr,nm,
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

      integer           np
      integer           ne
      integer           nt
      integer           na
      integer           nr
      integer           nm

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,na,nt,nr,nm,2)
      double precision  trSUMTAL(np,ne,na,nt,nr,nm,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           ip
      integer           ie
      integer           it
      integer           ia
      integer           ir
      integer           im

      integer,save   :: maxcasd
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

           if( ntf > 1 ) then

            ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
            !                 + F r2 sig_N2 xi wi
            trSUMTAL(ip,ie,ia,it,ir,im,1) =
     &      trSUMTAL(ip,ie,ia,it,ir,im,1) +
     &              ( sumfactor * weightRate(ntf) *
     &         trRES(ip,ie,ia,it,ir,im,1) )

            ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
            !                    + F**2 r2**2 sig_N2(xiwi)**2
            trSUMTAL(ip,ie,ia,it,ir,im,2) =
     &      trSUMTAL(ip,ie,ia,it,ir,im,2) +
     &              ( sumfactor**2 * weightRate(ntf)**2 *
     &         trRES(ip,ie,ia,it,ir,im,2) )

           else
            trSUMTAL(ip,ie,ia,it,ir,im,1) =
     &              ( sumfactor * weightRate(ntf) *
     &         trRES(ip,ie,ia,it,ir,im,1) )

            trSUMTAL(ip,ie,ia,it,ir,im,2) =
     &              ( sumfactor**2 * weightRate(ntf)**2 *
     &         trRES(ip,ie,ia,it,ir,im,2) )

           end if

           if( trSUMTAL(ip,ie,ia,it,ir,im,1) > cmax )
     &          cmax = trSUMTAL(ip,ie,ia,it,ir,im,1)
           if( trSUMTAL(ip,ie,ia,it,ir,im,1) < cmin )
     &          cmin = trSUMTAL(ip,ie,ia,it,ir,im,1)

          end do     ! ip loop end
         end do      ! ie loop end
        end do       ! it loop end
       end do        ! ia loop end
      end do         ! ir loop end
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
        ErrID = 'L:434/R:sumtal_calc_stdev_crsreg/F:sumtalcross.f' !E52_002_001


        write(*,'('' **** Error in sumtal calc stdev crsreg ****''
     & /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_crsreg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_crsrz(m,ntf,
     &                             np,ne,na,nt,nr,nz,nm,
     &                             sumfactor,weightRate,
     &                             trRES,trSUMTAL,
     &                             tzRES,tzSUMTAL,
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
      double precision  trSUMTAL(np,ne,na,nt,(nr+1)*nz,nm,2)
      double precision  tzRES   (np,ne,na,nt,nr*(nz+1),nm,2)
      double precision  tzSUMTAL(np,ne,na,nt,nr*(nz+1),nm,2)


      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           iz
      integer           ir
      integer           it
      integer           ia
      integer           ie
      integer           ip
      integer           im

      integer,save   :: maxcasd

*-----------------------------------------------------------------------
      integer           irf
      integer           izf

         irf(ir,iz) = ir + ( iz - 1 ) * ( nr + 1 )
         izf(ir,iz) = iz + ( ir - 1 ) * ( nz + 1 )

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

            if( ntf > 1 ) then

             ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
             !                 + F r2 sig_N2 xi wi
             tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) =
     &       tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) +
     &               ( sumfactor * weightRate(ntf) *
     &         tzRES(ip,ie,ia,it,izf(ir,iz),im,1) )

             ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
             !                    + F**2 r2**2 sig_N2(xiwi)**2
             tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,2) =
     &       tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,2) +
     &               ( sumfactor**2 * weightRate(ntf)**2 *
     &          tzRES(ip,ie,ia,it,izf(ir,iz),im,2) )

            else
             tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) =
     &               ( sumfactor * weightRate(ntf) *
     &          tzRES(ip,ie,ia,it,izf(ir,iz),im,1) )

             tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,2) =
     &               ( sumfactor**2 * weightRate(ntf)**2 *
     &          tzRES(ip,ie,ia,it,izf(ir,iz),im,2) )

            end if

            if( tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) > cmax )
     &           cmax = tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1)
            if( tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) < cmin )
     &           cmin = tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1)

           end do    ! ip loop end
          end do     ! ie loop end
         end do      ! ia loop end
        end do       ! it loop end
       end do        ! ir loop end
      end do         ! iz loop end
      end do         ! im loop end

! sumover
      calcfact1 =  sumfactor * weightRate(ntf)
      calcfact2 = sumfactor**2 * weightRate(ntf)**2
      call sumtal_calc_stdev_crsrz_trtz_sub(1,m, ntf,
     &                                 calcfact1,calcfact2)

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
!     resc2SUMTAL = resc2SUMTAL

      ! N = N1 + N2
!     resc3SUMTAL = resc3SUMTAL + resc3(m)

!     resc2(m) = resc2SUMTAL
!     resc3(m) = resc3SUMTAL

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

            if( ntf > 1 ) then

             ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
             !                 + F r2 sig_N2 xi wi
             trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,1) =
     &       trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,1) +
     &               ( sumfactor * weightRate(ntf) *
     &         trRES(ip,ie,ia,it,irf(ir,iz),im,1) )

             ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
             !                    + F**2 r2**2 sig_N2(xiwi)**2
             trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,2) =
     &       trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,2) +
     &               ( sumfactor**2 * weightRate(ntf)**2 *
     &          trRES(ip,ie,ia,it,irf(ir,iz),im,2) )

            else
             trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,1) =
     &               ( sumfactor * weightRate(ntf) *
     &          trRES(ip,ie,ia,it,irf(ir,iz),im,1) )

             trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,2) =
     &               ( sumfactor**2 * weightRate(ntf)**2 *
     &          trRES(ip,ie,ia,it,irf(ir,iz),im,2) )

            end if


           end do    ! ip loop end
          end do     ! ie loop end
         end do      ! ia loop end
        end do       ! it loop end
       end do        ! ir loop end
      end do         ! iz loop end
      end do         ! im loop end

! sumover
      calcfact1 =  sumfactor * weightRate(ntf)
      calcfact2 = sumfactor**2 * weightRate(ntf)**2
      call sumtal_calc_stdev_crsrz_trtz_sub(0,m, ntf,
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
        ErrID = 'L:671/R:sumtal_calc_stdev_crsrz/F:sumtalcross.f' !E52_002_002


        write(*,'(/'' ***** Error in sumtal calc stdev crsrz *****''/)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_crsrz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_crsxyz(m,ntf,
     &                             nt,np,na,ne,ny,nx,nz,nm,
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
      double precision  trSUMTAL(np,ne,na,nt,nx*ny*(nz+1),nm,2)

      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

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
      do iz = 1, nz + 1
       do ix = 1, nx
        do iy = 1, ny
         do ie = 1, ne
          do ia = 1, na
           do ip = 1, np
            do it = 1, nt

             if( ntf > 1 ) then

              ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
              !                 + F r2 sig_N2 xi wi
              trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) =
     &        trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) +
     &          ( sumfactor * weightRate(ntf) *
     &           trRES(ip,ie,ia,it,icf(ix,iy,iz),im,1) )

              ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
              !                    + F**2 r2**2 sig_N2(xiwi)**2
              trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,2) =
     &        trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,2) +
     &          ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie,ia,it,icf(ix,iy,iz),im,2) )

             else
              trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) =
     &          ( sumfactor * weightRate(ntf) *
     &           trRES(ip,ie,ia,it,icf(ix,iy,iz),im,1) )

              trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,2) =
     &          ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie,ia,it,icf(ix,iy,iz),im,2) )

             end if

             if( trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) > cmax )
     &             cmax = trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1)
             if( trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) < cmin )
     &             cmin = trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1)

            end do   ! it loop end
           end do    ! ip loop end
          end do     ! ia loop end
         end do      ! ie loop end
        end do       ! iy loop end
       end do        ! ix loop end
      end do         ! iz loop end
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
        ErrID = 'L:840/R:sumtal_calc_stdev_crsxyz/F:sumtalcross.f' !E52_002_003


        write(*,'(/'' ***** Error in sumtal calc stdev crsxyz *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_crsxyz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_crsreg(m,ntf,nfile,
     &                               np,ne,nt,na,nr,nm,
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

      integer           np
      integer           ne
      integer           nt
      integer           na
      integer           nr
      integer           nm

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (np,ne,na,nt,nr,nm,2)
      double precision  trSUMTAL(np,ne,na,nt,nr,nm,2)

      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir
      integer           ia
      integer           it
      integer           ie
      integer           ip
      integer           im

      double precision  Xa, sigx
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
     &                      trRES(ip,ie,ia,it,ir,im,1),
     &                      trRES(ip,ie,ia,it,ir,im,2),
     &                      1.0d+0)
            trRES(ip,ie,ia,it,ir,im,1) = Xa
            trRES(ip,ie,ia,it,ir,im,2) = sigx
           end if

           if( ntf > 1 ) then

            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(ip,ie,ia,it,ir,im,1) =
     &      trSUMTAL(ip,ie,ia,it,ir,im,1) +
     &               weightRate(ntf)
     &       * trRES(ip,ie,ia,it,ir,im,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(ip,ie,ia,it,ir,im,2) =
     &      trSUMTAL(ip,ie,ia,it,ir,im,2) +
     &      (  trRES(ip,ie,ia,it,ir,im,2)
     &       * trRES(ip,ie,ia,it,ir,im,1) )**2
     &       * weightRate(ntf)**2

           ! ntf = 1
           else
            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(ip,ie,ia,it,ir,im,1) =
     &               weightRate(ntf)
     &       * trRES(ip,ie,ia,it,ir,im,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(ip,ie,ia,it,ir,im,2) =
     &      (  trRES(ip,ie,ia,it,ir,im,2)
     &       * trRES(ip,ie,ia,it,ir,im,1) )**2
     &       * weightRate(ntf)**2
           end if

          end do     ! ip loop end
         end do      ! ie loop end
        end do       ! it loop end
       end do        ! ia loop end
      end do         ! ir loop end
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
        do ia = 1, na
         do it = 1, nt
          do ie = 1, ne
           do ip = 1, np

            ! X_bar
            trSUMTAL(ip,ie,ia,it,ir,im,1) =
     &      trSUMTAL(ip,ie,ia,it,ir,im,1) * sumfactor/sumWR

            ! sig_x
            trSUMTAL(ip,ie,ia,it,ir,im,2) = sqrt(
     &      trSUMTAL(ip,ie,ia,it,ir,im,2)) * sumfactor/sumWR


            ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
            trSUMTAL(ip,ie,ia,it,ir,im,2) =
     &     (trSUMTAL(ip,ie,ia,it,ir,im,2)**2 *
     &           resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &           resc3SUMTAL *
     &      trSUMTAL(ip,ie,ia,it,ir,im,1)**2) *
     &          (resc2SUMTAL/resc3SUMTAL)**2

            ! Sigma xi wi = X_bar W
            trSUMTAL(ip,ie,ia,it,ir,im,1) =
     &      trSUMTAL(ip,ie,ia,it,ir,im,1) * resc2SUMTAL

            if( trSUMTAL(ip,ie,ia,it,ir,im,1) > cmax )
     &           cmax = trSUMTAL(ip,ie,ia,it,ir,im,1)
            if( trSUMTAL(ip,ie,ia,it,ir,im,1) < cmin )
     &           cmin = trSUMTAL(ip,ie,ia,it,ir,im,1)

           end do     ! ip loop end
          end do      ! ie loop end
         end do       ! it loop end
        end do        ! ia loop end
       end do         ! ir loop end
       end do         ! im loop end

!sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_crsreg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_crsrz(m,ntf,nfile,
     &                               np,ne,na,nt,nr,nz,nm,
     &                               sumfactor,weightRate,sumWR,
     &                               trRES,trSUMTAL,
     &                               tzRES,tzSUMTAL,
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

      integer           np
      integer           ne
      integer           na
      integer           nt
      integer           nr
      integer           nz
      integer           nm

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (np,ne,na,nt,(nr+1)*nz,nm,2)
      double precision  trSUMTAL(np,ne,na,nt,(nr+1)*nz,nm,2)
      double precision  tzRES   (np,ne,na,nt,nr*(nz+1),nm,2)
      double precision  tzSUMTAL(np,ne,na,nt,nr*(nz+1),nm,2)


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

         irf(ir,iz) = ir + ( iz - 1 ) * ( nr + 1 )
         izf(ir,iz) = iz + ( ir - 1 ) * ( nz + 1 )

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
     &                       tzRES(ip,ie,ia,it,izf(ir,iz),im,1),
     &                       tzRES(ip,ie,ia,it,izf(ir,iz),im,2),
     &                       1.0d+0)
             tzRES(ip,ie,ia,it,izf(ir,iz),im,1) = Xa
             tzRES(ip,ie,ia,it,izf(ir,iz),im,2) = sigx
            end if

            if( ntf > 1 ) then

             ! X_bar = F Sigma rj/r Xj_bar
             tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) =
     &       tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) +
     &              weightRate(ntf)
     &        * tzRES(ip,ie,ia,it,izf(ir,iz),im,1)

             ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
             tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,2) =
     &       tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,2) +
     &       (  tzRES(ip,ie,ia,it,izf(ir,iz),im,2)
     &        * tzRES(ip,ie,ia,it,izf(ir,iz),im,1) )**2
     &            * weightRate(ntf)**2

            ! ntf = 1
            else
             ! X_bar = F Sigma rj/r Xj_bar
             tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) =
     &              weightRate(ntf)
     &        * tzRES(ip,ie,ia,it,izf(ir,iz),im,1)

             ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
             tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,2) =
     &       (  tzRES(ip,ie,ia,it,izf(ir,iz),im,2)
     &        * tzRES(ip,ie,ia,it,izf(ir,iz),im,1) )**2
     &        * weightRate(ntf)**2
            end if

           end do    ! ip loop end
          end do     ! ie loop end
         end do      ! ia loop end
        end do       ! it loop end
       end do        ! ir loop end
      end do         ! iz loop end
      end do         ! im loop end

!sumover  tr_sum -> trSUMTAL_sum
      call calc_stdev_trtz_sum(1, m, 1.0d0)
      call sumtal_calc_average_trtz_sub(1, m, ntf,weightRate(ntf))

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
       do iz = 1, nz + 1
        do ir = 1, nr
         do it = 1, nt
          do ia = 1, na
           do ie = 1, ne
            do ip = 1, np

             ! X_bar
             tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) =
     &       tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) * sumfactor/sumWR

             ! sig_x
             tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,2) = sqrt(
     &       tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,2)) * sumfactor/sumWR


             ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
             tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,2) =
     &      (tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,2)**2 *
     &            resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &            resc3SUMTAL *
     &       tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1)**2) *
     &           (resc2SUMTAL/resc3SUMTAL)**2

             ! Sigma xi wi = X_bar W
             tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) =
     &       tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) * resc2SUMTAL

             if( tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) > cmax )
     &            cmax = tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1)
             if( tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1) < cmin )
     &            cmin = tzSUMTAL(ip,ie,ia,it,izf(ir,iz),im,1)

            end do    ! ip loop end
           end do     ! ie loop end
          end do      ! ia loop end
         end do       ! it loop end
        end do        ! ir loop end
       end do         ! iz loop end
       end do         ! im loop end

!sumover
      call sumtal_calc_average_trtz_sub2(1,m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

      end if
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
     &                       trRES(ip,ie,ia,it,irf(ir,iz),im,1),
     &                       trRES(ip,ie,ia,it,irf(ir,iz),im,2),
     &                       1.0d+0)
             trRES(ip,ie,ia,it,irf(ir,iz),im,1) = Xa
             trRES(ip,ie,ia,it,irf(ir,iz),im,2) = sigx
            end if

            if( ntf > 1 ) then

             ! X_bar = F Sigma rj/r Xj_bar
             trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,1) =
     &       trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,1) +
     &              weightRate(ntf)
     &        * trRES(ip,ie,ia,it,irf(ir,iz),im,1)

             ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
             trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,2) =
     &       trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,2) +
     &       (  trRES(ip,ie,ia,it,irf(ir,iz),im,2)
     &        * trRES(ip,ie,ia,it,irf(ir,iz),im,1) )**2
     &            * weightRate(ntf)**2

            ! ntf = 1
            else
             ! X_bar = F Sigma rj/r Xj_bar
             trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,1) =
     &              weightRate(ntf)
     &        * trRES(ip,ie,ia,it,irf(ir,iz),im,1)

             ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
             trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,2) =
     &       (  trRES(ip,ie,ia,it,irf(ir,iz),im,2)
     &        * trRES(ip,ie,ia,it,irf(ir,iz),im,1) )**2
     &        * weightRate(ntf)**2
            end if

           end do    ! ip loop end
          end do     ! ie loop end
         end do      ! ia loop end
        end do       ! it loop end
       end do        ! ir loop end
      end do         ! iz loop end
      end do         ! im loop end

!sumover  tr_sum -> trSUMTAL_sum
      call calc_stdev_trtz_sum(0, m, 1.0d0)
      call sumtal_calc_average_trtz_sub(0,m, ntf,weightRate(ntf))

!-----------------------------------------------------------------------
  200 continue
      ! W = r1 N1 w1_bar + r2 N2 w2_bar
!     resc2SUMTAL = resc2SUMTAL

      ! N = N1 + N2
!     resc3SUMTAL = resc3SUMTAL + resc3(m)
!-----------------------------------------------------------------------
      if( ntf == nfile ) then
       do im = 1, nm
       do iz = 1, nz
        do ir = 1, nr + 1
         do it = 1, nt
          do ia = 1, na
           do ie = 1, ne
            do ip = 1, np

             ! X_bar
             trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,1) =
     &       trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,1) * sumfactor/sumWR

             ! sig_x
             trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,2) = sqrt(
     &       trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,2)) * sumfactor/sumWR


             ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
             trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,2) =
     &      (trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,2)**2 *
     &            resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &            resc3SUMTAL *
     &       trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,1)**2) *
     &           (resc2SUMTAL/resc3SUMTAL)**2

             ! Sigma xi wi = X_bar W
             trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,1) =
     &       trSUMTAL(ip,ie,ia,it,irf(ir,iz),im,1) * resc2SUMTAL


            end do    ! ip loop end
           end do     ! ie loop end
          end do      ! ia loop end
         end do       ! it loop end
        end do        ! ir loop end
       end do         ! iz loop end
       end do         ! im loop end

!sumover
      call sumtal_calc_average_trtz_sub2(0, m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_crsrz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_crsxyz(m,ntf,nfile,
     &                               nt,np,na,ne,ny,nx,nz,nm,
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
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (np,ne,na,nt,nx*ny*(nz+1),nm,2)
      double precision  trSUMTAL(np,ne,na,nt,nx*ny*(nz+1),nm,2)


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

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny
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

             if( ntf > 1 ) then

              ! X_bar = F Sigma rj/r Xj_bar
              trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) =
     &        trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) +
     &           weightRate(ntf)
     &         * trRES(ip,ie,ia,it,icf(ix,iy,iz),im,1)

              ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
              trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,2) =
     &        trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,2) +
     &        (  trRES(ip,ie,ia,it,icf(ix,iy,iz),im,2)
     &         * trRES(ip,ie,ia,it,icf(ix,iy,iz),im,1) )**2
     &         * weightRate(ntf)**2

             ! ntf = 1
             else
              ! X_bar = F Sigma rj/r Xj_bar
              trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) =
     &           weightRate(ntf)
     &         * trRES(ip,ie,ia,it,icf(ix,iy,iz),im,1)

              ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
              trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,2) =
     &        (  trRES(ip,ie,ia,it,icf(ix,iy,iz),im,2)
     &         * trRES(ip,ie,ia,it,icf(ix,iy,iz),im,1) )**2
     &         * weightRate(ntf)**2
             end if

            end do   ! it loop end
           end do    ! ip loop end
          end do     ! ia loop end
         end do      ! ie loop end
        end do       ! iy loop end
       end do        ! ix loop end
      end do         ! iz loop end
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
       do iz = 1, nz + 1
        do ix = 1, nx
         do iy = 1, ny
          do ie = 1, ne
           do ia = 1, na
            do ip = 1, np
             do it = 1, nt

              ! X_bar
              trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) =
     &        trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) * sumfactor/sumWR

              ! sig_x
              trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,2) = sqrt(
     &        trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,2))* sumfactor/sumWR


              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,2) =
     &       (trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,2)**2 *
     &        resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &        resc3SUMTAL *
     &        trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1)**2) *
     &        (resc2SUMTAL/resc3SUMTAL)**2

              ! Sigma xi wi = X_bar W
              trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) =
     &        trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) * resc2SUMTAL

              if( trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) > cmax )
     &              cmax = trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1)
              if( trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1) < cmin )
     &              cmin = trSUMTAL(ip,ie,ia,it,icf(ix,iy,iz),im,1)

             end do   ! it loop end
            end do    ! ip loop end
           end do     ! ia loop end
          end do      ! ie loop end
         end do       ! iy loop end
        end do        ! ix loop end
       end do         ! iz loop end
       end do         ! im loop end

!sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_crsxyz

!***********************************************************************
!                                                                      *
!  sumover subroutine                                                  *
!                                                                      *
!***********************************************************************

!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_crsrz_trtz_sub(kind,m, ntf,
     &                                      calcfact1,calcfact2)
!                                                                      *
!***********************************************************************

      use RESTALMOD, only: irestalm_sum,lrestalm_sum
      use TALMOD
!$      use TALMOD0
      use sumtallymod

      implicit none

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      integer :: kind, m, iax, ntf

      double precision :: calcfact(2),calcfact1,calcfact2

      real(8),pointer :: p_sum(:)
!      real(8),pointer :: resta_sum(:)
      real(8),pointer :: sumtal_sum(:)

      integer :: ln_sum, iskip

       calcfact(1) = calcfact1
       calcfact(2) = calcfact2

      do iax =1,6

       if(kind == 0) then
!       resta_sum => trRES_sum(irestalm_sum(m,iax):)

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(p_sum,m,iax)
!$          ln_sum = italsize0_2_sum(m,iax)/2
!$       else
            call GET_TR_HEAD_POINTER_SUM(p_sum,m,iax)
            ln_sum = italsize_2_sum(m,iax)/2
C for nonshared_tally option
!$       end if

         sumtal_sum => trSUMTAL_SUM(irestalm_sum(m,iax):)

       else
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

         sumtal_sum => trSUMTAL_SUM(irestalm_sum(m,iax)+iskip:)

       endif

        if(ln_sum > 0 ) then

          call calc_stdev_restal2sumtally(p_sum,sumtal_sum,
     &       ln_sum, ntf, calcfact)
        endif
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_trtz_sub(kind,m, ntf,
     &                                      calcfact)
!                                                                      *
!***********************************************************************

      use RESTALMOD, only: irestalm_sum,lrestalm_sum
      use TALMOD
!$      use TALMOD0
      use sumtallymod

      implicit none

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      integer :: kind, m, iax, ntf

      double precision :: calcfact

      real(8),pointer :: p_sum(:)
!      real(8),pointer :: resta_sum(:)
      real(8),pointer :: sumtal_sum(:)

      integer :: ln_sum, iskip

      do iax =1,6
       if(kind == 0) then
!       resta_sum => trRES_sum(irestalm_sum(m,iax):)

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(p_sum,m,iax)
!$          ln_sum = italsize0_2_sum(m,iax)/2
!$       else
            call GET_TR_HEAD_POINTER_SUM(p_sum,m,iax)
            ln_sum = italsize_2_sum(m,iax)/2
C for nonshared_tally option
!$       end if

         sumtal_sum => trSUMTAL_SUM(irestalm_sum(m,iax):)

       else
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

         sumtal_sum => trSUMTAL_SUM(irestalm_sum(m,iax)+iskip:)

       endif

        if(ln_sum > 0 ) then

          call calc_average_restal2sumtally(p_sum,sumtal_sum,
     &       ln_sum, ntf, calcfact)
        endif
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_trtz_sub2(kind, m, ntf,
     &           sumfactor1,sumWR1,resc2SUMTAL1,resc3SUMTAL1)
!                                                                      *
!***********************************************************************

      use TALMOD, only: italsize_sum,italsize_2_sum
      use RESTALMOD, only: irestalm_sum,lrestalm_sum
      use sumtallymod

      implicit none

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      integer :: kind, m, iax, ntf

      real(8) :: sumfactor1,sumWR1,resc2SUMTAL1,resc3SUMTAL1

      real(8),pointer :: sumtal_sum(:)

      integer :: ln_sum, iskip

      do iax =1,6

        if(kind == 0) then
          sumtal_sum => trSUMTAL_SUM(irestalm_sum(m,iax):)
          ln_sum = italsize_2_sum(m,iax)/2
        else

          iskip = italsize_2_sum(m,iax)
          sumtal_sum => trSUMTAL_SUM(irestalm_sum(m,iax)+iskip:)
          ln_sum = (italsize_sum(m,iax)-italsize_2_sum(m,iax))/2
        endif

        if(ln_sum > 0 ) then

          call calc_average_sumtally(sumtal_sum,
     &       ln_sum, ntf,
     &       sumfactor1,sumWR1,resc2SUMTAL1,resc3SUMTAL1)
        endif
      enddo

      return
      end

