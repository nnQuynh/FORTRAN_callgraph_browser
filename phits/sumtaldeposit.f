************************************************************************
*                                                                      *
      subroutine sumtal_read_tdeposit(m,iax,ntf,ierr)
*                                                                      *
*     original routine is read_tdeposit in restdeposit.f               *
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

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

        common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                  rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

c$$$        dimension     idas(1)
        dimension idas(mdas*2) !2025.2.26 S.H. temporary for fbounds-check
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

        if ( any( itaxs(m,iax) .eq. (/ 7, 8, 9, 10, 12, 13 /) )
     &      .and. ittwo(m) .ne. 4 ) noe = 2

*-----------------------------------------------------------------------
*   open restart file
*-----------------------------------------------------------------------
        call sumtal_open_resfile
     &       (m,noe,ntf,jsn,jsi,dsin,idsi,ill,ilf,newtall,ierr)

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
         call check_tdeposit(m,iax,jsn(ioe),jsi(ioe),
     &                   dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                   ierr)
        endif

        if (ierr .ne. 0 ) goto 900

*-----------------------------------------------------------------------

        rewind(jsi(ioe))
        ill(jsn(ioe),ioe) = 1

*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          idas0 = nmmax
          idas1 = lmmax
          idas2 = idas1 + ittnm(m)

          idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
          idas4 = idas2 + itrgn(m)
     &          + ( itrgn(m) + mod(itrgn(m),2) ) / 2
          idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
          idasa = idas3

        else if ( itmsh(m) .eq. 2 ) then

          idas0 = nmmax
          idas1 = lmmax
          idasa = idas1 + ittnm(m)

        else if ( itmsh(m) .eq. 3 ) then

          idas0 = nmmax
          idas1 = lmmax
          idas2 = ( idas1 + ittnm(m) - 1 ) * 2 + 1
          idasa = idas1 + ittnm(m)

        else if( itmsh(m) .eq. 4 ) then
cFURUTA20190121 ! copy from reg
          idas0 = nmmax
          idas1 = lmmax
          idas2 = idas1 + ittnm(m)

          idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
          idas4 = idas2 + itrgn(m)
     &          + ( itrgn(m) + mod(itrgn(m),2) ) / 2
          idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
          idasa = idas3

        end if

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_depstreg(m,iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itpan(m),itrgn(m),itrgm(m),itenm(m),
     &      ittnm(m),
     &      idas_itreg(itreg(m)),
     &      das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call read_depstrz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itpan(m),itrnm(m),itznm(m),itenm(m),
     &      ittnm(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call read_depstxyz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itpan(m),
     &      itxnm(m),itynm(m),itznm(m),itenm(m),
     &      ittnm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &      trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then

          call read_depsttet(m,iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itpan(m),itrgn(m),itenm(m),ittnm(m),
     &      das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &      trRES(irestalm(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_depstreg(m,
     &      itpan(m),itrgn(m),itrgm(m),itenm(m),
     &      ittnm(m),
     &      idas_itreg(itreg(m)),
     &      das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call restore_depstrz(m,
     &      itpan(m),itrnm(m),itznm(m),itenm(m),
     &      ittnm(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call restore_depstxyz(m,
     &      itpan(m),
     &      itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),itenm(m),
     &      ittnm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &      trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then

          call restore_depsttet(m,
     &      itpan(m),itrgn(m),itrgm(m),itenm(m),ittnm(m),
     &      das_iterg(iterg(m)),das_ittrg(ittrg(m)),
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

      end subroutine sumtal_read_tdeposit

!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_depstreg(m,ntf,
     &                             np,nei,ne,nt,nr,
     &                             sumfactor,weightRate,
     &                             trRES,trSUMTAL,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by S.Hashimoto on 2015/1/19                              *
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
      integer           ne, nei
      integer           nt
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,0:ne,nr,nt,2)
      double precision  trSUMTAL(np,0:ne,nr,nt,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir

      integer           ip
      integer           ie
      integer           it

      integer,save   :: maxcasd

      integer   nr0
      integer   itrwgtsum
      integer   itrgn1
      integer   itrgm1
      integer   itrncd
      common /tall70/ itrwgtsum(itlmax), itrgn1(itlmax),
     &                itrgm1(itlmax), itrncd(itlmax)

!-----------------------------------------------------------------------

      ierr = 0

      if( itrwgtsum(m) .eq. 1 )then
       nr0 = 1
      else
       nr0 = nr
      endif

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr0
       do it = 1, nt
        do ie = nei, ne
         do ip = 1, np

          if( ntf > 1 ) then

           ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
           !                 + F r2 sig_N2 xi wi
           trSUMTAL(ip,ie,ir,it,1) =
     &            trSUMTAL(ip,ie,ir,it,1) +
     &            ( sumfactor * weightRate(ntf) *
     &            trRES(ip,ie,ir,it,1) )

           ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
           !                    + F**2 r2**2 sig_N2(xiwi)**2
           trSUMTAL(ip,ie,ir,it,2) =
     &          trSUMTAL(ip,ie,ir,it,2) +
     &          ( sumfactor**2 * weightRate(ntf)**2 *
     &          trRES(ip,ie,ir,it,2) )

          else
           trSUMTAL(ip,ie,ir,it,1) =
     &            ( sumfactor * weightRate(ntf) *
     &            trRES(ip,ie,ir,it,1) )

           trSUMTAL(ip,ie,ir,it,2) =
     &          ( sumfactor**2 * weightRate(ntf)**2 *
     &          trRES(ip,ie,ir,it,2) )

          end if

          if( trSUMTAL(ip,ie,ir,it,1) > cmax )
     &         cmax = trSUMTAL(ip,ie,ir,it,1)
          if( trSUMTAL(ip,ie,ir,it,1) < cmin )
     &         cmin = trSUMTAL(ip,ie,ir,it,1)

         end do   ! ip loop end
        end do    ! ie loop end
       end do     ! it loop end
      end do      ! ir loop end

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
        ErrID = 'L:413/R:sumtal_calc_stdev_depstreg/F:sumtaldeposit.f' !E54_002_001


        write(*,'(/'' ***** Error in sumtal calc stdev depstreg *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_depstreg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_depstrz(m,ntf,
     &                             np,nei,ne,nt,nr,nz,
     &                             sumfactor,weightRate,
     &                             trRES,trSUMTAL,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by S.Hashimoto on 2015/1/16                              *
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
      integer           ne, nei
      integer           nt
      integer           nr
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,0:ne,nt,nr,nz,2)
      double precision  trSUMTAL(np,0:ne,nt,nr,nz,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir
      integer           iz

      integer           ip
      integer           ie
      integer           it

      integer,save   :: maxcasd
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do iz = 1, nz
       do ir = 1, nr
        do it = 1, nt
         do ie = nei, ne
          do ip = 1, np

           if( ntf > 1 ) then

            ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
            !                 + F r2 sig_N2 xi wi
            trSUMTAL(ip,ie,it,ir,iz,1) =
     &             trSUMTAL(ip,ie,it,ir,iz,1) +
     &             ( sumfactor * weightRate(ntf) *
     &             trRES(ip,ie,it,ir,iz,1) )

            ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
            !                    + F**2 r2**2 sig_N2(xiwi)**2
            trSUMTAL(ip,ie,it,ir,iz,2) =
     &           trSUMTAL(ip,ie,it,ir,iz,2) +
     &           ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie,it,ir,iz,2) )

           else
            trSUMTAL(ip,ie,it,ir,iz,1) =
     &             ( sumfactor * weightRate(ntf) *
     &             trRES(ip,ie,it,ir,iz,1) )

            trSUMTAL(ip,ie,it,ir,iz,2) =
     &           ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie,it,ir,iz,2) )

           end if

           if( trSUMTAL(ip,ie,it,ir,iz,1) > cmax )
     &          cmax = trSUMTAL(ip,ie,it,ir,iz,1)
           if( trSUMTAL(ip,ie,it,ir,iz,1) < cmin )
     &          cmin = trSUMTAL(ip,ie,it,ir,iz,1)

          end do   ! ip loop end
         end do    ! ie loop end
        end do     ! it loop end
       end do      ! ir loop end
      end do       ! iz loop end

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
        ErrID = 'L:565/R:sumtal_calc_stdev_depstrz/F:sumtaldeposit.f' !E54_002_002


        write(*,'(/'' ***** Error in sumtal calc stdev depstrz *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_depstrz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev_depstxyz(m,ntf,
     &                             np,nei,ne,nt,nx,ny,nz,
     &                             sumfactor,weightRate,
     &                             trRES,trSUMTAL,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by S.Hashimoto on 2015/1/16                              *
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
      integer           ne, nei
      integer           nt
      integer           nx
      integer           ny
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,0:ne,nt,nx*ny*nz,2)
      double precision  trSUMTAL(np,0:ne,nt,nx*ny*nz,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1,calcfact2

      integer           ierr
!-----------------------------------------------------------------------
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

      do iz = 1, nz
       do iy = 1, ny
        do ix = 1, nx
         do it = 1, nt
          do ie = nei, ne
           do ip = 1, np

            if( ntf > 1 ) then

              ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
              !                 + F r2 sig_N2 xi wi
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) =
     &        trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) +
     &          ( sumfactor * weightRate(ntf) *
     &           trRES(ip,ie,it,icf(ix,iy,iz),1) )

              ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
              !                    + F**2 r2**2 sig_N2(xiwi)**2
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) =
     &        trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) +
     &          ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie,it,icf(ix,iy,iz),2) )

             else
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) =
     &          ( sumfactor * weightRate(ntf) *
     &           trRES(ip,ie,it,icf(ix,iy,iz),1) )

              trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) =
     &          ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie,it,icf(ix,iy,iz),2) )

             end if

             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) > cmax )
     &             cmax = trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)
             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) < cmin )
     &             cmin = trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)

           end do   ! ip loop end
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
        ErrID = 'L:725/R:sumtal_calc_stdev_depstxyz/F:sumtaldeposit.f' !E54_002_003


        write(*,'(/'' ***** Error in sumtal calc stdev depstxyz *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev_depstxyz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_depstreg(m,ntf,nfile,
     &                               np,nei,ne,nt,nr,
     &                               sumfactor,weightRate,sumWR,
     &                               trRES,trSUMTAL,
     &                               resc2SUMTAL,resc3SUMTAL,
     &                               ierr)
!                                                                      *
!     created by S.Hashimoto on 2015/1/19                              *
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
      integer           ne, nei
      integer           nt
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (np,0:ne,nr,nt,2)
      double precision  trSUMTAL(np,0:ne,nr,nt,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir

      integer           ip
      integer           ie
      integer           it
      double precision  Xa, sigx

      integer   nr0
      integer   itrwgtsum
      integer   itrgn1
      integer   itrgm1
      integer   itrncd
      common /tall70/ itrwgtsum(itlmax), itrgn1(itlmax),
     &                itrgm1(itlmax), itrncd(itlmax)

!-----------------------------------------------------------------------

      ierr = 0

      if( itrwgtsum(m) .eq. 1 ) then
       nr0 = 1
      else
       nr0 = nr
      endif

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr0
       do it = 1, nt
        do ie = nei, ne
         do ip = 1, np

          if( trRES(ip,ie,ir,it,1) .gt. 0.d0 ) then

           call calc_deposit_stdev(m,Xa,sigx,
     &            trRES(ip,ie,ir,it,1),
     &            trRES(ip,ie,ir,it,2),
     &            1.0d+0,ip)
           trRES(ip,ie,ir,it,1) = Xa
           trRES(ip,ie,ir,it,2) = sigx
          end if

          if( ntf > 1 ) then
           ! X_bar = F Sigma rj/r Xj_bar
           trSUMTAL(ip,ie,ir,it,1) =
     &            trSUMTAL(ip,ie,ir,it,1) +
     &            weightRate(ntf)
     &            * trRES(ip,ie,ir,it,1)

           ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
           trSUMTAL(ip,ie,ir,it,2) =
     &          trSUMTAL(ip,ie,ir,it,2) +
     &          (  trRES(ip,ie,ir,it,2)
     &          * trRES(ip,ie,ir,it,1) )**2
     &          * weightRate(ntf)**2

           ! ntf = 1

          else
           ! X_bar = F Sigma rj/r Xj_bar
           trSUMTAL(ip,ie,ir,it,1) =
     &            weightRate(ntf)
     &            * trRES(ip,ie,ir,it,1)

           ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
           trSUMTAL(ip,ie,ir,it,2) =
     &          (  trRES(ip,ie,ir,it,2)
     &          * trRES(ip,ie,ir,it,1) )**2
     &          * weightRate(ntf)**2

          end if

         end do   ! ip loop end
        end do    ! ie loop end
       end do     ! it loop end
      end do      ! ir loop end

!sumover  tr_sum -> trSUMTAL_sum
      call calc_deposit_stdev_tr_sum(m, 1.0d0, np)
      call sumtal_calc_average_sub(m, ntf,weightRate(ntf))

!-----------------------------------------------------------------------

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = 1d0

      ! N = N1 + N2
      resc3SUMTAL = 2d0
!-----------------------------------------------------------------------

      if( ntf == nfile ) then
       do ir = 1, nr0
        do it = 1, nt
         do ie = nei, ne
          do ip = 1, np

           ! X_bar
           trSUMTAL(ip,ie,ir,it,1) =
     &            trSUMTAL(ip,ie,ir,it,1) * sumfactor/sumWR

           ! sig_x
           trSUMTAL(ip,ie,ir,it,2) = sqrt(
     &          trSUMTAL(ip,ie,ir,it,2)) * sumfactor/sumWR


           ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
           trSUMTAL(ip,ie,ir,it,2) =
     &          (trSUMTAL(ip,ie,ir,it,2)**2 *
     &          resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &          resc3SUMTAL *
     &          trSUMTAL(ip,ie,ir,it,1)**2) *
     &          (resc2SUMTAL/resc3SUMTAL)**2

           ! Sigma xi wi = X_bar W
           trSUMTAL(ip,ie,ir,it,1) =
     &          trSUMTAL(ip,ie,ir,it,1) * resc2SUMTAL

C S.H. added the following processes for negative values when isumtally=3.
C This is tentative. Correct calculation can be done if processes of
C subroutine calc_stdev is performed here, but preparation of the last
C parameter of calc_stdev is needed. (2020.12.16)
           if( trSUMTAL(ip,ie,ir,it,1) < 0d0 ) then
              trSUMTAL(ip,ie,ir,it,1) = 0d0
              trSUMTAL(ip,ie,ir,it,2) = 0d0
           end if

           if( trSUMTAL(ip,ie,ir,it,1) > cmax )
     &          cmax = trSUMTAL(ip,ie,ir,it,1)
           if( trSUMTAL(ip,ie,ir,it,1) < cmin )
     &          cmin = trSUMTAL(ip,ie,ir,it,1)

          end do   ! ip loop end
         end do    ! ie loop end
        end do     ! it loop end
       end do      ! ir loop end

!sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_depstreg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_depstrz(m,ntf,nfile,
     &                               np,nei,ne,nt,nr,nz,
     &                               sumfactor,weightRate,sumWR,
     &                               trRES,trSUMTAL,
     &                               resc2SUMTAL,resc3SUMTAL,
     &                               ierr)
!                                                                      *
!     created by S.Hashimoto on 2015/1/16                              *
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
      integer           ne, nei
      integer           nt
      integer           nr
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (np,0:ne,nt,nr,nz,2)
      double precision  trSUMTAL(np,0:ne,nt,nr,nz,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir
      integer           iz

      integer           ip
      integer           ie
      integer           it
      double precision  Xa, sigx
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do iz = 1, nz
       do ir = 1, nr
        do it = 1, nt
         do ie = nei, ne
          do ip = 1, np

           if( trRES(ip,ie,it,ir,iz,1) .gt. 0.d0 ) then
            call calc_deposit_stdev(m,Xa,sigx,
     &             trRES(ip,ie,it,ir,iz,1),
     &             trRES(ip,ie,it,ir,iz,2),
     &             1.0d+0,ip)
            trRES(ip,ie,it,ir,iz,1) = Xa
            trRES(ip,ie,it,ir,iz,2) = sigx
           end if

           if( ntf > 1 ) then
            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(ip,ie,it,ir,iz,1) =
     &             trSUMTAL(ip,ie,it,ir,iz,1) +
     &             weightRate(ntf)
     &             * trRES(ip,ie,it,ir,iz,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(ip,ie,it,ir,iz,2) =
     &           trSUMTAL(ip,ie,it,ir,iz,2) +
     &           (  trRES(ip,ie,it,ir,iz,2)
     &           * trRES(ip,ie,it,ir,iz,1) )**2
     &           * weightRate(ntf)**2

            ! ntf = 1
           else
            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL(ip,ie,it,ir,iz,1) =
     &             weightRate(ntf)
     &             * trRES(ip,ie,it,ir,iz,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL(ip,ie,it,ir,iz,2) =
     &           (  trRES(ip,ie,it,ir,iz,2)
     &           * trRES(ip,ie,it,ir,iz,1) )**2
     &           * weightRate(ntf)**2

           end if

          end do   ! ip loop end
         end do    ! ie loop end
        end do     ! it loop end
       end do      ! ir loop end
      end do       ! iz loop end

!sumover  tr_sum -> trSUMTAL_sum
      call calc_deposit_stdev_tr_sum(m, 1.0d0, np)
      call sumtal_calc_average_sub(m, ntf,weightRate(ntf))

!-----------------------------------------------------------------------

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = 1d0

      ! N = N1 + N2
      resc3SUMTAL = 2d0
!-----------------------------------------------------------------------

      if( ntf == nfile ) then
       do iz = 1, nz
        do ir = 1, nr
         do it = 1, nt
          do ie = nei, ne
           do ip = 1, np

            ! X_bar
            trSUMTAL(ip,ie,it,ir,iz,1) =
     &             trSUMTAL(ip,ie,it,ir,iz,1) * sumfactor/sumWR

            ! sig_x
            trSUMTAL(ip,ie,it,ir,iz,2) = sqrt(
     &           trSUMTAL(ip,ie,it,ir,iz,2)) * sumfactor/sumWR


            ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
            trSUMTAL(ip,ie,it,ir,iz,2) =
     &           (trSUMTAL(ip,ie,it,ir,iz,2)**2 *
     &           resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &           resc3SUMTAL *
     &           trSUMTAL(ip,ie,it,ir,iz,1)**2) *
     &           (resc2SUMTAL/resc3SUMTAL)**2

            ! Sigma xi wi = X_bar W
            trSUMTAL(ip,ie,it,ir,iz,1) =
     &           trSUMTAL(ip,ie,it,ir,iz,1) * resc2SUMTAL

C S.H. added the following processes for negative values when isumtally=3.
C This is tentative. Correct calculation can be done if processes of
C subroutine calc_stdev is performed here, but preparation of the last
C parameter of calc_stdev is needed. (2020.12.16)
            if( trSUMTAL(ip,ie,it,ir,iz,1) < 0d0 ) then
               trSUMTAL(ip,ie,it,ir,iz,1) = 0d0
               trSUMTAL(ip,ie,it,ir,iz,2) = 0d0
            end if

            if( trSUMTAL(ip,ie,it,ir,iz,1) > cmax )
     &           cmax = trSUMTAL(ip,ie,it,ir,iz,1)
            if( trSUMTAL(ip,ie,it,ir,iz,1) < cmin )
     &           cmin = trSUMTAL(ip,ie,it,ir,iz,1)

           end do   ! ip loop end
          end do    ! ie loop end
         end do     ! it loop end
        end do      ! ir loop end
       end do       ! iz loop end

!sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_depstrz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_depstxyz(m,ntf,nfile,
     &                               np,nei,ne,nt,nx,ny,nz,
     &                               sumfactor,weightRate,sumWR,
     &                               trRES,trSUMTAL,
     &                               resc2SUMTAL,resc3SUMTAL,
     &                               ierr)
!                                                                      *
!     created by S.Hashimoto on 2015/1/16                              *
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
      integer           ne, nei
      integer           nt
      integer           nx
      integer           ny
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (np,0:ne,nt,nx*ny*nz,2)
      double precision  trSUMTAL(np,0:ne,nt,nx*ny*nz,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
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

      do iz = 1, nz
       do iy = 1, ny
        do ix = 1, nx
         do it = 1, nt
          do ie = nei, ne
           do ip = 1, np

            if( trRES(ip,ie,it,icf(ix,iy,iz),1) .gt. 0.d0 ) then
             call calc_deposit_stdev(m,Xa,sigx,
     &              trRES(ip,ie,it,icf(ix,iy,iz),1),
     &              trRES(ip,ie,it,icf(ix,iy,iz),2),
     &              1.0d+0,ip)
             trRES(ip,ie,it,icf(ix,iy,iz),1) = Xa
             trRES(ip,ie,it,icf(ix,iy,iz),2) = sigx

            end if

            if( ntf > 1 ) then
             ! X_bar = F Sigma rj/r Xj_bar
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) =
     &              trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) +
     &              weightRate(ntf)
     &              * trRES(ip,ie,it,icf(ix,iy,iz),1)

             ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) =
     &            trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) +
     &            (  trRES(ip,ie,it,icf(ix,iy,iz),2)
     &            * trRES(ip,ie,it,icf(ix,iy,iz),1) )**2
     &            * weightRate(ntf)**2

             ! ntf = 1
            else
             ! X_bar = F Sigma rj/r Xj_bar
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) =
     &              weightRate(ntf)
     &              * trRES(ip,ie,it,icf(ix,iy,iz),1)

             ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) =
     &            (  trRES(ip,ie,it,icf(ix,iy,iz),2)
     &            * trRES(ip,ie,it,icf(ix,iy,iz),1) )**2
     &            * weightRate(ntf)**2

            end if

           end do   ! ip loop end
          end do    ! ie loop end
         end do     ! it loop end
        end do      ! ix loop end
       end do       ! iy loop end
      end do        ! iz loop end

!sumover  tr_sum -> trSUMTAL_sum
      call calc_deposit_stdev_tr_sum(m, 1.0d0, np)
      call sumtal_calc_average_sub(m, ntf,weightRate(ntf))

!-----------------------------------------------------------------------

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = 1d0

      ! N = N1 + N2
      resc3SUMTAL = 2d0
!-----------------------------------------------------------------------

      if( ntf == nfile ) then
       do iz = 1, nz
        do iy = 1, ny
         do ix = 1, nx
          do it = 1, nt
           do ie = nei, ne
            do ip = 1, np

             ! X_bar
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) =
     &              trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) * sumfactor/sumWR

             ! sig_x
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) = sqrt(
     &            trSUMTAL(ip,ie,it,icf(ix,iy,iz),2)) * sumfactor/sumWR


             ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) =
     &            (trSUMTAL(ip,ie,it,icf(ix,iy,iz),2)**2 *
     &            resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &            resc3SUMTAL *
     &            trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)**2) *
     &            (resc2SUMTAL/resc3SUMTAL)**2

             ! Sigma xi wi = X_bar W
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) =
     &            trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) * resc2SUMTAL

C S.H. added the following processes for negative values when isumtally=3.
C This is tentative. Correct calculation can be done if processes of
C subroutine calc_stdev is performed here, but preparation of the last
C parameter of calc_stdev is needed. (2020.12.16)
             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) < 0d0 ) then
                trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) = 0d0
                trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) = 0d0
             end if

             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) > cmax )
     &            cmax = trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)
             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) < cmin )
     &            cmin = trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)

            end do   ! ip loop end
           end do    ! ie loop end
          end do     ! it loop end
         end do      ! ix loop end
        end do       ! iy loop end
       end do        ! iz loop end

!sumover
      call sumtal_calc_average_sub2(m, ntf,
     &           sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average_depstxyz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev2_depstreg(m,ntf,nfile,
     &                             np,nei,ne,nt,nr,
     &                             sumfactor,weightRate,
     &                             trRES,trSUMTAL,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by T.Miura on 2015/07/31                                 *
!                                                                      *
!***********************************************************************
      use partmod, only: itpan, itpat, jtpat, itmxpt ! frtati 2021/10/05
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
      integer           m             ! tally number
      integer           ntf           ! Count of tally file
      integer           nfile         ! Number of tally files

      integer           np
      integer           ne, nei
      integer           nt
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,0:ne,nr,nt,2)
      double precision  trSUMTAL(np,0:ne,nr,nt,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1, calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir

      integer           ip
      integer           ie
      integer           it

      integer,save   :: maxcasd

      integer   nr0
      integer   itrwgtsum
      integer   itrgn1
      integer   itrgm1
      integer   itrncd
      common /tall70/ itrwgtsum(itlmax), itrgn1(itlmax),
     &                itrgm1(itlmax), itrncd(itlmax)

!-----------------------------------------------------------------------
      ierr = 0

      if( itrwgtsum(m) .eq. 1 ) then
       nr0 = 1
      else
       nr0 = nr
      endif

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      ! N = Sigma Nj
      resc3SUMTAL = resc3SUMTAL + resc3(m)

      do ir = 1, nr0
       do it = 1, nt
        do ie = nei, ne
         do ip = 1, np

          if( ntf > 1 ) then
           if ( itpat(m,ip,1) .eq. 20 ) then
           ! Sigma Nj Xk_bar_j
            trSUMTAL(ip,ie,ir,it,1) =
     &             trSUMTAL(ip,ie,ir,it,1) +
     &             trRES(ip,ie,ir,it,1)

           else
           ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
           !                 + F r2 sig_N2 xi wi
            trSUMTAL(ip,ie,ir,it,1) =
     &             trSUMTAL(ip,ie,ir,it,1) +
     &             ( sumfactor * weightRate(ntf) *
     &             trRES(ip,ie,ir,it,1) )

           ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
           !                    + F**2 r2**2 sig_N2(xiwi)**2
            trSUMTAL(ip,ie,ir,it,2) =
     &           trSUMTAL(ip,ie,ir,it,2) +
     &           ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie,ir,it,2) )

           end if

          ! ntf = 1
          else
           if ( itpat(m,ip,1) .eq. 20 ) then
           ! Sigma Nj Xk_bar_j
            trSUMTAL(ip,ie,ir,it,1) =
     &             trRES(ip,ie,ir,it,1)

           else
            trSUMTAL(ip,ie,ir,it,1) =
     &             ( sumfactor * weightRate(ntf) *
     &             trRES(ip,ie,ir,it,1) )

            trSUMTAL(ip,ie,ir,it,2) =
     &           ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie,ir,it,2) )

           end if

          end if

          if ( itpat(m,ip,1) .eq. 20 ) then
           if ( ntf == nfile ) then
            if ( trSUMTAL(ip,ie,ir,it,1) > 0.0d0 ) then

            ! 1 / sqrt( Sigma Nj Xk_bar_j )
             trSUMTAL(ip,ie,ir,it,2) = 1.0d0 / sqrt(
     &       trSUMTAL(ip,ie,ir,it,1))

            ! X_bar k = F / sigma Nj * Sigma Nj Xk_bar_j
             trSUMTAL(ip,ie,ir,it,1) = sumfactor *
     &            trSUMTAL(ip,ie,ir,it,1)

             if( trSUMTAL(ip,ie,ir,it,1) > cmax )
     &             cmax = trSUMTAL(ip,ie,ir,it,1)
             if( trSUMTAL(ip,ie,ir,it,1) < cmin )
     &             cmin = trSUMTAL(ip,ie,ir,it,1)

            else
             trSUMTAL(ip,ie,ir,it,2) = 0.0d0
            end if

           end if
          end if

         end do   ! ip loop end
        end do    ! ie loop end
       end do     ! it loop end
      end do      ! ir loop end

! sumover
      calcfact1 =  sumfactor * weightRate(ntf)
      calcfact2 = sumfactor**2 * weightRate(ntf)**2
      call sumtal_calc_stdev2_depstreg_sub(m, ntf,
     &            nfile, np, itlmax,itmxpt,itpat,
     &            calcfact1, calcfact2, sumfactor)

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = resc2SUMTAL
     &            + weightRate(ntf) * resc3(m) * resc2(m)/resc3(m)

      resc2(m) = resc2SUMTAL
      resc3(m) = resc3SUMTAL

!-----------------------------------------------------------------------
      if( ntf > 1 ) then
       if( maxcasd /= maxcasres ) then
        ierr = 1
        m_err = 'maxcas is inconsistent.'
        ErrCha = ''
        ErrID = 'L:1515/R:sumtal_calc_stdev2_depstreg/F:sumtaldeposit.f' !E54_002_004


        write(*,'(/'' ***** Error in sumtal calc stdev2 depstreg *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev2_depstreg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev2_depstrz(m,ntf,nfile,
     &                             np,nei,ne,nt,nr,nz,
     &                             sumfactor,weightRate,
     &                             trRES,trSUMTAL,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by T.Miura on 2015/07/31                                 *
!                                                                      *
!***********************************************************************
      use partmod, only: itpan, itpat, jtpat, itmxpt ! frtati 2021/10/05
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
      integer           m             ! tally number
      integer           ntf           ! Count of tally file
      integer           nfile         ! Number of tally files

      integer           np
      integer           ne, nei
      integer           nt
      integer           nr
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,0:ne,nt,nr,nz,2)
      double precision  trSUMTAL(np,0:ne,nt,nr,nz,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1, calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir
      integer           iz

      integer           ip
      integer           ie
      integer           it

      integer,save   :: maxcasd
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      ! N = N1 + N2
      resc3SUMTAL = resc3SUMTAL + resc3(m)

      do iz = 1, nz
       do ir = 1, nr
        do it = 1, nt
         do ie = nei, ne
          do ip = 1, np

           if( ntf > 1 ) then
            if ( itpat(m,ip,1) .eq. 20 ) then
            ! Sigma Nj Xk_bar_j
             trSUMTAL(ip,ie,it,ir,iz,1) =
     &              trSUMTAL(ip,ie,it,ir,iz,1) +
     &              trRES(ip,ie,it,ir,iz,1)

            else
            ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
            !                 + F r2 sig_N2 xi wi
             trSUMTAL(ip,ie,it,ir,iz,1) =
     &              trSUMTAL(ip,ie,it,ir,iz,1) +
     &              ( sumfactor * weightRate(ntf) *
     &              trRES(ip,ie,it,ir,iz,1) )

            ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
            !                    + F**2 r2**2 sig_N2(xiwi)**2
             trSUMTAL(ip,ie,it,ir,iz,2) =
     &            trSUMTAL(ip,ie,it,ir,iz,2) +
     &            ( sumfactor**2 * weightRate(ntf)**2 *
     &            trRES(ip,ie,it,ir,iz,2) )

            end if

           ! ntf = 1
           else
            if ( itpat(m,ip,1) .eq. 20 ) then
             trSUMTAL(ip,ie,it,ir,iz,1) =
     &              trRES(ip,ie,it,ir,iz,1)

            else
             trSUMTAL(ip,ie,it,ir,iz,1) =
     &              ( sumfactor * weightRate(ntf) *
     &              trRES(ip,ie,it,ir,iz,1) )

             trSUMTAL(ip,ie,it,ir,iz,2) =
     &            ( sumfactor**2 * weightRate(ntf)**2 *
     &            trRES(ip,ie,it,ir,iz,2) )

            end if

           end if

          if ( itpat(m,ip,1) .eq. 20 ) then
           if ( ntf == nfile ) then
            if ( trSUMTAL(ip,ie,it,ir,iz,1) > 0.0d0 ) then

             ! 1 / sqrt( Sigma Nj Xk_bar_j )
             trSUMTAL(ip,ie,it,ir,iz,2) = 1.0d0 / sqrt(
     &       trSUMTAL(ip,ie,it,ir,iz,1))

             ! X_bar k = F / sigma Nj * Sigma Nj Xk_bar_j
             trSUMTAL(ip,ie,it,ir,iz,1) = sumfactor *
     &       trSUMTAL(ip,ie,it,ir,iz,1)

             if( trSUMTAL(ip,ie,it,ir,iz,1) > cmax )
     &             cmax = trSUMTAL(ip,ie,it,ir,iz,1)
             if( trSUMTAL(ip,ie,it,ir,iz,1) < cmin )
     &             cmin = trSUMTAL(ip,ie,it,ir,iz,1)

            else
             trSUMTAL(ip,ie,it,ir,iz,2) = 0.0d0
            end if

           end if
          end if

          end do   ! ip loop end
         end do    ! ie loop end
        end do     ! it loop end
       end do      ! ir loop end
      end do       ! iz loop end

! sumover
      calcfact1 =  sumfactor * weightRate(ntf)
      calcfact2 = sumfactor**2 * weightRate(ntf)**2
      call sumtal_calc_stdev2_depstreg_sub(m, ntf,
     &            nfile, np, itlmax,itmxpt,itpat,
     &            calcfact1, calcfact2, sumfactor)

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = resc2SUMTAL
     &            + weightRate(ntf) * resc3(m) * resc2(m)/resc3(m)

      resc2(m) = resc2SUMTAL
      resc3(m) = resc3SUMTAL

!-----------------------------------------------------------------------
      if( ntf > 1 ) then
       if( maxcasd /= maxcasres ) then
        ierr = 1
        m_err = 'maxcas is inconsistent.'
        ErrCha = ''
        ErrID = 'L:1708/R:sumtal_calc_stdev2_depstrz/F:sumtaldeposit.f' !E54_002_005


        write(*,'(/'' ***** Error in sumtal calc stdev2 depstrz *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev2_depstrz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev2_depstxyz(m,ntf,nfile,
     &                             np,nei,ne,nt,nx,ny,nz,
     &                             sumfactor,weightRate,
     &                             trRES,trSUMTAL,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by T.Miura on 2015/07/31                                 *
!                                                                      *
!***********************************************************************
      use partmod, only: itpan, itpat, jtpat, itmxpt ! frtati 2021/10/05
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
      integer           m             ! tally number
      integer           ntf           ! Count of tally file
      integer           nfile         ! Number of tally files

      integer           np
      integer           ne, nei
      integer           nt
      integer           nx
      integer           ny
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,0:ne,nt,nx*ny*nz,2)
      double precision  trSUMTAL(np,0:ne,nt,nx*ny*nz,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1, calcfact2

      integer           ierr
!-----------------------------------------------------------------------
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

      ! N = N1 + N2
      resc3SUMTAL = resc3SUMTAL + resc3(m)

      do iz = 1, nz
       do iy = 1, ny
        do ix = 1, nx
         do it = 1, nt
          do ie = nei, ne
           do ip = 1, np

            if( ntf > 1 ) then
             if ( itpat(m,ip,1) .eq. 20 ) then
             ! Sigma Nj Xk_bar_j
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) =
     &               trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) +
     &               trRES(ip,ie,it,icf(ix,iy,iz),1)

            ! ntf = 1

             else
              ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
              !                 + F r2 sig_N2 xi wi
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) =
     &               trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) +
     &               ( sumfactor * weightRate(ntf) *
     &               trRES(ip,ie,it,icf(ix,iy,iz),1) )

              ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
              !                    + F**2 r2**2 sig_N2(xiwi)**2
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) =
     &             trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) +
     &             ( sumfactor**2 * weightRate(ntf)**2 *
     &             trRES(ip,ie,it,icf(ix,iy,iz),2) )

             end if

            else
             if ( itpat(m,ip,1) .eq. 20 ) then
             ! Sigma Nj Xk_bar_j
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) =
     &               trRES(ip,ie,it,icf(ix,iy,iz),1)


             else
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) =
     &               ( sumfactor * weightRate(ntf) *
     &               trRES(ip,ie,it,icf(ix,iy,iz),1) )

              trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) =
     &             ( sumfactor**2 * weightRate(ntf)**2 *
     &             trRES(ip,ie,it,icf(ix,iy,iz),2) )

             end if

            end if

          if ( itpat(m,ip,1) .eq. 20 ) then
           if ( ntf == nfile ) then
            if ( trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) > 0.0d0 ) then

              ! 1 / sqrt( Sigma Nj Xk_bar_j )
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) = 1.0d0 / sqrt(
     &              trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) )

              ! X_bar k = F / sigma Nj * Sigma Nj Xk_bar_j
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)=sumfactor*
     &            trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)

             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) > cmax )
     &            cmax = trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)
             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) < cmin )
     &            cmin = trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)

            else
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) = 0.0d0
            end if

            end if
          end if

           end do   ! ip loop end
          end do    ! ie loop end
         end do     ! it loop end
        end do      ! ix loop end
       end do       ! iy loop end
      end do        ! iz loop end

! sumover
      calcfact1 =  sumfactor * weightRate(ntf)
      calcfact2 = sumfactor**2 * weightRate(ntf)**2
      call sumtal_calc_stdev2_depstreg_sub(m, ntf,
     &            nfile, np, itlmax,itmxpt,itpat,
     &            calcfact1, calcfact2, sumfactor)

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = resc2SUMTAL
     &            + weightRate(ntf) * resc3(m) * resc2(m)/resc3(m)

      resc2(m) = resc2SUMTAL
      resc3(m) = resc3SUMTAL

!-----------------------------------------------------------------------
      if( ntf > 1 ) then
       if( maxcasd /= maxcasres ) then
        ierr = 1
        m_err = 'maxcas is inconsistent.'
        ErrCha = ''
        ErrID = 'L:1912/R:sumtal_calc_stdev2_depstxyz/F:sumtaldeposit.f' !E54_002_006


        write(*,'(/'' ***** Error in sumtal calc stdev2 depstxyz *****''
     &            /)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev2_depstxyz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average2_depstreg(m,ntf,nfile,
     &                               np,nei,ne,nt,nr,
     &                               sumfactor,weightRate,sumWR,
     &                               trRES,trSUMTAL,
     &                               resc2SUMTAL,resc3SUMTAL,
     &                               ierr)
!                                                                      *
!     created by T.Miura on 2015/07/31                                 *
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
      integer           m             ! tally number
      integer           ntf           ! Count of tally file
      integer           nfile         ! Number of tally files

      integer           np
      integer           ne, nei
      integer           nt
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (np,0:ne,nr,nt,2)
      double precision  trSUMTAL(np,0:ne,nr,nt,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir

      integer           ip
      integer           ie
      integer           it

      integer           nf

      double precision  Xa, sigx

      integer   nr0
      integer   itrwgtsum
      integer   itrgn1
      integer   itrgm1
      integer   itrncd
      common /tall70/ itrwgtsum(itlmax), itrgn1(itlmax),
     &                itrgm1(itlmax), itrncd(itlmax)

!-----------------------------------------------------------------------

      ierr = 0

      if( itrwgtsum(m) .eq. 1 )then
       nr0 = 1
      else
       nr0 = nr
      endif

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr0
       do it = 1, nt
        do ie = nei, ne
         do ip = 1, np

          if( trRES(ip,ie,ir,it,1) .gt. 0.d0 ) then
           call calc_deposit_stdev(m,Xa,sigx,
     &            trRES(ip,ie,ir,it,1),
     &            trRES(ip,ie,ir,it,2),
     &            1.0d+0,ip)
           trRES(ip,ie,ir,it,1) = Xa
           trRES(ip,ie,ir,it,2) = sigx
          end if

          if( ntf > 1 ) then
           ! Sigma (rj/r) Xk_bar_j
           trSUMTAL(ip,ie,ir,it,1) =
     &     trSUMTAL(ip,ie,ir,it,1) +
     &                 weightRate(ntf) / sumWR
     &      * trRES(ip,ie,ir,it,1)

           ! Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2
           trSUMTAL(ip,ie,ir,it,2) =
     &     trSUMTAL(ip,ie,ir,it,2) +
     &                  (weightRate(ntf) / sumWR)**2
     &      *(trRES(ip,ie,ir,it,2)
     &      * trRES(ip,ie,ir,it,1))**2

          ! ntf = 1
          else
           ! Sigma (rj/r) Xk_bar_j
           trSUMTAL(ip,ie,ir,it,1) =
     &                 weightRate(ntf) / sumWR
     &      * trRES(ip,ie,ir,it,1)

           ! Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2
           trSUMTAL(ip,ie,ir,it,2) =
     &                  (weightRate(ntf) / sumWR)**2
     &      *(trRES(ip,ie,ir,it,2)
     &      * trRES(ip,ie,ir,it,1))**2

          end if

         end do   ! ip loop end
        end do    ! ie loop end
       end do     ! it loop end
      end do      ! ir loop end

! sumover  tr_sum -> trSUMTAL_sum
      call calc_deposit_stdev_tr_sum(m, 1.0d0, np)
      call sumtal_calc_average_sub(m, ntf,weightRate(ntf)/sumWR)

!-----------------------------------------------------------------------

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = 1d0

      ! N = N1 + N2
      resc3SUMTAL = 2d0
!-----------------------------------------------------------------------

      if( ntf == nfile ) then

       do ir = 1, nr0
        do it = 1, nt
         do ie = nei, ne
          do ip = 1, np

           ! X_bar_k = F Sigma (rj/r) Xk_bar_j
           trSUMTAL(ip,ie,ir,it,1) = sumfactor *
     &     trSUMTAL(ip,ie,ir,it,1)

C S.H. added the following processes for negative values when isumtally=3.
C This is tentative. Correct calculation can be done if processes of
C subroutine calc_stdev is performed here, but preparation of the last
C parameter of calc_stdev is needed. (2020.12.16)
           if( trSUMTAL(ip,ie,ir,it,1) < 0d0 ) then
              trSUMTAL(ip,ie,ir,it,1) = 0d0
              trSUMTAL(ip,ie,ir,it,2) = 0d0
           end if

           if( trSUMTAL(ip,ie,ir,it,1) .gt. 0.d0 ) then
            ! F sqrt{ Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2} / X_bar_k
            trSUMTAL(ip,ie,ir,it,2) = sumfactor * sqrt(
     &             trSUMTAL(ip,ie,ir,it,2))/
     &             trSUMTAL(ip,ie,ir,it,1)
           end if

           if( trSUMTAL(ip,ie,ir,it,1) > cmax )
     &           cmax = trSUMTAL(ip,ie,ir,it,1)
           if( trSUMTAL(ip,ie,ir,it,1) < cmin )
     &           cmin = trSUMTAL(ip,ie,ir,it,1)

          end do   ! ip loop end
         end do    ! ie loop end
        end do     ! it loop end
       end do      ! ir loop end

!sumover
      call sumtal_calc_average2_sub2(m, ntf, sumfactor)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average2_depstreg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average2_depstrz(m,ntf,nfile,
     &                               np,nei,ne,nt,nr,nz,
     &                               sumfactor,weightRate,sumWR,
     &                               trRES,trSUMTAL,
     &                               resc2SUMTAL,resc3SUMTAL,
     &                               ierr)
!                                                                      *
!     created by T.Miura on 2015/07/31                                 *
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
      integer           m             ! tally number
      integer           ntf           ! Count of tally file
      integer           nfile         ! Number of tally files

      integer           np
      integer           ne, nei
      integer           nt
      integer           nr
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (np,0:ne,nt,nr,nz,2)
      double precision  trSUMTAL(np,0:ne,nt,nr,nz,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir
      integer           iz

      integer           ip
      integer           ie
      integer           it

      integer           nf

      double precision  Xa, sigx
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do iz = 1, nz
       do ir = 1, nr
        do it = 1, nt
         do ie = nei, ne
          do ip = 1, np

           if( trRES(ip,ie,it,ir,iz,1) .gt. 0.d0 ) then
            call calc_deposit_stdev(m,Xa,sigx,
     &             trRES(ip,ie,it,ir,iz,1),
     &             trRES(ip,ie,it,ir,iz,2),
     &             1.0d+0,ip)
            trRES(ip,ie,it,ir,iz,1) = Xa
            trRES(ip,ie,it,ir,iz,2) = sigx
           end if

           if( ntf > 1 ) then
            ! Sigma (rj/r) Xk_bar_j
            trSUMTAL(ip,ie,it,ir,iz,1) =
     &      trSUMTAL(ip,ie,it,ir,iz,1) +
     &                     weightRate(ntf) / sumWR
     &       * trRES(ip,ie,it,ir,iz,1)

            ! Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2
            trSUMTAL(ip,ie,it,ir,iz,2) =
     &      trSUMTAL(ip,ie,it,ir,iz,2) +
     &                    (weightRate(ntf) / sumWR)**2
     &       *(trRES(ip,ie,it,ir,iz,2)
     &       * trRES(ip,ie,it,ir,iz,1))**2

           ! ntf = 1
           else
            ! Sigma (rj/r) Xk_bar_j
            trSUMTAL(ip,ie,it,ir,iz,1) =
     &                     weightRate(ntf) / sumWR
     &       * trRES(ip,ie,it,ir,iz,1)

            ! Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2
            trSUMTAL(ip,ie,it,ir,iz,2) =
     &                    (weightRate(ntf) / sumWR)**2
     &       *(trRES(ip,ie,it,ir,iz,2)
     &       * trRES(ip,ie,it,ir,iz,1))**2

           end if

          end do   ! ip loop end
         end do    ! ie loop end
        end do     ! it loop end
       end do      ! ir loop end
      end do       ! iz loop end

! sumover  tr_sum -> trSUMTAL_sum
      call calc_deposit_stdev_tr_sum(m, 1.0d0, np)
      call sumtal_calc_average_sub(m, ntf,weightRate(ntf)/sumWR)

!-----------------------------------------------------------------------

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = 1d0

      ! N = N1 + N2
      resc3SUMTAL = 2d0
!-----------------------------------------------------------------------
      if( ntf == nfile ) then

       do iz = 1, nz
        do ir = 1, nr
         do it = 1, nt
          do ie = nei, ne
           do ip = 1, np

            ! X_bar_k = F Sigma (rj/r) Xk_bar_j
            trSUMTAL(ip,ie,it,ir,iz,1) = sumfactor *
     &      trSUMTAL(ip,ie,it,ir,iz,1)

C S.H. added the following processes for negative values when isumtally=3.
C This is tentative. Correct calculation can be done if processes of
C subroutine calc_stdev is performed here, but preparation of the last
C parameter of calc_stdev is needed. (2020.12.16)
            if( trSUMTAL(ip,ie,it,ir,iz,1) < 0d0 ) then
               trSUMTAL(ip,ie,it,ir,iz,1) = 0d0
               trSUMTAL(ip,ie,it,ir,iz,2) = 0d0
            end if

            if( trSUMTAL(ip,ie,it,ir,iz,1) .gt. 0.d0 ) then
            ! F sqrt{ Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2} / X_bar_k
             trSUMTAL(ip,ie,it,ir,iz,2) = sumfactor * sqrt(
     &              trSUMTAL(ip,ie,it,ir,iz,2))/
     &              trSUMTAL(ip,ie,it,ir,iz,1)
            end if

            if( trSUMTAL(ip,ie,it,ir,iz,1) > cmax )
     &           cmax = trSUMTAL(ip,ie,it,ir,iz,1)
            if( trSUMTAL(ip,ie,it,ir,iz,1) < cmin )
     &           cmin = trSUMTAL(ip,ie,it,ir,iz,1)

           end do   ! ip loop end
          end do    ! ie loop end
         end do     ! it loop end
        end do      ! ir loop end
       end do       ! iz loop end

!sumover
      call sumtal_calc_average2_sub2(m, ntf, sumfactor)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average2_depstrz



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average2_depstxyz(m,ntf,nfile,
     &                               np,nei,ne,nt,nx,ny,nz,
     &                               sumfactor,weightRate,sumWR,
     &                               trRES,trSUMTAL,
     &                               resc2SUMTAL,resc3SUMTAL,
     &                               ierr)
!                                                                      *
!     created by T.Miura on 2015/07/31                                 *
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
      integer           m             ! tally number
      integer           ntf           ! Count of tally file
      integer           nfile         ! Number of tally files

      integer           np
      integer           ne, nei
      integer           nt
      integer           nx
      integer           ny
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (np,0:ne,nt,nx*ny*nz,2)
      double precision  trSUMTAL(np,0:ne,nt,nx*ny*nz,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ix
      integer           iy
      integer           iz

      integer           ip
      integer           ie
      integer           it

      integer           nf

      double precision  Xa, sigx
!-----------------------------------------------------------------------
      integer           icf

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do iz = 1, nz
       do iy = 1, ny
        do ix = 1, nx
         do it = 1, nt
          do ie = nei, ne
           do ip = 1, np

            if( trRES(ip,ie,it,icf(ix,iy,iz),1) .gt. 0.d0 ) then
             call calc_deposit_stdev(m,Xa,sigx,
     &              trRES(ip,ie,it,icf(ix,iy,iz),1),
     &              trRES(ip,ie,it,icf(ix,iy,iz),2),
     &              1.0d+0,ip)
             trRES(ip,ie,it,icf(ix,iy,iz),1) = Xa
             trRES(ip,ie,it,icf(ix,iy,iz),2) = sigx
            end if

            if( ntf > 1 ) then
             ! Sigma (rj/r) Xk_bar_j
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) =
     &       trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) +
     &                          weightRate(ntf) / sumWR
     &        * trRES(ip,ie,it,icf(ix,iy,iz),1)

             ! Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) =
     &       trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) +
     &                         (weightRate(ntf) / sumWR)**2
     &        *(trRES(ip,ie,it,icf(ix,iy,iz),2)
     &        * trRES(ip,ie,it,icf(ix,iy,iz),1))**2

            ! ntf = 1
            else
             ! Sigma (rj/r) Xk_bar_j
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) =
     &                          weightRate(ntf) / sumWR
     &        * trRES(ip,ie,it,icf(ix,iy,iz),1)

             ! Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) =
     &                         (weightRate(ntf) / sumWR)**2
     &        *(trRES(ip,ie,it,icf(ix,iy,iz),2)
     &        * trRES(ip,ie,it,icf(ix,iy,iz),1))**2

            end if

           end do   ! ip loop end
          end do    ! ie loop end
         end do     ! it loop end
        end do      ! ix loop end
       end do       ! iy loop end
      end do        ! iz loop end

! sumover  tr_sum -> trSUMTAL_sum
      call calc_deposit_stdev_tr_sum(m, 1.0d0, np)
      call sumtal_calc_average_sub(m, ntf,weightRate(ntf)/sumWR)

!-----------------------------------------------------------------------

      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = 1d0

      ! N = N1 + N2
      resc3SUMTAL = 2d0
!-----------------------------------------------------------------------
      if( ntf == nfile ) then

       do iz = 1, nz
        do iy = 1, ny
         do ix = 1, nx
          do it = 1, nt
           do ie = nei, ne
            do ip = 1, np

             ! X_bar_k = F Sigma (rj/r) Xk_bar_j
             trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) = sumfactor *
     &       trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)

C S.H. added the following processes for negative values when isumtally=3.
C This is tentative. Correct calculation can be done if processes of
C subroutine calc_stdev is performed here, but preparation of the last
C parameter of calc_stdev is needed. (2020.12.16)
             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) < 0d0 ) then
                trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) = 0d0
                trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) = 0d0
             end if

             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) .gt. 0.d0 ) then
              ! F sqrt{ Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2} / X_bar_k
              trSUMTAL(ip,ie,it,icf(ix,iy,iz),2) = sumfactor * sqrt(
     &               trSUMTAL(ip,ie,it,icf(ix,iy,iz),2))/
     &               trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)
             end if

             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) > cmax )
     &             cmax = trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)
             if( trSUMTAL(ip,ie,it,icf(ix,iy,iz),1) < cmin )
     &             cmin = trSUMTAL(ip,ie,it,icf(ix,iy,iz),1)

            end do   ! ip loop end
           end do    ! ie loop end
          end do     ! it loop end
         end do      ! ix loop end
        end do       ! iy loop end
       end do        ! iz loop end

!sumover
      call sumtal_calc_average2_sub2(m, ntf, sumfactor)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average2_depstxyz

!***********************************************************************
!                                                                      *
! sumover sunroutine                                                   *
!                                                                      *
!***********************************************************************

!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev2_depstreg_sub(m, ntf,
     &           nfile1, np, itlmax,itmxpt,itpat,
     &           calcfact1,calcfact2,sumfactor1)
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

      integer :: m, iax, ntf, nfile1, np
      integer :: itlmax,itmxpt,itpat(itlmax,itmxpt,3)

      double precision :: calcfact(2),calcfact1,calcfact2,sumfactor1

      real(8),pointer :: p_sum(:)
!      real(8),pointer :: resta_sum(:)
      real(8),pointer :: sumtal_sum(:)

      integer :: ln_sum

       calcfact(1) = calcfact1
       calcfact(2) = calcfact2

      do iax =1,6
!       resta_sum => trRES_sum(irestalm_sum(m,iax):)

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(p_sum,m,iax)
!$       else
            call GET_TR_HEAD_POINTER_SUM(p_sum,m,iax)
C for nonshared_tally option
!$       end if

        sumtal_sum => trSUMTAL_SUM(irestalm_sum(m,iax):)
        ln_sum = lrestalm_sum(m,iax)/2

        if(ln_sum > 0 ) then

          call calc_stdev2_restal2sumtally(p_sum,sumtal_sum,
     &       m, nfile1, np, itlmax,itmxpt,itpat,
     &       ln_sum, ntf, calcfact, sumfactor1)
        endif
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine calc_stdev2_restal2sumtally(
     &       trRES_SUM,trSUMTAL_SUM,
     &       m, nfile, np,itlmax,itmxpt,itpat,
     &       ln_sum, ntf, calcfact, sumfactor)
!                                                                      *
!***********************************************************************

      implicit none

      integer :: m,nfile,np,ln_sum,ntf,i,isum,ip
      real(8) :: trRES_sum(ln_sum,2),trSUMTAL_SUM(ln_sum,2),calcfact(2)
      real(8) :: sumfactor
      integer :: itlmax,itmxpt,itpat(itlmax,itmxpt,3)

      ip = 0
      do isum=1,ln_sum
        ip = ip + 1
        if(ip > np) then
           ip= 1
        endif

        if ( itpat(m,ip,1) .eq. 20 ) then
            if( ntf > 1 ) then
              trSUMTAL_SUM(isum,1) =  trSUMTAL_SUM(isum,1)
     &                     +  trRES_SUM(isum,1)
            else
              trSUMTAL_SUM(isum,1) =  trRES_SUM(isum,1)
            endif

        else
           do i=1,2
             if( ntf > 1 ) then


                trSUMTAL_SUM(isum,i) =  trSUMTAL_SUM(isum,i)
     &                     +  calcfact(i) * trRES_SUM(isum,i)

             else
               trSUMTAL_SUM(isum,i) = calcfact(i) * trRES_SUM(isum,i)
             end if
           enddo
        endif

        if ( itpat(m,ip,1) .eq. 20 ) then
          if ( ntf == nfile ) then

            if ( trSUMTAL_SUM(isum,1) > 0.0d0 ) then

              trSUMTAL_SUM(isum,2) = 1.0d0 / sqrt(trSUMTAL_SUM(isum,1))

              trSUMTAL_SUM(isum,1) = sumfactor * trSUMTAL_SUM(isum,1)

            else
              trSUMTAL_SUM(isum,2) = 0.0d0

            end if

          end if
        end if
      enddo

      return
      end
