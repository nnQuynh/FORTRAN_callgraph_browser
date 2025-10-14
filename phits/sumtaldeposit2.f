************************************************************************
*                                                                      *
      subroutine sumtal_read_tdeposit2(m,iax,ntf,ierr)
!     original routine is read_tdeposit2 in restdeposit2.f             *
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
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
        common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                  rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)

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

        if ( any( itaxs(m,iax) .eq. (/ 3, 4, 5, 6, 7, 8 /) )
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
        call check_tdeposit2(m,iax,jsn(ioe),jsi(ioe),
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
          idas2 = idas1 + ittnm(m)
          idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
          idasa = idas2 + itrgn(m)
     &          + ( itrgn(m) + mod(itrgn(m),2) ) / 2

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

          call read_deposit2reg(m,iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itpan(m),itrgn(m),itrgm(m),itenm(m),itenm2(m),
     &      ittnm(m),idas_itreg(itreg(m)),
     &      das_iterg(iterg(m)),das_iterg2(iterg2(m)),
     &      das_ittrg(ittrg(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)))

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

          call restore_deposit2reg(m,
     &      itpan(m),itrgn(m),itrgm(m),itenm(m),itenm2(m),
     &      ittnm(m),idas_itreg(itreg(m)),
     &      das_iterg(iterg(m)),das_iterg2(iterg2(m)),
     &      das_ittrg(ittrg(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)))

  900   continue
*-----------------------------------------------------------------------
*   close restart file
*-----------------------------------------------------------------------

        do ioe = 1, noe

          close(jsi(ioe))

        end do

*-----------------------------------------------------------------------
  999   continue

      end subroutine sumtal_read_tdeposit2


!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev2_deposit2reg(m,ntf,nfile,
     &                             np,nt,ne1,ne2,
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
      integer           m
      integer           ntf
      integer           nfile

      integer           np
      integer           nt
      integer           ne1
      integer           ne2

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,0:ne1,0:ne2,nt,2)
      double precision  trSUMTAL(np,0:ne1,0:ne2,nt,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

! sumover
      double precision  calcfact1, calcfact2

      integer           ierr
!-----------------------------------------------------------------------
      integer           ie1
      integer           ie2
      integer           it
      integer           ip

      integer,save   :: maxcasd
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      ! N = Sigma Nj
      resc3SUMTAL = resc3SUMTAL + resc3(m)

      do ie1 = 1, ne1
       do ie2 = 1, ne2
        do it = 1, nt
         do ip = 1, np

          if( ntf > 1 ) then
           if ( itpat(m,ip,1) .eq. 20 ) then
           ! Sigma Nj Xk_bar_j
            trSUMTAL(ip,ie1,ie2,it,1) =
     &             trSUMTAL(ip,ie1,ie2,it,1) +
     &             trRES(ip,ie1,ie2,it,1)

           else
            ! sig_N1+N2 xi wi = F r1 sig_N1 xi wi
            !                 + F r2 sig_N2 xi wi
            trSUMTAL(ip,ie1,ie2,it,1) =
     &             trSUMTAL(ip,ie1,ie2,it,1) +
     &             ( sumfactor * weightRate(ntf) *
     &             trRES(ip,ie1,ie2,it,1) )

            ! sig_N1+N2(xiwi)**2 = F**2 r1**2 sig_N1(xiwi)**2
            !                    + F**2 r2**2 sig_N2(xiwi)**2
            trSUMTAL(ip,ie1,ie2,it,2) =
     &           trSUMTAL(ip,ie1,ie2,it,1) +
     &           ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie1,ie2,it,2) )

           end if

          ! ntf = 1
          else
           if ( itpat(m,ip,1) .eq. 20 ) then
            ! Sigma Nj Xk_bar_j
            trSUMTAL(ip,ie1,ie2,it,1) =
     &             trRES(ip,ie1,ie2,it,1)

           else
            trSUMTAL(ip,ie1,ie2,it,1) =
     &             ( sumfactor * weightRate(ntf) *
     &             trRES(ip,ie1,ie2,it,1) )

            trSUMTAL(ip,ie1,ie2,it,2) =
     &           ( sumfactor**2 * weightRate(ntf)**2 *
     &           trRES(ip,ie1,ie2,it,2) )

           end if

          end if

          if ( itpat(m,ip,1) .eq. 20 ) then
           if ( ntf == nfile ) then
            if ( trSUMTAL(ip,ie1,ie2,it,1) > 0.0d0 ) then

             ! 1 / sqrt( Sigma Nj Xk_bar_j )
             trSUMTAL(ip,ie1,ie2,it,2) = 1.0d0 / sqrt(
     &              trSUMTAL(ip,ie1,ie2,it,1))

             ! X_bar k = F / sigma Nj * Sigma Nj Xk_bar_j
             trSUMTAL(ip,ie1,ie2,it,1) = sumfactor *
     &            trSUMTAL(ip,ie1,ie2,it,1)

             if( trSUMTAL(ip,ie1,ie2,it,1) > cmax )
     &            cmax = trSUMTAL(ip,ie1,ie2,it,1)
             if( trSUMTAL(ip,ie1,ie2,it,1) < cmin )
     &            cmin = trSUMTAL(ip,ie1,ie2,it,1)

            else
             trSUMTAL(ip,ie1,ie2,it,2) = 0.0d0
            end if

           end if
          end if

         end do      ! ip  loop end
        end do       ! it  loop end
       end do        ! ie2 loop end
      end do         ! ie1 loop end

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
        ErrID = 'L:333/R:sumtal_calc_stdev2_deposit2reg/F:sumtaldepos'//
     &          'it2.f' !E55_001_001


        write(*,'(/
     &      '' ***** Error in sumtal calc stdev2 deposit2reg *****''/)')
        call ErrWrite(ErrID, ErrCha)
        write(*,'('' error = '',(a)/)') trim(m_err)
       end if
      end if
      maxcasd = maxcasres
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev2_deposit2reg



!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_averag2_deposit2reg(m,ntf,nfile,
     &                               np,nt,ne1,ne2,
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
      integer           m
      integer           ntf
      integer           nfile

      integer           np
      integer           nt
      integer           ne1
      integer           ne2

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  sumWR         ! Sum of the weight rates
      double precision  trRES   (np,0:ne1,0:ne2,nt,2)
      double precision  trSUMTAL(np,0:ne1,0:ne2,nt,2)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ie1
      integer           ie2
      integer           it
      integer           ip

      integer           nf

      double precision  Xa, sigx
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ie1 = 1, ne1
       do ie2 = 1, ne2
        do it = 1, nt
         do ip = 1, np

          if( trRES(ip,ie1,ie2,it,1) .gt. 0.d0 ) then
           call calc_deposit_stdev(m,Xa,sigx,
     &                     trRES(ip,ie1,ie2,it,1),
     &                     trRES(ip,ie1,ie2,it,2),
     &                     1.0d+0,ip)
           trRES(ip,ie1,ie2,it,1) = Xa
           trRES(ip,ie1,ie2,it,2) = sigx
          end if

          if( ntf > 1 ) then
           ! Sigma (rj/r) Xk_bar_j
           trSUMTAL(ip,ie1,ie2,it,1) =
     &     trSUMTAL(ip,ie1,ie2,it,1) +
     &                   weightRate(ntf) / sumWR
     &      * trRES(ip,ie1,ie2,it,1)

           ! Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2
           trSUMTAL(ip,ie1,ie2,it,2) =
     &     trSUMTAL(ip,ie1,ie2,it,2) +
     &                  (weightRate(ntf) / sumWR)**2
     &      *(trRES(ip,ie1,ie2,it,2)
     &      * trRES(ip,ie1,ie2,it,1))**2

          ! ntf = 1
          else
           ! Sigma (rj/r) Xk_bar_j
           trSUMTAL(ip,ie1,ie2,it,1) =
     &                   weightRate(ntf) / sumWR
     &      * trRES(ip,ie1,ie2,it,1)

           ! Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2
           trSUMTAL(ip,ie1,ie2,it,2) =
     &                  (weightRate(ntf) / sumWR)**2
     &      *(trRES(ip,ie1,ie2,it,2)
     &      * trRES(ip,ie1,ie2,it,1))**2

          end if

         end do      ! ip  loop end
        end do       ! it  loop end
       end do        ! ie2 loop end
      end do         ! ie1 loop end

! sumover  tr_sum -> trSUMTAL_sum
      call calc_deposit_stdev_tr_sum(m, 1.0d0, np)
      call sumtal_calc_average_sub(m, ntf,weightRate(ntf)/sumWR)

!-----------------------------------------------------------------------
  100 continue
      ! W = r1 N1 w1_bar + r2 N2 w2_bar
      resc2SUMTAL = resc2SUMTAL
     &            + weightRate(ntf) * resc3(m) * resc2(m)/resc3(m)

      ! N = N1 + N2
      resc3SUMTAL = resc3SUMTAL + resc3(m)
!-----------------------------------------------------------------------
      if( ntf == nfile ) then

       do ie1 = 1, ne1
        do ie2 = 1, ne2
         do it = 1, nt
          do ip = 1, np

           ! X_bar_k = F Sigma (rj/r) Xk_bar_j
           trSUMTAL(ip,ie1,ie2,it,1) = sumfactor *
     &     trSUMTAL(ip,ie1,ie2,it,1)

           if ( trSUMTAL(ip,ie1,ie2,it,1) > 0.0d0 ) then
            ! F sqrt{ Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2} / X_bar_k
            trSUMTAL(ip,ie1,ie2,it,2) = sumfactor * sqrt(
     &      trSUMTAL(ip,ie1,ie2,it,2))/
     &      trSUMTAL(ip,ie1,ie2,it,1)

            if( trSUMTAL(ip,ie1,ie2,it,1) > cmax )
     &            cmax = trSUMTAL(ip,ie1,ie2,it,1)
            if( trSUMTAL(ip,ie1,ie2,it,1) < cmin )
     &            cmin = trSUMTAL(ip,ie1,ie2,it,1)
           end if

          end do      ! ip  loop end
         end do       ! it  loop end
        end do        ! ie2 loop end
       end do         ! ie1 loop end

!sumover
      call sumtal_calc_average2_sub2(m, ntf, sumfactor)

       resc2(m) = resc2SUMTAL
       resc3(m) = resc3SUMTAL

      end if
!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_averag2_deposit2reg
