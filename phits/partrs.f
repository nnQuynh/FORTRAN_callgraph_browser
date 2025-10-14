************************************************************************
*                                                                      *
      subroutine partrs(mark,markp,nbeta,icge,itmak,emint)
                                                  !(emint:Takeshi Kai) *
*                                                                      *
*       particle transport                                             *
*       and region check                                               *
*       modified by K.Niita on 2003/10/12                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       icge   : error code                                            *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
************************************************************************
      use ion_track_structure, only: lflgTS
      use MMBANKMOD !FURUTA
      use MEMBANKMOD !FURUTA
      use ELEDATAMOD, only : ichem
!<-20220128murofushi add
      use moddas_region
      use moddas_ggs
      use ets_art, only: etstrn_art
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      common /eparm/  esmax, esmin, emin(20)
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      common /spsgn/  nspsgn
      common /tcntl/  icntl, inucr
      common /paraj/  mstz(300), parz(300)

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      common /inout/  in,io
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /cgerr/  nlost, ilost, igerr, icger, ncger, nrecover
      common /cgstr/  novp, nrovp(3,1000)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /gravit/ grav(3), igrav

*-----------------------------------------------------------------------

      common /elmgrg/ nelctf, nereg, inerc, inert, kelcs

      common /insprd/ initsp
!$OMP THREADPRIVATE(/insprd/)
      common /magreg/ nmreg, ingrc, ingrt, kmags, magtin, imgusr
      common /srsph/  pmphs
!$OMP THREADPRIVATE(/srsph/)
      common /argcns/ kcar, kczs, kcze

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      data initsp /0/

*-----------------------------------------------------------------------

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

      common /ndemax/dnmax(20)

      common /tdcyl/ tdcylife, itdcyw
!$OMP THREADPRIVATE(/tdcyl/)

      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax
      common / tsminmax  / tsmax

      common /wwin00/ iwwxyz
!$OMP THREADPRIVATE(/wwin00/)

      common /adjoint/ iadjnt ! T.Sato 2024/01/07

*-----------------------------------------------------------------------
*     initialization of spread
*-----------------------------------------------------------------------

         if( initsp .eq. 0 .and. nspred .ne. 0 ) then

            initsp = initsp + 1

            call sprdint

         end if

*-----------------------------------------------------------------------
*     initialization of time frag
*-----------------------------------------------------------------------

            itmak = 0

            itdcyw = 0     ! S.Abe 2015/08/10

*-----------------------------------------------------------------------
*     initialization of weight window frag
*-----------------------------------------------------------------------

            iwwxyz = 0

*-----------------------------------------------------------------------
*     initial jcoll = 0
*-----------------------------------------------------------------------

            jcoll = 0

*-----------------------------------------------------------------------
*     If lflgTS is 1, become 2 and energy recovery for restricted LET is conducted
*-----------------------------------------------------------------------

            lflgTS = lflgTS * 2

*-----------------------------------------------------------------------
*        off region particles ( for source particles )
*-----------------------------------------------------------------------

            if( mat .le. -1 ) then

               mark  = -1
               nbeta = 2
               icge  = 0

               ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1)
               xc(ibkxc+no,ipomp+1) = x(ibkx+no,ipomp+1)
               yc(ibkyc+no,ipomp+1) = y(ibky+no,ipomp+1)
               zc(ibkzc+no,ipomp+1) = z(ibkz+no,ipomp+1)

               return

            end if

*-----------------------------------------------------------------------
*        energy cutoff particles
*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

                  emint = emin(ityp)

*-----------------------------------------------------------------------
               if( mntsc .gt. 0 .and.
     &             ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.
     &           ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &           ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &             e(ibke+no,ipomp+1) .ge. etsmin ) )  emint = etsmin

*-----------------------------------------------------------------------

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------

            if( ( e(ibke+no,ipomp+1) .le. emint .and. mat .gt. 0 ) .or.
     &            e(ibke+no,ipomp+1) .eq. 0.0d0 ) then

               mark  = 1
               nbeta = 3
               icge  = 0

               ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1)
               xc(ibkxc+no,ipomp+1) = x(ibkx+no,ipomp+1)
               yc(ibkyc+no,ipomp+1) = y(ibky+no,ipomp+1)
               zc(ibkzc+no,ipomp+1) = z(ibkz+no,ipomp+1)

               return

            end if

*-----------------------------------------------------------------------
*        check initial cell
*-----------------------------------------------------------------------

               ici = 1

               if(mark.eq.1.and.markp.eq.0
     &              .and.itetpos(ibtetpos+no,ipomp+1).ne.0) !FURUTA20160902
     &              iii=itetpos(ibtetpos+no,ipomp+1) !FURUTA20160607

            call gomsor(x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  nmed(ibknmd+no,ipomp+1),
     &                  iblz(ibkblz+no,ipomp+1),
     &                  mark,markp,ici)

               if( mark .le. -2  ) goto 100

               markp = 1

*-----------------------------------------------------------------------
*     getflt : [ fpl ] sampling of flight path length
*-----------------------------------------------------------------------

            if( icntl .ne. 5 .and. icntl .ne. 14 ) then

              call getflt(fpl)

            else

               fpl = 1.0d+20

               goto 30

            end if

*-----------------------------------------------------------------------
*        electron transfer
*-----------------------------------------------------------------------

         if( ( ityp .eq. 12 .or. ityp .eq. 13 )  .and.
     &         mat .ne. 0 .and. fpl .lt. 1.0d+20 .and.
     &         iegsemi .eq. 0 ) then

               call eletr(mark,markp,fpl,nbeta,itmak)

               goto 100

         end if

*-----------------------------------------------------------------------
*        electron transfer by EGS
*-----------------------------------------------------------------------

            if( ( ityp .eq. 12 .or. ityp .eq. 13 .or.
     &            ityp .eq. 14 ) .and.
     &            iegsemi .ne. 0 ) then




              if( nwsprd .gt. 0 .and. nmreg .gt. 0 ) then

                  idsm = ingrt
                  jdsm = 0

                  kdsm = kmags
                  ldsm = 0

                  do i = 1, nmreg

                    jj  = 0

                    jdsm = jdsm + 1
                    ntrn = idas_ingrt(idsm+jdsm)
                    jdsm = jdsm + 1
                    mtrn = idas_ingrt(idsm+jdsm)

                    ldsm  = ldsm + 1
                    a_mag = das_kmags(kdsm+ldsm)
                    ldsm  = ldsm + 1
                    b_mag = das_kmags(kdsm+ldsm)
                    ldsm  = ldsm + 1
                    s_mag = das_kmags(kdsm+ldsm)
                    ldsm  = ldsm + 1
                    p_mag = das_kmags(kdsm+ldsm)
                    ldsm  = ldsm + 1
                    t_mag = das_kmags(kdsm+ldsm)
                    ldsm  = ldsm + 1
                    ldsm  = ldsm + 1
                    u_mag = das_kmags(kdsm+ldsm)

                    do j = 1, ntrn

                     call tregck(iblz1,ilev1,ilat1,
     &                        mtrn,idas_ingrt(idsm+jdsm+1),jj,icc)

                     if( icc .ne. 0 ) then

                       isp = nint( s_mag )

                       !ASTOM 2019/01/08
                       if( ( isp .eq. 2 .or. isp .eq. 4 .or.
     &                       isp .eq. 6 .or. isp .eq. 8 .or.
     &                       isp .eq. -1 .or. isp .eq. -2 .or.
     &                       isp .eq. -3 .or. isp .eq. -4) .and.
     &                       jtyp .ne. 0 ) then

*-----------------------------------------------------------------------
*           Wobbler magnet for charged particles for IHI
*-----------------------------------------------------------------------
                        if( p_mag .gt. -1000.0d0 ) then

                          b_mag = sin( pmphs + p_mag ) * b_mag

                        end if

*-----------------------------------------------------------------------

                        call egs_magfld(mark,markp,fpl,nbeta,itmak,
     &                             a_mag,b_mag,s_mag,p_mag,t_mag,u_mag)

                        goto 100

                       end if

                     endif

                    end do

                    jdsm = jdsm + mtrn

                  end do

              endif

*-----------------------------------------------------------------------
*        charge particle transfer under electro magnetic field
*           nereg    : number of region of electro magnetic field
*-----------------------------------------------------------------------

              if( nelctf .gt. 0 .and. nereg .gt. 0 .and.
     &            jtyp .ne. 0 ) then

                idsm = inert
                jdsm = 0

                kdsm = kelcs
                ldsm = 0

                do i = 1, nereg

                  jj  = 0

                  jdsm = jdsm + 1
                  ntrn = idas_inert(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_inert(idsm+jdsm)

                  ldsm  = ldsm + 1
                  s_elf = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  s_mgf = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  t_elf = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 2
                  t_mgf = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 2
                  emap_type = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  mmap_type = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  a_elmg = das_kelcs(kdsm+ldsm)

                  do j = 1, ntrn

                    call tregck(iblz1,ilev1,ilat1,
     &                          mtrn,idas_inert(idsm+jdsm+1),jj,icc)

                    if( icc .ne. 0 ) then

                      if(s_elf.ne.0.0) then

                         call clrcmn3()
                         call getflt(fpl)

                      endif

                      !AdvanceSoft Hasemi 2019/11/27 2029/12/20
                      call egs_elmgfd(mark,markp,fpl,nbeta,itmak,
     &                                s_elf,s_mgf,t_elf,t_mgf,e_chs,
     &                                emap_type, mmap_type)

                      goto 100

                    endif

                  end do

                  jdsm = jdsm + mtrn

                end do

              endif
*-----------------------------------------------------------------------

                   mat1 = 0
               if( mntsc .gt. 0 .and.
     &             ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.
     &           ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &           ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &             e(ibke+no,ipomp+1) .ge. etsmin ) )
     &              mat1 = ntscell( idgr(iblz(ibkblz+no,ipomp+1)) )

              if(mat1.eq.0)then
                if(ityp.eq.12) then
                 etsoff = emin(12)
                elseif(ityp.eq.13) then
                 etsoff = emin(13)
                else
                 etsoff = emin(14)
                endif
              else
                etsoff = etsmax
              endif

              if((e(ibke+no,ipomp+1).gt.etsoff.and.ityp.eq.12) .or.
     .           (e(ibke+no,ipomp+1).gt.etsoff.and.ityp.eq.13) .or.
     .                                   ityp.eq.14 .or. mat1.eq.0) then

                call egstr(mark,markp,fpl,nbeta,itmak)

              else

! hirata etsart 20230427
               if(mat1.eq.-1 .and.
     &        ichem(mat,1) .eq. 22014 .and.ichem(mat,2) .eq. 0) then
                 call etstrn02(mark,markp,fpl,nbeta,itmak)
               else if(mat1.eq.-1 .and.
     &        ichem(mat,1) .eq. 1 .and.ichem(mat,2) .eq. 0) then
                 call etstrn(mark,markp,fpl,nbeta,itmak)
               else if(mat1.eq.-1 ) then
                 call etstrn_art(mark,markp,fpl,nbeta,itmak)
               else
                call etstrn(mark,markp,fpl,nbeta,itmak)
               endif
              endif

             if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     .             mat1 .ne. 0                      .and.
     .           ( e(ibke+no,ipomp+1) .le. etsmax )         .and.
     .           ( e(ibke+no,ipomp+1) .ge. etsmin )         )then
               emint = etsmin
             else
               emint = emin(ityp) * dble(max(1,ibryf(ityp,ktyp)))
             endif


              goto 100

            end if

*-----------------------------------------------------------------------
*        charge particle or neutron transfer under magnetic field
*           nmreg    : number of region of magnetic field
*-----------------------------------------------------------------------

         if( nwsprd .gt. 0 .and. nmreg .gt. 0 ) then

                  idsm = ingrt
                  jdsm = 0

                  kdsm = kmags
                  ldsm = 0

            do i = 1, nmreg

                  jj  = 0

                  jdsm = jdsm + 1
                  ntrn = idas_ingrt(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_ingrt(idsm+jdsm)

                  ldsm  = ldsm + 1
                  a_mag = das_kmags(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  b_mag = das_kmags(kdsm+ldsm)
                  if(iadjnt.eq.2) b_mag=-b_mag ! T.Sato 2024/01/07, direction of magnetic field is reversed in the charged particle adjoint mode
                  ldsm  = ldsm + 1
                  s_mag = das_kmags(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  p_mag = das_kmags(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  t_mag = das_kmags(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  ldsm  = ldsm + 1
                  u_mag = das_kmags(kdsm+ldsm)

               do j = 1, ntrn

                  call tregck(iblz1,ilev1,ilat1,
     &                        mtrn,idas_ingrt(idsm+jdsm+1),jj,icc)

                  if( icc .ne. 0 ) goto 20

               end do

                  jdsm = jdsm + mtrn

            end do

               goto 25

   20       continue

*-----------------------------------------------------------------------

               isp = nint( s_mag )

            if( isp .ne.  2  .and. isp .ne.   4 .and.
     &          isp .ne. 60  .and. isp .ne.  61 .and.
     &          isp .ne. 62  .and.
     &          isp .ne. 100 .and. isp .ne. 101 .and.
     &          isp .ne. 102 .and. isp .ne. 103 .and.
     &          isp .ne. 104 .and. isp .ne. 106 .and.
     &          isp .ne. -1 .and. isp .ne. -2 .and.
     &          isp .ne. -3 .and. isp .ne. -4 .and.
     &          isp .ne. -101 .and. isp .ne. -102 .and.
     &          isp .ne. -103 .and. isp .ne. -104) goto 25
            if( ( isp .eq. 2 .or. isp .eq. 4 .or.
     &            isp .eq. -1 .or. isp .eq. -2 .or.
     &            isp .eq. -3 .or. isp .eq. -4) .and.
     &            jtyp .eq. 0 ) goto 25
            if( ( isp .eq.  60 .or. isp .eq.  61 .or.
     &            isp .eq.  62 .or.
     &            isp .eq. 100 .or. isp .eq. 101 .or.
     &            isp .eq. 102 .or. isp .eq. 103 .or.
     &            isp .eq. 104 .or. isp .eq. 106 .or.
     &            isp .eq. -101 .or. isp .eq. -102 .or.
     &            isp .eq. -103 .or. isp .eq. -104) .and.
     &            ityp .ne. 2 ) goto 25

*-----------------------------------------------------------------------
*           Wobbler magnet for charged particles for IHI
*-----------------------------------------------------------------------

            if( p_mag .gt. -1000.0d0 ) then

               if( jtyp .ne. 0 ) then

                  b_mag = sin( pmphs + p_mag ) * b_mag

               end if

            end if

*-----------------------------------------------------------------------

               call magfld(mark,markp,fpl,nbeta,itmak,
     &                     a_mag,b_mag,s_mag,p_mag,t_mag,u_mag)
               goto 100

         end if

   25    continue

*-----------------------------------------------------------------------
*        charge particle transfer under electro magnetic field
*           nereg    : number of region of electro magnetic field
*-----------------------------------------------------------------------

         if( nelctf .gt. 0 .and. nereg .gt. 0 .and.
     &       jtyp .ne. 0 ) then

                  idsm = inert
                  jdsm = 0

                  kdsm = kelcs
                  ldsm = 0

            do i = 1, nereg

                  jj  = 0

                  jdsm = jdsm + 1
                  ntrn = idas_inert(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_inert(idsm+jdsm)

                  ldsm  = ldsm + 1
                  s_elf = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  s_mgf = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  t_elf = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 2
                  t_mgf = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 2
                  emap_type = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  mmap_type = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  a_elmg = das_kelcs(kdsm+ldsm)

               do j = 1, ntrn

                  call tregck(iblz1,ilev1,ilat1,
     &                        mtrn,idas_inert(idsm+jdsm+1),jj,icc)

                  if( icc .ne. 0 ) goto 21

               end do

                  jdsm = jdsm + mtrn

            end do

               goto 26

   21       continue

*-----------------------------------------------------------------------

               !AdvanceSoft Hasemi 2019/11/27 2019/12/21
               call elmgfd(mark,markp,fpl,nbeta,itmak,
     &                     s_elf,s_mgf,t_elf,t_mgf,e_chs,
     &                     emap_type,mmap_type,a_elmg)

               goto 100

         end if

   26    continue

*-----------------------------------------------------------------------
*        proton beam transfer by fpl in medium
*        beam spread
*           only for nspred > 0, name(no) = 1, proton, no void
*           nfcs(no) =0,1 not 2 ; no spread for forced collisions
*-----------------------------------------------------------------------

             do i=1,mntsc
             if(ntsc(i).eq.iblz(ibkblz+no,ipomp+1))then
               mat1 = ntscell(idgr(ntsc(i)))
               goto 1122
             endif
             enddo
               mat1 = 0
 1122          continue

             if( ityp .eq. 1              .and.
     .           mat1 .eq. 1              .and.
     .         ( e(ibke+no,ipomp+1) .le. ptsmax ) .and.
     .         ( e(ibke+no,ipomp+1) .ge. ptsmin ) )then
               call ionts_trn(mark,markp,fpl,nbeta,itmak)
               return
             elseif( ktyp .eq. 6000012    .and.
     .           mat1 .eq. 1              .and.
     .         ( e(ibke+no,ipomp+1) .le. ctsmax ) .and.
     .         ( e(ibke+no,ipomp+1) .ge. ctsmin ) )then
               call ionts_trn(mark,markp,fpl,nbeta,itmak)
               return

             endif


         if( jtyp .ne. 0     .and.
     &       nspred .gt. 0   .and. mat .ne. 0 )then
          if((( nspsgn .eq. 1 .and. name(ibknam+no,ipomp+1) .eq. 1 ).or.
     &         nspsgn .eq. -1 ).and.
     &       nfcs(ibknfc+no,ipomp+1) .ne. 2 .and.
     &       arg(kcar+mat) .gt. 0.0     ) then

               call sprd(mark,markp,fpl,nbeta,itmak)

               goto 100

          endif
         end if

*-----------------------------------------------------------------------
*        particle transfer by fpl in free space or neutron gravity
*-----------------------------------------------------------------------

   30 continue

            if( ityp .eq. 2 .and. igrav .ne. 0 .and.
     &          e(ibke+no,ipomp+1) .le. 1.d-6 ) then

               call ngravt(mark,markp,fpl,nbeta,itmak)

            else

               call parfre(mark,markp,fpl,nbeta,itmak)

            end if

               goto 100

*-----------------------------------------------------------------------

  100 continue

      if(junf.ne.0)then !20220909frtati
      if(abs(lat(1,icl)).eq.3.and.iii.ne.
     &                      itetpos(ibtetpos+no,ipomp+1))then
         !FURUTA20190903
       itetpos(ibtetpos+no,ipomp+1)=iii  !FURUTA20160607
       iblz(ibkblz+no,ipomp+1)=idrg(icl) !FURUTA20171020
       call tetragetmat(matnum) !FURUTA20200612
       nmed(ibknmd+no,ipomp+1)=matnum   !FURUTA20200612
      endif
      endif !20220909frtati

*-----------------------------------------------------------------------
*        surface cross: reset the forced collision flag
*-----------------------------------------------------------------------

         if( mark .eq. 0 .or. mark .eq. 2 ) then

            nfcs(ibknfc+no,ipomp+1) = 0

         end if

*-----------------------------------------------------------------------
*        surface cross: keep importance of the previous non-void cell
*-----------------------------------------------------------------------

         if( ( mark .eq. 0 .or. mark .eq. 2 ) .and. mat .gt. 0 ) then

            wtnz(ibkwnz+no,ipomp+1) = aimp(ityp,iblz1,ilev1,ilat1,ii1)

         end if

*-----------------------------------------------------------------------
*        error in CG/GG ( mark = -2, icge = -2 ) : lost particles
*-----------------------------------------------------------------------

         if( mark .eq. -2 ) then

               icge = -1

               ilost = ilost + 1





        ErrCha = ''
        ErrID = 'L:756/R:partrs/F:partrs.f' !W06_002_001
        call ErrWrite(ErrID,ErrCha)

        write (*,'(''*** Lost particle ***'')')

        if( npe .gt. 0 ) then ! T.Sato 2019/01/07, for MPI
         write(*,'( '' my ip       = '',i3)') me
        endif

        write(*,'(''batch & history number:'',2i20)') nobch,nocas

        write(*,'(''initial & final cells :'',2i20)') iblz1, iblz2

        write(*,'( ''location (x,y,z)        :'',3es16.8)')
     &  ,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),z(ibkz+no,ipomp+1)

        write(*,'( ''direction vector (u,v,w):'',3es16.8)')
     &  ,u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),w(ibkw+no,ipomp+1)

c ---------------------------------------------------------------------

            if( ilost .ge. nlost ) then

             write(ErrCha,'(''Calculation is terminated!!!'')')
             ErrID = 'L:780/R:partrs/F:partrs.f' !E06_002_001
             call ErrWrite(ErrID,ErrCha)
             write(ErrCha,'(''number of'',
     &       '' lost particles exceeds nlost (='',i8,'')'')') nlost
             ErrID = 'L:784/R:partrs/F:partrs.f'
             call ErrWrite(ErrID,ErrCha)

             call parastop( 802 )

            end if

               return

         end if

*-----------------------------------------------------------------------
*        error in CG/GG ( mark = -3, -4, -5 ): region error
*        goto 81, start again after gomsor
*-----------------------------------------------------------------------

         if( mark .le. -3 ) then

                  icge  = icge + 1
                  icger = icger + 1

*-----------------------------------------------------------------------
*              booking the error and write information
*-----------------------------------------------------------------------

                  nb1 = iblz1
                  nb2 = iblz2

                  if( nb1 .lt. nb2 ) then

                     nb3 = nb1
                     nb1 = nb2
                     nb2 = nb3

                  end if
!$OMP CRITICAL (nrovp_crit)
                  do i = 1, novp

                     if( nrovp(1,i) .eq. nb1 .and.
     &                   nrovp(2,i) .eq. nb2 ) then

                        nrovp(3,i) = nrovp(3,i) + 1

                        goto 200

                     end if

                  end do

                  novp = novp + 1

                  if( novp .le. 1000 ) then

                     nrovp(1,novp) = nb1
                     nrovp(2,novp) = nb2
                     nrovp(3,novp) = 1

                  end if

  200             continue
!$OMP END CRITICAL (nrovp_crit)
*-----------------------------------------------------------------------
*           icge gt igerr
*-----------------------------------------------------------------------

            if( icge .gt. igerr ) then

                  if( icger .le. nlost ) then

                  write(6,'(/''*** error in region check *** no ='',
     &                  i6)') icger

                 if( npe .gt. 0 )
     &            write(6,'( '' my ip       = '',i3)') me

                  write(6,'( '' nbch ncs no = '',3i10)')
     &                  nobch, nocas, no
                  write(6,'( '' ityp, e(no) = '',i4,1x,e17.8)')
     &                  ityp, e(ibke+no,ipomp+1)
                  write(6,'( ''        mark ='',3i6)')
     &                  mark
                  write(6,'( '' reg. ini fin     ='',5i6)')
     &                  iblz1, iblz2
                  write(6,'( '' mat. ini fin     ='',5i6)')
     &                  idmn(mat), idmn(nmed(ibknmd+no,ipomp+1))
                  write(6,'( '' x, y, z :'',3e17.8)')
     &                  x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1)
                  write(6,'( '' xc,yc,zc:'',3e17.8)')
     &                  xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1)
                  write(6,'( '' u, v, w :'',3e17.8)')
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1)

                   write(6,'(''*** failed to recovering ***'')')
        write(6,'(''*** LOST particle! Geometry error occurrs ***'')')

                  endif

            if( ncger .ge. nlost ) then

               write(6,'(/''*** Number of unrecovered error exceeds ''
     &         ,''nlost : nlost ='',i8)') nlost
               call parastop( 802 )

            end if

                  ncger = ncger + 1


                  icge  = -2

                  return

*-----------------------------------------------------------------------
*           icge > 1
*-----------------------------------------------------------------------

            else if( icge .ge. 1 ) then

               if( mark .ne. -4 ) then

                  x(ibkx+no,ipomp+1) = x(ibkx+no,ipomp+1) +
     &                                  u(ibku+no,ipomp+1) * parz(28)
                  y(ibky+no,ipomp+1) = y(ibky+no,ipomp+1) +
     &                                  v(ibkv+no,ipomp+1) * parz(28)
                  z(ibkz+no,ipomp+1) = z(ibkz+no,ipomp+1) +
     &                                  w(ibkw+no,ipomp+1) * parz(28)

               end if

                  ici   = 0
                  mark  = 1
                  markp = 0

                  call gomsor(x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                        z(ibkz+no,ipomp+1),
     &                        u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                        w(ibkw+no,ipomp+1),
     &                        nmed(ibknmd+no,ipomp+1),
     &                        iblz(ibkblz+no,ipomp+1),
     &                        mark,markp,ici)

                     if( mark .le. -2 ) then

                        icge  = -2
                        return

                     end if

                  xc(ibkxc+no,ipomp+1) = x(ibkx+no,ipomp+1)
                  yc(ibkyc+no,ipomp+1) = y(ibky+no,ipomp+1)
                  zc(ibkzc+no,ipomp+1) = z(ibkz+no,ipomp+1)
                  ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1)

                  return

            end if

*-----------------------------------------------------------------------
*        recovered errors
*-----------------------------------------------------------------------

         else if( icge .gt. 0 ) then

            if( icger .le. nrecover ) then

                  write(6,'(/''*** warning in region check *** no ='',
     &                  i6)') icger

                  if( npe .gt. 0 )
     &            write(6,'( '' my ip       = '',i3)') me

                  write(6,'( '' nbch ncs no = '',3i10)')
     &                  nobch, nocas, no
                  write(6,'( '' ityp, e(no) = '',i4,1x,e17.8)')
     &                  ityp, e(ibke+no,ipomp+1)
                  write(6,'( ''        mark ='',3i6)')
     &                  mark
                  write(6,'( '' reg. ini fin     ='',5i6)')
     &                  iblz1, iblz2
                  write(6,'( '' mat. ini fin     ='',5i6)')
     &                  idmn(mat), idmn(nmed(ibknmd+no,ipomp+1))
                  write(6,'( '' x, y, z :'',3e17.8)')
     &                  x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1)
                  write(6,'( '' xc,yc,zc:'',3e17.8)')
     &                  xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1)
                  write(6,'( '' u, v, w :'',3e17.8)')
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1)

                  write(6,'(''*** succeeded in recovering *** icge ='',
     &                      i5)') icge

            end if

            icge = -3

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine parfre(mark,markp,dist,nbeta,itmak)
*                                                                      *
*                                                                      *
*       particle transport by d                                        *
*       and region check                                               *
*       modified by K.Niita on 2005/02/02                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distance                                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
************************************************************************

c--------------------------------
      use MMBANKMOD !FURUTA
      use moddas_mesh


      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /eparm/  esmax, esmin, emin(20)
      common /paraj/  mstz(300), parz(300)
      common /tcntl/  icntl, inucr
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx

*-----------------------------------------------------------------------
      common /celdg/  rhog(kvlmax)

      common /tdcyl/ tdcylife, itdcyw
!$OMP THREADPRIVATE(/tdcyl/)

      common /  regdc    / idrg(kvlmax), idgr(kvmmax)
      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax
      common / tsminmax  / tsmax

*-----------------------------------------------------------------------

      common /wwin00/ iwwxyz
!$OMP THREADPRIVATE(/wwin00/)
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww
      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------
*        dv = fpl : flight length
*-----------------------------------------------------------------------

                  iff = 0

                  dv  = dist

                  esav = e(ibke+no,ipomp+1)
                  tsav = t(ibkt+no,ipomp+1)
                  xsav = x(ibkx+no,ipomp+1)
                  ysav = y(ibky+no,ipomp+1)
                  zsav = z(ibkz+no,ipomp+1)

*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

                  emint = emin(ityp)

*-----------------------------------------------------------------------
               if( mntsc .gt. 0 .and.
     &             ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.
     &           ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &           ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &             e(ibke+no,ipomp+1) .ge. etsmin ) )  emint = etsmin

*-----------------------------------------------------------------------

            else

                  emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------
*        delt1 : maximum flight step
*        delt  : flight step
*            parz(27) : deltm
*            parz(163): deltc for nedisp
*-----------------------------------------------------------------------
*        dens = rhog(mat): density (g/cm3)
*        Correction factor for deltc and deltm 2015/08/13, Y. Iwamoto
*-----------------------------------------------------------------------

                  dens = 1.0

               if( mstz(111) .eq. 1 .and. mat .gt. 0 )
     &            dens = rhog(mat)

                  delt1 = min( parz(27)/dens, dist )

               if( mstz(111) .eq. 1 .and. mat .eq. 0 )
     &            delt1 = 1.d100                        ! 2020/3/27  Ogawa  max flight path is 1.d100 in vacuum

               if( jtyp .ne. 0 .and. mat .gt. 0 .and.
     &             (nedisp .ne. 0 .or. nspred .ne. 0) )
     &            delt1 = min( parz(163)/max(1.0,dens), delt1 ) ! T.Sato 2023/08/10 max value should be deltc

*-----------------------------------------------------------------------

                  delt  = delt1

*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
                  iwwxyz = 0
                  jwwxyz = 0
                  kwwxyz = 0

*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt, dv=dv-delt
*        nbeta = 2 : dv < delt,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

            call gomdis(dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                      z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  mark,markp,nmed(ibknmd+no,ipomp+1),
     &                  iblz(ibkblz+no,ipomp+1))

               if( mark .le. -2 ) goto 600

*-----------------------------------------------------------------------
*        WW of xyz mesh
*-----------------------------------------------------------------------

         if( iwwdp .gt. 0 .and. iwmsh .eq. 3 ) then

               dwwt = delt
               if( dwwt .gt. dpr ) dwwt = dpr
               iwctl = 1

            call wwxdis(iwctl,dpx,dwwt,
     &                  x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  iwxnm(1),iwynm(1),iwznm(1),
     &                  das_iwxrg(iwxrg(1)),das_iwyrg(iwyrg(1)),
     &                  das_iwzrg(iwzrg(1)))

            if( dpx .le. dpr ) then

               if( dpx .eq. dpr ) kwwxyz = 1

               dpr = dpx
               jwwxyz = 1

            end if

         end if

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

                  ee1 = e(ibke+no,ipomp+1)
                  ecc = ee1
                  delt2 = delt

            if( jtyp .ne. 0 .and. mat .ne. 0 .and.
     &          icntl .ne. 5 .and. icntl .ne. 14 ) then

                  call rainge(ee1,rng1,mat,ityp,ktyp,jtyp,rtyp)
                  call ecol(ecc,delt2,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

            end if

                  ec(ibkec+no,ipomp+1) = ecc

*-----------------------------------------------------------------------
*        next step or collision or stop
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

               if( delt2 .lt. delt ) then

                  nbeta = 3
                  delt = delt2

               end if

               if( ecc .le. emint )  nbeta = 3

*-----------------------------------------------------------------------

               if( nbeta .eq. 1 ) then

                  goto 800

               else

                  goto 700

               end if

*-----------------------------------------------------------------------
*        cross surface or stop
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

            if( jtyp .ne. 0 .and. mat .ne. 0 .and.
     &          icntl .ne. 5 .and. icntl .ne. 14 ) then

                  dpr1 = dpr

                  call ecol(ecc,dpr1,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

                  ec(ibkec+no,ipomp+1) = ecc

               if( ecc .le. emint .or. dpr1 .lt. dpr ) then

                  call rainge(emint,rngs,mat,ityp,ktyp,jtyp,rtyp)

                  nbeta = 3
                  delt  = rng1 - rngs
                  ec(ibkec+no,ipomp+1) = emint

                  goto 700

               end if

            end if

*-----------------------------------------------------------------------
*           WW of xyz mesh
*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 0 ) then

                  nbeta = 2
                  delt = dpr
                  iwwxyz = 1

                  goto 700

            end if

*-----------------------------------------------------------------------
*           search new cell
*-----------------------------------------------------------------------

               call gomnew(1,dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                     z(ibkz+no,ipomp+1),
     &                     xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                     zc(ibkzc+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     ec(ibkec+no,ipomp+1),
     &                     nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1),
     &                     mark,markp)

*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 1 ) then

                  iwwxyz = 2

            end if

*-----------------------------------------------------------------------

               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------
*              call tally if itstep = 1,  ncol = 15
*-----------------------------------------------------------------------

               if( itstep .ne. 0 .and. nedisp .ne. 0 ) then

                  call timtrs(itmak,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no,ipomp+1)
                  tsav = tc(ibktc+no,ipomp+1)
                  xsav = xc(ibkxc+no,ipomp+1)
                  ysav = yc(ibkyc+no,ipomp+1)
                  zsav = zc(ibkzc+no,ipomp+1)

                  ncol = 15

                  call analyz(ncol,mark)

               end if

*-----------------------------------------------------------------------

                  e(ibke+no,ipomp+1) = ec(ibkec+no,ipomp+1)
                  t(ibkt+no,ipomp+1) = tc(ibktc+no,ipomp+1)
                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------

  600 continue

                  e(ibke+no,ipomp+1) = esav
                  t(ibkt+no,ipomp+1) = tsav
                  x(ibkx+no,ipomp+1) = xsav
                  y(ibky+no,ipomp+1) = ysav
                  z(ibkz+no,ipomp+1) = zsav

      if( nbeta .eq. 3 .and. idcyc(nkf(ibknkf+no,ipomp+1)) .eq. 1 )
     &                                      itdcyw = 1

                  call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sprd(mark,markp,dist,nbeta,itmak)
*                                                                      *
*                                                                      *
*       proton beam spread by Coulomb                                  *
*       and region check                                               *
*       modified by K.Niita on 2010/12/23                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distance                                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use moddas_mesh

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /eparm/  esmax, esmin, emin(20)
      common /paraj/  mstz(300), parz(300)
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx

      common /tdcyl/ tdcylife, itdcyw
!$OMP THREADPRIVATE(/tdcyl/)

      common /celdg/  rhog(kvlmax)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax
      common / tsminmax  / tsmax

*-----------------------------------------------------------------------

      common /wwin00/ iwwxyz
!$OMP THREADPRIVATE(/wwin00/)
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww
      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------
*        dv = dist : flight length
*-----------------------------------------------------------------------

               iff = 0

               dv  = dist

               esav = e(ibke+no,ipomp+1)
               tsav = t(ibkt+no,ipomp+1)
               xsav = x(ibkx+no,ipomp+1)
               ysav = y(ibky+no,ipomp+1)
               zsav = z(ibkz+no,ipomp+1)

*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

                  emint = emin(ityp)

*-----------------------------------------------------------------------
               if( mntsc .gt. 0 .and.
     &             ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.
     &           ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &           ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &             e(ibke+no,ipomp+1) .ge. etsmin ) )  emint = etsmin

*-----------------------------------------------------------------------

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------
*        delt0 : minimum flight step
*        delt1 : maximum flight step
*-----------------------------------------------------------------------

                  delt0 = parz(26)

*-----------------------------------------------------------------------
* Correction factor for deltc and deltm 2015/08/13, Y. Iwamoto

            if( mstz(111) .eq. 0 ) then
                  delt1 = min( parz(27), dist )
               if( jtyp .ne. 0 .and. mat .gt. 0 .and.
     &             (nedisp .ne. 0 .or. nspred .ne. 0) )
     &             delt1 = min( parz(163), delt1 )
            end if
*-----------------------------------------------------------------------
*        delt1 : maximum flight step
*        delt  : flight step
*            parz(27) : deltm
*            parz(163): deltc for nedisp
*            dist: 1e10
*        dens = rhog(mat): density (g/cm3)
*-----------------------------------------------------------------------

            if( mstz(111) .eq. 1 ) then

                ddeltm = 0.0
                ddeltc = 0.0
               if( mat .ne .0 ) dens = rhog(mat) ! S.Abe 2024/11/27 avoid error when mat=0

               if( mat .eq. 0 ) dens = 1.0

                ddeltm = parz(27) / dens

                   delt1 = min( ddeltm, dist )

               if( jtyp .ne. 0 .and. mat .gt. 0 .and.
     &             nedisp .ne. 0 ) then

                    ddeltc = parz(163) / max(1.0d0,dens) ! T.Sato 2023/08/10 max value should be deltc

                    delt1 = min( ddeltc, delt1 )

               end if

            end if

*-----------------------------------------------------------------------

                  delt  = delt1

*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
                  iwwxyz = 0
                  jwwxyz = 0
                  kwwxyz = 0

*-----------------------------------------------------------------------
                  ee1 = e(ibke+no,ipomp+1)
                  call rainge(ee1,rngm,mat,ityp,ktyp,jtyp,rtyp)

                  delt = min(delt,1.2*rngm)

*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt,     dv = dv - delt
*        nbeta = 2 : delt = dv,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if

*-----------------------------------------------------------------------

               if( delt .le. delt0 ) then

                  idelt = 1

               else

                  idelt = 0

               end if

*-----------------------------------------------------------------------
*        transport particle by delt through sprdtrs(delt)
*-----------------------------------------------------------------------

                  uprv = u(ibku+no,ipomp+1)
                  vprv = v(ibkv+no,ipomp+1)
                  wprv = w(ibkw+no,ipomp+1)

               call sprdtrs(delt)

               call gomupp(mark,markp,
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1))

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

               call gomdis(dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                     z(ibkz+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     mark,markp,nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1))

                  if( mark .le. -2 ) goto 600

*-----------------------------------------------------------------------
*        back to the initial by reducing the delta
*-----------------------------------------------------------------------

         if( delt .ge. dpr .and. idelt .eq. 0 ) then

                  if( nbeta .eq. 1 ) dv = dv + delt

                  u(ibku+no,ipomp+1) = uprv
                  v(ibkv+no,ipomp+1) = vprv
                  w(ibkw+no,ipomp+1) = wprv

                  call gomupp(mark,markp,
     &                        u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                        w(ibkw+no,ipomp+1))

                  delt = max( delt0, dpr * 0.9 )

                  nmed(ibknmd+no,ipomp+1) = mat
                  iblz(ibkblz+no,ipomp+1) = iblz1

                  goto 500

         end if

*-----------------------------------------------------------------------
*        WW of xyz mesh
*-----------------------------------------------------------------------

         if( iwwdp .gt. 0 .and. iwmsh .eq. 3 ) then

               dwwt = delt
               if( dwwt .gt. dpr ) dwwt = dpr
               iwctl = 1

            call wwxdis(iwctl,dpx,dwwt,
     &                  x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  iwxnm(1),iwynm(1),iwznm(1),
     &                  das_iwxrg(iwxrg(1)),das_iwyrg(iwyrg(1)),
     &                  das_iwzrg(iwzrg(1)))

            if( dpx .le. dpr ) then

               if( dpx .eq. dpr ) kwwxyz = 1

               dpr = dpx
               jwwxyz = 1

            end if

         end if

*-----------------------------------------------------------------------
*        energy range and final energy
*-----------------------------------------------------------------------

                  ee1 = e(ibke+no,ipomp+1)
                  delt2 = delt

                  call rainge(ee1,rng1,mat,ityp,ktyp,jtyp,rtyp)
                  call ecol(ecc,delt2,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

*-----------------------------------------------------------------------
*        next step or collision or stop
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

               if( delt2 .lt. delt ) then

                  nbeta = 3
                  delt = delt2

               end if

               if( ecc .le. emint )  nbeta = 3

                  ec(ibkec+no,ipomp+1) = ecc

*-----------------------------------------------------------------------

               if( nbeta .eq. 1 ) then

                  goto 800

               else

                  goto 700

               end if

*-----------------------------------------------------------------------
*        cross surface or stop
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

                  dpr1 = dpr

                  call ecol(ecc,dpr1,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

                  ec(ibkec+no,ipomp+1) = ecc

               if( ecc .le. emint .or. dpr1 .lt. dpr ) then

                  call rainge(emint,rngs,mat,ityp,ktyp,jtyp,rtyp)

                  nbeta = 3
                  delt  = rng1 - rngs
                  ec(ibkec+no,ipomp+1) = emint

                  goto 700

               end if

*-----------------------------------------------------------------------
*           WW of xyz mesh
*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 0 ) then

                  nbeta = 2
                  delt = dpr
                  iwwxyz = 1

                  goto 700

            end if

*-----------------------------------------------------------------------
*           search new cell
*-----------------------------------------------------------------------

               call gomnew(1,dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                     z(ibkz+no,ipomp+1),
     &                     xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                     zc(ibkzc+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     ec(ibkec+no,ipomp+1),
     &                     nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1),
     &                     mark,markp)

*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 1 ) then

                  iwwxyz = 2

            end if

*-----------------------------------------------------------------------

               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------
*              call tally if itstep = 1, ncol = 15
*-----------------------------------------------------------------------

               if( itstep .ne. 0 ) then

                  call timtrs(itmak,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no,ipomp+1)
                  tsav = tc(ibktc+no,ipomp+1)
                  xsav = xc(ibkxc+no,ipomp+1)
                  ysav = yc(ibkyc+no,ipomp+1)
                  zsav = zc(ibkzc+no,ipomp+1)

                  ncol = 15

                  call analyz(ncol,mark)

               end if

*-----------------------------------------------------------------------

                  e(ibke+no,ipomp+1) = ec(ibkec+no,ipomp+1)
                  t(ibkt+no,ipomp+1) = tc(ibktc+no,ipomp+1)
                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------

  600 continue

                  e(ibke+no,ipomp+1) = esav
                  t(ibkt+no,ipomp+1) = tsav
                  x(ibkx+no,ipomp+1) = xsav
                  y(ibky+no,ipomp+1) = ysav
                  z(ibkz+no,ipomp+1) = zsav

      if( nbeta .eq. 3 .and. idcyc(nkf(ibknkf+no,ipomp+1)) .eq. 1 )
     &                                           itdcyw = 1

                  call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ngravt(mark,markp,dist,nbeta,itmak)
*                                                                      *
*                                                                      *
*       neutron with gravity                                           *
*       and region check                                               *
*       modified by K.Niita on 2010/12/23                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distance                                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use moddas_mesh

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( grvc = 980.665d-18 )
      parameter ( rlit = 29.97925d0 )

*-----------------------------------------------------------------------


      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /eparm/  esmax, esmin, emin(20)
      common /paraj/  mstz(300), parz(300)
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      common /gravit/ grav(3), igrav
*-----------------------------------------------------------------------
      common /celdg/  rhog(kvlmax)

*-----------------------------------------------------------------------

      common /wwin00/ iwwxyz
!$OMP THREADPRIVATE(/wwin00/)
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww
      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------
*        dv = dist : flight length
*-----------------------------------------------------------------------

               iff = 0

               dv  = dist

               esav = e(ibke+no,ipomp+1)
               tsav = t(ibkt+no,ipomp+1)
               xsav = x(ibkx+no,ipomp+1)
               ysav = y(ibky+no,ipomp+1)
               zsav = z(ibkz+no,ipomp+1)

*-----------------------------------------------------------------------
*        delt0 : minimum flight step
*        delt1 : maximum flight step
*-----------------------------------------------------------------------

                  delt0 = parz(26)


*-----------------------------------------------------------------------
* Correction factor for deltc and deltm 2015/08/13, Y. Iwamoto

            if( mstz(111) .eq. 0 ) then
                  delt1 = min( parz(27), dist )
                  delt1 = min( parz(163), delt1 )
            end if
*-----------------------------------------------------------------------
*        delt1 : maximum flight step
*        delt  : flight step
*            parz(27) : deltm
*            parz(163): deltc for nedisp
*            dist: 1e10
*        dens = rhog(mat): density (g/cm3)
*-----------------------------------------------------------------------

            if( mstz(111) .eq. 1 ) then

                ddeltm = 0.0
                ddeltc = 0.0
                dens = rhog(mat)

               if( mat .eq. 0 ) dens = 1.0

                ddeltm = parz(27) / dens

                   delt1 = min( ddeltm, dist )

                    ddeltc = parz(163) / max(1.0d0,dens) ! T.Sato 2023/08/10 max value should be deltc
                    delt1 = min( ddeltc, delt1 )

            end if
*-----------------------------------------------------------------------

                  delt  = delt1

*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
                  iwwxyz = 0
                  jwwxyz = 0
                  kwwxyz = 0

*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt,     dv = dv - delt
*        nbeta = 2 : delt = dv,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if

*-----------------------------------------------------------------------

               if( delt .le. delt0 ) then

                  idelt = 1

               else

                  idelt = 0

               end if

*-----------------------------------------------------------------------
*        transport particle by delt
*-----------------------------------------------------------------------

                  delts = delt

                  uprv = u(ibku+no,ipomp+1)
                  vprv = v(ibkv+no,ipomp+1)
                  wprv = w(ibkw+no,ipomp+1)

*-----------------------------------------------------------------------
*           gravity
*-----------------------------------------------------------------------

               xc(ibkxc+no,ipomp+1) = x(ibkx+no,ipomp+1) +
     &                               delt * u(ibku+no,ipomp+1)
               yc(ibkyc+no,ipomp+1) = y(ibky+no,ipomp+1) +
     &                               delt * v(ibkv+no,ipomp+1)
               zc(ibkzc+no,ipomp+1) = z(ibkz+no,ipomp+1) +
     &                               delt * w(ibkw+no,ipomp+1)

               vel = sqrt( 2.0 * e(ibke+no,ipomp+1) / rtyp ) * rlit
               timd = delt / vel

               vxx = u(ibku+no,ipomp+1) * vel - grvc * grav(1) * timd
               vyy = v(ibkv+no,ipomp+1) * vel - grvc * grav(2) * timd
               vzz = w(ibkw+no,ipomp+1) * vel - grvc * grav(3) * timd

               vab = sqrt( vxx**2 + vyy**2 + vzz**2 )

               ec(ibkec+no,ipomp+1) = 0.5 * rtyp * ( vab / rlit )**2

               xc(ibkxc+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
     &                      - 0.5 * grvc * grav(1) * timd**2
               yc(ibkyc+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
     &                      - 0.5 * grvc * grav(2) * timd**2
               zc(ibkzc+no,ipomp+1) = zc(ibkzc+no,ipomp+1)
     &                      - 0.5 * grvc * grav(3) * timd**2

               u(ibku+no,ipomp+1)   = vxx / vab
               v(ibkv+no,ipomp+1)   = vyy / vab
               w(ibkw+no,ipomp+1)   = vzz / vab

               delt = sqrt( ( xc(ibkxc+no,ipomp+1) -
     &                                  x(ibkx+no,ipomp+1) )**2
     &                    + ( yc(ibkyc+no,ipomp+1) -
     &                                  y(ibky+no,ipomp+1) )**2
     &                    + ( zc(ibkzc+no,ipomp+1) -
     &                                  z(ibkz+no,ipomp+1) )**2 )

*-----------------------------------------------------------------------

               call gomupp(mark,markp,
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1))

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

               call gomdis(dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                     z(ibkz+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     mark,markp,nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1))

                  if( mark .le. -2 ) goto 600

*-----------------------------------------------------------------------
*        back to the initial by reducing the delta
*-----------------------------------------------------------------------

         if( delt .ge. dpr .and. idelt .eq. 0 ) then

                  if( nbeta .eq. 1 ) dv = dv + delts

                  u(ibku+no,ipomp+1) = uprv
                  v(ibkv+no,ipomp+1) = vprv
                  w(ibkw+no,ipomp+1) = wprv

               call gomupp(mark,markp,
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1))

                  delt = max( delt0, dpr * 0.9 )

                  nmed(ibknmd+no,ipomp+1) = mat
                  iblz(ibkblz+no,ipomp+1) = iblz1

                  goto 500

         end if

*-----------------------------------------------------------------------
*        WW of xyz mesh
*-----------------------------------------------------------------------

         if( iwwdp .gt. 0 .and. iwmsh .eq. 3 ) then

               dwwt = delt
               if( dwwt .gt. dpr ) dwwt = dpr
               iwctl = 1

            call wwxdis(iwctl,dpx,dwwt,
     &                  x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  iwxnm(1),iwynm(1),iwznm(1),
     &                  das_iwxrg(iwxrg(1)),das_iwyrg(iwyrg(1)),
     &                  das_iwzrg(iwzrg(1)))

            if( dpx .le. dpr ) then

               if( dpx .eq. dpr ) kwwxyz = 1

               dpr = dpx
               jwwxyz = 1

            end if

         end if

*-----------------------------------------------------------------------
*        next step or collisions
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

            if( nbeta .eq. 1 ) then

               goto 800

            else

               goto 700

            end if

*-----------------------------------------------------------------------
*           WW of xyz mesh
*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 0 ) then

                  nbeta = 2
                  delt = dpr
                  iwwxyz = 1

                  goto 700

            end if

*-----------------------------------------------------------------------
*        cross surface and search new cell
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

               call gomnew(1,dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                     z(ibkz+no,ipomp+1),
     &                     xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                     zc(ibkzc+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     ec(ibkec+no,ipomp+1),
     &                     nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1),
     &                     mark,markp)

*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 1 ) then

                  iwwxyz = 2

            end if

*-----------------------------------------------------------------------

               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------
*              call tally if itstep = 1, ncol = 15
*-----------------------------------------------------------------------

               if( itstep .ne. 0 ) then

                  call timtrs(itmak,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no,ipomp+1)
                  tsav = tc(ibktc+no,ipomp+1)
                  xsav = xc(ibkxc+no,ipomp+1)
                  ysav = yc(ibkyc+no,ipomp+1)
                  zsav = zc(ibkzc+no,ipomp+1)

                  ncol = 15

                  call analyz(ncol,mark)

               end if

*-----------------------------------------------------------------------

                  e(ibke+no,ipomp+1) = ec(ibkec+no,ipomp+1)
                  t(ibkt+no,ipomp+1) = tc(ibktc+no,ipomp+1)
                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------

  600 continue

                  e(ibke+no,ipomp+1) = esav
                  t(ibkt+no,ipomp+1) = tsav
                  x(ibkx+no,ipomp+1) = xsav
                  y(ibky+no,ipomp+1) = ysav
                  z(ibkz+no,ipomp+1) = zsav

                  call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine magfld(mark,markp,dist,nbeta,itmak,
     &                  a_mag,b_mag,s_mag,p_mag,t_mag,u_mag)
*                                                                      *
*                                                                      *
*       charge particle or neutron transfer under magnetic field       *
*       and region check                                               *
*       modified by K.Niita on 2005/03/30                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distace                                               *
*       a_mag  : magnet gap(mm)                                        *
*       b_mag  : magnet field at pole tip  [kG]                        *
*       s_mag  : speicies of magnet dypole:2, quad:4, sext:6, oct:8    *
*       p_mag  : phase of magnetic field for charge particle           *
*                polarization of neuteron                              *
*       t_mag  : transform id                                          *
*       u_mag  : critical time of time dependent magnetic field        *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA

      use MEMBANKMOD !FURUTA
      use moddas_mesh

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      parameter ( rlit = 29.97925d0 )

      common /inout/  in,io
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /paraj/  mstz(300), parz(300)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /eparm/  esmax, esmin, emin(20)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      common /magreg/ nmreg, ingrc, ingrt, kmags, magtin, imgusr

      common /tdcyl/ tdcylife, itdcyw
!$OMP THREADPRIVATE(/tdcyl/)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax
      common / tsminmax  / tsmax

      common /spsgn/  nspsgn

      common /argcns/ kcar, kczs, kcze

*-----------------------------------------------------------------------

      common /wwin00/ iwwxyz
!$OMP THREADPRIVATE(/wwin00/)
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww
      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------
      common /celdg/  rhog(kvlmax)

*-----------------------------------------------------------------------
*        dv = dist : distance
*-----------------------------------------------------------------------

               iff = 0

               dv  = dist

               esav = e(ibke+no,ipomp+1)
               tsav = t(ibkt+no,ipomp+1)
               xsav = x(ibkx+no,ipomp+1)
               ysav = y(ibky+no,ipomp+1)
               zsav = z(ibkz+no,ipomp+1)

*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

                  emint = emin(ityp)

*-----------------------------------------------------------------------
               if( mntsc .gt. 0 .and.
     &             ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.
     &           ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &           ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &             e(ibke+no,ipomp+1) .ge. etsmin ) )  emint = etsmin

*-----------------------------------------------------------------------

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------
*        delt0 : minimum flight step
*        delt1 : maximum flight step ( parz(109) for deltg )
*                                    ( parz(27)  for normal )
*                                    ( parz(163) for deltc )
*                                    ( parz(171) for deltt msec)
*-----------------------------------------------------------------------

                  delt0 = parz(26)

*-----------------------------------------------------------------------
* Correction factor for deltc and deltm 2015/08/13, Y. Iwamoto
*        mstz(111) = idelt
*        dens = rhog(mat): density (g/cm3)

               if( mstz(111) .eq. 0 ) then
                  delt1 = min( parz(27), dist )
                  if( jtyp .ne. 0 .and. mat .gt. 0 .and.
     &             (nedisp .ne. 0 .or. nspred .ne. 0) )
     &                delt1 = min( parz(163), delt1 )
               else
                  ddeltm = 0.d0
                  ddeltc = 0.d0
                  if( mat .eq. 0 ) then ! T.Sato 2022/12/16, avoid access violation
                   dens = 1.d0
                  else
                   dens = rhog(mat)
                  endif

                  ddeltm = parz(27) / dens
                  delt1 = min( ddeltm, dist )

                  if( jtyp .ne. 0 .and. mat .gt. 0 .and.
     &             (nedisp .ne. 0 .or. nspred .ne. 0) ) then
                     ddeltc = parz(163) / max(1.0d0,dens) ! T.Sato 2023/08/10 max value should be deltc
                     delt1 = min( ddeltc, delt1 )
                  endif
               endif

               delt1 = min( parz(109), delt1 )

*-----------------------------------------------------------------------

               if( u_mag .gt. -1.0d+9 .and. esav .lt. 1.e-6 ) then

                     timd = delt1 * sqrt( rtyp / 2.0 / esav )
     &                    / rlit * 1.e-6

                  if( timd .gt. parz(171) ) then

                     delt1 = parz(171) / sqrt( rtyp / 2.0 / esav )
     &                     * rlit / 1.e-6

                  end if

               end if


                  delt  = delt1

*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
                  iwwxyz = 0
                  jwwxyz = 0
                  kwwxyz = 0

*-----------------------------------------------------------------------
         if( jtyp .ne. 0     .and.
     &       nspred .gt. 0   .and. mat .ne. 0 )then
          if(( ( nspsgn .eq. 1 .and.
     &            name(ibknam+no,ipomp+1) .eq. 1 ) .or.
     &         nspsgn .eq. -1 ).and.
     &       nfcs(ibknfc+no,ipomp+1) .ne. 2 .and.
     &       arg(kcar+mat) .gt. 0.0     ) then

                  ee1 = e(ibke+no,ipomp+1)
                  call rainge(ee1,rngm,mat,ityp,ktyp,jtyp,rtyp)

                  delt = min(delt,rngm)

          endif
         end if

*-----------------------------------------------------------------------
*           Time depedent magnet
*-----------------------------------------------------------------------

            if( u_mag .gt. -1.0d+9 ) then

                  ptime = abs(t(ibkt+no,ipomp+1)) * 1.e-6

               if( imgusr .eq. 1 ) then

                  call usrmgt1(b_mag,u_mag,ptime,o_mag)

               else if( imgusr .eq. 2 ) then

                  call usrmgt2(b_mag,u_mag,ptime,o_mag)

               end if

            else

                  o_mag = b_mag

            end if

*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt,     dv = dv - delt
*        nbeta = 2 : delt = dv,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if

*-----------------------------------------------------------------------

               if( delt .le. delt0 ) then

                  idelt = 1

               else

                  idelt = 0

               end if

*-----------------------------------------------------------------------
*        transport particle by delt(distance) under magnetic field
*-----------------------------------------------------------------------

                  delts = delt

                  uprv = u(ibku+no,ipomp+1)
                  vprv = v(ibkv+no,ipomp+1)
                  wprv = w(ibkw+no,ipomp+1)

                  sxsv = spx(ibkspx+no,ipomp+1)
                  sysv = spy(ibkspy+no,ipomp+1)
                  szsv = spz(ibkspz+no,ipomp+1)

               call magtrs(a_mag,o_mag,s_mag,p_mag,t_mag,
     &                     delt,dpr,udir,vdir,wdir,mark)
               if(mark.eq.-2) return ! T.Sato 2022/03/22

*-----------------------------------------------------------------------
         if( jtyp .ne. 0     .and.
     &       nspred .gt. 0   .and. mat .ne. 0 )then
          if(( ( nspsgn .eq. 1 .and.
     &          name(ibknam+no,ipomp+1) .eq. 1 ) .or.
     &         nspsgn .eq. -1 ).and.
     &       nfcs(ibknfc+no,ipomp+1) .ne. 2 .and.
     &       arg(kcar+mat) .gt. 0.0     ) then

                  xprv = x(ibkx+no,ipomp+1)
                  yprv = y(ibky+no,ipomp+1)
                  zprv = z(ibkz+no,ipomp+1)

                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1) - dpr * udir
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1) - dpr * vdir
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1) - dpr * wdir

                  u(ibku+no,ipomp+1) = udir
                  v(ibkv+no,ipomp+1) = vdir
                  w(ibkw+no,ipomp+1) = wdir

               call sprdtrs(dpr)

                  x(ibkx+no,ipomp+1) = xprv
                  y(ibky+no,ipomp+1) = yprv
                  z(ibkz+no,ipomp+1) = zprv

                  udir = u(ibku+no,ipomp+1)
                  vdir = v(ibkv+no,ipomp+1)
                  wdir = w(ibkw+no,ipomp+1)

                  u(ibku+no,ipomp+1) = uprv
                  v(ibkv+no,ipomp+1) = vprv
                  w(ibkw+no,ipomp+1) = wprv

          endif
         end if

*-----------------------------------------------------------------------

               call gomupp(mark,markp,
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1))

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

               call gomdis(dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                     z(ibkz+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     mark,markp,nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1))

                  if( mark .le. -2 ) goto 600

*-----------------------------------------------------------------------
*        back to the initial by reducing the delta
*-----------------------------------------------------------------------

         if( delt .ge. dpr .and. idelt .eq. 0 ) then

                  if( nbeta .eq. 1 ) dv = dv + delts

                  u(ibku+no,ipomp+1) = uprv
                  v(ibkv+no,ipomp+1) = vprv
                  w(ibkw+no,ipomp+1) = wprv

                  spx(ibkspx+no,ipomp+1) = sxsv
                  spy(ibkspy+no,ipomp+1) = sysv
                  spz(ibkspz+no,ipomp+1) = szsv

               call gomupp(mark,markp,
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1))

                  delt = max( delt0, dpr * 0.9 )

                  nmed(ibknmd+no,ipomp+1) = mat
                  iblz(ibkblz+no,ipomp+1) = iblz1

                  goto 500

         end if

*-----------------------------------------------------------------------
*        WW of xyz mesh
*-----------------------------------------------------------------------

         if( iwwdp .gt. 0 .and. iwmsh .eq. 3 ) then

               dwwt = delt
               if( dwwt .gt. dpr ) dwwt = dpr
               iwctl = 1

            call wwxdis(iwctl,dpx,dwwt,
     &                  x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  iwxnm(1),iwynm(1),iwznm(1),
     &                  das_iwxrg(iwxrg(1)),das_iwyrg(iwyrg(1)),
     &                  das_iwzrg(iwzrg(1)))

            if( dpx .le. dpr ) then

               if( dpx .eq. dpr ) kwwxyz = 1

               dpr = dpx
               jwwxyz = 1

            end if

         end if

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

                  ee1 = e(ibke+no,ipomp+1)
                  ecc = ee1
                  delt2 = delt

            if( jtyp .ne. 0 .and. mat .ne. 0 ) then

                  call rainge(ee1,rng1,mat,ityp,ktyp,jtyp,rtyp)
                  call ecol(ecc,delt2,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

            end if

                  ec(ibkec+no,ipomp+1) = ecc

*-----------------------------------------------------------------------
*        next step or collision or stop
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

               if( delt2 .lt. delt ) then

                  nbeta = 3
                  delt = delt2

               end if

               if( ecc .le. emint )  nbeta = 3

*-----------------------------------------------------------------------
*              for dipole, delt is sometimes changed ( small value )
*-----------------------------------------------------------------------

               if( nbeta .eq. 1 ) then

                  dv = dv + ( delts - delt )

                  goto 800

               else if( nbeta .eq. 2 ) then

                  dv = dv - delt

                  if( dv .le. 0.0d0 ) goto 700

                  nbeta = 1

                  goto 800

               else

                  goto 700

               end if

*-----------------------------------------------------------------------
*        cross surface or stop
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

            if( jtyp .ne. 0 .and. mat .ne. 0 ) then

                  dpr1 = dpr

                  call ecol(ecc,dpr1,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

                  ec(ibkec+no,ipomp+1) = ecc

               if( ecc .le. emint .or. dpr1 .lt. dpr ) then

                  call rainge(emint,rngs,mat,ityp,ktyp,jtyp,rtyp)

                  nbeta = 3
                  delt  = rng1 - rngs
                  ec(ibkec+no,ipomp+1) = emint

                  goto 700

               end if

            end if

*-----------------------------------------------------------------------
*           WW of xyz mesh
*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 0 ) then

                  nbeta = 2
                  delt = dpr
                  iwwxyz = 1

                  goto 700

            end if

*-----------------------------------------------------------------------
*           search new cell
*-----------------------------------------------------------------------

               call gomnew(1,dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                     z(ibkz+no,ipomp+1),
     &                     xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                     zc(ibkzc+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     ec(ibkec+no,ipomp+1),
     &                     nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1),
     &                     mark,markp)

*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 1 ) then

                  iwwxyz = 2

            end if

*-----------------------------------------------------------------------

               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

                  u(ibku+no,ipomp+1) = udir
                  v(ibkv+no,ipomp+1) = vdir
                  w(ibkw+no,ipomp+1) = wdir

               call gomupp(mark,markp,
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1))

*-----------------------------------------------------------------------
*              call tally if itstep = 1, ncol = 15
*-----------------------------------------------------------------------

               if( itstep .ne. 0 ) then

                  call timtrs(itmak,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no,ipomp+1)
                  tsav = tc(ibktc+no,ipomp+1)
                  xsav = xc(ibkxc+no,ipomp+1)
                  ysav = yc(ibkyc+no,ipomp+1)
                  zsav = zc(ibkzc+no,ipomp+1)

                  ncol = 15

                  call analyz(ncol,mark)

               end if

*-----------------------------------------------------------------------

                  e(ibke+no,ipomp+1) = ec(ibkec+no,ipomp+1)
                  t(ibkt+no,ipomp+1) = tc(ibktc+no,ipomp+1)
                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------

  600 continue

                  e(ibke+no,ipomp+1) = esav
                  t(ibkt+no,ipomp+1) = tsav
                  x(ibkx+no,ipomp+1) = xsav
                  y(ibky+no,ipomp+1) = ysav
                  z(ibkz+no,ipomp+1) = zsav

      if( nbeta .eq. 3 .and. idcyc(nkf(ibknkf+no,ipomp+1)) .eq. 1 )
     &               itdcyw = 1

                  call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine elmgfd(mark,markp,dist,nbeta,itmak,
     &                  s_elf,s_mgf,t_elf,t_mgf,e_chs,
     &                  emap_type, mmap_type, a_elmg)
*                                                                      *
*                                                                      *
*       charge particle under electro magnetic field                   *
*       and region check                                               *
*       modified by K.Niita on 2011/01/10                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distace                                               *
*       s_elf  : electric field (kV/cm)                                *
*       s_mgf  : dipole magnet field  [kG]                             *
*       t_elf  : transform for electirc field                          *
*       t_mgf  : transform for magnetic field                          *
*       e_chs  : charge state for this region                          *
*       emap_type : type of electric data map                          *
*                   -1:xyz-list                                        *
*                   -2:rz-list                                         *
*                   -3:xyz-map                                         *
*                   -4:rz-map                                          *
*       mmap_type : type of magnetic data map                          *
*                   -1:xyz-list                                        *
*                   -2:rz-list                                         *
*                   -3:xyz-map                                         *
*                   -4:rz-map                                          *
*       a_elmg  : gap                                                  *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
* AdvanceSoft Hasemi 2019/12/21 modified for electric map              *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA

      use MEMBANKMOD !FURUTA
      use moddas_mesh

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      parameter ( rlit = 29.97925d0 )

      common /inout/  in,io
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /paraj/  mstz(300), parz(300)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /eparm/  esmax, esmin, emin(20)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx

      common /elmgrg/ nelctf, nereg, inerc, inert, kelcs

      common /tdcyl/ tdcylife, itdcyw
!$OMP THREADPRIVATE(/tdcyl/)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax
      common / tsminmax  / tsmax

      common /spsgn/  nspsgn

*-----------------------------------------------------------------------

      common /wwin00/ iwwxyz
!$OMP THREADPRIVATE(/wwin00/)
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww
      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------
      common /argcns/ kcar, kczs, kcze

*-----------------------------------------------------------------------
      common /celdg/  rhog(kvlmax)

*-----------------------------------------------------------------------
*        transform ( electric field = x,  magnetic field = y )
*-----------------------------------------------------------------------
               itre = nint( t_elf )
               call trnsuv(1.d0,0.d0,0.d0,elcx,elcy,elcz,itre)

               itrm = nint( t_mgf )
               call trnsuv(0.d0,1.d0,0.d0,bmgx,bmgy,bmgz,itrm)

*-----------------------------------------------------------------------
*        charge state
*-----------------------------------------------------------------------

               chgp = ctyp

*-----------------------------------------------------------------------
*        constant
*-----------------------------------------------------------------------

               elct = chgp * s_elf * 1.d-3 * rlit
               bmgt = chgp * s_mgf / 3.3356 * rlit

*-----------------------------------------------------------------------
*        dv = dist : distance
*-----------------------------------------------------------------------

               iff = 0

               dv  = dist

               esav = e(ibke+no,ipomp+1)
               tsav = t(ibkt+no,ipomp+1)
               xsav = x(ibkx+no,ipomp+1)
               ysav = y(ibky+no,ipomp+1)
               zsav = z(ibkz+no,ipomp+1)

*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

                  emint = emin(ityp)

*-----------------------------------------------------------------------
               if( mntsc .gt. 0 .and.
     &             ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.
     &           ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &           ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &             e(ibke+no,ipomp+1) .ge. etsmin ) )  emint = etsmin

*-----------------------------------------------------------------------

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------
*        delt0 : minimum flight step
*        delt1 : maximum flight step ( parz(109) for deltg )
*                                    ( parz(27)  for normal )
*                                    ( parz(163) for deltc )
*-----------------------------------------------------------------------

                  delt0 = parz(26)

*-----------------------------------------------------------------------
* Correction factor for deltc and deltm 2015/08/13, Y. Iwamoto
*        mstz(111) = idelt
*        dens = rhog(mat): density (g/cm3)

               if( mstz(111) .eq. 0 ) then
                  delt1 = min( parz(27), dist )
                  if( jtyp .ne. 0 .and. mat .gt. 0 .and.
     &             (nedisp .ne. 0 .or. nspred .ne. 0) )
     &                delt1 = min( parz(163), delt1 )
               else
                  ddeltm = 0.d0
                  ddeltc = 0.d0
                  if( mat .eq. 0 ) then
                    dens = 1.d0
                  else
                    dens = rhog(mat)
                  end if

                  ddeltm = parz(27) / dens
                  delt1 = min( ddeltm, dist )

                  if( jtyp .ne. 0 .and. mat .gt. 0 .and.
     &             (nedisp .ne. 0 .or. nspred .ne. 0) ) then
                     ddeltc = parz(163) / max(1.0d0,dens) ! T.Sato 2023/08/10 max value should be deltc
                     delt1 = min( ddeltc, delt1 )
                  endif
               endif

               delt1 = min( parz(109), delt1 )

                  delt  = delt1

*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
                  iwwxyz = 0
                  jwwxyz = 0
                  kwwxyz = 0

*-----------------------------------------------------------------------
         if( jtyp .ne. 0     .and.
     &       nspred .gt. 0   .and. mat .ne. 0 )then
          if(( ( nspsgn .eq. 1 .and.
     &          name(ibknam+no,ipomp+1) .eq. 1 ) .or.
     &         nspsgn .eq. -1 ).and.
     &       nfcs(ibknfc+no,ipomp+1) .ne. 2 .and.
     &       arg(kcar+mat) .gt. 0.0     ) then

                  ee1 = e(ibke+no,ipomp+1)
                  call rainge(ee1,rngm,mat,ityp,ktyp,jtyp,rtyp)

                  delt = min(delt,rngm)

          endif
         end if

*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt,     dv = dv - delt
*        nbeta = 2 : delt = dv,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if

*-----------------------------------------------------------------------

               if( delt .le. delt0 ) then

                  idelt = 1

               else

                  idelt = 0

               end if

*-----------------------------------------------------------------------
*        transport particle by delt(distance) under electro magnetic field
*-----------------------------------------------------------------------

                  delts = delt

                  uprv = u(ibku+no,ipomp+1)
                  vprv = v(ibkv+no,ipomp+1)
                  wprv = w(ibkw+no,ipomp+1)

*-----------------------------------------------------------------------

               call elmtrs(elct,elcx,elcy,elcz,
     &                     bmgt,bmgx,bmgy,bmgz,chgp,
     &                     delt,dpr,udir,vdir,wdir,
     &                     emap_type,mmap_type,a_elmg,
     &                     t_elf,t_mgf)

*-----------------------------------------------------------------------
         if( jtyp .ne. 0     .and.
     &       nspred .gt. 0   .and. mat .ne. 0 )then
          if(( ( nspsgn .eq. 1 .and.
     &         name(ibknam+no,ipomp+1) .eq. 1 ) .or.
     &         nspsgn .eq. -1 ).and.
     &       nfcs(ibknfc+no,ipomp+1) .ne. 2 .and.
     &       arg(kcar+mat) .gt. 0.0     ) then

                  xprv = x(ibkx+no,ipomp+1)
                  yprv = y(ibky+no,ipomp+1)
                  zprv = z(ibkz+no,ipomp+1)

                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1) - dpr * udir
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1) - dpr * vdir
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1) - dpr * wdir

                  u(ibku+no,ipomp+1) = udir
                  v(ibkv+no,ipomp+1) = vdir
                  w(ibkw+no,ipomp+1) = wdir

               call sprdtrs(dpr)

                  x(ibkx+no,ipomp+1) = xprv
                  y(ibky+no,ipomp+1) = yprv
                  z(ibkz+no,ipomp+1) = zprv

                  udir = u(ibku+no,ipomp+1)
                  vdir = v(ibkv+no,ipomp+1)
                  wdir = w(ibkw+no,ipomp+1)

                  u(ibku+no,ipomp+1) = uprv
                  v(ibkv+no,ipomp+1) = vprv
                  w(ibkw+no,ipomp+1) = wprv

          endif
         end if

*-----------------------------------------------------------------------

               call gomupp(mark,markp,
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1))

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

               call gomdis(dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                     z(ibkz+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     mark,markp,nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1))

                  if( mark .le. -2 ) goto 600

*-----------------------------------------------------------------------
*        back to the initial by reducing the delta
*-----------------------------------------------------------------------

         if( delt .ge. dpr .and. idelt .eq. 0 ) then

                  if( nbeta .eq. 1 ) dv = dv + delts

                  u(ibku+no,ipomp+1) = uprv
                  v(ibkv+no,ipomp+1) = vprv
                  w(ibkw+no,ipomp+1) = wprv

               call gomupp(mark,markp,
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1))

                  delt = max( delt0, dpr * 0.9 )

                  nmed(ibknmd+no,ipomp+1) = mat
                  iblz(ibkblz+no,ipomp+1) = iblz1

                  goto 500

         end if

*-----------------------------------------------------------------------
*        WW of xyz mesh
*-----------------------------------------------------------------------

         if( iwwdp .gt. 0 .and. iwmsh .eq. 3 ) then

               dwwt = delt
               if( dwwt .gt. dpr ) dwwt = dpr
               iwctl = 1

            call wwxdis(iwctl,dpx,dwwt,
     &                  x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  iwxnm(1),iwynm(1),iwznm(1),
     &                  das_iwxrg(iwxrg(1)),das_iwyrg(iwyrg(1)),
     &                  das_iwzrg(iwzrg(1)))

            if( dpx .le. dpr ) then

               if( dpx .eq. dpr ) kwwxyz = 1

               dpr = dpx
               jwwxyz = 1

            end if

         end if

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

                  ee1 = e(ibke+no,ipomp+1)
                  ecc = ee1
                  delt2 = delt

            if( jtyp .ne. 0 .and. mat .ne. 0 ) then

                  call rainge(ee1,rng1,mat,ityp,ktyp,jtyp,rtyp)
                  call ecol(ecc,delt2,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

            end if

                  ec(ibkec+no,ipomp+1) = ecc

*-----------------------------------------------------------------------
*        next step or collision or stop
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

               if( delt2 .lt. delt ) then

                  nbeta = 3
                  delt = delt2

               end if

               if( ecc .le. emint )  nbeta = 3

*-----------------------------------------------------------------------
*              for dipole, delt is sometimes changed ( small value )
*-----------------------------------------------------------------------

               if( nbeta .eq. 1 ) then

                  dv = dv + ( delts - delt )

                  goto 800

               else if( nbeta .eq. 2 ) then

                  dv = dv - delt

                  if( dv .le. 0.0d0 ) goto 700

                  nbeta = 1

                  goto 800

               else

                  goto 700

               end if

*-----------------------------------------------------------------------
*        cross surface or stop
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

            if( jtyp .ne. 0 .and. mat .ne. 0 ) then

                  dpr1 = dpr

                  call ecol(ecc,dpr1,ee1,rng1,
     &                      mat,ityp,ktyp,jtyp,rtyp)

                  ec(ibkec+no,ipomp+1) = ecc

               if( ecc .le. emint .or. dpr1 .lt. dpr ) then

                  call rainge(emint,rngs,mat,ityp,ktyp,jtyp,rtyp)

                  nbeta = 3
                  delt  = rng1 - rngs
                  ec(ibkec+no,ipomp+1) = emint

                  goto 700

               end if

            end if

*-----------------------------------------------------------------------
*           WW of xyz mesh
*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 0 ) then

                  nbeta = 2
                  delt = dpr
                  iwwxyz = 1

                  goto 700

            end if

*-----------------------------------------------------------------------
*           search new cell
*-----------------------------------------------------------------------

               call gomnew(1,dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                     z(ibkz+no,ipomp+1),
     &                     xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                     zc(ibkzc+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     ec(ibkec+no,ipomp+1),
     &                     nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1),
     &                     mark,markp)

*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 1 ) then

                  iwwxyz = 2

            end if

*-----------------------------------------------------------------------

               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

                  u(ibku+no,ipomp+1) = udir
                  v(ibkv+no,ipomp+1) = vdir
                  w(ibkw+no,ipomp+1) = wdir

               call gomupp(mark,markp,
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1))

*-----------------------------------------------------------------------
*              call tally if itstep = 1, ncol = 15
*-----------------------------------------------------------------------

               if( itstep .ne. 0 ) then

                  call timtrs(itmak,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no,ipomp+1)
                  tsav = tc(ibktc+no,ipomp+1)
                  xsav = xc(ibkxc+no,ipomp+1)
                  ysav = yc(ibkyc+no,ipomp+1)
                  zsav = zc(ibkzc+no,ipomp+1)

                  ncol = 15

                  call analyz(ncol,mark)

               end if

*-----------------------------------------------------------------------

                  e(ibke+no,ipomp+1) = ec(ibkec+no,ipomp+1)
                  t(ibkt+no,ipomp+1) = tc(ibktc+no,ipomp+1)
                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------

  600 continue

                  e(ibke+no,ipomp+1) = esav
                  t(ibkt+no,ipomp+1) = tsav
                  x(ibkx+no,ipomp+1) = xsav
                  y(ibky+no,ipomp+1) = ysav
                  z(ibkz+no,ipomp+1) = zsav

      if( nbeta .eq. 3 .and. idcyc(nkf(ibknkf+no,ipomp+1)) .eq. 1 )
     &                                       itdcyw = 1

                  call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine timtrs(itmak,mark)
*                                                                      *
*       time transport of the particle                                 *
*       and check the time limit                                       *
*       modified by K.Niita on 2004/01/07                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       itmak  : =0, less than tmax                                    *
*                =1, greater than tmax                                 *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( rlit = 29.97925d0 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /tparm/  tmax(20)

      common /tdcyl/ tdcylife, itdcyw
!$OMP THREADPRIVATE(/tdcyl/)

*-----------------------------------------------------------------------
*     t(ibkt+no) < 0 : timer stopped the clock.
*-----------------------------------------------------------------------

               itmak = 0

           if( t(ibkt+no,ipomp+1) .lt. 0.d0 ) return

*-----------------------------------------------------------------------

               ekin = ( ec(ibkec+no,ipomp+1) + e(ibke+no,ipomp+1) )/ 2.0
               dist = sqrt( ( xc(ibkxc+no,ipomp+1)
     &                    - x(ibkx+no,ipomp+1) )**2
     &                    + ( yc(ibkyc+no,ipomp+1)
     &                    - y(ibky+no,ipomp+1) )**2
     &                    + ( zc(ibkzc+no,ipomp+1)
     &                    - z(ibkz+no,ipomp+1) )**2 )

                  timd = 0.0

               if( ekin .gt. 0.1 .or. rtyp .eq. 0.0d0 ) then

                  timd = dist * ( ekin + rtyp )
     &                 / sqrt( ekin * ( ekin + 2.0 * rtyp ) )
     &                 / rlit

               else if( ekin .gt. 0.0 ) then

                  timd = dist * sqrt( rtyp / 2.0 / ekin ) / rlit

               end if

                  tc(ibktc+no,ipomp+1) = t(ibkt+no,ipomp+1) + timd

*-----------------------------------------------------------------------
*           time over
*-----------------------------------------------------------------------

            if( tc(ibktc+no,ipomp+1) .gt. tmax(ityp) ) then

                  timd = ( tc(ibktc+no,ipomp+1) - tmax(ityp) )

                  dd = dist

               if( ekin .gt. 0.1 .or. rtyp .eq. 0.0d0 ) then

                  dd = dist - timd / ( ekin + rtyp )
     &                      * sqrt( ekin * ( ekin + 2.0 * rtyp ) )
     &                      * rlit

               else if( ekin .gt. 0.0 ) then

                  dd = dist - timd / sqrt( rtyp / 2.0 / ekin ) * rlit

               end if

                  xc(ibkxc+no,ipomp+1) = x(ibkx+no,ipomp+1) +
     &                                   u(ibku+no,ipomp+1) * dd
                  yc(ibkyc+no,ipomp+1) = y(ibky+no,ipomp+1) +
     &                                   v(ibkv+no,ipomp+1) * dd
                  zc(ibkzc+no,ipomp+1) = z(ibkz+no,ipomp+1) +
     &                                   w(ibkw+no,ipomp+1) * dd

                  itmak = 1
                  mark  = 1

            end if

*-----------------------------------------------------------------------
      if( itdcyw .eq. 1 .and. tdcylife .gt. 0.d0 ) then

       if( tdcylife .gt. tmax(ityp) ) then

        tc(ibktc+no,ipomp+1) = tmax(ityp)
        itmak = 1
        mark = 1

       else

        tc(ibktc+no,ipomp+1) = tdcylife

       endif

      endif

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine egstr(mark,markp,dist,nbeta,itmak)
*                                                                      *
*                                                                      *
*       particle transport by d                                        *
*       and region check                                               *
*       modified by K.Niita on 2014/08/10                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distance                                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use moddas_mesh

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /eparm/  esmax, esmin, emin(20)
      common /paraj/  mstz(300), parz(300)
      common /tcntl/  icntl, inucr
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx

*-----------------------------------------------------------------------
      common /ndemax/ dnmax(20)

      include 'include/egs5_h.f'
!     only destep is needed
      common/STACK/
     * egs5e(MXSTACK),      ! Total energy of particle (including rest
     * egs5x(MXSTACK),                                          ! X-po
     * egs5y(MXSTACK),                                          ! Y-po
     * egs5z(MXSTACK),                                          ! Z-po
     * egs5u(MXSTACK),                             ! X-axis direction
     * egs5v(MXSTACK),                             ! Y-axis direction
     * egs5w(MXSTACK),                             ! Z-axis direction
     * egs5uf(MXSTACK),         ! Electric field vectors of polarized
     * egs5vf(MXSTACK),
     * egs5wf(MXSTACK),
     * dnear(MXSTACK),          ! Estimated distance to nearest bounda
     * egs5wt(MXSTACK),                                    ! Particle
     * k1step(MXSTACK),      ! Scat stren to next hinge
     * k1rsd(MXSTACK),       ! Scat stren from hinge to end of step
     * k1init(MXSTACK),      ! Scat of prev hinge end of step
     * time(MXSTACK),
     * deinitial,
     * deresid,
     * denstep,
     * iq(MXSTACK),        ! Particle charge, -1(e-), 0(photons), +1(e
     * ir(MXSTACK),                                      ! Region numb
     * latch(MXSTACK),                               ! Latching variab
     * np,                                         ! Stack pointer ind
     * latchi                        ! Initialization for latch variab

c>ada.2014.12> only destep is needed ..?
!$OMP THREADPRIVATE(/STACK/)
      include 'include/egs5_misc.f'
      real*8 egs5x,egs5y,egs5z,egs5u,egs5v,egs5w,egs5uf,egs5vf,egs5wf,
     $       dnear,egs5wt
      real*8 egs5e,time
      real*8 k1step,k1rsd,k1init
      real*8 deinitial, deresid, denstep
      integer iq,ir,latch,np,latchi
      real*8 k1i,k1r,k1s

      common /egs5cmn3/thard,tinel,tmscat,hardstep,sig,scpow,dedx,sig0,
     $                 ams
!$OMP THREADPRIVATE(/egs5cmn3/)
      common /egs5cmn4/k1i,k1r,k1s
!$OMP THREADPRIVATE(/egs5cmn4/)
      real*8 pstep,dpmfp,gmfp
      common /egs5cmn6/pstep,dpmfp,gmfp
!$OMP THREADPRIVATE(/egs5cmn6/)

      include 'include/egs5_epcont.f'

      include 'include/egs5_useful.f'

      common /egs5cmn10/denstepold,denstepnew,deinit
      real*8 denstepold,denstepnew,deinit
!$OMP THREADPRIVATE(/egs5cmn10/)

*-----------------------------------------------------------------------
      common /celdg/  rhog(kvlmax)

      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax
      common / tsminmax  / tsmax

*-----------------------------------------------------------------------

      common /wwin00/ iwwxyz
!$OMP THREADPRIVATE(/wwin00/)
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww
      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)
      dimension     idas(1)
      equivalence ( das, idas )
      save infcheck,dprold,zcold ! T.Sat 2025/03/19 avoid compile error using gfortran OpenMP
      data infcheck/0/    ! infinite loop check, T.Sato 2023/08/11
      data dprold/1.0e20/ ! infinite loop check, T.Sato 2023/08/11
      data zcold/1.0e20/ ! infinite loop check, T.Sato 2023/08/11
!$OMP THREADPRIVATE(infcheck,dprold,zcold)
*-----------------------------------------------------------------------
*        dv = fpl : flight length
*-----------------------------------------------------------------------

                  iff = 0

                  dv  = dist

                  esav = e(ibke+no,ipomp+1)
                  tsav = t(ibkt+no,ipomp+1)
                  xsav = x(ibkx+no,ipomp+1)
                  ysav = y(ibky+no,ipomp+1)
                  zsav = z(ibkz+no,ipomp+1)

*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

                  emint = emin(ityp)

*-----------------------------------------------------------------------
               if( mntsc .gt. 0 .and.
     &             ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.
     &           ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &           ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &             e(ibke+no,ipomp+1) .ge. etsmin ) )  emint = etsmin

*-----------------------------------------------------------------------

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------
*        delt1 : maximum flight step
*        delt  : flight step
*-----------------------------------------------------------------------
* Correction factor for deltc and deltm 2015/08/13, Y. Iwamoto

            if( mstz(111) .eq. 0 ) then
                  delt1 = min( parz(27), dist )
            end if

*-----------------------------------------------------------------------
*        delt1 : maximum flight step
*        delt  : flight step
*            parz(27) : deltm
*            parz(163): deltc for nedisp
*            dist: 1e10
*        dens = rhog(mat): density (g/cm3)
*-----------------------------------------------------------------------

            if( mstz(111) .eq. 1 ) then

                ddeltm = 0.0
                ddeltc = 0.0
               if( mat .eq. 0 ) then
                 dens = 1.d-100
               else
                 dens = rhog(mat)
               end if

                ddeltm = parz(27) / dens

                   delt1 = min( ddeltm, dist )




            end if

*-----------------------------------------------------------------------

                  delt  = delt1
	if(ityp.eq.14) delt = 1.0e30 ! T.Sato 2016/03/25, avoid longer range of photon

*-----------------------------------------------------------------------

  500 continue

            iwwxyz = 0
            jwwxyz = 0
            kwwxyz = 0

*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt, dv=dv-delt
*        nbeta = 2 : dv < delt,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

            call gomdis(dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  mark,markp,nmed(ibknmd+no,ipomp+1),
     &                  iblz(ibkblz+no,ipomp+1))

               if( mark .le. -2 ) goto 600

*-----------------------------------------------------------------------
*        WW of xyz mesh
*-----------------------------------------------------------------------

         if( iwwdp .gt. 0 .and. iwmsh .eq. 3 ) then

               dwwt = delt
               if( dwwt .gt. dpr ) dwwt = dpr
               iwctl = 1

            call wwxdis(iwctl,dpx,dwwt,
     &                  x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  iwxnm(1),iwynm(1),iwznm(1),
     &                  das_iwxrg(iwxrg(1)),das_iwyrg(iwyrg(1)),
     &                  das_iwzrg(iwzrg(1)))

            if( dpx .le. dpr ) then

               if( dpx .eq. dpr ) kwwxyz = 1

               dpr = dpx
               jwwxyz = 1

            end if

         end if

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

               ecc = e(ibke+no,ipomp+1)
               delt2 = delt
               ec(ibkec+no,ipomp+1) = ecc

*---------------------------------------------------------------------
cc H.Iwase 2014/1/9  (should be called as a subroutine)
*        EGS5 electron transport
*---------------------------------------------------------------------

            if( jtyp .ne. 0 .and. mat .ne. 0 .and.
     &          icntl .ne. 5 .and. icntl .ne. 14 ) then

               tstep = min(delt,dpr)
               edep = 0d0
               if( mat .ne. 0 ) edep = tstep * dedx

               thard = thard  - tstep
               tinel = tinel  - tstep
               tmscat = tmscat - tstep

cc H.Iwase 2015/4/2 for keeping the old energy for [t-track]
               denstepold = denstep

               hardstep = thard  * sig
               k1s = tmscat * scpow
               denstep = tinel  * dedx

cc H.Iwase 2015/4/2 for keeping the old energy for [t-track]
               denstepnew = denstep
               deinit    = deinitial

               if(mstz(85).eq.99) then
                  write(93,*)'--------------------------------renewed'
                  write(93,'(a10,f15.9)')'(hardstep)',hardstep
                  write(93,'(a10,f15.9)')'(k1s)'     ,k1s
                  write(93,'(a10,f15.9)')'(denstep)' ,denstep
                  write(93,*)'--------------------------------'
               endif

           end if

*-----------------------------------------------------------------------
*        next step or collision or stop
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

            if( jtyp .ne. 0 .and. mat .ne. 0 .and.
     &          icntl .ne. 5 .and. icntl .ne. 14 ) then

               ein =  e(ibke+no,ipomp+1)
               eout = e(ibke+no,ipomp+1)

               irin = idgr(iblz(ibkblz+no,ipomp+1))
               call egs5ede(irin,ityp,ein,eout,nbeta)

               ec(ibkec+no,ipomp+1) = eout
               ecc = eout

            end if

            if( nbeta .eq. 1 ) then

               goto 800

            else

               goto 700

            end if

*-----------------------------------------------------------------------
*        cross surface or stop
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

*-----------------------------------------------------------------------
*           WW of xyz mesh
*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 0 ) then

                  nbeta = 2
                  delt = dpr
                  iwwxyz = 1

                  goto 700

            end if

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 1 ) iwwxyz = 2

*-----------------------------------------------------------------------
*           search new cell
*-----------------------------------------------------------------------

            irold = idgr(iblz(ibkblz+no,ipomp+1))

            call gomnew(1,dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  ec(ibkec+no,ipomp+1),
     &                  nmed(ibknmd+no,ipomp+1),iblz(ibkblz+no,ipomp+1),
     &                  mark,markp)

cc H.Iwase 2014/1/9  (photon transport: should be called as subroutine)
*---------------------------------------------------------------------
*           EGS5 photon transport
*---------------------------------------------------------------------

            if( ityp .eq. 14 ) then

        if(dprold.eq.dpr.and.dpr.lt.1.0e-3.and.
     &  zc(ibkzc+no,ipomp+1).eq.zcold) then ! infinite loop check, T.Sato 2023/08/11
         infcheck=infcheck+1
         if(infcheck.gt.1000) then
          mark=-2
          infcheck=0
          return
         endif
        else
         infcheck=0
         dprold=dpr
         zcold=zc(ibkzc+no,ipomp+1)
        endif

               pstep = pstep - dpr
               if( mat .ne. 0 .and. gmfp.ne.0) ! T.Sato 2022/09/08 to avoid NaN
     &         dpmfp = max(0.d0,dpmfp-dpr/gmfp)

               medold = medium
               irnew = idgr(iblz(ibkblz+no,ipomp+1))

               if (irnew .ne. irold) then ! Region has changed
                  irl = irnew
                  medium = med(irl)
               end if

               if( mark .eq. -1 ) pstep = 0d0

            endif

*-----------------------------------------------------------------------

               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------
*              call tally if itstep = 1,  ncol = 15
*-----------------------------------------------------------------------

               if( itstep .ne. 0 .and.
     &           ( nedisp .ne. 0 .or. nspred .ne. 0 ) ) then

                  call timtrs(itmak,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no,ipomp+1)
                  tsav = tc(ibktc+no,ipomp+1)
                  xsav = xc(ibkxc+no,ipomp+1)
                  ysav = yc(ibkyc+no,ipomp+1)
                  zsav = zc(ibkzc+no,ipomp+1)

                  ncol = 15

                  call analyz(ncol,mark)

               end if

*-----------------------------------------------------------------------

                  e(ibke+no,ipomp+1) = ec(ibkec+no,ipomp+1)
                  t(ibkt+no,ipomp+1) = tc(ibktc+no,ipomp+1)
                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------

  600 continue

                  e(ibke+no,ipomp+1) = esav
                  t(ibkt+no,ipomp+1) = tsav
                  x(ibkx+no,ipomp+1) = xsav
                  y(ibky+no,ipomp+1) = ysav
                  z(ibkz+no,ipomp+1) = zsav

                  call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine eletr(mark,markp,dist,nbeta,itmak)
*                                                                      *
*                                                                      *
*       electron transfer                                              *
*       and region check                                               *
*       modified by K.Niita on 2003/10/12                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distance                                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /eparm/  esmax, esmin, emin(20)

*-----------------------------------------------------------------------

      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, noz, mtel
!$OMP THREADPRIVATE(/elect/)


*-----------------------------------------------------------------------

      data delt0 /1.d-5/
      data delt1 /1.d-6/

*-----------------------------------------------------------------------

               nbeta = 2

*-----------------------------------------------------------------------
*        delt = dist : range of electron in ns substeop
*        delt0 : small distance before boundary
*-----------------------------------------------------------------------

               delt  = dist

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

               call gomdis(dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                     z(ibkz+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     mark,markp,nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1))

                  if( mark .le. -2 ) return

*-----------------------------------------------------------------------
*        cross the boundary
*-----------------------------------------------------------------------
*             mark: description
*              -1 : outgoing to the void region
*               0 : pass the forward surface
*               2 : reflect surface
*-----------------------------------------------------------------------

         if( delt .ge. 1.0d+10 .or.
     &     ( delt .ge. dpr .and. dpr .lt. delt0 ) ) then

               call gomnew(1,dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                     z(ibkz+no,ipomp+1),
     &                     xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                     zc(ibkzc+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     ec(ibkec+no,ipomp+1),
     &                     nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1),
     &                     mark,markp)

               ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1)

               ns = 0

               return

*-----------------------------------------------------------------------
*        collision at just before the boundary
*-----------------------------------------------------------------------

         else if( delt .ge. dpr .and. dpr .ge. delt1 ) then

               call gomupr(mark,markp,dpr-delt1)

               d = dpr - delt1

*-----------------------------------------------------------------------
*        collision at delt
*-----------------------------------------------------------------------

         else if( delt .lt. dpr ) then

               call gomupr(mark,markp,delt)

               d = delt

         end if

*-----------------------------------------------------------------------
*        change u and energy
*-----------------------------------------------------------------------

               uint(1) = u(ibku+no,ipomp+1)
               uint(2) = v(ibkv+no,ipomp+1)
               uint(3) = w(ibkw+no,ipomp+1)

               am = 1.
               call defelc(am,d,1.d0,1.d0,1,
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     uint(1),uint(2),uint(3))

               call gomupp(mark,markp,
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1))

               delc = d
               eint = e(ibke+no,ipomp+1)

               ec(ibkec+no,ipomp+1) =
     &                     max( 0.0d0, e(ibke+no,ipomp+1) - d * qs )
               noz = n1
               qo  = qs

*-----------------------------------------------------------------------
*        time evolution
*-----------------------------------------------------------------------

               call timtrs(itmak,mark)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sprdtrs(delt)
*                                                                      *
*       proton beam spread by Coulomb                                  *
*       modified by K.Niita on 2003/10/07                              *
*                                                                      *
*       transport particle by delt(distance)                           *
*                                                                      *
*         initial values                                               *
*                                                                      *
*              e(no), x(no), y(no), z(no), u(no), v(no), w(no)         *
*                                                                      *
*           final values                                               *
*                                                                      *
*              xc(no), yc(no), zc(no), u(no), v(no), w(no)             *
*                                                                      *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use MEMBANKMOD !FURUTA
      use moddas_material

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.1415926535898d0 )

*-----------------------------------------------------------------------

      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /argcns/ kcar, kczs, kcze

! T.Sato, 2015/03/10, introduce free parameter for Moliere 1st
      common /aspara/aspara1,aspara2

      common /kmat1g/ kmat(kvlmax)

*-----------------------------------------------------------------------

         if( delt .le. 0.0d0 ) then

               xc(ibkxc+no,ipomp+1) = x(ibkx+no,ipomp+1)
               yc(ibkyc+no,ipomp+1) = y(ibky+no,ipomp+1)
               zc(ibkzc+no,ipomp+1) = z(ibkz+no,ipomp+1)

               return

         end if

*-----------------------------------------------------------------------
*        nspred = 1 : NMTC original
*               = 2 : First order Moliere model ( by Meigo and Harada )
*   not comlete = 3 : Third order Moliere model ( by Meigo )
*-----------------------------------------------------------------------

         if( nspred .eq. 1 .or. nspred .eq. 2 .or.
     &     ( nspred .eq. 10 .and.
     &       ityp .ne. 1 .and. ityp .lt. 15 ) ) then

               sqxox = arg(kcar+mat) * dsqrt(delt)

            if( nspred .eq. 1 ) then

               sigth = 0.369 * abs(jtyp)
     &               / (e(ibke+no,ipomp+1)*
     &                 (1.+rtyp/(e(ibke+no,ipomp+1)+rtyp)))
     &               * sqxox / 0.866

               sigx = delt * sigth

            else if( nspred .eq. 2 .or. nspred .eq. 10 ) then

               sqxox1= sqxox / sqrt(431.273)
*                                   431.273 = 716.4 * 6.02e23 /1e24
! T.Sato, 2015/03/10, introduce free parameter, and fixed bug
               sigth = aspara1 * abs(jtyp)
     &               / (e(ibke+no,ipomp+1)*
     &               (1.+rtyp/(e(ibke+no,ipomp+1)+rtyp)))
     &               * sqxox1 * ( 1.0 + aspara2 * log( sqxox1**2 ) )

               sigx = delt * sigth

            end if

*-----------------------------------------------------------------------
               sigx = min( delt, sigx )
*-----------------------------------------------------------------------

  560          r1  = gaurn(x1)
               r2  = gaurn(x2)

               rm = sigx * sqrt( r1**2 + r2**2 )

               aa  = delt**2 - rm**2

               if( aa .lt. 0.0d0 ) goto 560

               thet = unirn(dummy) * 2.0 * pi
               cs   = cos(thet)
               sn   = sin(thet)

               xpr = rm * cs
               ypr = rm * sn

*-----------------------------------------------------------------------

         else if( nspred .eq. 3 ) then

               sqxox = arg(kcar+mat) * dsqrt(delt)

               p   = sqrt( e(ibke+no,ipomp+1)**2 +
     &                     2. * rtyp * e(ibke+no,ipomp+1) )
               bet = p / rtyp / sqrt( 1. + ( p / rtyp )**2 )

  570          continue

               sud = denh_das(kmat0+mat)
               sum = sud * 2.0
     &             * ( log( 1. + 3.34
     &             * ( 1.0 / 137.036 / bet )**2 ) )

               nel = nint( dnel_das(kmat0+mat) )

            do i = 1, nel

               zzik = zz_das(kmat(mat)+i)
                aik = a_das(kmat(mat)+i)
               dnik = den_das(kmat(mat)+i)

               sum = sum + dnik * zzik * ( zzik + 1.0 )
     &             * ( log( 1. + 3.34
     &             * ( zzik / 137.036 / bet )**2 ) )

               sud = sud + dnik

            end do

               zxmcs = sum / sud

               chiccmcs = 0.26057 * zsmcs(kczs+mat) * sud


               chick2 = chiccmcs / ( bet * p )**2 * delt

               chia2 = 2.016e-5 / p**2
     &               * exp( ( zxmcs - zemcs(kcze+mat) )
     &               / zsmcs(kczs+mat) )

               omega = chick2 / ( 1.167 * chia2 )

            if( omega .lt. 2.8 ) then

               call ruthscat(omega,chia2,xpr,ypr)

            else

               call moliere(omega,chick2,xpr,ypr)

            end if

               xpr = delt * xpr
               ypr = delt * ypr

               aa  = delt**2 - xpr**2 - ypr**2

            if( aa .lt. 0.0d0 ) goto 570

*-----------------------------------------------------------------------
*        angle straggling by ATIMA
*-----------------------------------------------------------------------

         else if( nspred .eq. 10 ) then

               ene = e(ibke+no,ipomp+1)

            if( ityp .ge. 15 ) then

               ap  = dble( ktyp - ktyp / 1000000 * 1000000 )
               zp  = dble( ktyp / 1000000 )

            else if( ityp .eq. 1 ) then

               ap  = 1.d0
               zp  = 1.d0

            else if( ityp .eq. 3 .or. ityp .eq. 5 .or.
     &               ityp .eq. 6 .or. ityp .eq. 7 .or.
     &               ityp .eq. 8 .or. ityp .eq. 10 .or.
     &               ityp .eq. 11 ) then

               ap = -1.d0 * dble(iabs(ktyp))
               zp = dble(iabs(jtyp))

            end if

               iway = 4

               call atima(ap,zp,ene,rtyp,mat,rng,delt,rm,iway)

               aa  = delt**2 - rm**2

               thet = unirn(dummy) * 2.0 * pi
               cs   = cos(thet)
               sn   = sin(thet)

               xpr = rm * cs
               ypr = rm * sn

         end if

*-----------------------------------------------------------------------

               zpr = dsqrt(aa)

               csth = w(ibkw+no,ipomp+1)
               snth = u(ibku+no,ipomp+1)**2 + v(ibkv+no,ipomp+1)**2

            if( snth .le. 0.0d0 ) then

               csphi = 1.0d0
               snphi = 0.0d0

            else

               snth  = dsqrt(snth)
               csphi = u(ibku+no,ipomp+1) / snth
               snphi = v(ibkv+no,ipomp+1) / snth

            end if

               cord = csth * xpr + snth * zpr

               xc(ibkxc+no,ipomp+1) = x(ibkx+no,ipomp+1) +
     &                                     csphi * cord - snphi * ypr
               yc(ibkyc+no,ipomp+1) = y(ibky+no,ipomp+1) +
     &                                     snphi * cord + csphi * ypr
               zc(ibkzc+no,ipomp+1) = z(ibkz+no,ipomp+1) -
     &                                     snth  * xpr  + csth  * zpr

               u(ibku+no,ipomp+1) =( csphi * cord - snphi * ypr ) / delt
               v(ibkv+no,ipomp+1) =( snphi * cord + csphi * ypr ) / delt
               w(ibkw+no,ipomp+1) =( -snth * xpr  + csth  * zpr ) / delt

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine sprdint
*                                                                      *
*                                                                      *
*       initialization of proton beam spread by Coulomb                *
*       modified by K.Niita on 10/02/2000                              *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       arg(kcar+i): out put in common                                 *
*                                                                      *
************************************************************************
      use MEMBANKMOD !FURUTA
      use moddas_material

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /spred/ nspred, nwsprd, nedisp, itstep, ndedx
      common /cparm/ maxbch,maxcas
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /argcns/ kcar, kczs, kcze

      common /kmat1g/ kmat(kvlmax)

*-----------------------------------------------------------------------

      if( nspred .eq. 0 ) return

*-----------------------------------------------------------------------
*     initialization of Coulomb spreading
*        nspred = 1 : NMTC original
*               = 2 : First order Moliere model ( by Meigo and Harada )
*   not comlete = 3 : Third order Moliere model ( by Meigo )
*-----------------------------------------------------------------------

      do 80 k = 1, mxmat

         if( nspred .eq. 1 ) then

            sumarg = denh_das(kmat0+k) * 21.132d+0 * 0.498d+0

         else if( nspred .eq. 2 .or. nspred .eq. 10 ) then

            sumarg = denh_das(kmat0+k) * 11.319d+0

         else if( nspred .eq. 3 ) then

            sumarg = denh_das(kmat0+k) * 11.319d+0
            sumzzs = denh_das(kmat0+k) * 2.
            sumzze= 0.0
            sumden= denh_das(kmat0+k)

         end if

            iii = nint( dnel_das(kmat0+k) )

         do 75 i = 1, iii

               zzik = zz_das(kmat(k)+i)
                aik = a_das(kmat(k)+i)
               dnik = den_das(kmat(k)+i)

            if( nspred .eq. 1 ) then

               sumarg = sumarg +
     &                ( dnik * zzik * (zzik+1.d+0)
     &                * ( 10.566d+0 - 0.333d+0 * log( zzik * aik ) )
     &                * 0.498d+0 )

            else if( nspred .eq. 2 .or. nspred .eq. 10 ) then

               sumarg = sumarg +
     &                  dnik * zzik * (zzik+1.d+0)
     &                * log( 287.d0 / sqrt( zzik ) )

            else if( nspred .eq. 3 ) then

               sumarg = sumarg +
     &                  dnik * zzik * (zzik+1.d+0)
     &                * log( 287.d0 / sqrt( zzik ) )

               sumzzs = sumzzs + dnik * zzik * (zzik+1.d+0)
               sumzze = sumzze + dnik * zzik * (zzik+1.d+0)
     &                * log(zzik)*(-2./3.)
               sumden = sumden + dnik

            end if

   75    continue

               arg(kcar+k) = sqrt(sumarg)

            if( nspred .eq. 3 ) then

               if( sumden .gt. 0.0 ) then

                  zsmcs(kczs+k) = sumzzs / sumden
                  zemcs(kcze+k) = sumzze / sumden

               else

                  zsmcs(kczs+k) = 0.0
                  zemcs(kcze+k) = 0.0

               end if

            end if

   80 continue

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine atten(x0,y0,z0,u0,v0,w0,
     &                 ipart,epart,distp,psigm,rsigm,ierr)
*                                                                      *
*                                                                      *
*       calculate the total mean free path = sgm * rho * distp         *
*       modified by K.Niita on 2015/03/18                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       x0,y0,z0 : initial position                                    *
*       u0,v0,w0 : unit vector                                         *
*       ipart    : type of particle                                    *
*       epart    : energy of particle                                  *
*       distp    : distance between r0 to r1                           *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       psigm  : total mean free path = sgm * rho * distp              *
*       rsigm  : atten. factor  = exp( - sgm * rho * distp )           *
*       ierr   : 0-> normal, 1-> error, 2-> reflect surface            *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /paraj/  mstz(300), parz(300)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

*-----------------------------------------------------------------------

      common /gg005/ xxx, yyy, zzz, uuu, vvv, www, tme, erg,
     &               dls, wgt, vel, dtc,
     &               icl, iii, jjj, kkk, jsu, iap, jgp, ipt,
     &               mtp, iexp, iex, idx,
     &               npa, ncp
!$OMP THREADPRIVATE(/gg005/)

*-----------------------------------------------------------------------
* for xstneu, xstgam

            iex0 = iex

*-----------------------------------------------------------------------
*        distance
*-----------------------------------------------------------------------

                  psigm = 0.d0
                  rsigm = 1.d0
                  ierr  = 0

*-----------------------------------------------------------------------

                  ici   = 0
                  mark  = 1
                  markp = 0

                  xx = x0
                  yy = y0
                  zz = z0

                  uu = u0
                  vv = v0
                  ww = w0

                  ein  = epart
                  delt = distp

*-----------------------------------------------------------------------

            call gomsor(xx,yy,zz,uu,vv,ww,
     &                  nmedi,iblzi,mark,markp,ici)

*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

            call gomdis(dpr,xx,yy,zz,uu,vv,ww,
     &                  mark,markp,nmedi,iblzi)

               if( mark .le. -2 ) goto 999

                  iclp = idgr(iblzi)
                  matp = nmedi
                  tmep = 0.d0

            if( matp .gt. 0 ) then
               rhop = denm(matp) ! T.Sato 2022/11/10, move to here

               if( ipart .eq. 2 ) then

                  call xstneu(0,sigtt,sigaa,iclp,ein,tmep,matp)

               else if( ipart .eq. 14 ) then

                  call xstgam(0,sigtt,sigaa,ein,matp)

               end if

            else

                  sigtt = 0.d0
                  rhop = 0.0

            end if

*-----------------------------------------------------------------------
*        final step
*-----------------------------------------------------------------------

         if( delt .le. dpr ) then

                psigm = psigm + sigtt * rhop * delt

             if( psigm .gt. -100 ) then
                rsigm = exp( - psigm )
             else
                rsigm = 0.0d0
             end if

                goto 888

*-----------------------------------------------------------------------
*        cross surface
*-----------------------------------------------------------------------

         else if( delt .gt. dpr ) then

                psigm = psigm + sigtt * rhop * dpr
                delt  = delt - dpr

*-----------------------------------------------------------------------
*           search new cell
*-----------------------------------------------------------------------

               call gomnew(0,dpr,xx,yy,zz,x2,y2,z2,uu,vv,ww,
     &                     epart,nmedi,iblzi,mark,markp)

                  if( mark .lt. 0 ) goto 999

                  xx = x2
                  yy = y2
                  zz = z2

               goto 500

         end if

*-----------------------------------------------------------------------

  999 continue
      ierr = 1

*-----------------------------------------------------------------------

  888 continue

                  ici   = 0
                  mark  = 1
                  markp = 0

                  iex   = iex0

      return
      end



*-----------------------------------------------------------------------

c  *********************************************************
      subroutine moliere(omega,chic2,dthxz,dthyz)
c
c     omega : mean number of scatters
c     chic2 : characteristic scattering angle (squared)
c.    *                                                                *
c.    * thred(na)=reduced angles of moliere theory                     *
c.    *                                                                *
c.    * f0i(na),f1i(na),f2i(na)= integrale of moliere functions        *
c.    *                                                                *
      implicit real*8 (a-h,o-z)

      parameter (eneper = 2.7182818d0) !,pi=3.141592654d0)
      dimension tint(40),arg(4),val(4),thred(40),f0i(40),f1i(40),f2i(40)
      data thred/
     +   0.00, 0.10, 0.20, 0.30
     +,  0.40, 0.50, 0.60, 0.70
     +,  0.80, 0.90, 1.00, 1.10
     +,  1.20, 1.30, 1.40, 1.50
     +,  1.60, 1.70, 1.80, 1.90
     +,  2.00, 2.20, 2.40, 2.60
     +,  2.80, 3.00, 3.20, 3.40
     +,  3.60, 3.80, 4.00, 5.00
     +,  6.00, 7.00, 8.00, 9.00
     +, 10.00,11.00,12.00,13.00/
      data f0i/
     +  0.000000e+00 ,0.995016e-02 ,0.392106e-01 ,0.860688e-01
     + ,0.147856e+00 ,0.221199e+00 ,0.302324e+00 ,0.387374e+00
     + ,0.472708e+00 ,0.555142e+00 ,0.632121e+00 ,0.701803e+00
     + ,0.763072e+00 ,0.815480e+00 ,0.859142e+00 ,0.894601e+00
     + ,0.922695e+00 ,0.944424e+00 ,0.960836e+00 ,0.972948e+00
     + ,0.981684e+00 ,0.992093e+00 ,0.996849e+00 ,0.998841e+00
     + ,0.999606e+00 ,0.999877e+00 ,0.999964e+00 ,0.999990e+00
     + ,0.999998e+00 ,0.999999e+00 ,0.100000e+01 ,0.100000e+01
     + ,0.100000e+01 ,0.100000e+01 ,0.100000e+01 ,0.100000e+01
     + ,1.,1.,1.,1./
      data f1i/
     +  0.000000e+00,0.414985e-02,0.154894e-01,0.310312e-01
     + ,0.464438e-01,0.569008e-01,0.580763e-01,0.468264e-01
     + ,0.217924e-01,-0.163419e-01,-0.651205e-01,-0.120503e+00
     + ,-0.178272e+00,-0.233580e+00,-0.282442e+00,-0.321901e+00
     + ,-0.350115e+00,-0.366534e+00,-0.371831e+00,-0.367378e+00
     + ,-0.354994e+00,-0.314803e+00,-0.266539e+00,-0.220551e+00
     + ,-0.181546e+00,-0.150427e+00,-0.126404e+00,-0.107830e+00
     + ,-0.933106e-01,-0.817375e-01,-0.723389e-01,-0.436650e-01
     + ,-0.294700e-01,-0.212940e-01,-0.161406e-01,-0.126604e-01
     + ,-0.102042e-01,-0.840465e-02,-0.704261e-02,-0.598886e-02/
      data f2i/
     +  0.0,0.121500e-01,0.454999e-01,0.913000e-01
     + ,0.137300e+00,0.171400e+00,0.183900e+00,0.170300e+00
     + ,0.132200e+00,0.763000e-01,0.126500e-01,-0.473500e-01
     + ,-0.936000e-01,-0.119750e+00,-0.123450e+00,-0.106300e+00
     + ,-0.732800e-01,-0.312400e-01,0.128450e-01,0.528800e-01
     + ,0.844100e-01,0.114710e+00,0.106200e+00,0.765830e-01
     + ,0.435800e-01,0.173950e-01,0.695001e-03,-0.809500e-02
     + ,-0.117355e-01,-0.125449e-01,-0.120280e-01,-0.686530e-02
     + ,-0.385275e-02,-0.231115e-02,-0.147056e-02,-0.982480e-03
     + ,-0.682440e-03,-0.489715e-03,-0.361190e-03,-0.272582e-03/
*
*     ------------------------------------------------------------------
*
*
      pi=4.*atan(1.0)
      twopi=2.*pi

* *** compute theta angle from moliere distribution
      chic  = sqrt(chic2)
      costh=1.
      sinth=0.
      th   =0.
      if(omega.le.eneper)go to 90
      cnst=log(omega)
      b=5.
c
      do 10 l=1,10
         if(abs(b).lt.1.e-10)then
            b=1.e-10
         endif
         db=-(b-log(abs(b))-cnst)/(1.-1./b)
         b=b+db
         if(abs(db).le.0.0001)go to 20
   10 continue
c
      go to 90
c
   20 continue
      if(b.le.0.)go to 90
      binv = 1./b
      tint(1) = 0.
c
      do 30 ja=2,4
         tint(ja)=f0i(ja)+(f1i(ja)+f2i(ja)*binv)*binv
   30 continue
c
      nmax = 4
   40 continue
      xint = unirn(dummy)
c
      do 50 na=3,40
         if(na.gt.nmax) then
            tint(na)=f0i(na)+(f1i(na)+f2i(na)*binv)*binv
            nmax=na
         endif
         if(xint.le.tint(na-1)) go to 60
   50 continue
c
      if(xint.le.tint(40)) then
         na=40
         goto 60
      else
         tmp=1.-(1.-b*(1.-xint))**5
         if(tmp.le.0.)go to 40
         thri=5./tmp
         go to 80
      endif
c
   60 continue
      na = max(na-1,3)
      na3 = na-3
c
      do 70 m=1,4
         na3m=na3+m
         arg(m)=tint(na3m)
         val(m)=thred(na3m)**2
   70 continue
c

      call polint(arg,val,4,xint,thri,err)    ! num rec
c
   80 continue

      th = chic * sqrt( abs( b * thri ) )


      if(th.gt.pi)go to 40

      sinth = sin(th)

      test=th*(unirn(dummy))**2
      if(test.gt.sinth)go to 40

      goto 100

   90 continue
*
* *** calculate sine and cosine of a random angle between 0 and 360 deg
*
  100 phi = unirn(dummy) * twopi


      dthxz = sinth * cos(phi)
      dthyz = sinth * sin(phi)

      return
      end

c  *********************************************************
      subroutine ruthscat(omega,chia2,dthxz,dthyz)
c
c     omega : number of scatters
c     chia2 : coulomb screening angle (squared)
c

      implicit real*8 (a-h,o-z)
c
      omega0 = 1.167 * omega
      nsc = poidev(omega0,idum)    ! num rec
      dthxz = 0.
      dthyz = 0.
      if( nsc .le. 0 ) go to 900
c
      do i=1,nsc
c
90     continue
       rn = unirn(dummy)
       thet = sqrt( chia2 * ((1./rn) - 1.) )

       tth = tan(thet)

       twopi=2.*4.*atan(1.0)
       phi = twopi* unirn(dummy)


        dthxz = dthxz + sin(tth) * cos(phi)
        dthyz = dthyz + sin(tth) * sin(phi)

      end do
c
900   continue
      return
      end

c  *********************************************************
      function poidev(xm,idum)

      implicit real*8 (a-h,o-z)

      integer idum
      parameter (pi=3.141592654)
      save alxm,g,oldm,sq
!$OMP THREADPRIVATE(alxm,g,oldm,sq)
      data oldm /-1./
      if (xm.lt.12.)then
        if (xm.ne.oldm) then
          oldm=xm
          g=exp(-xm)
        endif
        em=-1
        t=1.
2       em=em+1.
        t=t*unirn(dummy)
        if (t.gt.g) goto 2
      else
        if (xm.ne.oldm) then
          oldm=xm
          sq=sqrt(2.*xm)
          alxm=log(xm)
          g=xm*alxm-gammln(xm+1.)
        endif
1       y=tan(pi*unirn(dummy))
        em=sq*y+xm
        if (em.lt.0.) goto 1
        em=int(em)
        t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
        if (unirn(dummy).gt.t) goto 1
      endif
      poidev=em
      return
      end
c  ********************************************************
      subroutine polint(xa,ya,n,x,y,dy)

      implicit real*8 (a-h,o-z)

      integer n,nmax
      real*8 dy,x,y,xa(n),ya(n)
      parameter (nmax=10)
      integer i,m,ns
      real*8 den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) write(*,*) 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
c  ********************************************************
      FUNCTION gammln(xx)

      implicit real*8 (a-h,o-z)

      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
!$OMP THREADPRIVATE(cof,stp)
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
