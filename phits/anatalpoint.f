!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_point(m,ntf,
     &                             np,nr,ne,nm,nt,
     &                             sumfactor,weightRate,
     &                             trRES,tranatal,nfile,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by T.Miura on 2021/09/10                                 *
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
      integer           nr
      integer           ne
      integer           nm
      integer           nt

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,nt,nr,nm,2)
      integer           nfile
      double precision  tranatal(np,ne,nt,nr,nm,2*nfile)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           im
      integer           ir
      integer           it
      integer           ie
      integer           ip

      double precision  Xa, sigx
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
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
     &                        trRES(ip,ie,it,ir,im,1),
     &                        trRES(ip,ie,it,ir,im,2),
     &                        1.0d+0)
              trRES(ip,ie,it,ir,im,1) = Xa
              trRES(ip,ie,it,ir,im,2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(ip,ie,it,ir,im,iat(1,ntf)) =
     &           trRES(ip,ie,it,ir,im,1)

              ! (sig_xj)**2
              tranatal(ip,ie,it,ir,im,iat(2,ntf)) =
     &        (  trRES(ip,ie,it,ir,im,2)
     &         * trRES(ip,ie,it,ir,im,1) )**2

              ! sig_x
              tranatal(ip,ie,it,ir,im,iat(2,ntf)) = dsqrt(
     &        tranatal(ip,ie,it,ir,im,iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(ip,ie,it,ir,im,iat(2,ntf)) =
     &       (tranatal(ip,ie,it,ir,im,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(ip,ie,it,ir,im,iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(ip,ie,it,ir,im,iat(1,ntf)) =
     &        tranatal(ip,ie,it,ir,im,iat(1,ntf)) * resc2(m)

            if( tranatal(ip,ie,it,ir,im,iat(1,ntf)) > cmax )
     &            cmax = tranatal(ip,ie,it,ir,im,iat(1,ntf))
            if( tranatal(ip,ie,it,ir,im,iat(1,ntf)) < cmin )
     &            cmin = tranatal(ip,ie,it,ir,im,iat(1,ntf))

        end do     ! ip loop end
       end do      ! ie loop end
       end do      ! it loop end
      end do       ! ir loop end
      end do       ! im loop end

! sumover
      call anatal_calc_anova_sub(m,ntf,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return
      end subroutine anatal_calc_anova_point



************************************************************************
*                                                                      *
      subroutine anatal_ptpointp(m,np,nr,ne,nm,nt,eb,tb,tr,
     &                  nfile,weightRate,idasa,manatally)
*                                                                      *
*       output the point tally                                         *
*       created by T.Miura on 2021/09/13                               *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      common /stat / istdev, irestart, ireschk
      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      include 'param.inc'

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
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
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

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall35/ iteps(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall62/ itpon(itlmax), rtpon(itlmax,20,4)

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

      common /tall67/ itstd(itlmax), rtstd(itlmax)

*-----------------------------------------------------------------------
      common /tall78/ itism(itlmax), itist(10,itlmax), itjst(10,itlmax),
     &                itkst(10,itlmax), itstt(itlmax), itsdd(itlmax)

      dimension idas(1)
      equivalence ( das, idas )

      common /istcut/ ist_cut, ist_bat
*-----------------------------------------------------------------------

      dimension   eb(ne+1)
!      dimension   ew(ne)
      dimension   tb(nt+1)
!      dimension   tw(nt)
      dimension   tr(np,ne,nt,nr,nm,2*nfile)

      integer     nfile, manatally
!      dimension   weightRate(nfile), rdata(2,nfile)
      dimension   weightRate(nfile)

      integer irst,nrst, itmpdata
!      dimension   anatalrst(np,ne+1,1,nt+1,nr*1*1,nm+1,2*nfile+3)
! sumover
!      dimension   anatalrst(np,ne+1,1,nt+1,(nr+1)*1*1,nm+1,2*nfile+3)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)
      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9  ! S.H. 2021.8.9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*21,cijaxs*21,cijaxs2*21
      integer iDaxis
!      dimension   dlr(nr+1) !  integer lr(nr) ---> double dlr(nr)
      real(8),allocatable :: ew(:),tw(:),rdata(:,:),dlr(:)

*-----------------------------------------------------------------------

      dimension tott(6,2)

      character hsunit(14)*29

      data hsunit( 1) / '[1/cm^2/source]              '/
      data hsunit( 2) / '[1/cm^2/MeV/source]          '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]     '/
      data hsunit(11) / '[1/cm^2/nsec/source]         '/
      data hsunit(12) / '[1/cm^2/nsec/MeV/source]     '/
      data hsunit(13) / '[1/cm^2/nsec/Lethargy/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

      character rpa*1
      data rpa /'}'/
      character yen*1

      dimension dt_one(1)
      data dt_one /1.0d0/
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------
      yen  = char(92)

      nrst = 2*nfile+3 ! S.H. 2021.8.15
      allocate (anatalrst(np,ne+1,1,nt+1,(nr+1)*1*1,nm+1,nrst))

      allocate (ew(ne),tw(nt),rdata(2,nfile),dlr(nr+1))

*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 2) = '[1/cm^2/(MeV/n)/source]      '
         hsunit(12) = '[1/cm^2/nsec/(MeV/n)/source] '
      end if

*-----------------------------------------------------------------------

      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*           itunt(m) = 1,  11   : /cm^2/source
*                                 /cm^2/sec/source
*                    = 2   12   : /cm^2/MeV/source
*                                 /cm^2/MeV/nsec/source
*                    = 3   13   : /cm^2/Lethargy/source
*                                 /cm^2/Lethargy/nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq. 11 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq. 2 .or. itunt(m) .eq. 12 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq. 3 .or. itunt(m) .eq. 13 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 11, 12, 13 : /.../nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .gt. 10 ) then

               do i = 1, nt

                  tw(i) = tb(i+1) - tb(i)

               end do
               tw_sum = tb(nt+1) - tb(1)

            else

               do i = 1, nt

                  tw(i) = 1.d+0

               end do
               tw_sum = 1.0d0

            end if

*-----------------------------------------------------------------------

               do ir = 1, nr
                  dlr(ir) = dble(ir)
               end do
               dlr(nr+1) = dble(nr)

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

               if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

                  c1 = 1.0d+0 / rsouin

               else

                  c1 = 0d0

               end if

            do 100 im = 1, nm
            do 100 ir = 1, nr
            do 100 it = 1, nt
            do 100 ie = 1, ne
            do 100 ip = 1, np
            do ntf = 1, nfile

               if( tr(ip,ie,it,ir,im,iat(1,ntf)) .gt. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,ir,im,iat(1,ntf)),
     &                            tr(ip,ie,it,ir,im,iat(2,ntf)),
     &                            rtfac(m)/ew(ie)/tw(it))

                  tr(ip,ie,it,ir,im,iat(1,ntf)) = Xa
                  tr(ip,ie,it,ir,im,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,it,ir,im,iat(1,ntf)) .gt. cmax )
     &                        cmax = tr(ip,ie,it,ir,im,iat(1,ntf))

                  if( tr(ip,ie,it,ir,im,iat(1,ntf)) .lt. cmin )
     &                        cmin = tr(ip,ie,it,ir,im,iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,it,ir,im,iat(2,ntf)) = 0.0

               end if

! sumover
               call ptpointp_sumover_stdev_ntf(0,m,ntf,
     &               ip,ie,it,ir,im,
     &               rtfac(m),ew(ie),tw(it),
     &               ew_sum,tw_sum)

            end do
  100       continue

!$OMP PARALLEL
!$OMP& private(ipomp,ir,im,it,ie,ip,ntf,isdz,rdata,answer,rerr)
!$OMP& private(iax,ir_a,it_a,ie_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
            do ir = 1, nr
            do 101 im = 1, nm
            do 101 it = 1, nt
            do 101 ie = 1, ne
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,it,ir,im,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,it,ir,im,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(ip,ie,1,it,ir,im,1) = answer
              anatalrst(ip,ie,1,it,ir,im,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,it,ir,im,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,it,ir,im,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(ip,ie,1,it,ir,im,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,1,it,ir,im,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,1,it,ir,im,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,1,it,ir,im,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,1,it,ir,im,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,1,it,ir,im,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,ir,im,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,ir,im,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,ir,im,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(ip,ie,1,it,ir,im,iat(ioe,ntf))
     &              = tr(ip,ie,it,ir,im,iat(ioe,ntf))
                end do
               end do

            end if

!sumover
            do iax = 1,itaxn(m)

              ie_a = ie
              it_a = it
              ir_a = ir
              ido_ana = 0
              if(itaxs(m,iax) == 1 .and. ie == 1) then
                ie_a = ne + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 2 .and. ir == 1) then
                ir_a = nr + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 11 .and. it == 1) then
                it_a = nt + 1
                ido_ana = 1
              endif

              if(ido_ana == 1) then

                 call get_pointp_tr_sum_data(m,iax,nfile,
     &                 ip,ie,it,ir,im,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   anatalrst(ip,ie_a,1,it_a,ir_a,im,1) = answer
                   anatalrst(ip,ie_a,1,it_a,ir_a,im,2) = rerr

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  anatalrst(ip,ie_a,1,it_a,ir_a,im,1) ! mean
     &              = fmval
                 if ( fmval .gt. 0.0d0 ) then
                  anatalrst(ip,ie_a,1,it_a,ir_a,im,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                 else
                   anatalrst(ip,ie_a,1,it_a,ir_a,im,2) ! total uncertainty
     &                 = 0.0d0
                 end if

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      anatalrst(ip,ie_a,1,it_a,ir_a,im,iat(ioe,ntf))
     &                = rdata(ioe,ntf)
                    end do
                  end do

                end if

              endif
            enddo

  101       continue
      enddo  ! ir loop T.Sato 2021/04/15
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

      do 900 iax = 1, itfln(m)

            noe = 1

         if( itall .eq. 3 .and. itism(m) .gt. 0 ) then

            noe  = noe + 1

         end if

*-----------------------------------------------------------------------

        do 900 ioe = 1, noe

*-----------------------------------------------------------------------

             nstd = 0

         if ( ioe .eq. 1 ) then

            if( ( itall .eq. 2 .and. nobch .lt. maxbch )
     &           .or. itall .eq. 4 ) then

               call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &              nobch,maxbch,npe)

            else

               fname = ctfln(m,iax)

            end if

         else if ( ioe .eq. 2 ) then

            if( itall. eq. 3 .and. itism(m) .gt. 0 ) then

              call mk_2distfn(ctfln(m,iax),fname,itfll(m,iax))
              nstd = 1

            else

              goto 900

            end if

         end if

*-----------------------------------------------------------------------

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tponech(iot,m,iax,1)

*-----------------------------------------------------------------------
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        energy, time axis
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1
     &        .or. itaxs(m,iax) .eq. 2 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 1   ! energy
           case ( 2 )
            iDaxis = 9   ! time
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.15

          if ( itpon(m) .eq. 1 ) then ! point detector
             imesh = 5
          else ! ring detector
             imesh = 6
          end if

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,iMeVperu,imesh,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,       0,    nt,     nr,  0,  0,   nm, nrst,
     &         eb,ew,    1,1,  tb,tw, dlr,  1,  1,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

               inum = 0

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
               write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

           if( itaxs(m,iax) .eq. 1) then
               if( itunt(m) .eq. 3 .or. itunt(m) .eq. 13 .or.
     &             itety(m) .eq. 3 .or. itety(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( (itaxs(m,iax) .eq. 2) ) then
               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

           if ( manatally .eq. 0 ) then ! user defined analysis
            if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: n",12x,"x",12x,
     &                     1000(a10,"),hh0",a3," n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
            else
                write(iot,'( "h: n",12x,"x",12x,
     &                     1000(a1,i0,a9,"),hh0",a3
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if

               write(iot,'( "#  lower        upper  ",3x,
     &                    1000(a1,2x,a8,4x,"r.err "))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h: n",12x,"x",12x,
     &                     1000(a10,"),hh0",a3," n n n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
                write(iot,'( "h: n",12x,"x",12x,
     &                     1000(a1,i0,a9,"),hh0",a3
     &                ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if

               write(iot,'( "#  lower        upper  ",3x,
     &                    1000(a1,2x,a8,4x,"r.err(tot, syst, stat)"))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 2 ) then ! c-value dependence
            if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h:   x",10x,
     &                     1000(a10,"),",a4," n   "))')
     &                     ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
                write(iot,'( "h:   x",10x,
     &                     1000(a1,i0,a9,"),",a4
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if

               write(iot,'( "#  c-value   ",
     &                    1000(a1,2x,a8,4x,"r.err "))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
           end if

               do ip = 1, np
                     tott(ip,1) = 0.d+0
               end do

           voll = 0.0d0
           njaxs = nijaxs + 1

           do ijaxs = 1, nijaxs

            if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
                write(iot,'(1p2e13.4,1000(1pe13.4,0p3f8.4))')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,4)
     &               ,ip=1,np) ! frtati 2021/10/05

            else if ( manatally .eq. 2 ) then ! c-value dependence
                write(iot,'(1pe13.4,1000(1pe13.4,0pf8.4))')
     &               fgaxs(ijaxs),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05

            end if

            voll = voll + delvol(ianataldata,ijaxs)

                  do ip = 1, np
               vn = anataldata(ip,ianataldata,ijaxs,1)
     &              * delvol(ianataldata,ijaxs)
                     tott(ip,1) = tott(ip,1) + vn
                  end do

               end do

               do ip = 1, np
           if( itaxs(m,iax) .eq. 1 ) then
                     if( itunt(m) .eq.  2 .or. itunt(m) .eq.  3 .or.
     &                   itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then
                        tott(ip,1) = tott(ip,1) / ewtt
                     end if
           else if( (itaxs(m,iax) .eq. 2) ) then
                  if( itunt(m) .gt.  10 ) then
                        tott(ip,1) = tott(ip,1) / twtt
                  end if

           end if

               end do

! sumover
            do ip = 1, np
            do irst=1,2
              tott(ip,irst) = anataldata(ip,ianataldata,njaxs,irst)
            enddo
            enddo

               write(iot,'(/"#   sum over ",   13x ,
     &               1000(1pe13.4,0pf8.4))')
     &              (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

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
                write(iot,'(a15,"&=&",i5)')
     &               cij(ijtmp), idnint(fg(itmp))
                itmp = itmp - ij(ijtmp) + nij(ijtmp)
             else
                write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                     /1pe13.4,2x,"$--$",1pe13.4)')
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

*-----------------------------------------------------------------------

         end if

            call prestart(m,iot) !OBINATA(2012.6.13)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------
      deallocate (anatalrst)
      deallocate (ew,tw,rdata,dlr)

      return
      end subroutine anatal_ptpointp

!***********************************************************************
!                                                                      *
! sumover subroutine                                                   *
!                                                                      *
!***********************************************************************

!***********************************************************************
!                                                                      *
      subroutine get_pointp_tr_sum_data(m,iax,nfile,
     &           ip,ie,it,ir,im,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,it,ir,im
      real(8) :: rdata(2,nfile)

      integer :: ntf

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

      real(8),pointer :: p_sum(:)

        do ntf=1,nfile

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$         p_sum => tr00_sum(ianatalm_sum(m,iax,ntf):)
!$       else
           p_sum => tr0_sum(ianatalm_sum(m,iax,ntf):)
C for nonshared_tally option
!$       end if

        call get_pointp_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &     itmsh_sum(m,iax),itmst_sum(m,iax),
     &     ip,ie,it,ir,im,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_pointp_tr_sum_data_sub(tr_sum,
     &          np,ne,nt,nr,nm,ip,ie,it,ir,im,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,nt,nr,nm,ip,ie,it,ir,im
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,ne,nt,nr,nm,2)

      do i=1,2
        rdata(i)= tr_sum(ip,ie,it,ir,im,i)
      enddo

      end
