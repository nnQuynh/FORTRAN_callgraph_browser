************************************************************************
*                                                                      *
      subroutine check_tvlm(m,iax,jsn,jsi,dsin,idsi,ill,ilf,ierr)
*                                                                      *
*   m: the tally number, index of ital.                                *

************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall31/ iterl(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall68/ itstp(itlmax), itmth(itlmax), rtvr0(itlmax),
     &                rtvx0(itlmax), rtvy0(itlmax), rtvz0(itlmax),
     &                rtvx1(itlmax), rtvy1(itlmax), rtvz1(itlmax)

      character irfile*100

      dimension idas(mdas*2)
      equivalence ( das, idas )


*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character title*80

*-----------------------------------------------------------------------

      dimension lschn(12), ischn(12)
      character schan(12)*8

      data icsu / 12 /

      data ( schan(i), i = 1, 12 ) /
     &    'mesh    ','file    ','resfile ','stdcut  ','s-type  ',
     &    'x0      ','y0      ','z0      ','r0      ','x1      ',
     &    'y1      ','z1      '/

      data ( lschn(i), i = 1, 12 ) /
     &     4,         4,         6,         6,         6,
     &     2,         2,         2,         2,         2,
     &     2,         2/


*-----------------------------------------------------------------------

      dimension iptyp(6), ipnkf(6)
      dimension imtyp(6,6), imnkf(6,6)
      dimension jstyp(6), jnkf0(6)
      dimension icount(9)

      dimension iaxis(6)
      character ifile(6)*100
      dimension lfile(6)


*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      character dkam*9

      character tname*10
      data      tname /'[t-volume]'/

      logical deqn4

      character rglnrf*200, rglnech*200

*-----------------------------------------------------------------------

      dimension jptyp(6), jpnkf(6)
      dimension kptyp(6,6), kpnkf(6,6)
      dimension ktln(6), ktli(6), ktls(6), kmst(6), knpat(6)
      dimension dkmax(6)
      dimension imst(6), kimst(6,6)
      dimension vtrs(13)

*-----------------------------------------------------------------------

      character*10 cepn(100)
      dimension lepn(100)
      character*20 ctl(2)

*-----------------------------------------------------------------------

            ierr  = 0
            jpn   = 0

            infil = 0

            lrfile = 0 !OBINATA
            irfflg = 0 !OBINATA

            iechrl = 72

            x0 = 0.d0
            y0 = 0.d0
            z0 = 0.d0

            stdcut = -1.0

         do i = 1, icsu

            ischn(i) = 0

         end do

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( jpn  .ne. 0 ) goto 800
               if( jpn  .eq. 3 ) goto 800

               if( iskip .ne. 0 ) goto 140

  150 continue

               if( ierr .ne. 0 ) goto 800
               if( jpn  .eq. 3 ) goto 800

*-----------------------------------------------------------------------
*        head of the section
*-----------------------------------------------------------------------

            if( chlw(i1:i1) .eq. '[' .and.
     &          chcm(i1:10) .eq. tname ) then

cKN 2017/01/06 it alows [t-volulme] off

               goto 140

            end if

*-----------------------------------------------------------------------
*        end of the section
*-----------------------------------------------------------------------

            if( chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 800

            end if

*-----------------------------------------------------------------------
*        identify the parameters
*-----------------------------------------------------------------------

            icl = i1

  200    continue

            chlc = chlw
            call chcomp(chlc,icl,i3,i5)

         do i = 1, icsu

            il = icl + lschn(i) - 1

            if( chlc(icl:il) .eq. schan(i)(1:lschn(i)) ) goto 100

         end do

               goto 800

*-----------------------------------------------------------------------
*        read value of parameters
*-----------------------------------------------------------------------

  100    continue

               ipm = i

               ischn( ipm ) = ischn( ipm ) + 1

            if( ipm .ne.  2 .and. ipm .ne. 3 .and.
     &          ipm .ne.  5 .and.
     &          ischn( ipm ) .gt. 1 ) goto 989

               ic = inumc(chlw,il+1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 997

               icl = inumc(chlw,ic,i3,';') - 1

*-----------------------------------------------------------------------
*        mesh = region, r-z, or xyz
*-----------------------------------------------------------------------

         if( ipm .eq. 1 ) then

            if( chlw(ic:ic+2) .eq. 'reg' ) then

               imesh = 1

            else

               goto 998

            end if

            if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               call tregion0(1,jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      ntrn,mtrn,ndsm,nvol,ivl,irvl,0,
     &                      rglnrf)

                  if( jpn  .eq. 3 ) goto 800
                  if( ntrn .lt. -1 ) goto 996

                  goto 150

*-----------------------------------------------------------------------
*        file name
*-----------------------------------------------------------------------

         else if( ipm .eq. 2 ) then

  700          infil = infil + 1

               if( infil .gt. 6 ) goto 990

               icf = min( inumc(chlw,ic,icl,' ') - 1, icl )

               lfile(infil) = icf - ic + 1
               ifile(infil)(1:icf-ic+1) = chin(ic:icf)

               ic = jnumc(chlw,icf+2,icl)

               if( ic .le. icl ) goto 700

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        restart file name
*-----------------------------------------------------------------------

         else if( ipm .eq. 3 ) then

               icf = min( inumc(chlw,ic,icl,' ') - 1, icl )

               lrfile = icf-ic+1

               irfile(1:lrfile) = chin(ic:icf)


*-----------------------------------------------------------------------
*        stdcut
*-----------------------------------------------------------------------

         else if( ipm .eq. 4 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               stdcut = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        s-type
*-----------------------------------------------------------------------

         else if( ipm .eq. 5 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               istyp = nint(cvvv)

               if( istyp .ne. 1 .and. istyp .ne. 2 ) goto 996

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        x0
*-----------------------------------------------------------------------

         else if( ipm .eq. 6 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               x0  = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        y0
*-----------------------------------------------------------------------

         else if( ipm .eq. 7 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               y0  = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        z0
*-----------------------------------------------------------------------

         else if( ipm .eq. 8 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               z0  = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        r0
*-----------------------------------------------------------------------

         else if( ipm .eq. 9 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               r0  = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        x1
*-----------------------------------------------------------------------

         else if( ipm .eq. 10 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               x1  = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        y1
*-----------------------------------------------------------------------

         else if( ipm .eq. 11 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               y1  = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        z1
*-----------------------------------------------------------------------

         else if( ipm .eq. 12 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               z1  = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------

         end if

      goto 140

  800 continue

*-----------------------------------------------------------------------
*     check
*-----------------------------------------------------------------------

      iec = 0

      call check_mesh(m,imesh,iec,cepn,lepn,ierr)

      if ( itmsh(m) .eq. 1 ) then

      call check_reg(m,rglnrf,iec,cepn,lepn,ierr)

      end if

      if( itstp(m) .eq. 1 ) then

        call check_x0y0z0(m,x0,y0,z0,iec,cepn,lepn,ierr)
        call check_r0(m,r0,iec,cepn,lepn,ierr)

      else if( itstp(m) .eq. 2 ) then

        call check_x0y0z0(m,x0,y0,z0,iec,cepn,lepn,ierr)
        call check_x1y1z1(m,x1,y1,z1,iec,cepn,lepn,ierr)

      end if

*-----------------------------------------------------------------------

      call judge_tall_check(iec,cepn,lepn,dsin(jsn),idsi(jsn),ierr)

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  989    m_err = 'Double definition of the parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:481/R:check_tvlm/F:restvlm.f'
         goto 999

  990    m_err = 'Too many file name in tally '//tname
         ErrCha = ''
         ErrID = 'L:486/R:check_tvlm/F:restvlm.f'
         goto 999

  993    m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:491/R:check_tvlm/F:restvlm.f'
         goto 999

  996    m_err = 's-type should be 1 or 2 in tally '//tname
         ErrCha = ''
         ErrID = 'L:496/R:check_tvlm/F:restvlm.f'
         goto 999

  997    m_err = 'Description of parameter is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:501/R:check_tvlm/F:restvlm.f'
         goto 999

  998    m_err = 'Unknown mesh parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:506/R:check_tvlm/F:restvlm.f'
         goto 999

*-----------------------------------------------------------------------

  999 continue

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

         write(ErrCha,*) 'Error: ' // m_err
         ErrID = 'L:518/R:check_tvlm/F:restvlm.f' !E67_001_001
         call ErrWrite(ErrID,ErrCha)

*-----------------------------------------------------------------------

      return

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_tvlm(m,iax,ierr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        use RESTALMOD
      use moddas_region

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /mpi00/ npe, me

*-----------------------------------------------------------------------

        common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
        common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                  rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

        dimension     idas(mdas*2)
        equivalence ( das, idas )

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

        call open_resfile(m,noe,jsn,jsi,dsin,idsi,ill,ilf,newtall,ierr)

        if ( newtall .ne. 0 ) goto 900  !! it's new tally
        if ( ierr    .ne. 0 ) goto 900

*-----------------------------------------------------------------------
*   check tally
*-----------------------------------------------------------------------

        do 800 ioe = 1, noe

        if(ireschk.eq.0) then  ! T.Sato 2013/10/19
        call check_tvlm(m,iax,jsn(ioe),jsi(ioe),
     &                  dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                  ierr)
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
     &                    itrnv(m),idas(itriv(m)),das(itrrv(m)))

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

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_vlmreg(m,jsn,jsi,dsin,idsi,ill,ilf,
     &                        nr,mr,kr,tr,
     &                        nvl,ivl,rvl)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

        common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
        character ctfln*100

        common /cparm/  maxbch,maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

*-----------------------------------------------------------------------

        dimension   vl(nr)
        dimension   lr(nr)
        dimension   ivl(nvl)
        dimension   rvl(nvl)
        dimension   kr(mr)
        dimension   tr(nr,2)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200

*-----------------------------------------------------------------------
*        reg axis
*-----------------------------------------------------------------------

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '[volume]',8,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

               read(jsi,*)

             do ir = 1, nr

                  read(jsi,'(i5,2x,i7,1p1e13.4,0pf8.4)')
     &                     idmm0, idmm1,
     &                     (tr(ir,k),k=1,2)

             end do

*-----------------------------------------------------------------------

  900 continue

      end subroutine



************************************************************************
*                                                                      *
      subroutine restore_vlmreg(m,nr,mr,kr,tr,nvl,ivl,rvl)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

        dimension   lr(nr)
        dimension   vl(nr)
        dimension   ivl(nvl)
        dimension   rvl(nvl)
        dimension   kr(mr)
        dimension   tr(nr,2)

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

        call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)

*-----------------------------------------------------------------------

          uni = 1.d0

        do 100 ir = 1, nr

          call invert_stdev(m,A,B,
     &                      tr(ir,1),
     &                      tr(ir,2),
     &                      uni)

          tr(ir,1) = A
          tr(ir,2) = B

  100 continue

      end subroutine



