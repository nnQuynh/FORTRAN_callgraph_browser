************************************************************************
*                                                                      *
      subroutine seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        ctln,ltln,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character*(*) ctln
        character chin*200, chlw*200, chcm*200

*-----------------------------------------------------------------------

  100   continue

        ierr  = 0
        jpn   = 0
        iskip = 0

        call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &             jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

        if ( ierr  .ne. 0 ) goto 900
        if ( jpn   .eq. 3 ) goto 900
        if ( iskip .ne. 0 ) goto 100

        if ( chcm(i1:i1+ltln-1) .ne. ctln(1:ltln) ) goto 100

*-----------------------------------------------------------------------

  900   continue

      end subroutine

************************************************************************
*                                                                      *
      subroutine seek_resfile2(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        ctln,ltln,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character*(*) ctln
        character chin*200, chlw*200, chcm*200

*-----------------------------------------------------------------------

  100   continue

        ierr  = 0
        jpn   = 0
        iskip = 0

        call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &             jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

        if ( ierr  .ne. 0 ) goto 900
        if ( jpn   .eq. 3 ) goto 900
        if ( iskip .ne. 0 ) goto 100

        if ( chin(1:i2) .ne. ctln(1:ltln) ) goto 100

*-----------------------------------------------------------------------

  900   continue

      end subroutine

************************************************************************
*                                                                      *
      subroutine seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200

*-----------------------------------------------------------------------
*  seek to 'newpage:'
*-----------------------------------------------------------------------
  100   continue

        ierr  = 0
        jpn   = 0
        iskip = 0

        call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &             jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

        if ( ierr  .ne. 0 ) goto 900
        if ( jpn   .eq. 3 ) goto 900
        if ( iskip .ne. 0 ) goto 100

        if ( chcm(i1:i1+8) .eq. '#newpage:' ) goto 900
        if ( chcm(i1:i1+7) .eq. 'newpage:'  ) goto 900

        goto 100

  900   continue

      end subroutine



************************************************************************
*                                                                      *
      subroutine seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)
*                                                                      *
* created by Kitamura on 23/03/31                                      *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200

*-----------------------------------------------------------------------
*  seek to '#restart' : restart maximum for normalization
*-----------------------------------------------------------------------
  100   continue

        ierr  = 0
        jpn   = 0
        iskip = 0
        facmx = 0.d0

        call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &             jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

        if ( ierr  .ne. 0 ) goto 900
        if ( jpn   .eq. 3 ) goto 900
        if ( iskip .ne. 0 ) goto 100

        if ( chcm(i1:i1+7) .eq. '#restart' ) then

            read( chcm(33:45),'(e14.7)') facmx
            goto 900

        end if

        goto 100

  900   continue

      end subroutine



************************************************************************
*                                                                      *
      subroutine seek_resfile_fa2(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                            facmr,facmz,ierr)
*                                                                      *
* created by Kitamura on 23/03/31                                      *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200

*-----------------------------------------------------------------------
*  seek to '#restart' : restart maximum for normalization
*-----------------------------------------------------------------------
  100   continue

        ierr  = 0
        jpn   = 0
        iskip = 0
        facmr = 0.d0
        facmz = 0.d0

        call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &             jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

        if ( ierr  .ne. 0 ) goto 900
        if ( jpn   .eq. 3 ) goto 900
        if ( iskip .ne. 0 ) goto 100

        if ( chcm(i1:i1+7) .eq. '#restart' ) then

            read( chcm(33:45),'(e14.7)') facmr
            read( chcm(46:58),'(e14.7)') facmz
            goto 900

        end if

        goto 100

  900   continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_trline(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       chlw,igcl,ics,il,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200

*-----------------------------------------------------------------------

        dimension ic(0:igcl)

*-----------------------------------------------------------------------

  100   continue

        call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$#',
     &             jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

        if ( ierr  .ne. 0 ) goto 900
        if ( jpn   .eq. 3 ) goto 900
        if ( iskip .ne. 0 ) goto 100

        !! ignore columns
        ic(0) = i1
        do i = 1, igcl
          call snum(chlw,ic(i-1),i3,ic(i),prn,ierr)
        end do

        ics = ic(igcl)
        il  = i3

*-----------------------------------------------------------------------

  900   continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine sel_dc2(m,dc2,ldc2)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)

        character*(*) dc2

*-----------------------------------------------------------------------

        if( ittwo(m) .eq. 1 ) then
          dc2  = 'h2:'
          ldc2 = 3
        else if( ittwo(m) .eq. 2 ) then
          dc2 = 'hd:'
          ldc2 = 3
        else if( ittwo(m) .eq. 3 ) then
          dc2 = 'hc:'
          ldc2 = 3
        else if( ittwo(m) .eq. 6 ) then
          dc2 = 'hd2:'
          ldc2 = 4
        else if( ittwo(m) .eq. 7 ) then
          dc2 = 'hc2:'
          ldc2 = 4
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine mk_2distfn(cbfn,cefn,ibfll)
*                                                                      *
*   make a file name of 2d-type tally's error.                         *
*                                                                      *
*   cbfn  : base file name. ( file / resfile )                         *
*   cefn  : error file name.                                           *
*   ibfll : length of cbfn.                                            *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        character cbfn*100
        character cefn*100

*-----------------------------------------------------------------------

          cefn = ' '

        do i = ibfll, 1, -1
          if ( cbfn(i:i) .eq. '.' ) goto 100
        end do

        if ( i .eq. 0 ) i = ibfll+1

  100 continue

        cefn(1:i-1) = cbfn(1:i-1)
        cefn(i:i+3) = '_StD'
        if ( i .lt. ibfll ) then
          cefn(i+4:ibfll+4) = cbfn(i:ibfll)
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine mk_2dnumfn(cbfn,cefn,ibfll,nobch,maxbch,npe)
*                                                                      *
*   make a file name with batch number for tally.                      *
*                                                                      *
*   cbfn  : base file name. ( file / resfile )                         *
*   cefn  : file name with baatch number.                              *
*   ibfll : length of cbfn.                                            *
*   nobch : the current batch number.                                  *
*   maxbch : maxbch in input file.                                     *
*   npe-1 : the number of the MPI parallel calculation                 *
*           used for the transport calculation.                        *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        character cbfn*100
        character cefn*100
        character fnume*10
        character cform*8
        integer iorder
        character corder*2
        integer nobch, maxbch, npe
        integer nobchnumfn, maxbchnumfn

*-----------------------------------------------------------------------

          cefn = ' '

        do i = ibfll, 1, -1
          if ( cbfn(i:i) .eq. '.' ) goto 100
        end do

        if ( i .eq. 0 ) i = ibfll+1

  100 continue

        cefn(1:i-1) = cbfn(1:i-1)

        if ( npe .le. 1 ) then
           nobchnumfn = nobch
           maxbchnumfn = maxbch
        else
           nobchnumfn = nobch/(npe-1)
           maxbchnumfn = maxbch/(npe-1)
        end if
        if ( nobchnumfn .lt. 1 ) nobchnumfn = 1
        iorder = aint(log10(real(maxbchnumfn))) + 1
        if ( iorder .lt. 3) iorder = 3
        write(corder,'(i2)') iorder
        corder = trim(adjustl(corder))
        cform='(i'//corder//'.'//corder//')'
        write(fnume,cform) nobchnumfn

        cefn(i:i+iorder) = '_'//fnume(1:iorder)
        if ( i .lt. ibfll ) then
          cefn(i+1+iorder:ibfll+1+iorder) = cbfn(i:ibfll)
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine mk_2dnerfn(cbfn,cefn,ibfll,fnume)
*                                                                      *
*   make a file name with batch nmber and error for tally.                       *
*                                                                      *
*   cbfn  : base file name. ( file / resfile )                         *
*   cefn  : file name with baatch number.                              *
*   ibfll : length of cbfn.                                            *
*   fnume : batch number.                                              *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        character cbfn*100
        character cefn*100
        character fnume*3

*-----------------------------------------------------------------------

          cefn = ' '

        do i = ibfll, 1, -1
          if ( cbfn(i:i) .eq. '.' ) goto 100
        end do

        if ( i .eq. 0 ) i = ibfll+1

  100 continue

        cefn(1:i-1) = cbfn(1:i-1)
        cefn(i:i+7) = '_'//fnume//'_err'
        if ( i .lt. ibfll ) then
          cefn(i+8:ibfll+8) = cbfn(i:ibfll)
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine mk_2derrfn(cbfn,cefn,ibfll)
*                                                                      *
*   make a file name of 2d-type tally's error.                         *
*                                                                      *
*   cbfn  : base file name. ( file / resfile )                         *
*   cefn  : error file name.                                           *
*   ibfll : length of cbfn.                                            *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        character cbfn*100
        character cefn*100

*-----------------------------------------------------------------------

          cefn = ' '

        do i = ibfll, 1, -1
          if ( cbfn(i:i) .eq. '.' ) goto 100
        end do

        if ( i .eq. 0 ) i = ibfll+1

  100 continue

        cefn(1:i-1) = cbfn(1:i-1)
        cefn(i:i+3) = '_err'
        if ( i .lt. ibfll ) then
          cefn(i+4:ibfll+4) = cbfn(i:ibfll)
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_point(m,npont,itpnt,stpon,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)

        common /tall62/ itpon(itlmax), rtpon(itlmax,20,4)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        dimension stpon(20,4)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itmsh(m) .ne. npont ) then
          iec = iec + 1
          cepn(iec) = 'number'
          lepn(iec) = 6
          ierr = 1
        end if

        if ( itpon(m) .ne. itpnt ) then
          iec = iec + 1
          cepn(iec) = 'poin/ring'
          lepn(iec) = 9
          ierr = 1
        end if

        do i = 1, npont

           if( rtpon(m,i,1) .ne. stpon(i,1) .or.
     &         rtpon(m,i,2) .ne. stpon(i,2) .or.
     &         rtpon(m,i,3) .ne. stpon(i,3) .or.
     &         rtpon(m,i,4) .ne. stpon(i,4) ) then

             iec = iec + 1
             cepn(iec) = 'values'
             lepn(iec) = 6
             ierr = 1

           end if

        end do

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_type(ct,iec,cepn,lepn,
     &                      itty,rtma,rtmi,itnm,
     &                      itp, rmax,rmin,inm )
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)
        character ct

*-----------------------------------------------------------------------
        if(itp.eq.0) return ! T.Sato 2024/01/15 not necessary to check when type is not specified in the old tally file

        if ( itty .ne. itp ) then
          iec = iec + 1
          cepn(iec) = ct // '-type'
          lepn(iec) = 6
        else if ( itty .ge. 2 ) then

          if ( rtma .ne. rmax ) then
            iec = iec + 1
            cepn(iec) = ct // 'max'
            lepn(iec) = 4
          end if

          if ( rtmi .ne. rmin ) then
            iec = iec + 1
            cepn(iec) = ct // 'min'
            lepn(iec) = 4
          end if

          if ( itnm .ne. inm ) then
            iec = iec + 1
            cepn(iec) = 'n' // ct
            lepn(iec) = 2
          end if

        end if

*-----------------------------------------------------------------------

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_setype(m,iec,cepn,lepn,
     &                        itp, emax,emin,ine )
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                  rtemi(itlmax), rtema(itlmax), rtedl(itlmax)

        character*10 cepn(100)
        dimension lepn(100)
        character ct

*-----------------------------------------------------------------------

        if ( itety(m) .ne. itp ) then
          iec = iec + 1
          cepn(iec) = 'se-type'
          lepn(iec) = 7
        else if ( itety(m) .ge. 2 ) then

          if ( rtema(m) .ne. emax ) then
            iec = iec + 1
            cepn(iec) = 'emax'
            lepn(iec) = 4
          end if

          if ( rtemi(m) .ne. emin ) then
            iec = iec + 1
            cepn(iec) = 'emin'
            lepn(iec) = 4
          end if

          if ( itenm(m) .ne. ine ) then
            iec = iec + 1
            cepn(iec) = 'ne'
            lepn(iec) = 2
          end if

        end if

*-----------------------------------------------------------------------

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_e1type(m,iec,cepn,lepn,
     &                        itp, emax,emin,ine )
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                  rtemi(itlmax), rtema(itlmax), rtedl(itlmax)

        character*10 cepn(100)
        dimension lepn(100)
        character ct

*-----------------------------------------------------------------------

        if ( itety(m) .ne. itp ) then
          iec = iec + 1
          cepn(iec) = 'e1-type'
          lepn(iec) = 8
        else if ( itety(m) .ge. 2 ) then

          if ( rtema(m) .ne. emax ) then
            iec = iec + 1
            cepn(iec) = 'emax'
            lepn(iec) = 4
          end if

          if ( rtemi(m) .ne. emin ) then
            iec = iec + 1
            cepn(iec) = 'emin'
            lepn(iec) = 4
          end if

          if ( itenm(m) .ne. ine ) then
            iec = iec + 1
            cepn(iec) = 'ne'
            lepn(iec) = 2
          end if

        end if

*-----------------------------------------------------------------------

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_e2type(m,iec,cepn,lepn,
     &                        itp, emax,emin,ine )
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                  rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)

        character*10 cepn(100)
        dimension lepn(100)
        character ct

*-----------------------------------------------------------------------

        if ( itety2(m) .ne. itp ) then
          iec = iec + 1
          cepn(iec) = 'e1-type'
          lepn(iec) = 8
        else if ( itety2(m) .ge. 2 ) then

          if ( rtema2(m) .ne. emax ) then
            iec = iec + 1
            cepn(iec) = 'emax'
            lepn(iec) = 4
          end if

          if ( rtemi2(m) .ne. emin ) then
            iec = iec + 1
            cepn(iec) = 'emin'
            lepn(iec) = 4
          end if

          if ( itenm2(m) .ne. ine ) then
            iec = iec + 1
            cepn(iec) = 'ne'
            lepn(iec) = 2
          end if

        end if

*-----------------------------------------------------------------------

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_cdiam(m,cdiam,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

      common /tall57/ itsun(itlmax),rtdim(itlmax),rtucv(itlmax),
     &rtrho(itlmax),itmodel(itlmax) ! T.Sato 2022/08/14

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0
        if ( rtdim(m) .ne. cdiam ) then
          iec = iec + 1
          cepn(iec) = 'cdiam'
          lepn(iec) = 5
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_mesh(m,imesh,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0
        if ( itmsh(m) .ne. imesh ) then
          iec = iec + 1
          cepn(iec) = 'mesh'
          lepn(iec) = 4
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_reg(m,rglnrf,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************
      use moddas_region

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
        common /tall24/ itrcm(itlmax), itrcg(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
        common /tall31/ iterl(itlmax)
        common /tall70/ itrwgtsum(itlmax), itrgn1(itlmax),
     &                  itrgm1(itlmax), itrncd(itlmax)

        dimension idas(mdas*2) !2023.1.11 S.H. temporary for fbounds-check
        equivalence ( das, idas )

*-----------------------------------------------------------------------

        character rglnrf*200, rglnech*200

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        idas1 = mmmax
        idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
        idas3 = idas2 + itrgn(m)

        open( unit=34, status='scratch' )
        if( itrwgtsum(m) .eq. 1 ) then
          call echrg_wgtsum(m,1,iterl(m),0,1,34,itrcm(m),
     &               idas_itrcg(itrcg(m)),
     &               itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &               das(idas1),idas(idas2),
     &               itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &               idas3)
        else
          call echrg(1,iterl(m),0,1,34,itrcm(m),idas_itrcg(itrcg(m)),
     &               itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &               das(idas1),idas(idas2),
     &               itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &               idas3)
        endif
        rewind(34)
          read(34,'(a)') rglnech
        close(34)

        lrgrf = len_trim(rglnrf)
        lrgec = len_trim(rglnech)

        icn = inumc(rglnrf, 1,lrgrf,'=') + 1
        jcn = inumc(rglnech,1,lrgec,'=') + 1

  811 continue

        ic = icn
        jc = jcn

        call snum(rglnrf, ic,lrgrf,icn,prni,ierr)
        call snum(rglnech,jc,lrgec,jcn,prnj,jerr)

        irg = nint(prni)
        jrg = nint(prnj)

        if ( ierr .eq. jerr ) then
          if ( ierr .ne. 0  ) goto 812
          if ( irg .eq. jrg ) goto 811
        end if

        ! error in the parameter
        iec = iec + 1
        cepn(iec) = 'reg'
        lepn(iec) = 3

  812 continue

        return

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_tet(m,ndim_mtetreg,mtetreg,
     &     iec,cepn,lepn,ierr)
*                                                                      *
*     Last Modified by T.Furuta on 2025/01/16
*                                                                      *
************************************************************************
      use moddas_region

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)

*-----------------------------------------------------------------------

        integer,intent(in) :: m
        integer,intent(inout) :: iec
        character(10),intent(inout) :: cepn(100)
        integer,intent(inout) :: lepn(100)
        integer,intent(out) :: ierr
        integer,intent(in) :: ndim_mtetreg
        integer,intent(in) :: mtetreg(ndim_mtetreg)

        logical iflag
*-----------------------------------------------------------------------
        ierr=0
        iflag=.false.
        ii=itreg(m)
        if(idas_itreg(ii).eq.mtetreg(1))then
         nr=mtetreg(3)
         if(idas_itreg(ii+2).eq.nr)then
          iflag=.true.
          do i=1,nr
           if(idas_itreg(ii+2+i).ne.mtetreg(3+i))then
            iflag=.false.
            exit
           endif
          enddo
         endif
        endif

        if(.not.iflag)then
         ierr=1
        ! error in the parameter
         iec = iec + 1
         cepn(iec) = 'tet'
         lepn(iec) = 3
        endif

        return

      end subroutine

************************************************************************
*                                                                      *
      subroutine check_x0y0(m,rzx0,rzy0,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall08/ rtrx0(itlmax), rtry0(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( rtrx0(m) .ne. rzx0 ) then
          iec = iec + 1
          cepn(iec) = 'x0'
          lepn(iec) = 2
          ierr = 1
        end if

        if ( rtry0(m) .ne. rzy0 ) then
          iec = iec + 1
          cepn(iec) = 'y0'
          lepn(iec) = 2
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_x0y0z0(m,rzx0,rzy0,rzz0,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
      common /tall68/ itstp(itlmax), itmth(itlmax), rtvr0(itlmax),
     &                rtvx0(itlmax), rtvy0(itlmax), rtvz0(itlmax),
     &                rtvx1(itlmax), rtvy1(itlmax), rtvz1(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( rtvx0(m) .ne. rzx0 ) then
          iec = iec + 1
          cepn(iec) = 'x0'
          lepn(iec) = 2
          ierr = 1
        end if

        if ( rtvy0(m) .ne. rzy0 ) then
          iec = iec + 1
          cepn(iec) = 'y0'
          lepn(iec) = 2
          ierr = 1
        end if

        if ( rtvz0(m) .ne. rzz0 ) then
          iec = iec + 1
          cepn(iec) = 'z0'
          lepn(iec) = 2
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_r0(m,rzr0,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------
      common /tall68/ itstp(itlmax), itmth(itlmax), rtvr0(itlmax),
     &                rtvx0(itlmax), rtvy0(itlmax), rtvz0(itlmax),
     &                rtvx1(itlmax), rtvy1(itlmax), rtvz1(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( rtvr0(m) .ne. rzr0 ) then
          iec = iec + 1
          cepn(iec) = 'r0'
          lepn(iec) = 2
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_x1y1z1(m,rzx1,rzy1,rzz1,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------
      common /tall68/ itstp(itlmax), itmth(itlmax), rtvr0(itlmax),
     &                rtvx0(itlmax), rtvy0(itlmax), rtvz0(itlmax),
     &                rtvx1(itlmax), rtvy1(itlmax), rtvz1(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( rtvx1(m) .ne. rzx1 ) then
          iec = iec + 1
          cepn(iec) = 'x1'
          lepn(iec) = 2
          ierr = 1
        end if

        if ( rtvy1(m) .ne. rzy1 ) then
          iec = iec + 1
          cepn(iec) = 'y1'
          lepn(iec) = 2
          ierr = 1
        end if

        if ( rtvz1(m) .ne. rzz1 ) then
          iec = iec + 1
          cepn(iec) = 'z1'
          lepn(iec) = 2
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_unit(m,iunt,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itunt(m) .ne. iunt ) then
          iec = iec + 1
          cepn(iec) = 'unit'
          lepn(iec) = 4
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_axis(m,iaxis,inaxi,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

        dimension iaxis(6)

*-----------------------------------------------------------------------

        ierr = 0

        if ( all( itaxs(m,1:inaxi) .ne. iaxis(1:inaxi) ) ) then
          iec = iec + 1
          cepn(iec) = 'axis'
          lepn(iec) = 4
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_output(m,iout,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall52/ itprd(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( ital(m) .eq. 3 ) then !! [t-yield]
          if ( itprd(m) .ne. iout ) ierr = 1
        else
          if ( itout(m) .ne. iout ) ierr = 1
        end if

        if ( ierr .ne. 0 ) then
          iec = iec + 1
          cepn(iec) = 'output'
          lepn(iec) = 6
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_primary(m,ipri,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

      common /tall74/ iprim(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0
        if ( iprim(itlmax) .ne. ipri ) then
          iec = iec + 1
          cepn(iec) = 'primary'
          lepn(iec) = 7
          ierr = 1
        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine check_2dtype(m,idtyp,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( ittwo(m) .ne. idtyp ) then
          iec = iec + 1
          cepn(iec) = '2d-type'
          lepn(iec) = 7
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_factor(m,rfact,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall21/ rtfac(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

        character*20 ctl(2)

*-----------------------------------------------------------------------

        ierr = 0

        write(ctl(1), '(1p1g14.7)') rtfac(m)
        write(ctl(2), '(1p1g14.7)') rfact

        if ( ctl(1) .ne. ctl(2) ) then
          iec = iec + 1
          cepn(iec) = 'factor'
          lepn(iec) = 6
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_letmat(m,letmat,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itlmt(m) .ne. letmat ) then
          iec = iec + 1
          cepn(iec) = 'letmat'
          lepn(iec) = 6
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_letmat1(m,letmat,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itlmt(m) .ne. letmat ) then
          iec = iec + 1
          cepn(iec) = 'letmat1'
          lepn(iec) = 7
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_letmat2(m,letmat,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall55/ itlmt2(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itlmt2(m) .ne. letmat ) then
          iec = iec + 1
          cepn(iec) = 'letmat2'
          lepn(iec) = 7
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_dedxfnc(m,idxfn,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall53/ itdfn(itlmax,2)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itdfn(m,1) .ne. idxfn ) then
          iec = iec + 1
          cepn(iec) = 'dedxfnc'
          lepn(iec) = 7
          ierr = 1
        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine check_dresol(m,dreso,iec,cepn,lepn,ierr)
*     T.Sato 2014/8/19                                                 *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

      common /tall61/ rtdre(itlmax),rtdfa(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( rtdre(m) .ne. dreso ) then
          iec = iec + 1
          cepn(iec) = 'dresol'
          lepn(iec) = 6
          ierr = 1
        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine check_dfano(m,dfan,iec,cepn,lepn,ierr)
*     T.Sato 2014/8/28                                                 *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

      common /tall61/ rtdre(itlmax),rtdfa(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( rtdfa(m) .ne. dfan ) then
          iec = iec + 1
          cepn(iec) = 'dfano'
          lepn(iec) = 5
          ierr = 1
        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine check_deposit(m,idepo,iec,cepn,lepn,ierr)
*     S.Abe 2015/12/03                                                 *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

      common /tall36/ itdpo(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itdpo(m) .ne. idepo ) then
          iec = iec + 1
          cepn(iec) = 'deposit'
          lepn(iec) = 7
          ierr = 1
        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine check_dedxfnc1(m,idxfn,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall53/ itdfn(itlmax,2)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itdfn(m,1) .ne. idxfn ) then
          iec = iec + 1
          cepn(iec) = 'dedxfnc1'
          lepn(iec) = 8
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_dedxfnc2(m,idxfn,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall54/ itdfn2(itlmax,2)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itdfn2(m,1) .ne. idxfn ) then
          iec = iec + 1
          cepn(iec) = 'dedxfnc2'
          lepn(iec) = 8
          ierr = 1
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine check_part(m,iptyp,ipnkf,inpat,iec,cepn,lepn,ierr)
*                                                                      *
*                                                                      *
************************************************************************
        use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05

*-----------------------------------------------------------------------

        implicit double precision (a-h,o-z)

        include 'param.inc'
        include 'param01.inc' ! frtati 2021/10/05

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

        dimension iptyp(mxpart), ipnkf(mxpart) ! frtati 2021/10/05 6 -> mxpart

        common /subtra/ isubt, ipsub(mxpart) ! kitamura22/03/31

*-----------------------------------------------------------------------

        ierr = 0

        if ( itpan(m) .ne. inpat ) then
          iec = iec + 1
          cepn(iec) = 'part'
          lepn(iec) = 4
          ierr = 1
        else if ( inpat .gt. 0 ) then

          if ( any(itpat(m,1:inpat,1) .ne. iptyp(1:inpat)) .or.       ! kitamura22/03/31
     &         any(itpat(m,1:inpat,2) .ne. ipnkf(1:inpat)) .or.       ! kitamura22/03/31
     &         any(itpat(m,1:inpat,3) .ne. ipsub(1:inpat))) then      ! kitamura22/03/31
            iec = iec + 1
            cepn(iec) = 'part'
            lepn(iec) = 4
            ierr = 1
          end if

        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine check_morp(m,morp,iec,cepn,lepn,ierr)
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

      common /tall71/ itmorp(itlmax),itname(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itmorp(m) .ne. morp ) then
          iec = iec + 1
          cepn(iec) = 'morp'
          lepn(iec) = 4
          ierr = 1
        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine check_maxact(m,maxact,iec,cepn,lepn,ierr)
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

      common /tall72/ itactnm(itlmax), itactrg(itlmax), itactmax(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itactnm(m)-1 .ne. maxact ) then
          iec = iec + 1
          cepn(iec) = 'maxact'
          lepn(iec) = 6
          ierr = 1
        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine check_nlatcel(m,nlatcel,iec,cepn,lepn,ierr)
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

      common /tall75/ itnlatcel(itlmax), itnlatmem(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itnlatcel(m) .ne. nlatcel ) then
          iec = iec + 1
          cepn(iec) = 'nlatcel'
          lepn(iec) = 7
          ierr = 1
        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine check_nlatmem(m,nlatmem,iec,cepn,lepn,ierr)
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

      common /tall75/ itnlatcel(itlmax), itnlatmem(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itnlatmem(m) .ne. nlatmem ) then
          iec = iec + 1
          cepn(iec) = 'nlatmem'
          lepn(iec) = 7
          ierr = 1
        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine check_enclos(m,ienclo,iec,cepn,lepn,ierr)
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

      common /tall73/ itenclo(itlmax), itangform(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itenclo(m) .ne. ienclo ) then
          iec = iec + 1
          cepn(iec) = 'enclos'
          lepn(iec) = 7
          ierr = 1
        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine check_angform(m,iangform,iec,cepn,lepn,ierr)
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

      common /tall73/ itenclo(itlmax), itangform(itlmax)

*-----------------------------------------------------------------------

        character*10 cepn(100)
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( itangform(m) .ne. iangform ) then
          iec = iec + 1
          cepn(iec) = 'iangform'
          lepn(iec) = 8
          ierr = 1
        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine judge_tall_check(iec,cepn,lepn,ctfl,itfl,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'err.inc'

*-----------------------------------------------------------------------

        character*10 cepn(100)
        character(itfl) ctfl
        dimension lepn(100)

*-----------------------------------------------------------------------

        ierr = 0

        if ( iec .gt. 0 ) then
          ierr = 1
          ErrCha = ''
          ErrID = 'L:1981/R:judge_tall_check/F:resutl.f' !E83_003_001
          call ErrWrite(ErrID,ErrCha)

          write(*, fmt='(a,a)', advance='no')
     &      'Error: inconsistent tally parameters ',
     &      cepn(1)(1:lepn(1))
          do i = 2, iec
            write(*, fmt='('' & '', a )', advance='no')
     &      cepn(i)(1:lepn(i))
          end do
          write(*, fmt='(a,a)')
     &      ' in ', ctfl
          write(*,'("If you are sure that those parameters are ",
     &    "consistent, please add ireschk=1 in [parameters]")')

        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine open_resfile(m,noe,
     &                        jsn,jsi,dsin,idsi,ill,ilf,
     &                        newtall,ierr)
*                                                                      *
*  subroutine open resfiles *.out, *_err.out.                          *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

        include 'err.inc'

*-----------------------------------------------------------------------

        common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
        character crfln*100

        common /mpi00/ npe, me
        common /talout/ itall
        common /cparm/  maxbch,maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP   THREADPRIVATE(/jcomon/)
        character fname*100, fnume*3

*-----------------------------------------------------------------------

        character dsin(0:9,2)*200
        dimension idsi(0:9,2)
        dimension ill(0:9,2), ilf(0:9,2)
        dimension jsn(2), jsi(2)
        dimension ierrs(2)

*-----------------------------------------------------------------------

        ierr = 0
        newtall = 0

*-----------------------------------------------------------------------
*   open restart file
*-----------------------------------------------------------------------

        do ioe = 1, noe

          jsn(ioe)    = 0
          jsi(ioe)    = 31 + ioe
          idsi(:,ioe) = 0
          ill(:,ioe)  = 1
          ilf(:,ioe)  = 10000000
          ierrs(ioe)  = 0

          jsni = jsn(ioe)
          jsii = jsi(ioe)

          if ( itall.eq.4 ) then
            if ( nobch.lt.1 ) then
              write(fnume,'(i3.3)') 1
            else if ( npe.gt.1 ) then
              write(fnume,'(i3.3)') nobch/(npe-1)+1
            else
              write(fnume,'(i3.3)') nobch+1
            end if
            if ( ioe .eq. 1 ) then
              idsi(jsni,ioe) = irfll(m)+4
              call mk_2dnumfn(crfln(m),dsin(jsni,ioe),irfll(m),
     &              nobch,maxbch,npe)
            else
              idsi(jsni,ioe) = irfll(m)+8
              call mk_2dnumfn(crfln(m),fname,irfll(m),
     &              nobch,maxbch,npe)
              call mk_2derrfn(fname,dsin(jsni,ioe),irfll(m)+4)
            end if
          else

          if ( ioe .eq. 1 ) then
            idsi(jsni,ioe) = irfll(m)
            dsin(jsni,ioe)(1:irfll(m)) = crfln(m)(1:irfll(m))
          else
            idsi(jsni,ioe) = irfll(m) + 4
            call mk_2derrfn(crfln(m),dsin(jsni,ioe),irfll(m))
          end if

          end if

          open(jsii,
     &         file=dsin(jsni,ioe)(1:idsi(jsni,ioe)),
     &         status='old', action='read', form='formatted',
     &         iostat=ierrs(ioe))

        end do

*-----------------------------------------------------------------------
*   check dose restart file exit
*-----------------------------------------------------------------------

        if( all(ierrs(1:noe) .ne. 0) .and. itrff(m) .eq. 0 ) then
          !! it's newly tally
          ierr = 0
          newtall = 1
          goto 900
        end if

        if( noe .eq. 1 .and. ierrs(1) .ne. 0 .and. itrff(m) .ne. 0 )then
          ierr = 1
        write(ErrCha,'(''Error: resfile = '', a, '' does not exist.'')')
     &          dsin(jsn(1),1)(1:idsi(jsn(1),1))
          ErrID = 'L:2112/R:open_resfile/F:resutl.f' !E83_004_001
          call ErrWrite(ErrID,ErrCha)

          write(ErrCha,'(''If you want to create new tally'',
     &              '' in the restart calculation,'',
     &              '' you should not specify resfile parameter.'')')
          ErrID = 'L:2118/R:open_resfile/F:resutl.f'
          call ErrWrite(ErrID,ErrCha)

          goto 900
        end if

        if( noe .eq. 2 .and. any(ierrs(1:noe) .ne. 0 ) ) then
          ierr = 1
          if ( itall.eq.4 ) ierr = 0 ! frtati 2021/03/06
      write(ErrCha,'(''Error: Both '', a,'' and '', a,'' should exist'',
     &              '' when you want to restart the calculation.'')')
     &          dsin(jsn(1),1)(1:idsi(jsn(1),1)),
     &          dsin(jsn(2),2)(1:idsi(jsn(2),2))
          ErrID = 'L:2131/R:open_resfile/F:resutl.f' !E83_005_001
          call ErrWrite(ErrID,ErrCha)
          goto 900
        end if

*-----------------------------------------------------------------------
  900   continue

        return

      end subroutine


************************************************************************
*                                                                      *
      subroutine tregion0(icc,jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                   ntrn,mtrn,ndsm,nvol,ivl,irvl,icfl,
     &                   rglnrf)
*                                                                      *
*       read 'reg =' sub-section of input tally section                *
*       modified by K.Niita on 2011/02/03                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region_mtrg

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension lpar(0:20)

      data klnmax /1000000/


*-----------------------------------------------------------------------

      character rglnrf*200

*-----------------------------------------------------------------------

            ierr = 0

      if( icfl .eq. 1 ) goto 150

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

               rglnrf = chlw(1:i3) // ' '

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

*-----------------------------------------------------------------------
*        error not for reg =
*-----------------------------------------------------------------------

         if( icc .eq. 1 ) then

            if( chcm(i1:i1+3) .ne. 'reg=' ) goto 998

         else if( icc .eq. 2 ) then

            if( chcm(i1:i1+8) .ne. 'reginbox=' ) goto 998

         end if

               ic = inumc(chlw,i1,i3,'=') + 1

*-----------------------------------------------------------------------
*        read region number
*-----------------------------------------------------------------------

               ic = jnumc(chlw,ic,i3)

               ntrn = 0
               mtrn = 0
               ibra = 0
               ilat = 0
               klat = 0
               ifis = 0
               kpar = 0
               iuni = 0

               nvol = 0

               do i = 0, 20
                  ipar(i) = 0
                  jpar(i) = 0
                  lpar(i) = 0
               end do

                  igm = 0
                  call moddas_allocate_int2(0, 20, MAX_NUM_MTRG, mtrg)
               ipmax = 0

            do i = 1, klnmax

                  ic = jnumc(chlw,ic,i3)

                  call tregion3(icn,chlw,i1,i3,ic,
     &                          igm,ipmax,ipar,jpar,
     &                          ibra,ilat,klat,ifis,kpar,iuni
     &                          ,MAX_NUM_MTRG,mtrg)

                     if( icn .gt. 900 ) goto 900
                     if( icn .eq. 500 ) goto 500

               if( i .lt. klnmax .and. ic .gt. i3 ) then

  147             call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) goto 500

                     if( iskip .ne. 0 ) goto 147

                     ifis = ifis + 1

                     ic = i1

               end if

            end do

                  goto 996

  500 continue

*-----------------------------------------------------------------------
*        volume
*-----------------------------------------------------------------------

            if( chcm(i1:i1+5) .eq. 'volume' .or.
     &          chcm(i1:i1+4) .eq. 'value' ) then

                  call tregvol(jsn,jsi,dsin,idsi,ill,ilf,
     &                         jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                         nvol,ivl,irvl)

            end if

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

         return

  900 continue

         if( icn .eq. 980 ) goto 980
         if( icn .eq. 984 ) goto 984
         if( icn .eq. 985 ) goto 985
         if( icn .eq. 986 ) goto 986
         if( icn .eq. 987 ) goto 987
         if( icn .eq. 988 ) goto 988
         if( icn .eq. 989 ) goto 989
         if( icn .eq. 990 ) goto 990
         if( icn .eq. 991 ) goto 991
         if( icn .eq. 992 ) goto 992
         if( icn .eq. 993 ) goto 993
         if( icn .eq. 995 ) goto 995
         if( icn .eq. 997 ) goto 997

*-----------------------------------------------------------------------

  982 continue
         m_err = 'volume of value section is wrong'
         ErrCha = ''
         ErrID = 'L:2334/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  983 continue
         m_err = 'region is too many or memory is lack'
         ErrCha = ''
         ErrID = 'L:2345/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  980 continue

         m_err = 'n1-n2 should be used like {n1-n2}'
         ErrCha = ''
         ErrID = 'L:2357/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  984 continue

         m_err = 'maximum lattice elements by commas is 1000.'
         ErrCha = ''
         ErrID = 'L:2369/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  985 continue

         m_err = 'maximum level of < is 10.'
         ErrCha = ''
         ErrID = 'L:2381/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  986 continue

         m_err = 'usage of u=# is wrong.'
         ErrCha = ''
         ErrID = 'L:2393/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  987 continue

         m_err = '< should be inside ( ).'
         ErrCha = ''
         ErrID = 'L:2405/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  988 continue

         m_err = 'maximun level of ( ) is 10.'
         ErrCha = ''
         ErrID = 'L:2417/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  989 continue

         m_err = 'usage of latice [i1:i2 i3:i4 i5:i6] is wrong.'
         ErrCha = ''
         ErrID = 'L:2429/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  990 continue

         m_err = 'usage of latice [i1 i2 i3] is wrong.'
         ErrCha = ''
         ErrID = 'L:2441/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'all or (all) is available. (all<4), (all 3) are not.'
         ErrCha = ''
         ErrID = 'L:2453/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'Usage of Parenthesis { n1 - n2 } is wrong'
         ErrCha = ''
         ErrID = 'L:2465/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'Usage of Parenthesis ( ) is wrong'
         ErrCha = ''
         ErrID = 'L:2477/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:2490/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Number of region is too large. max(klnmax)=1000000'
         ErrCha = ''
         ErrID = 'L:2502/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of region number is wrong.'
         ErrCha = ''
         ErrID = 'L:2514/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'After { mesh = reg } line should be { reg = }.'
         ErrCha = ''
         ErrID = 'L:2526/R:tregion0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine


************************************************************************
*                                                                      *
      subroutine ttetmesh0(jsn,jsi,dsin,idsi,ill,ilf,
     &     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &     ndim_mtetreg,mtetreg)
*                                                                      *
*       read 'tet =' sub-section of input tally section                *
*       Created by T.Furuta on 2025/01/16                              *
*                                                                      *
************************************************************************
      use moddas

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character(200),intent(out) :: chin, chlw, chcm

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      integer ipar

      data klnmax /1000000/

*-----------------------------------------------------------------------

      integer, intent(in) :: ndim_mtetreg
      integer, intent(out) :: mtetreg(ndim_mtetreg)

*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

*-----------------------------------------------------------------------
*        error not for reg =
*-----------------------------------------------------------------------

               if( chcm(i1:i1+3) .ne. 'reg=' ) goto 998

               ic = inumc(chlw,i1,i3,'=') + 1

*-----------------------------------------------------------------------
*        read region number
*-----------------------------------------------------------------------

               ic = jnumc(chlw,ic,i3)

               ifis = 0
               ipar = 0
               jpar = 0
               kpar = 0

            do i = 1, klnmax

                  ic = jnumc(chlw,ic,i3)

                  call ttetmesh3(icn,chlw,i1,i3,ic,
     &                 ipar,jpar,kpar,ifis,
     &                 ndim_mtetreg,mtetreg)

                     if( icn .gt. 900 ) goto 900
                     if( icn .eq. 500 ) goto 500

               if( i .lt. klnmax .and. ic .gt. i3 ) then

  147             call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) goto 500

                     if( iskip .ne. 0 ) goto 147

                     ifis = ifis + 1

                     ic = i1

               end if

            end do

                  goto 996

*-----------------------------------------------------------------------
*        modify the input
*-----------------------------------------------------------------------

  500    continue

         if(ipar.gt.1)then
          icn=996
          goto 900
         endif

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

         return

  900 continue

         if( icn .eq. 996 ) goto 996
         if( icn .eq. 997 ) goto 997

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Only 1 region is allowed with mesh = tet'
         ErrCha = ''
         ErrID = 'L:2680/R:ttetmesh0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of region number is wrong.'
         ErrCha = ''
         ErrID = 'L:2692/R:ttetmesh0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue
         m_err = 'After { mesh = tet } line should be { reg = }.'
         ErrCha = ''
         ErrID = 'L:2703/R:ttetmesh0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

         return

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine tdepreg0(icc,jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                    ntrn,mtrn,ndsm,nvol,ivl,irvl,icfl,
     &                    rglnrf)
*                                                                      *
*       read 'reg =' sub-section of input tally section                *
*                                                                      *
*       modified subroutine 'tregion' to read                          *
*        "no cell operator ethres" sub-section and                     *
*        "cell cond0 cond1 ..." sub-section                            *
*       for 'reg=weightsum' in [T-Deposit]                             *
*                                                                      *
*       modified by S.Abe on 2016/11/15                                *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*       arguments for read input file:                                 *
*           jsn,jsi,dsin,idsi,ill,ilf,                                 *
*           jpn,chin,chlw,chcm,i1,i2,i3,i4                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*    add output :                                                      *
*          ntrn      total number of cell                              *
*          mtrn      total number of cell information x ... kr(x)      *
*        rglnrf                                                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     local variable :                                                 *
*           ign      start position of das memory                      *
*         ncond      number of condition                               *
*         mcond      counter for number of condition                   *
*          nadd      maximum number of additional condition in input   *
*          madd      counter for number of additional condition        *
*          nrsq      counter for number of reading index               *
*       irsq(5)      array for  reading index                          *
*                                                                      *
*         mcell      counter for number of efficiency list             *
*         nrsq2      counter for number of reading index               *
*        noflag      position of "cell" in index                       *
*      nefflist      number of "list" in index                         *
*       ncount1      counter for number of reading efficiency line     *
*       ncount2      counter for number of reading efficiency value    *
*         ntrn2      number of cell for efficiency list                *
*         mtrn2      number of cell information for efficiency list    *
*                                                                      *
************************************************************************
      use tdepwgtsum_local
      use moddas
      use moddas_region
      use moddas_region_mtrg

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension lpar(0:20)

      data klnmax /1000000/


      dimension irsq(5)

*-----------------------------------------------------------------------

      character rglnrf*200

*-----------------------------------------------------------------------

      ierr = 0

      if( icfl .eq. 1 ) goto 150

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

      call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) return

      if( iskip .ne. 0 ) goto 140

      rglnrf = chlw(1:i3) // ' '

  150 continue

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) return

*-----------------------------------------------------------------------
*     error not for reg =
*-----------------------------------------------------------------------

      if( icc .eq. 1 ) then

         if( chcm(i1:i1+3) .ne. 'reg=' ) goto 998

      else if( icc .eq. 2 ) then

         if( chcm(i1:i1+8) .ne. 'reginbox=' ) goto 998

      end if

      ic = inumc(chlw,i1,i3,'=') + 1

*-----------------------------------------------------------------------
*     read region number for weighted summation
*-----------------------------------------------------------------------

      ic = jnumc(chlw,ic,i3)

      if( chlw(ic:ic+10) .eq. 'weightsum' ) then

         iwgtsum = 1
         igm = 1
         call moddas_allocate_int(
     &           MAX_NUM_REGION_TEMPORARY, idas_region_temporary)

*-----------------------------------------------------------------------
*        initialize
*-----------------------------------------------------------------------

         ntrn = 0
         mtrn = 0
         ntrn1 = 0
         mtrn1 = 0
         ntrn2 = 0
         mtrn2 = 0

         ign = igm

         ncond = -1
         mcond = 0
         nadd = 1
         madd = 1
         nrsq = 0
         do i = 1, 5
            irsq(i) = i
         enddo

         ncell = -1
         mcell = 0
         nrsq2 = 0
         noflag = 0
         nefflist = 0
         ncount1 = 0

*-----------------------------------------------------------------------
*        read one line for conditon section
*-----------------------------------------------------------------------

  100    continue
         call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &              jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
         if( iskip .ne. 0 ) goto 100

*-----------------------------------------------------------------------
*        read "ncond=" value
*-----------------------------------------------------------------------

         if( ncond .lt. 0 ) then

            if( ierr .ne. 0 ) goto 961
            if( jpn  .eq. 3 ) goto 961

            if( chcm(i1:i1+5) .ne. 'ncond=' ) goto 961

            ic = inumc(chlw,i1,i3,'=') + 1
            ic = jnumc(chlw,ic,i3)

            call snum(chlw,ic,i3,ic2,cvvv,ierr)
            if( ierr .ne. 0 ) goto 962

            ncond = nint( cvvv )
            if( ncond .lt. 0 ) goto 962
            goto 100

         endif

*-----------------------------------------------------------------------
*        read index of condition
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

            if( ierr .ne. 0 ) goto 963
            if( jpn  .eq. 3 ) goto 963

            ic = i1

  110       if( ic .gt. i3 ) goto 120

            if(     chlw(ic:ic+1) .eq. 'no' ) then
               irsq(nrsq+1) = 1
               ic = ic + 2
            elseif( chlw(ic:ic+3) .eq. 'cell' ) then
               irsq(nrsq+1) = 2
               ic = ic + 4
            elseif( chlw(ic:ic+7) .eq. 'operator' ) then
               irsq(nrsq+1) = 3
               ic = ic + 8
            elseif( chlw(ic:ic+6) .eq. 'ethres' ) then
               irsq(nrsq+1) = 4
               ic = ic + 6
            elseif( chlw(ic:ic+6) .eq. 'list' ) then
               irsq(nrsq+1) = 5
               ic = ic + 4
            else
               goto 120
            endif

            nrsq = nrsq + 1
            if( nrsq .gt. 5 ) goto 963

            ic = jnumc(chlw,ic,i3)
            goto 110

*-----------------------------------------------------------------------
*        check index of condition
*-----------------------------------------------------------------------

  120       continue

            if( nrsq .gt. 0 ) then

               irno = 0
               ircel = 0
               irope = 0
               ireth = 0
               irlist = 0

               do k = 1, nrsq
                  if( irsq(k) .eq. 1 ) irno   = irno + 1
                  if( irsq(k) .eq. 2 ) ircel  = ircel + 1
                  if( irsq(k) .eq. 3 ) irope  = irope + 1
                  if( irsq(k) .eq. 4 ) ireth  = ireth + 1
                  if( irsq(k) .eq. 5 ) irlist = irlist + 1
               enddo

               if( irno .ne. 1 .or. ircel .ne. 1 .or.
     &             irope .ne. 1 .or. ireth .ne. 1 .or.
     &             irlist .ne. 1 ) goto 963

               call allocate_depwgtsum01(ncond,nadd)

               goto 100

            else

               goto 963

            endif

         endif

*-----------------------------------------------------------------------
*        read condition parameter
*-----------------------------------------------------------------------

         if( ierr .ne. 0 ) goto 964
         if( jpn  .eq. 3 ) goto 964

         ic = i1
         ic2 = i1
         ntrn = 0

         do k = 1, nrsq

            ic = jnumc(chlw,ic2,i3)

            if( irsq(k) .eq. 1 ) then   ! no

               if( chlw(ic:ic+2) .eq. 'and' ) then

                  if( mcond .eq. 0 ) goto 963

                  itmp1 = 1
                  madd = madd + 1
                  if( madd .gt. nadd ) then
                     call reallocate_depwgtsum01(ncond,nadd,madd)
                     nadd = madd
                  endif

                  ic = ic + 3
                  ic2 = ic + 1

               elseif( chlw(ic:ic+4) .eq. 'ncell' ) then

                  if( mcond .ne. ncond ) goto 964
                  goto 201   ! goto next index section

               else

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 964

                  itmp1 = nint( cvvv )
                  if( itmp1 .eq. 0 ) goto 965

                  do im = 1, mcond
                     if( itmp1 .eq. icond(1,im,1) ) goto 966
                  enddo

                  mcond = mcond + 1
                  if( mcond .gt. ncond ) goto 964

                  madd = 1

               endif

            elseif( irsq(k) .eq. 2 ) then   ! cell

               call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ign,ierr
     &                       ,MAX_NUM_REGION_TEMPORARY
     &                       ,idas_region_temporary)
               if( ierr .ne. 0 ) goto 964

               ign = ign + mtrn
               ntrn1 = ntrn1 + ntrn
               mtrn1 = mtrn1 + mtrn

            elseif( irsq(k) .eq. 3 ) then   ! operator

               if(     chlw(ic:ic+1) .eq. 'lt' ) then
                  itmp2 = 1
                  ic = ic + 2
                  ic2 = ic + 1
               elseif( chlw(ic:ic+1) .eq. 'le' ) then
                  itmp2 = 2
                  ic = ic + 2
                  ic2 = ic + 1
               elseif( chlw(ic:ic+1) .eq. 'eq' ) then
                  itmp2 = 3
                  ic = ic + 2
                  ic2 = ic + 1
               elseif( chlw(ic:ic+1) .eq. 'ge' ) then
                  itmp2 = 4
                  ic = ic + 2
                  ic2 = ic + 1
               elseif( chlw(ic:ic+1) .eq. 'gt' ) then
                  itmp2 = 5
                  ic = ic + 2
                  ic2 = ic + 1
               endif

            elseif( irsq(k) .eq. 4 ) then   ! ethres

               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 964

               tmp3 = dble( cvvv )
               if( tmp3 .lt. 0.d0 ) goto 967

            elseif( irsq(k) .eq. 5 ) then   ! list

               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 964

               itmp4 = nint( cvvv )
               if( itmp4 .eq. 0 ) goto 968

            else

               goto 963

            endif

         enddo

         icond(1,mcond,madd) = itmp1
         icond(2,mcond,madd) = itmp2
         icond(3,mcond,madd) = itmp4
         ethres(mcond,madd) = tmp3

         goto 100

*-----------------------------------------------------------------------
*        read one line for efficiency section
*-----------------------------------------------------------------------

  200    continue
         call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &              jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
         if( iskip .ne. 0 ) goto 200

*-----------------------------------------------------------------------
*        read "ncell=" value
*-----------------------------------------------------------------------

  201    continue

         if( ncell .lt. 0 ) then

            if( ierr .ne. 0 ) goto 964
            if( jpn  .eq. 3 ) goto 964

            if( chcm(i1:i1+5) .ne. 'ncell=' ) goto 971

            ic = inumc(chlw,i1,i3,'=') + 1
            ic = jnumc(chlw,ic,i3)

            call snum(chlw,ic,i3,ic2,cvvv,ierr)
            if( ierr .ne. 0 ) goto 972

            ncell = nint( cvvv )
            if( ncell .le. 0 ) goto 972

            goto 200

         endif

*-----------------------------------------------------------------------
*        read index of efficiency
*-----------------------------------------------------------------------

         if( nefflist .eq. 0 ) then

            if( ierr .ne. 0 ) goto 973
            if( jpn  .eq. 3 ) goto 973

            ic = i1

  210       if( ic .gt. i3 ) goto 220

            if(     chlw(ic:ic+3) .eq. 'cell' ) then
               ic = ic + 4
               ic2 = ic + 1
               if( noflag .ne. 0 ) goto 973
               noflag = nrsq2 + 1
            elseif( chlw(ic:ic+3) .eq. 'list' ) then
               ic = ic + 4
               ic2 = ic + 1
               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 973
               ic = ic2 - 1
               nefflist = nefflist + 1
            else
               goto 973
            endif

            nrsq2 = nrsq2 + 1

            ic = jnumc(chlw,ic,i3)
            goto 210

*-----------------------------------------------------------------------

  220       continue

            call allocate_depwgtsum02(ncell,nefflist)

            nefflist = 0
            ic = i1

  230       if( ic .gt. i3 ) goto 240

            if(     chlw(ic:ic+3) .eq. 'cell' ) then
               ic = ic + 4
               ic2 = ic + 1
            elseif( chlw(ic:ic+3) .eq. 'list' ) then
               ic = ic + 4
               ic2 = ic + 1
               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 973
               ic = ic2 - 1
               nefflist = nefflist + 1
               numlist(nefflist) = nint( cvvv )
               do i = 1, nefflist-1
                  if( numlist(nefflist) .eq. numlist(i) ) goto 973
               enddo
            else
               goto 973
            endif

            ic = jnumc(chlw,ic,i3)
            goto 230

*-----------------------------------------------------------------------
*        check correspondence of list number in condition
*-----------------------------------------------------------------------

  240       continue

            if( noflag .eq. 0 ) goto 973
            if( nefflist .le. 0 ) goto 973

            do i = 1, mcond

               nagree = 0

               do j = 1, nefflist
                  if( numlist(j) .eq. 0 ) cycle
                  if( icond(3,i,1) .eq. numlist(j) ) then
                     nagree = 1
                     exit
                  endif
               enddo

               if( nagree .eq. 0 ) goto 974

            enddo

            goto 200

         endif

*-----------------------------------------------------------------------
*        read efficiency parameter
*-----------------------------------------------------------------------

         if( ncount1 .eq. ncell ) goto 300

         if( ierr .ne. 0 ) goto 975
         if( jpn  .eq. 3 ) goto 975

         ic = i1
         ic2 = i1
         ntrn = 0

         ncount1 = ncount1 + 1
         ncount2 = 0

         do k = 1, nefflist + 1

            ic = jnumc(chlw,ic2,i3)

            if( k .eq. noflag ) then

               call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ign,ierr
     &                       ,MAX_NUM_REGION_TEMPORARY
     &                       ,idas_region_temporary)
               if( ierr .ne. 0 ) goto 975

               ign = ign + mtrn
               ntrn2 = ntrn2 + ntrn
               mtrn2 = mtrn2 + mtrn

            else

               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 975

               ncount2 = ncount2 + 1
               depeff(ncount1,ncount2) = cvvv

            endif

         enddo

         goto 200

*-----------------------------------------------------------------------

  300    continue

         ntrn = 0
         mtrn = 0

         call deallocate_depwgtsum
         call moddas_deallocate_int(idas_region_temporary)

*-----------------------------------------------------------------------
*     normal case
*-----------------------------------------------------------------------

      else

         iwgtsum = 0

         ntrn = 0
         mtrn = 0
         ibra = 0
         ilat = 0
         klat = 0
         ifis = 0
         kpar = 0
         iuni = 0

         nvol = 0

         do i = 0, 20
            ipar(i) = 0
            jpar(i) = 0
            lpar(i) = 0
         end do

         igm = 0
         call moddas_allocate_int2(0, 20, MAX_NUM_MTRG, mtrg)
         ipmax = 0

         do i = 1, klnmax

            ic = jnumc(chlw,ic,i3)

            call tregion3(icn,chlw,i1,i3,ic,
     &                    igm,ipmax,ipar,jpar,
     &                    ibra,ilat,klat,ifis,kpar,iuni
     &                    ,MAX_NUM_MTRG,mtrg)
            if( icn .gt. 900 ) goto 900
            if( icn .eq. 500 ) goto 500
            if( jpn .eq. 3 ) return

            if( i .lt. klnmax .and. ic .gt. i3 ) then

  147          call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 500

               if( iskip .ne. 0 ) goto 147

               ifis = ifis + 1

               ic = i1

            end if

         end do

         goto 996

*-----------------------------------------------------------------------
*        modify the input
*-----------------------------------------------------------------------

  500    continue

      endif

*-----------------------------------------------------------------------
*     volume
*-----------------------------------------------------------------------

      if( chcm(i1:i1+5) .eq. 'volume' .or.
     &    chcm(i1:i1+4) .eq. 'value' ) then

         call tregvol(jsn,jsi,dsin,idsi,ill,ilf,
     &                jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                nvol,ivl,irvl)

      end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*     error messages
*-----------------------------------------------------------------------

  950 continue
      m_err = 'Memory error: mmmax exceeds mdas '//
     &        ': Please extend mdas in param.inc'
      ErrCha = ''
      ErrID = 'L:3395/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  951 continue
      m_err = 'region is too many or memory is lack'
      ErrCha = ''
      ErrID = 'L:3403/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr  = 1
      return

  961 continue
      m_err = 'After {reg = weightsum} line should be {ncond =}.'
      ErrCha = ''
      ErrID = 'L:3412/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  962 continue
      m_err = 'Number of condition is wrong.'
      ErrCha = ''
      ErrID = 'L:3420/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  963 continue
      m_err = 'Index of condition is wrong.'
      ErrCha = ''
      ErrID = 'L:3428/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  964 continue
      m_err = 'Input value in condition is wrong.'
      ErrCha = ''
      ErrID = 'L:3436/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  965 continue
      m_err = 'Zero is not allowed for condition number.'
      ErrCha = ''
      ErrID = 'L:3444/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  966 continue
      m_err = 'Condition number is duplicated.'
      ErrCha = ''
      ErrID = 'L:3452/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  967 continue
      m_err = 'Threshold energy is less than zero.'
      ErrCha = ''
      ErrID = 'L:3460/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  968 continue
      m_err = 'Zero is not allowded for list number.'
      ErrCha = ''
      ErrID = 'L:3468/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return

  971 continue
      m_err = 'After condition of weightsum should be {ncell =}.'
      ErrCha = ''
      ErrID = 'L:3477/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  972 continue
      m_err = 'Number of ncell is wrong.'
      ErrCha = ''
      ErrID = 'L:3485/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  973 continue
      m_err = 'Index of ncell is wrong.'
      ErrCha = ''
      ErrID = 'L:3493/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  974 continue
      m_err = 'List number defined in condition does not exist.'
      ErrCha = ''
      ErrID = 'L:3501/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  975 continue
      m_err = 'Input value in ncell is wrong.'
      ErrCha = ''
      ErrID = 'L:3509/R:tdepreg0/F:resutl.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return

*-----------------------------------------------------------------------
*     errors for tregion3
*-----------------------------------------------------------------------

  900 continue

         if( icn .eq. 980 ) goto 980
         if( icn .eq. 984 ) goto 984
         if( icn .eq. 985 ) goto 985
         if( icn .eq. 986 ) goto 986
         if( icn .eq. 987 ) goto 987
         if( icn .eq. 988 ) goto 988
         if( icn .eq. 989 ) goto 989
         if( icn .eq. 990 ) goto 990
         if( icn .eq. 991 ) goto 991
         if( icn .eq. 992 ) goto 992
         if( icn .eq. 993 ) goto 993
         if( icn .eq. 995 ) goto 995
         if( icn .eq. 997 ) goto 997
         if( icn .eq. 999 ) goto 999

*-----------------------------------------------------------------------

  980 continue

         m_err = 'n1-n2 should be used like {n1-n2}'
         ErrCha = ''
         ErrID = 'L:3542/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  984 continue

         m_err = 'maximum lattice elements by commas is 1000.'
         ErrCha = ''
         ErrID = 'L:3554/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  985 continue

         m_err = 'maximum level of < is 10.'
         ErrCha = ''
         ErrID = 'L:3566/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  986 continue

         m_err = 'usage of u=# is wrong.'
         ErrCha = ''
         ErrID = 'L:3578/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  987 continue

         m_err = '< should be inside ( ).'
         ErrCha = ''
         ErrID = 'L:3590/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  988 continue

         m_err = 'maximun level of ( ) is 10.'
         ErrCha = ''
         ErrID = 'L:3602/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  989 continue

         m_err = 'usage of latice [i1:i2 i3:i4 i5:i6] is wrong.'
         ErrCha = ''
         ErrID = 'L:3614/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  990 continue

         m_err = 'usage of latice [i1 i2 i3] is wrong.'
         ErrCha = ''
         ErrID = 'L:3626/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'all or (all) is available. (all<4), (all 3) are not.'
         ErrCha = ''
         ErrID = 'L:3638/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'Usage of Parenthesis { n1 - n2 } is wrong'
         ErrCha = ''
         ErrID = 'L:3650/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'Usage of Parenthesis ( ) is wrong'
         ErrCha = ''
         ErrID = 'L:3662/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:3675/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Number of region is too large. max(klnmax)=1000000'
         ErrCha = ''
         ErrID = 'L:3687/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of region number is wrong.'
         ErrCha = ''
         ErrID = 'L:3699/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'After { mesh = reg } line should be { reg = }.'
         ErrCha = ''
         ErrID = 'L:3711/R:tdepreg0/F:resutl.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue
         l_err = ill(jsn)
         k_err = jsn
         ierr = 1
         return

*-----------------------------------------------------------------------

      end
