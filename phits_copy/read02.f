! Nais_2024 >>>

      subroutine read_dump_parameters(j,isdmp0)
      implicit real*8 (a-h,o-z)

      parameter(icsu0=1)

      include 'param.inc'
      include 'err.inc'

      common /isorsf/ isorf(isrc),lsfile(isrc), sfile(isrc)
      integer :: isorf,lsfile
      character :: sfile*100
      character :: filnm*100,filnmdmp*100

      integer, intent(in) :: j
      integer, intent(out) :: isdmp0(0:30)

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character chins*200, chlws*200, chcms*200
      character ctmp*200
      parameter ( icolms = 200 )

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension lschn0(icsu0), ischn0(icsu0), spava0(icsu0)
      character schan0(icsu0)*8

      logical exex

      data ( schan0(i), i = 1, icsu0 ) /
     &    'dump    '/

      data ( lschn0(i), i = 1, icsu0 ) /  ! character length of parameter
     &     4/

      ischn0(:) = 0
      spava0(:) = 0.0d0

      isdmp0(:) = 0

      jsn = 0

      ill(:) = 0
      ilf(:) = 0
      idsi(:) = 0
      dsin(:) = " "

      filnmdmp = sfile(j)
      inx=index(filnmdmp,"_dmp.out")
      if(inx >= 2) then
        filnm=filnmdmp(1:inx-1)//".out"
      else
        goto 900
      endif

      jpn = 0

      jsi = isorf(j)

      inquire( file = filnm, exist = exex )

      if( exex .eqv. .false. ) then
        write(ErrCha,'("dump parameter file error")')
        ErrID = 'L:70/R:read_dump_parameters/F:read02.f' !E03_015_001
        call ErrWrite(ErrID,ErrCha)
        stop
      endif

      open(jsi, file = filnm,
     &                        form='formatted',status = 'old' )
   10 continue
       read(jsi,'(a)',end=90) chin
       ilf(jsn)=ilf(jsn)+1
       goto 10
   90 rewind(jsi)

!      write(*,*) 'jsi,filnm=',jsi,filnm(1:40)

  100 continue

        call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$[',
     &             jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

         if( ierr .ne. 0 ) return
         if( jpn  .eq. 3 ) goto 900

         if( iskip .ne. 0 ) goto 100

*-----------------------------------------------------------------------
*        identify the parameters
*-----------------------------------------------------------------------

            icl = i1

  160       continue

            chlc = chlw
            call chcomp(chlc,icl,i3,i5)

         do i = 1, icsu0

            il = icl + lschn0(i) - 1

            if( chlc(icl:il) .eq. schan0(i)(1:lschn0(i)) ) goto 200

         end do

         goto 100

*-----------------------------------------------------------------------
*        read value of parameters
*-----------------------------------------------------------------------

  200    continue


      ic = inumc(chlw,icl+4,i3,'=') + 1

      icl = inumc(chlw,ic,i3,';') - 1

      call onum(chlw,ic,icl,cvvv,ierr)

      if( ierr .ne. 0 ) then
        write(ErrCha,'("dump error")')
        ErrID = 'L:131/R:read_dump_parameters/F:read02.f' !E03_015_001
        call ErrWrite(ErrID,ErrCha)
        goto 900
      endif

      isdmp0(0) = nint( cvvv )

  300 continue
      call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &          jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) then
        write(ErrCha,'("dump error")')
        ErrID = 'L:145/R:read_dump_parameters/F:read02.f' !E03_015_003
        call ErrWrite(ErrID,ErrCha)
        goto 900
      endif

      if( iskip .ne. 0 ) goto 300

      ic = i1

      do k = 1, abs( isdmp0(0) )

        if( ic .gt. i3 ) then

  320     continue
          call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

          if( ierr .ne. 0 ) return
          if( jpn  .eq. 3 ) then
            write(ErrCha,'("dump error")')
            ErrID = 'L:165/R:read_dump_parameters/F:read02.f' !E03_015_004
            call ErrWrite(ErrID,ErrCha)
            goto 900
          endif

          if( iskip .ne. 0 ) goto 320

          ic = i1

        end if

        ic = jnumc(chlw,ic,i3)

        call snum(chlw,ic,i3,ic2,cvvv,ierr)

        if( ierr .ne. 0 ) then
          write(ErrCha,'("dump error")')
          ErrID = 'L:182/R:read_dump_parameters/F:read02.f' !E03_015_005
          call ErrWrite(ErrID,ErrCha)
          goto 900
        endif

        isdmp0(k) = nint( cvvv )

        if( isdmp0(k) .gt. 20 .or.
     &      isdmp0(k) .le.  0 ) then
          write(ErrCha,'("dump error")')
          ErrID = 'L:192/R:read_dump_parameters/F:read02.f' !E03_015_006
          call ErrWrite(ErrID,ErrCha)
          goto 900
        endif

        ic = ic2

      end do

  900 continue

      close(jsi)


      if(isdmp0(0) == 0) then
        write(*,'(1x,a,i5)') 'tally=',j
        write(*,'(1x,a,i5)') 'dump =',isdmp0(0)
        write(ErrCha,'("dump error")')
        ErrID = 'L:210/R:read_dump_parameters/F:read02.f' !E03_015_007
        call ErrWrite(ErrID,ErrCha)
        stop
      else
        do i = 1, abs( isdmp0(0) )
           if(isdmp0(i).eq.0) then
              write(ErrCha,'("dump error")')
              ErrID = 'L:217/R:read_dump_parameters/F:read02.f' !E03_015_008
              call ErrWrite(ErrID,ErrCha)
             stop
           endif
        end do

      endif

      return
      end subroutine read_dump_parameters
! Nais_2024

************************************************************************
*                                                                      *
      subroutine ndatmax(jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [datamax] section of input files                          *
*       modified by K.Niita on 2016/07/26                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

c frtati 2021/12/17 added [data max] photo-nuclear, deuteron, alpha
      common /nntmax/ dmxdxx(6,500),indmm,ipdmm(6,6),ipdnn(6),ipdpt(6),
     &                nucdxx(6,500),matdxx(6,500)
      common /ndtmax/ indmp, indmn, indmu, indmd, indma,
     & nucdxp(500), nucdxn(500), nucdxu(500), nucdxd(500), nucdxa(500),
     & matdxp(500), matdxn(500), matdxu(500), matdxd(500), matdxa(500)
      common /ddtmax/ dmxdxp(500), dmxdxn(500),
     &                dmxdxu(500), dmxdxd(500), dmxdxa(500)

      dimension dmxdmm(500)
      dimension nucdpn(500)
      dimension matdpn(500)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 8 * 0 /

      dimension iptyp(20)
      dimension jstyp(6), jnkf0(6)

      common /subtra/ isubt, ipsub(mxpart)   ! kitamura22/03/31

*-----------------------------------------------------------------------

      character element(104)*3,cnuc*3

      data element/
     & 'h  ','he ','li ','be ','b  ','c  ','n  ','o  ',
     & 'f  ','ne ','na ','mg ','al ','si ','p  ','s  ',
     & 'cl ','ar ','k  ','ca ','sc ','ti ','v  ','cr ',
     & 'mn ','fe ','co ','ni ','cu ','zn ','ga ','ge ',
     & 'as ','se ','br ','kr ','rb ','sr ','y  ','zr ',
     & 'nb ','mo ','tc ','ru ','rh ','pd ','ag ','cd ',
     & 'in ','sn ','sb ','te ','i  ','xe ','cs ','ba ',
     & 'la ','ce ','pr ','nd ','pm ','sm ','eu ','gd ',
     & 'tb ','dy ','ho ','er ','tm ','yb ','lu ','hf ',
     & 'ta ','w  ','re ','os ','ir ','pt ','au ','hg ',
     & 'tl ','pb ','bi ','po ','at ','rn ','fr ','ra ',
     & 'ac ','th ','pa ','u  ','np ','pu ','am ','cm ',
     & 'bk ','cf ','es ','fm ','md ','no ','lr ','ku '/

      logical deqn1
      logical dcom2

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0
            inpat = 0
            nndmx = 0

            indmm = indmm + 1
            if( indmm .gt. 6 ) goto 983

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               ic  = i1
               icl = i1

               chlc = chlw
               call chcomp(chlc,icl,i3,i5)

            if(chlc(icl:icl+3) .eq. 'part' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = i3

  400          continue

                  ic = jnumc(chlw,ic,i3)

                  if( ic .gt. i3 ) goto 140

*-----------------------------------------------------------------------

                  call rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

                     if( ierr .eq. 994 ) goto 993
                     if( ierr .eq. 998 ) goto 993
                     if( istyp .le. 0 )  goto 993
                     if( isubt .eq. 1 )  goto 993   ! kitamura22/03/31

                        inpat = inpat + 1
                        if( inpat .gt. 2 ) goto 992

                        iptyp(inpat) = istyp

                     goto 400

            end if

         end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+6) .eq. 'nucleus' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 7

            else if( chlw(ic:ic+3) .eq. 'dmax' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 4

            else if( chlw(ic:ic+2) .eq. 'mat' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 3

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 3 ) goto 997

                  innuc = 0
                  indnn = 0
                  innon = 0
                  inmat = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) innuc = innuc + 1
                  if( imsq(k) .eq. 2 ) indnn = indnn + 1
                  if( imsq(k) .eq. 3 ) inmat = inmat + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( innuc .ne. 1 .or. indnn .ne. 1 .or.
     &                inmat .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  innuc = 1
                  indnn = 1
                  inmat = 1
                  innon = 0

                  nrsq = 3

            end if

         end if

*-----------------------------------------------------------------------

            if( inpat .gt. 0 ) then

               do k = 1, inpat

cfrtati 2021/12/17 added iptyp=14,15,18
                  if( iptyp(k) .ne. 1 .and. iptyp(k) .ne. 2  .and.
     &                iptyp(k) .ne.14 .and. iptyp(k) .ne.15  .and.
     &                iptyp(k) .ne.18 ) goto 981

               end do

            end if

            if( inpat .eq. 0 ) goto 980

*-----------------------------------------------------------------------
*        read informations : non
*-----------------------------------------------------------------------

               nndmx = nndmx + 1

               if( nndmx .gt. 500 ) goto 998

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

*-----------------------------------------------------------------------
*        read datamax informations
*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 2 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  dmxdmm(nndmx) = cvvv

*-----------------------------------------------------------------------
*        read mat informations
*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 3 ) then

               if( chlw(ic:ic+2) .eq. 'all' ) then

                     matdpn(nndmx) = 0

                     ic2 = ic + 3

               else

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                      matdpn(nndmx) = nint( cvvv )

               end if

*-----------------------------------------------------------------------
*        nucleus
*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 1 ) then

               if( chlw(ic:ic+2) .eq. 'all' ) then

                     nucdpn(nndmx) = 0

                     ic2 = ic + 3

               else

                        isa = 0
                        isn = 0

                        ica = 0
                        icb = 0
                        icm = 0
                        icn = 0

                     do i = ic, i3

                        if( chlw(i:i) .ge. 'a' .and.
     &                      chlw(i:i) .le. 'z' ) then

                           isa = isa + 1

                           if( isa .eq. 1 ) ica = i

                           icb = i

                        else if( deqn1( chlw(i:i) ) ) then

                           isn = isn + 1

                           if( isn .eq. 1 ) icm = i

                           icn = i

                        else if( dcom2( chlw(i:i) ) .or.
     &                           i .eq. i3 ) then

                           icd = i
                           goto 502

                        else

                           goto 982

                        end if

                     end do

                        icd = i3

  502                continue

                     if( ica .eq. 0 .or. icb .eq. 0 ) goto 983
                     if( icb-ica .lt. 0 .or. icb-ica .gt. 1 ) goto 982

                     cnuc = chlw(ica:icb)//'  '

                     do j = 1, 104

                        if( cnuc(1:3) .eq. element(j)(1:3) ) then

                           icha = j

                           goto 452

                        end if

                     end do

                           goto 982

  452                continue

                     if( icha .gt. 104 )  goto 982

                  if( icm .eq. 0 .or. icn .eq. 0 ) then

                     nucdpn(nndmx) = icha * 1000

                  else

                     if( isn .gt. 3 ) goto 982

                     read(chlw(icm:icn),'(i5)') masi

                     if( masi .lt. icha ) goto 982
                     if( masi-icha .gt. maxnt ) goto 982

                     nucdpn(nndmx) = icha * 1000 + masi

                  end if

                     ic2 = icd + 1

               end if

            end if

*-----------------------------------------------------------------------

         end do

         goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

            do i = 1, nndmx
            do j = 1, inpat

               if( iptyp(j) .eq. 1 ) then

                  indmp = indmp + 1
                  nucdxp(indmp) = nucdpn(i)
                  matdxp(indmp) = matdpn(i)
                  dmxdxp(indmp) = dmxdmm(i)

               else if( iptyp(j) .eq. 2 ) then

                  indmn = indmn + 1
                  nucdxn(indmn) = nucdpn(i)
                  matdxn(indmn) = matdpn(i)
                  dmxdxn(indmn) = dmxdmm(i)
               else if ( iptyp(j) .eq. 14 ) then

                  indmu = indmu + 1
                  nucdxu(indmu) = nucdpn(i)
                  matdxu(indmu) = matdpn(i)
                  dmxdxu(indmu) = dmxdmm(i)

               else if ( iptyp(j) .eq. 15 ) then

                  indmd = indmd + 1
                  nucdxd(indmd) = nucdpn(i)
                  matdxd(indmd) = matdpn(i)
                  dmxdxd(indmd) = dmxdmm(i)

               else if ( iptyp(j) .eq. 18 ) then

                  indma = indma + 1
                  nucdxa(indma) = nucdpn(i)
                  matdxa(indma) = matdpn(i)
                  dmxdxa(indma) = dmxdmm(i)

               end if

            end do
            end do

                  ipdnn(indmm) = nndmx
                  ipdpt(indmm) = inpat ! frtati 2021/12/17

            do i = 1, inpat
              ipdmm(indmm,i) = iptyp(i)
            end do

            do i = 1, nndmx

                  nucdxx(indmm,i) = nucdpn(i)
                  matdxx(indmm,i) = matdpn(i)
                  dmxdxx(indmm,i) = dmxdmm(i)

            end do

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  980 continue

         m_err = 'Particle is not defined in [datamax]'
         ErrCha = ''
         ErrID = 'L:740/R:ndatmax/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  982 continue

         m_err = 'Description of Nucleus is wrong'
         ErrCha = ''
         ErrID = 'L:752/R:ndatmax/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  983 continue

         m_err = 'Max [data max] section is 6'
         ErrCha = ''
         ErrID = 'L:764/R:ndatmax/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  981 continue

         m_err = 'Particle type is not available'
         ErrCha = ''
         ErrID = 'L:776/R:ndatmax/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'Number of particles is too large (<3)'
         ErrCha = ''
         ErrID = 'L:788/R:ndatmax/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'Name of particle is wrong'
         ErrCha = ''
         ErrID = 'L:800/R:ndatmax/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'Description of part = is wrong'
         ErrCha = ''
         ErrID = 'L:812/R:ndatmax/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[datamax] is wrong.'
         ErrCha = ''
         ErrID = 'L:825/R:ndatmax/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Number of datamax '//
     &           'exceeds 500'
         ErrCha = ''
         ErrID = 'L:838/R:ndatmax/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [datanax] is wrong.'
         ErrCha = ''
         ErrID = 'L:850/R:ndatmax/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine multip(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [multiplier] section of input files                       *
*       modified by K.Niita on 2010/01/31                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_multiplier

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc' ! frtati 2021/10/05
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /multipl/ imltp,iimlt(multmax),inmlt(multmax),
C MATSUDA 2024.11.25 (lagrange: ilmlt)
     &  ilmlt(multmax),
     &  idmlt(multmax),ismlt(multmax),impan(multmax),impat(multmax,6,3),
     &                 jmpat(multmax,6,6,2), imdfl(multmax)
C MATSUDA 2024.11.25 (multplf: file name)
      common /multplf/ imltf(multmax),lmltfile(multmax),mltfile(multmax)
      character mltfile*100

      dimension iptyp(mxpart), ipnkf(mxpart)
      dimension imtyp(mxpart,mxpart), imnkf(mxpart,mxpart)
      dimension jstyp(mxpart), jnkf0(mxpart)

      common /subtra/ isubt, ipsub(mxpart)   ! kitamura22/03/31

C MATSUDA 2024.12.09 (multpl27,35: x-, y-txt, and epsout)
      common /multpl27/ impxl(multmax),impxt(multmax),
     &                  impyl(multmax),impyt(multmax)
      character impxt*200, impyt*200
      common /multpl35/ impeps(multmax)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)
      dimension ill(0:9), ilf(0:9)

C MATSUDA 2024.12.09 (multpl27: x-, y-txt)
      character cxtxt*200
      character cytxt*200

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0
            nnum  = 0
            ninp  = -1
             icmp = 0
            inpat = 0

            imltp = imltp + 1

            if( imltp .gt. multmax ) goto 991

C MATSUDA 2024.11.25 (lagrange: ilmlt, default value)
            ilmlt(imltp) = 2
            idmlt(imltp) = 0
            imdfl(imltp) = 0
            imltf(imltp) = 0
            impxl(imltp) = 0
            impyl(imltp) = 0
            impeps(imltp) = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        number and interpolation
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               ic  = i1
               icl = i1

               chlc = chlw
               call chcomp(chlc,icl,i3,i5)

*-----------------------------------------------------------------------

            if( chlc(icl:icl+5) .eq. 'number' ) then

               ic = inumc(chlw,icl+6,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = inumc(chlw,ic,i3,';') - 1

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               idmlt(imltp) = nint( cvvv )

               if( idmlt(imltp) .lt. -299 .or.
     &             idmlt(imltp) .gt. -200 ) goto 992

                  nnum = 1

*-----------------------------------------------------------------------

            else if( chlc(icl:icl+12) .eq. 'interpolation' ) then

               icl = icl + 13
               if( chlc(icl:icl) .ne. '=' ) goto 994
               icl = icl + 1

               if( chlc(icl:icl+2) .eq. 'log' ) then

                  ninp = -1

               else if( chlc(icl:icl+2) .eq. 'lin' ) then

                  ninp = 1

               else if( chlc(icl:icl+3) .eq. 'glow' ) then

                  ninp = -2

               else if( chlc(icl:icl+4) .eq. 'ghigh' ) then

                  ninp = 2

C MATSUDA 2024.11.25 (options: xlin, ylog, xlog, and ylin)
               else if( chlc(icl:icl+4) .eq. 'xlin' ) then

                  ninp = 3

               else if( chlc(icl:icl+4) .eq. 'ylog' ) then

                  ninp = 4

               else if( chlc(icl:icl+4) .eq. 'xlog' ) then

                  ninp = -3

               else if( chlc(icl:icl+4) .eq. 'ylin' ) then

                  ninp = -4

               else

                  goto 994

               end if

*-----------------------------------------------------------------------
C MATSUDA 2024.11.25 (lagrange)

            else if( chlc(icl:icl+7) .eq. 'lagrange' ) then

               ic = inumc(chlw,icl+2,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = inumc(chlw,ic,i3,';') - 1

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               ilmlt(imltp) = nint( cvvv )

               if( ilmlt(imltp) .lt. 2 .and.
     &             ilmlt(imltp) .gt. 4 ) goto 994

*-----------------------------------------------------------------------
C MATSUDA 2024.12.09 (x-txt)

            else if( chlc(icl:icl+4) .eq. 'x-txt' ) then

               ic = inumc(chlw,icl+2,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 994

               icl = inumc(chlw,ic,i3,';') - 1

               ict = min( ic + 199, i2 )

               impxt(imltp) = chin(ic:ict)

               impxl(imltp) = ict - ic + 1

*-----------------------------------------------------------------------
C MATSUDA 2024.12.09 (y-txt)

            else if( chlc(icl:icl+4) .eq. 'y-txt' ) then

               ic = inumc(chlw,icl+2,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 994

               icl = inumc(chlw,ic,i3,';') - 1

               ict = min( ic + 199, i2 )

               impyt(imltp) = chin(ic:ict)

               impyl(imltp) = ict - ic + 1

*-----------------------------------------------------------------------
C MATSUDA 2024.12.09 (epsout)

            else if( chlc(icl:icl+5) .eq. 'epsout' ) then

               ic = inumc(chlw,icl+2,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = inumc(chlw,ic,i3,';') - 1

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               impeps(imltp) = nint( cvvv )

*-----------------------------------------------------------------------

            else if( chlc(icl:icl+1) .eq. 'ne' ) then

               ic = inumc(chlw,icl+2,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = inumc(chlw,ic,i3,';') - 1

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               inmlt(imltp) = nint( cvvv )

               if( inmlt(imltp) .le. 0 ) goto 994

                  nrsq = inmlt(imltp)

*-----------------------------------------------------------------------
*           particle name
*-----------------------------------------------------------------------

            else if( chlc(icl:icl+3) .eq. 'part' ) then

               ic = inumc(chlw,icl+4,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = i3

  400       continue

               ic = jnumc(chlw,ic,icl)

               if( ic .gt. icl ) then

                  icl = jnumc(chlw,icl+2,i3)


                  goto 140

               end if

*-----------------------------------------------------------------------

            call rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

               if( ierr .eq. 994 ) goto 990
               if( ierr .eq. 998 ) goto 990
               if( isubt .eq. 1 )  goto 990   ! kitamura22/03/31

*-----------------------------------------------------------------------

               inpat = inpat + 1

               if( inpat .gt. 6 ) goto 989

                  iptyp(inpat) = istyp
                  ipnkf(inpat) = inkf0

               if( istyp .lt. 0 ) then

                  do i = 1, -istyp

                     imtyp(inpat,i) = jstyp(i)
                     imnkf(inpat,i) = jnkf0(i)

                  end do

               end if

               goto 400

C MATSUDA 2024.11.25 (lagrange: file)
*-----------------------------------------------------------------------
*           input file name
*-----------------------------------------------------------------------

            else if( chlc(icl:icl+3) .eq. 'file' ) then

               ic = inumc(chlw,icl+4,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = i3

               ic1 = jnumc(chlw,ic,icl)
               ic2 = min( i3, inumc(chlw,ic+1,icl,' ') )

               lmltfile(imltp) = ic2 - ic1 + 1

               do i = 1, lmltfile(imltp)

                  mltfile(imltp)(i:i) = chin(ic1+i-1:ic1+i-1)

               end do

               do i = lmltfile(imltp) + 1, 100

                  mltfile(imltp)(i:i) = ' '

               end do

               imltf(imltp) = 1

*-----------------------------------------------------------------------

            else

               goto 994

            end if

               goto 140

         end if

*-----------------------------------------------------------------------
*        check multiplier parameters
*-----------------------------------------------------------------------

         if( nrsq .gt. 0 ) then

               iimlt(imltp) = ninp

               if( idmlt(imltp) .eq. 0 ) goto 993

*-----------------------------------------------------------------------
*           summary of particle
*-----------------------------------------------------------------------

            if( inpat .gt. 0 ) then

                  impan(imltp) = inpat

               do i = 1, inpat

                  impat(imltp,i,1) = iptyp(i)
                  impat(imltp,i,2) = ipnkf(i)
                  impat(imltp,i,3) = 0  ! T.Sato 2024/12/04, always 0 for multiplier part (has meaning for tally part)

                  if( impat(imltp,i,1) .lt. 0 ) then

                     do k = 1, -impat(imltp,i,1)

                        jmpat(imltp,i,k,1) = imtyp(i,k)
                        jmpat(imltp,i,k,2) = imnkf(i,k)

                     end do

                  end if

               end do

            else if( inpat .eq. 0 ) then

                  impan(imltp) = 1
                  impat(imltp,1,1) = 20
                  impat(imltp,1,2) = 0
                  impat(imltp,1,3) = 0  ! T.Sato 2024/12/04, always 0 for multiplier part (has meaning for tally part)

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        read multiplier informations
*-----------------------------------------------------------------------

               igm = ismlt(imltp)
               igr = nrsq
               igg = 2 * igr
               call moddas_reallocate_dbl(
     &                 multmax, imltp, igg, ismlt, gmsh_ismlt)

                  ic = i1

            do i = 1, igg

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 995

                  gmsh_ismlt(igm-1+i) = cvvv

C MATSUDA 2024.11.25 (Error 988)
                  if( mod(i,2) .eq. 1 .and. i .gt. 2 ) then
                     if( cvvv-gmsh_ismlt(igm-3+i) .le. 0.d0 ) goto 988
                  end if

               ic = ic2

               if( i .lt. igg .and. ic .gt. i3 ) then

  149             call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                  if( ierr .ne. 0 ) goto 995
                  if( jpn  .eq. 3 ) goto 995

                  if( iskip .ne. 0 ) goto 149

                  ic = i1

               end if

            end do

                  icmp = 1
                  nrsq = 0

               goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

               if( icmp .eq. 0 ) goto 995
C MATSUDA 2024.11.25 (Error 987)
               if( ilmlt(imltp) .gt.  2 .and.
     &             inmlt(imltp) .lt.  4 ) goto 987
C MATSUDA 2024.12.09 (Error 986)
!               igr = inmlt(imltp)
!               igm = ismlt(imltp) - 1
!               if( ( iimlt(imltp) .eq. -1 .or.
!     &               iimlt(imltp) .eq. -3 .or.
!     &               iimlt(imltp) .eq. -4 ) .and.
!     &             gmsh_ismlt(igm+1) .le. 0.d0 ) goto 986
!
!               if(   iimlt(imltp) .eq. -1 .or.
!     &               iimlt(imltp) .eq.  3 .or.
!     &               iimlt(imltp) .eq.  4 ) then
!                  do i = 1, igr
!                     if( gmsh_ismlt(igm+2*i) .le. 0.d0 ) goto 986
!                  end do
!               end if

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:1374/R:multip/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

C MATSUDA 2024.11.25 (Error 987, 988)
*-----------------------------------------------------------------------

  987 continue

         m_err = 'Lagrange interpolation requires 3 or more data'
         ErrCha = ''
         ErrID = 'L:1387/R:multip/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  988 continue

         m_err = 'The data is not in ascending oder'
         ErrCha = ''
         ErrID = 'L:1399/R:multip/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  989 continue

         m_err = 'Number of particles is too large (<6)'
         ErrCha = ''
         ErrID = 'L:1411/R:multip/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  990 continue

         m_err = 'Name of particle is wrong'
         ErrCha = ''
         ErrID = 'L:1423/R:multip/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'Number of [multiplier] section is too large, '//
     &           ' Please extend multmax in param.inc' ! T.Sato 2018/12/08
         ErrCha = ''
         ErrID = 'L:1436/R:multip/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'ID number should be -200 - -299'
         ErrCha = ''
         ErrID = 'L:1448/R:multip/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'ID number should be defined -200 - -299'
         ErrCha = ''
         ErrID = 'L:1460/R:multip/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'Description of parameters is wrong'
         ErrCha = ''
         ErrID = 'L:1472/R:multip/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Definition of data sequences in '//
     &           '[multiplier] is wrong.'
         ErrCha = ''
         ErrID = 'L:1485/R:multip/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine split(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [splitting] section of input files                        *
*       modified by K.Niita on 2004/12/05                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /splreg/ isptn, npreg(6), mnspt(6,0:20),
     &                ipgrc(6), ipgrt(6), ksplt(6), isplt(6),
     &                ispct(6,9), ispem(6,2)
      common /splrge/ espem(6,2)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 3, 7 * 0 /

      character dkam*6

      dimension s_fc(kvlmax)

      dimension iptyp(20)
      dimension jstyp(6), jnkf0(6)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0
            inpat = 0

            isptn = isptn + 1
            if( isptn .gt. 6 ) goto 991

            npreg(isptn) = 0

            ispct(isptn,1) = 0
            ispct(isptn,2) = 0
            ispct(isptn,3) = 0
            ispct(isptn,4) = -9999
            ispct(isptn,5) =  9999
            ispct(isptn,6) = -9999
            ispct(isptn,7) =  9999
            ispct(isptn,8) = -9999
            ispct(isptn,9) =  9999

            ispem(isptn,1) = 0
            ispem(isptn,2) = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               ic  = i1
               icl = i1

               chlc = chlw
               call chcomp(chlc,icl,i3,i5)

*-----------------------------------------------------------------------

            if(chlc(icl:icl+7) .eq. 'ctmin(1)' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               if( ic .gt. i3 ) goto 994
               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               ispct(isptn,1) = 1
               ispct(isptn,4) = nint( cvvv )

               goto 140

*-----------------------------------------------------------------------

            else if(chlc(icl:icl+7) .eq. 'ctmax(1)' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               if( ic .gt. i3 ) goto 994
               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               ispct(isptn,1) = 1
               ispct(isptn,5) = nint( cvvv )

               goto 140

*-----------------------------------------------------------------------

            else if(chlc(icl:icl+7) .eq. 'ctmin(2)' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               if( ic .gt. i3 ) goto 994
               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               ispct(isptn,2) = 1
               ispct(isptn,6) = nint( cvvv )

               goto 140

*-----------------------------------------------------------------------

            else if(chlc(icl:icl+7) .eq. 'ctmax(2)' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               if( ic .gt. i3 ) goto 994
               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               ispct(isptn,2) = 1
               ispct(isptn,7) = nint( cvvv )

               goto 140

*-----------------------------------------------------------------------

            else if(chlc(icl:icl+7) .eq. 'ctmin(3)' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               if( ic .gt. i3 ) goto 994
               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               ispct(isptn,3) = 1
               ispct(isptn,8) = nint( cvvv )

               goto 140

*-----------------------------------------------------------------------

            else if(chlc(icl:icl+7) .eq. 'ctmax(3)' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               if( ic .gt. i3 ) goto 994
               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               ispct(isptn,3) = 1
               ispct(isptn,9) = nint( cvvv )

               goto 140

*-----------------------------------------------------------------------

            else if(chlc(icl:icl+3) .eq. 'emin' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               if( ic .gt. i3 ) goto 994
               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               ispem(isptn,1) = 1
               espem(isptn,1) = cvvv

               goto 140

*-----------------------------------------------------------------------

            else if(chlc(icl:icl+3) .eq. 'emax' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               if( ic .gt. i3 ) goto 994
               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               ispem(isptn,2) = 1
               espem(isptn,2) = cvvv

               goto 140

*-----------------------------------------------------------------------

            else if(chlc(icl:icl+3) .eq. 'part' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = i3

  400          continue

                  ic = jnumc(chlw,ic,i3)

                  if( ic .gt. i3 ) goto 140

*-----------------------------------------------------------------------

                  call rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

                     if( ierr .eq. 994 ) goto 993
                     if( ierr .eq. 998 ) goto 993

                  if( istyp .gt. 0 ) then

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = istyp

                  else if( istyp .lt. 0 ) then

                     do i = 1, -istyp

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = jstyp(i)

                     end do

                  end if

                  goto 400

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'r-in' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 4

            else if( chlw(ic:ic+4) .eq. 'r-out' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 5

            else if( chlw(ic:ic+5) .eq. 'factor' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 6

            else

               goto 997

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 3 ) goto 997

                  irin = 0
                  irot = 0
                  imfc = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) irin = irin + 1
                  if( imsq(k) .eq. 2 ) irot = irot + 1
                  if( imsq(k) .eq. 3 ) imfc = imfc + 1
               end do

                  if( irin .ne. 1 .or. irot .ne. 1 .or.
     &                imfc .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  irin = 1
                  irot = 1
                  imfc = 1

                  nrsq = 3

            end if

         end if

*-----------------------------------------------------------------------

               isplt(isptn) = 3

            do k = 1, nrsq

               if( isplt(isptn) .eq. 3 .and. imsq(k) .eq. 1 )
     &             isplt(isptn) = 0
               if( isplt(isptn) .eq. 3 .and. imsq(k) .eq. 2 )
     &             isplt(isptn) = 1

            end do

*-----------------------------------------------------------------------

            if( inpat .gt. 0 ) then

               do k = 1, inpat

                  if( iptyp(k) .eq. 20 ) inpat = 0

               end do

            end if

            if( inpat .eq. 0 ) then

               do k = 1, 19

                  iptyp(k) = k

               end do

                  inpat = 19

            end if

*-----------------------------------------------------------------------
*        read splitting informations
*-----------------------------------------------------------------------

               npreg(isptn) = npreg(isptn) + 1

               if( npreg(isptn) .gt. kvlmax ) goto 998

            if( npreg(isptn) .eq. 1 ) then

                  call moddas_reallocate_int(
     &                    6, isptn, MAX_NUM_IPGRC, ipgrc, idas_ipgrc)

                  idsm  = ipgrc(isptn)
                  jdsm  = -1

            end if

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 1 .or. imsq(k) .eq. 2 ) then

                     ndsm = idsm + jdsm + 3

                  call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ndsm,ierr
     &                          ,MAX_NUM_IPGRC,idas_ipgrc)

                  if( ierr .ne. 0 ) goto 999

                     jdsm = jdsm + 1
                     idas_ipgrc(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                     idas_ipgrc(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn

            else if( imsq(k) .eq. 3 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               s_fc(npreg(isptn)) = cvvv

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

            if( npreg(isptn) .eq. 0 ) then

               isptn = isptn - 1
               return

            end if

               if( jdsm > MAX_NUM_IPGRC ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.split@read02.f ?dimension over idas_ipgrc?'
     &                    //' jdsm > MAX_NUM_IPGRC'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_IPGRC@moddas.f=',MAX_NUM_IPGRC,')'
                  ErrID = 'L:2033/R:split/F:read02.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 6, isptn, jdsm+1, ipgrc, idas_ipgrc)

               mnspt(isptn, 0) = npreg(isptn)
               mnspt(isptn,20) = inpat

            do j = 1, inpat

               mnspt(isptn,j) = iptyp(j)

            end do

               call moddas_reallocate_dbl(
     &                 6, isptn, npreg(isptn)*1+1, ksplt, das_ksplt)

               idsm  = ksplt(isptn)
               jdsm  = 0


            do i = 1, npreg(isptn)

               jdsm = jdsm + 1
               das_ksplt(idsm+jdsm) = s_fc(i)

            end do

            if( isplt(isptn) .ne. 0 ) npreg(isptn) = - npreg(isptn)

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:2076/R:split/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'Number of [splitting] section is'//
     &           ' too large (=<6)'
         ErrCha = ''
         ErrID = 'L:2089/R:split/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'Number of particles is too large (<19)'
         ErrCha = ''
         ErrID = 'L:2101/R:split/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'Name of particle is wrong'
         ErrCha = ''
         ErrID = 'L:2113/R:split/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'Description of part, ctmin, ctmax = is wrong'
         ErrCha = ''
         ErrID = 'L:2125/R:split/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[splitting] is wrong.'
         ErrCha = ''
         ErrID = 'L:2138/R:split/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of region for splitting '//
     &           'field exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:2152/R:split/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [splitting] is wrong.'
         ErrCha = ''
         ErrID = 'L:2164/R:split/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine supmir(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [super mirror] section of input files                     *
*       modified by K.Niita on 2004/12/01                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /smireg/ nsreg, isgrc, isgrt, ksmir, ismir

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 3, 4, 5, 6, 7, 3 * 0 /

      character dkam*6

      dimension s_mm(kvlmax), s_r0(kvlmax), s_qc(kvlmax)
      dimension s_am(kvlmax), s_wm(kvlmax)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 .and. nsreg .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'r-in' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 4

            else if( chlw(ic:ic+4) .eq. 'r-out' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 5

            else if( chlw(ic:ic+1) .eq. 'mm' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 2

            else if( chlw(ic:ic+1) .eq. 'r0' ) then

               imsq( mrsq + 1 ) = 4
               ic = ic + 2

            else if( chlw(ic:ic+2) .eq. 'qc' ) then

               imsq( mrsq + 1 ) = 5
               ic = ic + 2

            else if( chlw(ic:ic+2) .eq. 'am' ) then

               imsq( mrsq + 1 ) = 6
               ic = ic + 2

            else if( chlw(ic:ic+2) .eq. 'wm' ) then

               imsq( mrsq + 1 ) = 7
               ic = ic + 2

            else

               goto 997

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 7 ) goto 997

                  irin = 0
                  irot = 0
                  immm = 0
                  imr0 = 0
                  imqc = 0
                  imam = 0
                  imwm = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) irin = irin + 1
                  if( imsq(k) .eq. 2 ) irot = irot + 1
                  if( imsq(k) .eq. 3 ) immm = immm + 1
                  if( imsq(k) .eq. 4 ) imr0 = imr0 + 1
                  if( imsq(k) .eq. 5 ) imqc = imqc + 1
                  if( imsq(k) .eq. 6 ) imam = imam + 1
                  if( imsq(k) .eq. 7 ) imwm = imwm + 1
               end do

                  if( irin .ne. 1 .or. irot .ne. 1 .or.
     &                immm .ne. 1 .or. imr0 .ne. 1 .or.
     &                imqc .ne. 1 .or. imam .ne. 1 .or.
     &                imwm .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  irin = 1
                  irot = 1
                  immm = 1
                  imr0 = 1
                  imqc = 1
                  imam = 1
                  imwm = 1

                  nrsq = 7

            end if

         end if

            ismir = 3

         do k = 1, nrsq

            if( ismir .eq. 3 .and. imsq(k) .eq. 1 ) ismir = 0
            if( ismir .eq. 3 .and. imsq(k) .eq. 2 ) ismir = 1

         end do

*-----------------------------------------------------------------------
*        read super mirror informations
*-----------------------------------------------------------------------

               nsreg = nsreg + 1

               if( nsreg .gt. kvlmax ) goto 998

            if( nsreg .eq. 1 ) then

                  isgrc = 1
                  iaddress_region(:) = 1
                  call moddas_reallocate_int(
     &                    2, 1, MAX_NUM_ISGRC, iaddress_region
     &                    , idas_isgrc)
                  idsm  = isgrc
                  jdsm  = -1

            end if

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 1 .or. imsq(k) .eq. 2 ) then

                     ndsm = idsm + jdsm + 3

                  call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ndsm,ierr
     &                          ,MAX_NUM_ISGRC,idas_isgrc)

                  if( ierr .ne. 0 ) goto 999

                     jdsm = jdsm + 1
                     idas_isgrc(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                     idas_isgrc(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn

            else if( imsq(k) .eq. 3 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               s_mm(nsreg) = cvvv

            else if( imsq(k) .eq. 4 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               s_r0(nsreg) = cvvv

            else if( imsq(k) .eq. 5 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               s_qc(nsreg) = cvvv

            else if( imsq(k) .eq. 6 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               s_am(nsreg) = cvvv

            else if( imsq(k) .eq. 7 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               s_wm(nsreg) = cvvv

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

               if( jdsm > MAX_NUM_ISGRC ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.supmir@read02.f ?dimension over idas_isgrc?'
     &                    //' jdsm > MAX_NUM_ISGRC'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_ISGRC@moddas.f=',MAX_NUM_ISGRC,')'
                  ErrID = 'L:2495/R:supmir/F:read02.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 2, 1, jdsm+1, iaddress_region, idas_isgrc)

               ksmir = 0
               call moddas_allocate_dbl(nsreg*5+1, das_ksmir)
               idsm  = ksmir
               jdsm  = 0


            do i = 1, nsreg

               jdsm = jdsm + 1
               das_ksmir(idsm+jdsm) = s_mm(i)
               jdsm = jdsm + 1
               das_ksmir(idsm+jdsm) = s_r0(i)
               jdsm = jdsm + 1
               das_ksmir(idsm+jdsm) = s_qc(i)
               jdsm = jdsm + 1
               das_ksmir(idsm+jdsm) = s_am(i)
               jdsm = jdsm + 1
               das_ksmir(idsm+jdsm) = s_wm(i)

            end do

            if( ismir .ne. 0 ) nsreg = - nsreg

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:2536/R:supmir/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[super mirror] is wrong.'
         ErrCha = ''
         ErrID = 'L:2549/R:supmir/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of region for super mirror '//
     &           'field exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:2563/R:supmir/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [super mirror] is wrong.'
         ErrCha = ''
         ErrID = 'L:2575/R:supmir/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine wwbias(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [WW Bias] section of input files                          *
*       modified by K.Niita on 2016/12/04                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_mesh
      use moddas_region

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'
      include 'param01.inc'   ! kitamura22/03/31

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /wwbiasn/ eenwb(6,100), iwbdp, mnwbp(6,0:20), kfwbp(6),
     &                 inwbc(6), ienwb(6), inwbt(6), iswbp, maxwb

*-----------------------------------------------------------------------

      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wbxyz/ iwbxty(6), iwbxnm(6), iwbxrg(6),
     &               rwbxmi(6), rwbxma(6), rwbxdl(6),
     &               iwbyty(6), iwbynm(6), iwbyrg(6),
     &               rwbymi(6), rwbyma(6), rwbydl(6),
     &               iwbzty(6), iwbznm(6), iwbzrg(6),
     &               rwbzmi(6), rwbzma(6), rwbzdl(6)
      common /wbtrs/ iwbtr(6,4), rwbtr(6,13)

*-----------------------------------------------------------------------

      logical dnen2

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(100)
      dimension iwwm(100)

      data imsq / 1, 99 * 0 /

      character dkam*6

      dimension iptyp(20)
      dimension jstyp(6), jnkf0(6)
      dimension vtrs(13)

      common /subtra/ isubt, ipsub(mxpart)   ! kitamura22/03/31

      dimension     idas(1)
      equivalence ( das, idas )

      integer :: kfwbp_tmp(2),inwbc_tmp(2)

      integer, allocatable :: idas_inwbc_tmp(:),idas_temporary(:)
      double precision, allocatable :: das_kfwbp_tmp(:)

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0
            inpat = 0
            ieng  = 0
            itim  = 0
            ietb  = 0
            ndefl = 0
            iregd = 0
            nnww0 = 0

            igkst  = 0
            idtt   = 0
            ktrs   = 0
            imesh  = 0

            iwbdp = iwbdp + 1

            if( iwbdp .gt. 6 ) goto 991

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        another definition line
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. iregd .gt. 0 .and.
     &        ( chlw(i1:i1+2) .eq. 'reg' .or.
     &          chlw(i1:i1+2) .eq. 'non' .or.
     &          chlw(i1:i1+2) .eq. 'xyz' .or.
     &          chlw(i1:i1+2) .eq. 'tet' .or.
     &          chlw(i1:i1+2) .eq. 'wwb' ) ) then

               if( ieng .eq. 0 ) goto 999

               if( iregd .eq. 1 ) nnww0 = nnwwp
               if( iregd .gt. 2 .and. nnwwp .ne. nnww0 ) goto 983

               ndefl = 1
               goto 110

            end if

*-----------------------------------------------------------------------
*        definition of particle, energy, time bins
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               ic  = i1
               icl = i1

               chlc = chlw
               call chcomp(chlc,icl,i3,i5)

*-----------------------------------------------------------------------
*           mesh
*-----------------------------------------------------------------------

            if( chlc(ic:ic+3) .eq. 'mesh' ) then

               ic = inumc(chlw,icl+3,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 994

                  if( chlw(ic:ic+2) .eq. 'reg' ) then

                     imesh = 1
                     if( iwmsh .ne. 0 .and. imesh .ne. iwmsh ) goto 982

                     goto 140

                  else if( chlw(ic:ic+2) .eq. 'xyz' ) then

                     imesh = 3
                     if( iwmsh .ne. 0 .and. imesh .ne. iwmsh ) goto 982

                  else if( chlw(ic:ic+2) .eq. 'tet' ) then

                     imesh = 4
                     if( iwmsh .ne. 0 .and. imesh .ne. iwmsh ) goto 982

                  else

                     goto 981

                  end if

                  if( imesh .eq. 3 ) then

                     call txymesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                            jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                            ixtp,inx,xmin,xmax,xdel,istxg,
     &                            iytp,iny,ymin,ymax,ydel,istyg,
     &                            iztp,inz,zmin,zmax,zdel,istzg)

                        if( ierr .ne. 0 ) return
                        if( jpn  .eq. 3 ) goto 981

                        if( ixtp .lt. 0 ) goto 981
                        if( iytp .lt. 0 ) goto 981
                        if( iztp .lt. 0 ) goto 981

                        if( inx .ge. 10000 ) goto 977
                        if( iny .ge. 10000 ) goto 977
                        if( inz .ge. 10000 ) goto 977

           call avoidint(ixtp,inx,xmin,xmax,xdel) ! T.Sato 2024/03/23
                        iwbxty(iwbdp) = ixtp
                        iwbxnm(iwbdp) = inx
                        rwbxdl(iwbdp) = xdel
                        rwbxmi(iwbdp) = xmin
                        rwbxma(iwbdp) = xmax

                        call moddas_reallocate_dbl(
     &                          6, iwbdp, inx+1, iwbxrg, das_iwbxrg)
                        das_iwbxrg(iwbxrg(iwbdp):iwbxrg(iwbdp)+inx)
     &                           = gmsh(istxg:istxg+inx)

           call avoidint(iytp,iny,ymin,ymax,ydel) ! T.Sato 2024/03/23
                        iwbyty(iwbdp) = iytp
                        iwbynm(iwbdp) = iny
                        rwbydl(iwbdp) = ydel
                        rwbymi(iwbdp) = ymin
                        rwbyma(iwbdp) = ymax

                        call moddas_reallocate_dbl(
     &                          6, iwbdp, iny+1, iwbyrg, das_iwbyrg)
                        das_iwbyrg(iwbyrg(iwbdp):iwbyrg(iwbdp)+iny)
     &                           = gmsh(istyg:istyg+iny)

           call avoidint(iztp,inz,zmin,zmax,zdel) ! T.Sato 2024/03/23
                        iwbzty(iwbdp) = iztp
                        iwbznm(iwbdp) = inz
                        rwbzdl(iwbdp) = zdel
                        rwbzmi(iwbdp) = zmin
                        rwbzma(iwbdp) = zmax

                        call moddas_reallocate_dbl(
     &                          6, iwbdp, inz+1, iwbzrg, das_iwbzrg)
                        das_iwbzrg(iwbzrg(iwbdp):iwbzrg(iwbdp)+inz)
     &                           = gmsh(istzg:istzg+inz)

                   if( iwbdp .gt. 1 ) then

                     do kk = 1, iwbdp - 1

                        if( iwbxty(kk) .ne. ixtp ) goto 980
                        if( iwbxnm(kk) .ne. inx  ) goto 980
                        if( rwbxdl(kk) .ne. xdel ) goto 980
                        if( rwbxmi(kk) .ne. xmin ) goto 980
                        if( rwbxma(kk) .ne. xmax ) goto 980

                        if( iwbyty(kk) .ne. iytp ) goto 980
                        if( iwbynm(kk) .ne. iny  ) goto 980
                        if( rwbydl(kk) .ne. ydel ) goto 980
                        if( rwbymi(kk) .ne. ymin ) goto 980
                        if( rwbyma(kk) .ne. ymax ) goto 980

                        if( iwbzty(kk) .ne. iztp ) goto 980
                        if( iwbznm(kk) .ne. inz  ) goto 980
                        if( rwbzdl(kk) .ne. zdel ) goto 980
                        if( rwbzmi(kk) .ne. zmin ) goto 980
                        if( rwbzma(kk) .ne. zmax ) goto 980

                     end do

                   end if

                   goto 150

                  else if( imesh .eq. 4)then

                   call moddas_allocate_int(
     &                  MAX_NUM_INWBC, idas_temporary)

                   call ttetmesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                  mtrn,ntrn,
     &                  MAX_NUM_INWBC,idas_temporary)

                   ireg=idas_temporary(1)

                   if( mtrn > MAX_NUM_INWBC ) then
                    write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.wwbias@read02.f'
     &                    //' ?dimension over idas_temporary?'
     &                    //' mtrn > MAX_NUM_INWBC'
     &                 ,' (mtrn=',mtrn,')'
     &                 ,' (MAX_NUM_INWBC@moddas.f=',MAX_NUM_INWBC,')'
                    ErrID = 'L:2879/R:wwbias/F:read02.f'
                    call ErrWrite(ErrID,ErrCha)
                   endif

                   if( ierr .ne. 0 ) return
                   if( jpn  .eq. 3 ) goto 975

                   if( iwbdp .gt. 1 ) then

                    do kk = 1, iwbdp - 1
                     if( idas_inwbc(inwbc(kk)+1) .ne. ireg ) goto 974
                    enddo

                   endif

                   goto 150

                  end if

            end if

*-----------------------------------------------------------------------
*           transform
*-----------------------------------------------------------------------

            if( chlc(ic:ic+3) .eq. 'trcl' .or.
     &          chlc(ic:ic+4) .eq. '*trcl') then

              if( chlc(ic:ic+3) .eq. 'trcl') then
                 ic = ic + 4
              else if( chlc(ic:ic+4) .eq. '*trcl') then
                 ktrs = 1
                 ic = ic + 5
              end if

               ic = inumc(chlw,ic,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               call ttrans(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     ic,ktrs,igkst,idtt,vtrs)

               if( ierr .ne. 0 ) goto 999

                     iwbtr(iwbdp,1) = igkst
                     iwbtr(iwbdp,2) = ktrs
                     iwbtr(iwbdp,3) = idtt
                     iwbtr(iwbdp,4) = 0

               if( igkst .gt. 1 ) then

                  do k = 1, 13

                     rwbtr(iwbdp,k) = vtrs(k)

                  end do

               end if

               goto 150

            end if

*-----------------------------------------------------------------------
*           energy bin
*-----------------------------------------------------------------------

            if( chlc(ic:ic+2) .eq. 'eng' .or.
     &          chlc(ic:ic+2) .eq. 'tim' ) then

               if( ietb .ne. 0 ) goto 986

               if(chlc(icl:icl+2) .eq. 'eng' ) ietb =  1
               if(chlc(icl:icl+2) .eq. 'tim' ) ietb = -1

               ic = inumc(chlw,icl+3,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 989

                  ieng = nint( cvvv )

               if( ieng .lt. 1 .or. ieng .gt. 99 ) goto 988

*-----------------------------------------------------------------------

  151                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( ierr .ne. 0 ) goto 987
                        if( jpn  .eq. 3 ) goto 987
                        if( iskip .ne. 0 ) goto 151

                  ic = i1

               do k = 1, ieng

                  if( ic .gt. i3 ) then

  152                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) goto 987
                     if( jpn  .eq. 3 ) goto 987
                     if( iskip .ne. 0 ) goto 152

                     ic = i1

                  end if

                     ic = jnumc(chlw,ic,i3)

                     call snum(chlw,ic,i3,ic2,cvvv,ierr)

                     if( ierr .ne. 0 ) goto 987

                     eenwb(iwbdp,k) = cvvv

                     ic = ic2

               end do

                     goto 140

            end if

*-----------------------------------------------------------------------
*           particle name
*-----------------------------------------------------------------------

            if( chlc(icl:icl+3) .eq. 'part' ) then

               ic = inumc(chlw,icl+4,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = i3

  400          continue

                  ic = jnumc(chlw,ic,i3)

                  if( ic .gt. i3 ) goto 140

*-----------------------------------------------------------------------

                  call rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

                     if( ierr .eq. 994 ) goto 993
                     if( ierr .eq. 998 ) goto 993
                     if( isubt .eq. 1 )  goto 993   ! kitamura22/03/31

                  if( istyp .gt. 0 ) then

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = istyp

                  else if( istyp .lt. 0 ) then

                     do i = 1, -istyp

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = jstyp(i)

                     end do

                  end if

                  goto 400

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

  110    continue

         if( nrsq .eq. 0 .or. ndefl .eq. 1 ) then

                  iregd = iregd + 1
                  nnwwp = 0
                  nrsq  = 0
                  ndefl = 0

*-----------------------------------------------------------------------

               if( iregd .eq. 1 .and. ieng .eq. 0 ) then

                  ieng = 1

                  ienwb(iwbdp) = 0

                  if( imesh .eq. 1 ) then

                     call moddas_reallocate_dbl(
     &                       6, iwbdp, kvlmax, kfwbp, das_kfwbp)

                     do j = 1, kvlmax

                        das_kfwbp(kfwbp(iwbdp)+j-1) = 0.0d0

                     end do

                  else if( imesh .eq. 3 ) then

                     ixyz = inx * iny * inz

                     call moddas_reallocate_dbl(
     &                       6, iwbdp, ixyz, kfwbp, das_kfwbp)

                     do j = 1, ixyz

                           das_kfwbp(kfwbp(iwbdp)+j-1) = 0.0d0

                     end do

                  else if( imesh .eq. 4 ) then

                   if( ndefl .eq. 0 )then

                     call moddas_allocate_int(
     &                   kvlmax, idas_inwbc_tmp)

                     call moddas_allocate_dbl(
     &                       kvlmax, das_kfwbp_tmp)

                     inwbc_tmp(1) = kvlmax
                     inwbc_tmp(2) = kvlmax
                     kfwbp_tmp(1) = kvlmax
                     kfwbp_tmp(2) = kvlmax

                     do j=1,kvlmax

                           idas_inwbc_tmp(j) = 0
                           das_kfwbp_tmp(j) = 0.0d0

                     enddo

                   endif

                 end if

                     call moddas_reallocate_dbl(
     &                       6, iwbdp, kvlmax, kfwbp, das_kfwbp)
!--->

                  do j = 1, kvlmax

                        das_kfwbp(kfwbp(iwbdp)+j-1) = 0.0d0

                  end do

*-----------------------------------------------------------------------

               else if( iregd .eq. 1 .and. ieng .ne. 0 ) then

                     ienwb(iwbdp) = ieng * ietb

                  if( imesh .eq. 1 ) then

                     call moddas_reallocate_dbl(
     &                       6, iwbdp, kvlmax*ieng, kfwbp, das_kfwbp)

                     do k = 1, ieng
                     do j = 1, kvlmax

                        das_kfwbp(kfwbp(iwbdp)+(k-1)*kvlmax+j-1) = 0.0d0

                     end do
                     end do

                  else if( imesh .eq. 3 ) then

                     ixyz = inx * iny * inz

                     call moddas_reallocate_dbl(
     &                       6, iwbdp, ixyz*ieng, kfwbp, das_kfwbp)

                     do k = 1, ieng
                     do j = 1, ixyz

                        das_kfwbp(kfwbp(iwbdp)+(k-1)*ixyz+j-1) = 0.0d0

                     end do
                     end do

                  else if( imesh .eq. 4 ) then

                   if(ndefl.eq.0)then

                     call moddas_allocate_int(
     &                       kvlmax*ieng, idas_inwbc_tmp)

                     call moddas_allocate_dbl(
     &                       kvlmax*ieng, das_kfwbp_tmp)

                     inwbc_tmp(1) = kvlmax*ieng
                     inwbc_tmp(2) = kvlmax*ieng
                     kfwbp_tmp(1) = kvlmax*ieng
                     kfwbp_tmp(2) = kvlmax*ieng

                     do j = 1, kvlmax*ieng

                        idas_inwbc_tmp(j) = 0
                        das_kfwbp_tmp(j) = 0.0d0

                     end do

                   endif

                  end if

               end if

*-----------------------------------------------------------------------

               if( iregd .eq. 1 .and. imesh .eq. 1 ) then

                  call moddas_reallocate_int(
     &                    6, iwbdp, MAX_NUM_INWWC, inwbc, idas_inwbc)

                  idsm = inwbc(iwbdp)
                  jdsm = 0

               else if( iregd .eq. 1 .and. imesh .eq. 3 ) then

                  idsm = inwbc(iwbdp)
                  jdsm = inx * iny * inz * 3

                  call moddas_reallocate_int(
     &                    6, iwbdp, jdsm+1, inwbc, idas_inwbc)

                  do k = 1, jdsm

                     idas_inwbc(idsm+k) = 0

                  end do

               else if( iregd .eq. 1 .and. imesh .eq. 4 ) then

                  idsm = inwbc(iwbdp)

               end if

*-----------------------------------------------------------------------

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'wwb' ) then

               imsq( mrsq + 1 ) = 2

               ic = ic + 3
               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               ic = ic2

               if( ierr .ne. 0 ) goto 985

               iwwm( mrsq + 1 ) = nint( cvvv )

               if( iwwm(mrsq+1) .le. 0 .or.
     &             iwwm(mrsq+1) .gt. ieng ) goto 985

            else if( chlw(ic:ic+2) .eq. 'xyz' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'tet' ) then

               imsq( mrsq + 1 ) = 4
               ic = ic + 3

            else

               goto 997

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

*-----------------------------------------------------------------------

  200       continue

            if( mrsq .gt. 0 ) then

                  inreg = 0
                  inwwp = 0
                  innon = 0
                  inxyz = 0
                  intet = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) inwwp = inwwp + 1
                  if( imsq(k) .eq. 3 ) inxyz = inxyz + 1
                  if( imsq(k) .eq. 4 ) intet = intet + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( imesh .eq. 1 ) then
                     if( iregd .eq. 1 .and. inreg .ne. 1 ) goto 997
                     if( inreg .gt. 1 ) goto 997
                  else if( imesh .eq. 3 ) then
                     if( iregd .eq. 1 .and. inxyz .ne. 1 ) goto 997
                     if( inxyz .gt. 1 ) goto 997
                  else if( imesh .eq. 4 ) then
                     if( iregd .eq. 1 .and. intet .ne. 1 ) goto 997
                     if( intet .gt. 1 ) goto 997
                  end if

                  nrsq = mrsq

                  goto 140

            else

                  nrsq = ieng + 1

                  if( imesh .eq. 1 ) then
                     imsq(1) = 1
                  else if( imesh .eq. 3)then
                     imsq(1) = 3
                  else if( imesh .eq. 4)then
                     imsq(1) = 4
                  end if

               do k = 2, nrsq

                  imsq(k) = 2
                  iwwm(k) = k - 1

               end do

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        summary of particle
*-----------------------------------------------------------------------

            if( inpat .gt. 0 ) then

               do k = 1, inpat

                  if( iptyp(k) .eq. 20 ) inpat = 0

               end do

            end if

            if( inpat .eq. 0 ) then

               do k = 1, 19

                  iptyp(k) = k

               end do

                  inpat = 19

            end if

*-----------------------------------------------------------------------
*        read weight window lower bound
*-----------------------------------------------------------------------

               nnwwp = nnwwp + 1

               if( nnwwp .gt. kvlmax ) goto 998

               if( imesh .eq. 1 .and. nnwwp .gt. kvlmax ) goto 998

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 .or.
     &              (( imsq(k) .eq. 1 .or. imsq(k) .eq. 3
     &              .or. imsq(k) . eq. 4 ) .and.
     &          iregd .gt. 1 ) ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 2 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( imesh .eq. 1 ) then

                       kdsm = kfwbp(iwbdp)
                       ldsm = iwwm(k)

                     das_kfwbp(kdsm+(ldsm-1)*kvlmax+nnwwp-1) = cvvv

                  else if( imesh .eq. 3 ) then

                        kdsm = kfwbp(iwbdp)
                        ldsm = iwwm(k)

                     das_kfwbp(kdsm+(ldsm-1)*ixyz+nnwwp-1) = cvvv

                  else if( imesh .eq. 4 ) then

                     das_kfwbp_tmp((nnwwp-1)*ieng+iwwm(k)) = cvvv

                  end if

            else if( imsq(k) .eq. 1 .and. iregd .eq. 1 ) then

                     ndsm = idsm + jdsm + 3

                  call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ndsm,ierr
     &                          ,MAX_NUM_INWBC,idas_inwbc)

                  if( ierr .ne. 0 ) goto 999

                     jdsm = jdsm + 1
                     idas_inwbc(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                     idas_inwbc(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn

            else if( imsq(k) .eq. 3 .and. iregd .eq. 1 ) then

               if( chlw(ic:ic) .eq. '(' ) then

                     ic = ic + 1
                     call snum(chlw,ic,i3,ic2,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 999

                     idas_inwbc(idsm+(nnwwp-1)*3+1-1) = nint( cvvv )

                     ic = jnumc(chlw,ic2,i3)
                     call snum(chlw,ic,i3,ic2,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 999

                     idas_inwbc(idsm+(nnwwp-1)*3+2-1) = nint( cvvv )

                     ic = jnumc(chlw,ic2,i3)
                     ib = inumc(chlw,ic,i3,')') - 1

                     call snum(chlw,ic,ib,ic2,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 999

                     idas_inwbc(idsm+(nnwwp-1)*3+3-1) = nint( cvvv )

                     ic2 = inumc(chlw,ic,i3,')') + 1

               else

                  goto 978

               end if

            else if( imsq(k) .eq. 4 .and. iregd .eq. 1 ) then

             if(nnwwp*ieng.gt.kfwbp_tmp(2))then

              call moddas_reallocate_int(
     &             2, 1, ieng*kvlmax, inwbc_tmp, idas_inwbc_tmp)

              call moddas_reallocate_dbl(
     &             2, 1, ieng*kvlmax, kfwbp_tmp, das_kfwbp_tmp)

              inwbc_tmp(1)=inwbc_tmp(2)
              kfwbp_tmp(1)=kfwbp_tmp(2)

             endif

             call snum(chlw,ic,i3,ic2,cvvv,ierr)

             if( ierr .ne. 0 ) goto 999

             idas_inwbc_tmp(nnwwp) = nint( cvvv )

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

            if( iregd .eq. 1 ) nnww0 = nnwwp
            if( iregd .gt. 1 .and. nnwwp .ne. nnww0 ) goto 983

            if( nnww0 .eq. 0 ) goto 999

               mnwbp(iwbdp, 0) = nnww0
               mnwbp(iwbdp,20) = inpat

            do j = 1, inpat

               mnwbp(iwbdp,j) = iptyp(j)

            end do

            if( iwbdp .gt. 1 .and. iwmsh .ne. imesh ) goto 976

               iwmsh = imesh

            if( imesh .eq. 1 ) then

               if( jdsm > MAX_NUM_INWBC ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.wwbias@read02.f ?dimension over idas_inwbc?'
     &                    //' jdsm > MAX_NUM_INWBC'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INWBC@moddas.f=',MAX_NUM_INWBC,')'
                  ErrID = 'L:3551/R:wwbias/F:read02.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 6, iwbdp, jdsm+1, inwbc, idas_inwbc)

*-----------------------------------------------------------------------

            else if( imesh .eq. 4 ) then

             call moddas_reallocate_int(
     &            6, iwbdp, nnww0+mtrn+1, inwbc, idas_inwbc )

             call moddas_reallocate_dbl(
     &            6, iwbdp, nnww0*ieng, kfwbp, das_kfwbp )

             idas_inwbc(idsm+1:idsm+mtrn)
     &            =idas_temporary(1:mtrn)

             idas_inwwc(idsm+mtrn+1)=nnww0

               do j = 1, nnww0

                  idas_inwbc(idsm+mtrn+j+1)
     &                = idas_inwbc_tmp(j)

                  do k = 1, ieng

                   das_kfwbp(kfwbp(iwbdp)+(k-1)*nnww0+j-1)
     &                  = das_kfwbp_tmp((j-1)*ieng+k)

                  enddo

               enddo

               call moddas_deallocate_int( idas_temporary )

               call moddas_deallocate_int( idas_inwbc_tmp )

               call moddas_deallocate_dbl( das_kfwbp_tmp )

            end if

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

         return

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:3606/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

 974  continue

         m_err = 'ireg for tet mesh should be consistent'//
     &        ' for each WW section'
         ErrCha = ''
         ErrID = 'L:3619/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return


  975 continue

         m_err = 'Description of mesh = tet is wrong'
         ErrCha = ''
         ErrID = 'L:3630/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  976 continue

         m_err = 'mesh should be the same for each WW section'
         ErrCha = ''
         ErrID = 'L:3642/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  977 continue

         m_err = 'Each mesh of (ix iy iz) should be smaller than 10000'
         ErrCha = ''
         ErrID = 'L:3652/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  978 continue

         m_err = 'Description of (ix iy iz) is wrong'
         ErrCha = ''
         ErrID = 'L:3662/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  979 continue

         m_err = 'mesh should be defined before part or eng'
         ErrCha = ''
         ErrID = 'L:3672/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  980 continue

         m_err = 'xyz mesh should be the same for each WW section'
         ErrCha = ''
         ErrID = 'L:3682/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  981 continue

         m_err = 'Description of mesh = xyz is wrong'
         ErrCha = ''
         ErrID = 'L:3692/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  982 continue

         m_err = 'mesh = reg, xyz, tet are not allowed at the same time'
         ErrCha = ''
         ErrID = 'L:3702/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  983 continue

         m_err = 'region definition is not consistent'
         ErrCha = ''
         ErrID = 'L:3714/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  985 continue

         m_err = 'number after wwb is wrong'
         ErrCha = ''
         ErrID = 'L:3726/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  986 continue

         m_err = 'eng =, or tim = is defined twice.'
         ErrCha = ''
         ErrID = 'L:3738/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  987 continue

         m_err = 'energy or time bins is wrong'
         ErrCha = ''
         ErrID = 'L:3750/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  988 continue

         m_err = 'energy or time bins should be 0 < n < 99.'
         ErrCha = ''
         ErrID = 'L:3762/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  989 continue

         m_err = 'number of eng =, tim = is wrong'
         ErrCha = ''
         ErrID = 'L:3774/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'Number of [WW Bias] section is'//
     &           ' too large (=<6)'
         ErrCha = ''
         ErrID = 'L:3787/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'Number of particles is too large (<19)'
         ErrCha = ''
         ErrID = 'L:3799/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'Name of particle is wrong'
         ErrCha = ''
         ErrID = 'L:3811/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'Description of part, eng, time = is wrong'
         ErrCha = ''
         ErrID = 'L:3823/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[WW Bias] is wrong.'
         ErrCha = ''
         ErrID = 'L:3836/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of WW Bias region '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:3850/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [WW Bias] is wrong.'
         ErrCha = ''
         ErrID = 'L:3862/R:wwbias/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine wwind(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [weightwindow] section of input files                     *
*       modified by K.Niita on 2017/12/28                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_mesh
      use moddas_region

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'
      include 'param01.inc'   ! kitamura22/03/31

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww

*-----------------------------------------------------------------------

      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)
      common /wwtrs/ iwwtr(6,4), rwwtr(6,13)

*-----------------------------------------------------------------------

      logical dnen2

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(100)
      dimension iwwm(100)

      data imsq / 1, 99 * 0 /

      character dkam*6

      dimension iptyp(20)
      dimension jstyp(6), jnkf0(6)
      dimension vtrs(13)

      common /subtra/ isubt, ipsub(mxpart)   ! kitamura22/03/31

      dimension     idas(1)
      equivalence ( das, idas )

      integer :: kfwwp_tmp(2),inwwc_tmp(2)

      integer, allocatable :: idas_inwwc_tmp(:),idas_temporary(:)
      double precision, allocatable :: das_kfwwp_tmp(:)

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0
            inpat = 0
            ieng  = 0
            itim  = 0
            ietb  = 0
            ndefl = 0
            iregd = 0
            nnww0 = 0
            kecho = 0
            idval = 0
            dvals = 1.0d0

            igkst  = 0
            idtt   = 0
            ktrs   = 0

            imesh = 1

            iwwdp = iwwdp + 1

            if( iwwdp .gt. 6 ) goto 991

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        another definition line
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. iregd .gt. 0 .and.
     &        ( chlw(i1:i1+2) .eq. 'reg' .or.
     &          chlw(i1:i1+2) .eq. 'non' .or.
     &          chlw(i1:i1+2) .eq. 'xyz' .or.
     &          chlw(i1:i1+2) .eq. 'tet' .or.
     &          chlw(i1:i1+1) .eq. 'ww' ) ) then

               if( ieng .eq. 0 ) goto 999

               if( iregd .eq. 1 ) nnww0 = nnwwp
               if( iregd .gt. 1 .and. nnwwp .ne. nnww0 ) goto 983

               ndefl = 1
               goto 110

            end if

*-----------------------------------------------------------------------
*        definition of particle, energy, time bins
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               ic  = i1
               icl = i1

               chlc = chlw
               call chcomp(chlc,icl,i3,i5)

*-----------------------------------------------------------------------
*           mesh
*-----------------------------------------------------------------------

            if( chlc(ic:ic+3) .eq. 'mesh' ) then

               ic = inumc(chlw,icl+3,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 994

                  if( chlw(ic:ic+2) .eq. 'reg' ) then

                     imesh = 1
                     if( iwmsh .ne. 0 .and. imesh .ne. iwmsh ) goto 982

                     goto 140

                  else if( chlw(ic:ic+2) .eq. 'xyz' ) then

                     imesh = 3
                     if( iwmsh .ne. 0 .and. imesh .ne. iwmsh ) goto 982

                  else if( chlw(ic:ic+2) .eq. 'tet' ) then

                     imesh = 4
                     if( iwmsh .ne. 0 .and. imesh .ne. iwmsh ) goto 982

                  else

                     goto 981

                  end if

                  if( imesh .eq. 3 ) then

                     call txymesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                            jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                            ixtp,inx,xmin,xmax,xdel,istxg,
     &                            iytp,iny,ymin,ymax,ydel,istyg,
     &                            iztp,inz,zmin,zmax,zdel,istzg)

                        if( ierr .ne. 0 ) return
                        if( jpn  .eq. 3 ) goto 981

                        if( ixtp .lt. 0 ) goto 981
                        if( iytp .lt. 0 ) goto 981
                        if( iztp .lt. 0 ) goto 981

                        if( inx .ge. 10000 ) goto 977
                        if( iny .ge. 10000 ) goto 977
                        if( inz .ge. 10000 ) goto 977

           call avoidint(ixtp,inx,xmin,xmax,xdel) ! T.Sato 2024/03/23
                        iwxty(iwwdp) = ixtp
                        iwxnm(iwwdp) = inx
                        rwxdl(iwwdp) = xdel
                        rwxmi(iwwdp) = xmin
                        rwxma(iwwdp) = xmax
                        call moddas_reallocate_dbl(
     &                          6, iwwdp, inx+1, iwxrg, das_iwxrg)
                        das_iwxrg(iwxrg(iwwdp):iwxrg(iwwdp)+inx)
     &                           = gmsh(istxg:istxg+inx)

           call avoidint(iytp,iny,ymin,ymax,ydel) ! T.Sato 2024/03/23
                        iwyty(iwwdp) = iytp
                        iwynm(iwwdp) = iny
                        rwydl(iwwdp) = ydel
                        rwymi(iwwdp) = ymin
                        rwyma(iwwdp) = ymax
                        call moddas_reallocate_dbl(
     &                          6, iwwdp, iny+1, iwyrg, das_iwyrg)
                        das_iwyrg(iwyrg(iwwdp):iwyrg(iwwdp)+iny)
     &                           = gmsh(istyg:istyg+iny)

           call avoidint(iztp,inz,zmin,zmax,zdel) ! T.Sato 2024/03/23
                        iwzty(iwwdp) = iztp
                        iwznm(iwwdp) = inz
                        rwzdl(iwwdp) = zdel
                        rwzmi(iwwdp) = zmin
                        rwzma(iwwdp) = zmax
                        call moddas_reallocate_dbl(
     &                          6, iwwdp, inz+1, iwzrg, das_iwzrg)
                        das_iwzrg(iwzrg(iwwdp):iwzrg(iwwdp)+inz)
     &                           = gmsh(istzg:istzg+inz)

                    if( iwwdp .gt. 1 ) then

                     do kk = 1, iwwdp - 1

                        if( iwxty(kk) .ne. ixtp ) goto 980
                        if( iwxnm(kk) .ne. inx  ) goto 980
                        if( rwxdl(kk) .ne. xdel ) goto 980
                        if( rwxmi(kk) .ne. xmin ) goto 980
                        if( rwxma(kk) .ne. xmax ) goto 980

                        if( iwyty(kk) .ne. iytp ) goto 980
                        if( iwynm(kk) .ne. iny  ) goto 980
                        if( rwydl(kk) .ne. ydel ) goto 980
                        if( rwymi(kk) .ne. ymin ) goto 980
                        if( rwyma(kk) .ne. ymax ) goto 980

                        if( iwzty(kk) .ne. iztp ) goto 980
                        if( iwznm(kk) .ne. inz  ) goto 980
                        if( rwzdl(kk) .ne. zdel ) goto 980
                        if( rwzmi(kk) .ne. zmin ) goto 980
                        if( rwzma(kk) .ne. zmax ) goto 980

                     end do

                    end if

                     goto 150

                  else if( imesh .eq. 4)then


                   call moddas_allocate_int(
     &                  MAX_NUM_INWWC, idas_temporary)

                   call ttetmesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                  mtrn,ntrn,
     &                  MAX_NUM_INWWC,idas_temporary)

                   ireg=idas_temporary(1)

                   if( mtrn > MAX_NUM_INWWC ) then
                    write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.wwind@read02.f'
     &                    //' ?dimension over idas_temporary?'
     &                    //' mtrn > MAX_NUM_INWWC'
     &                 ,' (mtrn=',mtrn,')'
     &                 ,' (MAX_NUM_iNWWC@moddas.f=',MAX_NUM_INWWC,')'
                    ErrID = 'L:4168/R:wwind/F:read02.f'
                    call ErrWrite(ErrID,ErrCha)
                   endif

                   if( ierr .ne. 0 ) return
                   if( jpn  .eq. 3 ) goto 975

                   if( iwwdp .gt. 1 ) then

                    do kk = 1, iwwdp - 1
                     if( idas_inwwc(inwwc(kk)+1) .ne. ireg ) goto 974
                    enddo

                   endif

                   goto 150

                  end if

            end if

*-----------------------------------------------------------------------
*           transform
*-----------------------------------------------------------------------

            if( chlc(ic:ic+3) .eq. 'trcl' .or.
     &          chlc(ic:ic+4) .eq. '*trcl') then

              if( chlc(ic:ic+3) .eq. 'trcl') then
                 ic = ic + 4
              else if( chlc(ic:ic+4) .eq. '*trcl') then
                 ktrs = 1
                 ic = ic + 5
              end if

               ic = inumc(chlw,ic,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               call ttrans(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     ic,ktrs,igkst,idtt,vtrs)

               if( ierr .ne. 0 ) goto 999

                     iwwtr(iwwdp,1) = igkst
                     iwwtr(iwwdp,2) = ktrs
                     iwwtr(iwwdp,3) = idtt
                     iwwtr(iwwdp,4) = 0

               if( igkst .gt. 1 ) then

                  do k = 1, 13

                     rwwtr(iwwdp,k) = vtrs(k)

                  end do

               end if

               goto 150

            end if

*-----------------------------------------------------------------------
*           default value
*-----------------------------------------------------------------------

            if( chlc(ic:ic+3) .eq. 'dval' ) then

               ic = inumc(chlw,icl+3,i3,'=') + 1

               if( ic .gt. i3 ) goto 999

               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 989

                  dvinp = cvvv
                  idval = 1

               goto 140

            end if

*-----------------------------------------------------------------------
*           echo choice
*-----------------------------------------------------------------------

            if( chlc(ic:ic+3) .eq. 'echo' ) then

               ic = inumc(chlw,icl+3,i3,'=') + 1

               if( ic .gt. i3 ) goto 999

               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 989

                  kecho = nint( cvvv )

               goto 140

            end if

*-----------------------------------------------------------------------
*           energy bin
*-----------------------------------------------------------------------

            if( chlc(ic:ic+2) .eq. 'eng' .or.
     &          chlc(ic:ic+2) .eq. 'tim' ) then

               if( ietb .ne. 0 ) goto 986

               if(chlc(icl:icl+2) .eq. 'eng' ) ietb =  1
               if(chlc(icl:icl+2) .eq. 'tim' ) ietb = -1

               ic = inumc(chlw,icl+3,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 989

                  ieng = nint( cvvv )

               if( ieng .lt. 1 .or. ieng .gt. 99 ) goto 988

*-----------------------------------------------------------------------

  151                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( ierr .ne. 0 ) goto 987
                        if( jpn  .eq. 3 ) goto 987
                        if( iskip .ne. 0 ) goto 151

                  ic = i1

               do k = 1, ieng

                  if( ic .gt. i3 ) then

  152                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) goto 987
                     if( jpn  .eq. 3 ) goto 987
                     if( iskip .ne. 0 ) goto 152

                     ic = i1

                  end if

                     ic = jnumc(chlw,ic,i3)

                     call snum(chlw,ic,i3,ic2,cvvv,ierr)

                     if( ierr .ne. 0 ) goto 987

                     eenww(iwwdp,k) = cvvv

                     ic = ic2

               end do

                     goto 140

            end if

*-----------------------------------------------------------------------
*           particle name
*-----------------------------------------------------------------------

            if( chlc(icl:icl+3) .eq. 'part' ) then

               ic = inumc(chlw,icl+4,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = i3

  400          continue

                  ic = jnumc(chlw,ic,i3)

                  if( ic .gt. i3 ) goto 140

*-----------------------------------------------------------------------

                  call rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

                     if( ierr .eq. 994 ) goto 993
                     if( ierr .eq. 998 ) goto 993
                     if( isubt .eq. 1 )  goto 993   ! kitamura22/03/31

                  if( istyp .gt. 0 ) then

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = istyp

                  else if( istyp .lt. 0 ) then

                     do i = 1, -istyp

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = jstyp(i)

                     end do

                  end if

                  goto 400

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

  110    continue

         if( nrsq .eq. 0 .or. ndefl .eq. 1 ) then

                  iregd = iregd + 1
                  nnwwp = 0
                  nrsq  = 0
                  ndefl = 0

*-----------------------------------------------------------------------

               if( iregd .eq. 1 .and. ieng .eq. 0 ) then

                     ieng = 1

                     ienww(iwwdp) = 0

                  if( imesh .eq. 1 ) then

                     call moddas_reallocate_dbl(
     &                       6, iwwdp, kvlmax, kfwwp, das_kfwwp)

                     do j = 1, kvlmax

                           das_kfwwp(kfwwp(iwwdp)+j-1) = 0.0d0

                     end do

                  else if( imesh .eq. 3 ) then

                     ixyz = inx * iny * inz

                     call moddas_reallocate_dbl(
     &                       6, iwwdp, ixyz, kfwwp, das_kfwwp)

                     do j = 1, ixyz

                           das_kfwwp(kfwwp(iwwdp)+j-1) = 0.0d0

                     end do

                     call moddas_reallocate_dbl(
     &                       6, iwwdp, ixyz, kgwwp, das_kgwwp)

                     do j = 1, ixyz

                           das_kgwwp(kgwwp(iwwdp)+j-1) = -1.0d0

                     end do

                  else if( imesh .eq. 4 ) then

                   if( ndefl .eq. 0 )then

                     call moddas_allocate_int(
     &                   kvlmax, idas_inwwc_tmp)

                     call moddas_allocate_dbl(
     &                       kvlmax, das_kfwwp_tmp)

                     inwwc_tmp(1) = kvlmax
                     inwwc_tmp(2) = kvlmax
                     kfwwp_tmp(1) = kvlmax
                     kfwwp_tmp(2) = kvlmax

                     do j=1,kvlmax

                           idas_inwwc_tmp(j) = 0
                           das_kfwwp_tmp(j) = 0.0d0

                     enddo

                   endif

                  end if

*-----------------------------------------------------------------------

               else if( iregd .eq. 1 .and. ieng .ne. 0 ) then

                     ienww(iwwdp) = ieng * ietb

                  if( imesh .eq. 1 ) then

                     call moddas_reallocate_dbl(
     &                       6, iwwdp, kvlmax*ieng, kfwwp, das_kfwwp)

                     do k = 1, ieng
                     do j = 1, kvlmax

                        das_kfwwp(kfwwp(iwwdp)+(k-1)*kvlmax+j-1) = 0.0d0

                     end do
                     end do

                  else if( imesh .eq. 3 ) then

                     ixyz = inx * iny * inz

                     call moddas_reallocate_dbl(
     &                       6, iwwdp, ixyz*ieng, kfwwp, das_kfwwp)

                     do k = 1, ieng
                     do j = 1, ixyz

                        das_kfwwp(kfwwp(iwwdp)+(k-1)*ixyz+j-1) = 0.0d0

                     end do
                     end do

                     call moddas_reallocate_dbl(
     &                       6, iwwdp, ixyz*ieng, kgwwp, das_kgwwp)

                     do k = 1, ieng
                     do j = 1, ixyz

                        das_kgwwp(kgwwp(iwwdp)+(k-1)*ixyz+j-1) = -1.0d0

                     end do
                     end do

                  else if( imesh .eq. 4 ) then

                   if(ndefl.eq.0)then

                     call moddas_allocate_int(
     &                       kvlmax*ieng, idas_inwwc_tmp)

                     call moddas_allocate_dbl(
     &                       kvlmax*ieng, das_kfwwp_tmp)

                     inwwc_tmp(1) = kvlmax*ieng
                     inwwc_tmp(2) = kvlmax*ieng
                     kfwwp_tmp(1) = kvlmax*ieng
                     kfwwp_tmp(2) = kvlmax*ieng

                     do j = 1, kvlmax*ieng

                        idas_inwwc_tmp(j) = 0
                        das_kfwwp_tmp(j) = 0.0d0

                     end do

                   endif

                  end if

               end if

*-----------------------------------------------------------------------

               if( iregd .eq. 1 .and. imesh .eq. 1 ) then

                  call moddas_reallocate_int(
     &                    6, iwwdp, MAX_NUM_INWWC, inwwc, idas_inwwc)

                  idsm = inwwc(iwwdp)
                  jdsm = 0

               else if( iregd .eq. 1 .and. imesh .eq. 3 ) then

                  idsm = inwwc(iwwdp)
                  jdsm = inx * iny * inz * 3
                  call moddas_reallocate_int(
     &                    6, iwwdp, jdsm+1, inwwc, idas_inwwc)

                  do k = 1, jdsm

                     idas_inwwc(idsm+k) = 0

                  end do

               else if( iregd .eq. 1 .and. imesh .eq. 4 ) then

                  idsm = inwwc(iwwdp)

               end if

*-----------------------------------------------------------------------

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+1) .eq. 'ww' ) then

               imsq( mrsq + 1 ) = 2

               ic = ic + 2
               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               ic = ic2

               if( ierr .ne. 0 ) goto 985

               iwwm( mrsq + 1 ) = nint( cvvv )

               if( iwwm(mrsq+1) .le. 0 .or.
     &             iwwm(mrsq+1) .gt. ieng ) goto 985

            else if( chlw(ic:ic+2) .eq. 'xyz' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'tet' ) then

               imsq( mrsq + 1 ) = 4
               ic = ic + 3

            else

               goto 997

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)

               goto 100

*-----------------------------------------------------------------------

  200       continue

            if( mrsq .gt. 0 ) then

                  inreg = 0
                  inwwp = 0
                  innon = 0
                  inxyz = 0
                  intet = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) inwwp = inwwp + 1
                  if( imsq(k) .eq. 3 ) inxyz = inxyz + 1
                  if( imsq(k) .eq. 4 ) intet = intet + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( imesh .eq. 1 ) then
                     if( iregd .eq. 1 .and. inreg .ne. 1 ) goto 997
                     if( inreg .gt. 1 ) goto 997
                  else if( imesh .eq. 3 ) then
                     if( iregd .eq. 1 .and. inxyz .ne. 1 ) goto 997
                     if( inxyz .gt. 1 ) goto 997
                  else if( imesh .eq. 4 ) then
                     if( iregd .eq. 1 .and. intet .ne. 1 ) goto 997
                     if( intet .gt. 1 ) goto 997
                  end if

                  nrsq = mrsq

                  goto 140

            else

                  nrsq = ieng + 1

                  if( imesh .eq. 1 ) then
                     imsq(1) = 1
                  else if( imesh .eq. 3)then
                     imsq(1) = 3
                  else if( imesh .eq. 4)then
                     imsq(1) = 4
                  end if

               do k = 2, nrsq

                  imsq(k) = 2
                  iwwm(k) = k - 1

               end do

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        summary of particle
*-----------------------------------------------------------------------

            if( inpat .gt. 0 ) then

               do k = 1, inpat

                  if( iptyp(k) .eq. 20 ) inpat = 0

               end do

            end if

            if( inpat .eq. 0 ) then

               do k = 1, 19

                  iptyp(k) = k

               end do

                  inpat = 19

            end if

*-----------------------------------------------------------------------
*        read weight window lower bound
*-----------------------------------------------------------------------

               nnwwp = nnwwp + 1

            if( imesh .eq. 1 .and. nnwwp .gt. kvlmax ) goto 998

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 .or.
     &              (( imsq(k) .eq. 1 .or. imsq(k) .eq. 3
     &              .or. imsq(k) . eq. 4 ) .and.
     &          iregd .gt. 1 ) ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 2 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( imesh .eq. 1 ) then

                        kdsm = kfwwp(iwwdp)
                        ldsm = iwwm(k)

                     das_kfwwp(kdsm+(ldsm-1)*kvlmax+nnwwp-1) = cvvv

                  else if( imesh .eq. 3 ) then

                        kdsm = kfwwp(iwwdp)
                        ldsm = iwwm(k)

                     das_kfwwp(kdsm+(ldsm-1)*ixyz+nnwwp-1) = cvvv

                  else if( imesh .eq. 4 ) then

                     das_kfwwp_tmp((nnwwp-1)*ieng+iwwm(k)) = cvvv

                  end if

            else if( imsq(k) .eq. 1 .and. iregd .eq. 1 ) then

                     ndsm = idsm + jdsm + 3

                  call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ndsm,ierr
     &                          ,MAX_NUM_INWWC,idas_inwwc)

                  if( ierr .ne. 0 ) goto 999

                     jdsm = jdsm + 1
                     idas_inwwc(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                     idas_inwwc(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn

            else if( imsq(k) .eq. 3 .and. iregd .eq. 1 ) then

               if( chlw(ic:ic) .eq. '(' ) then

                     ic = ic + 1
                     call snum(chlw,ic,i3,ic2,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 999

                     idas_inwwc(idsm+(nnwwp-1)*3+1-1) = nint( cvvv )

                     ic = jnumc(chlw,ic2,i3)
                     call snum(chlw,ic,i3,ic2,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 999

                     idas_inwwc(idsm+(nnwwp-1)*3+2-1) = nint( cvvv )

                     ic = jnumc(chlw,ic2,i3)
                     ib = inumc(chlw,ic,i3,')') - 1

                     call snum(chlw,ic,ib,ic2,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 999

                     idas_inwwc(idsm+(nnwwp-1)*3+3-1) = nint( cvvv )

                     ic2 = inumc(chlw,ic,i3,')') + 1

               else

                  goto 978

               end if

            else if( imsq(k) .eq. 4 .and. iregd .eq. 1 ) then

             if(nnwwp*ieng.gt.kfwwp_tmp(2))then

              call moddas_reallocate_int(
     &             2, 1, ieng*kvlmax, inwwc_tmp, idas_inwwc_tmp)

              call moddas_reallocate_dbl(
     &             2, 1, ieng*kvlmax, kfwwp_tmp, das_kfwwp_tmp)

              inwwc_tmp(1)=inwwc_tmp(2)
              kfwwp_tmp(1)=kfwwp_tmp(2)

             endif

             call snum(chlw,ic,i3,ic2,cvvv,ierr)

             if( ierr .ne. 0 ) goto 999

             idas_inwwc_tmp(nnwwp) = nint( cvvv )

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

            if( iregd .eq. 1 ) nnww0 = nnwwp
            if( iregd .gt. 1 .and. nnwwp .ne. nnww0 ) goto 983

            if( nnww0 .eq. 0 ) goto 999

               mnwwp(iwwdp, 0) = nnww0
               mnwwp(iwwdp,20) = inpat

            do j = 1, inpat

               mnwwp(iwwdp,j) = iptyp(j)

            end do

            if( iwwdp .gt. 1 .and. iwmsh .ne. imesh ) goto 976

               iwmsh = imesh

*-----------------------------------------------------------------------

            if( imesh .eq. 1 ) then

               if( jdsm > MAX_NUM_INWWC ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.wwind@read02.f ?dimension over idas_inwwc?'
     &                    //' jdsm > MAX_NUM_INWWC'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INWWC@moddas.f=',MAX_NUM_INWWC,')'
                  ErrID = 'L:4895/R:wwind/F:read02.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 6, iwwdp, jdsm+1, inwwc, idas_inwwc)

*-----------------------------------------------------------------------

            else if( imesh .eq. 3 ) then

                        if( idval .ne. 0 ) dvals = dvinp
                        ixyz = inx * iny * inz

               do j = 1, nnww0

                        idsm = inwwc(iwwdp)
                        ix = idas_inwwc(idsm+(j-1)*3+1-1)
                        iy = idas_inwwc(idsm+(j-1)*3+2-1)
                        iz = idas_inwwc(idsm+(j-1)*3+3-1)
                     icf = iz + inz * ( iy - 1 )+ iny * inz * ( ix -1 )

                  do k = 1, ieng

                        das_kgwwp(kgwwp(iwwdp)+(k-1)*ixyz+icf-1)
     &                     = das_kfwwp(kfwwp(iwwdp)+(k-1)*ixyz+j-1)

                  end do

               end do

                  do jx = 1, inx
                  do jy = 1, iny
                  do jz = 1, inz

                     icf = jz + inz * ( jy - 1 )+ iny * inz * ( jx - 1 )

                  do k = 1, ieng

                     if( das_kgwwp(kgwwp(iwwdp)+(k-1)*ixyz+icf-1)
     &                   .lt. 0.0 )
     &                   das_kgwwp(kgwwp(iwwdp)+(k-1)*ixyz+icf-1)
     &                      = dvals

                  end do
                  end do
                  end do
                  end do

            else if( imesh .eq. 4 ) then

             call moddas_reallocate_int(
     &            6, iwwdp, nnww0+mtrn+1, inwwc, idas_inwwc )

             call moddas_reallocate_dbl(
     &            6, iwwdp, nnww0*ieng, kfwwp, das_kfwwp )

             idas_inwwc(idsm+1:idsm+mtrn)
     &            =idas_temporary(1:mtrn)

             idas_inwwc(idsm+mtrn+1)=nnww0

               do j = 1, nnww0

                  idas_inwwc(idsm+mtrn+j+1)
     &                = idas_inwwc_tmp(j)

                  do k = 1, ieng

                   das_kfwwp(kfwwp(iwwdp)+(k-1)*nnww0+j-1)
     &                  = das_kfwwp_tmp((j-1)*ieng+k)

                  enddo

               enddo

               call moddas_deallocate_int( idas_temporary )

               call moddas_deallocate_int( idas_inwwc_tmp )

               call moddas_deallocate_dbl( das_kfwwp_tmp )

            end if

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

         return

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:4990/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  974 continue

         m_err = 'ireg for tet mesh should be consistent'//
     &        ' for each WW section'
         ErrCha = ''
         ErrID = 'L:5003/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  975 continue

         m_err = 'Description of mesh = tet is wrong'
         ErrCha = ''
         ErrID = 'L:5013/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  976 continue

         m_err = 'mesh should be the same for each WW section'
         ErrCha = ''
         ErrID = 'L:5025/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  977 continue

         m_err = 'Each mesh of (ix iy iz) should be smaller than 10000'
         ErrCha = ''
         ErrID = 'L:5035/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  978 continue

         m_err = 'Description of (ix iy iz) is wrong'
         ErrCha = ''
         ErrID = 'L:5045/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  979 continue

         m_err = 'mesh should be defined before part or eng'
         ErrCha = ''
         ErrID = 'L:5055/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  980 continue

         m_err = 'xyz mesh should be the same for each WW section'
         ErrCha = ''
         ErrID = 'L:5065/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  981 continue

         m_err = 'Description of mesh = xyz is wrong'
         ErrCha = ''
         ErrID = 'L:5075/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  982 continue

         m_err = 'mesh = reg, xyz, tet are not allowed at the same time'
         ErrCha = ''
         ErrID = 'L:5085/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  983 continue

         m_err = 'region definition is not consistent'
         ErrCha = ''
         ErrID = 'L:5097/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  985 continue

         m_err = 'number after ww is wrong'
         ErrCha = ''
         ErrID = 'L:5109/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  986 continue

         m_err = 'eng =, or tim = is defined twice.'
         ErrCha = ''
         ErrID = 'L:5121/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  987 continue

         m_err = 'energy or time bins is wrong'
         ErrCha = ''
         ErrID = 'L:5133/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  988 continue

         m_err = 'energy or time bins should be 0 < n < 99.'
         ErrCha = ''
         ErrID = 'L:5145/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  989 continue

         m_err = 'number of eng =, tim = is wrong'
         ErrCha = ''
         ErrID = 'L:5157/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'Number of [weight window] section is'//
     &           ' too large (=<6)'
         ErrCha = ''
         ErrID = 'L:5170/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'Number of particles is too large (<19)'
         ErrCha = ''
         ErrID = 'L:5182/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'Name of particle is wrong'
         ErrCha = ''
         ErrID = 'L:5194/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'Description of part, eng, time = is wrong'
         ErrCha = ''
         ErrID = 'L:5206/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[weight window] is wrong.'
         ErrCha = ''
         ErrID = 'L:5219/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of weight window region '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:5233/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [weight window] is wrong.'
         ErrCha = ''
         ErrID = 'L:5245/R:wwind/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine frgdat(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [frag data] section of input files                        *
*       modified by K.Niita on 2003/10/15                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)
      character frgfl*200

*-----------------------------------------------------------------------

      dimension jstyp(6), jnkf0(6)
      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 3, 4, 6 * 0 /

      character dkam*6

      logical deqn4
      logical dnen1

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'opt' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'proj' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'targ' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'file' ) then

               imsq( mrsq + 1 ) = 4
               ic = ic + 4

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 4 ) goto 997

                  inopt = 0
                  inprj = 0
                  intag = 0
                  infil = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inopt = inopt + 1
                  if( imsq(k) .eq. 2 ) inprj = inprj + 1
                  if( imsq(k) .eq. 3 ) intag = intag + 1
                  if( imsq(k) .eq. 4 ) infil = infil + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inopt .ne. 1 .or. inprj .ne. 1 .or.
     &                intag .ne. 1 .or. infil .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inopt = 1
                  inprj = 1
                  intag = 1
                  infil = 1
                  innon = 0

                  nrsq = 4

            end if

         end if

*-----------------------------------------------------------------------
*        read frag data information
*-----------------------------------------------------------------------

               ifrgd = ifrgd + 1

               if( ifrgd .gt. kvlmax ) goto 998

                  ifgdf(ifrgd,1) = 0
                  ifgdf(ifrgd,2) = 0
                  ifgdf(ifrgd,3) = 0
                  ifgdf(ifrgd,4) = 0
                  ifgdf(ifrgd,5) = 0

                  ic2  = i1

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

*-----------------------------------------------------------------------

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 1 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  ifgdf(ifrgd,1) = nint( cvvv )

               if( ifgdf(ifrgd,1) .gt. 5 .or.
     &             ifgdf(ifrgd,1) .lt. 0 ) goto 996

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 2 .or. imsq(k) .eq. 3 ) then

                     imn = imsq(k)

*-----------------------------------------------------------------------

               call rdpname(ic,i3,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

                  if( ierr .eq. 994 ) goto 984
                  if( ierr .eq. 998 ) goto 999

                  if( istyp .lt. 0 ) then

                     istyp = jstyp(1)
                     inkf0 = jnkf0(1)

                  end if

*-----------------------------------------------------------------------

                  ic2 = ic

CS.Hashimoto revised for extension of kind of particles. (2014.9.24)
              if ( imn .eq. 3 ) then ! target should be proton or nucleus

               if( abs(inkf0) .le. 1000000 .and.
     &             abs(inkf0) .ne. 2212 ) goto 989

              end if

                  ifgdf(ifrgd,imn) = inkf0

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 4 ) then

                  ic1 = ic
                  ic2 = min( i3, inumc(chlw,ic1+1,i3,' ') )

                  iname = ic2 - ic1 + 1

               do i = 1, iname

                  frgfl(ifrgd)(i:i) = chin(ic1+i-1:ic1+i-1)

               end do

               do i = iname + 1, 200

                  frgfl(ifrgd)(i:i) = ' '

               end do

                  ifgdf(ifrgd,4) = iname

                  ic = ic2

*-----------------------------------------------------------------------

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  984    m_err = 'Name of proj or targ is wrong in [frag data]'
         ErrCha = ''
         ErrID = 'L:5547/R:frgdat/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  989    m_err = 'proj or targ shloud be nucleus in [frag data]'
         ErrCha = ''
         ErrID = 'L:5557/R:frgdat/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'opt should be 0, 1, 2, 3, 4, 5 in [frag data].'
         ErrCha = ''
         ErrID = 'L:5569/R:frgdat/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[frag data] is wrong.'
         ErrCha = ''
         ErrID = 'L:5582/R:frgdat/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of frag data '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:5596/R:frgdat/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [frag data] is wrong.'
         ErrCha = ''
         ErrID = 'L:5608/R:frgdat/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine mattc(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [mat time change] section of input files                  *
*       modified by K.Niita on 2004/06/22                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /mttmc/  smttc(kvlmax), mttcn, mttc1(kvlmax), mttc2(kvlmax)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 3, 7 * 0 /

      character dkam*6

      logical deqn4
      logical dnen1

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'mat' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'time' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 4

            else if( chlw(ic:ic+5) .eq. 'change' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 6

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inmat = 0
                  intim = 0
                  incmt = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inmat = inmat + 1
                  if( imsq(k) .eq. 2 ) intim = intim + 1
                  if( imsq(k) .eq. 3 ) incmt = incmt + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inmat .ne. 1 .or.
     &                intim .ne. 1 .or. incmt .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inmat = 1
                  intim = 1
                  incmt = 1
                  innon = 0

                  nrsq = 3

            end if

         end if

*-----------------------------------------------------------------------
*        read time and change information
*-----------------------------------------------------------------------

               mttcn = mttcn + 1

               if( mttcn .gt. kvlmax ) goto 998

                  mttc1(mttcn) = 0
                  mttc2(mttcn) = 0
                  smttc(mttcn) = -1.0d0

                  ic2  = i1

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 1 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  mttc1(mttcn) = nint( cvvv )

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 2 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  smttc(mttcn) = cvvv

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 3 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  mttc2(mttcn) = nint( cvvv )

*-----------------------------------------------------------------------

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[mat time change] is wrong.'
         ErrCha = ''
         ErrID = 'L:5857/R:mattc/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of mat time change '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:5871/R:mattc/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [mat time change] is wrong.'
         ErrCha = ''
         ErrID = 'L:5883/R:mattc/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine matnc(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [mat name color] section of input files                   *
*       modified by K.Niita on 2003/04/24                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /mtnmc/  smtnc(kvlmax), dmtnc(kvlmax,2),
     &                mtncn, mtnc(kvlmax,2), nmtnc(kvlmax,2)
      character dmtnc*80

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 3, 7 * 0 /

      character dkam*6

      logical deqn4
      logical dnen1

*-----------------------------------------------------------------------
      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'mat' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'name' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 4

            else if( chlw(ic:ic+4) .eq. 'color' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 5

            else if( chlw(ic:ic+3) .eq. 'size' ) then

               imsq( mrsq + 1 ) = 4
               ic = ic + 5

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inmat = 0
                  innam = 0
                  incol = 0
                  insiz = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inmat = inmat + 1
                  if( imsq(k) .eq. 2 ) innam = innam + 1
                  if( imsq(k) .eq. 3 ) incol = incol + 1
                  if( imsq(k) .eq. 4 ) insiz = insiz + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inmat .ne. 1 .or.
     &              ( innam .ne. 1 .and. incol .ne. 1 .and.
     &                insiz .ne. 1 ) .or.
     &                innam .gt. 1 .or. incol .gt. 1 .or.
     &                insiz .gt. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inmat = 1
                  innam = 1
                  incol = 1
                  insiz = 0
                  innon = 0

                  nrsq = 3

            end if

         end if

*-----------------------------------------------------------------------
*        read name, color and size information
*-----------------------------------------------------------------------

               mtncn = mtncn + 1

               if( mtncn .gt. kvlmax ) goto 998

                  nmtnc(mtncn,1) = 0
                  nmtnc(mtncn,2) = 0
                  smtnc(mtncn) = 1.0d0

                  ic2  = i1

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 1 ) then

                        ilev = 0
                        ic = ic - 1

  180                ic = ic + 1
                     if( ic .gt. i3 ) goto 190

                     if( chlw(ic:ic) .eq. '{' ) then

                        if( ilev .gt. 0 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '-' ) then

                        if( ilev .ne. 1 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '}' ) then

                        if( ilev .ne. 2 ) goto 995
                        ilev = 0

                        if( ntrf .le. ntri ) goto 995

                        mtnc(mtncn,1) = ntri
                        mtnc(mtncn,2) = ntrf

                        goto 190

                     else if( deqn4( chlw(ic:ic) ) ) then

                           ici = ic

                        do j = ic + 1, i3

                           if( dnen1( chlw(j:j) ) ) goto 167

                        end do

  167                      icf = j - 1

                           if( icf .lt. ici ) goto 996

                        call onum(chlw,ici,icf,cvvv,ierr)

                           if( ierr .ne. 0 ) goto 996

                           ic = icf

                           if( nint(cvvv) .lt. 0 .or.
     &                         nint(cvvv) .ge. kvmmax ) goto 993

*-----------------------------------------------------------------------

                        if( ilev .eq. 0 ) then

                           mtnc(mtncn,1) = nint( cvvv )
                           mtnc(mtncn,2) = nint( cvvv )

                           goto 190

                        else if( ilev .eq. 1 ) then

                           ntri = nint( cvvv )

                        else if( ilev .eq. 2 ) then

                           ntrf = nint( cvvv )

                        end if

*-----------------------------------------------------------------------

                     end if

                           goto 180

  190                continue

                     ic2 = ic + 1

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 2 .or. imsq(k) .eq. 3 ) then

                  if( imsq(k) .eq. 2 ) then
                     im = 1
                  else
                     im = 2
                  end if

                        ilev = 0
                        ntri = 0
                        ntrf = 0
                        ic = ic - 1

  181                ic = ic + 1
                     if( ic .gt. i3 ) goto 191

                     if( ntri .eq. 0 .and.
     &                   chlw(ic:ic) .eq. ' ' ) then

                        goto 181

                     else if( ntri .eq. 0 .and. ic .eq. i3 ) then

                        ntri = ic
                        ntrf = ic

                        goto 191

                     else if( ntri .gt. 0 .and. ilev .eq. 0 .and.
     &                      ( chlw(ic:ic) .eq. ' ' .or.
     &                        ic .eq. i3 ) ) then

                        ntrf = ic
                        if( ic .ne. i3 ) ntrf = ic - 1

                        goto 191

                     else if( ntri .gt. 0 .and. ilev .gt. 0 .and.
     &                        chlw(ic:ic) .eq. ' ' ) then

                        goto 181

                     else if( chlw(ic:ic) .eq. '{' .and.
     &                        chlw(ic-1:ic-1) .ne. yen ) then

                        if( ilev .gt. 0 ) goto 994
                        ilev = ilev + 1
                        ntri = ic + 1

                     else if( chlw(ic:ic) .eq. '}' .and.
     &                        chlw(ic-1:ic-1) .ne. yen ) then

                        if( ilev .ne. 1 ) goto 994
                        ilev = 0

                        ntrf = ic - 1

                        goto 191

                     else if( ntri .eq. 0 ) then

                        ntri = ic

*-----------------------------------------------------------------------

                     end if

                           goto 181

  191                continue

                        if( ntrf .lt. ntri ) goto 994
                        if( ntrf - ntri + 1 .gt. 80 ) goto 992

                     ii = 0
                  do 192 jj = ntri, ntrf
                     if( chlw(jj:jj) .eq. yen .and.
     &                   chlw(jj+1:jj+1) .eq. '{' ) goto 192
                     if( chlw(jj:jj) .eq. yen .and.
     &                   chlw(jj+1:jj+1) .eq. '}' ) goto 192
                        ii = ii + 1
                        dmtnc(mtncn,im)(ii:ii) = chin(jj:jj)
  192             continue
                        nmtnc(mtncn,im) = ii

                     if( im .eq. 2 .and.
     &                    ntrf - ntri + 1 .eq. 7 .and.
     &                    chlw(ntri:ntrf) .eq. 'default' ) then

                        nmtnc(mtncn,im) = 0

                     end if

                     ic2 = ic + 1

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 4 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  smtnc(mtncn) = cvvv

*-----------------------------------------------------------------------

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  992 continue

         m_err = 'definition of name or color is too long'
         ErrCha = ''
         ErrID = 'L:6309/R:matnc/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'material number should be 1 - kvmmax-1'
         ErrCha = ''
         ErrID = 'L:6321/R:matnc/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'Description of material name or color '//
     &           'is wrong.'
         ErrCha = ''
         ErrID = 'L:6334/R:matnc/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of material number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:6347/R:matnc/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Description of material number '//
     &           '{n1-n2} or n1 is wrong.'
         ErrCha = ''
         ErrID = 'L:6360/R:matnc/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[mat name color] is wrong.'
         ErrCha = ''
         ErrID = 'L:6373/R:matnc/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of mat name color '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:6387/R:matnc/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [mat name color] is wrong.'
         ErrCha = ''
         ErrID = 'L:6399/R:matnc/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine elasop(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [elastic option] section of input files                   *
*       modified by K.Niita on 2005/04/05                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /elreg/  elarg(kvlmax,4), mlrgn, melrg(kvlmax,2), ielusr

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 9 * 0 /

      character dkam*6

      logical deqn4
      logical dnen1

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+1) .eq. 'c1' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 2

            else if( chlw(ic:ic+1) .eq. 'c2' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 2

            else if( chlw(ic:ic+1) .eq. 'c3' ) then

               imsq( mrsq + 1 ) = 4
               ic = ic + 2

            else if( chlw(ic:ic+1) .eq. 'c4' ) then

               imsq( mrsq + 1 ) = 5
               ic = ic + 2

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 1 ) goto 997

                  inmat = 0
                  incc1 = 0
                  incc2 = 0
                  incc3 = 0
                  incc4 = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inmat = inmat + 1
                  if( imsq(k) .eq. 2 ) incc1 = incc1 + 1
                  if( imsq(k) .eq. 3 ) incc2 = incc2 + 1
                  if( imsq(k) .eq. 4 ) incc3 = incc3 + 1
                  if( imsq(k) .eq. 5 ) incc4 = incc4 + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inmat .ne. 1 .or.
     &                incc1 .gt. 1 .or. incc2 .gt. 1 .or.
     &                incc3 .gt. 1 .or. incc4 .gt. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inmat = 1
                  incc1 = 0
                  incc2 = 0
                  incc3 = 0
                  incc4 = 0
                  innon = 0

                  nrsq = 1

            end if

         end if

*-----------------------------------------------------------------------
*        read name and size information
*-----------------------------------------------------------------------

               mlrgn = mlrgn + 1

               if( mlrgn .gt. kvlmax ) goto 998

                  elarg(mlrgn,1) = 0.0d0
                  elarg(mlrgn,2) = 0.0d0
                  elarg(mlrgn,3) = 0.0d0
                  elarg(mlrgn,4) = 0.0d0

                  ic2  = i1

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 1 ) then

                        ilev = 0
                        ic = ic - 1

  180                ic = ic + 1
                     if( ic .gt. i3 ) goto 190

                     if( chlw(ic:ic) .eq. '{' ) then

                        if( ilev .gt. 0 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '-' ) then

                        if( ilev .ne. 1 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '}' ) then

                        if( ilev .ne. 2 ) goto 995
                        ilev = 0

                        if( ntrf .le. ntri ) goto 995

                        melrg(mlrgn,1) = ntri
                        melrg(mlrgn,2) = ntrf

                        goto 190

                     else if( deqn4( chlw(ic:ic) ) ) then

                           ici = ic

                        do j = ic + 1, i3

                           if( dnen1( chlw(j:j) ) ) goto 167

                        end do

  167                      icf = j - 1

                           if( icf .lt. ici ) goto 996

                        call onum(chlw,ici,icf,cvvv,ierr)

                           if( ierr .ne. 0 ) goto 996

                           ic = icf

                           if( nint(cvvv) .lt. 0 .or.
     &                         nint(cvvv) .ge. kvmmax ) goto 993

*-----------------------------------------------------------------------

                        if( ilev .eq. 0 ) then

                           melrg(mlrgn,1) = nint( cvvv )
                           melrg(mlrgn,2) = nint( cvvv )

                           goto 190

                        else if( ilev .eq. 1 ) then

                           ntri = nint( cvvv )

                        else if( ilev .eq. 2 ) then

                           ntrf = nint( cvvv )

                        end if

*-----------------------------------------------------------------------

                     end if

                           goto 180

  190                continue

                     ic2 = ic + 1

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 2 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  elarg(mlrgn,1) = cvvv

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 3 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  elarg(mlrgn,2) = cvvv

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 4 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  elarg(mlrgn,3) = cvvv

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 5 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  elarg(mlrgn,4) = cvvv

*-----------------------------------------------------------------------

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  993 continue

         m_err = 'region or cell number should be 1 - kvmmax-1'
         ErrCha = ''
         ErrID = 'L:6759/R:elasop/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region or cell number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:6772/R:elasop/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Description of reigon or cell number '//
     &           '{n1-n2} or n1 is wrong.'
         ErrCha = ''
         ErrID = 'L:6785/R:elasop/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[elatic option] is wrong.'
         ErrCha = ''
         ErrID = 'L:6798/R:elasop/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of elastic option '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:6812/R:elasop/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [elastic option] is wrong.'
         ErrCha = ''
         ErrID = 'L:6824/R:elasop/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine regnm(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [reg name] section of input files                         *
*       modified by K.Niita on 2003/04/24                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /mtreg/  smtrg(kvlmax), dmtrg(kvlmax),
     &                mtrgn, mtrg(kvlmax,2), nmtrg(kvlmax)
      character dmtrg*80

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 8 * 0 /

      character dkam*6

      logical deqn4
      logical dnen1

*-----------------------------------------------------------------------
      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'name' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'size' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 5

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inmat = 0
                  innam = 0
                  insiz = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inmat = inmat + 1
                  if( imsq(k) .eq. 2 ) innam = innam + 1
                  if( imsq(k) .eq. 3 ) insiz = insiz + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inmat .ne. 1 .or.
     &              ( innam .ne. 1 .and. insiz .ne. 1 ) .or.
     &                innam .gt. 1 .or.  insiz .gt. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inmat = 1
                  innam = 1
                  insiz = 0
                  innon = 0

                  nrsq = 2

            end if

         end if

*-----------------------------------------------------------------------
*        read name and size information
*-----------------------------------------------------------------------

               mtrgn = mtrgn + 1

               if( mtrgn .gt. kvlmax ) goto 998

                  nmtrg(mtrgn) = 0
                  smtrg(mtrgn) = 1.0d0

                  ic2  = i1

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 1 ) then

                        ilev = 0
                        ic = ic - 1

  180                ic = ic + 1
                     if( ic .gt. i3 ) goto 190

                     if( chlw(ic:ic) .eq. '{' ) then

                        if( ilev .gt. 0 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '-' ) then

                        if( ilev .ne. 1 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '}' ) then

                        if( ilev .ne. 2 ) goto 995
                        ilev = 0

                        if( ntrf .le. ntri ) goto 995

                        mtrg(mtrgn,1) = ntri
                        mtrg(mtrgn,2) = ntrf

                        goto 190

                     else if( deqn4( chlw(ic:ic) ) ) then

                           ici = ic

                        do j = ic + 1, i3

                           if( dnen1( chlw(j:j) ) ) goto 167

                        end do

  167                      icf = j - 1

                           if( icf .lt. ici ) goto 996

                        call onum(chlw,ici,icf,cvvv,ierr)

                           if( ierr .ne. 0 ) goto 996

                           ic = icf

                           if( nint(cvvv) .lt. 0 .or.
     &                         nint(cvvv) .ge. kvmmax ) goto 993

*-----------------------------------------------------------------------

                        if( ilev .eq. 0 ) then

                           mtrg(mtrgn,1) = nint( cvvv )
                           mtrg(mtrgn,2) = nint( cvvv )

                           goto 190

                        else if( ilev .eq. 1 ) then

                           ntri = nint( cvvv )

                        else if( ilev .eq. 2 ) then

                           ntrf = nint( cvvv )

                        end if

*-----------------------------------------------------------------------

                     end if

                           goto 180

  190                continue

                     ic2 = ic + 1

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 2 ) then

                        ilev = 0
                        ntri = 0
                        ntrf = 0
                        ic = ic - 1

  181                ic = ic + 1
                     if( ic .gt. i3 ) goto 191

                     if( ntri .eq. 0 .and.
     &                   chlw(ic:ic) .eq. ' ' ) then

                        goto 181

                     else if( ntri .eq. 0 .and. ic .eq. i3 ) then

                        ntri = ic
                        ntrf = ic

                        goto 191

                     else if( ntri .gt. 0 .and. ilev .eq. 0 .and.
     &                      ( chlw(ic:ic) .eq. ' ' .or.
     &                        ic .eq. i3 ) ) then

                        ntrf = ic
                        if( ic .ne. i3 ) ntrf = ic - 1

                        goto 191

                     else if( ntri .gt. 0 .and. ilev .gt. 0 .and.
     &                        chlw(ic:ic) .eq. ' ' ) then

                        goto 181

                     else if( chlw(ic:ic) .eq. '{' .and.
     &                        chlw(ic-1:ic-1) .ne. yen ) then

                        if( ilev .gt. 0 ) goto 994
                        ilev = ilev + 1
                        ntri = ic + 1

                     else if( chlw(ic:ic) .eq. '}' .and.
     &                        chlw(ic-1:ic-1) .ne. yen ) then

                        if( ilev .ne. 1 ) goto 994
                        ilev = 0

                        ntrf = ic - 1

                        goto 191

                     else if( ntri .eq. 0 ) then

                        ntri = ic

*-----------------------------------------------------------------------

                     end if

                           goto 181

  191                continue

                        if( ntrf .lt. ntri ) goto 994
                        if( ntrf - ntri + 1 .gt. 80 ) goto 992

                     ii = 0
                  do 192 jj = ntri, ntrf
                     if( chlw(jj:jj) .eq. yen .and.
     &                   chlw(jj+1:jj+1) .eq. '{' ) goto 192
                     if( chlw(jj:jj) .eq. yen .and.
     &                   chlw(jj+1:jj+1) .eq. '}' ) goto 192
                        ii = ii + 1
                        dmtrg(mtrgn)(ii:ii) = chin(jj:jj)
  192             continue
                        nmtrg(mtrgn) = ii

                     ic2 = ic + 1

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 3 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  smtrg(mtrgn) = cvvv

*-----------------------------------------------------------------------

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  992 continue

         m_err = 'definition of name is too long'
         ErrCha = ''
         ErrID = 'L:7225/R:regnm/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'region or cell number should be 1 - kvmmax-1'
         ErrCha = ''
         ErrID = 'L:7237/R:regnm/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'Description of reg name or size '//
     &           'is wrong.'
         ErrCha = ''
         ErrID = 'L:7250/R:regnm/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region or cell number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:7263/R:regnm/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Description of reigon or cell number '//
     &           '{n1-n2} or n1 is wrong.'
         ErrCha = ''
         ErrID = 'L:7276/R:regnm/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[reg name] is wrong.'
         ErrCha = ''
         ErrID = 'L:7289/R:regnm/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of reg name '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:7303/R:regnm/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [reg name] is wrong.'
         ErrCha = ''
         ErrID = 'L:7315/R:regnm/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine timers(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [timer] section of input files                            *
*       modified by K.Niita on 2005/12/02                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /tmtreg/ ntmrg, intmc, intmt, ktime

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      character dkam*6

      dimension itims(4,kvlmax)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0
            ntmrg = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 .and. ntmrg .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+1) .eq. 'in' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 2

            else if( chlw(ic:ic+2) .eq. 'out' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'coll' ) then

               imsq( mrsq + 1 ) = 4
               ic = ic + 4

            else if( chlw(ic:ic+2) .eq. 'ref' ) then

               imsq( mrsq + 1 ) = 5
               ic = ic + 3

            else

               goto 997

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inreg = 0
                  inin  = 0
                  inout = 0
                  incol = 0
                  inref = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) inin  = inin  + 1
                  if( imsq(k) .eq. 3 ) inout = inout + 1
                  if( imsq(k) .eq. 4 ) incol = incol + 1
                  if( imsq(k) .eq. 5 ) inref = inref + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inreg .ne. 1 ) goto 997
                  if( inin  .gt. 1 .or. inout .gt. 1 .or.
     &                incol .gt. 1 .or. inref .gt. 1 ) goto 997
                  if( inin  .eq. 0 .and. inout .eq. 0 .and.
     &                incol .eq. 0 .and. inref .eq. 0 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inreg = 1
                  inin  = 1
                  inout = 1
                  incol = 1
                  inref = 1
                  innon = 0

                  imsq(1) = 1
                  imsq(2) = 2
                  imsq(3) = 3
                  imsq(4) = 4
                  imsq(5) = 5

                  nrsq = 5

            end if

         end if

*-----------------------------------------------------------------------
*        read timer informations
*-----------------------------------------------------------------------

               ntmrg = ntmrg + 1

               if( ntmrg .gt. kvlmax ) goto 998

                  itims(1,ntmrg) = 0
                  itims(2,ntmrg) = 0
                  itims(3,ntmrg) = 0
                  itims(4,ntmrg) = 0

            if( ntmrg .eq. 1 ) then

                  intmc = 0
                  iaddress_region(:) = 0
                  call moddas_reallocate_int(
     &                    2, 1, MAX_NUM_INTMC, iaddress_region
     &                    , idas_intmc)
                  idsm  = intmc
                  jdsm  = 0

            end if

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 1 ) then

                     ndsm = idsm + jdsm + 3

                  call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ndsm,ierr
     &                          ,MAX_NUM_INTMC,idas_intmc)

                  if( ierr .ne. 0 ) goto 999

                     jdsm = jdsm + 1
                     idas_intmc(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                     idas_intmc(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn

            else if( imsq(k) .eq. 2 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999
               if( nint(cvvv) .ne. -1 .and.
     &             nint(cvvv) .ne.  0 .and.
     &             nint(cvvv) .ne.  1 ) goto 994

               itims(1,ntmrg) = nint(cvvv)

            else if( imsq(k) .eq. 3 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999
               if( nint(cvvv) .ne. -1 .and.
     &             nint(cvvv) .ne.  0 .and.
     &             nint(cvvv) .ne.  1 ) goto 994

               itims(2,ntmrg) = nint(cvvv)

            else if( imsq(k) .eq. 4 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999
               if( nint(cvvv) .ne. -1 .and.
     &             nint(cvvv) .ne.  0 .and.
     &             nint(cvvv) .ne.  1 ) goto 994

               itims(3,ntmrg) = nint(cvvv)

            else if( imsq(k) .eq. 5 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999
               if( nint(cvvv) .ne. -1 .and.
     &             nint(cvvv) .ne.  0 .and.
     &             nint(cvvv) .ne.  1 ) goto 994

               itims(4,ntmrg) = nint(cvvv)

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

               if( jdsm > MAX_NUM_INTMC ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.timers@read02.f ?dimension over idas_intmc?'
     &                    //' jdsm > MAX_NUM_INTMC'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INTMC@moddas.f=',MAX_NUM_INTMC,')'
                  ErrID = 'L:7638/R:timers/F:read02.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 2, 1, jdsm+1, iaddress_region, idas_intmc)

               ktime = 0
               call moddas_allocate_int(ntmrg*4+2, idas_ktime)
               idsm  = ktime
               jdsm  = 0


            do i = 1, ntmrg

               jdsm = jdsm + 1
               idas_ktime(idsm+jdsm) = itims(1,i)
               jdsm = jdsm + 1
               idas_ktime(idsm+jdsm) = itims(2,i)
               jdsm = jdsm + 1
               idas_ktime(idsm+jdsm) = itims(3,i)
               jdsm = jdsm + 1
               idas_ktime(idsm+jdsm) = itims(4,i)

            end do

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:7675/R:timers/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'timer id should be -1 or 0 or 1'
         ErrCha = ''
         ErrID = 'L:7687/R:timers/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[timer] is wrong.'
         ErrCha = ''
         ErrID = 'L:7700/R:timers/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of region for timer '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:7714/R:timers/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [timer] is wrong.'
         ErrCha = ''
         ErrID = 'L:7726/R:timers/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine counts(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [counter] section of input files                          *
*       modified by K.Niita on 2015/12/01                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,mxcntprt,2) ! S.H. set mxcntprt (2022.3.24)

      common /cntech/ icnech(3,0:maxcntr+1) ! S.H. set maxcntr+1(non case) (2022.3.24)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(maxcntr+1) ! S.H. set maxcntr+1(non case) (2022.3.24)

      character dkam*6

      dimension icnts(3,maxcntr-1,kvlmax)   ! S.H. set maxcntr-1(reg case) (2022.3.24)

      dimension     idas(1)
      equivalence ( das, idas )

      dimension iptyp(mxcntprt), inkfp(mxcntprt) ! S.H. set mxcntprt (2022.3.24)
      dimension jstyp(6), jnkf0(6)

      common /subtra/ isubt, ipsub(mxpart)   ! kitamura22/03/31

      dimension icontind(3), idasadress(3) ! frtati 2024/05/21

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 1
            ncsq  = 0
            nmreg = 0
            inpat = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        counter number
*-----------------------------------------------------------------------

            ic = i1

         if( chlw(ic:ic+6) .eq. 'counter' ) then

            if( ncsq .gt. 0 ) then

               if( nmreg .eq. 0 ) goto 993

               ncreg(icont) = nmreg

               mmmax = mmmax + ( jdsm + mod(jdsm,2) ) / 2 + 1
               if( mmmax .gt. mdas ) goto 970

!               call moddas_reallocate_int(
!     &                 3, icont, nmreg*(maxcntr-1)+2, kcont, idas_kcont) ! T.Sato 2024/03/21
               call moddas_reallocate_int(
     &                 3, ncsq, nmreg*(maxcntr-1)+2, kcont, idas_kcont) ! frtati 2024/05/21
!               idsm  = kcont(icont)
               idsm  = kcont(ncsq) ! frtati 2024/05/21
               jdsm  = 0


               do i = 1, nmreg

                  do j = 1, maxcntr-1   ! S.H. set maxcntr-1(reg case) (2022.3.24)
                     jdsm = jdsm + 1
                     idas_kcont(idsm+jdsm) = icnts(icont,j,i)
                  enddo

               end do

               inpat = 0

            end if

               ic = inumc(chlw,ic+7,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 999

               call onum(chlw,ic,i3,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               icont = nint( cvvv )

               if( icont .lt. 1 .or. icont .gt. 3 ) goto 996
               if( ncntc(icont) .ne. 0 ) goto 995

               ncntc(icont) = 1

               ncsq  = ncsq + 1
               nrsq  = 0
               nmreg = 0

               icontind(ncsq) = icont ! frtati 2024/05/21

               goto 140

         end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 .and. nmreg .eq. 0 .and. ncsq .gt. 0 ) then

               ic  = i1
               icl = i1

               chlc = chlw
               call chcomp(chlc,icl,i3,i5)

            if(chlc(icl:icl+3) .eq. 'part' .or.
     &         chlc(icl:icl+4) .eq. '*part' ) then

               if(chlc(icl:icl+3) .eq. 'part' ) then

                  inppn = 1
                  ic = inumc(chlw,icl+4,i3,'=') + 1

               else if(chlc(icl:icl+4) .eq. '*part' ) then

                  inppn = -1
                  ic = inumc(chlw,icl+5,i3,'=') + 1

               end if

               if( ic .gt. i3 ) goto 984

               icl = i3

  400          continue

                  ic = jnumc(chlw,ic,i3)

                  if( ic .gt. i3 ) goto 140

*-----------------------------------------------------------------------

                  call rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

                     if( ierr .eq. 994 ) goto 984
                     if( ierr .eq. 998 ) goto 984
                     if( isubt .eq. 1 )  goto 984   ! kitamura22/03/31

*-----------------------------------------------------------------------

                  if( istyp .gt. 0 ) then

                        inpat = inpat + 1
                        if( inpat .gt. mxcntprt ) goto 983 ! S.H. set mxcntprt (2022.3.24)

                        iptyp(inpat) = istyp
                        inkfp(inpat) = inkf0

                  else if( istyp .lt. 0 ) then

                     do i = 1, -istyp

                        inpat = inpat + 1
                        if( inpat .gt. mxcntprt ) goto 983 ! S.H. set mxcntprt (2022.3.24)

                        iptyp(inpat) = jstyp(i)
                        inkfp(inpat) = jnkf0(i)

                     end do

                  end if

*-----------------------------------------------------------------------

                  goto 400

            end if

*-----------------------------------------------------------------------
*           particles
*-----------------------------------------------------------------------

            if( inpat .gt. 0 ) then

                  icpan(icont) = inpat * inppn

               do i = 1, inpat

                  icpat(icont,i,1) = iptyp(i)
                  icpat(icont,i,2) = inkfp(i)

               end do

            else if( inpat .eq. 0 ) then

                  icpan(icont) = 1
                  icpat(icont,1,1) = mxcntprt ! S.H. set mxcntprt (2022.3.24)
                  icpat(icont,1,2) = 0

            end if

         end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 .and. nmreg .eq. 0 .and. ncsq .gt. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+1) .eq. 'in' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 2

            else if( chlw(ic:ic+2) .eq. 'out' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'coll' ) then

               imsq( mrsq + 1 ) = 4
               ic = ic + 4

            else if( chlw(ic:ic+2) .eq. 'ref' ) then

               imsq( mrsq + 1 ) = 5
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'fiss' ) then

               imsq( mrsq + 1 ) = 6
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'elst' ) then

               imsq( mrsq + 1 ) = 7
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'iels' ) then

               imsq( mrsq + 1 ) = 8
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'nucl' ) then

               imsq( mrsq + 1 ) = 9
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'dcay' ) then

               imsq( mrsq + 1 ) = 10
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'atom' ) then

               imsq( mrsq + 1 ) = 11
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'delt' ) then

               imsq( mrsq + 1 ) = 12
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'fluo' ) then

               imsq( mrsq + 1 ) = 13
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'auge' ) then

               imsq( mrsq + 1 ) = 14
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'brem' ) then

               imsq( mrsq + 1 ) = 15
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'phel' ) then

               imsq( mrsq + 1 ) = 16
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'cmpt' ) then

               imsq( mrsq + 1 ) = 17
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'pprd' ) then

               imsq( mrsq + 1 ) = 18
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'anih' ) then

               imsq( mrsq + 1 ) = 19
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'msct' ) then

               imsq( mrsq + 1 ) = 20
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'rayl' ) then

               imsq( mrsq + 1 ) = 21
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'knoe' ) then   ! S.Abe 2018/02/26

               imsq( mrsq + 1 ) = 22
               ic = ic + 4

            else if( chlw(ic:ic+4) .eq. 'ndata' ) then   ! S.H. 2021.12.16

               imsq( mrsq + 1 ) = 23
               ic = ic + 5

            else

               goto 997

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inreg = 0
                  inin  = 0
                  inout = 0
                  incol = 0
                  inref = 0
                  infis = 0
                  innon = 0
                  inels = 0
                  iniel = 0
                  inncr = 0
                  indcy = 0
                  inato = 0
                  indlr = 0
                  influ = 0
                  inaug = 0
                  inbrm = 0
                  inphe = 0
                  incmp = 0
                  inppd = 0
                  inanh = 0
                  inmst = 0
                  inray = 0
                  inkoe = 0   ! S.Abe 2018/02/26
                  indat = 0   ! S.H. 2021.12.16

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) inin  = inin  + 1
                  if( imsq(k) .eq. 3 ) inout = inout + 1
                  if( imsq(k) .eq. 4 ) incol = incol + 1
                  if( imsq(k) .eq. 5 ) inref = inref + 1
                  if( imsq(k) .eq. 6 ) infis = infis + 1
                  if( imsq(k) .eq. 7 ) inels = inels + 1
                  if( imsq(k) .eq. 8 ) iniel = iniel + 1
                  if( imsq(k) .eq. 9 ) inncr = inncr + 1
                  if( imsq(k) .eq. 10) indcy = indcy + 1
                  if( imsq(k) .eq. 11) inato = inato + 1
                  if( imsq(k) .eq. 12) indlr = indlr + 1
                  if( imsq(k) .eq. 13) influ = influ + 1
                  if( imsq(k) .eq. 14) inaug = inaug + 1
                  if( imsq(k) .eq. 15) inbrm = inbrm + 1
                  if( imsq(k) .eq. 16) inphe = inphe + 1
                  if( imsq(k) .eq. 17) incmp = incmp + 1
                  if( imsq(k) .eq. 18) inppd = inppd + 1
                  if( imsq(k) .eq. 19) inanh = inanh + 1
                  if( imsq(k) .eq. 20) inmst = inmst + 1
                  if( imsq(k) .eq. 21) inray = inray + 1
                  if( imsq(k) .eq. 22) inkoe = inkoe + 1   ! S.Abe 2018/02/26
                  if( imsq(k) .eq. 23) indat = indat + 1   ! S.H. 2021.12.16
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inreg .ne. 1 ) goto 997
                  if( inin  .gt. 1 .or. inout .gt. 1 .or.
     &                incol .gt. 1 .or. inref .gt. 1 .or.
     &                infis .gt. 1 .or. inels .gt. 1 .or.
     &                iniel .gt. 1 .or. inncr .gt. 1 .or.
     &                indcy .gt. 1 .or. inato .gt. 1 .or.
     &                indlr .gt. 1 .or. influ .gt. 1 .or.
     &                inaug .gt. 1 .or. inbrm .gt. 1 .or.
     &                inphe .gt. 1 .or. incmp .gt. 1 .or.
     &                inppd .gt. 1 .or. inanh .gt. 1 .or.
     &                inmst .gt. 1 .or. inray .gt. 1 .or.
     &                inkoe .gt. 1 .or. indat .gt. 1 ) goto 997   ! S.H. 2021.12.16
                  if( inin  .eq. 0 .and. inout .eq. 0 .and.
     &                incol .eq. 0 .and. inref .eq. 0 .and.
     &                infis .eq. 0 .and. inels .eq. 0 .and.
     &                iniel .eq. 0 .and. inncr .eq. 0 .and.
     &                indcy .eq. 0 .and. inato .eq. 0 .and.
     &                indlr .eq. 0 .and. influ .eq. 0 .and.
     &                inaug .eq. 0 .and. inbrm .eq. 0 .and.
     &                inphe .eq. 0 .and. incmp .eq. 0 .and.
     &                inppd .eq. 0 .and. inanh .eq. 0 .and.
     &                inmst .eq. 0 .and. inray .eq. 0 .and.
     &                inkoe .eq. 0 .and. indat .eq. 0 ) goto 997   ! S.H. 2021.12.16

                  nrsq = mrsq

               icnech(icont,0) = nrsq
               do k = 1, nrsq
                  icnech(icont,k) = imsq(k)
               enddo

                  goto 140

            else

                  inreg = 1
                  inin  = 1
                  inout = 1
                  incol = 1
                  inref = 1
                  infis = 1
                  innon = 0
                  inels = 1
                  iniel = 1
                  inncr = 1
                  indcy = 1
                  inato = 1
                  indlr = 1
                  influ = 1
                  inaug = 1
                  inbrm = 1
                  inphe = 1
                  incmp = 1
                  inppd = 1
                  inanh = 1
                  inmst = 1
                  inray = 1
                  inkoe = 1   ! S.Abe 2018/02/26
                  indat = 1   ! S.H. 2021.12.16

               do k = 1, maxcntr+1   ! S.H. set maxcntr+1(non case) (2022.3.24)
                  imsq(k) = k
               enddo

                  nrsq = maxcntr+1   ! S.H. set maxcntr+1(non case) (2022.3.24)

               icnech(icont,0) = 5
               do k = 1, 5
                  icnech(icont,k) = imsq(k)
               enddo

            end if

         end if

         if( ncsq .eq. 0 ) goto 999

*-----------------------------------------------------------------------
*        read counter informations
*-----------------------------------------------------------------------

               nmreg = nmreg + 1

               if( nmreg .gt. kvlmax ) goto 998

               do j = 1, maxcntr-1   ! S.H. set maxcntr-1(reg case) (2022.3.24)
                  icnts(icont,j,nmreg) = 0
               enddo

            if( nmreg .eq. 1 ) then

cfrtati 2024/05/21
!                  call moddas_reallocate_int(
!     &                    3, icont, MAX_NUM_INCRC, incrc, idas_incrc)
                  if( ncsq.gt.1 ) then
                    call moddas_reduce_int(
     &              3, ncsq-1, idasadress(ncsq-1)+1, incrc, idas_incrc)
                  end if
                  call moddas_reallocate_int(
     &            3, ncsq, MAX_NUM_INCRC, incrc, idas_incrc)
!                  idsm  = incrc(icont)
                  idsm  = incrc(ncsq) ! frtati 2024/05/21
                  jdsm  = 0

            end if

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 1 ) then

                     ndsm = idsm + jdsm + 3

                  call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ndsm,ierr
     &                          ,MAX_NUM_INCRC,idas_incrc)

                  if( ierr .ne. 0 ) goto 999

                     jdsm = jdsm + 1
                     idas_incrc(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                     idas_incrc(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn
                     idasadress(ncsq) = jdsm ! frtati 2024/05/21

            else

               do j = 1, maxcntr-1   ! S.H. set maxcntr-1(reg case) (2022.3.24)

                  if( imsq(k) .eq. j+1 ) then

                     call snum(chlw,ic,i3,ic2,cvvv,ierr)

                     if( ierr .ne. 0 ) goto 999
                     if( cvvv .gt. 10000.d0 .or.
     &                   cvvv .lt. -9999.d0 ) goto 994

                     icnts(icont,j,nmreg) = nint(cvvv)

                     exit

                  endif

               enddo

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

               if( nmreg .eq. 0 ) goto 993

               ncreg(icont) = nmreg

               if( jdsm > MAX_NUM_INCRC ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.counts@read02.f ?dimension over idas_incrc?'
     &                    //' jdsm > MAX_NUM_INCRC'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INCRC@moddas.f=',MAX_NUM_INCRC,')'
                  ErrID = 'L:8385/R:counts/F:read02.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

cfrtati 2024/05/21
!               call moddas_reduce_int(
!     &                 3, icont, jdsm+1, incrc, idas_incrc)
               call moddas_reduce_int(
     &                  3, ncsq, idasadress(ncsq)+1, incrc, idas_incrc)
               idasadress(:) = incrc(:)
               do i = 1, ncsq
                 incrc(icontind(i)) = idasadress(i)
               end do

!               call moddas_reallocate_int(
!     &                 3, icont, nmreg*(maxcntr-1)+2, kcont, idas_kcont)  ! T.Sato 2024/03/21
                call moddas_reallocate_int(
     &                  3, ncsq, nmreg*(maxcntr-1)+2, kcont, idas_kcont)  ! frtati 2024/05/21
!               idsm  = kcont(icont)
               idsm  = kcont(ncsq) ! frtati 2024/05/21
               jdsm  = 0


               do i = 1, nmreg

                  do j = 1, maxcntr-1   ! S.H. set maxcntr-1(reg case) (2022.3.24)
                     jdsm = jdsm + 1
                     idas_kcont(idsm+jdsm) = icnts(icont,j,i)
                  enddo

               end do
cfrtati 2024/05/21
               idasadress(:) = kcont(:)
               do i = 1, ncsq
                 kcont(icontind(i)) = idasadress(i)
               end do

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:8433/R:counts/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  983 continue

         m_err = 'maximum number of part is 20'
         ErrCha = ''
         ErrID = 'L:8445/R:counts/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  984 continue

         m_err = 'description of part = is wrong'
         ErrCha = ''
         ErrID = 'L:8457/R:counts/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'there is no region in this counter definition'
         ErrCha = ''
         ErrID = 'L:8469/R:counts/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'counter unit should be from -9999 to 10000'
         ErrCha = ''
         ErrID = 'L:8481/R:counts/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         write(dkam,'(i6)') icont
         m_err = 'definition of counter ='//dkam//' is duplicated'
         ErrCha = ''
         ErrID = 'L:8494/R:counts/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'counter number is 1-3.'
         ErrCha = ''
         ErrID = 'L:8506/R:counts/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[counter] is wrong.'
         ErrCha = ''
         ErrID = 'L:8519/R:counts/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of region for counter '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:8533/R:counts/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [counter] is wrong.'
         ErrCha = ''
         ErrID = 'L:8545/R:counts/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine repcl(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [repeated collisions] section of input files              *
*       modified by K.Niita on 2019/10/20                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region
      use moddas_repeated_collisions

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /rclmsg/ ircln, isrcl, maxrr,
     &                mnrcl(6,0:20), mrrcl(kvlmax),
     &                krcls(6), lrcls(6), inrlc(6), inrlt(6),
     &                irpem(6,2), irman(6), irmct(6), jsmat(6)
      common /rclemm/ erpem(6,2)

*-----------------------------------------------------------------------

      logical dnen2

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 8 * 0 /

      character dkam*6

      dimension irclt(kvlmax)
      data irclt / kvlmax*1 /
      dimension jrclt(kvlmax)
      data jrclt / kvlmax*1 /

      dimension iptyp(20)
      dimension jstyp(6), jnkf0(6)

      common /subtra/ isubt, ipsub(mxpart)   ! kitamura22/03/31

      dimension     idas(1)
      equivalence ( das, idas )


      logical deqn1
      logical deqn4
      logical dcom2

*-----------------------------------------------------------------------

      character element(104)*3,cnuc*3

      data element/
     & 'h  ','he ','li ','be ','b  ','c  ','n  ','o  ',
     & 'f  ','ne ','na ','mg ','al ','si ','p  ','s  ',
     & 'cl ','ar ','k  ','ca ','sc ','ti ','v  ','cr ',
     & 'mn ','fe ','co ','ni ','cu ','zn ','ga ','ge ',
     & 'as ','se ','br ','kr ','rb ','sr ','y  ','zr ',
     & 'nb ','mo ','tc ','ru ','rh ','pd ','ag ','cd ',
     & 'in ','sn ','sb ','te ','i  ','xe ','cs ','ba ',
     & 'la ','ce ','pr ','nd ','pm ','sm ','eu ','gd ',
     & 'tb ','dy ','ho ','er ','tm ','yb ','lu ','hf ',
     & 'ta ','w  ','re ','os ','ir ','pt ','au ','hg ',
     & 'tl ','pb ','bi ','po ','at ','rn ','fr ','ra ',
     & 'ac ','th ','pa ','u  ','np ','pu ','am ','cm ',
     & 'bk ','cf ','es ','fm ','md ','no ','lr ','ku '/

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0
            inpat = 0
            intar = 0
            nnrcl = 0
            imoth = 0
            jmoth = 1

            ircln = ircln + 1

            if( ircln .gt. 6 ) goto 991

            irpem(ircln,1) = 0
            irpem(ircln,2) = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        particle, mother, emin, emax
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               ic  = i1
               icl = i1

               chlc = chlw
               call chcomp(chlc,icl,i3,i5)

*-----------------------------------------------------------------------

         if(chlc(icl:icl+3) .eq. 'part' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = i3

  400          continue

                  ic = jnumc(chlw,ic,i3)

                  if( ic .gt. i3 ) goto 140

                  call rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

                     if( ierr .eq. 994 ) goto 993
                     if( ierr .eq. 998 ) goto 993
                     if( isubt .eq. 1 )  goto 993   ! kitamura22/03/31

                  if( istyp .gt. 0 ) then

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = istyp

                  else if( istyp .lt. 0 ) then

                     do i = 1, -istyp

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = jstyp(i)

                     end do

                  end if

                  goto 400

*-----------------------------------------------------------------------
*        mother
*-----------------------------------------------------------------------

         else if(chlc(icl:icl+5) .eq. 'mother' ) then

               ic = inumc(chlw,icl+6,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 994

               icl = inumc(chlw,ic,i3,';') - 1

            if( chlw(ic:ic+2) .eq. 'all' ) then

               imoth = 0

            else if( deqn4( chlw(ic:ic) ) ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994

               imoth = nint( cvvv )

               if( imoth .lt. 0 ) then

                  imoth = -imoth
                  jmoth = -1

               end if

               if( imoth .eq. 0 ) goto 994

                  nsmat = jsmat(ircln)
                  call moddas_reallocate_int(
     &                    6, ircln, imoth, jsmat, ismat_jsmat)

  145                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( ierr .ne. 0 ) return
                        if( jpn  .eq. 3 ) goto 994

                        if( iskip .ne. 0 ) goto 145

                  ic = i1

               do k = 1, imoth

                  if( ic .gt. i3 ) then

  146                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) goto 994

                     if( iskip .ne. 0 ) goto 146

                     ic = i1

                  end if

                        ic = jnumc(chlw,ic,i3)

                        isa = 0
                        isn = 0

                        ica = 0
                        icb = 0
                        icm = 0
                        icn = 0

                     do i = ic, i3

                        if( chlw(i:i) .ge. 'a' .and.
     &                      chlw(i:i) .le. 'z' ) then

                           isa = isa + 1

                           if( isa .eq. 1 ) ica = i

                           icb = i

                        else if( deqn1( chlw(i:i) ) ) then

                           isn = isn + 1

                           if( isn .eq. 1 ) icm = i

                           icn = i

                        else if( dcom2( chlw(i:i) ) .or.
     &                           i .eq. i3 ) then

                           icd = i
                           goto 502

                        else

                           goto 994

                        end if

                     end do

                        icd = i3

  502                continue

                     if( ica .eq. 0 .or. icb .eq. 0 ) goto 994
                     if( icb-ica .lt. 0 .or. icb-ica .gt. 1 ) goto 994

                     cnuc = chlw(ica:icb)//'  '

                     do j = 1, 104

                        if( cnuc(1:3) .eq. element(j)(1:3) ) then

                           icha = j

                           goto 452

                        end if

                     end do

                           goto 994

  452                continue

                     if( icha .gt. 104 )  goto 994

                  if( icm .eq. 0 .or. icn .eq. 0 ) then

                     ismat_jsmat(nsmat-1+k) = icha * 1000

                  else

                     if( isn .gt. 3 ) goto 994

                     read(chlw(icm:icn),'(i5)') masi

                     if( masi .lt. icha ) goto 994
                     if( masi-icha .gt. maxnt ) goto 994

                     ismat_jsmat(nsmat-1+k) = icha * 1000 + masi

                  end if

                     ic = icd + 1

               end do

            else

               goto 994

            end if

               goto 140

*-----------------------------------------------------------------------

         else if(chlc(icl:icl+3) .eq. 'emin' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               if( ic .gt. i3 ) goto 995
               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 995

               irpem(ircln,1) = 1
               erpem(ircln,1) = cvvv

               goto 140

*-----------------------------------------------------------------------

         else if(chlc(icl:icl+3) .eq. 'emax' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               if( ic .gt. i3 ) goto 995
               icl = i3

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 995

               irpem(ircln,2) = 1
               erpem(ircln,2) = cvvv

               goto 140

*-----------------------------------------------------------------------

         end if

         end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+5) .eq. 'n-coll' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 6

            else if( chlw(ic:ic+5) .eq. 'n-evap' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 6

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inreg = 0
                  incol = 0
                  inevp = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) incol = incol + 1
                  if( imsq(k) .eq. 3 ) inevp = inevp + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inreg .ne. 1 .or. incol .ne. 1 ) goto 997
                  if( inevp .gt. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inreg = 1
                  incol = 1
                  inevp = 0
                  innon = 0

                  nrsq = 2

            end if

         end if

*-----------------------------------------------------------------------

            if( inpat .gt. 0 ) then

               do k = 1, inpat

                  if( iptyp(k) .eq. 20 ) inpat = 0

               end do

            end if

            if( inpat .eq. 0 ) then

               do k = 1, 19

                  iptyp(k) = k

               end do

                  inpat = 19

            end if

*-----------------------------------------------------------------------
*        read repeated collisions informations
*-----------------------------------------------------------------------

               nnrcl = nnrcl + 1

               if( nnrcl .gt. kvlmax ) goto 998

            if( nnrcl .eq. 1 ) then

                  call moddas_reallocate_int(
     &                    6, ircln, MAX_NUM_INRLC, inrlc, idas_inrlc)

                  idsm = inrlc(ircln)
                  jdsm = 0

            end if

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 2 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  irclt(nnrcl) = nint( cvvv )

            else if( imsq(k) .eq. 3 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  jrclt(nnrcl) = nint( cvvv )

            else if( imsq(k) .eq. 1 ) then

                     ndsm = idsm + jdsm + 3

                  call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ndsm,ierr
     &                          ,MAX_NUM_INRLC,idas_inrlc)

                  if( ierr .ne. 0 ) goto 999

                     jdsm = jdsm + 1
                     idas_inrlc(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                     idas_inrlc(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

            if( nnrcl .eq. 0 ) then

               ircln = ircln - 1
               return

            end if

               if( jdsm > MAX_NUM_INRLC ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.repcl@read02.f ?dimension over idas_inrlc?'
     &                    //' jdsm > MAX_NUM_INRLC'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INRLC@moddas.f=',MAX_NUM_INRLC,')'
                  ErrID = 'L:9159/R:repcl/F:read02.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 6, ircln, jdsm+1, inrlc, idas_inrlc)

               mnrcl(ircln, 0) = nnrcl
               mnrcl(ircln,20) = inpat

            do j = 1, inpat

               mnrcl(ircln,j) = iptyp(j)

            end do

*-----------------------------------------------------------------------

               call moddas_reallocate_int(
     &                 6, ircln, nnrcl, krcls, idas_krcls)
               idsm = krcls(ircln)


            do i = 1, nnrcl

               if( abs( irclt(i) ) * jrclt(i) .le. 1 ) goto 989
               idas_krcls(idsm-1+i) = irclt(i)

            end do

*-----------------------------------------------------------------------

         if( inevp .gt. 0 ) then

               call moddas_reallocate_int(
     &                 6, ircln, nnrcl, lrcls, idas_lrcls)
               idsm = lrcls(ircln)


            do i = 1, nnrcl

               if( jrclt(i) .le. 0 ) goto 988
               idas_lrcls(idsm-1+i) = jrclt(i)

            end do

         end if

*-----------------------------------------------------------------------

                  irman(ircln) = imoth
                  irmct(ircln) = jmoth

               if( imoth .gt. 0 ) then


               end if

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:9228/R:repcl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  988 continue

         m_err = 'n-evap of [repeated collisions] section is'//
     &           ' greater than zero'
         ErrCha = ''
         ErrID = 'L:9241/R:repcl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  989 continue

         m_err = 'factor of [repeated collisions] section is'//
     &           ' greater than one'
         ErrCha = ''
         ErrID = 'L:9254/R:repcl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'Number of [repeated collisions] section is'//
     &           ' too large (=<6)'
         ErrCha = ''
         ErrID = 'L:9267/R:repcl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'Number of particles is too large (<19)'
         ErrCha = ''
         ErrID = 'L:9279/R:repcl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'Name of particle is wrong'
         ErrCha = ''
         ErrID = 'L:9291/R:repcl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'Description of part or mother= is wrong'
         ErrCha = ''
         ErrID = 'L:9303/R:repcl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of emin or emax= is wrong'
         ErrCha = ''
         ErrID = 'L:9315/R:repcl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[repeated collisions] is wrong.'
         ErrCha = ''
         ErrID = 'L:9328/R:repcl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of repeated collisions '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:9342/R:repcl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [repeated collisions] is wrong.'
         ErrCha = ''
         ErrID = 'L:9354/R:repcl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine foccl(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [forced collisions] section of input files                *
*       modified by K.Niita on 2001/12/06                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'
      include 'param01.inc'   ! kitamura22/03/31

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      common /fclmsg/ ifcln, isfcl, maxrg,
     &                mnfcl(6,0:20), mrfcl(kvlmax),
     &                kfcls(6), inflc(6), inflt(6)

*-----------------------------------------------------------------------

      logical dnen2

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 8 * 0 /

      character dkam*6

      dimension vfclt(kvlmax)
      data vfclt / kvlmax*0.0d0 /

      dimension iptyp(20)
      dimension jstyp(6), jnkf0(6)

      common /subtra/ isubt, ipsub(mxpart)   ! kitamura22/03/31

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0
            inpat = 0
            nnfcl = 0

            ifcln = ifcln + 1

            if( ifcln .gt. 6 ) goto 991

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               ic  = i1
               icl = i1

               chlc = chlw
               call chcomp(chlc,icl,i3,i5)

            if(chlc(icl:icl+3) .eq. 'part' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = i3

  400          continue

                  ic = jnumc(chlw,ic,i3)

                  if( ic .gt. i3 ) goto 140

*-----------------------------------------------------------------------

                  call rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

cfrtati 2023/11/28 forced collision with dmax(15) and dmax(18) enabled
!                     if ( istyp.eq.15 .and. parz(145).gt.0d0 ) then
!                        goto 961
!                     else if ( istyp.eq.18 .and. parz(148).gt.0d0 ) then
!                        goto 962
!                     end if

                     if( ierr .eq. 994 ) goto 993
                     if( ierr .eq. 998 ) goto 993
                     if( isubt .eq. 1 )  goto 993   ! kitamura22/03/31

                  if( istyp .gt. 0 ) then

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = istyp

                  else if( istyp .lt. 0 ) then

                     do i = 1, -istyp

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = jstyp(i)

                     end do

                  end if

                  goto 400

            end if

         end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1

            else if( chlw(ic:ic+2) .eq. 'fcl' ) then

               imsq( mrsq + 1 ) = 2

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = ic + 3
               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inreg = 0
                  infcl = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) infcl = infcl + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inreg .ne. 1 .or. infcl .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inreg = 1
                  infcl = 1
                  innon = 0

                  nrsq = 2

            end if

         end if

*-----------------------------------------------------------------------

            if( inpat .gt. 0 ) then

               do k = 1, inpat

                  if( iptyp(k) .eq. 20 ) inpat = 0

               end do

            end if

            if( inpat .eq. 0 ) then

               do k = 1, 19

                  iptyp(k) = k

               end do

                  inpat = 19

            end if

*-----------------------------------------------------------------------
*        read forced collisions informations
*-----------------------------------------------------------------------

               nnfcl = nnfcl + 1

               if( nnfcl .gt. kvlmax ) goto 998

            if( nnfcl .eq. 1 ) then

                  call moddas_reallocate_int(
     &                    6, ifcln, MAX_NUM_INFLC, inflc, idas_inflc)

                  idsm = inflc(ifcln)
                  jdsm = 0

            end if

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 2 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  vfclt(nnfcl) = cvvv

            else if( imsq(k) .eq. 1 ) then

                     ndsm = idsm + jdsm + 3

                  call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ndsm,ierr
     &                          ,MAX_NUM_INFLC,idas_inflc)

                  if( ierr .ne. 0 ) goto 999

                     jdsm = jdsm + 1
                     idas_inflc(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                     idas_inflc(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

            if( nnfcl .eq. 0 ) then

               ifcln = ifcln - 1
               return

            end if

               if( jdsm > MAX_NUM_INRLC ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.foccl@read02.f ?dimension over idas_inflc?'
     &                    //' jdsm > MAX_NUM_INFLC'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INFLC@moddas.f=',MAX_NUM_INFLC,')'
                  ErrID = 'L:9721/R:foccl/F:read02.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 6, ifcln, jdsm+1, inflc, idas_inflc)

               mnfcl(ifcln, 0) = nnfcl
               mnfcl(ifcln,20) = inpat

            do j = 1, inpat

               mnfcl(ifcln,j) = iptyp(j)

            end do

               call moddas_reallocate_dbl(
     &                 6, ifcln, nnfcl+1, kfcls, das_kfcls)
               idsm = kfcls(ifcln)


            do i = 1, nnfcl

               if( abs( vfclt(i) ) .gt. 1.0d0 ) goto 990

               das_kfcls(idsm-1+i) = vfclt(i)

            end do

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

 961     continue

         m_err = '[forced collisions] cannot be used in tandem with '//
     &           'deuteron data library'
         ErrCha = ''
         ErrID = 'L:9761/R:foccl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

 962     continue

         m_err = '[forced collisions] cannot be used in tandem with '//
     &           'alpha data library'
         ErrCha = ''
         ErrID = 'L:9772/R:foccl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:9785/R:foccl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  990 continue

         m_err = 'factor of [forced collisions] section is'//
     &           ' greater than one'
         ErrCha = ''
         ErrID = 'L:9798/R:foccl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'Number of [forced collisions] section is'//
     &           ' too large (=<6)'
         ErrCha = ''
         ErrID = 'L:9811/R:foccl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'Number of particles is too large (<19)'
         ErrCha = ''
         ErrID = 'L:9823/R:foccl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'Name of particle is wrong'
         ErrCha = ''
         ErrID = 'L:9835/R:foccl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'Description of part = is wrong'
         ErrCha = ''
         ErrID = 'L:9847/R:foccl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[forced collisions] is wrong.'
         ErrCha = ''
         ErrID = 'L:9860/R:foccl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of forced collisions '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:9874/R:foccl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [forced collisions] is wrong.'
         ErrCha = ''
         ErrID = 'L:9886/R:foccl/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine impot(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [importance] section of input files                       *
*       modified by K.Niita on 2001/12/06                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'
      include 'param01.inc'   ! kitamura22/03/31

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /impmsg/ iimpn, isimp, iswct, isstr, maxip,
     &                mnimp(7,0:20),
     &                kfimp(7), inimc(7), inimt(7)

*-----------------------------------------------------------------------

      logical dnen2

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 8 * 0 /

      character dkam*6

      dimension vimpt(kvlmax)
      data vimpt / kvlmax*0.0d0 /

      dimension iptyp(20)
      dimension jstyp(6), jnkf0(6)

      common /subtra/ isubt, ipsub(mxpart)   ! kitamura22/03/31

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0
            inpat = 0
            nnimp = 0

            iimpn = iimpn + 1

            if( iimpn .gt. 6 ) goto 991

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               ic  = i1
               icl = i1

               chlc = chlw
               call chcomp(chlc,icl,i3,i5)

            if(chlc(icl:icl+3) .eq. 'part' ) then

               nrsq = 0

               ic = inumc(chlw,icl+4,i3,'=') + 1

               if( ic .gt. i3 ) goto 994

               icl = i3

  400          continue

                  ic = jnumc(chlw,ic,i3)

                  if( ic .gt. i3 ) goto 140

*-----------------------------------------------------------------------

                  call rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

                     if( ierr .eq. 994 ) goto 993
                     if( ierr .eq. 998 ) goto 993
                     if( isubt .eq. 1 )  goto 993   ! kitamura22/03/31

                  if( istyp .gt. 0 ) then

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = istyp

                  else if( istyp .lt. 0 ) then

                     do i = 1, -istyp

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = jstyp(i)

                     end do

                  end if

                  goto 400

            end if

         end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1

            else if( chlw(ic:ic+2) .eq. 'imp' ) then

               imsq( mrsq + 1 ) = 2

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = ic + 3
               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inreg = 0
                  inimp = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) inimp = inimp + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inreg .ne. 1 .or. inimp .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inreg = 1
                  inimp = 1
                  innon = 0

                  nrsq = 2

            end if

         end if

*-----------------------------------------------------------------------

            if( inpat .gt. 0 ) then

               do k = 1, inpat

                  if( iptyp(k) .eq. 20 ) inpat = 0

               end do

            end if

            if( inpat .eq. 0 ) then

               do k = 1, 19

                  iptyp(k) = k

               end do

                  inpat = 19

            end if

*-----------------------------------------------------------------------
*        read importance
*-----------------------------------------------------------------------

               nnimp = nnimp + 1

               if( nnimp .gt. kvlmax ) goto 998

            if( nnimp .eq. 1 ) then

                  call moddas_reallocate_int(
     &                    7, iimpn, MAX_NUM_INIMC, inimc, idas_inimc)

                  idsm = inimc(iimpn)
                  jdsm = 0

            end if

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 2 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  vimpt(nnimp) = cvvv

            else if( imsq(k) .eq. 1 ) then

                     ndsm = idsm + jdsm + 3

                  call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ndsm,ierr
     &                          ,MAX_NUM_INIMC,idas_inimc)

                  if( ierr .ne. 0 ) goto 999

                     jdsm = jdsm + 1
                     idas_inimc(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                     idas_inimc(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

            if( nnimp .eq. 0 ) then

               iimpn = iimpn - 1
               return

            end if

               if( jdsm > MAX_NUM_INIMC ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.impot@read02.f ?dimension over idas_inimc?'
     &                    //' jdsm > MAX_NUM_INIMC'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INIMC@moddas.f=',MAX_NUM_INIMC,')'
                  ErrID = 'L:10242/R:impot/F:read02.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 7, iimpn, jdsm+1, inimc, idas_inimc)

               mnimp(iimpn, 0) = nnimp
               mnimp(iimpn,20) = inpat

            do j = 1, inpat

               mnimp(iimpn,j) = iptyp(j)

            end do

               call moddas_reallocate_dbl(
     &                 7, iimpn, nnimp+1, kfimp, das_kfimp)
               idsm = kfimp(iimpn)


            do i = 1, nnimp

               if( vimpt(i) .lt. 0.0d0 ) goto 990

               das_kfimp(idsm-1+i) = vimpt(i)

            end do

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:10282/R:impot/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  990 continue

         m_err = 'factor of [importance] section is'//
     &           ' less than zero'
         ErrCha = ''
         ErrID = 'L:10295/R:impot/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'Number of [importance] section is'//
     &           ' too large (=<6)'
         ErrCha = ''
         ErrID = 'L:10308/R:impot/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'Number of particles is too large (<19)'
         ErrCha = ''
         ErrID = 'L:10320/R:impot/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'Name of particle is wrong'
         ErrCha = ''
         ErrID = 'L:10332/R:impot/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'Description of part = is wrong'
         ErrCha = ''
         ErrID = 'L:10344/R:impot/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[importance] is wrong.'
         ErrCha = ''
         ErrID = 'L:10357/R:impot/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of importance '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:10371/R:impot/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [importance] is wrong.'
         ErrCha = ''
         ErrID = 'L:10383/R:impot/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine delte(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [Delta Ray] section of input files                        *
*       modified by K.Niita on 2010/07/27                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /delesg/ rdels(kvlmax), mndel, ndels(kvlmax)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 8 * 0 /

      character dkam*6

      dimension ntrg(kvlmax)

      logical deqn4
      logical dnen1

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1

            else if( chlw(ic:ic+2) .eq. 'del' ) then

               imsq( mrsq + 1 ) = 2

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = ic + 3
               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inreg = 0
                  indel = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) indel = indel + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inreg .ne. 1 .or. indel .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inreg = 1
                  indel = 1
                  innon = 0

                  nrsq = 2

            end if

         end if

*-----------------------------------------------------------------------
*        read delt ray informations
*-----------------------------------------------------------------------

               mndel = mndel + 1

               if( mndel .gt. kvlmax ) goto 998

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 2 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               rdels(mndel) = cvvv

               if( cvvv .lt. 1.d-3 ) goto 994

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 1 ) then

                        ibra = 0
                        ilev = 0
                        ic = ic - 1

  180                ic = ic + 1

                     if( chlw(ic:ic) .eq. '(' ) then

                        if( ilev .gt. 0 ) goto 996
                        ibra = ibra + 1

                     else if( chlw(ic:ic) .eq. ')' ) then

                        ibra = ibra - 1
                        if( ibra .lt. 0 ) goto 996
                        if( ibra .eq. 0 ) goto 190

                     else if( chlw(ic:ic) .eq. '{' ) then

                        if( ilev .gt. 0 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '-' ) then

                        if( ilev .ne. 1 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '}' ) then

                        if( ilev .ne. 2 ) goto 995
                        ilev = 0

                        if( ntrf .le. ntri ) goto 995

                        do i = 0, ntrf - ntri

                           ntrn = ntrn + 1
                           ntrg(ntrn) = ntri + i

                        end do

                        if( ibra .eq. 0 ) goto 190

                     else if( deqn4( chlw(ic:ic) ) ) then

                           ici = ic

                        do j = ic + 1, i3

                           if( dnen1( chlw(j:j) ) ) goto 167

                        end do

  167                      icf = j - 1

                           if( icf .lt. ici ) goto 996

                        call onum(chlw,ici,icf,cvvv,ierr)

                           if( ierr .ne. 0 ) goto 996

                           ic = icf

*-----------------------------------------------------------------------

                        if( ilev .eq. 1 ) then

                           ntri = nint( cvvv )

                        else if( ilev .eq. 2 ) then

                           ntrf = nint( cvvv )

                        else

                           ntrn = ntrn + 1
                           ntrg(ntrn) = nint( cvvv )

                        end if

                        if( ilev .eq. 0 .and. ibra .eq. 0 ) goto 190

*-----------------------------------------------------------------------

                     end if

                           goto 180

  190                continue

                     ic2 = ic + 1

*-----------------------------------------------------------------------

            end if

         end do

               if( mndel + ntrn .gt. kvlmax ) goto 998

            do i = 1, ntrn

               mn = mndel + i - 1

               ndels(mn) = ntrg(i)
               rdels(mn) = rdels(mndel)

            end do

               mndel = mndel + ntrn - 1

         goto 140

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  994 continue

         m_err = 'Min. energy of delta should be '//
     &           'greater than 1 keV.'
         ErrCha = ''
         ErrID = 'L:10709/R:delte/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:10722/R:delte/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Description of region number '//
     &           '( {n1-n2} n3 n4 ) is wrong.'
         ErrCha = ''
         ErrID = 'L:10735/R:delte/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[Delta Ray] is wrong.'
         ErrCha = ''
         ErrID = 'L:10748/R:delte/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of Delta Ray '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:10762/R:delte/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [Delta Ray] is wrong.'
         ErrCha = ''
         ErrID = 'L:10774/R:delte/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine trackst(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [Track Structure] section of input files                  *
*       written by T.Sato based on subroutine delte                    *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /tscmsg/ ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /tscreg/ ntscell(kvlmax)
      common /etsart/ bgets(kvlmax), ebgets(kvlmax),
     &                wvets(kvlmax), ewvets(kvlmax)
*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 3, 4, 6 * 0 /

      character dkam*6

      dimension ntrg(kvlmax)

      logical deqn4
      logical dnen1

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1

            else if( chlw(ic:ic+2) .eq. 'mid' ) then

               imsq( mrsq + 1 ) = 2
! Hirata 20230210 ETSART Band gap inputs
            else if( chlw(ic:ic+2) .eq. 'ebg' ) then

               imsq( mrsq + 1 ) = 3

            else if( chlw(ic:ic+2) .eq. 'wvl' ) then

               imsq( mrsq + 1 ) = 4

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = ic + 3
               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inreg = 0
                  indel = 0
                  innon = 0
                  inband= 0
                  inwvl= 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg  = inreg + 1
                  if( imsq(k) .eq. 2 ) indel  = indel + 1
                  if( imsq(k) .eq. 3 ) inband = inband + 1
                  if( imsq(k) .eq. 4 ) inwvl  = inwvl + 1
                  if( imsq(k) .eq. 0 ) innon  = innon + 1
               end do

                  if( inreg .ne. 1 .or. indel .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inreg = 1
                  indel = 1
                  inband= 1
                  inwvl = 1
                  innon = 0

                  nrsq = 4

            end if

         end if

*-----------------------------------------------------------------------
*        read track structure informations
*-----------------------------------------------------------------------

               mntsc = mntsc + 1

               if( mntsc .gt. kvlmax ) goto 998

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 2 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               ktsc(mntsc) = int(cvvv)

*-----------------------------------------------------------------------
! Hirata 20240111 etsart w value (w-value)
            else if( imsq(k) .eq. 4 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) then
                wvets(mntsc) = 0.d0
               else
                wvets(mntsc) = cvvv
               endif


! Hirata 20230210 etsart BG
            else if( imsq(k) .eq. 3 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) then
                bgets(mntsc) = 0.d0
               else
                bgets(mntsc) = cvvv
               endif

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 1 ) then

                        ibra = 0
                        ilev = 0
                        ic = ic - 1

  180                ic = ic + 1

                     if( chlw(ic:ic) .eq. '(' ) then

                        if( ilev .gt. 0 ) goto 996
                        ibra = ibra + 1

                     else if( chlw(ic:ic) .eq. ')' ) then

                        ibra = ibra - 1
                        if( ibra .lt. 0 ) goto 996
                        if( ibra .eq. 0 ) goto 190

                     else if( chlw(ic:ic) .eq. '{' ) then

                        if( ilev .gt. 0 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '-' ) then

                        if( ilev .ne. 1 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '}' ) then

                        if( ilev .ne. 2 ) goto 995
                        ilev = 0

                        if( ntrf .le. ntri ) goto 995

                        do i = 0, ntrf - ntri

                           ntrn = ntrn + 1
                           ntrg(ntrn) = ntri + i

                        end do

                        if( ibra .eq. 0 ) goto 190

                     else if( deqn4( chlw(ic:ic) ) ) then

                           ici = ic

                        do j = ic + 1, i3

                           if( dnen1( chlw(j:j) ) ) goto 167

                        end do

  167                      icf = j - 1

                           if( icf .lt. ici ) goto 996

                        call onum(chlw,ici,icf,cvvv,ierr)

                           if( ierr .ne. 0 ) goto 996

                           ic = icf

*-----------------------------------------------------------------------

                        if( ilev .eq. 1 ) then

                           ntri = nint( cvvv )

                        else if( ilev .eq. 2 ) then

                           ntrf = nint( cvvv )

                        else

                           ntrn = ntrn + 1
                           ntrg(ntrn) = nint( cvvv )

                        end if

                        if( ilev .eq. 0 .and. ibra .eq. 0 ) goto 190

*-----------------------------------------------------------------------

                     end if

                           goto 180

  190                continue

                     ic2 = ic + 1

*-----------------------------------------------------------------------

            end if

         end do

               if( mntsc + ntrn .gt. kvlmax ) goto 998

            do i = 1, ntrn

               mn = mntsc + i - 1

               ntsc(mn) = ntrg(i)
               ktsc(mn) = ktsc(mntsc)

            end do

               mntsc = mntsc + ntrn - 1

         goto 140

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:11137/R:trackst/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Description of region number '//
     &           '( {n1-n2} n3 n4 ) is wrong.'
         ErrCha = ''
         ErrID = 'L:11150/R:trackst/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[Track Structure] is wrong.'
         ErrCha = ''
         ErrID = 'L:11163/R:trackst/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of Track Structure '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:11177/R:trackst/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [Track Structure] is wrong.'
         ErrCha = ''
         ErrID = 'L:11189/R:trackst/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine volum(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [volume] section of input files                           *
*       modified by K.Niita on 12/06/2001                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /volmsg/ rvols(kvlmax), mnvol, nvols(kvlmax)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 8 * 0 /

      character dkam*6

      dimension ntrg(kvlmax)

      logical deqn4
      logical dnen1

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1

            else if( chlw(ic:ic+2) .eq. 'vol' ) then

               imsq( mrsq + 1 ) = 2

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = ic + 3
               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inreg = 0
                  invol = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) invol = invol + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inreg .ne. 1 .or. invol .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inreg = 1
                  invol = 1
                  innon = 0

                  nrsq = 2

            end if

         end if

*-----------------------------------------------------------------------
*        read volume informations
*-----------------------------------------------------------------------

               mnvol = mnvol + 1

               if( mnvol .gt. kvlmax ) goto 998

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 2 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               rvols(mnvol) = cvvv

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 1 ) then

                        ibra = 0
                        ilev = 0
                        ic = ic - 1

  180                ic = ic + 1

                     if( chlw(ic:ic) .eq. '(' ) then

                        if( ilev .gt. 0 ) goto 996
                        ibra = ibra + 1

                     else if( chlw(ic:ic) .eq. ')' ) then

                        ibra = ibra - 1
                        if( ibra .lt. 0 ) goto 996
                        if( ibra .eq. 0 ) goto 190

                     else if( chlw(ic:ic) .eq. '{' ) then

                        if( ilev .gt. 0 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '-' ) then

                        if( ilev .ne. 1 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '}' ) then

                        if( ilev .ne. 2 ) goto 995
                        ilev = 0

                        if( ntrf .le. ntri ) goto 995

                        do i = 0, ntrf - ntri

                           ntrn = ntrn + 1
                           ntrg(ntrn) = ntri + i

                        end do

                        if( ibra .eq. 0 ) goto 190

                     else if( deqn4( chlw(ic:ic) ) ) then

                           ici = ic

                        do j = ic + 1, i3

                           if( dnen1( chlw(j:j) ) ) goto 167

                        end do

  167                      icf = j - 1

                           if( icf .lt. ici ) goto 996

                        call onum(chlw,ici,icf,cvvv,ierr)

                           if( ierr .ne. 0 ) goto 996

                           ic = icf

*-----------------------------------------------------------------------

                        if( ilev .eq. 1 ) then

                           ntri = nint( cvvv )

                        else if( ilev .eq. 2 ) then

                           ntrf = nint( cvvv )

                        else

                           ntrn = ntrn + 1
                           ntrg(ntrn) = nint( cvvv )

                        end if

                        if( ilev .eq. 0 .and. ibra .eq. 0 ) goto 190

*-----------------------------------------------------------------------

                     end if

                           goto 180

  190                continue

                     ic2 = ic + 1

*-----------------------------------------------------------------------

            end if

         end do

               if( mnvol + ntrn .gt. kvlmax ) goto 998

            do i = 1, ntrn

               mn = mnvol + i - 1

               nvols(mn) = ntrg(i)
               rvols(mn) = rvols(mnvol)

            end do

               mnvol = mnvol + ntrn - 1

         goto 140

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:11513/R:volum/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Description of region number '//
     &           '( {n1-n2} n3 n4 ) is wrong.'
         ErrCha = ''
         ErrID = 'L:11526/R:volum/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[volume] is wrong.'
         ErrCha = ''
         ErrID = 'L:11539/R:volum/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of volume '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:11553/R:volum/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [volume] is wrong.'
         ErrCha = ''
         ErrID = 'L:11565/R:volum/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine tempe(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [temperature] section of input files                      *
*       modified by K.Niita on 2002/03/11                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /tmpmsg/ rtmps(kvlmax), mntmp, ntmps(kvlmax)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 8 * 0 /

      character dkam*6

      dimension ntrg(kvlmax)

      logical deqn4
      logical dnen1

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1

            else if( chlw(ic:ic+2) .eq. 'tmp' ) then

               imsq( mrsq + 1 ) = 2

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = ic + 3
               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inreg = 0
                  intmp = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) intmp = intmp + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inreg .ne. 1 .or. intmp .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inreg = 1
                  intmp = 1
                  innon = 0

                  nrsq = 2

            end if

         end if

*-----------------------------------------------------------------------
*        read temperature informations
*-----------------------------------------------------------------------

               mntmp = mntmp + 1

               if( mntmp .gt. kvlmax ) goto 998

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 2 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               rtmps(mntmp) = cvvv

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 1 ) then

                        ibra = 0
                        ilev = 0
                        ic = ic - 1

  180                ic = ic + 1

                     if( chlw(ic:ic) .eq. '(' ) then

                        if( ilev .gt. 0 ) goto 996
                        ibra = ibra + 1

                     else if( chlw(ic:ic) .eq. ')' ) then

                        ibra = ibra - 1
                        if( ibra .lt. 0 ) goto 996
                        if( ibra .eq. 0 ) goto 190

                     else if( chlw(ic:ic) .eq. '{' ) then

                        if( ilev .gt. 0 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '-' ) then

                        if( ilev .ne. 1 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '}' ) then

                        if( ilev .ne. 2 ) goto 995
                        ilev = 0

                        if( ntrf .le. ntri ) goto 995

                        do i = 0, ntrf - ntri

                           ntrn = ntrn + 1
                           ntrg(ntrn) = ntri + i

                        end do

                        if( ibra .eq. 0 ) goto 190

                     else if( deqn4( chlw(ic:ic) ) ) then

                           ici = ic

                        do j = ic + 1, i3

                           if( dnen1( chlw(j:j) ) ) goto 167

                        end do

  167                      icf = j - 1

                           if( icf .lt. ici ) goto 996

                        call onum(chlw,ici,icf,cvvv,ierr)

                           if( ierr .ne. 0 ) goto 996

                           ic = icf

*-----------------------------------------------------------------------

                        if( ilev .eq. 1 ) then

                           ntri = nint( cvvv )

                        else if( ilev .eq. 2 ) then

                           ntrf = nint( cvvv )

                        else

                           ntrn = ntrn + 1
                           ntrg(ntrn) = nint( cvvv )

                        end if

                        if( ilev .eq. 0 .and. ibra .eq. 0 ) goto 190

*-----------------------------------------------------------------------

                     end if

                           goto 180

  190                continue

                     ic2 = ic + 1

*-----------------------------------------------------------------------

            end if

         end do

               if( mntmp + ntrn .gt. kvlmax ) goto 998

            do i = 1, ntrn

               mn = mntmp + i - 1

               ntmps(mn) = ntrg(i)
               rtmps(mn) = rtmps(mntmp)

            end do

               mntmp = mntmp + ntrn - 1

         goto 140

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:11889/R:tempe/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Description of region number '//
     &           '( {n1-n2} n3 n4 ) is wrong.'
         ErrCha = ''
         ErrID = 'L:11902/R:tempe/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[temperature] is wrong.'
         ErrCha = ''
         ErrID = 'L:11915/R:tempe/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of temperature '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:11929/R:tempe/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [temperature] is wrong.'
         ErrCha = ''
         ErrID = 'L:11941/R:tempe/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine bremb(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [brems bias] section of input files                       *
*       modified by K.Niita on 2002/03/13                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /brsmsg/ mnbrs, icbrs, mbbrs(kvlmax), cbrem(49)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)
      dimension ill(0:9), ilf(0:9)
      dimension imsq(10)
      data imsq / 1, 2, 8 * 0 /
      character dkam*6

      dimension ntrg(49)

      logical deqn4
      logical dnen1

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

         do i = 1, 49

            cbrem(i) = 1.0

         end do

            icbrs = 1
            imate = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        material definition
*-----------------------------------------------------------------------

         if( chcm(i1:i1+7) .eq. 'material' ) then

               imate = imate + 1
               if( imate .gt. 1 ) goto 993

               ic = inumc(chlw,i1+8,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 993

            if( chlw(ic:ic+2) .eq. 'all' ) then

                  mnbrs = 0
                  icbrs = 0

            else if( deqn4( chlw(ic:ic) ) ) then

               call onum(chlw,ic,i3,cvvv,ierr)

               if( ierr .ne. 0 ) goto 993

               mnbrs = nint( cvvv )

               if( mnbrs .lt. 0 ) then

                  mnbrs = -mnbrs
                  icbrs = -1

               end if

               if( mnbrs .eq. 0 ) goto 993

  151                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( ierr .ne. 0 ) return
                        if( jpn  .eq. 3 ) goto 993

                        if( iskip .ne. 0 ) goto 151

                  ic = i1

               do k = 1, mnbrs

                  if( ic .gt. i3 ) then

  152                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) goto 993

                     if( iskip .ne. 0 ) goto 152

                     ic = i1

                  end if

                     ic = jnumc(chlw,ic,i3)

                     call snum(chlw,ic,i3,ic2,cvvv,ierr)

                     if( ierr .ne. 0 ) goto 993

                     matei = nint( cvvv )

                     ic = ic2

                     mbbrs(k) = matei

               end do

            else

               goto 993

            end if

               goto 140

         end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0

            else if( chlw(ic:ic+2) .eq. 'num' ) then

               imsq( mrsq + 1 ) = 1

            else if( chlw(ic:ic+3) .eq. 'bias' ) then

               imsq( mrsq + 1 ) = 2

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = ic + 3
               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  innum = 0
                  inbrs = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) innum = innum + 1
                  if( imsq(k) .eq. 2 ) inbrs = inbrs + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( innum .ne. 1 .or. inbrs .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  innum = 1
                  inbrs = 1
                  innon = 0

                  nrsq = 2

            end if

         end if

*-----------------------------------------------------------------------
*        read brems bias informations
*-----------------------------------------------------------------------

               ic2  = i1
               ntrn = 0
               bbrms = 0.0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 2 ) then

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  bbrms = cvvv

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 1 ) then

                        ibra = 0
                        ilev = 0
                        ic = ic - 1

  180                ic = ic + 1

                     if( chlw(ic:ic) .eq. '(' ) then

                        if( ilev .gt. 0 ) goto 996
                        ibra = ibra + 1

                     else if( chlw(ic:ic) .eq. ')' ) then

                        ibra = ibra - 1
                        if( ibra .lt. 0 ) goto 996
                        if( ibra .eq. 0 ) goto 190

                     else if( chlw(ic:ic) .eq. '{' ) then

                        if( ilev .gt. 0 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '-' ) then

                        if( ilev .ne. 1 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '}' ) then

                        if( ilev .ne. 2 ) goto 995
                        ilev = 0

                        if( ntrf .le. ntri ) goto 995
                        if( ntri .le. 0 .or.
     &                      ntrf .gt. 49 ) goto 998

                        do i = 0, ntrf - ntri

                           ntrn = ntrn + 1
                           if( ntrn .gt. 49 ) goto 998
                           ntrg(ntrn) = ntri + i

                        end do

                        if( ibra .eq. 0 ) goto 190

                     else if( deqn4( chlw(ic:ic) ) ) then

                           ici = ic

                        do j = ic + 1, i3

                           if( dnen1( chlw(j:j) ) ) goto 167

                        end do

  167                      icf = j - 1

                           if( icf .lt. ici ) goto 996

                        call onum(chlw,ici,icf,cvvv,ierr)

                           if( ierr .ne. 0 ) goto 996

                           ic = icf

*-----------------------------------------------------------------------

                        if( ilev .eq. 1 ) then

                           ntri = nint( cvvv )

                        else if( ilev .eq. 2 ) then

                           ntrf = nint( cvvv )

                        else

                           ntrn = ntrn + 1
                           if( ntrn .gt. 49 ) goto 998
                           ntrg(ntrn) = nint( cvvv )

                           if( ntrg(ntrn) .le. 0 .or.
     &                         ntrg(ntrn) .gt. 49 ) goto 998

                        end if

                        if( ilev .eq. 0 .and. ibra .eq. 0 ) goto 190

*-----------------------------------------------------------------------

                     end if

                           goto 180

  190                continue

                     ic2 = ic + 1

*-----------------------------------------------------------------------

            end if

         end do

               if( ntrn .eq. 0 ) goto 999
               if( bbrms .le. 0.0 ) goto 994

            do i = 1, ntrn

               cbrem( ntrg(i) ) = bbrms

            end do

         goto 140

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  993 continue

         m_err = 'material description is wrong.'
         ErrCha = ''
         ErrID = 'L:12355/R:bremb/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'bias factor should be positive.'
         ErrCha = ''
         ErrID = 'L:12367/R:bremb/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:12380/R:bremb/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Description of number '//
     &           '( {n1-n2} n3 n4 ) is wrong.'
         ErrCha = ''
         ErrID = 'L:12393/R:bremb/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[brems bias] is wrong.'
         ErrCha = ''
         ErrID = 'L:12406/R:bremb/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') 49
         m_err = 'Number of bias factors '//
     &           'should be witin 1 and 49 = '// dkam
         ErrCha = ''
         ErrID = 'L:12420/R:bremb/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [brems bias] is wrong.'
         ErrCha = ''
         ErrID = 'L:12432/R:bremb/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine photw(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [photon weight] section of input files                    *
*       modified by K.Niita on 2002/03/19                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /pwtmsg/ rpwts(kvlmax), mnpwt, npwts(kvlmax)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 8 * 0 /

      character dkam*6

      dimension ntrg(kvlmax)

      logical deqn4
      logical dnen1

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1

            else if( chlw(ic:ic+2) .eq. 'pwt' ) then

               imsq( mrsq + 1 ) = 2

            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = ic + 3
               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 2 ) goto 997

                  inreg = 0
                  inpwt = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) inpwt = inpwt + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inreg .ne. 1 .or. inpwt .ne. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inreg = 1
                  inpwt = 1
                  innon = 0

                  nrsq = 2

            end if

         end if

*-----------------------------------------------------------------------
*        read temperature informations
*-----------------------------------------------------------------------

               mnpwt = mnpwt + 1

               if( mnpwt .gt. kvlmax ) goto 998

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 2 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               rpwts(mnpwt) = cvvv

*-----------------------------------------------------------------------

            else if( imsq(k) .eq. 1 ) then

                        ibra = 0
                        ilev = 0
                        ic = ic - 1

  180                ic = ic + 1

                     if( chlw(ic:ic) .eq. '(' ) then

                        if( ilev .gt. 0 ) goto 996
                        ibra = ibra + 1

                     else if( chlw(ic:ic) .eq. ')' ) then

                        ibra = ibra - 1
                        if( ibra .lt. 0 ) goto 996
                        if( ibra .eq. 0 ) goto 190

                     else if( chlw(ic:ic) .eq. '{' ) then

                        if( ilev .gt. 0 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '-' ) then

                        if( ilev .ne. 1 ) goto 995
                        ilev = ilev + 1

                     else if( chlw(ic:ic) .eq. '}' ) then

                        if( ilev .ne. 2 ) goto 995
                        ilev = 0

                        if( ntrf .le. ntri ) goto 995

                        do i = 0, ntrf - ntri

                           ntrn = ntrn + 1
                           ntrg(ntrn) = ntri + i

                        end do

                        if( ibra .eq. 0 ) goto 190

                     else if( deqn4( chlw(ic:ic) ) ) then

                           ici = ic

                        do j = ic + 1, i3

                           if( dnen1( chlw(j:j) ) ) goto 167

                        end do

  167                      icf = j - 1

                           if( icf .lt. ici ) goto 996

                        call onum(chlw,ici,icf,cvvv,ierr)

                           if( ierr .ne. 0 ) goto 996

                           ic = icf

*-----------------------------------------------------------------------

                        if( ilev .eq. 1 ) then

                           ntri = nint( cvvv )

                        else if( ilev .eq. 2 ) then

                           ntrf = nint( cvvv )

                        else

                           ntrn = ntrn + 1
                           ntrg(ntrn) = nint( cvvv )

                        end if

                        if( ilev .eq. 0 .and. ibra .eq. 0 ) goto 190

*-----------------------------------------------------------------------

                     end if

                           goto 180

  190                continue

                     ic2 = ic + 1

*-----------------------------------------------------------------------

            end if

         end do

               if( mnpwt + ntrn .gt. kvlmax ) goto 998

            do i = 1, ntrn

               mn = mnpwt + i - 1

               npwts(mn) = ntrg(i)
               rpwts(mn) = rpwts(mnpwt)

            end do

               mnpwt = mnpwt + ntrn - 1

         goto 140

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:12756/R:photw/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Description of region number '//
     &           '( {n1-n2} n3 n4 ) is wrong.'
         ErrCha = ''
         ErrID = 'L:12769/R:photw/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[photon weight] is wrong.'
         ErrCha = ''
         ErrID = 'L:12782/R:photw/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of photon weight '//
     &           'exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:12796/R:photw/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [photon weight] is wrong.'
         ErrCha = ''
         ErrID = 'L:12808/R:photw/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine mgnet(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [magnetic field] section of input files                   *
*       modified by K.Niita on 2003/09/03                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /magreg/ nmreg, ingrc, ingrt, kmags, magtin, imgusr

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 3, 4, 6 * 0 /

      character dkam*6

      dimension a_mag(kvlmax), b_mag(kvlmax), s_mag(kvlmax)
      dimension t_mag(kvlmax), p_mag(kvlmax), u_mag(kvlmax)

      character map1(kvlmax)*200 ! T.Sato 2019/01/14
      dimension icmap1(kvlmax)   ! T.Sato 2019/01/19

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 .and. nmreg .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'typ' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'gap' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'mgf' ) then

               imsq( mrsq + 1 ) = 4
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'trcl' ) then

               imsq( mrsq + 1 ) = 5
               ic = ic + 4

            else if( chlw(ic:ic+10) .eq. 'phase/polar' .or.
     &               chlw(ic:ic+10) .eq. 'polar/phase' ) then

               imsq( mrsq + 1 ) = 6
               ic = ic + 11

            else if( chlw(ic:ic+4) .eq. 'phase' .or.
     &               chlw(ic:ic+4) .eq. 'polar' ) then

               imsq( mrsq + 1 ) = 6
               ic = ic + 5

            else if( chlw(ic:ic+3) .eq. 'time' ) then

               imsq( mrsq + 1 ) = 7
               ic = ic + 4

*           ASTOM 2018/11/22
            else if( chlw(ic:ic+3) .eq. 'file' ) then

               imsq( mrsq + 1 ) = 8
               ic = ic + 4

            else

               goto 997

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 4 ) goto 997

                  inreg = 0
                  intyp = 0
                  ingap = 0
                  inmgf = 0
                  innon = 0
                  itrcl = 0
                  iphas = 0
                  itime = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) intyp = intyp + 1
                  if( imsq(k) .eq. 3 ) ingap = ingap + 1
                  if( imsq(k) .eq. 4 ) inmgf = inmgf + 1
                  if( imsq(k) .eq. 5 ) itrcl = itrcl + 1
                  if( imsq(k) .eq. 6 ) iphas = iphas + 1
                  if( imsq(k) .eq. 7 ) itime = itime + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inreg .ne. 1 .or. intyp .ne. 1 .or.
     &                ingap .ne. 1 .or. inmgf .ne. 1 .or.
     &                itrcl .gt. 1 .or. iphas .gt. 1 .or.
     &                itime .gt. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inreg = 1
                  intyp = 1
                  ingap = 1
                  inmgf = 1
                  itrcl = 0
                  iphas = 0
                  itime = 0
                  innon = 0

                  nrsq = 4

            end if

         end if

*-----------------------------------------------------------------------
*        read magnetic field informations
*-----------------------------------------------------------------------

               nmreg = nmreg + 1

               if( nmreg .gt. kvlmax ) goto 998

                  a_mag(nmreg)  = 0.0d0
                  b_mag(nmreg)  = 0.0d0
                  s_mag(nmreg)  = 0.0d0
                  t_mag(nmreg)  = 0.0d0
                  u_mag(nmreg)  = -1.0d+9
                  p_mag(nmreg)  = -10000.0d0

            if( nmreg .eq. 1 ) then

                  ingrc = 0
                  iaddress_region(:) = 0
                  call moddas_reallocate_int(
     &                    2, 1, MAX_NUM_INGRC, iaddress_region
     &                    , idas_ingrc)
                  idsm  = ingrc
                  jdsm  = 0

            end if

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 1 ) then

                     ndsm = idsm + jdsm + 3

                  call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ndsm,ierr
     &                          ,MAX_NUM_INGRC,idas_ingrc)

                  if( ierr .ne. 0 ) goto 999

                     jdsm = jdsm + 1
                     idas_ingrc(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                    idas_ingrc(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn

            else if( imsq(k) .eq. 2 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               s_mag(nmreg) = cvvv

            else if( imsq(k) .eq. 3 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               a_mag(nmreg) = cvvv

            else if( imsq(k) .eq. 4 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               b_mag(nmreg) = cvvv

            else if( imsq(k) .eq. 5 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               t_mag(nmreg) = cvvv

            else if( imsq(k) .eq. 7 ) then

               if( chlw(ic:ic+2) .eq. 'non' ) then

                  u_mag(nmreg) = -1.0d+9
                  ic2 = ic + 3

               else

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  u_mag(nmreg) = -1.0d+9
                  if( cvvv  .gt. -1.0d+9 )
     &            u_mag(nmreg) = cvvv

               end if

            else if( imsq(k) .eq. 6 ) then

               if( chlw(ic:ic+2) .eq. 'non' ) then

                  p_mag(nmreg) = -10000.0d0
                  ic2 = ic + 3

               else

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  p_mag(nmreg) = -10000.0d0
                  if( cvvv  .gt. -10000.0d0 )
     &            p_mag(nmreg) = cvvv

               end if

*           T.Sato 2019/01/18
            else if( imsq(k).eq.8 ) then
               call mapnum(chin,200,ic,ic1,ic2,ierr)
               if( ierr .ne. 0 ) goto 999
               map1(nmreg) = chin(ic1:ic2)
               icmap1(nmreg) = ic2-ic1+1
            end if
         end do

         goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

               if( jdsm > MAX_NUM_INGRC ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.mgnet@read02.f ?dimension over idas_ingrc?'
     &                    //' jdsm > MAX_NUM_INGRC'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INGRC@moddas.f=',MAX_NUM_INGRC,')'
                  ErrID = 'L:13191/R:mgnet/F:read02.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 2, 1, jdsm+1, iaddress_region, idas_ingrc)

               kmags = 0
               call moddas_allocate_dbl(nmreg*7+1, das_kmags)
               idsm  = kmags
               jdsm  = 0


            ixyzlist=0 ! number of xyz list
            irzlist=0  ! number of r-z list
            ixyzmap=0  ! number of xyz map
            irzmap=0   ! number of r-z map

            do i = 1, nmreg

               jdsm = jdsm + 1
               das_kmags(idsm+jdsm) = a_mag(i)
               jdsm = jdsm + 1
               das_kmags(idsm+jdsm) = b_mag(i)
               jdsm = jdsm + 1
               das_kmags(idsm+jdsm) = s_mag(i)
               jdsm = jdsm + 1
               das_kmags(idsm+jdsm) = p_mag(i)
               jdsm = jdsm + 1
               das_kmags(idsm+jdsm) = 0
               jdsm = jdsm + 1
               das_kmags(idsm+jdsm) = t_mag(i)
               jdsm = jdsm + 1
               das_kmags(idsm+jdsm) = u_mag(i)

               if(s_mag(i).eq.-1.or.s_mag(i).eq.-101) then
                if(ixyzlist.eq.0) then
                 call readmagxyzlist(map1(i),icmap1(i))
                 ixyzlist=1
                else
                 goto 991 ! only 1 xyzlist is allowed
                endif
               elseif(s_mag(i).eq.-2.or.s_mag(i).eq.-102) then
                if(irzlist.eq.0) then
                 call readmagrzlist(map1(i),icmap1(i))
                 irzlist=1
                else
                 goto 992 ! only 1 rzlist is allowed
                endif
               elseif(s_mag(i).eq.-3.or.s_mag(i).eq.-103) then
                if(ixyzmap.eq.0) then
                 call readmagxyzmap(map1(i),icmap1(i))
                 ixyzmap=1
                else
                 goto 993 ! only 1 rxyzmap is allowed
                endif
               elseif(s_mag(i).eq.-4.or.s_mag(i).eq.-104) then
                if(irzmap.eq.0) then
                 call readmagrzmap(map1(i),icmap1(i))
                 irzmap=1
                else
                 goto 994 ! only 1 rzmap is allowed
                endif
               endif
            end do

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:13268/R:mgnet/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'Only 1 xyz magnetic field list is allowed in '//
     &           '[magnetic field]'
         ErrCha = ''
         ErrID = 'L:13281/R:mgnet/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  992 continue

         m_err = 'Only 1 r-z magnetic field list is allowed in '//
     &           '[magnetic field]'
         ErrCha = ''
         ErrID = 'L:13292/R:mgnet/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  993 continue

         m_err = 'Only 1 xyz magnetic field map is allowed in '//
     &           '[magnetic field]'
         ErrCha = ''
         ErrID = 'L:13303/R:mgnet/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  994 continue

         m_err = 'Only 1 r-z magnetic field map is allowed in '//
     &           '[magnetic field]'
         ErrCha = ''
         ErrID = 'L:13314/R:mgnet/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[magnetic field] is wrong.'
         ErrCha = ''
         ErrID = 'L:13327/R:mgnet/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of region for magnetic '//
     &           'field exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:13341/R:mgnet/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [magnetic field] is wrong.'
         ErrCha = ''
         ErrID = 'L:13353/R:mgnet/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine elmgf(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [electro magnetic field] section of input files           *
*       modified by K.Niita on 2011/01/11                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /elmgrg/ nelctf, nereg, inerc, inert, kelcs

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension imsq(10)

      data imsq / 1, 2, 3, 7 * 0 /

      character dkam*6

      dimension s_elf(kvlmax), s_mgf(kvlmax),
     &          t_elf(kvlmax), t_mgf(kvlmax)

      dimension emap_type(kvlmax), mmap_type(kvlmax)
      dimension a_elmg(kvlmax)
      character map1(kvlmax)*200
      dimension icmap1(kvlmax)
      character map2(kvlmax)*200
      dimension icmap2(kvlmax)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               imsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               imsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'elf' ) then

               imsq( mrsq + 1 ) = 2
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'mgf' ) then

               imsq( mrsq + 1 ) = 3
               ic = ic + 3

            else if( chlw(ic:ic+4) .eq. 'trcle' ) then

               imsq( mrsq + 1 ) = 4
               ic = ic + 5

            else if( chlw(ic:ic+4) .eq. 'trclm' ) then

               imsq( mrsq + 1 ) = 5
               ic = ic + 5

            else if( chlw(ic:ic+3) .eq. 'type') then

               imsq( mrsq + 1 ) = 6
               ic = ic + 4

            else if( chlw(ic:ic+3) .eq. 'typm') then

               imsq( mrsq + 1 ) = 7
               ic = ic + 4

            else if( chlw(ic:ic+2) .eq. 'gap') then

               imsq( mrsq + 1 ) = 8
               ic = ic + 3

            else if( chlw(ic:ic+4) .eq. 'filee') then

               imsq( mrsq + 1 ) = 9
               ic = ic + 5

            else if( chlw(ic:ic+4) .eq. 'filem') then

               imsq( mrsq + 1 ) = 10
               ic = ic + 5

            else

               goto 997

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .lt. 3 ) goto 997

                  inreg = 0
                  inelf = 0
                  inmgf = 0
                  itrce = 0
                  itrcm = 0
                  innon = 0

               do k = 1, mrsq
                  if( imsq(k) .eq. 1 ) inreg = inreg + 1
                  if( imsq(k) .eq. 2 ) inelf = inelf + 1
                  if( imsq(k) .eq. 3 ) inmgf = inmgf + 1
                  if( imsq(k) .eq. 4 ) itrce = itrce + 1
                  if( imsq(k) .eq. 5 ) itrcm = itrcm + 1
                  if( imsq(k) .eq. 0 ) innon = innon + 1
               end do

                  if( inreg .ne. 1 .or. inelf .ne. 1 .or.
     &                inmgf .ne. 1 .or. itrce .gt. 1 .or.
     &                itrcm .gt. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inreg = 1
                  inelf = 1
                  inmgf = 1
                  itrce = 0
                  itrcm = 0
                  innon = 0

                  nrsq = 3

            end if

         end if

*-----------------------------------------------------------------------
*        read magnetic field informations
*-----------------------------------------------------------------------

               nereg = nereg + 1

               if( nereg .gt. kvlmax ) goto 998

                  s_elf(nereg)  = 0.0d0
                  s_mgf(nereg)  = 0.0d0
                  t_elf(nereg)  = 0.0d0
                  t_mgf(nereg)  = 0.0d0
                  ! AdvanceSoft Hasemi 2019/11/18
                  emap_type(nereg) = 0.0d0
                  mmap_type(nereg) = 0.0d0
                  a_elmg(nereg)  = 0.0d0

            if( nereg .eq. 1 ) then

                  inerc = 0
                  iaddress_region(:) = 0
                  call moddas_reallocate_int(
     &                    2, 1, MAX_NUM_INERC, iaddress_region
     &                    , idas_inerc)
                  idsm  = inerc
                  jdsm  = 0

            end if

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( imsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( imsq(k) .eq. 1 ) then

                     ndsm = idsm + jdsm + 3

                  call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ndsm,ierr
     &                          ,MAX_NUM_INERC,idas_inerc)

                  if( ierr .ne. 0 ) goto 999

                     jdsm = jdsm + 1
                     idas_inerc(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                     idas_inerc(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn

            else if( imsq(k) .eq. 2 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               s_elf(nereg) = cvvv

            else if( imsq(k) .eq. 3 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               s_mgf(nereg) = cvvv

            else if( imsq(k) .eq. 4 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               t_elf(nereg) = cvvv

            else if( imsq(k) .eq. 5 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               t_mgf(nereg) = cvvv

            else if( imsq(k) .eq. 6 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               emap_type(nereg) = cvvv

            else if( imsq(k) .eq. 7 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               mmap_type(nereg) = cvvv

            else if( imsq(k) .eq. 8 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               a_elmg(nereg) = cvvv

            else if( imsq(k) .eq. 9 ) then

               call mapnum(chin,200,ic,ic1,ic2,ierr)

               if( ierr .ne. 0 ) goto 999

               map1(nereg) = chin(ic1:ic2)
               icmap1(nereg) = ic2-ic1+1

            else if( imsq(k) .eq. 10 ) then

               call mapnum(chin,200,ic,ic1,ic2,ierr)

               if( ierr .ne. 0 ) goto 999

               map2(nereg) = chin(ic1:ic2)
               icmap2(nereg) = ic2-ic1+1

            end if

         end do

         goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

               if( jdsm > MAX_NUM_INERC ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.elmgf@read02.f ?dimension over idas_inerc?'
     &                    //' jdsm > MAX_NUM_INERC'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INERC@moddas.f=',MAX_NUM_INERC,')'
                  ErrID = 'L:13736/R:elmgf/F:read02.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 2, 1, jdsm+1, iaddress_region, idas_inerc)

               kelcs = 0
               call moddas_allocate_dbl(nereg*9+1, das_kelcs)
               idsm  = kelcs
               jdsm  = 0


               ielfxyzlist=0 ! number of elf xyz list
               ielfrzlist=0  ! number of elf r-z list
               ielfxyzmap=0  ! number of elf xyz map
               ielfrzmap=0   ! number of elf r-z map
               imgfxyzlist=0 ! number of mgf xyz list
               imgfrzlist=0  ! number of mgf r-z list
               imgfxyzmap=0  ! number of mgf xyz map
               imgfrzmap=0   ! number of mgf r-z map

            do i = 1, nereg

               jdsm = jdsm + 1
               das_kelcs(idsm+jdsm) = s_elf(i)
               jdsm = jdsm + 1
               das_kelcs(idsm+jdsm) = s_mgf(i)
               jdsm = jdsm + 1
               das_kelcs(idsm+jdsm) = 0
               jdsm = jdsm + 1
               das_kelcs(idsm+jdsm) = t_elf(i)
               jdsm = jdsm + 1
               das_kelcs(idsm+jdsm) = 0
               jdsm = jdsm + 1
               das_kelcs(idsm+jdsm) = t_mgf(i)
               jdsm = jdsm + 1
               das_kelcs(idsm+jdsm) = emap_type(i)
               jdsm = jdsm + 1
               das_kelcs(idsm+jdsm) = mmap_type(i)
               jdsm = jdsm + 1
               das_kelcs(idsm+jdsm) = a_elmg(i)

               ! read electric map
               if(icmap1(i).gt.0) then
                  if(emap_type(i) .eq. -1) then
                     if(ielfxyzlist .eq. 0) then
                        call readelcxyzlist(map1(i),icmap1(i))
                        ielfxyzlist = 1
                     else
                        goto 991
                     endif
                  elseif(emap_type(i) .eq. -2) then
                     if(ielfrzlist .eq. 0) then
                        call readelcrzlist(map1(i),icmap1(i))
                        ielfrzlist = 1
                     else
                        goto 992
                     endif
                  elseif(emap_type(i) .eq. -3) then
                     if(ielfxyzmap .eq. 0) then
                        call readelcxyzmap(map1(i),icmap1(i))
                        ielfxyzmap = 1
                     else
                        goto 993
                     endif
                  elseif(emap_type(i) .eq. -4) then
                     if(ielfrzmap .eq. 0) then
                        call readelcrzmap(map1(i),icmap1(i))
                        ielfrzmap = 1
                     else
                        goto 994
                     endif
                  endif
               endif

               ! read magnetic map
               if(icmap2(i).gt.0) then
                  if(mmap_type(i) .eq. -1) then
                     if(imgfxyzlist .eq. 0) then
                        call readmagxyzlist(map2(i),icmap2(i))
                        imgfxyzlist = 1
                     else
                        goto 991
                     endif
                  elseif(mmap_type(i) .eq. -2) then
                     if(imgfrzlist .eq. 0) then
                        call readmagrzlist(map2(i),icmap2(i))
                        imgfrzlist = 1
                     else
                        goto 992
                     endif
                  elseif(mmap_type(i) .eq. -3) then
                     if(imgfxyzmap .eq. 0) then
                        call readmagxyzmap(map2(i),icmap2(i))
                        imgfxyzmap = 1
                     else
                        goto 993
                     endif
                  elseif(mmap_type(i) .eq. -4) then
                     if(imgfrzmap .eq. 0) then
                        call readmagrzmap(map2(i),icmap2(i))
                        imgfrzmap = 1
                     else
                        goto 994
                     endif
                  endif
               endif

            end do


         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:13859/R:elmgf/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

!AdvanceSoft Hasemi 2019/11/18
  991 continue

         m_err = 'Only 1 xyz electric field list is allowed in '//
     &           '[electro magnetic field]'
         ErrCha = ''
         ErrID = 'L:13871/R:elmgf/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  992 continue

         m_err = 'Only 1 r-z electric field list is allowed in '//
     &           '[electro magnetic field]'
         ErrCha = ''
         ErrID = 'L:13882/R:elmgf/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  993 continue

         m_err = 'Only 1 xyz electric field map is allowed in '//
     &           '[electro magnetic field]'
         ErrCha = ''
         ErrID = 'L:13893/R:elmgf/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  994 continue

         m_err = 'Only 1 r-z electric field map is allowed in '//
     &           '[electro magnetic field]'
         ErrCha = ''
         ErrID = 'L:13904/R:elmgf/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in '//
     &           '[electro magnetic field] is wrong.'
         ErrCha = ''
         ErrID = 'L:13917/R:elmgf/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         write(dkam,'(i6)') kvlmax
         m_err = 'Number of region for electro magnetic '//
     &           'field exceeds kvlmax = '// dkam
         ErrCha = ''
         ErrID = 'L:13931/R:elmgf/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [electro magnetic field] is wrong.'
         ErrCha = ''
         ErrID = 'L:13943/R:elmgf/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine sours(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,npcunt,ierr)
*                                                                      *
*       read [source] section of input files                           *
*       modified by K.Niita on 2012/12/13                              *
*                                                                      *
************************************************************************

      use COSMICMOD ! T.Sato 2021/05/25
      use moddas
      use moddas_mesh
      use moddas_region
      use moddas_source

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'angel01.inc'
C MATSUDA 2017.05.29
      include 'risrcparam.inc'

      include 'param-physcnst.inc' ! T.Sato 2023/07/23 for using electron mass

      include 'err.inc'

*     risrcparam.inc: include parameters
*     integer nuclinmax, nuclpumax, maxchain, chainmax
*     parameter ( nuclinmax =    100 ) nuclides including daughter
*     parameter ( nuclpumax =  1,000 ) Not used
*     parameter ( maxchain  =     23 ) from Number of chain (.NDX)
*     parameter ( chainmax  =  4,050 ) from Number of linear chain (.NDX)
*
*     integer rimax1, rimax2
*     parameter ( rimax1    =  4,000 ) from Number of auger electron (.ACK)
*     parameter ( rimax2    = 20,000 ) rimax1 * 5

*-----------------------------------------------------------------------
*     from angel01.inc for constants
*-----------------------------------------------------------------------

      common /rval1/ cval(mxcval), aval(mxcval)

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /mpi00/  npe, me
      common /paraj/ mstz(300), parz(300)

      common /isocor/ iscorr, itcorr, imlwt(isrc)
      common /isomul/ smlwt(isrc), totfact, imsrc
      common /isorfs/ ispfs(isrc), rspfn, rspfz, ispfn
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)
      common /isorori/ jstypori(isrc) ! T.Sato, original jstyp written in input file
      common /isorbia/ isbias(isrc)   ! T.Sato, source bias method for xyz-mesh source
      common /isbeam/ sx2(isrc),sy2(isrc),
     &       sxmrad1(isrc),sxmrad2(isrc),symrad1(isrc),symrad2(isrc)   ! T.Sato, beam emittance source

      common /isorsp/ sx0(isrc), sy0(isrc), sz0(isrc), sx1(isrc),
     &                sy1(isrc), sz1(isrc), sr0(isrc), se0(isrc),
     &                sdir(isrc), srx(isrc), sry(isrc), swem(isrc),
     &                sphi(isrc), sdom(isrc), swt0(isrc)
      common /isorsf/ isorf(isrc),lsfile(isrc), sfile(isrc)
      character sfile*100

      common /isorsn/ sr1(isrc), sr2(isrc), isrn(isrc)
      common /isorsr/ nsmx(isrc), nsrn(isrc), nsrc(isrc)
      common /isorsc/ isort(isrc,4), rsort(isrc,13)
      common /isorsd/ isdmp(isrc,0:30), jsdmp(isrc,0:30)
      common /isorsa/ sfactor(isrc)
      common /isorpn/ ssx(isrc), ssy(isrc), ssz(isrc)
!$OMP THREADPRIVATE(/isorpn/)
      common /isodct/ sdl0(isrc), sdl1(isrc), sdl2(isrc), sdpf(isrc),
     &                sdxw(isrc), sdyw(isrc), sdrd(isrc), sdebg(isrc),
     &                sdsxp(isrc), sdsxn(isrc), sdsyp(isrc),
     &                sdsyn(isrc), dnorm(isrc), prs1(isrc), prs2(isrc),
     &                prs3(isrc), prs4(isrc), psxp(isrc), psxn(isrc),
     &                psyp(isrc), psyn(isrc), sdls(isrc), sdrs(isrc),
     &                sdxs(isrc), sdys(isrc)

      common /isorsl/ nglp(isrc), ngli(isrc), ngla(isrc), ngfl(isrc),
     &                nglc(isrc), nglw(isrc)

      common /isousr/ jsusr(isrc,0:30)

*-----------------------------------------------------------------------

      common /isorse/ ngrp(isrc), ngei(isrc), ngea(isrc), ngfe(isrc),
     &                ngft(isrc), ngll(isrc), ngpi(isrc), ngpw(isrc)

      common /isorfx/ isnm(isrc), lsfx(isrc), srfx(isrc)
      character srfx*200


      common /isorsg/ seg0(isrc),seg1(isrc),seg2(isrc),seg3(isrc),
     &                set0(isrc),set1(isrc),set2(isrc),set3(isrc),
     &                jetyp(isrc),jptyp(isrc)

*-----------------------------------------------------------------------

      common /isoraa/ narp(isrc), naei(isrc), naea(isrc), nafe(isrc),
     &                naft(isrc), nall(isrc), napi(isrc), napw(isrc)

      common /isorfa/ isnn(isrc), lsfy(isrc), srfy(isrc)
      character srfy*200


      common /isoras/ sag1(isrc),sag2(isrc),jatyp(isrc),jqtyp(isrc)
      common /isospg/ spg1(isrc),spg2(isrc) ! T.Sato 2021/02/26

*-----------------------------------------------------------------------

      common /isorsm/ stm0(isrc),stmw(isrc),stmc(isrc),stmd(isrc),
     &                jttyp(isrc),jttpn(isrc)

      common /isortt/ ntrp(isrc), ntei(isrc), ntea(isrc), ntfe(isrc),
     &                ntft(isrc), ntll(isrc), ntpi(isrc), ntpw(isrc)

      common /isorft/ isll(isrc), lsfz(isrc), srfz(isrc)
      character srfz*200


      common /isorts/ stg1(isrc),stg2(isrc),jotyp(isrc)

*-----------------------------------------------------------------------

      common /isormd/ istdd(isrc),istcc(isrc),
     &                isxtp(isrc),isinx(isrc),istxx(isrc),
     &                isytp(isrc),isiny(isrc),istyy(isrc),
     &                isztp(isrc),isinz(isrc),istzz(isrc),
     &                sxmin(isrc),sxmax(isrc),sxdel(isrc),
     &                symin(isrc),symax(isrc),sydel(isrc),
     &                szmin(isrc),szmax(isrc),szdel(isrc)

*-----------------------------------------------------------------------
      common /itetsor/ ksoutnode,itetreg(isrc)

*-----------------------------------------------------------------------

      common /isosuf/ issuf(isrc), iscut(isrc), isvct(isrc,8),
     &                issfd(isrc), isvfd(isrc,8),
     &                ivsfd(isrc), ivvfd(isrc,8),
     &                dvsfd(isrc,4), dvvfd(isrc,8,4)

*-----------------------------------------------------------------------
      common /ibchsor/ ibcsn(isrc), ibsor(isrc,isrc)
      integer :: ibcsn = 0
      logical deqn1

*-----------------------------------------------------------------------

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------
cFURUTA20150515 78=>80  idmpmode and dmpmulti are added
cKN 2015/11/26 add 79-91, sift idmpmode 92 and dmpmulti 93 after those
C MATSUDA 2016.07.31 94=>97  for RI source

cKN 2016/12/29 97=>98
C MATSUDA 2017.05.29 98=>99  for RI source (Auger)
cFURUTA20180115 99=>100 for tetra source
C MATSUDA 2018.08.15 100=>101  for RI source (character X-rays)
C MATSUDA 2019.05.08 101=>102  for RI source (annihilation photons)
! T.Sato 2019.05.19 102=>105  for counter
! T.Sato 2020/09/21, 109 -> 113 for beam source
cFURUTA20201110 113->114
! T.Sato 2020/12/07, cosmic-ray source 114->125
! T.Sato 2021/02/26, cosmic-ray source (phi limitation) 125 -> 127
cfrtati 2021/03/06, 127 -> 128 for ibatch
*-----------------------------------------------------------------------
      parameter(icsu=128)

      dimension lschn(icsu), ischn(icsu), spava(icsu)
      character schan(icsu)*8

      data ( schan(i), i = 1, icsu ) /
     &    's-type  ','proj    ','x0      ','x1      ','y0      ',
     &    'y1      ','z0      ','z1      ','r0      ','e0      ',
     &    'dir     ','rx      ','ry      ','wem     ','rn      ',
     &    'rb      ','ne      ','e-type  ','eg0     ','eg1     ',
     &    'eg2     ','eg3     ','et0     ','et1     ','et2     ',
     &    'file    ','phi     ','dom     ','t-type  ','t0      ',
     &    'tw      ','tc      ','tn      ','td      ','r1      ',
     &    'r2      ','reg     ','ntmax   ','wgt     ','trcl    ',
     &    '*trcl   ','dump    ','factor  ','sx      ','sy      ',
     &    'sz      ','dl1     ','dl2     ','dpf     ','dxw     ',
     &    'dyw     ','drd     ','dl0     ','debg    ','dsxp    ',
     &    'dsxn    ','dsyp    ','dsyn    ','dls     ','drs     ',
     &    'dxs     ','dys     ','p-type  ','nl      ','<source>',
     &    'totfact ','f(x)    ','nm      ','a-type  ','na      ',
     &    'ag1     ','ag2     ','nn      ','g(x)    ','q-type  ',
     &    'ispfs   ','iscorr  ','izst    ','ntt     ','tg1     ',
     &    'tg2     ','ll      ','h(x)    ','o-type  ','x2      ',
     &    'y2      ','z2      ','x3      ','y3      ','z3      ',
     &    'exa     ','idmpmode','dmpmulti','ni      ','actlow  ',
     &    'dtime   ','norm    ','mesh    ','iaugers ','tetreg  ',
     &    'icharctx','iannih  ','cnt(1)  ','cnt(2)  ','cnt(3)  ',
     &    'suf     ','cut     ','et3     ','isbias  ','xmrad1  ',
     &    'xmrad2  ','ymrad1  ','ymrad2  ','jpsf    ','icyear  ',
     &    'icmonth ','icday   ','glat    ','glong   ','alti    ',
     &    'environ ','icenv   ','solarmod','rigid   ','depatom ',
     &    'pg1     ','pg2     ','ibatch  '/

      data ( lschn(i), i = 1, icsu ) /  ! character length of parameter
     &     6,         4,         2,         2,         2,
     &     2,         2,         2,         2,         2,
     &     3,         2,         2,         3,         2,
     &     2,         2,         6,         3,         3,
     &     3,         3,         3,         3,         3,
     &     4,         3,         3,         6,         2,
     &     2,         2,         2,         2,         2,
     &     2,         3,         5,         3,         4,
     &     5,         4,         6,         2,         2,
     &     2,         3,         3,         3,         3,
     &     3,         3,         3,         4,         4,
     &     4,         4,         4,         3,         3,
     &     3,         3,         6,         2,         8,
     &     7,         4,         2,         6,         2,
     &     3,         3,         2,         4,         6,
     &     5,         6,         4,         3,         3,
     &     3,         2,         4,         6,         2,
     &     2,         2,         2,         2,         2,
     &     3,         8,         8,         2,         6,
     &     5,         4,         4,         7,         6,
     &     8,         6,         6,         6,         6,
     &     3,         3,         3,         6,         6,
     &     6,         6,         6,         4,         6,
     &     7,         5,         4,         5,         4,
     &     7,         5,         8,         5,         7,
     &     3,         3,         6/

*-----------------------------------------------------------------------

      dimension kstyp(6), knkf0(6)

      character chin*200, chlw*200, chcm*200
      character chlc*200
      character chlg*200
      character chla*200

      character chlt*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)
      character dkam*9

      logical exex

      dimension ssort(13)

      data imsrc  / 0 /
      data ispfn  / 0 /
      data rspfn  / 0.0d0 /
      data rspfz  / 0.0d0 /

      data nsrc   /isrc*0/ ! hirata 2023/06/17
      dimension cval0(mxcval),cval1(mxcval),cval2(mxcval),cval3(mxcval)

*-----------------------------------------------------------------------

      character chme*5, filnm*100

*-----------------------------------------------------------------------

      data bign / -1.d10 /

*-----------------------------------------------------------------------
      common /stat / istdev, irestart, ireschk
      real(8) dmpmulti
      integer idmpmode,ibchjmp,idmpjmp
      common /stat2/ dmpmulti,idmpmode,ibchjmp,idmpjmp(2)
      integer idmpsors
      integer jpsf
      common /stat3/ jpsf

*-----------------------------------------------------------------------
C MATSUDA 2016.07.31  for RI source
      common /risource/ chalct(nuclinmax,isrc), act000(nuclinmax,isrc),
     &                  activy(nuclinmax,isrc), thf(nuclinmax,isrc),
     &                  decayt(isrc)
      real  chalct
      double precision  act000, activy, thf, decayt

C MATSUDA 2017.05.29 (iaugers)
C MATSUDA 2018.08.15 (icharacterx)
      common /risrc00/ niorg(isrc), nicur(isrc), norm(isrc),
     &                 iaugers(isrc), icharacterx(isrc), iannih(isrc),
     &                 aclow(isrc)
      integer  niorg, nicur, norm, iaugers, icharacterx, iannih
      double precision  aclow
      common /risrc01/ normfact(isrc), asfsum(isrc)
      double precision  normfact, asfsum

C T.Sato 2021/09/01 necessary for cosmic-ray source
      common /risrc02/ smlwt2(isrc), totfact2

      real chza
      character chaa*7
      integer j33, nk0, nk1
      double precision karival
      double precision, allocatable :: e_lowk(:), e_upk(:), ratkari(:)
      double precision, allocatable :: prbkari(:)

      common /cntsrc/icntsrc(3,isrc) ! initial counter

! T.Sato 2020/12/07 for cosmic-ray source
      common /cosmicint/ipcosmic(isrc),icenv(isrc)
      common /cosmicreal/solarmod(isrc),rigid(isrc),
     &        depatom(isrc),ground(isrc),environ(isrc)
      common /cosmicani/ie511(isrc),annihratio(isrc) ! T.Sato 2023/07/23

CCSE set proj=all flag (2022.03.02) >>>>>

! Nais_2024 >>>
! T-cross Dump
      dimension isdmp0(0:30)
! Nais_2024 <<<

      iprojall = npcunt

      if( iprojall .eq. -1 ) then

         open(niws,form='formatted',status='scratch')   ! opened in source subroutine


         npcunt = 0   ! total incident particles
         iripmode = 0   ! process mode for ripmake sub. 1:<> subsection, 2:[] section
         ipaflg = 0   ! flag of "proj=all" in <source> sub-section
         issflg = 0   ! flag of <source> subsection input

      endif
CCSE set proj=all flag (2022.03.02) <<<<<
*-----------------------------------------------------------------------

            ierr  = 0
            imrd  = 0
            imlcnt  = 0
            totfact = 1.0d0
            smlwts  = 1.0d0
            iscorr = 0
            idmpsors = 0 !FURUTA20150515
            ksoutnode = 0 !FURUTA20180115
            ierrMSG = 0 ! T.Sato 2017/07/13
            icntsrc(:,:)=0 ! T.Sato 2019/05/19
cABE 2023/05/25, change the default value of iannih = 1
            iannih(:) = 1

*-----------------------------------------------------------------------
*     multi-source: j
*-----------------------------------------------------------------------

 5000 continue

               imsrc = imsrc + 1

               if( imsrc .gt. isrc ) then
                write(ErrCha,'("Number of multi-source is larger than",
     &          " isrc. Please increase isrc written in param.inc and",
     &          " recompile PHITS")')
                ErrID = 'L:14322/R:sours/F:read02.f' !E03_006_001
                call ErrWrite(ErrID,ErrCha)
                goto 962
               endif

               j     = imsrc
               smlwt(j) = smlwts
               ie511(j)=0 ! T.Sato 2023/07/23

               imlcnt   = 0

*-----------------------------------------------------------------------

         do i = 1, icsu

            ischn(i) = 0
            spava(i) = 0.0d0

         end do

            nsrn(j) = 0

            spava(11) = 1.0 ! T.Sato 2017/07/13, default value for dir

            spava(27) = -10000.0
            spava(28) = -10000.0
            spava(29) = 0.0
            spava(30) = 0.0
            spava(33) = 1.0
            spava(34) = 0.0
            spava(38) = 1000.0
            spava(39) = 1.0d0
            spava(43) = 1.0d0
            spava(44) = 0.0d0
            spava(45) = 0.0d0
            spava(46) = 0.0d0
            spava(64) = -1.0d0
            spava(69) =  0.0d0
            spava(70) = -1.0d0
            spava(78) = -10000.0d0
            spava(95) = 1.0d-10
            spava(96) = -10.0d0  ! T.Sato 2016/08/22
            spava(97) = 0.0d0
            spava(99) = 0.0
           spava(101) = 0.0
cABE 2023/05/25, change the default value of iannih = 1
           spava(102) = 1.0

            igkst  = 0
            idtt   = 0
            ktrs   = 0
            istyp(j)  = -10000

         do i = 0, 30
            isdmp(j,i) = 0
            jsdmp(j,i) = 0
            jsusr(j,i) = 0
         end do

            nglp(j) = spava(64)
           ispfs(j) = spava(76)

*-----------------------------------------------------------------------

            spava(106) = -1.0
            spava(107) =  0.0
            spava(108) =  0.5 ! et3 is the power index of Maxwellian, default is 0.5, T.Sato 2020/02/05
            incut = 0

            issuf(j) = -1
            iscut(j) = 0
         do i = 1, 8
            isvct(j,i) = 0
         end do

            jpsf = 0

*-------- Cosmic-ray source, T.Sato 2020/12/07 -------------------------
            spava(115)=2009 ! icyear
            spava(116)=10   ! icmonth
            spava(117)=20   ! icday
            spava(118)=90.0 ! glat (latitude in degree)
            spava(119)=0.0 ! glong (longitude in degree)
            spava(120)=0.0  ! alti (altitude in km)
            spava(121)=0.15  ! ground (water density or aircraft mass)
            spava(122)=0     ! environment, <0:SEP, 0:GCR, 1:ideal atmosphere, 2:on ground, 3:aircraft (pilot), 4:aircraft(cabin), 5:blackhole
            spava(123)=0     ! solarmod (FFP in MV)
            spava(124)=0     ! rigid (cut-off rigidity in GV)
            spava(125)=0     ! depatom (atmospheric depth in g/cm2)

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800

               if( iskip .ne. 0 ) goto 140

  150 continue

*-----------------------------------------------------------------------
*        end of source section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               if( iprojall .eq. -1 ) then

                  iripmode = 2
                  call ripmake(imsrc,npcunt,iripmode,ipaflg,issflg,chin,
     &                         irptyp,ierr)

                  select case(irptyp)
                     case (1)   ! photon
                        istyp(j) = 14
                        inkf0(j) = 22
                        kstyp(1) = 14
                        knkf0(1) = 22
                     case (2)   ! electron
                        istyp(j) = 12
                        inkf0(j) = 11
                        kstyp(1) = 12
                        knkf0(1) = 11
                     case (3)   ! positron
                        istyp(j) = 13
                        inkf0(j) = -11
                        kstyp(1) = 13
                        knkf0(1) = -11
                     case (4)   ! alpha
                        istyp(j) = 18
                        kstyp(1) = 18
                        inkf0(j) = 2000004
                        knkf0(1) = 2000004
                     case (5)   ! frtati 2022/12/26 SF neutron
                        istyp(j) = 2
                        inkf0(j) = 2112
                        kstyp(1) = 2
                        knkf0(1) = 2112
                  end select
                  ipaflg = 0

                  if( ierr .ne. 0 ) return

               endif

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

               goto 999

*-----------------------------------------------------------------------
*        read value of parameters
*-----------------------------------------------------------------------

  100    continue
               if( iprojall .eq. -1 ) then
                  if ( (i.ne.65) .and. (i.ne.66) ) then   ! not <source>, totfact
                     write(niws,'(a)') chin(1:i2)         ! work souce section
                  end if
               endif

               ipm = i

               ischn( ipm ) = 1

               ic = inumc(chlw,il+1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) then
                write(ErrCha,'("read value of parameters error")')
                ErrID = 'L:14519/R:sours/F:read02.f' !E03_007_001
                call ErrWrite(ErrID,ErrCha)
                goto 997
              endif

               icl = inumc(chlw,ic,i3,';') - 1

*-----------------------------------------------------------------------
*        total factor
*-----------------------------------------------------------------------

         if( ipm .eq. 66 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,'("total factor error")')
                ErrID = 'L:14536/R:sours/F:read02.f' !E03_008_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               totfact = cvvv

*-----------------------------------------------------------------------
*        multi correlation source
*-----------------------------------------------------------------------

         else if( ipm .eq. 77 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,'("multi correlation source error")')
                ErrID = 'L:14553/R:sours/F:read02.f' !E03_009_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               if( nint( cvvv ) .ne. 0 ) iscorr = 1

               if( nint( cvvv ) .eq. 2 ) iscorr = 2  ! T.Sato 2013/11/15 for correlation

               if( nint( cvvv ) .eq. 3 ) iscorr = 3  ! T.Sato 2013/11/15 for correlation

*-----------------------------------------------------------------------
*        multi-source <source> subsection
*-----------------------------------------------------------------------

         else if( ipm .eq. 65 ) then

            if( iprojall .eq. -1 ) then

               iripmode = 1
               call ripmake(imsrc,npcunt,iripmode,ipaflg,issflg,chin,
     &                      irptyp,ierr)

               select case(irptyp)
                  case (1)   ! photon
                     istyp(j) = 14
                     inkf0(j) = 22
                     kstyp(1) = 14
                     knkf0(1) = 22
                  case (2)   ! electron
                     istyp(j) = 12
                     inkf0(j) = 11
                     kstyp(1) = 12
                     knkf0(1) = 11
                  case (3)   ! positron
                     istyp(j) = 13
                     inkf0(j) = -11
                     kstyp(1) = 13
                     knkf0(1) = -11
                  case (4)   ! alpha
                     istyp(j) = 18
                     kstyp(1) = 18
                  case (5)   ! frtati 2022/12/26 SF neutron
                     istyp(j) = 2
                     inkf0(j) = 2112
                     kstyp(1) = 2
                     knkf0(1) = 2112
               end select
               ipaflg = 0

               if( ierr .ne. 0 ) return
               issflg = 1
               write(niws,'(a)') chin(1:i2)      ! work source sub-section

            endif

               call onum(chlw,ic,icl,cvvv,ierr)

            if( ierr .ne. 0 ) then
             write(ErrCha,'("multi-source <source> subsection error")')
             ErrID = 'L:14613/R:sours/F:read02.f' !E03_010_001
             call ErrWrite(ErrID,ErrCha)
             goto 998
            endif

            if( ischn(1) .eq. 0 .and.
     &        ( imsrc .eq. 1 .or. imrd .eq. 0 ) ) then

               smlwt(imsrc) = cvvv
               imrd = imrd + 1

            else

               smlwts = cvvv
               imlcnt = 1
               goto 800

            end if

*-----------------------------------------------------------------------
*        input file name
*-----------------------------------------------------------------------

         else if( ipm .eq. 26 ) then

               ic1 = jnumc(chlw,ic,icl)
               ic2 = min( i3, inumc(chlw,ic+1,icl,' ') )

               lsfile(j) = ic2 - ic1 + 1

               do i = 1, lsfile(j)

                  sfile(j)(i:i) = chin(ic1+i-1:ic1+i-1)

               end do

               do i = lsfile(j) + 1, 100

                  sfile(j)(i:i) = ' '

               end do

*-----------------------------------------------------------------------
*        projectile name
*-----------------------------------------------------------------------

         else if( ipm .eq. 2 ) then

               ic = jnumc(chlw,ic,icl)

            call rdpname(ic,icl,chlw,istyp(j),inkf0(j),kstyp,knkf0,ierr)

               if( ierr .eq. 994 ) then
                write(ErrCha,'("projectile name error")')
                ErrID = 'L:14667/R:sours/F:read02.f' !E03_011_001
                call ErrWrite(ErrID,ErrCha)
                goto 983
               endif
               if( ierr .eq. 998 ) then
                write(ErrCha,'("projectile value error")')
                ErrID = 'L:14673/R:sours/F:read02.f' !E03_011_002
                call ErrWrite(ErrID,ErrCha)
                goto 983
               endif
               if( istyp(j) .lt. -1 ) then
                write(ErrCha,*)
     & 'Error of s-type. The setting is less than -1.'
                ErrID = 'L:14680/R:sours/F:read02.f' !E03_012_001
                call ErrWrite(ErrID,ErrCha)
                goto 964
               endif

               if( istyp(j) .lt. 0 ) then

                  istyp(j) = kstyp(1)
                  inkf0(j) = knkf0(1)

               end if

               ipaflg = 1
               if( iprojall .eq. -1 .and. istyp(j) .eq. 20  ) then  ! proj=all
                  ipaflg = 2
               else if( istyp(j) .ge. 19 .and. inkf0(j) .eq. 0 ) then
                  write(ErrCha,*)
     & 'Error of s-type. The setting is 19 or more.'
                  ErrID = 'L:14698/R:sours/F:read02.f' !E03_012_002
                  call ErrWrite(ErrID,ErrCha)
                  goto 983
               endif

*-----------------------------------------------------------------------
*        direction
*-----------------------------------------------------------------------

         else if( ipm .eq. 11 ) then

               ic = jnumc(chlw,ic,icl)

            if( chlw(ic:ic+2) .eq. 'all' ) then

               spava(11) = 200.0

            else if( chlw(ic:ic+3) .eq. '+all' ) then

               spava(11) = 2.0

            else if( chlw(ic:ic+3) .eq. '-all' ) then

               spava(11) = -2.0

            else if( chlw(ic:ic+2) .eq. 'iso' ) then

               spava(11) = -3.0

            else if( chlw(ic:ic+1) .eq. 'cr' ) then  ! T.Sato 2024/06/18, cosmic-ray angular distribution

               spava(11) = -4.0d0

            else if( chlw(ic:ic+3) .eq. 'data' ) then

               spava(11) = 300.0

            else

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,'("direction error")')
                ErrID = 'L:14741/R:sours/F:read02.f' !E03_013_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               if( abs( cvvv ) .gt. 1.0d0 ) then
                write(ErrCha,'("dir is greater than 1.0.")')
                ErrID = 'L:14748/R:sours/F:read02.f' !E03_013_002
                call ErrWrite(ErrID,ErrCha)
                goto 993
               endif

               spava(11) = cvvv

            end if

C S.H. added all option for phi (2018.10.2)
*-----------------------------------------------------------------------
*        phi (azimuthal angle)
*-----------------------------------------------------------------------

         else if( ipm .eq. 27 ) then

               ic = jnumc(chlw,ic,icl)

            if( chlw(ic:ic+2) .eq. 'all' ) then

               spava(27) = -1000.0

            else

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,'("phi (azimuthal angle) error")')
                ErrID = 'L:14776/R:sours/F:read02.f' !E03_014_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               if( cvvv  .le. -1000.0 ) then
                write(ErrCha,*)
     & 'phi should be greater than -1000 in [source].'
                ErrID = 'L:14784/R:sours/F:read02.f' !E03_014_002
                call ErrWrite(ErrID,ErrCha)
                goto 930
               endif

               spava(27) = cvvv

            end if

*-----------------------------------------------------------------------
*        region
*-----------------------------------------------------------------------

         else if( ipm .eq. 37 ) then

               ndsm = 1
               call moddas_allocate_int(
     &                 MAX_NUM_NSRC, idas_nsrc_temporary)

               call tregion(1,jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      nsrn(j),msrn,ndsm,nvol,ivl,irvl,1
     &                      ,MAX_NUM_NSRC,idas_nsrc_temporary)

                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 800

                  goto 150

*-----------------------------------------------------------------------
*        transform
*-----------------------------------------------------------------------

         else if( ipm .eq. 40 .or. ipm .eq. 41 ) then

                  if( ipm .eq. 41 ) ktrs = 1

               call ttrans(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     ic,ktrs,igkst,idtt,ssort)

                  do k = 1, 13

                     rsort(j,k) = ssort(k)

                  end do

                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 800

                  goto 150

*-----------------------------------------------------------------------
*        dump
*-----------------------------------------------------------------------

         else if( ipm .eq. 42 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,'("dump error")')
                ErrID = 'L:14846/R:sours/F:read02.f' !E03_015_001
                call ErrWrite(ErrID,ErrCha)
                goto 982
               endif
               if( jpsf .gt. 0)then
                write(ErrCha,'("dump error")')
                ErrID = 'L:14852/R:sours/F:read02.f' !E03_015_007
                call ErrWrite(ErrID,ErrCha)
                goto 936
               endif

               isdmp(j,0) = nint( cvvv )

            if( isdmp(j,0) .ne. 0 ) then

  151                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( ierr .ne. 0 ) return
                        if( jpn  .eq. 3 ) then
                         write(ErrCha,'("dump error")')
                         ErrID = 'L:14867/R:sours/F:read02.f' !E03_015_002
                         call ErrWrite(ErrID,ErrCha)
                         goto 982
                        endif

                        if( iskip .ne. 0 ) goto 151

                        if( iprojall .eq. -1 ) then
                           write(niws,'(a)') chin(1:i2)      ! work souce section
                        endif

                  ic = i1

               do k = 1, abs( isdmp(j,0) )

                  if( ic .gt. i3 ) then

  152                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                      write(ErrCha,'("dump error")')
                      ErrID = 'L:14890/R:sours/F:read02.f' !E03_015_003
                      call ErrWrite(ErrID,ErrCha)
                      goto 982
                     endif

                     if( iskip .ne. 0 ) goto 152

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

                  end if

                     ic = jnumc(chlw,ic,i3)

                     call snum(chlw,ic,i3,ic2,cvvv,ierr)

                     if( ierr .ne. 0 ) then
                      write(ErrCha,'("dump error")')
                      ErrID = 'L:14911/R:sours/F:read02.f' !E03_015_004
                      call ErrWrite(ErrID,ErrCha)
                      goto 982
                     endif

                     isdmp(j,k) = nint( cvvv )

                     if( isdmp(j,k) .gt. 20 .or.
     &                   isdmp(j,k) .le.  0 ) then
                        write(ErrCha,'("dump error")')
                        ErrID = 'L:14921/R:sours/F:read02.f' !E03_015_005
                        call ErrWrite(ErrID,ErrCha)
                        goto 982
                     endif

                     jsdmp(j,isdmp(j,k)) = k

                     ic = ic2

               end do

! Nais_2024 >>>
            else

               call read_dump_parameters(j,isdmp0)

               isdmp(j,0) = isdmp0(0)
               do k = 1, abs( isdmp(j,0) )
                  isdmp(j,k) = isdmp0(k)
                  jsdmp(j,isdmp(j,k)) = k
               end do
! Nais_2024 <<<

            end if
            if(idmpjmp(1).eq.0)then
             if(jsdmp(j,18).gt.0.and.jsdmp(j,19).gt.0)then
              idmpmode=1
             else
              idmpmode=0
             endif
            endif
            idmpsors=1
c------------------
               goto 140

*-----------------------------------------------------------------------
*        jpsf
*-----------------------------------------------------------------------
         else if( ipm .eq. 114 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,'("dump error")')
                ErrID = 'L:14965/R:sours/F:read02.f' !E03_015_006
                call ErrWrite(ErrID,ErrCha)
                goto 937
               endif

               jpsf = nint( cvvv )

               if(jpsf.gt.0)then
                if( ischn(42) .eq. 1)then
                 write(ErrCha,'("dump error")')
                 ErrID = 'L:14975/R:sours/F:read02.f' !E03_015_007
                 call ErrWrite(ErrID,ErrCha)
                 goto 936
                endif
                isdmp(j,0) = 6 + jpsf
                isdmp(j,1) = 1
                isdmp(j,2) = 2
                isdmp(j,3) = 3
                isdmp(j,4) = 5
                isdmp(j,5) = 6
                isdmp(j,6) = 7
                isdmp(j,7) = 8
                if(jpsf.eq.2) isdmp(j,8) = 9
                jsdmp(j,1) = 1
                jsdmp(j,2) = 2
                jsdmp(j,3) = 3
                jsdmp(j,5) = 4
                jsdmp(j,6) = 5
                jsdmp(j,7) = 6
                jsdmp(j,8) = 7
                if(jpsf.eq.2) jsdmp(j,9) = 8
                idmpmode=0
                idmpsors=1
               endif
               goto 140

*-----------------------------------------------------------------------
*        f(x) information for energy function
*-----------------------------------------------------------------------

         else if( ipm .eq. 67 ) then

               lsfx(j) = icl - ic + 1

            do i = 1, lsfx(j)

               srfx(j)(i:i) = chin(ic+i-1:ic+i-1)
                  chlg(i:i) = chlw(ic+i-1:ic+i-1)

            end do

            do i = 1, mxcval

               cval1(i) = cval(i)

            end do

*-----------------------------------------------------------------------
*        g(x) information for angle function
*-----------------------------------------------------------------------

         else if( ipm .eq. 74 ) then

               lsfy(j) = icl - ic + 1

            do i = 1, lsfy(j)

               srfy(j)(i:i) = chin(ic+i-1:ic+i-1)
                  chla(i:i) = chlw(ic+i-1:ic+i-1)

            end do

            do i = 1, mxcval

               cval2(i) = cval(i)

            end do

*-----------------------------------------------------------------------
*        h(x) information for time function
*-----------------------------------------------------------------------

         else if( ipm .eq. 83 ) then

               lsfz(j) = icl - ic + 1

            do i = 1, lsfz(j)

               srfz(j)(i:i) = chin(ic+i-1:ic+i-1)
                  chlt(i:i) = chlw(ic+i-1:ic+i-1)

            end do

            do i = 1, mxcval

               cval3(i) = cval(i)

            end do

*-----------------------------------------------------------------------
*        nm : energy mesh for e-type = 3, 13, 5, 15, 6, 16
*-----------------------------------------------------------------------

         else if( ipm .eq. 68 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,*)
     & 'There is an error in e-type = 3, 13, 5, 15, 6, 16.'
                ErrID = 'L:15075/R:sours/F:read02.f' !E03_016_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               spava( ipm ) = cvvv

               if( ischn(18) .eq. 0 ) then
                write(ErrCha,*)
     & 'nm : energy mesh for e-type = 3, 13, 5, 15, 6, 16 error'
                ErrID = 'L:15085/R:sours/F:read02.f' !E03_016_002
                call ErrWrite(ErrID,ErrCha)
                goto 967
               endif

               ngrp(j)   = abs( nint( spava(68) ) )

               ngll(j)  = 1
               if( spava(68) .lt. 0.d0 ) ngll(j) = -1

               call moddas_reallocate_dbl(isrc, j, ngrp(j), ngei, egmin)
               call moddas_reallocate_dbl(isrc, j, ngrp(j), ngea, egmax)
               call moddas_reallocate_dbl(isrc, j, ngrp(j), ngfe, fegrp)
               call moddas_reallocate_dbl(isrc, j, ngrp(j), ngft, rfe)
               call moddas_reallocate_dbl(isrc, j, ngrp(j), ngpi, prw)
               call moddas_reallocate_dbl(isrc, j, ngrp(j), ngpw, pwt)

*-----------------------------------------------------------------------
*        nn : angle or cos mesh for a-type = 5, 6, 15, 16
*-----------------------------------------------------------------------

         else if( ipm .eq. 73 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,*)
     & 'angle or cos mesh for a-type = 5, 6, 15, 16 error.'
                ErrID = 'L:15113/R:sours/F:read02.f' !E03_017_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               spava( ipm ) = cvvv

               if( ischn(69) .eq. 0 ) then
                write(ErrCha,*)
     & 'angle or cos mesh for a-type = 5, 6, 15, 16 error.'
                ErrID = 'L:15123/R:sours/F:read02.f' !E03_017_002
                call ErrWrite(ErrID,ErrCha)
                goto 967
               endif

               if( abs( nint( spava(73) ) ) .le. 1 ) then
                write(*,*) ' *** Warning: '//
     &                   'Angular mesh of multisource No.', j,
     &                   ' is 1. nn is Replaced with 100.'
                narp(j) = 100

               else

                narp(j)   = abs( nint( spava(73) ) )

               endif

               nall(j)  = 1
               if( spava(73) .lt. 0.d0 ) nall(j) = -1

               call moddas_reallocate_dbl(isrc, j, narp(j), naei, agmin)
               call moddas_reallocate_dbl(isrc, j, narp(j), naea, agmax)
               call moddas_reallocate_dbl(isrc, j, narp(j), nafe, fagrp)
               call moddas_reallocate_dbl(isrc, j, narp(j), naft, rfa)
               call moddas_reallocate_dbl(isrc, j, narp(j), napi, paw)
               call moddas_reallocate_dbl(isrc, j, narp(j), napw, pat)


*-----------------------------------------------------------------------
*        ll : time for t-type = 5, 6
*-----------------------------------------------------------------------

         else if( ipm .eq. 82 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,'("ll : time for t-type = 5, 6 error")')
                ErrID = 'L:15161/R:sours/F:read02.f' !E03_018_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               spava( ipm ) = cvvv

               if( ischn(29) .eq. 0 ) then
                write(ErrCha,'("ll : time for t-type = 5, 6 error")')
                ErrID = 'L:15170/R:sours/F:read02.f' !E03_018_002
                call ErrWrite(ErrID,ErrCha)
                goto 960
               endif

               if( abs( nint( spava(82) ) ) .le. 1 ) then
                write(*,*) ' *** Warning: '//
     &                   'Time mesh of multisource No.', j,
     &                   ' is 1. ll is Replaced with 100.'
                ntrp(j) = 100

               else

                ntrp(j)   = abs( nint( spava(82) ) )

               endif

               ntll(j)  = 1
               if( spava(82) .lt. 0.d0 ) ntll(j) = -1

               call moddas_reallocate_dbl(isrc, j, ntrp(j), ntei, tgmin)
               call moddas_reallocate_dbl(isrc, j, ntrp(j), ntea, tgmax)
               call moddas_reallocate_dbl(isrc, j, ntrp(j), ntfe, ftgrp)
               call moddas_reallocate_dbl(isrc, j, ntrp(j), ntft, rft)
               call moddas_reallocate_dbl(isrc, j, ntrp(j), ntpi, ptw)
               call moddas_reallocate_dbl(isrc, j, ntrp(j), ntpw, ptt)

*-----------------------------------------------------------------------
*        idmpmode : event by event statistical processing for dump source
*                 (D=1 when dump source contains NOCAS and NOBCH)
*
*                 =1: On
*                 =0: Off
* FURUTA20150515
*-----------------------------------------------------------------------

         else if( ipm .eq. 92 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,'("idmpmode error")')
                ErrID = 'L:15212/R:sours/F:read02.f' !E03_019_001
                call ErrWrite(ErrID,ErrCha)
                goto 982
               endif

                idmpmode = nint( cvvv )

                idmpjmp(1)=1

                if(idmpmode.gt.0)then
                 if(imsrc.gt.1) then
                  write(ErrCha,'("idmpmode error")')
                  ErrID = 'L:15224/R:sours/F:read02.f' !E03_019_002
                  call ErrWrite(ErrID,ErrCha)
                  goto 956
                 endif

                 if(irestart.eq.1) then
                 write(ErrCha,*)
     & 'istdev < 0 is not allowed when idmpmode = 1 ',
     & 'or dmpmulti not= 0.0.'
                  ErrID = 'L:15233/R:sours/F:read02.f' !E03_019_003
                  call ErrWrite(ErrID,ErrCha)
                  goto 957
                 endif

                 if(idmpjmp(2).eq.1.and.dmpmulti.eq.0.0) then
                  write(ErrCha,*)
     & 'dmpmul = 0.0 is not allowed when idmpmode = 1.'
                  ErrID = 'L:15241/R:sours/F:read02.f' !E03_019_004
                  call ErrWrite(ErrID,ErrCha)
                  goto 960
                 endif
                endif

*-----------------------------------------------------------------------
*        dmpmulti : Multiplication factor for # of dump souce with
*                   proper statistical treatment for two-step calculation
*                   Multiplication of number of decimal places is realised
*                   by russian roulette
*
*                 = 0.0 special option
*                       dump source is recursively used till the end of
*                       calculation
* FURUTA20150515
*-----------------------------------------------------------------------

         else if( ipm .eq. 93 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,'("dmpmulti error")')
                ErrID = 'L:15265/R:sours/F:read02.f' !E03_020_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               dmpmulti =  cvvv

               idmpjmp(2)=1

               if(dmpmulti.lt.0.0) then
                write(ErrCha,'("dmpmulti should be positive or =0.")')
                ErrID = 'L:15276/R:sours/F:read02.f' !E03_020_002
                call ErrWrite(ErrID,ErrCha)
                goto 952
               endif

               if(dmpmulti.ne.0.0.and.irestart.eq.1) then
                write(ErrCha,*)
     & 'istdev < 0 is not allowed when idmpmode = 1 or ',
     & 'dmpmulti not= 0.0.'
                ErrID = 'L:15285/R:sours/F:read02.f' !E03_020_003
                call ErrWrite(ErrID,ErrCha)
                goto 957
               endif

               if(dmpmulti.eq.0.0.and.idmpjmp(1).eq.1
     &              .and.idmpmode.gt.0) then
        write(ErrCha,*)
     & 'dmpmul = 0.0 is not allowed when idmpmode = 1.'
                ErrID = 'L:15294/R:sours/F:read02.f' !E03_020_004
                call ErrWrite(ErrID,ErrCha)
               goto 960
               endif

*-----------------------------------------------------------------------
*        Lower limit of activity (Bq)         by MATSUDA 2016.07.31
*-----------------------------------------------------------------------
         else if( ipm .eq. 95 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,'("Lower limit of activity (Bq) error")')
                ErrID = 'L:15308/R:sours/F:read02.f' !E03_021_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               aclow(j)  =   cvvv

*-----------------------------------------------------------------------
*        decay time (cooling time)            by MATSUDA 2016.07.31
*-----------------------------------------------------------------------
         else if( ipm .eq. 96 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,'("decay time (cooling time) error")')
                ErrID = 'L:15324/R:sours/F:read02.f' !E03_022_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               decayt(j)  =  cvvv

*-----------------------------------------------------------------------
*        normalization per Bq(0) / photon(1)  by MATSUDA 2016.07.31
*-----------------------------------------------------------------------
         else if( ipm .eq. 97 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,*)
     & 'normalization per Bq(0) / photon(1) error'
                ErrID = 'L:15341/R:sours/F:read02.f' !E03_023_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               if( nint( cvvv ) .ne. 0 .and.
     &             nint( cvvv ) .ne. 1 ) then
                write(ErrCha,*)
     & 'normalization per Bq(0) / photon(1) error'
                ErrID = 'L:15350/R:sours/F:read02.f' !E03_023_002
                call ErrWrite(ErrID,ErrCha)
               goto 998
               endif

               norm(j)  = nint( cvvv )

*-----------------------------------------------------------------------
*        Auger electron with(0), w/o(1), only(2)  by MATSUDA 2017.05.29
*-----------------------------------------------------------------------
         else if( ipm .eq. 99 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,*)
     & 'Auger electron with(0) , w/o(1), only(2) error'
                ErrID = 'L:15367/R:sours/F:read02.f' !E03_024_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif
               if( nint( cvvv ) .ne. 0 .and.
     &             nint( cvvv ) .ne. 1 .and.
     &             nint( cvvv ) .ne. 2 ) then
               write(ErrCha,*)
     & 'Auger electron with(0), w/o(1), only(2) error'
                ErrID = 'L:15376/R:sours/F:read02.f' !E03_024_002
                call ErrWrite(ErrID,ErrCha)
               goto 998
               endif

               iaugers(j)  = nint( cvvv )

*-----------------------------------------------------------------------
*        Character X-rays with(0), w/o(1), only(2) by MATSUDA 2018.08.15
*-----------------------------------------------------------------------
         else if( ipm .eq. 101 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,*)
     & 'Character X-rays with(0) , w/o(1), only(2) error'
                ErrID = 'L:15393/R:sours/F:read02.f' !E03_025_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif
               if( nint( cvvv ) .ne. 0 .and.
     &             nint( cvvv ) .ne. 1 .and.
     &             nint( cvvv ) .ne. 2 ) then
                write(ErrCha,*)
     & 'Character X-rays with(0), w/o(1), only(2) error'
                ErrID = 'L:15402/R:sours/F:read02.f' !E03_025_002
                call ErrWrite(ErrID,ErrCha)
               goto 998
             endif

               icharacterx(j)  = nint( cvvv )

*-----------------------------------------------------------------------
*        Annihilation photons with(0), w/o(1) by MATSUDA 2019.05.08
*-----------------------------------------------------------------------
         else if( ipm .eq. 102 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,*)
     & 'Annihilation photons with(0) , w/o(1) error'
                ErrID = 'L:15419/R:sours/F:read02.f' !E03_025_003
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif
               if( nint( cvvv ) .ne. 0 .and.
     &             nint( cvvv ) .ne. 1 ) then
                write(ErrCha,*)
     & 'Annihilation photons with(0), w/o(1) error'
                ErrID = 'L:15427/R:sours/F:read02.f' !E03_025_004
                call ErrWrite(ErrID,ErrCha)
               goto 998
             endif

               iannih(j)  = nint( cvvv )

*-----------------------------------------------------------------------
*        Set initial counter by T.Sato 2019/05/19
*-----------------------------------------------------------------------
         else if(ipm.ge.103.and.ipm.le.105) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,*)
     & 'Something wrong in source counter setting'
                ErrID = 'L:15444/R:sours/F:read02.f' !E03_111_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               icntsrc(ipm-102,j)  = nint( cvvv )

*-----------------------------------------------------------------------
*        ibatch by frtati 2021/03/06
*-----------------------------------------------------------------------
         else if( ipm.eq.128 ) then

           ibcsn(j) = 0
           if( chlw(ic:ic+2) .eq. 'all' ) then
             ibcsn(j) = 0
             goto 157
           end if

            do itt = 1, nbchmax
             ic = jnumc(chlw,ic,i3)
             if( deqn1( chlw(ic:ic) ) ) then
               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 998
               ibcsn(j) = ibcsn(j) + 1
               ibsor(ibcsn(j),j) = nint( cvvv )
               ic = ic2
             else if ( chlw(ic:ic).eq.'{' ) then
               itc1 = ic+index(chlw(ic:),'-')-1
               itc2 = ic+index(chlw(ic:),'}')-1
               if ( itc1.eq.0 .or. itc2.eq.0 ) goto 998
               call onum(chlw,ic+1,itc1-1,cvvv,ierr)
               if( ierr.ne.0 ) goto 998
               itn1 = int(cvvv)
               call onum(chlw,itc1+1,itc2-1,cvvv,ierr)
               if ( ierr.ne.0 ) goto 998
               itn2 = int(cvvv)
               if ( itn2.lt.itn1 ) goto 998
               do ittt = 0, itn2-itn1
                 ibcsn(j) = ibcsn(j) + 1
                 ibsor(ibcsn(j),j) = itn1+ittt
               end do
               ic = itc2 + 1
             else if ( ibcsn(j).gt.0 ) then
               goto 157
             else
               goto 998
             end if
           end do

  157      continue

*-----------------------------------------------------------------------
*        mesh = xyz
*-----------------------------------------------------------------------

         else if( ipm .eq. 98 ) then

               if( chlw(ic:ic+2) .eq. 'reg' ) then

                  imesh = 1

               else if( chlw(ic:ic+2) .eq. 'r-z' ) then

                  imesh = 2

               else if( chlw(ic:ic+2) .eq. 'xyz' ) then

                  imesh = 3

               else
                  write(ErrCha,'("mesh card error")')
                  ErrID = 'L:15515/R:sours/F:read02.f' !E03_026_001
                  call ErrWrite(ErrID,ErrCha)
                  goto 998

               end if

               if( imesh .eq. 3 ) then

                  call txymesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                         jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                         ixtp,inx,xmin,xmax,xdel,istxg,
     &                         iytp,iny,ymin,ymax,ydel,istyg,
     &                         iztp,inz,zmin,zmax,zdel,istzg)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                      write(ErrCha,'("Error with mesh = xyz")')
                      ErrID = 'L:15532/R:sours/F:read02.f' !E03_026_002
                      call ErrWrite(ErrID,ErrCha)
                      goto 997
                     endif

                     if( ixtp .lt. 0 ) then
           write(ErrCha,'("Error with mesh = xyz. It is x mesh <0.")')
           ErrID = 'L:15539/R:sours/F:read02.f' !E03_026_003
           call ErrWrite(ErrID,ErrCha)
                      goto 997
                     endif

                     if( iytp .lt. 0 ) then
           write(ErrCha,'("Error with mesh = xyz. It is y mesh <0.")')
           ErrID = 'L:15546/R:sours/F:read02.f' !E03_026_004
           call ErrWrite(ErrID,ErrCha)
                      goto 997
                     endif

                     if( iztp .lt. 0 ) then
           write(ErrCha,'("Error with mesh = xyz. It is z mesh <0.")')
           ErrID = 'L:15553/R:sours/F:read02.f' !E03_026_005
           call ErrWrite(ErrID,ErrCha)
                      goto 997
                     endif

               else

cKN now, only mesh=xyz
                     ierrMSG=1 ! T.Sato 2017/07/13
                     write(ErrCha,'("mesh = xyz error")')
                     ErrID = 'L:15563/R:sours/F:read02.f' !E03_026_006
                     call ErrWrite(ErrID,ErrCha)
                     goto 997

               end if

                     jnx = iabs( inx )
                     jny = iabs( iny )
                     jnz = iabs( inz )

                  istdg = mmmax
                  inxyz = jnx * jny * jnz
                  istdc = istdg + 2 * inxyz + 1
                  mmmax = istdc + inxyz + 1

*-----------------------------------------------------------------------

                  ic = i1

                  ll = 0

               do ka = 1, jnz
               do ja = 1, jny
               do ia = 1, jnx

                  ii = ia
                  jj = ja
                  kk = ka
                  ll = ll + 1

                  if( inx .lt. 0 ) ii = jnx - ia + 1
                  if( iny .lt. 0 ) jj = jny - ja + 1
                  if( inz .lt. 0 ) kk = jnz - ka + 1

                  l = ii + jnx * ( jj - 1 )  + jnx * jny * ( kk - 1 )

                  if( ic .gt. i3 ) then

  154                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                      write(ErrCha,'("mesh = xyz error")')
                      ErrID = 'L:15607/R:sours/F:read02.f' !E03_026_007
                      call ErrWrite(ErrID,ErrCha)
                      goto 997
                     endif

                     if( iskip .ne. 0 ) goto 154

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

                  end if

                     ic = jnumc(chlw,ic,i3)

                     call snum(chlw,ic,i3,ic2,cvvv,ierr)

                     if( ierr .ne. 0 ) then
                      write(ErrCha,'("mesh = xyz error")')
                      ErrID = 'L:15628/R:sours/F:read02.f' !E03_026_008
                      call ErrWrite(ErrID,ErrCha)
                      goto 998
                     endif

                     das(istdg + l ) = cvvv
                     das(istdg + inxyz + ll ) = cvvv

                     ic = ic2

               end do
               end do
               end do

*-----------------------------------------------------------------------

                  goto 140

*-----------------------------------------------------------------------
*        s-type = 26 for surface source
*-----------------------------------------------------------------------

         else if( ipm .eq. 107 ) then

  160    continue

            if( ic .gt. i3 ) goto 140

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  incut = incut + 1

                  if( incut .gt. 8 ) goto 939

                  isvct(j,incut) = nint(cvvv)

                  ic = jnumc(chlw,ic2,i3)

                  goto 160

*-----------------------------------------------------------------------
*        set values
*-----------------------------------------------------------------------

         else

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) then
                write(ErrCha,'("set values error")')
                ErrID = 'L:15680/R:sours/F:read02.f' !E03_027_001
                call ErrWrite(ErrID,ErrCha)
                goto 998
               endif

               spava( ipm ) = cvvv

         end if

*-----------------------------------------------------------------------
*        next parameters
*-----------------------------------------------------------------------

         if( ipm .ne. 17 .and. ipm .ne. 63 .and. ipm .ne. 64 .and.
     &       ipm .ne. 70 .and. ipm .ne. 75 .and.
     &       ipm .ne. 79 .and. ipm .ne. 84 .and. ipm .ne. 94 ) then

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

               goto 140

*-----------------------------------------------------------------------
*        energy bin
*-----------------------------------------------------------------------

         else if( ipm .eq. 17 ) then

               if( chlw(icl+1:icl+1) .eq. ';' ) then
                write(ErrCha,'("energy bin error")')
                ErrID = 'L:15711/R:sours/F:read02.f' !E03_028_001
                call ErrWrite(ErrID,ErrCha)
                goto 994
               endif
               if( ischn(94) .ne. 0 ) then
        write(ErrCha,*)
     & 'Parameters ni is already specified in [source].'
        ErrID = 'L:15718/R:sours/F:read02.f' !E03_028_002
        call ErrWrite(ErrID,ErrCha)
                goto 941
               endif

               ngrp(j)  = abs( nint( spava(17) ) )

               ngll(j)  = 1
               if( spava(17) .lt. 0.d0 ) ngll(j) = -1

               if( nint(spava(18)) .eq. 28 .or.
     &             nint(spava(18)) .eq. 29 ) then
                 allocate(e_lowk(ngrp(j)),e_upk(ngrp(j)),
     &                    ratkari(ngrp(j)))

               else
                 call moddas_reallocate_dbl(
     &                   isrc, j, ngrp(j), ngei, egmin)
                 call moddas_reallocate_dbl(
     &                   isrc, j, ngrp(j), ngea, egmax)
                 call moddas_reallocate_dbl(
     &                   isrc, j, ngrp(j), ngfe, fegrp)
                 call moddas_reallocate_dbl(isrc, j, ngrp(j), ngft, rfe)
                 call moddas_reallocate_dbl(isrc, j, ngrp(j), ngpi, prw)
                 call moddas_reallocate_dbl(isrc, j, ngrp(j), ngpw, pwt)
               end if

               if( nint(spava(18)) .eq.  8 .or.
     &             nint(spava(18)) .eq. 18 .or.
     &             nint(spava(18)) .eq.  9 .or.
     &             nint(spava(18)) .eq. 19 ) then

                  nmva = 2 * ngrp(j)


               else if( nint(spava(18)) .eq. 22 .or.
     &                  nint(spava(18)) .eq. 32 .or.
     &                  nint(spava(18)) .eq. 23 .or.
     &                  nint(spava(18)) .eq. 33 .or.
     &                  nint(spava(18)) .eq. 28 .or.
     &                  nint(spava(18)) .eq. 29 ) then

                  nmva = 3 * ngrp(j)

c
               else

                  nmva = 2 * ngrp(j) + 1

               end if

  141                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                      write(ErrCha,'("energy bin error")')
                      ErrID = 'L:15775/R:sours/F:read02.f' !E03_028_003
                      call ErrWrite(ErrID,ErrCha)
                      goto 996
                     endif

                     if( iskip .ne. 0 ) goto 141

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

               do i = 1, nmva

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) then
                   m_err = 'Energy bin is wrong or insufficient'//
     &             ' in this or above line(s)'
                   ErrCha = ''
                   ErrID = 'L:15796/R:sours/F:read02.f'
                   goto 998
                  endif

                if( nint(spava(18)) .eq. 22 .or.
     &              nint(spava(18)) .eq. 32 .or.
     &              nint(spava(18)) .eq. 23 .or.
     &              nint(spava(18)) .eq. 33 .or.
     &              nint(spava(18)) .eq. 28 .or.
     &              nint(spava(18)) .eq. 29 ) then

                  if( i/3*3 .eq. i ) then

                    if( nint(spava(18)) .eq. 28 .or.
     &                  nint(spava(18)) .eq. 29 ) then
                      ratkari(i/3) = cvvv
                    else
                      fegrp(ngfe(j)+(i-1)/3+1) = cvvv
                    end if

                  else if( i/3*3+1 .eq. i ) then

                    if( nint(spava(18)) .eq. 28 .or.
     &                  nint(spava(18)) .eq. 29 ) then
                      e_lowk(1+i/3) = cvvv
                    else
                      egmin(ngei(j)+(i-1)/3+1) = cvvv
                    end if

                  else if( i/3*3+2 .eq. i ) then

                    if( nint(spava(18)) .eq. 28 .or.
     &                  nint(spava(18)) .eq. 29 ) then
                      e_upk(1+i/3) = cvvv
                    else
                      egmax(ngea(j)+(i-1)/3+1) = cvvv
                    end if

                  end if

                else

                  if( i .eq. 2 * ngrp(j) + 1 ) then

                      egmax(ngea(j)+ngrp(j)) = cvvv

                  else if( i/2*2 .ne. i ) then

                      egmin(ngei(j)+1+i/2) = cvvv

                  else if( i/2*2 .eq. i ) then

                      fegrp(ngfe(j)+i/2) = cvvv

                  end if

                end if

                  ic = ic2

                  if( i .lt. nmva .and. ic .gt. i3 ) then

  142                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                      write(ErrCha,'("energy bin error")')
                      ErrID = 'L:15864/R:sours/F:read02.f' !E03_028_005
                      call ErrWrite(ErrID,ErrCha)
                      goto 996
                     endif

                     if( iskip .ne. 0 ) goto 142

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

                  end if

               end do

               if( nint(spava(18)) .eq. 28 .or.
     &             nint(spava(18)) .eq. 29 ) then
               else
                 if( ngll(j) .lt. 0 .and.
     &               egmin(ngei(j)+1) .le. 0.0d0 ) then
           write(ErrCha,'("Error in energy bin. It is emin <0.")')
           ErrID = 'L:15887/R:sours/F:read02.f' !E03_028_006
           call ErrWrite(ErrID,ErrCha)
                  goto 996
                 endif

                 if( nint(spava(18)) .eq. 22 .or.
     &               nint(spava(18)) .eq. 32 .or.
     &               nint(spava(18)) .eq. 23 .or.
     &               nint(spava(18)) .eq. 33 ) then
                 else if( nint(spava(18)) .eq.  8 .or.
     &               nint(spava(18)) .eq. 18 .or.
     &               nint(spava(18)) .eq.  9 .or.
     &               nint(spava(18)) .eq. 19 ) then

                   do i = 1, ngrp(j)

                      egmax(ngea(j)+i) = egmin(ngei(j)+i)

                   end do

                 else

                   do i = 1, ngrp(j) - 1

                      egmax(ngea(j)+i) = egmin(ngei(j)+i+1)

                   end do

                 end if

               end if

*-----------------------------------------------------------------------

            if( nint(spava(18)) .eq. 21 .or.
     &          nint(spava(18)) .eq. 31 .or.
     &          nint(spava(18)) .eq. 24 .or.
     &          nint(spava(18)) .eq. 34 ) then

               do i = 1, ngrp(j)

                  fegrp(ngfe(j)+i) = fegrp(ngfe(j)+i)
     &            * ( egmax(ngea(j)+i) - egmin(ngei(j)+i) )

               end do

            end if

            if( nint(spava(18)) .eq. 28 .or.
     &          nint(spava(18)) .eq. 29 ) then
            else

              do i = 1, ngrp(j)

                 prw(ngpi(j)+i) = fegrp(ngfe(j)+i)

              end do

            end if

              goto 140

*-----------------------------------------------------------------------
*        angle bin
*-----------------------------------------------------------------------

         else if( ipm .eq. 70 ) then

               if( chlw(icl+1:icl+1) .eq. ';' ) then
                write(ErrCha,'("angle bin error")')
                ErrID = 'L:15957/R:sours/F:read02.f' !E03_029_001
                call ErrWrite(ErrID,ErrCha)
                goto 994
               endif

               narp(j)  = abs( nint( spava(70) ) )

               nall(j)  = 1
               if( spava(70) .lt. 0.d0 ) nall(j) = -1

               call moddas_reallocate_dbl(isrc, j, narp(j), naei, agmin)
               call moddas_reallocate_dbl(isrc, j, narp(j), naea, agmax)
               call moddas_reallocate_dbl(isrc, j, narp(j), nafe, fagrp)
               call moddas_reallocate_dbl(isrc, j, narp(j), naft, rfa)
               call moddas_reallocate_dbl(isrc, j, narp(j), napi, paw)
               call moddas_reallocate_dbl(isrc, j, narp(j), napw, pat)

               nmva = 2 * narp(j) + 1

  241                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                      write(ErrCha,'("angle bin error")')
                      ErrID = 'L:15982/R:sours/F:read02.f' !E03_029_002
                      call ErrWrite(ErrID,ErrCha)
                      goto 996
                     endif

                     if( iskip .ne. 0 ) goto 241

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

               do i = 1, nmva

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) then
                   m_err = 'Angle bin is wrong or insufficient'//
     &             ' in this or above line(s)'
                   ErrCha = ''
                   ErrID = 'L:16003/R:sours/F:read02.f'
                   goto 998
                  endif

                  if( i .eq. nmva ) then

                     agmax(naea(j)+narp(j)) = cvvv

                  else if( i/2*2 .ne. i ) then

                     agmin(naei(j)+1+i/2) = cvvv

                  else if( i/2*2 .eq. i ) then

                     fagrp(nafe(j)+i/2) = cvvv

                  end if

                  ic = ic2

                  if( i .lt. nmva .and. ic .gt. i3 ) then

  242                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                      write(ErrCha,'("angle bin error")')
                      ErrID = 'L:16031/R:sours/F:read02.f' !E03_029_004
                      call ErrWrite(ErrID,ErrCha)
                      goto 996
                     endif

                     if( iskip .ne. 0 ) goto 242

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

                  end if

               end do

               do i = 1, narp(j) - 1

                  agmax(naea(j)+i) = agmin(naei(j)+i+1)

               end do

               do i = 1, narp(j)

                     paw(napi(j)+i) = fagrp(nafe(j)+i)

               end do

               goto 140

*-----------------------------------------------------------------------
*        time bin
*-----------------------------------------------------------------------

         else if( ipm .eq. 79 ) then

               if( chlw(icl+1:icl+1) .eq. ';' ) then
                write(ErrCha,'("time bin error")')
                ErrID = 'L:16070/R:sours/F:read02.f' !E03_030_001
                call ErrWrite(ErrID,ErrCha)
                goto 994
               endif

               ntrp(j)  = abs( nint( spava(79) ) )

               ntll(j)  = 1
               if( spava(79) .lt. 0.d0 ) ntll(j) = -1

               call moddas_reallocate_dbl(isrc, j, ntrp(j), ntei, tgmin)
               call moddas_reallocate_dbl(isrc, j, ntrp(j), ntea, tgmax)
               call moddas_reallocate_dbl(isrc, j, ntrp(j), ntfe, ftgrp)
               call moddas_reallocate_dbl(isrc, j, ntrp(j), ntft, rft)
               call moddas_reallocate_dbl(isrc, j, ntrp(j), ntpi, ptw)
               call moddas_reallocate_dbl(isrc, j, ntrp(j), ntpw, ptt)

               nmva = 2 * ntrp(j) + 1

  251                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                      write(ErrCha,'("time bin error")')
                      ErrID = 'L:16095/R:sours/F:read02.f' !E03_030_002
                      call ErrWrite(ErrID,ErrCha)
                      goto 996
                     endif

                     if( iskip .ne. 0 ) goto 251

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

               do i = 1, nmva

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) then
                   m_err = 'Time bin is wrong or insufficient'//
     &             ' in this or above line(s)'
                   ErrCha = ''
                   ErrID = 'L:16116/R:sours/F:read02.f'
                   goto 998
                  endif

                  if( i .eq. nmva ) then

                     tgmax(ntea(j)+ntrp(j)) = cvvv

                  else if( i/2*2 .ne. i ) then

                     tgmin(ntei(j)+1+i/2) = cvvv

                  else if( i/2*2 .eq. i ) then

                     ftgrp(ntfe(j)+i/2) = cvvv

                  end if

                  ic = ic2

                  if( i .lt. nmva .and. ic .gt. i3 ) then

  252                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                      write(ErrCha,'("time bin error")')
                      ErrID = 'L:16144/R:sours/F:read02.f' !E03_030_004
                      call ErrWrite(ErrID,ErrCha)
                      goto 996
                     endif

                     if( iskip .ne. 0 ) goto 252

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

                  end if

               end do

               do i = 1, ntrp(j) - 1

                  tgmax(ntea(j)+i) = tgmin(ntei(j)+i+1)

               end do

               do i = 1, ntrp(j)

                     ptw(ntpi(j)+i) = ftgrp(ntfe(j)+i)

               end do

               goto 140

*-----------------------------------------------------------------------
*        energy probability bin
*-----------------------------------------------------------------------

         else if( ipm .eq. 63 ) then

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 984

               jptyp(j)  = abs( nint( spava(63) ) )

               if( ischn(17) .eq. 0 .and. ischn(68) .eq. 0 ) then
                write(ErrCha,'("energy probability bin error")')
                ErrID = 'L:16187/R:sours/F:read02.f' !E03_031_001
                call ErrWrite(ErrID,ErrCha)
                goto 972
               endif

               if( ischn(18) .eq. 0 ) goto 996

         if( nint(spava(18)) .eq.  1 .or. nint(spava(18)) .eq. 11 .or.
     &       nint(spava(18)) .eq. 21 .or. nint(spava(18)) .eq. 31 .or.
     &       nint(spava(18)) .eq.  8 .or. nint(spava(18)) .eq. 18 .or.
     &       nint(spava(18)) .eq.  5 .or. nint(spava(18)) .eq. 15 .or.
     &       nint(spava(18)) .eq.  3 ) goto 961

            if( nint(spava(18)) .eq. 29 ) allocate(prbkari(ngrp(j)))

            if( jptyp(j) .eq. 1 ) then

  161                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                      write(ErrCha,'("energy probability bin error")')
                      ErrID = 'L:16210/R:sours/F:read02.f' !E03_031_002
                      call ErrWrite(ErrID,ErrCha)
                      goto 996
                     endif

                     if( iskip .ne. 0 ) goto 161

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

               do i = 1, ngrp(j)

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) then
                   m_err = 'Energy probability is wrong'//
     &             ' in this or above line(s)'
                   ErrCha = ''
                   ErrID = 'L:16231/R:sours/F:read02.f'
                   goto 998
                  endif

                  if( nint(spava(18)) .eq. 29 ) then
                     prbkari(i) = cvvv
                  else
                     prw(ngpi(j)+i) = cvvv
                  end if

                  ic = ic2

                  if( i .lt. ngrp(j) .and. ic .gt. i3 ) then

  162                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                      write(ErrCha,'("energy probability bin error")')
                      ErrID = 'L:16251/R:sours/F:read02.f' !E03_031_004
                      call ErrWrite(ErrID,ErrCha)
                      goto 996
                     endif

                     if( iskip .ne. 0 ) goto 162

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

                  end if

               end do

            else if( jptyp(j) .eq. 0 ) then

               do i = 1, ngrp(j)

                     prw(ngpi(j)+i) = 1.0d0
                     prbkari(i) = 1.0d0

               end do

            end if

               goto 140

*-----------------------------------------------------------------------
*        angle probability bin
*-----------------------------------------------------------------------

         else if( ipm .eq. 75 ) then

               if( chlw(icl+1:icl+1) .eq. ';' ) then
                write(ErrCha,'("angle probability bin error")')
                ErrID = 'L:16289/R:sours/F:read02.f' !E03_032_001
                call ErrWrite(ErrID,ErrCha)
                goto 984
               endif

               jqtyp(j)  = abs( nint( spava(75) ) )

               if( ischn(70) .eq. 0 .and. ischn(73) .eq. 0 ) goto 972
               if( ischn(69) .eq. 0 ) then
                write(ErrCha,'("angle probability bin error")')
                ErrID = 'L:16299/R:sours/F:read02.f' !E03_032_002
                call ErrWrite(ErrID,ErrCha)
                goto 996
               endif

            if( jqtyp(j) .eq. 1 ) then

  261                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                     write(ErrCha,'("angle probability bin error")')
                     ErrID = 'L:16312/R:sours/F:read02.f' !E03_032_003
                     call ErrWrite(ErrID,ErrCha)
                     goto 996
                     endif

                     if( iskip .ne. 0 ) goto 261

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

               do i = 1, narp(j)

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) then
                   m_err = 'Angle probability is wrong'//
     &             ' in this or above line(s)'
                   ErrCha = ''
                   ErrID = 'L:16333/R:sours/F:read02.f'
                   goto 998
                  endif

                     paw(napi(j)+i) = cvvv

                  ic = ic2

                  if( i .lt. narp(j) .and. ic .gt. i3 ) then

  262                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                     write(ErrCha,'("angle probability bin error")')
                     ErrID = 'L:16349/R:sours/F:read02.f' !E03_032_005
                     call ErrWrite(ErrID,ErrCha)
                     goto 996
                     endif

                     if( iskip .ne. 0 ) goto 262

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

                  end if

               end do

            else if( jqtyp(j) .eq. 0 ) then

               do i = 1, narp(j)

                     paw(napi(j)+i) = 1.0d0

               end do

            end if

               goto 140

*-----------------------------------------------------------------------
*        time probability bin
*-----------------------------------------------------------------------

         else if( ipm .eq. 84 ) then

               if( chlw(icl+1:icl+1) .eq. ';' ) then
                write(ErrCha,'("time probability bin error")')
                ErrID = 'L:16386/R:sours/F:read02.f' !E03_033_001
                call ErrWrite(ErrID,ErrCha)
                goto 984
               endif

               jotyp(j)  = abs( nint( spava(84) ) )

               if( ischn(79) .eq. 0 .and. ischn(82) .eq. 0 ) then
                write(ErrCha,'("time probability bin error")')
                ErrID = 'L:16395/R:sours/F:read02.f' !E03_033_002
                call ErrWrite(ErrID,ErrCha)
                goto 872
               endif
               if( ischn(29) .eq. 0 ) then
                write(ErrCha,'("time probability bin error")')
                ErrID = 'L:16401/R:sours/F:read02.f' !E03_033_003
                call ErrWrite(ErrID,ErrCha)
                goto 996
               endif

            if( jotyp(j) .eq. 1 ) then

  271                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                     write(ErrCha,'("time probability bin error")')
                     ErrID = 'L:16414/R:sours/F:read02.f' !E03_033_004
                     call ErrWrite(ErrID,ErrCha)
                     goto 996
                     endif

                     if( iskip .ne. 0 ) goto 271

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

               do i = 1, ntrp(j)

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) then
                   m_err = 'Time probability is wrong'//
     &             ' in this or above line(s)'
                   ErrCha = ''
                   ErrID = 'L:16435/R:sours/F:read02.f'
                   goto 998
                  endif

                     ptw(napi(j)+i) = cvvv

                  ic = ic2

                  if( i .lt. narp(j) .and. ic .gt. i3 ) then

  272                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                     write(ErrCha,'("time probability bin error")')
                     ErrID = 'L:16451/R:sours/F:read02.f' !E03_033_006
                     call ErrWrite(ErrID,ErrCha)
                     goto 996
                     endif

                     if( iskip .ne. 0 ) goto 272

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

                  end if

               end do

            else if( jotyp(j) .eq. 0 ) then

               do i = 1, ntrp(j)

                     ptw(ntpi(j)+i) = 1.0d0

               end do

            end if

               goto 140

*-----------------------------------------------------------------------
*        duct length probability bin
*-----------------------------------------------------------------------

         else if( ipm .eq. 64 ) then

               if( chlw(icl+1:icl+1) .eq. ';' ) then
                write(ErrCha,'("duct length probability bin error")')
                ErrID = 'L:16488/R:sours/F:read02.f' !E03_034_001
                call ErrWrite(ErrID,ErrCha)
                goto 971
               endif

               nglp(j)  = nint( spava(64) )

               if( nglp(j) < 0 ) then
                  write(ErrCha,'(a,i5,a,i5)')
     &                 'duct length probability bin: nglp(',j
     &                 ,')=',nglp(j)
                  ErrID = 'L:16499/R:sours/F:read02.f'
                  call ErrWrite(ErrID,ErrCha)
               endif
               call moddas_reallocate_dbl(isrc, j, nglp(j), ngli, slmin)
               call moddas_reallocate_dbl(isrc, j, nglp(j), ngla, slmax)
               call moddas_reallocate_dbl(isrc, j, nglp(j), ngfl, flgrp)
               call moddas_reallocate_dbl(isrc, j, nglp(j), nglc, rgl)
               call moddas_reallocate_dbl(isrc, j, nglp(j), nglw, rgw)

               nmva = 2 * nglp(j) + 1

  171                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                if( ierr .ne. 0 ) return
                if( jpn  .eq. 3 ) then
                write(ErrCha,'("duct length probability bin error")')
                ErrID = 'L:16516/R:sours/F:read02.f' !E03_034_002
                call ErrWrite(ErrID,ErrCha)
                goto 996
                endif

                     if( iskip .ne. 0 ) goto 171

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

               do i = 1, nmva

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                if( ierr .ne. 0 ) then
                 m_err = 'Duct length probability is wrong'//
     &           ' in this or above line(s)'
                 ErrCha = ''
                 ErrID = 'L:16537/R:sours/F:read02.f'
                 goto 998
                endif

                  if( i .eq. nmva ) then

                     slmax(ngla(j)+nglp(j)) = cvvv

                  else if( i/2*2 .ne. i ) then

                     slmin(ngli(j)+1+i/2) = cvvv

                  else if( i/2*2 .eq. i ) then

                     flgrp(ngfl(j)+i/2) = cvvv

                  end if

                  ic = ic2

                  if( i .lt. nmva .and. ic .gt. i3 ) then

  172                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) then
                  write(ErrCha,'("duct length probability bin error")')
                  ErrID = 'L:16565/R:sours/F:read02.f' !E03_034_004
                  call ErrWrite(ErrID,ErrCha)
                  goto 996
                  endif

                     if( iskip .ne. 0 ) goto 172

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

                  end if

               end do

               do i = 1, nglp(j) - 1

                  slmax(ngla(j)+i) = slmin(ngli(j)+i+1)

               end do

               goto 140

*-----------------------------------------------------------------------
*        nuclide and intensity bin      by MATSUDA 2016.07.31
*-----------------------------------------------------------------------

         else if( ipm .eq. 94 ) then

               if( chlw(icl+1:icl+1) .eq. ';' ) then
                write(ErrCha,'("nuclide and intensity bin error")')
                ErrID = 'L:16598/R:sours/F:read02.f' !E03_035_001
                call ErrWrite(ErrID,ErrCha)
                goto 942
               endif
               if( ischn(17) .ne. 0 ) then
                write(ErrCha,'("nuclide and intensity bin error")')
                ErrID = 'L:16604/R:sours/F:read02.f' !E03_035_002
                call ErrWrite(ErrID,ErrCha)
               goto 943
               endif

               niorg(j)  = abs( nint( spava(94) ) )

  341                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                write(ErrCha,'("nuclide and intensity bin error")')
                ErrID = 'L:16617/R:sours/F:read02.f' !E03_035_003
                call ErrWrite(ErrID,ErrCha)
                      goto 944
                     endif

                     if( iskip .ne. 0 ) goto 341

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1

               do i = 1, niorg(j)
                     icl = inumc(chlw, ic, i3, ' ') - 1
                     if( icl .gt. i3 ) then
                write(ErrCha,'("nuclide and intensity bin error")')
                ErrID = 'L:16634/R:sours/F:read02.f' !E03_035_004
                call ErrWrite(ErrID,ErrCha)
                      goto 997
                     endif

                     if( (icl-ic+1) .gt. 7 ) then
                write(ErrCha,'("nuclide and intensity bin error")')
                ErrID = 'L:16641/R:sours/F:read02.f' !E03_035_005
                call ErrWrite(ErrID,ErrCha)
                      goto 944
                     endif

                 chaa = '       '
                 chaa = chlw(ic:icl)
                 if( (icl-ic+1) .eq. 7 ) then
                 else
                   do ii = (icl-ic+2), 7
                     chaa(ii:ii) = ' '
                   end do
                 end if

                 call getelm(2,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*       you can choose INPUT or OUTPUT style
*         INPUT:
*           2: chaa(ANY)   -> chza,8digit(1001.0(Zaid+))
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )
*          12: chza(Zaid+) -> chaa,7digit(H1     )
*          13: chza(Zaid+) -> chaa,7digit(1H     )
*          14: chza(Zaid+) -> chaa,7digit(%H   1 )

                 if( ierr .eq. 1 ) then
                   l_err = ill(jsn)
                   k_err = jsn

! T.Sato 2019/01/29 output error ID when something wrong in reading RI name
       write(ErrCha,*) 'Error in reading RI name in [source]'
       ErrID = 'L:16673/R:sours/F:read02.f' !E03_109_001
       call ErrWrite(ErrID,ErrCha)

                     return
                 end if

                 chalct(i,j) = chza

                 icl = icl + 1
                 ic = jnumc(chlw, icl, i3)
                 if( i .le. niorg(j) .and. ic .gt. i3 ) then

  342              call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                   if( ierr .ne. 0 ) return
                   if( jpn  .eq. 3 ) then
                write(ErrCha,'("nuclide and intensity bin error")')
                ErrID = 'L:16691/R:sours/F:read02.f' !E03_035_006
                call ErrWrite(ErrID,ErrCha)
                    goto 944
                   endif

                   if( iskip .ne. 0 ) goto 342

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                     ic = i1
                   call snum(chlw,ic,i3,ic2,cvvv,ierr)

                   if( ierr .ne. 0 ) then
                write(ErrCha,'("nuclide and intensity bin error")')
                ErrID = 'L:16707/R:sours/F:read02.f' !E03_035_007
                call ErrWrite(ErrID,ErrCha)
                    goto 944
                   endif

                 else

                   call snum(chlw,ic,i3,ic2,cvvv,ierr)

                   if( ierr .ne. 0 ) then
                write(ErrCha,'("nuclide and intensity bin error")')
                ErrID = 'L:16718/R:sours/F:read02.f' !E03_035_008
                call ErrWrite(ErrID,ErrCha)
                    goto 944
                   endif

                 end if

                   act000(i,j) = cvvv
                 if( i .eq. niorg(j) ) then
                 else

                     ic = ic2

                   if( i .lt. niorg(j) .and. ic .gt. i3 ) then

  343                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) then
                write(ErrCha,'("nuclide and intensity bin error")')
                ErrID = 'L:16739/R:sours/F:read02.f' !E03_035_009
                call ErrWrite(ErrID,ErrCha)
                      goto 944
                     endif

                     if( iskip .ne. 0 ) goto 343

                     if( iprojall .eq. -1 ) then
                        write(niws,'(a)') chin(1:i2)      ! work souce section
                     endif

                       ic = i1

                   else
                     ic = jnumc(chlw, ic, i3)

                   end if
                 end if
               end do

               goto 140

         end if

*-----------------------------------------------------------------------
*     summary and check
*-----------------------------------------------------------------------

  800 continue

         imrd = imrd + 1

*-----------------------------------------------------------------------
*        compativility check for event dump source
* FURUTA20150515
*-----------------------------------------------------------------------

         if(idmpsors.eq.0)idmpmode=0
         if(idmpmode.eq.1)then
          if(imsrc .gt. 1) then
          write(ErrCha,*)
     & 'compativility check for event dump source error'
          ErrID = 'L:16781/R:sours/F:read02.f' !E03_036_001
          call ErrWrite(ErrID,ErrCha)
          goto 958
          endif
          if(irestart .eq. 1) goto 959
          if(jsdmp(j,18).le.0.or.jsdmp(j,19).le.0)then
           write(ErrCha,'(''Either NOCAS or NOBCH'',
     &             '' does not exist in dump'')')
           ErrID = 'L:16789/R:sours/F:read02.f' !E03_036_002
           call ErrWrite(ErrID,ErrCha)
           write(ErrCha,'(''idmpmode changed to 0'')')
           ErrID = 'L:16792/R:sours/F:read02.f'
           call ErrWrite(ErrID,ErrCha)
           idmpmode=0
          endif
          if(idmpjmp(2).eq.0)then
           dmpmulti=1.0d0
          else
           if(dmpmulti.eq.0.0) then
          write(ErrCha,*)
     & 'compativility check for event dump source error'
          ErrID = 'L:16802/R:sours/F:read02.f' !E03_036_003
          call ErrWrite(ErrID,ErrCha)
            goto 960
           endif
          endif
         else
          if(idmpjmp(2).eq.0)dmpmulti=0.0d0
         endif

*-----------------------------------------------------------------------

         if( ischn(1) .eq. 0 ) then
          ierrMSG=2 ! T.Sato 2017/07/13
          goto 997
         endif
               jstyp(j) = nint( spava(1) )
               jstypori(j) = nint( spava(1) ) ! T.Sato, original jstyp written in input file

      if(ischn(10).eq.1.and.ischn(18).eq.1) then ! warning for both e0 and e-type are specified
       if(jstyp(j).eq.4.or.jstyp(j).eq.5.or.jstyp(j).eq.6.or.
     & jstyp(j).eq.8.or.jstyp(j).eq.10.or.jstyp(j).eq.12.or.
     & jstyp(j).eq.14.or.jstyp(j).eq.16.or.jstyp(j).eq.19.or.
     & jstyp(j).eq.21.or.jstyp(j).eq.23.or.jstyp(j).eq.25) then
        write(ErrCha,'("Warning: Both e0 and e-type are specified,",
     &  " but only e-type is used in the simulation")')
       else
        write(ErrCha,'("Warning: both e0 and e-type are specified,",
     &  " but only e0 is used in the simulation")')
       endif
       ErrID = 'L:16831/R:sours/F:read02.f' !E03_037_001
       call ErrWrite(ErrID,ErrCha)
      endif
! Check mono-energy and energy distribution, T.Sato 2017/07/13
      if(jstyp(j).eq.1.and.ischn(10).eq.0) then ! mono-energy type, but no e0
       if(ischn(18).ne.0) then
        jstyp(j)=4
       else
        ierrMSG=5
        write(ErrCha,'("mono-energy type, but no e0")')
        ErrID = 'L:16841/R:sours/F:read02.f' !E03_037_001
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.4.and.ischn(18).eq.0) then ! energy distribution type, but no e-type
       if(ischn(10).ne.0) then
        jstyp(j)=1
       else
        ierrMSG=5
        write(ErrCha,'("energy distribution type, but no e-type")')
        ErrID = 'L:16851/R:sours/F:read02.f' !E03_038_001
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.2.and.ischn(10).eq.0) then ! mono-energy type, but no e0
       if(ischn(18).ne.0) then
        jstyp(j)=5
       else
        ierrMSG=5
        write(ErrCha,'("mono-energy type, but no e0")')
        ErrID = 'L:16861/R:sours/F:read02.f' !E03_037_002
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.5.and.ischn(18).eq.0) then ! energy distribution type, but no e-type
       if(ischn(10).ne.0) then
        jstyp(j)=2
       else
        ierrMSG=5
        write(ErrCha,'("energy distribution type, but no e-type")')
        ErrID = 'L:16871/R:sours/F:read02.f' !E03_038_002
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.3.and.ischn(10).eq.0) then ! mono-energy type, but no e0
       if(ischn(18).ne.0) then
        jstyp(j)=6
       else
        ierrMSG=5
        write(ErrCha,'("mono-energy type, but no e0")')
        ErrID = 'L:16881/R:sours/F:read02.f' !E03_037_003
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.6.and.ischn(18).eq.0) then ! energy distribution type, but no e-type
       if(ischn(10).ne.0) then
        jstyp(j)=3
       else
        ierrMSG=5
        write(ErrCha,'("energy distribution type, but no e-type")')
        ErrID = 'L:16891/R:sours/F:read02.f' !E03_038_003
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.7.and.ischn(10).eq.0) then ! mono-energy type, but no e0
       if(ischn(18).ne.0) then
        jstyp(j)=8
       else
        ierrMSG=5
        write(ErrCha,'("mono-energy type, but no e0")')
        ErrID = 'L:16901/R:sours/F:read02.f' !E03_037_004
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.8.and.ischn(18).eq.0) then ! energy distribution type, but no e-type
       if(ischn(10).ne.0) then
        jstyp(j)=7
       else
        ierrMSG=5
        write(ErrCha,'("energy distribution type, but no e-type")')
        ErrID = 'L:16911/R:sours/F:read02.f' !E03_038_004
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.9.and.ischn(10).eq.0) then ! mono-energy type, but no e0
       if(ischn(18).ne.0) then
        jstyp(j)=10
       else
        ierrMSG=5
        write(ErrCha,'("mono-energy type, but no e0")')
        ErrID = 'L:16921/R:sours/F:read02.f' !E03_037_005
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.10.and.ischn(18).eq.0) then ! energy distribution type, but no e-type
       if(ischn(10).ne.0) then
        jstyp(j)=9
       else
        ierrMSG=5
        write(ErrCha,'("energy distribution type, but no e-type")')
        ErrID = 'L:16931/R:sours/F:read02.f' !E03_038_005
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.11.and.ischn(10).eq.0) then ! mono-energy type, but no e0, T.Sato 2020/11/29
       if(ischn(18).ne.0) then
        jstyp(j)=-11  ! T.Sato 2020/11/29, s-type=-11 is used for s-type=11 with e-type
       else
        ierrMSG=5
        write(ErrCha,'("mono-energy type, but no e0")')
        ErrID = 'L:16941/R:sours/F:read02.f' !E03_037_006
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.13.and.ischn(10).eq.0) then ! mono-energy type, but no e0
       if(ischn(18).ne.0) then
        jstyp(j)=14
       else
        ierrMSG=5
        write(ErrCha,'("mono-energy type, but no e0")')
        ErrID = 'L:16951/R:sours/F:read02.f' !E03_037_006
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.14.and.ischn(18).eq.0) then ! energy distribution type, but no e-type
       if(ischn(10).ne.0) then
        jstyp(j)=13
       else
        ierrMSG=5
        write(ErrCha,'("energy distribution type, but no e-type")')
        ErrID = 'L:16961/R:sours/F:read02.f' !E03_038_006
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.15.and.ischn(10).eq.0) then ! mono-energy type, but no e0
       if(ischn(18).ne.0) then
        jstyp(j)=16
       else
        ierrMSG=5
        write(ErrCha,'("mono-energy type, but no e0")')
        ErrID = 'L:16971/R:sours/F:read02.f' !E03_037_007
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.16.and.ischn(18).eq.0) then ! energy distribution type, but no e-type
       if(ischn(10).ne.0) then
        jstyp(j)=15
       else
        ierrMSG=5
        write(ErrCha,'("energy distribution type, but no e-type")')
        ErrID = 'L:16981/R:sours/F:read02.f' !E03_038_007
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.18.and.ischn(10).eq.0) then ! mono-energy type, but no e0
       if(ischn(18).ne.0) then
        jstyp(j)=19
       else
        ierrMSG=5
        write(ErrCha,'("mono-energy type, but no e0")')
        ErrID = 'L:16991/R:sours/F:read02.f' !E03_037_008
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.19.and.ischn(18).eq.0) then ! energy distribution type, but no e-type
       if(ischn(10).ne.0) then
        jstyp(j)=18
       else
        ierrMSG=5
        write(ErrCha,'("energy distribution type, but no e-type")')
        ErrID = 'L:17001/R:sours/F:read02.f' !E03_038_008
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.20.and.ischn(10).eq.0) then ! mono-energy type, but no e0
       if(ischn(18).ne.0) then
        jstyp(j)=21
       else
        ierrMSG=5
        write(ErrCha,'("mono-energy type, but no e0")')
        ErrID = 'L:17011/R:sours/F:read02.f' !E03_037_009
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.21.and.ischn(18).eq.0) then ! energy distribution type, but no e-type
       if(ischn(10).ne.0) then
        jstyp(j)=20
       else
        ierrMSG=5
        write(ErrCha,'("energy distribution type, but no e-type")')
        ErrID = 'L:17021/R:sours/F:read02.f' !E03_038_009
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.22.and.ischn(10).eq.0) then ! mono-energy type, but no e0
       if(ischn(18).ne.0) then
        jstyp(j)=23
       else
        ierrMSG=5
        write(ErrCha,'("mono-energy type, but no e0")')
        ErrID = 'L:17031/R:sours/F:read02.f' !E03_037_010
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.23.and.ischn(18).eq.0) then ! energy distribution type, but no e-type
       if(ischn(10).ne.0) then
        jstyp(j)=22
       else
        ierrMSG=5
        write(ErrCha,'("energy distribution type, but no e-type")')
        ErrID = 'L:17041/R:sours/F:read02.f' !E03_038_010
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.24.and.ischn(10).eq.0) then ! mono-energy type, but no e0
       if(ischn(18).ne.0) then
        jstyp(j)=25
       else
        ierrMSG=5
        write(ErrCha,'("mono-energy type, but no e0")')
        ErrID = 'L:17051/R:sours/F:read02.f' !E03_037_011
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      elseif(jstyp(j).eq.25.and.ischn(18).eq.0) then ! energy distribution type, but no e-type
       if(ischn(10).ne.0) then
        jstyp(j)=24
       else
        ierrMSG=5
        write(ErrCha,'("energy distribution type, but no e-type")')
        ErrID = 'L:17061/R:sours/F:read02.f' !E03_038_011
        call ErrWrite(ErrID,ErrCha)
        goto 997
       endif
      endif

         if( ischn(2) .eq. 0 .and. jstyp(j) .le.  16 ) then
          ierrMSG=3 ! T.Sato 2017/07/13
          write(ErrCha,'("mono-energy and energy distribution")')
          ErrID = 'L:17070/R:sours/F:read02.f' !E03_039_001
          call ErrWrite(ErrID,ErrCha)
          goto 997
         endif
         if( ischn(2) .eq. 0 .and.
     &       (jstyp(j) .ge.  18 .and. jstyp(j) .ne. 100) ) then
          ierrMSG=3 ! T.Sato 2017/07/13
          write(ErrCha,'("mono-energy and energy distribution")')
          ErrID = 'L:17078/R:sours/F:read02.f' !E03_039_002
          call ErrWrite(ErrID,ErrCha)
          goto 997
         endif

*-----------------------------------------------------------------------
*        global factor
*-----------------------------------------------------------------------

               sfactor(j) = spava(43)

*-----------------------------------------------------------------------
*        charge state
*-----------------------------------------------------------------------

               lstyp(j) = nint( spava(78) )

*-----------------------------------------------------------------------
*        spin
*-----------------------------------------------------------------------

               ssx(j) = spava(44)
               ssy(j) = spava(45)
               ssz(j) = spava(46)

*-----------------------------------------------------------------------
*        transformation
*-----------------------------------------------------------------------

               isort(j,1) = igkst
               isort(j,2) = ktrs
               isort(j,3) = idtt
               isort(j,4) = 0

*-----------------------------------------------------------------------
*        region parameter
*-----------------------------------------------------------------------

         if( nsrn(j) .gt. 0 ) then

                  call moddas_reallocate_int(
     &                    isrc, j, MAX_NUM_NSRC, nsrc, idas_nsrc)
                  idsm = nsrc(j)
                  jdsm = 0

                  jdsm = jdsm + 1
                  idas_nsrc(idsm+jdsm) = nsrn(j)

                  jdsm = jdsm + 1
                  idas_nsrc(idsm+jdsm) = msrn

                  do i = 1, msrn

                     jdsm = jdsm + 1
                     idas_nsrc(idsm+jdsm)
     &                  = idas_nsrc_temporary(ndsm+i-1)

                  end do

                  if( jdsm > MAX_NUM_NSRC ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.sours@read02.f'
     &                       //' ?dimension over idas_nsrc?'
     &                       //' jdsm > MAX_NUM_NSRC'
     &                    ,' (jdsm=',jdsm,')'
     &                    ,' (MAX_NUM_NSRC@moddas.f=',MAX_NUM_NSRC,')'
                     ErrID = 'L:17144/R:sours/F:read02.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                  call moddas_reduce_int(
     &                    isrc, j, jdsm+1, nsrc, idas_nsrc)
                  call moddas_deallocate_int(idas_nsrc_temporary)

         end if

*-----------------------------------------------------------------------
*        the other parameters
*-----------------------------------------------------------------------

               nsmx(j) = nint( spava(38) )
               swt0(j) = spava(39)
              ispfs(j) = nint( spava(76) )

*-----------------------------------------------------------------------

         if( jstyp(j) .eq. 1 ) then

               sx0(j) = 0.d0
               sy0(j) = 0.d0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)

               sz0(j)    = spava( 7)
               sz1(j)    = spava( 8)
               sr0(j)    = spava( 9)
               se0(j)    = spava(10)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)
               sr1(j)    = spava(35)
               if( sr1(j) .gt. sr0(j) ) then
                write(ErrCha,*)
     & 'inner radius r1 is greater than r0 in [source]'
                ErrID = 'L:17183/R:sours/F:read02.f' !E03_041_001
                call ErrWrite(ErrID,ErrCha)
                goto 968
               endif

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 4 ) then


               sx0(j) = 0.d0
               sy0(j) = 0.d0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)

               sz0(j)    = spava( 7)
               sz1(j)    = spava( 8)
               sr0(j)    = spava( 9)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)
               sr1(j)    = spava(35)
               if( sr1(j) .gt. sr0(j) ) then
                write(ErrCha,*)
     & 'At s-type = 4. The radius r1 of [source] is larger than r0.'
                ErrID = 'L:17209/R:sours/F:read02.f' !E03_042_001
                call ErrWrite(ErrID,ErrCha)
                goto 968
               endif

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 2 ) then


               sx0(j)    = spava( 3)
               sx1(j)    = spava( 4)
               sy0(j)    = spava( 5)
               sy1(j)    = spava( 6)
               sz0(j)    = spava( 7)
               sz1(j)    = spava( 8)
               se0(j)    = spava(10)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 5 ) then


               sx0(j)    = spava( 3)
               sx1(j)    = spava( 4)
               sy0(j)    = spava( 5)
               sy1(j)    = spava( 6)
               sz0(j)    = spava( 7)
               sz1(j)    = spava( 8)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 3 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) sz0(j) = spava(7)

               sx1(j)    = spava( 4)
               sy1(j)    = spava( 6)
               sz1(j)    = spava( 8)
!  T.Sato 2024/07/14, add warning
         if(sx1(j).lt.0.0) then
          sx1(j) =abs(sx1(j))
          write(ErrCha,'("Warning: x1 in [source] is changed to",es12.4,
     &    " because it should be positive for s-type = 3")') sx1(j)
          ErrID = 'L:17266/R:sours/F:read02.f' !E03_042_001
          call ErrWrite(ErrID,ErrCha)
         endif
         if(sy1(j).lt.0.0) then
          sy1(j) =abs(sy1(j))
          write(ErrCha,'("Warning: y1 in [source] is changed to",es12.4,
     &    " because it should be positive for s-type = 3")') sy1(j)
          ErrID = 'L:17273/R:sours/F:read02.f' !E03_042_001
          call ErrWrite(ErrID,ErrCha)
         endif
         if(sz1(j).lt.0.0) then
          sz1(j) =abs(sz1(j))
          write(ErrCha,'("Warning: z1 in [source] is changed to",es12.4,
     &    " because it should be positive for s-type = 3")') sz1(j)
          ErrID = 'L:17280/R:sours/F:read02.f' !E03_042_001
          call ErrWrite(ErrID,ErrCha)
         endif
               se0(j)    = spava(10)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 6 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) sz0(j) = spava(7)

               sx1(j)    = spava( 4)
               sy1(j)    = spava( 6)
               sz1(j)    = spava( 8)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 7 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) sz0(j) = spava(7)

               sx1(j)    = spava( 4)
               sy1(j)    = spava( 6)
               sz1(j)    = spava( 8)
               se0(j)    = spava(10)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

               isrn(j) = 2

               if( ischn(15) .ne. 0 ) isrn(j)   = nint( spava(15) )
               if( isrn(j) .lt. 2 .or. isrn(j)/2*2 .ne. isrn(j) ) then
                write(ErrCha,'("s-type = 7 error")')
                ErrID = 'L:17334/R:sours/F:read02.f' !E03_043_001
                call ErrWrite(ErrID,ErrCha)
                goto 986
               endif

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 8 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) sz0(j) = spava(7)

               sx1(j)    = spava( 4)
               sy1(j)    = spava( 6)
               sz1(j)    = spava( 8)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

               isrn(j) = 2

               if( ischn(15) .ne. 0 ) isrn(j)   = nint( spava(15) )
               if( isrn(j) .lt. 2 .or. isrn(j)/2*2 .ne. isrn(j) ) then
                write(ErrCha,'("s-type = 8 error")')
                ErrID = 'L:17364/R:sours/F:read02.f' !E03_044_001
                call ErrWrite(ErrID,ErrCha)
                goto 986
               endif

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 9 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) sz0(j) = spava(7)

               se0(j)    = spava(10)
               sdir(j)   = spava(11)
               sr1(j)    = spava(35)
               sr2(j)    = spava(36)

            if(sdir(j).lt.-3.5) sdir(j)=-3.0  ! for s-type = 9, dir = cr and iso are equivalent

            if(sdir(j).gt.-3.5.and.sdir(j).lt.-2.5) then ! for 'iso' source, T.Sato 2017/11/21
             if( sr2(j) .gt. sr1(j) ) then
              write(ErrCha,'("r1 < r2 is not allowed for dir = iso")')
              ErrID = 'L:17392/R:sours/F:read02.f' !E03_045_001
              call ErrWrite(ErrID,ErrCha)
              goto 885
             endif
             ! T.Sato 2020/12/13, limitation of zenith angle
             if(ischn(71).eq.0) then
              sag1(j)=-1.0d0
             else
              sag1(j)=spava(71)
             endif
             if(ischn(72).eq.0) then
              sag2(j)=1.0d0
             else
              sag2(j)=spava(72)
             endif
             if(sag1(j).ge.sag2(j)) then
              goto 925
             endif
             ! T.Sato 2021/02/26, limitation of azimuth angle
             if(ischn(126).eq.0) then
              spg1(j)=0.0d0
             else
              spg1(j)=spava(126)
             endif
             if(ischn(127).eq.0) then
              spg2(j)=360.0d0
             else
              spg2(j)=spava(127)
             endif
! T.Sato 2024/09/07
             if(spg1(j).lt.-360d0.or.spg1(j).gt.360d0) goto 923
             if(spg2(j).lt.-360d0.or.spg2(j).gt.360d0) goto 924
             if(spg1(j).lt.0.0) spg1(j)=spg1(j)+360.0d0
             if(spg2(j).lt.0.0) spg2(j)=spg2(j)+360.0d0
            else
             if( sr1(j) .gt. sr2(j) ) then
              write(ErrCha,'("s-type = 9 error")')
              ErrID = 'L:17429/R:sours/F:read02.f' !E03_045_002
              call ErrWrite(ErrID,ErrCha)
              goto 985
             endif
            endif

cFURUTA20191129 bugfix
            if(sr2(j).eq.0.0)then
             write(6,*)'Warning: r2=0 is not allowed in s-type=9'
             write(6,*)'         r2=1.0E-10 is used instead'
             sr2(j)=1.0d-10
            endif

            if( sdir(j) .eq. 200.0 .or. sdir(j) .eq. 0.0 ) then

               sdir(j) =  200.0

            else if( sdir(j) .gt. 0.9999.and. sdir(j) .lt. 1.0001 ) then  ! T.Sato 2024/07/20 allow intermediate value

               sdir(j) =  1.0d0

            else if( sdir(j) .lt. -0.9999.and. sdir(j) .gt.-1.0001) then

               sdir(j) = -1.0d0

            else if( sdir(j) .ge. 1.5 .and. sdir(j) .lt. 2.5 ) then

               sdir(j) = 2.0

            else if( sdir(j) .le. -1.5 .and. sdir(j) .gt. -2.5 ) then

               sdir(j) = -2.0

            end if

            isbias(j) = nint(spava(109)) ! T.Sato, 0:normal, 1:dead particle is not produced, 2020/03/03

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 10 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) sz0(j) = spava(7)

               sdir(j)   = spava(11)
               sr1(j)    = spava(35)
               sr2(j)    = spava(36)

            if(sdir(j).lt.-3.5) sdir(j)=-3.0  ! for s-type = 10, dir = cr and iso are equivalent

            if(sdir(j).gt.-3.5.and.sdir(j).lt.-2.5) then ! for 'iso' source, T.Sato 2017/11/21
             if( sr2(j) .gt. sr1(j) ) then
              write(ErrCha,'("s-type = 10 error")')
              ErrID = 'L:17488/R:sours/F:read02.f' !E03_046_001
              call ErrWrite(ErrID,ErrCha)
              goto 885
             endif
             ! T.Sato 2020/12/13, limitation of zenith angle
             if(ischn(71).eq.0) then
              sag1(j)=-1.0d0
             else
              sag1(j)=spava(71)
             endif
             if(ischn(72).eq.0) then
              sag2(j)=1.0d0
             else
              sag2(j)=spava(72)
             endif
             if(sag1(j).ge.sag2(j)) then
              goto 926
             endif
             ! T.Sato 2021/02/26, limitation of azimuth angle
             if(ischn(126).eq.0) then
              spg1(j)=0.0d0
             else
              spg1(j)=spava(126)
             endif
             if(ischn(127).eq.0) then
              spg2(j)=360.0d0
             else
              spg2(j)=spava(127)
             endif
! T.Sato 2024/09/07
             if(spg1(j).lt.-360d0.or.spg1(j).gt.360d0) goto 923
             if(spg2(j).lt.-360d0.or.spg2(j).gt.360d0) goto 924
             if(spg1(j).lt.0.0) spg1(j)=spg1(j)+360.0d0
             if(spg2(j).lt.0.0) spg2(j)=spg2(j)+360.0d0
            else
             if( sr1(j) .gt. sr2(j) ) then
              write(ErrCha,'("s-type = 10 error")')
              ErrID = 'L:17525/R:sours/F:read02.f' !E03_046_002
              call ErrWrite(ErrID,ErrCha)
              goto 985
             endif
            endif

cFURUTA20191129 bugfix
            if(sr2(j).eq.0.0)then
             write(6,*)'Warning: r2=0 is not allowed in s-type=9'
             write(6,*)'         r2=1.0E-10 is used instead'
             sr2(j)=1.0d-10
            endif

            if( sdir(j) .eq. 200.0 .or. sdir(j) .eq. 0.0 ) then

               sdir(j) =  200.0

            else if( sdir(j) .gt. 0.9999.and. sdir(j) .lt. 1.0001 ) then  ! T.Sato 2024/07/20 allow intermediate value

               sdir(j) =  1.0d0

            else if( sdir(j) .lt. -0.9999.and. sdir(j) .gt.-1.0001) then

               sdir(j) = -1.0d0

            else if( sdir(j) .ge. 1.5 .and. sdir(j) .lt. 2.5 ) then

               sdir(j) = 2.0

            else if( sdir(j) .le. -1.5 .and. sdir(j) .gt. -2.5 ) then

               sdir(j) = -2.0

            end if

            isbias(j) = nint(spava(109)) ! T.Sato, 0:normal, 1:dead particle is not produced, 2020/03/03

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 11 .or. jstyp(j) .eq. -11) then ! T.Sato 2020/11/29 s-type=-11 is used for e-type case


               sx0(j) = 0.0
               sy0(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)

               sx1(j)    = spava( 4)
               sy1(j)    = spava( 6)
               sz0(j)    = spava( 7)
               sz1(j)    = spava( 8)
               se0(j)    = spava(10)
               sdir(j)   = spava(11)
               srx(j)    = spava(12)
               sry(j)    = spava(13)
               swem(j)   = spava(14)

! T.Sato 2020/09/21, Gaussian mode
               sx2(j)    = spava(85)
               sy2(j)    = spava(86)
               sxmrad1(j)= spava(110)
               sxmrad2(j)= spava(111)
               symrad1(j)= spava(112)
               symrad2(j)= spava(113)

               if( sdir(j) .eq. 200.0 .or. sdir(j) .eq. 0.0 ) then
                write(ErrCha,'("s-type = 11 error")')
                ErrID = 'L:17593/R:sours/F:read02.f' !E03_047_001
                call ErrWrite(ErrID,ErrCha)
                  goto 992

               else if( sdir(j) .gt. 0.0 ) then

                  sdir(j) =  1.0

               else if( sdir(j) .lt. 0.0 ) then

                  sdir(j) = -1.0

               end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 12 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) sz0(j) = spava(7)

               sdir(j)   = spava(11)

               if( sdir(j) .eq. 200.0 .or. sdir(j) .eq. 0.0 ) then
              write(ErrCha,'("s-type = 12 error")')
              ErrID = 'L:17624/R:sours/F:read02.f' !E03_047_002
              call ErrWrite(ErrID,ErrCha)
                  goto 991

               else if( sdir(j) .gt. 0.0 ) then

                  sdir(j) =  1.0

               else if( sdir(j) .lt. 0.0 ) then

                  sdir(j) = -1.0

               end if

            inquire( file = sfile(j), exist = exex )

            if( exex .eqv. .false. ) then
              write(ErrCha,'("s-type = 12 error")')
              ErrID = 'L:17642/R:sours/F:read02.f' !E03_048_001
              call ErrWrite(ErrID,ErrCha)
              goto 995
            endif

               isorf(j) = 4000 + j

               open(isorf(j), file = sfile(j), status = 'old' )

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 13 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0
               sz1(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) then
                                     sz0(j) = spava(7)
                                     sz1(j) = sz0(j)
               end if
               if( ischn(8) .ne. 0 ) sz1(j) = spava(8)

               se0(j)    = spava(10)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)
               sr1(j)    = spava(35)

! T.Sato 2024/07/14, add warning
         if(sr1(j).lt.0.0) then
          sr1(j) =abs(sr1(j))
          write(ErrCha,'("Warning: r1 in [source] is changed to",es12.4,
     &    " because it should be positive for s-type = 13")') sr1(j)
          ErrID = 'L:17680/R:sours/F:read02.f' !E03_042_001
          call ErrWrite(ErrID,ErrCha)
         endif

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 14 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0
               sz1(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) then
                                     sz0(j) = spava(7)
                                     sz1(j) = sz0(j)
               end if
               if( ischn(8) .ne. 0 ) sz1(j) = spava(8)

               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)
               sr1(j)    = spava(35)

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 15 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0
               sz1(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) then
                                     sz0(j) = spava(7)
                                     sz1(j) = sz0(j)
               end if
               if( ischn(8) .ne. 0 ) sz1(j) = spava(8)

               se0(j)    = spava(10)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)
               sr1(j)    = spava(35)

               isrn(j) = 2

               if( ischn(15) .ne. 0 ) isrn(j)   = nint( spava(15) )
               if( isrn(j) .lt. 2 .or. isrn(j)/2*2 .ne. isrn(j) ) then
                write(ErrCha,'("s-type = 15 error")')
                ErrID = 'L:17736/R:sours/F:read02.f' !E03_049_001
                call ErrWrite(ErrID,ErrCha)
                goto 986
               endif

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 16 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0
               sz1(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) then
                                     sz0(j) = spava(7)
                                     sz1(j) = sz0(j)
               end if
               if( ischn(8) .ne. 0 ) sz1(j) = spava(8)

               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)
               sr1(j)    = spava(35)

               isrn(j) = 2

               if( ischn(15) .ne. 0 ) isrn(j)   = nint( spava(15) )
               if( isrn(j) .lt. 2 .or. isrn(j)/2*2 .ne. isrn(j) ) then
                write(ErrCha,'("s-type = 16 error")')
                ErrID = 'L:17769/R:sours/F:read02.f' !E03_050_001
                call ErrWrite(ErrID,ErrCha)
                goto 986
               endif

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 18 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0
               sr0(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) sz0(j) = spava(7)
               if( ischn(9) .ne. 0 ) sr0(j) = spava(9)

               sx1(j) = spava(4)
               sy1(j) = spava(6)
               sz1(j) = spava(8)
               sr1(j) = spava(35)
               sr2(j) = spava(36)

               se0(j)    = spava(10)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 19 ) then


               sx0(j) = 0.0
               sy0(j) = 0.0
               sz0(j) = 0.0
               sr0(j) = 0.0

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) sz0(j) = spava(7)
               if( ischn(9) .ne. 0 ) sr0(j) = spava(9)

               sx1(j) = spava(4)
               sy1(j) = spava(6)
               sz1(j) = spava(8)
               sr1(j) = spava(35)
               sr2(j) = spava(36)

               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 20 ) then


               sx0(j) = spava(3)
               sy0(j) = spava(5)
               sz0(j) = spava(7)

               sx1(j) = spava(4)
               sy1(j) = spava(6)
               sz1(j) = spava(8)

               sr0(j) = spava(85)
               sr1(j) = spava(86)
               sr2(j) = spava(87)

               srx(j) = spava(88)
               sry(j) = spava(89)
               swem(j)= spava(90)

               sdl0(j) = spava(91)

               if( spava(91) .lt. 0.0d0 ) then
                write(ErrCha,'("s-type = 20 error")')
                ErrID = 'L:17850/R:sours/F:read02.f' !E03_051_001
                call ErrWrite(ErrID,ErrCha)
                goto 940
               endif

               se0(j)    = spava(10)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 21 ) then


               sx0(j) = spava(3)
               sy0(j) = spava(5)
               sz0(j) = spava(7)

               sx1(j) = spava(4)
               sy1(j) = spava(6)
               sz1(j) = spava(8)

               sr0(j) = spava(85)
               sr1(j) = spava(86)
               sr2(j) = spava(87)

               srx(j) = spava(88)
               sry(j) = spava(89)
               swem(j)= spava(90)

               sdl0(j) = spava(91)

               if( spava(91) .lt. 0.0d0 ) then
                write(ErrCha,'("s-type = 21 error")')
                ErrID = 'L:17885/R:sours/F:read02.f' !E03_052_001
                call ErrWrite(ErrID,ErrCha)
                goto 959
               endif

               sdir(j) = spava(11)
               sphi(j) = spava(27)
               sdom(j) = spava(28)

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 22 .or. jstyp(j) .eq. 23 ) then

            if( jstyp(j) .eq. 22 ) then


               se0(j)    = spava(10)
               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

            else if( jstyp(j) .eq. 23 ) then


               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

            end if

               isbias(j) = nint(spava(109)) ! T.Sato, 0:normal, 1:weight bias mode

               istdd(j) = istdg
               istcc(j) = istdc

               call moddas_reallocate_dbl(
     &                 isrc, j, jnx+1, istxx, das_istxx)
               das_istxx(istxx(j):istxx(j)+jnx) = gmsh(istxg:istxg+jnx)

               call moddas_reallocate_dbl(
     &                 isrc, j, jny+1, istyy, das_istyy)
               das_istyy(istyy(j):istyy(j)+jny) = gmsh(istyg:istyg+jny)

               call moddas_reallocate_dbl(
     &                 isrc, j, jnz+1, istzz, das_istzz)
               das_istzz(istzz(j):istzz(j)+jnz) = gmsh(istzg:istzg+jnz)

               isxtp(j) = ixtp
               isinx(j) = jnx
               sxmin(j) = xmin
               sxmax(j) = xmax
               sxdel(j) = xdel

               isytp(j) = iytp
               isiny(j) = jny
               symin(j) = ymin
               symax(j) = ymax
               sydel(j) = ydel

               isztp(j) = iztp
               isinz(j) = jnz
               szmin(j) = zmin
               szmax(j) = zmax
               szdel(j) = zdel

*-----------------------------------------------------------------------

                  sek = 0.d0
                  isbmax = 0

               do k = 1, inxyz

                  sek = sek + das(istdg + k )
                  if(das(istdg + k ).ne.0.0) isbmax=isbmax+1

               end do

               if(isbias(j).eq.0) then  ! T.Sato, normal bias mode

                  das(istdc + 1 ) = das(istdg + 1 ) / sek

                  do k = 2, inxyz

                     das(istdc + k ) = das(istdc + k - 1 )
     &                               + das(istdg + k ) / sek

                  end do

               elseif(isbias(j).ge.1) then  ! weight bias mode

                  do k = 1, inxyz

                     das(istdg + k ) = das(istdg + k )*isbmax / sek

                  end do

                  if(isbias(j).ge.2) then  ! equal distribution mode

                     if(das(istdg+1).eq.0.0) then
                        das(istdc+1) = 0.0
                     else
                        das(istdc+1) = 1.0d0
                     endif

                     do k = 2, inxyz

                        if(das(istdg+k).eq.0.0) then
                           das(istdc+k)=das(istdc+k-1)
                        else
                           das(istdc+k)=das(istdc+k-1)+1.0d0
                        endif

                     end do

                  endif

               endif

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 17 ) then

            if( ischn(26) .eq. 0 ) then
             write(ErrCha,'("s-type = 17 error")')
             ErrID = 'L:18009/R:sours/F:read02.f' !E03_053_001
             call ErrWrite(ErrID,ErrCha)
             goto 997
            endif

            if( jpsf.gt.0 .and. ischn(7) .eq. 0) then
             write(ErrCha,'("jpsf >0 requires z0 specification")')
             ErrID = 'L:18016/R:sours/F:read02.f' !E03_053_001
             call ErrWrite(ErrID,ErrCha)
             goto 997
            endif

               sx0(j) =  2.0 * bign
               sy0(j) =  2.0 * bign
               sz0(j) =  2.0 * bign

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) sz0(j) = spava(7)

               sx1(j) = sx0(j)
               sy1(j) = sy0(j)
               sz1(j) = sz0(j)

               if( ischn(4) .ne. 0 ) sx1(j) = spava(4)
               if( ischn(6) .ne. 0 ) sy1(j) = spava(6)
               if( ischn(8) .ne. 0 ) sz1(j) = spava(8)

               sdir(j)   = 10000.0

            if( ischn(11) .ne. 0 ) then

               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

            end if

               swt0(j)   = -10000.0
               if( ischn(39) .ne. 0 ) swt0(j) = spava(39)

               se0(j)    = -10000.0
               jetyp(j)  = -10000

               if( ischn(10) .ne. 0 ) se0(j) = spava(10)

            if( ischn(18) .ne. 0 ) then

               se0(j)    = -10000.0
               jetyp(j)  = nint( spava(18) )

            end if

               jttyp(j)  = -10000

            if( ischn(29) .ne. 0 ) then

               jttyp(j)  = nint( spava(29) )

            end if

*-----------------------------------------------------------------------

            if( jsdmp(j,1) .eq. 0 .and. istyp(j) .le. -1000 ) then
             write(ErrCha,'("the other parameters error")')
             ErrID = 'L:18074/R:sours/F:read02.f' !E03_054_001
             call ErrWrite(ErrID,ErrCha)
             goto 981
            endif
            if( jsdmp(j,2) .eq. 0 .and. sx0(j) .le. bign ) then
             write(ErrCha,'("the other parameters error")')
             ErrID = 'L:18080/R:sours/F:read02.f' !E03_054_002
             call ErrWrite(ErrID,ErrCha)
             goto 981
            endif

            if( jsdmp(j,3) .eq. 0 .and. sy0(j) .le. bign ) then
             write(ErrCha,'("the other parameters error")')
             ErrID = 'L:18087/R:sours/F:read02.f' !E03_054_003
             call ErrWrite(ErrID,ErrCha)
             goto 981
            endif

            if( jsdmp(j,4) .eq. 0 .and. sz0(j) .le. bign ) then
             write(ErrCha,'("the other parameters error")')
             ErrID = 'L:18094/R:sours/F:read02.f' !E03_054_004
             call ErrWrite(ErrID,ErrCha)
             goto 981
            endif
            if( jsdmp(j,5) .eq. 0 .and. sdir(j) .gt. 1000.0 ) then
             write(ErrCha,'("the other parameters error")')
             ErrID = 'L:18100/R:sours/F:read02.f' !E03_054_005
             call ErrWrite(ErrID,ErrCha)
             goto 981
            endif

            if( jsdmp(j,6) .eq. 0 .and. sdir(j) .gt. 1000.0 ) then
             write(ErrCha,'("the other parameters error")')
             ErrID = 'L:18107/R:sours/F:read02.f' !E03_054_006
             call ErrWrite(ErrID,ErrCha)
             goto 981
            endif

            if( jsdmp(j,7) .eq. 0 .and. sdir(j) .gt. 1000.0 ) then
             write(ErrCha,'("the other parameters error")')
             ErrID = 'L:18114/R:sours/F:read02.f' !E03_054_007
             call ErrWrite(ErrID,ErrCha)
             goto 981
            endif

            if( jsdmp(j,8) .eq. 0 .and. se0(j) .le. -1000.0 .and.
     &          jetyp(j) .le. -1000 ) then
             write(ErrCha,'("the other parameters error")')
             ErrID = 'L:18122/R:sours/F:read02.f' !E03_054_008
             call ErrWrite(ErrID,ErrCha)
             goto 981
            endif


            if( jsdmp(j,9) .eq. 0 .and. swt0(j) .le. -1000.0 )
     &          swt0(j) = 1.d0

            if( jsdmp(j,1) .gt. 0 .and. istyp(j) .gt. -1000 )
     &          jsdmp(j,1) = 0
            if( jsdmp(j,2) .gt. 0 .and. sx0(j) .gt. bign )
     &          jsdmp(j,2) = 0
            if( jsdmp(j,3) .gt. 0 .and. sy0(j) .gt. bign )
     &          jsdmp(j,3) = 0
            if( jsdmp(j,4) .gt. 0 .and. sz0(j) .gt. bign )
     &          jsdmp(j,4) = 0

            if( sdir(j) .lt. 1000.0 ) then
                 jsdmp(j,5) = 0
                 jsdmp(j,6) = 0
                 jsdmp(j,7) = 0
            end if

            if( jsdmp(j,8) .gt. 0 .and.
     &        ( se0(j) .gt. -1000.0 .or. jetyp(j) .gt. -1000 ) )
     &          jsdmp(j,8) = 0
            if( jsdmp(j,9) .gt. 0 .and. swt0(j) .gt. -10000.0 )
     &          jsdmp(j,9) = 0
            if( jsdmp(j,10) .gt. 0 .and. jttyp(j) .gt. -1000 )
     &          jsdmp(j,10) = 0

               ssr = ssx(j)**2 + ssy(j)**2 + ssz(j)**2

            if( ssr .gt. 1.d-8 ) then
               jsdmp(j,14) = 0
               jsdmp(j,15) = 0
               jsdmp(j,16) = 0
            end if

*-----------------------------------------------------------------------

            if( npe .le. 1 ) then

                     filnm = sfile(j)(1:lsfile(j))
                     lfiln = lsfile(j)

            else if( me .gt. 0 ) then

                     iorder = aint(log10(real(npe)))+1
                     if ( iorder .lt. 3) iorder = 3

                     write(chme,'(i5.5)') me

                     filnm = sfile(j)(1:lsfile(j))
     &                       //'.' // chme(6-iorder:5)
                     lfiln = lsfile(j)+1+iorder

            else if(me.eq.0.and.npe.gt.1)then

             iorder = aint(log10(real(npe)))+1
             if ( iorder .lt. 3) iorder = 3

             lfiln = lsfile(j)+1+iorder
             do k=1,10**iorder-1
              write(chme,'(i5.5)') k
              filnm = sfile(j)(1:lsfile(j))
     &             //'.' // chme(6-iorder:5)
              inquire( file = filnm, exist = exex )
              if( exex .eqv. .false. ) then
               exit
              endif
             enddo
             ndumpmax=k-1

             if(npe-1.gt.ndumpmax)then
              write(ErrCha,'("dump error")')
              ErrID = 'L:18199/R:sours/F:read02.f' !E03_055_001
              call ErrWrite(ErrID,ErrCha)
              goto 855
             elseif(npe-1.lt.ndumpmax)then
              write(6,'("Warning: # of prepared dump file is ",i5)')
     &             ndumpmax
              write(6,'("         But # of MPI specified is (npe-1)",
     &             " = ",i5)') npe-1
              write(6,*)'You may need to re-distribute dump files. ',
     &             'Ref. Sec. 8 of the PHITS manual.'
             endif

            end if

            if( me .gt. 0 .or. npe .le. 1 ) then

                  inquire( file = filnm, exist = exex )
                  if( exex .eqv. .false. ) then
                  write(ErrCha,'("dump error")')
                  ErrID = 'L:18218/R:sours/F:read02.f' !E03_055_001
                  call ErrWrite(ErrID,ErrCha)
                  goto 955
                  endif

                  isorf(j) = 4000 + j

               if( isdmp(j,0) .gt. 0 ) then

                  open(isorf(j), file = filnm,
     &                        form='unformatted',status = 'old' )

                  read(isorf(j),iostat=ios,err=954) dummy

               else

                  open(isorf(j), file = filnm,
     &                        form='formatted',status = 'old' )

                  read(isorf(j),*,iostat=ios,err=954) dummy

               end if

                  if( ios .eq. -1 ) then
             write(ErrCha,'("dump error")')
             ErrID = 'L:18243/R:sours/F:read02.f' !E03_055_002
             call ErrWrite(ErrID,ErrCha)
                    goto 953
                  endif

                  rewind isorf(j)

            end if
            if(me.eq.0)then
             if(idmpmode.eq.1)then
              do i=lsfile(j),1,-1
               if(sfile(j)(i:i).eq.'.')exit
              enddo
              if(i.eq.0)i=lsfile(j)+1
              filnm(1:i-5) = sfile(j)(1:i-5)
              filnm(i-4:lsfile(j)-4) = sfile(j)(i:lsfile(j))
              lfiln = lsfile(j)-4
              inquire( file = filnm(1:lfiln), exist = exex )
              if( exex .eqv. .false. ) then
             write(ErrCha,'("dump error")')
             ErrID = 'L:18263/R:sours/F:read02.f' !E03_055_003
             call ErrWrite(ErrID,ErrCha)
                goto 955
              endif
             endif
            endif
c------------------
*-----------------------------------------------------------------------
         else if( jstyp(j) .eq. 24 .or. jstyp(j) .eq. 25 ) then
               if(ischn(100).eq.0) then
             write(ErrCha,'("s-type = 24 or 25 error")')
             ErrID = 'L:18274/R:sours/F:read02.f' !E03_056_001
             call ErrWrite(ErrID,ErrCha)
                goto 951
               endif
               itetreg(j) = spava(100)
               sdir(j)    = spava(11)
               sphi(j)    = spava(27)
               sdom(j)    = spava(28)
               if( jstyp(j) .eq. 24 ) se0(j)    = spava(10)

*-----------------------------------------------------------------------
         else if( jstyp(j) .eq. 100 ) then

               sx0(j) = 2.0 * bign
               sy0(j) = 2.0 * bign
               sz0(j) = 2.0 * bign

               if( ischn(3) .ne. 0 ) sx0(j) = spava(3)
               if( ischn(5) .ne. 0 ) sy0(j) = spava(5)
               if( ischn(7) .ne. 0 ) sz0(j) = spava(7)

               sx1(j) = sx0(j)
               sy1(j) = sy0(j)
               sz1(j) = sz0(j)

               if( ischn(4) .ne. 0 ) sx1(j) = spava(4)
               if( ischn(6) .ne. 0 ) sy1(j) = spava(6)
               if( ischn(8) .ne. 0 ) sz1(j) = spava(8)

               sdir(j)   = 10000.0

            if( ischn(11) .ne. 0 ) then

               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

            end if

               swt0(j)   = -10000.0
               if( ischn(39) .ne. 0 ) swt0(j) = spava(39)

               se0(j)    = -10000.0
               jetyp(j)  = -10000

               if( ischn(10) .ne. 0 ) se0(j) = spava(10)

            if( ischn(18) .ne. 0 ) then

               se0(j)    = -10000.0
               jetyp(j)  = nint( spava(18) )

            end if

               jttyp(j)  = -10000

            if( ischn(29) .ne. 0 ) then

               jttyp(j)  = nint( spava(29) )

            end if

*-----------------------------------------------------------------------

            if( istyp(j) .gt. -1000  ) jsusr(j,1) = 1

            if( sx0(j) .gt. bign ) jsusr(j,2) = 1
            if( sy0(j) .gt. bign ) jsusr(j,3) = 1
            if( sz0(j) .gt. bign ) jsusr(j,4) = 1

            if( sdir(j) .lt. 1000.0 ) then
                jsusr(j,5) = 1
                jsusr(j,6) = 1
                jsusr(j,7) = 1
            end if

            if( se0(j) .gt. -1000.0 .or. jetyp(j) .gt. -1000 )
     &         jsusr(j,8) = 1

            if( swt0(j) .gt. -10000.0 ) jsusr(j,9) = 1

            if(  jttyp(j) .gt. -1000 ) jsusr(j,10) = 1

               ssr = ssx(j)**2 + ssy(j)**2 + ssz(j)**2

            if( ssr .gt. 1.d-8 ) then
               jsusr(j,14) = 1
               jsusr(j,15) = 1
               jsusr(j,16) = 1
            end if

            if( ischn(37) .ne. 0 ) jsusr(j,17) = 1
            if( ischn(38) .ne. 0 ) jsusr(j,18) = 1
            if( ischn(40) .ne. 0 ) jsusr(j,19) = 1
            if( ischn(41) .ne. 0 ) jsusr(j,19) = 1
            if( ischn(43) .ne. 0 ) jsusr(j,20) = 1

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 26 ) then

               se0(j)    = -10000.0
               jetyp(j)  = -1
            if( ischn(18) .ne. 0 ) then
               jetyp(j)  = nint( spava(18) )
            else if( ischn(10) .ne. 0 ) then
               se0(j) = spava(10)
            end if

               sdir(j)   = spava(11)
               sphi(j)   = spava(27)
               sdom(j)   = spava(28)

               issuf(j)  = nint( spava(106) )
               if( issuf(j) .le. 0 ) goto 938

               iscut(j)  = incut

*-----------------------------------------------------------------------

         else

          ierrMSG=4 ! T.Sato 2017/07/13
           write(ErrCha,'("s-type = 100 error.")')
           ErrID = 'L:18398/R:sours/F:read02.f' !E03_057_001
           call ErrWrite(ErrID,ErrCha)
           goto 997

         end if

*-----------------------------------------------------------------------
*     source energy group
*-----------------------------------------------------------------------
ccse 2022.03.28 proj=all, e-type not 28,29
         jwetyp = nint( spava(18) )
         if(jwetyp.ne.28 .and. jwetyp.ne.29 .and.istyp(j).eq.20) then
           write(ErrCha,'("proj=all but e-type is not 28 nor 29")')
           ErrID = 'L:18411/R:sours/F:read02.f' !E03_058_001
           call ErrWrite(ErrID,ErrCha)
           goto 997
         endif


         if( jstyp(j) .eq.   4 .or.
     &       jstyp(j) .eq.   5 .or.
     &       jstyp(j) .eq.   6 .or.
     &       jstyp(j) .eq.   8 .or.
     &       jstyp(j) .eq.  10 .or.
     &       jstyp(j) .eq. -11 .or.  ! T.Sato 2020/11/29
     &       jstyp(j) .eq.  14 .or.
     &       jstyp(j) .eq.  16 .or.
     &       jstyp(j) .eq.  19 .or.
     &       jstyp(j) .eq.  21 .or.
     &       jstyp(j) .eq.  23 .or.
     &       ( jstyp(j) .eq.  26 .and. jetyp(j) .gt. 0 ) .or.
     &       jstyp(j) .eq.  25 .or.
     &        ( jstyp(j) .eq.  17 .and. jetyp(j) .gt. 0 ) .or.
     &     ( jstyp(j) .eq. 100 .and. jetyp(j) .gt. 0 ) ) then

            if( ischn(18) .eq. 0 ) then
             write(ErrCha,'("e-type error.")')
             ErrID = 'L:18435/R:sours/F:read02.f' !E03_058_001
             call ErrWrite(ErrID,ErrCha)
              goto 997
            endif

               jetyp(j)  = nint( spava(18) )

*-----------------------------------------------------------------------
            if( jetyp(j) .eq. 1 .or. jetyp(j) .eq. 11 ) then

               if( ischn(17) .eq. 0 ) then
             write(ErrCha,
     & '("There is an error in ne in e-type = 1 or e-type = 11.")')
             ErrID = 'L:18448/R:sours/F:read02.f' !E03_059_001
             call ErrWrite(ErrID,ErrCha)
                goto 997
               endif

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 21 .or. jetyp(j) .eq. 31 ) then

               if( ischn(17) .eq. 0 ) then
             write(ErrCha,
     & '("There is an error in ne in e-type = 21 or e-type = 31.")')
             ErrID = 'L:18460/R:sours/F:read02.f' !E03_060_001
             call ErrWrite(ErrID,ErrCha)
                goto 997
               endif

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 8 .or. jetyp(j) .eq. 18 ) then

               if( ischn(17) .eq. 0 ) then
             write(ErrCha,
     & '("There is an error in ne in e-type = 8 or e-type = 18.")')
             ErrID = 'L:18472/R:sours/F:read02.f' !E03_061_001
             call ErrWrite(ErrID,ErrCha)
                goto 997
               endif

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 22 .or. jetyp(j) .eq. 32 ) then

               if( ischn(17) .eq. 0 ) then
             write(ErrCha,
     & '("There is an error in ne in e-type = 22 or e-type = 32.")')
             ErrID = 'L:18484/R:sours/F:read02.f' !E03_062_001
             call ErrWrite(ErrID,ErrCha)
                goto 997
              endif

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 5 .or. jetyp(j) .eq. 15 ) then

               if( ischn(68) .eq. 0 .or. ischn(20) .eq. 0 .or.
     &             ischn(21) .eq. 0 .or. ischn(67) .eq. 0 ) goto 997

                  seg1(j)   = spava(20)
                  seg2(j)   = spava(21)
                  isnm(j)   = spava(68)

                  eming = seg1(j)
                  emaxg = seg2(j)

                  if( eming .gt. emaxg ) then
             write(ErrCha,*)
     & 'There is an error in nm or eg1 or eg2 or',
     & ' f (x) in e-type = 5 or e-type = 15.'
             ErrID = 'L:18507/R:sours/F:read02.f' !E03_063_001
             call ErrWrite(ErrID,ErrCha)
                    goto 997
                  endif

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 4 .or. jetyp(j) .eq. 14 ) then

               if( ischn(17) .eq. 0 ) then
             write(ErrCha,
     & '("In e-type = 4 or e-type = 14, there is an error in p-type.")')
             ErrID = 'L:18519/R:sours/F:read02.f' !E03_064_001
             call ErrWrite(ErrID,ErrCha)
                goto 997
               endif

               if( ischn(63) .eq. 0 ) then

                  do i = 1, ngrp(j)

                     prw(ngpi(j)+i) = 1.0d0

                  end do

               end if

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 24 .or. jetyp(j) .eq. 34 ) then

               if( ischn(17) .eq. 0 ) then
             write(ErrCha,
     & '("There is an error in p-type in e-type = 24 or e-type = 34.")')
             ErrID = 'L:18541/R:sours/F:read02.f' !E03_065_001
             call ErrWrite(ErrID,ErrCha)
                goto 997
               endif

               if( ischn(63) .eq. 0 ) then

                  do i = 1, ngrp(j)

                     prw(ngpi(j)+i) = 1.0d0

                  end do

               end if

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 9 .or. jetyp(j) .eq. 19 ) then

               if( ischn(17) .eq. 0 ) then
             write(ErrCha,
     & '("In e-type = 9 or e-type = 19, there is an error in ne.")')
             ErrID = 'L:18563/R:sours/F:read02.f' !E03_066_001
             call ErrWrite(ErrID,ErrCha)
                goto 997
               endif


               if( ischn(63) .eq. 0 ) then

                  do i = 1, ngrp(j)

                     prw(ngpi(j)+i) = 1.0d0

                  end do

               end if

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 23 .or. jetyp(j) .eq. 33 ) then

               if( ischn(17) .eq. 0 ) then
             write(ErrCha,
     & '("There is an error in ne in e-type = 23 or e-type = 33.")')
             ErrID = 'L:18586/R:sours/F:read02.f' !E03_067_001
             call ErrWrite(ErrID,ErrCha)
                goto 997
               endif

               if( ischn(63) .eq. 0 ) then

                  do i = 1, ngrp(j)

                     prw(ngpi(j)+i) = 1.0d0

                  end do

               end if

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 6 .or. jetyp(j) .eq. 16 ) then

               if( ischn(68) .eq. 0 .or. ischn(20) .eq. 0 .or.
     &             ischn(21) .eq. 0 .or. ischn(67) .eq. 0 ) then
             write(ErrCha,*)
     & 'There is an error in nm or eg1 or eg2 or ',
     & 'f (x) in e-type = 6 or e-type = 16.'
             ErrID = 'L:18610/R:sours/F:read02.f' !E03_068_001
             call ErrWrite(ErrID,ErrCha)
               goto 997
             endif

                  seg1(j)   = spava(20)
                  seg2(j)   = spava(21)
                  isnm(j)   = spava(68)

               if( ischn(63) .eq. 0 ) then

                  do i = 1, ngrp(j)

                     prw(ngpi(j)+i) = 1.0d0

                  end do

               end if

                  eming = seg1(j)
                  emaxg = seg2(j)

                  if( eming .gt. emaxg ) then
             write(ErrCha,'("eg1 > eg2")')
             ErrID = 'L:18634/R:sours/F:read02.f' !E03_069_001
             call ErrWrite(ErrID,ErrCha)
                    goto 997
                  endif

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 2 .or. jetyp(j) .eq. 12 ) then

               if( ischn(19) .eq. 0 .or. ischn(20) .eq. 0 .or.
     &             ischn(21) .eq. 0 .or. ischn(22) .eq. 0 ) then
             write(ErrCha,*)
     & 'There is an error in eg 0 or eg 1 or eg 2 or eg 3 ',
     & 'in e-type = 2 or e-type = 12.'
             ErrID = 'L:18648/R:sours/F:read02.f' !E03_070_001
             call ErrWrite(ErrID,ErrCha)
               goto 997
             endif

                  seg0(j)   = spava(19)
                  seg1(j)   = spava(20)
                  seg2(j)   = spava(21)
                  seg3(j)   = spava(22)

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 3 .or. jetyp(j) .eq. 7 ) then

               if( ischn(23) .eq. 0 .or. ischn(24) .eq. 0 .or.
     &             ischn(25) .eq. 0 ) then
             write(ErrCha,*)
     & 'In e-type = 3 or e-type = 7 there is an error ',
     & 'in et 0, et 1 or et 2.'
             ErrID = 'L:18667/R:sours/F:read02.f' !E03_071_001
             call ErrWrite(ErrID,ErrCha)
               goto 997
               endif

                  set0(j)   = spava(23)
                  set1(j)   = spava(24)
                  set2(j)   = spava(25)
                  set3(j)   = spava(108) ! power index of Maxwellian, T.Sato 2020/02/05

               if( ischn(68) .eq. 0 ) then

                  ischn(68) = 1
                  spava(68) = -200

                  ngrp(j)  = abs( nint( spava(68) ) )

                  ngll(j)  = 1
                  if( spava(68) .lt. 0.d0 ) ngll(j) = -1

                  call moddas_reallocate_dbl(
     &                    isrc, j, ngrp(j), ngei, egmin)
                  call moddas_reallocate_dbl(
     &                    isrc, j, ngrp(j), ngea, egmax)
                  call moddas_reallocate_dbl(
     &                    isrc, j, ngrp(j), ngfe, fegrp)
                  call moddas_reallocate_dbl(
     &                    isrc, j, ngrp(j), ngft, rfe)
                  call moddas_reallocate_dbl(
     &                    isrc, j, ngrp(j), ngpi, prw)
                  call moddas_reallocate_dbl(
     &                    isrc, j, ngrp(j), ngpw, pwt)

               end if

                  eming = set1(j)
                  emaxg = set2(j)

                  if( eming .gt. emaxg ) then
             write(ErrCha,'("et1> et2 for e-type = 3 or e-type = 7.")')
             ErrID = 'L:18707/R:sours/F:read02.f' !E03_071_002
             call ErrWrite(ErrID,ErrCha)
                    goto 997
                  endif

               if( jetyp(j) .eq. 7 .and. ischn(63) .eq. 0 ) then

                  do i = 1, ngrp(j)

                     prw(ngpi(j)+i) = 1.0d0

                  end do

               end if

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 28 .or. jetyp(j) .eq. 29 ) then

               if( ischn(17) .eq. 0 .and. ischn(94) .eq. 0 ) then
             write(ErrCha,'("e-type=28, 29 error")')
             ErrID = 'L:18728/R:sours/F:read02.f' !E03_072_001
             call ErrWrite(ErrID,ErrCha)
                goto 945
               endif
               if( ischn(95) .eq. 0 ) then
c                set default value
               end if
               if( ischn(96) .eq. 0 ) then
                 if( ischn(17) .eq. 0 .and. ischn(94) .eq. 1 ) then
*                  go to 946
                   decayt(j) = spava(96)
                 else
                   decayt(j) = spava(96)
                 end if
               end if
               if( ischn(97) .eq. 0 ) then
c                set default value
                   norm(j) = nint( spava(97) )
               end if
C MATSUDA 2017.05.29 (iaugers)
               if( ischn(99) .eq. 0 ) then
c                set default value
                   iaugers(j) = nint( spava(99) )
               end if
C MATSUDA 2018.08.15 (icharacterx)
               if( ischn(101) .eq. 0 ) then
c                set default value
                   icharacterx(j) = nint( spava(101) )
               end if

C MATSUDA 2019.05.08 (iannih)
               if( ischn(102) .eq. 0 ) then
c                set default value
                   iannih(j) = nint( spava(102) )
               end if

               if( ischn(17) .eq. 1 .and. ischn(94) .eq. 0 ) then

! MATSUDA 2017.05.29
                 if( inkf0(j) .eq. 11 .or. inkf0(j) .eq. -11 ) then
*                  ngrp(j) = nint( spava(17) )
                   ngll(j) =  1
c                       ! if negative, emin(1) must change from Zero.
                 else
*                  ngrp(j) = nint( spava(17) )
*                  ngll(j) =  1
                 end if

                 call moddas_reallocate_dbl(
     &                   isrc, j, ngrp(j), ngei, egmin)
                 call moddas_reallocate_dbl(
     &                   isrc, j, ngrp(j), ngea, egmax)
                 call moddas_reallocate_dbl(
     &                   isrc, j, ngrp(j), ngfe, fegrp)
                 call moddas_reallocate_dbl(isrc, j, ngrp(j), ngft, rfe)
                 call moddas_reallocate_dbl(isrc, j, ngrp(j), ngpi, prw)
                 call moddas_reallocate_dbl(isrc, j, ngrp(j), ngpw, pwt)

                 do n = 1, ngrp(j)

                   egmin(ngei(j)+n) = e_lowk(n)
                   egmax(ngea(j)+n) = e_upk(n)
                   fegrp(ngfe(j)+n) = ratkari(n)
                 end do

                 if( jetyp(j) .eq. 28 ) then
                   do n = 1, ngrp(j)
                     prw(ngpi(j)+n) = fegrp(ngfe(j)+n)
                   end do
                 else if( ischn(63) .eq. 1 ) then  ! jetyp(j) = 29
                   do n = 1, ngrp(j)
                     prw(ngpi(j)+n) = prbkari(n)
                   end do
                   deallocate(prbkari)
                 else if( ischn(63) .eq. 0 ) then  ! jetyp(j) = 29
                   do n = 1, ngrp(j)
                     prw(ngpi(j)+n) = 1.0d+00
                   end do
                 end if

                 niorg(j) = 0
                 normfact(j) = 1.0d+00

                 deallocate(e_lowk,e_upk,ratkari)

               else
                 j33 = j
                 nk0 = niorg(j)
*                karival = 1.0d+00 * normfact(j)

                 call setrisrc(j33,nk0,nk1,ierr)

                   if( ierr .ne. 0 ) then
                     l_err = ill(jsn)
                     k_err = jsn
                     ierr  = 1
                     return
                   end if

                   nicur(j) = nk1

                 if( jetyp(j) .eq. 28 ) then
                   do n = 1, ngrp(j)
                     prw(ngpi(j)+n) = fegrp(ngfe(j)+n)
                   end do
                 else  ! jetyp(j) = 29
                   do n = 1, ngrp(j)
                     prw(ngpi(j)+n) = 1.0d+00
                   end do
                 end if

               end if

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 20 ) then

             inquire( file = sfile(j), exist = exex )
             if( exex .eqv. .false. ) then
             write(ErrCha,'("e-type = 20 error.")')
             ErrID = 'L:18848/R:sours/F:read02.f' !E03_073_001
             call ErrWrite(ErrID,ErrCha)
              goto 995
             endif
             isorf(j) = 4000 + j
             open(isorf(j), file = sfile(j), status = 'old' )

             call etalfile(j,isorf(j),ndata,'#!$',ierr)

             close(isorf(j))

             if ( ierr .eq. 1 ) then
                m_err = 'This file cannot be used for source'
                ErrCha = ''
                ErrID = 'L:18862/R:sours/F:read02.f'
                l_err = ill(jsn)
                k_err = jsn
                return
             end if

*-----------------------------------------------------------------------

            else if( jetyp(j) .eq. 25 .or. jetyp(j) .eq. 26) then  ! T.Sato 2020/12/03  Cosmic-ray source

              ipcosmic(j)=10000  ! Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
              if(istyp(j).eq.1) then ! proton
               ipcosmic(j)=1
              elseif(istyp(j).eq.2) then ! neutron
               ipcosmic(j)=0
              elseif(istyp(j).eq.6) then ! muon+
               ipcosmic(j)=29
              elseif(istyp(j).eq.7) then ! muon-
               ipcosmic(j)=30
              elseif(istyp(j).eq.12) then ! electron
               ipcosmic(j)=31
              elseif(istyp(j).eq.13) then ! positron
               ipcosmic(j)=32
              elseif(istyp(j).eq.14) then ! photon
               ipcosmic(j)=33
              elseif(istyp(j).eq.18) then ! alpha
               ipcosmic(j)=2
              elseif(istyp(j).eq.19) then ! heavy ion
               ipcosmic(j)=inkf0(j)/1000000
              endif
              if(ipcosmic(j).ge.34) goto 929 ! no source information

              icyear=nint(spava(115))   ! Year
              icmonth=nint(spava(116))  ! Month
              icday=nint(spava(117))    ! Date
              glat=spava(118)     ! Latitude (deg), -90 =< glat =< 90
              glong=spava(119)    ! Longitude (deg), -180 =< glong =< 180
              alti=spava(120)     ! Altitude (km)
              icenv(j)=nint(spava(122)) ! icenv, <0:SEP, 0:GCR, 1:ideal atmosphere, 2:on ground, 3:aircraft (pilot), 4:aircraft(cabin), 5:blackhole
! ground: Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
              environ(j)=spava(121) ! remember original environ parameter
              if(icenv(j).eq.2) then
               ground(j)=spava(121) ! water density
              elseif(icenv(j).eq.3) then
               ground(j)=-spava(121) ! aircraft mass in 100 ton (-10 < ground(j) < 0)
              elseif(icenv(j).eq.4) then
               ground(j)=-spava(121)-10 ! aircraft mass in 100 ton (ground(j) < -10)
              elseif(icenv(j).eq.5) then
               ground(j)=100d0         ! blackhole
              else
               ground(j)=10d0    ! ideal atmosphere or in space
              endif

! ground: SPE event ID for icenv < 0
              if(icenv(j).eq.-1) then
               if(ischn(121).eq.1) then
                ground(j)=spava(121)
               else
                ground(j)=1.0d0 ! Feb.1956 event is the default
               endif
              endif

              idummy=0
              if(icenv(j).ge.0.and.ischn(123).eq.0) then
               solarmod(j)=getHP(icyear,icmonth,icday,idummy) ! Solar activity (W index) of the day
               if(idummy.ge.4) goto 925
              else
               solarmod(j)=spava(123)
              endif

              if(ischn(124).eq.0) then
               rigid(j)=getr(glat,glong)    ! Vertical cut-off rigidity (GV)
              else
               rigid(j)=spava(124)
              endif

              if(icenv(j).ge.1.and.ischn(125).eq.0) then
               glat=100.0 ! use US standard atmosphere 1986
               depatom(j)=getd(alti,glat)   ! Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.
              elseif(icenv(j).eq.0.or.icenv(j).eq.-2) then ! LEO mode, depatom should be altitude in km
               depatom(j)=alti
              else
               depatom(j)=spava(125)
              endif

       if((icenv(j).eq.0.and.depatom(j).ne.0.0).or.icenv(j).eq.-2) then ! LEO mode
        if(rigid(j).ne.0)
     &  write(*,'("Warning: rigid is ignored in the LEO mode in ",i3,
     &  "th multisource")') j
        if(ischn(125).eq.1)
     &  write(*,'("Warning: depatom is ignored in the LEO mode in ",i3,
     &  "th multisource")') j
       endif

       if((ischn(120).eq.1.or.ischn(125).eq.1).and.icenv(j).eq.-1)
     &  write(*,'("Warning: alti/depatom is ignored",
     &  " for icenv = -1 in ",i3,"th multisource")')j

              if(j.ge.2) then
               if(icenv(j).ne.icenv(j-1)) goto 928 ! cosmic-ray environment should be same for all multisource
              endif

              eval=1000.0
              tmp=getCosSpec(icenv(j),ipcosmic(j),solarmod(j),rigid(j),
     &        depatom(j),eval,ground(j)) ! dummy call to setup parameter

              if(ischn(20).eq.0) then ! default minimum energy
               if(icenv(j).le.0) then ! GCR or SEP mode
                seg1(j)=1.0
               elseif(ipcosmic(j).eq.0) then ! neutron
                seg1(j)=1.0e-8
               else
                seg1(j)=1.0e-2
               endif
              else
               seg1(j)   = spava(20)
              endif

              if(ischn(21).eq.0) then ! default maximum energy
               if(icenv(j).le.-1) then ! SEP or TP mode
                seg2(j)=1.0e5 ! upto 10 GeV
               elseif(ipcosmic(j).eq.0.or.ipcosmic(j).ge.31) then
                seg2(j) = 1.0e4 ! for neutron, photon, electron, positoron, up to 10 GeV
               elseif(ipcosmic(j).eq.29.or.ipcosmic(j).eq.30) then
                seg2(j) = 1.0e8 ! for muon, up to 100 TeV
               else
                seg2(j) = 1.0e6 ! up to 1 TeV
               endif
              else
               seg2(j)   = spava(21)
              endif

             ! T.Sato 2021/11/29, limitation of zenith angle is effective for all s-type
             if(ischn(71).eq.0) then
              sag1(j)=-1.0d0
             else
              sag1(j)=spava(71)
             endif
             if(ischn(72).eq.0) then
              sag2(j)=1.0d0
             else
              sag2(j)=spava(72)
             endif
             if(sag1(j).ge.sag2(j)) then
              goto 925
             endif

            isbias(j) = nint(spava(109)) ! T.Sato, 0:normal, 1:dead particle is not produced, 2020/03/03

              if( ischn(68) .eq. 0 ) then
               ischn(68) = 1
               spava(68) = -200

               ngrp(j)  = abs( nint( spava(68) ) )

               ngll(j)  = 1
               if( spava(68) .lt. 0.d0 ) ngll(j) = -1

               call moddas_reallocate_dbl(isrc, j, ngrp(j), ngei, egmin)
               call moddas_reallocate_dbl(isrc, j, ngrp(j), ngea, egmax)
               call moddas_reallocate_dbl(isrc, j, ngrp(j), ngfe, fegrp)
               call moddas_reallocate_dbl(isrc, j, ngrp(j), ngft, rfe)
               call moddas_reallocate_dbl(isrc, j, ngrp(j), ngpi, prw)
               call moddas_reallocate_dbl(isrc, j, ngrp(j), ngpw, pwt)

              endif

              eming = seg1(j)
              emaxg = seg2(j)

              if( eming .gt. emaxg ) then
               write(ErrCha,'("et1> et2 for e-type = 25 or 26.")')
               ErrID = 'L:19034/R:sours/F:read02.f' !E03_071_002
               call ErrWrite(ErrID,ErrCha)
               goto 997
              endif

              if( jetyp(j) .eq. 26 .and. ischn(63) .eq. 0 ) then

               do i = 1, ngrp(j)
                prw(ngpi(j)+i) = 1.0d0
               end do

              end if

**** No such e-type *************************
            else

             write(ErrCha,'("unknown e-type",i5)') jetyp(j)
             ErrID = 'L:19051/R:sours/F:read02.f' !E03_073_002
             call ErrWrite(ErrID,ErrCha)
               goto 997

            end if

         end if

*-----------------------------------------------------------------------
*     source angle group
*-----------------------------------------------------------------------

               jatyp(j)  = nint( spava(69) )

         if( jatyp(j) .eq.   1 .or.
     &       jatyp(j) .eq.   4 .or.
     &       jatyp(j) .eq.   5 .or.
     &       jatyp(j) .eq.   6 .or.
     &       jatyp(j) .eq.  11 .or.
     &       jatyp(j) .eq.  14 .or.
     &       jatyp(j) .eq.  15 .or.
     &       jatyp(j) .eq.  16 .or.
     &     ( jstyp(j) .eq.  17 .and. jatyp(j) .gt. 0 ) .or.
     &     ( jstyp(j) .eq. 100 .and. jatyp(j) .gt. 0 ) ) then

          if( sdir(j) .le. 250.0 .or. sdir(j) .gt. 350.0 ) then
           write(ErrCha,'("It is dir <= 250.0 or dir> 350.0.")')
           ErrID = 'L:19078/R:sours/F:read02.f' !E03_074_001
           call ErrWrite(ErrID,ErrCha)
           goto 965
          endif

*-----------------------------------------------------------------------

            if( jatyp(j) .eq. 1 .or. jatyp(j) .eq. 11 ) then

               if( ischn(70) .eq. 0 ) then
             write(ErrCha,'("a-type = 1 .or. a-type = 11 error")')
             ErrID = 'L:19089/R:sours/F:read02.f' !E03_075_001
             call ErrWrite(ErrID,ErrCha)
                goto 997
               endif

*-----------------------------------------------------------------------

            else if( jatyp(j) .eq. 5 .or. jatyp(j) .eq. 15 ) then

               if( ischn(73) .eq. 0 .or. ischn(71) .eq. 0 .or.
     &             ischn(72) .eq. 0 .or. ischn(74) .eq. 0 ) then
             write(ErrCha,*)
     & 'It is nn, ag1, ag2, g (x) error of ',
     & 'a-type = 5 or a-type = 15.'
             ErrID = 'L:19103/R:sours/F:read02.f' !E03_076_001
             call ErrWrite(ErrID,ErrCha)
               goto 997
             endif

                  sag1(j)   = spava(71)
                  sag2(j)   = spava(72)
                  isnn(j)   = spava(73)

                  aming = sag1(j)
                  amaxg = sag2(j)

                  if( aming .gt. amaxg ) then
             write(ErrCha,'("a-type = 5 or a-type = 15 ag1> ag2.")')
             ErrID = 'L:19117/R:sours/F:read02.f' !E03_076_002
             call ErrWrite(ErrID,ErrCha)
                    goto 997
                  endif

                  if( jatyp(j) .eq. 5 .and. aming .lt. -1.d0 ) then
             write(ErrCha,'("a-type = 5 or a-type = 15 ag1 <-1.")')
             ErrID = 'L:19124/R:sours/F:read02.f' !E03_076_003
             call ErrWrite(ErrID,ErrCha)
                    goto 997
                  endif

                  if( jatyp(j) .eq. 5 .and. amaxg .gt.  1.d0 ) then
             write(ErrCha,'("a-type = 5 or a-type = 15 ag2> 1.")')
             ErrID = 'L:19131/R:sours/F:read02.f' !E03_076_004
             call ErrWrite(ErrID,ErrCha)
                    goto 997
                  endif

*-----------------------------------------------------------------------

            else if( jatyp(j) .eq. 4 .or. jatyp(j) .eq. 14 ) then

               if( ischn(70) .eq. 0 ) then
             write(ErrCha,
     & '("a-type = 4 or a-type = 14 ag1 is not set.")')
             ErrID = 'L:19143/R:sours/F:read02.f' !E03_077_001
             call ErrWrite(ErrID,ErrCha)
                goto 997
              endif

               if( ischn(75) .eq. 0 ) then

                  do i = 1, narp(j)

                     paw(napi(j)+i) = 1.0d0

                  end do

               end if

*-----------------------------------------------------------------------

            else if( jatyp(j) .eq. 6 .or. jatyp(j) .eq. 16 ) then

               if( ischn(73) .eq. 0 .or. ischn(71) .eq. 0 .or.
     &             ischn(72) .eq. 0 .or. ischn(74) .eq. 0 ) then

             write(ErrCha,
     & '("Error of a-type = 6 or a-type = 16 nn, ag1, ag2, g (x).")')
             ErrID = 'L:19167/R:sours/F:read02.f' !E03_078_001
             call ErrWrite(ErrID,ErrCha)
               goto 997
             endif

                  sag1(j)   = spava(71)
                  sag2(j)   = spava(72)
                  isnn(j)   = spava(73)

               if( ischn(75) .eq. 0 ) then

                  do i = 1, narp(j)

                     paw(napi(j)+i) = 1.0d0

                  end do

               end if

                  aming = sag1(j)
                  amaxg = sag2(j)

                  if( aming .gt. amaxg ) then
             write(ErrCha,'("a-type = 6 or a-type = 16 ag1> ag2.")')
             ErrID = 'L:19191/R:sours/F:read02.f' !E03_078_002
             call ErrWrite(ErrID,ErrCha)
                    goto 997
                  endif

                  if( jatyp(j) .eq. 6 .and. aming .lt. -1.d0 ) then
             write(ErrCha,'("a-type = 6 or a-type = 16 ag1 <-1.")')
             ErrID = 'L:19198/R:sours/F:read02.f' !E03_078_003
             call ErrWrite(ErrID,ErrCha)
                    goto 997
                  endif

                  if( jatyp(j) .eq. 6 .and. amaxg .gt.  1.d0 ) then
             write(ErrCha,'("a-type = 6 or a-type = 16 ag2> 1.")')
             ErrID = 'L:19205/R:sours/F:read02.f' !E03_078_004
             call ErrWrite(ErrID,ErrCha)
                    goto 997
                  endif

*-----------------------------------------------------------------------

            else
             write(ErrCha,'("Error with a-type = 6 or a-type = 16.")')
             ErrID = 'L:19214/R:sours/F:read02.f' !E03_078_005
             call ErrWrite(ErrID,ErrCha)
               goto 997

            end if

         end if

*-----------------------------------------------------------------------
*     source time group
*-----------------------------------------------------------------------

               jttyp(j)  = nint( spava(29) )

         if( jttyp(j) .lt. 0 .or.
     &     ( jttyp(j) .gt. 6 .and. jttyp(j) .ne. 100 ) ) then
             write(ErrCha,'("It is an error of tw.")')
             ErrID = 'L:19231/R:sours/F:read02.f' !E03_079_001
             call ErrWrite(ErrID,ErrCha)
         goto 987
       endif

         if( jttyp(j) .eq.   1 .or.
     &       jttyp(j) .eq.   2 .or.
     &       jttyp(j) .eq.   3 .or.
     &       jttyp(j) .eq.   4 .or.
     &       jttyp(j) .eq.   5 .or.
     &       jttyp(j) .eq.   6 .or.
     &       jttyp(j) .eq. 100 .or.
     &     ( jstyp(j) .eq.  17 .and. jttyp(j) .gt. 0 ) .or.
     &     ( jstyp(j) .eq. 100 .and. jttyp(j) .gt. 0 ) ) then

*-----------------------------------------------------------------------

            if( jttyp(j) .eq. 1 .or. jttyp(j) .eq. 2 ) then

                  stm0(j) = spava(30)

               if( ischn(31) .eq. 0 ) then
             write(ErrCha,'("It is an error of tc.")')
             ErrID = 'L:19254/R:sours/F:read02.f' !E03_080_001
             call ErrWrite(ErrID,ErrCha)
                goto 990
               endif

                  stmw(j) = spava(31)

               if( ischn(32) .eq. 0 .and. jttyp(j) .eq. 2 ) then

                  stmc(j) = stmw(j) * 10.0

               else if( ischn(32) .ne. 0 ) then

                  stmc(j) = spava(32)

               end if

                  jttpn(j) = nint( spava(33) )

               if( jttpn(j) .lt. 1 ) then
             write(ErrCha,'("It is an error of tn.")')
             ErrID = 'L:19275/R:sours/F:read02.f' !E03_081_001
             call ErrWrite(ErrID,ErrCha)
                goto 989
               endif

               if( jttpn(j) .gt. 1 .and. ischn(34) .eq. 0 ) then
             write(ErrCha,'("It is an error of td.")')
             ErrID = 'L:19282/R:sours/F:read02.f' !E03_082_001
             call ErrWrite(ErrID,ErrCha)
                goto 988
               endif

                  stmd(j) = spava(34)

*-----------------------------------------------------------------------

            else if( jttyp(j) .eq. 3 ) then

               if( ischn(79) .eq. 0 ) then
             write(ErrCha,'("It is an error of tg1.")')
             ErrID = 'L:19295/R:sours/F:read02.f' !E03_083_001
             call ErrWrite(ErrID,ErrCha)
                goto 997
               endif

*-----------------------------------------------------------------------

            else if( jttyp(j) .eq. 5 ) then

               if( ischn(82) .eq. 0 .or. ischn(80) .eq. 0 .or.
     &             ischn(81) .eq. 0 .or. ischn(83) .eq. 0 ) then
             write(ErrCha,'("tg2, ll, h (x), o-type error.")')
             ErrID = 'L:19307/R:sours/F:read02.f' !E03_084_001
             call ErrWrite(ErrID,ErrCha)
               goto 997
             endif

                  stg1(j)   = spava(80)
                  stg2(j)   = spava(81)
                  isll(j)   = spava(82)

                  tming = stg1(j)
                  tmaxg = stg2(j)

                  if( tming .gt. tmaxg ) then
             write(ErrCha,'("Error of tg1> tg2.")')
             ErrID = 'L:19321/R:sours/F:read02.f' !E03_085_001
             call ErrWrite(ErrID,ErrCha)
                    goto 997
                  endif

*-----------------------------------------------------------------------

            else if( jttyp(j) .eq. 4 ) then

               if( ischn(79) .eq. 0 ) then
             write(ErrCha,'("Error of tg1.")')
             ErrID = 'L:19332/R:sours/F:read02.f' !E03_086_001
             call ErrWrite(ErrID,ErrCha)
                goto 997
               endif

               if( ischn(84) .eq. 0 ) then

                  do i = 1, ntrp(j)

                     ptw(ntpi(j)+i) = 1.0d0

                  end do

               end if

*-----------------------------------------------------------------------

            else if( jatyp(j) .eq. 6 ) then

               if( ischn(82) .eq. 0 .or. ischn(80) .eq. 0 .or.
     &             ischn(81) .eq. 0 .or. ischn(83) .eq. 0 ) then
             write(ErrCha,'("tg2, ll, h (x), o-type error.")')
             ErrID = 'L:19354/R:sours/F:read02.f' !E03_087_001
             call ErrWrite(ErrID,ErrCha)
               goto 997
             endif

                  stg1(j) = spava(80)
                  stg2(j) = spava(81)
                  isll(j) = spava(82)

               if( ischn(84) .eq. 0 ) then

                  do i = 1, ntrp(j)

                     ptw(ntpi(j)+i) = 1.0d0

                  end do

               end if

                  tming = stg1(j)
                  tmaxg = stg2(j)

                  if( tming .gt. tmaxg ) then
             write(ErrCha,'("Error of tg1> tg2.")')
             ErrID = 'L:19378/R:sours/F:read02.f' !E03_088_001
             call ErrWrite(ErrID,ErrCha)
                    goto 997
                  endif

*-----------------------------------------------------------------------

            else if( jttyp(j) .eq. 100 ) then

               if( ischn(80) .eq. 0 .or. ischn(81) .eq. 0 ) then
             write(ErrCha,
     & '("t-type 100 tg1 = 0 or tg2 = 0 error error.")')
             ErrID = 'L:19390/R:sours/F:read02.f' !E03_089_001
             call ErrWrite(ErrID,ErrCha)
               goto 997
               endif

                  stg1(j) = spava(80)
                  stg2(j) = spava(81)

                  tming = stg1(j)
                  tmaxg = stg2(j)

                  if( tming .gt. tmaxg ) then
             write(ErrCha,'("Error of tg1> tg2.")')
             ErrID = 'L:19403/R:sours/F:read02.f' !E03_090_001
             call ErrWrite(ErrID,ErrCha)
                    goto 997
                  endif

*-----------------------------------------------------------------------

            else
             write(ErrCha,'("Error of t-type 100.")')
             ErrID = 'L:19412/R:sours/F:read02.f' !E03_091_001
             call ErrWrite(ErrID,ErrCha)
               goto 997

            end if

         end if

*-----------------------------------------------------------------------
*     special for duct
*-----------------------------------------------------------------------

         if( ( jstyp(j) .eq. 1 .or.
     &         jstyp(j) .eq. 4 .or.
     &         jstyp(j) .eq. 2 .or.
     &         jstyp(j) .eq. 5 ) .and.
     &       ( sdom(j) .gt. -10.5 .and. sdom(j) .lt. -9.5 ) ) then

               if( ischn(47) .eq. 0 .or. ischn(48) .eq. 0 .or.
     &             ischn(53) .eq. 0 ) goto 980

               if(   ischn(52) .eq. 0  .and.
     &             ( jstyp(j) .eq. 1 .or. jstyp(j) .eq. 4 ) ) then
             write(ErrCha,'("special for duct error")')
             ErrID = 'L:19436/R:sours/F:read02.f' !E03_092_001
             call ErrWrite(ErrID,ErrCha)
               goto 979
               endif

               if( ( ischn(50) .eq. 0 .or. ischn(51) .eq. 0 ) .and.
     &             ( jstyp(j) .eq. 2 .or. jstyp(j) .eq. 5 ) ) then
             write(ErrCha,'("special for duct error")')
             ErrID = 'L:19444/R:sours/F:read02.f' !E03_093_001
             call ErrWrite(ErrID,ErrCha)
               goto 978
               endif

               sdir(j) = 1.0d0
               sz1(j)  = sz0(j)

               sdl0(j) = spava(53)
               sdl1(j) = spava(47)
               sdl2(j) = spava(48)

               if( sdl0(j) .le. 0.0d0 .or. sdl1(j) .le. 0.0d0 .or.
     &             sdl2(j) .le. 0.0d0 .or. sdl0(j) .gt. sdl1(j) .or.
     &             sdl1(j) .gt. sdl2(j) ) then
             write(ErrCha,*)
     & 'It is an error of dl0 <= 0 or dl1 <= 0 or ',
     & 'dl2 <= 0 or dl0> dl1 or dl1> dl2.'
             ErrID = 'L:19462/R:sours/F:read02.f' !E03_094_001
             call ErrWrite(ErrID,ErrCha)
               goto 996
               endif

            if( ischn(49) .ne. 0 ) then

               sdpf(j) = spava(49)

               if( sdpf(j) .ge. 1.0d0 .or. sdpf(j) .le. 0.0d0 ) then
             write(ErrCha,'("special for duct error")')
             ErrID = 'L:19473/R:sours/F:read02.f' !E03_095_001
             call ErrWrite(ErrID,ErrCha)
               goto 977
               endif

            else

               sdpf(j) = 0.2d0

            end if

            if( jstyp(j) .eq. 1 .or. jstyp(j) .eq. 4 ) then

               sdrd(j) = spava(52)

            else if( jstyp(j) .eq. 2 .or. jstyp(j) .eq. 5 ) then

               sdxw(j) = spava(50)
               sdyw(j) = spava(51)

               sdrd(j) = sqrt( sdxw(j)**2 + sdyw(j)**2 ) / 2.0d0

            end if

            if( ischn(54) .ne. 0 ) then

               sdebg(j) = 1.0d0

            else

               sdebg(j) = 0.0d0

            end if

*-----------------------------------------------------------------------

               sdsxp(j) = 1.0d0
               sdsxn(j) = 1.0d0
               sdsyp(j) = 1.0d0
               sdsyn(j) = 1.0d0

            if( ischn(55) .ne. 0 ) sdsxp(j) = spava(55)
            if( ischn(56) .ne. 0 ) sdsxn(j) = spava(56)
            if( ischn(57) .ne. 0 ) sdsyp(j) = spava(57)
            if( ischn(58) .ne. 0 ) sdsyn(j) = spava(58)

            if( sdsxp(j) .le. 0.0d0 ) then
             write(ErrCha,'("special for duct error")')
             ErrID = 'L:19521/R:sours/F:read02.f' !E03_096_001
             call ErrWrite(ErrID,ErrCha)
             goto 975
            endif
            if( sdsxn(j) .le. 0.0d0 ) then
             write(ErrCha,'("special for duct error")')
             ErrID = 'L:19527/R:sours/F:read02.f' !E03_096_002
             call ErrWrite(ErrID,ErrCha)
             goto 975
            endif
            if( sdsyp(j) .le. 0.0d0 ) then
             write(ErrCha,'("special for duct error")')
             ErrID = 'L:19533/R:sours/F:read02.f' !E03_096_003
             call ErrWrite(ErrID,ErrCha)
             goto 975
            endif
            if( sdsyn(j) .le. 0.0d0 ) then
             write(ErrCha,'("special for duct error")')
             ErrID = 'L:19539/R:sours/F:read02.f' !E03_096_004
             call ErrWrite(ErrID,ErrCha)
             goto 975
            endif

            if( jstyp(j) .eq. 2 .or. jstyp(j) .eq. 5 ) then

               prs1(j) = 1.d0 / sdsxp(j)
               prs2(j) = 1.d0 / sdsyp(j)
               prs3(j) = 1.d0 / sdsxn(j)
               prs4(j) = 1.d0 / sdsyn(j)

               psum = 4.d0 / ( prs1(j) + prs2(j) + prs3(j) + prs4(j) )

               prs1(j) = prs1(j) * psum
               prs2(j) = prs2(j) * psum
               prs3(j) = prs3(j) * psum
               prs4(j) = prs4(j) * psum

               totslg =  sdyw(j) * sdsxp(j) + sdxw(j) * sdsyp(j)
     &                +  sdyw(j) * sdsxn(j) + sdxw(j) * sdsyn(j)

               psxp(j) = ( sdyw(j) * sdsxp(j) ) / totslg
               psyp(j) = ( sdxw(j) * sdsyp(j) ) / totslg
               psxn(j) = ( sdyw(j) * sdsxn(j) ) / totslg
               psyn(j) = ( sdxw(j) * sdsyn(j) ) / totslg

               pnorm = prs1(j) * psxp(j) + prs2(j) * psyp(j)
     &               + prs3(j) * psxn(j) + prs4(j) * psyn(j)

               prs1(j) = prs1(j) / pnorm
               prs2(j) = prs2(j) / pnorm
               prs3(j) = prs3(j) / pnorm
               prs4(j) = prs4(j) / pnorm

               psxp(j) = psxp(j)
               psyp(j) = psxp(j) + psyp(j)
               psxn(j) = psyp(j) + psxn(j)
               psyn(j) = psxn(j) + psyn(j)

            end if

*-----------------------------------------------------------------------

               sdls(j) = -1.0d0
               sdrs(j) = -1.0d0
               sdxs(j) = -1.0d0
               sdys(j) = -1.0d0

            if( ischn(59) .ne. 0 ) sdls(j) = spava(59)
            if( ischn(60) .ne. 0 ) sdrs(j) = spava(60)
            if( ischn(61) .ne. 0 ) sdxs(j) = spava(61)
            if( ischn(62) .ne. 0 ) sdys(j) = spava(62)

            if( sdls(j) .gt. 0.0d0 ) then

               if( sdls(j) .gt. sdl0(j) ) then
             write(ErrCha,'("(drs <= 0 and dxs <= 0) or dys <= 0.")')
             ErrID = 'L:19597/R:sours/F:read02.f' !E03_097_001
             call ErrWrite(ErrID,ErrCha)
                goto 974
               endif

                  if( sdrs(j) .le. 0.0d0 .and.
     &              ( sdxs(j) .le. 0.0d0 .or.
     &                sdys(j) .le. 0.0d0 ) ) then
             write(ErrCha,'("(drs> 0 and dxs> 0) or dys> 0 error.")')
             ErrID = 'L:19606/R:sours/F:read02.f' !E03_098_001
             call ErrWrite(ErrID,ErrCha)
                  goto 973
                endif

                  if( sdrs(j) .gt. 0.0d0 .and.
     &              ( sdxs(j) .gt. 0.0d0 .or.
     &                sdys(j) .gt. 0.0d0 ) ) then
             write(ErrCha,'("special for duct error")')
             ErrID = 'L:19615/R:sours/F:read02.f' !E03_099_001
             call ErrWrite(ErrID,ErrCha)
             goto 973
             endif

            end if

*-----------------------------------------------------------------------

            if( nglp(j) .gt. 0 ) then

               if( slmin(ngli(j)+1) .ne. sdl1(j) ) then
               write(ErrCha,'("special for duct error")')
               ErrID = 'L:19628/R:sours/F:read02.f' !E03_100_001
               call ErrWrite(ErrID,ErrCha)
               goto 970
               endif
               if( slmax(ngla(j)+nglp(j)) .ne. sdl2(j) ) then
               write(ErrCha,'("special for duct error")')
               ErrID = 'L:19634/R:sours/F:read02.f' !E03_101_001
               call ErrWrite(ErrID,ErrCha)
                goto 969
               endif

                     srpe = 0.0

               do ig = 1, nglp(j)

                  if( flgrp(ngfl(j)+ig) .gt. 0.0d0 ) then

                     rgw(nglw(j)+ig) = 1.d0 / flgrp(ngfl(j)+ig)
                     srpe = srpe + rgw(nglw(j)+ig)

                  else

                     rgw(nglw(j)+ig) = 0.d0

                  end if

               end do

               do ig = 1, nglp(j)

                     rgw(nglw(j)+ig) = rgw(nglw(j)+ig)
     &                               / srpe * nglp(j)

               end do

                     totol = 0.0

               do ig = 1, nglp(j)

                     totol = totol
     &                     + ( slmax(ngla(j)+ig)
     &                     - slmin(ngli(j)+ig) )
     &                     * flgrp(ngfl(j)+ig)

               end do

               do ig = 1, nglp(j)

                     rgl(nglc(j)+ig) = ( slmax(ngla(j)+ig)
     &                                 - slmin(ngli(j)+ig) )
     &                                 * flgrp(ngfl(j)+ig) / totol

               end do

                     pnorm = 0.0

               do ig = 1, nglp(j)

                     pnorm = pnorm + rgw(nglw(j)+ig)
     &                             * rgl(nglc(j)+ig)

               end do

               do ig = 1, nglp(j)

                     rgw(nglw(j)+ig) = rgw(nglw(j)+ig) / pnorm

               end do

               do ig = 2, nglp(j)

                     rgl(nglc(j)+ig) = rgl(nglc(j)+ig-1)
     &                               + rgl(nglc(j)+ig)

               end do

            end if

*-----------------------------------------------------------------------

               dnorm(j)  = ( sdl0(j) / sdl1(j) )**2

*-----------------------------------------------------------------------

         else if( sdom(j) .gt. -10.5 .and. sdom(j) .lt. -9.5 ) then
            write(ErrCha,'("special for duct error")')
            ErrID = 'L:19714/R:sours/F:read02.f' !E03_102_001
            call ErrWrite(ErrID,ErrCha)
            goto 966

         end if

*-----------------------------------------------------------------------
*        energy mesh preparation
*-----------------------------------------------------------------------

      if( jetyp(j) .eq.  3 .or. jetyp(j) .eq.  7 .or.
     &    jetyp(j) .eq. 25 .or. jetyp(j) .eq. 26 .or. ! T.Sato 2020/12/06
     &    jetyp(j) .eq.  5 .or. jetyp(j) .eq. 15 .or.
     &    jetyp(j) .eq.  6 .or. jetyp(j) .eq. 16 ) then

         if( jetyp(j) .ne. 3 .and. jetyp(j) .ne. 7 ) then

            do i = 1, mxcval

               cval0(i) = cval(i)
               cval(i) = cval1(i)

            end do

         end if

         if( ngll(j) .gt. 0 ) then

               edif = ( emaxg - eming ) / dble( ngrp(j) )

         else

               edif = log( emaxg / eming ) / dble( ngrp(j) )

         end if

         if( ngll(j) .gt. 0 ) then

            do ig = 1, ngrp(j)

               egmin(ngei(j)+ig) = eming + dble( ig - 1 ) * edif
               egmax(ngea(j)+ig) = eming + dble( ig ) * edif

               eval = eming +  ( ig - 0.5d0  ) * edif

               if( jetyp(j) .eq. 3 .or. jetyp(j) .eq. 7 ) then

cABE 2017/08/17, further revised by T.Sato 2020/02/05
                  fegrp(ngfe(j)+ig) = eval**set3(j)
     &                              * exp( -eval / set0(j) )

               elseif( jetyp(j) .eq. 25 .or. jetyp(j) .eq. 26 ) then ! T.Sato 2020/12/06

                  fegrp(ngfe(j)+ig) =
     &            getCosSpec(icenv(j),ipcosmic(j),solarmod(j),rigid(j),
     &            depatom(j),eval,ground(j))

       if(ipcosmic(j).eq.33) then ! T.Sato 2023/07/23 annihilation gamma correction
        if(egmin(ngei(j)+ig).lt.rstms(12)*1.0d3.and.
     &     egmax(ngei(j)+ig).ge.rstms(12)*1.0d3) then
         ie511(j)=ig
         flux511bin=fegrp(ngfe(j)+ig)
     &   *(egmax(ngei(j)+ig)-egmin(ngei(j)+ig))
         annihratio(j)=(flux511bin
     &   +get511flux(solarmod(j),rigid(j),depatom(j)))/flux511bin
         fegrp(ngfe(j)+ig) = fegrp(ngfe(j)+ig) * annihratio(j)
        endif
       endif

               else

                  call gnum(chlg,1,lsfx(j),cvvv,eval,ierr)

                  fegrp(ngfe(j)+ig) = cvvv

               end if

! T.Sato 2021/12/01, fegrp should be integral fluence in the bin in the same as ngll < 0
               fegrp(ngfe(j)+ig) = fegrp(ngfe(j)+ig)
     &         * ( egmax(ngea(j)+ig) - egmin(ngei(j)+ig) )

            end do

         else

            do ig = 1, ngrp(j)

               egmin(ngei(j)+ig) = eming
     &                           * exp( dble( ig - 1 ) * edif )
               egmax(ngea(j)+ig) = eming
     &                           * exp( dble( ig ) * edif )

               eval = eming * exp( ( ig - 0.5d0  ) * edif )

               if( jetyp(j) .eq. 3 .or. jetyp(j) .eq. 7 ) then

cABE 2017/08/17, further revised by T.Sato 2020/02/05
                  fegrp(ngfe(j)+ig) = eval**set3(j)
     &                              * exp( -eval / set0(j) )

               elseif( jetyp(j) .eq. 25 .or. jetyp(j) .eq. 26 ) then  ! T.Sato 2020/12/06

                  fegrp(ngfe(j)+ig) =
     &            getCosSpec(icenv(j),ipcosmic(j),solarmod(j),rigid(j),
     &            depatom(j),eval,ground(j))

       if(ipcosmic(j).eq.33) then ! T.Sato 2023/07/23 annihilation gamma correction
        if(egmin(ngei(j)+ig).lt.rstms(12)*1.0d3.and.
     &     egmax(ngei(j)+ig).ge.rstms(12)*1.0d3) then
         ie511(j)=ig
         flux511bin=fegrp(ngfe(j)+ig)
     &   *(egmax(ngei(j)+ig)-egmin(ngei(j)+ig))
         annihratio(j)=(flux511bin
     &   +get511flux(solarmod(j),rigid(j),depatom(j)))/flux511bin
         fegrp(ngfe(j)+ig) = fegrp(ngfe(j)+ig) * annihratio(j)
        endif
       endif

               else

                  call gnum(chlg,1,lsfx(j),cvvv,eval,ierr)

                  fegrp(ngfe(j)+ig) = cvvv

               end if

               fegrp(ngfe(j)+ig) = fegrp(ngfe(j)+ig)
     &         * ( egmax(ngea(j)+ig) - egmin(ngei(j)+ig) )

            end do

         end if

! cosmic-ray ground level correction
         if( (jetyp(j).eq.25.or.jetyp(j).eq.26).and.icenv(j).ge.1) then  ! T.Sato 2021/12/01
          call CosmicAngSetup(j) ! set up angular table for terrestrial mode
          if((icenv(j).eq.2.and.abs(inkf0(j)).eq.13).or. ! ground level muon
     &    (icenv(j).eq.5.and.inkf0(j).eq.2112)) then          ! neutron for black-hole mode
           do ig = 1, ngrp(j)
            fegrp(ngfe(j)+ig) = fegrp(ngfe(j)+ig)* angtable(-1,ig,j)
           enddo
          endif
         endif

         if( jetyp(j) .ne. 3 .and. jetyp(j) .ne. 7 ) then

            do i = 1, mxcval

               cval(i) = cval0(i)

            end do

         end if

      end if

*-----------------------------------------------------------------------

      if( jetyp(j) .eq.  3 .or. jetyp(j) .eq. 25 .or.  ! T.Sato 2020/12/06
     &    jetyp(j) .eq.  5 .or. jetyp(j) .eq. 15 ) then

               do i = 1, ngrp(j)

                  prw(ngpi(j)+i) = fegrp(ngfe(j)+i)

               end do

      end if

*-----------------------------------------------------------------------
      if( jetyp(j) .eq.  1 .or. jetyp(j) .eq. 11 .or.
     &    jetyp(j) .eq.  4 .or. jetyp(j) .eq. 14 .or.
     &    jetyp(j) .eq. 21 .or. jetyp(j) .eq. 31 .or.
     &    jetyp(j) .eq. 24 .or. jetyp(j) .eq. 34 .or.
     &    jetyp(j) .eq.  8 .or. jetyp(j) .eq. 18 .or.
     &    jetyp(j) .eq.  9 .or. jetyp(j) .eq. 19 .or.
     &    jetyp(j) .eq.  5 .or. jetyp(j) .eq. 15 .or.
     &    jetyp(j) .eq.  6 .or. jetyp(j) .eq. 16 .or.
     &    jetyp(j) .eq.  3 .or. jetyp(j) .eq.  7 .or.
     &    jetyp(j) .eq. 22 .or. jetyp(j) .eq. 32 .or.
     &    jetyp(j) .eq. 23 .or. jetyp(j) .eq. 33 .or.
c
     &    jetyp(j) .eq. 20 .or.                      ! S.H. (2016.12.30)
     &    jetyp(j) .eq. 25 .or. jetyp(j) .eq. 26 .or. ! T.Sato 2020/12/06
     &    jetyp(j) .eq. 28 .or. jetyp(j) .eq. 29 ) then

*-----------------------------------------------------------------------

               srfe = 0.0

            do ig = 1, ngrp(j)

               srfe = srfe + fegrp(ngfe(j)+ig)

            end do

c            if(srfe.eq.0) goto 927
            if(srfe.eq.0) then
                imsrc = imsrc - 1
                write(*,*) 'Warning. Integral fluence of one multisource
     &  is zero. This source component is skipped'
                goto 6000
            endif

            if(jetyp(j) .eq. 25 .or. jetyp(j) .eq. 26) then ! T.Sato 2020/12/06
             normfact(j)=srfe ! for cosmic-ray source, absolute flux (/cm2/s) is important
            endif

            do ig = 1, ngrp(j)

               fegrp(ngfe(j)+ig) = fegrp(ngfe(j)+ig) / srfe

            end do

*-----------------------------------------------------------------------

               srpp = 0.0

            do ig = 1, ngrp(j)

               srpp = srpp + prw(ngpi(j)+ig)

            end do

            do ig = 1, ngrp(j)

               prw(ngpi(j)+ig) = prw(ngpi(j)+ig) / srpp

            end do

*-----------------------------------------------------------------------

       if ( jetyp(j) .ne. 20 ) then ! S.H. (2016.12.30)

                  srpe = 0.0
                  ngrx = 0

            do ig = 1, ngrp(j)

               if( prw(ngpi(j)+ig) .gt. 0.d0 ) then

                  pwt(ngpw(j)+ig) = fegrp(ngfe(j)+ig)
     &                              / prw(ngpi(j)+ig)
                  srpe = srpe + pwt(ngpw(j)+ig)

                  ngrx = ngrx + 1

               else

                  pwt(ngpw(j)+ig) = 0.0d0

               end if

            end do

            do ig = 1, ngrp(j)

               pwt(ngpw(j)+ig) = pwt(ngpw(j)+ig) / srpe * ngrx

            end do

       else ! S.H. (2016.12.30)
        do ig = 1, ngrp(j)
           pwt(ngpw(j)+ig) = 1d0 ! weight=1.0 for all energy bins in e-type=20
        end do

       end if

*-----------------------------------------------------------------------
cKN 2013/08/17 Bug was fixed, no affect on the result, except for echo.

               rfe(ngft(j)+1) = prw(ngpi(j)+1)

            do ig = 2, ngrp(j)

               rfe(ngft(j)+ig) = rfe(ngft(j)+ig-1)
     &                         + prw(ngpi(j)+ig)

            end do

*-----------------------------------------------------------------------
cKN 2013/08/17, differential type

         if( jetyp(j) .eq. 21 .or. jetyp(j) .eq. 31 .or.
     &       jetyp(j) .eq. 24 .or. jetyp(j) .eq. 34 ) then

            do ig = 1, ngrp(j)

               fegrp(ngfe(j)+ig) = fegrp(ngfe(j)+ig) * srfe
     &                     / ( egmax(ngea(j)+ig) - egmin(ngei(j)+ig) )

            end do

         else

            do ig = 1, ngrp(j)

               fegrp(ngfe(j)+ig) = fegrp(ngfe(j)+ig) * srfe

            end do

         end if

            do ig = 1, ngrp(j)

               prw(ngpi(j)+ig) = prw(ngpi(j)+ig) * srpp

            end do

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*        angle mesh preparation
*-----------------------------------------------------------------------

      if( jatyp(j) .eq.  5 .or. jatyp(j) .eq. 15 .or.
     &    jatyp(j) .eq.  6 .or. jatyp(j) .eq. 16 ) then

            do i = 1, mxcval

               cval0(i) = cval(i)
               cval(i) = cval2(i)

            end do

         if( nall(j) .gt. 0 ) then

               adif = ( amaxg - aming ) / dble( narp(j) )

         else

               adif = log( amaxg / aming ) / dble( narp(j) )

         end if

         if( nall(j) .gt. 0 ) then

            do ig = 1, narp(j)

               agmin(naei(j)+ig) = aming + dble( ig - 1 ) * adif
               agmax(naea(j)+ig) = aming + dble( ig ) * adif

               angl = aming +  ( ig - 0.5d0  ) * adif

                  call gnum(chla,1,lsfy(j),cvvv,angl,ierr)

                  fagrp(nafe(j)+ig) = cvvv

               if( jatyp(j) .gt. 10 ) then

                  fagrp(nafe(j)+ig) = fagrp(nafe(j)+ig)
     &                              * sin( angl / 180.0 * pi )

               end if

            end do

         else

            do ig = 1, narp(j)

               agmin(naei(j)+ig) = aming
     &                           * exp( dble( ig - 1 ) * adif )
               agmax(naea(j)+ig) = aming
     &                           * exp( dble( ig ) * adif )

               angl = aming * exp( ( ig - 0.5d0  ) * adif )

                  call gnum(chla,1,lsfy(j),cvvv,angl,ierr)

                  fagrp(nafe(j)+ig) = cvvv

               if( jatyp(j) .gt. 10 ) then

                  fagrp(nafe(j)+ig) = fagrp(nafe(j)+ig)
     &                              * sin( angl / 180.0 * pi )

               end if

                  fagrp(nafe(j)+ig) = fagrp(nafe(j)+ig)
     &            * ( agmax(naea(j)+ig) - agmin(naei(j)+ig) )

            end do

         end if

            do i = 1, mxcval

               cval(i) = cval0(i)

            end do

      end if

*-----------------------------------------------------------------------

      if( jatyp(j) .eq.  5 .or. jatyp(j) .eq. 15 ) then

               do i = 1, narp(j)

                  paw(napi(j)+i) = fagrp(nafe(j)+i)

               end do

      end if

*-----------------------------------------------------------------------

      if( jatyp(j) .eq.  1 .or. jatyp(j) .eq. 11 .or.
     &    jatyp(j) .eq.  4 .or. jatyp(j) .eq. 14 .or.
     &    jatyp(j) .eq.  5 .or. jatyp(j) .eq. 15 .or.
     &    jatyp(j) .eq.  6 .or. jatyp(j) .eq. 16 ) then

*-----------------------------------------------------------------------

               srfe = 0.0

            do ig = 1, narp(j)

               srfe = srfe + fagrp(nafe(j)+ig)

            end do

            if(srfe.eq.0) goto 927

            do ig = 1, narp(j)

               fagrp(nafe(j)+ig) = fagrp(nafe(j)+ig) / srfe

            end do

               srpe = 0.0

            do ig = 1, narp(j)

               srpe = srpe + paw(napi(j)+ig)

            end do

            do ig = 1, narp(j)

               paw(napi(j)+ig) = paw(napi(j)+ig) / srpe

            end do

                  srpe = 0.0
                  ngrx = 0

            do ig = 1, narp(j)

               if( paw(napi(j)+ig) .gt. 0.d0 ) then

                  pat(napw(j)+ig) = fagrp(nafe(j)+ig)
     &                              / paw(napi(j)+ig)
                  srpe = srpe + pat(napw(j)+ig)

                  ngrx = ngrx + 1

               else

                  pat(napw(j)+ig) = 0.0d0

               end if

            end do

            do ig = 1, narp(j)

               pat(napw(j)+ig) = pat(napw(j)+ig) / srpe * ngrx

            end do

            do ig = 1, narp(j)

               fagrp(nafe(j)+ig) = paw(napi(j)+ig)

            end do

               rfa(naft(j)+1) = fagrp(nafe(j)+1)

            do ig = 2, narp(j)

               rfa(naft(j)+ig) = rfa(naft(j)+ig-1)
     &                         + fagrp(nafe(j)+ig)

            end do

      end if

*-----------------------------------------------------------------------
*        time mesh preparation
*-----------------------------------------------------------------------

      if( jttyp(j) .eq. 5 .or. jttyp(j) .eq. 6 ) then

            do i = 1, mxcval

               cval0(i) = cval(i)
               cval(i) = cval3(i)

            end do

         if( ntll(j) .gt. 0 ) then

               tdif = ( tmaxg - tming ) / dble( ntrp(j) )

         else

               tdif = log( tmaxg / tming ) / dble( ntrp(j) )

         end if

         if( ntll(j) .gt. 0 ) then

            do ig = 1, ntrp(j)

               tgmin(ntei(j)+ig) = tming + dble( ig - 1 ) * tdif
               tgmax(ntea(j)+ig) = tming + dble( ig ) * tdif

               tngl = tming +  ( ig - 0.5d0  ) * tdif

                  call gnum(chlt,1,lsfz(j),cvvv,tngl,ierr)

                  ftgrp(ntfe(j)+ig) = cvvv

            end do

         else

            do ig = 1, ntrp(j)

               tgmin(ntei(j)+ig) = tming
     &                           * exp( dble( ig - 1 ) * tdif )
               tgmax(ntea(j)+ig) = tming
     &                           * exp( dble( ig ) * tdif )

               tngl = tming * exp( ( ig - 0.5d0  ) * tdif )

                  call gnum(chlt,1,lsfz(j),cvvv,tngl,ierr)

                  ftgrp(ntfe(j)+ig) = cvvv

                  ftgrp(ntfe(j)+ig) = ftgrp(ntfe(j)+ig)
     &            * ( tgmax(ntea(j)+ig) - tgmin(ntei(j)+ig) )

            end do

         end if

            do i = 1, mxcval

               cval(i) = cval0(i)

            end do

      end if

*-----------------------------------------------------------------------

      if( jttyp(j) .eq.  5 ) then

               do i = 1, ntrp(j)

                  ptw(ntpi(j)+i) = ftgrp(ntfe(j)+i)

               end do

      end if

*-----------------------------------------------------------------------

      if( jttyp(j) .eq.  3 .or. jttyp(j) .eq.  4 .or.
     &    jttyp(j) .eq.  5 .or. jttyp(j) .eq.  6 ) then

*-----------------------------------------------------------------------

               srfe = 0.0

            do ig = 1, ntrp(j)

               srfe = srfe + ftgrp(ntfe(j)+ig)

            end do

            if(srfe.eq.0) goto 927

            do ig = 1, ntrp(j)

               ftgrp(ntfe(j)+ig) = ftgrp(ntfe(j)+ig) / srfe

            end do

               srpe = 0.0

            do ig = 1, ntrp(j)

               srpe = srpe + ptw(ntpi(j)+ig)

            end do

            do ig = 1, ntrp(j)

               ptw(ntpi(j)+ig) = ptw(ntpi(j)+ig) / srpe

            end do

                  srpe = 0.0
                  ngrx = 0

            do ig = 1, ntrp(j)

               if( ptw(ntpi(j)+ig) .gt. 0.d0 ) then

                  ptt(ntpw(j)+ig) = ftgrp(ntfe(j)+ig)
     &                              / ptw(ntpi(j)+ig)
                  srpe = srpe + ptt(ntpw(j)+ig)

                  ngrx = ngrx + 1

               else

                  ptt(ntpw(j)+ig) = 0.0d0

               end if

            end do

            do ig = 1, ntrp(j)

               ptt(ntpw(j)+ig) = ptt(ntpw(j)+ig) / srpe * ngrx

            end do

            do ig = 1, ntrp(j)

               ftgrp(ntfe(j)+ig) = ptw(ntpi(j)+ig)

            end do

               rft(ntft(j)+1) = ftgrp(ntfe(j)+1)

            do ig = 2, ntrp(j)

               rft(ntft(j)+ig) = rft(ntft(j)+ig-1)
     &                         + ftgrp(ntfe(j)+ig)

            end do

      end if

*-----------------------------------------------------------------------
*        multi-source:  <source> subsection
*-----------------------------------------------------------------------

         if( ispfs(j) .ne. 0 ) ispfn = 1

 6000    if( imlcnt .eq. 1 ) goto 5000

*-----------------------------------------------------------------------
*        multi correlation source
*-----------------------------------------------------------------------

         if( iscorr .ne. 0 .and. imsrc .eq. 1 ) then
          write(ErrCha,'("multi correlation source error")')
          ErrID = 'L:20380/R:sours/F:read02.f' !E03_103_001
          call ErrWrite(ErrID,ErrCha)
          goto 963
         endif

         if( iscorr .ne. 0 .and. imsrc .gt. 1 ) then

               isek = 0

            do k = 1, imsrc

               imlwt(k) = nint( smlwt(k) )
               isek = isek + imlwt(k)

            end do

               itcorr = isek

         end if

*-----------------------------------------------------------------------
! T. Sato 2021/09/01, setup memory for cosmic-ray
      if(jetyp(1).eq.25.or.jetyp(1).eq.26) then ! cosmic-ray source mode, update smlwt and totfact here
       totfact2=totfact ! original totfact
       tmp=0.0
       do j=1,imsrc
        smlwt2(j)=smlwt(j)            ! original <source>, should be 1
        smlwt(j)=smlwt(j)*normfact(j) ! updated <source>: +- total flux (/cm2/s)
        tmp=tmp+normfact(j)
       enddo
       totfact=totfact*tmp  ! (/s) when original totfact is source area (cm2)
      endif
*-----------------------------------------------------------------------

 8888 continue

      if( iprojall .eq. -1 ) then
         write(niws,'(a)') chlw(1:i2)
         close(niws)   ! close & delete, unit niws is scratch file

         if ( npcunt .gt. 0 ) then
           imsrc = 0     ! reset multi-source counter.  common /isomul/imsrc

           ispfn = 0     ! common /isorfs/
           rspfn = 0.0d0 ! common /isorfs/
           rspfz = 0.0d0 ! common /isorfs/

           ibcsn = 0     ! common /ibchsor/
         end if

      endif

*-----------------------------------------------------------------------
      return

*-----------------------------------------------------------------------
! T.Sato 2020/12/06 for cosmic-ray
  923 continue

         m_err = 'pg1 should be between -360 and 360 for dir = iso'
         ErrCha = ''
         ErrID = 'L:20441/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  924 continue

         m_err = 'pg2 should be between -360 and 360 for dir = iso'
         ErrCha = ''
         ErrID = 'L:20451/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  925 continue

         m_err = 'ag1 > ag2 is not allowed for dir = iso'
         ErrCha = ''
         ErrID = 'L:20461/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  926 continue

         m_err = 'Specified date is wrong (too past or too recent)'
         ErrCha = ''
         ErrID = 'L:20471/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  927 continue

         m_err = 'integral fluence is zero, check previous multisource'
         ErrCha = ''
         ErrID = 'L:20481/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  928 continue

         m_err = 'icenv should be same for all multi-source'
         ErrCha = ''
         ErrID = 'L:20491/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  929 continue

         m_err = 'There is no cosmic-ray source for the projectile'
         ErrCha = ''
         ErrID = 'L:20501/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return


*-----------------------------------------------------------------------
C S.H. added all option for phi (2018.10.2)
  930 continue

         m_err = 'phi should be greater than -1000 in [source]'
         ErrCha = ''
         ErrID = 'L:20514/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  936 continue

         m_err = 'jpsf can not be used with dump in [source]'
         ErrCha = ''
         ErrID = 'L:20525/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  937 continue

         m_err = 'jpsf discription is wrong in [source]'
         ErrCha = ''
         ErrID = 'L:20535/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  938 continue

         m_err = 'suf = is missing for s-type=26.'
         ErrCha = ''
         ErrID = 'L:20546/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  939 continue

         m_err = '# of cut should be less than 9 for s-type=26.'
         ErrCha = ''
         ErrID = 'L:20556/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  940 continue

         m_err = 'exa should not be negative.'
         ErrCha = ''
         ErrID = 'L:20567/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  941 continue

         m_err = 'Parameters ni is already specified in [source]'
         ErrCha = ''
         ErrID = 'L:20579/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  942 continue

         m_err = 'After [ ni = ], [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:20591/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  943 continue

         m_err = 'Parameters ne is already specified in [source]'
         ErrCha = ''
         ErrID = 'L:20603/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  944 continue

         m_err = 'Data of RI source (nuclide and activity) is wrong'
         ErrCha = ''
         ErrID = 'L:20615/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  945 continue

         m_err = 'Parameters ni or ne must be set under e-type=28, 29.'
         ErrCha = ''
         ErrID = 'L:20627/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  946 continue

         m_err = 'Parameters dtime must be set under e-type=28, 29.'
         ErrCha = ''
         ErrID = 'L:20639/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  950    write(dkam,'(i9)') mdas
         m_err = 'Memory over in [source] : mdas ='//dkam
         ErrCha = ''
         ErrID = 'L:20650/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  961 continue

         m_err = 'p-type is not necessary in this e-type'
         ErrCha = ''
         ErrID = 'L:20660/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  962 continue

         m_err = 'number of multi source is exceed isrc'
         ErrCha = ''
         ErrID = 'L:20670/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  963 continue

         m_err = 'multi correlation source should be with multi source'
         ErrCha = ''
         ErrID = 'L:20680/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  964 continue

         m_err = 'projectile should be one in <source>'
         ErrCha = ''
         ErrID = 'L:20690/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  965 continue

         m_err = 'a-type should be set under dir = data.'
         ErrCha = ''
         ErrID = 'L:20700/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  966 continue

         m_err = 'dom=-10 should be set under s-type=1, 2, 4, 5.'
         ErrCha = ''
         ErrID = 'L:20710/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  967 continue

         m_err = 'nm should be set under e-type=3, 5, 6, 7, 15, 16.'
         ErrCha = ''
         ErrID = 'L:20720/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  968 continue

         m_err = 'inner radius r1 is greater than r0 in [source]'
         ErrCha = ''
         ErrID = 'L:20730/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  969 continue

         m_err = 'lmax should be dl2.'
         ErrCha = ''
         ErrID = 'L:20740/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  970 continue

         m_err = 'lmin should be dl1.'
         ErrCha = ''
         ErrID = 'L:20750/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  971 continue

         m_err = 'After nl = , [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:20760/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  972 continue

         m_err = 'ne(na) or nm(nn) should be defined before p(q)-type'
         ErrCha = ''
         ErrID = 'L:20770/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  973 continue

         m_err = 'slit size should be defined'
         ErrCha = ''
         ErrID = 'L:20780/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  974 continue

         m_err = 'slit position should be < dl0'
         ErrCha = ''
         ErrID = 'L:20790/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  975 continue

         m_err = 'dsxp,dsxn,dsyp,dsyn should be > 0'
         ErrCha = ''
         ErrID = 'L:20800/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  976 continue

         m_err = 'dl1 and dl2 should be dl2 > dl1 > 0'
         ErrCha = ''
         ErrID = 'L:20810/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  977 continue

         m_err = 'dpf should be 0 < dpf < 1'
         ErrCha = ''
         ErrID = 'L:20820/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  978 continue

         m_err = 'dxw or dyw is missing for square duct source'
         ErrCha = ''
         ErrID = 'L:20830/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  979 continue

         m_err = 'drd is missing for cylinder duct source'
         ErrCha = ''
         ErrID = 'L:20840/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  980 continue

         m_err = 'dl0, dl1, dl2 or dpf is missing for duct source'
         ErrCha = ''
         ErrID = 'L:20850/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  981 continue

         m_err = 'dump data is lack in [source]'
         ErrCha = ''
         ErrID = 'L:20860/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  982 continue

         m_err = 'dump discription is wrong in [source]'
         ErrCha = ''
         ErrID = 'L:20870/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  983 continue

         m_err = 'projectile name is wrong in [source]'
         ErrCha = ''
         ErrID = 'L:20880/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  984 continue

         m_err = 'After [ p-type = ], [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:20890/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  885 continue

         m_err = 'sphere radius r1 should be larger than r2 for dir=iso'
         ErrCha = ''
         ErrID = 'L:20900/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  985 continue

         m_err = 'inner radius r1 is greater than r2 in [source]'
         ErrCha = ''
         ErrID = 'L:20910/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  986 continue

         m_err = 'dimension of parabola should be even in [source]'
         ErrCha = ''
         ErrID = 'L:20920/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  987 continue

         m_err = 't-type is wrong in [source]'
         ErrCha = ''
         ErrID = 'L:20930/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  988 continue

         m_err = 'distance between time peaks is not '//
     &           'defined in [source]'
         ErrCha = ''
         ErrID = 'L:20941/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  989 continue

         m_err = 'number of peaks is wrong in [source]'
         ErrCha = ''
         ErrID = 'L:20951/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  990 continue

         m_err = 't-type is not 0, '//
     &           'but width is not defined in [source].'
         ErrCha = ''
         ErrID = 'L:20962/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  991 continue

         m_err = 'dir = all is not available in s-type = 12'
         ErrCha = ''
         ErrID = 'L:20972/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'dir = all is not available in s-type = 11'
         ErrCha = ''
         ErrID = 'L:20984/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'abs( dir ) .gt. 1.0 in [source]'
         ErrCha = ''
         ErrID = 'L:20996/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  872 continue

         m_err = 'nt or ll should be defined before o-type'
         ErrCha = ''
         ErrID = 'L:21007/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'After [ ne = ], [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:21019/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Input file in [source] does not exist. file name = '//
     &            sfile(j)(1:lsfile(j))
         ErrCha = ''
         ErrID = 'L:21032/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  960 continue

         m_err = 'dmpmul = 0.0 is not allowed'//
     &        ' when idmpmode = 1.'
         ErrCha = ''
         ErrID = 'L:21044/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  959 continue

         m_err = 'istdev < 0 is not allowed'//
     &        ' when idmpmode = 1 (default).'
         ErrCha = ''
         ErrID = 'L:21056/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  958 continue

         m_err = '[idmpmode = 1] (default) is not compatible'//
     &        ' with multi source.'
         ErrCha = ''
         ErrID = 'L:21068/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  957 continue

         m_err = 'istdev < 0 is not allowed'//
     &        ' when idmpmode = 1 or dmpmulti not= 0.0.'
         ErrCha = ''
         ErrID = 'L:21080/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  956 continue

         m_err = '[idmpmode = 1] is not compatible with multi source.'
         ErrCha = ''
         ErrID = 'L:21091/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  855 continue ! Fail to MPI dump file check -------------

         write(*,'(/" ***** Error Message from Dump File *****"/)')
         call MPI_dumpFileCheck(ndumpmax)
         m_err = '# of dump files is inconsistent.'
         ErrCha = ''
         ErrID = 'L:21105/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

         return

*-----------------------------------------------------------------------

  955 continue

         m_err = 'Dump file in [source] does not exist. file name = '//
     &            filnm(1:lfiln)
         ErrCha = ''
         ErrID = 'L:21119/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

         return

*-----------------------------------------------------------------------

  954 continue

         m_err = 'Dump file in [source] is wrong. file name = '//
     &            filnm(1:lfiln)
         ErrCha = ''
         ErrID = 'L:21133/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  953 continue

         m_err = 'Dump file in [source] is empty. file name = '//
     &            filnm(1:lfiln)
         ErrCha = ''
         ErrID = 'L:21146/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  952 continue

         m_err = 'dmpmulti should be positive or =0.'
         ErrCha = ''
         ErrID = 'L:21157/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
 951     continue

         m_err = 'tetreg is required for s-type = 24 or 25.'
         ErrCha = ''
         ErrID = 'L:21168/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Data of energy bin in [source] is wrong'
         ErrCha = ''
         ErrID = 'L:21180/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

      if(ierrMSG.eq.1) then ! T.Sato 2017/07/14
       m_err = 'mesh should be xyz in [source]'
       ErrCha = ''
       ErrID = 'L:21193/R:sours/F:read02.f'
      elseif(ierrMSG.eq.2) then
       m_err='s-type should be defined in [source] or previous <source>'
       ErrCha = ''
       ErrID = 'L:21197/R:sours/F:read02.f'
      elseif(ierrMSG.eq.3) then
       m_err = 'proj should be defined in [source] or previous <source>'
       ErrCha = ''
       ErrID = 'L:21201/R:sours/F:read02.f'
      elseif(ierrMSG.eq.4) then
       m_err = 'Unknown s-type in [source] or previous <source>'
       ErrCha = ''
       ErrID = 'L:21205/R:sours/F:read02.f'
      elseif(ierrMSG.eq.5) then
       m_err='Either e0 or e-type is necessary'
     & //' in [source] or previous <source>'
       ErrCha = ''
       ErrID = 'L:21210/R:sours/F:read02.f'
      else
       m_err = 'Something is wrong in [source]'
       ErrCha = ''
       ErrID = 'L:21214/R:sours/F:read02.f'
      endif

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

         return

*-----------------------------------------------------------------------

  998 continue  ! T.Sato 2024/09/23 error message has been already assigned

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Unknown parameter in [source]'
         ErrCha = ''
         ErrID = 'L:21238/R:sours/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
      subroutine ripmake(j,npcunt,iripmode,ipaflg,issflg,chout,
     &                   irptyp,ierr)
*                                                                      *
*       expansion to multi-source subsection when "proj=all"           *
*       create by T.Miura on 2021/11/30                                *
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'risrcparam.inc'
      include 'err.inc'

*     risrcparam.inc: include parameters
*     integer nuclinmax, nuclpumax, maxchain, chainmax
*     parameter ( nuclinmax =    100 ) nuclides including daughter
*     parameter ( nuclpumax =  1,000 ) Not used
*     parameter ( maxchain  =     23 ) from Number of chain (.NDX)
*     parameter ( chainmax  =  4,050 ) from Number of linear chain (.NDX)
*
*     integer rimax1, rimax2
*     parameter ( rimax1    =  4,000 ) from Number of auger electron (.ACK)
*     parameter ( rimax2    = 20,000 ) rimax1 * 5
!----------------------------------------------------------------------
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

      common /risource/ chalct(nuclinmax,isrc), act000(nuclinmax,isrc),
     &                  activy(nuclinmax,isrc), thf(nuclinmax,isrc),
     &                  decayt(isrc)
      real  chalct  ! ex) 55137.0 ( = Cs-137 )
      double precision  act000, activy, thf, decayt

      common /risrc00/ niorg(isrc), nicur(isrc), norm(isrc),
     &                 iaugers(isrc), icharacterx(isrc), iannih(isrc),
     &                 aclow(isrc)
      integer  niorg, nicur, norm, iaugers, icharacterx, iannih
      double precision  aclow
!----------------------------------------------------------------------
      character m_err*200
      common /error/ m_err, l_err, k_err
!----------------------------------------------------------------------
      character chin*200, chlw*200, chcm*200
      character chwrk*200

      logical   exex
      character filnm*100
      character NUCL*7

      real chza
      character chaa*7

      integer,intent(in) :: j       ! multi-source number
      integer,intent(inout) :: npcunt  ! total incident particles
      integer,intent(in) :: iripmode  ! process mode for risearch sub. 1:<> sub-section, 2:[] section
      integer,intent(in) :: ipaflg  ! flag of proj=all
      integer,intent(in) :: issflg  ! flag of <source> subsection input
      character,intent(in) :: chout*200

      integer,intent(out):: irptyp  ! (output) particle type
      integer,intent(out):: ierr

      common /isomul/ smlwt(isrc), totfact, imsrc
      real*8,save :: totfact_dum = 0.d0

!----------------------------------------------------------------------

      character*2   nu                  ! nuclide
      character*4   ms                  ! mass
      character*7   RI
      integer       irismax
      parameter(irismax=5)   ! Maximum number of RI source type
      integer       irflag(irismax,nuclinmax) ! photon/electorn/positoron/alpha/neutron
      integer       nirps (irismax)

      common /mpi00/ npe, me
      character chprojall*200
      character chme*5


!======================================================================

      ierr   = 0
      irptyp = 0

      if( npe .le. 1 ) then
         chprojall = 'risrc-extract.tmp'
      else
         iorder = aint(log10(real(npe))) + 1
         if ( iorder .lt. 3) iorder = 3
         write(chme,'(i5.5)') me
         chprojall = 'risrc-extract.tmp'//chme(6-iorder:5)
      endif
      open(noms,file=chprojall,
     &          form='formatted',status='unknown',access='append')

      if ( (j.gt.0) .and. (ipaflg.gt.1) ) then  ! check proj=all ?

!----------------------------------------------------------------------
!     search incident particles
!----------------------------------------------------------------------

        filnm = chfn(21)(1:ilfn(21)) // 'phits_risors.dat'   ! /dchain-sp/data/
        inquire( file = filnm(1:ilfn(21)+16), exist = exex )
        if( exex .eqv. .false. ) then
          write(ErrCha,'("incident particles file error")')
          ErrID = 'L:21357/R:ripmake/F:read02.f'
          call ErrWrite(ErrID,ErrCha)
          goto 995
        endif
        open(nirp,file=filnm,status='old')

        nirps(:) = 0
        do in = 1, niorg(j)      ! loop for number of nuclide
          rewind(nirp)
          read(nirp,*)           ! skip 1st-line comment

          chza = chalct(in,j)
          call getelm(12,chaa,chza,ierr)
*                     12: chza(Zaid+) -> chaa,7digit(H1     )
          NUCL = chaa

  100     continue
          read(nirp,'(2x,a2,1x,a4,5i3)',end=996)
     &                   nu,   ms,(irflag(ip,in),ip=1,irismax)

cABE 2023/05/25, if iannih=0, positron frag is turned off to avoid double count the annihilation photon
          if( iannih(j) .eq. 0 ) irflag(3,in) = 0
cABE 2023/05/25, if iannih=0, photon flag is turned on for C-11, N-13, O-15 and F-18 to consider the annihilation photon
          if( iannih(j) .eq. 0 .and.
     &        (( nu .eq. " C" .and. ms .eq. " 11 " ) .or.
     &        ( nu .eq. " N" .and. ms .eq. " 13 " ) .or.
     &        ( nu .eq. " O" .and. ms .eq. " 15 " ) .or.
     &        ( nu .eq. " F" .and. ms .eq. " 18 " )) ) then ! revised T.Sato 2023/09/09
             irflag(1,in) = 1
          endif

          RI = trim(adjustL(nu)) // trim(adjustL(ms)) // " "
          if ( NUCL .ne. RI ) go to 100   ! not fund, next data

          do ip = 1, irismax
            nirps(ip) = nirps(ip) + irflag(ip,in)
          end do
        end do      ! in = 1, niorg(j)

        close (nirp)   ! phits_risors.dat

!----------------------------------------------------------------------
!       write multi-source subsections
!----------------------------------------------------------------------


        do ip = 1, irismax     ! photon/electron/positoron/alpha/neutron
          if ( nirps(ip) .gt. 0 ) then
            rewind(niws)   ! work souce section file
            totfact_dum = totfact_dum + smlwt(imsrc)

            if( issflg.eq.0 ) write(noms,'(" <Source> =   1.000")')

  200       continue
            ! read temp. source subsection data
            read(niws,'(a200)', iostat=ios ) chin
            if( ios .eq. -1 ) cycle   ! EOF
            call chlngt(chin,200,i1,i2)
            chlw = chin
            call chcaps(chlw,i1,i2,i3,'#!$')
            chcm = chlw
            call chcomp(chcm,i1,i3,i4)

!......................................................................
            if ( index(chcm(i1:i4),'ni=') .gt. 0 ) then    ! nuclide and activity
              write(noms,'(7x,"ni = ",i3)') nirps(ip)      ! ni=xxx
              chwrk = ' '                                  ! buf clear

              read(niws,'(a200)', iostat=ios ) chin
              if( ios .eq. -1 ) go to 997   ! EOF
              call chlngt(chin,200,i1,i2)
              chlw = chin
              call chcaps(chlw,i1,i2,i3,'#!$')
              chcm = chlw
              call chcomp(chcm,i1,i3,i4)

              ic = i1

              do i = 1, niorg(j)
                icl = inumc(chlw, ic, i3, ' ') - 1

                if ( irflag(ip,i) .gt. 0 )
     &            chwrk(ic:icl)=chin(ic:icl)               ! set nuclide

                icl = icl + 1
                ic = jnumc(chlw, icl, i3)
                if( i .le. niorg(j) .and. ic .gt. i3 ) then
                  write(noms,'(a)') trim(chwrk)            ! write buf
                  chwrk = ' '                              ! buf clear

                  read(niws,'(a200)', iostat=ios ) chin
                  if( ios .eq. -1 ) go to 997   ! EOF
                  call chlngt(chin,200,i1,i2)
                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,'#!$')
                  chcm = chlw
                  call chcomp(chcm,i1,i3,i4)

                  ic = i1
                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                else
                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                end if

                if ( irflag(ip,i) .gt. 0 )
     &            chwrk(ic:ic2-1)=chin(ic:ic2-1)           ! set activity

                if( i .eq. niorg(j) ) then
                  write(noms,'(a)') trim(chwrk)            ! write buf
                  chwrk = ' '                              ! buf clear
                else
                  ic = ic2
                  if( i .lt. niorg(j) .and. ic .gt. i3 ) then
                    write(noms,'(a)') trim(chwrk)          ! write buf
                    chwrk = ' '                            ! buf clear

                    read(niws,'(a200)', iostat=ios ) chin
                    if( ios .eq. -1 ) go to 997   ! EOF
                    call chlngt(chin,200,i1,i2)
                    chlw = chin
                    call chcaps(chlw,i1,i2,i3,'#!$')
                    chcm = chlw
                    call chcomp(chcm,i1,i3,i4)

                    ic = i1
                  else
                    ic = jnumc(chlw, ic, i3)
                  end if
                end if

              end do   ! i = 1, niorg(j)

              if ( chwrk .ne. ' ' ) write(noms,'(a)') trim(chwrk) ! write buf

!......................................................................
            else  if ( chcm(i1:i4) .eq. 'proj=all' ) then  ! replace proj=xxx
              npcunt = npcunt + 1
              select case (ip)
                case (1)
                  write(noms,'(5x,"proj = photon")')
                  if ( irptyp .eq. 0 ) irptyp = 1
                case (2)
                  write(noms,'(5x,"proj = electron")')
                  if ( irptyp .eq. 0 ) irptyp = 2
                case (3)
                  write(noms,'(5x,"proj = positron")')
                  if ( irptyp .eq. 0 ) irptyp = 3
                case (4)
                  write(noms,'(5x,"proj = alpha")')
                  if ( irptyp .eq. 0 ) irptyp = 4
                case (5) ! frtati 2022/12/26 for SF neutron
                  write(noms,'(5x,"proj = neutron")')
                  if ( irptyp .eq. 0 ) irptyp = 5
              end select

!......................................................................
            else
              write(noms,'(a)') chin(1:i2)
            end if

            go to 200   ! next line form scrach file

          end if        ! nirps(ip) > 0
        end do          ! i = 1, irismax (photon/electron/positoron/alpha/neutron)

!----------------------------------------------------------------------
!     Do not expand the input data.
!     Proj /= "all" or [source] section common parameter.
!----------------------------------------------------------------------

      else               ! no assignment nuclide
        rewind(niws)     ! work souce section file

        do
          read (niws,'(a)',end=899) chin
          write(noms,'(a)') trim(chin)
        end do

  899   continue
        if ( ipaflg.ge.1 ) totfact_dum = totfact_dum + smlwt(imsrc)
      end if

!======================================================================
!     normal return
      select case (iripmode)
        case (1)      ! <source> sub-section

        case (2)      ! [] section
          totfact_dum = totfact_dum * dsign(1.d0, totfact)
          write(noms,'(/2x,"totfact = ",f0.1/)') totfact_dum
          write(noms,'(a)') trim(chout)

      end select

  900 continue
      close (noms)   ! extend multi-source sub-section (risrc-extract.tmp)

      rewind(niws)   ! work souce section file (scratch file)

      return

!----------------------------------------------------------------------
!     file end (phits_risors.dat)
  995 continue
      m_err = 'file does not exist. file name = '//
     &         filnm(1:ilfn(21)+16)
      ErrCha = ''
      ErrID = 'L:21564/R:ripmake/F:read02.f'
      l_err = 1
      k_err = 1
      ierr  = 1

      close(noms)
      return

  996 continue
      m_err = NUCL//' is not in the list of RI source'
      ErrCha = ''
      ErrID = 'L:21575/R:ripmake/F:read02.f'
      l_err = 1
      k_err = 1
      ierr  = 1

      close(noms)
      close(nirp)
      return

  997 continue
      m_err = 'nuclide and activity error.'
      ErrCha = ''
      ErrID = 'L:21587/R:ripmake/F:read02.f'
      l_err = 1
      k_err = 1
      ierr  = 1

      close(noms)
      close(nirp)
      return
      end subroutine ripmake


************************************************************************
*                                                                      *
      subroutine etalfile(j,isorf,ndata,cmmt,ierr)
*                                                                      *
*       read energy distribution from tally output file                *
*       modified by S.Hashimoto on 2016.12.30                          *
*                                                                      *
*        Output files of [t-track], [t-cross], [t-point], [t-product], *
*      [t-time], and [t-star] tallies with axis=eng can be used.       *
*                                                                      *
************************************************************************
      use moddas
      use moddas_source

      implicit none
      include 'param.inc'

*-----------------------------------------------------------------------

      integer j, isorf, ierr
      character*(*) cmmt
      integer iii, ios, i1,i2,i3,i4, il, ic, icl
      character chin*200, chlw*200, chcm*200
      double precision cvvv

*-----------------------------------------------------------------------

      integer ipart, ipage
      integer itally, iunit, ideriv, ndata, ifrag, idata, icolumn
      double precision el, eu, ycolumn(12), ysum

*-----------------------------------------------------------------------

      common /isorse/ ngrp(isrc), ngei(isrc), ngea(isrc), ngfe(isrc),
     &                ngft(isrc), ngll(isrc), ngpi(isrc), ngpw(isrc)
      integer ngrp, ngei, ngea, ngfe, ngft, ngll, ngpi, ngpw

*-----------------------------------------------------------------------

      integer navtal, iavtal
      parameter ( navtal = 6 ) ! number of available tallies
      character cavtal(navtal)*12
      integer lcavtal(navtal)
      data ( cavtal(iavtal), iavtal = 1, navtal ) /
     &    '[t-track]   ','[t-cross]   ','[t-point]   ','[t-product] ',
     &    '[t-time]    ','[t-star]    '/
      data ( lcavtal(iavtal), iavtal = 1, navtal ) /
     &     9            ,9             ,9             ,11            ,
     &     8            ,8             /

*-----------------------------------------------------------------------
      integer    inumc, jnumc
      external   inumc, jnumc
*-----------------------------------------------------------------------
* initialization
      do iii=1,200
         chin(iii:iii)=' '
      end do

      ipart = 1
      ipage = 1
      ngll(j) = 1

*-----------------------------------------------------------------------
* to check a kind of tally
      read(isorf,'(a200)',iostat=ios) chin
      if( ios .eq. -1 ) goto 900
      call chlngt(chin,200,i1,i2)
      chlw = chin
      call chcaps(chlw,i1,i2,i3,cmmt)
      chcm = chlw
      call chcomp(chcm,i1,i3,i4)
      itally = 0
      do iavtal = 1, navtal
       if( chcm(i1:i4) .eq. cavtal(iavtal)(1:lcavtal(iavtal)) ) then
          itally = iavtal
       end if
      end do
      if( itally .eq. 0 ) goto 900


* to check unit and axis of tally output file
      ifrag = 0
      do while ( ifrag .ne. 1 )
       read(isorf,'(a200)',iostat=ios) chin
       if( ios .eq. -1 ) goto 900

       call chlngt(chin,200,i1,i2)
       chlw = chin
       call chcaps(chlw,i1,i2,i3,cmmt)
       chcm = chlw
       call chcomp(chcm,i1,i3,i4)

       if(i1 .eq. 0) cycle


       il = i1 + 4 - 1
       if( chcm(i1:il) .eq. 'unit' ) then ! start of checking unit
        ic = inumc(chlw,il+1,i3,'=') + 1
        ic = jnumc(chlw,ic,i3)
        if( ic .gt. i3 ) goto 900
        call onum(chlw,ic,i3,cvvv,ierr)
        iunit = nint( cvvv )

        ideriv = 0 ! ideriv=0: non-derivative, =1: derivative
* [t-track]
        if ( itally .eq. 1 ) then
         if ( iunit.eq.2 .or. iunit.eq.12 ) then
            ideriv = 1
            ngll(j) = 1
         else if ( iunit.eq.3 .or. iunit.eq.13 ) then
            ideriv = 1
            ngll(j) = -1
         end if

* [t-cross]
        else if ( itally .eq. 2 ) then
         if ( iunit.eq.2 .or. iunit.eq.5 .or.
     &          iunit.eq.12 .or. iunit.eq.15 ) then
            ideriv = 1
            ngll(j) = 1
         else if ( iunit.eq.3 .or. iunit.eq.6 .or.
     &          iunit.eq.13 .or. iunit.eq.16 ) then
            ideriv = 1
            ngll(j) = -1
         end if

* [t-point]
        else if ( itally .eq. 3 ) then
         if ( iunit.eq.2 .or. iunit.eq.12 ) then
            ideriv = 1
            ngll(j) = 1
         else if ( iunit.eq.3 .or. iunit.eq.13 ) then
            ideriv = 1
            ngll(j) = -1
         end if

* [t-product]
        else if ( itally .eq. 4 ) then
         if ( iunit.eq.3 .or. iunit.eq.4 .or.
     &          iunit.eq.13 .or. iunit.eq.14 .or.
     &          iunit.eq.23 .or. iunit.eq.24.or.
     &          iunit.eq.33 .or. iunit.eq.34 ) then
            ideriv = 1
            ngll(j) = 1
         else if ( iunit.eq.5 .or. iunit.eq.6 .or.
     &          iunit.eq.15 .or. iunit.eq.16 .or.
     &          iunit.eq.25 .or. iunit.eq.26 .or.
     &          iunit.eq.35 .or. iunit.eq.36 ) then
            ideriv = 1
            ngll(j) = -1
         end if

* [t-time]
        else if ( itally .eq. 5 ) then
         if ( iunit.eq.4 ) then
            ideriv = 1
            ngll(j) = 1
         end if

* [t-star]
        else if ( itally .eq. 6 ) then
         if ( iunit.eq.2 .or. iunit.eq.12 ) then
            ideriv = 1
            ngll(j) = 1
         end if
        end if

       end if ! end of checking unit


       il = i1 + 4 - 1
       if( chcm(i1:il) .eq. 'axis' ) then
        ic = inumc(chlw,il+1,i3,'=') + 1
        ic = jnumc(chlw,ic,i3)
        if( ic .gt. i3 ) goto 900
        if ( chlw(ic:i3) .ne. 'eng' ) goto 900
       end if


       if( chin(1:4) .eq. 'h: n' ) ifrag = 1

      end do

*-----------------------------------------------------------------------
* to count the number of energy bin

      read(isorf,*)
      ndata = 0
      ifrag = 0
      do while ( ifrag .ne. 1 )

       read(isorf,'(a200)',iostat=ios) chin
       if( ios .eq. -1 ) goto 900
       call chlngt(chin,200,i1,i2)
       if( i1 .eq. 0 .and. i2 .eq. 0 ) ifrag = 1
       ndata = ndata + 1

      end do
      ndata = ndata - 1

      do idata = 1, ndata+1
       backspace(isorf)
      end do

*-----------------------------------------------------------------------
* to set energy distribution

      ngrp(j)  = ndata
      call moddas_reallocate_dbl(isrc, j, ngrp(j), ngei, egmin)
      call moddas_reallocate_dbl(isrc, j, ngrp(j), ngea, egmax)
      call moddas_reallocate_dbl(isrc, j, ngrp(j), ngfe, fegrp)
      call moddas_reallocate_dbl(isrc, j, ngrp(j), ngft, rfe)
      call moddas_reallocate_dbl(isrc, j, ngrp(j), ngpi, prw)
      call moddas_reallocate_dbl(isrc, j, ngrp(j), ngpw, pwt)

      ysum = 0d0
      do idata = 1, ndata

       read(isorf,'(1p2e13.4,6(1pe13.4,0pf8.4))')
     &        el,eu,(ycolumn(icolumn), icolumn=1,12)

       egmin(ngei(j)+idata) = el
       egmax(ngea(j)+idata) = eu
       fegrp(ngfe(j)+idata) = ycolumn(2*ipart-1)

       if ( ideriv .eq. 0 ) then
          prw(ngpi(j)+idata) = ycolumn(2*ipart-1)

       else if ( ideriv.eq.1 .and. ngll(j).eq.1 ) then
          prw(ngpi(j)+idata) = ycolumn(2*ipart-1) * (eu-el)

       else if ( ideriv.eq.1 .and. ngll(j).eq.-1 ) then
          prw(ngpi(j)+idata) = ycolumn(2*ipart-1) * dlog(eu/el)

       end if

       ysum = ysum + ycolumn(2*ipart-1)

      end do
      if ( ysum .le. 0d0 ) goto 900

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

 900  continue
      ierr  = 1

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine region(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [region] section of input files                           *
*       modified by K.Niita on 12/06/2000                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_character

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar

      common /regdu/  iuni(kvlmax)
      common /regdm/  idmg(kvlmax)
      common /regda/  ichl(kvlmax), chsm(kvlmax), ichmx, iod
      character       chsm*10
      common /regdb/  nrsq, irsq(10)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /regdd/  ivolm, iimpo
      common /regdn/  irden

      common /impreg/ dimp(kvlmax)
      common /volreg/ dvol(kvlmax)
      common /celda/  deng(kvlmax)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chtm2*2

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      data irsq / 1, 2, 4, 100, 6 * 0 /

      character dkam*6

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

            if( iregn .ne. 0 ) goto 989

            if( igcel .ne. 0 .or. igsuf .ne. 0 ) goto 988

*-----------------------------------------------------------------------

            mci  = ( mmmax - 1 ) * 8 + 1
            mcmx = ( mdas  - 1 ) * 8 + 1 - mci

            iod = 26

            open(iod,status='scratch',form='unformatted')

            ichmx = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

               call moddas_allocate_cha(MAX_NUM_CHRG, chrg)
               mci = 0
  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 .and. iregn .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               irsq( mrsq + 1 ) = 0

            else if( chlw(ic:ic+2) .eq. 'num' ) then

               irsq( mrsq + 1 ) = 1

            else if( chlw(ic:ic+2) .eq. 'mat' ) then

               irsq( mrsq + 1 ) = 2

            else if( chlw(ic:ic+2) .eq. 'imp' ) then

               irsq( mrsq + 1 ) = 3

            else if( chlw(ic:ic+2) .eq. 'sym' ) then

               irsq( mrsq + 1 ) = 4

            else if( chlw(ic:ic+2) .eq. 'uni' ) then

               irsq( mrsq + 1 ) = 5

            else if( chlw(ic:ic+2) .eq. 'vol' ) then

               irsq( mrsq + 1 ) = 6

            else if( chlw(ic:ic+2) .eq. 'den' ) then

               irsq( mrsq + 1 ) = 7

            else if( chlw(ic:ic+2) .eq. 'def' ) then

               irsq( mrsq + 1 ) = 100

            else

               goto 997

            end if

               mrsq = mrsq + 1

               ic = ic + 3
               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

               if( mrsq .eq. 1 ) goto 997
               if( irsq(mrsq) .ne. 100 ) goto 996

                  inumb = 0
                  imatt = 0
                  iimpo = 0
                  isymb = 0
                  iuniv = 0
                  ivolm = 0
                  irden = 0

               do k = 1, mrsq - 1
                  if( irsq(k) .eq. 1 ) inumb = inumb + 1
                  if( irsq(k) .eq. 2 ) imatt = imatt + 1
                  if( irsq(k) .eq. 3 ) iimpo = iimpo + 1
                  if( irsq(k) .eq. 4 ) isymb = isymb + 1
                  if( irsq(k) .eq. 5 ) iuniv = iuniv + 1
                  if( irsq(k) .eq. 6 ) ivolm = ivolm + 1
                  if( irsq(k) .eq. 7 ) irden = irden + 1
               end do

                  if( imatt .eq. 0 ) goto 995

                  if( inumb .gt. 1 .or. imatt .gt. 1 .or.
     &                iimpo .gt. 1 .or. isymb .gt. 1 .or.
     &                iuniv .gt. 1 .or. ivolm .gt. 1 .or.
     &                irden .gt. 1 ) goto 997

                  nrsq = mrsq

                  goto 140

            else

                  inumb = 1
                  imatt = 1
                  iimpo = 0
                  isymb = 1
                  iuniv = 0
                  ivolm = 0
                  irden = 0

                  nrsq = 4

            end if

         end if

*-----------------------------------------------------------------------
*        sequential region definition lines
*-----------------------------------------------------------------------

         if( chlw(i1:i1) .eq. '+' .or.
     &       chlw(i1:i1) .eq. '-' .or.
     &       chlw(i1:i1+2) .eq. 'or ' ) then

               if( ichl(iregn) .le. 2 ) goto 991

               if( ichl(iregn) + i3 - i1 + 2 .gt. mcmx ) goto 990

                     chrg(mci+ichl(iregn)+1:mci+ichl(iregn)+1) = ' '
                     ichl(iregn) = ichl(iregn) + 1

                  do i = i1, i3

                     chrg(mci+ichl(iregn)+i-i1+1:
     &                    mci+ichl(iregn)+i-i1+1) =
     &               chin(i:i)

                  end do

                  ichl(iregn) = ichl(iregn) + i3 - i1 + 1

            goto 140

         end if

*-----------------------------------------------------------------------
*        read region informations
*-----------------------------------------------------------------------

               iregn = iregn + 1

               if( iregn .gt. 1 ) then

                  ichmx = max(ichmx,ichl(iregn-1))

                  write(iod) (chrg(mci+k:mci+k),k=1,ichl(iregn-1))

               end if

               if( iregn .gt. kvlmax ) goto 998

               ic2  = i1

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( irsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( irsq(k) .eq. 1 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               idrg(iregn) = nint( cvvv )


            else if( irsq(k) .eq. 2 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               idmg(iregn) = nint( cvvv )


            else if( irsq(k) .eq. 3 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               dimp(iregn) = cvvv


            else if( irsq(k) .eq. 4 ) then

               if( chlw(ic:ic) .lt. 'a' .or.
     &             chlw(ic:ic) .gt. 'z' ) goto 999

               ic2 = inumc(chlw,ic,i3,' ')
               chsm(iregn) = chin(ic:ic2-1)


            else if( irsq(k) .eq. 5 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               iuni(iregn) = nint( cvvv )


            else if( irsq(k) .eq. 6 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               dvol(iregn) = cvvv


            else if( irsq(k) .eq. 7 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               deng(iregn) = cvvv


            else if( irsq(k) .eq. 100 ) then

                  ichl(iregn) = i3 - ic + 1

                  do i = 1, ichl(iregn)

                     chrg(mci+i:mci+i) = chin(ic+i-1:ic+i-1)

                  end do

            end if

         end do

            if( inumb .eq. 0 ) idrg(iregn) = iregn
            if( iuniv .eq. 0 ) iuni(iregn) = 0

            if( isymb .eq. 0 ) then

                  irch = iregn / 100
                  irnm = iregn - irch * 100
                  write(chtm2,'(i2.2)') irnm
                  chsm(iregn) = char(ichar('a')+irch)//chtm2

            end if

*-----------------------------------------------------------------------
*        check region number and matter
*-----------------------------------------------------------------------

            if( idrg(iregn) .le. 0 .or.
     &          idrg(iregn) .ge. 9999 ) goto 994

            if( idgr(idrg(iregn)) .ne. 0 ) goto 993

               idgr(idrg(iregn)) = iregn

               if( idmg(iregn) .ge. 9999 ) goto 992
               if( idmg(iregn) .lt. -1 )   goto 992

*-----------------------------------------------------------------------

         goto 140

*-----------------------------------------------------------------------

 1000 continue

               ichmx = max(ichmx,ichl(iregn))

               write(iod) (chrg(mci+k:mci+k),k=1,ichl(iregn))

         if( ichmx > MAX_NUM_CHRG ) then
            write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &           'sub.region@read02.f'
     &              //' ?dimension over chrg?'
     &              //' ichmx > MAX_NUM_CHRG'
     &           ,' (ichmx=',ichmx,')'
     &           ,' (MAX_NUM_CHRG@moddas.f=',MAX_NUM_CHRG,')'
            ErrID = 'L:22259/R:region/F:read02.f'
            call ErrWrite(ErrID,ErrCha)
         endif
         call moddas_deallocate_cha(chrg)
      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  988 continue

         m_err = 'GG section is already appeared.'
         ErrCha = ''
         ErrID = 'L:22273/R:region/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  989 continue

         m_err = '[region] section is duplicated'
         ErrCha = ''
         ErrID = 'L:22285/R:region/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  990 continue

         m_err = 'Sequential line of definitin is too long '
         ErrCha = ''
         ErrID = 'L:22297/R:region/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'Sequential line of definitin is wrong'
         ErrCha = ''
         ErrID = 'L:22309/R:region/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         write(dkam,'(i6)') idmg(iregn)
         m_err = 'ID number of material should be'//
     &           ' -1, 0, 1 - 9998 = '// dkam
         ErrCha = ''
         ErrID = 'L:22323/R:region/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         write(dkam,'(i6)') idrg(iregn)
         m_err = 'ID number of region is duplicated. = '// dkam
         ErrCha = ''
         ErrID = 'L:22336/R:region/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         write(dkam,'(i6)') idrg(iregn)
         m_err = 'ID number of region should be'//
     &           ' 1 - 9998 = '// dkam
         ErrCha = ''
         ErrID = 'L:22350/R:region/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = '<mat> is missing in definition of [region].'
         ErrCha = ''
         ErrID = 'L:22362/R:region/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = '<def> should be the last in definition of [region].'
         ErrCha = ''
         ErrID = 'L:22374/R:region/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Definition of data sequences in [region] is wrong.'
         ErrCha = ''
         ErrID = 'L:22386/R:region/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Number of region exceeds kvlmax in param.inc.'
         ErrCha = ''
         ErrID = 'L:22398/R:region/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [region] is wrong.'
         ErrCha = ''
         ErrID = 'L:22410/R:region/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine bodyin(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [body] section of input files                             *
*       modified by K.Niita on 26/05/2000                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

      parameter ( ibmt = 23 )

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

      common /inggs/ iog, igcel, ioa, igsuf, iob, igtrs
      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------

      dimension bval(50)

      dimension ibva(50), iblg(50)
      character chbd(50)*3

      dimension iplg(4), ipva(4)
      character chpa(4)*5

      character chtit*60

      dimension irby(10)
      data irby / 2, 1, 100, 7*0 /

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension ibck(kvmmax)
      data ibck / kvmmax*0 /

      character dkam*6

*-----------------------------------------------------------------------

      data ( chbd(i), i = 1, ibmt ) /
     &         'arb','sph','rcc','rec','trc',
     &         'ell','box','wed','rpp','gel',
     &         'tor','qua','bpp','wpp','p  ',
     &         'px ','py ','pz ','ps ','c  ',
     &         'cx ','cy ','cz '/

      data ( iblg(i), i = 1, ibmt ) /
     &           3,    3,    3,    3,    3,
     &           3,    3,    3,    3,    3,
     &           3,    3,    3,    3,    2,
     &           2,    2,    2,    2,    2,
     &           2,    2,    2/

      data ( ibva(i), i = 1, ibmt ) /
     &          30,    4,    7,   12,    8,
     &           7,   12,   12,    6,   12,
     &           9,   10,    9,    9,    4,
     &           1,    1,    1,    7,    7,
     &           4,    4,    4/

      data chpa / 'ivopt', 'idbg ', 'ibod ', 'naz  '/
      data iplg /  5,       4,       4,       3     /
      data ipva /  0,       0,       1,       0     /

      data chtit /'Body Information'/

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

            ierr = 0
            nbdy = 0

         if( ibody .gt. 0 ) goto 995

         if( igcel .ne. 0 .or. igsuf .ne. 0 ) goto 991

*-----------------------------------------------------------------------
*        open temporary binary file : unit 18
*-----------------------------------------------------------------------

            itby = 18

            open(itby,form='unformatted',status='scratch')

*-----------------------------------------------------------------------
*        read title
*-----------------------------------------------------------------------

                  is = inumc(chin,i1,i3,']')
                  if( chin(is+1:is+1) .eq. ' ' .or.
     &                chlw(is+1:is+1) .eq. tub ) is = is + 1

            if( is .le. i3 ) then

               do i = 1, 60
                     j = is + i
                  if( j .le. i3 ) then
                     chtit(i:i) = chin(j:j)
                  else
                     chtit(i:i) = ' '
                  end if
               end do

            end if

            write(itby) chtit

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of body section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        read parameters
*-----------------------------------------------------------------------

            if( ibody .eq. 0 .and. nbdy .eq. 0 ) then

               icl = i1

  100          continue

               do i = 1, 4

                  il = icl + iplg(i) - 1

                  if( chlw(icl:il) .eq. chpa(i)(1:iplg(i)) ) then

                     ic = inumc(chlw,il+1,i3,'=') + 1

                     if( ic .gt. i3 ) goto 997

                     icl = inumc(chlw,ic,i3,';') - 1

                     ipm = i

                     call onum(chlw,ic,icl,cvvv,ierr)

                     if( ierr .ne. 0 ) goto 998

                     ipva(ipm) = nint( cvvv )

                     icl = jnumc(chlw,icl+2,i3)

                     if( icl .le. i3 ) goto 100

                     goto 140

                  end if

               end do

                  write(itby) ipva

            end if

*-----------------------------------------------------------------------
*        definition of the data sequences
*-----------------------------------------------------------------------

         if( ibody .eq. 0 .and. nbdy .eq. 0 ) then

               mbdy = 0

               ic = i1

  200       if( ic .gt. i3 ) goto 400

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               irby( mbdy + 1 ) = 0

            else if( chlw(ic:ic+2) .eq. 'sym' ) then

               irby( mbdy + 1 ) = 1

            else if( chlw(ic:ic+2) .eq. 'num' ) then

               irby( mbdy + 1 ) = 2

            else if( chlw(ic:ic+2) .eq. 'def' ) then

               irby( mbdy + 1 ) = 100

            else

               goto 400

            end if

               mbdy = mbdy + 1

               ic = ic + 3
               ic = jnumc(chlw,ic,i3)
               goto 200

  400       continue

            if( mbdy .gt. 0 ) then

               if( mbdy .eq. 1 ) goto 994
               if( irby(mbdy) .ne. 100 ) goto 994

                  inumb = 0
                  isymb = 0

               do k = 1, mbdy - 1
                  if( irby(k) .eq. 1 ) isymb = isymb + 1
                  if( irby(k) .eq. 2 ) inumb = inumb + 1
               end do

                  if( inumb .gt. 1 ) goto 994
                  if( isymb .eq. 0 ) goto 994
                  if( isymb .gt. 1 ) goto 994

                  nbdy  = mbdy

                  goto 140

            else

                  inumb = 1

                  nbdy  = 3

            end if

         end if

*-----------------------------------------------------------------------
*        read body informations
*-----------------------------------------------------------------------

               ibody = ibody + 1

               ic2 = i1

         do k = 1, nbdy

               ic = jnumc(chlw,ic2,i3)

            if( irby(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( irby(k) .eq. 1 ) then

               do i = 1, ibmt

                  il = ic + iblg(i) - 1

                  if( chlw(ic:il) .eq. chbd(i)(1:iblg(i)) ) goto 300

               end do

                  goto 999

  300          ipm = i
               ic2 = il + 1

            else if( irby(k) .eq. 2 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ibnum = nint( cvvv )


            else if( irby(k) .eq. 100 ) then

               do i = 1, ibva(ipm)

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 997

                  bval(i) = cvvv

                  ic = ic2

                  if( i .lt. ibva(ipm) .and. ic .gt. i3 ) then

  143                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) goto 997

                     if( iskip .ne. 0 ) goto 143

                     ic = 1

                  end if

               end do

            end if

         end do

               if( ipva(3) .eq. 0 ) ibnum = ibody

               if( ibnum .le. 0 .or. ibnum .ge. 9999 ) goto 992

               if( ibck(ibnum) .ne. 0 ) goto 993
               ibck(ibnum) = ibody

               write(itby) chbd(ipm), ibnum,
     &                    ibva(ipm), ( bval(i), i = 1, ibva(ipm) )


               goto 140

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  991 continue

         m_err = 'GG section already appeared.'
         ErrCha = ''
         ErrID = 'L:22794/R:bodyin/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         write(dkam,'(i6)') ibnum
         m_err = 'ID number of body should be'//
     &           ' 1 - 9998 = '// dkam
         ErrCha = ''
         ErrID = 'L:22808/R:bodyin/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         write(dkam,'(i6)') ibnum
         m_err = 'ID number of body is duplicated. = '// dkam
         ErrCha = ''
         ErrID = 'L:22821/R:bodyin/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'Definition of [body] is wrong .'
         ErrCha = ''
         ErrID = 'L:22833/R:bodyin/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = '[body] section already appeared.'
         ErrCha = ''
         ErrID = 'L:22845/R:bodyin/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of Body Parameters is wrong.'
         ErrCha = ''
         ErrID = 'L:22857/R:bodyin/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Value of parameter is wrong in [body].'
         ErrCha = ''
         ErrID = 'L:22869/R:bodyin/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Unknown Body Name.'
         ErrCha = ''
         ErrID = 'L:22881/R:bodyin/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine mater(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [material] section of input files                         *
*       modified by K.Niita on 30/05/2000                              *
*                                                                      *
************************************************************************
      use dedx_file
      use ELEDATAMOD, only : allocate_frac_ichem, chemi_form, ichem,
     & frac  ! fraction

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /kmat1o/ iom1, iom2, iom3
      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      dimension lschn(20), ischn(20)
      character schan(20)*10

      data icsu / 10 /

      data ( schan(i), i = 1, 10 ) /
     &    'gas     ','estep   ','nlib    ','plib    ','elib    ',
     &    'cond    ','pnlib   ','hlib    ','dedxfile','chem    '/

      data ( lschn(i), i = 1, 10 ) /
     &     3,         5,         4,         4,         4,
     &     4,         5,         4,         8,         4/

*-----------------------------------------------------------------------

      logical deqn1
      logical dcom2

      dimension iseq(10)
      character chin*200, chlw*200, chcm*200
      character ch60*60,ch200*200
      character dsin(0:9)*200
      dimension idsi(0:9)
      dimension ill(0:9), ilf(0:9)
      character ckam*3
      character dkam*6
      character hs*10, ht*10
      dimension ix(3)

      data iseq / 1, 2, 8*0 /

*-----------------------------------------------------------------------

         if( mxmat .gt. 0 ) goto 994

            ks   = 0 ! S.H. added (2016.12.16)
            ierr  = 0
            mxmat = 0
            kseq  = 0
            icelp = 0
            idmat = 0
            imts  = -1

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue
! T.Sato 2016/12/28, c is comment signal or not
           if(mstz(124).eq.1) then ! c is comment, former default
            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
           else ! c is not comment, current default
            call readm(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
           endif
               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of material section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        parameter or s(a,b) sequential lines
*-----------------------------------------------------------------------

         ic = i1

         if( icelp .ne. 0 ) goto 300

*-----------------------------------------------------------------------
*        definition of data sequence
*        nucleus = 1, density = 2
*-----------------------------------------------------------------------

         if( chlw(i1:i1+2) .eq. 'nuc' .or.
     &       chlw(i1:i1+2) .eq. 'den' ) then

               if( kseq  .ne. 0 ) goto 998
               if( mxmat .gt. 0 ) goto 998

  200       continue

            if( ic .lt. i3 ) then

               if( chlw(ic:ic+2) .eq. 'nuc' ) then

                  kseq = kseq + 1
                  iseq(kseq) = 1
                  ic = jnumc(chlw,ic+3,i3)
                  goto 200

               else if( chlw(ic:ic+2) .eq. 'den' ) then

                  kseq = kseq + 1
                  iseq(kseq) = 2
                  ic = jnumc(chlw,ic+3,i3)
                  goto 200

               end if

            end if

               if( kseq .gt. 2 ) goto 998

               goto 140

         end if

*-----------------------------------------------------------------------
*        read material information
*-----------------------------------------------------------------------

         if( chlw(i1:i1+2) .eq. 'mat' .or.
     &     ( chlw(i1:i1) .eq. 'm' .and. deqn1(chlw(i1+1:i1+1)) ) ) then

            if( chlw(i1:i1+2) .eq. 'mat' .and.
     &          chlw(i1+3:i1+3) .ne. '(' .and.
     &          chlw(i1+3:i1+3) .ne. '[' .and.
     &          chlw(i1+3:i1+3) .ne. '{' ) goto 997

            if( chlw(i1:i1+2) .eq. 'mat' ) then

                  ic = i1 + 3

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 997

            else

                  ic = i1 + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 996

            end if

                  idmat = nint( cvvv )

                  ic = ic2

*-----------------------------------------------------------------------
*        new material
*-----------------------------------------------------------------------

            mxmat = mxmat + 1

            if( mxmat .eq. 1 ) then

               iom1 = 23
               iom2 = 24

               open(iom1,form='unformatted',status='scratch')
               open(iom2,form='unformatted',status='scratch')

            end if

            if( mxmat .gt. 1 ) then

                     if( imts .eq. -1 ) imts = 0

                     write(iom1) nel
                     write(iom1) denh, libh
                     write(iom1) igas, istp, inlb, iplb, ielb, icnd
                     write(iom1) iulb, ihlb
                     write(iom1) imts
                     write(iom1) (idedx(ii),ii=1,20)

               if( nel .gt. 0 .or. imts .gt. 0 ) rewind iom2

               if( nel .gt. 0 ) then

                  do i = 1, nel

                      read(iom2) icha, masi, denst, libi
                     write(iom1) icha, masi, denst, libi

                  end do

               end if

               if( imts .gt. 0 ) then

                  do i = 1, imts

                      read(iom2) ix(1), ix(2), ix(3)
                     write(iom1) ix(1), ix(2), ix(3)

                  end do

               end if

               if( nel .gt. 0 .or. imts .gt. 0 ) rewind iom2

            end if


*-----------------------------------------------------------------------

                nel = 0
               denh = 0.0d0

               libh = ichar(' ') + 256 * ( ichar(' ')
     &                           + 256 * ichar(' ') )

               ichn = 0

               isdn = 0
               ks   = 0

               igas = 0
               istp = 0
               icnd = 0
               imts = -1

               inlb = 0
               iplb = 0
               ielb = 0
               iulb = 0
               ihlb = 0
               idedx = 0

            do i = 1, icsu
               ischn(i) = 0
            end do

*-----------------------------------------------------------------------

               if( mxmat .ge. kvlmax ) then

                  write(ckam,'(i3)') kvlmax
                  m_err = 'Number of material exceeds kvlmax = '//
     &                    ckam
                  ErrCha = ''
                  ErrID = 'L:23172/R:mater/F:read02.f'
                  l_err = ill(jsn)
                  k_err = jsn
                  ierr  = 1
                  return

               end if

               if( idmat .ge. kvmmax .or. idmat .le. 0 ) then

                  write(dkam,'(i6)') idmat
                  m_err ='ID number of material should '//
     &                   'be 1 - kvmmax-1, = '//dkam
                  ErrCha = ''
                  ErrID = 'L:23186/R:mater/F:read02.f'
                  l_err = ill(jsn)
                  k_err = jsn
                  ierr  = 1
                  return

               end if

               if( idnm(idmat) .gt. 0 ) then

                  write(dkam,'(i6)') idmat
                  m_err = 'ID number of material is duplicated = '//
     &                    dkam
                  ErrCha = ''
                  ErrID = 'L:23200/R:mater/F:read02.f'
                  l_err = ill(jsn)
                  k_err = jsn
                  ierr  = 1
                  return

               end if

*-----------------------------------------------------------------------

               idnm(idmat) = mxmat
               idmn(mxmat) = idmat

         end if

*-----------------------------------------------------------------------
*        read s(a,b) information
*-----------------------------------------------------------------------

         if( chlw(i1:i1+1) .eq. 'mt' .and.
     &       deqn1(chlw(i1+2:i1+2)) ) then

                  ic = i1 + 2

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 996

                  ismat = nint( cvvv )

                  if( ismat .ne. idmat ) goto 988

                  ic = ic2

                  imts = imts + 1

         end if

*-----------------------------------------------------------------------
*        read data
*-----------------------------------------------------------------------

            ic = ic - 1

  250       continue

            ic = ic + 1

            if( ic .gt. i3 ) goto 140

            ic = jnumc(chlw,ic,i3)

*-----------------------------------------------------------------------
*        keyward parameters ( gas, estep, nlib, plib, elib, cond )
*-----------------------------------------------------------------------

            do i = 1, icsu

               il = ic + lschn(i) - 1

               if( chlw(ic:il) .eq. schan(i)(1:lschn(i)) ) goto 270

            end do

               goto 280

  270       continue

               icelp = i
               ieqc  = 0

               ischn( icelp ) = ischn( icelp ) + 1
               if( ischn(icelp) .gt. 1 ) goto 992

*-----------------------------------------------------------------------

  300       continue

            if( ieqc .eq. 0 ) then

               ic = inumc(chlw,ic,i3,'=')
               if( ic .gt. i3 ) goto 991

               ieqc = 1
               ic = jnumc(chlw,ic+1,i3)
               if( ic .gt. i3 )  goto 140

            end if

*-----------------------------------------------------------------------

            if( icelp .eq. 1 ) then

               call snum(chlw,ic,i3,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 998

               igas = nint( cvvv )

            else if( icelp .eq. 2 ) then

               call snum(chlw,ic,i3,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 998

               istp = nint( cvvv )

            else if( ( icelp .ge. 3 .and. icelp .le. 5 ) .or.
     &                 icelp .eq. 7 .or.  icelp .eq. 8 ) then

                  icm = 0
                  icn = 0
                  ica = 0
                  icb = 0
                  isa = 0
                  isn = 0

               i = ic - 1

  360          i = i + 1

               if( dcom2( chlw(i:i) ) .or. i .gt. i3 ) goto 350

               if( deqn1( chlw(i:i) ) ) then

                     isn = isn + 1
                     if( isn .eq. 1 ) icm = i
                     icn = i

               else if( chlw(i:i) .ge. 'a' .and.
     &                  chlw(i:i) .le. 'z' ) then

                     isa = isa + 1
                     if( isa .eq. 1 ) ica = i
                     icb = i

               end if

               goto 360

  350          continue

               icl = i

               if( ica .gt. 0 .and. icb - ica .gt. 0 ) goto 990
               if( icn .gt. 0 .and. icn - icm .gt. 1 ) goto 990
               if( ica .gt. 0 .and. icn .eq. 0 ) goto 990

               if( icn .eq. 0 ) then

                  ckam(1:1) = '0'
                  ckam(2:2) = '0'

               else if( icn .gt. 0 .and. icn .eq. icm ) then

                  ckam(1:1) = chlw(icn:icn)
                  ckam(2:2) = '0'

               else

                  ckam(1:1) = chlw(icm:icm)
                  ckam(2:2) = chlw(icn:icn)

               end if

               if( ica .gt. 0 ) then

                  ckam(3:3) = chlw(ica:ica)

               else if( icelp .eq. 3 ) then

                  ckam(3:3) = 'c'

               else if( icelp .eq. 4 ) then

                  ckam(3:3) = 'p'

               else if( icelp .eq. 5 ) then

                  ckam(3:3) = 'e'

               else if( icelp .eq. 7 ) then

                  ckam(3:3) = 'u'

               else if( icelp .eq. 8 ) then

                  ckam(3:3) = 'h'

               end if

               if( icelp .eq. 3 ) then

                  inlb = ichar(ckam(1:1)) + 256 * ( ichar(ckam(2:2))
     &                                    + 256 * ichar(ckam(3:3)) )

               else if( icelp .eq. 4 ) then

                  iplb = ichar(ckam(1:1)) + 256 * ( ichar(ckam(2:2))
     &                                    + 256 * ichar(ckam(3:3)) )

               else if( icelp .eq. 5 ) then

                  ielb = ichar(ckam(1:1)) + 256 * ( ichar(ckam(2:2))
     &                                    + 256 * ichar(ckam(3:3)) )

               else if( icelp .eq. 7 ) then

                  iulb = ichar(ckam(1:1)) + 256 * ( ichar(ckam(2:2))
     &                                    + 256 * ichar(ckam(3:3)) )

               else if( icelp .eq. 8 ) then

                  ihlb = ichar(ckam(1:1)) + 256 * ( ichar(ckam(2:2))
     &                                    + 256 * ichar(ckam(3:3)) )

               end if

            else if( icelp .eq. 6 ) then

               call snum(chlw,ic,i3,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 998

               icnd = nint( cvvv )

C ----------------------------------------------------------------------

            else if( icelp .eq. 9 ) then

               ic1 = ic

               do k = ic1, i3

                 if( chlw(k:k) .eq. ' ' .or. k .eq. i3 ) goto 8290

              end do

 8290         continue

              ic2 = k

              ch60 = chlw(ic1:ic2)
              call char_int(1,ch60,idedx)

               ic = jnumc(chlw,ic2+1,i3)
               icelp = 0

                goto 250
C ----------------------------------------------------------------------
C read chemical composition

            else if( icelp .eq. 10 ) then

              if(.not.allocated(ichem)) then
               call allocate_frac_ichem
              endif

              i   = 0

              ic1 = ic

              do while(ic1 .lt. i3)

  400         do k = ic1, i3

                 if( chlw(k+1:k+1) .eq. ' ' .or. k .eq. i3 ) goto 410

              end do

  410         continue

              ic2 = k

              ch200=chlw(ic1:ic2) ! T.Sato 2022/11/20
              call chemi_form(ch200,ic2-ic1+1,idchem,ierr)

              if(ierr .eq. 1) then

              write(ErrCha,*) 'This material does not exist in track-
     &structure molecule/atom library ', chlw(ic1:ic2)
              ErrID = 'L:23478/R:mater/F:read02.f'
              call ErrWrite(ErrID,ErrCha)
              call parastop( 500 )
              endif
              i = i + 1
              ichem(mxmat,i) = idchem

              ic1 = ic2

              do k = ic1, i3

                 if( chlw(k:k) .eq. ' ' .or. k .eq. i3 ) goto 420

              end do

  420         continue

              if(i .eq. 1 .and. k .eq. i3) then
                  frac(mxmat,1) = 1.d0 ! if material is 1 and density not specified, assume it is 100%
                  icl           = k + 1
                  exit
              endif

              ic = k

              call snum(chlw,ic,i3,icl,cvvv,ierr)
              if( ierr .ne. 0 ) goto 998

              frac(mxmat,i)  = cvvv

              ic1 = jnumc(chlw,icl,i3)

              enddo

            end if

               ic = jnumc(chlw,icl,i3) - 1
               icelp = 0

            goto 250

*-----------------------------------------------------------------------

  280    continue

*-----------------------------------------------------------------------
*        read s(a,b) information which is last part of one material
*-----------------------------------------------------------------------

         if( imts .ge. 0 ) then

  291          ic1 = ic

            do k = ic1, i3

               if( chlw(k:k) .eq. ' ' .or. k .eq. i3 ) goto 290

            end do

  290       continue

               ic2 = k
               if( ic2 - ic1 .gt. 10 ) goto 987

               imts = imts + 1

               hs = ' '
               ht = ' '
               hs = chlw(ic1:ic2)

               if( index(hs,'.') .eq. 0 .and. ic2 .lt. 9 ) then
                  ic2 = ic2 + 1
                  hs(ic2:ic2) = '.'
               else if( index(hs,'.') .eq. 0 ) then
                  goto 987
               end if

               ht(8-index(hs,'.'):10) = hs
               if( ht(8:8) .ne. ' ' .and. ht(9:9) .eq. ' ' )
     &         ht(9:9) = '0'

               call zaid(1,ht,ix)

               write(iom2) ix(1), ix(2), ix(3)

               ic = jnumc(chlw,ic2+1,i3)

               if( ic .le. i3 ) goto 291

               goto 140

         end if

*-----------------------------------------------------------------------
*        nucleus and density ( fraction )
*-----------------------------------------------------------------------

            ks = ks + 1

*-----------------------------------------------------------------------
*        read nucleus
*-----------------------------------------------------------------------

         if( iseq(ks) .eq. 1 ) then

               if( mxmat .eq. 0 ) goto 998

               call readnc(chin,chlw,ic,i3,icha,masi,libi,ierr)

                  if( ierr .ne. 0 ) goto 999
                  if( icha .gt. 104 )  goto 999

                  if( ks .eq. 2 ) goto 260
                  goto 250

*-----------------------------------------------------------------------
*        read density
*-----------------------------------------------------------------------

         else if( iseq(ks) .eq. 2 ) then

               if( mxmat .eq. 0 ) goto 998

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  denst = cvvv

                  ic = ic2 - 1

                  if( ks .eq. 2 ) goto 260
                  goto 250

         end if

*-----------------------------------------------------------------------
*        store the data
*-----------------------------------------------------------------------

  260    continue

               if( denst .ne. 0.0d0 ) then

                  if( isdn .eq. 0 ) isdn = nint(sign(1.0d0,denst))
                  if( isdn .ne. nint(sign(1.0d0,denst)) ) goto 995

               end if

            if( icha .eq. 1 .and. masi .eq. 1 ) then

!               ichn = ichn + 1
!               if( ichn .gt. 1 ) goto 989
               denh = denh + denst ! allow duplicate definition of 1H

               libh = libi

            else

               nel = nel + 1

               write(iom2) icha, masi, denst, libi

            end if

*-----------------------------------------------------------------------

            ks = 0
            goto 250

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

            if( ks .ne. 0 ) goto 993

            if( mxmat .gt. 0 ) then

                     if( imts .eq. -1 ) imts = 0

                     write(iom1) nel
                     write(iom1) denh, libh
                     write(iom1) igas, istp, inlb, iplb, ielb, icnd
                     write(iom1) iulb, ihlb
                     write(iom1) imts
                     write(iom1) (idedx(ii),ii=1,20)

               if( nel .gt. 0 .or. imts .gt. 0 ) rewind iom2

               if( nel .gt. 0 ) then

                  do i = 1, nel

                      read(iom2) icha, masi, denst, libi
                     write(iom1) icha, masi, denst, libi

                  end do

               end if

               if( imts .gt. 0 ) then

                  do i = 1, imts

                      read(iom2) ix(1), ix(2), ix(3)
                     write(iom1) ix(1), ix(2), ix(3)

                  end do

               end if

                     close( iom2 )

            end if

      return

*-----------------------------------------------------------------------

  987 continue

         m_err = 's(a,b) data name is wrong'
         ErrCha = ''
         ErrID = 'L:23703/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  988 continue

         m_err = 'number of mt is different from material number'
         ErrCha = ''
         ErrID = 'L:23715/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  989 continue

         m_err = '1H is appeared twice'
         ErrCha = ''
         ErrID = 'L:23727/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  990 continue

         m_err = 'value of parameter is wrong'
         ErrCha = ''
         ErrID = 'L:23739/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'value of parameter is missing'
         ErrCha = ''
         ErrID = 'L:23751/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'duplicate parameter in one material'
         ErrCha = ''
         ErrID = 'L:23763/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'odd number entries in [material]'
         ErrCha = ''
         ErrID = 'L:23775/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = '[material] section is duplicated'
         ErrCha = ''
         ErrID = 'L:23787/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'density should be positive or negative in a material'
         ErrCha = ''
         ErrID = 'L:23799/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Description of m### or mt### is wrong.'
         ErrCha = ''
         ErrID = 'L:23811/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of mat[###] is wrong.'
         ErrCha = ''
         ErrID = 'L:23823/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Description of [material] is wrong.'
         ErrCha = ''
         ErrID = 'L:23835/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Specification of Nucleus or Density is wrong.'
         ErrCha = ''
         ErrID = 'L:23847/R:mater/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine readnc(chin,chlw,ic,i3,icha,masi,libi,ierr)
*                                                                      *
*       read nucleus data                                              *
*       modified by K.Niita on 31/07/2000                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      character chin*200
      character chlw*200

      character ckam*3

      logical deqn1
      logical dcom2

*-----------------------------------------------------------------------

      character element(104)*3,cnuc*3

      data element/
     & 'h  ','he ','li ','be ','b  ','c  ','n  ','o  ',
     & 'f  ','ne ','na ','mg ','al ','si ','p  ','s  ',
     & 'cl ','ar ','k  ','ca ','sc ','ti ','v  ','cr ',
     & 'mn ','fe ','co ','ni ','cu ','zn ','ga ','ge ',
     & 'as ','se ','br ','kr ','rb ','sr ','y  ','zr ',
     & 'nb ','mo ','tc ','ru ','rh ','pd ','ag ','cd ',
     & 'in ','sn ','sb ','te ','i  ','xe ','cs ','ba ',
     & 'la ','ce ','pr ','nd ','pm ','sm ','eu ','gd ',
     & 'tb ','dy ','ho ','er ','tm ','yb ','lu ','hf ',
     & 'ta ','w  ','re ','os ','ir ','pt ','au ','hg ',
     & 'tl ','pb ','bi ','po ','at ','rn ','fr ','ra ',
     & 'ac ','th ','pa ','u  ','np ','pu ','am ','cm ',
     & 'bk ','cf ','es ','fm ','md ','no ','lr ','ku '/

*-----------------------------------------------------------------------

            ierr = 0

            ckam = ' '
            libi = ichar(' ') + 256 * ( ichar(' ')
     &                        + 256 * ichar(' ') )

*-----------------------------------------------------------------------
*        208Pb.24c  type ( 208-Pb is also OK )
*        Pb-208.24c type ( Pb208 is also OK )
*        82208.24c  type : iz * 1000 + ia + im * 0.01
*-----------------------------------------------------------------------

  300       continue

                  isa = 0
                  isn = 0
                  ipr = 0

                  icm = 0
                  icn = 0
                  ick = 0
                  icl = 0
                  ica = 0
                  icb = 0
                  icc = 0

               i = ic - 1

  360          i = i + 1

               if( dcom2( chlw(i:i) ) .or. i .gt. i3 ) goto 350

               if( chlw(i:i) .eq. '.' ) then

                     ipr = ipr + 1
                     if( ipr .gt. 1 ) goto 999

                     ick = i + 1

               else if( deqn1( chlw(i:i) ) ) then

                  if( ipr .eq. 0 ) then

                     isn = isn + 1
                     if( isn .eq. 1 ) icm = i

                     icn = i

                  else

                     icl = i

                  end if

               else if( chlw(i:i) .ge. 'a' .and.
     &                  chlw(i:i) .le. 'z' ) then

                  if( ipr .eq. 0 ) then

                     isa = isa + 1
                     if( isa .eq. 1 ) ica = i

                     icb = i

                  else

                     icc = icc + 1
                     if( icc .gt. 1 ) goto 999

                     ckam(3:3) = chlw(i:i)

                  end if

               else if( chlw(i:i) .ne. '-' ) then

                     goto 999

               end if

               goto 360

  350       continue

               ic = i

*-----------------------------------------------------------------------

            if( ica .gt. 0 .and. icb .lt. ica ) goto 999
            if( ica .gt. 0 .and. icb .gt. ica+1 ) goto 999
            if( icm .gt. 0 .and. icn .lt. icm ) goto 999
            if( icc .gt. 0 .and. ick .eq. 0 ) goto 999
            if( ick .gt. 0 .and. icl - ick .gt. 1 ) goto 999

*-----------------------------------------------------------------------
*           nucleus from character
*-----------------------------------------------------------------------

            if( isa .gt. 0 ) then

                  cnuc = chlw(ica:icb)//'  '

                  do j = 1, 104

                     if( cnuc(1:3) .eq. element(j)(1:3) ) goto 370

                  end do

                     goto 999

  370             icha = j

               if( icm .gt. 0 ) then

                  read(chlw(icm:icn),'(i5)') masi

               else

                  masi = 0

               end if

*-----------------------------------------------------------------------
*           nucleus from  82208
*-----------------------------------------------------------------------

            else

                  if( icm .eq. 0 ) goto 999

                  call onum(chlw,icm,icn,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  icha = int( cvvv / 1000.0 ) ! T.Sato 2023/04/30 change nint to int because mass number is over 500 for meta stable
                  masi = nint( cvvv - icha * 1000.0d0 )

            end if

c$$$            if(masi.gt.400) masi = masi - 400  ! T.Sato 2023/04/30 for meta-stable, mass number is +400
            if(masi.gt.400) then
              write(*,'(/"Warning :: Only neutrons and photons should ",
     &       "be transported when meta-stable target nuclei are ",
     &       "specified in [material].")')
            end if

*-----------------------------------------------------------------------
*           library number
*-----------------------------------------------------------------------

            if( ipr .gt. 0 ) then

               if( ick .gt. 0 .and. icl .eq. ick ) then

                  ckam(1:1) = chlw(ick:ick)
                  ckam(2:2) = '0'

               else

                  ckam(1:1) = chlw(ick:ick)
                  if(icl.ge.1) then
                   ckam(2:2) = chlw(icl:icl)
                  else
                   ckam = ' '  ! T.Sato 2024/12/11, avoid error when '.' is specified without library name
                  endif

               end if

               if( ckam(1:3) .eq. '00 ' .or.
     &             ckam(1:3) .eq. '000' )
     &             ckam = ' '

               libi = ichar(ckam(1:1)) + 256 * ( ichar(ckam(2:2))
     &                                 + 256 * ichar(ckam(3:3)) )


            end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

  999 continue

      ierr = 1

      return
      end

************************************************************************
*                                                                      *
      subroutine param(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [parameters] section of input files                       *
*       modified by K.Niita on 23/05/2000                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar
      common /parai/ ipsq(400)
      common /paraj/ mstz(300), parz(300)
      common /parak/ icnu, icdf(400), icdl(400), chnm(400)
      character chnm*8
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200

      common /userp/ idam(100), rdam(100)

      common /stat / istdev, irestart, ireschk
      common /stat2/ dmpmulti,idmpmode,ibchjmp,idmpjmp(2) !FURUTA20150515

*-----------------------------------------------------------------------

      integer*8 :: ibitrseed ! S.H. xorshift (2020.2.6)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension ihan(100) ! T.Sato 2016/06/04, =1: file parameter has been changed
      ihan(:)=0
*-----------------------------------------------------------------------

            ierr  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 235  ! goto 235 instead of return T.Sato 2018/08/28

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of parameters section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 235 ! goto 235 instead of return T.Sato 2018/08/28

            end if

*-----------------------------------------------------------------------
*        identify the parameters
*-----------------------------------------------------------------------

            icl = i1

  200    continue

            chlc = chlw
            call chcomp(chlc,icl,i3,i5)

         if( chlc(icl:icl+3) .eq. 'mstz' ) then

               impz = 1
               goto 300

         else if( chlc(icl:icl+3) .eq. 'parz' ) then

               impz = -1
               goto 300

         else if( chlc(icl:icl+3) .eq. 'idam' ) then

               imus = 1
               goto 500

         else if( chlc(icl:icl+3) .eq. 'rdam' ) then

               imus = -1
               goto 500

         else if( chlc(icl:icl+3) .eq. 'file' ) then

               goto 400

         else if( chlc(icl:icl+5) .eq. 'istdev' ) then

               goto 600

c S.H. xorshift (2020.2.6)
         else if( chlc(icl:icl+7) .eq. 'bitrseed' ) then

               goto 700

C Zenkaku Space Replace Log Flag
         else if( chlc(icl:icl+3) .eq. 'zspc' ) then
               ZSPC_Eflag = 1

         else if( chlc(icl:icl+2) .eq. 'lib' ) then

               goto 800

         else

            do i = 1, icnu

               il = icl + icdl(i) - 1

               if( chlc(icl:il) .eq. chnm(i)(1:icdl(i)) ) goto 100

            end do

               goto 999

         end if

*-----------------------------------------------------------------------
*        read value of named parameters
*-----------------------------------------------------------------------

  100    continue

               ic = inumc(chlw,il+1,i3,'=') + 1

               if( ic .gt. i3 ) goto 997

               icl = inumc(chlw,ic,i3,';') - 1

               ipm = i

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 998

*-----------------------------------------------------------------------
*           set values
*-----------------------------------------------------------------------

               if( icdf(ipm) .gt. 0 ) then

                  mstz( icdf(ipm) ) = nint( cvvv )

               else

                  parz( -icdf(ipm) ) = cvvv

               end if

                  ipara = ipara + 1
                  ipsq(ipara) = icdf(ipm)

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

               goto 140

*-----------------------------------------------------------------------
*        read value of non-named parameters
*-----------------------------------------------------------------------

  300    continue

               ic1 = inumc(chlw,icl+4,i3,'(')
               ic2 = inumc(chlw,ic1,i3,')')

               if( ic1 .gt. i3 .or. ic2 .gt. i3 ) goto 997

               call onum(chlw,ic1,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               imnm = nint( cvvv )

               ic = inumc(chlw,ic2+1,i3,'=') + 1

               if( ic .gt. i3 ) goto 997

               icl = inumc(chlw,ic,i3,';') - 1

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 998

*-----------------------------------------------------------------------
*           set values
*-----------------------------------------------------------------------

               if( impz .gt. 0 ) then

                  mstz( imnm ) = nint( cvvv )

               else

                  parz( imnm ) = cvvv

               end if

                  ipara = ipara + 1
                  ipsq(ipara) = impz * imnm

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

               goto 140

*-----------------------------------------------------------------------
*        read value of user defined parameters
*-----------------------------------------------------------------------

  500    continue

               ic1 = inumc(chlw,icl+4,i3,'(')
               ic2 = inumc(chlw,ic1,i3,')')

               if( ic1 .gt. i3 .or. ic2 .gt. i3 ) goto 997

               call onum(chlw,ic1,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               imnm = nint( cvvv )

               ic = inumc(chlw,ic2+1,i3,'=') + 1

               if( ic .gt. i3 ) goto 997

               icl = inumc(chlw,ic,i3,';') - 1

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 998

*-----------------------------------------------------------------------
*           set values
*-----------------------------------------------------------------------

               if( imus .gt. 0 ) then

                  idam( imnm ) = nint( cvvv )

               else

                  rdam( imnm ) = cvvv

               end if

                  ipara = ipara + 1
                  ipsq(ipara) = imus * ( imnm + 1000 )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

               goto 140

*-----------------------------------------------------------------------
*        read file name of file( ) = parameter
*-----------------------------------------------------------------------

  400    continue

               ic1 = inumc(chlw,icl+4,i3,'(')
               ic2 = inumc(chlw,ic1,i3,')')

               if( ic1 .gt. i3 .or. ic2 .gt. i3 ) goto 997

               call onum(chlw,ic1,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               imnm = nint( cvvv )
               ihan(imnm) = 1 ! T.Sato 2016/06/04, file parameter has been changed

               ic = inumc(chlw,ic2+1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 997

               ic1 = inumc(chlw,ic,i3,' ') - 1
               ic2 = inumc(chlw,ic,i3,';') - 1
               icl = min( ic1, ic2, i3 )

                  icfn( imnm ) = 1
                  ilfn( imnm ) = icl - ic + 1

               do k = ic, icl

                  chfn( imnm )(k-ic+1:k-ic+1) = chin(k:k)

               end do

               do k = icl-ic+2, 100

                  chfn( imnm )(k:k) = ' '

               end do

               icl = inumc(chlw,ic,i3,';') + 1
               icl = jnumc(chlw,icl,i3)

               if( icl .le. i3 ) goto 200

               goto 140

*-----------------------------------------------------------------------
*        read istdev
*-----------------------------------------------------------------------

  600    continue

               ic = inumc(chlw,icl,i3,'=') + 1

               if( ic .gt. i3 ) goto 997

               call onum(chlw,ic,i3,cvvv,ierr)

               istdev = nint( cvvv )
               if ( istdev .lt. 0 )then
                irestart   = 1
                if(idmpmode.gt.0.or.dmpmulti.ne.0.0) goto 996 !FURUTA20150515
               endif

               mstz(73) = istdev
               ipara = ipara + 1
               ipsq(ipara) = 73

               if ( all(istdev .ne. (/ -2, -1, 0, 1, 2 /)) ) goto 998
               if ( istdev .lt. 0 ) istdev = -istdev

               goto 140

*-----------------------------------------------------------------------
*        read bitrseed
*-----------------------------------------------------------------------

  700    continue

               ic = inumc(chlw,icl,i3,'=') + 1

               if( ic .gt. i3 ) goto 997

               read(chlw(ic+1:i3),'(b64)') ibitrseed
               parz(195) = transfer(ibitrseed,parz(195))

               ipara = ipara + 1
               ipsq(ipara) = -195

               goto 140

*-----------------------------------------------------------------------
*        read and set lib values ! frtati 2021/12/17
*-----------------------------------------------------------------------

  800    continue

            do ii = 1, icnu
               il = icl + icdl(ii) - 1
               if( chlc(icl:il) .eq. chnm(ii)(1:icdl(ii)) ) then
                 ipm = ii
                 exit
               end if
            end do

               ic = inumc(chlw,il+1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3) ! frtati 2022/03/24

               if( ic .gt. i3 ) goto 997

               icl = inumc(chlw,ic,i3,';') - 1


               mstz(icdf(ipm)) = ichar(chlw(ic:ic))
     &                          +256*( ichar(chlw(ic+1:ic+1))
     &                                +256*ichar(chlw(ic+2:ic+2)) )

               ipara = ipara + 1
               ipsq(ipara) = icdf(ipm)

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

               goto 140

*-----------------------------------------------------------------------
 996  continue

         m_err = 'istdev < 0 is not allowed'//
     &        ' when idmpmode = 1 or dmpmulti != 0.0.'
         ErrCha = ''
         ErrID = 'L:24511/R:param/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of parameter is wrong'
         ErrCha = ''
         ErrID = 'L:24523/R:param/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Value of parameter is wrong'
         ErrCha = ''
         ErrID = 'L:24535/R:param/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Unknown parameter'
         ErrCha = ''
         ErrID = 'L:24547/R:param/F:read02.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

! T.Sato 2016/06/04, revised 2018/08/28
  235 continue
      if(ihan(1).eq.0) then  ! file(1) is not defined, read from environmental variable
       call GET_ENVIRONMENT_VARIABLE('PHITSPATH', status=ios,
     & length = ilength)
       if(ios.eq.0) then ! PhitsPath exist
        ilfn(1)=ilength
        call GET_ENVIRONMENT_VARIABLE('PHITSPATH', value=chfn(1))
       endif
      endif

      if(ihan(7).ne.1) then
       chfn(7)=chfn(1)(1:ilfn(1))//'/data/xsdir.jnd'
       ilfn(7)=ilfn(1)+15
      endif
      if(ihan(20).ne.1) then
       chfn(20)=chfn(1)(1:ilfn(1))//'/XS/egs/'
       ilfn(20)=ilfn(1)+8
      endif
      if(ihan(21).ne.1) then
       chfn(21)=chfn(1)(1:ilfn(1))//'/dchain-sp/data/'
       ilfn(21)=ilfn(1)+16
      endif
      if(ihan(24).ne.1) then
       chfn(24)=chfn(1)(1:ilfn(1))//'/data/'
       ilfn(24)=ilfn(1)+6
      endif
      if(ihan(25).ne.1) then
       chfn(25)=chfn(1)(1:ilfn(1))//'/XS/tra/'
       ilfn(25)=ilfn(1)+8
      endif

      if(ihan(26).ne.1) then
       chfn(26)=chfn(1)(1:ilfn(1))//'/data/multiplier'
       ilfn(26)=ilfn(1)+16
      endif

      if(ihan(27).ne.1) then
       chfn(27)=chfn(1)(1:ilfn(1))//'/XS/yield/'
       ilfn(27)=ilfn(1)+10
      endif

      if(ihan(28).ne.1) then
       chfn(28)=chfn(1)(1:ilfn(1))//'/data/aama.dat'
       ilfn(28)=ilfn(1)+14
      endif

c dedxfile
      if(ihan(29).ne.1) then
       chfn(29)=chfn(1)(1:ilfn(1))//'/data/dedx'
       ilfn(29)=ilfn(1)+10
      endif

      return

      end


************************************************************************
*                                                                      *
      block data nmjamin
*                                                                      *
*       block data for input of PHITS                                  *
*       last modified by K.Niita and S. Hashimoto on 2011/08/02        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

      common /paraj/ mstz(300), parz(300)
      common /parak/ icnu, icdf(400), icdl(400), chnm(400)
      character chnm*8
      common /paral/ lpcn(400), lpcr(400), pcmn(400), pcmr(400)
      character pcmn*70, pcmr*70  !OBINATA(2012.6.13)
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*     Set file name
*-----------------------------------------------------------------------

      data icfn /100*0/

      data chfn( 1) /'c:/phits'/
      data ilfn( 1) / 8 /

      data chfn( 2) /'cgview.in'/
      data ilfn( 2) / 9 /

      data chfn( 3) /'cgview.set'/
      data ilfn( 3) / 10 /

      data chfn( 4) /'marspf.in'/
      data ilfn( 4) / 9 /

      data chfn( 6) /'phits.out'/
      data ilfn( 6) / 9 /

      data chfn( 7) /'c:/phits/data/xsdir.jnd'/
      data ilfn( 7) / 23 /

      data chfn(10) /'fort.10'/
      data ilfn(10) / 7 /

      data chfn(11) /'nuclcal.out'/
      data ilfn(11) / 11 /

      data chfn(12) /'fort.12'/
      data ilfn(12) / 7 /

      data chfn(13) /'fort.13'/
      data ilfn(13) / 7 /


      data chfn(15) /'dumpall.dat'/
      data ilfn(15) / 11 /

      data chfn(16) /'checkpt.dat'/
      data ilfn(16) / 11 /

      data chfn(17) /'restart.dat'/
      data ilfn(17) / 11 /

      data chfn(18) /'voxel.bin'/
      data ilfn(18) / 9 /

      data chfn(19) /'gcell.bin'/
      data ilfn(19) / 9 /

      data chfn(20) /'c:/phits/XS/egs'/
      data ilfn(20) /15/

      data chfn(21) /'c:/phits/dchain-sp/data'/
      data ilfn(21) /23/

      data chfn(22) /'batch.out'/
      data ilfn(22) /9/

      data chfn(23) /'pegs5'/
      data ilfn(23) /5/

      data chfn(24) /'c:/phits/data'/
      data ilfn(24) /13/

C T.Sato 2017/05/12, for track structure database
      data chfn(25) /'c:/phits/XS/tra'/
      data ilfn(25) /15/

      data chfn(26) /'c:/phits/data/multiplier'/
      data ilfn(26) /24/

      data chfn(27) /'c:/phits/XS/yield/'/
      data ilfn(27) /18/

      data chfn(28) /'c:/phits/data/aama.dat'/
      data ilfn(28) / 22 /

C dedxfile
      data chfn(29) /'c:/phits/data/dedx'/
      data ilfn(29) / 18 /

      data chfn(30) /'tetra.bin'/
      data ilfn(30) / 9 /

*-----------------------------------------------------------------------
*     Set character parameters
*-----------------------------------------------------------------------
cKN 2024/12/09

      data icnu /373/

*-----------------------------------------------------------------------
! T.Sato geometry memory extension factor
      data chnm(373)/'igeomem'/
      data icdf(373)/167/
      data icdl(373)/7/

*-----------------------------------------------------------------------
cFURUTA20241209 Option to dump rijk into rijkdmp.inp at each histroy
      data chnm(372)/'idmprijk'/
      data icdf(372)/166/
      data icdl(372)/8/

*-----------------------------------------------------------------------
c S.H. 2024.11.28 Option to display messages about ndata file in [t-yield]
      data chnm(371)/'iyldfile'/
      data icdf(371)/165/
      data icdl(371)/8/

*-----------------------------------------------------------------------
c Ogawa 2024/10/11  Rutherford bias angle upper boundary
      data chnm(370)/'ruthmax'/
      data icdf(370)/-206/
      data icdl(370)/7/

c Ogawa 2024/10/11  Rutherford bias angle lower boundary
      data chnm(369)/'ruthmin'/
      data icdf(369)/-205/
      data icdl(369)/7/

*-----------------------------------------------------------------------
c Ogawa 2024/10/11  T-product scores projectile as secondary ?
      data chnm(368)/'iprjprd'/
      data icdf(368)/164/
      data icdl(368)/7/

*-----------------------------------------------------------------------
cKN 2024/03/26  restriction WW parameters for huge bank particle
      data chnm(367)/'iwwbnk'/
      data icdf(367)/163/
      data icdl(367)/6/

*-----------------------------------------------------------------------
c T.Sato 2024/03/24, small gap for avoid integer value in weight window mesh
      data chnm(366)/'deltxyz'/
      data icdf(366)/-204/
      data icdl(366)/7/

*-----------------------------------------------------------------------
c T.Sato 2024/02/18 number of segmentation in charged-particle forced collisions
      data chnm(365)/'nfcseg'/
      data icdf(365)/162/
      data icdl(365)/6/

*-----------------------------------------------------------------------
c S.H. 2023.10.27 number of c-values output
      data chnm(364)/'ncvalout'/
      data icdf(364)/161/
      data icdl(364)/8/

*-----------------------------------------------------------------------
      data chnm(363)/'tsxcl'/
      data icdf(363)/160/
      data icdl(363)/5/

*-----------------------------------------------------------------------
      data chnm(362)/'iemdebug'/
      data icdf(362)/159/
      data icdl(362)/8/

*-----------------------------------------------------------------------
      data chnm(361)/'itgchk'/
      data icdf(361)/158/
      data icdl(361)/6/

*-----------------------------------------------------------------------
      data chnm(360)/'mdbpseud'/
      data icdf(360)/157/
      data icdl(360)/8/

! JI Marquez Damian 2021/10/21
      data chnm(359)/'isans'/
      data icdf(359)/156/
      data icdl(359)/5/

      data chnm(358)/'dpnmax'/
      data icdf(358)/-203/
      data icdl(358)/6/
      data (chnm(i),i=351,357)
     &  /'lib(1)','lib(2)','lib(14)','lib(15)','lib(16)'
     &  ,'lib(17)','lib(18)'/
      data (icdf(i),i=351,357)/149,150,151,152,153,154,155/
      data (icdl(i),i=351,357)/6,6,7,7,7,7,7/

      data chnm(350)/'rtrckflp'/
      data icdf(350)/-202/
      data icdl(350)/8/

      data chnm(349)/'xsmemory'/
      data icdf(349)/-201/
      data icdl(349)/8/

      data chnm(348)/'icxnp'/
      data icdf(348)/148/
      data icdl(348)/5/

      data chnm(347)/'ierrout'/
      data icdf(347)/147/
      data icdl(347)/7/

      data chnm(346)/'scm_rcls'/
      data icdf(346)/-200/
      data icdl(346)/8/

      data chnm(345)/'scm_d'/
      data icdf(345)/-199/
      data icdl(345)/5/

      data chnm(344)/'scm_h0'/
      data icdf(344)/-198/
      data icdl(344)/6/

      data chnm(343)/'iqmdscm'/
      data icdf(343)/146/
      data icdl(343)/7/

      data chnm(342)/'epseudo'/
      data icdf(342)/-197/
      data icdl(342)/7/

      data chnm(341)/'ichkmat'/
      data icdf(341)/145/
      data icdl(341)/7/

      data chnm(340)/'italsh'/
      data icdf(340)/144/
      data icdl(340)/6/

      data chnm(339)/'ifission'/
      data icdf(339)/143/
      data icdl(339)/8/

      data chnm(338)/'icells'/
      data icdf(338)/142/
      data icdl(338)/6/

      data chnm(337)/'tsmax'/
      data icdf(337)/-196/
      data icdl(337)/5/

      data chnm(336)/'iscinful'/
      data icdf(336)/141/
      data icdl(336)/8/

c S.H. xorshift (2020.2.6)
      data chnm(335)/'itimrand'/
      data icdf(335)/140/
      data icdl(335)/8/

      data chnm(334)/'bitrseed'/
      data icdf(334)/-195/
      data icdl(334)/8/

      data chnm(333)/'nrandgen'/
      data icdf(333)/139/
      data icdl(333)/8/

*-----------------------------------------------------------------------
* added by nais 2019.11.21
      data chnm(332)/'gsline'/
      data icdf(332)/138/
      data icdl(332)/6/

*-----------------------------------------------------------------------
      data chnm(331)/'gmumul'/
      data icdf(331)/-194/
      data icdl(331)/6/

*-----------------------------------------------------------------------
      data chnm(330)/'igmuppd'/
      data icdf(330)/137/
      data icdl(330)/7/

*-----------------------------------------------------------------------
      data chnm(329)/'istdbat'/
      data icdf(329)/136/
      data icdl(329)/7/

*-----------------------------------------------------------------------
      data chnm(328)/'italecho'/
      data icdf(328)/135/
      data icdl(328)/8/

*-----------------------------------------------------------------------
      data chnm(327)/'istdcut'/
      data icdf(327)/134/
      data icdl(327)/7/

*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      data chnm(326)/'itetauto'/
      data icdf(326)/133/
      data icdl(326)/8/

*-----------------------------------------------------------------------
      data chnm(325)/'ifixchg'/
      data icdf(325)/132/
      data icdl(325)/7/

*-----------------------------------------------------------------------
      data chnm(324)/'ntrnore'/
      data icdf(324)/131/
      data icdl(324)/7/

*-----------------------------------------------------------------------
! T.Sato 2018/03/06 for Kerma control
      data (chnm(i),i=322,323)/'ikerman','ikermap'/
      data (icdf(i),i=322,323)/129,130/
      data (icdl(i),i=322,323)/7,7/

*-----------------------------------------------------------------------
      data chnm(321)/'em-emode'/
      data icdf(321)/-193/
      data icdl(321)/8/

*-----------------------------------------------------------------------
      data chnm(320)/'nucdata'/
      data icdf(320)/128/
      data icdl(320)/7/

*-----------------------------------------------------------------------

      data chnm(319)/'irlet'/
      data icdf(319)/127/
      data icdl(319)/5/

*-----------------------------------------------------------------------
! T.Sato 2017/05/11 for Track Structure
      data (chnm(i),i=317,318)/'etsmin','etsmax'/
      data (icdf(i),i=317,318)/-191,-192/
      data (icdl(i),i=317,318)/6,6/

*-----------------------------------------------------------------------

      data chnm(316)/'ngem'/
      data icdf(316)/126/
      data icdl(316)/4/

*-----------------------------------------------------------------------

      data chnm(315)/'timeout'/
      data icdf(315)/-190/
      data icdl(315)/7/

*-----------------------------------------------------------------------

      data chnm(314)/'iwwbias'/
      data icdf(314)/125/
      data icdl(314)/7/

*-----------------------------------------------------------------------
c T.Sato 2016.12.28 commnet in [material] section
      data chnm(313)/'icommat'/
      data icdf(313)/124/
      data icdl(313)/7/

*-----------------------------------------------------------------------
cMIURA 2016.09.30 tally output unit, the energy per nucleon. [MeV/n]
      data chnm(312)/'imevperu'/
      data icdf(312)/123/
      data icdl(312)/8/

*-----------------------------------------------------------------------
      data (chnm(i),i=309,311)/'itetra','ntetsurf','ntetelem'/
      data (icdf(i),i=309,311)/120,121,122/
      data (icdl(i),i=309,311)/6,8,8/

*-----------------------------------------------------------------------

      data chnm(308)/'icrdm'/
      data icdf(308)/119/
      data icdl(308)/5/

*-----------------------------------------------------------------------

      data chnm(307)/'natural'/
      data icdf(307)/118/
      data icdl(307)/7/

*-----------------------------------------------------------------------

      data chnm(306)/'adjemax'/
      data icdf(306)/-189/
      data icdl(306)/7/

*-----------------------------------------------------------------------
      data chnm(305)/'icxspi  '/
      data icdf(305)/117/
      data icdl(305)/6/

*-----------------------------------------------------------------------
      data chnm(304)/'prmui1  '/
      data icdf(304)/-188/
      data icdl(304)/6/

*-----------------------------------------------------------------------
      data chnm(303)/'ismode  '/
      data icdf(303)/116/
      data icdl(303)/6/

*-----------------------------------------------------------------------
cFURUTA20160126 ATIMA database by Wada

      data (chnm(i),i=301,302)/'mdbatima','dbcutoff'/
      data (icdf(i),i=301,302)/115,-187/
      data (icdl(i),i=301,302)/8,8/

*-----------------------------------------------------------------------

      data chnm(300)/'iadjoint'/
      data icdf(300)/114/
      data icdl(300)/8/

      data chnm(299)/'nonu    '/
      data icdf(299)/113/
      data icdl(299)/4/

*-----------------------------------------------------------------------
! T.Sato 2015/08/31, consideration of multiple scattering in EGS5
      data chnm(298)/'imsegs  '/
      data icdf(298)/112/
      data icdl(298)/6/

*-----------------------------------------------------------------------
      data chnm(297)/'idelt  '/
      data icdf(297)/111/
      data icdl(297)/5/

*-----------------------------------------------------------------------
      data chnm(296)/'pnimul  '/
      data icdf(296)/-186/
      data icdl(296)/6/

*-----------------------------------------------------------------------
      data (chnm(i),i=294,295)/'imubrm  ','imuppd  '/
      data (icdf(i),i=294,295)/109,110/
      data (icdl(i),i=294,295)/6,6/

*-----------------------------------------------------------------------
cSato PEGS mode, chard
      data (chnm(i),i=290,293)/'ipegs','chard','iunrst','epstfl'/
      data (icdf(i),i=290,293)/106,-185,107,108/
      data (icdl(i),i=290,293)/5,5,6,6/

*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      data chnm(289)/'itetvol '/
      data icdf(289)/105/
      data icdl(289)/7/

*-----------------------------------------------------------------------
cEGS
      data (chnm(i),i=281,288)/
     &    'ibound  ','iaprim  ','iegsout ','gasegs  ','nrecover',
     &    'ascat1  ','ascat2  ','imucap  '/

      data (icdf(i),i=281,288)/
     &    100,       101,       102,      -182,       103,
     &   -183,      -184,       104/

      data (icdl(i),i=281,288)/
     &      6,         6,         7,         6,       8,
     &      6,         6,         6/

*-----------------------------------------------------------------------

      data (chnm(i),i=261,280)/
     &    'infout  ','mvoww   ','imuint  ','emumin  ','emumax  ',
     &    'negs    ','iedgfl  ','iauger  ','iraylr  ','lpolar  ',
     &    'incohr  ','iprofr  ','impacr  ','iegsrand','luxlev  ',
     &    'ieispl  ','neispl  ','ibrdst  ','iprdst  ','iphter  '/

      data (icdf(i),i=261,280)/
     &     82,        83,        84,      -180,      -181,
     &     85,        86,        87,        88,        89,
     &     90,        91,        92,        93,        94,
     &     95,        96,        97,        98,        99/

      data (icdl(i),i=261,280)/
     &      6,         5,         6,         6,         6,
     &      4,         6,         6,         6,         6,
     &      6,         6,         6,         8,         6,
     &      6,         6,         6,         6,         6/

*-----------------------------------------------------------------------

      data (chnm(i),i=241,260)/
     &    'isaba   ','ivoxel  ','ipnint  ','ielctf  ','inclg   ',
     &    'inclv   ','istdev  ','ismm    ','dsck    ','icxsni  ',
     &    'incelf  ','idwba   ','einclmin','einclmax','eielfmin',
     &    'eielfmax','irndmode','irqmd   ','ireschk ','ifbm    '/

      data (icdf(i),i=241,260)/
     &     67,        68,        69,        70,        71,
     &     72,        73,        74,      -175,        75,
     &     76,        77,      -176,      -177,      -178,
     &   -179,        78,        79,        80,        81/

      data (icdl(i),i=241,260)/
     &      5,         6,         6,         6,         5,
     &      5,         6,         4,         4,         6,
     &      6,         5,         8,         8,         8,
     &      8,         8,         5,         7,         4/

*-----------------------------------------------------------------------

      data (chnm(i),i=221,240)/
     &    'icrhi   ','bplus   ','nwsors  ','judge   ','imadj   ',
     &    'ejamqmd ','ndedx   ','esmax   ','deltt   ','e-mode  ',
     &    'ge1     ','ge2     ','ih2o    ','dumpall ','idpara  ',
     &    'usrmgt  ','usrelst ','nqtmax  ','c9decay ','iidfs   '/

      data (icdf(i),i=221,240)/
     &     54,      -168,        55,        56,        57,
     &   -169,        58,      -170,      -171,        59,
     &   -172,      -173,      -174,        60,        61,
     &     62,        63,        64,        65,        66/

      data (icdl(i),i=221,240)/
     &      5,         5,         6,         5,         5,
     &      7,         5,         5,         5,         6,
     &      3,         3,         4,         7,         6,
     &      6,         7,         6,         7,         5/

*-----------------------------------------------------------------------

      data (chnm(i),i=201,220)/
     &    'dbcn(17)','dbcn(18)','icput   ','ipara   ','ieleh   ',
     &    'wupn    ','wsurvn  ','mxspln  ','mwhere  ','dircha  ',
     &    'pngdr   ','esmin   ','nedisp  ','deltc   ','eqmdmin ',
     &    'itstep  ','lionprd ','gravx   ','gravy   ','gravz   '/
cKN 2016/08/08 nkerma->kerma, 2025.3.6. S.H. changed kerma->lionprd

      data (icdf(i),i=201,220)/
     &   -158,      -159,        44,        45,        46,
     &   -160,      -161,        47,        48,        49,
     &     50,      -162,        51,      -163,      -164,
     &     52,        53,      -165,      -166,      -167/

      data (icdl(i),i=201,220)/
     &      8,         8,         5,         5,         5,
     &      4,         6,         6,         6,         6,
     &      5,         5,         6,         5,         7,
     &      6,         7,         5,         5,         5/

*-----------------------------------------------------------------------

      data (chnm(i),i=181,200)/
     &    'dmax(18)','dmax(19)','dmax(20)','maxbnk  ','jmout   ',
     &    'emcnf   ','iunr    ','dnb     ','ides    ','nocoh   ',
     &    'iphot   ','ibad    ','istrg   ','bnum    ','xnum    ',
     &    'rnok    ','enum    ','numb    ','emcpf   ','kmout   '/

      data (icdf(i),i=181,200)/
     &   -148,      -149,      -150,        34,        35,
     &   -151,        36,      -152,        37,        38,
     &     39,        40,        41,      -153,      -154,
     &   -155,      -156,        42,      -157,        43/

      data (icdl(i),i=181,200)/
     &      8,         8,         8,         6,         5,
     &      5,         4,         3,         4,         5,
     &      5,         4,         5,         4,         4,
     &      4,         4,         4,         5,         5/

*-----------------------------------------------------------------------

      data (chnm(i),i=161,180)/
     &    'cmin(18)','cmin(19)','cmin(20)','dmax(1) ','dmax(2) ',
     &    'dmax(3) ','dmax(4) ','dmax(5) ','dmax(6) ','dmax(7) ',
     &    'dmax(8) ','dmax(9) ','dmax(10)','dmax(11)','dmax(12)',
     &    'dmax(13)','dmax(14)','dmax(15)','dmax(16)','dmax(17)'/

      data (icdf(i),i=161,180)/
     &   -128,      -129,      -130,      -131,      -132,
     &   -133,      -134,      -135,      -136,      -137,
     &   -138,      -139,      -140,      -141,      -142,
     &   -143,      -144,      -145,      -146,      -147/

      data (icdl(i),i=161,180)/
     &      8,         8,         8,         7,         7,
     &      7,         7,         7,         7,         7,
     &      7,         7,         8,         8,         8,
     &      8,         8,         8,         8,         8/

*-----------------------------------------------------------------------

      data (chnm(i),i=141,160)/
     &    'iptall  ','deltg   ','eqmdnu  ','cmin(1) ','cmin(2) ',
     &    'cmin(3) ','cmin(4) ','cmin(5) ','cmin(6) ','cmin(7) ',
     &    'cmin(8) ','cmin(9) ','cmin(10)','cmin(11)','cmin(12)',
     &    'cmin(13)','cmin(14)','cmin(15)','cmin(16)','cmin(17)'/

      data (icdf(i),i=141,160)/
     &     33,      -109,      -110,      -111,      -112,
     &   -113,      -114,      -115,      -116,      -117,
     &   -118,      -119,      -120,      -121,      -122,
     &   -123,      -124,      -125,      -126,      -127/

      data (icdl(i),i=141,160)/
     &      6,         5,         6,         7,         7,
     &      7,         7,         7,         7,         7,
     &      7,         7,         8,         8,         8,
     &      8,         8,         8,         8,         8/

*-----------------------------------------------------------------------

      data (chnm(i),i=121,140)/
     &    'wc2(7)  ','wc2(8)  ','wc2(9)  ','wc2(10) ','wc2(11) ',
     &    'wc2(12) ','wc2(13) ','wc2(14) ','wc2(15) ','wc2(16) ',
     &    'wc2(17) ','wc2(18) ','wc2(19) ','wc2(20) ','matadd  ',
     &    'inpara  ','igpara  ','ipcut   ','ippara  ','iggcm'/

      data (icdf(i),i=121,140)/
     &    -95,       -96,       -97,       -98,       -99,
     &   -100,      -101,      -102,      -103,      -104,
     &   -105,      -106,      -107,      -108,        27,
     &     28,        29,        30,        31,        32/

      data (icdl(i),i=121,140)/
     &      6,         6,         6,         7,         7,
     &      7,         7,         7,         7,         7,
     &      7,         7,         7,         7,         6,
     &      6,         6,         5,         6,         5/

*-----------------------------------------------------------------------

      data (chnm(i),i=101,120)/
     &    'wc1(7)  ','wc1(8)  ','wc1(9)  ','wc1(10) ','wc1(11) ',
     &    'wc1(12) ','wc1(13) ','wc1(14) ','wc1(15) ','wc1(16) ',
     &    'wc1(17) ','wc1(18) ','wc1(19) ','wc1(20) ','wc2(1)  ',
     &    'wc2(2)  ','wc2(3)  ','wc2(4)  ','wc2(5)  ','wc2(6)  '/

      data (icdf(i),i=101,120)/
     &    -75,       -76,       -77,       -78,       -79,
     &    -80,       -81,       -82,       -83,       -84,
     &    -85,       -86,       -87,       -88,       -89,
     &    -90,       -91,       -92,       -93,       -94/

      data (icdl(i),i=101,120)/
     &      6,         6,         6,         7,         7,
     &      7,         7,         7,         7,         7,
     &      7,         7,         7,         7,         6,
     &      6,         6,         6,         6,         6/

*-----------------------------------------------------------------------

      data (chnm(i),i=81,100)/
     &    'swtm(7) ','swtm(8) ','swtm(9) ','swtm(10)','swtm(11)',
     &    'swtm(12)','swtm(13)','swtm(14)','swtm(15)','swtm(16)',
     &    'swtm(17)','swtm(18)','swtm(19)','swtm(20)','wc1(1)  ',
     &    'wc1(2)  ','wc1(3)  ','wc1(4)  ','wc1(5)  ','wc1(6)  '/

      data (icdf(i),i=81,100)/
     &    -55,       -56,       -57,       -58,       -59,
     &    -60,       -61,       -62,       -63,       -64,
     &    -65,       -66,       -67,       -68,       -69,
     &    -70,       -71,       -72,       -73,       -74/

      data (icdl(i),i=81,100)/
     &      7,         7,         7,         8,         8,
     &      8,         8,         8,         8,         8,
     &      8,         8,         8,         8,         6,
     &      6,         6,         6,         6,         6/

*-----------------------------------------------------------------------

      data (chnm(i),i=61,80)/
     &    'tmax(7) ','tmax(8) ','tmax(9) ','tmax(10)','tmax(11)',
     &    'tmax(12)','tmax(13)','tmax(14)','tmax(15)','tmax(16)',
     &    'tmax(17)','tmax(18)','tmax(19)','tmax(20)','swtm(1) ',
     &    'swtm(2) ','swtm(3) ','swtm(4) ','swtm(5) ','swtm(6) '/

      data (icdf(i),i=61,80)/
     &    -35,       -36,       -37,       -38,       -39,
     &    -40,       -41,       -42,       -43,       -44,
     &    -45,       -46,       -47,       -48,       -49,
     &    -50,       -51,       -52,       -53,       -54/

      data (icdl(i),i=61,80)/
     &      7,         7,         7,         8,         8,
     &      8,         8,         8,         8,         8,
     &      8,         8,         8,         8,         7,
     &      7,         7,         7,         7,         7/

*-----------------------------------------------------------------------

      data (chnm(i),i=41,60)/
     &    'imout   ','incut   ','igcut   ','ivout   ','ipout   ',
     &    'delt0   ','nlost   ','igerr   ','deltm   ','igchk   ',
     &    'deltb   ','ielms   ','inucl   ','itall   ','tmax(1) ',
     &    'tmax(2) ','tmax(3) ','tmax(4) ','tmax(5) ','tmax(6) '/

      data (icdf(i),i=41,60)/
     &     16,        17,        18,        19,        20,
     &    -26,        21,        22,       -27,        23,
     &    -28,        24,        25,        26,       -29,
     &    -30,       -31,       -32,       -33,       -34/

      data (icdl(i),i=41,60)/
     &      5,         5,         5,         5,         5,
     &      5,         5,         5,         5,         5,
     &      5,         5,         5,         5,         7,
     &      7,         7,         7,         7,         7/

*-----------------------------------------------------------------------

      data (chnm(i),i=21,40)/
     &    'emin(19)','emin(20)','maxcas  ','maxbch  ','npidk   ',
     &    'andit   ','nevap   ','nspred  ','ielas   ','isobar  ',
     &    'ipreeq  ','level   ','igamma  ','inmed   ','imagnf  ',
     &    'ejamnu  ','ejampi  ','eisobar ','icntl   ','inucr   '/

      data (icdf(i),i=21,40)/
     &    -19,       -20,         3,         4,         5,
     &    -22,         6,         7,         8,         9,
     &     10,        11,        12,        13,        14,
     &    -23,       -24,       -25,         1,        15/

      data (icdl(i),i=21,40)/
     &      8,         8,         6,         6,         5,
     &      5,         5,         6,         5,         6,
     &      6,         5,         6,         5,         6,
     &      6,         6,         7,         5,         5/

*-----------------------------------------------------------------------

      data (chnm(i),i=1,20)/
     &    'rseed   ','irskip  ','emin(1) ','emin(2) ','emin(3) ',
     &    'emin(4) ','emin(5) ','emin(6) ','emin(7) ','emin(8) ',
     &    'emin(9) ','emin(10)','emin(11)','emin(12)','emin(13)',
     &    'emin(14)','emin(15)','emin(16)','emin(17)','emin(18)'/

      data (icdf(i),i=1,20)/
     &    -21,         2,        -1,        -2,        -3,
     &     -4,        -5,        -6,        -7,        -8,
     &     -9,       -10,       -11,       -12,       -13,
     &    -14,       -15,       -16,       -17,       -18/

      data (icdl(i),i=1,20)/
     &      5,         6,         7,         7,         7,
     &      7,         7,         7,         7,         7,
     &      7,         8,         8,         8,         8,
     &      8,         8,         8,         8,         8/

*-----------------------------------------------------------------------
*     Set defaults in the mstz array
*-----------------------------------------------------------------------
*        1     2     3     4     5     6     7     8     9    10
      data (mstz(i),i=101,200)/
     $   1,    0,    0,    1,    0,    0,    0,    0,    1,    1,
     1   1,    1,    1,    0,  500,    0,    1,    1,    0,    0,
     2 100,  200,    0,    0,    0,    0,    1,    1,    0,    0,
     3   0,    0,    0,    0,    1,    0,    0,    2,    1,    0,
     4   0,    3,    0,    0,    0,    0,    0,    1,
     & 6828082, 6500402, 7680050, 7286834, 7483701, 7549237, 6369330, ! frtati 2023/08/17 JENDL5 extensions
     5                                 0,  200,    0,    0,    0,
     6   0,   50,    1,    0,    0,    0,10000,    0,    0,    0,
     7   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     8   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     9   0,    0,    0,    0,    0,    0,    0,    0,    0,    0/

*-----------------------------------------------------------------------
cABE 2018/01/26: change 0 to 1 at index 52(itstep).
*OBINATA(2012.10.2): change 1 to 0 at index 73(istdev).
*        1     2     3     4     5     6     7     8     9    10
      data (mstz(i),i=1,100)/
     $   0,    0,   10,   10,    0,    3,    0,    2,    0,    0,
     1   3,    2,    1,    0,    1,    0,    0,    0,    0,    1,
     2  10,    1,    0,  720,    0,    0,    1,    0,    0,    0,
     3   0,    0,    0,10000,    0,    0,    1,    0,    0,    0,
     4   0,    0,    0,    0,    0,    0,    5,    0,    0,    0,
     5   0,    1,   -1,    2,    0,    1,    1,    3,    0,    0,
     6   3,    1,    1,  150,    0,    0,    0,    0,    0,    0,
     7   1,    0,    0,    0,    0,    0,    0,    0,    0,    1,
     8   0,    7,    0,    1,   -1,    1,    1,    1,    0,    1,
     9   1,    1,   -1,    1,    0,    0,    1,    1,    1,    1/

*-----------------------------------------------------------------------
*     Set defaults in the parz array
*-----------------------------------------------------------------------

      data parz(204) /1.d-4/ ! T.Sato 2024/03/24, deltxyz
      data parz(203) /1.d-3/ ! frtati 2021/12/17 for dpnmax
      data parz(202) /0.2d0/
      data parz(201) /1.0d0/ ! after 201, set default parameter by each parameter

*-----------------------------------------------------------------------
*         1       2       3       4       5       6       7       8
*         9      10
      data (parz(i),i=101,200)/
     & -100.d0,-100.d0,-100.d0,-100.d0,-100.d0,-100.d0,-100.d0,-100.d0,
     &1.012345d0,20.d0,
     1    0.d0,   0.d0,   0.d0,   0.d0,   0.d0,   0.d0,   0.d0,   0.d0,
     &    0.d0,   0.d0,
     2    0.d0,   0.d0,   0.d0,   0.d0,   0.d0,   0.d0,   0.d0,   0.d0,
     &    0.d0,   0.d0,
     3    0.d0,  20.d0,   0.d0,   0.d0,   0.d0,   0.d0,   0.d0,   0.d0,
     &    0.d0,   0.d0,
     4    0.d0,   0.d0,   0.d0,   1.d3,   0.d0,   0.d0,   0.d0,   0.d0,
     &    0.d0,   0.d0,
     5    0.d0,  -1.d0,   1.d0,   1.d0,   1.d0,   1.d0, 100.d0,   0.d0,
     &    0.d0,   5.d0,
     6 -100.d0,  1.d-3,2.012345d0,10.d0,  0.d0,   0.d0,   0.d0,   0.d0,
     & 3000.d0,   1.d6,
     7    1.d0,   1.d0,   1.d0,  -1.d0,   0.d0,   1.d0,3000.d0,   1.d0,
     & 3500.d0, 200.d0,
     8    1.d6, 0.03d0, 13.6d0, 3.8d-2,  1.d-1,   1.d0,   0.d0,  5.d-2,
cKN 2016/12/31, T.Sato 2017/1/31 aspara2 = 0.088 -> 0.038
cKN 2017/02/15, parz(189)=3.0 for adjoint mode
     &    3.d0,  -1.d0,
     9   1.d-6,  1.d-2,  20.d0, 100.d0,   0.d0,  1.d-3,  10.d0, 0.25d0,
cKN 2018/02/12           em-emode=20MeV
! T.Sato 2020/03/18 tsmax=1d-3, will be changed to 1d8 in future
     &  1.75d0,  -1.124d0 /

*-----------------------------------------------------------------------
*         1       2       3       4       5       6       7       8
*         9      10
      data (parz(i),i=1,100)/
     $  1.0d-3, 1.0d-11, 1.0d-3, 1.0d-3, 1.0d0,  1.0d-3, 1.0d-3, 1.0d-3,
     &  1.0d-3, 1.0d-3,
     1   1.0d0,  1.0d9,  1.0d9,  1.0d-3, 1.0d-3, 1.0d-3, 1.0d-3, 1.0d-3,
     &  1.0d-3,  1.0d9,
     2   0.0d0,  0.0d0,  20.d0,  20.d0,  0.0d0,  0.1d0,  20.12345d0,
     &   1.d-09, 1.0d9,  1.0d9,
     3   1.0d9,  1.0d9,  1.0d9,  1.0d9,  1.0d9,  1.0d9,  1.0d9,  1.0d9,
     &   1.0d9,  1.0d9,
     4   1.0d9,  1.0d9,  1.0d9,  1.0d9,  1.0d9,  1.0d9,  1.0d9,  1.0d9,
     &    1.d0,   1.d0,
     5    1.d0,   1.d0,   1.d0,   1.d0,   1.d0,   1.d0,   1.d0,   1.d0,
     &    1.d0,   1.d0,
     6    1.d0,   1.d0,   1.d0,   1.d0,   1.d0,   1.d0,   1.d0,   1.d0,
     &    -0.5,   -0.5,
     7    -0.5,   -0.5,   -0.5,   -0.5,   -0.5,   -0.5,   -0.5,   -0.5,
     &    -0.5,    0.0,
     8     0.0,   -0.5,   -0.5,   -0.5,   -0.5,   -0.5,   -0.5,   -0.5,
     & -100.d0,-100.d0,
     9 -100.d0,-100.d0,-100.d0,-100.d0,-100.d0,-100.d0,-100.d0,-100.d0,
     & -100.d0,-100.d0/

*-----------------------------------------------------------------------
*     Set comments for mstz
*-----------------------------------------------------------------------

*           5    1    5    2    5    3    5    4    5    5    5    6
! T.Sato 2024/12/25
      data pcmn(167)
     &/'(D=10000) Geometry memory extension size'/
      data lpcn(167)/ 40 /

cFURUTA Option to dump rijk into rijkdmp.inp at each history
      data pcmn(166)
     &/'(D=0) Option to dump rijk into rijkdmp.inp at each history'/
      data lpcn(166)/ 58 /

c S.H. 2024.11.28 Option to display messages about ndata file in [t-yield]
      data pcmn(165)
     &/'(D=0) Option to display messages about ndata file in [t-yield]'/
      data lpcn(165)/ 62 /

c Ogawa 2024/10/11  flag to score elastic reaction products in Track structure by t-product
      data pcmn(164)
     &/'(D=0)  flag to score projectile Rutherford recoil by t-product'/
      data lpcn(164)/ 62 /

cKN 2024/03/26  restriction WW parameters for huge bank particle
      data pcmn(163)
     &/'(D=1) Control WW parameter based on maxbnk'/
      data lpcn(163)/ 42 /

c T.Sato 2024/02/18 number of segmentation in charged-particle forced collision
      data pcmn(162)
     &/'(D=50) Number of segmentation in forced collision'/
      data lpcn(162)/ 49 /

c S.H. 2023.10.27 number of c-values output
      data pcmn(161)
     &/'(D=0) Number of c-values output'/
      data lpcn(161)/ 31 /

      data pcmn(160)
     &/'(D=1) Track structure exclusive mode. 0:on, 1:off'/ ! In the path where KURBUC cross section is non-zero, other cross sections are nullfied.
      data lpcn(160)/ 49 /
      data pcmn(159)
     &/'(D=0) mark & markp reset for egs magnetic field, 0:yes 1: no'/
      data lpcn(159)/ 60 /

      data pcmn(158)
     &/'(D=0) Tetra-mesh geometry check 0: off, 1: on'/
      data lpcn(158)/ 50 /
      data pcmn(157)
     &/'(D=200) Size of cross-section database for pseudo collision'/
      data lpcn(157)/ 59 /
      data pcmn(156)
     &/'(D=0) SANS reaction 0: off, 1: on'/
      data lpcn(156)/ 33 /
      data pcmn(155)
     &/'(D=51a) Extension for alpha library'/
      data lpcn(155)/ 35 /
      data pcmn(154)
     &/'(D=51s) Extension for He-3 library'/
      data lpcn(154)/ 34 /
      data pcmn(153)
     &/'(D=51r) Extension for triton library'/
      data lpcn(153)/ 36 /
      data pcmn(152)
     &/'(D=51o) Extension for dueteron library'/
      data lpcn(152)/ 38 /
      data pcmn(151)
     &/'(D=51u) Extension for photo-nuclear library'/
      data lpcn(151)/ 43 /
      data pcmn(150)
     &/'(D=51c) Extension for high-energy neutron library'/
      data lpcn(150)/ 49 /
      data pcmn(149)
     &/'(D=51h) Extension for proton library'/
      data lpcn(149)/ 36 /

      data pcmn(148)
     &/'(D=1) Cross section for n-p scattering 0: JAM, 1: JENDL/HE'/
      data lpcn(148)/ 58 /

      data pcmn(147)
     &/'(D=0) Error line output 0: off, 1: on'/
      data lpcn(147)/ 37 /

      data pcmn(146)
     &/'(D=0) Surface Coalescense Model for JQMD. 0: off, 1: on'/
      data lpcn(146)/ 55 /

      data pcmn(145)
     &/'(D=0) Same material check option for gshow 0:Check, 1:No check'/
      data lpcn(145)/ 62 /

      data pcmn(144)
     &/'(D=0) Shared tally for multi threads. 0:Non-shared, 1:shared'/
      data lpcn(144)/ 60 /

      data pcmn(143)
     &/'(D=0) Fission model. 0: Default(Ver.1), 1:Legacy, 2:Iwamoto'/
      data lpcn(143)/ 59 /

      data pcmn(142)
     &/'(D=3) mode. 0:no echo, 1:read, 2:write, 3:nomal'/
      data lpcn(142)/ 47 /

      data pcmn(141)
     &/'(D=0) SCINFUL mode. 0:off, 1:on'/
      data lpcn(141)/ 31 /

c S.H. xorshift (2020.2.6)
      data pcmn(140)
     &/'(D=0) option for time dependent initial random number'/
      data lpcn(140)/ 53 /

      data pcmn(139)
     &/'(D=1) pseudo-random number generator, 0:LCG, 1: xorshift'/
      data lpcn(139)/ 56 /

      data pcmn(138)
     &/'(D=2) 0: no for lat, 1: with, 2: no between the same cell'/
      data lpcn(138)/ 57 /

      data pcmn(137)
     &/'(D=0) 0: no, 1: cosider photon-induced muon pair production'/
      data lpcn(137)/ 59 /

      data pcmn(136)
     &/'(D=1) minimum batch number for stdcut'/
      data lpcn(136)/ 37 /

      data pcmn(135)
     &/'(D=1) =0 no input echo for tally output file'/
      data lpcn(135)/ 44 /

      data pcmn(134)
     &/'(D=0) stdcut for each tally  0: off, 1: on'/
      data lpcn(134)/ 42 /

      data pcmn(133)
     &/'(D=0) =1 Automatic mode for tetrahedrons'/
      data lpcn(133)/ 40 /

      data pcmn(132)
     &/'(D=0) ATIMA charge state, 0:effective, 1:fixed'/
      data lpcn(132)/ 46 /

      data pcmn(131)
     &/'(D=0) Neutrino reactions, 0: off, 1: on'/
      data lpcn(131)/ 39 /

      data pcmn(130)
     &/'(D=0) 0: Auto, 1:do not use, 2:use kerma for photon'/
      data lpcn(130)/ 51 /

      data pcmn(129)
     &/'(D=0) 0: Auto, 1:do not use, 2:use kerma for neutron'/
      data lpcn(129)/ 52 /

      data pcmn(128)
     &/'(D=1) 0: do not use, 1: use normal nuclear data'/
      data lpcn(128)/ 47 /

      data pcmn(127)
     &/'(D=1) 0: conventional mode, 1: restricted LET mode'/
      data lpcn(127)/ 50 /

c Ogawa 2017/04/28 Option to activate Advanced GEM
      data pcmn(126)
     &/'(D=0) 0: Default (Ver.1), 1: GEM Ver.1, 2:GEM Ver.2,
     & 100:GEM-Furihata'/
      data lpcn(126)/ 69 /

      data pcmn(125)
     &/'(D=0) WW Bias on(1) off(0)'/
      data lpcn(125)/ 26 /

c T.Sato 2016.12.28 commnet in [material] section
      data pcmn(124)
     &/'(D=0) c is a comment signal or not in [material], 0:NO, 1:YES'/
      data lpcn(124)/ 61 /

cMIURA 2016.09.30 tally output unit, the energy per nucleon. [MeV/n]
      data pcmn(123)
     &/'(D=0) 0:[MeV] or 1:[MeV/n] is unit of tally output'/
      data lpcn(123)/ 50 /

      data pcmn(122)
     &/'(D=200) Number of tetra element allowed in sub-block'/
      data lpcn(122)/ 52 /

      data pcmn(121)
     &/'(D=100) Number of tetra outer surface allowed in sub-section'/
      data lpcn(121)/ 60 /

      data pcmn(120)
     &/'(D=0) tetra data is read(=1)/write(=2) on binary'/
      data lpcn(120)/ 48 /


      data pcmn(119)
     &/'(D=0) 0: =icrhi, 1: MWO formula for Deuteron'/
      data lpcn(119)/ 44 /

      data pcmn(118)
     &/'(D=1) 0:  No, 1: Expand Nat Nucleus, 2: Expand and Echo'/
      data lpcn(118)/ 55 /

      data pcmn(117)
     &/'(D=1) 0:geometrical,1:PHITS original, pion cross section model'/
      data lpcn(117)/ 62 /

c T.Sato 2016/2/17 Shielding distribution calculation mode
      data pcmn(116)
     &/'(D=0) =1, shielding distribution calculation mode'/
      data lpcn(116)/ 49 /

cFURUTA20160126 ATIMA database by Wada
      data pcmn(115)
     &/'(D=500) max database size of ATIMA'/
      data lpcn(115)/ 34 /

      data pcmn(114)
     &/'(D=0) adjoint mode for photon (=1) or CP (=2)'/
      data lpcn(114)/ 45 /

      data pcmn(113)
     &/'(D=1) real fission, =0; Fission turnoff, treated as capture'/
      data lpcn(113)/ 64 /

      data pcmn(112)
     &/'(D=1) multiple scattering option, 0:EGS original, 1:PHITS'/
      data lpcn(112)/ 57 /

      data pcmn(111)
     &/'(D=1) 0: no, 1: deltm and deltc divided by density'/
      data lpcn(111)/ 50 /

      data pcmn(110)
     &/'(D=1) 0: no, 1: muon-induced pair production'/
      data lpcn(110)/ 44 /

      data pcmn(109)
     &/'(D=1) 0: no, 1: muon-induced bremsstrahlung'/
      data lpcn(109)/ 43 /

      data pcmn(108)
     &/'[EGS](D=0) 0:no, 1:consider density effect based on ICRU90'/
      data lpcn(108)/ 58 /

      data pcmn(107)
     &/'[EGS](D=0):No output, =1:Output information in PEGS5'/
      data lpcn(107)/ 52 /

      data pcmn(106)
     &/'[EGS](D=0) -1:PEGSonly, 0:Full, 1:PEGSrun+PHITS, 2:skipPEGS'/
      data lpcn(106)/ 59 /

      data pcmn(105)
     &/'(D=0) =1 Volume calculation for tetrahedrons'/
      data lpcn(105)/ 44 /

cABE 2018/10/30, add description for =2
      data pcmn(104)
     &/'(D=1) 0: no, 1: muon capture, 2: muon capture reading file(28)'/
      data lpcn(104)/ 62 /

      data pcmn(103)
     &/'(D=0) number of output warnings for geometry recovering'/
      data lpcn(103)/ 55 /

      data pcmn(102)
     &/'[EGS](D=0) no EGS file, 1=>pegs5.inp,pegs5.dat, 2=> all'/
      data lpcn(102)/ 60 /

      data pcmn(101)
     &/'[EGS](D=1) Brems. XS correction, 0=>Motz 1=>ICRU37 2=>No'/
      data lpcn(101)/ 56 /

      data pcmn(100)
     &/'[EGS](D=1) 0=>Free Compton XS, 1=>Bound total Compton XS'/
      data lpcn(100)/ 56 /

      data pcmn(99)
     &/'[EGS](D=1) Photoelectron angle, 0=>Default 1=>Sampling'/
      data lpcn(99)/ 54 /

      data pcmn(98)
     &/'[EGS](D=1) Order of sampling of polar angles of pair electrons'/
      data lpcn(98)/ 67 /

      data pcmn(97)
     &/'[EGS](D=1) Brems. polar angle, 0=>Default 1=>Sampling'/
      data lpcn(97)/ 53 /

      data pcmn(96)
     &/'[EGS](D=0) Number of EII x-rays for splitting when IEISPL=1'/
      data lpcn(96)/ 64 /

      data pcmn(95)
     &/'[EGS](D=0) X-rays by EII, 0=>No splitting 1=>Splitting'/
      data lpcn(95)/ 54 /

      data pcmn(94)
     &/'[EGS](D=1) Luxury level of random number generator'/
      data lpcn(94)/ 55 /

      data pcmn(93)
     &/'[EGS](D=-1) initial random number for EGS, <0 use PHITS'/
      data lpcn(93)/ 55 /

      data pcmn(92)
     &/'[EGS](D=1) Electron Impact Ionization, 0=>Ignore 1=>Consider'/
      data lpcn(92)/ 60 /

      data pcmn(91)
     &/'[EGS](D=1) Doppler broadening, 0=>Ignore 1=>Consider'/
      data lpcn(91)/ 52/

      data pcmn(90)
     &/'[EGS](D=1) Incoherent scattering function, 1=>Consider'/
      data lpcn(90)/ 54 /

      data pcmn(89)
     &/'[EGS](D=0) Linearly polarized photon scattering, 1=>Consider'/
      data lpcn(89)/ 60 /

      data pcmn(88)
     &/'[EGS](D=1) Rayleigh scattering, 0=>Ignore 1=>Consider'/
      data lpcn(88)/ 53 /

      data pcmn(87)
     &/'[EGS](D=1) K and L-Auger electrons, 1=>Consider'/
      data lpcn(87)/ 47 /

      data pcmn(86)
     &/'[EGS](D=1) K and L-edge fluorescent photons, 1=>Consider'/
      data lpcn(86)/ 56 /

      data pcmn(85)/'(D=-1) =-1:original, =0:No, =1:EGS, =2:EGS(HE)'/  ! T.Sato 2022/12/16
      data lpcn(85)/ 46 /

*                        5    1    5    2    5    3    5    4    5    5
      data pcmn(84)/'(D=1) 0: no, 1: muon interaction'/
      data lpcn(84)/ 32 /

      data pcmn(83)/'(D=0) =1 no WW for void-void boundary'/
      data lpcn(83)/ 37 /

      data pcmn(82)/'(D=7) Print out info. -1,0,1,2,3,4,5,6,7,8'/
      data lpcn(82)/ 42 /

      data pcmn(81)/'(D=0) 0: no, 1: Activate Fermi breakup Model'/
      data lpcn(81)/ 44 /

      data pcmn(80)/'(D=1) Restart, 0:Check consistency, 1:No check'/
      data lpcn(80)/ 46 /

      data pcmn(79)/'(D=0) 0: JQMD legacy version, 1: JQMD-2.0'/
      data lpcn(79)/ 41 /

      data pcmn(78)/'(D=0) random mode'/
      data lpcn(78)/ 17 /

      data pcmn(77)
     &  /'(D=0) additional option of discrete levels using DWBA'/
      data lpcn(77)/ 53 /

      data pcmn(76)/'(D=0) 0:no, 1:INC-ELF for p,n,alpha, 2:p,n only'/
      data lpcn(76)/ 47 /

      data pcmn(75)/'(D=0) 0: Pearlstein-Niita, 1: KUROTAMA, 2: Sato'/
      data lpcn(75)/ 47 /

      data pcmn(74)
     &  /'(D=0) 0: no, 1: Activate statistical multi-fragmentation'/
      data lpcn(74)/ 56 /

      data pcmn(73)
     &  /'(D=0) 0:Auto, 1,2:Batch or History variance, <0:Restart mode'/
      data lpcn(73)/ 60 /

      data pcmn(72)/'(D=0) version control option of INCL model'/
      data lpcn(72)/ 42 /

      data pcmn(71)
     &  /'(D=1) 0:no, 1:INCL for p,n,pi,d,t,3He,alpha, 2:p,n,pi only'/
      data lpcn(71)/ 58 /

      data pcmn(70)/'(D=0) electro magnetic field'/
      data lpcn(70)/ 29 /

      data pcmn(69)/'(D=0) 0: no, 1: consider photo-nuclear reaction'/
      data lpcn(69)/ 47 /

      data pcmn(68)/'(D=0) voxel data is read(=1)/write(=2) on binary'/
      data lpcn(68)/ 48 /

      data pcmn(67)/'(D=0) S(a,b) angular dis. option'/
      data lpcn(67)/ 32 /

      data pcmn(66)/'(D=0) Induced Fission option'/
      data lpcn(66)/ 28 /

      data pcmn(65)/'(D=0) 9C decay option'/
      data lpcn(65)/ 21 /

      data pcmn(64)/'(D=150) max time step for JQMD'/
      data lpcn(64)/ 30 /

      data pcmn(63)/'(D=1) usrelst subroutine option, 1 or 2'/
      data lpcn(63)/ 39 /

      data pcmn(62)/'(D=1) usrmgt subroutine option, 1 or 2'/
      data lpcn(62)/ 38 /

      data pcmn(61)/'(D=3) dump file name option'/
      data lpcn(61)/ 27 /

      data pcmn(60)/"(D=0) 0: No, 1: dump all information on file(15)"/
      data lpcn(60)/ 48 /

      data pcmn(59)/"(D=0) 0: Normal, 1: Event generator mode Ver.1,
     & 2: Ver.2"/
      data lpcn(59)/ 56 /

      data pcmn(58)/"(D=3) 0: SPAR, 1: ATIMA, 2: SPAR all for dE/dx,
     & 3: ATIMA all for dE/dx"/
      data lpcn(58)/ 70 /

      data pcmn(57)/"(D=1) 0: no, 1: moving frame adjust"/
      data lpcn(57)/ 35 /

      data pcmn(56)/"(D=1) 0: normal, 1: new judge of ielst"/
      data lpcn(56)/ 38 /

      data pcmn(55)/"(D=0) # of source written in output"/
      data lpcn(55)/ 35 /

      data pcmn(54)/"(D=2) 0: Shen, 1: NASA, 2: KUROTAMA, 3: Hashimoto"/
      data lpcn(54)/ 49 /

      data pcmn(53)/"(D=-1) Neutron data library opt. for light ions."/
      data lpcn(53)/ 48 /

      data pcmn(52)/'(D=1) 0: no, 1: tally by a step in mag and strg'/
      data lpcn(52)/ 47 /

      data pcmn(51)/'(D=0) 0: no, 1: with e straggling, 10: ATIMA'/
      data lpcn(51)/ 44 /

      data pcmn(50)/'(D=0) 0: no photonuclear, 1: analog, -1: biased'/
      data lpcn(50)/ 47 /

      data pcmn(49)/'(D=0) 0: slash for unix, 1: back slash for PCDOS'/
      data lpcn(49)/ 48 /

      data pcmn(48)/'(D=0) action point in Weight Window'/
      data lpcn(48)/ 35 /

      data pcmn(47)/'(D=5) max split number in Weight Window'/
      data lpcn(47)/ 39 /

      data pcmn(46)/'(D=0) electron above dmax(12), 1: e=dmax(12)'/
      data lpcn(46)/ 44 /

      data pcmn(45)/'(D=0) 0: input, 1: all parameter out'/
      data lpcn(45)/ 36 /

      data pcmn(44)/'(D=0) 0: no cpu time for each process'/
      data lpcn(44)/ 37 /

      data pcmn(43)/'(D=0) 1: Output library information in detail'/
      data lpcn(43)/ 45 /

      data pcmn(42)/'(D=0) >0: on each substep 0: nominal brems.'/
      data lpcn(42)/ 43 /

      data pcmn(41)/'(D=0) 0: sampled straggl. 1: expected value.'/
      data lpcn(41)/ 44 /

      data pcmn(40)/'(D=0) 0: full brems. 1: simple brems. ang. dist.'/
      data lpcn(40)/ 48 /

      data pcmn(39)/'(D=0) 0/1 electrons will/not produce photons'/
      data lpcn(39)/ 44 /

      data pcmn(38)/'(D=0) 0/1 photon coherent scatt. will/not occure'/
      data lpcn(38)/ 48 /

      data pcmn(37)/'(D=1) 0/1 photons will/not produce electrons'/
      data lpcn(37)/ 44 /

      data pcmn(36)/'(D=0) 0/1=on/off unresolved res. range prob. tab.'/
      data lpcn(36)/ 49 /

      data pcmn(35)/'(D=0) Den.echo, 0:input, 1:number, 2:weight'/
      data lpcn(35)/ 43 /

      data pcmn(34)/'(D=10000) maximum bank memory length'/
      data lpcn(34)/ 36 /

      data pcmn(33)/'(D=0) 0:no 1:save memory for itall=1,2 or para'/
      data lpcn(33)/ 46 /

      data pcmn(32)/'(D=0) 0:no 1:print GG message and warning'/
      data lpcn(32)/ 41 /

      data pcmn(31)/'(D=0) pcut file name option'/
      data lpcn(31)/ 27 /

      data pcmn(30)/'(D=0) cutoff proton on file, 0:No, 2:with Time'/
      data lpcn(30)/ 46 /

      data pcmn(29)/'(D=0) gcut file name option'/
      data lpcn(29)/ 27 /

      data pcmn(28)/'(D=0) ncut file name option'/
      data lpcn(28)/ 27 /

      data pcmn(27)/'(D=1) 0:no 1:add matt with different den for GG'/
      data lpcn(27)/ 47 /

      data pcmn(26)/'(D=0) 0:no tally at batch, 1:same, 2:different'/
      data lpcn(26)/ 46 /

      data pcmn(25)/'(D=0) 0:linear, 1:log for energy mesh of inucr'/
      data lpcn(25)/ 46 /

      data pcmn(24)/'(D=720) mesh of elastic angular distribution'/
      data lpcn(24)/ 44 /

      data pcmn(23)/'(D=0) Check overlap region errors in CG and GG'/
      data lpcn(23)/ 46 /

      data pcmn(22)/'(D=1) Recovery number of CG and GG errors'/
      data lpcn(22)/ 41 /

      data pcmn(21)/'(D=10) Lost part: allowance of CG and GG errors'/
      data lpcn(21)/ 47 /

      data pcmn(20)/'(D=1) import. echo 0: [importance], 1: [region]'/
      data lpcn(20)/ 47 /

      data pcmn(19)/'(D=0) volume echo 0: [volume], 1: [region]'/
      data lpcn(19)/ 42 /

      data pcmn(18)/'(D=0) cutoff photon on file, 0:No, 2:with Time'/
      data lpcn(18)/ 46 /

      data pcmn(17)/'(D=0) cutoff neutron on file, 0:No, 2:with Time'/
      data lpcn(17)/ 47 /

      data pcmn(16)/'(D=0) MAT format, 0: 208Pb, 1: Pb-208, 2: MCNP'/
      data lpcn(16)/ 46 /

      data pcmn(15)/'(D=1) options of nuclear reaction for icntl = 1'/
      data lpcn(15)/ 47 /

      data pcmn(14)/'(D=0) magnetic field'/
      data lpcn(14)/ 20 /

      data pcmn(13)/'(D=1) 0:free, 1:Cugnon, 2:Cugnon new, for NN'/
      data lpcn(13)/ 44 /

      data pcmn(12)/'(D=2) 0:No, 1:Old, 2:EBITEM, 3:EBITEM+Isomer'/
      data lpcn(12)/ 44 /

      data pcmn(11)/'(D=3) 1:A/8, 2:Baba, 3:Ignatyuk for nevap=1'/
      data lpcn(11)/ 43 /

      data pcmn(10)/'(D=0) pre-equilibrium'/
      data lpcn(10)/ 21 /

      data pcmn( 9)/'(D=0) Isobar model, with <eisobar>'/
      data lpcn( 9)/ 34 /

      data pcmn( 8)/'(D=2) 0:non, 1:neutron, 2:neut. and proton'/
      data lpcn( 8)/ 42 /

      data pcmn( 7)/'(D=0) 0:no, 1:original, 2:Lynch(Moliere), 10:ATIMA'
     &/
      data lpcn( 7)/ 50 /

      data pcmn( 6)/'(D=3) 0:non, 1:DRES, 2:SDM, 3:GEM, for evap.'/
      data lpcn( 6)/ 44 /

      data pcmn( 5)/'(D=0) 0:reaction, 1:stopped, for negative mesons'/
      data lpcn( 5)/ 48 /

      data pcmn( 4)/'(D=10) number of batches'/
      data lpcn( 4)/ 24 /

      data pcmn( 3)/'(D=10) number of particles per one batch'/
      data lpcn( 3)/ 40 /

      data pcmn( 2)/'(D=0) >0:skip events, <0:skip random number'/
      data lpcn( 2)/ 43 /

      data pcmn( 1)/'(D=0) 3:ECH 5:NOR 6:SRC 7,8:GSH 11:DSH 12:DUMP'/
      data lpcn( 1)/ 46 /

*                        5    1    5    2    5    3    5    4    5    5
*-----------------------------------------------------------------------
*     Set comments for parz
*-----------------------------------------------------------------------
*           5    1    5    2    5    3    5    4    5    5    5    6
! T.Sato 2024/10/11
      data pcmr(206)
     &/'(D=0.0) Rutherford bias angle upper boundary'/
      data lpcr(206)/ 44 /

      data pcmr(205)
     &/'(D=0.0) Rutherford bias angle lower boundary'/
      data lpcr(205)/ 44 /

! T.Sato 2024/03/24
      data pcmr(204)
     &/'(D=1.0e-4) small gap to avoid integer value in xyz mesh'/
      data lpcr(204)/ 55 /

cfrtati 2021/12/17
      data pcmr(203)
     &/'(D=emin(14)) data max. energy for photo-nuclear lib. (MeV)'/
      data lpcr(203)/ 58 /

      data pcmr(202)
     &/'(D=0.2) Track banked particles if maxbank*rtrkflp is occupied'/
      data lpcr(202)/ 61 /

      data pcmr(201)
     &/'(D=5.0) Relative memory space used for nuclear data'/
      data lpcr(201)/ 51 /

      data pcmr(200)
     &/'(D=-1.124) parameter Rc for JQMD surface coalescence model'/
      data lpcr(200)/ 58 /

      data pcmr(199)
     &/'(D=1.75) parameter D for JQMD surface coalescence model'/
      data lpcr(199)/ 55 /

      data pcmr(198)
     &/'(D=0.25) parameter h0 for JQMD surface coalescence model'/
      data lpcr(198)/ 56 /

      data pcmr(197)
     &/'(D=1.d1) energy for pseudo reaction for charged particles'/
      data lpcr(197)/ 57 /

      data pcmr(196)
     &/'(D=1.d-3) maximum energy of proton and ion track structure'/
      data lpcr(196)/ 58 /

c S.H. xorshift (2020.2.6)
      data pcmr(195)
     &/'(D=[omitted]) bit data of initial random seed'/
      data lpcr(195)/ 45 /

      data pcmr(194)
     &/'(D=1.0) multiplying factor for photon-induced mu-pair prod. CS'/
      data lpcr(194)/ 62 /

      data pcmr(193)
     &/'(D=20) maximum energy of e-mode for neutron'/
      data lpcr(193)/ 43 /

      data pcmr(192)
     &/'(D=1.0D-2) Switch energy of electron for track structure mode'/
      data lpcr(192)/ 61 /

      data pcmr(191)
     &/'(D=1.0D-6) Cutoff energy of electron for track structure mode'/
      data lpcr(191)/ 61 /

      data pcmr(190)
     &/'(D=-1.0) Cutoff CPU time (sec)'/
      data lpcr(190)/ 30 /

      data pcmr(189)
     &/'(D=3.0) Maxmum energy for adjoint photon transport'/
      data lpcr(189)/ 50 /

      data pcmr(188)
     &/'(D=0.05) threshold of Bjorken variable for muon interaction CS'/
      data lpcr(188)/ 62 /

cFURUTA20160126 ATIMA database by Wada
      data pcmr(187)
     &/'(D=0.0) energy cutoff of ATIMA database (MeV/n)'/
      data lpcr(187)/ 47 /

      data pcmr(186)
     &/'(D=1.0) multiplying factor for photonuclear interaction CS'/
      data lpcr(186)/ 58 /

      data pcmr(185)
     &/'(D=0.1) Minimum length of a region (used in EGS5)'/
      data lpcr(185)/ 49 /

      data pcmr(184)/
     &'(D=0.038) e parameter in Lynch formula for nspred = 2'/
      data lpcr(184)/ 53/

      data pcmr(183)/
     &'(D=13.6) S2 parameter in Lynch formula for nspred = 2'/
      data lpcr(183)/ 53/

      data pcmr(182)/'(D=0.03) maxium density (g/cm3) of gas in EGS'/
      data lpcr(182)/ 45/

      data pcmr(181)/'(D=1.0d6) max energy of muon interaction (MeV)'/
      data lpcr(181)/ 46 /

      data pcmr(180)/'(D=200.0) min energy of muon interaction (MeV)'/
      data lpcr(180)/ 46 /

      data pcmr(179)/'(D=3500.0) max energy of INC-ELF (MeV)'/
      data lpcr(179)/ 38 /

      data pcmr(178)/'(D=1.0) min energy of INC-ELF (MeV)'/
      data lpcr(178)/ 35 /

      data pcmr(177)/'(D=3000.0) max energy of INCL (MeV/n)'/
      data lpcr(177)/ 37 /

      data pcmr(176)/'(D=1.0) min energy of INCL (MeV/n)'/
      data lpcr(176)/ 34 /

      data pcmr(175)
     &  /"(D=0.0) Density Symmetry Coefficient of KUROTAMA (MeV)"/
      data lpcr(175)/ 54 /

      data pcmr(174)/"(D=-1) Ion.Pot.for H2O, <0:75, >0 IP eV"/
      data lpcr(174)/ 39 /

      data pcmr(173)/'(D=1MeV) ge2; e-mode parameter'/
      data lpcr(173)/ 30 /

      data pcmr(172)/'(D=1MeV) ge1; e-mode parameter'/
      data lpcr(172)/ 30 /

      data pcmr(171)/'(D=1msec) max flight time for T-dep. Mag. Field'/
      data lpcr(171)/ 47 /

      data pcmr(170)
     &     /'(D=1.0e6) max energy/u for charged particle range'/
      data lpcr(170)/ 49 /

      data pcmr(169)/'(D=3000.0) switch energy from QMD to JAMQMD'/
      data lpcr(169)/ 43 /

      data pcmr(168)/'(D=0.0) extra impact parameter for QMD'/
      data lpcr(168)/ 38 /

      data pcmr(167)/'(D=0.0) z component of gravity'/
      data lpcr(167)/ 30 /

      data pcmr(166)/'(D=0.0) y component of gravity'/
      data lpcr(166)/ 30 /

      data pcmr(165)/'(D=0.0) x component of gravity'/
      data lpcr(165)/ 30 /

      data pcmr(164)/'(D=10.0) min energy / u for QMD reaction'/
      data lpcr(164)/ 40 /

      data pcmr(163)/'(D=2.012345) max flight mesh for nedisp'/
      data lpcr(163)/ 41 /

      data pcmr(162)
     &     /'(D=0.001) min energy/u for charged particle range'/
      data lpcr(162)/ 50 /

      data pcmr(161)
     &    /'(D=-100->0.6*wupn) survive factor in Weight Window'/
      data lpcr(161)/ 50 /

      data pcmr(160)/'(D=5) upper weight factor in Weight Window'/
      data lpcr(160)/ 42 /

      data pcmr(159)/'(D=0.) dbcn(18): energy indexing algorithm'/
      data lpcr(159)/ 42 /

      data pcmr(158)/'(D=0.) dbcn(17): angular treatment'/
      data lpcr(158)/ 34 /

      data pcmr(157)/'(D=100) upper limit of photons (MeV)'/
      data lpcr(157)/ 36 /

      data pcmr(156)/'(D=1) for photon-induced secondary electrons'/
      data lpcr(156)/ 49 /

      data pcmr(155)/'(D=1) # of knock-on electrons'/
      data lpcr(155)/ 29 /

      data pcmr(154)/'(D=1) for smapling of x-ray photons'/
      data lpcr(154)/ 35 /

      data pcmr(153)/'(D=1) for smapling of brems. photons'/
      data lpcr(153)/ 36 /

      data pcmr(152)/'(D=-1) # of delayed neutron from fission'/
      data lpcr(152)/ 40 /

      data pcmr(151)/'(D=0) analog vs implicit capture (MeV)'/
      data lpcr(151)/ 38 /

      data pcmr(150)/'(D=emin(20)) data max. energy of Non (MeV)'/
      data lpcr(150)/ 42 /

      data pcmr(149)/'(D=emin(19)) data max. energy of Nucleus (MeV/n)'/
      data lpcr(149)/ 48 /

      data pcmr(148)/'(D=emin(18)) data max. energy of Alpha (MeV/n)'/
      data lpcr(148)/ 46 /

      data pcmr(147)/'(D=emin(17)) data max. energy of 3He (MeV/n)'/
      data lpcr(147)/ 44 /

      data pcmr(146)/'(D=emin(16)) data max. energy of triton (MeV/n)'/
      data lpcr(146)/ 47 /

      data pcmr(145)/'(D=emin(15)) data max. energy of deuteron (MeV/n)'
     &/
      data lpcr(145)/ 49 /

      data pcmr(144)/'(D=1.0e+3) data max. energy of photon (MeV)'/
      data lpcr(144)/ 43 /

      data pcmr(143)/'(D=emin(13)) data max. energy of positron (MeV)'/
      data lpcr(143)/ 47 /

      data pcmr(142)/'(D=emin(12)) data max. energy of electron (MeV)'/
      data lpcr(142)/ 47 /

      data pcmr(141)/'(D=emin(11)) data max. energy of other (MeV)'/
      data lpcr(141)/ 43 /

      data pcmr(140)/'(D=emin(10)) data max. energy of kaon- (MeV)'/
      data lpcr(140)/ 44 /

*                         5    1    5    2    5    3    5    4    5    5
      data pcmr(139)/'(D=emin(9)) data max. energy of kaon0 (MeV)'/
      data lpcr(139)/ 43 /

      data pcmr(138)/'(D=emin(8)) data max. energy of kaon+ (MeV)'/
      data lpcr(138)/ 43 /

      data pcmr(137)/'(D=emin(7)) data max. energy of muon- (MeV)'/
      data lpcr(137)/ 43 /

      data pcmr(136)/'(D=emin(6)) data max. energy of muon+ (MeV)'/
      data lpcr(136)/ 43 /

      data pcmr(135)/'(D=emin(5)) data max. energy of pion- (MeV)'/
      data lpcr(135)/ 43 /

      data pcmr(134)/'(D=emin(4)) data max. energy of pion0 (MeV)'/
      data lpcr(134)/ 43 /

      data pcmr(133)/'(D=emin(3)) data max. energy of pion+ (MeV)'/
      data lpcr(133)/ 43 /

      data pcmr(132)/'(D=20.0) data max. energy of neutron (MeV)'/
      data lpcr(132)/ 42 /

      data pcmr(131)/'(D=emin(1)) data max. energy of proton (MeV)'/
      data lpcr(131)/ 44 /


*                         5    1    5    2    5    3    5    4    5    5
      data pcmr(130)/'(D=emin(20)) reaction cut-off of Non (MeV)'/
      data lpcr(130)/ 42 /

      data pcmr(129)/'(D=emin(19)) reaction cut-off of Nucleus (MeV/n)'/
      data lpcr(129)/ 48 /

      data pcmr(128)/'(D=emin(18)) reaction cut-off of Alpha (MeV/n)'/
      data lpcr(128)/ 46 /

      data pcmr(127)/'(D=emin(17)) reaction cut-off of 3He (MeV/n)'/
      data lpcr(127)/ 44 /

      data pcmr(126)/'(D=emin(16)) reaction cut-off of triton (MeV/n)'/
      data lpcr(126)/ 47 /

      data pcmr(125)/'(D=emin(15)) reaction cut-off of deuteron (MeV/n)'
     &/
      data lpcr(125)/ 49 /

      data pcmr(124)/'(D=emin(14)) reaction cut-off of photon (MeV)'/
      data lpcr(124)/ 45 /

      data pcmr(123)/'(D=emin(13)) reaction cut-off of positron (MeV)'/
      data lpcr(123)/ 47 /

      data pcmr(122)/'(D=emin(12)) reaction cut-off of electron (MeV)'/
      data lpcr(122)/ 47 /

      data pcmr(121)/'(D=emin(11)) reaction cut-off of other (MeV)'/
      data lpcr(121)/ 43 /

      data pcmr(120)/'(D=emin(10)) reaction cut-off of kaon- (MeV)'/
      data lpcr(120)/ 44 /

*                         5    1    5    2    5    3    5    4    5    5
      data pcmr(119)/'(D=emin(9)) reaction cut-off of kaon0 (MeV)'/
      data lpcr(119)/ 43 /

      data pcmr(118)/'(D=emin(8)) reaction cut-off of kaon+ (MeV)'/
      data lpcr(118)/ 43 /

      data pcmr(117)/'(D=emin(7)) reaction cut-off of muon- (MeV)'/
      data lpcr(117)/ 43 /

      data pcmr(116)/'(D=emin(6)) reaction cut-off of muon+ (MeV)'/
      data lpcr(116)/ 43 /

      data pcmr(115)/'(D=emin(5)) reaction cut-off of pion- (MeV)'/
      data lpcr(115)/ 43 /

      data pcmr(114)/'(D=emin(4)) reaction cut-off of pion0 (MeV)'/
      data lpcr(114)/ 43 /

      data pcmr(113)/'(D=emin(3)) reaction cut-off of pion+ (MeV)'/
      data lpcr(113)/ 43 /

      data pcmr(112)/'(D=emin(2)) reaction cut-off of neutron (MeV)'/
      data lpcr(112)/ 45 /

      data pcmr(111)/'(D=emin(1)) reaction cut-off of proton (MeV)'/
      data lpcr(111)/ 44 /

*                         5    1    5    2    5    3    5    4    5    5
      data pcmr(110)/'(D=20.0) energy of QMD for nucleon (MeV)'/
      data lpcr(110)/ 40 /

      data pcmr(109)
     &    /'(D=1.012345) max. flight mesh for magnet field (cm)'/
      data lpcr(109)/ 51 /

*                         5    1    5    2    5    3    5    4    5    5
      data pcmr(108)/'(D=-100->wc1/2) weight cutoff of Non'/
      data lpcr(108)/ 36 /

      data pcmr(107)/'(D=-100->wc1/2) weight cutoff of Nucleus'/
      data lpcr(107)/ 40 /

      data pcmr(106)/'(D=-100->wc1/2) weight cutoff of Alpha'/
      data lpcr(106)/ 38 /

      data pcmr(105)/'(D=-100->wc1/2) weight cutoff of 3He'/
      data lpcr(105)/ 36 /

      data pcmr(104)/'(D=-100->wc1/2) weight cutoff of triton'/
      data lpcr(104)/ 39 /

      data pcmr(103)/'(D=-100->wc1/2) weight cutoff of deuteron'/
      data lpcr(103)/ 41 /

      data pcmr(102)/'(D=-100->wc1/2) weight cutoff of photon'/
      data lpcr(102)/ 39 /

      data pcmr(101)/'(D=-100->wc1/2) weight cutoff of positron'/
      data lpcr(101)/ 41 /

      data pcmr(100)/'(D=-100->wc1/2) weight cutoff of electron'/
      data lpcr(100)/ 41 /
*                         5    1    5    2    5    3    5    4    5    5

*                        5    1    5    2    5    3    5    4    5    5
      data pcmr(99)/'(D=-100->wc1/2) weight cutoff of other particles'/
      data lpcr(99)/ 48 /

      data pcmr(98)/'(D=-100->wc1/2) weight cutoff of kaon-'/
      data lpcr(98)/ 38 /

      data pcmr(97)/'(D=-100->wc1/2) weight cutoff of kaon0'/
      data lpcr(97)/ 38 /

      data pcmr(96)/'(D=-100->wc1/2) weight cutoff of kaon+'/
      data lpcr(96)/ 38 /

      data pcmr(95)/'(D=-100->wc1/2) weight cutoff of muon-'/
      data lpcr(95)/ 38 /

      data pcmr(94)/'(D=-100->wc1/2) weight cutoff of muon+'/
      data lpcr(94)/ 38 /

      data pcmr(93)/'(D=-100->wc1/2) weight cutoff of pion-'/
      data lpcr(93)/ 38 /

      data pcmr(92)/'(D=-100->wc1/2) weight cutoff of pion0'/
      data lpcr(92)/ 38 /

      data pcmr(91)/'(D=-100->wc1/2) weight cutoff of pion+'/
      data lpcr(91)/ 38 /

      data pcmr(90)/'(D=-100->wc1/2) weight cutoff of neutron'/
      data lpcr(90)/ 40 /

      data pcmr(89)/'(D=-100->wc1/2) weight cutoff of proton'/
      data lpcr(89)/ 39 /

*                        5    1    5    2    5    3    5    4    5    5
      data pcmr(88)/'(D=-0.5) minimum weight of Non'/
      data lpcr(88)/ 30 /

      data pcmr(87)/'(D=-0.5) minimum weight of Nucleus'/
      data lpcr(87)/ 34 /

      data pcmr(86)/'(D=-0.5) minimum weight of Alpha'/
      data lpcr(86)/ 32 /

      data pcmr(85)/'(D=-0.5) minimum weight of 3He'/
      data lpcr(85)/ 30 /

      data pcmr(84)/'(D=-0.5) minimum weight of triton'/
      data lpcr(84)/ 33 /

      data pcmr(83)/'(D=-0.5) minimum weight of deuteron'/
      data lpcr(83)/ 35 /

      data pcmr(82)/'(D=-0.5) minimum weight of photon'/
      data lpcr(82)/ 33 /

      data pcmr(81)/'(D=0.0) minimum weight of positron'/
      data lpcr(81)/ 34 /

      data pcmr(80)/'(D=0.0) minimum weight of electron'/
      data lpcr(80)/ 34 /

*                        5    1    5    2    5    3    5    4    5    5
      data pcmr(79)/'(D=-0.5) minimum weight of other particles'/
      data lpcr(79)/ 42 /

      data pcmr(78)/'(D=-0.5) minimum weight of kaon-'/
      data lpcr(78)/ 32 /

      data pcmr(77)/'(D=-0.5) minimum weight of kaon0'/
      data lpcr(77)/ 32 /

      data pcmr(76)/'(D=-0.5) minimum weight of kaon+'/
      data lpcr(76)/ 32 /

      data pcmr(75)/'(D=-0.5) minimum weight of muon-'/
      data lpcr(75)/ 32 /

      data pcmr(74)/'(D=-0.5) minimum weight of muon+'/
      data lpcr(74)/ 32 /

      data pcmr(73)/'(D=-0.5) minimum weight of pion-'/
      data lpcr(73)/ 32 /

      data pcmr(72)/'(D=-0.5) minimum weight of pion0'/
      data lpcr(72)/ 32 /

      data pcmr(71)/'(D=-0.5) minimum weight of pion+'/
      data lpcr(71)/ 32 /

      data pcmr(70)/'(D=-0.5) minimum weight of neutron'/
      data lpcr(70)/ 34 /

      data pcmr(69)/'(D=-0.5) minimum weight of proton'/
      data lpcr(69)/ 33 /

*                        5    1    5    2    5    3    5    4    5    5
      data pcmr(68)/'(D=1.0) minimum source weight of Non'/
      data lpcr(68)/ 36 /

      data pcmr(67)/'(D=1.0) minimum source weight of Nucleus'/
      data lpcr(67)/ 40 /

      data pcmr(66)/'(D=1.0) minimum source weight of Alpha'/
      data lpcr(66)/ 38 /

      data pcmr(65)/'(D=1.0) minimum source weight of 3He'/
      data lpcr(65)/ 36 /

      data pcmr(64)/'(D=1.0) minimum source weight of triton'/
      data lpcr(64)/ 39 /

      data pcmr(63)/'(D=1.0) minimum source weight of deuteron'/
      data lpcr(63)/ 41 /

      data pcmr(62)/'(D=1.0) minimum source weight of photon'/
      data lpcr(62)/ 39 /

      data pcmr(61)/'(D=1.0) minimum source weight of positron'/
      data lpcr(61)/ 41 /

      data pcmr(60)/'(D=1.0) minimum source weight of electron'/
      data lpcr(60)/ 41 /

*                        5    1    5    2    5    3    5    4    5    5
      data pcmr(59)/'(D=1.0) minimum source weight of other particles'/
      data lpcr(59)/ 48 /

      data pcmr(58)/'(D=1.0) minimum source weight of kaon-'/
      data lpcr(58)/ 38 /

      data pcmr(57)/'(D=1.0) minimum source weight of kaon0'/
      data lpcr(57)/ 38 /

      data pcmr(56)/'(D=1.0) minimum source weight of kaon+'/
      data lpcr(56)/ 38 /

      data pcmr(55)/'(D=1.0) minimum source weight of muon-'/
      data lpcr(55)/ 38 /

      data pcmr(54)/'(D=1.0) minimum source weight of muon+'/
      data lpcr(54)/ 38 /

      data pcmr(53)/'(D=1.0) minimum source weight of pion-'/
      data lpcr(53)/ 38 /

      data pcmr(52)/'(D=1.0) minimum source weight of pion0'/
      data lpcr(52)/ 38 /

      data pcmr(51)/'(D=1.0) minimum source weight of pion+'/
      data lpcr(51)/ 38 /

      data pcmr(50)/'(D=1.0) minimum source weight of neutron'/
      data lpcr(50)/ 40 /

      data pcmr(49)/'(D=1.0) minimum source weight of proton'/
      data lpcr(49)/ 39 /

*                        5    1    5    2    5    3    5    4    5    5
      data pcmr(48)/'(D=emin()) cut-off time of Non (ns)'/
      data lpcr(48)/ 35 /

      data pcmr(47)/'(D=emin()) cut-off time of Nucleus (ns)'/
      data lpcr(47)/ 39 /

      data pcmr(46)/'(D=emin()) cut-off time of Alpha (ns)'/
      data lpcr(46)/ 37 /

      data pcmr(45)/'(D=emin()) cut-off time of 3He (ns)'/
      data lpcr(45)/ 35 /

      data pcmr(44)/'(D=emin()) cut-off time of triton (ns)'/
      data lpcr(44)/ 38 /

      data pcmr(43)/'(D=emin()) cut-off time of deuteron (ns)'/
      data lpcr(43)/ 40 /

      data pcmr(42)/'(D=emin()) cut-off time of photon (ns)'/
      data lpcr(42)/ 38 /

      data pcmr(41)/'(D=emin()) cut-off time of positron (ns)'/
      data lpcr(41)/ 40 /

      data pcmr(40)/'(D=emin()) cut-off time of electron (ns)'/
      data lpcr(40)/ 40 /

*                        5    1    5    2    5    3    5    4    5    5
      data pcmr(39)/'(D=1.0e+9) cut-off time of other particles (ns)'/
      data lpcr(39)/ 47 /

      data pcmr(38)/'(D=1.0e+9) cut-off time of kaon- (ns)'/
      data lpcr(38)/ 37 /

      data pcmr(37)/'(D=1.0e+9) cut-off time of kaon0 (ns)'/
      data lpcr(37)/ 37 /

      data pcmr(36)/'(D=1.0e+9) cut-off time of kaon+ (ns)'/
      data lpcr(36)/ 37 /

      data pcmr(35)/'(D=1.0e+9) cut-off time of muon- (ns)'/
      data lpcr(35)/ 37 /

      data pcmr(34)/'(D=1.0e+9) cut-off time of muon+ (ns)'/
      data lpcr(34)/ 37 /

      data pcmr(33)/'(D=1.0e+9) cut-off time of pion- (ns)'/
      data lpcr(33)/ 37 /

      data pcmr(32)/'(D=1.0e+9) cut-off time of pion0 (ns)'/
      data lpcr(32)/ 37 /

      data pcmr(31)/'(D=1.0e+9) cut-off time of pion+ (ns)'/
      data lpcr(31)/ 37 /

*                        5    1    5    2    5    3    5    4    5    5
      data pcmr(30)/'(D=1.0e+9) cut-off time of neutron (ns)'/
      data lpcr(30)/ 39 /

      data pcmr(29)/'(D=1.0e+9) cut-off time of proton (ns)'/
      data lpcr(29)/ 38 /

      data pcmr(28)/'(D=1.e-9) flight mesh after boundary (cm)'/
      data lpcr(28)/ 41 /

      data pcmr(27)/'(D=20.12345) max. flight mesh for transport (cm)'/
      data lpcr(27)/ 48 /

      data pcmr(26)/'(D=0.1) flight mesh for spread and magfield (cm)'/
      data lpcr(26)/ 48 /

      data pcmr(25)/'(D=0.0) switching energy of ISOBAR (MeV)'/
      data lpcr(25)/ 37 /

      data pcmr(24)/'(D=20.0) energy of JAM for pion (MeV)'/
      data lpcr(24)/ 37 /

      data pcmr(23)/'(D=20.0) energy of JAM for nucleon (MeV)'/
      data lpcr(23)/ 40 /

      data pcmr(22)/'(D=0.0) 0:iso+forward, 1:iso, 2:forward'/
      data lpcr(22)/ 39 /

      data pcmr(21)/'(D=0.0) 0:original, /=0:this value'/
      data lpcr(21)/ 34 /

      data pcmr(20)/'(D=1.e+9) cut-off energy of Non (MeV)'/
      data lpcr(20)/ 37 /

      data pcmr(19)/'(D=1.0e-3) cut-off energy of Nucleus (MeV/n)'/
      data lpcr(19)/ 44 /

      data pcmr(18)/'(D=1.0e-3) cut-off energy of Alpha (MeV/n)'/
      data lpcr(18)/ 42 /

      data pcmr(17)/'(D=1.0e-3) cut-off energy of 3He (MeV/n)'/
      data lpcr(17)/ 40 /

      data pcmr(16)/'(D=1.0e-3) cut-off energy of triton (MeV/n)'/
      data lpcr(16)/ 43 /

      data pcmr(15)/'(D=1.0e-3) cut-off energy of deuteron (MeV/n)'/
      data lpcr(15)/ 45 /

      data pcmr(14)/'(D=1.0e-3) cut-off energy of photon (MeV)'/
      data lpcr(14)/ 41 /

      data pcmr(13)/'(D=1.e+9) cut-off energy of positron (MeV)'/
      data lpcr(13)/ 42 /

      data pcmr(12)/'(D=1.e+9) cut-off energy of electron (MeV)'/
      data lpcr(12)/ 42 /

      data pcmr(11)/'(D=1.0) cut-off energy of other particles (MeV)'/
      data lpcr(11)/ 47 /

      data pcmr(10)/'(D=1.0e-3) cut-off energy of kaon- (MeV)'/
      data lpcr(10)/ 40 /

*                        5    1    5    2    5    3    5    4    5    5
      data pcmr( 9)/'(D=1.0e-3) cut-off energy of kaon0 (MeV)'/
      data lpcr( 9)/ 40 /

      data pcmr( 8)/'(D=1.0e-3) cut-off energy of kaon+ (MeV)'/
      data lpcr( 8)/ 40 /

      data pcmr( 7)/'(D=1.0e-3) cut-off energy of muon- (MeV)'/
      data lpcr( 7)/ 40 /

      data pcmr( 6)/'(D=1.0e-3) cut-off energy of muon+ (MeV)'/
      data lpcr( 6)/ 40 /

      data pcmr( 5)/'(D=1.0e-3) cut-off energy of pion- (MeV)'/
      data lpcr( 5)/ 40 /

      data pcmr( 4)/'(D=1.0e-3) cut-off energy of pion0 (MeV)'/
      data lpcr( 4)/ 40 /

      data pcmr( 3)/'(D=1.0e-3) cut-off energy of pion+ (MeV)'/
      data lpcr( 3)/ 40 /

      data pcmr( 2)/'(D=1.0e-11) cut-off energy of neutron (MeV)'/
      data lpcr( 2)/ 43 /

      data pcmr( 1)/'(D=1.0e-3) cut-off energy of proton (MeV)'/
      data lpcr( 1)/ 41 /

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine title(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [title] section of input files                            *
*       modified by K.Niita on 22/05/2000                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar
      common /htitl/ iclgt(100), ctitl(100)
      character ctitl*200

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

            ierr  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of title section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        store title messages
*-----------------------------------------------------------------------

               ititl = ititl + 1

               iclgt(ititl) = i3
               ctitl(ititl)(1:i3) = chin(1:i3)

*-----------------------------------------------------------------------

         goto 140

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine arrayg(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [arrayg] section of input files                           *
*       modified by K.Niita on 12/06/2000                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------
*        temporary file number
*-----------------------------------------------------------------------

            itar = 19

*-----------------------------------------------------------------------

            ierr  = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               return

            end if

*-----------------------------------------------------------------------
*        write array infomation on the temporary file 19
*-----------------------------------------------------------------------

            if( llarr .eq. 0 ) then

               open(itar,form='formatted',status='scratch')

            end if


            llarr = llarr + 1

            write(itar,'(i6,200a1)') i2, ( chin(i:i), i = 1, i2 )


         goto 140

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine cgview
*                                                                      *
*       output cgview.in file                                          *
*       modified by K.Niita on 26/07/2000                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_character

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar

      common /regdu/  iuni(kvlmax)
      common /regdm/  idmg(kvlmax)
      common /regda/  ichl(kvlmax), chsm(kvlmax), ichmx, iod
      character       chsm*10
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /regdb/  nrsq, irsq(10)

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      character cblan*200

      character chdf*80

      character chtit*60
      character chbd*3
      dimension bval(50)

      character chsc*10

      logical   exex

      character chin*200

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

         if( icntl .ne. 2 ) return
         if( iregn .le. 0 ) return
         if( ibody .le. 0 ) return

*-----------------------------------------------------------------------

         do i = 1, 200

            cblan(i:i) = ' '

         end do

*-----------------------------------------------------------------------
*        Files
*-----------------------------------------------------------------------

            iot = 21

            open(iot, file = chfn(2), status = 'unknown' )

*-----------------------------------------------------------------------

            rewind iod

            mcmx = ( mdas / 2 - 1 ) * 8 + 1
               call moddas_allocate_cha(MAX_NUM_CHRG, chrg)
            mci = 0

*-----------------------------------------------------------------------
*     Header
*-----------------------------------------------------------------------

         write(iot,'( ''$CROSS SECTION'')')
         write(iot,'( ''$END CROSS SECTION'')')
         write(iot,'(/''NMEMO(100)'')')
         write(iot,'(/''$GEOMETRY'')')

*-----------------------------------------------------------------------
*        body
*-----------------------------------------------------------------------

               rewind itby
               read(itby) chtit
               read(itby) ipva

            do k = 1, ibody

                  read(itby) chbd, ibnum, ibva, ( bval(i), i = 1, ibva )

                  call chcptl(chbd,1,3)

*-----------------------------------------------------------------------

               if( chbd .eq. 'REC' ) then

                  chbd     = 'TEC'
                  ibva     = 13
                  bval(13) = 1.0

               end if

               if( chbd .eq. 'GEL' ) then

                  r1 = sqrt( bval( 4)**2 + bval( 5)**2 + bval( 6)**2 )
                  r2 = sqrt( bval( 7)**2 + bval( 8)**2 + bval( 9)**2 )
                  r3 = sqrt( bval(10)**2 + bval(11)**2 + bval(12)**2 )

                  bval(10) = r1
                  bval(11) = r2
                  bval(12) = r3

               end if

               if( chbd .eq. 'TOR' ) then

                  chbd     = 'ELT'
                  ibva     = 10

                     d1 = bval(1)
                     d2 = bval(2)
                     d3 = bval(3)

                     d7 = bval(4)
                     d8 = bval(6)
                     d9 = bval(5)
                    d10 = 0.0

                  if( nint(bval(7)) .eq. 1 ) then

                     d4 = 1.0
                     d5 = 0.0
                     d6 = 0.0

                  else if( nint(bval(7)) .eq. 2 ) then

                     d4 = 0.0
                     d5 = 1.0
                     d6 = 0.0

                  else if( nint(bval(7)) .eq. 3 ) then

                     d4 = 0.0
                     d5 = 0.0
                     d6 = 1.0

                  end if

                     bval( 1) = d1
                     bval( 2) = d2
                     bval( 3) = d3
                     bval( 4) = d4
                     bval( 5) = d5
                     bval( 6) = d6
                     bval( 7) = d7
                     bval( 8) = d8
                     bval( 9) = d9
                     bval(10) = d10

               end if

*-----------------------------------------------------------------------

                  if( ibva .eq. 4 ) then

                     write(iot,'(a3,'' ('',i4,1x,
     &                          4(1p1e15.7),'' )'')')
     &                    chbd, ibnum, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 6 ) then

                     write(iot,'(a3,'' ('',i4,1x,
     &                          4(1p1e15.7)/10x,2(1p1e15.7),'' )'')')
     &                    chbd, ibnum, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 7 ) then

                     write(iot,'(a3,'' ('',i4,1x,
     &                          4(1p1e15.7)/10x,3(1p1e15.7),'' )'')')
     &                    chbd, ibnum, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 8 ) then

                     write(iot,'(a3,'' ('',i4,1x,
     &                          4(1p1e15.7)/10x,4(1p1e15.7),'' )'')')
     &                    chbd, ibnum, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 9 ) then

                     write(iot,'(a3,'' ('',i4,1x,
     &                          4(1p1e15.7)/10x,4(1p1e15.7),
     &                                     /10x,1(1p1e15.7),'' )'')')
     &                    chbd, ibnum, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 10 ) then

                     write(iot,'(a3,'' ('',i4,1x,
     &                          4(1p1e15.7)/10x,4(1p1e15.7),
     &                                     /10x,2(1p1e15.7),'' )'')')
     &                    chbd, ibnum, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 12 ) then

                     write(iot,'(a3,'' ('',i4,1x,
     &                          4(1p1e15.7)/10x,4(1p1e15.7),
     &                                     /10x,4(1p1e15.7),'' )'')')
     &                    chbd, ibnum, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 13 ) then

                     write(iot,'(a3,'' ('',i4,1x,
     &                          4(1p1e15.7)/10x,4(1p1e15.7),
     &                                     /10x,4(1p1e15.7),
     &                                     /10x,1(1p1e15.7),'' )'')')
     &                    chbd, ibnum, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 30 ) then

                     write(iot,'(a3,'' ('',i4,1x,
     &                          4(1p1e15.7),5(/10x,4(1p1e15.7)))')
     &                    chbd, ibnum, ( bval(i), i = 1, 24 )

                     write(iot,'(10x,6i10,'' )'')')
     &                    ( nint( bval(i) ), i = 25, 30 )

                  end if

            end do

         write(iot,'(/''END''/)')

*-----------------------------------------------------------------------
*        region
*-----------------------------------------------------------------------

            ild1 = 22
            ild0 = 72 - ild1

            ioe = 22

            open(ioe,status='scratch',form='unformatted')

*-----------------------------------------------------------------------

         do i = 1, iregn

            rewind ioe
            read(iod) (chrg(mci+k:mci+k),k=1,ichl(i))

*-----------------------------------------------------------------------

                  ilrm = ichl(i)
                  isqd = 0
                  isrm = 0

  430             isqd = isqd + 1

                  if( ilrm .le. ild0 ) then

                     if( isqd .eq. 1 ) then

                           ildf = ilrm

                        do k = 1, ilrm

                           chdf(k:k) = chrg(mci+k+isrm:mci+k+isrm)

                        end do

                     else

                           write(ioe) ilrm
                           write(ioe)
     &                     (chrg(mci+k+isrm:mci+k+isrm),k=1,ilrm)

                     end if

                  else

                     do k = ild0, 1, -1

                        if( chrg(mci+k+isrm:mci+k+isrm) .eq. ' ' )
     &                  goto 420

                     end do

  420                k2 = k

                     if( isqd .eq. 1 ) then

                           ildf = k2 - 1

                        do k = 1, k2 - 1

                           chdf(k:k) = chrg(mci+k+isrm:mci+k+isrm)

                        end do

                     else

                           write(ioe) k2 - 1
                           write(ioe)
     &                     (chrg(mci+k+isrm:mci+k+isrm),k=1,k2-1)

                     end if

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 430

                  end if

*-----------------------------------------------------------------------

                  chsc = chsm(i)
                  call chcptl(chsc,1,3)
                  call chcptl(chdf,1,80)

               if( idmg(i) .eq. -1 ) then

                  idcg = -1000

               else

                  idcg = idmg(i)

               end if

               write(iot,'(1x,a3,'':  '',a3,'':''i6'' :   '',200a1)')
     &                     chsc, chsc, idcg,
     &                     (chdf(j:j),j=1,ildf)


            if( isqd .gt. 1 ) then

                  rewind ioe

               do k = 2, isqd

                  read(ioe) ildf
                  read(ioe) (chdf(j:j),j=1,ildf)

                  call chcptl(chdf,1,80)

                  write(iot,'(200a1)') (cblan(j:j),j=1,ild1),
     &                    (chdf(j:j),j=1,ildf)

               end do

            end if

*-----------------------------------------------------------------------

         end do

            close(ioe)

         write(iot,'(/''$END GEOMETRY''/)')

*-----------------------------------------------------------------------
*     read and write cgview setting file
*-----------------------------------------------------------------------

            inquire( file = chfn(3), exist = exex )

            if( exex .eqv. .false. ) then

               if( icfn(3) .eq. 1 ) then

                  write( 6,'(/'' Error : cgview setting file'',
     &                         '' does not exist.''/
     &                         '' file name = '',100a1)')
     &                  ( chfn(3)(i:i), i = 1, ilfn(3) )

               end if

                  goto 600

            end if

            ios = 22

            open(ios, file = chfn(3), status = 'old' )

*-----------------------------------------------------------------------

  400    continue

            read(ios,'(a200)',end = 500 ) chin

            call chlngt(chin,200,i1,i2)

            if( i1 .eq. 0 .and. i2 .eq. 0 ) then

               write(iot,'()')

            else

               if( chin(i1:i2) .eq. 'QUIT' .or.
     &             chin(i1:i2) .eq. 'quit' ) then

                  write(iot,'(''QUIT'')')
                  write(iot,'(''/'')')

                  goto 500

               else

                  write(iot,'(200a1)') ( chin(i:i), i = 1, i2 )

               end if

            end if

            goto 400

  500    continue

         close(ios)

*-----------------------------------------------------------------------

  600    continue

         close(iot)
         close(itby)

*-----------------------------------------------------------------------

      call moddas_deallocate_cha(chrg)
      return
      end

