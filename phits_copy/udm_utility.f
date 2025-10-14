************************************************************************
      module udm_Utility
************************************************************************

      use udm_Parameter
      implicit real*8(a-h,o-z)
      contains

************************************************************************
      subroutine udm_initialize_Utility
************************************************************************
      common /paran/ icfn(100), ilfn(100), chfn(100) ! phits paths
      common /wparm/  swtm(20), wc1(20), wc2(20) !FURUTA
      character chfn*200
C     ------------------------------------------------------------------
      integer*8 :: iranji64 ! S.H. xorshift (2020.2.6)
      common /randm4/ rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,
     &                rnrtc,nstrid,inif, iranji64
      integer*8 :: iransb64 ! S.H. xorshift (2020.2.6)
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)
C     ------------------------------------------------------------------
      udm_phits_path=chfn(1)(1:ilfn(1))
C     ------------------------------------------------------------------
C     To make random numbers available at the time of 'initialize' of the udm module.
C     If this is not done, random numbers cannot be used until a non-zero value is assigned to iransb64 in ovly12.f.
      iransb64 = iranji64
      do i=1,1000
        tmp=unirn(dummy)
      enddo
C     ------------------------------------------------------------------
      if(iudmodel .gt. 0) then
        i_udm_tmp = 0
        do i=1,19
          if(i .eq. 12 .or. i .eq. 13) cycle ! 12 and 13 are 0 by default.
          if(wc1(i) .ne. 0.5d0) cycle ! Skip for set values.
          if(i .eq. 2) then
            ! If neutron's weight cutoff is set to 0, calculation time increase so much.
            bias_max=1d0
            if(pnimul > bias_max) bias_max=pnimul
            if(gmumul > bias_max) bias_max=gmumul
            do iUDM = 1, udm_int_num
              if(udm_bias(iUDM) > bias_max) bias_max=udm_bias(iUDM)
            enddo
            wc1(i) = 0.1d0/bias_max
          else
            wc1(i) = 0.0d0
          endif
          wc2(i) = 0.5d0 * wc1(i) * swtm(i)
          i_udm_tmp = i_udm_tmp + 1
        enddo
        if(i_udm_tmp .ne. 0)
     & write(*,*) "*** wc1(i) values are modified to alomost zero."
      endif

C     ------------------------------------------------------------------
      end subroutine udm_initialize_Utility

************************************************************************
      function udm_phits_path_fun()
************************************************************************
      common /paran/ icfn(100), ilfn(100), chfn(100) ! phits paths
      character chfn*100
      character udm_phits_path_fun*200
      udm_phits_path_fun=chfn(1)(1:ilfn(1))
      return
      end function


************************************************************************
*                                                                      *
      subroutine allocate_udinteract
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      allocate(udm_int_name     (udm_module_nMax))
      allocate(udm_bias         (udm_module_nMax))
      allocate(udm_int_param    (udm_module_nMax, udm_int_param_nMax))
      allocate(udm_int_param_str(udm_module_nMax, udm_int_param_nMax))

      do i=1,udm_module_nMax
        udm_int_name(i) = " "
        udm_bias    (i) = 0d0
        do j=1,udm_int_param_nMax
          udm_int_param    (i,j) = 0d0
          udm_int_param_str(i,j) = " "
        enddo
      enddo

      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine read_udinteract(jsn,jsi,dsin,idsi,ill,ilf,
     &                           jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [User Defined Interaction] section of input files         *
*       modified by S.Abe on 2021/11/04                                *
*                                                                      *
************************************************************************

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

      character(len=200) :: tokens(42)
      integer ntokens

*-----------------------------------------------------------------------

      ierr  = 0

*-----------------------------------------------------------------------

      call allocate_udinteract
      i_udm_tmp = 0

*-----------------------------------------------------------------------

  140 continue

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

      call readl(jsn,jsi,dsin,idsi,ill,ilf,'#%!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) return
      if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------

      if( chlw(i1:i1+4-1) .eq. 'name' ) goto 140      ! skip data sequence line

      if( chlw(i1:i1+5-1)=='n_int' .or. chlw(i1:i1+6-1)=='number' )then
        print*,"'n_int' and 'number' have been deprecated. The number"
     &       //" of interactions used is now automatically determined."
        goto 140
      endif

*-----------------------------------------------------------------------
*     end of section
*-----------------------------------------------------------------------

      if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

         jpn = 1
         return

      end if

*-----------------------------------------------------------------------

      iudmodel=1

*-----------------------------------------------------------------------

      i_udm_tmp = i_udm_tmp + 1
      udm_int_num=i_udm_tmp

*-----------------------------------------------------------------------
*     read name
*-----------------------------------------------------------------------

      ic2  = i1
      ic = jnumc(chlw,ic2,i3)

      ntri = ic
      ntrf = 0
      ic = ic - 1

  151 ic = ic + 1

      if( chlw(ic:ic) .eq. ' ' ) then
         ntrf = ic - 1
         goto 161
      else if( ic .eq. i3 ) then
         goto 995
      else
         goto 151
      endif

  161 continue

      if( ntrf .lt. ntri ) goto 995
      if( ntrf - ntri + 1 .gt. 200 ) goto 994

      ii = 0
      do jj = ntri, ntrf
         ii = ii + 1
         udm_int_name(i_udm_tmp)(ii:ii) = chin(jj:jj)
      enddo

*-----------------------------------------------------------------------
*     read bias
*-----------------------------------------------------------------------

      call snum(chlw,ic,i3,ic2,cvvv,ierr)

      if( ierr .ne. 0 ) goto 999

      udm_bias(i_udm_tmp) = cvvv

      if( ic2 .ge. i3 ) goto 181

*-----------------------------------------------------------------------
*     read parameters (as a string format)
*-----------------------------------------------------------------------

      call tokenize_chlw(chlw, tokens, ntokens)
C     tokens(1): name
C     tokens(2): bias
      if(ntokens > 2) then
        do iTok=3, ntokens
          ipar=iTok-2
          udm_int_param_str(i_udm_tmp,ipar)=tokens(iTok)
        enddo
      endif

*-----------------------------------------------------------------------
*     read parameters (in double-precision format)
*-----------------------------------------------------------------------

      ipar = 0

  171 ipar = ipar + 1

      ic = jnumc(chlw,ic2,i3)

      call snum(chlw,ic,i3,ic2,cvvv,ierr)

C     if( ierr .ne. 0 ) goto 999
      if( ierr .eq. 0 ) then
        udm_int_param(i_udm_tmp,ipar) = cvvv
      else ! It is a string.
C       Manually update the index of the next character's start position
        ic2=ic+len_trim(udm_int_param_str(i_udm_tmp,ipar))+1
      endif

C     if( ic2 .ge. i3 ) goto 181
      if( ic2-1 .ge. i3 ) goto 181
C     "-1" prevents it from being ignored when the last number is a single-digit
C     natural number and there is exactly one preceding space.

      goto 171

  181 continue

*-----------------------------------------------------------------------

      goto 140

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  993 continue

         m_err = 'Description of bias is wrong'
         ErrCha = ''
         ErrID = 'L:287/R:read_udinteract/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  994 continue

         m_err = 'Name should be less than 200 characters'//
     &           ' in [user defined interaction]'
         ErrCha = ''
         ErrID = 'L:298/R:read_udinteract/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Only Name is defined in [user defined interaction]'
         ErrCha = ''
         ErrID = 'L:310/R:read_udinteract/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'number is not defined in the first line'
         ErrCha = ''
         ErrID = 'L:322/R:read_udinteract/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return


*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of parameter is wrong'
         ErrCha = ''
         ErrID = 'L:335/R:read_udinteract/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Value of parameter is wrong'
         ErrCha = ''
         ErrID = 'L:347/R:read_udinteract/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [user defined interaction] is wrong.'
         ErrCha = ''
         ErrID = 'L:359/R:read_udinteract/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine

*-----------------------------------------------------------------------



************************************************************************
      subroutine tokenize_chlw(chlw, tokens, ntokens)
*     The string chlw consists of multiple substrings separated by one or more spaces.
*     A substring enclosed in double or single quotation marks is treated as a single token, even if it contains spaces.
*     chlw is split into individual tokens and stored in the array tokens.
*     The number of obtained tokens is set to ntokens.
************************************************************************
      implicit none
      character(*), intent(in)  :: chlw
      character(len=*), dimension(42) :: tokens
      integer, intent(out)      :: ntokens

      integer :: i, lench, start
      character(len=1) :: ch, delim
      character(len=nStrMax) :: token  ! Adjust the maximum length of each token as needed.

      ntokens = 0
      lench = len_trim(chlw)
      i = 1

      do while (i <= lench)
         ! Skip spaces.
         do while (i <= lench .and. chlw(i:i) == ' ')
            i = i + 1
         end do
         if (i > lench) exit

         ch = chlw(i:i)
         if (ch == '"' .or. ch == "'") then
            ! Tokens enclosed in quotation marks.
            delim = ch
            start = i + 1
            i = start
            do while (i <= lench .and. chlw(i:i) /= delim)
               i = i + 1
            end do
            if (i <= lench) then
               token = chlw(start:i-1)
               i = i + 1  ! Skip the closing quotation mark.
            else
               token = chlw(start:lench)  ! If no matching quotation mark is found.
            end if
         else
            ! A normal token up to the next space.
            start = i
            do while (i <= lench .and. chlw(i:i) /= ' ')
               i = i + 1
            end do
            token = chlw(start:i-1)
         end if

         ntokens = ntokens + 1
         if (ntokens > 42) exit
         tokens(ntokens) = token
      end do

      end subroutine tokenize_chlw
************************************************************************





************************************************************************
*                                                                      *
      subroutine allocate_udparticle
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      allocate(udm_part_name     (udm_module_nMax))
      allocate(udm_part_kf       (udm_module_nMax))
      allocate(udm_part_param    (udm_module_nMax, udm_part_param_nMax))
      allocate(udm_part_param_str(udm_module_nMax, udm_part_param_nMax))

      do i=1,udm_module_nMax
        udm_part_name(i) = " "
        udm_part_kf  (i) = 0
        do j=1,udm_part_param_nMax
          udm_part_param    (i,j) = 0d0
          udm_part_param_str(i,j) = " "
        enddo
      enddo

      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine read_udpart(jsn,jsi,dsin,idsi,ill,ilf,
     &                           jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [User Defined Particle] section of input files            *
*       modified by S.Abe on 2021/11/04                                *
*                                                                      *
************************************************************************

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

      character(len=200) :: tokens(42)
      integer ntokens

*-----------------------------------------------------------------------

      ierr  = 0

*-----------------------------------------------------------------------

      call allocate_udparticle
      i_udm_tmp = 0

*-----------------------------------------------------------------------

  140 continue

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

      call readl(jsn,jsi,dsin,idsi,ill,ilf,'#%!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) return
      if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------

      if( chlw(i1:i1+4-1) .eq. 'name' ) goto 140      ! skip data sequence line

      if( chlw(i1:i1+6-1)=='n_part' .or. chlw(i1:i1+6-1)=='number' )then
        print*,"'n_part' and 'number' have been deprecated. The number"
     &       //" of interactions used is now automatically determined."
        goto 140
      endif

*-----------------------------------------------------------------------
*     end of section
*-----------------------------------------------------------------------

      if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

         jpn = 1
         return

      end if

*-----------------------------------------------------------------------

      iudmodel=1

*-----------------------------------------------------------------------

      i_udm_tmp = i_udm_tmp + 1
      udm_part_num=i_udm_tmp

*-----------------------------------------------------------------------
*     read name
*-----------------------------------------------------------------------

      ic2  = i1
      ic = jnumc(chlw,ic2,i3)

      ntri = ic
      ntrf = 0
      ic = ic - 1

  151 ic = ic + 1

      if( chlw(ic:ic) .eq. ' ' ) then
         ntrf = ic - 1
         goto 161
      else if( ic .eq. i3 ) then
         goto 995
      else
         goto 151
      endif

  161 continue

      if( ntrf .lt. ntri ) goto 995
      if( ntrf - ntri + 1 .gt. 200 ) goto 994

      ii = 0
      do jj = ntri, ntrf
         ii = ii + 1
         udm_part_name(i_udm_tmp)(ii:ii) = chin(jj:jj)
      enddo

*-----------------------------------------------------------------------
*     read kf-code
*-----------------------------------------------------------------------

      call snum(chlw,ic,i3,ic2,cvvv,ierr)

      if( ierr .ne. 0 ) goto 999

      udm_part_kf(i_udm_tmp) = nint( cvvv )

      if( ic2 .ge. i3 ) goto 181

*-----------------------------------------------------------------------
*     read parameters (as a string format)
*-----------------------------------------------------------------------

      call tokenize_chlw(chlw, tokens, ntokens)
C     tokens(1): name
C     tokens(2): kf-code
      if(ntokens > 2) then
        do iTok=3, ntokens
          ipar=iTok-2
          udm_part_param_str(i_udm_tmp,ipar)=tokens(iTok)
        enddo
      endif

*-----------------------------------------------------------------------
*     read parameters (in double-precision format)
*-----------------------------------------------------------------------

      ipar = 0

  171 ipar = ipar + 1

      ic = jnumc(chlw,ic2,i3)

      call snum(chlw,ic,i3,ic2,cvvv,ierr)

C     if( ierr .ne. 0 ) goto 999
      if( ierr .eq. 0 ) then
        udm_part_param(i_udm_tmp,ipar) = cvvv
      else ! It is a string.
C       Manually update the index of the next character's start position
        ic2=ic+len_trim(udm_part_param_str(i_udm_tmp,ipar))+1
      endif

C     if( ic2 .ge. i3 ) goto 181
      if( ic2-1 .ge. i3 ) goto 181
C     "-1" prevents it from being ignored when the last number is a single-digit
C     natural number and there is exactly one preceding space.

      goto 171

  181 continue

*-----------------------------------------------------------------------

      goto 140

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  994 continue

         m_err = 'Name should be less than 200 characters'//
     &           ' in [user defined particle]'
         ErrCha = ''
         ErrID = 'L:654/R:read_udpart/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Only Name is defined in [user defined particle]'
         ErrCha = ''
         ErrID = 'L:666/R:read_udpart/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'number is not defined in the first line'
         ErrCha = ''
         ErrID = 'L:678/R:read_udpart/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return


*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of parameter is wrong'
         ErrCha = ''
         ErrID = 'L:691/R:read_udpart/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Value of parameter is wrong'
         ErrCha = ''
         ErrID = 'L:703/R:read_udpart/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [user defined particle] is wrong.'
         ErrCha = ''
         ErrID = 'L:715/R:read_udpart/F:udm_utility.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine

*-----------------------------------------------------------------------











************************************************************************
      subroutine fill_mat_info(mk)
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)
*-----------------------------------------------------------------------
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      dimension dnel(1), denh(1), zz(1), a(1), dens(1)
      equivalence ( das, dnel, denh, zz, a, dens )
      common /kmat1g/ kmat(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

*-----------------------------------------------------------------------

      imat = mk

      num_nuclide=1

      do m = jmd(1+mk), jmd(1+mk+1) - 1
        ie = lme(2,m)
        izm = iza(m) / 1000       ! Z: proton number
        ims = iza(m) - 1000 * izm ! A: mass number
        icntm = m - jmd(1+mk) + 1
C         write(*,"(I5,I5,f12.6)") izm, ims, fme(m)

        ! -------------------------
        if( ims .eq. 0 ) then
          !--------------------------------------------------------------
          lemm = nint( dnel(kmat(imat)+1) )
          unnorm_AR_sum = 0.0d0 ! averaged cross section with the Unnormalized Abundance Ratio
          do lem = 1, lemm
            izt = zz(kmat(imat)+(lem-1)*3+3)
            if( izt .eq. izm ) then
              unnorm_AR = dens(kmat(imat)+(lem-1)*3+5) ! Unnormalized Abundance Ratio
              unnorm_AR_sum = unnorm_AR_sum + unnorm_AR
            end if
          end do
          !--------------------------------------------------------------
          do lem = 1, lemm
            izt = zz(kmat(imat)+(lem-1)*3+3)
            if( izt .eq. izm ) then
              unnorm_AR = dens(kmat(imat)+(lem-1)*3+5)
              sekt = sekt + unnorm_AR
              ims = nint( a(kmat(imat)+(lem-1)*3+4) ) ! mass number
C               write(*,"(A,I4,I5,f12.6)") "-",izt,ims,
C      &          unnorm_AR/unnorm_AR_sum*fme(m)
              mat_Z(num_nuclide)=izt
              mat_A(num_nuclide)=ims
              mat_ratio(num_nuclide)=unnorm_AR/unnorm_AR_sum*fme(m)
              num_nuclide=num_nuclide+1
              if(num_nuclide .gt. num_nuclide_max) then
                print*,
     &               "[CAUTION] Increase 'num_nuclide_max' 200 or more."
                return
              endif
            end if
          end do
        !--------------------------------------------------------------
        else
          mat_Z(num_nuclide)=izm
          mat_A(num_nuclide)=ims
          mat_ratio(num_nuclide)=fme(m)
          num_nuclide=num_nuclide+1
          if(num_nuclide .gt. num_nuclide_max) then
            print*,"[CAUTION] Increase 'num_nuclide_max' 200 or more."
            return
          endif
          endif
      end do
      num_nuclide=num_nuclide-1

      end subroutine fill_mat_info

************************************************************************
      function get_random_0to1()
************************************************************************
      implicit real*8 (a-h,o-z)
      get_random_0to1 = unirn(dummy)
      return
      end function

************************************************************************
      function get_random(rmin,rmax)
************************************************************************
      implicit real*8 (a-h,o-z)
      get_random = rmin+(rmax-rmin)*unirn(dummy)
      return
      end function

************************************************************************
      function get_random_int(imin,imax)
************************************************************************
      implicit real*8 (a-h,o-z)
      integer get_random_int
      do
        get_random_int = imin+int((imax-imin+1)*unirn(dummy))
        if(imin<=get_random_int .and. get_random_int<=imax) return
      enddo
      end function

************************************************************************
      subroutine two_body_decay_uniform(kf1,kf2)
************************************************************************
      implicit real*8 (a-h,o-z)
      integer kf1, kf2
      double precision m0, m1, m2
      m0=get_mass(udm_kf_for_21)
      m1=get_mass(kf1)
      m2=get_mass(kf2)
      if(m0 < m1+m2) then
        print*,"*** [CAUTION] two_body_decay_uniform failed"
        print*,"*** ",udm_kf_for_21,"-->",kf1,kf2
        return
      endif
      call dklos2(
     &      get_ityp(udm_kf_for_21),
     &               udm_kf_for_21,
     &               udm_Kin,
     &      get_ityp(kf1), kf1,
     &      get_ityp(kf2), kf2)
*-----------------------------------------------------------------------
      set_decay_success=.true.
*-----------------------------------------------------------------------
      end subroutine two_body_decay_uniform

************************************************************************
      subroutine three_body_decay_uniform(kf1,kf2,kf3)
************************************************************************
      implicit real*8 (a-h,o-z)
      integer kf1, kf2, kf3
      double precision m0, m1, m2, m3
      m0=get_mass(udm_kf_for_21)
      m1=get_mass(kf1)
      m2=get_mass(kf2)
      m3=get_mass(kf3)
      if(m0 < m1+m2+m3) then
        print*,"*** [CAUTION] three_body_decay_uniform failed"
        print*,"*** ",udm_kf_for_21,"-->",kf1,kf2,kf3
        return
      endif
      call dklos3(
     &      get_ityp(udm_kf_for_21),
     &               udm_kf_for_21,
     &               udm_Kin,
     &      get_ityp(kf1), kf1,
     &      get_ityp(kf2), kf2,
     &      get_ityp(kf3), kf3)
*-----------------------------------------------------------------------
      set_decay_success=.true.
*-----------------------------------------------------------------------
      end subroutine three_body_decay_uniform

************************************************************************
      function get_mass(kf)
************************************************************************
      double precision get_mass
      integer ityp, kf
      ityp=get_ityp(kf)
      get_mass=rmtyp(ityp,kf)
      return
      end function

************************************************************************
      subroutine fill_final_state
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      real*8 event(20)

*-----------------------------------------------------------------------

      include 'param00.inc'
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      integer iZN_i(2)
      double precision mass_i_MeV
      double precision absp_i
      double precision total_energy_MeV
      double precision kinetic_energy_MeV

*-----------------------------------------------------------------------

      numpal = 0
      rumpal = 0.d0

      do i=1,set_final_state_number

        nclsts = i

        kf_i  =set_kf(i)
        ityp_i=get_ityp(kf_i)
        iZN_i =get_proton_neutron_number(kf_i)

        iclusts(nclsts) = ipatf(ityp_i, kf_i)

        jclusts(0,nclsts) = 0
        jclusts(1,nclsts) = iZN_i(1)
        jclusts(2,nclsts) = iZN_i(2)
        jclusts(3,nclsts) = ityp_i
        jclusts(4,nclsts) = 0
        jclusts(5,nclsts) = ichgf(ityp_i, kf_i)
        jclusts(6,nclsts) = ibryf(ityp_i, kf_i)
        jclusts(7,nclsts) = kf_i
        jclusts(8,nclsts) = set_isomer_level(i)

        mass_i_MeV = rmtyp(ityp_i, kf_i)
        total_energy_MeV = set_Total_Energy_in_MeV(i)
        kinetic_energy_MeV = total_energy_MeV - mass_i_MeV
        absp_i = sqrt( set_Px_in_MeV(i)**2
     &               + set_Py_in_MeV(i)**2
     &               + set_Pz_in_MeV(i)**2 )
        qclusts(0,nclsts)  = 0.0
        qclusts(1,nclsts)  = set_Px_in_MeV(i) / absp_i
        qclusts(2,nclsts)  = set_Py_in_MeV(i) / absp_i
        qclusts(3,nclsts)  = set_Pz_in_MeV(i) / absp_i
        qclusts(4,nclsts)  = total_energy_MeV / 1000d0 ! GeV
        qclusts(5,nclsts)  = mass_i_MeV       / 1000d0 ! GeV
        qclusts(6,nclsts)  = set_excitation_energy_in_MeV(i)
        qclusts(7,nclsts)  = kinetic_energy_MeV
        qclusts(8,nclsts)  = 1.0 ! wgt / wg0
        qclusts(9,nclsts)  = 0.0
        qclusts(10,nclsts) = 0.0d0
        qclusts(11,nclsts) = 0.0d0
        qclusts(12,nclsts) = 0.0d0

        numpal(ityp_i) = numpal(ityp_i) + 1
        rumpal(ityp_i) = rumpal(ityp_i) + qclusts(8,nclsts)

      enddo
*-----------------------------------------------------------------------
      set_decay_success=.true.
*-----------------------------------------------------------------------
      end subroutine fill_final_state

************************************************************************
      subroutine fill_final_state_for_decay
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      real*8 event(20)

*-----------------------------------------------------------------------

      include 'param00.inc'
      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)
      integer iZN_i(2)
      double precision mass_i_MeV
      double precision absp_i
      double precision total_energy_MeV
      double precision kinetic_energy_MeV

*-----------------------------------------------------------------------

      nclst=set_final_state_number

      do i=1,set_final_state_number

        kf_i  =set_kf(i)
        ityp_i=get_ityp(kf_i)
        iZN_i =get_proton_neutron_number(kf_i)

        numpat(ityp_i)  = numpat(ityp_i) + 1
        rumpat(ityp_i)  = rumpat(ityp_i) + 1

        iclust(i) = ipatf(ityp_i, kf_i)

        jclust(0,i) = 0
        jclust(1,i) = iZN_i(1)
        jclust(2,i) = iZN_i(2)
        jclust(3,i) = ityp_i
        jclust(4,i) = 0
        jclust(5,i) = ichgf(ityp_i, kf_i)
        jclust(6,i) = ibryf(ityp_i, kf_i)
        jclust(7,i) = kf_i
        jclust(8,i) = set_isomer_level(i)

        mass_i_MeV = rmtyp(ityp_i, kf_i)
        total_energy_MeV = set_Total_Energy_in_MeV(i)
        kinetic_energy_MeV = total_energy_MeV - mass_i_MeV
        absp_i = sqrt( set_Px_in_MeV(i)**2
     &               + set_Py_in_MeV(i)**2
     &               + set_Pz_in_MeV(i)**2 )
        qclust(0,i)  = 0.0
        qclust(1,i)  = set_Px_in_MeV(i) / absp_i
        qclust(2,i)  = set_Py_in_MeV(i) / absp_i
        qclust(3,i)  = set_Pz_in_MeV(i) / absp_i
        qclust(4,i)  = total_energy_MeV / 1000d0 ! GeV
        qclust(5,i)  = mass_i_MeV       / 1000d0 ! GeV
        qclust(6,i)  = set_excitation_energy_in_MeV(i)
        qclust(7,i)  = kinetic_energy_MeV
        qclust(8,i)  = 1.0 ! wgt / wg0
        qclust(9,i)  = 0.0
        qclust(10,i) = 0.0d0
        qclust(11,i) = 0.0d0
        qclust(12,i) = 0.0d0

      enddo
*-----------------------------------------------------------------------
      set_decay_success=.true.
*-----------------------------------------------------------------------
      end subroutine fill_final_state_for_decay

************************************************************************
      subroutine multiply_weight_by(ratio)
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
*-----------------------------------------------------------------------
      wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * ratio
      end subroutine multiply_weight_by

************************************************************************
      function udm_photonuclear_Xsec(Z,N,Kin)
      ! Z: atomic number
      ! N: neutron number
      ! Kin: incident photon energy (MeV)
************************************************************************
      implicit none
      double precision udm_photonuclear_Xsec
      double precision phxs
      integer Z,N
      double precision Kin
      udm_photonuclear_Xsec=phxs(Z,N,Kin)
      return
      end function

************************************************************************
      function udm_int2char(i)
************************************************************************
      implicit none
      character(len=50) :: udm_int2char
      character(len=10) :: format
      integer i
C     ---------
      if    (0<=i .and. i <10) then; format="(I1)"
      elseif(10<=i .and. i <100) then; format="(I2)"
      elseif(100<=i .and. i <1000) then; format="(I3)"
      elseif(1000<=i .and. i <10000) then; format="(I4)"
      elseif(10000<=i .and. i <100000) then; format="(I5)"
      elseif(100000<=i .and. i <1000000) then; format="(I6)"
      elseif(1000000<=i .and. i <10000000) then; format="(I7)"
      elseif(10000000<=i .and. i <100000000) then; format="(I8)"
      elseif(100000000<=i .and. i <1000000000) then; format="(I9)"
C     ---------
      elseif(-9<=i .and. i <0) then; format="(I2)"
      elseif(-99<=i .and. i <-9) then; format="(I3)"
      elseif(-999<=i .and. i <-99) then; format="(I4)"
      elseif(-9999<=i .and. i <-999) then; format="(I5)"
      elseif(-99999<=i .and. i <-9999) then; format="(I6)"
      elseif(-999999<=i .and. i <-99999) then; format="(I7)"
      elseif(-9999999<=i .and. i <-999999) then; format="(I8)"
      elseif(-99999999<=i .and. i <-9999999) then; format="(I9)"
      elseif(-999999999<=i .and. i <-99999999) then; format="(I10)"
C     ---------
      else
        print*,"error: function udm_int2char"
        stop
      endif
C     ---------
      write (udm_int2char,format) i
C     ---------
      return
      end function

************************************************************************
      function udm_double2char(format,x)
************************************************************************
      implicit none
      character(len=50) :: udm_double2char
      character(len=*)  :: format
      double precision x
C     ---------
      write (udm_double2char,format) x
C     ---------
      return
      end function

************************************************************************
      subroutine udm_get_atomic_information(Z,A,c,c2)
************************************************************************
      implicit none
      integer Z
      double precision A
      character(len=3 ) :: c
      character(len=13) :: c2
C     https://www.angelo.edu/faculty/kboudrea/periodic/structure_mass.htm
      if(Z==1  )then;A=1.00797  ;c="H"  ;c2="Hydrogen"     ;return;endif
      if(Z==2  )then;A=4.00260  ;c="He" ;c2="Helium"       ;return;endif
      if(Z==3  )then;A=6.941    ;c="Li" ;c2="Lithium"      ;return;endif
      if(Z==4  )then;A=9.01218  ;c="Be" ;c2="Beryllium"    ;return;endif
      if(Z==5  )then;A=10.81    ;c="B"  ;c2="Boron"        ;return;endif
      if(Z==6  )then;A=12.011   ;c="C"  ;c2="Carbon"       ;return;endif
      if(Z==7  )then;A=14.0067  ;c="N"  ;c2="Nitrogen"     ;return;endif
      if(Z==8  )then;A=15.9994  ;c="O"  ;c2="Oxygen"       ;return;endif
      if(Z==9  )then;A=18.998403;c="F"  ;c2="Fluorine"     ;return;endif
      if(Z==10 )then;A=20.179   ;c="Ne" ;c2="Neon"         ;return;endif
      if(Z==11 )then;A=22.98977 ;c="Na" ;c2="Sodium"       ;return;endif
      if(Z==12 )then;A=24.305   ;c="Mg" ;c2="Magnesium"    ;return;endif
      if(Z==13 )then;A=26.98154 ;c="Al" ;c2="Aluminum"     ;return;endif
      if(Z==14 )then;A=28.0855  ;c="Si" ;c2="Silicon"      ;return;endif
      if(Z==15 )then;A=30.97376 ;c="P"  ;c2="Phosphorus"   ;return;endif
      if(Z==16 )then;A=32.06    ;c="S"  ;c2="Sulfur"       ;return;endif
      if(Z==17 )then;A=35.453   ;c="Cl" ;c2="Chlorine"     ;return;endif
      if(Z==18 )then;A=39.948   ;c="Ar" ;c2="Argon"        ;return;endif
      if(Z==19 )then;A=39.0983  ;c="K"  ;c2="Potassium"    ;return;endif
      if(Z==20 )then;A=40.08    ;c="Ca" ;c2="Calcium"      ;return;endif
      if(Z==21 )then;A=44.9559  ;c="Sc" ;c2="Scandium"     ;return;endif
      if(Z==22 )then;A=47.90    ;c="Ti" ;c2="Titanium"     ;return;endif
      if(Z==23 )then;A=50.9415  ;c="V"  ;c2="Vanadium"     ;return;endif
      if(Z==24 )then;A=51.996   ;c="Cr" ;c2="Chromium"     ;return;endif
      if(Z==25 )then;A=54.9380  ;c="Mn" ;c2="Manganese"    ;return;endif
      if(Z==26 )then;A=55.847   ;c="Fe" ;c2="Iron"         ;return;endif
      if(Z==27 )then;A=58.9332  ;c="Co" ;c2="Cobalt"       ;return;endif
      if(Z==28 )then;A=58.70    ;c="Ni" ;c2="Nickel"       ;return;endif
      if(Z==29 )then;A=63.546   ;c="Cu" ;c2="Copper"       ;return;endif
      if(Z==30 )then;A=65.38    ;c="Zn" ;c2="Zinc"         ;return;endif
      if(Z==31 )then;A=69.72    ;c="Ga" ;c2="Gallium"      ;return;endif
      if(Z==32 )then;A=72.59    ;c="Ge" ;c2="Germanium"    ;return;endif
      if(Z==33 )then;A=74.9216  ;c="As" ;c2="Arsenic"      ;return;endif
      if(Z==34 )then;A=78.96    ;c="Se" ;c2="Selenium"     ;return;endif
      if(Z==35 )then;A=79.904   ;c="Br" ;c2="Bromine"      ;return;endif
      if(Z==36 )then;A=83.80    ;c="Kr" ;c2="Krypton"      ;return;endif
      if(Z==37 )then;A=85.4678  ;c="Rb" ;c2="Rubidium"     ;return;endif
      if(Z==38 )then;A=87.62    ;c="Sr" ;c2="Strontium"    ;return;endif
      if(Z==39 )then;A=88.9059  ;c="Y"  ;c2="Yttrium"      ;return;endif
      if(Z==40 )then;A=91.22    ;c="Zr" ;c2="Zirconium"    ;return;endif
      if(Z==41 )then;A=92.9064  ;c="Nb" ;c2="Niobium"      ;return;endif
      if(Z==42 )then;A=95.94    ;c="Mo" ;c2="Molybdenum"   ;return;endif
      if(Z==43 )then;A=98       ;c="Tc" ;c2="Technetium"   ;return;endif
      if(Z==44 )then;A=101.07   ;c="Ru" ;c2="Ruthenium"    ;return;endif
      if(Z==45 )then;A=102.9055 ;c="Rh" ;c2="Rhodium"      ;return;endif
      if(Z==46 )then;A=106.4    ;c="Pd" ;c2="Palladium"    ;return;endif
      if(Z==47 )then;A=107.868  ;c="Ag" ;c2="Silver"       ;return;endif
      if(Z==48 )then;A=112.41   ;c="Cd" ;c2="Cadmium"      ;return;endif
      if(Z==49 )then;A=114.82   ;c="In" ;c2="Indium"       ;return;endif
      if(Z==50 )then;A=118.69   ;c="Sn" ;c2="Tin"          ;return;endif
      if(Z==51 )then;A=121.75   ;c="Sb" ;c2="Antimony"     ;return;endif
      if(Z==52 )then;A=127.60   ;c="Te" ;c2="Tellurium"    ;return;endif
      if(Z==53 )then;A=126.9045 ;c="I"  ;c2="Iodine"       ;return;endif
      if(Z==54 )then;A=131.30   ;c="Xe" ;c2="Xenon"        ;return;endif
      if(Z==55 )then;A=132.9054 ;c="Cs" ;c2="Cesium"       ;return;endif
      if(Z==56 )then;A=137.33   ;c="Ba" ;c2="Barium"       ;return;endif
      if(Z==57 )then;A=138.9055 ;c="La" ;c2="Lanthanum"    ;return;endif
      if(Z==58 )then;A=140.12   ;c="Ce" ;c2="Cerium"       ;return;endif
      if(Z==59 )then;A=140.9077 ;c="Pr" ;c2="Praseodymium" ;return;endif
      if(Z==60 )then;A=144.24   ;c="Nd" ;c2="Neodymium"    ;return;endif
      if(Z==61 )then;A=145      ;c="Pm" ;c2="Promethium"   ;return;endif
      if(Z==62 )then;A=150.4    ;c="Sm" ;c2="Samarium"     ;return;endif
      if(Z==63 )then;A=151.96   ;c="Eu" ;c2="Europium"     ;return;endif
      if(Z==64 )then;A=157.25   ;c="Gd" ;c2="Gadolinium"   ;return;endif
      if(Z==65 )then;A=158.9254 ;c="Tb" ;c2="Terbium"      ;return;endif
      if(Z==66 )then;A=162.50   ;c="Dy" ;c2="Dysprosium"   ;return;endif
      if(Z==67 )then;A=164.9304 ;c="Ho" ;c2="Holmium"      ;return;endif
      if(Z==68 )then;A=167.26   ;c="Er" ;c2="Erbium"       ;return;endif
      if(Z==69 )then;A=168.9342 ;c="Tm" ;c2="Thulium"      ;return;endif
      if(Z==70 )then;A=173.04   ;c="Yb" ;c2="Ytterbium"    ;return;endif
      if(Z==71 )then;A=174.967  ;c="Lu" ;c2="Lutetium"     ;return;endif
      if(Z==72 )then;A=178.49   ;c="Hf" ;c2="Hafnium"      ;return;endif
      if(Z==73 )then;A=180.9479 ;c="Ta" ;c2="Tantalum"     ;return;endif
      if(Z==74 )then;A=183.85   ;c="W"  ;c2="Tungsten"     ;return;endif
      if(Z==75 )then;A=186.207  ;c="Re" ;c2="Rhenium"      ;return;endif
      if(Z==76 )then;A=190.2    ;c="Os" ;c2="Osmium"       ;return;endif
      if(Z==77 )then;A=192.22   ;c="Ir" ;c2="Iridium"      ;return;endif
      if(Z==78 )then;A=195.09   ;c="Pt" ;c2="Platinum"     ;return;endif
      if(Z==79 )then;A=196.9665 ;c="Au" ;c2="Gold"         ;return;endif
      if(Z==80 )then;A=200.59   ;c="Hg" ;c2="Mercury"      ;return;endif
      if(Z==81 )then;A=204.37   ;c="Tl" ;c2="Thallium"     ;return;endif
      if(Z==82 )then;A=207.2    ;c="Pb" ;c2="Lead"         ;return;endif
      if(Z==83 )then;A=208.9804 ;c="Bi" ;c2="Bismuth"      ;return;endif
      if(Z==84 )then;A=209      ;c="Po" ;c2="Polonium"     ;return;endif
      if(Z==85 )then;A=210      ;c="At" ;c2="Astatine"     ;return;endif
      if(Z==86 )then;A=222      ;c="Rn" ;c2="Radon"        ;return;endif
      if(Z==87 )then;A=223      ;c="Fr" ;c2="Francium"     ;return;endif
      if(Z==88 )then;A=226.0254 ;c="Ra" ;c2="Radium"       ;return;endif
      if(Z==89 )then;A=227.0278 ;c="Ac" ;c2="Actinium"     ;return;endif
      if(Z==90 )then;A=232.0381 ;c="Th" ;c2="Thorium"      ;return;endif
      if(Z==91 )then;A=231.0359 ;c="Pa" ;c2="Protactinium" ;return;endif
      if(Z==92 )then;A=238.029  ;c="U"  ;c2="Uranium"      ;return;endif
      if(Z==93 )then;A=237.0482 ;c="Np" ;c2="Neptunium"    ;return;endif
      if(Z==94 )then;A=242      ;c="Pu" ;c2="Plutonium"    ;return;endif
      if(Z==95 )then;A=243      ;c="Am" ;c2="Americium"    ;return;endif
      if(Z==96 )then;A=247      ;c="Cm" ;c2="Curium"       ;return;endif
      if(Z==97 )then;A=247      ;c="Bk" ;c2="Berkelium"    ;return;endif
      if(Z==98 )then;A=251      ;c="Cf" ;c2="Californium"  ;return;endif
      if(Z==99 )then;A=252      ;c="Es" ;c2="Einsteinium"  ;return;endif
      if(Z==100)then;A=257      ;c="Fm" ;c2="Fermium"      ;return;endif
      if(Z==101)then;A=258      ;c="Md" ;c2="Mendelevium"  ;return;endif
      if(Z==102)then;A=250      ;c="No" ;c2="Nobelium"     ;return;endif
      if(Z==103)then;A=260      ;c="Lr" ;c2="Lawrencium"   ;return;endif
      if(Z==104)then;A=261      ;c="Rf" ;c2="Rutherfordium";return;endif
      if(Z==105)then;A=262      ;c="Db" ;c2="Dubnium"      ;return;endif
      if(Z==106)then;A=263      ;c="Sg" ;c2="Seaborgium"   ;return;endif
      if(Z==107)then;A=262      ;c="Bh" ;c2="Bohrium"      ;return;endif
      if(Z==108)then;A=255      ;c="Hs" ;c2="Hassium"      ;return;endif
      if(Z==109)then;A=256      ;c="Mt" ;c2="Meitnerium"   ;return;endif
      if(Z==110)then;A=269      ;c="Ds" ;c2="Darmstadtium" ;return;endif
      if(Z==111)then;A=272      ;c="Rg" ;c2="Roentgenium"  ;return;endif
      if(Z==112)then;A=277      ;c="Uub";c2="Ununbiium"    ;return;endif
      end subroutine udm_get_atomic_information

************************************************************************
      integer function udm_thread_id()
************************************************************************
      common /ipomp0/ipomp,npomp
!$OMP THREADPRIVATE(/ipomp0/)
      udm_thread_id=ipomp
      return
      end function

************************************************************************
      function udm_tolower(s)
C     A function that converts uppercase letters in the argument s to lowercase and returns the result.
************************************************************************
      character*(*) s
      character*(len(s)) udm_tolower
      character*26 alph, lower
      data alph  /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data lower /'abcdefghijklmnopqrstuvwxyz'/
      integer i, pos, n
      n = len(s)
      do 10 i = 1, n
         pos = index(alph, s(i:i))
         if (pos .gt. 0) then
            udm_tolower(i:i) = lower(pos:pos)
         else
            udm_tolower(i:i) = s(i:i)
         end if
10    continue
      return
      end

************************************************************************
      end module udm_Utility
************************************************************************



