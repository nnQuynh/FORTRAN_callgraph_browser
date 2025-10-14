!------------------------------egs5init.f-----------------------------
      subroutine egs5init
      use EGS5_MS_MOD !FURUTA20140825
      use EGS5_BCOMP_MOD !OGAWA20140830
      USE egs5_brempr_mod ,only: ibrdst, iprdst !<- 2015.08xx allocatable 1
      USE egs5_eiicom_mod ,only: ieispl, neispl !<- 2015.08xx allocatable 2
      USE egs5_media_mod   !<- 2015.08xx allocatable
      USE egs5_edge_mod    !<- 2015.08xx allocatable
      USE egs5_thresh_mod  !<- 2015.08xx allocatable
      implicit none
      save
      include 'include/egs5_h.f'
      include 'include/egs5_bounds.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_usersc.f'
      include 'include/randomm.f'

      include 'include/egs5_userxt.f'

      real*8 ein,xin,yin,zin,uin,vin,win,wtin,emax
      integer iqin, irin, i, j

*-----------------------------------------------------------------------

! T.Sato 2015/07/25, add ipegs parameter
      integer iegsemi, iegsout, ipegs
      common /egsemi/ iegsemi, iegsout, ipegs
      integer npeme
      common /mpiegs/ npeme
      integer iegsrand
      common /egsrand/ iegsrand

      integer npe, me
      common /mpi00/ npe, me
*-----------------------------------------------------------------------
c
c> add "mxmat" for  egs5 parameter "MXMED" . 2015/08/xx. Local / Common.
      integer  mxmat, mxmat0, mxnel
      common /kmat1a/ mxmat, mxmat0, mxnel
      integer imxmed , flg_local !<- debug1
      character(len=50) cfname   !<- debug2

      include 'param.inc'
      integer idrg,idgr
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      integer ititl,ipara,ibody,iregn,llarr,itby,itar
      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar

c> ada.allocatable
c...
c> ada.allocatable
      integer n, nn, m, mat, kmat, nmat
      integer ,allocatable :: nel(:), zz(:,:), a(:,:)
      real*8 ,allocatable :: den(:),stoich(:,:)
      character*24 ,allocatable :: medarr(:)
      character*24 ,allocatable :: medname(:)
c< ada.allocatable

      real*8  pden, tstoich, rmass
      character*1 matchar
      integer n1,n2,n3,n4,r10,r100,r1000
      character*1 c1,c2,c3,c4
c
c for asymt(z)
      include './pegscommons/elmtbc.f'

c for watbl
      include './pegscommons/elemtb.f'

c for emin
      common /eparm/  esmax, esmin, emin(20)
      real*8 esmax, esmin, emin, emine, eming

c for emax_end
      real*8 emaxend, emaxg    ! emaxe defined egs5_usersc
      common /egs5cmn1/emaxend

c for idmg
      integer idmg,idmn,idnm,imat
      common /regdm/idmg(kvlmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

c a new parameter for EGS5 user discard function
      integer egs5disc
      common /egs5cmn2/egs5disc(kvlmax)

c material name !>ada.allocatable / character*24 medname(MXMED)
      character*1 asymt1,asymt2
      integer k
c
      real*8 thard, tinel, tmscat, hardstep, sig, scpow, dedx, sig0, ams
      common /egs5cmn3/thard,tinel,tmscat,hardstep,sig,scpow,dedx,sig0,
     $                 ams
!$OMP THREADPRIVATE(/egs5cmn3/)

      real*8 k1i,k1r,k1s
      common /egs5cmn4/k1i,k1r,k1s
!$OMP THREADPRIVATE(/egs5cmn4/)

      integer IRAYL,IBOUND,INCOH,ICPROF1,IMPACT,IAPRIM,IUNRSTinp,iepstfl
      real*8 gasdens,GASP

      integer med_p2e(kvlmax)
      common /egs5cmn7/med_p2e

      integer mstz
      real*8 parz
      common /paraj/  mstz(300), parz(300)

! 2015/5/28 T.Sato, Datapath common
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn

*-----------------------------------------------------------------------
c> exhange material size.  from EGS5's parameter "MXMED" to Input "mxmat" of common/kmat1a/
      imxmed = mxmat
      call ALLOCATE_EGS5_MS(NFIT,NFIT1,NEXFIT1,NMSE,NK1,imxmed) ! call ALLOCATE_EGS5_MS(NFIT,NFIT1,NEXFIT1,NMSE,NK1,MXMED)
      call ALLOCATE_EGS5_BCOMP(imxmed, MXSCTFF, MXCP, MXNS)     ! call ALLOCATE_EGS5_BCOMP(MXMED, MXSCTFF, MXCP, MXNS)
c
      call ALLOCATES_EGS5
c>ada.allocatable local
! T.Sato, extend MXEL to MXEL*10 because isotopes are distinguished when PHITS reads material
       allocate( nel(imxmed), zz(imxmed,MXEL*10), a(imxmed,MXEL*10),
     + den(imxmed), stoich(imxmed,MXEL*10),
     + medarr(imxmed), medname(imxmed), stat=flg_local)
       if( flg_local .ne. 0 ) then
         write(*,*)'# egs5init local array can not allocate !! '
         STOP'# egs5init local array can not allocate !! '
       endif
c<ada.allocatable local

      npeme = me

      emine = min(emin(12),emin(13))+RM
      emaxe = emaxend+RM
      eming = emin(14)
      emaxg = emaxend

      if(emine.gt.emaxe)emine=emaxe-1
      if(eming.gt.emaxg)eming=emaxg-1

*-----------------------------------------------------------------------
cc H.Iwase 2014/1/31  (set EGS5 parameters)

         if( mstz(93) .ge. 0 ) then
            inseed = mstz(93)
            iegsrand = 0
         else
            inseed = 0
            iegsrand = 1
         end if
! -------------------------------------------------------------
! EGS5 parameters (should be set by PHITS parameters in future)
! -------------------------------------------------------------
! set pegs parameters

cc$$$cc H.Iwase 2014/4/4 (set pegs5 parameters)
cc not needed; this is already set in pegs5.f

! T.Sato 2014/8/29 Consider correlation of parameters
      gasdens = parz( 182 )
      if(mstz(91).eq.1) mstz(90) = 1  ! incohr should be 1 when iprofr = 1
      if(mstz(90).eq.1) mstz(100) = 1 ! ibound should be 1 when incohr = 1

cc H.Iwase 2015/4/2 set PEGS parameter
cc                  and the setting of EGS parameters moved bottom
cc                  so to be set from PHITS parameters
cc                  after the block_set routine, which inizialize them
         if( mstz(88) .eq. 1 )  IRAYL   = 1
         if( mstz(90) .eq. 1 )  INCOH   = 1
         if( mstz(91) .eq. 1)   ICPROF1 = -3
         if( mstz(92) .eq. 1)   IMPACT  = 1
         if( mstz(100).eq. 1)   IBOUND  = 1
         if( mstz(101).eq. 1)   IAPRIM  = 1
         IUNRSTinp = mstz(107)

         GASP   = 1.0

! not needed
! egs parameters

! -------------------------
! initialize parameters
! -------------------------
      do i = 1, mxmat  !<- correct src.20150808 .
         nel(i) = 0
      enddo

      thard  = 0d0
      tinel  = 0d0
      tmscat = 0d0

      k1i = 0d0
      k1r = 0d0
      k1s = 0d0


c (new) reading material data from PHITS and making indexes of
c converting mat of PHITS to EGS and EGS to PHITS.

c====================================
c phits-ir phits-mat | egs-ir egs-mat
c====================================
c  1          -1     |   1       0
c 10           1     |   2       1
c  3           3     |   3       2
c  2           2     |   4       3
c  5           0     |   5       0

c---------------------------------------------------------

*-----------------------------------------------------------------------
         open(567,file=chfn(23)(1:ilfn(23))//'.tmp',status='unknown')

*-----------------------------------------------------------------------

c egs flow
c (0) med(ir) = j is defined by usercode
c (1) medarr(j) = TA is defined by usercode
c (2) media <= medarr(i)
c (3) HATCH looks TA (medarr(j)) in pegs.inp
c
c memo: egs can treat different density inputs in one material


*-----------------------------------------------------------------------

c reading PHITS-mat
      i = 1
 199  read(567,*,end=200)mat,den(mat),nel(mat),zz(mat,nel(mat)),
     $                      a(mat,nel(mat)), stoich(mat,nel(mat))
c
! T.Sato 2015/09/30, debug
      if( (mat .gt. mxmat))then
       write(*,*)'Error in egs5init, mat > mxmat'
       write(*,*)'mat=',mat,' mxmat=',mxmat
       STOP
      endif
c

      i = i + 1
      go to 199
 200  close(567)

      if( npeme .eq. 0 . and. ipegs .le. 0)
     &   open(566,file=chfn(23)(1:ilfn(23))//'.inp',status='unknown')

*-----------------------------------------------------------------------

      nmed = mat

      do i = 1, nmed         ! creat a medarr
         medarr(i)  = '                        '
         medname(i) = '                        '
      enddo

*-----------------------------------------------------------------------

c PHITS-mat --> EGS-mat (medarr)
      do mat = 1, nmed

cc H.Iwase 2014/8/12 moved from below
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! calc density
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if( den(mat) .lt. 0 ) then
               den(mat) = -den(mat) ! density was given by (g/cm3)

            elseif( den(mat) .eq. 0d0 ) then
               den(mat) = 1d0       ! set dummy density for PEGS

            endif

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! sum up for unique Z for nel(mat)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 811        continue

            do i = 1, nel(mat)-1

               do j = i+1, nel(mat)

                  if( zz(mat,i) .eq. zz(mat,j) ) then ! find same Z
                     stoich(mat,i) = stoich(mat,i) + stoich(mat,j) ! sumup

                     if     ( j .eq. nel(mat) )then   ! in the case j = nel(mat)
                        nel(mat) = nel(mat) - 1
                        goto 811
                     elseif ( j .lt. nel(mat) )then ! in the case j < nel(mat)
                        do k = j, nel(mat)-1
                           stoich(mat,k) = stoich(mat,k+1)
                           zz(mat,k) = zz(mat,k+1)
                            a(mat,k) =  a(mat,k+1)
                        enddo
                        nel(mat) = nel(mat) - 1
                        goto 811
                     else
                        write(*,*)'error in egs5init.f, ask developer'
                        stop
                     endif

                  endif
               enddo

            enddo
! T.Sato 2015/09/30, check element / material
       if(nel(mat).gt.MXEPERMED) then
        write(*,*)'#No. of elements in a material should be '//
     &         'less than MXEPERMED in EGS mode'
        write(*,*)'#In material no.',idmn(mat),
     &       ', no. of elements is',nel(mat)
        stop
       endif

cc H.Iwase 2014/8/12 add ELEM-gas treatment
            if(nel(mat).eq.1) then ! ELEM
! T.Sato 2015/3/2, H, N and O should be composit material irrespective of their density
            if(abs(den(mat)) .le. gasdens .or.  ! gas
     $      zz(mat,1).eq.1.or.zz(mat,1).eq.7.or.zz(mat,1).eq.8)then ! these element should be composit irrespective of their density

            if(      zz(mat,1).ne.2  .and. zz(mat,1).ne.10
     $         .and. zz(mat,1).ne.18 .and. zz(mat,1).ne.36
     $         .and. zz(mat,1).ne.54 .and. zz(mat,1).ne.86 )then  ! except noble gas

               nel(mat)      = 2
               zz(mat,2)     = zz(mat,1)
                a(mat,2)     =  a(mat,1)
               stoich(mat,1) = stoich(mat,1)/2.0
               stoich(mat,2) = stoich(mat,1)

            endif
            endif
         endif
         i=1

cc H.Iwase 2014/8/12 change rule on matname, MAT0001

               do n = 1, nel(mat)

                  asymt1=asymt(zz(mat,n))(1:1)
                  asymt2=asymt(zz(mat,n))(2:2)

                  if( asymt2.ne.' ') then               ! no blank in name like "Fe"
                   if(i.le.23) then  ! T.Sato 2015/10/17, reduced from 24 to 23 on 2024/09/22
                     medname(mat)(i:i)     = asymt1     ! then both characters are
                     medname(mat)(i+1:i+1) = asymt2     ! put in the medname
                   endif
                   i=i+2
                  else                                  ! the name includes blank
                   if(i.le.24) then  ! T.Sato 2015/10/17
                     medname(mat)(i:i) = asymt1         ! like " H" then only the
                   endif
                   i=i+1                              ! latter is used
                  endif


               enddo
               call getmedname(medname(mat),i,abs(den(mat)),gasdens,
     &         nel(mat))
               medarr(mat)  = 'MAT                     '
               i=4

               n4 = mat/1000
               n3 = ( mat-n4*1000 )/100
               n2 = ( mat-n4*1000 -n3*100 )/10
               n1 = mat-n4*1000-n3*100-n2*10

               write(c4(1:1),'(i1)')n4
               write(c3(1:1),'(i1)')n3
               write(c2(1:1),'(i1)')n2
               write(c1(1:1),'(i1)')n1

               medarr(mat)(i:i)     = c4
               medarr(mat)(i+1:i+1) = c3
               medarr(mat)(i+2:i+2) = c2
               medarr(mat)(i+3:i+3) = c1

      end do

*-----------------------------------------------------------------------

      do j=1,nmed    ! initialize media
        do i=1,24
          media(i,j)='    '
        end do
      end do

      do j=1,nmed
        do i=1,24
          media(i,j)(1:1)=medarr(j)(i:i)
        end do
      end do

! =====================================================================
! write pegs5.inp for all mat
! =====================================================================

*-----------------------------------------------------------------------
      if( npeme .eq. 0 .and. ipegs. le. 0) then
*-----------------------------------------------------------------------

      do mat = 1, nmed

         if(mat .gt. 0)then

!  T.Sato 2020/07/02, ICRU90 mode
         if(mstz(108).eq.1.and.(medname(mat)(1:3).eq.'H2O'.or.
     &   medname(mat)(1:2).eq.'C-'.or.
     &   medname(mat)(1:12).eq.'AIR-GAS-ICRU')) then
          iepstfl = 1
         else
          iepstfl = 0
         endif

c--------------
c     ELEM
c--------------


            if ( nel(mat) .le. 1) then ! 1 material --> ELEM

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! write pegs5.inp
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               write(566,'(a4)')'ELEM'

               if(abs(den(mat)) .gt. gasdens)then ! non-gas

c T.Sato 2014/8/29 Output all parameters, 2015/7/25, add IUNRSTinp,iepstfl
                write(566,302)den(mat),IRAYL,IBOUND,INCOH,ICPROF1,
     $          IMPACT,IAPRIM,IUNRSTinp,iepstfl
               else             ! gas
                  write(566,301)IRAYL,IBOUND,INCOH,ICPROF1,IMPACT,
     $                 IAPRIM,IUNRSTinp,iepstfl,den(mat),GASP

               endif
               write(566,'(a24,6x,a24)')medarr(mat),medname(mat)

               write(566,'(a2)')asymt(zz(mat,1))


 3             write(566,'(a4)')'ENER'
               write(566,500)' &INP AE=',emine,',AP=',eming,
     $                        ',UE=',emaxe,',UP=',emaxg,
     $                   ' &END'

               write(566,'(a4)')'TEST'
               write(566,*)'&INP  &END'
               write(566,'(a4)')'PWLF'
               write(566,*)'&INP  &END'
               write(566,'(a4)')'DECK'
               write(566,*)'&INP  &END'

 301           format('&INP ','IRAYL=',i1,',IBOUND=',i1,',INCOH=',i1,
     $      ',ICPROF=',i3,',IMPACT=',i1,',IAPRIM=',i1,',IUNRST=',i1,
     $      ',EPSTFL=',i1,',RHO=',1pe12.6,',GASP=',1pe12.6,' &END')
 302           format('&INP ','RHO=',1pe12.6,',IRAYL=',i1,',IBOUND=',i1,
     $         ',INCOH=',i1,',ICPROF=',i3,',IMPACT=',i1,',IAPRIM=',i1,
     $      ',IUNRST=',i1,',EPSTFL=',i1,' &END')

c--------------
c or  COMP
c--------------
            elseif ( nel(mat) .ge. 2 ) then ! more than 2 materials --> COMP

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! write pegs5.inp
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

               if(abs(den(mat)) .gt. gasdens)then ! non-GAS

                  write(566,'(a4)')'COMP'

                  write(566,480)nel(mat),den(mat)
 480              format('&INP NE=',i2,',RHO=',1pe12.6,',PZ=',$)

                  do n = 1, nel(mat)
                     write(566,481)stoich(mat,n)
                     if(n.lt.nel(mat)) write(566,482)', '
                  enddo


c T.Sato 2014/8/29 Output all parameters, 2015/7/25, add IUNRSTinp,iepstfl
                  write(566,483)IRAYL,IBOUND,INCOH,ICPROF1,IMPACT,
     $          IAPRIM,IUNRSTinp,iepstfl
 483           format(',IRAYL=',i1,',IBOUND=',i1,',INCOH=',i1,
     $         ',ICPROF=',i3,',IMPACT=',i1,',IAPRIM=',i1,
     $         ',IUNRST=',i1,',EPSTFL=',i1,' &END')

                  write(566,'(a24,6x,a24)')medarr(mat),medname(mat)

                  do n = 1, nel(mat)
                     write(566,'(a3$)')asymt(zz(mat,n))
                  enddo


                  write(566,*)''
                  write(566,'(a4)')'ENER'

                  write(566,500)' &INP AE=',emine,',AP=',eming,
     $                 ',UE=',emaxe,',UP=',emaxg,
     $                 ' &END'

                  write(566,'(a4)')'TEST'
                  write(566,*)'&INP  &END'
                  write(566,'(a4)')'PWLF'
                  write(566,*)'&INP  &END'
                  write(566,'(a4)')'DECK'
                  write(566,*)'&INP  &END'

               else             ! GAS

                  write(566,'(a4)')'COMP'

                  write(566,484)nel(mat),den(mat),GASP
 484              format('&INP NE=',i2,',RHO=',1pe12.6,',GASP=',1pe12.6,
     $                 ',PZ=',$)

                  do n = 1, nel(mat)
                     write(566,481)stoich(mat,n)
                     if(n.lt.nel(mat)) write(566,482)', '
                  enddo


c T. Sato 2014/8/29 Output all parameters, 2015/7/25, add IUNRSTinp,iepstfl
                  write(566,483)IRAYL,IBOUND,INCOH,ICPROF1,IMPACT,
     $          IAPRIM,IUNRSTinp,iepstfl

                  write(566,'(a24,6x,a24)')medarr(mat),medname(mat)


                  do n = 1, nel(mat)
                     write(566,'(a3$)')asymt(zz(mat,n))
                  enddo


                  write(566,*)''
                  write(566,'(a4)')'ENER'
                  write(566,500)' &INP AE=',emine,',AP=',eming,
     $                 ',UE=',emaxe,',UP=',emaxg,
     $                 ' &END'
                  write(566,'(a4)')'TEST'
                  write(566,*)'&INP  &END'
                  write(566,'(a4)')'PWLF'
                  write(566,*)'&INP  &END'
                  write(566,'(a4)')'DECK'
                  write(566,*)'&INP  &END'

               endif

            endif

         endif

      enddo

      close(566)

*-----------------------------------------------------------------------
      end if
*-----------------------------------------------------------------------


 481  format(g13.5$)
 482  format(a2$)
 500  format(a9,es12.5,a4,es12.5,a4,es12.5,a4,es12.5,a5)


*-----------------------------------------------------------------------

      if( iegsout .eq. 2 .and. npeme .eq. 0 ) then
         open(UNIT= 506,FILE=chfn(23)(1:ilfn(23))//'-egs5.out',
     &   STATUS='unknown')
      else
         open(UNIT= 506,form='formatted',STATUS='scratch')
      end if

      call counters_out(0)
      call block_set                 ! Initialize some general variables

      nreg = iregn

      do i = 1, nmed !<- correct 20150808, normalized to density T.Sato 2016/03/14
       if(parz(185).gt.0.0) then
         chard(i) = parz(185)
       else
         chard(i) = abs(parz(185))/den(i)
       endif
      enddo

*-----------------------------------------------------------------------

      if( npeme .eq. 0 .and. ipegs. le. 1) call pegs5

*-----------------------------------------------------------------------

      if(ipegs.lt.0) then
       write(*,*) 'Only PEGS is executed when ipegs = -1'
       write(*,*) 'PEGS is finished successfully'
       stop
      endif

      j = 0
      do i = 1, iregn


         ecut(i) = min(emin(12), emin(13))
         ecut(i) = ecut(i) + RM
         pcut(i) = emin(14)

*-----------------------------------------------------------------------


! procedure in egs5ede.f
! (1) if e < ecut, then set e = ecut with deresid ne 0
! (2) the next set deresid = 0 and e = emin

c iwase, convert PHITS-mat to EGS-mat
         if( idmg(i).gt.0 ) then
           imat = idnm(idmg(i))
         else
           imat = 0
         end if
          if( imat .gt. 0 ) then
             j = j + 1
c iwase check it
             med(i) = imat          ! med = 1,2,3,...
             egs5disc(i) = 0        ! idisc = 0 for materials

          elseif ( imat .eq. 0 ) then
             med(i) = 0             ! med = 0
             egs5disc(i) = 0        ! idisc = 0 for inner void

          elseif ( imat .eq. -1 )then
             med(i) = 0             ! med = 0
             egs5disc(i) = 1        ! idisc = 1 for outer void

          endif

      enddo

*-----------------------------------------------------------------------

120   FORMAT(/,' inseed=',I12,5X,
     *         ' (seed for generating unique sequences of Ranlux)')


      call rluxinit  ! Initialize the Ranlux random-number generator

*-----------------------------------------------------------------------

      emaxe = emaxend

      write(506,130)

      open(UNIT=KMPI,FILE=chfn(23)(1:ilfn(23))//'.dat',STATUS='unknown')

      if( iegsout .eq. 2 .and. npeme .eq. 0 ) then
         open(UNIT=KMPO,FILE=chfn(23)(1:ilfn(23))//'-egs5.dummy',
     &   STATUS='unknown')
      else
         open(UNIT=KMPO,form='formatted',status='scratch')
      end if

      write(506,140)

130   format(/' Start tutor1'/' Call hatch to get cross-section data')
140   FORMAT(/,' HATCH-call comes next',/)

*-----------------------------------------------------------------------
cc H.Iwase 2015/4/2 set EGS parameters (moved here)

      do  i = 1, iregn

         iedgfl(i) = mstz(86) ! EGS parameter
         iauger(i) = mstz(87) ! EGS parameter
         iraylr(i) = mstz(88) ! EGS parameter
         lpolar(i) = mstz(89) ! EGS parameter
         incohr(i) = mstz(90) ! EGS parameter
         impacr(i) = mstz(92) ! EGS parameter
         iphter(i) = mstz(99) ! EGS parameter
         iprofr(i) = mstz(91) ! EGS parameter

      enddo

         ieispl = mstz(95) ! EGS parameter
         neispl = mstz(96) ! EGS parameter
         ibrdst = mstz(97) ! EGS parameter
         iprdst = mstz(98) ! EGS parameter

*-----------------------------------------------------------------------

      call hatch

*-----------------------------------------------------------------------

         close(UNIT=KMPI)
         close(UNIT=KMPO)

*-----------------------------------------------------------------------

      write(506,150) ae(1)-RM, ap(1)
      write(506,160)

150   format(/' Knock-on electrons can be created and any electron ',
     *'followed down to' /T40,F8.3,' MeV kinetic energy'/
     *' Brem photons can be created and any photon followed down to',
     */T40,F8.3,' MeV')
160   format(/T19,'Kinetic energy(MeV)',T40,'charge',T48,
     *'angle w.r.t. z axis-degrees')

*-----------------------------------------------------------------------
       deallocate( nel, zz, a, den, stoich, medarr, medname)
      return
      end
!------------------------------ end--------------------------
!======= ADD ALLOCATE EGS5 COMMONS =====================================!
! allocatable array of EGS5
      SUBROUTINE ALLOCATES_EGS5
      USE egs5_brempr_mod !<- 2015.08xx allocatable 1
      USE egs5_eiicom_mod !<- 2015.08xx allocatable 2
      USE egs5_elecin_mod !<- 2015.08xx allocatable 3
      USE pegs_dcsstr_mod !<- 2015.08xx allocatable 4
      USE egs5_media_mod  !<- 2015.08xx allocatable 5
      USE egs5_edge_mod     !<- 2015.08xx allocatable
          USE egs5_photin_mod   !<- 2015.08xx allocatable
          USE egs5_mscon_mod    !<- 2015.08xx allocatable
          USE egs5_scpw_mod     !<- 2015.08xx allocatable
          USE egs5_thresh_mod   !<- 2015.08xx allocatable
      implicit none
      include 'include/egs5_h.f'
!
!> add "mxmat" for  egs5 parameter "MXMED" . 2015/08/xx. Local / Common.
!
      integer  mxmat, mxmat0, mxnel
      common /kmat1a/ mxmat, mxmat0, mxnel
!> work and debug
      integer imxmed !<- work of MXMED
!.....
      imxmed = mxmat ! material size from PHITS input.  !     imxmed = 109 !  lec01 dummy
!
      call ALLOCATE_EGS5_BREMPR(imxmed, MXEL, MXPWR2I)
      call ALLOCATE_EGS5_EIICOM(imxmed, MXEKE,MXEPERMED)
      call ALLOCATE_EGS5_ELECIN
     + ( imxmed, MXEKE, MSSTEPS, MXVRT1, MXVRT2, MXJREFF)
      call  ALLOCATE_PEGS_DCSSTR
     + ( imxmed, NEGRID, NREDA)
      call  ALLOCATE_EGS5_MEDIA( imxmed)
!....
      call  ALLOCATE_EGS5_EDGE( imxmed, MXEPERMED, MXREG)
      call ALLOCATE_EGS5_PHOTIN( imxmed, MXGE, MXRAYFF)
      call ALLOCATE_EGS5_MSCON
     + ( imxmed, NESCPW, NMSE, NK1, NFIT, NFIT1, NEXFIT)
      call ALLOCATE_EGS5_SCPW(imxmed, MXEKE)
      call ALLOCATE_EGS5_THRESH(imxmed)
!
      END SUBROUTINE ALLOCATES_EGS5
!=======================================================================!
      SUBROUTINE DEALLOCATES_EGS5
      use egs5_ms_mod
      use egs5_bcomp_mod
      USE egs5_brempr_mod !<- 2015.08xx allocatable 1
      USE egs5_eiicom_mod !<- 2015.08xx allocatable 2
      USE egs5_elecin_mod !<- 2015.08xx allocatable 3
      USE pegs_dcsstr_mod !<- 2015.08xx allocatable 4
      USE egs5_media_mod  !<- 2015.08xx allocatable 5
      USE egs5_edge_mod       !<- 2015.08xx allocatable
          USE egs5_photin_mod     !<- 2015.08xx allocatable
          USE egs5_mscon_mod      !<- 2015.08xx allocatable
          USE egs5_scpw_mod       !<- 2015.08xx allocatable
          USE egs5_thresh_mod     !<- 2015.08xx allocatable
         call DEALLOCATE_EGS5_MS
         call DEALLOCATE_EGS5_BCOMP
         call DEALLOCATE_EGS5_BREMPR
         call DEALLOCATE_EGS5_EIICOM
c
         call DEALLOCATE_EGS5_ELECIN
         call DEALLOCATE_PEGS_DCSSTR  !<- use only in egs5init
         call DEALLOCATE_EGS5_MEDIA
c..
         call DEALLOCATE_EGS5_EDGE  !<- NG over Material.
c..
         call DEALLOCATE_EGS5_PHOTIN
         call DEALLOCATE_EGS5_MSCON
         call DEALLOCATE_EGS5_SCPW
         call DEALLOCATE_EGS5_THRESH
      END SUBROUTINE DEALLOCATES_EGS5
!
!======= END ALLOCATE EGS5 COMMONS =====================================!

! Change proper material name which is listed in Table 2.3 of SLAC-R-730/KEK-2005-8
! for Density correction
      subroutine getmedname(medname,i,density,gasdens,nel)
      implicit none
      character medname*24
      integer i,k,nel
      real*8 density,gasdens
      k=min(24,i-1) ! real character length, should be less tha 24, T.Sato 2024/09/22
      if(density.lt.gasdens) then ! Gas
       if(medname(1:k).eq.'H' ) medname='H2-GAS                  '
       if(medname(1:k).eq.'HE') medname='HE-GAS                  '
       if(medname(1:k).eq.'N' ) medname='N2-GAS                  '
       if(medname(1:k).eq.'O' ) medname='O2-GAS                  '
       if(medname(1:k).eq.'NE') medname='NE-GAS                  '
       if(medname(1:k).eq.'AR') medname='AR-GAS                  '
       if(medname(1:k).eq.'XE') medname='XE-GAS                  '
       if(medname(1:k).eq.'RN') medname='RN-GAS                  '
       if(medname(1:k).eq.'NO') medname='AIR-GAS                 '
       if(medname(1:k).eq.'ON') medname='AIR-GAS                 '
       if(medname(1:k).eq.'NOAR') medname='AIR-GAS                 '
       if(medname(1:k).eq.'ONAR') medname='AIR-GAS                 '
       if(medname(1:k).eq.'CO') medname='CO2-GAS                 '
       if(medname(1:k).eq.'OC') medname='CO2-GAS                 '
       if(medname(1:k).eq.'CNOAR') medname='AIR-GAS-ICRU            '
      else ! Non-Gas
       if(medname(1:k).eq.'H' ) medname='H2-LIQUID               '
       if(density.ge.1.55.and.density.lt.1.85.and.  ! T.Sato 2020/07/02, accept for wider density range
     1    medname(1:k).eq.'C' ) medname='C-1.70G/CM**2           '
       if(density.ge.1.85.and.density.lt.2.14.and.  ! T.Sato 2020/07/02, accept for wider density range
     1    medname(1:k).eq.'C' ) medname='C-2.00G/CM**2           '
       if(density.ge.2.14.and.density.lt.2.40.and.  ! T.Sato 2020/07/02, accept for wider density range
     1    medname(1:k).eq.'C' ) medname='C-2.265G/CM**2          '
       if(density.gt.0.9.and.density.lt.1.1.and.
     1    medname(1:k).eq.'HO') medname='H2O                     '
       if(density.gt.0.9.and.density.lt.1.1.and.
     1    medname(1:k).eq.'OH') medname='H2O                     '
       if(medname(1:k).eq.'INA') medname='NAI                     '
       if(medname(1:k).eq.'ICS') medname='CSI                     '
       if(medname(1:k).eq.'SIO') medname='SIO2                    '
       if(medname(1:k).eq.'OSI') medname='SIO2                    '
      endif
      if(nel.eq.2) then ! avoid mismatch by combination of two element
       if(medname(1:k).eq.'SI' ) medname='SI-comp                 '
       if(medname(1:k).eq.'SC' ) medname='SC-comp                 '
       if(medname(1:k).eq.'CO' ) medname='CO-comp                 '
       if(medname(1:k).eq.'NI' ) medname='NI-comp                 '
       if(medname(1:k).eq.'NB' ) medname='NB-comp                 '
       if(medname(1:k).eq.'IN' ) medname='IN-comp                 '
       if(medname(1:k).eq.'SN' ) medname='SN-comp                 '
       if(medname(1:k).eq.'SB' ) medname='SB-comp                 '
       if(medname(1:k).eq.'CS' ) medname='CS-comp                 '
       if(medname(1:k).eq.'HF' ) medname='HF-comp                 '
       if(medname(1:k).eq.'HO' ) medname='HO-comp                 '
       if(medname(1:k).eq.'YB' ) medname='YB-comp                 '
       if(medname(1:k).eq.'OS' ) medname='OS-comp                 '
       if(medname(1:k).eq.'PB' ) medname='PB-comp                 '
       if(medname(1:k).eq.'BI' ) medname='BI-comp                 '
       if(medname(1:k).eq.'PO' ) medname='PO-comp                 '
       if(medname(1:k).eq.'NP' ) medname='NP-comp                 '
       if(medname(1:k).eq.'PU' ) medname='PU-comp                 '
       if(medname(1:k).eq.'CF' ) medname='CF-comp                 '
       if(medname(1:k).eq.'NO' ) medname='NO-comp                 '
      endif
      return
      end
