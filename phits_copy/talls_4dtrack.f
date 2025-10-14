************************************************************************
      module t4dtrack_mod
************************************************************************
      implicit real*8(a-h,o-z)
      include 'param.inc'
C     ==================================================================
      integer, parameter :: nBufNucMax = 100 ! Maximum number of unstable nuclei buffered per source.
      integer, parameter :: nEventsMax = 10000 ! Maximum number of accepted events.
C     ==================================================================
C     Universal variables
      integer, save :: it4dtrack = 0
      integer, save :: nSection = 0
      integer, save :: isec
!$OMP THREADPRIVATE(isec)
C     ==================================================================
C     Variables for each section
      integer, save :: no_pre(itlmax)
      double precision, save :: x_start(itlmax)
      double precision, save :: y_start(itlmax)
      double precision, save :: z_start(itlmax)
!$OMP THREADPRIVATE(no_pre,x_start,y_start,z_start)
      logical, save :: buff_exist(itlmax) ! buffered start point exist or not
      character*30, save :: StartPoint(itlmax,20) ! Buffered start point
      integer, save :: n_line_of_this_source(itlmax)
      integer, save :: n_line_of_this_track(itlmax)
      logical, save :: lastLineWasEndTrack(itlmax)
!$OMP THREADPRIVATE(n_line_of_this_source,n_line_of_this_track)
!$OMP THREADPRIVATE(buff_exist,StartPoint,lastLineWasEndTrack)
      integer, save :: n_line_of_this_file(itlmax)
C     ----------
      character*200, save :: filename(itlmax)
C     ----------
      double precision, save :: stepminEM(itlmax)
C     ----------
C     AcceptedHistories
      integer, save :: n_AcceptedHistories(itlmax)
      integer, allocatable,save:: AcceptedHistories(:,:,:) ! accepted (nocas , nobch)
      logical, save :: is_RejectedEvent(itlmax)
      integer, save :: historyMax(itlmax)
      integer, save :: n_WrittenHistories(itlmax)
!$OMP THREADPRIVATE(is_RejectedEvent)
C     ----------
C     PinnedOption
      integer, save :: PinnedOption(itlmax)
      integer, save :: PinnedTimeOption(itlmax)
      integer, allocatable,save:: PinnedParticles(:,:)
      double precision, save :: PinnedTimeFix(itlmax) ! ns
      double precision, save :: HalflifeMinSEC(itlmax) ! second
      double precision, save :: HalflifeMaxSEC(itlmax) ! second
      double precision, save :: aHalf(itlmax)
      double precision, save :: bHalf(itlmax)
      integer, save :: nBufNuc(itlmax)
      character*30, save :: BufNuc(itlmax,nBufNucMax,20) ! buffer_nuclei
!$OMP THREADPRIVATE(nBufNuc,BufNuc)
C     ----------
C     checkCellMode
      logical, save :: checkCellMode(itlmax)
      logical, allocatable,save:: acceptedCells(:,:)
      logical, save :: TrackStart_No_Yes(itlmax)
!$OMP THREADPRIVATE(TrackStart_No_Yes)
C     ----------
C     checkParticleMode
      integer, save :: checkParticleMode(itlmax)
      integer, allocatable,save:: acceptedParticles(:,:)
      integer, allocatable,save:: rejectedParticles(:,:)
      logical, save :: NucleusIsAccepted(itlmax)
C     ----------
      character*10, save :: format(itlmax)
      integer, save :: verbose(itlmax)
      character*10, save :: fmtW(itlmax)
C     ----------
      double precision, save :: SkipRate(itlmax)
      integer, parameter :: n_kfSkipMax = 20 ! Maximum number of kf-code for SkipRate
      integer, save ::      n_kfSkip(itlmax)
      integer, save ::        kfSkip(itlmax,n_kfSkipMax)
      double precision, save :: SkipRate_of(itlmax,n_kfSkipMax)
      integer :: seed_T4D = 123456789
C     ----------
      logical, save :: comment_for_track_structure = .true.
!$OMP THREADPRIVATE(seed_T4D)
C     ==================================================================
C     Temporary variables
      integer, save :: iSave ! Used to buffer unstable nuclei. Like a global variable.
      integer, save :: kfSave
!$OMP THREADPRIVATE(iSave,kfSave)
      integer, save :: nt4domp = -1
C     ==================================================================
      contains
************************************************************************


************************************************************************
      subroutine write_header(i)
************************************************************************
      if(trim(format(isec))=="t4d") then
C       version info
        write(i,'(A)') "#$ t4dv1"

C       history info
        write(i,'(A)',advance='no') "# [history info] h: history-ID"
        if(verbose(isec) >= 1) then
        write(i,'(A)',advance='no') " historyID-in-batch batch-ID"
        endif
        write(i,*)

C       track info
        write(i,'(A)',advance='no') "# [  track info] t: kf-code"
        if(verbose(isec) >= 1) then
        write(i,'(A)',advance='no') " name"
        endif
        if(verbose(isec) >= 2) then
        write(i,'(A)',advance='no') " track-ID"
        endif
        write(i,*)

C       point info
        write(i,'(A)',advance='no') "# [  point info]"
        write(i,'(A)',advance='no') " x[cm] y[cm] z[cm] time[ns]"
        write(i,'(A)',advance='no') " kinetic-energy[MeV] weight"
        if(verbose(isec) >= 1) then
        write(i,'(A)',advance='no') " cell cell_before"
        endif
        if(verbose(isec) >= 2) then
        write(i,'(A)',advance='no') " ncol jcoll kcoll nclsts"
        write(i,'(A)',advance='no') " counter1 counter2 counter3"
        endif
        write(i,*)
        write(i,'(A)') "# ----------------------------------------"

      elseif(trim(format(isec))=="pict") then
        write(i,'(A)') "GSTA-FREE-TIME"
        write(i,'(A)') "MEND"
      endif
      end subroutine


************************************************************************
      integer function numberOfElementsIn(str)
************************************************************************
      character*(*) str
      logical tmp
      numberOfElementsIn=0
      tmp = .false.
      do i=1,len_trim(str)
         if(str(i:i) /= ' ' .and. .not.tmp) then
            numberOfElementsIn=numberOfElementsIn+1 ! Count up when a non-space character is found.
            tmp = .true.
         elseif (str(i:i) == ' ') then
            tmp = .false.
         endif
      end do
      return
      end function


************************************************************************
      subroutine init_t4dtrack()
************************************************************************
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      common /tscmsg/ ktsc(kvlmax), mntsc, ntsc(kvlmax)
      integer i

      do i=1,itlmax
C      ----------
C      Variables for each section
       n_line_of_this_source(i) = 0
       n_line_of_this_track(i) = 0
       n_line_of_this_file(i) = 0
       lastLineWasEndTrack(i) = .true.
C      ----------
       is_RejectedEvent(i) = .false.
C      ----------
       TrackStart_No_Yes(i) = .false.
C      ----------
       ! when track structure is on
       if(mntsc.ne.0) then
        if(historyMax(i)==-1) historyMax(i) = 1
        if(SkipRate(i) < 0d0) SkipRate(i)=0.99d0

        if(comment_for_track_structure) then
         comment_for_track_structure=.false.
         print*,
     &   "When [track structure] is used, the following [t-4Dtrack] par"
     &   //"ameters are automatically set to prevent the output file fr"
     &   //"om becoming excessively large:"
         print*,"  [t-4Dtrack]"
         print*,"  HistoryMax = 1   $ Output tracks for only 1 history"
         print*,"  SkipRate = 0.99  $ Skip 99% of track step output"
         print*,"You can manually change these values if needed."
         print*
        endif

       endif
      enddo

      return
      end subroutine init_t4dtrack

************************************************************************
      subroutine read_t4dtrack(jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
************************************************************************

      use moddas
      use moddas_region

C       use usrtalmod
      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err
      common /subtra/ isubt, ipsub(mxpart)  ! kitamura22/03/31

*-----------------------------------------------------------------------

      integer, parameter :: nparam = 19
      character schan(nparam)*20
      dimension lschn(nparam)
      data ( schan(i) , lschn(i), i = 1, nparam ) /
     &      'file'            ,  4, ! ipm==1
     &      'stepminem'       ,  9, ! ipm==2
     &      'historyfile'     , 11, ! ipm==3
     &      'pinnedoption'    , 12, ! ipm==4
     &      'pinnedparticles' , 15, ! ipm==5
     &      'pinnedtimeoption', 16, ! ipm==6
     &      'pinnedtimefix'   , 13, ! ipm==7
     &      'halflifemin'     , 11, ! ipm==8
     &      'halflifemax'     , 11, ! ipm==9
     &      'ahalf'           ,  5, ! ipm==10
     &      'bhalf'           ,  5, ! ipm==11
     &      'reg'             ,  3, ! ipm==12
     &      'part'            ,  4, ! ipm==13
     &      'historymax'      , 10, ! ipm==14
     &      'format'          ,  6, ! ipm==15
     &      'verbose'         ,  7, ! ipm==16
     &      'precision'       ,  9, ! ipm==17
     &      'skiprate-of'     , 11, ! ipm==18
     &      'skiprate'        ,  8/ ! ipm==19
C     * 'schan' should be written in lowercase letters.
C     * If skiprate-of is placed before skiprate, it cannot be read. When
C       using multiple strings with a common prefix, the longer one must
C       be written first.
*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

C     character dkam*9

      character tname*12
      data      tname /'[t-4Dtrack]'/

      character Param*500
      character ParamLow*500
      character charTMP*100
      integer iTMP
      integer iTmpAcc(mxpart)
      integer iTmpRej(mxpart)

*-----------------------------------------------------------------------

      ierr  = 0

*-----------------------------------------------------------------------

      if(npe.gt.2)goto 987
      call chkt4domp()
C     Initilize [t-4Dtrack]
      if(it4dtrack==0) then
        it4dtrack = 1
        do i=1,itlmax
C         ----------
          stepminEM(i) = -1.0d0 ! cm
C         ----------
          n_AcceptedHistories(i) = 0
          historyMax(i) = -1
          n_WrittenHistories(i) = 0
C         ----------
          PinnedOption(i) = 0
          PinnedTimeOption(i) = 0
          PinnedTimeFix(i) = 10.0d0 ! ns
          HalflifeMinSEC(i) = 1e-22
          HalflifeMaxSEC(i) = 3e+32 ! ~1e+25 year
          aHalf(i) = 1d0
          bHalf(i) = 0d0
C         ----------
          checkCellMode(i)=.false.
C         ----------
          checkParticleMode(i) = 0 ! .false.
          NucleusIsAccepted(i) = .false.
C         ----------
          format(i)="t4d"
          verbose(i)=1
          fmtW(i)="(es13.5)" ! <=> precsion = 6
C         ----------
          SkipRate(i) = -1d0
          n_kfSkip(i) = 0
C         ----------
        enddo
      endif

*-----------------------------------------------------------------------

      nSection = nSection + 1
      isec = nSection

      if(nSection*nt4domp.gt.200)then
       goto 988
      endif

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

      call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
      !     12345678901234567890123456789012
      ! ex:"  file = track.out   $ file name"
      ! i1=3 : Column number (index) of the beginning of non-space text.
      ! i2=32: Column number (index) of the end of non-space text.
      ! i3=18: Column number (index) of the end of non-space text, excluding commented-out sections.
      ! i4=16: ?
      ! chlw(i1:i3): Non-space text excluding the commented-out parts.

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) goto 800

      if( iskip .ne. 0 ) goto 140

  150 continue

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) goto 800

*-----------------------------------------------------------------------
*     end of the section
*-----------------------------------------------------------------------
      if( chlw(i1:i1) .eq. '[' ) then

         jpn = 1

C        ----------------------------------------------

         goto 800 ! retuen

      end if

*-----------------------------------------------------------------------
*     identify the parameters
*-----------------------------------------------------------------------
  200 continue

      do i = 1, nparam
         iTextBegin = i1
         iTextEnd   = i1+lschn(i)-1
         if( chlw(iTextBegin:iTextEnd)==trim(schan(i)) ) then
            ipm = i
            goto 100
         endif
      end do

      goto 998 ! Unknown parameter is found

*-----------------------------------------------------------------------
*        read value of parameters
*-----------------------------------------------------------------------
  100 continue

C     -------------------
      !     12345678901234567890123456789012
      ! ex:"  file = track.out   $ file name"
      ic = inumc(chlw,i1,i3,'=') + 1 ! inumc(chin,i1,i2,cc): get column number of character cc
      ! -> ic=9
      ic = jnumc(chlw,ic,i3) ! jnumc(chin,i1,i2): get column number of non space character
      ! -> ic=10
      iParamBegin = ic
      iParamEnd = i3
      if( iParamBegin > iParamEnd ) goto 997 ! Description of parameter is wrong
      icl = inumc(chlw,ic,i3,';') - 1

C     -------------------
C     Get Parameter from input file.
      Param   =chin(iParamBegin:iParamEnd)
      ParamLow=chlw(iParamBegin:iParamEnd) ! lowercase

C     Count the number of elements (nElem) in Param.
      nElem=numberOfElementsIn(Param)

*-----------------------------------------------------------------------
      if( ipm==1 ) then ! file
         filename(isec)=trim(Param)
*-----------------------------------------------------------------------
      else if( ipm==2 ) then
         read(Param,*,iostat=ios) stepminEM(isec)
         if(ios/=0) goto 997
*-----------------------------------------------------------------------
      else if( ipm==3 ) then ! HistoryFile (AcceptedHistories)
C        -------------------
C        Count the number of valid lines.
         n_AcceptedHistories(isec) = 0
         ! 2000 is for temporary
         open(2000,file=trim(Param),status="old",err=996)
         do
            read(2000, '(A)', iostat=ios) charTMP
            if(ios /= 0) exit
            read(charTMP, *, iostat=ios) dtmp
            if(ios == 0) n_AcceptedHistories(isec) =
     &                   n_AcceptedHistories(isec) + 1
         enddo
C        -------------------
         if(n_AcceptedHistories(isec)==0) then
           goto 989
         endif
C        -------------------
C        Store data in an array.
         rewind(2000)
         if(.not. allocated(AcceptedHistories))
     &     allocate( AcceptedHistories(itlmax,nEventsMax,2) )
         i=0
         do
            read(2000, '(A)', iostat=ios) charTMP
            if(ios /= 0) exit
C           ----------------
C           if(numberOfElementsIn(charTMP)==2) then
            if(.true.) then
               read(charTMP, *, iostat=ios) dNOCAS, dNOBCH
               if(ios == 0) then
                  i=i+1
                  if(i>nEventsMax) then
            print*,"Caution: [t-4Dtrack]. The number of events is up to"
     &      , nEventsMax, ". Please increase nEventsMax."
                     exit
                  endif
                  AcceptedHistories(isec,i,1)=int(dNOCAS)
                  AcceptedHistories(isec,i,2)=int(dNOBCH)
               endif
            endif
C           ----------------
            if(i==n_AcceptedHistories(isec)) exit
         enddo
C        -------------------
         close(2000)
*-----------------------------------------------------------------------
      else if( ipm==4 ) then
         read(Param,*,iostat=ios) PinnedOption(isec)
         if(ios/=0) goto 997
*-----------------------------------------------------------------------
      else if( ipm==5 ) then
         read(Param,*,iostat=ios) PinnedParticles
         if(ios/=0) goto 997
*-----------------------------------------------------------------------
      else if( ipm==6 ) then
         read(Param,*,iostat=ios) PinnedTimeOption(isec)
         if(ios/=0) goto 997
*-----------------------------------------------------------------------
      else if( ipm==7 ) then
         read(Param,*,iostat=ios) PinnedTimeFix(isec)
         if(ios/=0) goto 997
*-----------------------------------------------------------------------
      else if( ipm==8 .or. ipm==9 ) then
         read(Param,*,iostat=ios) tmp, charTMP
         if(ios/=0) goto 997
         if     (trim(charTMP)=="s"  ) then; tmp2=1d0;
         else if(trim(charTMP)=="min") then; tmp2=60d0;
         else if(trim(charTMP)=="h"  ) then; tmp2=3600d0;
         else if(trim(charTMP)=="d"  ) then; tmp2=3600d0*24d0;
         else if(trim(charTMP)=="y"  ) then; tmp2=3600d0*24d0*365d0;
         else                              ; goto 994;
         endif
         if(ipm==8) HalflifeMinSEC(isec) = tmp*tmp2
         if(ipm==9) HalflifeMaxSEC(isec) = tmp*tmp2
*-----------------------------------------------------------------------
      else if( ipm==10 ) then
         read(Param,*,iostat=ios) aHalf(isec)
         if(ios/=0) goto 997
*-----------------------------------------------------------------------
      else if( ipm==11 ) then
         read(Param,*,iostat=ios) bHalf(isec)
         if(ios/=0) goto 997
*-----------------------------------------------------------------------
      else if( ipm==12 ) then ! reg
C        ----------------------------
C        read 'reg'
         call moddas_allocate_int(MAX_NUM_ITREG, idas_itreg_temporary)

         ndsm = 1
         call tregion(1,jsn,jsi,dsin,idsi,ill,ilf,
     &                jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                ntrn,mtrn,ndsm,nvol,ivl,irvl,1
     &                ,MAX_NUM_ITREG,idas_itreg_temporary)

         if( mtrn > MAX_NUM_ITREG ) then
            write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &           'sub.tprodct@tallsm2.f'
     &              //' ?dimension over idas_itreg_temporary?'
     &              //' mtrn > MAX_NUM_ITREG'
     &           ,' (mtrn=',mtrn,')'
     &           ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
            ErrID = 'L:522/R:read_t4dtrack/F:talls_4dtrack.f'
            call ErrWrite(ErrID,ErrCha)
         endif

         if( ierr .ne. 0 ) return
         if( jpn  .eq. 3 ) goto 800
         if( ntrn .lt. -1 ) goto 995

C        ----------------------------
C        acceptedCells
         if(.not. allocated(acceptedCells)) then
            allocate( acceptedCells(itlmax,kvmmax) )
            do i=1,itlmax
            do j=1,kvmmax
               acceptedCells(i,j)=.false.
            enddo
            enddo
         endif

         do i=1,ntrn ! ntrn = Number of cell in 'reg'
            icell=idas_itreg_temporary(i)
            if(icell==5000000) then ! reg = all
               checkCellMode(isec)=.false.
               exit
            endif
            acceptedCells(isec,icell)=.true.
            checkCellMode(isec)=.true.
         enddo

         goto 150

*-----------------------------------------------------------------------
      else if( ipm==13 ) then ! part

         inpat = 0
         nAcc = 0
         nRej = 0
         do i=1,mxpart
            iTmpAcc(i)=0 ! kf-code of accepted particles
            iTmpRej(i)=0 ! kf-code of rejected particles
         enddo

  400    continue

*--------------------------------------------------------------------

         ic = jnumc(chlw,ic,icl)
         if( ic .gt. icl ) then
            icl = jnumc(chlw,icl+2,i3)
            if( icl .le. i3 ) goto 200
C           -----------------
            if(nAcc > 0 .and. nRej ==0) checkParticleMode(isec)=1
            if(nAcc ==0 .and. nRej > 0) checkParticleMode(isec)=2
            if(nAcc > 0 .and. nRej > 0) checkParticleMode(isec)=3
C           -----------------
            if( .not. allocated(acceptedParticles) ) then
               allocate (acceptedParticles(itlmax,mxpart))
               allocate (rejectedParticles(itlmax,mxpart))
               do i=1,itlmax
               do j=1,mxpart
                 acceptedParticles(i,j)=-1
                 rejectedParticles(i,j)=-1
               enddo
               enddo
            endif
C           -----------------
            do i=1,nAcc
               if(iTmpAcc(i)==0) then ! part = all
                  checkParticleMode(isec)=0
                  goto 140
               endif
               acceptedParticles(isec,i)=iTmpAcc(i)
            enddo
            do i=1,nRej
               rejectedParticles(isec,i)=iTmpRej(i)
            enddo
C           -----------------
            goto 140
         end if

*--------------------------------------------------------------------

         call rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

         if( ierr .eq. 994 ) goto 992
         if( ierr .eq. 998 ) goto 991

*-----------------------------------------------------------------

         inpat = inpat + 1
         if( inpat .gt. mxpart ) goto 990

*-----------------------------------------------------------------

         if    (isubt==0) then ! Accepted particle
            nAcc = nAcc + 1
            if(istyp==19 .and. inkf0==0) then ! part = nucleus
              inkf0=1
              NucleusIsAccepted(isec)=.true.
            endif
            iTmpAcc(nAcc)=inkf0
         elseif(isubt==1) then ! Rejected particle
            nRej = nRej + 1
            if(istyp==19 .and. inkf0==0) inkf0=1
            iTmpRej(nRej)=inkf0
         endif
         goto 400

*-----------------------------------------------------------------------
      else if( ipm==14 ) then ! HistoryMax
         read(Param,*,iostat=ios) historyMax(isec)
         if(ios/=0) goto 997

*-----------------------------------------------------------------------
      else if( ipm==15 ) then ! format
         if(trim(ParamLow)/="t4d" .and. trim(ParamLow)/="pict") then
           goto 986
         endif
         format(isec)=trim(ParamLow)

*-----------------------------------------------------------------------
      else if( ipm==16 ) then ! verbose
         read(Param,*,iostat=ios) verbose(isec)
         if(ios/=0) goto 985

*-----------------------------------------------------------------------
      else if( ipm==17 ) then ! precision
         read(Param,*,iostat=ios) iTMP
         if(ios/=0) goto 984
         if(iTMP<1) goto 983
         write(fmtW(isec),'(a,i0,a,i0,a)') '(es',iTMP+7,'.',iTMP-1,')'

*-----------------------------------------------------------------------
      else if( ipm==18 ) then
         itmp1=inumc(chlw,1,len(chlw),'(')
         itmp2=inumc(chlw,1,len(chlw),')')
         if(itmp1==len(chlw)+1) goto 982
         if(itmp2==len(chlw)+1) goto 982
         n_kfSkip(isec)=n_kfSkip(isec)+1
         ikf = n_kfSkip(isec)
C        -----
         read(chlw(itmp1+1:itmp2-1),*,iostat=ios) iTMP
         if(ios/=0) goto 982
         read(Param                ,*,iostat=ios) tmp
         if(ios/=0) goto 982
         kfSkip     (isec,ikf) = iTMP
         SkipRate_of(isec,ikf) = tmp

*-----------------------------------------------------------------------
      else if( ipm==19 ) then
         read(Param,*,iostat=ios) SkipRate(isec)
         if(ios/=0) goto 981


*-----------------------------------------------------------------------
      end if ! end ipm
      goto 140

*-----------------------------------------------------------------------
  800 continue
      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------
  998    m_err = 'Unknown parameter is found in '//tname
         ErrCha = ''
         ErrID = 'L:689/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  997    m_err = 'Description of parameter is wrong in '//tname
         ErrCha = ''
         ErrID = 'L:694/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  996    m_err = 'AcceptedHistories file is not found in '//tname
         ErrCha = ''
         ErrID = 'L:699/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  995    m_err = 'Description of reg is wrong in '//tname
         ErrCha = ''
         ErrID = 'L:704/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  994    m_err = 'Time unit of HalflifeMinSEC,HalflifeMaxSEC '//
     &           'should be s,min,h,d,y in '//tname
         ErrCha = ''
         ErrID = 'L:710/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  993    m_err = 'AcceptedHistories file is wrong in '//tname
         ErrCha = ''
         ErrID = 'L:715/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  992    m_err = 'Name of particle is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:720/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  991    m_err = 'Description of parameter is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:725/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  990    m_err = 'Number of particles is larger than mxpart in '//tname
         ErrCha = ''
         ErrID = 'L:730/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  989    m_err = 'HistoryFile description is incorrect or empty in '
     &   //tname
         ErrCha = ''
         ErrID = 'L:736/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  988    m_err = '# of T-4Dtrack * # of OpenMP exceeds 100. '//
     &           'Reduce # of T-4Dtrack or # of OpenMP'
         ErrCha = ''
         ErrID = 'L:742/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  987    write(ErrCha,'("warning: '//
     &           'T-4Dtrack is not compatible '//
     &        'with MPI parallelization. '//
     &        'T-4Dtrack is ignored.")')
         ErrID = 'L:749/R:read_t4dtrack/F:talls_4dtrack.f'
         if(me.eq.0)call ErrWrite(ErrID,ErrCha)
         ierr=900
         return

  986    m_err = 'Description of format is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:756/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  985    m_err = 'Description of verbose is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:761/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  984    m_err = 'Description of precision is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:766/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  983    m_err = 'precision should be greater than 0 in tally '//tname
         ErrCha = ''
         ErrID = 'L:771/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  982    m_err = 'Description of SkipRate-of is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:776/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999

  981    m_err = 'Description of SkipRate is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:781/R:read_t4dtrack/F:talls_4dtrack.f'
         goto 999


*-----------------------------------------------------------------------
  999 continue

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

*-----------------------------------------------------------------------
      return
      end subroutine

************************************************************************
      character*10 function pName(kf)
************************************************************************
      integer kf
C     ------------------------------------------------------------------
C     Normal particles
      if(kf==   2212) then; pName="    proton"; return; endif
      if(kf==   2112) then; pName="   neutron"; return; endif
      if(kf==    211) then; pName="     pion+"; return; endif
      if(kf==    111) then; pName="     pion0"; return; endif
      if(kf==   -211) then; pName="     pion-"; return; endif
      if(kf==    -13) then; pName="     muon+"; return; endif
      if(kf==     13) then; pName="     muon-"; return; endif
      if(kf==    321) then; pName="     kaon+"; return; endif
      if(kf==    311) then; pName="     kaon0"; return; endif
      if(kf==   -321) then; pName="     kaon-"; return; endif
      if(kf==     11) then; pName="  electron"; return; endif
      if(kf==    -11) then; pName="  positron"; return; endif
      if(kf==     22) then; pName="    photon"; return; endif
      if(kf==1000002) then; pName="  deuteron"; return; endif
      if(kf==1000003) then; pName="    triton"; return; endif
      if(kf==2000003) then; pName="       3He"; return; endif
      if(kf==2000004) then; pName="     alpha"; return; endif
      if(isNucleus(kf))then;pName="   nucleus"; return; endif
      if(kf==     12) then; pName="      nu_e"; return; endif
      if(kf==    -12) then; pName="  nu_e_bar"; return; endif
      if(kf==     14) then; pName="     nu_mu"; return; endif
      if(kf==    -14) then; pName=" nu_mu_bar"; return; endif
      if(kf==     16) then; pName="    nu_tau"; return; endif
      if(kf==    -16) then; pName="nu_tau_bar"; return; endif
      if(kf==    130) then; pName="   K0-long"; return; endif
      if(kf==    310) then; pName="  K0-short"; return; endif
C     ------------------------------------------------------------------
C     User defined particles
      if(900000<=abs(kf).and.abs(kf)<=999999) then
        write(pName, '(I10)') kf
        return
      endif
C     ------------------------------------------------------------------
C     Handling of Other Particles
      if(.true.) then
        write(pName, '(I10)') kf
        return
      else
        pName="    others"
        return
      endif
C     ------------------------------------------------------------------
      end function

************************************************************************
      integer function kf_of(x)
************************************************************************
      character*(*) x
      read(x, *, iostat=ios) kf_of

      ! x is an integer.
      if(ios==0) return

      ! x is a string.
      if(trim(x)=="proton"  )then; kf_of=   2212; return; endif
      if(trim(x)=="neutron" )then; kf_of=   2112; return; endif
      if(trim(x)=="pion+"   )then; kf_of=    211; return; endif
      if(trim(x)=="pion0"   )then; kf_of=    111; return; endif
      if(trim(x)=="pion-"   )then; kf_of=   -211; return; endif
      if(trim(x)=="muon+"   )then; kf_of=    -13; return; endif
      if(trim(x)=="muon-"   )then; kf_of=     13; return; endif
      if(trim(x)=="kaon+"   )then; kf_of=    321; return; endif
      if(trim(x)=="kaon0"   )then; kf_of=    311; return; endif
      if(trim(x)=="kaon-"   )then; kf_of=   -321; return; endif
      if(trim(x)=="electron")then; kf_of=     11; return; endif
      if(trim(x)=="positron")then; kf_of=    -11; return; endif
      if(trim(x)=="photon"  )then; kf_of=     22; return; endif
      if(trim(x)=="deuteron")then; kf_of=1000002; return; endif
      if(trim(x)=="triton"  )then; kf_of=1000003; return; endif
      if(trim(x)=="3he"     )then; kf_of=2000003; return; endif
      if(trim(x)=="alpha"   )then; kf_of=2000004; return; endif
      if(trim(x)=="nucleus" )then; kf_of=      1; return; endif
      if(trim(x)=="all"     )then; kf_of=      0; return; endif
      print*,"error: Description of part is wrong."
      stop
      end function


************************************************************************
      subroutine t4dtrack(ncol)
      do i=1,nSection
        isec=i
        call t4dtrack_of(isec,ncol)
      enddo
      end subroutine


************************************************************************
*                                                                      *
      subroutine t4dtrack_of(isec,ncol)
*                                                                      *
*        sample subroutine for t4dtrack.                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*       ncol = 1 : start of calculation                                *
*              2 : end of calculation                                  *
*              3 : end of a batch                                      *
*              4 : source                                              *
*              5 : detection of geometry error                         *
*              6 : recovery of geometry error                          *
*              7 : termination by geometry error                       *
*              8 : termination by weight cut-off                       *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : (n,x) reaction                                      *
*             14 : (n,n'x) reaction                                    *
*             15 : sequential transport only for tally                 *
*             16 : surface cross for WW of xyz mesh                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        In the distributed-memory parallel computing,                 *
*          npe : total number of used Processor Elements               *
*          me : ID number of each processor                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        In the shared memory parallel computing,                      *
*          ipomp : ID number of each core                              *
*          npomp : total number of used cores                          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        nocas : current event number in this batch                    :
*        nobch : current batch number                                  *
*        rcasc : real number of NOCAS+maxcas*(NOBCH-1)                 *
*        rsouin : sum of the weight of source particle                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        no : cascade id in this event                                 *
*        idmn(mat) : material id                                       *
*        ityp : particle type                                          *
*        ktyp : particle kf-code                                       *
*        jtyp : charge number of the particle                          *
*        mtyp : baryon number of the particle                          *
*        rtyp : rest mass of the particle (MeV)                        *
*        oldwt : weight of the particle at (x,y,z)                     *
*        qs : dE/dx of electron at (x,y,z)                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        iblz1 : cell id at (x,y,z)                                    *
*        iblz2 : cell id after crossing                                *
*        ilev1 : level structure id of the cell at (x,y,z)             *
*          ilat1(i,j)                                                  *
*        ilev2 : level structure id of the cell after crossing         *
*          ilat2(i,j)                                                  *
*        costha : cosine of theta on surface crossing                  *
*        uang(1) : x of position(?) at surface crossing                *
*        uang(2) : y of position(?) at surface crossing                *
*        uang(3) : z of position(?) at surface crossing                *
*        nsurf : surface number                                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        name(no,ipomp+1) : collision number of the particle                   *
*        ncnt(1,no,ipomp+1) : values of counter 1                              *
*        ncnt(2,no,ipomp+1) : values of counter 2                              *
*        ncnt(3,no,ipomp+1) : values of counter 3                              *
*        wt(no,ipomp+1) : weight of the particle at (xc,yc,zc)                 *
*        u(no,ipomp+1) : x, y, z-components of unit vector of                  *
*        v(no,ipomp+1) :momentum of the particle                               *
*        w(no,ipomp+1) :                                                       *
*        e(no,ipomp+1) : energy of the particle at (x,y,z) (MeV)               *
*        t(no,ipomp+1) : time of the particle at (x,y,z) (nsec)                *
*        x(no,ipomp+1) : x, y, z-position coordinates of                       *
*        y(no,ipomp+1) :the preceding event point (cm)                         *
*        z(no,ipomp+1) :                                                       *
*        ec(no,ipomp+1) : energy of the particle at (xc,yc,zc) (MeV)           *
*        tc(no,ipomp+1) : time of the particle at (xc,yc,zc) (nsec)            *
*        xc(no,ipomp+1) : x, y, z-position coordinates of                      *
*        yc(no,ipomp+1) :the particle (cm)                                     *
*        zc(no,ipomp+1) :                                                      *
*        spx(no,ipomp+1) : x, y, z-components of unit vector of                *
*        spy(no,ipomp+1) :spin direction of the particle                       *
*        spz(no,ipomp+1) :                                                     *
*        nzst(no,ipomp+1)                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        nclsts : the number of produced particle and nucleus          *
*        mathz : Z number of the mother nucleus                        *
*        mathn : N number of the mother nucleus                        *
*        jcoll : reaction type id1                                     *
*        kcoll : reaction type id2                                     *
*                                                                      *
*----------------------------------------------------------------------*
*        jcoll reaction type identifier                                *
*----------------------------------------------------------------------*
*                                                                      *
*        jcoll : =  0, nothing happen                                  *
*                =  1, Hydrogen collisions                             *
*                =  2, Particle Decays                                 *
*                =  3, Elastic collisions                              *
*                =  4, High Energy Nuclear collisions                  *
*                =  5, Heavy Ion reactions                             *
*                =  6, Neutron reactions by data                       *
*                =  7, Photon reactions by data                        *
*                =  8, Electron reactions by data                      *
*                =  9, P,d,a, and photo-nuclear reactions by data      *
*                = 10, Neutron event mode                              *
*                = 11, Delta Ray production                            *
*                = 12, Muon atomic interaction                         *
*                = 13, Photon by EGS5                                  *
*                = 14, Electron by EGS5                                *
*                = 15, Photon photonuclear interaction                 *
*                = 16, Negative muon captured by nucleon               *
*                = 17, Muon photonuclear interaction                   *
*                                                                      *
*----------------------------------------------------------------------*
*        kcoll reaction type identifier                                *
*----------------------------------------------------------------------*
*                                                                      *
*        kcoll : =  0, normal                                          *
*                =  1, high energy fission                             *
*                =  2, high energy absorption                          *
*                =  3, low energy n elastic                            *
*                =  4, low energy n non-elastic                        *
*                =  5, low energy n fission                            *
*                =  6, low energy n absorption                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        nclsts   : total number of out going particles and nuclei     *
*                                                                      *
*        iclusts(nclsts)                                               *
*                                                                      *
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclusts(i,nclsts)                                             *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4, status of the particle 0: real, <0 : dead      *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        qclusts(i,nclsts)                                             *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, x-component of unit vector of momentum         *
*                  = 2, y-component of unit vector of momentum         *
*                  = 3, z-component of unit vector of momentum         *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight                                         *
*                  = 9, time                                           *
*                  = 10, x                                             *
*                  = 11, y                                             *
*                  = 12, z                                             *
*                                                                      *
*        jcount(i,nclsts)                                              *
*                                                                      *
************************************************************************
*                                                                      *
*        kf code table                                                 *
*                                                                      *
*           kf-code: ityp :  description                               *
*                                                                      *
*             2212 :   1  :  proton                                    *
*             2112 :   2  :  neutron                                   *
*              211 :   3  :  pion (+)                                  *
*              111 :   4  :  pion (0)                                  *
*             -211 :   5  :  pion (-)                                  *
*              -13 :   6  :  muon (+)                                  *
*               13 :   7  :  muon (-)                                  *
*              321 :   8  :  kaon (+)                                  *
*              311 :   9  :  kaon (0)                                  *
*             -321 :  10  :  kaon (-)                                  *
*                                                                      *
*               11 :  12  :  electron                                  *
*              -11 :  13  :  positron                                  *
*               22 :  14  :  photon                                    *
*          1000002 :  15  :  deuteron                                  *
*          1000003 :  16  :  triton                                    *
*          2000003 :  17  :  3he                                       *
*          2000004 :  18  :  alpha                                     *
*      Z*1000000+A :  19  :  nucleus                                   *
*                                                                      *
*           kf-code of the other transport particles (ityp=11)         *
*               12 :         nu_e                                      *
*               14 :         nu_mu                                     *
*              221 :         eta                                       *
*              331 :         eta'                                      *
*             -311 :         k0bar                                     *
*              130 :         K_L0                                      *
*              310 :         K_S0                                      *
*            -2112 :         nbar                                      *
*            -2212 :         pbar                                      *
*             3122 :         Lambda0                                   *
*             3222 :         Sigma+                                    *
*             3212 :         Sigma0                                    *
*             3112 :         Sigma-                                    *
*             3322 :         Xi0                                       *
*             3312 :         Xi-                                       *
*             3334 :         Omega-                                    *
*                                                                      *
************************************************************************

      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      common /cusrtally/ iusrtally, iudtf(50)
      common /cudtpara/ udtpara(0:9)

*-----------------------------------------------------------------------

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /sors0/ rcasc00,rsouin00,initsor
!$OMP THREADPRIVATE(/sors0/)

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /wtsav2/ oldwt2
!$OMP THREADPRIVATE(/wtsav2/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /tlcost/ costha, uang(3), nsurf
!$OMP THREADPRIVATE(/tlcost/)
      common /celdb/  idsn(kvlmax), idtn(kvlmax)

cKN 2012/10/31
      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, noz, mtel
!$OMP THREADPRIVATE(/elect/)

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

      integer*8 :: iransb64 ! S.H. xorshift (2020.2.6)
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)

      common /taliin/ rsouin, nzztin, nrgnin
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

      integer unit_number

      logical isAccepted
      logical TrackStart
      logical TrackEnd_Yes_No
      logical l_tmp
      logical :: isFound = .false.

*-----------------------------------------------------------------------
*       ncol = 1 : start of calculation                                *
*              2 : end of calculation                                  *
*              3 : end of a batch                                      *
*              4 : source                                              *
*              5 : detection of geometry error                         *
*              6 : recovery of geometry error                          *
*              7 : termination by geometry error                       *
*              8 : termination by weight cut-off                       *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : (n,x) reaction                                      *
*             14 : (n,n'x) reaction                                    *
*             15 : sequential transport only for tally                 *
*             16 : surface cross for WW of xyz mesh                    *
*-----------------------------------------------------------------------

      unit_number=UnitNumber(isec)

C     ==================================================================
      if( ncol==4 ) then ! 4 : source
        no_pre(isec)=0 ! initialize

C       First, close the previous track, nuclei buffer and the source.
        if(.not. is_RejectedEvent(isec)) then
          call line("end track",unit_number,0)
          call line("end source",unit_number,0)
        endif
        nBufNuc(isec)=0 ! initialize buffered nuclei

C       ------------
C       Event veto
        if(n_AcceptedHistories(isec) > 0) then
          is_RejectedEvent(isec)=.true.
          do i=1, n_AcceptedHistories(isec)
            if( nocas==AcceptedHistories(isec,i,1) .and.
     &          nobch==AcceptedHistories(isec,i,2) ) then
              is_RejectedEvent(isec)=.false. ! This is accepted event.
              exit
            endif
          enddo
        endif
        if(is_RejectedEvent(isec)) return

C       ------------
C       Event veto
        if(historyMax(isec) >= 0) then

          if(nt4domp.gt.0)then
            l_tmp = my_rcasc() <= historyMax(isec) ! OMP mode
          else
            l_tmp = n_WrittenHistories(isec) < historyMax(isec) ! SINGLE mode
          endif

          if(l_tmp) then
            is_RejectedEvent(isec)=.false. ! This is accepted event.
          else
            is_RejectedEvent(isec)=.true.  ! This is rejected event.
          endif
        endif

C       ------------
C       Let's start new source!
        call line("new source",unit_number,0)
      endif

C     ==================================================================
      if( 8<=ncol .and. ncol<=15 ) then ! write step

        if(is_RejectedEvent(isec)) return

C       ---------------
        TrackStart = (no/=no_pre(isec)) ! no = track number
C       ----
        if(checkCellMode(isec)) then
          TrackStart_No_Yes(isec) = ncol==10 .and.
     &      ( (.not. acceptedCells(isec,iblz1))
     &         .and. acceptedCells(isec,iblz2) )
          ! 10 : geometry boundary crossing
          ! From "NoCell" to "YesCell"
          TrackStart = TrackStart .or. TrackStart_No_Yes(isec)
        endif

C       ---------------
C       buffer_nuclei within the specified half-life.
        if(PinnedOption(isec) > 0) then
          if((ncol==13 .or. ncol==14) .and. nclsts>=2) then
            do i=1,nclsts
              if(checkCellMode(isec)) then
                if(.not. acceptedCells(isec,iblz2)) return
              endif
              kfSave=jclusts(7,i)
              if(nucleiInHalfLife(kfSave)) then ! Find the target nucleus from the final state.
                nBufNuc(isec) = nBufNuc(isec) + 1
                iSave=i
                if(nBufNuc(isec) > nBufNucMax) cycle
                call line("buffer_nuclei",unit_number,ncol)
              endif
            enddo
          endif
        endif

C       ---------------
        if(TrackStart) then ! Start point found.
          no_pre(isec)=no ! update track number 'no'
C         ---------------
C         First, close the previous track,
          call line("end track",unit_number,ncol)
C         ---------------
C         (If buffered nuclear information exists, write it out at the end of the track.)
          call line("write nuclei",unit_number,ncol)
C         ---------------
C         and, buffer the starting point of the track.
          call line("buffer_start_point",0,ncol)
          call update_StartPoint(ncol)
        endif

C       ---------------
        if(checkCellMode(isec)) then
          TrackEnd_Yes_No = (
     &      ncol==10 .and.
     &      ( acceptedCells(isec,iblz1) .and.
     & (.not. acceptedCells(isec,iblz2)) )
     &    )
          if(TrackEnd_Yes_No) then
C           First, writing out the buffer,
            if(buff_exist(isec))
     &        call line("write_start_point",unit_number,ncol)
C           then writing out the step.
            call line("write_end_point",unit_number,ncol)
            call line("end track",unit_number,ncol)
            return
          endif
        endif
        ! From "YesCell" to "NoCell"


C       ---------------
C       Conditions for vetoing the step.
C       ---------------
        if(stepminEM(isec) > 0d0) then
          if(StepLength() < stepminEM(isec)) then
            if(abs(ktyp)==11 .or. ktyp==22) then
              if(isSkip(ktyp,ncol,jcoll,kcoll,nclsts)) return
            endif
          endif
        endif
C       ---------------
        if(checkCellMode(isec)) then
          if(.not. acceptedCells(isec,iblz2)) return
          if((.not. acceptedCells(isec,iblz1)) .and. ncol==12) return ! This happens when going from the rejected cell to the outer void.
        endif
C       ---------------
        if(checkParticleMode(isec) > 0) then
          isAccepted=.false.
C         condition 1 ----------
          if     (checkParticleMode(isec)==1) then
            isAccepted = isPartContained("Acc",ktyp)
          else if(checkParticleMode(isec)==2) then
            isAccepted = (.not.isPartContained("Rej",ktyp))
          else if(checkParticleMode(isec)==3) then
            isAccepted = isPartContained("Acc",ktyp)
     &      .and. ( .not.isPartContained("Rej",ktyp) )
          endif
C         condition 2 ----------
          isAccepted = isAccepted
     &       .or. (NucleusIsAccepted(isec) .and. nucleiInHalfLife(ktyp))
          if(.not.isAccepted) return
        endif
C       ---------------
        isFound = .false.
        if(n_kfSkip(isec) > 0) then
C         Search for the target KF-code
          do ikf=1,n_kfSkip(isec)
            if(ktyp==kfSkip(isec, ikf)) then
              isFound = .true.
              exit
            endif
          enddo
C         ---
          if(isFound) then
            if(SkipRate_of(isec,ikf) > random_T4D()) then
              if(isSkip(ktyp,ncol,jcoll,kcoll,nclsts)) return
            endif
          endif
        endif
C       ---------------
        if(SkipRate(isec) > 0d0) then
          if(.not. isFound) then
            if(SkipRate(isec) > random_T4D()) then
              if(isSkip(ktyp,ncol,jcoll,kcoll,nclsts)) return
            endif
          endif
        endif
C       ---------------


C       ---------------
C       First, writing out the buffer,
        if(buff_exist(isec))
     &    call line("write_start_point",unit_number,ncol)

C       then writing out the step.
        if(.not.TrackStart_No_Yes(isec))
     &    call line("write_end_point",unit_number,ncol)
        call update_StartPoint(ncol)

        if(PinnedOption(isec) /= 0) then
          if(isStop(ktyp,ncol)) then
            if(isPinnedParticle(ktyp))
     &        call line("write_end_point(Pinned)",unit_number,ncol)
          endif
        endif
C       ---------------

C     ==================================================================
      else if( ncol==-1 ) then ! -1 : start of batch calculation
        no_pre(isec)=0           ! initialize

        call opent4dompfile(unit_number,isec)

C     ==================================================================
      else if( ncol==1 ) then ! 1 : start of calculation
        no_pre(isec)=0 ! initialize
        open(2000+isec,file=filename(isec),
     &                                  status='replace',action='write')
        call write_header(2000+isec)




C     ==================================================================
      else if( ncol==2 ) then ! 2: end of calculation
        close(2000+isec)

C     ==================================================================
      else if( ncol==-3 ) then ! 3: end of batch
        call line("end track",unit_number,0)
        call line("end source",unit_number,0)
        call cpt4dompfile(unit_number,isec,2000+isec)
        call line("new source",unit_number,0)

C     ==================================================================
      endif
      return
      end subroutine


************************************************************************
      double precision function StepLength()
************************************************************************
      use MMBANKMOD
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
      StepLength=sqrt(
     &  (x_start(isec)-xc(no,ipomp+1))**2
     & +(y_start(isec)-yc(no,ipomp+1))**2
     & +(z_start(isec)-zc(no,ipomp+1))**2)
      return
      end function
************************************************************************


************************************************************************
      subroutine update_StartPoint(ncol)
************************************************************************
      use MMBANKMOD
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      if(ncol==10) then ! 10 : geometry boundary crossing
        x_start(isec)=xc(no,ipomp+1)
        y_start(isec)=yc(no,ipomp+1)
        z_start(isec)=zc(no,ipomp+1)
      else
        x_start(isec)=x(no,ipomp+1)
        y_start(isec)=y(no,ipomp+1)
        z_start(isec)=z(no,ipomp+1)
      endif
      end subroutine


************************************************************************
      double precision function PinnedTime(kf)
************************************************************************
      PinnedTime=0d0
      if(PinnedTimeOption(isec)==0) then
        !0: Fix
        PinnedTime=PinnedTimeFix(isec)
      elseif( PinnedTimeOption(isec)==1 .or.
     &        PinnedTimeOption(isec)==2 ) then
        !1: Half-life of unstable nuclei
        !2: aHalf * log10(Half-life/s) + bHalf
        iZ=kf / 1000000
        iA=mod(kf, 1000000)
        dT_SEC=halflife(iZ,iA) ! Change to a sampling method later.
        dT_ns=dT_SEC*1e+9 ! s -> ns
        if(PinnedTimeOption(isec)==1)
     &                 PinnedTime=aHalf(isec)*dT_ns
        if(PinnedTimeOption(isec)==2)
     &                 PinnedTime=aHalf(isec)*log10(dT_ns)+bHalf(isec)
        if(PinnedTime < 0d0) PinnedTime=0d0
      endif
      return
      end function

************************************************************************
      integer function my_rcasc()
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /cparm/  maxbch,maxcas
      my_rcasc=nocas+maxcas*(nobch-1)
      return
      end function

************************************************************************
      subroutine line(type,unit_number,ncol)
C     write line
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)
*-----------------------------------------------------------------------
      include 'param00.inc'
      include 'param.inc'

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /sors0/ rcasc00,rsouin00,initsor
!$OMP THREADPRIVATE(/sors0/)

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /wtsav2/ oldwt2
!$OMP THREADPRIVATE(/wtsav2/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)

      character*(*) type
      integer unit_number
      character*99 cTMP

C     ------------------------------------------------------------------
      if(type=="write_end_point".or.type=="write_end_point(Pinned)")then
        if(n_line_of_this_track(isec)==0) return
        n_line_of_this_track(isec) = n_line_of_this_track(isec) + 1
        n_line_of_this_file(isec) = n_line_of_this_file(isec) + 1
        time=tc(no,ipomp+1)
        if(type=="write_end_point(Pinned)")
     &    time = time + PinnedTime(ktyp)
C       ~~~~~~~~~~
        if(trim(format(isec))=="t4d") then
          write(unit_number,fmtW(isec),advance='no') xc(no,ipomp+1)
          write(unit_number,fmtW(isec),advance='no') yc(no,ipomp+1)
          write(unit_number,fmtW(isec),advance='no') zc(no,ipomp+1)
          write(unit_number,fmtW(isec),advance='no') time
          write(unit_number,fmtW(isec),advance='no') ec(no,ipomp+1)
          write(unit_number,fmtW(isec),advance='no') wt(no,ipomp+1)
          if(verbose(isec) >= 1) then
          write(unit_number,'(i8)',advance='no') iblz2
          write(unit_number,'(i8)',advance='no') iblz1
          endif
          if(verbose(isec) >= 2) then
          write(unit_number,'(i3)',advance='no') ncol
          write(unit_number,'(i3)',advance='no') jcoll
          write(unit_number,'(i3)',advance='no') kcoll
          write(unit_number,'(x,i5)',advance='no') nclsts
          write(unit_number,'(x,i0)',advance='no') ncnt(1,no,ipomp+1)
          write(unit_number,'(x,i0)',advance='no') ncnt(2,no,ipomp+1)
          write(unit_number,'(x,i0)',advance='no') ncnt(3,no,ipomp+1)
          endif
C       ~~~~~~~~~~
        elseif(trim(format(isec))=="pict") then
          write(unit_number,'(i10)'   ,advance='no') ktyp
          write(unit_number,fmtW(isec),advance='no') xc(no,ipomp+1)
          write(unit_number,fmtW(isec),advance='no') yc(no,ipomp+1)
          write(unit_number,fmtW(isec),advance='no') zc(no,ipomp+1)
          write(unit_number,fmtW(isec),advance='no') ec(no,ipomp+1)
          write(unit_number,'(i8)'    ,advance='no') iblz2
          write(unit_number,fmtW(isec),advance='no') wt(no,ipomp+1)
          write(unit_number,fmtW(isec),advance='no') time
          if(verbose(isec) >= 1) then
          write(unit_number,'(i8)',advance='no') iblz1
          endif
          if(verbose(isec) >= 2) then
          write(unit_number,'(a)',advance='no') pName(ktyp)
          write(unit_number,'(i6)',advance='no') no
          write(unit_number,'(i3)',advance='no') ncol
          write(unit_number,'(i3)',advance='no') jcoll
          write(unit_number,'(i3)',advance='no') kcoll
          write(unit_number,'(x,i5)',advance='no') nclsts
          write(unit_number,'(x,i0)',advance='no') ncnt(1,no,ipomp+1)
          write(unit_number,'(x,i0)',advance='no') ncnt(2,no,ipomp+1)
          write(unit_number,'(x,i0)',advance='no') ncnt(3,no,ipomp+1)
          endif
        endif
C       ~~~~~~~~~~
        write(unit_number,*) ! new line
C     ------------------------------------------------------------------
      else if(type=="buffer_start_point") then
        if(TrackStart_No_Yes(isec)) then
C         ~~~~~~~~~~
          if(trim(format(isec))=="t4d") then
            write(StartPoint(isec, 1),fmtW(isec)) xc(no,ipomp+1)
            write(StartPoint(isec, 2),fmtW(isec)) yc(no,ipomp+1)
            write(StartPoint(isec, 3),fmtW(isec)) zc(no,ipomp+1)
            write(StartPoint(isec, 4),fmtW(isec)) tc(no,ipomp+1)
            write(StartPoint(isec, 5),fmtW(isec)) ec(no,ipomp+1)
            write(StartPoint(isec, 6),fmtW(isec)) wt(no,ipomp+1)
            write(StartPoint(isec, 7),'(i8)') iblz1
            write(StartPoint(isec, 8),'(i8)') iblz2
            write(StartPoint(isec, 9),'(i3)') ncol
            write(StartPoint(isec,10),'(i3)') jcoll
            write(StartPoint(isec,11),'(i3)') kcoll
            write(StartPoint(isec,12),'(x,i5)') nclsts
            write(StartPoint(isec,13),'(x,i0)') ncnt(1,no,ipomp+1)
            write(StartPoint(isec,14),'(x,i0)') ncnt(2,no,ipomp+1)
            write(StartPoint(isec,15),'(x,i0)') ncnt(3,no,ipomp+1)
C           save for track header
            write(StartPoint(isec,20),'(i10)') ktyp
            write(StartPoint(isec,19),'(a)'  ) pName(ktyp)
            write(StartPoint(isec,18),'(i6)' ) no
C         ~~~~~~~~~~
          elseif(trim(format(isec))=="pict") then
            write(StartPoint(isec, 1),'(i10)'   ) ktyp
            write(StartPoint(isec, 2),fmtW(isec)) xc(no,ipomp+1)
            write(StartPoint(isec, 3),fmtW(isec)) yc(no,ipomp+1)
            write(StartPoint(isec, 4),fmtW(isec)) zc(no,ipomp+1)
            write(StartPoint(isec, 5),fmtW(isec)) ec(no,ipomp+1)
            write(StartPoint(isec, 6),'(i8)'    ) iblz1
            write(StartPoint(isec, 7),fmtW(isec)) wt(no,ipomp+1)
            write(StartPoint(isec, 8),fmtW(isec)) tc(no,ipomp+1)
            write(StartPoint(isec, 9),'(i8)') iblz2
            write(StartPoint(isec,10),'(a)') pName(ktyp)
            write(StartPoint(isec,11),'(i6)') no
            write(StartPoint(isec,12),'(i3)') ncol
            write(StartPoint(isec,13),'(i3)') jcoll
            write(StartPoint(isec,14),'(i3)') kcoll
            write(StartPoint(isec,15),'(x,i5)') nclsts
            write(StartPoint(isec,16),'(x,i0)') ncnt(1,no,ipomp+1)
            write(StartPoint(isec,17),'(x,i0)') ncnt(2,no,ipomp+1)
            write(StartPoint(isec,18),'(x,i0)') ncnt(3,no,ipomp+1)
          endif
C         ~~~~~~~~~~
        else
C         ~~~~~~~~~~
          if(trim(format(isec))=="t4d") then
            write(StartPoint(isec, 1),fmtW(isec)) x(no,ipomp+1)
            write(StartPoint(isec, 2),fmtW(isec)) y(no,ipomp+1)
            write(StartPoint(isec, 3),fmtW(isec)) z(no,ipomp+1)
            write(StartPoint(isec, 4),fmtW(isec)) t(no,ipomp+1)
            write(StartPoint(isec, 5),fmtW(isec)) e(no,ipomp+1)
            write(StartPoint(isec, 6),fmtW(isec)) oldwt
            write(StartPoint(isec, 7),'(i8)') iblz1
            write(StartPoint(isec, 8),'(i8)') iblz2
            write(StartPoint(isec, 9),'(i3)') ncol
            write(StartPoint(isec,10),'(i3)') jcoll
            write(StartPoint(isec,11),'(i3)') kcoll
            write(StartPoint(isec,12),'(x,i5)') nclsts
            write(StartPoint(isec,13),'(x,i0)') ncnt(1,no,ipomp+1)
            write(StartPoint(isec,14),'(x,i0)') ncnt(2,no,ipomp+1)
            write(StartPoint(isec,15),'(x,i0)') ncnt(3,no,ipomp+1)
C           save for track header
            write(StartPoint(isec,20),'(i10)') ktyp
            write(StartPoint(isec,19),'(a)'  ) pName(ktyp)
            write(StartPoint(isec,18),'(i6)' ) no
C         ~~~~~~~~~~
          elseif(trim(format(isec))=="pict") then
            write(StartPoint(isec, 1),'(i10)'   ) ktyp
            write(StartPoint(isec, 2),fmtW(isec)) x(no,ipomp+1)
            write(StartPoint(isec, 3),fmtW(isec)) y(no,ipomp+1)
            write(StartPoint(isec, 4),fmtW(isec)) z(no,ipomp+1)
            write(StartPoint(isec, 5),fmtW(isec)) e(no,ipomp+1)
            write(StartPoint(isec, 6),'(i8)'    ) iblz1
            write(StartPoint(isec, 7),fmtW(isec)) oldwt
            write(StartPoint(isec, 8),fmtW(isec)) t(no,ipomp+1)
            write(StartPoint(isec, 9),'(i8)') iblz2
            write(StartPoint(isec,10),'(a)') pName(ktyp)
            write(StartPoint(isec,11),'(i6)') no
            write(StartPoint(isec,12),'(i3)') ncol
            write(StartPoint(isec,13),'(i3)') jcoll
            write(StartPoint(isec,14),'(i3)') kcoll
            write(StartPoint(isec,15),'(x,i5)') nclsts
            write(StartPoint(isec,16),'(x,i0)') ncnt(1,no,ipomp+1)
            write(StartPoint(isec,17),'(x,i0)') ncnt(2,no,ipomp+1)
            write(StartPoint(isec,18),'(x,i0)') ncnt(3,no,ipomp+1)
          endif
C         ~~~~~~~~~~
        endif
        buff_exist(isec)=.true.
C     ------------------------------------------------------------------
      else if(type=="buffer_nuclei") then
C       ~~~~~~~~~~
        if(trim(format(isec))=="t4d") then
        write(BufNuc(isec,nBufNuc(isec), 1),fmtW(isec))qclusts(10,iSave) ! x
        write(BufNuc(isec,nBufNuc(isec), 2),fmtW(isec))qclusts(11,iSave) ! y
        write(BufNuc(isec,nBufNuc(isec), 3),fmtW(isec))qclusts(12,iSave) ! z
        write(BufNuc(isec,nBufNuc(isec), 4),fmtW(isec))qclusts( 9,iSave) ! time   ! If the order of 'time' changes, modify 'j==2 .and. i==4'
        write(BufNuc(isec,nBufNuc(isec), 5),fmtW(isec))qclusts( 7,iSave) ! kinetic energy (MeV)
        write(BufNuc(isec,nBufNuc(isec), 6),fmtW(isec))qclusts( 8,iSave) ! weight
        write(BufNuc(isec,nBufNuc(isec), 7),'(i8)') iblz1
        write(BufNuc(isec,nBufNuc(isec), 8),'(i8)') iblz2
        write(BufNuc(isec,nBufNuc(isec), 9),'(i3)') ncol
        write(BufNuc(isec,nBufNuc(isec),10),'(i3)') jcoll
        write(BufNuc(isec,nBufNuc(isec),11),'(i3)') kcoll
        write(BufNuc(isec,nBufNuc(isec),12),'(x,i5)') nclsts
        write(BufNuc(isec,nBufNuc(isec),13),'(x,i0)') ncnt(1,no,ipomp+1)
        write(BufNuc(isec,nBufNuc(isec),14),'(x,i0)') ncnt(2,no,ipomp+1)
        write(BufNuc(isec,nBufNuc(isec),15),'(x,i0)') ncnt(3,no,ipomp+1)
C       save for track header
        write(BufNuc(isec,nBufNuc(isec),20),'(i10)') kfSave
        write(BufNuc(isec,nBufNuc(isec),19),'(a)'  ) pName(kfSave)
        write(BufNuc(isec,nBufNuc(isec),18),'(i6)' ) no
C       ~~~~~~~~~~
        elseif(trim(format(isec))=="pict") then
        write(BufNuc(isec,nBufNuc(isec), 1),'(i10)'   )kfSave   ! If the order of 'kf' changes, modify 'read(BufNuc(isec,iB,1),*) kf'.
        write(BufNuc(isec,nBufNuc(isec), 2),fmtW(isec))qclusts(10,iSave) ! x
        write(BufNuc(isec,nBufNuc(isec), 3),fmtW(isec))qclusts(11,iSave) ! y
        write(BufNuc(isec,nBufNuc(isec), 4),fmtW(isec))qclusts(12,iSave) ! z
        write(BufNuc(isec,nBufNuc(isec), 5),fmtW(isec))qclusts( 7,iSave) ! kinetic energy (MeV)
        write(BufNuc(isec,nBufNuc(isec), 6),'(i8)'    )iblz1
        write(BufNuc(isec,nBufNuc(isec), 7),fmtW(isec))qclusts(8,iSave) ! weight
        write(BufNuc(isec,nBufNuc(isec), 8),fmtW(isec))qclusts(9,iSave) ! time   ! If the order of 'time' changes, modify 'j==2 .and. i==8'
        write(BufNuc(isec,nBufNuc(isec), 9),'(i8)') iblz2
        write(BufNuc(isec,nBufNuc(isec),10),'(a)') pName(kfSave)
        write(BufNuc(isec,nBufNuc(isec),11),'(i6)') no ! no
        write(BufNuc(isec,nBufNuc(isec),12),'(i3)') ncol
        write(BufNuc(isec,nBufNuc(isec),13),'(i3)') jcoll
        write(BufNuc(isec,nBufNuc(isec),14),'(i3)') kcoll
        write(BufNuc(isec,nBufNuc(isec),15),'(x,i5)') nclsts
        write(BufNuc(isec,nBufNuc(isec),16),'(x,i0)') ncnt(1,no,ipomp+1)
        write(BufNuc(isec,nBufNuc(isec),17),'(x,i0)') ncnt(2,no,ipomp+1)
        write(BufNuc(isec,nBufNuc(isec),18),'(x,i0)') ncnt(3,no,ipomp+1)
        endif
C       ~~~~~~~~~~
C     ------------------------------------------------------------------
      else if(type=="write_start_point") then
        n_line_of_this_source(isec) = n_line_of_this_source(isec) + 1
        n_line_of_this_track(isec) = n_line_of_this_track(isec) + 1
        n_line_of_this_file(isec) = n_line_of_this_file(isec) + 1
        if(n_line_of_this_source(isec)==1) then
          ! Write out the source number only for sources where a track exists.
C         ~~~~~~~~~~
          if(trim(format(isec))=="t4d") then
            write(unit_number,'(a)'   ,advance='no') "h:"
            write(unit_number,'(x,i0)',advance='no') my_rcasc()
            if(verbose(isec) >= 1) then
            write(unit_number,'(x,i0)',advance='no') int(nocas)
            write(unit_number,'(x,i0)',advance='no') int(nobch)
            endif
            write(unit_number,*) ! new line
C         ~~~~~~~~~~
          elseif(trim(format(isec))=="pict") then
            write(unit_number,'(i1,x,i9)')
     &        0, my_rcasc()
          endif
C         ~~~~~~~~~~
        endif
C       ~~~~~~~~~~
        if(trim(format(isec))=="t4d") then
          imax=6
          if(verbose(isec)==1) imax=8
          if(verbose(isec)==2) imax=15
          read(StartPoint(isec,20),*) i1 ! ktyp
          read(StartPoint(isec,18),*) i2 ! no
          write(unit_number,'(a)'   ,advance='no') "t:"
          write(unit_number,'(x,i0)',advance='no') i1 ! ktyp
          if(verbose(isec) >= 1) then
          cTMP=trim(StartPoint(isec,19)) ! pName(ktyp)
          write(unit_number,'(x,a)' ,advance='no') trim(cTMP) ! pName(ktyp)
          endif
          if(verbose(isec) >= 2) then
          write(unit_number,'(x,i0)',advance='no') i2 ! no
          endif
          write(unit_number,*) ! new line
C       ~~~~~~~~~~
        elseif(trim(format(isec))=="pict") then
          imax=8
          if(verbose(isec)==1) imax=9
          if(verbose(isec)==2) imax=18
        endif
C       ~~~~~~~~~~
        do i=1, imax
          write(unit_number,'(a)',advance='no') trim(StartPoint(isec,i))
        enddo
        write(unit_number,*) ! new line
        buff_exist(isec)=.false.
        lastLineWasEndTrack(isec)=.false.
C     ------------------------------------------------------------------
      else if(type=="write nuclei") then
        do iB=1,nBufNuc(isec)
          n_line_of_this_source(isec) = n_line_of_this_source(isec) + 1
          n_line_of_this_track(isec) = n_line_of_this_track(isec) + 1
          n_line_of_this_file(isec) = n_line_of_this_file(isec) + 1
C         ~~~~~~~~~~
          if(trim(format(isec))=="t4d") then
            do j=1,2 ! Write out two lines. Add pinned time to the second line.
              imax=6
              if(verbose(isec)==1) imax=8
              if(verbose(isec)==2) imax=15
C             -----
              if(j==1) then
                read(BufNuc(isec,iB,20),*) i1 ! ktyp
                read(BufNuc(isec,iB,18),*) i2 ! no
                write(unit_number,'(a)'   ,advance='no') "t:"
                write(unit_number,'(x,i0)',advance='no') i1 ! ktyp
                if(verbose(isec) >= 1) then
                cTMP=trim(BufNuc(isec,iB,19))//"(Pinned)" ! pName(ktyp)
                write(unit_number,'(x,a)' ,advance='no') trim(cTMP)
     &
                endif
                if(verbose(isec) >= 2) then
                write(unit_number,'(x,i0)',advance='no') i2 ! no
                endif
                write(unit_number,*) ! new line
              endif
C             -----
              do i=1, imax
                if(j==2 .and. i==4) then
                  ! Pin the second instance of time information (i=4).
                  read(BufNuc(isec,iB, i),*) time
                  read(BufNuc(isec,iB,20),*) kf
                  time=time+PinnedTime(kf)
                  write(unit_number,fmtW(isec),advance='no') time
                else
                  ! Normal write-out.
                  write(unit_number,'(a)',advance='no')
     &              trim(BufNuc(isec,iB,i))
                endif
              enddo
              write(unit_number,*) ! new line
            enddo
C         ~~~~~~~~~~
          elseif(trim(format(isec))=="pict") then
            do j=1,2 ! Write out two lines. Add pinned time to the second line.
              imax=8
              if(verbose(isec)==1) imax=9
              if(verbose(isec)==2) imax=18
              do i=1, imax
                if(j==2 .and. i==8) then
                  ! Pin the second instance of time information (i=8).
                  read(BufNuc(isec,iB,i),*) time
                  read(BufNuc(isec,iB,1),*) kf
                  time=time+PinnedTime(kf)
                  write(unit_number,fmtW(isec),advance='no') time
                else
                  ! Normal write-out.
                  write(unit_number,'(a)',advance='no')
     &              trim(BufNuc(isec,iB,i))
                endif
              enddo
              write(unit_number,*) ! new line
            enddo
          endif
C         ~~~~~~~~~~
          if(trim(format(isec))=="t4d") then
            continue
C         ~~~~~~~~~~
          elseif(trim(format(isec))=="pict") then
            write(unit_number,'(i2)') -1
          endif
C         ~~~~~~~~~~
          lastLineWasEndTrack(isec)=.false.
          nBufNuc(isec)=0 ! initialize buffered nuclei
        enddo
C     ------------------------------------------------------------------
      else if(type=="end track") then
        if(.not. lastLineWasEndTrack(isec)) then
C         ~~~~~~~~~~
          if(trim(format(isec))=="t4d") then
            continue
C         ~~~~~~~~~~
          elseif(trim(format(isec))=="pict") then
            write(unit_number,'(i2)') -1
          endif
C         ~~~~~~~~~~
          lastLineWasEndTrack(isec)=.true.
          n_line_of_this_track(isec)=0
        endif
        buff_exist(isec)=.false.
C     ------------------------------------------------------------------
      else if(type=="new source") then
C       This might be a source without any tracks, so the source number is not written here. Wait until 'write_start_point'.
        n_line_of_this_source(isec)=0
C     ------------------------------------------------------------------
      else if(type=="end source") then
        if(n_line_of_this_source(isec) > 0) then
C         ~~~~~~~~~~
          if(trim(format(isec))=="t4d") then
            continue
C         ~~~~~~~~~~
          elseif(trim(format(isec))=="pict") then
            write(unit_number,'(a)') "9 -1"
          endif
C         ~~~~~~~~~~
          n_WrittenHistories(isec) = n_WrittenHistories(isec) + 1
        endif
        lastLineWasEndTrack(isec)=.true.
C     ------------------------------------------------------------------
      else
        print*,"error: line"
        stop
      endif
      end subroutine


************************************************************************
      logical function isStop(ktyp,ncol)
************************************************************************
      if(isNucleus(ktyp)) then
        isStop=(ncol==11) ! 11 : termination by energy cut-off
      else
        isStop=(ncol==11) ! tentative. There seem to be moments to pin even when ncol != 11.
      endif
      return
      end function


************************************************************************
      logical function isBoundary(ncol)
************************************************************************
*             10 : geometry boundary crossing                          *
*             12 : termination by escape or leakage                    *
      isBoundary = ncol==10 .or. ncol==12
      return
      end function


************************************************************************
      logical function isSkip(ktyp,ncol,jcoll,kcoll,nclsts)
************************************************************************
      logical isImportranReaction
      isImportranReaction = (nclsts > 1)
      if(ktyp==2112) then
        isImportranReaction = isImportranReaction .and.
     &  .not.(ncol==14 .and.  jcoll==10
     &                 .and. (kcoll==3 .or. kcoll==4)
     &                 .and.  nclsts==2)
      endif
      isSkip = (.not. isStop(ktyp,ncol)) .and.
     &         (.not. isBoundary(ncol)) .and.
     &         (.not. isImportranReaction) ! not increasing particles
      return
      end function


************************************************************************
      logical function isPinnedParticle(kf)
************************************************************************

      isPinnedParticle=.false.
C     ------------------------
      if(PinnedOption(isec)==1) then ! Display unstable nuclei.
        if(nucleiInHalfLife(kf)) isPinnedParticle=.true.

C     ------------------------
      else if(PinnedOption(isec)==2) then ! Display specified particles (by PinnedParticles).
        do i=1, size(PinnedParticles)
          if(kf==PinnedParticles(isec,i)) then
            isPinnedParticle=.true.
            return
          endif
        enddo

C     ------------------------
      endif

      return
      end function


************************************************************************
      logical function isNucleus(kf)
************************************************************************
C     Z*1000000+A
      isNucleus = 1*1000000+2 <= kf .and. kf <= 119*1000000+1
      return
      end function


************************************************************************
      logical function isPartContained(type,kf)
************************************************************************
      character*(*) type
      integer kf
      integer itmp
      include 'param01.inc'

      isPartContained=.false.

      if(type=="Acc") then
        do i=1,mxpart
          itmp=acceptedParticles(isec,i)
          if(itmp == -1) return
          if(kf==itmp) then
            isPartContained=.true.
            return
          endif
        enddo
      endif

      if(type=="Rej") then
        do i=1,mxpart
          itmp=rejectedParticles(isec,i)
          if(itmp == -1) return
          if(kf==itmp) then
            isPartContained=.true.
            return
          endif
        enddo
      endif

      return
      end function


************************************************************************
      logical function nucleiInHalfLife(kf)
************************************************************************
      if(.not. isNucleus(kf)) then
        nucleiInHalfLife=.false.
        return
      endif
      iZ=kf / 1000000
      iA=mod(kf, 1000000)
      nucleiInHalfLife=(HalflifeMinSEC(isec) < halflife(iZ,iA) .and.
     &                           halflife(iZ,iA) < HalflifeMaxSEC(isec))
      return
      end function


************************************************************************
      double precision function halflife(Z,A)
*     Unit: second
************************************************************************
      integer Z,A
C Mass Excess from BNL (Nucl Wallet Cards 2010)
      INTEGER ATB(3070),ZTB(3070)
      REAL*8 MEXCES(3070)
      COMMON/WC/ATB,ZTB,MEXCES
      REAL*8 DUMMY(3070)
      COMMON/DUM/DUMMY
C     ------------
      if     (Z < 27) then; istart=1;    !  1<=Z<27
      else if(Z < 43) then; istart=501;  ! 27<=Z<43
      else if(Z < 57) then; istart=1001; ! 43<=Z<57
      else if(Z < 71) then; istart=1501; ! 57<=Z<71
      else if(Z < 85) then; istart=2001; ! 71<=Z<85
      else                ; istart=2501; ! 85<=Z
      endif

C     ------------
      do i=istart,3070
        if(Z==ZTB(i)) then
          if(A==ATB(i)) then
            halflife=DUMMY(i)
            return
          endif
        endif
      enddo

C     ------------
      if(Z==73 .and. A==178) then
        halflife=2.36d0*3600d0
        return
      endif

C     ------------
      if(Z==81 .and. A==190) then
        halflife=2.6d0*60d0
        return
      endif

C     ------------
      print*,"Warning: double precision function halflife(Z,A)",Z ,A
      return
      end function

************************************************************************
      integer function UnitNumber(isec)
************************************************************************
      implicit none
      integer,intent(in) :: isec
      integer ipomp,npomp
      common /ipomp0/ipomp,npomp
!$OMP THREADPRIVATE(/ipomp0/)

      if(nt4domp.gt.0)then
       UnitNumber=800+ipomp+npomp*(isec-1)
      else
       UnitNumber=2000+isec
      endif
      return
      end function

************************************************************************
      subroutine opent4dompfile(unit_number,isec)
************************************************************************
      implicit none
      integer,intent(in) :: unit_number,isec
      integer :: ipomp,npomp
      common /ipomp0/ipomp,npomp
!$OMP THREADPRIVATE(/ipomp0/)
      character(100) :: filnm
      character(8) chipomp,chme
      integer :: iord

      if(nt4domp.gt.0)then
       iord=aint(log10(real(ipomp+1+npomp*(isec-1)))) + 1
       write(chipomp,'(i8.8)') ipomp+1+npomp*(isec-1)
       filnm='ompt4dtmp-'//chipomp(9-iord:8)
!$OMP CRITICAL (ompfile_crit) ! To avoid Intel Compiler bug 20220523
       open(unit_number,file=filnm,status='unknown')
!$OMP END CRITICAL (ompfile_crit)
      endif
      return
      end subroutine opent4dompfile

************************************************************************
      subroutine cpt4dompfile(unit_number,isec,unit0)
************************************************************************
      implicit none
      integer,intent(in) :: unit_number,isec,unit0
      character(200)dummy
      integer :: ios

      if(nt4domp.gt.0)then
       rewind(unit_number)
       do
        read(unit_number,'(a)',iostat=ios)dummy
        if(ios.lt.0)exit
        write(unit0,'(a)')trim(dummy)
       enddo
       close(unit_number,status='DELETE')
      endif
      return
      end subroutine cpt4dompfile

************************************************************************
      subroutine chkt4domp()
************************************************************************
!$    use omp_lib
      implicit none
      if(nt4domp.eq.-1)then
       nt4domp=0
!$OMP PARALLEL
!$OMP SINGLE
!$    nt4domp=OMP_GET_NUM_THREADS()
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
       if(nt4domp.eq.1)nt4domp=0
      endif
      return
      end subroutine chkt4domp

************************************************************************
      double precision function random_T4D()
************************************************************************
      implicit none
      ! Linear congruential method
      integer k, t
      integer a, m, q, r
      parameter (a = 16807, m = 2147483647, q = 127773, r = 2836)
      k = seed_T4D / q
      t = a * (seed_T4D - k * q) - r * k
      if (t .gt. 0) then
         seed_T4D = t
      else
         seed_T4D = t + m
      endif
      random_T4D = dble(seed_T4D) / dble(m)
      return
      end

************************************************************************
      end module t4dtrack_mod
************************************************************************
