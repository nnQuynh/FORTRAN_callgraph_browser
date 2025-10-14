************************************************************************
      module liboutmod
*     created by N.Furutachi on 2022/12/28                             *
************************************************************************
      implicit real*8 (a-h,o-z)
      integer :: mxlibnuc, mxmt
      parameter( mxlibnuc = 100 )
      parameter( mxmt = 100, mxmtg = 500 )
      parameter( pi=3.141592653589d0 )

      integer :: ilibtype, icsout, iddxout, igamout, iinfout
      integer :: nlibnuc, nmtout, mxtrial, iapair
      integer :: ilibpart, iemittype
      integer :: iplotlib, iangdisp, icompare, ilibepsout ! frtati 2023/12/07
      integer :: nepx, napx, nepg, napg, nepi
      integer :: ilibnuc(mxlibnuc)
      integer :: mtout(mxmt)
      real*8 :: eneinc, eincsmin, eincsmax, elibreadmax, fcmevpu

      integer :: isettype, isetnuc, isettrial

      integer,allocatable :: iddx(:,:,:), iddxcp(:,:,:)
      integer,allocatable :: iddxg(:,:,:)
      integer,allocatable :: ipmarr(:)
      real*8,allocatable :: emshx(:), amshx(:)
      real*8,allocatable :: emshg(:), amshg(:), emshi(:)
      real*8,allocatable :: ddx(:,:,:), ddxall(:,:), ddxpart(:,:)
      real*8,allocatable :: ddxg(:,:), ddxc(:,:)

      real*8,allocatable :: threshold_mt(:)
      real*8,allocatable :: sig_part(:), sig_partcp(:)
      real*8,allocatable :: sig_partg(:), discgam(:)

      real*8 :: angtem(2)
      character*100 :: cfilen, cfiletem, cfilesig, cfileddx, cfilegam
      character*100 :: cfileddxd
      character*20 :: cangelkey, cangcom
      character*40 :: cangeltitle ! frtati 2023/12/13
      character*200 :: cangelsig, cangelddx, cangelgam

*-----------------------------------------------------------------------
      parameter(icsumx=30)
      integer :: icsu
      dimension lschn(icsumx), ischn(icsumx)
      character schan(icsumx)*8

*-----------------------------------------------------------------------
      contains

************************************************************************
      subroutine allocate_liboutarray(mxnepx,mxangp,mxnepg,mxnepi)
************************************************************************
      mxenep = max(mxnepx,mxnepg)

      allocate(emshx(mxenep+1),amshx(mxangp+1))
      allocate(emshg(mxenep+1),amshg(2),emshi(mxnepi+1))
      allocate(ipmarr(max(mxenep,mxangp)*2+3))
      allocate(iddx(mxangp,mxenep,mxmt),iddxcp(mxangp,mxenep,mxmt))
      allocate(iddxg(1,mxenep,mxmtg))
      allocate(ddx(mxangp,mxenep,mxmt))
      allocate(ddxall(mxangp,mxenep),ddxpart(mxangp,mxenep))
      allocate(ddxg(1,mxenep))
      if( icompare.eq.1 ) allocate(ddxc(mxangp,mxenep))

      allocate(threshold_mt(mxmt))
      allocate(sig_part(mxmt),sig_partcp(mxmt))
      allocate(sig_partg(mxmtg),discgam(mxmtg))

      emshx = 0.d0
      amshx = 0.d0
      emshg = 0.d0
      amshg = 0.d0
      emshi = 0.d0

      iddx = 0
      iddxcp = 0
      iddxg = 0
      ddx = 0.d0
      ddxall = 0.d0
      ddxpart = 0.d0
      ddxg = 0.d0

      return
      end subroutine allocate_liboutarray

************************************************************************
      subroutine acewrite
************************************************************************
      use MEMBANKMOD
      use GGMARRAYMOD
      use GGMBANKMOD
      use moddas_ggs

*-----------------------------------------------------------------------
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------
      real*8 :: yd(2)
      character:: hs*10
      character:: ha*80,hc*181,hd*10,hf*64,hk*70,hm*10,hp*10,hr*70,ht*10
      character:: emitpartname*8
      character:: fname_angel*200 ! frtati 2023/12/07

*-----------------------------------------------------------------------
! set libout parameters using settings common for icntl=1
      call setlibparam

*-----------------------------------------------------------------------
! output files
      ioa = 101
      ioc = 102
      iod = 103
      iog = 104

      if( iemittype.eq.1 ) emitpartname = "neutron"
      if( iemittype.eq.9 ) emitpartname = "proton"
      if( iemittype.eq.31 ) emitpartname = "deuteron"
      if( iemittype.eq.32 ) emitpartname = "triton"
      if( iemittype.eq.33 ) emitpartname = "3He"
      if( iemittype.eq.34 ) emitpartname = "alpha"

      if( iinfout.ne.0 ) then
        open(ioa,file="libout.out",status='replace')
        call echolibout(ioa)
        write(ioa,'("#",78("-"))')
        write(ioa,'("# Infomation of nuclear data library")')
      end if
      if( icsout.ne.0 ) then
        open(ioc,file=trim(cfilesig),status='replace')
        write(ioc,'("#",78("-"))')
        write(ioc,'("# Cross-sections (b)")')
      end if
      if( iddxout.ne.0 ) then
        open(iod,file=trim(cfileddx),status='replace')
        if( iddxout.eq.1 ) then
          write(iod,'("#",78("-"))')
          write(iod,'("# Differential-cross sections ",
     &          "with respect to energy (mb/MeV)")')
        else if( iddxout.eq.2 ) then
          write(iod,'("#",78("-"))')
          write(iod,'("# Double-differential ",
     &          "cross-sections (mb/sr/MeV)")')
        end if
        write(iod,'("# E_in =",ES10.2," (MeV/u)")') fcmevpu*eneinc
        write(iod,'("# emission particle: ",A)') emitpartname
      end if
      if( igamout.ne.0 ) then
        open(iog,file=trim(cfilegam),status='replace')
        write(iog,'("#",78("-"))')
        write(iog,'("# Gamma-spectrum (mb/MeV)")')
        write(iog,'("# E_in =",ES10.2," (MeV/u)")') fcmevpu*eneinc
      end if


*-----------------------------------------------------------------------
! initial setteing

      erg = eneinc

! for bidirectional search from previous energy/angle point
      iepos = 1
      iapos = 1
      do i = 1, max(napx,nepx)*2+3, 2
        ipmarr(i) = i/2
      end do
      do i = 2, max(napx,nepx)*2+2, 2
        ipmarr(i) = -i/2
      end do

*-----------------------------------------------------------------------
      nucloop: do inuc = 1, nlibnuc ! output for each nucleus

*-----------------------------------------------------------------------
      mn = 1
      do km = 1, mix
        if( km .ge. jmd(1+mn+1) ) mn = mn + 1
        if( iza(km).eq.ilibnuc(inuc) ) then
          if( lme(ilibtype,km).eq.0 ) then
            write(ErrCha,'(I6," : Library is not loaded")')ilibnuc(inuc)
            ErrID = 'L:183/R:acewrite/F:liboutmod.f'
            call ErrWrite(ErrID,ErrCha)
            cycle nucloop
          end if
          iex = lme(ilibtype,km)
          exit
        end if
      end do
      if( km.eq.mix+1 ) then
        write(ErrCha,'(I6," : Library is not loaded")') ilibnuc(inuc)
        ErrID = 'L:193/R:acewrite/F:liboutmod.f'
        call ErrWrite(ErrID,ErrCha)
        cycle nucloop
      end if
      elibmx = xss(jxs(1,iex)+nxs(3,iex)-1)
      if( erg.gt.elibmx ) then
        write(ErrCha,'("Incident energy of",ES10.2," is greater than ",
     &  "maximum energy of library:",ES10.2)') erg, elibmx
        ErrID = 'L:201/R:acewrite/F:liboutmod.f'
        call ErrWrite(ErrID,ErrCha)
        cycle nucloop
      end if

      if( ilibtype.eq.9 ) then
        itypein = 1
      else if( ilibtype.eq.31 ) then
        itypein = 15
      else if( ilibtype.eq.34 ) then
        itypein = 18
      end if

cfrtati 2023/12/13
      icsmtout = 0
      do i = 1, nmtout
        if( mtout(i).gt.0 .and. mtout(i).lt.4 ) then
          icsmtout = 1
        end if
      end do

! save iex and set dummy cell and temperature
      iex_save = iex ! save to avoid change of iex in xstneu
      icli_dummy = 1
      tmei_dummy = 0.d0 ! frtati 2023/12/13

*--------
      if( icsout.ne.0 ) then
        write(ioc,'("#",78("-"))')
        write(ioc,'("# ZA = ",I6)') ilibnuc(inuc)
        write(ioc,'("#",78("-"))')
      end if
      if( iddxout.ne.0 ) then
        write(iod,'("#",78("-"))')
        write(iod,'("# ZA = ",I6)') ilibnuc(inuc)
        write(iod,'("#",78("-"))')
      end if
      if( igamout.ne.0 ) then
        write(iog,'("#",78("-"))')
        write(iog,'("# ZA = ",I6)') ilibnuc(inuc)
        write(iog,'("#",78("-"))')
      end if

      if( iinfout.ne.0 ) then
        hc = ' '
        call zaid(2,hc,ixc(1,iex))
        write(ioa,'("#",78("-"))')
        write(ioa,'("  Library: ",A)') trim(hc)
        cangcom = trim(hc)
        do 10 j = 11, 181
          k = mod(ixc((j+1)/3,iex)/256**mod(j+1,3),256)
          if( k .eq. 0 ) goto 20
   10     hc(j:j) = char(k)
   20     iu = 10
        do 30 j = 2, 7
   30     call nxtsym(hc,' ',iu+1,it,iu,0)
          ht = hc(it:iu)
        write(ioa,'("  Table length: ",A10)') ht
        call nxtsym(hc,' ',11,it,iu,0)
        call nxtsym(hc,' ',iu+1,it,iu,0)
        hf = hc(it:iu)
        write(ioa,'("  Tables from file: ",A,/)') hf
        write(ioa,'("#",78("-"))')

        write(ioa,'("# Maximum energy of library (MeV)")')
        write(ioa,'(ES10.2)') elibmx
      end if

*--------
*   jxs(1): ESZ=location of energy table
*   nxs(3): NES=number of energies
      nesz = jxs(1,iex)
      nnes = nxs(3,iex)
      if( iinfout.ne.0 ) then
        write(ioa,'("# Number of incident energy points")')
        write(ioa,'(I10)') nnes
      end if

      if( icsout.ne.0 .and. icsmtout.eq.0 ) then ! frtati 2023/12/13
        if( abs(icsout).eq.1 ) then
          write(ioc,'("# E_in(MeV/u)  total(b)")')
          if( iplotlib.eq.1 ) then
            cangelkey = "total"
            call setlibangel(ioc,1)
          end if
        else if( abs(icsout).eq.2 ) then
         write(ioc,'("# E_in(MeV/u) total(b) abs.(b) ela.(b) non.(b)")')
          if( iplotlib.eq.1 ) call setlibangel(ioc,2)
        end if
        if( icsout.lt.0.d0 ) then
          do i = 1, nnes
            ii0 = nesz+i-1
            ii1 = nesz+nnes+i-1
            ii2 = nesz+2*nnes+i-1
            ii3 = nesz+3*nnes+i-1
            if( xss(ii0).ge.eincsmin .and. xss(ii0).le.eincsmax ) then
              if( abs(icsout).eq.1 ) then
                write(ioc,'(2ES11.3)') fcmevpu*xss(ii0), xss(ii1)
              else if( abs(icsout).eq.2 ) then
                write(ioc,'(5ES11.3)') fcmevpu*xss(ii0), xss(ii1),
     &          xss(ii2), xss(ii3), xss(ii1)-xss(ii3)
              end if
            end if
          end do
        else
cfrtati 2023/12/13
          ieit = 0
          do j = 1, nepi + 1
            do i = ieit, nnes - 1
              ii0 = nesz+i-1
              if( emshi(j).ge.xss(ii0) .and. emshi(j).le.xss(ii0+1))then
                ieit = i
                exit
              end if
            end do
            if( i.eq.nnes ) cycle
            ii1 = nesz+nnes+ieit-1
            ii2 = nesz+2*nnes+ieit-1
            ii3 = nesz+3*nnes+ieit-1
            if( emshi(j).eq.xss(ii0) .or. emshi(j).eq.xss(ii0+1) ) then
              tot = xss(ii1)
              abo = xss(ii2)
              ela = xss(ii3)
            else
              frc = (emshi(j)-xss(ii0))/(xss(ii0+1)-xss(ii0))
              tot = xss(ii1) + frc*(xss(ii1+1)-xss(ii1))
              abo = xss(ii2) + frc*(xss(ii2+1)-xss(ii2))
              ela = xss(ii3) + frc*(xss(ii3+1)-xss(ii3))
            end if
            rea = tot - ela
            if( icsout.eq.1 ) then
              write(ioc,'(2ES11.3)') fcmevpu*emshi(j), tot
            else if( icsout.eq.2 ) then
              write(ioc,'(5ES11.3)') fcmevpu*emshi(j), tot,
     &        abo, ela, rea
            end if
          end do
        end if
        write(ioc,*)
        write(ioc,*)
      end if

*--------
*   jxs(3): MTR=location of MT array
*   nxs(4): NTR=number of reactions excluding elastic
      nmtr = jxs(3,iex)
      nntr = nxs(4,iex)

      if( iinfout.ge.1 ) then
        write(ioa,'("# MT")')
        write(ioa,'(10I5)') ( nint(xss(nmtr+j-1)),j=1,nntr )
      end if

*--------
*   jxs(5): TYR=location of reaction type array
*   nxs(4): NTR=number of reactions excluding elastic
*   note: The sign indicates the system for scattering: negative = CM system; positive = LAB system.
      ntyr = jxs(5,iex)
      nnr = nxs(4,iex)

      if( iinfout.ge.1 ) then
        write(ioa,'("# Neutron release number of MT")')
        write(ioa,'(10I8)') ( int(xss(ntyr+j-1)),j=1,nntr )
      end if

*--------
*   jxs(6): LSIG=location of table of cross-section locators
      nlsig = jxs(6,iex)
      nsig = jxs(7,iex)

      if( icsout.ne.0 ) then
        do i = 1, nntr
          iloca = int( xss(nlsig+i-1) )
          ii0 = int(xss(nmtr+i-1))
          threshold_mt(i) = xss( nesz+int(xss(nsig+iloca-1))+1-2 )
          do j = 1, nmtout
            if ( ii0.eq.mtout(j) ) then
              write(ioc,'("# MT = ",I5)') mtout(j)
              write(ioc,'("# E_in(MeV/u)   CS(b)")')
              if( iplotlib.eq.1 ) then
                write(cangelkey,'("MT=",I3)') ii0
                call setlibangel(ioc,1)
              end if
              if( icsout.lt.0.d0 ) then
                do k = 1, int( xss(nsig+iloca) )
                  if( xss( nesz+int(xss(nsig+iloca-1))+k-2 )
     &               .ge.eincsmin .and.
     &                xss( nesz+int(xss(nsig+iloca-1))+k-2 )
     &               .le.eincsmax ) then
                    write(ioc,'(2ES12.4)')
     &              fcmevpu*xss( nesz+int(xss(nsig+iloca-1))+k-2 ),
     &              xss( nsig+iloca+k )
                  end if
                end do
              else
cfrtati 2023/12/13
                ieit = 0
                do l = 1, nepi + 1
                  if( emshi(l).lt.threshold_mt(i) ) cycle
                  do k = ieit, int( xss(nsig+iloca) ) - 1
                    ii1 = nesz+int(xss(nsig+iloca-1))+k-2
                    if( emshi(l).ge.xss(ii1) .and.
     &                  emshi(l).le.xss(ii1+1) ) then
                      ieit = k
                      exit
                    end if
                  end do
                  if( k.eq.int(xss(nsig+iloca)) ) cycle
                  ii2 = nsig+iloca+ieit
                  if(emshi(l).eq.xss(ii1).or.emshi(l).eq.xss(ii1+1))then
                    csmt = xss(ii2)
                  else
                    frc = (emshi(l)-xss(ii1))/(xss(ii1+1)-xss(ii1))
                    csmt = xss(ii2)+frc*(xss(ii2+1)-xss(ii2))
                  end if
                  write(ioc,'(2ES12.4)') fcmevpu*emshi(l), csmt
                end do
              end if
              write(ioc,*)
              write(ioc,*)
              exit
            end if
          end do
        end do
      end if

cfrtati 2023/12/13
      do j = 1, nmtout
        if( mtout(j).lt.4 .or. mtout(j).gt.1000 ) then
          if( mtout(j).lt.0 ) then
            write(ioc,'("# special number = ",I3)') mtout(j)
            if( mtout(j).gt.-4.or.mtout(j).eq.-5.or.mtout(j).eq.-6 )then
              write(ioc,'("# E_in(MeV/u)   CS(b)")')
              cangeltitle = "cross-section [b]"
            else if( mtout(j).eq.-4 ) then
              write(ioc,'("# E_in(MeV/u)   heating(MeV/collision)")')
              cangeltitle = "average heating [MeV/collision]"
            else if( mtout(j).eq.-7 ) then
              write(ioc,'("# E_in(MeV/u)   number(n/fission)")')
              cangeltitle = "neutron emission number [n/fission]"
            else if( mtout(j).eq.-8 ) then
              write(ioc,'("# E_in(MeV/u)   heating(MeV/fission)")')
              cangeltitle = "average heating [MeV/fission]"
            end if
            if( iplotlib.eq.1 ) then
              write(cangelkey,'("number=",I3)') mtout(j)
              call setlibangel(ioc,0)
            end if
          else
            write(ioc,'("# MT = ",I5)') mtout(j)
            write(ioc,'("# E_in(MeV/u)   CS(b)")')
            if( iplotlib.eq.1 ) then
              write(cangelkey,'("MT=",I6)') mtout(j)
              call setlibangel(ioc,1)
            end if
          end if
          do l = 1, nepi + 1
            erg = emshi(l)
            if( ilibtype.eq.1 ) then
              call xstneu(0,sigt,sigaa,icli_dummy,erg,tmei_dummy,mn)
              iex = iex_save
              csmt = getxs(mtout(j))
            else if( ilibtype.eq.9 .or. ilibtype.eq.31 .or.
     &               ilibtype.eq.34 ) then
              einn = erg
              call sig_tot(itypein,sigt,sigaa,einn,mn)
              iex = iex_save
              if( mtout(j).eq.1 ) then
                csmt = rtc(5,iex)
              else if( mtout(j).eq.2 ) then
                csmt = rtcel(iex)
              else if( mtout(j).eq.3 ) then
                csmt = rtc(5,iex) - rtcel(iex)
              else if( mtout(j).eq.-4 ) then
                ii3 = jxs(1,iex)+4*nxs(3,iex)+ktc(1,iex)
                csmt = xss(ii3-1) + rtc(1,iex)*( xss(ii3)-xss(ii3-1) )
              end if
            end if
            if( csmt.eq.0.d0 ) cycle
            write(ioc,'(2ES12.4)') fcmevpu*emshi(l), csmt
          end do
          erg = eneinc
          write(ioc,*)
          write(ioc,*)
        end if
      end do

*--------
*   jxs(13): MTRP=location of photon production MT array
*   nxs(6): NTRP=number of photon production reactions excluding elastic
      nmtrp = jxs(13,iex)
      nntrp = nxs(6,iex)
      if( iinfout.ge.1 ) then
        write(ioa,'("# Photon production MT")')
        write(ioa,'(10I8)') ( int(xss(nmtrp+j-1)),j=1,nntrp )
      end if

*--------
*   jxs(12): GPD=location of photon production data
*   nxs(6): NTRP=number of photon production reactions excluding elastic
      ngpd = jxs(12,iex)

      if( icsout.ne.0 .and. igamout.ge.1 ) then
        write(ioc,'("# Total photon production cross section")')
        write(ioc,'("# E_in(MeV/u)   CS(b)")')
        if( icsout.lt.0.d0 ) then
          do i = 1, nnes
            write(ioc,'(2E12.4)') fcmevpu*xss(nesz+i-1), xss(ngpd+i-1)
          end do
        else
          ieit = 1
          do j = 1, nepi
            do i = 1, nnes
              if( xss(nesz+i-1).ge.emshi(j) .and.
     &        xss(nesz+i-1).lt.emshi(j+1) ) then
                write(ioc,'(2ES12.4)') fcmevpu*xss(nesz+i-1),
     &          xss(ngpd+i-1)
                ieit = i + 1
                exit
              end if
            end do
          end do
        end if
        write(ioc,*)
        write(ioc,*)
      end if

*--------
      if( iinfout.ne.0 ) then
        write(ioa,'("# Num. of other emission particles")')
        write(ioa,'(I4)') nxs(7,iex)
        if( nxs(7,iex).gt.0 ) then
          write(ioa,'("# Emission particle type other than ",
     &    "incident particle type")')
          write(ioa,'(100I4)')
     &    (nint(xss(jxs(30,iex)+i-1)),i=1,nxs(7,iex))
        end if
        if( iemittype.ne.ilibtype ) then
          do i = 1, nxs(7,iex)
            if( nint(xss(jxs(30,iex)+i-1)).eq.iemittype ) exit
          end do
          if( i.eq.nxs(7,iex)+1 ) then
            write(ErrCha,'("Emission particle type is ",
     &      "not foud for ",I6)') ilibnuc(inuc)
            ErrID = 'L:465/R:acewrite/F:liboutmod.f'
            call ErrWrite(ErrID,ErrCha)
            cycle nucloop
          end if
        end if
      end if

*-----------------------------------------------------------------------
! set rtc and ktc array for nucleus in material mn
      if( ilibtype.eq.1 ) then
        call xstneu(0,sigt,sigaa,icli_dummy,erg,tmei_dummy,mn)
      else if( ilibtype.eq.9 .or. ilibtype.eq.31 .or.
     &         ilibtype.eq.34 ) then
        einn = erg
        call sig_tot(itypein,sigt,sigaa,einn,mn)
      end if
      iex = iex_save

*-----------------------------------------------------------------------
*  ddx calculation start

! ixre and ntyn are used in xstcas
*-------- neutron production ddx
      if ( iddxout.ge.1 ) then
        iddx = 0
        iddxcp = 0
        ddx = 0.d0
        ddxall = 0.d0
        ddxpart = 0.d0
        sig_part = 0.d0
        sig_partcp = 0.d0
        sig_nonela = 0.d0 ! frtati 2022/02/01

        esave = erg

! sampling of particle type same with incident particle type
! notice: only in this process, MT is randomly selected.
! in other process, partial cross-section of each MT is used
        if( ilibtype.eq.iemittype ) then
          ipt = iemittype

          do ixre = 1, nxs(5,iex)
            is = jxs(7,iex)+nint(xss(jxs(6,iex)+ixre-1))
            ic = ktc(1,iex)+1-nint(xss(is-1))
            if(ic.lt.1.or.ic.gt.nint(xss(is)).or.
     &         ic.eq.nint(xss(is)).and.rtc(1,iex).ne.0.) cycle
            sig_part(ixre) = xss(ic+is)
            if ( rtc(1,iex).ne.0.d0 ) then
              sig_part(ixre) = sig_part(ixre)+rtc(1,iex)
     &                        *(xss(ic+is+1)-xss(ic+is))
            end if
            nswr = abs(nint(xss(jxs(5,iex)+ixre-1)))
            if ( nswr.gt.100 ) then
              l = jxs(11,iex)+nswr-101
              cmult = acefcn(l,erg,ln)
              sig_part(ixre) = cmult * sig_part(ixre)
            end if
            sig_nonela = sig_nonela + sig_part(ixre)
          end do

          trialloop: do itrial = 1, mxtrial
            r = 0.d0
            rand = rang()
            do ixre = 1, nxs(5,iex)
              r = r + sig_part(ixre) / sig_nonela
              if( rand.le.r ) exit
            end do
            ntyn = nint(xss(jxs(5,iex)+ixre-1)) ! used in xstcas through common
            nswr = abs(ntyn)
            qwr = xss(jxs(4,iex)+ixre-1) ! used
            mtpwr = nint(xss(jxs(3,iex)+ixre-1))
            iawr = jxs(9,iex)
            kawr = nint(xss(jxs(8,iex)+ixre))
            idwr = jxs(11,iex)
            kdwr = nint(xss(jxs(10,iex)+ixre-1))
            if ( abs(ntyn).gt.100 ) then ! neutron multiplicity is given by yield data
              nswr = 1
            end if
            do ins = 1, nswr
              call xstcas(1,ilibtype,qwr,iawr,kawr,idwr,kdwr)
              do iebin = 1, nepx*2
                iet = iepos + ipmarr(iebin)
                if ( iet.lt.1 .or. iet.gt.nepx ) cycle
                if ( colout(1,1).gt.emshx(iet) .and.
     &               colout(1,1).le.emshx(iet+1) ) then
                  do iabin = 1, napx*2
                    iat = iapos + ipmarr(iabin)
                    if ( iat.lt.1 .or. iat.gt.napx ) cycle
                    if ( colout(2,1).gt.amshx(iat+1) .and.
     &                   colout(2,1).le.amshx(iat) ) then
                      iddx(iat,iet,ixre) =
     &                iddx(iat,iet,ixre) + 1
                      iepos = iet
                      iapos = iat
                      cycle trialloop
                    end if
                  end do ! iabin
                end if
              end do ! iebin
            end do ! ins
          end do trialloop
        end if

! sampling of particle type different from incident particle type
        do i = 1, nxs(7,iex)
          ip = nint(xss(jxs(30,iex)+i-1))
          if( ip.ne.iemittype ) cycle
          ipt = ip

          is = ixs(1,i,iex)+1
          ic = ktc(1,iex)+1-nint(xss(is-1))
          if(ic.lt.1.or.ic.gt.nint(xss(is)).or.
     &       ic.eq.nint(xss(is)).and.rtc(1,iex).ne.0.) cycle
          t = xss(ic+is)
          if(rtc(1,iex).ne.0.) t = t+rtc(1,iex)*
     &    (xss(ic+is+1)-xss(ic+is))
          if(t.eq.0.) cycle

          ixloopcp: do ixre = 1, nint(xss(jxs(31,iex)+i-1))
            if(nint(xss(jxs(31,iex)+i-1)).gt.1) then
              l = ixs(5,i,iex)+nint(xss(ixs(4,i,iex)+ixre-1))+1
              do jy = 1, 2
                yd(jy) = acefcn(l,xss(jxs(1,iex)
     &          + ktc(1,iex)+jy-2)*1.00000000001d0,ln)
                if(yd(1).eq.0..and.yd(2).eq.0.) cycle ixloopcp
                ix = nint(xss(l-1))
                is = jxs(7,iex)+nint(xss(jxs(6,iex)+ix-1))
                x1 = 0.
                x2 = 0.
              end do
              ic = ktc(1,iex)-nint(xss(is-1))+1
              if(ic.gt.0.and.ic.le.nint(xss(is))) x1 = xss(is+ic)
              ic = ic+1
              if(ic.gt.0.and.ic.le.nint(xss(is))) x2 = xss(is+ic)
              s = x1*yd(1)+rtc(1,iex)*(x2*yd(2)-x1*yd(1))
              sig_partcp(ixre) = s
            else
              sig_partcp(ixre) = t
            end if

            ntyn = nint(xss(ixs(3,i,iex)+ixre-1)) ! used in xstcas through common
            mtpwr = nint(xss(ixs(2,i,iex)+ixre-1))
            iawr = ixs(7,i,iex)
            kawr = nint(xss(ixs(6,i,iex)+ixre-1))
            idwr = ixs(9,i,iex)
            kdwr = nint(xss(ixs(8,i,iex)+ixre-1))

            trialloopcp: do itrial = 1, mxtrial
              call xstcas(1,ilibtype,0.d0,iawr,kawr,idwr,kdwr)
              do iebin = 1, nepx*2
                iet = iepos + ipmarr(iebin)
                if ( iet.lt.1 .or. iet.gt.nepx ) cycle
                if ( colout(1,1).gt.emshx(iet) .and.
     &               colout(1,1).le.emshx(iet+1) ) then
                  do iabin = 1, napx*2
                    iat = iapos + ipmarr(iabin)
                    if ( iat.lt.1 .or. iat.gt.napx ) cycle
                    if ( colout(2,1).gt.amshx(iat+1) .and.
     &                   colout(2,1).le.amshx(iat) ) then
                      iddxcp(iat,iet,ixre) =
     &                iddxcp(iat,iet,ixre) + 1
                      iepos = iet
                      iapos = iat
                      cycle trialloopcp
                    end if
                  end do ! iabin
                end if
              end do ! iebin
            end do trialloopcp
          end do ixloopcp
        end do ! mtloop

        erg = esave

        do ixre = 1, nxs(5,iex)
          do iebin = 1, nepx
            do iabin = 1, napx
              ddx(iabin,iebin,ixre) =
     &        (iddx(iabin,iebin,ixre)*sig_nonela +
     &         iddxcp(iabin,iebin,ixre)*sig_partcp(ixre))
     &        /dble(mxtrial)/(emshx(iebin+1)-emshx(iebin))
     &        /(2.d0*pi*(amshx(iabin)-amshx(iabin+1)))
            end do
          end do
          ddxall(:,:) = ddxall(:,:) + ddx(:,:,ixre)
        end do

        ddxpart = 0.d0
        do ixre = 1, nxs(5,iex)
          if ( ixre.ge.4 .and. ixre.le.30 ) then
            ddxpart(:,:) = ddxpart(:,:) + ddx(:,:,ixre)
          end if
        end do

! partial cross-sections
        if( iinfout.ne.0 ) then
          if( nxs(5,iex).gt.1 .or. ilibtype.eq.iemittype ) then
            if( nxs(5,iex).gt.1 ) then
              write(ioa,'("# Partial cross-sections of MTs")')
            else
              write(ioa,'("# Production cross-section of particle ",
     &        "type same with incident particle")')
            end if
            write(ioa,'("# MT  CS(b)  multiplicity")')
            do ixre = 1, nxs(5,iex)
              nswr = abs(nint(xss(jxs(5,iex)+ixre-1)))
              if ( nswr.gt.100 ) then
                l = jxs(11,iex)+nswr-101
                cmult = acefcn(l,erg,ln)
              else
                cmult = dble(nswr)
              end if
              write(ioa,'(I4,ES11.3,F6.2)')
     &        nint(xss(jxs(3,iex)+ixre-1)), sig_part(ixre), cmult
            end do
          else
            do i = 1, nxs(7,iex)
              ip = nint(xss(jxs(30,iex)+i-1))
              if( ip.eq.1 ) then
                write(ioa,'("# Particle production cross-sction ",
     &          "(multiplicity included)")')
                write(ioa,'("# MT  CS(b)")')
                do ixre = 1, nint(xss(jxs(31,iex)+i-1))
                  write(ioa,'(I4,ES11.3)')
     &            nint(xss(ixs(2,i,iex)+ixre-1)), sig_partcp(ixre)
                end do
              end if
            end do
          end if
          write(ioa,*)
        end if

! ddx(mb)
        if( iddxout.eq.2 ) then ! ddx
          if( iangdisp.eq.0 ) then
            do iabin = 1, napx
              if( iapair.eq.1 .and. mod(iabin,2).eq.0 ) cycle
              angtem(1) = 180.d0/pi*acos(amshx(iabin))
              angtem(2) = 180.d0/pi*acos(amshx(iabin+1))
              write(iod,'("# angle(deg): ",F5.1," - ",F5.1)')
     &        angtem(1), angtem(2)
              write(iod,'("#",A10,2A11)')'emin(MeV)','emax(MeV)','CS'
              if( iplotlib.eq.1 ) then
                cangelkey = trim(emitpartname)
                if( icompare.eq.0 ) then
                  call setlibangel(iod,5)
                else if( icompare.eq.1 ) then
                  call setlibangel(iod,7)
                end if
              end if
              if( icompare.eq.0 ) then
                do iebin = 1, nepx
                  write(iod,'(3ES11.3)') emshx(iebin), emshx(iebin+1),
     &            1.d3*ddxall(iabin,iebin)
                end do
              else if( icompare.eq.1 ) then
                do iebin = 1, nepx
                  write(iod,'(4ES11.3)') emshx(iebin), emshx(iebin+1),
     &            1.d3*ddxall(iabin,iebin), ddxc(iabin,iebin)
                end do
              end if
              write(iod,*)
              write(iod,*)
              if( nmtout.gt.0 ) then
                do i = 1, nmtout
                  do ixre = 1, nxs(5,iex)
                    if( nint(xss(jxs(3,iex)+ixre-1)).eq.mtout(i) ) then
                      write(iod,'("# angle(deg): ",F5.1," - ",F5.1)')
     &                angtem(1), angtem(2)
                      write(iod,'("# MT = ",I5)') mtout(i)
                      write(iod,'("#",A10,2A11)')
     &                'emin(MeV)','emax(MeV)','CS'
                      if( iplotlib.eq.1 ) then
                        write(cangelkey,'("MT=",I3" neutron")')mtout(i)
                        call setlibangel(iod,5)
                      end if
                      do iebin = 1, nepx
                        write(iod,'(3ES11.3)') emshx(iebin),
     &                  emshx(iebin+1), 1.d3*ddx(iabin,iebin,ixre)
                      end do
                      write(iod,*)
                      write(iod,*)
                      exit
                    end if
                  end do
                end do
              end if
            end do
          else if( iangdisp.eq.1 ) then
            iastep = 1
            if( iapair.eq.1 ) iastep = 2
            write(iod,'("#",A10,2A)')'emin(MeV)',' emax(MeV)',
     &      ' (CS(i),i=1,na)'
            if( iplotlib.eq.1 ) then
              cangelkey = trim(emitpartname)
              if( icompare.eq.0 ) then
                call setlibangel(iod,6)
              else if( icompare.eq.1 ) then
                call setlibangel(iod,8)
              end if
            end if
            if( icompare.eq.0 ) then
              do iebin = 1, nepx
                write(iod,'(100ES11.3)') emshx(iebin), emshx(iebin+1),
     &          (1.d3*ddxall(iabin,iebin),iabin=1,napx,iastep)
              end do
            else if( icompare.eq.1 ) then
              do iebin = 1, nepx
                write(iod,'(200ES11.3)') emshx(iebin), emshx(iebin+1),
     &          (1.d3*ddxall(iabin,iebin),iabin=1,napx,iastep),
     &          (ddxc(iabin,iebin),iabin=1,napx/iastep+1-mod(iastep,2))
              end do
            end if
            write(iod,*)
            write(iod,*)
            if( nmtout.gt.0 ) then
              do i = 1, nmtout
                do ixre = 1, nxs(5,iex)
                  if( nint(xss(jxs(3,iex)+ixre-1)).eq.mtout(i) ) then
                    write(iod,'("# MT = ",I5)') mtout(i)
                    write(iod,'("#",A10,2A)')
     &              'emin(MeV)','emax(MeV)','(CS(i),i=1,na)'
                    if( iplotlib.eq.1 ) then
                      write(cangelkey,'("MT=",I3" neutron")')mtout(i)
                      call setlibangel(iod,6)
                    end if
                    do iebin = 1, nepx
                      write(iod,'(100ES11.3)')
     &                emshx(iebin), emshx(iebin+1),
     &                (1.d3*ddx(iabin,iebin,ixre),iabin=1,napx,iastep)
                    end do
                    write(iod,*)
                    write(iod,*)
                    exit
                  end if
                end do
              end do
            end if
          end if
        else if( iddxout.eq.1 ) then! energy distribution
          write(iod,'("#",A10,2A11)')'emin(MeV)','emax(MeV)','CS'
          if( iplotlib.eq.1 ) then
            cangelkey = trim(emitpartname)
            call setlibangel(iod,3)
          end if
          do iebin = 1, nepx
            write(iod,'(3ES11.3)') emshx(iebin), emshx(iebin+1),
     &      1.d3*ddxall(1,iebin)*(4.d0*pi)
          end do
          write(iod,*)
          write(iod,*)
          if( nmtout.gt.0 ) then
            do i = 1, nmtout
              do ixre = 1, nxs(5,iex)
                if( nint(xss(jxs(3,iex)+ixre-1)).eq.mtout(i) ) then
                  write(iod,'("# MT = ",I5)') mtout(i)
                  write(iod,'("#",A10,2A11)')
     &            'emin(MeV)','emax(MeV)','CS'
                  if( iplotlib.eq.1 ) then
                    write(cangelkey,'("MT=",I3," neutron")')mtout(i)
                    call setlibangel(iod,3)
                  end if
                  do iebin = 1, nepx
                    write(iod,'(3ES11.3)') emshx(iebin),
     &              emshx(iebin+1), 1.d3*ddx(1,iebin,ixre)*(4.d0*pi)
                  end do
                  write(iod,*)
                  write(iod,*)
                end if
              end do
            end do
          end if
        else if( iddxout.eq.4 ) then ! partial ddx
          do iabin = 1, napx
            if( iapair.eq.1 .and. mod(iabin,2).eq.0 ) cycle
            write(iod,'("# angle(deg): ",F5.1," - ",F5.1)')
     &      180.d0/pi*acos(amshx(iabin)),
     &      180.d0/pi*acos(amshx(iabin+1))
            write(iod,'("#",A10,2A11)')'emin(MeV)','emax(MeV)','CS'
            do iebin = 1, nepx
              write(iod,'(3ES11.3)') emshx(iebin), emshx(iebin+1),
     &        1.d3*ddxpart(iabin,iebin)
            end do
            write(iod,*)
            write(iod,*)
          end do
        else if( iddxout.eq.3 ) then
          write(iod,'("#",A10,2A11)')'emin(MeV)','emax(MeV)','CS'
          do iebin = 1, nepx
            write(iod,'(3ES11.3)') emshx(iebin), emshx(iebin+1),
     &      1.d3*ddxpart(1,iebin)*(4.d0*pi)
          end do
          write(iod,*)
          write(iod,*)
        end if

      end if ! iddxout

*-------- photon production ddx
      if ( igamout.ge.1 ) then
        iddxg = 0
        ddxg = 0.d0

        ipt = 2 ! gamma emission
        esave = erg

        mtloopg: do ixre = 1, nxs(6,iex)
          l = jxs(15,iex)+nint(xss(jxs(14,iex)+ixre-1))+1
          if( nint(xss(l-2)).eq.13 ) then ! mf=13
            ic = ktc(1,iex)-nint(xss(l-1))+1
            if(ic.lt.1.or.ic.gt.nint(xss(l)).or.
     &         ic.eq.nint(xss(l)).and.rtc(1,iex).ne.0.) cycle
            t = xss(ic+l)+rtc(1,iex)*(xss(ic+l+1)-xss(ic+l))
            if( t.eq.0.d0 ) cycle
            sig_partg(ixre) = t
          else if( nint(xss(jxs(13,iex)+ixre-1)).eq.5001 ) then ! use total photon production
            lt = jxs(12,iex)+ktc(1,iex)
            sig_partg(ixre) = (xss(lt-1)+rtc(1,iex)*(xss(lt)-xss(lt-1)))
          else ! mf=12
            do jy=1,2
              yd(jy) = acefcn(l,xss(jxs(1,iex)
     &                + ktc(1,iex)+jy-2)*1.00000000001d0,ln)
            end do
            if(yd(1).eq.0..and.yd(2).eq.0.) cycle
            ix = nint(xss(l-1)) ! reaction MT that is multiplied
            x1 = 0.d0 ! partial neutron CS
            x2 = 0.d0
            is = jxs(21,iex)+1
            if(ix.gt.0) then
              is = jxs(7,iex) + nint(xss(jxs(6,iex)+ix-1))
            end if
            ic = min(ktc(1,iex)-nint(xss(is-1))+1,nint(xss(is)))
            if(ic.gt.0) then
              x1 = xss(is+ic)
            end if
            ic = min(ic+1,nint(xss(is)))
            if(ic.gt.0) then
              x2 = xss(is+ic)
            end if
            sig_partg(ixre) = x1*yd(1)+(x2*yd(2)-x1*yd(1))*rtc(1,iex)
          end if

          idiscflag = 0 ! cycle if discrete gamma
          ipprodmt = nint(1.d-3*xss(jxs(13,iex)+ixre-1))
          if( ipprodmt.gt.50 .and. ipprodmt.lt.91 ) idiscflag = 1

          ntyn = 0
          qwr = 0.d0
          mtpwr = nint(xss(jxs(13,iex)+ixre-1))
          iawr = jxs(17,iex)
          kawr = 0
          if(jxs(16,iex).ne.0) then
            kawr = nint(xss(jxs(16,iex)+ixre-1))
          end if
          idwr = jxs(19,iex)
          kdwr = nint(xss(jxs(18,iex)+ixre-1))

          trialloopg: do itrial = 1, mxtrial
            call xstcas(1,ilibtype,qwr,iawr,kawr,idwr,kdwr)
            do iebin = 1, nepg
              if ( colout(1,1).gt.emshg(iebin) .and.
     &             colout(1,1).le.emshg(iebin+1) ) then
                do iabin = 1, napg
                  if ( colout(2,1).gt.amshg(iabin+1) .and.
     &                 colout(2,1).le.amshg(iabin) ) then
                    if( idiscflag.eq.1 ) then
                      discgam(ixre) = colout(1,1)
                      iddxg(1,iebin,ixre) = mxtrial
                      cycle mtloopg
                    end if
                    iddxg(1,iebin,ixre) = iddxg(1,iebin,ixre)
     &                                  + 1
                    cycle trialloopg
                  end if
                end do ! iabin
              end if
            end do ! iebin
          end do trialloopg

        end do mtloopg

        erg = esave

        do ixre = 1, nxs(6,iex)
          nswr = abs(nint(xss(jxs(5,iex)+ixre-1)))
          do iebin = 1, nepg
            ddxg(1,iebin) = ddxg(1,iebin) + iddxg(1,iebin,ixre)
     &      /dble(mxtrial)*sig_partg(ixre)
     &      /(emshg(iebin+1)-emshg(iebin))
     &      /(4.d0*pi)
          end do
        end do

! discrete gamma
        if( iinfout.ne.0 ) then
          write(ioa,'("# Partial photon production CS")')
          write(ioa,'("# MT  CS(b)  discrete-gamma(MeV)")')
          do ixre = 1, nxs(6,iex)
            if( sig_partg(ixre).ne.0.d0 ) then
              write(ioa,'(I7,2ES11.3)') nint(xss(jxs(13,iex)+ixre-1)),
     &        sig_partg(ixre), discgam(ixre)
            end if
          end do
          write(ioa,*)
        end if

! gamma-spectrum
        write(iog,'("#",A10,2A11)')'emin(MeV)','emax(MeV)','CS'
        if( iplotlib.eq.1 ) then
          cangelkey = "gamma"
          if( icompare.eq.0 ) then
            call setlibangel(iog,3)
          else if( icompare.eq.1 ) then
            call setlibangel(iog,4)
          end if
        end if
        if( icompare.eq.0 ) then
          do iebin = 1, nepg
            write(iog,'(3ES11.3)') emshg(iebin), emshg(iebin+1),
     &      1.d3*ddxg(1,iebin)*(4.d0*pi)
          end do
        else if( icompare.eq.1 ) then
          do iebin = 1, nepg
            write(iog,'(4ES11.3)') emshg(iebin), emshg(iebin+1),
     &      1.d3*ddxg(1,iebin)*(4.d0*pi), ddxc(1,iebin)
          end do
        end if

      end if ! igamout

*-----------------------------------------------------------------------
      end do nucloop

*-----------------------------------------------------------------------
      close(ioa)
      close(ioc)
      close(iod)
      close(iog)

c frtati 2023/12/07
      if( ilibepsout.eq.1 ) then
        iot = 31
        if( icsout.ne.0 ) then
          fname_angel = cfilesig
          open(iot, file = fname_angel, status = 'unknown' )
          call a_angel(0,fname_angel)
        end if
        if( iddxout.ne.0 ) then
          fname_angel = cfileddx
          open(iot, file = fname_angel, status = 'unknown' )
          call a_angel(0,fname_angel)
        end if
        if( igamout.ne.0 ) then
          fname_angel = cfilegam
          open(iot, file = fname_angel, status = 'unknown' )
          call a_angel(0,fname_angel)
        end if
      end if

*-----------------------------------------------------------------------
      return
      end subroutine acewrite

************************************************************************
      subroutine libout(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
************************************************************************
      use moddas_mesh
      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------
      character chin*200, chlw*200, chcm*200
      character chlc*200
      character dsin(0:9)*200
      dimension idsi(0:9)
      dimension ill(0:9), ilf(0:9)
      character cht*200

      character m_err*200
      common /error/ m_err, l_err, k_err

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
cfrtati 2023/12/08 epsout added
      icsu = 26
      schan(1:icsu) =  (/
     &    'libtype ','nucleus ','eneinc  ','sigout  ','ddxout  ',
     &    'gamout  ','mxsample','e-type  ','a-type  ','file    ',
     &    'plotlib ','info    ','eincsmin','eincsmax','mtout   ',
     &    'delcos  ','eg-type ','ei-type ','part    ','angelsig',
     &    'angelddx','angelgam','elibread','angdisp ','ddxdata ',
     &    'epsout  '/)
      lschn(1:icsu) =  (/
     &     7,         7,         6,         6,         6,
     &     6,         8,         6,         6,         4,
     &     7,         4,         8,         8,         5,
     &     6,         7,         7,         4,         8,
     &     8,         8,         8,         7,         7,
     &     6/)

*-----------------------------------------------------------------------
      ierr  = 0

! default values
      ilibtype = 1
      ilibpart = 2
      nlibnuc = 0
      icsout = 0
      iddxout = 2
      igamout = 0
      iinfout = 1
      nmtout = 0
      iemittype = 1
      iplotlib = 1
      ielibread = 0
      iangdisp = 0
      icompare = 0
      ilibepsout = 0 ! frtati 2023/12/07

      isettype = 0
      isetnuc = 0
      isettrial = 0

      ine = 0
      ina = 0
      ineg = 0
      inei = 0
      iatp = 0 ! frtati 2022/02/01

      eneinc = 15.d0
      eincsmin = 0.d0
      eincsmax = 1.d3
      elibreadmax = 0.d0
      fcmevpu = 1.d0

      mxtrial = 1000000
      iapair = 0
      isetgam = 0

      delcos = 1.d-3

      cfilen(1:100) = ' '
      cangelsig(1:200) = ' '
      cangelddx(1:200) = ' '
      cangelgam(1:200) = ' '
      cfileddxd(1:100) = ' '

      do i = 1, icsu
        ischn(i) = 0
      end do

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------
  140 continue

      call readl(jsn,jsi,dsin,idsi,ill,ilf,'#%!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) goto 800
      if( iskip .ne. 0 ) goto 140

  150 continue

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) goto 800

*-----------------------------------------------------------------------
*        end of libout section
*-----------------------------------------------------------------------
      if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then
        jpn = 1
        goto 800
      end if

*-----------------------------------------------------------------------
*        identify the parameters
*-----------------------------------------------------------------------
      icl = i1

  200 continue

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
  100 continue

      ipm = i
      ischn( ipm ) = ischn( ipm ) + 1

      ic = inumc(chlw,il+1,i3,'=') + 1
      ic = jnumc(chlw,ic,i3)
      if( ic .gt. i3 ) goto 997

      icl = inumc(chlw,ic,i3,';') - 1

*-----------------------------------------------------------------------
      if( ipm.eq.1 ) then ! library type
        isettype = 1

        if( chlw(ic:ic+6).eq.'neutron' ) then
          ilibtype = 1
          ilibpart = 2
        else if( chlw(ic:ic+5).eq.'proton' ) then
          ilibtype = 9
          ilibpart = 1
        else if( chlw(ic:ic+7).eq.'deuteron' ) then
          ilibtype = 31
          ilibpart = 15
        else if( chlw(ic:ic+4).eq.'alpha' ) then
          ilibtype = 34
          ilibpart = 18
        else if( chlw(ic:ic+5).eq.'photon' ) then
          ilibtype = 99
        else
          goto 950
        end if

*-----------------------------------------------------------------------
      else if( ipm.eq.19 ) then ! emission particle type

        if( chlw(ic:ic+6).eq.'neutron' ) then
          iemittype = 1
        else if( chlw(ic:ic+5).eq.'proton' ) then
          iemittype = 9
        else if( chlw(ic:ic+7).eq.'deuteron' ) then
          iemittype = 31
        else if( chlw(ic:ic+5).eq.'triton' ) then
          iemittype = 32
        else if( chlw(ic:ic+2).eq.'3he' ) then
          iemittype = 33
        else if( chlw(ic:ic+4).eq.'alpha' ) then
          iemittype = 34
        else if( chlw(ic:ic+5).eq.'photon' ) then
          isetgam = 1
        else
          goto 951
        end if

*-----------------------------------------------------------------------
      else if( ipm.eq.2 ) then ! nucleus
        isetnuc = 1

        if( deqn1( chlw(ic:ic) ) ) then

          call onum(chlw,ic,icl,cvvv,ierr)
          if( ierr .ne. 0 ) goto 997

          nlibnuc = nint( cvvv )

  143     call readl(jsn,jsi,dsin,idsi,ill,ilf,'#%!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
          if( ierr .ne. 0 ) return
          if( jpn  .eq. 3 ) goto 981
          if( iskip .ne. 0 ) goto 143

          ic = i1

          do k = 1, nlibnuc
            if( ic .gt. i3 ) then
  144         call readl(jsn,jsi,dsin,idsi,ill,ilf,'#%!$',
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
              if( ierr .ne. 0 ) return
              if( jpn  .eq. 3 ) goto 981
              if( iskip .ne. 0 ) goto 144

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
              if( chlw(i:i).ge.'a' .and. chlw(i:i).le.'z' ) then
                isa = isa + 1
                if( isa .eq. 1 ) ica = i
                icb = i
              else if( deqn1( chlw(i:i) ) ) then
                isn = isn + 1
                if( isn .eq. 1 ) icm = i
                icn = i
              else if( dcom2( chlw(i:i) ) .or. i .eq. i3 ) then
                icd = i
                goto 500
              else
                goto 981
              end if
            end do

            icd = i3

  500       continue

            if( ica .eq. 0 .or. icb .eq. 0 ) goto 981
            if( icb-ica .lt. 0 .or. icb-ica .gt. 1 ) goto 981

            cnuc = chlw(ica:icb)//'  '

            do j = 1, 104
              if( cnuc(1:3) .eq. element(j)(1:3) ) then
                icha = j
                goto 450
              end if
            end do
            goto 981

  450       continue

            if( icm .eq. 0 .or. icn .eq. 0 ) then
              ilibnuc(k) = icha * 1000
            else
              if( isn .gt. 3 ) goto 981
              read(chlw(icm:icn),'(i5)') masi
              if( masi .lt. icha ) goto 981

              ilibnuc(k) = icha * 1000 + masi
            end if

            ic = icd + 1
          end do

        else
          goto 981
        end if

*-----------------------------------------------------------------------
      else if( ipm.eq.3 ) then ! eneinc

        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        eneinc = cvvv

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.4 ) then ! csout

        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        icsout = nint( cvvv )

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.5 ) then ! ddxout

        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        iddxout = nint( cvvv )

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.6 ) then ! gamout

        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        igamout = nint( cvvv )

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.7 ) then ! mxsample
        isettrial = 1

        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        mxtrial = nint( cvvv )

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.8 ) then ! e-type

        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997
        if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

        ietp = nint( cvvv )

        call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &              jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &              'e',ietp,ine,emin,emax,edel,isteg)
        if( ierr .ne. 0 ) return

        goto 150

*-----------------------------------------------------------------------
      else if( ipm.eq.9 ) then ! a-type

        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997
        if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

        iatp = nint( cvvv )
        if( abs(iatp).eq.3 .or. abs(iatp).eq.4 ) then
          jatp = 1
          iapair = 1
        else
          jatp = abs( iatp )
        end if

        call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &              jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &              'a',jatp,ina,amin,amax,adel,istag)
        if( ierr .ne. 0 ) return

        goto 150

*-----------------------------------------------------------------------
      else if( ipm.eq.10 ) then ! file
        icf = min( inumc(chlw,ic,icl,' ') - 1, icl )
        cfilen(1:icf-ic+1) = chin(ic:icf)

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.11 ) then ! plotlib
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        iplotlib = nint( cvvv )

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.12 ) then ! iinfout
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        iinfout = nint( cvvv )

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.13 ) then ! eincsmin
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        eincsmin = cvvv

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.14 ) then ! eincsmax
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        eincsmax = cvvv

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.15 ) then ! mtout
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        nmtout = nint( cvvv )

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

        if( deqn4( chlw(ic:ic) ) ) then
          call onum(chlw,ic,icl,cvvv,ierr)
          if( ierr .ne. 0 ) goto 979

          nmtout = nint( cvvv )
          if( nmtout.eq.0 ) goto 979

  151     call readl(jsn,jsi,dsin,idsi,ill,ilf,'#%!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
          if( ierr .ne. 0 ) return
          if( jpn  .eq. 3 ) goto 979
          if( iskip .ne. 0 ) goto 151

          ic = i1

          do k = 1, nmtout
            if( ic .gt. i3 ) then
  152              call readl(jsn,jsi,dsin,idsi,ill,ilf,'#%!$',
     &             jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
              if( ierr .ne. 0 ) return
              if( jpn  .eq. 3 ) goto 979
              if( iskip .ne. 0 ) goto 152

              ic = i1
            end if
              ic = jnumc(chlw,ic,i3)
              call snum(chlw,ic,i3,ic2,cvvv,ierr)
              if( ierr .ne. 0 ) goto 979

              mtout(k) = nint( cvvv )

              ic = ic2
          end do

        else
          goto 979
        end if

*-----------------------------------------------------------------------
      else if( ipm.eq.16 ) then ! delcos
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        delcos = cvvv

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.17 ) then ! eg-type

        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997
        if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

        iegtp = nint( cvvv )

        call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &              jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &              'e',iegtp,ineg,emin,emax,edel,istegg)
        if( ierr .ne. 0 ) return

        goto 150

*-----------------------------------------------------------------------
      else if( ipm.eq.18 ) then ! ei-type

        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997
        if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

        ieitp = nint( cvvv )

        call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &              jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &              'e',ieitp,inei,emin,emax,edel,istegi)
        if( ierr .ne. 0 ) return

        goto 150

*-----------------------------------------------------------------------
      else if( ipm.eq.20 ) then ! angel parameter
        ict = min( ic + 199, i2 )
        cangelsig = chin(ic:ict)

*-----------------------------------------------------------------------
      else if( ipm.eq.21 ) then ! angel parameter
        ict = min( ic + 199, i2 )
        cangelddx = chin(ic:ict)

*-----------------------------------------------------------------------
      else if( ipm.eq.22 ) then ! angel parameter
        ict = min( ic + 199, i2 )
        cangelgam = chin(ic:ict)

*-----------------------------------------------------------------------
      else if( ipm.eq.23 ) then ! elibread
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        elibreadmax = cvvv
        ielibread = 1

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.24 ) then ! angdisp
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        iangdisp = nint( cvvv )

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      else if( ipm.eq.25 ) then ! ddxdata
        icf = min( inumc(chlw,ic,icl,' ') - 1, icl )
        cfileddxd(1:icf-ic+1) = chin(ic:icf)
        icompare = 1

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
cfrtati 2023/12/07
      else if( ipm.eq.26 ) then ! epsout
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        ilibepsout = nint( cvvv )

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
      end if

      goto 140

*-----------------------------------------------------------------------
  800 continue ! read completed

! check
      if( iddxout.eq.2 .and. (ine.eq.0 .or. ina.eq.0) ) then
        m_err = 'Both e-type and a-type must be defined for ddxout = 2'
        ErrCha = ''
        ErrID = 'L:1709/R:libout/F:liboutmod.f'
        goto 999
      else if( iddxout.eq.1 .and. ine.eq.0 ) then
        m_err = 'e-type must be defined for ddxout = 1'
        ErrCha = ''
        ErrID = 'L:1714/R:libout/F:liboutmod.f'
        goto 999
      end if
      if( igamout.eq.1 .and. ineg.eq.0 ) then
        m_err = 'eg-type must be defined for igamout = 1'
        ErrCha = ''
        ErrID = 'L:1720/R:libout/F:liboutmod.f'
        goto 999
      end if
      if( icsout.gt.0 .and. inei.eq.0 ) then
        m_err = 'ei-type must be defined for sigout > 1'
        ErrCha = ''
        ErrID = 'L:1726/R:libout/F:liboutmod.f'
        goto 999
      end if
      if( abs(iatp).eq.3 .and. mod(ina,2).eq.0 ) then
        m_err = 'na must be odd for a-type = 3 of libout section'
        ErrCha = ''
        ErrID = 'L:1732/R:libout/F:liboutmod.f'
        goto 999
      end if
      if( igamout.eq.1 .and. isetgam.eq.1 ) then
        write(*,'("Both part = photon and gamout are defiend. "/
     &  "e-type is used for gamma-spectrum instead of eg-type")')
      end if
      if( iddxout.eq.2 .and. isetgam.eq.1 ) then
        m_err = 'ddxout=2 is not available for part = photon'
        ErrCha = ''
        ErrID = 'L:1742/R:libout/F:liboutmod.f'
        goto 999
      end if
cfrtati 2023/12/13
      if(icsout.lt.0 .and. nmtout.gt.0 .and. minval(mtout(:)).lt.4) then
        m_err = 'sigout must be positive for mt < 4'
        ErrCha = ''
        ErrID = 'L:1749/R:libout/F:liboutmod.f'
        goto 999
      end if

*-----------------------------------------------------------------------
! set factor for MeV -> MeV/u
      if( ilibtype.eq.31 ) then
        fcmevpu = 0.5d0
      else if( ilibtype.eq.34 ) then
        fcmevpu = 0.25d0
      end if

! set igamout if part = photon
      if( isetgam.eq.1 ) then
        iddxout = 0
        igamout = 1
        ineg = ine
        nepg = ineg
        istegg = isteg
        ina = 0
        if( adjustl(cangelgam).eq.' ' ) then
          cangelgam = cangelddx
        end if
      end if

! allocate libout array
      if( abs(iatp).eq.4 ) ina = ina * 2 + 1
      if( iddxout.eq.1 ) ina = 1
      if( iddxout.eq.0 .and. igamout.eq.1 ) ina = 1
      nepx = ine
      napx = ina
      nepg = ineg
      napg = 1
      nepi = inei

      if( ine.eq.0 ) nepx = 1
      if( ina.eq.0 ) napx = 1
      if( ineg.eq.0 ) nepg = 1
      if( inei.eq.0 ) nepi = 1
      call allocate_liboutarray( nepx,napx,nepg,nepi )

! set mesh
      if( iddxout.ge.1 ) then
        emshx(1:1+ine) = gmsh(isteg:isteg+ine)
        if( iddxout.eq.1 ) then
          amshx(1) =  1.d0
          amshx(2) = -1.d0
        else
          if( iatp.lt.0 ) then
            if( abs(iatp).ne.4 ) then
              amshx(1:1+ina) = cos(pi/180.d0*gmsh(istag:istag+ina))
            else
              do i = 1, ina, 2
                amshx(i) = cos(pi/180.d0*gmsh(istag+(i-1)/2))
                if( amshx(i).gt.1.d0-delcos ) then
                  amshx(i+1) = amshx(i) - delcos
                else if( amshx(i).lt.-1.d0+delcos ) then
                  amshx(i) = amshx(i) + delcos
                  amshx(i+1) = amshx(i) - delcos
                else
                  amshx(i) = amshx(i) + 0.5d0 * delcos
                  amshx(i+1) = amshx(i) - delcos
                end if
              end do
            end if
          else
            if( abs(iatp).ne.4 ) then
              do i = 1, ina + 1
                ii = ina - i + 1
                amshx(i) = gmsh(istag+ii)
              end do
            else
              do i = 1, ina, 2
                ii = ina - i + 1
                amshx(i) = gmsh(istag+(ii-1)/2)
                if( amshx(i).gt.1.d0-delcos ) then
                  amshx(i+1) = amshx(i) - delcos
                else if( amshx(i).lt.-1.d0+delcos ) then
                  amshx(i) = amshx(i) + delcos
                  amshx(i+1) = amshx(i) - delcos
                else
                  amshx(i) = amshx(i) + 0.5d0 * delcos
                  amshx(i+1) = amshx(i) - delcos
                end if
              end do
            end if
          end if
        end if
      end if
      if( igamout.ge.1 ) then
        emshg(1:1+ineg) = gmsh(istegg:istegg+ineg)
        amshg(1) =  1.d0
        amshg(2) = -1.d0
      end if
      if( icsout.ge.1 ) then
        emshi(1:1+inei) = gmsh(istegi:istegi+inei)
      end if

! set elibreadmax, and MeV/u -> MeV
      if( ielibread.eq.0 ) then
        if( ilibtype.eq.1 ) then
          if( eneinc.gt.20.d0 .or. emshi(1+inei).gt.20.d0 ) then
            elibreadmax = 200.d0
          else
            elibreadmax = 20.d0
          end if
        else
          elibreadmax = 200.d0 * fcmevpu
        end if
      end if

! MeV/u -> MeV
      eneinc = eneinc / fcmevpu
      emshi(:) = emshi(:) / fcmevpu

! set eincsmin
      if( inei.ne.0 ) then
        eincsmin = emshi(1)
        eincsmax = emshi(inei+1)
      end if

! set filename
      if( cfilen.eq.' ' ) then
        cfilesig = "sig.out"
        cfileddx = "ddx.out"
        cfilegam = "gam.out"
      else
        indlst = lnblnk(cfilen)
        indexp = 0
        indexp = index(cfilen,'.')
        iocount = 0
        if( icsout.ne.0 ) iocount = iocount + 1
        if( iddxout.ne.0 ) iocount = iocount + 1
        if( igamout.ne.0 ) iocount = iocount + 1
        if( iocount.eq.1 .and. indexp.eq.0 ) then
          cfilen = cfilen//".out"
        end if
        if( iocount.eq.1 ) then
          cfiletem = cfilen
          cfilesig = cfiletem
          cfileddx = cfiletem
          cfilegam = cfiletem
        else if( iocount.gt.1 ) then
          cfiletem = cfilen
          if( indexp.eq.0 ) then
            cfiletem(indlst+5:indlst+8) = ".out"
            indexp = indlst+1
          else
            cfiletem(indexp+4:indlst+4) = cfilen(indexp:indlst)
          end if
          cfilesig = cfiletem
          cfileddx = cfiletem
          cfilegam = cfiletem
          cfilesig(indexp:indexp+3) = "_sig"
          cfileddx(indexp:indexp+3) = "_ddx"
          cfilegam(indexp:indexp+3) = "_gam"
        end if
      end if

! read file to compare ddx
      if( icompare.eq.1 ) then
        io = 100
        open(io,file=cfileddxd,status='old')
        ireadcount = 0
        nar = napx
        if( iapair.eq.1 ) nar = napx/2 + 1
        do i = 1, 1000000
          read(io,'(A)') cht
          if( cht(1:2).eq.'h:' ) then
            ireadcount = ireadcount + 1
            read(io,*)
            do iebin = 1, nepx
              read(io,*) blk, blk, ddxc(ireadcount,iebin)
            end do
            if( ireadcount.eq.nar ) exit
          end if
        end do
        close(io)
      end if

*-----------------------------------------------------------------------
      return

*-----------------------------------------------------------------------
  950    m_err = 'Selected libtype is not availavle'
         ErrCha = ''
         ErrID = 'L:1935/R:libout/F:liboutmod.f'
         goto 999
  951    m_err = 'Selected emission particle is not availavle'
         ErrCha = ''
         ErrID = 'L:1939/R:libout/F:liboutmod.f'
         goto 999

  979    m_err = 'Description of mtout parameter is wrong'
         ErrCha = ''
         ErrID = 'L:1944/R:libout/F:liboutmod.f'
         goto 999
  993    m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:1948/R:libout/F:liboutmod.f'
         goto 999
  997    m_err = 'Description of parameter is wrong in libout'
         ErrCha = ''
         ErrID = 'L:1952/R:libout/F:liboutmod.f'
         goto 999
  981    m_err = 'Description of nucleus parameter is wrong in libout'
         ErrCha = ''
         ErrID = 'L:1956/R:libout/F:liboutmod.f'
         goto 999

*-----------------------------------------------------------------------
  999 continue

      l_err = ill(jsn)
      k_err = jsn
      ierr  = 1

*-----------------------------------------------------------------------
      return
      end subroutine libout

************************************************************************
      subroutine echolibout(ioa)
************************************************************************
      include 'err.inc'

*-----------------------------------------------------------------------
      real*8 :: dec(icsu)
      integer :: iec(icsu)
      character*100 :: fmti, fmtd, fmtc
      character*100 :: cec(icsu), ca(icsu)

*-----------------------------------------------------------------------
      write(ioa,'("#",78("-"))')
      write(ioa,'("# Input echo of [libout] section")')
      write(ioa,'("#",78("-"))')

      iec(1) = ilibtype
      iec(2) = nlibnuc
      iec(3) = 1
      iec(4) = icsout
      iec(5) = iddxout
      iec(6) = igamout
      iec(7) = mxtrial
      iec(8) = -99
      iec(9) = -99
      if( cfilen.ne.' ' ) then
        iec(10) = 1
      else
        iec(10) = -99
      end if
      iec(11) = iplotlib
      iec(12) = iinfout
      iec(13) = 1
      iec(14) = 1
      iec(15) = nmtout
      iec(16) = -99
      iec(17) = -99
      iec(18) = -99
      iec(19) = iemittype
      do i = 1, 3
        if( i.eq.1 .and. adjustl(cangelsig).ne.' ' .or.
     &      i.eq.2 .and. adjustl(cangelddx).ne.' ' .or.
     &      i.eq.3 .and. adjustl(cangelgam).ne.' ' ) then
          iec(20+i-1) = 1
        else
          iec(20+i-1) = -99
        end if
      end do
      iec(23) = -99

      dec(3) = eneinc
      dec(13) = eincsmin
      dec(14) = eincsmax
      dec(23) = elibreadmax

      ca(10) = trim(cfilen)
      ca(20) = cangelsig
      ca(21) = cangelddx
      ca(22) = cangelgam

      cec(1) = "(D=proj) Incident particle type of library. "//
     &"1:n, 9:p, 31:d, 34:a"
      cec(2) = "Target nucleus of library"
      cec(3) = "(D=15.d0) Incident energy to derive ddx (MeV)"
      cec(4) = "(D=0) 2:Derive all integrated cross-sections (b), "//
     &"1:derive only total"
      cec(5) = "(D=2) 1:Derive energy-distribution (mb/MeV), "//
     &"2:derive ddx (mb/sr/MeV)"
      cec(6) = "(D=0) 1:Derive gamma-spectrum (mb/MeV)"
      cec(7) = "(D=maxcas) Number of sampling to derive ddx"
      cec(10) = "Output file name"
      cec(11) = "(D=1) 1:Set ANGEL parameter to plot cross-sections"
      cec(12) = "(D=1) 1:Output general information of library"
      cec(13) = "(D=0.d0) Minimum of E_inc range of integrated CS"
      cec(14) = "(D=1.d3) Maximum of E_inc range of integrated CS"
      cec(15) = "(D=0) Specify MTs to derive cross-sections"
      cec(19) = "(D=1) Emission particle type. "//
     &"1:n, 9:p, 31:d, 32:t, 33:3he, 34:a"
      cec(20) = "Angel parameter for integrated cross-section"
      cec(21) = "Angel parameter for DDX"
      cec(22) = "Angel parameter for gamma-spectrum"
      cec(23) = "If dmax is not defined, This value is used "//
     &          "except for dmax(2)"

      fmti = '("  ",A8," =",I10,"  # ",A)'
      fmtd = '("  ",A8," =",ES10.2,"  # ",A)'
      fmtc = '("  ",A8," =  ",A)'

      write(ioa,'("[ libout ]")')
      do ipm = 1, 23
        if( iec(ipm).eq.-99 ) cycle
        if( ipm.eq.3.or.ipm.eq.13.or.ipm.eq.14.or.ipm.eq.23 ) then
          write(ioa,fmtd) schan(ipm), dec(ipm), trim(cec(ipm))
        else if( ipm.eq.10 .or. ipm.ge.20 .and. ipm.le.22 ) then
          write(ioa,fmtc) schan(ipm), trim(ca(ipm))
        else
          write(ioa,fmti) schan(ipm), iec(ipm), trim(cec(ipm))
        end if
        if( ipm.eq.2 ) then
          write(ioa,'(10I8)') (ilibnuc(i),i=1,nlibnuc)
        end if
        if( ipm.eq.15 .and. nmtout.ne.0 ) then
          write(ioa,'(10I5)') (mtout(i),i=1,nmtout)
        end if
      end do
      write(ioa,*)

      write(ioa,'("# Mesh information")')
      if( nepi.gt.1 ) then
        write(ioa,'("# Incident energy mesh for integrated CS")')
        write(ioa,'("  ne = ",I8)') nepi
        write(ioa,'(10ES10.2)') (emshi(i),i=1,nepi+1)
      end if
      if( iddxout.gt.0 ) then
        write(ioa,'("# Energy mesh for differential cross-section")')
        write(ioa,'("  ne = ",I8)') nepx
        write(ioa,'(10ES10.2)') (emshx(i),i=1,nepx+1)
        if( iddxout.eq.2 ) then
          write(ioa,'("# Angle mesh for differential cross-section")')
          write(ioa,'("  na = ",I8)') napx
          write(ioa,'(10ES10.2)') (amshx(i),i=1,napx+1)
        end if
      end if
      if( igamout.gt.0 ) then
        write(ioa,'("# Energy mesh for gamma-emission diff. CS")')
        write(ioa,'("  ne = ",I8)') nepg
        write(ioa,'(10ES10.2)') (emshg(i),i=1,nepg+1)
      end if

      write(ioa,*)

      return
      end subroutine echolibout

************************************************************************
      subroutine setlibangel(io,mode)
************************************************************************
      logical,save :: firstcall(4) = .true.
      character,dimension(50) :: cangfac*14, cangfacd*15

*-----------------------------------------------------------------------
      ioind = io - 100

      if( firstcall(ioind) ) then
        firstcall(ioind) = .false.
      else
        write(io,'("newpage: ")')
      end if

      if( mode.eq.6 .or. mode.eq.8 ) then
        naw = napx
        if( iapair.eq.1 ) naw = napx/2
        do i = 1, naw
         write(cangfac(i),'(" y*",ES7.1,",h0l")') 10.d0**(-i)
        end do
        if( mode.eq.8 ) then
          do i = 1, naw + 1
           write(cangfacd(i),'(" y*",ES7.1,",h0dr")') 10.d0**(-i+1)
          end do
        end if
      end if

      write(io,'("wt:")')
      write(io,'("library = ",A)') trim(cangcom)
      if( mode.gt.2 ) then
        write(io,'("E_{in} =",ES10.2,"[MeV/u]")') fcmevpu*eneinc
      end if
      if( mode.eq.5 .or. mode.eq.7 ) then
        write(io,'("angle:")')
        write(io,'("lower =",ES10.2," [deg]")') angtem(1)
        write(io,'("upper =",ES10.2," [deg]")') angtem(2)
      end if
      write(io,'("e:")')

      write(io,'("p: xlin ylog")')
      if( mode.eq.6 .or. mode.eq.8 ) then
        write(io,'("p: xfac(0.8) form(1.1) scal(1.1) afac(0.8)")')
      end if
      if( ioind.eq.2 .and. adjustl(cangelsig).ne.' ' ) then
        write(io,'("p: ",A)') trim(cangelsig)
      else if( ioind.eq.3 .and. adjustl(cangelddx).ne.' ' ) then
        write(io,'("p: ",A)') trim(cangelddx)
      else if( ioind.eq.4 .and. adjustl(cangelgam).ne.' ' ) then
        write(io,'("p: ",A)') trim(cangelgam)
      end if
      if( mode.eq.0 ) then ! frtati 2023/12/13
        write(io,'("x: incident energy [MeV/u]")')
        write(io,'(A)') "y: " // trim(cangeltitle)
        write(io,'("h: x y(",A,"),3l")') trim(cangelkey)
      else if( mode.eq.1 ) then
        write(io,'("x: incident energy [MeV/u]")')
        write(io,'("y: cross-section [b]")')
        write(io,'("h: x y(",A,"),3l")') trim(cangelkey)
      else if( mode.eq.2 ) then
        write(io,'("x: incident energy [MeV/u]")')
        write(io,'("y: cross-section [b]")')
        write(io,'("h: x y1(total),3l n y2(elastic),3r ",
     &  "y3(non-ela),3b")')
      else if( mode.eq.3 .or. mode.eq.4 ) then
        write(io,'("x: outgoing particle energy [MeV]")')
        write(io,'("y: d\sigma/dE [mb/MeV]")')
        if( mode.eq.3 ) then
          write(io,'("h: x n y(",A,"),h0l")') trim(cangelkey)
        else if( mode.eq.4 ) then
          write(io,'("h: x n y(",A,"),h0l y(data),h0dr")')
     &    trim(cangelkey)
        end if
      else if( mode.ge.5 .and. mode.le.8 ) then
        write(io,'("x: outgoing particle energy [MeV]")')
        write(io,'("y: d^2\sigma/dE/dA [mb/sr/MeV]")')
        if( mode.eq.5 ) then
          write(io,'("h: x n y(",A,"),h0l")') trim(cangelkey)
        else if( mode.eq.6 ) then
          write(io,'("h: x n y(",A,"),h0l",50A)')
     &    trim(cangelkey), (cangfac(i),i=1,naw)
        else if( mode.eq.7 ) then
          write(io,'("h: x n y(",A,"),h0l y(data),h0dr")')
     &    trim(cangelkey)
        else if( mode.eq.8 ) then
          write(io,'("h: x n y(",A,"),h0l",100A)')
     &    trim(cangelkey),(cangfac(i),i=1,naw),(cangfacd(i),i=1,naw+1)
        end if
      end if

      return
      end subroutine setlibangel

************************************************************************
      subroutine setlibparam
************************************************************************
      use moddas_material
      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------
      common /kmat1g/ kmat(kvlmax)

      common /cparm/  maxbch,maxcas
      common /isomul/ smlwt(isrc), totfact, imsrc
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)

*-----------------------------------------------------------------------
      if( isettype.eq.0 ) then
        if( istyp(imsrc).eq.1 ) then
          ilibtype = 9
        else if( istyp(imsrc).eq.2 ) then
          ilibtype = 1
        else if( istyp(imsrc).eq.15 ) then
          ilibtype = 31
          fcmevpu = 0.5d0
        else if( istyp(imsrc).eq.18 ) then
          ilibtype = 34
          fcmevpu = 0.25d0
          write(*,'("Incident particle type defined in [source] ",
     &    "section is not available for inucr=17. ",
     &    "This setting is neglected.")')
        end if
! MeV/u -> MeV
        eneinc = eneinc / fcmevpu
        emshi(:) = emshi(:) / fcmevpu
      end if

      if( isetnuc.eq.0 ) then
        mat = 1
        lem = nint( dnel_das(kmat0+mat) )
        nlibnuc = lem
        do i = 1, nlibnuc
          itz = idnint( zz_das(kmat(mat)+i) )
          ita = idnint( a_das(kmat(mat)+i) )
          ilibnuc(i) = 1000 * itz + ita
        end do
      end if

      if( isettrial.eq.0 ) then
        mxtrial = maxcas
      end if

      return
      end subroutine setlibparam

************************************************************************
      end module liboutmod
