!ggm09.f
!
      logical function checkpncross(eein, ireg,imat)
      implicit real*8 (a-h,o-z)
      common /dpnmaxcom/ dpnmax ! photon nuclear library maximum energy
      checkpncross= .false.
      if(eein.le.dpnmax) then
      checkpncross= .true.
      endif
      return
      end
*
************************************************************************
*                                                                      *
      subroutine xstpni_pn(mm,sigt,ein,mk)
*       calculate photonuclear cross sections in cell icl.             *
*       Last modified by nais                                          *                                                                    *
************************************************************************
      use GGMBANKMOD
      use GGMARRAYMOD
      use moddas_material
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'
      common /pnint/ ipnint
      logical idokind
      common /kmat1g/ kmat(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)
      common /pnixs/ totmpni,pnixsm(pnlmax),ipnikind(pnlmax),egypni,
     &               dlibmax
!$OMP THREADPRIVATE(/pnixs/)
      common /dpnmaxcom/ dpnmax ! photon nuclear library maximum energy
*-----------------------------------------------------------------------
*     initialization of this material and particle
*-----------------------------------------------------------------------

      erg  = ein

        egypni = erg
!
! clear
        totmpni = 0.0d0
        dlibmax = 0.0d0
!
        imat = mk
        do m = jmd(1+mk), jmd(1+mk+1) - 1
          ie = lme(2,m)
          rtc(7,ie) = 0.0
          icntm = m - jmd(1+mk) + 1
          pnixsm(icntm) = 0.0d0
          ipnikind(icntm) = 0
        enddo
        sigt  = 0.0d0
!
*-----------------------------------------------------------------------
*     photonuclear part
*-----------------------------------------------------------------------
        if( ipnint .ne. 0 ) then
*
*!!!     call pnctot(mk,totmpni,ein)
*        calculate the photonuclear cross section in material mk.      *
*        local variables: m - isotope index, it - table index,         *
*          ib - upper energy index, ic - lower energy index,           *
*        global variables:                                             *
*          jmd - material isotope indices                              *
*          lmn - photonuclear isotope table indices                    *
*          erg - current photon energy                                 *
*          pnt(mk) - the lowest photonuclear threshold in material mk  *
*          ktc(1,table) - table index above current energy             *
*          rtc(1,table) - current linear table interpolation factor    *
*          rtc(2,table) - current total cross section                  *
*          rtc(6,table) - current energy being interpolated            *
*-----------------------------------------------------------------------
*        return if below photonuclear threshold for material.
*-----------------------------------------------------------------------
!
!
*-----------------------------------------------------------------------
*        accumulate the cross sections by nuclide in material.
*-----------------------------------------------------------------------

          do m = jmd(1+mk), jmd(1+mk+1) - 1

            idokind = .false.

            if(erg.gt.pnt(mk)) then

              it = lmn(m)
              icntm = m - jmd(1+mk) + 1
              if(it /= 0) then
                if(dpnmax == 0.0d0 .or. erg <= dpnmax) then
                  if(erg <= xss(jxs(1,it)+nxs(3,it)-1)) then
                    dlibmax = max(dlibmax,xss(jxs(1,it)+nxs(3,it)-1))
                    idokind = .true.
                  endif
                end if
              end if
            endif

! max energy
            if(idokind) then
              ic = jxs(1,it)
              ib = ic+nxs(3,it)-1
*
* new setup
              rtc(6,it) = erg
*-----------------------------------------------------------------------
*        find index of current energy.
*        if below photonuclear threshold, use zero values.
*-----------------------------------------------------------------------
              if(erg < xss(ic)) then
                ktc(1,it) = 0
                rtc(1,it) = 0.
                rtc(2,it) = 0.
              else
*-----------------------------------------------------------------------
*        if above last energy value, use boundary values.
*-----------------------------------------------------------------------
                if(erg < xss(ib)) then
*-----------------------------------------------------------------------
*          handle normal value within table.
*-----------------------------------------------------------------------
                  do while(ib-ic /= 1)
                    ih = (ic+ib)/2
                    if(erg <  xss(ih)) then
                     ib = ih
                    else
                      ic = ih
                    end if
                  end do
                  ktc(1,it) = ib-jxs(1,it)
                  rtc(1,it) = (erg-xss(ic))/(xss(ib)-xss(ic))
                  k = ktc(1,it)+jxs(2,it)
                  rtc(2,it) = (xss(k)-xss(k-1))*rtc(1,it)+xss(k-1)
                else
                  ktc(1,it) = ic-jxs(1,it)
                  rtc(1,it) = 0.
                  rtc(2,it) = xss(jxs(2,it)+ktc(1,it))
                end if
              end if
*
*-----------------------------------------------------------------------
*        accumulate the photonuclear microscopic xs for the material.
*-----------------------------------------------------------------------
              pnixsm(icntm) = rtc(2,it)
              ipnikind(icntm) = 1
              totmpni = totmpni+pnixsm(icntm)*fme(m)
*
            else
*
              ie = lme(2,m)
              izm = iza(m) / 1000
              ims = iza(m) - 1000 * izm
              icntm = m - jmd(1+mk) + 1
*-----------------------------------------------------------------------
*     expand natural nucleus and calculate average x-section
*-----------------------------------------------------------------------
              if( ims == 0 ) then
                lemm  = nint( dnel_das(kmat0+imat) )
                sekt = 0.0d0
                sekc = 0.0d0
                do lem = 1, lemm
                  izt = zz_das(kmat(imat)+lem)
                  if( izt  == izm ) then
                    sekt = sekt + den_das(kmat(imat)+lem)
                    inm = nint( a_das(kmat(imat)+lem) ) - izm
                    xsc = phxs(izm,inm,erg)
                    sekc = sekc + xsc * dens(kmat(imat)+lem)
                  end if
                end do
                xsc = sekc / sekt
              else
                inm = ims - izm
                xsc = phxs(izm,inm,erg)
              end if
!
              rtc(7,ie) = xsc
              pnixsm(icntm) = xsc
              ipnikind(icntm) = 2
              totmpni = totmpni + pnixsm(icntm) * fme(m)
!
            end if
          end do
*-----------------------------------------------------------------------
        end if
*-----------------------------------------------------------------------

      sigt  = totmpni
*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
!      subroutine getpnitype(mk,icntm,ikind,sig)                        *
*       calculate photonuclear cross sections in cell icl.             *
*       Last modified by nais                                          *                                                                    *
*       Not used anymore because [data max] section is introduced      *                                                                    *
************************************************************************


************************************************************************
        subroutine cg_coll(eein,wgti,ireg,imat,icntm)
*       calculate a collision of a photon with an atom.               *
************************************************************************
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
     &                      ,iaevt, ibevt, jaevt, jbevt
     &                      ,aevts,aevtr,bevts,bevtr
      use GGMBANKMOD
      use NGSDATAMOD, only : bindeg
      use GGMARRAYMOD
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'
      parameter ( rpmass = 938.27d0, rnmass = 939.58d0 )
      include 'param00.inc'
      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn),  qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /pnint/  ipnint
      common /pionflag/ iqdflag, ipionflag, istrflag, iphreturn
!$OMP THREADPRIVATE(/pionflag/)
      common /kmat1g/ kmat(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)
      common /wparm/  swtm(20), wc1(20), wc2(20)
      common /wparm0/ wc01(20), wc02(20)
!$OMP THREADPRIVATE(/wparm0/)
      common /qmdflg/ irqmd
      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)
      common /tarzmn/ itarz, itarm, itarn
!$OMP THREADPRIVATE(/tarzmn/)
      common /pnixs/ totmpni,pnixsm(pnlmax),ipnikind(pnlmax),egypni,
     &               dlibmax
!$OMP THREADPRIVATE(/pnixs/)
      common /nrfmem/ spis, spgr, levabs, lflgnrf, mpole
!$OMP THREADPRIVATE(/nrfmem/)
      common /nrfdir/ usave(3)
!$OMP THREADPRIVATE(/nrfdir/)
      common /cascid/ jcasc
!$OMP THREADPRIVATE(/cascid/)
!
      dimension suu(3,2)
      dimension xv(3)
      integer igdrqd
      integer jgdrqd
      dimension uvw(3)
      character(len=10) :: ht
      integer, save :: iwng
      data iwng /0/
!$OMP THREADPRIVATE(iwng)
!
      if(mcal .ne. 0) return
!
      nter = 0
      ipt  = 2
      idx  = 0
      erg  = eein
      eg0  = eein
      icl  = ireg
      wgt  = wgti
      wg0  = wgti
      wga  = 0.0d0
      tme  = 0.0d0
      ecmp = 0.0d0
      usave(1) = uuu
      usave(2) = vvv
      usave(3) = www
      uuu  = 0.0
      vvv  = 0.0
      www  = 1.0
      numpal(0:20) = 0
      rumpal(0:20) = 0.0d0
      dxd = 0.0
      xv(1) = 0.0
      xv(2) = 0.0
      xv(3) = 0.0
      uold(1) = 0.0
      uold(2) = 0.0
      uold(3) = 1.0
      ntyn = 0
      jsu  = 0
      mk   = imat

      m = icntm + jmd(1+mk) - 1
      m1 = m
      iexp = m

      iex = lmn(m)

      mathz = iza(iexp) / 1000 ! S.H. 2024.10.15
      mathn = iza(iexp) - 1000 * mathz - mathz ! S.H. 2024.10.15
!
! sample cumulative pn x-section to find the collision nuclide.
!
! photonuclear collision was with nnth nuclide, iex.
!
      ncp = 0
      jsu = 0
      ei  = erg
!
! sample from each species of emitted particle.
! presently, ipt=1 or 2  (neutrons/photons).
! presently, ipt=3 (electron).
      DO_230: do jp=1,nxs(5,iex)
        ipt = ixs(1,jp,iex)
      if( ipt <=0 .or. (ipt >= 4 .and. ipt <= 8) .or.
     &   (ipt >= 10 .and. ipt <= 30 ) .or. ipt >= 35 )  cycle DO_230
!
! determine the total photonuclear production cross-section, tp.
        it = ixs(3,jp,iex)
        if( ktc(1,iex) < nint(xss(it)) )  cycle DO_230
        if( ktc(1,iex) >=nint(xss(it))+nint(xss(it+1))-1 )  cycle DO_230
        ij = it+ktc(1,iex)-nint(xss(it))+2
        tp = xss(ij)+rtc(1,iex)*(xss(ij+1)-xss(ij))
        fp=tp/rtc(2,iex)
        np=int(rang()+fp)
        if(np == 0) cycle DO_230

110     continue

! sample emission parameters and bank each particle.
        DO_220: do ii=1,np
          wgt = wgti
          ixre = 1
          nr = ixs(2,jp,iex)
          if( nr==1 )  go to 180

! sample which reaction, ixre, produces photonuclear particle
! by subtracting partial cross sections, px, from rx=tp*rang().
! return here up to mp=100 times if partials do not add to tp.
          mp = 0
          is = ixs(8,jp,iex)
120       continue
          rx = tp*rang()

! loop over the available reactions until xs fulfilled.
! ks=the next yield block offset locator.
! jx=mftype, the type of data to check.
          DO_170: do ixre_tmp=1,nr
            ixre = ixre_tmp
            ks = nint(xss(ixs(7,jp,iex)+ixre-1))
            jx = nint(xss(is+ks-1))

! handle mftype = 6, 12 or 16.
!   energy dependent yield xs multiplier.
!   data is in format:
!     xss(is+ks-1) = mftype - reaction type indicator
!     xss(is+ks)   = mtmult - mt reaction applied to
!     xss(is+ks+...) std. energy/data table entries
!       nr, [nbt(i:i=1..nr), int(i:i=1..nr),]
!       nv, e(i:i=1..nv), y(i:i=1..nv)
            if( jx==6 .or. jx==12 .or. jx==16 ) then

! get yield value for incident photon energy.
              yd = acefcn_pn(is+ks+1,ei,ln)

! look up pointer to mt reaction xs data.
              do im=1,nxs(4,iex)
                if( nint(xss(is+ks)) == nint(xss(jxs(6,iex)+im-1)) )
     &           go to 140
              end do
              call zaid(2,ht,ixl(1,iex))
              write(*,*) ' *** Warning: '//
     &         'Cound not find mt for yield multiplier. zaid = '//ht//
     &         ' in cg_coll'
              write(*,'(''  ixr = '',i5,'' mt = '',i5)')
     &                  zero+ixre, nint(xss(is+ks))
              write(*,*) ' *** Warning: '//
     &         'Photonuclear yield multiplier mt wrong'//
     &         ' in cg_coll'
              return
140           continue
! get partial production cross section, px=yd*xs from
! yield data, yd, and reaction cross section, xs.
! look up reaction xs using mt index im found above.
              xs = 0.
              it = jxs(9,iex)+nint(xss(jxs(8,iex)+im-1))-1
              ij = it+ktc(1,iex)-nint(xss(it))+2
              if( ktc(1,iex)>=nint(xss(it)) .and.
     &          ktc(1,iex)<nint(xss(it))+nint(xss(it+1))-1 ) then
                xs = xss(ij)+rtc(1,iex)*(xss(ij+1)-xss(ij))
              endif
              px = yd*xs

            else
! handle mftype = 13.
!   get partial production cross section, px, directly.
!   data is in format:
!     xss(is+ks)   = mftype - reaction type indicator
!     xss(is+ks+1) = ie - index of first value
!     xss(is+ks+2) = nv - number of values listed
!     xss(is+ks+3..is+ks+2+nv-1) = pxs(ie..ie+nv-1)

              if( nint(xss(is+ks))==13 ) then
                px = 0.
                it = is+ks+1
                ij = it+ktc(1,iex)-nint(xss(it))+2
                if(ktc(1,iex) >= nint(xss(it)).and.ktc(1,iex)
     &           < nint(xss(it))+nint(xss(it+1))-1)px=xss(ij)+
     &           rtc(1,iex)*(xss(ij+1)-xss(ij))
              else
! bad trouble if not mftype 6, 12, 13 or 16.
                call zaid(2,ht,ixl(1,iex))
                write(*,*) ' *** Warning: '//
     &         'Cound not find reaction type mftype. zaid = '//ht//
     &         ' in cg_coll'
                write(*,'(''  ixr = '',i5,'' mft = '',i5)')
     &             zero+ixre, nint(xss(is+ks-1))
                write(*,*) ' *** Warning: '//
     &                     'Photonuclear reaction type not found'//
     &                     ' in cg_coll'
                return
              endif
            endif

! if the partial xs has been fulfilled, use this reaction, ixre.
            rx = rx-px
            if( rx<0. )  go to 180
          end do DO_170
          ixre = ixre_tmp

! resample partial xs if they do not sum to total, tp.
          mp   = mp+1
          npum = npum+1
          if( mp<100 )  go to 120
          if(iwng .le. 2) 
     &    write(*,*) ' *** Warning: '//
     &   'Photonuclear reaction type not found'//
     &   ' in cg_coll'
          if(iwng .eq. 2) 
     &    write(*,*) ' *** This message is not printed anymore'
          iwng = iwng + 1
          return

180       continue
! ixre (sampled reaction) now found.  sample a new particle.
          erg  = ei
          ipsc = 0
          npa  = 1
          mtp  = nint(xss(ixs(5,jp,iex)+ixre-1))
          ntyn = nint(xss(ixs(6,jp,iex)+ixre-1))
          ia   = ixs(10,jp,iex)
          ka   = nint(xss(ixs(9,jp,iex)+ixre-1))
          id   = ixs(12,jp,iex)
          kd   = nint(xss(ixs(11,jp,iex)+ixre-1))
          call acecas(1,2,zero,ia,ka,id,kd)
          if( kdb/=0 )  return
          erg  = colout(1,1)
          vel  = slite*sqrt(erg*(erg+2.*gpt(ipt)))/(erg+gpt(ipt))
          call rotas(colout(2,1),uold,uvw,lev,irt)
          uuu = uvw(1)
          vvv = uvw(2)
          www = uvw(3)
!
          if( kdb/=0 )  return
!
! neutron
      if(ipt == 1)then
        if(nclst <= 0) then
          nclst =  1
        else
          nclst =  nclst + 1
        end if
!
        iclust(nclst) = 2
        jclust(0,nclst) = 0
        jclust(1,nclst) = 0
        jclust(2,nclst) = 1
        jclust(3,nclst) = 2
        jclust(4,nclst) = 0
        jclust(5,nclst) = 0
        jclust(6,nclst) = 1
        jclust(7,nclst) = 2112
        jclust(8,nclst) = 0
        rms = gpt(1)
        rmg = rms / 1000.
        pr  = sqrt( erg * ( erg + 2.0 * rms ) ) / 1000.
        ett = sqrt( pr**2 + rmg**2 )
        pxl = pr * uuu
        pyl = pr * vvv
        pzl = pr * www

        qclust(0,nclst)  = 0.0d0
        qclust(1,nclst)  = pxl
        qclust(2,nclst)  = pyl
        qclust(3,nclst)  = pzl
        qclust(4,nclst)  = ett
        qclust(5,nclst)  = rmg
        qclust(6,nclst)  = 0.0d0
        qclust(7,nclst)  = erg
        qclust(8,nclst)  = wgt / wg0
        qclust(9,nclst)  = tme
        qclust(10,nclst) = 0.0d0
        qclust(11,nclst) = 0.0d0
        qclust(12,nclst) = 0.0d0

! photon
      else if(ipt == 2) then

        if(nclst <= 0) then
          nclst =  1
        else
          nclst =  nclst + 1
        end if

       iclust(nclst) = 4
       jclust(0,nclst) = 0
       jclust(1,nclst) = 0
       jclust(2,nclst) = 0
       jclust(3,nclst) = 14
       jclust(4,nclst) = 0
       jclust(5,nclst) = 0
       jclust(6,nclst) = 0
       jclust(7,nclst) = 22
       jclust(8,nclst) = 0
       pr   = erg / 1000.0d0
       pxl = pr * uuu
       pyl = pr * vvv
       pzl = pr * www
       qclust(0,nclst)  = 0.0
       qclust(1,nclst)  = pxl
       qclust(2,nclst)  = pyl
       qclust(3,nclst)  = pzl
       qclust(4,nclst)  = pr
       qclust(5,nclst)  = 0.0
       qclust(6,nclst)  = 0.0
       qclust(7,nclst)  = erg
       qclust(8,nclst)  = wgt / wg0
       qclust(9,nclst)  = 0.0
       qclust(10,nclst) = 0.0d0
       qclust(11,nclst) = 0.0d0
       qclust(12,nclst) = 0.0d0
!
! electron
      else if(ipt == 3)then
!
        if(nclst <= 0) then
          nclst =  1
        else
          nclst =  nclst + 1
        end if

        iclust(nclst) = 7      ! 7:e-,e+, 4:photon
        jclust(0,nclst) = 0
        jclust(1,nclst) = 0
        jclust(2,nclst) = 0
        jclust(3,nclst) = 12   ! 12:e-, 13:e+, 14:photon
        jclust(4,nclst) = 0
        jclust(5,nclst) = -1    ! charge
        jclust(6,nclst) = 0
        jclust(7,nclst) = 11  ! kf code
        jclust(8,nclst) = 0

        rms = gpt(3)
        rmg = rms / 1000.
        pr  = sqrt( erg * ( erg + 2.0 * rms ) ) / 1000.
        ett = sqrt( pr**2 + rmg**2 )
        pxl = pr * uuu
        pyl = pr * vvv
        pzl = pr * www

        qclust(0,nclst)  = 0.0
        qclust(1,nclst)  = pxl
        qclust(2,nclst)  = pyl
        qclust(3,nclst)  = pzl
        qclust(4,nclst)  = ett
        qclust(5,nclst)  = rmg
        qclust(6,nclst)  = 0.0
        qclust(7,nclst)  = erg
        qclust(8,nclst)  = wgt / wg0
        qclust(9,nclst)  = 0.0
        qclust(10,nclst) = 0.0
        qclust(11,nclst) = 0.0
        qclust(12,nclst) = 0.0
!
! proton deuteron triton helium3 and alpha
      else if(ipt == 9 .or.  (ipt >= 31 .and. ipt <= 34) )then

           if(ipt == 9) then ! proton
             itypin = 1
             iptin = 9
             kfin = 2212
             izin = 1
             inin = 0
           elseif(ipt == 31) then ! deuteron
              itypin = 15
              iptin = 31
              kfin = 1000002
              izin = 1
              inin = 1
            elseif(ipt == 32) then ! triton
              itypin = 16
              iptin = 32
              kfin = 1000003
              izin = 1
              inin = 2
            elseif(ipt == 33 ) then ! helium3
              itypin = 17
              iptin = 33
              kfin = 2000003
              izin = 2
              inin = 1
            elseif(ipt  == 34 ) then ! alpha
              itypin = 18
              iptin = 34
              kfin = 2000004
              izin = 2
              inin = 2
            endif

            if(nclst <= 0) then
              nclst =  1
            else
              nclst =  nclst + 1
            end if

            iclust(nclst) = 1

            jclust(0,nclst) = 0
            jclust(1,nclst) = izin
            jclust(2,nclst) = inin
            jclust(3,nclst) = itypin
            jclust(4,nclst) = 0
            jclust(5,nclst) = izin
            jclust(6,nclst) = izin+inin
            jclust(7,nclst) = kfin
            jclust(8,nclst) = 0
            rms = rmtyp(itypin,kfin)
            rmg = rms / 1000.
            pr  = sqrt( erg * ( erg + 2.0 * rms ) ) / 1000.
            ett = sqrt( pr**2 + rmg**2 )
            pxl = pr * uuu
            pyl = pr * vvv
            pzl = pr * www
            qclust(0,nclst)  = 0.0
            qclust(1,nclst)  = pxl
            qclust(2,nclst)  = pyl
            qclust(3,nclst)  = pzl
            qclust(4,nclst)  = ett
            qclust(5,nclst)  = rmg
            qclust(6,nclst)  = 0.0
            qclust(7,nclst)  = erg
            qclust(8,nclst)  = wgt / wg0
            qclust(9,nclst)  = tme
            qclust(10,nclst) = 0.0d0
            qclust(11,nclst) = 0.0d0
            qclust(12,nclst) = 0.0d0
          endif
!

        end do DO_220
!
      end do DO_230
!
      call nevap(0)
      rnpnt(5) = rnpnt(5) + 1.0d0
      rnpntr(5) = rnpntr(5) + 1.0d0
      kdecay(4) = 0
!
      uuu = usave(1)
      vvv = usave(2)
      www = usave(3)
      rncnt(33) = rncnt(33) + 1.0d0
      call cputime(33)
!
      return
      end subroutine cg_coll
!
      subroutine acecas(ls,ip,q,ia,ka,id,kd)
! sample the emission energy and scattering angle from the
! appropriate law data in the laboratory coordinate system.
!
! ls - the current particle index in the colout array
! ip - the ipt particle type for the incident particle
!  q - the q value for the reaction being sampled
! ia - the first word of the relevant AND block in the xss array
! ka - the offset to the first word of the table in AND block
! id - the first word of the relevant DLW block in the xss array
! kd - the offset to the first word of the table in DLW block
!
! returns the sample emission parameters in the lab system:
!   colout(1,ls) - the emission energy
!   colout(2,ls) - the emission scattering angle
!
! Law 4/44/61 makes use of a biased distribution which
! can affect the outgoing particle weight, wgt.
!
! ipsc, kdb, and tpd also modified.
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      real*8, parameter :: ep = 0.000001
      character(len=10)     :: ht

      integer, save :: iwng1
      data iwng1 /0/
!$OMP THREADPRIVATE(iwng1)

! ***********************************************************************
! select the law.
      nx = 0
      colout(3,ls) = 0.
10    continue
      t1 = rang()
      n  = id+kd
      do
        if( nint(xss(n-1))==0 )  exit
        t1 = t1-acefcn_pn(n+2,erg,ln)
        if( t1<0. )  exit
        n = id+nint(xss(n-1))
      enddo

! ***********************************************************************
! use the selected law to sample the energy (and possibly angle).
! if law samples without error, go to sample angle or coordinate
! transform as appropriate.
      colout(1,ls) = -huge
      lw = nint(xss(n))
      iw = id-1+nint(xss(n+1))

      select case( lw )

        case( 1 )
! law 1 (from endf Law 1) -- tabular equiprobable energy bins.
          call acetbl(iw,ic,r,ln)
          nt = nint(xss(iw+ln))
          iw = iw+ln+nt*(ic-1)
          k = int(rang()*(nt-1)+1)
          fr = rang()
          if( r==0. ) then
! sample from single table.
            colout(1,ls) = xss(iw+k)+fr*(xss(iw+k+1)-xss(iw+k))
          else
! sample by scaled interpolation between tables.
            t1 = xss(iw+1)+r*(xss(iw+1+nt)-xss(iw+1))
            i  =  iw
            if( rang()<=r )  i=i+nt
            colout(1,ls) =
     & t1+(xss(i+k)+fr*(xss(i+k+1)-xss(i+k))-xss(i+1))*
     & (xss(iw+nt)+r*(xss(iw+2*nt)-xss(iw+nt))-t1)/(xss(i+nt)-xss(i+1))
          endif
          go to 260

        case( 2 )
! law 2 -- discrete photon lines.
! photon production from neutrons only.
! (see endf-102 rev. 2/97 manual page 12.3).
          colout(1,ls) = xss(iw+1)
          if( nint(xss(iw))==2 ) colout(1,ls)=colout(1,ls)+
     &     erg*awn(iex)/(awn(iex)+1.)
          go to 260

        case( 3, 33 )
! law 3 & 33 -- level scattering.
! law 3 applies only for neutron in-neutron out scattering.
! law 33 for a more general combination of particle types.
          colout(1,ls) = xss(iw+1)*(erg-xss(iw))
          go to 260

        case( 4, 44, 61 )
!  law 4 (from endf law 1) -- continuous erg tabular distribution.
!  law 44 (from endf law 1) -- kalbach-87 correlated formalism.
!  law 61 (from endf law 1) -- correlated tab energy-angle dist.
          call acetbl(iw,ic,r,ln)
          nr = nint(xss(iw))
          lb = iw-1+ln+ic
          lc = id+nint(xss(lb))
          ld = lc
          lf = lc
!
! xss(lc,ld,lf) is an overloaded variable.  it contains the
! number of points in the emission distribution and if part
! of the continuum has been expunged, it contains 0.5 times
! the cumulative probability of the portion expunged.
          np = int(xss(lc)+ep)
          mp = np
          jj = nint(xss(lc-1))
          nd = 0
          if(jj < 9999)nd=jj/10
          if( r/=0. ) then
            ld = id+nint(xss(lb+1))
            mp = int(xss(ld)+ep)
            t1 = xss(lc+nd+1)+r*(xss(ld+nd+1)-xss(lc+nd+1))
            t2 = xss(lc+np)+r*(xss(ld+mp)-xss(lc+np))
            ra = rang()
            if( ra<r ) then
              lf = ld
              jj = nint(xss(ld-1))
              if( jj<9999 .and. jj/10/=nd )
     &         write(*,*) ' *** Warning: '//
     &        'Wrong number of discrete lines for law 4/44'//
     &        ' in acecas'
              if( jj>10000 )  nd=0
            endif
          endif
          if( jj==9999 )  return
          jj = jj-nd*10
          r1 = rang()
!
! check for a hit in the discrete part.
! use histogram or corresponding-point interpolation.
          if( nd/=0 ) then
            do ih=1,nd
              cc = xss(lc+2*np+ih)+r*(xss(ld+2*mp+ih)-xss(lc+2*np+ih))
              if( cc>=r1 ) then
                t = xss(lc+ih)+r*(xss(ld+ih)-xss(lc+ih))
                if( t<0. )  t=erg*awn(iex)/(awn(iex)+1.)-t
                if( lw==44 ) then
                  tpd(1) =
     &        xss(lc+ih+3*np)+r*(xss(ld+ih+3*mp)-xss(lc+ih+3*np))
                  tpd(2) =
     &        xss(lc+ih+4*np)+r*(xss(ld+ih+4*mp)-xss(lc+ih+4*np))
                endif
                go to 152
              endif
            end do
            if( np==nd .and. mp==nd )  go to 152
          endif
!
! handle a hit in the continuous part.
! use histogram or unit-base interpolation inside table.
! use scaled interpolation between tables.
          np = int(xss(lf)+ep)
          ns = nd+1
          if( nd/=0 )  r1=1.-(1.-r1)*(1.-xss(lf+2*np+nd))/(1.-cc)

! adjust for the energy cutoff if necessary.
! see subroutine expung for law 4/44 distribution. 0 < wc < 1
          if( jj>=9999 ) then
            jj=jj/10000
            ns=nint(xss(lf-1))-jj*10000
            wc=2.*(xss(lf)-np)
            r1=1.-wc*(1.-r1)
            wgt=wgt*wc
          endif
!
! find the energy bin in the continuous part.
          ic = lf+2*np+ns
          ib = lf+3*np
          do
            if( ib-ic==1 )  exit
            ih = (ic+ib)/2
            if( r1>=xss(ih) ) then
              ic = ih
            else
              ib = ih
            endif
          enddo
          ln = ic-2*np
          fa = xss(ln+np)
          ea = xss(ln)
!
! sample from linear-linear interpolated bin.
          if( jj==1 )  go to 140
          bb = (xss(ln+np+1)-fa)/(xss(ln+1)-ea)
          if( bb==0 )  go to 140
          t = ea+(sqrt(max(zero,fa**2+2.*bb*(r1-xss(ic))))-fa)/bb
          if( lw==44 ) then
            fb = (t-xss(ln))/(xss(ln+1)-xss(ln))
            tpd(1) = xss(ic+np)+fb*(xss(ic+1+np)-xss(ic+np))
            tpd(2) = xss(ic+2*np)+fb*(xss(ic+1+2*np)-xss(ic+2*np))
          endif
          go to 150
!
! sample from histogram bin.
140       continue
          t = ea+(r1-xss(ic))/fa
          if( lw==44 ) then
            tpd(1) = xss(ic+np)
            tpd(2) = xss(ic+2*np)
          endif
!
! use scaled interpolation between energies.
150       continue
          if(r /= 0.)t=t1+(t-xss(lf+nd+1))*(t2-t1)
     &                /(xss(lf+np)-xss(lf+nd+1))
152       continue
          colout(1,ls) = t
          if( lw==4 ) then
            go to 260
          elseif( lw==44 ) then
! sample law 44 -- kalbach-87 angular systematics.
! tpd(1)=r, tpd(2)=a
            if( ka/=-1 .or. ntyn>=0 )  go to 295
            ipsc = 14
            if( rang()>tpd(1) ) then
              t1 = (2.*rang()-1.)*sinh(tpd(2))
              colout(2,ls) = log(t1+sqrt(t1**2+1.))/tpd(2)
            else
              r2 = rang()
              colout(2,ls) = log(r2*exp(tpd(2))+(1.-r2)
     %                     *exp(-tpd(2)))/tpd(2)
            endif
            go to 280
          elseif( lw==61 ) then
! sample law 61 -- tabulated angular distribution.
            ipsc = 16
!
! if jj=1 (i.e., histogram on e-primes) always use "ic."
! if jj=2 (lin-lin on e-primes), use the distribution for
! the e-prime closest to "r1" (in cdf space).
            lb = ic+np
            if( jj/=1 .and. (xss(ib)-r1 < r1-xss(ic)) )  lb=lb+1
!
! sample from appropriate distribution.  unlike the and block, in
! law 61 only isotropic or tabular angular information is passed.
            ixcos = 0
            lm = nint(xss(lb))
            if( lm==0 )  colout(2,ls)=2.*rang()-1.
            if( lm>0  )  colout(2,ls)=acecos_pn(id,-lm)
            if( lm<0  )  write(*,*) ' *** Warning: '//
     &                  'Data: law 61 contained a bad table pointer'//
     &                  ' in acecas'
            go to 280
          else
            go to 295
          endif

        case( 5 )
!  law 5 (from endf law 5) -- general evaporation spectrum.
          t1 = acefcn_pn(iw,erg,ln)
          i = iw+ln+1+int(rang()*(nint(xss(iw+ln))-1))
          colout(1,ls) = t1*(xss(i)+rang()*(xss(i+1)-xss(i)))
          go to 260

        case( 7 )
!  law 7 (from endf law 7) -- simple maxwell fission spectrum.
          t1 = acefcn_pn(iw,erg,ln)
          t3 = erg-xss(iw+ln)
          if( t3>0. ) then
            do
              t4 = rang()**2
              t2 = t4+rang()**2
              if( t2<=1. ) then
                t2 = log(rang())*t4/t2
                colout(1,ls) = -t1*(t2+log(rang()))
! reject if outside range 0 ... e-u
                if( colout(1,ls)<=t3 )  exit
              endif
            enddo
          endif
          go to 260

        case( 9 )
!  law 9 (from endf law 9) -- evaporation spectrum.
          t1 = acefcn_pn(iw,erg,ln)
          t2 = erg-xss(iw+ln)
          if( t2>0. ) then
            do
              fr = rang()
              colout(1,ls) = -t1*log(fr*rang())
! reject if outside range 0 ... e-u
              if( colout(1,ls)<=t2 )  exit
            enddo
          endif
          go to 260

        case( 11 )
! law 11 (from endf law 11) -- energy dependent watt spectrum.
          t1=acefcn_pn(iw,erg,ln)
          t2 = acefcn_pn(iw+ln,erg,lb)
          if( erg>xss(iw+ln+lb) ) then
            t5=sqrt((1.+.125*t1*t2)**2-1.)+1.+.125*t1*t2
            do
              t = -log(rang())
              colout(1,ls) = t1*t5*t
              if( ((1.-t5)*(1.+t)-log(rang()))**2 <= t2*colout(1,ls) )
     &        exit
            enddo
          endif
          go to 260

        case( 22 )
!  law 22 (from uk law 2) -- tabular linear functions.
          call acetbl(iw,ic,r,ln)
          ie = id-1+nint(xss(iw+ln+ic-1))
          nf = nint(xss(ie))
          do
            iw = ie
            fr = rang()
            do
              iw = iw+1
              fr = fr-xss(iw)
              if( fr<0. )  exit
             enddo
          if( iw<=ie+nf )  exit
          enddo
          colout(1,ls) = xss(iw+2*nf)*(erg-xss(iw+nf))
          go to 260

        case( 24 )
!  law 24 (from uk law 6) -- tabular energy multipliers.
          call acetbl(iw,ic,r,ln)
          i = iw+ln+1+nint(xss(iw+ln))*(ic-1)+int(rang()
     &      *(nint(xss(iw+ln))-1))
          colout(1,ls) = erg*(xss(i)+rang()*(xss(i+1)-xss(i)))
          go to 260

        case( 66 )
!  law 66 (from endf law 6) -- n-body phase space distribution.
          nb = nint(xss(iw))
          ap = xss(iw+1)
          if(ipt > 2)ap=ap*gpt(1)/gpt(ipt)
          do
            r = rang()**2
            s = r+rang()**2
            if( s<=1. )  exit
          enddo
          x = -r*log(s)/s-log(rang())
          do
             r = rang()**2
             s = r+rang()**2
            if( s<=1. )  exit
          enddo
          p = rang()
          if( nb==3 ) then
            y = -r*log(s)/s-log(p)
          elseif( nb==4 ) then
            p = p*rang()
            y = -r*log(s)/s-log(p)
          else
            p = p*rang()*rang()
            p = p*rang()
            y = -r*log(s)/s-log(p)
          endif
          t = x/(x+y)
          aw = awn(iex)
          colout(1,ls) = t*((ap-1.)/ap)*(erg*aw/(aw+1.)+q)
          colout(2,ls) = 2.*rang()-1.
          ixcos = 0
          if( ntyn>=0 )  go to 295
          go to 280

        case( 67 )
!  law 67 (endf/b-vi law 7) -- correleted energy-angle scatter.
          call acetbl(iw,ic,r,ln)
          cs = acecos_pn(ia,ka)
          ipsc = 15
          colout(2,ls) = cs
          colout(1,ls) = acecs6_pn(0,id,iw,ic,r,cs)
          if( ntyn<=0 )  go to 295
          if( colout(1,ls)<0. )  go to 300
          go to 290

        case default
          go to 300
      end select

! ***********************************************************************
! if not correlated energy-angle, calculate the cosine.
260   continue
      colout(2,ls) = acecos_pn(ia,ka)

! ***********************************************************************
! adjust if energy and cosine are given in center-of-mass system.
280   continue
      if( colout(1,ls)<0. )  go to 300
      ergace = colout(1,ls)
      if( ntyn<0 ) then
! formulas below are from p. 2 of x-6:res-93-68.
! these formulas assume two-body kinematics.
! seamon's formulas specify atomic weight ratios (to neutron)
! for incident particle (a), awr=gpt(ipt_incident)/gpt(1)
! for exiting particle (b), awr=gpt(ipt)/gpt(1)
! for target (a), awr=awn(iex)
c        if(.false.) then  ! classical conversion
c            a1 = gpt(ip)/gpt(1)
c            a2 = gpt(ipt)/gpt(1)
c            a3 = awn(iex)
c            t1  =  colout(1,ls)
c            t2 = erg*a1*a2/(a3+a1)**2
c            t3 = 2.*sqrt(a2*a1*erg*colout(1,ls))*colout(2,ls)/(a3+a1)
c            t4 = t1+t2+t3
c            s1 = colout(2,ls)*sqrt(colout(1,ls)/t4)
c            s2 = sqrt(a1*a2*erg/t4)/(a3+a1)
c            colout(1,ls) = t4
c            colout(2,ls) = s1+s2
c        elseif(.true.) then  ! Relativistic conversion
         call relphocvt(colout(1,ls),colout(2,ls),erg,gpt(ipt),awn(iex)
     &    ,lw)
c        endif
      endif

! ***********************************************************************
! resample energy up to 100 times if > emx(1).
290   continue
      if(ipt >= 31 .and. ipt <= 34) then
        if(emx(ipt) > 0.0d0) then
          if( colout(1,ls)<=emx(ipt) )  return
        else
          if( colout(1,ls)<=emx(9) )  return
        endif
      end if
      if( colout(1,ls)<=emx(ipt) )  return
      

      if(iwng1 .lt. 10) then

       write(*,*) ' *** Warning: '//
     &         'Energy of particle from inelastic collision > emx'//
     &         ' in acecas'
       write(*,'(''  egr = '',1pe12.5)') erg
       iwng1 = iwng1 + 1
      elseif( iwng1 .eq. 10 )  then
      write(*,*) ' *** Warning: '//
     &           'Erg > emx happened 100 times'//
     &           ' in acecas. This message is not printed anymore'
       iwng1 = iwng1 + 1
      endif
      
      nx = nx + 1
      if(nx .lt. 100) goto 10 ! return to find energy. 
      
      return

! ***********************************************************************
! print debug information for cross-section table errors.
295   continue
      colout(1,ls) = huge

300   continue
      call zaid(2,ht,ixl(1,iex))

      write(*,310)  ht,erg,ixre,mtp,ntyn,lw,colout(1,ls)
310   format(/  ' error in cross-section table ',a10/ ' energy in =',
     & 1pe12.4,5x,  'reaction index =',i3,5x, 'mt =',i4,5x, 'ty =',i4,
     & 5x, 'law =',i3,5x,  'energy out =',1pe12.4)

      if(colout(1,ls)==-huge ) then
        write(*,*) ' *** Warning: '//
     &      'An inappropriate or non-existent law was selected'//
     &      ' in acecas'
      elseif( colout(1,ls)<0.     ) then
        write(*,*) ' *** Warning: '//
     &              'Emission energy was negative'//
     &              ' in acecas'
      elseif( colout(1,ls)==huge  ) then
        write(*,*) ' *** Warning: '//
     &              'Faulty cross-section data'//
     &              ' in acecas'
      elseif( colout(1,ls)>erg    ) then
        write(*,*) ' *** Warning: '//
     &              'Emission energy exceeds incident energy'//
     &              ' in acecas'
      endif
      return
      end subroutine acecas
!
      subroutine relphocvt(eout,aout,erg,rmn,rmtr,law)
! Convert energy in CM frame to Lab frame in relativistic manners
! Input and output
! eout : outgoing neutron energy (converted from CM to Lab)
! aout : outgoing neutron angular cosine (in Lab all the way)
!
! Input
! erg  : Incoming photon Energy (MeV)
! rmn  : neutron mass used in nuclear data (MeV/c**2)
! rmtr : target mass in unit of rmn
!
! Created on 2023/6/30 by Ogawa
      implicit real*8 (a-h,o-z)

      rmt  = rmtr * rmn ! target mass by unit of neutron mass
      
! frame conversion from Lab to CM. 
      if(law .ne. 3 .and. law .ne. 33) then ! frame conversion from Lab to CM. Deuterium is exempted because it is sampled in CM. This might be wrong for X-sec other than JENDL5
          pnlb = sqrt(eout**2 + 2.d0 * rmn * eout) ! neutron momentum in lab
          preslb = sqrt(erg**2 - 2.d0 * erg * pnlb * aout + pnlb**2) ! residue momentum in lab
          eout = abs(eout - (sqrt(preslb**2 + rmt**2) - rmt)) ! 2024/6/28  take abs because eout is negative for low initial eout. Think some day. 
      endif

! Doppler broadening 
      vf   = erg/sqrt((erg/2.d0)**2+rmt**2) ! velocity of CM frame from lab frame
      vncm = sqrt(1.d0-1.d0/(eout/rmn + 1.d0)**2) ! velocity of neutron in CM frame
      vncmlon = vncm * aout ! longitudinal velocity
      vncmtra = vncm * sqrt(1.d0-min(1.d0,aout**2)) ! transverse velocity
      vnlblon = (vf+vncmlon)/(1.d0+vf*vncmlon)! longitudinal velocity of neutron in lab frame
      vnlb    = sqrt(vnlblon**2+vncmtra**2)   ! total velocity of neutron in lab frame
      eout    = rmn * (1.d0/sqrt(1.d0-vnlb**2)-1.d0)
c      aout    = vncmtra/sqrt(vncmtra**2+vnlblon**2)

      end subroutine relphocvt

      subroutine acetbl(l,ii,r,ln)
! get interpolation parameters r and ii of the energy in the
! table at xss(l).
! data:  nr,nbttc(i=1,nr),int(i=1,nr),nf,e(i=1,nf)
      use GGMARRAYMOD !2020ASTOM

      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'
!
      nr = nint(xss(l))
      ie = l+2*nr+1
      nf = nint(xss(ie))
      ln = 2*(nr+1)+nf
      r = 0.
!
! use extreme value if energy is off either end of the table.
      if( erg>=xss(ie+nf) )  go to 70
      if( erg<=xss(ie+1) )  go to 80
!
! binary search for the location of the energy in the table.
      ic = ie+1
      ib = ie+nf
10    continue
      if( ib-ic==1 )  go to 30
      ih = (ic+ib)/2
      if( erg<xss(ih) )  go to 20
      ic = ih
      go to 10
20    continue
      ib = ih
      go to 10
30    continue
      ii = ic-ie
!
! calculate interpolation fraction r unless int=1 (histogram).
      if( nr==0 )  go to 60
      do n = 1,nr
        if( ib-ie<=nint(xss(l+n)) )  go to 50
      enddo
      n = nr
50    continue
      if( nint(xss(l+nr+n))==1 )  return
60    continue
      if( erg-xss(ic)<1.e-6*(xss(ib)-xss(ic)) )  return
      r = (erg-xss(ic))/(xss(ib)-xss(ic))
      return
!
70    continue
      ii = nf
      return
80    continue
      ii = 1
      return
      end subroutine acetbl
!
      function acefcn_pn(l,eg,ln)
! evaluate a function of eg from the table at xss(l).
! also return ln, the length of the table.
! data:  nr,nbttc(i=1,nr),int(i=1,nr),nf,e(i=1,nf),f(i=1,nf)
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'
!
       nr = nint(xss(l))
      ie = l+2*nr+1
      nf = nint(xss(ie))
      ln = 2*(nr+nf+1)

! use the extreme value if eg is off either end of the table.
      if( eg>=xss(ie+nf) ) then
        acefcn_pn = xss(ie+2*nf)
        return
      elseif( eg<=xss(ie+1)  ) then
        acefcn_pn = xss(ie+nf+1)
        return
      endif

! binary search for the location of eg in the table.
      ic = ie+1
      ib = ie+nf
10    continue
      if( ib-ic==1 )  go to 30
      ih = (ic+ib)/2
      if( eg>=xss(ih) ) then
        ic = ih
      else
        ib = ih
      endif
      go to 10
30    continue
      ea = xss(ic)
      eb = xss(ib)
      fa = xss(ic+nf)
      fb = xss(ib+nf)

! find out which kind of interpolation should be used.
      if( nr==0 ) then
        acefcn_pn = fa+(fb-fa)*(eg-ea)/(eb-ea)
       return
      endif
      do n = 1,nr
        if( ib-ie<=nint(xss(l+n)) )  go to 50
      enddo
      n = nr

! interpolate between table entries.
50    continue
      select case( nint(xss(l+nr+n)) )
        case( 1 )
          acefcn_pn = fa
        case( 2 )
          acefcn_pn = fa+(fb-fa)*(eg-ea)/(eb-ea)
        case( 3 )
          acefcn_pn = fa+(fb-fa)*log(eg/ea)/log(eb/ea)
        case( 4 )
          acefcn_pn = fa*(fb/fa)**((eg-ea)/(eb-ea))
        case( 5 )
          acefcn_pn = fa*(fb/fa)**(log(eg/ea)/log(eb/ea))
      end select
      return
      end function acefcn_pn
!
      function acecos_pn(ia,ka)
! sample a cosine from the distribution at ia,ka.  return the
! cosine (acecos) and the index in xss of the first word of the
! cosine table used (ixcos). ka=0: isotropic. ka<0: tabular only.
! 32-equiprobable bin cosine data:
!    ne,(e(i),i=1,ne),(lmu(i),i=1,ne),(mu(lmu+j),j=1,33)
! tabular probability angular distribution data:
!    ne,(e(i),i=1,ne),(lmu(i),i=1,ne),jj,np,
!    (mu(lmu+j),j=1,np),(pdf(lmu+j),j=1,np),(cdf(lmu+j),j=1,np)
!    where jj is the interpolation flag at lmu(i).
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

! find the cosine table by binary search on the energy table.
      if( ka==0 ) then
! isotropic case
        acecos_pn = 2.*rang()-1.
        ixcos = 0
        return
      endif
      lm = ka
      if( ka>0 ) then
        ic = ia+ka
        n  = nint(xss(ic-1))
        ib = ic-1+n
        do
          if( ib-ic==1 )  exit
          ih = (ic+ib)/2
          if( erg>=xss(ih) ) then
            ic = ih
          else
            ib = ih
          endif
        enddo
!
! sample between adjoining tables by interpolation fraction.
        if( rang()*(xss(ib)-xss(ic)) < erg-xss(ic) )  ic=ib
        lm = nint(xss(ic+n))
        if( lm==0 ) then
! isotropic case
          acecos_pn = 2.*rang()-1.
          ixcos = 0
          return
        elseif( lm>0 ) then
! sample from table of 32 equiprobable cosine groups.
          t1 = rang()*32.
          kr = t1
          ixcos = ia-1+lm
          acecos_pn = xss(ixcos+kr)+(t1-kr)
     &              *(xss(ixcos+kr+1)-xss(ixcos+kr))
          return
        endif
      endif

! tabular probability angular distribution.
      k = ia-1-lm
      ixcos = -k
      jj = nint(xss(k))
      np = nint(xss(k+1))
      rn = rang()
!
! binary search of cumulative density function.
      ic = k+2*np+2
      ib = k+3*np+1
      do
        if( ib-ic==1 ) exit
        ih = (ic+ib)/2
        if( rn>=xss(ih) ) then
          ic = ih
        else
          ib = ih
        endif
      enddo
!
      fa = xss(ic-np)
      ca = xss(ic-2*np)
      acecos_pn = ca+(rn-xss(ic))/fa
      if( jj==1 )  return
      bb = (xss(ic-np+1)-fa)/(xss(ic-2*np+1)-ca)
      if( bb/=0. )
     &  acecos_pn=ca+(sqrt(max(zero,fa**2+2.*bb*(rn-xss(ic))))-fa)/bb
      return
      end function acecos_pn
!
      function acecs6_pn(ii,id,iw,jc,r,cs)
! sample energy acecs6_pn given angle cs, law 67 (endf/b-vi law 7).
! ii=0 acecs6_pn called from acecas (transport): save id,iw,jc,r
!      and random numbers generated in sampling energy;
! ii=1 acecs6_pn called from calcps (next-event estimator):
!      use id,iw,jc,r, and random numbers saved from transport.
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      dimension el(2), eh(2)
      integer    :: lx(2)

      t = 0.
      if( ii==0 ) then
        id0 = id
        iw0 = iw
        ic0 = jc
        rr0 = r
      endif

! go to appropriate table for incident energy
      jw = iw0+2*nint(xss(iw0))+1
      ne = nint(xss(jw))
      lx(1) = id0+nint(xss(jw+ne+ic0))-1
      lx(2) = 0
      ix = 1
      ir = 0
      if( rr0/=0. ) then
        lx(2) = id0+nint(xss(jw+ne+ic0+1))-1
        if( ii==0 )  rnb(1) = rang()
        ir = 1
        if( rnb(1)<=rr0 )  ix = 2
      endif
      DO_110: do m = 1,2
        if( lx(m)==0 )  cycle DO_110
        le = lx(m)
        mu = nint(xss(le))
        nm = nint(xss(le+1))
!
! find appropriate table for sampled cosine
        do iq = 1,nm-1
          if( xss(le+iq+2)>=cs )  exit
        enddo
        if( mu/=1 ) then
          ir = ir+1
          if( ii==0 )  rnb(ir) = rang()
          if( rnb(ir)<=(cs-xss(le+iq+1))/(xss(le+iq+2)
     &                 - xss(le+iq+1)))  iq = iq+1
        endif

! sample from tabulated energy distribution
        lb = id0+nint(xss(le+nm+iq+1))
        jj = nint(xss(lb-1))
        np = nint(xss(lb))
        el(m) = xss(lb+1)
        eh(m) = xss(lb+np)
        if( m/=ix )  cycle DO_110
        ir = ir+1
        if( ii==0 )  rnb(ir) = rang()
        ic = lb+2*np+1
        ib = lb+3*np
60      continue
        if( ib-ic==1 )  go to 80
        ih = (ic+ib)/2
        if( rnb(ir)>=xss(ih) ) then
          ic = ih
        else
         ib = ih
        endif
        go to 60
80      continue
        l = ic-2*np
        fa = xss(l+np)
        ea = xss(l)
        if( jj==1 )  go to 90
        bb = (xss(l+np+1)-fa)/(xss(l+1)-ea)
        if( bb==0. )  go to 90
        t = ea+(sqrt(max(zero,fa**2+2.*bb*(rnb(ir)-xss(ic))))-fa)/bb
        go to 100
90      continue
        t = ea+(rnb(ir)-xss(ic))/fa
100     continue
        acecs6_pn = t
      enddo DO_110
!
! use scaled interpolation between energies
      if( rr0==0. )  return
      t1 = el(1)+rr0*(el(2)-el(1))
      t2 = eh(1)+rr0*(eh(2)-eh(1))
      acecs6_pn = t1+(t-el(ix))*(t2-t1)/(eh(ix)-el(ix))
      return
      end function acecs6_pn
!
      subroutine rotas( c, a, d, l, ir )
! Description:
! Sample a direction d at an angle arccos(c) from axis a and
! at an azimuthal angle sampled uniformly.

! Modules used:
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'
      real*8, intent(in)    ::  a(3), c
      real*8, intent(inout) ::  d(3)
      integer,    intent(in)    ::  l
      integer,    intent(inout) ::  ir

      character(len=16) :: h
      real*8        :: r, s, t, t1, t2

      if( abs(c) > 1.0001d0) then
! Cosine too large, abort.
        write(h,'(f16.10)') c
        write(*,*) ' *** Warning: '//'Cosine ='//h//' in rotas'
        return
!
      else if ( abs(c) >= 1.0d0) then
! special handling for 1.0 <= abs(c) <= 1.0001 or more.
        d(1) = c*a(1)
        d(2) = c*a(2)
        d(3) = c*a(3)

      else
! Normal case.
        do
          t1 = 2.0d0*rang()-1.0d0
          t2 = 2.0d0*rang()-1.0d0
          r = t1**2+t2**2
          if( r<=1.0d0 ) exit
        enddo
        r = sqrt((1.0d0-c**2)/r)
        t1 = t1*r
        t2 = t2*r
!
! Two ranges for abs(c) < 1.0 .
        if( abs(a(3))<=0.9d0 ) then
          s = sqrt(a(1)**2+a(2)**2)
          t = 1.0d0/s
          d(1) = a(1)*c+(t1*a(1)*a(3)-t2*a(2))*t
          d(2) = a(2)*c+(t1*a(2)*a(3)+t2*a(1))*t
          d(3) = a(3)*c-t1*s
! renormalize every 50 calls to prevent error buildup.
          ir = ir-1
          if( ir==0 ) then
            ir = 50
            s = 1.0d0/sqrt(d(1)**2+d(2)**2+d(3)**2)
            d(1) = d(1)*s
            d(2) = d(2)*s
            d(3) = d(3)*s
          endif
!
        else
! special handling for the case of exceptionally large a(3).
          s = sqrt(a(1)**2+a(3)**2)
          t = 1.0d0/s
          d(1) = a(1)*c+(t1*a(1)*a(2)+t2*a(3))*t
          d(2) = a(2)*c-t1*s
          d(3) = a(3)*c+(t1*a(3)*a(2)-t2*a(1))*t
        endif
      endif
      return
      end subroutine rotas
