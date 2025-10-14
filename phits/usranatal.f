************************************************************************
*                                                                      *
      subroutine usranatal(jmax,rdata,answer,rerr)
*                                                                      *
*        sample subroutine for user defined anatally                   *
*        input :                                                       *
*           jmax : number of tally output files                        *
*           rdata(2,jmax) : 1 for tally data, 2 for relative error     *
*                                                                      *
*        output :                                                      *
*           answer : output data                                       *
*           rerr   : relative error of output data                     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'angel01.inc'

      integer jmax
      real(8) :: rdata(2,jmax) 

      common /rval1/ cval(mxcval), aval(mxcval) ! you can use parameter defined in PHITS input file
      if(cval(99).gt.1.0e-10) then
       call smk_bnct(jmax,rdata,answer,rerr) ! c99 is defined
      else
       call weighted_sum(jmax,rdata,answer,rerr) ! c99[0] or not defined
      endif

      return
      end

************************************************************************
*                                                                      *
      subroutine weighted_sum(jmax,rdata,answer,rerr)
*                                                                      *
*        input :                                                       *
*           jmax : number of tally output files                        *
*           rdata(2,jmax) : 1 for tally data, 2 for relative error     *
*                                                                      *
*        output :                                                      *
*           answer : weighted sum                                      *
*           rerr   : relative error of weighted sum                    *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'angel01.inc'

      integer jmax
      real(8) :: rdata(2,jmax)

      common /rval1/ cval(mxcval), aval(mxcval) ! you can use parameter defined in PHITS input file

      sum1=0.0
      sum2=0.0  
      do j=1,jmax
       if(j.gt.98) then
        write(*,'(''only up to 98 files can be analyzed in the '',
     &  '' default setting in usranatal.f'')')
       endif 
       sum1=sum1+cval(j)*rdata(1,j)
       sum2=sum2+(cval(j)*rdata(1,j)*rdata(2,j))**2
      enddo     
      if(sum1.gt.0) then
       answer=sum1
       rerr=sqrt(sum2)/sum1
      else
       answer=0.0
       rerr=0.0
      endif

      return
      end

************************************************************************
*                                                                      *
      subroutine smk_bnct(jmax,rdata,answer,rerr)
*                                                                      *
*        outout photo isoeffective dose for BNCT based on SMK model    *
*        input :                                                       *
*           jmax : number of tally output files                        *
*           rdata(2,jmax) : 1 for tally data, 2 for relative error     *
*                                                                      *
*        output :                                                      *
*           answer : photon isoeffective dose                          *
*           rerr   : relative error of output data                     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'angel01.inc'

      parameter (ncomp = 4) ! number of dose component (should be equal to jmax)
      parameter (nGauss=200) ! number of bin for drug density Gaussian

      parameter (nmax = 32768)   ! number of data for FFT (should be 2*nmaxsqrt**2)
      parameter (nzbin= 16384)   ! number of zbin (should be nmax/2)
      parameter (nmaxsqrt = 128) ! for FFT

      common /rval1/ cval(mxcval), aval(mxcval) ! you can use parameter defined in PHITS input file

      integer jmax
      real(8) :: rdata(2,jmax) ! input tally data

      real(8) :: zdeach(ncomp+1)  ! dose-mean domain specific energy, ncomp+1 for boron data depending on drug heterogeneity
      real(8) :: zseach(ncomp+1)  ! saturation corrected specific energy
      real(8) :: zneach(ncomp+1)  ! dose-mean nucleus specific energy
      real(8) :: zfeach(ncomp+1)  ! frequency-mean nucleus specific energy = (8/9)*zn for mu-randomness, (45/64)*zn for i-randomness
      real(8) :: DnRatio(ncomp+1) ! nuclear dose ratio (D_nuclear / D_kerma
      integer :: irandomness(ncomp+1) ! -1: do not consider heterogeneity, 0: mu-randomness, 1: i-randomness
      real(8) :: xGauss(0:nGauss),GaussDist(0:nGauss) ! Drug density Gaussian distribution
      real(8) :: zhig(0:nzbin) ! zbin
      real(8) :: Fzdz(0:nzbin) ! probability density of cell-nucleus dose
      real(8) :: ram(0:ncomp+1)   ! mean number of hit for each compoenent, ram(0) is total number of hit
      save xGauss,GaussDist,zhig,zfeach

! for FFT
      integer n, ip(0 : nmaxsqrt + 1)
      real*8 a(0 : nmax), w(0 : nmax * 5 / 4 - 1), t(0 : nmax / 2) 
      real*8 Fz1(0 : nmax) ! f_1
      real*8 Fzn(0 : nmax) ! f_n
      real*8 FT1(0 : nmax) ! F{F_1}
      real*8 FTn(0 : nmax) ! F{F_n}

! BPA data in the order of photon, hydro, nitro, and boron extracellular and intracellular doses.
      data irandomness/-1,0,0,0,0/ ! -1: do not consider heterogeneity, 0: mu-randomness, 1: i-randomness
      data zdeach/2.61,42.2,51.3,178.2,178.2/  ! dose-mean domain specific energy
      data zseach/2.61,42.2,51.3,178.2,178.2/  ! same as zd because y0 is too high
      data zneach/2.86E-03,1.28E-01,2.70E-01,0.787,0.787/    ! dose mean nucleus specific energy
      data DnRatio/1.0,1.0,1.0,0.0901,0.693/ ! microscopic dose correction factor (D_nuclear / D_kerma)*N for BPA
! BSH data in the order of photon, hydro, nitro, and boron extracellular and intracellular doses.
!      data irandomness/-1,0,0,0,0/ ! -1: do not consider heterogeneity, 0: mu-randomness, 1: i-randomness
!      data zdeach/2.61,42.2,51.3,167.1,167.1/  ! dose-mean domain specific energy
!      data zseach/2.61,42.2,51.3,167.1,167.1/  ! same as zd because y0 is too high
!      data zneach/2.86E-03,1.28E-01,2.70E-01,0.651,0.651/    ! dose mean nucleus specific energy
!      data DnRatio/1.0,1.0,1.0,0.213,0.295/ ! microscopic dose correction factor (D_nuclear / D_kerma)*N for BPA

      data a0/0.0422/    ! Alpha0 parameter (Gy-1)
      data b0/0.00822/   ! Beta0 parameter (Gy-2)
      data Aref/0.0633/  ! Alpha for reference treatment
      data Xfrac/1000.0/ ! Gy/Fraction for reference treatment, more than 100 indicate RBE-weighted dose
      data scale/1.0/    ! scaling factor of dose
      data g0/4.33/       ! dose rate effect parameter (h-1)
      data Tr/0.0/       ! Irradiation time (h)
      data sigma/0.0/    ! sigma for drug density
      data imode/1/      ! FFT SMK (=1), Taylor-expansion SMK (=2), z*-based MK (=3)
      data zwid/5.0d-3/  ! zbin width

      save Gfact         ! G parameter for considering dose rate effect
      save Bref          ! beta for reference treatment, generally equal to b0 so no default value

      data ifirst/0/

!$OMP CRITICAL
      if(ifirst.eq.0) then ! initialization
       ifirst=1
       if(jmax.ne.ncomp) then
        write(*,*) 'Number of tally output files should be same as ',
     &  'number of dose component =',ncomp
        stop
       endif
       do j=1,jmax+1
        if(cval(j*10+1).gt.1e-10.and.cval(j*10+1).lt.1e5)
     &  zdeach(j)=cval(j*10+1)     ! zd (Gy) for c11,c21,c31,c41,c51
        if(cval(j*10+2).gt.1e-10.and.cval(j*10+2).lt.1e5)
     &  zseach(j)=cval(j*10+2)     ! z* (Gy) for c12,c22,c32,c42,c52
        if(cval(j*10+3).gt.1e-10.and.cval(j*10+3).lt.1e5)
     &  zneach(j)=cval(j*10+3)     ! zn (Gy) for c13,c23,c33,c43,c53
        if(irandomness(j).le.0) then   ! mu-randomness (or randomness is not considered)
         zfeach(j)=zneach(j)*8.0d0/9.0d0 ! frequency-mean nucleus specific energy
        else  ! i-randomness
         zneach(j)=zneach(j)*32.0/45.0 ! convert (zn) for u-randomness to (zn) for i-randomness
         zfeach(j)=zneach(j)*45.0/64.0 ! frequency-mean nucleus specific energy
        endif
        if(cval(j*10+4).gt.1e-10.and.cval(j*10+4).lt.1e5)
     &  DnRatio(j)=cval(j*10+4)     ! DnRatio (Gy) for c14,c24,c34,c44,c54
       enddo
       if(abs(cval(90)).gt.1e-10.and.abs(cval(90)).lt.1e20) then
        scale=abs(cval(90)) ! Scaling factor of dose
        if(cval(90).lt.0.0.and.rdata(1,1).ne.0) scale=scale/rdata(1,1)  ! negative value indicates use 1st data for normalization
       endif
       if(cval(91).gt.1e-10.and.cval(91).lt.1e5) a0=cval(91)     ! a0 (Gy-1)
       if(cval(92).gt.1e-10.and.cval(92).lt.1e5) b0=cval(92)     ! b0 (Gy-2)
       if(cval(93).gt.1e-10.and.cval(93).lt.1e5) Aref=cval(93)   ! Alpha for reference treatment (Gy-1)
       if(cval(94).eq.0.or.(cval(94).gt.1e-10.and.cval(94).lt.1e5)) then
        Bref=cval(94)   ! beta for reference treatment (Gy-1)
       else
        Bref=b0 ! basically same as b0
       endif
       if(cval(95).gt.1e-10.and.cval(95).lt.1e5) Xfrac=cval(95)  ! Gy/Fraction for reference treatment
       if(cval(96).gt.1e-10.and.cval(96).lt.1e5) sigma=cval(96)  ! sigma for drug density
       if(cval(97).gt.1e-10.and.cval(97).lt.1e5) g0=cval(97)     ! gamma0 (h-1)
       if(cval(98).gt.1e-10.and.cval(98).lt.1e5) Tr=cval(98)     ! Irradiation time (h)
       if(cval(99).gt.1e-10.and.cval(99).lt.1e5) imode=nint(cval(99))  ! FFT SMK (=1), Taylor-expansion SMK (=2), z*-based MK (=3)
       if(Tr.eq.0) then
        Gfact=1.0d0
       else
        Gfact=2.0d0*(exp(-g0*Tr)-1+g0*Tr)/((g0*Tr)**2) 
       endif

!  set up for drug distribution
       GaussDist(:)=0.0
       do ig=0,nGauss
        xGauss(ig)=ig*2.0d0/nGauss  ! x from 0 to 2.0
       enddo
       if(sigma.ne.0.0) then ! Gaussian distribution
        sum=0.0
        do ig=0,nGauss
         tmp=(xGauss(ig)-1.0d0)**2/(2*sigma**2)
         if(tmp.lt.-log(1.0e-3)) GaussDist(ig)=exp(-tmp)
         sum=sum+GaussDist(ig)
        enddo
        do ig=0,nGauss
         GaussDist(ig)=GaussDist(ig)/sum
        enddo
       else
        GaussDist(nGauss/2)=1.0 ! Delta function
       endif

! set up for zbin
       do iz=0,nzbin
        zhig(iz)=iz*zwid ! 0.01 Gy step, up to 81.92 Gy
       enddo

      endif     

!$OMP END CRITICAL

*---------- end of initialization --------------------------------------------

      answer=0.0
      rerr=0.0

      D1=0.0  ! total nucleus dose
      D2=0.0 ! square of total dose error
      do j=1,jmax
       Ratio=DnRatio(j)
       if(j.eq.jmax) Ratio=Ratio+DnRatio(j+1) ! for boron dose, sum of inner and extra-cellular dose ratio
       D1=D1+rdata(1,j)*Ratio
       D2=D2+(rdata(1,j)*rdata(2,j)*Ratio)**2
      enddo

      if(D1.eq.0.0) return ! no dose at all
      
      rerr=sqrt(D2)/D1  ! relative error is assumed to be same as that of total dose

      SF=0.0 ! total survival fraction
      do ig=0,nGauss
       if(GaussDist(ig).ne.0) then
        D=0.0 ! dose in nucleus (not kerma dose)
        Dhetero=0.0 ! heterogenous dose in nucleus
        zd=0.0
        zs=0.0
        zn=0.0
        zf=0.0
        ram(:)=0.0 ! number of hit
        do j=1,jmax+1
         if(j.ne.jmax+1) then ! photon, hydrogen, nitrogen, boron extra-cellular dose
          Ratio=DnRatio(j)
         else ! boron dose
          Ratio=DnRatio(j)*xGauss(ig)
         endif 
         zd=zd+rdata(1,min(jmax,j))*zdeach(j)*Ratio
         zs=zs+rdata(1,min(jmax,j))*zseach(j)*Ratio
         zn=zn+rdata(1,min(jmax,j))*zneach(j)*Ratio
         D=D+rdata(1,min(jmax,j))*Ratio
         if(irandomness(j).ge.0) then ! dose heterogeneity must be considered         
          Dhetero=Dhetero+rdata(1,min(jmax,j))*Ratio
          ram(j)=rdata(1,min(jmax,j))/zfeach(j)*Ratio
          ram(0)=ram(0)+ram(j)
         endif
        enddo
        if(D.eq.0.0) then
         Slog=0.0 !
        else
         zd=zd/D  ! dose-mean domain specific energy
         zs=zs/D  ! saturation-corrected speicific energy
         zn=zn/D  ! dose-mean nuclear speicific energy
         asmk=a0+b0*zs ! Alpha_smk
         bsmk=b0*(zs/zd)*Gfact ! Beta_smk
         if(imode.ge.2.or.ram(0).eq.0.0) then ! Taylor-expansion or z*-based mode
          D=D*scale ! consider scaling factor
          Slog=-asmk*D-bsmk*D**2
          if(imode.eq.2) Slog=Slog+ ! Taylor-expansion mode
     &    log(1.0+D*(-bsmk+0.5*(asmk+2*bsmk*D)**2)*zn) ! additional term of SMK model
         else ! FFT mode
          Fzdz(:)=0.0 ! set up f1(z)
          do j=1,jmax+1
           if(irandomness(j).eq.0) then ! mu-randomness
            izmax=int(zneach(j)/zwid*4.0/3.0)+1
            sum=0.0
            do iz=1,izmax
             sum=sum+iz*1.0d0/izmax
            enddo
!            tmp=ram(j)/ram(0)/((izmax+1)*izmax/2)
            do iz=1,izmax
             Fzdz(iz)=Fzdz(iz)+ram(j)/ram(0)*(iz*1.0d0/izmax)/sum
            enddo
           elseif(irandomness(j).eq.1) then ! i-randomness
            izmax=int(zneach(j)/zwid*15.0/8.0)+1
            sum=0.0
            do iz=1,izmax
             sum=sum+(1.0d0-((iz*1.0d0/izmax)**2))
            enddo
            do iz=1,izmax
             Fzdz(iz)=Fzdz(iz)+ram(j)/ram(0)*
     &       (1.0d0-(iz*1.0d0/izmax)**2)/sum
            enddo
           endif
          enddo
          sum1=0.0
          sum2=0.0
          do iz=0,nzbin ! normalized and re-evaluate ram(0)
           sum1=sum1+Fzdz(iz)
           sum2=sum2+Fzdz(iz)*zhig(iz)
          enddo
          Fzdz(:)=Fzdz(:)/sum1
          zf=sum2/sum1
          ram(0)=Dhetero/zf ! re-evaluate ram(0) because zf could be changed slightly
          D=D*scale         ! consider scaling factor
          Dhetero=Dhetero*scale ! consider scaling factor
          ram(0)=ram(0)*scale ! mean number of events is also proportional to scale factor
! Poisson and FFT analysis
          nPmax=10000
          n=nmax ! always max
          ip(0) = 0
          Fzn(:)=0.0
          do m=0,n-1,2
           Fzn(m)=Fzdz(m/2)
          enddo
          Fzdz(:)=0.0  ! from now, Fzdz is used for integral probability density
          Fzdz(0)=Poisson(ram(0),0) 
          do k=1,nPmax
           P=Poisson(ram(0),k)
           do m=0,n-1,2
            Fzdz(m/2)=Fzdz(m/2)+P*Fzn(m) ! add f_k(z) 
           enddo
           if(k.gt.ram(0).and.P.lt.1.0e-10) exit  ! not necessary anymore
           a(:)=Fzn(:)
           call rdft(n, 1, a, ip, w) ! FFT
           if(k.eq.1) FT1(:)=a(:)
           do m=0,n-1,2
            FTn(m)=a(m)*FT1(m)-a(m+1)*FT1(m+1)
            FTn(m+1)=a(m)*FT1(m+1)+a(m+1)*FT1(m)
           enddo
           a(:)=FTn(:)
           call rdft(n, -1, a, ip, w) ! Inverse FFT
           sum1=0.0
           do m=0,n-1,2
            Fzn(m)=a(m)/(n/2) 
            if(Fzn(m).lt.1.0e-10) Fzn(m)=0.0 ! ignore too small value
            sum1=sum1+Fzn(m)
            Fzn(m+1)=0.0  ! ignore imaginary part for next calculation
           enddo
           Fzn(:)=Fzn(:)/sum1 ! f_k+1(z), normalized to 1.0
           if(Fzn(n-1).ne.0.0) then
            write(*,*) 'nzbin in usranatal.f is not enough at k=',k
            stop
           endif
          enddo 
          if(k.eq.nPmax+1) then
           write(*,*) 'nPmax in usranatal.f is not enough'
           stop
          endif
          SFcell=0.0
          sum1=0.0
          sum2=0.0
          do iz=0,nzbin
           z=zhig(iz)+(D-Dhetero) ! add fixed photon component
           if(Fzdz(iz).ne.0) then
            SFcell=SFcell+Fzdz(iz)*exp(-asmk*z-bsmk*z**2)
            sum1=sum1+Fzdz(iz)
            sum2=sum2+Fzdz(iz)*z
           endif
          enddo
          SFcell=SFcell/sum1  ! normalized to 1.0
          Slog=log(SFcell) 
         endif
        endif
        SF=SF+exp(Slog)*GaussDist(ig)  ! Total survival fraction
       endif
      enddo
      if(Xfrac.ge.100) then ! RBE-weighted dose
       answer=(-Aref+sqrt(Aref**2-4*Bref*log(SF)))/(2*Bref)
      else ! EQDX
       answer=-log(SF)/(Aref+Bref*Xfrac)
      endif

      return

      end


c ****************************************************
      function Poisson(ram,k)  ! Possion distribution 
c ****************************************************
      implicit integer (i-m)
	implicit real*8 (a-h, o-z)
	
	data rammax/50.0/
	data kmax/50/
	
      if((ram.lt.k.and.ram.gt.50.0).
     1or.(ram.gt.rammax.or.k.gt.kmax)) then ! Gauss distribution approximation 
       x=k+0.5
       Poisson=exp(-((x-ram)**2/(2*ram)))/sqrt(2.0*acos(-1.0)*ram)
      else
       dk=1
       do i=1,k
        dk=dk*i
       enddo
       Poisson=ram**k*exp(-ram)/dk
      endif
      
      if(Poisson.lt.1.0e-5) Poisson=0.0
      
      return
      
      end


********************************************************************************
*    Fast Fourier Transform
*    Reference:
*    * Masatake MORI, Makoto NATORI, Tatuo TORII: Suchikeisan, 
*      Iwanamikouzajyouhoukagaku18, Iwanami, 1982 (Japanese)
*    * Henri J. Nussbaumer: Fast Fourier Transform and Convolution 
*      Algorithms, Springer Verlag, 1982
*    * C. S. Burrus, Notes on the FFT (with large FFT paper list)
*      http://www-dsp.rice.edu/research/fft/fftnote.asc

*    Copyright(C) 1996-2001 Takuya OOURA
*    download: http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html
*    You may use, copy, modify this code for any purpose and 
*    without fee. You may distribute this ORIGINAL package.
********************************************************************************

      subroutine rdft(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), nw, nc
      real*8 a(0 : n - 1), w(0 : *), xi
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. 4 * nc) then
          nc = n / 4
          call makect(nc, ip, w(nw))
      end if
      if (isgn .ge. 0) then
          if (n .gt. 4) then
              call bitrv2(n, ip(2), a)
              call cftfsub(n, a, w)
              call rftfsub(n, a, nc, w(nw))
          else if (n .eq. 4) then
              call cftfsub(n, a, w)
          end if
          xi = a(0) - a(1)
          a(0) = a(0) + a(1)
          a(1) = xi
      else
          a(1) = 0.5d0 * (a(0) - a(1))
          a(0) = a(0) - a(1)
          if (n .gt. 4) then
              call rftbsub(n, a, nc, w(nw))
              call bitrv2(n, ip(2), a)
              call cftbsub(n, a, w)
          else if (n .eq. 4) then
              call cftfsub(n, a, w)
          end if
      end if
      end

!
! -------- initializing routines --------
!
      subroutine makewt(nw, ip, w)
      integer nw, ip(0 : *), j, nwh
      real*8 w(0 : nw - 1), delta, x, y
      ip(0) = nw
      ip(1) = 1
      if (nw .gt. 2) then
          nwh = nw / 2
          delta = atan(1.0d0) / nwh
          w(0) = 1
          w(1) = 0
          w(nwh) = cos(delta * nwh)
          w(nwh + 1) = w(nwh)
          if (nwh .gt. 2) then
              do j = 2, nwh - 2, 2
                  x = cos(delta * j)
                  y = sin(delta * j)
                  w(j) = x
                  w(j + 1) = y
                  w(nw - j) = y
                  w(nw - j + 1) = x
              end do
              call bitrv2(nw, ip(2), w)
          end if
      end if
      end
!
      subroutine makect(nc, ip, c)
      integer nc, ip(0 : *), j, nch
      real*8 c(0 : nc - 1), delta
      ip(1) = nc
      if (nc .gt. 1) then
          nch = nc / 2
          delta = atan(1.0d0) / nch
          c(0) = cos(delta * nch)
          c(nch) = 0.5d0 * c(0)
          do j = 1, nch - 1
              c(j) = 0.5d0 * cos(delta * j)
              c(nc - j) = 0.5d0 * sin(delta * j)
          end do
      end if
      end
!
! -------- child routines --------
!
      subroutine bitrv2(n, ip, a)
      integer n, ip(0 : *), j, j1, k, k1, l, m, m2
      real*8 a(0 : n - 1), xr, xi, yr, yi
      ip(0) = 0
      l = n
      m = 1
      do while (8 * m .lt. l)
          l = l / 2
          do j = 0, m - 1
              ip(m + j) = ip(j) + l
          end do
          m = m * 2
      end do
      m2 = 2 * m
      if (8 * m .eq. l) then
          do k = 0, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + 2 * m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 - m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + 2 * m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              j1 = 2 * k + m2 + ip(k)
              k1 = j1 + m2
              xr = a(j1)
              xi = a(j1 + 1)
              yr = a(k1)
              yi = a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
          end do
      else
          do k = 1, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
          end do
      end if
      end
!
      subroutine bitrv2conj(n, ip, a)
      integer n, ip(0 : *), j, j1, k, k1, l, m, m2
      real*8 a(0 : n - 1), xr, xi, yr, yi
      ip(0) = 0
      l = n
      m = 1
      do while (8 * m .lt. l)
          l = l / 2
          do j = 0, m - 1
              ip(m + j) = ip(j) + l
          end do
          m = m * 2
      end do
      m2 = 2 * m
      if (8 * m .eq. l) then
          do k = 0, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + 2 * m2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 - m2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + 2 * m2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              k1 = 2 * k + ip(k)
              a(k1 + 1) = -a(k1 + 1)
              j1 = k1 + m2
              k1 = j1 + m2
              xr = a(j1)
              xi = -a(j1 + 1)
              yr = a(k1)
              yi = -a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              k1 = k1 + m2
              a(k1 + 1) = -a(k1 + 1)
          end do
      else
          a(1) = -a(1)
          a(m2 + 1) = -a(m2 + 1)
          do k = 1, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + m2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              k1 = 2 * k + ip(k)
              a(k1 + 1) = -a(k1 + 1)
              a(k1 + m2 + 1) = -a(k1 + m2 + 1)
          end do
      end if
      end
!
      subroutine cftfsub(n, a, w)
      integer n, j, j1, j2, j3, l
      real*8 a(0 : n - 1), w(0 : *)
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      l = 2
      if (n .gt. 8) then
          call cft1st(n, a, w)
          l = 8
          do while (4 * l .lt. n)
              call cftmdl(n, l, a, w)
              l = 4 * l
          end do
      end if
      if (4 * l .eq. n) then
          do j = 0, l - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = a(j + 1) + a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = a(j + 1) - a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i + x2i
              a(j2) = x0r - x2r
              a(j2 + 1) = x0i - x2i
              a(j1) = x1r - x3i
              a(j1 + 1) = x1i + x3r
              a(j3) = x1r + x3i
              a(j3 + 1) = x1i - x3r
          end do
      else
          do j = 0, l - 2, 2
              j1 = j + l
              x0r = a(j) - a(j1)
              x0i = a(j + 1) - a(j1 + 1)
              a(j) = a(j) + a(j1)
              a(j + 1) = a(j + 1) + a(j1 + 1)
              a(j1) = x0r
              a(j1 + 1) = x0i
          end do
      end if
      end
!
      subroutine cftbsub(n, a, w)
      integer n, j, j1, j2, j3, l
      real*8 a(0 : n - 1), w(0 : *)
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      l = 2
      if (n .gt. 8) then
          call cft1st(n, a, w)
          l = 8
          do while (4 * l .lt. n)
              call cftmdl(n, l, a, w)
              l = 4 * l
          end do
      end if
      if (4 * l .eq. n) then
          do j = 0, l - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = -a(j + 1) - a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = -a(j + 1) + a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i - x2i
              a(j2) = x0r - x2r
              a(j2 + 1) = x0i + x2i
              a(j1) = x1r - x3i
              a(j1 + 1) = x1i - x3r
              a(j3) = x1r + x3i
              a(j3 + 1) = x1i + x3r
          end do
      else
          do j = 0, l - 2, 2
              j1 = j + l
              x0r = a(j) - a(j1)
              x0i = -a(j + 1) + a(j1 + 1)
              a(j) = a(j) + a(j1)
              a(j + 1) = -a(j + 1) - a(j1 + 1)
              a(j1) = x0r
              a(j1 + 1) = x0i
          end do
      end if
      end
!
      subroutine cft1st(n, a, w)
      integer n, j, k1, k2
      real*8 a(0 : n - 1), w(0 : *)
      real*8 wk1r, wk1i, wk2r, wk2i, wk3r, wk3i
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      x0r = a(0) + a(2)
      x0i = a(1) + a(3)
      x1r = a(0) - a(2)
      x1i = a(1) - a(3)
      x2r = a(4) + a(6)
      x2i = a(5) + a(7)
      x3r = a(4) - a(6)
      x3i = a(5) - a(7)
      a(0) = x0r + x2r
      a(1) = x0i + x2i
      a(4) = x0r - x2r
      a(5) = x0i - x2i
      a(2) = x1r - x3i
      a(3) = x1i + x3r
      a(6) = x1r + x3i
      a(7) = x1i - x3r
      wk1r = w(2)
      x0r = a(8) + a(10)
      x0i = a(9) + a(11)
      x1r = a(8) - a(10)
      x1i = a(9) - a(11)
      x2r = a(12) + a(14)
      x2i = a(13) + a(15)
      x3r = a(12) - a(14)
      x3i = a(13) - a(15)
      a(8) = x0r + x2r
      a(9) = x0i + x2i
      a(12) = x2i - x0i
      a(13) = x0r - x2r
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(10) = wk1r * (x0r - x0i)
      a(11) = wk1r * (x0r + x0i)
      x0r = x3i + x1r
      x0i = x3r - x1i
      a(14) = wk1r * (x0i - x0r)
      a(15) = wk1r * (x0i + x0r)
      k1 = 0
      do j = 16, n - 16, 16
          k1 = k1 + 2
          k2 = 2 * k1
          wk2r = w(k1)
          wk2i = w(k1 + 1)
          wk1r = w(k2)
          wk1i = w(k2 + 1)
          wk3r = wk1r - 2 * wk2i * wk1i
          wk3i = 2 * wk2i * wk1r - wk1i
          x0r = a(j) + a(j + 2)
          x0i = a(j + 1) + a(j + 3)
          x1r = a(j) - a(j + 2)
          x1i = a(j + 1) - a(j + 3)
          x2r = a(j + 4) + a(j + 6)
          x2i = a(j + 5) + a(j + 7)
          x3r = a(j + 4) - a(j + 6)
          x3i = a(j + 5) - a(j + 7)
          a(j) = x0r + x2r
          a(j + 1) = x0i + x2i
          x0r = x0r - x2r
          x0i = x0i - x2i
          a(j + 4) = wk2r * x0r - wk2i * x0i
          a(j + 5) = wk2r * x0i + wk2i * x0r
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j + 2) = wk1r * x0r - wk1i * x0i
          a(j + 3) = wk1r * x0i + wk1i * x0r
          x0r = x1r + x3i
          x0i = x1i - x3r
          a(j + 6) = wk3r * x0r - wk3i * x0i
          a(j + 7) = wk3r * x0i + wk3i * x0r
          wk1r = w(k2 + 2)
          wk1i = w(k2 + 3)
          wk3r = wk1r - 2 * wk2r * wk1i
          wk3i = 2 * wk2r * wk1r - wk1i
          x0r = a(j + 8) + a(j + 10)
          x0i = a(j + 9) + a(j + 11)
          x1r = a(j + 8) - a(j + 10)
          x1i = a(j + 9) - a(j + 11)
          x2r = a(j + 12) + a(j + 14)
          x2i = a(j + 13) + a(j + 15)
          x3r = a(j + 12) - a(j + 14)
          x3i = a(j + 13) - a(j + 15)
          a(j + 8) = x0r + x2r
          a(j + 9) = x0i + x2i
          x0r = x0r - x2r
          x0i = x0i - x2i
          a(j + 12) = -wk2i * x0r - wk2r * x0i
          a(j + 13) = -wk2i * x0i + wk2r * x0r
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j + 10) = wk1r * x0r - wk1i * x0i
          a(j + 11) = wk1r * x0i + wk1i * x0r
          x0r = x1r + x3i
          x0i = x1i - x3r
          a(j + 14) = wk3r * x0r - wk3i * x0i
          a(j + 15) = wk3r * x0i + wk3i * x0r
      end do
      end
!
      subroutine cftmdl(n, l, a, w)
      integer n, l, j, j1, j2, j3, k, k1, k2, m, m2
      real*8 a(0 : n - 1), w(0 : *)
      real*8 wk1r, wk1i, wk2r, wk2i, wk3r, wk3i
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      m = 4 * l
      do j = 0, l - 2, 2
          j1 = j + l
          j2 = j1 + l
          j3 = j2 + l
          x0r = a(j) + a(j1)
          x0i = a(j + 1) + a(j1 + 1)
          x1r = a(j) - a(j1)
          x1i = a(j + 1) - a(j1 + 1)
          x2r = a(j2) + a(j3)
          x2i = a(j2 + 1) + a(j3 + 1)
          x3r = a(j2) - a(j3)
          x3i = a(j2 + 1) - a(j3 + 1)
          a(j) = x0r + x2r
          a(j + 1) = x0i + x2i
          a(j2) = x0r - x2r
          a(j2 + 1) = x0i - x2i
          a(j1) = x1r - x3i
          a(j1 + 1) = x1i + x3r
          a(j3) = x1r + x3i
          a(j3 + 1) = x1i - x3r
      end do
      wk1r = w(2)
      do j = m, l + m - 2, 2
          j1 = j + l
          j2 = j1 + l
          j3 = j2 + l
          x0r = a(j) + a(j1)
          x0i = a(j + 1) + a(j1 + 1)
          x1r = a(j) - a(j1)
          x1i = a(j + 1) - a(j1 + 1)
          x2r = a(j2) + a(j3)
          x2i = a(j2 + 1) + a(j3 + 1)
          x3r = a(j2) - a(j3)
          x3i = a(j2 + 1) - a(j3 + 1)
          a(j) = x0r + x2r
          a(j + 1) = x0i + x2i
          a(j2) = x2i - x0i
          a(j2 + 1) = x0r - x2r
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j1) = wk1r * (x0r - x0i)
          a(j1 + 1) = wk1r * (x0r + x0i)
          x0r = x3i + x1r
          x0i = x3r - x1i
          a(j3) = wk1r * (x0i - x0r)
          a(j3 + 1) = wk1r * (x0i + x0r)
      end do
      k1 = 0
      m2 = 2 * m
      do k = m2, n - m2, m2
          k1 = k1 + 2
          k2 = 2 * k1
          wk2r = w(k1)
          wk2i = w(k1 + 1)
          wk1r = w(k2)
          wk1i = w(k2 + 1)
          wk3r = wk1r - 2 * wk2i * wk1i
          wk3i = 2 * wk2i * wk1r - wk1i
          do j = k, l + k - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = a(j + 1) + a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = a(j + 1) - a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i + x2i
              x0r = x0r - x2r
              x0i = x0i - x2i
              a(j2) = wk2r * x0r - wk2i * x0i
              a(j2 + 1) = wk2r * x0i + wk2i * x0r
              x0r = x1r - x3i
              x0i = x1i + x3r
              a(j1) = wk1r * x0r - wk1i * x0i
              a(j1 + 1) = wk1r * x0i + wk1i * x0r
              x0r = x1r + x3i
              x0i = x1i - x3r
              a(j3) = wk3r * x0r - wk3i * x0i
              a(j3 + 1) = wk3r * x0i + wk3i * x0r
          end do
          wk1r = w(k2 + 2)
          wk1i = w(k2 + 3)
          wk3r = wk1r - 2 * wk2r * wk1i
          wk3i = 2 * wk2r * wk1r - wk1i
          do j = k + m, l + (k + m) - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = a(j + 1) + a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = a(j + 1) - a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i + x2i
              x0r = x0r - x2r
              x0i = x0i - x2i
              a(j2) = -wk2i * x0r - wk2r * x0i
              a(j2 + 1) = -wk2i * x0i + wk2r * x0r
              x0r = x1r - x3i
              x0i = x1i + x3r
              a(j1) = wk1r * x0r - wk1i * x0i
              a(j1 + 1) = wk1r * x0i + wk1i * x0r
              x0r = x1r + x3i
              x0i = x1i - x3r
              a(j3) = wk3r * x0r - wk3i * x0i
              a(j3 + 1) = wk3r * x0i + wk3i * x0r
          end do
      end do
      end
!
      subroutine rftfsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
      m = n / 2
      ks = 2 * nc / m
      kk = 0
      do j = 2, m - 2, 2
          k = n - j
          kk = kk + ks
          wkr = 0.5d0 - c(nc - kk)
          wki = c(kk)
          xr = a(j) - a(k)
          xi = a(j + 1) + a(k + 1)
          yr = wkr * xr - wki * xi
          yi = wkr * xi + wki * xr
          a(j) = a(j) - yr
          a(j + 1) = a(j + 1) - yi
          a(k) = a(k) + yr
          a(k + 1) = a(k + 1) - yi
      end do
      end
!
      subroutine rftbsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
      a(1) = -a(1)
      m = n / 2
      ks = 2 * nc / m
      kk = 0
      do j = 2, m - 2, 2
          k = n - j
          kk = kk + ks
          wkr = 0.5d0 - c(nc - kk)
          wki = c(kk)
          xr = a(j) - a(k)
          xi = a(j + 1) + a(k + 1)
          yr = wkr * xr + wki * xi
          yi = wkr * xi - wki * xr
          a(j) = a(j) - yr
          a(j + 1) = yi - a(j + 1)
          a(k) = a(k) + yr
          a(k + 1) = yi - a(k + 1)
      end do
      a(m + 1) = -a(m + 1)
      end
!
      subroutine dctsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
      m = n / 2
      ks = nc / n
      kk = 0
      do j = 1, m - 1
          k = n - j
          kk = kk + ks
          wkr = c(kk) - c(nc - kk)
          wki = c(kk) + c(nc - kk)
          xr = wki * a(j) - wkr * a(k)
          a(j) = wkr * a(j) + wki * a(k)
          a(k) = xr
      end do
      a(m) = c(0) * a(m)
      end
!
      subroutine dstsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
      m = n / 2
      ks = nc / n
      kk = 0
      do j = 1, m - 1
          k = n - j
          kk = kk + ks
          wkr = c(kk) - c(nc - kk)
          wki = c(kk) + c(nc - kk)
          xr = wki * a(k) - wkr * a(j)
          a(k) = wkr * a(k) + wki * a(j)
          a(j) = xr
      end do
      a(m) = c(0) * a(m)
      end
!





