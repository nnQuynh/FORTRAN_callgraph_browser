************************************************************************
*                                                                      *
      subroutine fissin(ifssev,apr,zpr,ex,px,py,pz,
     &                  afrg,zfrg,exfrg,erfrg,
     &                  alp0,bet0,gam0,ptt0)
*                                                                      *
*                                                                      *
*       main control routine of Nakahara fission model                 *
*       last modified by K.Niita on 04/02/2000                         *
*                                                                      *
*     input:                                                           *
*                                                                      *
*        apr,zpr    : mass and proton number of mather nucleus         *
*        ex         : excitation energy of mather (MeV)                *
*        px,py,pz   : momentum vector of mather (GeV)                  *
*                                                                      *
*     output:                                                          *
*                                                                      *
*        afrg(2),zfrg(2) : mass and proton number of daughters         *
*        exfrg(2)        : excitation energies of daughters            *
*        erfrg(2)        : recoil energies of daughters                *
*        alp0,bet0,gam0  : unit vectors of recoil daughters            *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------
*     common for fission
*-----------------------------------------------------------------------

      common/barior/efb,smlan,qn,efba,efbb,smlaf
!$OMP THREADPRIVATE(/barior/)

*-----------------------------------------------------------------------

      dimension erffg(2),affg(2),zffg(2),exffg(2)
      dimension erfrg(2),afrg(2),zfrg(2),exfrg(2)
      dimension alp0(2),bet0(2),gam0(2),ptt0(2)

*-----------------------------------------------------------------------
*        fission
*-----------------------------------------------------------------------

         if( apr .lt. 50.0 ) return

*-----------------------------------------------------------------------

                  call fisbar(apr,zpr,ex)

               if( apr .lt. 225.0 ) exxx = efb
               if( apr .ge. 225.0 ) exxx = dmax1(efba,efbb)

*-----------------------------------------------------------------------

         if( ex .lt. exxx + 0.000001 ) return

*-----------------------------------------------------------------------

                  call fispr(pf,apr,zpr,ex)
                  r = unirn(dummy)

*-----------------------------------------------------------------------

         if( r .gt. pf ) return

*-----------------------------------------------------------------------

               if( zpr .lt. 90.0 ) then

                  call subfsm(affg,apr,zpr,ex)

               else

                  call acnfsm(affg,apr,zpr,ex)

               end if

*-----------------------------------------------------------------------

         if( min(affg(1),affg(2)) .lt. 6.0 ) return

*-----------------------------------------------------------------------

                  call fischg(affg,zffg,apr,zpr,ex)

*-----------------------------------------------------------------------

         if( min(zffg(1),zffg(2)) .lt. 0.0 ) return
         if( affg(1) .le. zffg(1) .or. affg(2) .le. zffg(2) ) return

*-----------------------------------------------------------------------

                  call frgen(affg,zffg,exffg,erffg,apr,zpr,ex)

*-----------------------------------------------------------------------

         if( min(exffg(1),exffg(2)) .lt. -10.0 ) return

*-----------------------------------------------------------------------
*        fission is occured
*-----------------------------------------------------------------------

                  ifssev = 2

               do k = 1, 2

                  afrg(k)  = affg(k)
                  zfrg(k)  = zffg(k)
                  exfrg(k) = exffg(k)
                  erfrg(k) = erffg(k)

               end do

                  am1 = afrg(1) * 0.93895
                  am2 = afrg(2) * 0.93895

                  pa1 = sqrt( ( erfrg(1) / 1000.0 )**2
     &                      + 2.0 * am1 * erfrg(1) / 1000.0 )

                  call gtiso(ux,uy,uz)

                  px1 = px + pa1 * ux
                  py1 = py + pa1 * uy
                  pz1 = pz + pa1 * uz

                  psq1 = sqrt( px1**2 + py1**2 + pz1**2 )

               if( psq1 .gt. 0.0d0 ) then

                  alp0(1) = px1 / psq1
                  bet0(1) = py1 / psq1
                  gam0(1) = pz1 / psq1
                  ptt0(1) = psq1

               else

                  alp0(1) = 0.0
                  bet0(1) = 0.0
                  gam0(1) = 1.0
                  ptt0(1) = 0.0

               end if

                  px2 = px - pa1 * ux
                  py2 = py - pa1 * uy
                  pz2 = pz - pa1 * uz

                  psq2 = sqrt( px2**2 + py2**2 + pz2**2 )

               if( psq2 .gt. 0.0d0 ) then

                  alp0(2) = px2 / psq2
                  bet0(2) = py2 / psq2
                  gam0(2) = pz2 / psq2
                  ptt0(2) = psq2

               else

                  alp0(2) = 0.0
                  bet0(2) = 0.0
                  gam0(2) = 1.0
                  ptt0(2) = 0.0

               end if

                  erfrg(1) = ( sqrt( psq1**2 + am1**2 ) - am1 )
     &                     * 1000.0
                  erfrg(2) = ( sqrt( psq2**2 + am2**2 ) - am2 )
     &                     * 1000.0

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine fisbar(apr,zpr,ex)
c ----------------------------------------------------------------------
c     computes fission barrior
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)

      common/barior/efb,smlan,qn,efba,efbb,smlaf
!$OMP THREADPRIVATE(/barior/)

      dimension fba(64),fbb(64)
      data fba/5.69,5.50,5.08,5.07,5.69,5.42,5.26,5.45,5.59,5.42,5.48,
     1         5.45,5.41,5.29,5.24,5.57,5.45,5.65,5.48,5.67,5.79,5.97,
     2         5.92,5.91,5.83,5.89,5.65,5.63,5.42,5.69,6.10,6.02,6.17,
     3         5.96,6.17,5.82,6.25,6.22,6.40,6.16,6.17,5.94,5.92,5.71,
     4         5.67,6.40,6.59,6.34,6.44,6.09,6.26,5.82,5.92,5.37,6.56,
     5         6.45,6.53,6.41,6.54,6.32,6.32,6.10,5.89,5.48/
      data fbb/7.89,7.68,7.30,7.34,7.35,7.14,7.04,6.58,6.79,6.68,6.80,
     1         6.84,6.86,6.79,6.80,6.26,6.21,6.48,6.38,5.75,5.95,6.20,
     2         6.23,6.29,6.28,6.40,6.23,6.26,6.12,5.21,5.69,5.68,5.93,
     3         5.79,6.08,5.79,5.34,5.39,5.65,5.48,5.46,5.41,5.52,5.32,
     4         5.34,4.87,5.15,4.98,5.16,4.89,5.13,4.77,4.94,4.45,4.50,
     5         4.38,4.54,4.50,4.72,4.57,4.65,4.50,4.36,4.02/
c ----------------------------------------------------------------------
c
    1 if(apr.gt.90.)go to 2
      efb=52.* dexp(-0.000139*(apr-90.0)**2)
      go to 100
    2 if(apr.gt.200.)go to 4
      efb=52.* dexp(-0.00008243*(apr-90.0)**2)
      go to 100
    4 if(apr.ge.225.)go to 5
      efb=23.* dexp(-0.005728*(apr -210.)**2)
      go to 100
    5 if(zpr.gt.87.9.and.zpr.lt.88.1)ifis=1
      if(zpr.gt.88.9.and.zpr.lt.89.1)ifis=5
      if(zpr.gt.89.9.and.zpr.lt.90.1)ifis=8
      if(zpr.gt.90.9.and.zpr.lt.91.1)ifis=16
      if(zpr.gt.91.9.and.zpr.lt.92.1)ifis=20
      if(zpr.gt.92.9.and.zpr.lt.93.1)ifis=30
      if(zpr.gt.93.9.and.zpr.lt.94.1)ifis=37
      if(zpr.gt.94.9.and.zpr.lt.95.1)ifis=46
      if(zpr.gt.95.1.and.zpr.lt.96.1)ifis=55
      if(ifis.ne.1)go to 6
      if(apr.gt.228.)go to 50
      ifiss=ifis+idint(apr-225.)
      go to 40
    6 if(ifis.ne.5)go to 7
      if(apr.lt.226..or.apr.gt.228)go to 50
      ifiss=ifis+idint(apr-226.)
      go to 40
    7 if(ifis.ne.8)go to 8
      if(apr.lt.227..or.apr.gt.234.)go to 50
      ifiss=ifis+idint(apr-227.)
      go to 40
    8 if(ifis.ne.16)go to 9
      if(apr.lt.230..or.apr.gt.233.)go to 50
      ifiss=ifis+idint(apr-230.)
      go to 40
    9 if(ifis.ne.20)go to 10
      if(apr.lt.231..or.apr.gt.240.)go to 50
      ifiss=ifis+idint(apr-231.)
      go to 40
   10 if(ifis.ne.30)go to 11
      if(apr.lt.233..or.apr.gt.239.)go to 50
      ifiss=ifis+idint(apr-233.)
      go to 40
   11 if(ifis.ne.37)go to 12
      if(apr.lt.237..or.apr.gt.245.)go to 50
      ifiss=ifis+idint(apr-237.)
      go to 40
   12 if(ifis.ne.46)go to 13
      if(apr.lt.239..or.apr.gt.247.)go to 50
      ifiss=ifis+idint(apr-239.)
      go to 40
   13 if(ifis.ne.55)go to 50
      if(apr.lt.241..or.apr.gt.250.)go to 50
      ifiss=ifis+idint(apr-241.)
   40 efba=fba(ifiss)
      efbb=fbb(ifiss)
      go to 100
   50 efba=6.0
      efbb=6.0
  100 return
      end


************************************************************************
*                                                                      *
      subroutine fispr(pf,apr,zpr,ex)
c ----------------------------------------------------------------------
c     computes fission probability
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)

      common/barior/efb,smlan,qn,efba,efbb,smlaf
!$OMP THREADPRIVATE(/barior/)

c ----------------------------------------------------------------------
c
    1 smlan=apr/10.
      if(apr.ge.225.) efb=dmax1(efba,efbb)
      sa=1.10107e-07
      sb=-2.12189e-04
      sc=1.11208
      smlaf=(sa*(ex**2)+sb*ex+sc)*smlan
      if(smlaf.lt.0.) smlaf=smlan
      dmino1=1.
      dmino0=0.
      qn=energy(apr-1.0,zpr)+energy(dmino1,dmino0)-energy(apr,zpr)
      xxxxf=ex-efb
      xxxxn=ex-qn
      if (xxxxf.lt.0.) then
      ex=efb
      endif
      if (xxxxn.lt.0.) then
      ex=qn
      endif
      rendf=dsqrt(smlaf*(ex-efb))
    5 rendn=dsqrt(smlan*(ex-qn))
      gratio= cbrt(apr**2)*smlaf*(ex-qn)* dexp(2.0*(rendn-rendf))/
     1  (3.05*smlan*(2.0*rendf-1.0))
   20 pf=1.0/(1.0+gratio)

  100 format(' *** warnig ex.le.efb in fispr ****')
  200 format(' *** warnig ex.le.qn in fispr ****')
      return
      end


************************************************************************
*                                                                      *
      subroutine subfsm(afrg,apr,zpr,ex)
c ----------------------------------------------------------------------
c     computes masses of fission fragments of subactinides
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)

      common/barior/efb,smlan,qn,efba,efbb,smlaf
!$OMP THREADPRIVATE(/barior/)

      dimension afrg(2)
c ----------------------------------------------------------------------
      if(apr.ge.225.)efb=0.5*(efba+efbb)
csa   changed by sasa (suggessed by nishida) 95.08.24
      extemp=ex
      if(ex.gt.50) extemp=50.
    1 width=extemp-efb+7.0
      sigma=0.8493218*width
    3 i=0
    4 xrand=anrmrn(dummy)
      xmu=0.5*apr
      xrand=dabs(xmu+sigma*xrand)
      if(xrand.ge.apr)i=i+1
      if(i.eq.1)go to 4
    8 if(xrand.ge.apr)go to 12
      afrg(1)=xrand
      afrg(2)=apr-afrg(1)
      return
   12 afrg(1)=apr
      afrg(2)=0.0
      return
      end


************************************************************************
*                                                                      *
      subroutine acnfsm(afrg,apr,zpr,ex)
c ----------------------------------------------------------------------
c     computes masses of fission fragments of actinides
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)

      common/barior/efb,smlan,qn,efba,efbb,smlaf
!$OMP THREADPRIVATE(/barior/)

      dimension afrg(2)
c ----------------------------------------------------------------------
    1 bara1=0.4*apr
      bara2=0.5*apr
      bara3=0.6*apr
      ex=ex+6.0
      if(ex.gt.25.)go to 8
      alpaa=dexp(0.5991*ex-13.1869)
      betaa=dexp(0.7013*ex-17.5325)
      go to 18
    8 if(ex.gt.40.)go to 12
      alpaa=dexp(0.2008*ex**0.8-0.8451)
      betaa=dexp(2.2672*dsqrt(ex)-11.3431)
      go to 18
   12 if(ex.gt.48.)go to 16
      alpaa=19.98159
      betaa=dexp(2.2672*dsqrt(ex)-11.3431)
      go to 18
   16 alpaa=19.98159
      betaa=78.61181
   18 ex=ex-6.0
      if(apr.ge.225.)efb=0.5*(efba+efbb)
csa   modified by sasa (suggessed by nishida)  95.08.24
      extemp=ex
      if(ex.gt.15.) extemp=15
      width=extemp-efb+7.0
      sigma=0.849321*width
      if(alpaa.le.betaa)go to 31
      xmu=0.5*(bara1-bara3)
   22 i=0
   23 xrand=anrmrn(dummy)
      xrand=dabs(xmu+sigma*xrand)
      afrg(1)=xrand+0.5*(bara1+bara3)
      if(afrg(1).ge.apr)i=i+1
      if(i.eq.1)go to 23
      if(afrg(1).ge.apr)afrg(1)=apr
      afrg(2)=apr-afrg(1)
      return
   31 xmu=bara2
   32 i=0
   33 xrand=anrmrn(dummy)
      xrand=dabs(xmu+sigma*xrand)
      afrg(1)=xrand
      if(afrg(1).ge.apr)i=i+1
      if(i.eq.1)go to 33
      if(afrg(1).ge.apr)afrg(1)=apr
      afrg(2)=apr-afrg(1)
      return
      end


************************************************************************
*                                                                      *
      subroutine fischg(afrg,zfrg,apr,zpr,ex)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)

      common/barior/efb,smlan,qn,efba,efbb,smlaf
!$OMP THREADPRIVATE(/barior/)

      dimension afrg(2),zfrg(2)
c ----------------------------------------------------------------------
c
    1 rho=1.1
      a2=0.5*apr
      zf=0.00269348*(1.0-0.625/rho)*cbrt(a2**2)
      zmp1=0.5*zpr*zf+afrg(1)*zpr*(1.0-zf)/apr
      temp=dsqrt(smlaf*(ex-efb  ))/smlaf
      ddz=503.2096*(1.0-1.4065073/cbrt(a2)+0.0148214*cbrt(apr**2))/
     1    (apr*temp)
   11 xmu=zmp1
      sigma=0.7071067/dsqrt(ddz)
      i=0
   14 xrand=anrmrn(dummy)
      xrand=dabs(xmu+sigma*xrand)
      if(xrand.ge.zpr)i=i+1
      if(i.eq.1)go to 14
      zfrg(1)=xrand
      if(zfrg(1).ge.zpr)zfrg(1)=zpr
      zfrg(2)=zpr-zfrg(1)
      return
      end


************************************************************************
*                                                                      *
      subroutine frgen(afrg,zfrg,exfrg,erfrg,
     &                 apr,zpr,ex)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)

      common/barior/efb,smlan,qn,efba,efbb,smlaf
!$OMP THREADPRIVATE(/barior/)

      dimension afrg(2),zfrg(2),exfrg(2),erfrg(2)
c ----------------------------------------------------------------------
c
    1 ekin=22.2+0.1071*zpr**2/cbrt(apr)
      erfrg(1)=afrg(2)*ekin/apr
      erfrg(2)=afrg(1)*ekin/apr
      a1=afrg(1)
    5 a2=afrg(2)
      z1=zfrg(1)
      z2=zfrg(2)
      etot=energy(apr,zpr)+ex-energy(a1,z1)-energy(a2,z2)
      exx=etot-ekin
   10 r=a1/(a1+a2)
      exfrg(1)=r*exx
      exfrg(2)=(1.0-r)*exx
      return
      end


************************************************************************
*                                                                      *
      function cbrt(xx)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c ----------------------------------------------------------------------
c
      cbrt = xx**0.33333333333d0
      return
      end


************************************************************************
*                                                                      *
      function anrmrn(dummy)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c ----------------------------------------------------------------------
c
      r = unirn(r)
      if(r.gt.0.8638) go to 10
      anrmrn=2.*(unirn(x)+unirn(y)+unirn(z)-1.5)
      return
10    if(r.gt.0.9745)go to 20
      anrmrn=1.5*(unirn(x)+unirn(y)-1.0)
      return
20    if(r.gt.0.997302039)go to 100
25    x=6.*unirn(x)-3.0
      y=0.358*unirn(x)
      xsq=x*x
      gx=17.49731196*dexp(-xsq*.5)
      ax=dabs(x)
      if(ax.gt.1.0) go to 30
      if(y.gt.(gx-17.44392294+4.73570326*xsq+2.15787544*ax))
     1go to 25
      anrmrn=x
      return
30    ax3=2.36785163*(3-ax)**2
      if(ax.gt.1.5) go to 40
      if(y.gt. (gx-ax3-2.15787544*(1.5-ax))) go to 25
      anrmrn=x
      return
40    if(y.gt.(gx-ax3)) go to 25
      anrmrn=x
      return
100   x=dsqrt(9+2*exprnf(x))
      if(unirn(x).gt.3/x) go to 100
      if(unirn(x).gt.0.5)x=-x
      anrmrn=x
      return
      end


************************************************************************
*                                                                      *
