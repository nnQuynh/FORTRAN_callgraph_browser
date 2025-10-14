
************************************************************************
*                                                                      *
      subroutine dpal3(mat,ityp,ktyp,jtyp,rtyp,
     &                 z1,a1,z2,a2,dpab,e1,e2,idpa,ieth)
*                                                                      *
*                                                                      *
*        Last Revised:     2020/12/21                                  *
*        Author:     Yosuke Iwamoto                                    *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              calculation of average dpa cross sections               *
*                                                                      *
*        Valuable:                                                     *
*        dpab: displacement cross section cm2 * traveling length cm    *
*                                                                      *
*        Projectile : z1, a1                                           *
*        Target     : z2, a2                                           *
*                                                                      *
*        References :						       *
*        Y. Iwamoto et al., Improvement of radiation damage calculation*
*	 in PHITS and tests for copper and tungsten irradiated 	       *
*	 with protons and heavy-ions over a wide energy range 	       *
*	 NIMB 274 (2012) 57-64. 				       *
*        Y. Iwamoto et al., Calculation of displacement cross-sections *
*	 for structural materials in accelerators using PHITS event    *
*	 generator and its applications to radiation damage, 	       *
*        JNST 51 (2013) 98-117. 				       *
*	 Y. Iwamoto, Implementing displacement damage  		       *
*	 calculations for electrons and gamma rays  		       *
*	 in the Particle and Heavy-Ion Transport code System,	       *
*	 NIMB 419 (2018) 32-37. 				       *
*	 Y. Iwamoto et al., Estimation of reliable 		       *
*	 displacements-per-atom based 	       			       *
*	 on athermal-recombination-corrected model in radiation        *
*	 environments at nuclear fission, fusion, 	       	       *
*	 and accelerator facilities	       			       *
*	 JNM 538 (2020) 152261.					       *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
*-----------------------------------------------------------------------
		ext2   = 0.0
		dedx1  = 0.0
		dedx2  = 0.0
		delts  = 0.0
		deltss = 0.0
		dpab   = 0.0
		dpad   = 0.0
		dpadd  = 0.0
		ext2 = (e1*e2)**0.5 + (e1-e2)*0.35

		if (ext2 .eq. 0.) then
			ext2 = 1.0d-10
		endif
*-----------------------------------------------------------------------
*	 dpa cross section for electron
*-----------------------------------------------------------------------
		if (ityp .eq. 12 .or. ityp .eq. 13) then

			call egs5ededx(e2,mat,ityp,dedx2)
        	        call egs5ededx(e1,mat,ityp,dedx1)

			if (dedx2 .le. 0. or.dedx1 .le. 0.   ) then
                     		dedx2 = 1.0d-20
                     		dedx1 = 1.0d-20
                 	endif

                      	deltss2 = e2 / dedx2
                      	deltss1 = e1 / dedx1
                      	deltss = deltss1 - deltss2

                	if (deltss.gt.1.0d+10 .or. deltss .le. 0.) then
                 		deltss=0.0
                	endif

                	call dpaele(z2,a2,ext2,dpad,idpa,ieth)
                	dpab = dpad * deltss * 1.0d+24
	                if (dpad .lt. 0. ) then
        	             dpab = 0.0
               		endif

		return
		endif

*-----------------------------------------------------------------------
*	 dpa cross section for charged particles
*-----------------------------------------------------------------------

		call rainge(e2,delts2,mat,ityp,ktyp,jtyp,rtyp)
		call rainge(e1,delts1,mat,ityp,ktyp,jtyp,rtyp)

                delts = delts1 - delts2

              	if (delts.gt.1.0d+20) then
                	delts=0.0
              	endif

		if (ityp .eq. 3 .or. ityp .eq. 5) then
                	a1 = 139.57 / 931.494
              	endif

              	if (ityp .eq. 6 .or. ityp .eq. 7) then
                	a1 = 105.659 / 931.494
              	endif

              	if (ityp .eq. 8 .or. ityp .eq. 10) then
                	a1 = 494.0 / 931.494
              	endif

               	call dpaeme(z1,a1,z2,a2,ext2,dpad,idpa,ieth)
                dpab = dpad * delts * 1.0d+24

* ----------------------------------------------------------------------

	return
	end


************************************************************************
*                                                                      *
      subroutine dpaeme(z1,am1,z2,am2,enrg,dpn,idpa,ieth)
*                                                                      *
*                                                                      *
*        Last Revised:     2020/12/21                                  *
*        Author:     Yosuke Iwamoto                                    *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              evaluate dpa cross sections through Rutherford Scatt.   *
*              by integration                                          *
*                                                                      *
*        Valuable:                                                     *
*              A: minimum value for integration                        *
*              B: maximum value for integration                        *
*              M:  number of division                                  *
*                                                                      *
*        Projectile : z1, am1                                          *
*        Target     : z2, am2                                          *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

	implicit double precision (a-h,o-z)
	data pi/3.14159265d+0/

*-----------------------------------------------------------------------
*        meshes of integration
*-----------------------------------------------------------------------

	parameter( M = 50 )

*-----------------------------------------------------------------------
*        Bohr radius [cm]
*-----------------------------------------------------------------------

	parameter(a0 = 5.3d-9)

*-----------------------------------------------------------------------
*        unit charge [esu]
*-----------------------------------------------------------------------

	parameter(esu = 4.80325d-10)

*-----------------------------------------------------------------------
*        1eV = 1.60219e-12 [erg]
*-----------------------------------------------------------------------

	parameter(eerg = 1.60219d-12)

*-----------------------------------------------------------------------

	dpn = 0.0

*-----------------------------------------------------------------------

	if( z1   .eq. 0.0 .or.
     &          am1  .le. 0.0 .or.
     &          z2   .le. 0.0 .or.
     &          am2  .le. 0.0 .or.
     &          enrg .le. 0.0 ) return

	if( z1 .le. 0.0) z1=1.0

*-----------------------------------------------------------------------
*        incident energy (eV)
*-----------------------------------------------------------------------

	e1 = enrg * 1.d6

*-----------------------------------------------------------------------
*        screening Coulomb potential
*        V(r)=z1z2/r*e^2*phi*(r/a)
*        Lindhard; Thomas Fermi model for neutral molecule [cm]
*-----------------------------------------------------------------------

	aa  = 0.8856 * a0 / ( z1**(2./3.) + z2**(2./3.) )**0.5

*-----------------------------------------------------------------------
*        for non dimensional energy [1/erg]
*-----------------------------------------------------------------------

	berg = aa * am2 / ( z1 * z2 * esu**2. * ( am1 + am2 ) )

*-----------------------------------------------------------------------
*        [1/eV]
*-----------------------------------------------------------------------

	bem = berg * eerg

*-----------------------------------------------------------------------
*        non dimensional
*-----------------------------------------------------------------------

	em = e1 * bem

*-----------------------------------------------------------------------
*        Threshold energy
*-----------------------------------------------------------------------

	call dpath(z2, ttd,ieth)

*-----------------------------------------------------------------------
*        Tmax, t=em^2 *(sin(q/2))^2=em^2 * (T/Tmax)
*-----------------------------------------------------------------------

	tmax = 4. * am1 * am2 * e1 / ( am1 + am2 )**2.

*-----------------------------------------------------------------------
*        when T = Ed,   t = td is assumed
*-----------------------------------------------------------------------

	td = em**2. * ttd / tmax

*-----------------------------------------------------------------------
*        max of t ; tdmax at theta = pai
*-----------------------------------------------------------------------

	tdmax = em**2.
	if(td.ge.tdmax) return

*-----------------------------------------------------------------------

	SSS = 0.0d0

	dx = log(tdmax / td)/ dble (M)

	A = td

	do 111 ii = 1, M

*-----------------------------------------------------------------------
*          max of integration
*-----------------------------------------------------------------------

	B = td * exp( dble(ii)*dx )

* -----------------------------------------
*      integration
* -----------------------------------------

	S = fnf(A,z2,am2,tmax,bem,em,ttd,z1,enrg,idpa,ieth)
     &     + fnf(B,z2,am2,tmax,bem,em,ttd,z1,enrg,idpa,ieth)

	area = S * (B-A) / 2.0

	SSS = SSS + area

	A = B

111     continue

* -------------------------------------------
*       displacement cross section [cm2]
* -------------------------------------------

	if (z1 .eq. 1 .and. enrg .ge. 1e1) then
		dpn = SSS
	else
		dpn = pi * aa**2. * SSS
	end if



	SSS = 0.0

	return
        end

************************************************************************
*                                                                      *
	function fnf(x,z,am,tmax,bem,em,td,z1,enrg,idpa,ieth)
*                                                                      *
*                                                                      *
*        Last Revised:     2020/12/21                                  *
*        Author:     Yosuke Iwamoto                                    *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate Coulomb scat cross section                 *
*        Valuable:                                                     *
*              x: energy                                               *
*              z: proton number of target                              *
*             am: mass number of target                                *
*             tmax:Incidence-PKA maximum energy (eV)                   *
*             bem: [1/eV]                                              *
*             em: non dimensional energy                               *
*             td: threshold energy (eV)                                *
*                                                                      *
************************************************************************

	implicit real*8(a-h,o-z)
	data pi/3.14159265d+0/

c Fitting parameter for the scattering function f(t^(1/2))
c     Thomas-Fermi
c	parameter(alamda = 1.309)
c	parameter(amm = 0.333)
c	parameter(aqq = 0.667)

c     Kr-C
c      alamda = 3.35
c      amm = 0.233
c      aqq = 0.445

c     Moliere
c      alamda = 3.07
c      amm = 0.216
c      aqq = 0.530

c     ZBL
      parameter(alamda = 5.01)
      parameter(amm = 0.203)
      parameter(aqq = 0.413)

*-----------------------------------------------------------------------
*       Selection of NRT-dpa or arc-dpa
*       idpa = 0: NRT-dpa, idpa=1: arc-dpa


*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*       Analytical approximation by Winterbon et. al.
c        a1 = alamda*x**0.1667*(1.+(2.*alamda*(x**0.667))**0.667)**(-1.5)

c       General form of f(t^(1/2))
        f1 = (1.+(2.*alamda*x**(1.-amm))**aqq)**(-1./aqq)
        a1=alamda*x**(1./2.-amm)* f1
        a2=2.*x*sqrt(x)

*-----------------------------------------------------------------------

	ak=0.1334*z**(2./3.)*am**(-0.5)

        tt = tmax/em**2.*x

c iwamoto, 2021/07/01 IEEE TRANSACTIONS ON NUCLEAR SCIENCE, VOL. 50, NO. 6, DECEMBER 2003
c        eg = tt*bem ! This is not correct!
         eg = 0.01014*z**(-7./3.)*tt

        eg34 = sqrt(sqrt(eg*eg*eg)) ! eg**(3./4.)
        g=3.4008*eg**(1./6.)+0.40244*eg34+eg
        gg=1./(1.+ak*g)

*-----------------------------------------------------------------------
*       Selection of NRT-dpa or arc-dpa
*       idpa = 0: NRT-dpa, idpa=1: arc-dpa
*-----------------------------------------------------------------------

	c = 1.0
	if(idpa.ne.0) then
		ttdam = gg*tt
		call arc(z, td, ttdam, c)
	end if

*----------------------------------------------

	fnf = a1/a2*gg*0.8*tt/(2.*td)*c

	if( tt .ge. td .and. tt .lt. 2.5*td) fnf = a1/a2
	if ( tt .lt. td)    fnf = 0.0

c iwamoto 2020/01/09, Current Opinion in Solid State
c  & Materials Science 23 (2019) 100757

*--------------------------------------------
* Coulomb relativistic for proton incidence
*--------------------------------------------
	if (z1 .eq. 1 .and. enrg .ge. 1e1) then

		rmp = 938.279
		rmn = 939.571
		av = 15.6
		as = 17.2
		ac = 0.6
		asym = 23.3
		am1 = 1.0
		am2 = am
		z1 = 1.0
		z2 = z
		bm1 = av*am2-as*am2**(2./3.)-ac*z2**2./am2**(1./3.)
		bm2 = -asym*(am2-z2-z2)**2/am2+12./am2**0.5
		bmt = bm1 + bm2
		rmtar = z2*rmp+(am2-z2)*rmn-bmt

c beta for proton
		beta = (1.-(rmp/(rmp+enrg)))**0.5

c gamma for proton
		gamma = 1./(1.-beta**2.)**0.5

c b-radius
		brad = z1*z2*197.323/137.036/rmp/beta**2*1.d-15*100.

c Tmax/4/gamma2
		pib2 = tmax/4./gamma**2.*pi*brad**2.

c alpha
		alpha = z2 / 137.036

		relfac1 = 1.-beta**2.*x/em**2.
		relfac2 = pi*alpha*beta*((x/em**2.)**0.5-x/em**2.)
		relfac = relfac1 + relfac2
		dmg = gg*0.8*tt/(2.*td)*c
		fac = tmax/em**2.

		fnf= (pib2*relfac/x**2./fac)*dmg

*------------------------------------------------------------------
		if( tt .ge. td .and. tt .lt. 2.5*td) then
			fnf = (pib2*relfac/x**2./fac)
		endif
		if ( tt .lt. td)    fnf = 0.0

c iwamoto 2020/01/09, Current Opinion in Solid State
c  & Materials Science 23 (2019) 100757

	endif

*-----------------------------------------------------------------------

	return
	end

************************************************************************
*                                                                      *
      subroutine dpaele(z2,am2,enrg,dpee,idpa,ieth)
*                                                                      *
*                                                                      *
*        Last Revised:     2021/08/31                                  *
*        Author:     Yosuke Iwamoto                                    *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*         calculate dpa cross sections for electron                    *
*              through Mckinley and Feshbach.                          *
*                                                                      *
*        Target     : z2, am2                                          *
*                                                                      *
*   James W. Corbett, Electron radiation damage                        *
*   in Semiconductors and Metals (1966)                                *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

	implicit double precision (a-h,o-z)

	data pi/3.14159265d+0/
*-----------------------------------------------------------------------
*        Clasic electron radius  [cm]
*-----------------------------------------------------------------------

	parameter(a0 = 2.81777d-13)

*-----------------------------------------------------------------------
*        electron mass [MeV]
*-----------------------------------------------------------------------

	parameter(ems = 0.511)
		dpee =0.0
		dpe=0.0

*-----------------------------------------------------------------------
*       Selection of NRT-dpa or arc-dpa
*       idpa = 0: NRT-dpa, idpa=1: arc-dpa
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        beta, gamma, alpha
*-----------------------------------------------------------------------

	beta  = (1.0 - (ems/(enrg+ems))**2.)**0.5
	gamma = (1.0 - beta**2.)**(-0.5)
	alpha = z2/137.

*-----------------------------------------------------------------------
*        Threshold energy
*-----------------------------------------------------------------------

	call dpath(z2, ttd,ieth)

*-----------------------------------------------------------------------
*        Maximum recoil energy
*-----------------------------------------------------------------------

	tmax = (560.8*(enrg/ems+2.0)*enrg/ems)/am2
	tdmx = tmax / ttd

*-----------------------------------------------------------------------
*        dpa cross section
*-----------------------------------------------------------------------

	piab = pi * alpha * beta
	bdas = 2. * z2 * a0 / beta**2. / gamma
	pib  = pi / 4. * bdas**2.

	dpe1 = tdmx - 1. - beta**2. * log(tdmx)
	dpe2 = piab * (2. * tdmx**0.5 - log(tdmx) - 2.)

	dpe = pib * (dpe1+dpe2)

        if(dpe .lt. 0.0) then
	        dpee = 0.0
	        return
        endif

*-----------------------------------------------------------------------
*        average kinetic energy of PKA and KP mpdel
*-----------------------------------------------------------------------
	tave1 = tmax*log(tdmx)-beta**2.*(tmax-ttd)
        tave2 = piab*(tmax+ttd-2.*(tmax*ttd)**0.5)
        tave  = pib * (tave1 + tave2)/dpe

c   electron energy loss limit
c   mean excitation potential
	if (z2 .lt. 13) ep = (12.+7./z2)*z2
        if (z2 .ge. 13) ep = (9.76+58.8*z2**(-1.19))*z2

        if (tave .lt. ttd) fact = 0.0
        if (tave .ge. ttd .and. tave .lt. 2.0*ttd) fact =1.0
        if (tave .ge. 2.0*ttd .and. tave .le. ep) fact = tave/2.0/ttd
        if (tave .ge. ep) fact = ep/2.0/ttd

*-----------------------------------------------------------------------
*       Selection of NRT-dpa or arc-dpa
*       idpa = 0: NRT-dpa, idpa=1: arc-dpa
*-----------------------------------------------------------------------

	c = 1.0
	if(idpa.ne.0) then
		ttdam = tave
		call arc(z2, ttd, ttdam, c)
	end if

*----------------------------------------------

	dpee = fact * dpe *c

*----------------------------------------------
        return
        end

************************************************************************
*                                                                      *
	subroutine dpaevl(z1,a1,z2,a2,enrg,dpn)
*                                                                      *
*       evaluate dpa crosssection by Lindhard-Robinson model           *
*       for anal-08.f and icntl=1                                      *
*            modified by K.Niita on 2001/04/17                         *
*       last modified by Y. Iwamoto on 2021/06/28                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
	implicit double precision (a-h,o-z)
	include 'param.inc'
	common /tall80/ itdpa(itlmax)
        common /tall81/ iteth(itlmax)
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*       Selection of NRT-dpa or arc-dpa
*       idpa = 0: NRT-dpa, idpa=1: arc-dpa
	idpa = itdpa(itnm)
        ieth = iteth(itnm)
*-----------------------------------------------------------------------

	dpn = 0.0

	if( z1 .le. 0.0 .or.
     &          a1 .le. 0.0 .or.
     &          z2 .le. 0.0 .or.
     &          a2 .le. 0.0 .or.
     &          enrg .le. 0.0 ) return

*-----------------------------------------------------------------------

	td  = 40.0
        bt  = 0.8

*-----------------------------------------------------------------------
*        Threshold energy
*-----------------------------------------------------------------------

        call dpath(z2, td,ieth)

*-----------------------------------------------------------------------

	dpn = 0.0
	ext = enrg * 1.d+6

	if( ext .gt. 0.0 ) then

		ep = 0.8853 * a2 * ext
     &               / ( 27.2 * z1 * z2
     &                 * ( z1**(2./3.) + z2**(2./3.) )**(1./2.)
     &                 * ( a1 + a2 ) )

		rk = 0.0793 * z1**(2./3.) * z2**(1./2.)
     &               * ( a1 + a2 )**(3./2.)
     &               / ( ( z1**(2./3.) + z2**(2./3.) )**(3./4.)
     &                   * a1**(3./2.) * a2**(1./2.) )

		ge = 3.4008 * ep**(1./6.)
     &               + 0.40244 * ep**(3./4.) + ep

		tdam = ext / ( 1. + rk * ge )

*-----------------------------------------------------------------------
*       Selection of NRT-dpa or arc-dpa
*       idpa = 0: NRT-dpa, idpa=1: arc-dpa
*-----------------------------------------------------------------------

		c = 1.0
		if(idpa.ne.0) then
			ttdam = tdam
			call arc(z2, td, ttdam, c)
		end if

*----------------------------------------------

		dpn = bt / 2.0 / td
     &                * tdam * c

c iwamoto 2017/05/19 for NJOY2012 comparison
		if ( tdam .lt. td) then
			dpn = 0.0
		endif

		if ( tdam .ge. td .and. tdam .lt. 2.5*td ) then
			dpn = 1.0
		endif

	end if
*-----------------------------------------------------------------------

	return
	end


************************************************************************
*                                                                      *
	subroutine dpath(z2,ttd,ieth)
*                                                                      *
*                                                                      *
*        Last Revised:     2020/12/21                                  *
*        Author:     Yosuke Iwamoto                                    *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*            determination of displacement threshold energy            *
*                                                                      *
*        Valuable:                                                     *
*           td: dpa threshold energy                                   *
*                                                                      *
*        threshold energy : ttd                                        *
*        Target     : z2                                               *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

	implicit double precision (a-h,o-z)
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*       Selection of displacement threshold energy
*       ieth = 0: before phits3.22, ieth=1: later phits3.23

*-----------------------------------------------------------------------
*        Displacement threshold energy
*        default value of previous version (before 3.22)
*-----------------------------------------------------------------------

        if (ieth .eq. 0) then

		ttd = 25.0

                if (nint(z2) .eq. 4)  ttd = 31.0
                if (nint(z2) .eq. 6)  ttd = 31.0
                if (nint(z2) .eq. 12) ttd = 25.0
                if (nint(z2) .eq. 13) ttd = 27.0
                if (nint(z2) .eq. 14) ttd = 25.0
                if (nint(z2) .eq. 20) ttd = 40.0
                if (nint(z2) .eq. 22) ttd = 40.0
                if (nint(z2) .eq. 23) ttd = 40.0
                if (nint(z2) .eq. 24) ttd = 40.0
                if (nint(z2) .eq. 25) ttd = 40.0
                if (nint(z2) .eq. 26) ttd = 40.0
                if (nint(z2) .eq. 27) ttd = 40.0
                if (nint(z2) .eq. 28) ttd = 40.0
                if (nint(z2) .eq. 29) ttd = 30.0
                if (nint(z2) .eq. 40) ttd = 40.0
                if (nint(z2) .eq. 41) ttd = 40.0
                if (nint(z2) .eq. 42) ttd = 60.0
                if (nint(z2) .eq. 47) ttd = 60.0
                if (nint(z2) .eq. 74) ttd = 90.0
                if (nint(z2) .eq. 79) ttd = 30.0
                if (nint(z2) .eq. 82) ttd = 25.0
		return
	endif

*-----------------------------------------------------------------------
*        Displacement threshold energy of PHITS version 3.23 later.
*         Fe, Cu, Ni, Pd, Pt, Ag, Au and W for the arc-dpa
*        See Table 5 on Journal of Nuclear Materials 512 (2018) 450-479
*        and Table 8 on Nuclear Energy and Technology 3 (2017) 169-175.
*        Others are obtained from Table 4 on Nuclear Energy and Technology 3 (2017) 169-175.
*-----------------------------------------------------------------------

        if (ieth .ne. 0) then

		ttd = 25.0

                if (nint(z2) .eq. 26) ttd = 40.0
                if (nint(z2) .eq. 28) ttd = 39.0
                if (nint(z2) .eq. 29) ttd = 33.0
                if (nint(z2) .eq. 46) ttd = 41.0
                if (nint(z2) .eq. 47) ttd = 39.0
                if (nint(z2) .eq. 74) ttd = 70.0
                if (nint(z2) .eq. 78) ttd = 42.0
                if (nint(z2) .eq. 79) ttd = 43.0

                if (nint(z2) .eq. 3)  ttd = 19.0
                if (nint(z2) .eq. 4)  ttd = 31.0
                if (nint(z2) .eq. 5)  ttd = 46.0
                if (nint(z2) .eq. 6)  ttd = 69.0
                if (nint(z2) .eq. 11) ttd = 17.0
                if (nint(z2) .eq. 12) ttd = 20.0
                if (nint(z2) .eq. 13) ttd = 27.0
                if (nint(z2) .eq. 14) ttd = 37.0
                if (nint(z2) .eq. 15) ttd = 20.0
                if (nint(z2) .eq. 16) ttd = 20.0
                if (nint(z2) .eq. 19) ttd = 16.0
                if (nint(z2) .eq. 20) ttd = 23.0
                if (nint(z2) .eq. 21) ttd = 33.0
                if (nint(z2) .eq. 22) ttd = 30.0
                if (nint(z2) .eq. 23) ttd = 57.0
                if (nint(z2) .eq. 24) ttd = 40.0
                if (nint(z2) .eq. 25) ttd = 33.0
                if (nint(z2) .eq. 27) ttd = 36.0
                if (nint(z2) .eq. 30) ttd = 29.0
                if (nint(z2) .eq. 31) ttd = 23.0
                if (nint(z2) .eq. 32) ttd = 35.0
                if (nint(z2) .eq. 33) ttd = 31.0
                if (nint(z2) .eq. 34) ttd = 23.0
                if (nint(z2) .eq. 35) ttd = 19.0

                if (nint(z2) .eq. 37)  ttd = 17.0
                if (nint(z2) .eq. 38)  ttd = 24.0
                if (nint(z2) .eq. 39)  ttd = 36.0
                if (nint(z2) .eq. 40)  ttd = 40.0
                if (nint(z2) .eq. 41) ttd = 78.0
                if (nint(z2) .eq. 42) ttd = 65.0
                if (nint(z2) .eq. 43) ttd = 58.0
                if (nint(z2) .eq. 44) ttd = 60.0
                if (nint(z2) .eq. 45) ttd = 51.0
                if (nint(z2) .eq. 48) ttd = 30.0
                if (nint(z2) .eq. 49) ttd = 12.0
                if (nint(z2) .eq. 50) ttd = 20.0
                if (nint(z2) .eq. 51) ttd = 22.0
                if (nint(z2) .eq. 52) ttd = 20.0
                if (nint(z2) .eq. 53) ttd = 16.0
                if (nint(z2) .eq. 55) ttd = 15.0
                if (nint(z2) .eq. 56) ttd = 22.0
                if (nint(z2) .eq. 57) ttd = 29.0
                if (nint(z2) .eq. 58) ttd = 28.0
                if (nint(z2) .eq. 59) ttd = 27.0
                if (nint(z2) .eq. 60) ttd = 28.0
                if (nint(z2) .eq. 61) ttd = 30.0
                if (nint(z2) .eq. 62) ttd = 27.0

                if (nint(z2) .eq. 63)  ttd = 17.0
                if (nint(z2) .eq. 64)  ttd = 24.0
                if (nint(z2) .eq. 65)  ttd = 36.0
                if (nint(z2) .eq. 66)  ttd = 34.0
                if (nint(z2) .eq. 67) ttd = 36.0
                if (nint(z2) .eq. 68) ttd = 37.0
                if (nint(z2) .eq. 69) ttd = 36.0
                if (nint(z2) .eq. 70) ttd = 27.0
                if (nint(z2) .eq. 71) ttd = 44.0
                if (nint(z2) .eq. 72) ttd = 61.0
                if (nint(z2) .eq. 73) ttd = 90.0
                if (nint(z2) .eq. 75) ttd = 60.0
                if (nint(z2) .eq. 76) ttd = 69.0
                if (nint(z2) .eq. 77) ttd = 58.0
                if (nint(z2) .eq. 80) ttd = 20.0
                if (nint(z2) .eq. 81) ttd = 24.0
                if (nint(z2) .eq. 82) ttd = 25.0
                if (nint(z2) .eq. 83) ttd = 23.0
                if (nint(z2) .eq. 84) ttd = 22.0
                if (nint(z2) .eq. 85) ttd = 22.0
                if (nint(z2) .eq. 87) ttd = 34.0
                if (nint(z2) .eq. 88) ttd = 24.0
                if (nint(z2) .eq. 89) ttd = 33.0
                if (nint(z2) .eq. 90) ttd = 44.0
                if (nint(z2) .eq. 91) ttd = 43.0
                if (nint(z2) .eq. 92) ttd = 39.0

		return
	endif
* ----------------------------------------------------------------------

	return
	end

************************************************************************
*                                                                      *
	subroutine arc(z, td, tdam, c)
*                                                                      *
*                                                                      *
*        Last Revised:     2020/12/21                                  *
*        Author:     Yosuke Iwamoto                                    *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*            determination of defection production efficiency          *
*                                                                      *
*        Valuable:                                                     *
*        Defection production efficiency   : c                         *
*        Target     : z                                                *
*        Threshold energy : td                                         *
*        Damage energy     : tdam                                      *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

	implicit double precision (a-h,o-z)

	barc = 0.0
	carc = 0.0

*-----------------------------------------------------------------------
* Material constants for damage production in Fe, Cu, Ni, Pd, Pt and W for the arc-dpa
* K. Nordlund et al., Journal of Nuclear Materials 512 (2018) 450-479.
*----------------------------------------------
c Fe
	if (z.eq.26) then
		barc = -0.568
		carc =  0.286
	end if
*----------------------------------------------
c Ni
	if (z.eq.28) then
		barc = -1.01
		carc =  0.23
        end if
*----------------------------------------------
c Cu
        if (z.eq.29) then
		barc = -0.68
		carc =  0.16
        end if
*----------------------------------------------
c Pd
        if (z.eq.46) then
		barc = -0.88
		carc =  0.15
        end if
*----------------------------------------------
c W
        if (z.eq.74) then
		barc = -0.56
		carc =  0.12
        end if
*----------------------------------------------
c Pt
        if (z.eq.78) then
		barc = -1.12
		carc =  0.11
        end if
*----------------------------------------------
* Material constants for damage production in Ag and Au for the arc-dpa
* Table 8 on A.Yu. Konobeyev et al. / Nuclear Energy and Technology 3 (2017) 169-175.
* Data were obtained using the method of molecular dynamics by K. Nordlund et al.
* Ag and Au data are Not indicated in JNM and Nature communication.
*----------------------------------------------
c Ag
        if (z.eq.47) then
		barc = -1.063006
		carc =  0.257325
        end if
*----------------------------------------------
c Au
        if (z.eq.79) then
		barc = -0.788966
		carc =  0.1304146
        end if
*----------------------------------------------
* Material constants for damage production for others obtained by available experimental data.
* Note that they were Not obtained by the method of molecular dynamics.
* A.Yu. Konobeyev et al. / Nuclear Energy and Technology 3 (2017) 169-175.
* the value of c-arcdpa parameter is taken from Table 7, the parameter b arcdpa is equal to -1.
*----------------------------------------------
c Li
        if (z.eq.3) then
		carc =  0.34
		barc = -1.0
	endif
c Be
        if (z.eq.4) then
		carc =  0.46
		barc = -1.0
	endif
c B
        if (z.eq.5) then
		carc =  0.58
        	barc = -1.0
	endif
c C
        if (z.eq.6) then
		carc =  0.71
        	barc = -1.0
	endif
c Na
        if (z.eq.11) then
		carc =  0.32
        	barc = -1.0
	endif
c Mg
        if (z.eq.12) then
		carc =  0.44
        	barc = -1.0
	endif
c Al
        if (z.eq.13) then
		carc =  0.44
		barc = -1.0
	endif
c Si
        if (z.eq.14) then
		carc =  0.50
        	barc = -1.0
	endif
c P
        if (z.eq.15) then
		carc =  0.36
        	barc = -1.0
	endif
c S
        if (z.eq.16) then
		carc =  0.36
        	barc = -1.0
	endif
c K
        if (z.eq.19) then
		carc =  0.33
        	barc = -1.0
	endif
c Ca
        if (z.eq.20) then
		carc =  0.41
        	barc = -1.0
	endif
c Sc
        if (z.eq.21) then
		carc =  0.53
		barc = -1.0
	endif
c Ti
        if (z.eq.22) then
		carc =  0.83
        	barc = -1.0
	endif
c V
        if (z.eq.23) then
		carc =  0.51
        	barc = -1.0
	endif
c Cr
        if (z.eq.24) then
		carc =  0.37
		barc = -1.0
	endif
c Mn
        if (z.eq.25) then
		carc =  0.33
		barc = -1.0
	endif
c Co
        if (z.eq.27) then
        	carc =  0.26
		barc = -1.0
	endif
c Zn
        if (z.eq.30) then
		carc =  0.37
		barc = -1.0
	endif
c Ga
        if (z.eq.31) then
		carc =  0.33
		barc = -1.0
	endif
c Ge
        if (z.eq.32) then
		carc =  0.43
		barc = -1.0
	endif
c As
        if (z.eq.33) then
		carc =  0.40
		barc = -1.0
	endif
c Se
        if (z.eq.34) then
		carc =  0.35
		barc = -1.0
	endif
c Br
        if (z.eq.35) then
		carc =  0.31
		barc = -1.0
	endif
c Rb
        if (z.eq.37) then
		carc =  0.31
		barc = -1.0
	endif

c Sr
        if (z.eq.38) then
		carc =  0.38
		barc = -1.0
	endif
c Y
        if (z.eq.39) then
		carc =  0.50
		barc = -1.0
	endif
c Zr
        if (z.eq.40) then
		carc =  0.70
		barc = -1.0
	endif
c Nb
        if (z.eq.41) then
		carc =  0.63
		barc = -1.0
	endif
c Mo
        if (z.eq.42) then
		carc =  0.46
		barc = -1.0
	endif
c Tc
        if (z.eq.43) then
		carc =  0.47
		barc = -1.0
	endif
c Ru
        if (z.eq.44) then
		carc =  0.44
		barc = -1.0
	endif
c Rh
        if (z.eq.45) then
		carc =  0.42
		barc = -1.0
	endif
c Cd
        if (z.eq.48) then
		carc =  0.45
		barc = -1.0
	endif
c In
        if (z.eq.49) then
		carc =  0.23
		barc = -1.0
	endif
c Sn
        if (z.eq.50) then
		carc =  0.70
		barc = -1.0
	endif
c Sb
        if (z.eq.51) then
		carc =  0.40
		barc = -1.0
	endif
c Te
        if (z.eq.52) then
		carc =  0.38
		barc = -1.0
	endif
c I
        if (z.eq.53) then
		carc =  0.33
		barc = -1.0
	endif
c Cs
        if (z.eq.55) then
		carc =  0.32
		barc = -1.0
	endif
c Ba
        if (z.eq.56) then
		carc =  0.40
		barc = -1.0
	endif
c La
        if (z.eq.57) then
		carc =  0.47
		barc = -1.0
	endif
c Ce
        if (z.eq.58) then
		carc =  0.46
		barc = -1.0
	endif
c Pr
        if (z.eq.59) then
		carc =  0.46
		barc = -1.0
	endif
c Nd
        if (z.eq.60) then
		carc =  0.46
		barc = -1.0
	endif
c Pm
        if (z.eq.61) then
		carc =  0.47
		barc = -1.0
	endif
c Sm
        if (z.eq.62) then
		carc =  0.42
		barc = -1.0
	endif
c Eu
        if (z.eq.63) then
		carc =  0.40
		barc = -1.0
	endif

c Gd
        if (z.eq.64) then
		carc =  0.49
		barc = -1.0
	endif
c Tb
        if (z.eq.65) then
		carc =  0.49
		barc = -1.0
	endif
c Dy
        if (z.eq.66) then
		carc =  0.46
		barc = -1.0
	endif
c Ho
        if (z.eq.67) then
		carc =  0.46
		barc = -1.0
	endif
c Er
        if (z.eq.68) then
		carc =  0.47
		barc = -1.0
	endif
c Tm
        if (z.eq.69) then
		carc =  0.44
		barc = -1.0
	endif
c Yb
        if (z.eq.70) then
		carc =  0.38
		barc = -1.0
	endif
c Lu
        if (z.eq.71) then
		carc =  0.50
		barc = -1.0
	endif
c Hf
        if (z.eq.72) then
		carc =  0.57
		barc = -1.0
	endif
c Ta
        if (z.eq.73) then
		carc =  0.72
		barc = -1.0
	endif
c Re
        if (z.eq.75) then
		carc =  0.87
		barc = -1.0
	endif
c Os
        if (z.eq.76) then
		carc =  0.62
		barc = -1.0
	endif
c Ir
        if (z.eq.77) then
		carc =  0.50
		barc = -1.0
	endif
c Hg
        if (z.eq.80) then
		carc =  0.28
		barc = -1.0
	endif
c Tl
        if (z.eq.81) then
		carc =  0.32
		barc = -1.0
	endif
c Pb
        if (z.eq.82) then
		carc =  0.33
		barc = -1.0
	endif
c Bi
        if (z.eq.83) then
		carc =  0.31
		barc = -1.0
	endif
c Po
        if (z.eq.84) then
		carc =  0.29
		barc = -1.0
	endif
c At
        if (z.eq.85) then
		carc =  0.30
		barc = -1.0
	endif
c Fr
        if (z.eq.87) then
		carc =  0.39
		barc = -1.0
	endif
c Ra
        if (z.eq.88) then
		carc =  0.31
		barc = -1.0
	endif
c Ac
        if (z.eq.89) then
		carc =  0.37
		barc = -1.0
	endif
c Th
        if (z.eq.90) then
		carc =  0.42
		barc = -1.0
	endif
c Pa
        if (z.eq.91) then
		carc =  0.42
		barc = -1.0
	endif
c U
        if (z.eq.92) then
		carc =  0.39
		barc = -1.0
	endif
*----------------------------------------------

	c = (1.0-carc)/(2.5*td)**barc*tdam**barc+carc
	if (c.gt. 1) then
		c = 1.0
	end if


*----------------------------------------------
	return
	end

