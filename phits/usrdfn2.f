************************************************************************
*                                                                      *
      subroutine usrdfn2(ityp,ktyp,jtyp,rtyp,
     &                    dedx,dtrk,eini,efin,dhet)
*                                                                      *
*        sample subroutine for user defined let factor for heat        *
*                                                                      *
*        input :                                                       *
*           ityp : type of particle                                    *
*           ktyp : kf code of particle                                 *
*           jtyp : charge of particle                                  *
*           rtyp : mass of particle (MeV)                              *
*           dedx : dE/dx (keV/um) at sqrt(eini*efin)                   *
*           dtrk : track (cm)                                          *
*           eini : initial energy (MeV)                                *
*           efin : final energy (MeV)                                  *
*                                                                      *
*        output :                                                      *
*           dhet : dE (default) * user defined factor                  *
*                                                                      *
*                                                                      *
************************************************************************
! Calculate RBE10-weighted dose for charged particle therapy
! Reference: T.Sato et al. Radiat. Res. 171, 107-117 (2009)

      parameter (ndiv=100)
      implicit real*8 (a-h,o-z)
      include 'zdist.inc'
      real*8, save:: yb(ndiv+1)
      dimension ratio(ndiv+1)

! Parameter set for T.Sato et al. Radiat. Res. 171, 107-117 (2009)
      data rd/0.300/      ! radius of domain in um
      data A0/0.0777/      ! Alpha0 parameter
      data y0/108/    ! y0 parameter
      data beta/0.05/      ! Beta parameter
      data AlphaX/0.192/   ! Alpha for X-ray
! Parameter set for T.Sato et al. Radiat. Prot. Dosim. 143, 491-496 (2010)

      data rho/1.0/        ! density g/cm^3
      data CelDiam/0.0/    ! Cell Diameter in um, calculated in this program
      data D10Rnew/0.0/    ! D10 for X ray, calculated in this program
      data unitchange/0.0/ ! convert g/cm^3*Gy to keV/um
      data ifirst/0/  ! check first time or not

       if(ifirst.eq.0) then  ! define e-boundary
        ifirst = 1
        CelDiam=rd*2.0    ! cell diameter
        unitchange=1.0e-15/1.6e-16  ! convert g/cm^3*Gy to keV/um3
        power10=-2.0
        do ie=1,ndiv+1
         yb(ie)=10**power10
         if(power10.ge.1.0.and.power10.lt.4.0) then
          power10=power10+0.05
         else
          power10=power10+0.1
         endif
        enddo
        D10Rnew=(-AlphaX+sqrt(AlphaX**2-4*beta*log(0.1)))/(2*beta)
       endif

       unitconv  = 1.0e-3 * (2.0 / 3.0 * Celdiam) ! Q(y) is defined for keV/um
       erg  = sqrt( eini * efin )

       if(dedx.le.1.0e-5) then ! low dE/dx, RBE=1.0
        dhet = (eini - efin)
        return
       endif

       if( ityp .eq. 12 .or. ityp .eq. 13 .or.
     &     ityp .eq.  6 .or. ityp .eq.  7 ) then

          call getfitele(CelDiam,erg)

       else

          izz = abs(jtyp)
          if( izz .gt. 26 ) izz = 26

          AA = 1.0
          if( ityp .ge. 15 )
     &    AA = dble( ktyp - ktyp / 1000000 * 1000000 )

          call getfit(CelDiam,izz,erg/AA,dedx*10.0)

       end if

          do i = 1, ndiv + 1
             ratio(i) = 0.0
          end do

       do i = 1, ndiv

          zlow = unitconv *   yb(i)
          zhig = unitconv *   yb(i+1)
          zwid = unitconv * ( yb(i+1) - yb(i) )
          zmid = unitconv * ( yb(i+1) + yb(i) ) / 2.d0

          if( zmid .gt. max( 5.0e-3, Ccurr(1)*100.0) )
     &    goto 777

          if( zhig .lt. Pmid(iponly) ) then

             ratio(i) = Pcurr(1) / Pmid(iponly) * zwid

          else

             do ip = iponly, npeak

                if(zlow .le. Pmid(ip) .and.
     &             zhig .gt. Pmid(ip) ) then

                   ratio(i) = ratio(i) + Pcurr(ip)

                end if

             end do

             ratio(i) = ratio(i) + zdistfunc(zmid)*zwid

          end if

          ratio(ndiv+1) = ratio(ndiv+1) + ratio(i)

       end do

  777  if( ratio(ndiv+1) .lt. 0.5) then

          write(*,*) 'Warning: Z bin is not enough !',ratio(ndiv+1),
     &    dedx,izz,aa,erg/aa

          do i=1,ndiv+1
           write(*,*) i,ratio(i),yb(i)
          enddo
          stop

          if( ratio(ndiv+1) .eq. 0.0) ratio(ne+1)=1.0

          izpoor = izpoor + 1
          if( izpoor .gt. 100000 ) then
             write(*,*) 'Too poor of Z bin !!!'
             call parastop( 666 )
          end if

       end if

       imax=i

       iq=1

c **** Calculation of y* **************
c **** ratio(i) = d(y)dy, so yfy*wid in kase-analy.for ****
       totys = 0.0
       do i=1,min(imax,ndiv)  ! y*

          ymid = ( yb(i+1) + yb(i) )/2.0
          ywid = yb(i+1)-yb(i)

          totys = totys +
     &    (1-exp(-ymid**2/y0**2))*ratio(i)/ymid
       enddo
       totys = totys*y0**2/ratio(ndiv+1)

c ****  Calculate Alpha and RBE  ******************

       A = a0+beta/(rho*acos(-1.0)*rd**2*unitchange)*totys ! Alpha used in Bio Dose Y*

       RBE = D10Rnew*2*Beta/(-A+sqrt(A**2-4*Beta*log(0.1)))  ! Bio Dose based on Y

       dhet = (eini - efin)*RBE

! for checking purpose

*-----------------------------------------------------------------------

      return
      end

