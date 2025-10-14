************************************************************************
*                                                                      *
      subroutine usrmgf3(xx,yy,zz,dmg,cmg,
     &                   bbx,bby,bbz,
     &                   dxx,dyx,dzx,
     &                   dxy,dyy,dzy,
     &                   dxz,dyz,dzz)
*                                                                      *
*        sample subroutine for user defined magnetic field.            *
*                                                                      *
*        Strength is given by mgf[T/m^2]                               *
*        Strength of additional field is given by gap[T]               *
*                                                                      *
*        input :                                                       *
*           xx, yy, zz    : position [cm]                              *
*           cmg           : strength of magnetic field [T/cm^2]        *
*           dmg           : additinal magnetic field [cm^2]            *
*                           dmg[cm^2] = gap[T] / cmg[T/cm^2]           *
*                                                                      *
*        output :                                                      *
*           bbx, bby, bbz : magnetic field [cm^2]                      *
*                           M[T] = cmg[T/cm^2] * bbx[cm^2]             *
*           dxx, dyx, dzx : dxx = d bbx / d x, ...                     *
*           dxy, dyy, dzy : dxy = d bbx / d y, ...                     *
*           dxz, dyz, dzz : derivative of magnetic field [cm]          *
*                           dMx/dx[T/cm] = cmg[T/cm^2] * dxx[cm]       *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------
*        example for sextapole
*-----------------------------------------------------------------------

            bbx = ( yy**2 - xx**2 ) / 2.d0
            bby = xx * yy
            bbz = dmg

            dxx = - xx
            dyx =   yy
            dzx = 0.d0

            dxy =   yy
            dyy =   xx
            dzy = 0.d0

            dxz = 0.d0
            dyz = 0.d0
            dzz = 0.d0

*-----------------------------------------------------------------------

      return
      end

