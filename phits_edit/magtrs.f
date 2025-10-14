*********************************************************************************
*                                                                      *
      subroutine elmtrs(elct,elcx,elcy,elcz,
     &                  bmgt,bmgx,bmgy,bmgz,chgp,
     &                  delt,dpr,udir,vdir,wdir,
     &                  emap_type,mmap_type,a_elmg,
     &                  t_elf,t_mgf)
*                                                                      *
*                                                                      *
*       particle transfer in void under electro magnetic field         *
*       uniform electric field and dipole magnet                       *
*       last modified by K.Niita on 2011/01/10                         *
*                                                                      *
*     input  :                                                         *
*       delt   : distance                                              *
*       elct   : electric field [MeV/nsec]                             *
*       elcx,elcy,elcz  : unit vector of electric field                *
*       bmgt   : magneric field [MeV/nsec]                             *
*       elcx,elcy,elcz  : unit vector of magnetic field                *
*       chgp   : charge state of the particle                          *
*       emap_type : type of electric data map                          *
*                   -1:xyz-list                                        *
*                   -2:rz-list                                         *
*                   -3:xyz-map                                         *
*                   -4:rz-map                                          *
*       mmap_type : type of magnetic data map                          *
*                   -1:xyz-list                                        *
*                   -2:rz-list                                         *
*                   -3:xyz-map                                         *
*                   -4:rz-map                                          *
*       a_elmg  : gap                                                  *
*       t_elf  : transform id for electoric field                      *
*       t_mgf  : transform id for magnetic field                       *
*                                                                      *
*     in common                                                        *
*                                                                      *
*       e(no)    : initial energy                                      *
*       x(no),y(no),z(no)                                              *
*                : initial position                                    *
*       u(no),v(no),w(no)                                              *
*                : initial direction                                   *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       ec(no)   : final energy                                        *
*       xc(no),yc(no),zc(no)                                           *
*                : final position                                      *
*                                                                      *
*       u(no),v(no),w(no)                                              *
*                : initial direction                                   *
*                                                                      *
*       udir,vdir,wdir                                                 *
*                : final direction                                     *
*                                                                      *
*       dpr      : final distance                                      *
*                                                                      *
* AdvanceSoft Hasemi 2019/12/21 modified for electric map                             *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA 2013/11/06
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      parameter ( pi   = 3.1415926535898d0 )
      parameter ( rlit = 29.97925d0 )

*-----------------------------------------------------------------------
*     zero electro magnetic field
*-----------------------------------------------------------------------

            if( elct .eq. 0.0d0 .and. bmgt .eq. 0.0d0 ) then

               xc(ibkxc+no,ipomp+1) = x(ibkx+no,ipomp+1) +
     &                               delt * u(ibku+no,ipomp+1)
               yc(ibkyc+no,ipomp+1) = y(ibky+no,ipomp+1) +
     &                               delt * v(ibkv+no,ipomp+1)
               zc(ibkzc+no,ipomp+1) = z(ibkz+no,ipomp+1) +
     &                               delt * w(ibkw+no,ipomp+1)

               udir = u(ibku+no,ipomp+1)
               vdir = v(ibkv+no,ipomp+1)
               wdir = w(ibkw+no,ipomp+1)

               dpr = delt

               return

            end if

*-----------------------------------------------------------------------
*     initial coordinate and direction
*-----------------------------------------------------------------------

               xx0 = x(ibkx+no,ipomp+1)
               yy0 = y(ibky+no,ipomp+1)
               zz0 = z(ibkz+no,ipomp+1)

               xx = xx0
               yy = yy0
               zz = zz0

               xp = u(ibku+no,ipomp+1)
               yp = v(ibkv+no,ipomp+1)
               zp = w(ibkw+no,ipomp+1)

               ecc = e(ibke+no,ipomp+1)
               rm0 = rtyp

*-----------------------------------------------------------------------
*      electro magnetic field
*-----------------------------------------------------------------------

            call elmgt1(xx,yy,zz,xp,yp,zp,ecc,rm0,delt,dpr,
     &                  elct,elcx,elcy,elcz,
     &                  bmgt,bmgx,bmgy,bmgz,
     &                  emap_type,mmap_type,a_elmg,
     &                  t_elf,t_mgf)

*-----------------------------------------------------------------------

               e(ibke+no,ipomp+1)   = ecc
               ec(ibkec+no,ipomp+1) = ecc

               xc(ibkxc+no,ipomp+1) = xx
               yc(ibkyc+no,ipomp+1) = yy
               zc(ibkzc+no,ipomp+1) = zz

               u(ibku+no,ipomp+1) = ( xx - xx0 ) / dpr
               v(ibkv+no,ipomp+1) = ( yy - yy0 ) / dpr
               w(ibkw+no,ipomp+1) = ( zz - zz0 ) / dpr

               udir = xp
               vdir = yp
               wdir = zp

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine elmgt1(xx,yy,zz,xp,yp,zp,ecc,rm0,delt,dpr,
     &                  elct,elcx,elcy,elcz,
     &                  bmgt,bmgx,bmgy,bmgz,
     &                  emap_type,mmap_type,a_elmg,
     &                  t_elf,t_mgf)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( rlit = 29.97925d0 )

*-----------------------------------------------------------------------

               itre = nint( t_elf )
               itrm = nint( t_mgf )

               pa0 = sqrt( ecc * ( ecc + 2.d0 * rm0 ) )

               rx0 = xx
               ry0 = yy
               rz0 = zz

               px0 = xp * pa0
               py0 = yp * pa0
               pz0 = zp * pa0

               et0 = sqrt( px0**2 + py0**2 + pz0**2 + rm0**2 )

               bx0 = px0 / et0
               by0 = py0 / et0
               bz0 = pz0 / et0

*-----------------------------------------------------------------------
               if(a_elmg .le. 0.0) then
                  nt = 10
               else
                  nt = a_elmg
               endif

               bt0 = pa0 / et0
               deltt = delt / bt0 / rlit
               dt = deltt / dble( nt )

*-----------------------------------------------------------------------

               x1 = 0
               x2 = 0
               y1 = 0
               y2 = 0
               z1 = 0
               z2 = 0
               r1 = 0
               r2 = 0
               elct0 = elct
               elcx0 = elcx
               elcy0 = elcy
               elcz0 = elcz
               bmgt0 = bmgt
               bmgx0 = bmgx
               bmgy0 = bmgy
               bmgz0 = bmgz

  100    continue

*-----------------------------------------------------------------------
*           RKG Second order
*-----------------------------------------------------------------------

               rx1 = rx0 + 0.5 * dt * bx0 * rlit
               ry1 = ry0 + 0.5 * dt * by0 * rlit
               rz1 = rz0 + 0.5 * dt * bz0 * rlit


               elct = elct0
               elcx = elcx0
               elcy = elcy0
               elcz = elcz0

               if(emap_type .ge. 0) then
                  ! elct_t = elct
               else
                  xxc = rx1
                  yyc = ry1
                  zzc = rz1
                  if(emap_type .eq. -1) then
                     call readelcmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,
     &                    elct,elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                    bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                    bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                    bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                    bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                    bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                    bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                    bxx2y2z2, byx2y2z2, bzx2y2z2)
                  else if(emap_type .eq. -2) then
                     call readelcmap2(xxc,yyc,zzc,r1,r2,z1,z2,elct,
     &                    elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                    brr2z1, bzr2z1, brr2z2, bzr2z2)
                  else if(emap_type .eq. -3) then
                     call readelcmap3(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,
     &                    elct,elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                    bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                    bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                    bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                    bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                    bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                    bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                    bxx2y2z2, byx2y2z2, bzx2y2z2)
                  else if(emap_type .eq. -4) then
                     call readelcmap4(xxc,yyc,zzc,r1,r2,z1,z2,elct,
     &                    elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                    brr2z1, bzr2z1, brr2z2, bzr2z2)
                  endif
                  elct = sqrt(elcx**2 + elcy**2 + elcz**2)
                  elct = SIGN(elct, elct0)
                  if(elct.ne.0) then ! T.Sato 2020/11/29
                   elcx1 = elcx/elct
                   elcy1 = elcy/elct
                   elcz1 = elcz/elct
                   call trnsuv(elcx1,elcy1,elcz1,elcx,elcy,elcz,itre)
                  endif
               endif


               bmgt = bmgt0
               bmgx = bmgx0
               bmgy = bmgy0
               bmgz = bmgz0

               if(mmap_type .ge. 0) then
                  ! bmgt_t = bmgt
               else
                  xxc = rx1
                  yyc = ry1
                  zzc = rz1
                  if(mmap_type .eq. -1) then
                     call readmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,
     &                    bmgt,bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                    bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                    bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                    bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                    bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                    bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                    bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                    bxx2y2z2, byx2y2z2, bzx2y2z2)
                  else if(mmap_type .eq. -2) then
                     call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,bmgt,
     &                    bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                    brr2z1, bzr2z1, brr2z2, bzr2z2)
                  else if(mmap_type .eq. -3) then
                     call readmap3(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,
     &                    bmgt,bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                    bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                    bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                    bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                    bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                    bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                    bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                    bxx2y2z2, byx2y2z2, bzx2y2z2)
                  else if(mmap_type .eq. -4) then
                     call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,bmgt,
     &                    bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                    brr2z1, bzr2z1, brr2z2, bzr2z2)
                  endif
                  bmgt = sqrt(bmgx**2 + bmgy**2 + bmgz**2)
                  bmgt = SIGN(bmgt, bmgt0)
                  if(bmgt.ne.0) then ! T.Sato 2020/11/29
                   bmgx1 = bmgx/bmgt
                   bmgy1 = bmgy/bmgt
                   bmgz1 = bmgz/bmgt
                   call trnsuv(bmgx1,bmgy1,bmgz1,bmgx,bmgy,bmgz,itrm)
                  endif
               endif

               fx0 = elct * elcx
     &             + bmgt * ( + bmgz * by0 - bmgy * bz0 )
               fy0 = elct * elcy
     &             + bmgt * ( - bmgz * bx0 + bmgx * bz0 )
               fz0 = elct * elcz
     &             + bmgt * ( + bmgy * bx0 - bmgx * by0 )

               px1 = px0 + 0.5 * dt * fx0
               py1 = py0 + 0.5 * dt * fy0
               pz1 = pz0 + 0.5 * dt * fz0

               et1 = sqrt( px1**2 + py1**2 + pz1**2 + rm0**2 )
               bx1 = px1 / et1
               by1 = py1 / et1
               bz1 = pz1 / et1

               rxf = rx0 + dt * bx1 * rlit
               ryf = ry0 + dt * by1 * rlit
               rzf = rz0 + dt * bz1 * rlit


               elct = elct0
               elcx = elcx0
               elcy = elcy0
               elcz = elcz0

               if(emap_type .ge. 0) then
                  ! elct_t = elct
               else
                  xxc = rxf
                  yyc = ryf
                  zzc = rzf
                  if(emap_type .eq. -1) then
                     call readelcmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,
     &                    elct,elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                    bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                    bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                    bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                    bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                    bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                    bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                    bxx2y2z2, byx2y2z2, bzx2y2z2)
                  else if(emap_type .eq. -2) then
                     call readelcmap2(xxc,yyc,zzc,r1,r2,z1,z2,elct,
     &                    elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                    brr2z1, bzr2z1, brr2z2, bzr2z2)
                  else if(emap_type .eq. -3) then
                     call readelcmap3(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,
     &                    elct,elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                    bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                    bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                    bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                    bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                    bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                    bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                    bxx2y2z2, byx2y2z2, bzx2y2z2)
                  else if(emap_type .eq. -4) then
                     call readelcmap4(xxc,yyc,zzc,r1,r2,z1,z2,elct,
     &                    elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                    brr2z1, bzr2z1, brr2z2, bzr2z2)
                  endif
                  elct = sqrt(elcx**2 + elcy**2 + elcz**2)
                  elct = SIGN(elct, elct0)
                  if(elct.ne.0) then ! T.Sato 2020/11/29
                   elcx1 = elcx/elct
                   elcy1 = elcy/elct
                   elcz1 = elcz/elct
                   call trnsuv(elcx1,elcy1,elcz1,elcx,elcy,elcz,itre)
                  endif
               endif


               bmgt = bmgt0
               bmgx = bmgx0
               bmgy = bmgy0
               bmgz = bmgz0

               if(mmap_type .ge. 0) then
                  ! bmgt_t = bmgt
               else
                  xxc = rxf
                  yyc = ryf
                  zzc = rzf
                  if(mmap_type .eq. -1) then
                     call readmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,
     &                    bmgt,bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                    bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                    bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                    bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                    bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                    bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                    bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                    bxx2y2z2, byx2y2z2, bzx2y2z2)
                  else if(mmap_type .eq. -2) then
                     call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,bmgt,
     &                    bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                    brr2z1, bzr2z1, brr2z2, bzr2z2)
                  else if(mmap_type .eq. -3) then
                     call readmap3(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,
     &                    bmgt,bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                    bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                    bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                    bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                    bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                    bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                    bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                    bxx2y2z2, byx2y2z2, bzx2y2z2)
                  else if(mmap_type .eq. -4) then
                     call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,bmgt,
     &                    bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                    dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                    brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                    brr2z1, bzr2z1, brr2z2, bzr2z2)
                  endif
                  bmgt = sqrt(bmgx**2 + bmgy**2 + bmgz**2)
                  bmgt = SIGN(bmgt, bmgt0)
                  if(bmgt.ne.0) then ! T.Sato 2020/11/29
                   bmgx1 = bmgx/bmgt
                   bmgy1 = bmgy/bmgt
                   bmgz1 = bmgz/bmgt
                   call trnsuv(bmgx1,bmgy1,bmgz1,bmgx,bmgy,bmgz,itrm)
                  endif
               endif

               fx1 = elct * elcx
     &             + bmgt * ( + bmgz * by1 - bmgy * bz1 )
               fy1 = elct * elcy
     &             + bmgt * ( - bmgz * bx1 + bmgx * bz1 )
               fz1 = elct * elcz
     &             + bmgt * ( + bmgy * bx1 - bmgx * by1 )

               pxf = px0 + dt * fx1
               pyf = py0 + dt * fy1
               pzf = pz0 + dt * fz1

*-----------------------------------------------------------------------

               dpr = dsqrt( ( rxf - xx )**2
     &                    + ( ryf - yy )**2
     &                    + ( rzf - zz )**2 )

            if( abs( dpr - delt ) / delt .lt. 0.001 ) goto 200

            if( dpr .gt. delt ) then

               dt = dt / 2.0

            else if( dpr .lt. delt ) then

               rx0 = rxf
               ry0 = ryf
               rz0 = rzf

               px0 = pxf
               py0 = pyf
               pz0 = pzf

               et0 = sqrt( px0**2 + py0**2 + pz0**2 + rm0**2 )

               bx0 = px0 / et0
               by0 = py0 / et0
               bz0 = pz0 / et0

            end if

               goto 100

*-----------------------------------------------------------------------

  200    continue

               elct = elct0
               elcx = elcx0
               elcy = elcy0
               elcz = elcz0
               bmgt = bmgt0
               bmgx = bmgx0
               bmgy = bmgy0
               bmgz = bmgz0

               xx = rxf
               yy = ryf
               zz = rzf

               pf2 = pxf**2 + pyf**2 + pzf**2
               pf0 = sqrt( pf2 )
               et0 = sqrt( pf2 + rm0**2 )

               ecc = et0 - rm0

               xp = pxf / pf0
               yp = pyf / pf0
               zp = pzf / pf0

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine magtrs(a_mag,b_mag,s_mag,p_mag,t_mag,
     &                  delt,dpr,udir,vdir,wdir,mark)
*                                                                      *
*                                                                      *
*       particle transfer in void under magnetic field                 *
*       created  by S.Meigo on 2000/07/01                              *
*       last modified by K.Niita on 2010/12/23                         *
*                                                                      *
*     input  :                                                         *
*       delt   : distance                                              *
*       a_mag  : magnet gap(mm)                                        *
*       b_mag  : magnet field at pole tip  [kG]                        *
*       s_mag  : speicies of magnet dypole:2, quad:4, sext:6, oct:8    *
*       p_mag  : phase of magnet                                       *
*       t_mag  : transform id                                          *
*                                                                      *
*     in common                                                        *
*                                                                      *
*       e(no)    : initial energy                                      *
*       x(no),y(no),z(no)                                              *
*                : initial position                                    *
*       u(no),v(no),w(no)                                              *
*                : initial direction                                   *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       ec(no)   : final energy                                        *
*       xc(no),yc(no),zc(no)                                           *
*                : final position                                      *
*                                                                      *
*       u(no),v(no),w(no)                                              *
*                : initial direction                                   *
*                                                                      *
*       udir,vdir,wdir                                                 *
*                : final direction                                     *
*                                                                      *
*       dpr      : final distance                                      *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /gravit/ grav(3), igrav

      parameter ( alph = 5.7688252d0 )
      parameter ( gamn = 1.8324712d+8 )
      parameter ( grvt = 980.665d0 )
      parameter ( cvel = 2.997925d+10 )
      parameter ( pi   = 3.1415926535898d0 )
      parameter ( rlit = 29.97925d0 )
      parameter ( grvc = 980.665d-18 )

      integer nn

      vxx = 0
      vyy = 0
      vzz = 0
      vel = 0

      sx = 0d0 ! S.H. initialization (2020.2.17)
      sy = 0d0
      sz = 0d0

*-----------------------------------------------------------------------
*     transform
*-----------------------------------------------------------------------

            itrs = nint( t_mag )

            call trnsxx(x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  xx0,yy0,zz0,itrs)

            call trnsuu(u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  uu0,vv0,ww0,itrs)

*-----------------------------------------------------------------------
*     gravity
*-----------------------------------------------------------------------

               jgrav = 0

         if( igrav .ne. 0 .and. ityp .eq. 2 .and.
     &       e(ibke+no,ipomp+1) .gt. 0.0 .and.
     &                    e(ibke+no,ipomp+1) .le. 1.d-6 ) then

               call trnsuu(grav(1),grav(2),grav(3),
     &                     gravx,gravy,gravz,itrs)

               jgrav = 1

         end if

*-----------------------------------------------------------------------
*     zero magnet field
*-----------------------------------------------------------------------

            if( b_mag .eq. 0.0d0  ) goto 1000

*-----------------------------------------------------------------------
*     which magnetic field and method
*-----------------------------------------------------------------------

            isp = nint( s_mag )

            if( isp .eq. 60 ) then

               isp = 6
               is6 = 0

            else if( isp .eq. 61 ) then

               isp = 6
               is6 = 1

            else if( isp .eq. 62 ) then

               isp = 100
               is6 = 6

            else if( isp .eq. 101 ) then

               isp = 100
               is6 = 1

            else if( isp .eq. 102 ) then

               isp = 100
               is6 = 2

            else if( isp .eq. 103 ) then

               isp = 100
               is6 = 3

            else if( isp .eq. 104 ) then

               isp = 100
               is6 = 4

            else if( isp .eq. 106 ) then

               isp = 100
               is6 = 6

            else if( isp .eq. -1 ) then
               isp = -1
               ineu = 0
               imap = 1

            else if( isp .eq. -101 ) then
               isp = -1
               ineu = 1
               imap = 1

            else if( isp .eq. -2 ) then
               isp = -1
               ineu = 0
               imap = 2

            else if( isp .eq. -102 ) then
               isp = -1
               ineu = 1
               imap = 2

            else if( isp .eq. -3 ) then
               isp = -1
               ineu = 0
               imap = 3

            else if( isp .eq. -103 ) then
               isp = -1
               ineu = 1
               imap = 3

            else if( isp .eq. -4 ) then
               isp = -1
               ineu = 0
               imap = 4

            else if( isp .eq. -104 ) then
               isp = -1
               ineu = 1
               imap = 4

            else
            end if

*-----------------------------------------------------------------------
*     drift
*-----------------------------------------------------------------------

            if( ( isp .eq. 4 .or. isp .eq. 6 ) .and.
     &            abs(ww0) .lt. 1.0d-5 ) goto 1000

            if( isp .eq. 2 .and.
     &          uu0**2 + ww0**2 .lt. 1.0d-5 ) goto 1000

*-----------------------------------------------------------------------
*      magnet field
*-----------------------------------------------------------------------

            if( isp .eq. 2 ) then

                  b = b_mag
                  p = dsqrt(e(ibke+no,ipomp+1)**2+
     &                             2.*rtyp*e(ibke+no,ipomp+1))/1000.

                  cg = ctyp
                  r0 = 33.356 * p / abs( cg * b ) * 100.d0

                  aw0 = dsqrt( uu0**2 + ww0**2 )

                  icm = -1
                  if( cg * b .lt. 0.0d0 ) icm = 1

               if( delt .gt. 2.0 * r0 / aw0 ) then

                  delt = 2.0 * r0 / aw0

               end if

               if( abs(ww0) .gt. 0.0d0 ) then

                  the2 = atan( uu0 / ww0 )

                  if( ww0 .lt. 0.0 ) the2 = the2 + pi

               else

                  the2 = pi / 2.0d0

                  if( uu0 .lt. 0.0 ) the2 = the2 + pi

               end if

                  cst = cos( the2 )
                  snt = sin( the2 )

*-----------------------------------------------------------------------

            else if( isp .eq. 4 ) then

                  a = a_mag
                  b = b_mag
                  p = dsqrt(e(ibke+no,ipomp+1)**2+
     &                        2.*rtyp*e(ibke+no,ipomp+1))/1000.

                  cg = ctyp

                  uu = uu0 / ww0
                  vv = vv0 / ww0

                  xp0 = atan(uu)
                  yp0 = atan(vv)

*-----------------------------------------------------------------------

            else if( isp .eq. 6 ) then

                  omg = sqrt( alph * abs( b_mag ) )
                  sfc = gamn * omg / alph / 10000.0
                  grv = grvt / omg**2
                  vel = sqrt( 2. * e(ibke+no,ipomp+1) / rtyp ) * cvel

                  xp0 = vel / omg * uu0
                  yp0 = vel / omg * vv0
                  zp0 = vel / omg * ww0

*-----------------------------------------------------------------------

            else if( isp .eq. 100 ) then

                  omg = sqrt( alph * abs( b_mag ) )
                  sfc = gamn * omg / alph / 10000.0
                  grv = grvt / omg**2
                  vel = sqrt( 2. * e(ibke+no,ipomp+1) / rtyp ) * cvel

                  cmg = abs( b_mag ) / 10000.d0
                  dmg = a_mag / cmg

                  xp0 = vel / omg * uu0
                  yp0 = vel / omg * vv0
                  zp0 = vel / omg * ww0

*-----------------------------------------------------------------------
*           ASTOM 2018/10/25
            else if( isp .eq. -1 ) then
*                 verocity of particle
! T.Sato 2020/10/04 consider relativistic
                  gmm = (rtyp + e(ibke+no,ipomp+1)) / rtyp
                  vel = sqrt( gmm**2 - 1. ) / gmm
                  vxx = vel * uu0
                  vyy = vel * vv0
                  vzz = vel * ww0
                  if( ineu .eq. 0 ) then
*                     charge
                      cg = ctyp
                      qm = cg / rtyp / gmm ! T.Sato 2020/10/04 consider relativistic
                      ssxa = 0
                      ssya = 0
                      ssza = 0
                  else if( ineu .eq. 1 ) then
*                     spin
                      ssxa = 0
                      ssya = 0
                      ssza = 0
                      call trnsuu(spx(ibkspx+no,ipomp+1),
     &                            spy(ibkspy+no,ipomp+1),
     &                            spz(ibkspz+no,ipomp+1),
     &                            ssxa,ssya,ssza,itrs)

                  else
                  end if
            end if

*-----------------------------------------------------------------------

         deltv = delt

         icc = 0
 5000    icc = icc + 1

*-----------------------------------------------------------------------

            if( isp .eq. 2 ) then

                  dll = deltv * aw0
                  al1 = dll / r0

                  dpr = sqrt( ( 2.0 * r0 * sin(al1/2.0) )**2
     &                      + ( vv0 / aw0 * dll )**2 )

*-----------------------------------------------------------------------

            else if( isp .eq. 4 ) then

                  xx = xx0
                  yy = yy0
                  xp = xp0
                  yp = yp0
                  dz = deltv * ww0

                  call quad(p,b,a,dz,xx,xp,yy,yp,cg)

                  xxc = xx
                  yyc = yy
                  zzc = zz0 + dz

                  dpr = dsqrt( ( xx0 - xxc )**2
     &                       + ( yy0 - yyc )**2
     &                       + ( zz0 - zzc )**2 )

*-----------------------------------------------------------------------
            else if( isp .eq. -1 ) then
                  xx = xx0
                  yy = yy0
                  zz = zz0
                  vx = vxx
                  vy = vyy
                  vz = vzz
                  nn = nint(100 * a_mag)
                  tt = (deltv / vel) / nn
                  x1 = 0
                  x2 = 0
                  y1 = 0
                  y2 = 0
                  z1 = 0
                  z2 = 0
                  r1 = 0
                  r2 = 0
                  do i = 1, nn
                  call magmap(xx,yy,zz,vx,vy,vz,tt,qm,
     &                   ssxa,ssya,ssza,alph,gamn,
     &                   ineu,imap,b_mag,
     &                   x1,x2,y1,y2,z1,z2,r1,r2,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2,mark)
                   if(mark.eq.-2) return ! T.Sato 2022/03/22
                  end do


                  xxc = xx
                  yyc = yy
                  zzc = zz

                  dpr = dsqrt( ( xx0 - xxc )**2
     &                       + ( yy0 - yyc )**2
     &                       + ( zz0 - zzc )**2 )
*-----------------------------------------------------------------------

            else if( isp .eq. 6 ) then

                     xx = xx0
                     xp = xp0
                     yy = yy0
                     yp = yp0
                     zz = zz0
                     zp = zp0

                     dz = deltv * ww0
                     thet = dz / zp0

*-----------------------------------------------------------------------
*                 initial spin
*-----------------------------------------------------------------------

                     ssr = spx(ibkspx+no,ipomp+1)**2
     &                   + spy(ibkspy+no,ipomp+1)**2
     &                   + spz(ibkspz+no,ipomp+1)**2

                  if( ssr .lt. 1.d-8 ) then

                        pol = p_mag
                        sn = 0.d0

                     if( pol .lt. -1.d0 ) then
                        pol  =  0.d0
                     else if( pol .gt. 1.d0 ) then
                        pol  =  100.d0
                        ipal = 0
                     end if

                     if( pol .ge. -1.d0 .and. pol .le. 1.d0 ) then

                        pol = ( pol + 1.d0 ) / 2.d0
                        ipal = 1
                        if( unirn(dummy) .gt. pol ) ipal = -1

                        sn = ( xx**2 + yy**2 ) / 2.d0

                     end if

                     if( sn .gt. 0.0d0 .and. ipal .ne. 0 ) then

                        sx = ipal * ( yy**2 - xx**2 ) / 2.d0 / sn
                        sy = ipal *   xx * yy  / sn
                        sz = 0.0d0

                     else

                        th = 2.0d0 * pi * unirn(dummy)
                        cs = 2.d0 *  unirn(dummy) - 1.d0
                        sn = sqrt( 1.d0 - cs**2 )
                        sx = sn * cos( th )
                        sy = sn * sin( th )
                        sz = cs

                     end if

                  else

                       call trnsuu(spx(ibkspx+no,ipomp+1),
     &                             spy(ibkspy+no,ipomp+1),
     &                             spz(ibkspz+no,ipomp+1),sx,sy,sz,itrs)

                  end if

*-----------------------------------------------------------------------

               if( is6 .eq. 0 ) then

                     call sexts0(dz,xx,xp,yy,yp,sx,sy,sz,thet)

                     xxc = xx
                     yyc = yy
                     zzc = zz + dz

                     dpr = dsqrt( ( xx0 - xxc )**2
     &                          + ( yy0 - yyc )**2
     &                          + ( zz0 - zzc )**2 )

               else if( is6 .eq. 1 ) then

                  call sexts1(xx,xp,yy,yp,zz,zp,sx,sy,sz,
     &                        thet,delt,dpr,sfc,
     &                        grv,gravx,gravy,gravz)

                     xxc = xx
                     yyc = yy
                     zzc = zz

               end if

*-----------------------------------------------------------------------
*        special version
*-----------------------------------------------------------------------

            else if( isp .eq. 100 ) then

                     xx = xx0
                     xp = xp0
                     yy = yy0
                     yp = yp0
                     zz = zz0
                     zp = zp0

                     dz = deltv * ww0
                     thet = dz / zp0

*-----------------------------------------------------------------------
*                 initial spin
*-----------------------------------------------------------------------

                     ssr = spx(ibkspx+no,ipomp+1)**2
     &                   + spy(ibkspy+no,ipomp+1)**2
     &                   + spz(ibkspz+no,ipomp+1)**2

                  if( ssr .lt. 1.d-8 ) then

                        pol = p_mag
                        bba = 0.d0

                     if( pol .lt. -1.d0 ) then
                        pol  =  0.d0
                     else if( pol .gt. 1.d0 ) then
                        pol  =  100.d0
                        ipal = 0
                     end if

                     if( pol .ge. -1.d0 .and. pol .le. 1.d0 ) then

                        pol = ( pol + 1.d0 ) / 2.d0
                        ipal = 1
                        if( unirn(dummy) .gt. pol ) ipal = -1

                        call magbdx(is6,xx,yy,zz,dmg,cmg,
     &                              bbx,bby,bbz,bba,
     &                              dxx,dyx,dzx,
     &                              dxy,dyy,dzy,
     &                              dxz,dyz,dzz)

                     end if

                     if( bba .gt. 0.0d0 .and. ipal .ne. 0 ) then

                        sx = ipal * bbx / bba
                        sy = ipal * bby / bba
                        sz = ipal * bbz / bba

                     else

                        th = 2.0d0 * pi * unirn(dummy)
                        cs = 2.d0 *  unirn(dummy) - 1.d0
                        sn = sqrt( 1.d0 - cs**2 )
                        sx = sn * cos( th )
                        sy = sn * sin( th )
                        sz = cs

                     end if

                  else

                        call trnsuu(spx(ibkspx+no,ipomp+1),
     &                              spy(ibkspy+no,ipomp+1),
     &                             spz(ibkspz+no,ipomp+1),sx,sy,sz,itrs)

                  end if

*-----------------------------------------------------------------------

                  call sexts2(xx,xp,yy,yp,zz,zp,sx,sy,sz,
     &                        thet,delt,dpr,sfc,dmg,cmg,
     &                        grv,gravx,gravy,gravz,is6)

                     xxc = xx
                     yyc = yy
                     zzc = zz

            end if

*-----------------------------------------------------------------------
*           adjust the distanc to delt
*-----------------------------------------------------------------------

            if( icc .le. 100 .and.
     &          abs( dpr - delt ) / delt .gt. 0.001 ) then

               factv = ( delt - dpr ) / delt
               factv = sign( min( 1.0d0, abs( factv ) ), factv )

               deltv = deltv + deltv * factv * 0.9

               goto 5000

            end if

*-----------------------------------------------------------------------

            if( isp .eq. 2 ) then

               cs = cos( al1 )
               sn = sin( al1 )
               zz = r0 * sn
               xx = ( r0 - r0 * cs ) * icm

               zzc = zz0 + zz * cst - xx * snt
               xxc = xx0 + zz * snt + xx * cst
               yyc = yy0 + dll * vv0 / aw0

               zz = cs
               xx = sn * icm

               vz = ( zz * cst - xx * snt ) * aw0
               vx = ( zz * snt + xx * cst ) * aw0
               vy = vv0

               av = dsqrt( vx**2 + vy**2 + vz**2 )

               uur = vx / av
               vvr = vy / av
               wwr = vz / av

            else if( isp .eq. 4 ) then

               vx = ww0 * tan( xp )
               vy = ww0 * tan( yp )
               vz = ww0 * 1.0d0

               av = dsqrt( vx**2 + vy**2 + vz**2 )

               uur = vx / av
               vvr = vy / av
               wwr = vz / av

            else if( isp .eq. 6 ) then

               av = dsqrt( xp**2 + yp**2 + zp**2 )
               vl = av * omg

               ec(ibkec+no,ipomp+1) = 0.5 * rtyp * ( vl / cvel )**2

               uur = xp / av
               vvr = yp / av
               wwr = zp / av

            else if( isp .eq. 100 ) then

               av = dsqrt( xp**2 + yp**2 + zp**2 )
               vl = av * omg

               ec(ibkec+no,ipomp+1) = 0.5 * rtyp * ( vl / cvel )**2

               uur = xp / av
               vvr = yp / av
               wwr = zp / av

               sv = dsqrt( sx**2 + sy**2 + sz**2 )

               sx = sx / sv
               sy = sy / sv
               sz = sz / sv

*           ASTOM 2018/10/25
            else if( isp .eq. -1 ) then

                av = dsqrt( vx**2 + vy**2 + vz**2 )

                uur = vx / av
                vvr = vy / av
                wwr = vz / av

                if( ineu .eq. 1 ) then
                     ssr = dsqrt(ssxa**2 + ssya**2 + ssza**2)
                     if( ssr .gt. 1.d-8 ) then
                         sx = ssxa / ssr
                         sy = ssya / ssr
                         sz = ssza / ssr
                     else
                         sx = 0.d0
                         sy = 0.d0
                         sz = 0.d0
                     end if

                 else
                 end if

            end if

               uuc = ( xxc - xx0 ) / dpr
               vvc = ( yyc - yy0 ) / dpr
               wwc = ( zzc - zz0 ) / dpr

*-----------------------------------------------------------------------
*        transform
*-----------------------------------------------------------------------

            call trnsxv(xxc,yyc,zzc,
     &                  xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1),itrs)

            call trnsuv(uuc,vvc,wwc,
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),itrs)

            call trnsuv(uur,vvr,wwr,
     &                  udir,vdir,wdir,itrs)

            call trnsuv(sx,sy,sz,
     &           spx(ibkspx+no,ipomp+1),spy(ibkspy+no,ipomp+1),
     &           spz(ibkspz+no,ipomp+1),itrs)


            return

*-----------------------------------------------------------------------

 1000 continue

               xc(ibkxc+no,ipomp+1) = x(ibkx+no,ipomp+1) +
     &                                 delt * u(ibku+no,ipomp+1)
               yc(ibkyc+no,ipomp+1) = y(ibky+no,ipomp+1) +
     &                                 delt * v(ibkv+no,ipomp+1)
               zc(ibkzc+no,ipomp+1) = z(ibkz+no,ipomp+1) +
     &                                 delt * w(ibkw+no,ipomp+1)

            if( jgrav .ne. 0 ) then

               ekin = e(ibke+no,ipomp+1)
               vel  = sqrt( 2.0 * ekin / rtyp ) * rlit
               timd = delt / vel

               vxx = u(ibku+no,ipomp+1) * vel - grvc * grav(1) * timd
               vyy = v(ibkv+no,ipomp+1) * vel - grvc * grav(2) * timd
               vzz = w(ibkw+no,ipomp+1) * vel - grvc * grav(3) * timd

               vab = sqrt( vxx**2 + vyy**2 + vzz**2 )

               ec(ibkec+no,ipomp+1) = 0.5 * rtyp * ( vab / rlit )**2

               xc(ibkxc+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
     &                      - 0.5 * grvc * grav(1) * timd**2
               yc(ibkyc+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
     &                      - 0.5 * grvc * grav(2) * timd**2
               zc(ibkzc+no,ipomp+1) = zc(ibkzc+no,ipomp+1)
     &                      - 0.5 * grvc * grav(3) * timd**2

               u(ibku+no,ipomp+1)   = vxx / vab
               v(ibkv+no,ipomp+1)   = vyy / vab
               w(ibkw+no,ipomp+1)   = vzz / vab

               delt = sqrt( ( xc(ibkxc+no,ipomp+1) -
     &                                          x(ibkx+no,ipomp+1) )**2
     &                    + ( yc(ibkyc+no,ipomp+1) -
     &                                          y(ibky+no,ipomp+1) )**2
     &                    + ( zc(ibkzc+no,ipomp+1) -
     &                                         z(ibkz+no,ipomp+1) )**2 )

            end if

*-----------------------------------------------------------------------

               udir = u(ibku+no,ipomp+1)
               vdir = v(ibkv+no,ipomp+1)
               wdir = w(ibkw+no,ipomp+1)

               dpr = delt

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine magbdx(is6,xx,yy,zz,dmg,cmg,
     &                  bbx,bby,bbz,bba,
     &                  dxx,dyx,dzx,
     &                  dxy,dyy,dzy,
     &                  dxz,dyz,dzz)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------
*     diple only
*-----------------------------------------------------------------------

      if( is6 .eq. 2 ) then

            bbx =  1.d0
            bby =  0.d0
            bbz =  0.d0

            bba = sqrt( bbx**2 + bby**2 + bbz**2 )

            dxx =  0.d0
            dyx =  0.d0
            dzx =  0.d0

            dxy =  0.d0
            dyy =  0.d0
            dzy =  0.d0

            dxz =  0.d0
            dyz =  0.d0
            dzz =  0.d0

*-----------------------------------------------------------------------
*     quadrapole and dipole
*-----------------------------------------------------------------------

      else if( is6 .eq. 4 ) then

            bbx =  xx
            bby = -yy
            bbz = dmg

            bba = sqrt( bbx**2 + bby**2 + bbz**2 )

            dxx =  1.d0
            dyx =  0.d0
            dzx =  0.d0

            dxy =  0.d0
            dyy = -1.d0
            dzy =  0.d0

            dxz =  0.d0
            dyz =  0.d0
            dzz =  0.d0

*-----------------------------------------------------------------------
*     sextapole and dipole
*-----------------------------------------------------------------------

      else if( is6 .eq. 6 ) then

            bbx = ( yy**2 - xx**2 ) / 2.d0
            bby = xx * yy
            bbz = dmg

            bba = sqrt( bbx**2 + bby**2 + bbz**2 )

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
*     from user subroutine 1
*-----------------------------------------------------------------------

      else if( is6 .eq. 1 ) then

            call usrmgf1(xx,yy,zz,dmg,cmg,
     &                   bbx,bby,bbz,
     &                   dxx,dyx,dzx,
     &                   dxy,dyy,dzy,
     &                   dxz,dyz,dzz)

            bba = sqrt( bbx**2 + bby**2 + bbz**2 )

*-----------------------------------------------------------------------
*     from user subroutine 3
*-----------------------------------------------------------------------

      else if( is6 .eq. 3 ) then

            call usrmgf3(xx,yy,zz,dmg,cmg,
     &                   bbx,bby,bbz,
     &                   dxx,dyx,dzx,
     &                   dxy,dyy,dzy,
     &                   dxz,dyz,dzz)

            bba = sqrt( bbx**2 + bby**2 + bbz**2 )

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sexts2(xx,xp,yy,yp,zz,zp,sx,sy,sz,
     &                  thet,delt,dpr,sfc,dmg,cmg,
     &                  grv,gravx,gravy,gravz,is6)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( pi   = 3.1415926535898d0 )

*-----------------------------------------------------------------------

            xx0 = xx
            vx0 = xp
            yy0 = yy
            vy0 = yp
            zz0 = zz
            vz0 = zp

            sx0 = sx
            sy0 = sy
            sz0 = sz

            sxf = sx
            syf = sy
            szf = sz

*-----------------------------------------------------------------------

            nt  = 10
            dt0 = thet / dble( nt )
            ispm = 0
            dt = dt0
            idt = 0

*-----------------------------------------------------------------------

  100    continue

               call magbdx(is6,xx0,yy0,zz0,dmg,cmg,
     &                     bbx0,bby0,bbz0,bba0,
     &                     dxx0,dyx0,dzx0,
     &                     dxy0,dyy0,dzy0,
     &                     dxz0,dyz0,dzz0)

*-----------------------------------------------------------------------

               sn = bba0
               dts = dt0
               if( idt .eq. 1 ) dts = dt

*-----------------------------------------------------------------------

            if( sn * sfc .lt. 1000.0 ) then

               ispm = 1

               dt = dts
               if( sn * sfc .gt. 0.d0 )
     &         dt = min( dts, 0.05 / ( sfc * sn ) )

            else

               ispm = 0
               dt = dts

            end if

*-----------------------------------------------------------------------
*           RKG Second order
*-----------------------------------------------------------------------

         if( ispm .eq. 0 ) then

*-----------------------------------------------------------------------
*           spin rotation
*-----------------------------------------------------------------------

               sba = ( sx0 * bbx0 + sy0 * bby0 + sz0 * bbz0 ) / bba0
               sba = max( -1.d0, sba )
               sba = min(  1.d0, sba )

               costh = bbz0 / bba0
               rt2 = ( bbx0 / bba0 )**2 + ( bby0 / bba0 )**2

            if( rt2 .eq. 0.0d0 ) then

               sinth = 0.0d0
               cosphi= 1.0d0
               sinphi= 0.0d0

            else

               rt = sqrt(rt2)
               sinth  = rt
               cosphi =   ( bbx0 / bba0 ) / rt
               sinphi = - ( bby0 / bba0 ) / rt

            end if

               adc = sx0
               bdc = sy0
               gdc = sz0

               t1  = costh  * adc + sinth  * gdc
               sxb = cosphi * t1  - sinphi * bdc
               syb = sinphi * t1  + cosphi * bdc
               szb = costh  * gdc - sinth  * adc

            if( syb .ge. 0.d0 .and. sxb .eq. 0.d0 ) then
               sbt = pi / 2.d0
            else if( syb .lt. 0.d0 .and. sxb .eq. 0.d0 ) then
               sbt = pi / 2.d0 * 3.d0
            else if( syb .ge. 0.d0 .and. sxb .gt. 0.d0 ) then
               sbt = atan( syb / sxb )
            else if( syb .ge. 0.d0 .and. sxb .lt. 0.d0 ) then
               sbt = atan( syb / sxb ) + pi
            else if( syb .lt. 0.d0 .and. sxb .gt. 0.d0 ) then
               sbt = atan( syb / sxb ) + 2.d0 * pi
            else if( syb .lt. 0.d0 .and. sxb .lt. 0.d0 ) then
               sbt = atan( syb / sxb ) + pi
            end if

*-----------------------------------------------------------------------

               xx1 = xx0 + 0.5 * dt * vx0
               yy1 = yy0 + 0.5 * dt * vy0
               zz1 = zz0 + 0.5 * dt * vz0

               fx0 = - (  bbx0 / bba0 * dxx0
     &                 +  bby0 / bba0 * dyx0
     &                 +  bbz0 / bba0 * dzx0 ) * sba
     &               - grv * gravx

               fy0 = - (  bbx0 / bba0 * dxy0
     &                 +  bby0 / bba0 * dyy0
     &                 +  bbz0 / bba0 * dzy0 ) * sba
     &               - grv * gravy

               fz0 = - (  bbx0 / bba0 * dxz0
     &                 +  bby0 / bba0 * dyz0
     &                 +  bbz0 / bba0 * dzz0 ) * sba
     &               - grv * gravz

               vx1 = vx0 + 0.5 * dt * fx0
               vy1 = vy0 + 0.5 * dt * fy0
               vz1 = vz0 + 0.5 * dt * fz0

               xxf = xx0 + dt * vx1
               yyf = yy0 + dt * vy1
               zzf = zz0 + dt * vz1

               call magbdx(is6,xx1,yy1,zz1,dmg,cmg,
     &                     bbx1,bby1,bbz1,bba1,
     &                     dxx1,dyx1,dzx1,
     &                     dxy1,dyy1,dzy1,
     &                     dxz1,dyz1,dzz1)

               fx1 = - (  bbx1 / bba1 * dxx1
     &                 +  bby1 / bba1 * dyx1
     &                 +  bbz1 / bba1 * dzx1 ) * sba
     &               - grv * gravx

               fy1 = - (  bbx1 / bba1 * dxy1
     &                 +  bby1 / bba1 * dyy1
     &                 +  bbz1 / bba1 * dzy1 ) * sba
     &               - grv * gravy

               fz1 = - (  bbx1 / bba1 * dxz1
     &                 +  bby1 / bba1 * dyz1
     &                 +  bbz1 / bba1 * dzz1 ) * sba
     &               - grv * gravz

               vxf = vx0 + dt * fx1
               vyf = vy0 + dt * fy1
               vzf = vz0 + dt * fz1

               call magbdx(is6,xxf,yyf,zzf,dmg,cmg,
     &                     bbxf,bbyf,bbzf,bbaf,
     &                     dxxf,dyxf,dzxf,
     &                     dxyf,dyyf,dzyf,
     &                     dxzf,dyzf,dzzf)

*-----------------------------------------------------------------------
*           spin rotation
*-----------------------------------------------------------------------

               costh = bbzf / bbaf
               rt2 = ( bbxf / bbaf )**2 + ( bbyf / bbaf )**2

            if( rt2 .eq. 0.0d0 ) then

               sinth = 0.0d0
               cosphi= 1.0d0
               sinphi= 0.0d0

            else

               rt = sqrt(rt2)
               sinth  = rt
               cosphi = ( bbxf / bbaf ) / rt
               sinphi = ( bbyf / bbaf ) / rt

            end if

               cosh  = sba
               sinh  = sqrt( 1.d0 - cosh**2 )

               thet = sbt + bba1 * sfc * dt

               cost = cos( thet )
               sint = sin( thet )

               adc = sinh * cost
               bdc = sinh * sint
               gdc = cosh

               t1  = costh  * adc + sinth  * gdc
               sxf = cosphi * t1  - sinphi * bdc
               syf = sinphi * t1  + cosphi * bdc
               szf = costh  * gdc - sinth  * adc

*-----------------------------------------------------------------------

         else

*-----------------------------------------------------------------------

               xx1 = xx0 + 0.5 * dt * vx0
               yy1 = yy0 + 0.5 * dt * vy0
               zz1 = zz0 + 0.5 * dt * vz0

               fx0 = - (  sx0 * dxx0
     &                 +  sy0 * dyx0
     &                 +  sz0 * dzx0 )
     &               - grv * gravx

               fy0 = - (  sx0 * dxy0
     &                 +  sy0 * dyy0
     &                 +  sz0 * dzy0 )
     &               - grv * gravy

               fz0 = - (  sx0 * dxz0
     &                 +  sy0 * dyz0
     &                 +  sz0 * dzz0 )
     &               - grv * gravz

               vx1 = vx0 + 0.5 * dt * fx0
               vy1 = vy0 + 0.5 * dt * fy0
               vz1 = vz0 + 0.5 * dt * fz0

               gx0 = sfc * ( sy0 * bbz0 - sz0 * bby0 )
               gy0 = sfc * ( sz0 * bbx0 - sx0 * bbz0 )
               gz0 = sfc * ( sx0 * bby0 - sy0 * bbx0 )

               sx1 = sx0 + 0.5 * dt * gx0
               sy1 = sy0 + 0.5 * dt * gy0
               sz1 = sz0 + 0.5 * dt * gz0

               xxf = xx0 + dt * vx1
               yyf = yy0 + dt * vy1
               zzf = zz0 + dt * vz1

               call magbdx(is6,xx1,yy1,zz1,dmg,cmg,
     &                     bbx1,bby1,bbz1,bba1,
     &                     dxx1,dyx1,dzx1,
     &                     dxy1,dyy1,dzy1,
     &                     dxz1,dyz1,dzz1)

               fx1 = - (  sx1 * dxx1
     &                 +  sy1 * dyx1
     &                 +  sz1 * dzx1 )
     &               - grv * gravx

               fy1 = - (  sx1 * dxy1
     &                 +  sy1 * dyy1
     &                 +  sz1 * dzy1 )
     &               - grv * gravy

               fz1 = - (  sx1 * dxz1
     &                 +  sy1 * dyz1
     &                 +  sz1 * dzz1 )
     &               - grv * gravz

               vxf = vx0 + dt * fx1
               vyf = vy0 + dt * fy1
               vzf = vz0 + dt * fz1

               gx1 = sfc * ( sy1 * bbz1 - sz1 * bby1 )
               gy1 = sfc * ( sz1 * bbx1 - sx1 * bbz1 )
               gz1 = sfc * ( sx1 * bby1 - sy1 * bbx1 )

               sxf = sx0 + dt * gx1
               syf = sy0 + dt * gy1
               szf = sz0 + dt * gz1

         end if

*-----------------------------------------------------------------------

            dpr = dsqrt( ( xxf - xx )**2
     &                 + ( yyf - yy )**2
     &                 + ( zzf - zz )**2 )

         if( abs( dpr - delt ) / delt .lt. 0.001 ) goto 200

         if( dpr .gt. delt ) then

            dt = dt / 2.0
            idt = 1

         else if( dpr .lt. delt ) then

            xx0 = xxf
            yy0 = yyf
            zz0 = zzf

            vx0 = vxf
            vy0 = vyf
            vz0 = vzf

            sx0 = sxf
            sy0 = syf
            sz0 = szf

         end if

             ssr = sqrt( sx0**2 + sy0**2 + sz0**2 )

           if( ssr .gt. 1.05 .or. ssr .lt. 0.95 )  then
             write(*,*) 'ERROR: integration is diverged !!', ssr
             write(*,*) '     ispm =', ispm
             write(*,*) '   sn*sfc =', sn * sfc
             write(*,*) '       dt =', dt
             stop 777
           end if

            goto 100

*-----------------------------------------------------------------------

  200    continue

            xx = xxf
            yy = yyf
            zz = zzf

            xp = vxf
            yp = vyf
            zp = vzf

            sx = sxf
            sy = syf
            sz = szf

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sexts1(xx,xp,yy,yp,zz,zp,sx,sy,sz,
     &                  thet,delt,dpr,sfc,
     &                  grv,gravx,gravy,gravz)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( pi   = 3.1415926535898d0 )

*-----------------------------------------------------------------------

            xx0 = xx
            vx0 = xp
            yy0 = yy
            vy0 = yp
            zz0 = zz
            vz0 = zp

            sx0 = sx
            sy0 = sy
            sz0 = sz

            sxf = sx
            syf = sy
            szf = sz

*-----------------------------------------------------------------------

            nt  = 10
            dt0 = thet / dble( nt )
            ispm = 0
            dt = dt0
            idt = 0

*-----------------------------------------------------------------------

  100    continue

               sn = ( xx0**2 + yy0**2 ) / 2.d0
               dts = dt0
               if( idt .eq. 1 ) dts = dt

*-----------------------------------------------------------------------

            if( sn * sfc .lt. 1000.0 ) then

               ispm = 1
               dt = dts
               if( sn * sfc .gt. 0.d0 )
     &         dt = min( dts, 0.05 / ( sfc * sn ) )

            else

               ispm = 0
               dt = dts

            end if

*-----------------------------------------------------------------------
*           RKG Second order
*-----------------------------------------------------------------------

         if( ispm .eq. 0 ) then

*-----------------------------------------------------------------------
*           spin rotation
*-----------------------------------------------------------------------

               sba = ( ( yy0**2 - xx0**2 ) / 2.d0 * sx0
     &                 + xx0 * yy0 * sy0 ) / sn
               sba = max( -1.d0, sba )
               sba = min(  1.d0, sba )

               costh  = 0.d0
               sinth  = 1.d0
               cosphi =   ( yy0**2 - xx0**2 ) / 2.d0 / sn
               sinphi = - xx0 * yy0 / sn

               adc = sx0
               bdc = sy0
               gdc = sz0

               t1  = costh  * adc + sinth  * gdc
               sxb = cosphi * t1  - sinphi * bdc
               syb = sinphi * t1  + cosphi * bdc
               szb = costh  * gdc - sinth  * adc

            if( syb .ge. 0.d0 .and. sxb .eq. 0.d0 ) then
               sbt = pi / 2.d0
            else if( syb .lt. 0.d0 .and. sxb .eq. 0.d0 ) then
               sbt = pi / 2.d0 * 3.d0
            else if( syb .ge. 0.d0 .and. sxb .gt. 0.d0 ) then
               sbt = atan( syb / sxb )
            else if( syb .ge. 0.d0 .and. sxb .lt. 0.d0 ) then
               sbt = atan( syb / sxb ) + pi
            else if( syb .lt. 0.d0 .and. sxb .gt. 0.d0 ) then
               sbt = atan( syb / sxb ) + 2.d0 * pi
            else if( syb .lt. 0.d0 .and. sxb .lt. 0.d0 ) then
               sbt = atan( syb / sxb ) + pi
            end if

*-----------------------------------------------------------------------

               xx1 = xx0 + 0.5 * dt * vx0
               yy1 = yy0 + 0.5 * dt * vy0
               zz1 = zz0 + 0.5 * dt * vz0

               vx1 = vx0 + 0.5 * dt * ( -xx0 * sba - grv * gravx )
               vy1 = vy0 + 0.5 * dt * ( -yy0 * sba - grv * gravy )
               vz1 = vz0 + 0.5 * dt * (            - grv * gravz )

               xxf = xx0 + dt * vx1
               yyf = yy0 + dt * vy1
               zzf = zz0 + dt * vz1

               vxf = vx0 + dt * ( -xx1 * sba - grv * gravx )
               vyf = vy0 + dt * ( -yy1 * sba - grv * gravy )
               vzf = vz0 + dt * (            - grv * gravz )

*-----------------------------------------------------------------------
*           spin rotation
*-----------------------------------------------------------------------

               costh  = 0.d0
               sinth  = 1.d0
               cosphi = ( yyf**2 - xxf**2 )
     &                / ( xxf**2 + yyf**2 )
               sinphi = 2.0 * xxf * yyf
     &                / ( xxf**2 + yyf**2 )

               cosh  = sba
               sinh  = sqrt( 1.d0 - cosh**2 )

               sn1  = ( xx0**1 + yy1**2 ) / 2.d0
               thet = sbt + sn1 * sfc * dt

               cost = cos( thet )
               sint = sin( thet )

               adc = sinh * cost
               bdc = sinh * sint
               gdc = cosh

               t1  = costh  * adc + sinth  * gdc
               sxf = cosphi * t1  - sinphi * bdc
               syf = sinphi * t1  + cosphi * bdc
               szf = costh  * gdc - sinth  * adc

*-----------------------------------------------------------------------

         else

*-----------------------------------------------------------------------

               xx1 = xx0 + 0.5 * dt * vx0
               yy1 = yy0 + 0.5 * dt * vy0
               zz1 = zz0 + 0.5 * dt * vz0

               vx1 = vx0 + 0.5 * dt
     &             * (  xx0 * sx0 - yy0 * sy0 - grv * gravx )
               vy1 = vy0 + 0.5 * dt
     &             * ( -yy0 * sx0 - xx0 * sy0 - grv * gravy )
               vz1 = vz0 + 0.5 * dt
     &             * (                        - grv * gravz )

               sx1 = sx0 + 0.5 * dt * sfc
     &             * ( - xx0 * yy0 * sz0 )
               sy1 = sy0 + 0.5 * dt * sfc
     &             * ( yy0**2 - xx0**2 ) / 2.d0 * sz0
               sz1 = sz0 + 0.5 * dt * sfc
     &             * ( xx0 * yy0 * sx0
     &             - ( yy0**2 - xx0**2 ) / 2.d0 * sy0 )

               xxf = xx0 + dt * vx1
               yyf = yy0 + dt * vy1
               zzf = zz0 + dt * vz1

               vxf = vx0 + dt
     &             * (  xx1 * sx1 - yy1 * sy1 - grv * gravx )
               vyf = vy0 + dt
     &             * ( -yy1 * sx1 - xx1 * sy1 - grv * gravy )
               vzf = vz0 + dt
     &             * (                        - grv * gravz )

               sxf = sx0 + dt * sfc
     &             * ( - xx1 * yy1 * sz1 )
               syf = sy0 + dt * sfc
     &             * ( yy1**2 - xx1**2 ) / 2.d0 * sz1
               szf = sz0 + dt * sfc
     &             * ( xx1 * yy1 * sx1
     &             - ( yy1**2 - xx1**2 ) / 2.d0 * sy1 )

         end if

*-----------------------------------------------------------------------

            dpr = dsqrt( ( xxf - xx )**2
     &                 + ( yyf - yy )**2
     &                 + ( zzf - zz )**2 )

         if( abs( dpr - delt ) / delt .lt. 0.001 ) goto 200

         if( dpr .gt. delt ) then

            dt = dt / 2.0
            idt = 1

         else if( dpr .lt. delt ) then

            xx0 = xxf
            yy0 = yyf
            vx0 = vxf
            vy0 = vyf
            zz0 = zzf

            sx0 = sxf
            sy0 = syf
            sz0 = szf

         end if

             ssr = sqrt( sx0**2 + sy0**2 + sz0**2 )

           if( ssr .gt. 1.05 .or. ssr .lt. 0.95 )  then
             write(*,*) 'ERROR: integration is diverged !!', ssr
             write(*,*) '     ispm =', ispm
             write(*,*) '   sn*sfc =', sn * sfc
             write(*,*) '       dt =', dt
             stop 777
           end if

            goto 100

*-----------------------------------------------------------------------

  200    continue

            xx = xxf
            yy = yyf
            zz = zzf
            xp = vxf
            yp = vyf

            sx = sxf
            sy = syf
            sz = szf

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine sexts0(dz,xx,xp,yy,yp,sx,sy,sz,thet)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

         sba = ( yy**2 - xx**2 ) * sx
     &       + 2.0 * xx * yy     * sy
         ipal = 1
         if( sba .lt. 0.0d0 ) ipal = -1
         if( sba .eq. 0.0d0 .and.
     &       unirn(dummy) .le. 0.5d0 ) ipal = -1

      if( ipal .gt. 0.0 ) then

         cs = cos( thet )
         sn = sin( thet )

         xn = xx * cs + xp * sn
         yn = yy * cs + yp * sn

         vx = -xx * sn + xp * cs
         vy = -yy * sn + yp * cs

      else

         cs = cosh( thet )
         sn = sinh( thet )

         xn = xx * cs + xp * sn
         yn = yy * cs + yp * sn

         vx = xx * sn + xp * cs
         vy = yy * sn + yp * cs

      end if

         xx = xn
         yy = yn
         xp = vx
         yp = vy

         sn = xx**2 + yy**2
         sx = ipal * ( yy**2 - xx**2 ) / sn
         sy = ipal * 2.d0 * xx * yy    / sn
         sz = 0.0d0

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine quad(p,b,a0,xl0,x,xp,y,yp,cg)
*                                                                      *
*       quadropole magnet                                              *
*       modified by K.Niita on 2003/08/19                              *
*                                                                      *
*       p    : momentum [GeV/c]
*       b    : field [kG]
*       a    : gap   [cm]
*       xl   : length [cm]
*       x,y  : position [cm]
*       x',y': divergence [rad]
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension coeff1(6,6)

*-----------------------------------------------------------------------

         a  = a0  / 100.d0
         al = xl0 / 100.d0
         xp = xp  * 1000.d0
         yp = yp  * 1000.d0

*-----------------------------------------------------------------------
*       a    : gap   [m]
*       xl   : length [m]
*       x,y  : position [cm]
*       x',y': divergence [mrad]
*-----------------------------------------------------------------------

      if( cg .gt. 0.0 ) then
         bb =  b
      else
         bb = -b
      end if

         akq = dsqrt( abs( b / ( a * p / cg * 33.356 ) ) )

*-----------------------------------------------------------------------
*     define terms for horizontal and vertical focus
*-----------------------------------------------------------------------

      if( bb .ge. 0.0 ) then

         c1 = cos(akq*al)
         s1 = sin(akq*al)
         c2 = cosh(akq*al)
         s2 = sinh(akq*al)

         sign = 1.0

      else

         c1 = cosh(akq*al)
         s1 = sinh(akq*al)
         c2 = cos(akq*al)
         s2 = sin(akq*al)

         sign = -1.0

      end if

*-----------------------------------------------------------------------
*     1st order terms
*-----------------------------------------------------------------------

      coeff1(1,1) = c1
      coeff1(2,2) = c1
      coeff1(3,3) = c2
      coeff1(4,4) = c2

      coeff1(1,2) = s1 *0.1 / akq
      coeff1(2,1) = -sign * akq * s1 * 10.0
      coeff1(3,4) = s2 * 0.1 / akq
      coeff1(4,3) = sign * akq * s2 * 10.0

      xn  = coeff1(1,1) * x + coeff1(1,2) * xp
      xpn = coeff1(2,1) * x + coeff1(2,2) * xp
      yn  = coeff1(3,3) * y + coeff1(3,4) * yp
      ypn = coeff1(4,3) * y + coeff1(4,4) * yp

*-----------------------------------------------------------------------

      x  = xn
      y  = yn
      xp = xpn / 1000.d0
      yp = ypn / 1000.d0

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine magmap(xx,yy,zz,vx,vy,vz,tt,qm,
     &                  ssx,ssy,ssz,alph,gamn,
     &                  ineu,imap,b_mag,
     &                  x1,x2,y1,y2,z1,z2,r1,r2,
     &                  bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                  bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                  bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                  bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                  bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                  bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                  bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                  bxx2y2z2, byx2y2z2, bzx2y2z2,
     &                  brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                  brr2z1, bzr2z1, brr2z2, bzr2z2,mark)
*                                                                      *
*       magnet                                                         *
*       modified by ASTOM                                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)
      parameter ( vc = 0.29979239d0 )
      parameter ( xc = 0.14989619d0 )
      parameter ( cvel = 2.997925d+8 )
      parameter ( cneu = 0.898755d+11 )
      include 'err.inc'

*-----------------------------------------------------------------------

      xxc = xx
      yyc = yy
      zzc = zz
      ssxc = ssx
      ssyc = ssy
      sszc = ssz
      ssrc = sqrt(ssxc**2 + ssyc**2 + sszc**2)
      if(ssrc.ne.0) then
       ssxc = ssxc /ssrc
       ssyc = ssyc /ssrc
       sszc = sszc /ssrc
      endif
      dssx = 0
      dssy = 0
      dssz = 0

*     read cubic 8 paramaters then interpolate
*     xx:particle location, x1:lower end of the cube, x2:upper
*     y:, z:, as same
*     subroutin readmap1 returns bi-linear interprated magnetic
*     field and x-y-z differencials of it
      if (imap .eq. 1) then
            call readmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
      else if (imap .eq. 2) then
            call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
      else if (imap .eq. 3) then
            call readmap3(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
      else if (imap .eq. 4) then
            call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
      else
      end if

*    Runge Kutta
      ckx0 = tt * vx
      cky0 = tt * vy
      ckz0 = tt * vz

      if (ineu .eq. 0) then
            ax = qm * ( vy * bfbz - vz * bfby )
            ay = qm * ( vz * bfbx - vx * bfbz )
            az = qm * ( vx * bfby - vy * bfbx )
      else
            ax = -1 * alph * (ssxc*dxbfbx+ssyc*dxbfby+sszc*dxbfbz)
     &            * 100 / (cneu * cvel)
            ay = -1 * alph * (ssxc*dybfbx+ssyc*dybfby+sszc*dybfbz)
     &            * 100 / (cneu * cvel)
            az = -1 * alph * (ssxc*dzbfbx+ssyc*dzbfby+sszc*dzbfbz)
     &            * 100 / (cneu * cvel)

            dssx = -1 * gamn * ( ssyc*bfbz - sszc*bfby )
            dssy = -1 * gamn * ( sszc*bfbx - ssxc*bfbz )
            dssz = -1 * gamn * ( ssxc*bfby - ssyc*bfbx )

      end if


      clx0 = tt * ax
      cly0 = tt * ay
      clz0 = tt * az

      cmx0 = (0.001/cvel) * tt * dssx
      cmy0 = (0.001/cvel) * tt * dssy
      cmz0 = (0.001/cvel) * tt * dssz

*-----------------------------------------------------------------------
      xxc = xx+ckx0/2
      yyc = yy+cky0/2
      zzc = zz+ckz0/2

      ssxc = ssx+cmx0/2
      ssyc = ssy+cmy0/2
      sszc = ssz+cmz0/2
      ssrc = sqrt(ssxc**2 + ssyc**2 + sszc**2)
      if(ssrc.ne.0) then
       ssxc = ssxc /ssrc
       ssyc = ssyc /ssrc
       sszc = sszc /ssrc
      endif

      if (imap .eq. 1) then
            call readmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
      else if (imap .eq. 2) then
            call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
      else if (imap .eq. 3) then
            call readmap3(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
      else if (imap .eq. 4) then
            call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
      else
      end if

*    Runge Kutta
      ckx1 = tt * (vx+xc*clx0/2)
      cky1 = tt * (vy+xc*cly0/2)
      ckz1 = tt * (vz+xc*clz0/2)

      if (ineu .eq. 0) then
            ax = qm * ((vy+cly0/2) * bfbz - (vz+clz0/2) * bfby)
            ay = qm * ((vz+clz0/2) * bfbx - (vx+clx0/2) * bfbz)
            az = qm * ((vx+clx0/2) * bfby - (vy+cly0/2) * bfbx)
      else
            ax = -1 * alph * (ssxc*dxbfbx+ssyc*dxbfby+sszc*dxbfbz)
     &            * 100 / (cneu * cvel)
            ay = -1 * alph * (ssxc*dybfbx+ssyc*dybfby+sszc*dybfbz)
     &            * 100 / (cneu * cvel)
            az = -1 * alph * (ssxc*dzbfbx+ssyc*dzbfby+sszc*dzbfbz)
     &            * 100 / (cneu * cvel)

            dssx = -1 * gamn * (ssyc*bfbz - ssyc*bfby)
            dssy = -1 * gamn * (sszc*bfbx - ssxc*bfbz)
            dssz = -1 * gamn * (ssxc*bfby - ssyc*bfbx)

      end if


      clx1 = tt * ax
      cly1 = tt * ay
      clz1 = tt * az

      cmx1 = (0.001/cvel) * tt * dssx
      cmy1 = (0.001/cvel) * tt * dssy
      cmz1 = (0.001/cvel) * tt * dssz

*-----------------------------------------------------------------------

      xxc = xx+ckx1/2
      yyc = yy+cky1/2
      zzc = zz+ckz1/2

      ssxc = ssx+cmx1/2
      ssyc = ssy+cmy1/2
      sszc = ssz+cmz1/2
      ssrc = sqrt(ssxc**2 + ssyc**2 + sszc**2)
      if(ssrc.ne.0) then
       ssxc = ssxc /ssrc
       ssyc = ssyc /ssrc
       sszc = sszc /ssrc
      endif

      if (imap .eq. 1) then
            call readmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
      else if (imap .eq. 2) then
            call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
      else if (imap .eq. 3) then
            call readmap3(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
      else if (imap .eq. 4) then
            call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
      else
      end if

*    Runge Kutta
      ckx2 = tt * (vx+xc*clx1/2)
      cky2 = tt * (vy+xc*cly1/2)
      ckz2 = tt * (vz+xc*clz1/2)

      if(ineu .eq. 0) then
            ax = qm * ( (vy+cly1/2) * bfbz - (vz+clz1/2) * bfby )
            ay = qm * ( (vz+clz1/2) * bfbx - (vx+clx1/2) * bfbz )
            az = qm * ( (vx+clx1/2) * bfby - (vy+cly1/2) * bfbx )
      else
            ax = -1 * alph * (ssxc*dxbfbx+ssyc*dxbfby+sszc*dxbfbz)
     &            * 100 / (cneu * cvel)
            ay = -1 * alph * (ssxc*dybfbx+ssyc*dybfby+sszc*dybfbz)
     &            * 100 / (cneu * cvel)
            az = -1 * alph * (ssxc*dzbfbx+ssyc*dzbfby+sszc*dzbfbz)
     &            * 100 / (cneu * cvel)

            dssx = -1 * gamn * ( ssyc*bfbz - ssyc*bfby )
            dssy = -1 * gamn * ( sszc*bfbx - ssxc*bfbz )
            dssz = -1 * gamn * ( ssxc*bfby - ssyc*bfbx )
      end if


      clx2 = tt * ax
      cly2 = tt * ay
      clz2 = tt * az

      cmx2 = (0.001/cvel) * tt * dssx
      cmy2 = (0.001/cvel) * tt * dssy
      cmz2 = (0.001/cvel) * tt * dssz

*-----------------------------------------------------------------------

      xxc = xx+ckx2
      yyc = yy+cky2
      zzc = zz+ckz2

      ssxc = ssx+cmx2
      ssyc = ssy+cmy2
      sszc = ssz+cmz2
      ssrc = sqrt(ssxc**2 + ssyc**2 + sszc**2)
      if(ssrc.ne.0) then
       ssxc = ssxc /ssrc
       ssyc = ssyc /ssrc
       sszc = sszc /ssrc
      endif

      if (imap .eq. 1) then
            call readmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
      else if (imap .eq. 2) then
            call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
      else if (imap .eq. 3) then
            call readmap3(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
      else if (imap .eq. 4) then
            call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
      else
      end if

*    Runge Kutta
      ckx3 = tt * (vx+xc*clx2)
      cky3 = tt * (vy+xc*cly2)
      ckz3 = tt * (vz+xc*clz2)

      if(ineu .eq. 0) then
            ax = qm * ( (vy+cly2) * bfbz - (vz+clz2) * bfby )
            ay = qm * ( (vz+clz2) * bfbx - (vx+clx2) * bfbz )
            az = qm * ( (vx+clx2) * bfby - (vy+cly2) * bfbx )
      else
            ax = -1 * alph * (ssxc*dxbfbx+ssyc*dxbfby+sszc*dxbfbz)
     &            * 100 / (cneu * cvel)
            ay = -1 * alph * (ssxc*dybfbx+ssyc*dybfby+sszc*dybfbz)
     &            * 100 / (cneu * cvel)
            az = -1 * alph * (ssxc*dzbfbx+ssyc*dzbfby+sszc*dzbfbz)
     &            * 100 / (cneu * cvel)

            dssx = -1 * gamn * ( ssyc*bfbz - ssyc*bfby )
            dssy = -1 * gamn * ( sszc*bfbx - ssxc*bfbz )
            dssz = -1 * gamn * ( ssxc*bfby - ssyc*bfbx )
      end if


      clx3 = tt * ax
      cly3 = tt * ay
      clz3 = tt * az

      cmx3 = (0.001/cvel) * tt * dssx
      cmy3 = (0.001/cvel) * tt * dssy
      cmz3 = (0.001/cvel) * tt * dssz

*-----------------------------------------------------------------------

      ckx = (ckx0 + 2*ckx1 + 2*ckx2 + ckx3) / 6
      cky = (cky0 + 2*cky1 + 2*cky2 + cky3) / 6
      ckz = (ckz0 + 2*ckz1 + 2*ckz2 + ckz3) / 6

      clx = (clx0 + 2*clx1 + 2*clx2 + clx3) / 6
      cly = (cly0 + 2*cly1 + 2*cly2 + cly3) / 6
      clz = (clz0 + 2*clz1 + 2*clz2 + clz3) / 6

      cmx = (cmx0 + 2*cmx1 + 2*cmx2 + cmx3) / 6
      cmy = (cmy0 + 2*cmy1 + 2*cmy2 + cmy3) / 6
      cmz = (cmz0 + 2*cmz1 + 2*cmz2 + cmz3) / 6

      xn = xx + ckx
      yn = yy + cky
      zn = zz + ckz
      vxn = vx + vc * clx
      vyn = vy + vc * cly
      vzn = vz + vc * clz

      ssxn = ssx + cmx
      ssyn = ssy + cmy
      sszn = ssz + cmz


*-----------------------------------------------------------------------

      vx = vxn
      vy = vyn
      vz = vzn
      xx = xn
      yy = yn
      zz = zn
      ssx = ssxn
      ssy = ssyn
      ssz = sszn

      if(vx.ge.1.0e30.or.vy.ge.1.0e30.or.vz.ge.1.0e30) then
       write(ErrCha,*) 'magnetic field is too strong',
     & ' to be handled by PHITS. This problem may be solved',
     & ' by increasing gap parameter.'
       ErrID = 'L:2660/R:magmap/F:magtrs.f' !E10_010_001
       call ErrWrite(ErrID,ErrCha)
       mark=-2 ! T.Sato 2022/03/22
      endif

*-----------------------------------------------------------------------

      return
      end

*****************************************************************
      subroutine readmagxyzlist(map1,icmap1) ! Read xyz magnetic field list
      use magmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /bxyzlistint/nx,ny,nz,ibxyzflip(3,3),ixyztype(3)

      character chin*200,cmmt*2 ! T.Sato 2019/01/13
      character map1*200

      cmmt="#$"  ! only two characters are used as comment

      open(790,file=map1(1:icmap1),status='old')
      iline=0
! initialization
      nx=1
      ny=1
      nz=1
      ibxyzflip(:,:)=0 ! flip index

      do iheader=1,1000
       iline=iline+1
       read(790, '(a)', err=989) chin
       call chlngt(chin,200,i1,i2)  ! calculate character length
       call chcaps(chin,i1,i2,i3,cmmt) ! CAPITAL to lower character
       i1=1                         ! delete space even for 1st character
       call chcomp(chin,i1,i2,i4)   ! delete space
       if(chin(1:4).eq.'data') then
        exit
       elseif(chin(1:3).eq.'nx=') then
        read(chin(4:200),*,err=990) nx
       elseif(chin(1:3).eq.'ny=') then
        read(chin(4:200),*,err=990) ny
       elseif(chin(1:3).eq.'nz=') then
        read(chin(4:200),*,err=990) nz
       elseif(chin(1:7).eq.'extendx') then
        ibxyzflip(1,:)=1
        if(chin(8:12).eq.'flip=') then
         if(chin(13:14).eq.'bx'.or.chin(16:17).eq.'bx'
     &   .or.chin(19:20).eq.'bx') ibxyzflip(1,1)=-1
         if(chin(13:14).eq.'by'.or.chin(16:17).eq.'by'
     &   .or.chin(19:20).eq.'by') ibxyzflip(1,2)=-1
         if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz'
     &   .or.chin(19:20).eq.'bz') ibxyzflip(1,3)=-1
        endif
       elseif(chin(1:7).eq.'extendy') then
        ibxyzflip(2,:)=1
        if(chin(8:12).eq.'flip=') then
         if(chin(13:14).eq.'bx'.or.chin(16:17).eq.'bx'
     &   .or.chin(19:20).eq.'bx') ibxyzflip(2,1)=-1
         if(chin(13:14).eq.'by'.or.chin(16:17).eq.'by'
     &   .or.chin(19:20).eq.'by') ibxyzflip(2,2)=-1
         if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz'
     &   .or.chin(19:20).eq.'bz') ibxyzflip(2,3)=-1
        endif
       elseif(chin(1:7).eq.'extendz') then
        ibxyzflip(3,:)=1
        if(chin(8:12).eq.'flip=') then
         if(chin(13:14).eq.'bx'.or.chin(16:17).eq.'bx'
     &   .or.chin(19:20).eq.'bx') ibxyzflip(3,1)=-1
         if(chin(13:14).eq.'by'.or.chin(16:17).eq.'by'
     &   .or.chin(19:20).eq.'by') ibxyzflip(3,2)=-1
         if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz'
     &   .or.chin(19:20).eq.'bz') ibxyzflip(3,3)=-1
        endif
       elseif(chin(1:1).ne.' ') then ! not comment line
        write(ErrCha,*) 'Skip unknown parameter during reading',
     &  ' magnetic field file at line =',iheader
        ErrID = 'L:2740/R:readmagxyzlist/F:magtrs.f' !W10_001_001
        call ErrWrite(ErrID,ErrCha)
       endif
      enddo
!     read header end

      if(nx.le.1) goto 991
      if(ny.le.1) goto 992
      if(nz.le.1) goto 993

      call ALLOCATE_bxyzlist(nx,ny,nz)

      do ix=1,nx
       do iy=1,ny
        do iz=1,nZ
         iline=iline+1
         read(790, *, err=996) xtmp,ytmp,ztmp,
     &   bxlist(ix,iy,iz),bylist(ix,iy,iz),bzlist(ix,iy,iz)
         if(iy.eq.1.and.iz.eq.1) then
          bxmesh(ix)=xtmp
          if(ix.ne.1) then
           if(bxmesh(ix).le.bxmesh(ix-1)) goto 1010
          endif
         elseif(bxmesh(ix).ne.xtmp) then
          goto 1000
         endif
         if(ix.eq.1.and.iz.eq.1) then
          bymesh(iy)=ytmp
          if(iy.ne.1) then
           if(bymesh(iy).le.bymesh(iy-1)) goto 1011
          endif
         elseif(bymesh(iy).ne.ytmp) then
          goto 1001
         endif
         if(ix.eq.1.and.iy.eq.1) then
          bzmesh(iz)=ztmp
          if(iz.ne.1) then
           if(bzmesh(iz).le.bzmesh(iz-1)) goto 1012
          endif
         elseif(bzmesh(iz).ne.ztmp) then
          goto 1002
         endif
        enddo
       enddo
      enddo

! check free or fixed mesh
      dX=bxmesh(2)-bxmesh(1)
      do ix=2,nx-1
       if(bxmesh(ix+1)-bxmesh(ix).ne.dX) exit ! different mesh
      enddo
      if(ix.eq.nx) then
       ixyztype(1)=2  ! fixed x-mesh
      else
       ixyztype(1)=1  ! free x-mesh
      endif

      dY=bymesh(2)-bymesh(1)
      do iy=2,ny-1
       if(bymesh(iy+1)-bymesh(iy).ne.dy) exit ! different mesh
      enddo
      if(iy.eq.ny) then
       ixyztype(2)=2  ! fixed y-mesh
      else
       ixyztype(2)=1  ! free y-mesh
      endif

      dZ=bzmesh(2)-bzmesh(1)
      do iz=2,nZ-1
       if(bzmesh(iz+1)-bzmesh(iz).ne.dZ) exit ! different mesh
      enddo
      if(iz.eq.nZ) then
       ixyztype(3)=2  ! fixed z-mesh
      else
       ixyztype(3)=1  ! free z-mesh
      endif

      close(790)

      return

  989  write(ErrCha,*) 'End of file during reading ',map1(1:icmap1)
       ErrID = 'L:2822/R:readmagxyzlist/F:magtrs.f' !E10_002_001
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  990  write(ErrCha,*) 'Error in reading header of ',map1(1:icmap1),
     & ' at line =',iheader
       ErrID = 'L:2828/R:readmagxyzlist/F:magtrs.f' !E10_003_001
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  991  write(ErrCha,*) 'nX should be greater than 1 in ',map1(1:icmap1),
     & ' ; nX=',nx
       ErrID = 'L:2834/R:readmagxyzlist/F:magtrs.f' !E10_004_001
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  992  write(ErrCha,*) 'nY should be greater than 1 in ',map1(1:icmap1),
     & ' ; nY=',ny
       ErrID = 'L:2840/R:readmagxyzlist/F:magtrs.f' !E10_004_002
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  993  write(ErrCha,*) 'nZ should be greater than 1 in ',map1(1:icmap1),
     & ' ; nZ=',nz
       ErrID = 'L:2846/R:readmagxyzlist/F:magtrs.f' !E10_004_003
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  996  write(ErrCha,*) 'Error in reading magnetic fields from ',
     & map1(1:icmap1),' at ix,iy,iz =',ix,iy,iz
       ErrID = 'L:2852/R:readmagxyzlist/F:magtrs.f' !E10_005_001
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1000  write(ErrCha,*) 'Inconsistent x-mesh in ',
     & map1(1:icmap1),' at ix,iy,iz =',ix,iy,iz
       ErrID = 'L:2858/R:readmagxyzlist/F:magtrs.f' !E10_006_001
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1001  write(ErrCha,*) 'Inconsistent y-mesh in ',
     & map1(1:icmap1),' at ix,iy,iz =',ix,iy,iz
       ErrID = 'L:2864/R:readmagxyzlist/F:magtrs.f' !E10_006_002
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1002  write(ErrCha,*) 'Inconsistent z-mesh in ',
     & map1(1:icmap1),' at ix,iy,iz =',ix,iy,iz
       ErrID = 'L:2870/R:readmagxyzlist/F:magtrs.f' !E10_006_003
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1010  write(ErrCha,*) 'X-mesh should be in accending order in ',
     & map1(1:icmap1),' at x(ix) & x(ix-1)',bxmesh(ix),bxmesh(ix-1)
       ErrID = 'L:2876/R:readmagxyzlist/F:magtrs.f' !E10_007_001
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1011  write(ErrCha,*) 'Y-mesh should be in accending order in ',
     & map1(1:icmap1),' at y(iy) & y(iy-1)',bymesh(iy),bymesh(iy-1)
       ErrID = 'L:2882/R:readmagxyzlist/F:magtrs.f' !E10_007_002
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1012  write(ErrCha,*) 'Z-mesh should be in accending order in ',
     & map1(1:icmap1),' at z(iz) & z(iz-1)',bzmesh(iz),bzmesh(iz-1)
       ErrID = 'L:2888/R:readmagxyzlist/F:magtrs.f' !E10_007_003
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1900  write(ErrCha,'('' Please check line '',i8)') iline
       ErrID = 'L:2893/R:readmagxyzlist/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)

 1999  close(790)
       call parastop( 120 )

      end

************************************************************************
*                                                                      *
      subroutine readmap1(xx,yy,zz,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
*                                                                      *
************************************************************************

      use magmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /bxyzlistint/nx,ny,nz,ibxyzflip(3,3),ixyztype(3)

*-----------------------------------------------------------------------
      smul = b_mag

      x0=bxmesh(1)
      x99=bxmesh(nx)
      y0=bymesh(1)
      y99=bymesh(ny)
      z0=bzmesh(1)
      z99=bzmesh(nz)

      if(ixyztype(1).ne.1) dX=(x99-x0)/(nx-1)  ! fixed mesh
      if(ixyztype(2).ne.1) dY=(y99-y0)/(ny-1)  ! fixed mesh
      if(ixyztype(3).ne.1) dZ=(z99-z0)/(nz-1)  ! fixed mesh

      if(ibxyzflip(1,1).eq.0.or.xx.ge.0.0) then ! no extended field, or particle in positive x
       xc = min(x99,max(x0,xx))
       xxr=smul
       yyr=smul
       zzr=smul
      else ! extended field
       xc = min(x99,max(x0,dabs(xx)))
       xxr=smul*ibxyzflip(1,1)
       yyr=smul*ibxyzflip(1,2)
       zzr=smul*ibxyzflip(1,3)
      endif

      if(ibxyzflip(2,1).eq.0.or.yy.ge.0.0) then ! no extended field, or particle in positive y
       yc = min(y99,max(y0,yy))
      else ! extended field, flip again if necessary
       yc = min(y99,max(y0,dabs(yy)))
       xxr=xxr*ibxyzflip(2,1)
       yyr=yyr*ibxyzflip(2,2)
       zzr=zzr*ibxyzflip(2,3)
      endif

      if(ibxyzflip(3,1).eq.0.or.zz.ge.0.0) then ! no extended field, or particle in positive z
       zc = min(z99,max(z0,zz))
      else ! extended field, flip again if necessary
       zc = min(z99,max(z0,dabs(zz)))
       xxr=xxr*ibxyzflip(3,1)
       yyr=yyr*ibxyzflip(3,2)
       zzr=zzr*ibxyzflip(3,3)
      endif

      if(ixyztype(1).eq.1) then ! free mesh, T.Sato 2020/12/20 bug fix
       do nx2=2,nx-1
        if(xc.lt.bxmesh(nx2)) exit
       enddo
       nx1=nx2-1
       x1=bxmesh(nx1)
       x2=bxmesh(nx2)
      else
       nX1 = min(nx-1,int((xc-x0)/dX) + 1)
       nX2 = nX1 + 1
       x1 = x0 + (nx1-1) * dx
       x2 = x1 + dx
      endif

      if(ixyztype(2).eq.1) then ! free mesh, T.Sato 2020/12/20 bug fix
       do ny2=2,ny-1
        if(yc.lt.bymesh(ny2)) exit
       enddo
       ny1=ny2-1
       y1=bymesh(ny1)
       y2=bymesh(ny2)
      else
       nY1 = min(ny-1,int((yc-Y0)/dY) + 1)
       nY2 = nY1 + 1
       y1 = y0 + (ny1-1) * dy
       y2 = y1 + dy
      endif

      if(ixyztype(3).eq.1) then ! free mesh, T.Sato 2020/12/20 bug fix
       do nz2=2,nz-1
        if(zc.lt.bzmesh(nz2)) exit
       enddo
       nz1=nz2-1
       z1=bzmesh(nz1)
       z2=bzmesh(nz2)
      else
       nZ1 = min(nz-1,int((zc-Z0)/dZ) + 1)
       nZ2 = nZ1 + 1
       z1 = z0 + (nz1-1) * dz
       z2 = z1 + dz
      endif

*-----------------------------------------------------------------------

      bxx1y1z1=bxlist(nx1,ny1,nz1)*xxr
      byx1y1z1=bylist(nx1,ny1,nz1)*yyr
      bzx1y1z1=bzlist(nx1,ny1,nz1)*zzr

      bxx1y1z2=bxlist(nx1,ny1,nz2)*xxr
      byx1y1z2=bylist(nx1,ny1,nz2)*yyr
      bzx1y1z2=bzlist(nx1,ny1,nz2)*zzr

      bxx1y2z1=bxlist(nx1,ny2,nz1)*xxr
      byx1y2z1=bylist(nx1,ny2,nz1)*yyr
      bzx1y2z1=bzlist(nx1,ny2,nz1)*zzr

      bxx1y2z2=bxlist(nx1,ny2,nz2)*xxr
      byx1y2z2=bylist(nx1,ny2,nz2)*yyr
      bzx1y2z2=bzlist(nx1,ny2,nz2)*zzr

      bxx2y1z1=bxlist(nx2,ny1,nz1)*xxr
      byx2y1z1=bylist(nx2,ny1,nz1)*yyr
      bzx2y1z1=bzlist(nx2,ny1,nz1)*zzr

      bxx2y1z2=bxlist(nx2,ny1,nz2)*xxr
      byx2y1z2=bylist(nx2,ny1,nz2)*yyr
      bzx2y1z2=bzlist(nx2,ny1,nz2)*zzr

      bxx2y2z1=bxlist(nx2,ny2,nz1)*xxr
      byx2y2z1=bylist(nx2,ny2,nz1)*yyr
      bzx2y2z1=bzlist(nx2,ny2,nz1)*zzr

      bxx2y2z2=bxlist(nx2,ny2,nz2)*xxr
      byx2y2z2=bylist(nx2,ny2,nz2)*yyr
      bzx2y2z2=bzlist(nx2,ny2,nz2)*zzr

      px = (xc-x1) / (x2-x1)
      py = (yc-y1) / (y2-y1)
      pz = (zc-z1) / (z2-z1)

*    interpolation
      bfbx = (1-py)*(1-px)*(1-pz)*bxx1y1z1
     &     + (1-py)*px*(1-pz)*bxx2y1z1
     &     + py*(1-px)*(1-pz)*bxx1y2z1
     &     + py*px*(1-pz)*bxx2y2z1
     &     + (1-py)*(1-px)*pz*bxx1y1z2
     &     + (1-py)*px*pz*bxx2y1z2
     &     + py*(1-px)*pz*bxx1y2z2
     &     + py*px*pz*bxx2y2z2

      bfby = (1-py)*(1-px)*(1-pz)*byx1y1z1
     &     + (1-py)*px*(1-pz)*byx2y1z1
     &     + py*(1-px)*(1-pz)*byx1y2z1
     &     + py*px*(1-pz)*byx2y2z1
     &     + (1-py)*(1-px)*pz*byx1y1z2
     &     + (1-py)*px*pz*byx2y1z2
     &     + py*(1-px)*pz*byx1y2z2
     &     + py*px*pz*byx2y2z2

      bfbz = (1-py)*(1-px)*(1-pz)*bzx1y1z1
     &     + (1-py)*px*(1-pz)*bzx2y1z1
     &     + py*(1-px)*(1-pz)*bzx1y2z1
     &     + py*px*(1-pz)*bzx2y2z1
     &     + (1-py)*(1-px)*pz*bzx1y1z2
     &     + (1-py)*px*pz*bzx2y1z2
     &     + py*(1-px)*pz*bzx1y2z2
     &     + py*px*pz*bzx2y2z2

*     d(bfbx)/dx
      dxbfbx = ((1-py)*(1-pz)*bxx1y1z1 + (1-py)*pz*bxx1y1z2
     &       + py*(1-pz)*bxx1y2z1 + py*pz*bxx1y2z2
     &       - (1-py)*(1-pz)*bxx2y1z1 - (1-py)*pz*bxx2y1z2
     &       - py*(1-pz)*bxx2y2z1 - py*pz*bxx2y2z2) / (x1-x2)
*     d(bfby)/dx
      dxbfby = ((1-py)*(1-pz)*byx1y1z1 + (1-py)*pz*byx1y1z2
     &       + py*(1-pz)*byx1y2z1 + py*pz*byx1y2z2
     &       - (1-py)*(1-pz)*byx2y1z1 - (1-py)*pz*byx2y1z2
     &       - py*(1-pz)*byx2y2z1 - py*pz*byx2y2z2) / (x1-x2)
*     d(bfbz)/dx
      dxbfbz = ((1-py)*(1-pz)*bzx1y1z1 + (1-py)*pz*bzx1y1z2
     &       + py*(1-pz)*bzx1y2z1 + py*pz*bzx1y2z2
     &       - (1-py)*(1-pz)*bzx2y1z1 - (1-py)*pz*bzx2y1z2
     &       - py*(1-pz)*bzx2y2z1 - py*pz*bzx2y2z2) / (x1-x2)

*     d(bfbx)/dy
      dybfbx = ((1-px)*(1-pz)*bxx1y1z1 + (1-px)*pz*bxx1y1z2
     &       + px*(1-pz)*bxx2y1z1 + px*pz*bxx2y1z2
     &       - (1-px)*(1-pz)*bxx1y2z1 - (1-px)*pz*bxx1y2z2
     &       - px*(1-pz)*bxx2y2z1 - px*pz*bxx2y2z2) / (y1-y2)
*     d(bfby)/dy
      dybfby = ((1-px)*(1-pz)*byx1y1z1 + (1-px)*pz*byx1y1z2
     &       + px*(1-pz)*byx2y1z1 + px*pz*byx2y1z2
     &       - (1-px)*(1-pz)*byx1y2z1 - (1-px)*pz*byx1y2z2
     &       - px*(1-pz)*byx2y2z1 - px*pz*byx2y2z2) / (y1-y2)
*     d(bfbz)/dy
      dybfbz = ((1-px)*(1-pz)*bzx1y1z1 + (1-px)*pz*bzx1y1z2
     &       + px*(1-pz)*bzx2y1z1 + px*pz*bzx2y1z2
     &       - (1-px)*(1-pz)*bzx1y2z1 - (1-px)*pz*bzx1y2z2
     &       - px*(1-pz)*bzx2y2z1 - px*pz*bzx2y2z2) / (y1-y2)

*     d(bfbx)/dz
      dzbfbx = ((1-px)*(1-py)*bxx1y1z1 + (1-px)*py*bxx1y2z1
     &       + px*(1-py)*bxx2y1z1 + px*py*bxx2y2z1
     &       - (1-px)*(1-py)*bxx1y1z2 - (1-px)*py*bxx1y2z2
     &       - px*(1-py)*bxx2y1z2 - px*py*bxx2y2z2) / (z1-z2)
*     d(bfby)/dz
      dzbfby = ((1-px)*(1-py)*byx1y1z1 + (1-px)*py*byx1y2z1
     &       + px*(1-py)*byx2y1z1 + px*py*byx2y2z1
     &       - (1-px)*(1-py)*byx1y1z2 - (1-px)*py*byx1y2z2
     &       - px*(1-py)*byx2y1z2 - px*py*byx2y2z2) / (z1-z2)
*     d(bfbz)/dz
      dzbfbz = ((1-px)*(1-py)*bzx1y1z1 + (1-px)*py*bzx1y2z1
     &       + px*(1-py)*bzx2y1z1 + px*py*bzx2y2z1
     &       - (1-px)*(1-py)*bzx1y1z2 - (1-px)*py*bzx1y2z2
     &       - px*(1-py)*bzx2y1z2 - px*py*bzx2y2z2) / (z1-z2)

      return
      end


*****************************************************************
      subroutine readmagrzlist(map1,icmap1) ! Read r-z magnetic field list
      use magmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /brzlistint/ nr,nz,ibzrflip,ibzzflip,irztype(2)

      character chin*200,cmmt*2 ! T.Sato 2019/01/13
      character map1*200

      cmmt="#$"  ! only two characters are used as comment

      open(790,file=map1(1:icmap1),status='old')
      iline=0
! initialization
      nr=1
      nz=1
      ibzrflip=0
      ibzzflip=0
      irztype(:)=1  ! free mesh (=1) or fixed mesh (<>1)

      do iheader=1,1000
       iline=iline+1
       read(790, '(a)', err=989) chin
       call chlngt(chin,200,i1,i2)  ! calculate character length
       call chcaps(chin,i1,i2,i3,cmmt) ! CAPITAL to lower character
       i1=1                         ! delete space even for 1st character
       call chcomp(chin,i1,i2,i4)   ! delete space
       if(chin(1:4).eq.'data') then
        exit
       elseif(chin(1:3).eq.'nr=') then
        read(chin(4:200),*,err=990) nr
       elseif(chin(1:3).eq.'nz=') then
        read(chin(4:200),*,err=990) nz
       elseif(chin(1:7).eq.'extendz') then
        ibzzflip=1
        ibzrflip=1
        if(chin(8:12).eq.'flip=') then
         if(chin(13:14).eq.'br'.or.chin(16:17).eq.'br') ibzrflip=-1
         if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz') ibzzflip=-1
        endif
       elseif(chin(1:1).ne.' ') then ! not comment line
        write(ErrCha,*) 'Skip unknown parameter during reading',
     &  ' magnetic field file at line =',iheader
        ErrID = 'L:3174/R:readmagrzlist/F:magtrs.f' !W10_001_002
        call ErrWrite(ErrID,ErrCha)
       endif
      enddo
!     read header end

      if(nr.le.1) goto 991
      if(nz.le.1) goto 992

      call ALLOCATE_brzlist(nr,nz)

      do ir=1,nR
       do iz=1,nZ
        iline=iline+1
        read(790, *, err=996) rtmp,ztmp,brzrlist(ir,iz),brzzlist(ir,iz)
        if(iz.eq.1) then
         brzrmesh(ir)=rtmp
         if(ir.ne.1) then
          if(brzrmesh(ir).le.brzrmesh(ir-1)) goto 999
         endif
        elseif(brzrmesh(ir).ne.rtmp) then
         goto 997
        endif
        if(ir.eq.1) then
         brzzmesh(iz)=ztmp
         if(iz.ne.1) then
          if(brzzmesh(iz).le.brzzmesh(iz-1)) goto 1000
         endif
        elseif(brzzmesh(iz).ne.ztmp) then
         goto 998
        endif
       enddo
      enddo

! check free or fixed mesh
      dR=brzrmesh(2)-brzrmesh(1)
      do ir=2,nR-1
       if(brzrmesh(ir+1)-brzrmesh(ir).ne.dR) exit ! different mesh
      enddo
      if(ir.eq.nR) then
       irztype(1)=2  ! fixed r-mesh
      else
       irztype(1)=1  ! free r-mesh
      endif

      dZ=brzzmesh(2)-brzzmesh(1)
      do iz=2,nZ-1
       if(brzzmesh(iz+1)-brzzmesh(iz).ne.dZ) exit ! different mesh
      enddo
      if(iz.eq.nZ) then
       irztype(2)=2  ! fixed r-mesh
      else
       irztype(2)=1  ! free r-mesh
      endif

      close(790)

      return

  989  write(ErrCha,*) 'End of file during reading ',map1(1:icmap1)
       ErrID = 'L:3234/R:readmagrzlist/F:magtrs.f' !E10_002_002
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  990  write(ErrCha,*) 'Error in reading header of ',map1(1:icmap1),
     & ' at line =',iheader
       ErrID = 'L:3240/R:readmagrzlist/F:magtrs.f' !E10_003_002
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  991  write(ErrCha,*) 'nR should be greater than 1 in ',map1(1:icmap1),
     & ' ; nR=',nr
       ErrID = 'L:3246/R:readmagrzlist/F:magtrs.f' !E10_004_004
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  992  write(ErrCha,*) 'nZ should be greater than 1 in ',map1(1:icmap1),
     & ' ; nZ=',nz
       ErrID = 'L:3252/R:readmagrzlist/F:magtrs.f' !E10_004_005
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  996  write(ErrCha,*) 'Error in reading Br and Bz from ',
     & map1(1:icmap1),' at ir & iz =',ir,iz
       ErrID = 'L:3258/R:readmagrzlist/F:magtrs.f' !E10_005_002
       call ErrWrite(ErrID,ErrCha)
       goto 1900

  997  write(ErrCha,*) 'Inconsistent r-mesh in ',
     & map1(1:icmap1),' at ir & iz =',ir,iz
       ErrID = 'L:3264/R:readmagrzlist/F:magtrs.f' !E10_006_004
       call ErrWrite(ErrID,ErrCha)
       goto 1900

  998  write(ErrCha,*) 'Inconsistent z-mesh in ',
     & map1(1:icmap1),' at ir & iz =',ir,iz
       ErrID = 'L:3270/R:readmagrzlist/F:magtrs.f' !E10_006_005
       call ErrWrite(ErrID,ErrCha)
       goto 1900

  999  write(ErrCha,*) 'R-mesh should be in accending order in ',
     & map1(1:icmap1),' at r(ir) & r(ir-1)',brzrmesh(ir),brzrmesh(ir-1)
       ErrID = 'L:3276/R:readmagrzlist/F:magtrs.f' !E10_007_004
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1000  write(ErrCha,*) 'Z-mesh should be in accending order in ',
     & map1(1:icmap1),' at z(iz) & z(iz-1)',brzzmesh(iz),brzzmesh(iz-1)
       ErrID = 'L:3282/R:readmagrzlist/F:magtrs.f' !E10_007_005
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1900  write(ErrCha,'('' Please check line '',i8)') iline
       ErrID = 'L:3287/R:readmagrzlist/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)

 1999  close(790)
       call parastop( 120 )

      end

************************************************************************
*                                                                      *
      subroutine readmap2(xx,yy,zz,r1,r2,z1,z2,b_mag,
     &               bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &               dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &               brr1z1, bzr1z1, brr1z2, bzr1z2,
     &               brr2z1, bzr2z1, brr2z2, bzr2z2)
*                                                                      *
************************************************************************

      use magmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /brzlistint/ nr,nz,ibzrflip,ibzzflip,irztype(2)

      smul = b_mag
      r0=brzrmesh(1)
      r99=brzrmesh(nr)
      z0=brzzmesh(1)
      z99=brzzmesh(nz)

      if(irztype(1).ne.1) dR=(r99-r0)/(nr-1)  ! fixed mesh
      if(irztype(2).ne.1) dZ=(z99-z0)/(nz-1)  ! fixed mesh

      rtmp=dsqrt(xx**2+yy**2)
      rc = min(r99,max(r0,rtmp))
      if(rtmp.gt.0.0) then
       xc = xx*rc/rtmp
       yc = yy*rc/rtmp
      else
       xc = xx
       yc = yy
      endif

      if(rc.eq.0.0) rc=1.0e-10 ! avoid NaN
      if(ibzzflip.eq.0.or.zz.ge.0.0) then ! no extended field, or particle in positive z
       zc = min(z99,max(z0,zz))
       rrr=smul
       zzr=smul
      else ! extended field
       zc = min(z99,max(z0,dabs(zz)))
       rrr=smul*ibzrflip  ! Br is flipped or not
       zzr=smul*ibzzflip  ! Bz is flipped or not
      endif

      if(irztype(1).eq.1) then ! free mesh, T.Sato 2020/12/20 bug fix
       do nr2=2,nr-1
        if(rc.lt.brzrmesh(nr2)) exit
       enddo
       nr1=nr2-1
       r1=brzrmesh(nr1)
       r2=brzrmesh(nr2)
      else
       nR1 = min(nr-1,int((rc-R0)/dR) + 1)
       nR2 = nR1 + 1
       r1 = r0 + (nr1-1) * dr
       r2 = r1 + dr
      endif

      if(irztype(2).eq.1) then ! free mesh,  T.Sato 2020/12/20 bug fix
       do nz2=2,nZ-1
        if(zc.lt.brzzmesh(nz2)) exit
       enddo
       nz1=nz2-1
       z1=brzzmesh(nz1)
       z2=brzzmesh(nz2)
      else
       nZ1 = min(nz-1,int((zc-Z0)/dZ) + 1)
       nZ2 = nZ1 + 1
       z1 = z0 + (nz1-1) * dz
       z2 = z1 + dz
      endif



*-----------------------------------------------------------------------

      !ASTOM 2019/03/27 modified
      brr1z1 = rrr * brzrlist(nR1,nz1)
      brr2z1 = rrr * brzrlist(nR2,nz1)
      brr1z2 = rrr * brzrlist(nR1,nz2)
      brr2z2 = rrr * brzrlist(nR2,nz2)

      bzr1z1 = zzr * brzzlist(nR1,nz1)
      bzr2z1 = zzr * brzzlist(nR2,nz1)
      bzr1z2 = zzr * brzzlist(nR1,nz2)
      bzr2z2 = zzr * brzzlist(nR2,nz2)


      pr = (rc -r1) / (r2 -r1)
      pz = (zc -z1) / (z2 -z1)

*    interpolation
      bfbr = (1-pz)*(1-pr)*brr1z1
     &     + (1-pz)*pr*brr2z1
     &     + pz*(1-pr)*brr1z2
     &     + pz*pr*brr2z2
      bfbz = (1-pz)*(1-pr)*bzr1z1
     &     + (1-pz)*pr*bzr2z1
     &     + pz*(1-pr)*bzr1z2
     &     + pz*pr*bzr2z2

      bfbx = bfbr * xc/rc
      bfby = bfbr * yc/rc
      bfbz = bfbz

      drbfbr = ((1-pz)*brr1z1 + pz*brr1z2
     &         - (1-pz)*brr2z1 - pz*brr2z2) / (r1-r2)
      dzbfbr = ((1-pr)*brr1z1 + pr*brr2z1
     &         - (1-pr)*brr1z2 - pr*brr2z2) / (z1-z2)

      drbfbz = ((1-pz)*bzr1z1 + pz*bzr1z2
     &         - (1-pz)*bzr2z1 - pz*bzr2z2) / (r1-r2)
      dzbfbz = ((1-pr)*bzr1z1 + pr*bzr2z1
     &         - (1-pr)*bzr1z2 - pr*bzr2z2) / (z1-z2)

*     transform coordinate
      cost = xc / rc
      sint = yc / rc

      dxbfbx = 2*xc*cost*drbfbr
      dxbfby = 2*xc*sint*drbfbr
      dxbfbz = 2*xc*drbfbz

      dybfbx = 2*yc*cost*drbfbr
      dybfby = 2*yc*sint*drbfbr
      dybfbz = 2*yc*drbfbz

      dzbfbx = cost*dzbfbr
      dzbfby = sint*dzbfbr
      dzbfbz = dzbfbz

      return
      end

*****************************************************************
      subroutine readmagxyzmap(map1,icmap1) ! Read xyz magnetic field map
      use magmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /bxyzmapint/ nx,ny,nz,ibxyzflip(3,3)
      common /bxyzmapreal/ x0,y0,z0,x99,y99,z99,dx,dy,dz

      character chin*200,cmmt*2 ! T.Sato 2019/01/13
      character map1*200

      cmmt="#$"  ! only two characters are used as comment

      open(790,file=map1(1:icmap1),status='old')
      iline=0

! initialization
      x0=0.0
      y0=0.0
      z0=0.0
      x99=0.0
      y99=0.0
      z99=0.0
      nx=0
      ny=0
      nz=0
      ixmap=1
      iymap=1
      izmap=1
      ibxyzflip(:,:)=0

      do iheader=1,1000
       iline=iline+1
       read(790, '(a)', err=989) chin
       call chlngt(chin,200,i1,i2)  ! calculate character length
       call chcaps(chin,i1,i2,i3,cmmt) ! CAPITAL to lower character
       i1=1                         ! delete space even for 1st character
       call chcomp(chin,i1,i2,i4)   ! delete space
       if(chin(1:4).eq.'data') then
        exit
       elseif(chin(1:5).eq.'xmin=') then
        read(chin(6:200),*,err=990) x0
       elseif(chin(1:5).eq.'ymin=') then
        read(chin(6:200),*,err=990) y0
       elseif(chin(1:5).eq.'zmin=') then
        read(chin(6:200),*,err=990) z0
       elseif(chin(1:3).eq.'nx=') then
        read(chin(4:200),*,err=990) nx
       elseif(chin(1:3).eq.'ny=') then
        read(chin(4:200),*,err=990) ny
       elseif(chin(1:3).eq.'nz=') then
        read(chin(4:200),*,err=990) nz
       elseif(chin(1:5).eq.'xmax=') then
        read(chin(6:200),*,err=990) x99
       elseif(chin(1:5).eq.'ymax=') then
        read(chin(6:200),*,err=990) y99
       elseif(chin(1:5).eq.'zmax=') then
        read(chin(6:200),*,err=990) z99
       elseif(chin(1:7).eq.'ibxmap=') then
        read(chin(8:200),*,err=990) ixmap
       elseif(chin(1:7).eq.'ibymap=') then
        read(chin(8:200),*,err=990) iymap
       elseif(chin(1:7).eq.'ibzmap=') then
        read(chin(8:200),*,err=990) izmap
       elseif(chin(1:7).eq.'extendx') then
        ibxyzflip(1,:)=1
        if(chin(8:12).eq.'flip=') then
         if(chin(13:14).eq.'bx'.or.chin(16:17).eq.'bx'
     &   .or.chin(19:20).eq.'bx') ibxyzflip(1,1)=-1
         if(chin(13:14).eq.'by'.or.chin(16:17).eq.'by'
     &   .or.chin(19:20).eq.'by') ibxyzflip(1,2)=-1
         if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz'
     &   .or.chin(19:20).eq.'bz') ibxyzflip(1,3)=-1
        endif
       elseif(chin(1:7).eq.'extendy') then
        ibxyzflip(2,:)=1
        if(chin(8:12).eq.'flip=') then
         if(chin(13:14).eq.'bx'.or.chin(16:17).eq.'bx'
     &   .or.chin(19:20).eq.'bx') ibxyzflip(2,1)=-1
         if(chin(13:14).eq.'by'.or.chin(16:17).eq.'by'
     &   .or.chin(19:20).eq.'by') ibxyzflip(2,2)=-1
         if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz'
     &   .or.chin(19:20).eq.'bz') ibxyzflip(2,3)=-1
        endif
       elseif(chin(1:7).eq.'extendz') then
        ibxyzflip(3,:)=1
        if(chin(8:12).eq.'flip=') then
         if(chin(13:14).eq.'bx'.or.chin(16:17).eq.'bx'
     &   .or.chin(19:20).eq.'bx') ibxyzflip(3,1)=-1
         if(chin(13:14).eq.'by'.or.chin(16:17).eq.'by'
     &   .or.chin(19:20).eq.'by') ibxyzflip(3,2)=-1
         if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz'
     &   .or.chin(19:20).eq.'bz') ibxyzflip(3,3)=-1
        endif
       elseif(chin(1:1).ne.' ') then ! not comment line
        write(ErrCha,*) 'Skip unknown parameter during reading',
     &  map1(1:icmap1),' at line =',iheader
        ErrID = 'L:3531/R:readmagxyzmap/F:magtrs.f' !W10_001_003
        call ErrWrite(ErrID,ErrCha)
       endif
      end do
!     read header end

      if(nx.le.1) goto 900
      if(ny.le.1) goto 901
      if(nz.le.1) goto 902

      if(x99.le.x0) goto 910
      if(y99.le.y0) goto 911
      if(z99.le.z0) goto 912

      dX = (X99-X0) / (nX-1)
      dY = (Y99-Y0) / (nY-1)
      dZ = (Z99-Z0) / (nZ-1)

      call ALLOCATE_bxyzmap(nx,ny,nz)

      if(ixmap.ne.0) then
       do iz=1, nZ
        do iy=1, nY
         iline=iline+1
         read(790, *, err=920) (bxmap(ix,iy,iz), ix=1,nX)  ! read Bx
        enddo
       enddo
      endif

      if(iymap.ne.0) then
       do iz=1, nZ
        do iy=1, nY
         iline=iline+1
         read(790, *, err=921) (bymap(ix,iy,iz), ix=1,nX)  ! read By
        enddo
       enddo
      endif

      if(izmap.ne.0) then
       do iz=1, nZ
        do iy=1, nY
         iline=iline+1
         read(790, *, err=922) (bzmap(ix,iy,iz), ix=1,nX)  ! read Bz
        enddo
       enddo
      endif

      close(790)

      return

  989  write(ErrCha,*) 'End of file during reading ',map1(1:icmap1)
       ErrID = 'L:3583/R:readmagxyzmap/F:magtrs.f' !E10_002_003
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  990  write(ErrCha,*) 'Error in reading header of ',map1(1:icmap1),
     & ' at line =',iheader
       ErrID = 'L:3589/R:readmagxyzmap/F:magtrs.f' !E10_003_003
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  900  write(ErrCha,*) 'nX should be greater than 1 in ',map1(1:icmap1),
     & ' ; nX=',nx
       ErrID = 'L:3595/R:readmagxyzmap/F:magtrs.f' !E10_004_006
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  901  write(ErrCha,*) 'nY should be greater than 1 in ',map1(1:icmap1),
     & ' ; nY=',ny
       ErrID = 'L:3601/R:readmagxyzmap/F:magtrs.f' !E10_004_007
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  902  write(ErrCha,*) 'nZ should be greater than 1 in ',map1(1:icmap1),
     & ' ; nZ=',nz
       ErrID = 'L:3607/R:readmagxyzmap/F:magtrs.f' !E10_004_008
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  910  write(ErrCha,*) 'x1 should be greater than x0 in ',map1(1:icmap1)
     & ,' ; x1,x0=',x99,x0
       ErrID = 'L:3613/R:readmagxyzmap/F:magtrs.f' !E10_008_001
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  911  write(ErrCha,*) 'y1 should be greater than y0 in ',map1(1:icmap1)
     & ,' ; y1,y0=',y99,y0
       ErrID = 'L:3619/R:readmagxyzmap/F:magtrs.f' !E10_008_002
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  912  write(ErrCha,*) 'z1 should be greater than z0 in ',map1(1:icmap1)
     & ,' ; z1,z0=',z99,z0
       ErrID = 'L:3625/R:readmagxyzmap/F:magtrs.f' !E10_008_003
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  920  write(ErrCha,*) 'Error in reading Bx from ',map1(1:icmap1),
     & ' at iz=',iz
       ErrID = 'L:3631/R:readmagxyzmap/F:magtrs.f' !E10_009_001
       call ErrWrite(ErrID,ErrCha)
       goto 1900

  921  write(ErrCha,*) 'Error in reading By from ',map1(1:icmap1),
     & ' at iz=',iz
       ErrID = 'L:3637/R:readmagxyzmap/F:magtrs.f' !E10_009_002
       call ErrWrite(ErrID,ErrCha)
       goto 1900

  922  write(ErrCha,*) 'Error in reading Bz from ',map1(1:icmap1),
     & ' at ir=',ir
       ErrID = 'L:3643/R:readmagxyzmap/F:magtrs.f' !E10_009_003
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1900  write(ErrCha,'('' Please check line '',i8)') iline
       ErrID = 'L:3648/R:readmagxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)

 1999  close(790)
       call parastop( 120 )

      end

************************************************************************
*                                                                      *
      subroutine readmap3(xx,yy,zz,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
*                                                                      *
************************************************************************

      use magmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /bxyzmapint/ nx,ny,nz,ibxyzflip(3,3)
      common /bxyzmapreal/ x0,y0,z0,x99,y99,z99,dx,dy,dz


*-----------------------------------------------------------------------
      smul = b_mag

      if(ibxyzflip(1,1).eq.0.or.xx.ge.0.0) then ! no extended field, or particle in positive x
       xc = min(x99,max(x0,xx))
       xxr=smul
       yyr=smul
       zzr=smul
      else ! extended field
       xc = min(x99,max(x0,dabs(xx)))
       xxr=smul*ibxyzflip(1,1)
       yyr=smul*ibxyzflip(1,2)
       zzr=smul*ibxyzflip(1,3)
      endif
      if(ibxyzflip(2,1).eq.0.or.yy.ge.0.0) then ! no extended field, or particle in positive y
       yc = min(y99,max(y0,yy))
      else ! extended field, flip again if necessary
       yc = min(y99,max(y0,dabs(yy)))
       xxr=xxr*ibxyzflip(2,1)
       yyr=yyr*ibxyzflip(2,2)
       zzr=zzr*ibxyzflip(2,3)
      endif

      if(ibxyzflip(3,1).eq.0.or.zz.ge.0.0) then ! no extended field, or particle in positive z
       zc = min(z99,max(z0,zz))
      else ! extended field, flip again if necessary
       zc = min(z99,max(z0,dabs(zz)))
       xxr=xxr*ibxyzflip(3,1)
       yyr=yyr*ibxyzflip(3,2)
       zzr=zzr*ibxyzflip(3,3)
      endif

      nX1 = min(nx-1,int((xc-x0)/dX) + 1)
      nX2 = nX1 + 1
      nY1 = min(ny-1,int((yc-Y0)/dY) + 1)
      nY2 = nY1 + 1
      nZ1 = min(nz-1,int((zc-Z0)/dZ) + 1)
      nZ2 = nZ1 + 1

      x1 = x0 + (nx1-1) * dx
      x2 = x1 + dx
      y1 = y0 + (ny1-1) * dy
      y2 = y1 + dy
      z1 = z0 + (nz1-1) * dz
      z2 = z1 + dz

*-----------------------------------------------------------------------

      bxx1y1z1=bxmap(nx1,ny1,nz1)*xxr
      byx1y1z1=bymap(nx1,ny1,nz1)*yyr
      bzx1y1z1=bzmap(nx1,ny1,nz1)*zzr

      bxx1y1z2=bxmap(nx1,ny1,nz2)*xxr
      byx1y1z2=bymap(nx1,ny1,nz2)*yyr
      bzx1y1z2=bzmap(nx1,ny1,nz2)*zzr

      bxx1y2z1=bxmap(nx1,ny2,nz1)*xxr
      byx1y2z1=bymap(nx1,ny2,nz1)*yyr
      bzx1y2z1=bzmap(nx1,ny2,nz1)*zzr

      bxx1y2z2=bxmap(nx1,ny2,nz2)*xxr
      byx1y2z2=bymap(nx1,ny2,nz2)*yyr
      bzx1y2z2=bzmap(nx1,ny2,nz2)*zzr

      bxx2y1z1=bxmap(nx2,ny1,nz1)*xxr
      byx2y1z1=bymap(nx2,ny1,nz1)*yyr
      bzx2y1z1=bzmap(nx2,ny1,nz1)*zzr

      bxx2y1z2=bxmap(nx2,ny1,nz2)*xxr
      byx2y1z2=bymap(nx2,ny1,nz2)*yyr
      bzx2y1z2=bzmap(nx2,ny1,nz2)*zzr

      bxx2y2z1=bxmap(nx2,ny2,nz1)*xxr
      byx2y2z1=bymap(nx2,ny2,nz1)*yyr
      bzx2y2z1=bzmap(nx2,ny2,nz1)*zzr

      bxx2y2z2=bxmap(nx2,ny2,nz2)*xxr
      byx2y2z2=bymap(nx2,ny2,nz2)*yyr
      bzx2y2z2=bzmap(nx2,ny2,nz2)*zzr

      px = (xc-x1) / (x2-x1)
      py = (yc-y1) / (y2-y1)
      pz = (zc-z1) / (z2-z1)

*    interpolation
      bfbx = (1-py)*(1-px)*(1-pz)*bxx1y1z1
     &     + (1-py)*px*(1-pz)*bxx2y1z1
     &     + py*(1-px)*(1-pz)*bxx1y2z1
     &     + py*px*(1-pz)*bxx2y2z1
     &     + (1-py)*(1-px)*pz*bxx1y1z2
     &     + (1-py)*px*pz*bxx2y1z2
     &     + py*(1-px)*pz*bxx1y2z2
     &     + py*px*pz*bxx2y2z2

      bfby = (1-py)*(1-px)*(1-pz)*byx1y1z1
     &     + (1-py)*px*(1-pz)*byx2y1z1
     &     + py*(1-px)*(1-pz)*byx1y2z1
     &     + py*px*(1-pz)*byx2y2z1
     &     + (1-py)*(1-px)*pz*byx1y1z2
     &     + (1-py)*px*pz*byx2y1z2
     &     + py*(1-px)*pz*byx1y2z2
     &     + py*px*pz*byx2y2z2

      bfbz = (1-py)*(1-px)*(1-pz)*bzx1y1z1
     &     + (1-py)*px*(1-pz)*bzx2y1z1
     &     + py*(1-px)*(1-pz)*bzx1y2z1
     &     + py*px*(1-pz)*bzx2y2z1
     &     + (1-py)*(1-px)*pz*bzx1y1z2
     &     + (1-py)*px*pz*bzx2y1z2
     &     + py*(1-px)*pz*bzx1y2z2
     &     + py*px*pz*bzx2y2z2

*     d(bfbx)/dx
      dxbfbx = ((1-py)*(1-pz)*bxx1y1z1 + (1-py)*pz*bxx1y1z2
     &       + py*(1-pz)*bxx1y2z1 + py*pz*bxx1y2z2
     &       - (1-py)*(1-pz)*bxx2y1z1 - (1-py)*pz*bxx2y1z2
     &       - py*(1-pz)*bxx2y2z1 - py*pz*bxx2y2z2) / (x1-x2)
*     d(bfby)/dx
      dxbfby = ((1-py)*(1-pz)*byx1y1z1 + (1-py)*pz*byx1y1z2
     &       + py*(1-pz)*byx1y2z1 + py*pz*byx1y2z2
     &       - (1-py)*(1-pz)*byx2y1z1 - (1-py)*pz*byx2y1z2
     &       - py*(1-pz)*byx2y2z1 - py*pz*byx2y2z2) / (x1-x2)
*     d(bfbz)/dx
      dxbfbz = ((1-py)*(1-pz)*bzx1y1z1 + (1-py)*pz*bzx1y1z2
     &       + py*(1-pz)*bzx1y2z1 + py*pz*bzx1y2z2
     &       - (1-py)*(1-pz)*bzx2y1z1 - (1-py)*pz*bzx2y1z2
     &       - py*(1-pz)*bzx2y2z1 - py*pz*bzx2y2z2) / (x1-x2)

*     d(bfbx)/dy
      dybfbx = ((1-px)*(1-pz)*bxx1y1z1 + (1-px)*pz*bxx1y1z2
     &       + px*(1-pz)*bxx2y1z1 + px*pz*bxx2y1z2
     &       - (1-px)*(1-pz)*bxx1y2z1 - (1-px)*pz*bxx1y2z2
     &       - px*(1-pz)*bxx2y2z1 - px*pz*bxx2y2z2) / (y1-y2)
*     d(bfby)/dy
      dybfby = ((1-px)*(1-pz)*byx1y1z1 + (1-px)*pz*byx1y1z2
     &       + px*(1-pz)*byx2y1z1 + px*pz*byx2y1z2
     &       - (1-px)*(1-pz)*byx1y2z1 - (1-px)*pz*byx1y2z2
     &       - px*(1-pz)*byx2y2z1 - px*pz*byx2y2z2) / (y1-y2)
*     d(bfbz)/dy
      dybfbz = ((1-px)*(1-pz)*bzx1y1z1 + (1-px)*pz*bzx1y1z2
     &       + px*(1-pz)*bzx2y1z1 + px*pz*bzx2y1z2
     &       - (1-px)*(1-pz)*bzx1y2z1 - (1-px)*pz*bzx1y2z2
     &       - px*(1-pz)*bzx2y2z1 - px*pz*bzx2y2z2) / (y1-y2)

*     d(bfbx)/dz
      dzbfbx = ((1-px)*(1-py)*bxx1y1z1 + (1-px)*py*bxx1y2z1
     &       + px*(1-py)*bxx2y1z1 + px*py*bxx2y2z1
     &       - (1-px)*(1-py)*bxx1y1z2 - (1-px)*py*bxx1y2z2
     &       - px*(1-py)*bxx2y1z2 - px*py*bxx2y2z2) / (z1-z2)
*     d(bfby)/dz
      dzbfby = ((1-px)*(1-py)*byx1y1z1 + (1-px)*py*byx1y2z1
     &       + px*(1-py)*byx2y1z1 + px*py*byx2y2z1
     &       - (1-px)*(1-py)*byx1y1z2 - (1-px)*py*byx1y2z2
     &       - px*(1-py)*byx2y1z2 - px*py*byx2y2z2) / (z1-z2)
*     d(bfbz)/dz
      dzbfbz = ((1-px)*(1-py)*bzx1y1z1 + (1-px)*py*bzx1y2z1
     &       + px*(1-py)*bzx2y1z1 + px*py*bzx2y2z1
     &       - (1-px)*(1-py)*bzx1y1z2 - (1-px)*py*bzx1y2z2
     &       - px*(1-py)*bzx2y1z2 - px*py*bzx2y2z2) / (z1-z2)

      return
      end




*****************************************************************
      subroutine readmagrzmap(map1,icmap1) ! Read r-z magnetic field map
      use magmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /brzmapint/ nr,nz,ibzrflip,ibzzflip
      common /brzmapreal/ r0,z0,r99,z99,dr,dz

      character chin*200,cmmt*2 ! T.Sato 2019/01/13
      character map1*200

      cmmt="#$"  ! only two characters are used as comment

      open(790,file=map1(1:icmap1),status='old')
      iline=0
! initialization
      r0=0.0
      z0=0.0
      nr=0
      nz=0
      r99=0.0
      z99=0.0
      ibzrflip=0
      ibzzflip=0
      irmap=1
      izmap=1

      do iheader=1,1000
       iline=iline+1
       read(790, '(a)', err=989) chin
       call chlngt(chin,200,i1,i2)  ! calculate character length
       call chcaps(chin,i1,i2,i3,cmmt) ! CAPITAL to lower character
       i1=1                         ! delete space even for 1st character
       call chcomp(chin,i1,i2,i4)   ! delete space
       if(chin(1:4).eq.'data') then
        exit
       elseif(chin(1:5).eq.'rmin=') then
        read(chin(6:200),*,err=990) r0
       elseif(chin(1:5).eq.'zmin=') then
        read(chin(6:200),*,err=990) z0
       elseif(chin(1:3).eq.'nr=') then
        read(chin(4:200),*,err=990) nr
       elseif(chin(1:3).eq.'nz=') then
        read(chin(4:200),*,err=990) nz
       elseif(chin(1:5).eq.'rmax=') then
        read(chin(6:200),*,err=990) r99
       elseif(chin(1:5).eq.'zmax=') then
        read(chin(6:200),*,err=990) z99
       elseif(chin(1:7).eq.'extendz') then
        ibzzflip=1
        ibzrflip=1
        if(chin(8:12).eq.'flip=') then
         if(chin(13:14).eq.'br'.or.chin(16:17).eq.'br') ibzrflip=-1
         if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz') ibzzflip=-1
        endif
       elseif(chin(1:7).eq.'ibrmap=') then
        read(chin(8:200),*,err=990) irmap
       elseif(chin(1:7).eq.'ibzmap=') then
        read(chin(8:200),*,err=990) izmap
       elseif(chin(1:1).ne.' ') then ! not comment line
        write(ErrCha,*) 'Skip unknown parameter during reading',
     &  map1(1:icmap1),' at line =',iheader
        ErrID = 'L:3911/R:readmagrzmap/F:magtrs.f' !W10_001_004
        call ErrWrite(ErrID,ErrCha)
       endif
      end do
!     read header end

      if(nr.le.1) goto 991
      if(nz.le.1) goto 992
      if(r99.le.0.0.or.r99.le.r0) goto 993
      if(z99.le.z0) goto 994

      dR = (R99-R0) / (nR-1)
      dZ = (Z99-Z0) / (nZ-1)

      call ALLOCATE_brzmap(nr,nz)

      if(irmap.eq.1) then
       do iz=1, nZ
        iline=iline+1
        read(790, *, err=996) (brzrmap(ir,iz), ir=1,nR)
       enddo
      endif

      if(izmap.eq.1) then
       do iz=1, nZ
        iline=iline+1
        read(790, *, err=997) (brzzmap(ir,iz), ir=1,nR)
       enddo
      endif

      close(790)

      return

  989  write(ErrCha,*) 'End of file during reading ',map1(1:icmap1)
       ErrID = 'L:3946/R:readmagrzmap/F:magtrs.f' !E10_002_004
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  990  write(ErrCha,*) 'Error in reading header of ',map1(1:icmap1),
     & ' at line =',iheader
       ErrID = 'L:3952/R:readmagrzmap/F:magtrs.f' !E10_003_004
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  991  write(ErrCha,*) 'nR should be greater than 1 in ',map1(1:icmap1)
     & ,' ; nR=',nr
       ErrID = 'L:3958/R:readmagrzmap/F:magtrs.f' !E10_004_009
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  992  write(ErrCha,*) 'nZ should be greater than 1 in ',map1(1:icmap1)
     & ,' ; nZ=',nz
       ErrID = 'L:3964/R:readmagrzmap/F:magtrs.f' !E10_004_010
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  993  write(ErrCha,*) 'r1 should be greater than 0 (or r0) in '
     & ,map1(1:icmap1),' ; r1,r0=',r99,r0
       ErrID = 'L:3970/R:readmagrzmap/F:magtrs.f' !E10_008_004
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  994  write(ErrCha,*) 'z1 should be greater than z0 in ',map1(1:icmap1)
     & ,' ; z0,z1=',z99,z0
       ErrID = 'L:3976/R:readmagrzmap/F:magtrs.f' !E10_008_005
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  996  write(ErrCha,*) 'Error in reading Br from ',map1(1:icmap1),
     & ' at iz=',iz
       ErrID = 'L:3982/R:readmagrzmap/F:magtrs.f' !E10_009_004
       call ErrWrite(ErrID,ErrCha)
       goto 1900

  997  write(ErrCha,*) 'Error in reading Bz from ',map1(1:icmap1),
     & ' at iz=',iz
       ErrID = 'L:3988/R:readmagrzmap/F:magtrs.f' !E10_009_005
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1900  write(ErrCha,'('' Please check line '',i8)') iline
       ErrID = 'L:3993/R:readmagrzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)

 1999  close(790)
       call parastop( 120 )

      end

************************************************************************
*                                                                      *
      subroutine readmap4(xx,yy,zz,r1,r2,z1,z2,b_mag,
     &               bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &               dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &               brr1z1, bzr1z1, brr1z2, bzr1z2,
     &               brr2z1, bzr2z1, brr2z2, bzr2z2)
*                                                                      *
************************************************************************

      use magmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'
      common /brzmapint/ nr,nz,ibzrflip,ibzzflip
      common /brzmapreal/ r0,z0,r99,z99,dr,dz

*-----------------------------------------------------------------------
      smul = b_mag
      rtmp=dsqrt(xx**2+yy**2)
      rc = min(r99,max(r0,rtmp))
      if(rtmp.gt.0.0) then
       xc = xx*rc/rtmp
       yc = yy*rc/rtmp
      else
       xc = xx
       yc = yy
      endif
      if(rc.eq.0.0) rc=1.0e-10 ! avoid NaN

      if(ibzzflip.eq.0.or.zz.ge.0.0) then ! no extended field, or particle in positive z
       zc = min(z99,max(z0,zz))
       rrr=smul
       zzr=smul
      else ! extended field
       zc = min(z99,max(z0,dabs(zz)))
       rrr=smul*ibzrflip  ! Br is flipped or not
       zzr=smul*ibzzflip  ! Bz is flipped or not
      endif

      nZ1 = min(nz-1,int((zc-Z0)/dZ) + 1)
      nZ2 = nZ1 + 1
      nR1 = min(nr-1,int((rc-R0)/dR) + 1)
      nR2 = nR1 + 1

      z1 = Z0 + (nZ1-1) * dZ
      z2 = z1 + dZ
      r1 = 0.0 + (nR1-1) * dR
      r2 = r1 + dR

*-----------------------------------------------------------------------

      brr1z1 = rrr * brzrmap(nR1,nz1)
      brr2z1 = rrr * brzrmap(nR2,nz1)
*      brr2z1 = zzr * brzrmap(nR2,nz1)
      brr1z2 = rrr * brzrmap(nR1,nz2)
      brr2z2 = rrr * brzrmap(nR2,nz2)
*      brr2z2 = zzr * brzrmap(nR2,nz2)

      bzr1z1 = zzr * brzzmap(nR1,nz1)
*      bzr1z1 = rrr * brzzmap(nR1,nz1)
      bzr2z1 = zzr * brzzmap(nR2,nz1)
      bzr1z2 = zzr * brzzmap(nR1,nz2)
*      bzr1z2 = rrr * brzzmap(nR1,nz2)
      bzr2z2 = zzr * brzzmap(nR2,nz2)


      pr = (rc -r1) / (r2 -r1)
      pz = (zc -z1) / (z2 -z1)

*    interpolation
      bfbr = (1-pz)*(1-pr)*brr1z1
     &     + (1-pz)*pr*brr2z1
     &     + pz*(1-pr)*brr1z2
     &     + pz*pr*brr2z2
      bfbz = (1-pz)*(1-pr)*bzr1z1
     &     + (1-pz)*pr*bzr2z1
     &     + pz*(1-pr)*bzr1z2
     &     + pz*pr*bzr2z2

      bfbx = bfbr * xc/rc
      bfby = bfbr * yc/rc
      bfbz = bfbz

      drbfbr = ((1-pz)*brr1z1 + pz*brr1z2
     &         - (1-pz)*brr2z1 - pz*brr2z2) / (r1-r2)
      dzbfbr = ((1-pr)*brr1z1 + pr*brr2z1
     &         - (1-pr)*brr1z2 - pr*brr2z2) / (z1-z2)

      drbfbz = ((1-pz)*bzr1z1 + pz*bzr1z2
     &         - (1-pz)*bzr2z1 - pz*bzr2z2) / (r1-r2)
      dzbfbz = ((1-pr)*bzr1z1 + pr*bzr2z1
     &         - (1-pr)*bzr1z2 - pr*bzr2z2) / (z1-z2)

*     transform coordinate
      cost = xc / rc
      sint = yc / rc

      dxbfbx = 2*xc*cost*drbfbr
      dxbfby = 2*xc*sint*drbfbr
      dxbfbz = 2*xc*drbfbz

      dybfbx = 2*yc*cost*drbfbr
      dybfby = 2*yc*sint*drbfbr
      dybfbz = 2*yc*drbfbz

      dzbfbx = cost*dzbfbr
      dzbfby = sint*dzbfbr
      dzbfbz = dzbfbz

      return

      end

*****************************************************************
      subroutine readelcxyzlist(map1,icmap1) ! Read xyz electric field list
      use elcmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /exyzlistint/nx,ny,nz,iexyzflip(3,3),ixyztype(3)

      character chin*200,cmmt*2
      character map1*200

      cmmt="#$"  ! only two characters are used as comment

      open(790,file=map1(1:icmap1),status='old')
      iline=0
! initialization
      nx=1
      ny=1
      nz=1
      iexyzflip(:,:)=0 ! flip index

      do iheader=1,1000
            iline=iline+1
            read(790, '(a)', err=989) chin
            call chlngt(chin,200,i1,i2)  ! calculate character length
            call chcaps(chin,i1,i2,i3,cmmt) ! CAPITAL to lower character
            i1=1                         ! delete space even for 1st character
            call chcomp(chin,i1,i2,i4)   ! delete space
            if(chin(1:4).eq.'data') then
            exit
            elseif(chin(1:3).eq.'nx=') then
            read(chin(4:200),*,err=990) nx
            elseif(chin(1:3).eq.'ny=') then
            read(chin(4:200),*,err=990) ny
            elseif(chin(1:3).eq.'nz=') then
            read(chin(4:200),*,err=990) nz
            elseif(chin(1:7).eq.'extendx') then
            iexyzflip(1,:)=1
            if(chin(8:12).eq.'flip=') then
            if(chin(13:14).eq.'bx'.or.chin(16:17).eq.'bx'
     &   .or.chin(19:20).eq.'bx') iexyzflip(1,1)=-1
            if(chin(13:14).eq.'by'.or.chin(16:17).eq.'by'
     &   .or.chin(19:20).eq.'by') iexyzflip(1,2)=-1
            if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz'
     &   .or.chin(19:20).eq.'bz') iexyzflip(1,3)=-1
            endif
            elseif(chin(1:7).eq.'extendy') then
            iexyzflip(2,:)=1
            if(chin(8:12).eq.'flip=') then
            if(chin(13:14).eq.'bx'.or.chin(16:17).eq.'bx'
     &   .or.chin(19:20).eq.'bx') iexyzflip(2,1)=-1
            if(chin(13:14).eq.'by'.or.chin(16:17).eq.'by'
     &   .or.chin(19:20).eq.'by') iexyzflip(2,2)=-1
            if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz'
     &   .or.chin(19:20).eq.'bz') iexyzflip(2,3)=-1
            endif
            elseif(chin(1:7).eq.'extendz') then
            iexyzflip(3,:)=1
            if(chin(8:12).eq.'flip=') then
            if(chin(13:14).eq.'bx'.or.chin(16:17).eq.'bx'
     &   .or.chin(19:20).eq.'bx') iexyzflip(3,1)=-1
            if(chin(13:14).eq.'by'.or.chin(16:17).eq.'by'
     &   .or.chin(19:20).eq.'by') iexyzflip(3,2)=-1
            if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz'
     &   .or.chin(19:20).eq.'bz') iexyzflip(3,3)=-1
            endif
            elseif(chin(1:1).ne.' ') then ! not comment line
            write(ErrCha,*) 'Skip unknown parameter during reading',
     &  ' electric field file at line =',iheader
            ErrID = 'L:4185/R:readelcxyzlist/F:magtrs.f'
            call ErrWrite(ErrID,ErrCha)
            endif
      enddo
!     read header end

      if(nx.le.1) goto 991
      if(ny.le.1) goto 992
      if(nz.le.1) goto 993

      call ALLOCATE_exyzlist(nx,ny,nz)

      do ix=1,nx
            do iy=1,ny
            do iz=1,nz
            iline=iline+1
            read(790, *, err=996) xtmp,ytmp,ztmp,
     &   exlist(ix,iy,iz),eylist(ix,iy,iz),ezlist(ix,iy,iz)
            if(iy.eq.1.and.iz.eq.1) then
            exmesh(ix)=xtmp
            if(ix.ne.1) then
            if(exmesh(ix).le.exmesh(ix-1)) goto 1010
            endif
            elseif(exmesh(ix).ne.xtmp) then
            goto 1000
            endif
            if(ix.eq.1.and.iz.eq.1) then
            eymesh(iy)=ytmp
            if(iy.ne.1) then
            if(eymesh(iy).le.eymesh(iy-1)) goto 1011
            endif
            elseif(eymesh(iy).ne.ytmp) then
            goto 1001
            endif
            if(ix.eq.1.and.iy.eq.1) then
            ezmesh(iz)=ztmp
            if(iz.ne.1) then
            if(ezmesh(iz).le.ezmesh(iz-1)) goto 1012
            endif
            elseif(ezmesh(iz).ne.ztmp) then
            goto 1002
            endif
            enddo
            enddo
      enddo

! check free or fixed mesh
      dx=exmesh(2)-exmesh(1)
      do ix=2,nx-1
            if(exmesh(ix+1)-exmesh(ix).ne.dx) exit ! different mesh
      enddo
      if(ix.eq.nx) then
            ixyztype(1)=2  ! fixed x-mesh
      else
            ixyztype(1)=1  ! free x-mesh
      endif

      dy=eymesh(2)-eymesh(1)
      do iy=2,ny-1
            if(eymesh(iy+1)-eymesh(iy).ne.dy) exit ! different mesh
      enddo
      if(iy.eq.ny) then
            ixyztype(2)=2  ! fixed y-mesh
      else
            ixyztype(2)=1  ! free y-mesh
      endif

      dz=ezmesh(2)-ezmesh(1)
      do iz=2,nZ-1
            if(ezmesh(iz+1)-ezmesh(iz).ne.dz) exit ! different mesh
      enddo
      if(iz.eq.nz) then
            ixyztype(3)=2  ! fixed z-mesh
      else
            ixyztype(3)=1  ! free z-mesh
      endif

      close(790)

      return

  989  write(ErrCha,*) 'End of file during reading ',map1(1:icmap1)
      ErrID = 'L:4267/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)
      goto 1999

  990  write(ErrCha,*) 'Error in reading header of ',map1(1:icmap1),
     & ' at line =',iheader
      ErrID = 'L:4273/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)
      goto 1999

  991  write(ErrCha,*) 'nX should be greater than 1 in ',map1(1:icmap1),
     & ' ; nX=',nx
      ErrID = 'L:4279/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)
      goto 1999

  992  write(ErrCha,*) 'nY should be greater than 1 in ',map1(1:icmap1),
     & ' ; nY=',ny
      ErrID = 'L:4285/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)
      goto 1999

  993  write(ErrCha,*) 'nZ should be greater than 1 in ',map1(1:icmap1),
     & ' ; nZ=',nz
      ErrID = 'L:4291/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)
      goto 1999

  996  write(ErrCha,*) 'Error in reading electric fields from ',
     & map1(1:icmap1),' at ix,iy,iz =',ix,iy,iz
      ErrID = 'L:4297/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)
      goto 1900

 1000  write(ErrCha,*) 'Inconsistent x-mesh in ',
     & map1(1:icmap1),' at ix,iy,iz =',ix,iy,iz
      ErrID = 'L:4303/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)
      goto 1900

 1001  write(ErrCha,*) 'Inconsistent y-mesh in ',
     & map1(1:icmap1),' at ix,iy,iz =',ix,iy,iz
      ErrID = 'L:4309/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)
      goto 1900

 1002  write(ErrCha,*) 'Inconsistent z-mesh in ',
     & map1(1:icmap1),' at ix,iy,iz =',ix,iy,iz
      ErrID = 'L:4315/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)
      goto 1900

 1010  write(ErrCha,*) 'X-mesh should be in accending order in ',
     & map1(1:icmap1),' at x(ix) & x(ix-1)',exmesh(ix),exmesh(ix-1)
      ErrID = 'L:4321/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)
      goto 1900

 1011  write(ErrCha,*) 'Y-mesh should be in accending order in ',
     & map1(1:icmap1),' at y(iy) & y(iy-1)',eymesh(iy),eymesh(iy-1)
      ErrID = 'L:4327/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)
      goto 1900

 1012  write(ErrCha,*) 'Z-mesh should be in accending order in ',
     & map1(1:icmap1),' at z(iz) & z(iz-1)',ezmesh(iz),ezmesh(iz-1)
      ErrID = 'L:4333/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)
      goto 1900

 1900  write(ErrCha,'('' Please check line '',i8)') iline
      ErrID = 'L:4338/R:readelcxyzlist/F:magtrs.f'
      call ErrWrite(ErrID,ErrCha)

 1999  close(790)
       call parastop( 120 )

      end

************************************************************************
*                                                                      *
      subroutine readelcmap1(xx,yy,zz,x1,x2,y1,y2,z1,z2,b_elc,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
*                                                                      *
************************************************************************
      use elcmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /exyzlistint/nx,ny,nz,iexyzflip(3,3),ixyztype(3)

*-----------------------------------------------------------------------
      smul = b_elc

      x0=exmesh(1)
      x99=exmesh(nx)
      y0=eymesh(1)
      y99=eymesh(ny)
      z0=ezmesh(1)
      z99=ezmesh(nz)

      if(ixyztype(1).ne.1) dX=(x99-x0)/(nx-1)  ! fixed mesh
      if(ixyztype(2).ne.1) dY=(y99-y0)/(ny-1)  ! fixed mesh
      if(ixyztype(3).ne.1) dZ=(z99-z0)/(nz-1)  ! fixed mesh

      if(iexyzflip(1,1).eq.0.or.xx.ge.0.0) then ! no extended field, or particle in positive x
       xc = min(x99,max(x0,xx))
       xxr=smul
       yyr=smul
       zzr=smul
      else ! extended field
       xc = min(x99,max(x0,dabs(xx)))
       xxr=smul*iexyzflip(1,1)
       yyr=smul*iexyzflip(1,2)
       zzr=smul*iexyzflip(1,3)
      endif

      if(iexyzflip(2,1).eq.0.or.yy.ge.0.0) then ! no extended field, or particle in positive y
       yc = min(y99,max(y0,yy))
      else ! extended field, flip again if necessary
       yc = min(y99,max(y0,dabs(yy)))
       xxr=xxr*iexyzflip(2,1)
       yyr=yyr*iexyzflip(2,2)
       zzr=zzr*iexyzflip(2,3)
      endif

      if(iexyzflip(3,1).eq.0.or.zz.ge.0.0) then ! no extended field, or particle in positive z
       zc = min(z99,max(z0,zz))
      else ! extended field, flip again if necessary
       zc = min(z99,max(z0,dabs(zz)))
       xxr=xxr*iexyzflip(3,1)
       yyr=yyr*iexyzflip(3,2)
       zzr=zzr*iexyzflip(3,3)
      endif

      if(ixyztype(1).eq.1) then ! free mesh, T.Sato bug fix on 2020/12/20
       do nx2=2,nx-1
        if(xc.lt.exmesh(nx2)) exit
       enddo
       nx1=nx2-1
       x1=exmesh(nx1)
       x2=exmesh(nx2)
      else
       nX1 = min(nx-1,int((xc-x0)/dX) + 1)
       nX2 = nX1 + 1
       x1 = x0 + (nx1-1) * dx
       x2 = x1 + dx
      endif

      if(ixyztype(2).eq.1) then ! free mesh, T.Sato bug fix on 2020/12/20
       do ny2=2,ny-1
        if(yc.lt.eymesh(ny2)) exit
       enddo
       ny1=ny2-1
       y1=eymesh(ny1)
       y2=eymesh(ny2)
      else
       nY1 = min(ny-1,int((yc-Y0)/dY) + 1)
       nY2 = nY1 + 1
       y1 = y0 + (ny1-1) * dy
       y2 = y1 + dy
      endif

      if(ixyztype(3).eq.1) then ! free mesh, T.Sato bug fix on 2020/12/20
       do nz2=2,nz-1
        if(zc.lt.ezmesh(nz2)) exit
       enddo
       nz1=nz2-1
       z1=ezmesh(nz1)
       z2=ezmesh(nz2)
      else
       nZ1 = min(nz-1,int((zc-Z0)/dZ) + 1)
       nZ2 = nZ1 + 1
       z1 = z0 + (nz1-1) * dz
       z2 = z1 + dz
      endif

*-----------------------------------------------------------------------

      bxx1y1z1=exlist(nx1,ny1,nz1)*xxr
      byx1y1z1=eylist(nx1,ny1,nz1)*yyr
      bzx1y1z1=ezlist(nx1,ny1,nz1)*zzr

      bxx1y1z2=exlist(nx1,ny1,nz2)*xxr
      byx1y1z2=eylist(nx1,ny1,nz2)*yyr
      bzx1y1z2=ezlist(nx1,ny1,nz2)*zzr

      bxx1y2z1=exlist(nx1,ny2,nz1)*xxr
      byx1y2z1=eylist(nx1,ny2,nz1)*yyr
      bzx1y2z1=ezlist(nx1,ny2,nz1)*zzr

      bxx1y2z2=exlist(nx1,ny2,nz2)*xxr
      byx1y2z2=eylist(nx1,ny2,nz2)*yyr
      bzx1y2z2=ezlist(nx1,ny2,nz2)*zzr

      bxx2y1z1=exlist(nx2,ny1,nz1)*xxr
      byx2y1z1=eylist(nx2,ny1,nz1)*yyr
      bzx2y1z1=ezlist(nx2,ny1,nz1)*zzr

      bxx2y1z2=exlist(nx2,ny1,nz2)*xxr
      byx2y1z2=eylist(nx2,ny1,nz2)*yyr
      bzx2y1z2=ezlist(nx2,ny1,nz2)*zzr

      bxx2y2z1=exlist(nx2,ny2,nz1)*xxr
      byx2y2z1=eylist(nx2,ny2,nz1)*yyr
      bzx2y2z1=ezlist(nx2,ny2,nz1)*zzr

      bxx2y2z2=exlist(nx2,ny2,nz2)*xxr
      byx2y2z2=eylist(nx2,ny2,nz2)*yyr
      bzx2y2z2=ezlist(nx2,ny2,nz2)*zzr

      px = (xc-x1) / (x2-x1)
      py = (yc-y1) / (y2-y1)
      pz = (zc-z1) / (z2-z1)

*    interpolation
      bfbx = (1-py)*(1-px)*(1-pz)*bxx1y1z1
     &     + (1-py)*px*(1-pz)*bxx2y1z1
     &     + py*(1-px)*(1-pz)*bxx1y2z1
     &     + py*px*(1-pz)*bxx2y2z1
     &     + (1-py)*(1-px)*pz*bxx1y1z2
     &     + (1-py)*px*pz*bxx2y1z2
     &     + py*(1-px)*pz*bxx1y2z2
     &     + py*px*pz*bxx2y2z2

      bfby = (1-py)*(1-px)*(1-pz)*byx1y1z1
     &     + (1-py)*px*(1-pz)*byx2y1z1
     &     + py*(1-px)*(1-pz)*byx1y2z1
     &     + py*px*(1-pz)*byx2y2z1
     &     + (1-py)*(1-px)*pz*byx1y1z2
     &     + (1-py)*px*pz*byx2y1z2
     &     + py*(1-px)*pz*byx1y2z2
     &     + py*px*pz*byx2y2z2

      bfbz = (1-py)*(1-px)*(1-pz)*bzx1y1z1
     &     + (1-py)*px*(1-pz)*bzx2y1z1
     &     + py*(1-px)*(1-pz)*bzx1y2z1
     &     + py*px*(1-pz)*bzx2y2z1
     &     + (1-py)*(1-px)*pz*bzx1y1z2
     &     + (1-py)*px*pz*bzx2y1z2
     &     + py*(1-px)*pz*bzx1y2z2
     &     + py*px*pz*bzx2y2z2

*     d(bfbx)/dx
      dxbfbx = ((1-py)*(1-pz)*bxx1y1z1 + (1-py)*pz*bxx1y1z2
     &       + py*(1-pz)*bxx1y2z1 + py*pz*bxx1y2z2
     &       - (1-py)*(1-pz)*bxx2y1z1 - (1-py)*pz*bxx2y1z2
     &       - py*(1-pz)*bxx2y2z1 - py*pz*bxx2y2z2) / (x1-x2)
*     d(bfby)/dx
      dxbfby = ((1-py)*(1-pz)*byx1y1z1 + (1-py)*pz*byx1y1z2
     &       + py*(1-pz)*byx1y2z1 + py*pz*byx1y2z2
     &       - (1-py)*(1-pz)*byx2y1z1 - (1-py)*pz*byx2y1z2
     &       - py*(1-pz)*byx2y2z1 - py*pz*byx2y2z2) / (x1-x2)
*     d(bfbz)/dx
      dxbfbz = ((1-py)*(1-pz)*bzx1y1z1 + (1-py)*pz*bzx1y1z2
     &       + py*(1-pz)*bzx1y2z1 + py*pz*bzx1y2z2
     &       - (1-py)*(1-pz)*bzx2y1z1 - (1-py)*pz*bzx2y1z2
     &       - py*(1-pz)*bzx2y2z1 - py*pz*bzx2y2z2) / (x1-x2)

*     d(bfbx)/dy
      dybfbx = ((1-px)*(1-pz)*bxx1y1z1 + (1-px)*pz*bxx1y1z2
     &       + px*(1-pz)*bxx2y1z1 + px*pz*bxx2y1z2
     &       - (1-px)*(1-pz)*bxx1y2z1 - (1-px)*pz*bxx1y2z2
     &       - px*(1-pz)*bxx2y2z1 - px*pz*bxx2y2z2) / (y1-y2)
*     d(bfby)/dy
      dybfby = ((1-px)*(1-pz)*byx1y1z1 + (1-px)*pz*byx1y1z2
     &       + px*(1-pz)*byx2y1z1 + px*pz*byx2y1z2
     &       - (1-px)*(1-pz)*byx1y2z1 - (1-px)*pz*byx1y2z2
     &       - px*(1-pz)*byx2y2z1 - px*pz*byx2y2z2) / (y1-y2)
*     d(bfbz)/dy
      dybfbz = ((1-px)*(1-pz)*bzx1y1z1 + (1-px)*pz*bzx1y1z2
     &       + px*(1-pz)*bzx2y1z1 + px*pz*bzx2y1z2
     &       - (1-px)*(1-pz)*bzx1y2z1 - (1-px)*pz*bzx1y2z2
     &       - px*(1-pz)*bzx2y2z1 - px*pz*bzx2y2z2) / (y1-y2)

*     d(bfbx)/dz
      dzbfbx = ((1-px)*(1-py)*bxx1y1z1 + (1-px)*py*bxx1y2z1
     &       + px*(1-py)*bxx2y1z1 + px*py*bxx2y2z1
     &       - (1-px)*(1-py)*bxx1y1z2 - (1-px)*py*bxx1y2z2
     &       - px*(1-py)*bxx2y1z2 - px*py*bxx2y2z2) / (z1-z2)
*     d(bfby)/dz
      dzbfby = ((1-px)*(1-py)*byx1y1z1 + (1-px)*py*byx1y2z1
     &       + px*(1-py)*byx2y1z1 + px*py*byx2y2z1
     &       - (1-px)*(1-py)*byx1y1z2 - (1-px)*py*byx1y2z2
     &       - px*(1-py)*byx2y1z2 - px*py*byx2y2z2) / (z1-z2)
*     d(bfbz)/dz
      dzbfbz = ((1-px)*(1-py)*bzx1y1z1 + (1-px)*py*bzx1y2z1
     &       + px*(1-py)*bzx2y1z1 + px*py*bzx2y2z1
     &       - (1-px)*(1-py)*bzx1y1z2 - (1-px)*py*bzx1y2z2
     &       - px*(1-py)*bzx2y1z2 - px*py*bzx2y2z2) / (z1-z2)

      return
      end

*****************************************************************
      subroutine readelcrzlist(map1,icmap1) ! Read r-z electric field list
      use elcmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /erzlistint/ nr,nz,iezrflip,iezzflip,irztype(2)

      character chin*200,cmmt*2
      character map1*200

      cmmt="#$"  ! only two characters are used as comment

      open(790,file=map1(1:icmap1),status='old')
      iline=0
! initialization
      nr=1
      nz=1
      iezrflip=0
      iezzflip=0
      irztype(:)=1  ! free mesh (=1) or fixed mesh (<>1)

      do iheader=1,1000
       iline=iline+1
       read(790, '(a)', err=989) chin
       call chlngt(chin,200,i1,i2)  ! calculate character length
       call chcaps(chin,i1,i2,i3,cmmt) ! CAPITAL to lower character
       i1=1                         ! delete space even for 1st character
       call chcomp(chin,i1,i2,i4)   ! delete space
       if(chin(1:4).eq.'data') then
        exit
       elseif(chin(1:3).eq.'nr=') then
        read(chin(4:200),*,err=990) nr
       elseif(chin(1:3).eq.'nz=') then
        read(chin(4:200),*,err=990) nz
       elseif(chin(1:7).eq.'extendz') then
         iezzflip=1
         iezrflip=1
         if(chin(8:12).eq.'flip=') then
           if(chin(13:14).eq.'br'.or.chin(16:17).eq.'br') iezrflip=-1
           if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz') iezzflip=-1
         endif
       elseif(chin(1:1).ne.' ') then ! not comment line
        write(ErrCha,*) 'Skip unknown parameter during reading',
     &  ' electric field file at line =',iheader
        ErrID = 'L:4617/R:readelcrzlist/F:magtrs.f'
        call ErrWrite(ErrID,ErrCha)
       endif
      enddo
!     read header end

      if(nr.le.1) goto 991
      if(nz.le.1) goto 992

      call ALLOCATE_erzlist(nr,nz)

      do ir=1,nR
       do iz=1,nZ
        iline=iline+1
        read(790, *, err=996) rtmp,ztmp,erzrlist(ir,iz),erzzlist(ir,iz)
        if(iz.eq.1) then
         erzrmesh(ir)=rtmp
         if(ir.ne.1) then
          if(erzrmesh(ir).le.erzrmesh(ir-1)) goto 999
         endif
        elseif(erzrmesh(ir).ne.rtmp) then
         goto 997
        endif
        if(ir.eq.1) then
         erzzmesh(iz)=ztmp
         if(iz.ne.1) then
          if(erzzmesh(iz).le.erzzmesh(iz-1)) goto 1000
         endif
        elseif(erzzmesh(iz).ne.ztmp) then
         goto 998
        endif
       enddo
      enddo

! check free or fixed mesh
      dR=erzrmesh(2)-erzrmesh(1)
      do ir=2,nR-1
       if(erzrmesh(ir+1)-erzrmesh(ir).ne.dR) exit ! different mesh
      enddo
      if(ir.eq.nR) then
       irztype(1)=2  ! fixed r-mesh
      else
       irztype(1)=1  ! free r-mesh
      endif

      dZ=erzzmesh(2)-erzzmesh(1)
      do iz=2,nZ-1
       if(erzzmesh(iz+1)-erzzmesh(iz).ne.dZ) exit ! different mesh
      enddo
      if(iz.eq.nZ) then
       irztype(2)=2  ! fixed r-mesh
      else
       irztype(2)=1  ! free r-mesh
      endif

      close(790)

      return

  989  write(ErrCha,*) 'End of file during reading ',map1(1:icmap1)
       ErrID = 'L:4677/R:readelcrzlist/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  990  write(ErrCha,*) 'Error in reading header of ',map1(1:icmap1),
     & ' at line =',iheader
       ErrID = 'L:4683/R:readelcrzlist/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  991  write(ErrCha,*) 'nR should be greater than 1 in ',map1(1:icmap1),
     & ' ; nR=',nr
       ErrID = 'L:4689/R:readelcrzlist/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  992  write(ErrCha,*) 'nZ should be greater than 1 in ',map1(1:icmap1),
     & ' ; nZ=',nz
       ErrID = 'L:4695/R:readelcrzlist/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  996  write(ErrCha,*) 'Error in reading Br and Bz from ',
     & map1(1:icmap1),' at ir & iz =',ir,iz
       ErrID = 'L:4701/R:readelcrzlist/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1900

  997  write(ErrCha,*) 'Inconsistent r-mesh in ',
     & map1(1:icmap1),' at ir & iz =',ir,iz
       ErrID = 'L:4707/R:readelcrzlist/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1900

  998  write(ErrCha,*) 'Inconsistent z-mesh in ',
     & map1(1:icmap1),' at ir & iz =',ir,iz
       ErrID = 'L:4713/R:readelcrzlist/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1900

  999  write(ErrCha,*) 'R-mesh should be in accending order in ',
     & map1(1:icmap1),' at r(ir) & r(ir-1)',erzrmesh(ir),erzrmesh(ir-1)
       ErrID = 'L:4719/R:readelcrzlist/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1000  write(ErrCha,*) 'Z-mesh should be in accending order in ',
     & map1(1:icmap1),' at z(iz) & z(iz-1)',erzzmesh(iz),erzzmesh(iz-1)
       ErrID = 'L:4725/R:readelcrzlist/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1900  write(ErrCha,'('' Please check line '',i8)') iline
       ErrID = 'L:4730/R:readelcrzlist/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)

 1999  close(790)
       call parastop( 120 )

      end

************************************************************************
*                                                                      *
      subroutine readelcmap2(xx,yy,zz,r1,r2,z1,z2,b_elc,
     &               bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &               dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &               brr1z1, bzr1z1, brr1z2, bzr1z2,
     &               brr2z1, bzr2z1, brr2z2, bzr2z2)
*                                                                      *
************************************************************************
      use elcmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /erzlistint/ nr,nz,iezrflip,iezzflip,irztype(2)

      smul = b_elc
      r0=erzrmesh(1)
      r99=erzrmesh(nr)
      z0=erzzmesh(1)
      z99=erzzmesh(nz)

      if(irztype(1).ne.1) dR=(r99-r0)/(nr-1)  ! fixed mesh
      if(irztype(2).ne.1) dZ=(z99-z0)/(nz-1)  ! fixed mesh

      rtmp=dsqrt(xx**2+yy**2)
      rc = min(r99,max(r0,rtmp))
      if(rtmp.gt.0.0) then
       xc = xx*rc/rtmp
       yc = yy*rc/rtmp
      else
       xc = xx
       yc = yy
      endif

      if(rc.eq.0.0) rc=1.0e-10 ! avoid NaN
      if(iezzflip.eq.0.or.zz.ge.0.0) then ! no extended field, or particle in positive z
       zc = min(z99,max(z0,zz))
       rrr=smul
       zzr=smul
      else ! extended field
       zc = min(z99,max(z0,dabs(zz)))
       rrr=smul*iezrflip  ! Br is flipped or not
       zzr=smul*iezzflip  ! Bz is flipped or not
      endif

      if(irztype(1).eq.1) then ! free mesh, T.Sato 2020/12/20 bug fix
       do nr2=2,nr-1
        if(rc.lt.erzrmesh(nr2)) exit
       enddo
       nr1=nr2-1
       r1=erzrmesh(nr1)
       r2=erzrmesh(nr2)
      else
       nR1 = min(nr-1,int((rc-R0)/dR) + 1)
       nR2 = nR1 + 1
       r1 = r0 + (nr1-1) * dr
       r2 = r1 + dr
      endif

      if(irztype(2).eq.1) then ! free mesh, T.Sato 2020/12/20 bug fix
       do nz2=2,nZ-1
        if(zc.lt.erzzmesh(nz2)) exit
       enddo
       nz1=nz2-1
       z1=erzzmesh(nz1)
       z2=erzzmesh(nz2)
      else
       nZ1 = min(nz-1,int((zc-Z0)/dZ) + 1)
       nZ2 = nZ1 + 1
       z1 = z0 + (nz1-1) * dz
       z2 = z1 + dz
      endif

*-----------------------------------------------------------------------

      brr1z1 = rrr * erzrlist(nR1,nz1)
      brr2z1 = rrr * erzrlist(nR2,nz1)
      brr1z2 = rrr * erzrlist(nR1,nz2)
      brr2z2 = rrr * erzrlist(nR2,nz2)

      bzr1z1 = zzr * erzzlist(nR1,nz1)
      bzr2z1 = zzr * erzzlist(nR2,nz1)
      bzr1z2 = zzr * erzzlist(nR1,nz2)
      bzr2z2 = zzr * erzzlist(nR2,nz2)

      pr = (rc -r1) / (r2 -r1)
      pz = (zc -z1) / (z2 -z1)

*    interpolation
      bfbr = (1-pz)*(1-pr)*brr1z1
     &     + (1-pz)*pr*brr2z1
     &     + pz*(1-pr)*brr1z2
     &     + pz*pr*brr2z2
      bfbz = (1-pz)*(1-pr)*bzr1z1
     &     + (1-pz)*pr*bzr2z1
     &     + pz*(1-pr)*bzr1z2
     &     + pz*pr*bzr2z2

      bfbx = bfbr * xc/rc
      bfby = bfbr * yc/rc
      bfbz = bfbz

      drbfbr = ((1-pz)*brr1z1 + pz*brr1z2
     &         - (1-pz)*brr2z1 - pz*brr2z2) / (r1-r2)
      dzbfbr = ((1-pr)*brr1z1 + pr*brr2z1
     &         - (1-pr)*brr1z2 - pr*brr2z2) / (z1-z2)

      drbfbz = ((1-pz)*bzr1z1 + pz*bzr1z2
     &         - (1-pz)*bzr2z1 - pz*bzr2z2) / (r1-r2)
      dzbfbz = ((1-pr)*bzr1z1 + pr*bzr2z1
     &         - (1-pr)*bzr1z2 - pr*bzr2z2) / (z1-z2)

*     transform coordinate
      cost = xc / rc
      sint = yc / rc

      dxbfbx = 2*xc*cost*drbfbr
      dxbfby = 2*xc*sint*drbfbr
      dxbfbz = 2*xc*drbfbz

      dybfbx = 2*yc*cost*drbfbr
      dybfby = 2*yc*sint*drbfbr
      dybfbz = 2*yc*drbfbz

      dzbfbx = cost*dzbfbr
      dzbfby = sint*dzbfbr
      dzbfbz = dzbfbz

      return
      end

*****************************************************************
      subroutine readelcxyzmap(map1,icmap1) ! Read xyz electric field map
      use elcmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /exyzmapint/ nx,ny,nz,iexyzflip(3,3)
      common /exyzmapreal/ x0,y0,z0,x99,y99,z99,dx,dy,dz

      character chin*200,cmmt*2
      character map1*200

      cmmt="#$"  ! only two characters are used as comment

      open(790,file=map1(1:icmap1),status='old')
      iline=0

! initialization
      x0=0.0
      y0=0.0
      z0=0.0
      x99=0.0
      y99=0.0
      z99=0.0
      nx=0
      ny=0
      nz=0
      ixmap=1
      iymap=1
      izmap=1
      iexyzflip(:,:)=0

      do iheader=1,1000
       iline=iline+1
       read(790, '(a)', err=989) chin
       call chlngt(chin,200,i1,i2)  ! calculate character length
       call chcaps(chin,i1,i2,i3,cmmt) ! CAPITAL to lower character
       i1=1                         ! delete space even for 1st character
       call chcomp(chin,i1,i2,i4)   ! delete space
       if(chin(1:4).eq.'data') then
        exit
       elseif(chin(1:5).eq.'xmin=') then
        read(chin(6:200),*,err=990) x0
       elseif(chin(1:5).eq.'ymin=') then
        read(chin(6:200),*,err=990) y0
       elseif(chin(1:5).eq.'zmin=') then
        read(chin(6:200),*,err=990) z0
       elseif(chin(1:3).eq.'nx=') then
        read(chin(4:200),*,err=990) nx
       elseif(chin(1:3).eq.'ny=') then
        read(chin(4:200),*,err=990) ny
       elseif(chin(1:3).eq.'nz=') then
        read(chin(4:200),*,err=990) nz
       elseif(chin(1:5).eq.'xmax=') then
        read(chin(6:200),*,err=990) x99
       elseif(chin(1:5).eq.'ymax=') then
        read(chin(6:200),*,err=990) y99
       elseif(chin(1:5).eq.'zmax=') then
        read(chin(6:200),*,err=990) z99
       elseif(chin(1:7).eq.'ibxmap=') then
        read(chin(8:200),*,err=990) ixmap
       elseif(chin(1:7).eq.'ibymap=') then
        read(chin(8:200),*,err=990) iymap
       elseif(chin(1:7).eq.'ibzmap=') then
        read(chin(8:200),*,err=990) izmap
       elseif(chin(1:7).eq.'extendx') then
        iexyzflip(1,:)=1
        if(chin(8:12).eq.'flip=') then
         if(chin(13:14).eq.'bx'.or.chin(16:17).eq.'bx'
     &   .or.chin(19:20).eq.'bx') iexyzflip(1,1)=-1
         if(chin(13:14).eq.'by'.or.chin(16:17).eq.'by'
     &   .or.chin(19:20).eq.'by') iexyzflip(1,2)=-1
         if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz'
     &   .or.chin(19:20).eq.'bz') iexyzflip(1,3)=-1
        endif
       elseif(chin(1:7).eq.'extendy') then
         iexyzflip(2,:)=1
        if(chin(8:12).eq.'flip=') then
         if(chin(13:14).eq.'bx'.or.chin(16:17).eq.'bx'
     &   .or.chin(19:20).eq.'bx') iexyzflip(2,1)=-1
         if(chin(13:14).eq.'by'.or.chin(16:17).eq.'by'
     &   .or.chin(19:20).eq.'by') iexyzflip(2,2)=-1
         if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz'
     &   .or.chin(19:20).eq.'bz') iexyzflip(2,3)=-1
        endif
       elseif(chin(1:7).eq.'extendz') then
         iexyzflip(3,:)=1
        if(chin(8:12).eq.'flip=') then
         if(chin(13:14).eq.'bx'.or.chin(16:17).eq.'bx'
     &   .or.chin(19:20).eq.'bx') iexyzflip(3,1)=-1
         if(chin(13:14).eq.'by'.or.chin(16:17).eq.'by'
     &   .or.chin(19:20).eq.'by') iexyzflip(3,2)=-1
         if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz'
     &   .or.chin(19:20).eq.'bz') iexyzflip(3,3)=-1
        endif
       elseif(chin(1:1).ne.' ') then ! not comment line
        write(ErrCha,*) 'Skip unknown parameter during reading',
     &  map1(1:icmap1),' at line =',iheader
        ErrID = 'L:4969/R:readelcxyzmap/F:magtrs.f'
        call ErrWrite(ErrID,ErrCha)
       endif
      end do
!     read header end

      if(nx.le.1) goto 900
      if(ny.le.1) goto 901
      if(nz.le.1) goto 902

      if(x99.le.x0) goto 910
      if(y99.le.y0) goto 911
      if(z99.le.z0) goto 912

      dX = (X99-X0) / (nX-1)
      dY = (Y99-Y0) / (nY-1)
      dZ = (Z99-Z0) / (nZ-1)

      call ALLOCATE_exyzmap(nx,ny,nz)

      if(ixmap.ne.0) then
       do iz=1, nZ
        do iy=1, nY
         iline=iline+1
         read(790, *, err=920) (exmap(ix,iy,iz), ix=1,nX)  ! read Ex
        enddo
       enddo
      endif

      if(iymap.ne.0) then
       do iz=1, nZ
        do iy=1, nY
         iline=iline+1
         read(790, *, err=921) (eymap(ix,iy,iz), ix=1,nX)  ! read Ey
        enddo
       enddo
      endif

      if(izmap.ne.0) then
       do iz=1, nZ
        do iy=1, nY
         iline=iline+1
         read(790, *, err=922) (ezmap(ix,iy,iz), ix=1,nX)  ! read Ez
        enddo
       enddo
      endif

      close(790)

      return

  989  write(ErrCha,*) 'End of file during reading ',map1(1:icmap1)
       ErrID = 'L:5021/R:readelcxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  990  write(ErrCha,*) 'Error in reading header of ',map1(1:icmap1),
     & ' at line =',iheader
       ErrID = 'L:5027/R:readelcxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  900  write(ErrCha,*) 'nX should be greater than 1 in ',map1(1:icmap1),
     & ' ; nX=',nx
       ErrID = 'L:5033/R:readelcxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  901  write(ErrCha,*) 'nY should be greater than 1 in ',map1(1:icmap1),
     & ' ; nY=',ny
       ErrID = 'L:5039/R:readelcxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  902  write(ErrCha,*) 'nZ should be greater than 1 in ',map1(1:icmap1),
     & ' ; nZ=',nz
       ErrID = 'L:5045/R:readelcxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  910  write(ErrCha,*) 'x1 should be greater than x0 in ',map1(1:icmap1)
     & ,' ; x1,x0=',x99,x0
       ErrID = 'L:5051/R:readelcxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  911  write(ErrCha,*) 'y1 should be greater than y0 in ',map1(1:icmap1)
     & ,' ; y1,y0=',y99,y0
       ErrID = 'L:5057/R:readelcxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  912  write(ErrCha,*) 'z1 should be greater than z0 in ',map1(1:icmap1)
     & ,' ; z1,z0=',z99,z0
       ErrID = 'L:5063/R:readelcxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  920  write(ErrCha,*) 'Error in reading Ex from ',map1(1:icmap1),
     & ' at iz=',iz
       ErrID = 'L:5069/R:readelcxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1900

  921  write(ErrCha,*) 'Error in reading Ey from ',map1(1:icmap1),
     & ' at iz=',iz
       ErrID = 'L:5075/R:readelcxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1900

  922  write(ErrCha,*) 'Error in reading Ez from ',map1(1:icmap1),
     & ' at ir=',ir
       ErrID = 'L:5081/R:readelcxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1900  write(ErrCha,'('' Please check line '',i8)') iline
       ErrID = 'L:5086/R:readelcxyzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)

 1999  close(790)
       call parastop( 120 )

      end

************************************************************************
*                                                                      *
      subroutine readelcmap3(xx,yy,zz,x1,x2,y1,y2,z1,z2,b_elc,
     &                   bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
*                                                                      *
************************************************************************
      use elcmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /exyzmapint/ nx,ny,nz,iexyzflip(3,3)
      common /exyzmapreal/ x0,y0,z0,x99,y99,z99,dx,dy,dz


*-----------------------------------------------------------------------
      smul = b_elc

      if(iexyzflip(1,1).eq.0.or.xx.ge.0.0) then ! no extended field, or particle in positive x
       xc = min(x99,max(x0,xx))
       xxr=smul
       yyr=smul
       zzr=smul
      else ! extended field
       xc = min(x99,max(x0,dabs(xx)))
       xxr=smul*iexyzflip(1,1)
       yyr=smul*iexyzflip(1,2)
       zzr=smul*iexyzflip(1,3)
      endif
      if(iexyzflip(2,1).eq.0.or.yy.ge.0.0) then ! no extended field, or particle in positive y
       yc = min(y99,max(y0,yy))
      else ! extended field, flip again if necessary
       yc = min(y99,max(y0,dabs(yy)))
       xxr=xxr*iexyzflip(2,1)
       yyr=yyr*iexyzflip(2,2)
       zzr=zzr*iexyzflip(2,3)
      endif

      if(iexyzflip(3,1).eq.0.or.zz.ge.0.0) then ! no extended field, or particle in positive z
       zc = min(z99,max(z0,zz))
      else ! extended field, flip again if necessary
       zc = min(z99,max(z0,dabs(zz)))
       xxr=xxr*iexyzflip(3,1)
       yyr=yyr*iexyzflip(3,2)
       zzr=zzr*iexyzflip(3,3)
      endif

      nX1 = min(nx-1,int((xc-x0)/dX) + 1)
      nX2 = nX1 + 1
      nY1 = min(ny-1,int((yc-Y0)/dY) + 1)
      nY2 = nY1 + 1
      nZ1 = min(nz-1,int((zc-Z0)/dZ) + 1)
      nZ2 = nZ1 + 1

      x1 = x0 + (nx1-1) * dx
      x2 = x1 + dx
      y1 = y0 + (ny1-1) * dy
      y2 = y1 + dy
      z1 = z0 + (nz1-1) * dz
      z2 = z1 + dz

*-----------------------------------------------------------------------

      bxx1y1z1=exmap(nx1,ny1,nz1)*xxr
      byx1y1z1=eymap(nx1,ny1,nz1)*yyr
      bzx1y1z1=ezmap(nx1,ny1,nz1)*zzr

      bxx1y1z2=exmap(nx1,ny1,nz2)*xxr
      byx1y1z2=eymap(nx1,ny1,nz2)*yyr
      bzx1y1z2=ezmap(nx1,ny1,nz2)*zzr

      bxx1y2z1=exmap(nx1,ny2,nz1)*xxr
      byx1y2z1=eymap(nx1,ny2,nz1)*yyr
      bzx1y2z1=ezmap(nx1,ny2,nz1)*zzr

      bxx1y2z2=exmap(nx1,ny2,nz2)*xxr
      byx1y2z2=eymap(nx1,ny2,nz2)*yyr
      bzx1y2z2=ezmap(nx1,ny2,nz2)*zzr

      bxx2y1z1=exmap(nx2,ny1,nz1)*xxr
      byx2y1z1=eymap(nx2,ny1,nz1)*yyr
      bzx2y1z1=ezmap(nx2,ny1,nz1)*zzr

      bxx2y1z2=exmap(nx2,ny1,nz2)*xxr
      byx2y1z2=eymap(nx2,ny1,nz2)*yyr
      bzx2y1z2=ezmap(nx2,ny1,nz2)*zzr

      bxx2y2z1=exmap(nx2,ny2,nz1)*xxr
      byx2y2z1=eymap(nx2,ny2,nz1)*yyr
      bzx2y2z1=ezmap(nx2,ny2,nz1)*zzr

      bxx2y2z2=exmap(nx2,ny2,nz2)*xxr
      byx2y2z2=eymap(nx2,ny2,nz2)*yyr
      bzx2y2z2=ezmap(nx2,ny2,nz2)*zzr

      px = (xc-x1) / (x2-x1)
      py = (yc-y1) / (y2-y1)
      pz = (zc-z1) / (z2-z1)

*    interpolation
      bfbx = (1-py)*(1-px)*(1-pz)*bxx1y1z1
     &     + (1-py)*px*(1-pz)*bxx2y1z1
     &     + py*(1-px)*(1-pz)*bxx1y2z1
     &     + py*px*(1-pz)*bxx2y2z1
     &     + (1-py)*(1-px)*pz*bxx1y1z2
     &     + (1-py)*px*pz*bxx2y1z2
     &     + py*(1-px)*pz*bxx1y2z2
     &     + py*px*pz*bxx2y2z2

      bfby = (1-py)*(1-px)*(1-pz)*byx1y1z1
     &     + (1-py)*px*(1-pz)*byx2y1z1
     &     + py*(1-px)*(1-pz)*byx1y2z1
     &     + py*px*(1-pz)*byx2y2z1
     &     + (1-py)*(1-px)*pz*byx1y1z2
     &     + (1-py)*px*pz*byx2y1z2
     &     + py*(1-px)*pz*byx1y2z2
     &     + py*px*pz*byx2y2z2

      bfbz = (1-py)*(1-px)*(1-pz)*bzx1y1z1
     &     + (1-py)*px*(1-pz)*bzx2y1z1
     &     + py*(1-px)*(1-pz)*bzx1y2z1
     &     + py*px*(1-pz)*bzx2y2z1
     &     + (1-py)*(1-px)*pz*bzx1y1z2
     &     + (1-py)*px*pz*bzx2y1z2
     &     + py*(1-px)*pz*bzx1y2z2
     &     + py*px*pz*bzx2y2z2

*     d(bfbx)/dx
      dxbfbx = ((1-py)*(1-pz)*bxx1y1z1 + (1-py)*pz*bxx1y1z2
     &       + py*(1-pz)*bxx1y2z1 + py*pz*bxx1y2z2
     &       - (1-py)*(1-pz)*bxx2y1z1 - (1-py)*pz*bxx2y1z2
     &       - py*(1-pz)*bxx2y2z1 - py*pz*bxx2y2z2) / (x1-x2)
*     d(bfby)/dx
      dxbfby = ((1-py)*(1-pz)*byx1y1z1 + (1-py)*pz*byx1y1z2
     &       + py*(1-pz)*byx1y2z1 + py*pz*byx1y2z2
     &       - (1-py)*(1-pz)*byx2y1z1 - (1-py)*pz*byx2y1z2
     &       - py*(1-pz)*byx2y2z1 - py*pz*byx2y2z2) / (x1-x2)
*     d(bfbz)/dx
      dxbfbz = ((1-py)*(1-pz)*bzx1y1z1 + (1-py)*pz*bzx1y1z2
     &       + py*(1-pz)*bzx1y2z1 + py*pz*bzx1y2z2
     &       - (1-py)*(1-pz)*bzx2y1z1 - (1-py)*pz*bzx2y1z2
     &       - py*(1-pz)*bzx2y2z1 - py*pz*bzx2y2z2) / (x1-x2)

*     d(bfbx)/dy
      dybfbx = ((1-px)*(1-pz)*bxx1y1z1 + (1-px)*pz*bxx1y1z2
     &       + px*(1-pz)*bxx2y1z1 + px*pz*bxx2y1z2
     &       - (1-px)*(1-pz)*bxx1y2z1 - (1-px)*pz*bxx1y2z2
     &       - px*(1-pz)*bxx2y2z1 - px*pz*bxx2y2z2) / (y1-y2)
*     d(bfby)/dy
      dybfby = ((1-px)*(1-pz)*byx1y1z1 + (1-px)*pz*byx1y1z2
     &       + px*(1-pz)*byx2y1z1 + px*pz*byx2y1z2
     &       - (1-px)*(1-pz)*byx1y2z1 - (1-px)*pz*byx1y2z2
     &       - px*(1-pz)*byx2y2z1 - px*pz*byx2y2z2) / (y1-y2)
*     d(bfbz)/dy
      dybfbz = ((1-px)*(1-pz)*bzx1y1z1 + (1-px)*pz*bzx1y1z2
     &       + px*(1-pz)*bzx2y1z1 + px*pz*bzx2y1z2
     &       - (1-px)*(1-pz)*bzx1y2z1 - (1-px)*pz*bzx1y2z2
     &       - px*(1-pz)*bzx2y2z1 - px*pz*bzx2y2z2) / (y1-y2)

*     d(bfbx)/dz
      dzbfbx = ((1-px)*(1-py)*bxx1y1z1 + (1-px)*py*bxx1y2z1
     &       + px*(1-py)*bxx2y1z1 + px*py*bxx2y2z1
     &       - (1-px)*(1-py)*bxx1y1z2 - (1-px)*py*bxx1y2z2
     &       - px*(1-py)*bxx2y1z2 - px*py*bxx2y2z2) / (z1-z2)
*     d(bfby)/dz
      dzbfby = ((1-px)*(1-py)*byx1y1z1 + (1-px)*py*byx1y2z1
     &       + px*(1-py)*byx2y1z1 + px*py*byx2y2z1
     &       - (1-px)*(1-py)*byx1y1z2 - (1-px)*py*byx1y2z2
     &       - px*(1-py)*byx2y1z2 - px*py*byx2y2z2) / (z1-z2)
*     d(bfbz)/dz
      dzbfbz = ((1-px)*(1-py)*bzx1y1z1 + (1-px)*py*bzx1y2z1
     &       + px*(1-py)*bzx2y1z1 + px*py*bzx2y2z1
     &       - (1-px)*(1-py)*bzx1y1z2 - (1-px)*py*bzx1y2z2
     &       - px*(1-py)*bzx2y1z2 - px*py*bzx2y2z2) / (z1-z2)

      return
      end

*****************************************************************
      subroutine readelcrzmap(map1,icmap1) ! Read r-z electric field map
      use elcmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      common /erzmapint/ nr,nz,iezrflip,iezzflip
      common /erzmapreal/ r0,z0,r99,z99,dr,dz

      character chin*200,cmmt*2 ! T.Sato 2019/01/13
      character map1*200

      cmmt="#$"  ! only two characters are used as comment

      open(790,file=map1(1:icmap1),status='old')
      iline=0
! initialization
      r0=0.0
      z0=0.0
      nr=0
      nz=0
      r99=0.0
      z99=0.0
      iezrflip=0
      iezzflip=0
      irmap=1
      izmap=1

      do iheader=1,1000
       iline=iline+1
       read(790, '(a)', err=989) chin
       call chlngt(chin,200,i1,i2)  ! calculate character length
       call chcaps(chin,i1,i2,i3,cmmt) ! CAPITAL to lower character
       i1=1                         ! delete space even for 1st character
       call chcomp(chin,i1,i2,i4)   ! delete space
       if(chin(1:4).eq.'data') then
        exit
       elseif(chin(1:5).eq.'rmin=') then
        read(chin(6:200),*,err=990) r0
       elseif(chin(1:5).eq.'zmin=') then
        read(chin(6:200),*,err=990) z0
       elseif(chin(1:3).eq.'nr=') then
        read(chin(4:200),*,err=990) nr
       elseif(chin(1:3).eq.'nz=') then
        read(chin(4:200),*,err=990) nz
       elseif(chin(1:5).eq.'rmax=') then
        read(chin(6:200),*,err=990) r99
       elseif(chin(1:5).eq.'zmax=') then
        read(chin(6:200),*,err=990) z99
       elseif(chin(1:7).eq.'extendz') then
        iezzflip=1
        iezrflip=1
        if(chin(8:12).eq.'flip=') then
         if(chin(13:14).eq.'br'.or.chin(16:17).eq.'br') iezrflip=-1
         if(chin(13:14).eq.'bz'.or.chin(16:17).eq.'bz') iezzflip=-1
        endif
       elseif(chin(1:7).eq.'ibrmap=') then
        read(chin(8:200),*,err=990) irmap
       elseif(chin(1:7).eq.'ibzmap=') then
        read(chin(8:200),*,err=990) izmap
       elseif(chin(1:1).ne.' ') then ! not comment line
        write(ErrCha,*) 'Skip unknown parameter during reading',
     &  map1(1:icmap1),' at line =',iheader
        ErrID = 'L:5345/R:readelcrzmap/F:magtrs.f'
        call ErrWrite(ErrID,ErrCha)
       endif
      end do
!     read header end

      if(nr.le.1) goto 991
      if(nz.le.1) goto 992
      if(r99.le.0.0.or.r99.le.r0) goto 993
      if(z99.le.z0) goto 994

      dR = (R99-R0) / (nR-1)
      dZ = (Z99-Z0) / (nZ-1)

      call ALLOCATE_erzmap(nr,nz)

      if(irmap.eq.1) then
       do iz=1, nZ
        iline=iline+1
        read(790, *, err=996) (erzrmap(ir,iz), ir=1,nR)
       enddo
      endif

      if(izmap.eq.1) then
       do iz=1, nZ
        iline=iline+1
        read(790, *, err=997) (erzzmap(ir,iz), ir=1,nR)
       enddo
      endif

      close(790)

      return

  989  write(ErrCha,*) 'End of file during reading ',map1(1:icmap1)
       ErrID = 'L:5380/R:readelcrzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  990  write(ErrCha,*) 'Error in reading header of ',map1(1:icmap1),
     & ' at line =',iheader
       ErrID = 'L:5386/R:readelcrzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  991  write(ErrCha,*) 'nR should be greater than 1 in ',map1(1:icmap1)
     & ,' ; nR=',nr
       ErrID = 'L:5392/R:readelcrzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  992  write(ErrCha,*) 'nZ should be greater than 1 in ',map1(1:icmap1)
     & ,' ; nZ=',nz
       ErrID = 'L:5398/R:readelcrzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  993  write(ErrCha,*) 'r1 should be greater than 0 (or r0) in '
     & ,map1(1:icmap1),' ; r1,r0=',r99,r0
       ErrID = 'L:5404/R:readelcrzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  994  write(ErrCha,*) 'z1 should be greater than z0 in ',map1(1:icmap1)
     & ,' ; z0,z1=',z99,z0
       ErrID = 'L:5410/R:readelcrzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1999

  996  write(ErrCha,*) 'Error in reading Er from ',map1(1:icmap1),
     & ' at iz=',iz
       ErrID = 'L:5416/R:readelcrzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1900

  997  write(ErrCha,*) 'Error in reading Ez from ',map1(1:icmap1),
     & ' at iz=',iz
       ErrID = 'L:5422/R:readelcrzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)
       goto 1900

 1900  write(ErrCha,'('' Please check line '',i8)') iline
       ErrID = 'L:5427/R:readelcrzmap/F:magtrs.f'
       call ErrWrite(ErrID,ErrCha)

 1999  close(790)
       call parastop( 120 )

      end

************************************************************************
*                                                                      *
      subroutine readelcmap4(xx,yy,zz,r1,r2,z1,z2,b_elc,
     &               bfbx,bfby,bfbz,dxbfbx,dxbfby,dxbfbz,
     &               dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &               brr1z1, bzr1z1, brr1z2, bzr1z2,
     &               brr2z1, bzr2z1, brr2z2, bzr2z2)
*                                                                      *
************************************************************************
      use elcmod
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'
      common /erzmapint/ nr,nz,iezrflip,iezzflip
      common /erzmapreal/ r0,z0,r99,z99,dr,dz

*-----------------------------------------------------------------------
      smul = b_elc
      rtmp=dsqrt(xx**2+yy**2)
      rc = min(r99,max(r0,rtmp))
      if(rtmp.gt.0.0) then
       xc = xx*rc/rtmp
       yc = yy*rc/rtmp
      else
       xc = xx
       yc = yy
      endif
      if(rc.eq.0.0) rc=1.0e-10 ! avoid NaN

      if(iezzflip.eq.0.or.zz.ge.0.0) then ! no extended field, or particle in positive z
       zc = min(z99,max(z0,zz))
       rrr=smul
       zzr=smul
      else ! extended field
       zc = min(z99,max(z0,dabs(zz)))
       rrr=smul*iezrflip  ! Er is flipped or not
       zzr=smul*iezzflip  ! Ez is flipped or not
      endif

      nZ1 = min(nz-1,int((zc-Z0)/dZ) + 1)
      nZ2 = nZ1 + 1
      nR1 = min(nr-1,int((rc-R0)/dR) + 1)
      nR2 = nR1 + 1

      z1 = Z0 + (nZ1-1) * dZ
      z2 = z1 + dZ
      r1 = 0.0 + (nR1-1) * dR
      r2 = r1 + dR

*-----------------------------------------------------------------------

      brr1z1 = rrr * erzrmap(nR1,nz1)
      brr2z1 = rrr * erzrmap(nR2,nz1)
      brr1z2 = rrr * erzrmap(nR1,nz2)
      brr2z2 = rrr * erzrmap(nR2,nz2)

      bzr1z1 = zzr * erzzmap(nR1,nz1)
      bzr2z1 = zzr * erzzmap(nR2,nz1)
      bzr1z2 = zzr * erzzmap(nR1,nz2)
      bzr2z2 = zzr * erzzmap(nR2,nz2)

      pr = (rc -r1) / (r2 -r1)
      pz = (zc -z1) / (z2 -z1)

*    interpolation
      bfbr = (1-pz)*(1-pr)*brr1z1
     &     + (1-pz)*pr*brr2z1
     &     + pz*(1-pr)*brr1z2
     &     + pz*pr*brr2z2
      bfbz = (1-pz)*(1-pr)*bzr1z1
     &     + (1-pz)*pr*bzr2z1
     &     + pz*(1-pr)*bzr1z2
     &     + pz*pr*bzr2z2

      bfbx = bfbr * xc/rc
      bfby = bfbr * yc/rc
      bfbz = bfbz

      drbfbr = ((1-pz)*brr1z1 + pz*brr1z2
     &         - (1-pz)*brr2z1 - pz*brr2z2) / (r1-r2)
      dzbfbr = ((1-pr)*brr1z1 + pr*brr2z1
     &         - (1-pr)*brr1z2 - pr*brr2z2) / (z1-z2)

      drbfbz = ((1-pz)*bzr1z1 + pz*bzr1z2
     &         - (1-pz)*bzr2z1 - pz*bzr2z2) / (r1-r2)
      dzbfbz = ((1-pr)*bzr1z1 + pr*bzr2z1
     &         - (1-pr)*bzr1z2 - pr*bzr2z2) / (z1-z2)

*     transform coordinate
      cost = xc / rc
      sint = yc / rc

      dxbfbx = 2*xc*cost*drbfbr
      dxbfby = 2*xc*sint*drbfbr
      dxbfbz = 2*xc*drbfbz

      dybfbx = 2*yc*cost*drbfbr
      dybfby = 2*yc*sint*drbfbr
      dybfbz = 2*yc*drbfbz

      dzbfbx = cost*dzbfbr
      dzbfby = sint*dzbfbr
      dzbfbz = dzbfbz

      return

      end

************************************************************************
*                                                                      *
      subroutine chkmaperr1(x, y, z, xp, yp, zp,
     &               X0, Y0, Z0, dx, dy, dz, i)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)
      include 'err.inc'

*----------------------------------------------------------------------

      !check 0 point
      if (i .eq. 1) then
           if (x .ne. X0 .or. y .ne. Y0 .or. z .ne. Z0) then
                  !0point error
                  write(ErrCha,*) 'Starting point of magnetic ',
     &                         'map data inconsistent. ',
     &            'Check header and first line of data section.'
                  ErrID = 'L:5562/R:chkmaperr1/F:magtrs.f' !E10_003_001
                  call ErrWrite(ErrID,ErrCha)
                  close(790)
                  call parastop( 120 )
            else
            end if
      !check dx, dy, dz
      else
            dxp = abs(x - xp)
            dyp = abs(y - yp)
            dzp = abs(z - zp)

            if ((dxp .ne. dX .and. dxp .ne. 0) .or.
     &          (dyp .ne. dY .and. dyp .ne. 0) .or.
     &          (dzp .ne. dZ .and. dzp .ne. 0)) then
                 !dx size error
                   write(ErrCha,*) 'Change width is inconsistent',
     &                  'in magnetic map data coordinate at line',
     &                  i, 'under the "data" tag. ',
     &           'Check header and line ', i, ' of data section.'
                  ErrID = 'L:5582/R:chkmaperr1/F:magtrs.f' !E10_004_001
                  call ErrWrite(ErrID,ErrCha)
                  close(790)
                  call parastop( 120 )
            else
            end if
      end if
      return
      end

************************************************************************
*                                                                      *
      subroutine chkmaperr2(r, z, rp, zp, R0, Z0, dr, dz, i)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)
      include 'err.inc'

*----------------------------------------------------------------------

      !check 0 point
      if (i .eq. 1) then
           if (r .ne. R0 .or. z .ne. Z0) then
                  !0point error
                  write(ErrCha,*) 'Starting point of magnetic ',
     &                         'map data inconsistent. ',
     &            'Check header and first line of data section.'
                  ErrID = 'L:5610/R:chkmaperr2/F:magtrs.f' !E10_005_001
                  call ErrWrite(ErrID,ErrCha)
                  close(790)
                  call parastop( 120 )
            else
            end if
      !check dr, dz
      else
            drp = abs(r - rp)
            dzp = abs(z - zp)

            if ((drp .ne. dR .and. drp .ne. 0) .or.
     &          (dzp .ne. dZ .and. dzp .ne. 0)) then
                 !dx size error
                   write(ErrCha,*) 'Change width is inconsistent',
     &                  'in magnetic map data coordinate at line',
     &                  i, 'under the "data" tag. ',
     &           'Check header and line ', i, ' of data section.'
                  ErrID = 'L:5628/R:chkmaperr2/F:magtrs.f' !E10_006_001
                  call ErrWrite(ErrID,ErrCha)
                  close(790)
                  call parastop( 120 )
            else
            end if
      end if
      return
      end
