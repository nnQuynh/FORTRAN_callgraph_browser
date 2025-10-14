************************************************************************
*                                                                      *
      subroutine vtk_gshow(
     &                 isoutg, isoutgmeta,
     &                 ior,icod,iout,ires,
     &                 igser,
     &                 nx,ny,nz,vx,vy,vz,imat,
     &                 idxyreg,idxymat,
     &                 nr,mr,kr,vl,mtrns)
*                                                                      *
*       output geometry boundary for vtk.                              *
*       this is modified with gshow subroutine.                        *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit none

      include 'param.inc'

*-----------------------------------------------------------------------

      double precision :: smin
      parameter ( smin = 0.2 ) ! min is blue

*-----------------------------------------------------------------------

      integer :: isoutg, isoutgmeta
      integer :: ior, icod, iout, ires
      integer :: igser
      integer :: nx, ny, nz
      double precision :: vx, vy, vz, vl
      integer :: imat
      integer :: idxyreg(nx,ny,nz), idxymat(nx,ny,nz)
      integer :: nr, mr, kr
      integer :: mtrns

*-----------------------------------------------------------------------

      integer :: iblz1, iblz2
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      integer :: ilev1, ilev2, ilat1, ilat2
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)

      integer :: mxmat, mxmat0, mxnel
      common /kmat1a/ mxmat, mxmat0, mxnel

      integer :: idrg, idgr
      common /regdc/ idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------
      integer :: igsher, icl01, icl02
      common /igsherr/ igsher, icl01, icl02
      character i3iop*3
      character iopfile*100

      integer :: icntl, inucr
      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      dimension   vx(nx)
      dimension   vy(ny)
      dimension   vz(nz)
      dimension   imat(nx,ny,nz)
      dimension   kr(mr)
      dimension   vl(nr)

*-----------------------------------------------------------------------

      integer :: itcl
      dimension itcl(kvlmax)

      double precision :: x, xc, xd
      dimension x(3)
      dimension xc(3)
      dimension xd(3)

      double precision :: rcd
      integer :: ics
      dimension rcd(8,3)
      dimension ics(8,3)

      integer :: icdx, icdy, icd0, icd1, icd2
      dimension icdx(8)
      dimension icdy(8)

      dimension icd0(3)
      dimension icd1(3)
      dimension icd2(3)

      data icd0 / 3, 1, 2 /
      data icd1 / 1, 3, 3 /
      data icd2 / 2, 2, 1 /

      data icdx / 1,  1,  0, -1, -1, -1,  0,  1/
      data icdy / 0,  1,  1,  1,  0, -1, -1, -1/

*-----------------------------------------------------------------------

      double precision :: dlg, epsb, dps1, dps2
      data dlg / 1.0d+10 /
      data epsb / 1.0d-08 /
      data dps1 / 0.12345678d-10 /
      data dps2 / 0.23456789d-10 /

*-----------------------------------------------------------------------

      character huni4*80
      character huni6*30
      character huni5*3

      data huni5 /'),i'/


      character nlat*20
      character mlat*20

      character coln*8

      common /mtnmcl/ dmhsb(-1:kvlmax,4),
     &                dmtnm(-1:kvlmax), dmtcl(-1:kvlmax),
     &                nmtnm(-1:kvlmax), nmtcl(-1:kvlmax)
      double precision :: dmhsb
      character dmtnm*80, dmtcl*30
      integer :: nmtnm, nmtcl

      common /mtnreg/ dmhsg(-1:kvlmax),
     &                nmtng(-1:kvlmax), dmtng(-1:kvlmax)
      double precision :: dmhsg
      character dmtng*80
      integer :: nmtng

      double precision :: ercol
      dimension ercol(3,3)
      data  ercol / 1.0, 0.133, 1.0, 1.0, 1.0, 1.0, 0.6, 1.0, 0.0/

*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk

      integer :: ioh
      integer :: i, i0, ik, il, ilatx, ilatc, ii, if4, if6, iffmn
      integer :: ieid, ic, iop, ioe, inum, inrg, ici, iblz, ibond
      integer :: iicl01, iicl02, ilaty, ilatz, inne, inne1
      integer :: ivod, ixys, is, ire, ire1, ire2, ires2
      integer :: ir, ir1, ip, ipl, inne2, inne3, iblz0, icr
      integer :: j, j0, jj, jk, jl
      integer :: k, k0, kk
      integer :: l, mark, markp, mm, nmed, nmed0
      integer :: idreg, idmat
      double precision :: ddel, rval
      double precision :: sizm, sizr, srt
      double precision :: u, v, w
      double precision :: x0, x1, x2, x3, xdel, xdel0, xdel1, xdeld
      double precision :: xp, xpma, xpmn, xv, xmin, xmin0, xminf, xmaxf
      double precision :: xvmax, xvmin
      double precision :: ymin, ymin0, yminf, yp, ypma, ymaxf
      double precision :: ydel, ydel0, ydel1, ydeld, y0, y1, y2, y3
      double precision :: yvmax, yvmin, yv, ypmn
      double precision :: zmin0, zp

      integer :: idmn, idnm
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------

      if(nlat3.gt.0) call tetranowgshow !FURUTA20221209

*-----------------------------------------------------------------------
*     open temporary file
*-----------------------------------------------------------------------

            ioh = 18
            open(ioh,status='scratch',form='unformatted')

*-----------------------------------------------------------------------

         if( igser .gt. 0 .or.
     &       icntl .eq. 7 .or. icntl .eq. 8 .or.
     &       icntl .eq. 9 .or. icntl .eq. 10 ) then

            if( igser .gt. 0 ) igsher = igser

            if( icntl .eq. 7 .or. icntl .eq. 8 .or.
     &          icntl .eq. 9 .or. icntl .eq. 10 ) igsher = 2

            ioe = 19
            open(ioe,status='scratch',form='unformatted')

         end if

            inne1 = 0
            inne2 = 0
            inne3 = 0

*-----------------------------------------------------------------------
*           color of material
*-----------------------------------------------------------------------

            k = 0

         do i = 1, mxmat

            itcl(i) = 0

         end do

            ivod = 0

cFURUTA20221209 move to below
*-----------------------------------------------------------------------
*        initialization
*        icod = 1(xy), 2(yz), 3(xz)
*-----------------------------------------------------------------------

            do i = 1, 8

               rcd(i,icd0(icod)) = 0.0d0
               rcd(i,icd1(icod)) = dble( icdx(i) )
               rcd(i,icd2(icod)) = dble( icdy(i) )

               srt = sqrt( dble(icdx(i))**2 + dble(icdy(i))**2 )

               rcd(i,icd1(icod)) = dble(icdx(i)) / srt
               rcd(i,icd2(icod)) = dble(icdy(i)) / srt

               ics(i,1) = icdx(i)
               ics(i,2) = icdy(i)
               ics(i,3) = 0

            end do

*-----------------------------------------------------------------------
*        save initial mesh
*-----------------------------------------------------------------------

               ires2 = ires**2

               xmin0 = vx(1)
               ymin0 = vy(1)
               zmin0 = vz(1)
               vz(1) = vz(1) + dps1

               xdel0 = ( vx(nx) - vx(1) ) / dble( nx - 1 )
               ydel0 = ( vy(ny) - vy(1) ) / dble( ny - 1 )

               xdel1 = ( vx(nx) - vx(1) ) / dble( ires )
               ydel1 = ( vy(ny) - vy(1) ) / dble( ires )

               xdeld = ( vx(nx) - vx(1) ) / dble( nx )
               ydeld = ( vy(ny) - vy(1) ) / dble( ny )

               xdel = ( xdel1 + 2.0 * xdeld ) / dble( nx - 1 )
               ydel = ( ydel1 + 2.0 * ydeld ) / dble( ny - 1 )

               ddel = min( xdel/dble(nx), ydel/dble(ny) )**2

!FURUTA2022 initialization move to here
*-----------------------------------------------------------------------
*        initialization
*        icod = 1(xy), 2(yz), 3(xz)
*-----------------------------------------------------------------------

            do i = 1, 8

               srt = sqrt( (xdel*icdx(i))**2 + (ydel*icdy(i))**2 )

               rcd(i,icd0(icod)) = 0.0d0
               rcd(i,icd1(icod)) = xdel/srt*icdx(i)
               rcd(i,icd2(icod)) = ydel/srt*icdy(i)

               ics(i,1) = icdx(i)
               ics(i,2) = icdy(i)
               ics(i,3) = 0

            end do

*-----------------------------------------------------------------------
*     do loop 900 for fine resolution
*-----------------------------------------------------------------------

      do 900 ire = 1, ires2

               ire1 = ire - ( ire - 1 ) / ires * ires
               ire2 = ( ire - 1 ) / ires + 1

               xmin = xmin0 + xdel1 * dble( ire1 - 1 ) - xdeld + dps1
               ymin = ymin0 + ydel1 * dble( ire2 - 1 ) - ydeld + dps2

               do i = 1, nx
                  vx(i) = xmin + xdel * dble( i - 1 )
               end do

               do i = 1, ny
                  vy(i) = ymin + ydel * dble( i - 1 )
               end do

*-----------------------------------------------------------------------
*        min and max of x, y and frame
*-----------------------------------------------------------------------

               call gval(x,1,1,1,icod,nx,ny,nz,vx,vy,vz)

                  xvmin = x(icd1(icod))
                  yvmin = x(icd2(icod))

               call gval(x,nx,ny,1,icod,nx,ny,nz,vx,vy,vz)

                  xvmax = x(icd1(icod))
                  yvmax = x(icd2(icod))

*-----------------------------------------------------------------------

                  xminf = xmin + xdeld - dps1
                  yminf = ymin + ydeld - dps2
                  xmaxf = xminf + xdel1
                  ymaxf = yminf + ydel1

*-----------------------------------------------------------------------
*        initialization of imat
*-----------------------------------------------------------------------

               idreg = 0
               idmat = 0

               do i = 1, nx
               do j = 1, ny
               do k = 1, nz

                  imat(i,j,k) = 0
                  idxyreg(i,j,k) = idreg
                  idxymat(i,j,k) = idmat

               end do
               end do
               end do

*-----------------------------------------------------------------------
*     start check each cell
*-----------------------------------------------------------------------

               inum = 0

  100 continue

                     ic = 1

               do k = 1, nz
               do j = 1, ny
               do i = 1, nx

                  if( imat(i,j,k) .eq. 0 ) then

                     call gval(x,i,j,k,icod,nx,ny,nz,vx,vy,vz)

                     u  = rcd(ic,1)
                     v  = rcd(ic,2)
                     w  = rcd(ic,3)

                     mark  =  1
                     markp =  0
                     ici   = -1

                     call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)

                     if( nmed .ne. -1 .and. mark .gt. -2 ) then

                        goto 200

                     else

                        imat(i,j,k) = 90000

                     end if

                  end if

               end do
               end do
               end do

*-----------------------------------------------------------------------

               goto 900

*-----------------------------------------------------------------------
*        start new region
*-----------------------------------------------------------------------

  200    continue

               rewind ioh

               iblz0 = iblz
               nmed0 = nmed

               i0 = i
               j0 = j
               k0 = k

               ip = 0

            if( ilev1 .gt. 1 .and. ilat1(2,1) .ne. 0 ) then

               ilatc = 1

               ilatx = ilat1(3,1)
               ilaty = ilat1(4,1)
               ilatz = ilat1(5,1)

            else

               ilatc = 0

            end if

*-----------------------------------------------------------------------
*        region values : ior = 1
*-----------------------------------------------------------------------

         if( ior .eq. 1 ) then

               mm = 0

            do ir = 1, nr

               call tregck(iblz1,ilev1,ilat1,mr,kr,mm,icr)

               if( icr .ne. 0 ) then

                  ir1 = ir
                  goto 400

               end if

            end do

                  ir1 = 0

  400       continue

*-----------------------------------------------------------------------
*     region value is scaled as smin - 1.0 ( blue to red )
*-----------------------------------------------------------------------

               if( ir1 .ne. 0 ) then

                  rval = vl(ir1) * ( 1.0 - smin ) + smin

               else

                  rval = -1.0d0

               end if

         end if

*-----------------------------------------------------------------------
*           ic = 5 : first left check
*-----------------------------------------------------------------------

            ic = 5

                  ii = i + ics(ic,1)
                  jj = j + ics(ic,2)
                  kk = k + ics(ic,3)

            if( ii .gt. 0 ) then

                  call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                  u  = rcd(ic,1)
                  v  = rcd(ic,2)
                  w  = rcd(ic,3)

                  mark  = 1
                  markp = 0
                  iblz  = iblz0

               call tetra0iii !FURUTA20160701 Bugfix
               call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),u,v,w,
     &                      nmed,iblz,mark,markp,mtrns)

               if( mark .eq. 1 ) then

                  imat(i,j,k) = imat(ii,jj,kk)
                  idxyreg(i,j,k) = idxyreg(ii,jj,kk)
                  idxymat(i,j,k) = idxymat(ii,jj,kk)
                  goto 100

               else if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                  x0 = xc(icd1(icod))
                  y0 = xc(icd2(icod))

               else

                  imat(i,j,k) = -90000
                  goto 100

               end if

            else

                  x0 = x(icd1(icod))
                  y0 = x(icd2(icod))

            end if

*-----------------------------------------------------------------------
*           ic = 6 : first left-down check
*-----------------------------------------------------------------------

            ic = 6

                  ii = i + ics(ic,1)
                  jj = j + ics(ic,2)
                  kk = k + ics(ic,3)

            if( ii .gt. 0 .and. jj .gt. 0 ) then

                  call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                  u  = rcd(ic,1)
                  v  = rcd(ic,2)
                  w  = rcd(ic,3)

                  mark  = 1
                  markp = 0
                  iblz  = iblz0

               call tetra0iii !FURUTA20160701 Bugfix
               call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),u,v,w,
     &                      nmed,iblz,mark,markp,mtrns)

               if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                  x1 = xc(icd1(icod))
                  y1 = xc(icd2(icod))

                  ip = ip + 1
                  write(ioh) x1, y1

               end if

            end if

*-----------------------------------------------------------------------
*           ic = 7 : first down check
*-----------------------------------------------------------------------

            ic = 7

                  ii = i + ics(ic,1)
                  jj = j + ics(ic,2)
                  kk = k + ics(ic,3)

            if( jj .gt. 0 ) then

                  call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                  u  = rcd(ic,1)
                  v  = rcd(ic,2)
                  w  = rcd(ic,3)

                  mark  = 1
                  markp = 0
                  iblz  = iblz0

               call tetra0iii !FURUTA20160701 Bugfix
               call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),u,v,w,
     &                      nmed,iblz,mark,markp,mtrns)

               if( mark .eq. 1 ) then

                  imat(i,j,k) = imat(ii,jj,kk)
                  idxyreg(i,j,k) = idxyreg(ii,jj,kk)
                  idxymat(i,j,k) = idxymat(ii,jj,kk)
                  goto 100

               else if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                  x2 = xc(icd1(icod))
                  y2 = xc(icd2(icod))

                  ip = ip + 1
                  write(ioh) x2, y2

               else

                  imat(i,j,k) = -90000
                  goto 100

               end if

            else

                  x2 = x(icd1(icod))
                  y2 = x(icd2(icod))

                  ip = ip + 1
                  write(ioh) x2, y2

            end if

*-----------------------------------------------------------------------
*        write new region
*-----------------------------------------------------------------------

               inum = inum + 1
               inrg = inum * 10000 + iblz0
               idreg = iblz0
               idmat = idmn(nmed0)

               imat(i,j,k) = inrg
               idxyreg(i,j,k) = idreg
               idxymat(i,j,k) = idmat

               ic = 8

*-----------------------------------------------------------------------
*        eight directions movement and final write information
*-----------------------------------------------------------------------

  300    continue

*-----------------------------------------------------------------------
*           iout = 1 : boundary
*                = 2 : boundary + material
*                = 3 : boundary + material number
*                = 4 : boundary + material + material number
*                = 5 : boundary + region number
*                = 6 : boundary + material + region number
*                = 7 : boundary + lattice number
*                = 8 : boundary + material + lattice number
*-----------------------------------------------------------------------

         if( i .eq. i0 .and. j .eq. j0 .and. k .eq. k0 .and.
     &     ( ic .eq. 5 .or. ic .eq. 6 .or. ic .eq. 7 ) ) then

                  ip = ip + 1
                  write(ioh) x0, y0

            if( iout .eq. 3 .or. iout .eq. 4 .or.
     &          iout .eq. 5 .or. iout .eq. 6 .or.
     &          iout .eq. 7 .or. iout .eq. 8 ) then

                     rewind ioh

                     xpmn =  1.e33
                     xpma = -1.e33
                     ypmn =  1.e33
                     ypma = -1.e33

                  do l = 1, ip

                     read(ioh) xp, yp

                     if( xp .gt. xpma ) xpma = xp
                     if( xp .lt. xpmn ) xpmn = xp
                     if( yp .gt. ypma ) ypma = yp
                     if( yp .lt. ypmn ) ypmn = yp

                  end do

                     xp = ( xpma + xpmn ) / 2.0
                     yp = ( ypma + ypmn ) / 2.0

                     ixys = 2

*-----------------------------------------------------------------------

                     zp = vz(k0)

                     call gxval(xc,icod,xp,yp,zp)

                        ii = i0
                        jj = j0

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  =  1
                        markp =  0
                        ici   = -1

                        call gomsort(xc(1),xc(2),xc(3),u,v,w,
     &                              nmed,iblz,mark,markp,ici,mtrns)

                  if( nmed .ne. -1 .and. mark .gt. -2 .and.
     &              ( nmed .ne. nmed0 .or. iblz .ne. iblz0 ) ) then

                        call gval(xc,ii,jj,k0,icod,nx,ny,nz,vx,vy,vz)

                        xp = xc(icd1(icod))
                        yp = xc(icd2(icod))
                        ixys = 1

  120                continue

                        ii = ii + 1
                        jj = jj + 1

                        mark  =  1
                        markp =  0
                        ici   = -1

                        call gomsort(xc(1),xc(2),xc(3),u,v,w,
     &                               nmed,iblz,mark,markp,ici,mtrns)

                     if( nmed .ne. -1 .and. mark .gt. -2 .and.
     &                   nmed .eq. nmed0 .and. iblz .eq. iblz0 ) then

                        call gval(xc,ii,jj,k0,icod,nx,ny,nz,vx,vy,vz)

                        xp = xc(icd1(icod))
                        yp = xc(icd2(icod))

                        goto 120

                     else

                        ii = ( ii - 1 + i0 ) / 2
                        jj = ( jj - 1 + j0 ) / 2

                        call gval(xc,ii,jj,k0,icod,nx,ny,nz,vx,vy,vz)

                        xp = xc(icd1(icod))
                        yp = xc(icd2(icod))
                        ixys = 2

                     end if

                  end if

*-----------------------------------------------------------------------

               if( iout .eq. 3 .or. iout .eq. 4 ) then

                     if4   = nmtnm(nmed0)
                     huni4 = dmtnm(nmed0)


               else if(   iout .eq. 5 .or. iout .eq. 6 .or.
     &                ( ( iout .eq. 7 .or. iout .eq. 8 ) .and.
     &                    ilatc .eq. 0 ) ) then

                     if4   = nmtng(idgr(iblz0))
                     huni4 = dmtng(idgr(iblz0))

               else if( ( iout .eq. 7 .or. iout .eq. 8 ) .and.
     &                    ilatc .ne. 0 ) then

                     write(nlat,'(a1,i4,a1,i4,a1,i4,a1)')
     &               '(', ilatx, ',', ilaty, ',', ilatz, ')'

                        jk = 0

                     do ik = 1, 20

                        if( nlat(ik:ik) .ne. ' ') then

                           jk = jk + 1
                           mlat(jk:jk) = nlat(ik:ik)

                        end if

                     end do

               end if

            end if

*-----------------------------------------------------------------------
*           boundary is in it or not
*-----------------------------------------------------------------------

               ibond = 0

                  rewind ioh

               do l = 1, ip

                  read(ioh) x0, y0

                     if( abs( x0 - xvmin ) .le. epsb .or.
     &                   abs( x0 - xvmax ) .le. epsb .or.
     &                   abs( y0 - yvmin ) .le. epsb .or.
     &                   abs( y0 - yvmax ) .le. epsb ) then

                        ibond = 1
                        goto 230

                     end if

               end do

  230          continue


*-----------------------------------------------------------------------
*           clip:
*-----------------------------------------------------------------------

               if ( isoutg.gt.0 .and. isoutgmeta.gt.0 ) then

                  rewind(ioh)
                  call vtk_clip_to_geometry(
     &                    isoutg, isoutgmeta,
     &                    ioh, ip,
     &                    icod, vz(1), ddel)

               end if

*-----------------------------------------------------------------------

               goto 100

         end if

*-----------------------------------------------------------------------
*        four directions movement
*-----------------------------------------------------------------------

  500    continue

               if( ic .eq. 1 ) goto 501
               if( ic .eq. 2 ) goto 502
               if( ic .eq. 3 ) goto 503
               if( ic .eq. 4 ) goto 504
               if( ic .eq. 5 ) goto 505
               if( ic .eq. 6 ) goto 506
               if( ic .eq. 7 ) goto 507
               if( ic .eq. 8 ) goto 508

*-----------------------------------------------------------------------
*        ic = 1
*-----------------------------------------------------------------------

  501    ic = 1

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( ii .le. nx ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                     call tetra0iii !FURUTA20160701 Bugfix
                     call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)

               if( mark .eq. 1 ) then

                  if( imat(ii,jj,kk) .lt. 0 ) goto 110

                  if( j + 1 .le. ny ) then

                        is = 3

                        xd(1) = x(1) + dlg * rcd(is,1)
                        xd(2) = x(2) + dlg * rcd(is,2)
                        xd(3) = x(3) + dlg * rcd(is,3)

                        u  = rcd(is,1)
                        v  = rcd(is,2)
                        w  = rcd(is,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                     call tetra0iii !FURUTA20160701 Bugfix
                     call gomprpt(0,x(1),x(2),x(3),xd(1),xd(2),xd(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)

                     if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                   mark .eq. 2 ) then

                        do jl = j + 1, ny

                           if( vy(jl) .gt. xd(icd2(icod)) ) goto 601

                           if( imat(i,jl,k) .gt. 0 .and.
     &                         imat(i,jl,k) .ne. inrg ) then

                              imat(i,j,k) = -90000
                              goto 110

                           else

                              imat(i,jl,k) = inrg
                              idxyreg(i,jl,k) = idreg
                              idxymat(i,jl,k) = idmat

                           end if

                        end do

  601                   continue

                     else

                        imat(i,j,k) = -90000
                        goto 110

                     end if

                     iblz=iblz0 !FURUTA20150714 TETRA

                  end if

                        i = ii
                        j = jj
                        k = kk

                        imat(i,j,k) = inrg
                        idxyreg(i,j,k) = idreg
                        idxymat(i,j,k) = idmat

                        x(1) = xc(1)
                        x(2) = xc(2)
                        x(3) = xc(3)

                        ic = 7
                        goto 300

               else if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                  mark .eq. 2 ) then

                        ic = 2
                        goto 310

               else

                        imat(i,j,k) = -90000
                        goto 110

               end if

            else

                        xc(1) = x(1)
                        xc(2) = x(2)
                        xc(3) = x(3)

                        ic = 3
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 2
*-----------------------------------------------------------------------

  502    ic = 2

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( ii .le. nx .and. jj .le. ny ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                     call tetra0iii !FURUTA20160701 Bugfix
                     call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)

               if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                        ic = 3
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 3
*-----------------------------------------------------------------------

  503    ic = 3

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( jj .le. ny ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                     call tetra0iii !FURUTA20160701 Bugfix
                     call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)

               if( mark .eq. 1 ) then

                  if( imat(ii,jj,kk) .lt. 0 ) goto 110

                  if( i - 1 .ge. 1 ) then

                        is = 5

                        xd(1) = x(1) + dlg * rcd(is,1)
                        xd(2) = x(2) + dlg * rcd(is,2)
                        xd(3) = x(3) + dlg * rcd(is,3)

                        u  = rcd(is,1)
                        v  = rcd(is,2)
                        w  = rcd(is,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                     call tetra0iii !FURUTA20160701 Bugfix
                     call gomprpt(0,x(1),x(2),x(3),xd(1),xd(2),xd(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)

                     if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                   mark .eq. 2 ) then

                        do il = i - 1, 1, -1

                           if( vx(il) .lt. xd(icd1(icod)) ) goto 603

                           if( imat(il,j,k) .gt. 0 .and.
     &                         imat(il,j,k) .ne. inrg ) then

                              imat(i,j,k) = -90000
                              goto 110

                           else

                              imat(il,j,k) = inrg
                              idxyreg(il,j,k) = idreg
                              idxymat(il,j,k) = idmat

                           end if

                        end do

  603                   continue

                     else

                        imat(i,j,k) = -90000
                        goto 110

                     end if

                  end if

                        i = ii
                        j = jj
                        k = kk

                        imat(i,j,k) = inrg
                        idxyreg(i,j,k) = idreg
                        idxymat(i,j,k) = idmat

                        x(1) = xc(1)
                        x(2) = xc(2)
                        x(3) = xc(3)

                        ic = 1
                        goto 300

               else if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                  mark .eq. 2 ) then

                        ic = 4
                        goto 310

               else

                        imat(i,j,k) = -90000
                        goto 110

               end if

            else

                        xc(1) = x(1)
                        xc(2) = x(2)
                        xc(3) = x(3)

                        ic = 5
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 4
*-----------------------------------------------------------------------

  504    ic = 4

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( ii .gt. 0 .and. jj .le. ny ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                     call tetra0iii !FURUTA20160701 Bugfix
                     call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)

               if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                        ic = 5
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 5
*-----------------------------------------------------------------------

  505    ic = 5

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( ii .gt. 0 ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                     call tetra0iii !FURUTA20160701 Bugfix
                     call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)

               if( mark .eq. 1 ) then

                  if( imat(ii,jj,kk) .lt. 0 ) goto 110

                  if( j - 1 .ge. 1 ) then

                        is = 7

                        xd(1) = x(1) + dlg * rcd(is,1)
                        xd(2) = x(2) + dlg * rcd(is,2)
                        xd(3) = x(3) + dlg * rcd(is,3)

                        u  = rcd(is,1)
                        v  = rcd(is,2)
                        w  = rcd(is,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                     call tetra0iii !FURUTA20160701 Bugfix
                     call gomprpt(0,x(1),x(2),x(3),xd(1),xd(2),xd(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)

                     if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                   mark .eq. 2 ) then

                        do jl = j - 1, 1, -1

                           if( vy(jl) .lt. xd(icd2(icod)) ) goto 605

                           if( imat(i,jl,k) .gt. 0 .and.
     &                         imat(i,jl,k) .ne. inrg ) then

                              imat(i,j,k) = -90000
                              goto 110

                           else

                              imat(i,jl,k) = inrg
                              idxyreg(i,jl,k) = idreg
                              idxymat(i,jl,k) = idmat

                           end if

                        end do

  605                   continue

                     else

                        imat(i,j,k) = -90000
                        goto 110

                     end if

                  end if

                        i = ii
                        j = jj
                        k = kk

                        imat(i,j,k) = inrg
                        idxyreg(i,j,k) = idreg
                        idxymat(i,j,k) = idmat

                        x(1) = xc(1)
                        x(2) = xc(2)
                        x(3) = xc(3)

                        ic = 3
                        goto 300

               else if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                  mark .eq. 2 ) then

                        ic = 6
                        goto 310

               else

                        imat(i,j,k) = -90000
                        goto 110

               end if

            else

                        xc(1) = x(1)
                        xc(2) = x(2)
                        xc(3) = x(3)

                        ic = 7
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 6
*-----------------------------------------------------------------------

  506    ic = 6

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( ii .gt. 0 .and. jj .gt. 0 ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                     call tetra0iii !FURUTA20160701 Bugfix
                     call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)

               if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                        ic = 7
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 7
*-----------------------------------------------------------------------

  507    ic = 7

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( jj .gt. 0 ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                     call tetra0iii !FURUTA20160701 Bugfix
                     call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)

               if( mark .eq. 1 ) then

                  if( imat(ii,jj,kk) .lt. 0 ) goto 110

                  if( i + 1 .le. nx ) then

                        is = 1

                        xd(1) = x(1) + dlg * rcd(is,1)
                        xd(2) = x(2) + dlg * rcd(is,2)
                        xd(3) = x(3) + dlg * rcd(is,3)

                        u  = rcd(is,1)
                        v  = rcd(is,2)
                        w  = rcd(is,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                     call tetra0iii !FURUTA20160701 Bugfix
                     call gomprpt(0,x(1),x(2),x(3),xd(1),xd(2),xd(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)

                     if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                   mark .eq. 2 ) then

                        do il = i + 1, nx

                           if( vx(il) .gt. xd(icd1(icod)) ) goto 607

                           if( imat(il,j,k) .gt. 0 .and.
     &                         imat(il,j,k) .ne. inrg ) then

                              imat(i,j,k) = -90000
                              goto 110

                           else

                              imat(il,j,k) = inrg
                              idxyreg(il,j,k) = idreg
                              idxymat(il,j,k) = idmat

                           end if

                        end do

  607                   continue

                     else

                        imat(i,j,k) = -90000
                        goto 110

                     end if

                  end if

                        i = ii
                        j = jj
                        k = kk

                        imat(i,j,k) = inrg
                        idxyreg(i,j,k) = idreg
                        idxymat(i,j,k) = idmat

                        x(1) = xc(1)
                        x(2) = xc(2)
                        x(3) = xc(3)

                        ic = 5
                        goto 300

               else if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                  mark .eq. 2 ) then

                        ic = 8
                        goto 310

               else

                        imat(i,j,k) = -90000
                        goto 110

               end if

            else

                        xc(1) = x(1)
                        xc(2) = x(2)
                        xc(3) = x(3)

                        ic = 1
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 8
*-----------------------------------------------------------------------

  508    ic = 8

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( ii .le. nx .and. jj .gt. 0 ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                     call tetra0iii !FURUTA20160701 Bugfix
                     call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)

               if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                        ic = 1
                        goto 310

               end if

            end if

                  goto 501

*-----------------------------------------------------------------------

  310       continue

                  x3 = xc(icd1(icod))
                  y3 = xc(icd2(icod))

                  ip = ip + 1
                  write(ioh) x3, y3

            goto 300

*-----------------------------------------------------------------------

  110 continue

               do k = 1, nz
               do j = 1, ny
               do i = 1, nx

                  if( imat(i,j,k) .eq. inrg ) then

                     imat(i,j,k) = -90000

                  end if

               end do
               end do
               end do

               goto 100

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------
*        restore initial mesh
*-----------------------------------------------------------------------

               do i = 1, nx
                  vx(i) = xmin0 + xdel0 * dble( i - 1 )
               end do

               do i = 1, ny
                  vy(i) = ymin0 + ydel0 * dble( i - 1 )
               end do

                  vz(1) = zmin0

*-----------------------------------------------------------------------

      close( ioh )

      if( igsher .ne. 0 ) then
         close( ioe )
      end if
         igsher = 0

      if(nlat3.gt.0) call tetrafingshow !FURUTA20221209

*-----------------------------------------------------------------------

      return
      end

