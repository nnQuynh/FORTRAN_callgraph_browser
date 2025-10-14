************************************************************************
*                                                                      *
      subroutine usrmgf1(xx,yy,zz,dmgc,cmgc,
     &                   bbx,bby,bbz,
     &                   dxx,dyx,dzx,
     &                   dxy,dyy,dzy,
     &                   dxz,dyz,dzz)
*                                                                      *
*        sample subroutine for user defined magnetic field.            *
*                                                                      *
*        Strength of additional field is given by gap[T]               *
*                                                                      *
*        input :                                                       *
*           xx, yy, zz    : position [cm]                              *
*           cmg           : strength of magnetic field [T/m^2]         *
*           dmg           : additinal z magnetic field [T]             *
*                           dmg[T] = gap[T]                            *
*                                                                      *
*        output :                                                      *
*           bbx, bby, bbz : magnetic field [T]                         *
*           dxx, dyx, dzx : dxx = d bbx / d x, ...                     *
*           dxy, dyy, dzy : dxy = d bbx / d y, ...                     *
*           dxz, dyz, dzz : derivative of magnetic field [T/cm]        *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( pi  = 3.1415926535898d0 )
      parameter ( t30 = 30.d0/180.d0*pi )

*-----------------------------------------------------------------------
*     example for theta0, theta10, theta20, theta30 data
*-----------------------------------------------------------------------

      data iread / 0 /
      save iread
!$OMP THREADPRIVATE(iread)
      data dx / 0.025 /
      data dy / 0.025 /
      data dz / 0.025 /

*-----------------------------------------------------------------------
*     read data from file
*-----------------------------------------------------------------------

         if( iread .eq. 0 ) then

            call readmgf

            iread = iread + 1

         end if

*-----------------------------------------------------------------------
*     change variables
*-----------------------------------------------------------------------

            dmg = dmgc * cmgc
            cmg = cmgc * 10000.0

*-----------------------------------------------------------------------
*     give magnetic field
*-----------------------------------------------------------------------

            call givemgf(xx,yy,zz,bbx,bby,bbz,dmg,cmg)

*-----------------------------------------------------------------------
*     give derivative of magnetic field
*-----------------------------------------------------------------------

            call givemgf(xx-dx,yy,zz,bx1,by1,bz1,dmg,cmg)
            call givemgf(xx+dx,yy,zz,bx2,by2,bz2,dmg,cmg)
            call givemgf(xx,yy-dy,zz,bx3,by3,bz3,dmg,cmg)
            call givemgf(xx,yy+dy,zz,bx4,by4,bz4,dmg,cmg)
            call givemgf(xx,yy,zz-dz,bx5,by5,bz5,dmg,cmg)
            call givemgf(xx,yy,zz+dz,bx6,by6,bz6,dmg,cmg)

            dxx = ( bx2 - bx1 ) / 2.d0 / dx
            dxy = ( bx4 - bx3 ) / 2.d0 / dy
            dxz = ( bx6 - bx5 ) / 2.d0 / dz

            dyx = ( by2 - by1 ) / 2.d0 / dx
            dyy = ( by4 - by3 ) / 2.d0 / dy
            dyz = ( by6 - by5 ) / 2.d0 / dz

            dzx = ( bz2 - bz1 ) / 2.d0 / dx
            dzy = ( bz4 - bz3 ) / 2.d0 / dy
            dzz = ( bz6 - bz5 ) / 2.d0 / dz

*-----------------------------------------------------------------------
*     change variables
*-----------------------------------------------------------------------

            bbx =  bbx / cmgc
            bby =  bby / cmgc
            bbz =  bbz / cmgc

            dxx =  dxx / cmgc
            dyx =  dyx / cmgc
            dzx =  dzx / cmgc

            dxy =  dxy / cmgc
            dyy =  dyy / cmgc
            dzy =  dzy / cmgc

            dxz =  dxz / cmgc
            dyz =  dyz / cmgc
            dzz =  dzz / cmgc

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine givemgf(xx0,yy0,zz0,
     &                   bbx,bby,bbz,dmg,cmg)
*                                                                      *
*        give magnetic field by interpolation                          *
*                                                                      *
*        input :                                                       *
*           xx, yy, zz    : position [cm]                              *
*           cmg           : strength of magnetic field [T/m^2]         *
*           dmg           : additinal z magnetic field [T]             *
*                           dmg[T] = gap[T]                            *
*                                                                      *
*        output :                                                      *
*           bbx, bby, bbz : magnetic field [T]                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'err.inc'

*-----------------------------------------------------------------------

      parameter ( pi  = 3.1415926535898d0 )

      parameter ( t10 = 10.d0/180.d0*pi )
      parameter ( t20 = 20.d0/180.d0*pi )
      parameter ( t30 = 30.d0/180.d0*pi )

      parameter ( r05 = 0.5d0 )
      parameter ( r25 = 2.5d0 )

      parameter ( zdd = 0.05d0 )
      parameter ( z00 = 0.d0 )
      parameter ( z01 = 60.d0 )
      parameter ( z02 = 200.d0 )
      parameter ( z03 = 200.d0 + z01 )

*-----------------------------------------------------------------------

      parameter ( smg = 12800.d0 )

*-----------------------------------------------------------------------

      parameter ( s00  = 0.d0 )
      parameter ( s10  = 0.173648177d0 )
      parameter ( s20  = 0.342020143d0 )
      parameter ( s30  = 0.5d0 )
      parameter ( s40  = 0.642787609d0 )
      parameter ( s60  = 0.866025403d0 )
      parameter ( s90  = 1.d0 )
      parameter ( s30m = -0.5d0 )
      parameter ( s40m = -0.642787609d0 )
      parameter ( s60m = -0.866025403d0 )
      parameter ( s90m = -1.d0 )
      parameter ( c10  = 0.984807753d0 )
      parameter ( c20  = 0.93969262d0 )
      parameter ( c40  = 0.766044443d0 )

*-----------------------------------------------------------------------

      dimension sv0(0:12), cv0(0:12)
      dimension sv1(0: 3), cv1(0: 3)
      dimension sv2(0: 3), cv2(0: 3)
      dimension sv3(0: 3)

      data sv1 / s00,s10,s20,s30 /
      data cv1 / s90,c10,c20,s60 /
      data sv2 / s00,s20,s40,s60 /
      data cv2 / s90,c20,c40,s30 /
      data sv3 / s00,s30,s60,s90 /
      data sv0 / s00,s30,s60,s90,s60,s30,s00,s30m,s60m,s90m,
     &           s60m,s30m,s00 /
      data cv0 / s90,s60,s30,s00,s30m,s60m,s90m,s60m,s30m,s00,
     &           s30,s60,s90 /

*-----------------------------------------------------------------------

      logical   exex

      character file(0:3)*20

      data file(0) / 'theta0.dat'  /
      data file(1) / 'theta10.dat' /
      data file(2) / 'theta20.dat' /
      data file(3) / 'theta30.dat' /

*-----------------------------------------------------------------------

      dimension b(0:1200,0:3,0:5,3)
      save b
!$OMP THREADPRIVATE(b)
      dimension bi(6,3), bf(3)

*-----------------------------------------------------------------------

            bbx = 0.d0
            bby = 0.d0
            bbz = 0.d0

            xx = xx0
            yy = yy0
            zz = zz0

*-----------------------------------------------------------------------
*        outer region
*-----------------------------------------------------------------------

         if( zz .lt. z00 ) zz = z00
         if( zz .gt. z03 ) zz = z03

            rr0 = sqrt( xx**2 + yy**2 )

         if( rr0 .gt. r25 ) then

            rrm = r25 / rr0
            xx  = xx * rrm
            yy  = yy * rrm

         end if

*-----------------------------------------------------------------------
*     r and theta values from xx,yy
*-----------------------------------------------------------------------

         if( yy .ge. 0.d0 .and. xx .eq. 0.d0 ) then
            th0 = pi / 2.d0
         else if( yy .lt. 0.d0 .and. xx .eq. 0.d0 ) then
            th0 = pi / 2.d0 * 3.d0
         else if( yy .ge. 0.d0 .and. xx .gt. 0.d0 ) then
            th0 = atan( yy / xx )
         else if( yy .ge. 0.d0 .and. xx .lt. 0.d0 ) then
            th0 = atan( yy / xx ) + pi
         else if( yy .lt. 0.d0 .and. xx .gt. 0.d0 ) then
            th0 = atan( yy / xx ) + 2.d0 * pi
         else if( yy .lt. 0.d0 .and. xx .lt. 0.d0 ) then
            th0 = atan( yy / xx ) + pi
         end if

*-----------------------------------------------------------------------

            ith = int( th0 / t30 )

         if( mod(ith,2) .eq. 0 ) then

            th1 = th0 - t30 * dble(ith)

         else

            th1 = t30 - ( th0 - t30 * dble(ith) )

         end if

*-----------------------------------------------------------------------

            jth = int( th1 / t10 )
            kth = jth + 1
            th2 = th1 - t10 * dble(jth)

            sin2t = sin( 2.d0 * th1 )
            cos2t = cos( 2.d0 * th1 )
            sin3t = sin( 3.d0 * th1 )

            irr = int( rr0 / r05 )
            rr1 = r05 * dble(irr)
            rr2 = rr0 - rr1
            rr3 = rr1 + r05

*-----------------------------------------------------------------------

         if( zz .ge. z01 .and. zz .le. z02 ) then

            izz =  0
            zz1 =  0.d0
            izd =  1

         else if( zz .ge. z00 .and. zz .lt. z01 ) then

            izz = int( ( z01 - zz ) / zdd ) + 1
            zz1 = -( z01 - zz ) + zdd * dble(izz)
            izd = -1

         else if( zz .gt. z02 .and. zz .le. z03 ) then

            izz = int( ( zz - z02 ) / zdd )
            zz1 = ( zz - z02 ) - zdd * dble(izz)
            izd =  1

         end if

*-----------------------------------------------------------------------

         do k = 1, 3

               if( irr .eq. 0 ) then

                  bi(1,k) = b(izz,jth,irr+1,k)
     &                      * ( rr0 / rr3 )**2

                  bi(2,k) = b(izz,jth+1,irr+1,k)
     &                      * ( rr0 / rr3 )**2

               else

                  bi(1,k) = ( b(izz,jth,irr+1,k)
     &                      * ( rr0 / rr3 )**2
     &                      * rr2
     &                      + b(izz,jth,irr  ,k)
     &                      * ( rr0 / rr1 )**2
     &                      * ( r05 - rr2 ) ) / r05

                  bi(2,k) = ( b(izz,jth+1,irr+1,k)
     &                      * ( rr0 / rr3 )**2
     &                      * rr2
     &                      + b(izz,jth+1,irr  ,k)
     &                      * ( rr0 / rr1 )**2
     &                      * ( r05 - rr2 ) ) / r05

               end if

               if( k .eq. 1 .and. jth .eq. 0 ) then

                  bi(3,k) = bi(2,k)
     &                      * sin2t / sv2(kth)

               else if( k .eq. 1 .and. jth .gt. 0 ) then

                  bi(3,k) = ( bi(2,k)
     &                      * sin2t / sv2(kth)
     &                      * th2
     &                      + bi(1,k)
     &                      * sin2t / sv2(jth)
     &                      * ( t10 - th2 ) ) / t10

               else if( k .eq. 2 ) then

                  bi(3,k) = ( bi(2,k)
     &                      * cos2t / cv2(kth)
     &                      * th2
     &                      + bi(1,k)
     &                      * cos2t / cv2(jth)
     &                      * ( t10 - th2 ) ) / t10

               else if( k .eq. 3 .and. jth .eq. 0 ) then

                  bi(3,k) = bi(2,k)
     &                      * sin3t / sv3(kth)

               else if( k .eq. 3 .and. jth .gt. 0 ) then

                  bi(3,k) = ( bi(2,k)
     &                      * sin3t / sv3(kth)
     &                      * th2
     &                      + bi(1,k)
     &                      * sin3t / sv3(jth)
     &                      * ( t10 - th2 ) ) / t10

               end if

*-----------------------------------------------------------------------

            if( izz .eq. 0 ) then

                  bf(k)   = bi(3,k)

            else

*-----------------------------------------------------------------------

               if( irr .eq. 0 ) then

                  bi(4,k) = b(izz+izd,jth,irr+1,k)
     &                      * ( rr0 / rr3 )**2

                  bi(5,k) = b(izz+izd,jth+1,irr+1,k)
     &                      * ( rr0 / rr3 )**2

               else

                  bi(4,k) = ( b(izz+izd,jth,irr+1,k)
     &                      * ( rr0 / rr3 )**2
     &                      * rr2
     &                      + b(izz+izd,jth,irr  ,k)
     &                      * ( rr0 / rr1 )**2
     &                      * ( r05 - rr2 ) ) / r05

                  bi(5,k) = ( b(izz+izd,jth+1,irr+1,k)
     &                      * ( rr0 / rr3 )**2
     &                      * rr2
     &                      + b(izz+izd,jth+1,irr  ,k)
     &                      * ( rr0 / rr1 )**2
     &                      * ( r05 - rr2 ) ) / r05

               end if

               if( k .eq. 1 .and. jth .eq. 0 ) then

                  bi(6,k) = bi(5,k)
     &                      * sin2t / sv2(kth)

               else if( k .eq. 1 .and. jth .gt. 0 ) then

                  bi(6,k) = ( bi(5,k)
     &                      * sin2t / sv2(kth)
     &                      * th2
     &                      + bi(4,k)
     &                      * sin2t / sv2(jth)
     &                      * ( t10 - th2 ) ) / t10

               else if( k .eq. 2 ) then

                  bi(6,k) = ( bi(5,k)
     &                      * cos2t / cv2(kth)
     &                      * th2
     &                      + bi(4,k)
     &                      * cos2t / cv2(jth)
     &                      * ( t10 - th2 ) ) / t10

               else if( k .eq. 3 .and. jth .eq. 0 ) then

                  bi(6,k) = bi(5,k)
     &                      * sin3t / sv3(kth)

               else if( k .eq. 3 .and. jth .gt. 0 ) then

                  bi(6,k) = ( bi(5,k)
     &                      * sin3t / sv3(kth)
     &                      * th2
     &                      + bi(4,k)
     &                      * sin3t / sv3(jth)
     &                      * ( t10 - th2 ) ) / t10

               end if

                  bf(k) = ( bi(6,k) * zz1
     &                    + bi(3,k) * ( zdd - zz1 ) ) / zdd

            end if

         end do

*-----------------------------------------------------------------------

            ixy = 1

         if( ith .eq.  2 .or. ith .eq.  3 .or.
     &       ith .eq.  6 .or. ith .eq.  7 .or.
     &       ith .eq. 10 .or. ith .eq. 11 ) ixy = -1

            izs = 1

         if( ith .eq.  2 .or. ith .eq.  3 .or.
     &       ith .eq.  6 .or. ith .eq.  6 .or.
     &       ith .eq. 10 .or. ith .eq. 11 ) izs = -1

            bf(1) = bf(1) * ixy
            bf(2) = bf(2) * ixy
            bf(3) = bf(3) * izs * izd

            if( mod(ith,2) .ne. 0 )  bf(2) = -bf(2)

            iang = ith
            if( mod(ith,2) .ne. 0 )
     &      iang = iang + 1

            coa = cv0(iang)
            sia = sv0(iang)
            bxf = bf(1)
            byf = bf(2)

            bf(1) =  bxf * coa - byf * sia
            bf(2) =  bxf * sia + byf * coa

*-----------------------------------------------------------------------

            bbx = bf(1) * cmg / smg
            bby = bf(2) * cmg / smg
            bbz = bf(3) * cmg / smg + dmg

*-----------------------------------------------------------------------
*        for debug
*-----------------------------------------------------------------------



*-----------------------------------------------------------------------

      return


************************************************************************
*                                                                      *
      entry readmgf
*                                                                      *
*        read magnetic field data                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

         do i = 0, 3

               inquire( file = file(i), exist = exex )
               if( exex .eqv. .false. ) then
                  write(ErrCha,*) 'file does not exist =>  ', file(i)
                  ErrID = 'L:524/R:givemgf/F:usrmgf1.f' !E10_001_001
                  call ErrWrite(ErrID,ErrCha)
                  stop 888
                  end if

               open(71, file = file(i), status = 'old' )

            do m = 0, 1200

                  read(71,*) ( ( b(m,i,j,k), k = 1, 3 ), j = 1, 5 )

                  do k = 1, 3
                     b(m,i,0,k) = 0.d0
                  end do

            end do

               close(71)

         end do

*-----------------------------------------------------------------------

      return
      end

