!***********************************************************************
!                                                                      *
      subroutine anatal_rearrange(itaty,iMeVperu,imesh,
     &     iDaxis,ntaxis,anatalrst,
     &     nd1,    nd2,    nd3,    nd4,    nd5,nd6,nd7,   nd8,nd9,
     &             db2,dw2,db3,dw3,db4,dw4,db5,db6,db7,dw7,
     &     vl,
     &     nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &     cij,cijaxs,cijaxs2,nanataldata,ierr)
!                                                                      *
!     created by S.Hashimoto on 2021/4/27                              *
!     modified by T.Miura on 2021/8/6                                   *
!                                                                      *
!     anataldata: rearranged data of anatalrst                         *
!     fg: rearranged mesh data of all meshes                           *
!     fgaxs: mesh data of specified axis                               *
!     delvol: delta volume, increments used in calculating sum over    *
!                                                                      *
!***********************************************************************

      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol

      implicit none
!-----------------------------------------------------------------------
      include 'err.inc'
      include 'param-physcnst.inc'
!-----------------------------------------------------------------------

      integer itaty, iMeVperu, imesh
      integer ierr

      integer iDaxis,ntaxis
      integer nd1,nd2,nd3,nd4,nd5,nd6,nd7,nd8,nd9
      real(8) db2(nd2+1), db3(nd3+1), db4(nd4+1),
     &        db5(nd5+1), db6(nd6+1), db7(nd7+1)
      real(8) dw2(nd2), dw3(nd3),   dw4(nd4),   dw7(nd7)
      real(8) db3tmp(nd3+1)
      real(8) vl(nd5,nd6,nd7)
      real(8) anatalrst(nd1,nd2+1,nd3+1,nd4+1,nd5*nd6*nd7,nd8+1,nd9) ! S.H. 2021.9.1
      integer nij(ntaxis),ibin(ntaxis),nijaxs,nijaxs2,nijaxs3,nfgmax
      character(*) :: cij(ntaxis),cijaxs,cijaxs2
      integer ijaxs,irst,ianataldata,nanataldata

      integer :: njaxs

      integer np, nrst
      integer itmp
      integer id1,id2,id3,id4,id5,id6,id7,id8,id9
      integer itmpdata,itmpxx,itmpyy,itmpzz

*-----------------------------------------------------------------------
      integer maxnvar, maxcvar
      parameter (maxnvar=9, maxcvar=10000)
      integer nvar,iduc(maxnvar), ivar
      real*8 cvalue(maxcvar,maxnvar)
      common /cnvar/ nvar,iduc
      common /ccvar/ cvalue

*-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! input : iDaxis  (axis type)
!       1 energy(energy1) axis
!       2 reg axis
!       3 x axis
!       4 y axis
!       5 z axis (xyz)
!       6 r axis
!       7 z axis (rz)
!       8 tet axis
!       9 time axis
!      10 angle_r axis  (12,13)
!      11 angle_p axis
!      12 cos (cos theta)
!      13 the (theta deg)
!      14 energy2 axis
!      15 mass axis
!      16 reg axis of t-yield
!      17 charge axis
!      18 chart  axis
!      19 kind axis
!      20 let axis
!      21 sed axis
!      22 act axis
!      23 wwg axis
!      31 xy axis     (matrix)
!      32 yz axis     (matrix)
!      33 xz axis     (matrix)
!      34 rz axis     (matrix)
!      35 t-eng axis  (matrix)
!      36 eng-t axis  (matrix)
!      37 t-e1 axis   (matrix)
!      38 e1-t axis   (matrix)
!      39 t-e2 axis   (matrix)
!      40 e2-t axis   (matrix)
!      41 e12 axis    (matrix)
!      42 e21 axis    (matrix)
!      51 c1 axis ! S.H. 2021.8.15
!-----------------------------------------------------------------------
!   dimension
!     nd1:(np): particle / level(t-yield)
!     nd2:(ne): energy / mz(t-yield)
!     nd3:      angle_p(t-cross) / energy2 / mn(t-yield)
!     nd4:(nt): time
!     nd5:(nx): x(xyz) / r(rz) / reg / tetpla
!     nd6:(ny): y(xyz) / z(rz)
!     nd7:(nz): z(xyz) / angle_r(t-track_rz, t-adjoint_rz)
!     nd8:(nm): Multiplier/kind/action
!     nd9:      nrst
!-----------------------------------------------------------------------
!
      integer icf,ix,iy,iz
      integer iat, iad, iaf
      icf(ix,iy,iz) = ix + ( iy - 1 ) * nd5 + ( iz - 1 ) * nd5 * nd6
      iat(iad,iaf) = iad + (iaf-1) * 2

      character yen*1
      yen  = char(92)

!-----------------------------------------------------------------------
C default setting of nij, ibin, cij

      nij (1) = nd1 ! np
      ibin(1) = 0
      nij (2) = nd2 ! ne
      if ( nd2 .eq. 0 ) nij (2) = 1 ! ne of t-deposit with output=dose
      ibin(2) = 1
      nij (3) = nd3 ! angle_p / energy2
      if ( nd3 .eq. 0 ) nij (3) = 1
      ibin(3) = 1
      nij (4) = nd4 ! nt
      if ( nd4 .eq. 0 ) nij (4) = 1
      ibin(4) = 1
      if ( imesh .eq. 1 ) then ! mesh = reg
       nij (5) = nd5 ! reg
       ibin(5) = 0
       nij (6) = 1
       ibin(6) = 0
       nij (7) = 1
       ibin(7) = 0
      else if ( imesh .eq. 2 ) then ! mesh = r-z
       nij (5) = nd5 ! r
       ibin(5) = 1
       nij (6) = nd6 ! z
       ibin(6) = 1
       nij (7) = nd7 ! angle_r
       ibin(7) = 1
      else if ( imesh .eq. 3 ) then ! mesh = xyz
       nij (5) = nd5 ! nx
       ibin(5) = 1
       nij (6) = nd6 ! ny
       ibin(6) = 1
       nij (7) = nd7 ! nz
       ibin(7) = 1
      else if ( imesh .eq. 4 ) then ! mesh = tet
       nij (5) = nd5 ! tet
       ibin(5) = 0
       nij (6) = 1
       ibin(6) = 0
       nij (7) = 1
       ibin(7) = 0
      else if ( imesh.eq.5 .or. imesh.eq.6 ) then ! point or ring detector of [t-point]
       nij (5) = nd5 ! number of point or ring detector
       ibin(5) = 0
       nij (6) = 1
       ibin(6) = 0
       nij (7) = 1
       ibin(7) = 0
      end if
      nij (8) = nd8 ! nm / kind/ action
      if ( nd8 .eq. 0 ) nij (8) = 1
      ibin(8) = 0
      nij (9) = nd9 ! nrst
      ibin(9) = 0

      cij(1) = 'Particle'
      if( iMeVperu.eq.0 ) then
         cij(2) = 'Energy [MeV]'
      else
         cij(2) = 'Energy [MeV/n]'
      end if
      if ( nd2 .eq. 0 ) cij(2) = 'F' ! ne of t-deposit with output=dose
      cij(3) = 'F'
      cij(4) = 'Time [nsec]'
      if ( nd4 .eq. 0 ) cij(4) = 'F' ! case that nt is not defined
      if ( imesh .eq. 1 ) then ! mesh = reg
       cij(5) = 'Region'
       cij(6) = 'F'
       cij(7) = 'F'
      else if ( imesh .eq. 2 ) then ! mesh = r-z
       cij(5) = 'r [cm]'
       cij(6) = 'z [cm]'
       cij(7) = 'F'
      else if ( imesh .eq. 3 ) then ! mesh = xyz
       cij(5) = 'x [cm]'
       cij(6) = 'y [cm]'
       cij(7) = 'z [cm]'
      else if ( imesh .eq. 4 ) then ! mesh = tet
       cij(5) = 'Tetla'
       cij(6) = 'F'
       cij(7) = 'F'
      else if ( imesh .eq. 5 ) then ! point detector of [t-point]
       cij(5) = 'Point'
       cij(6) = 'F'
       cij(7) = 'F'
      else if ( imesh .eq. 6 ) then ! ring detector of [t-point]
       cij(5) = 'Ring'
       cij(6) = 'F'
       cij(7) = 'F'
      end if
      cij(8) = 'Multiplier set'
      if ( nd8 .eq. 0 ) cij(8) = 'F' ! tallies without using multiplier
      cij(9) = 'nrst'
      ivar = 1 ! temporary specification in the case of manatally=2

!-----------------------------------------------------------------------

C make data set of axis blank
      nijaxs2 = 1
                                     ! [t-track]
      if ( iDaxis.eq.1 ) then        ! energy(energy1) axis
         nij(2) = 1
         nijaxs = nd2
         cijaxs = cij(2)
         cij(2) = 'F'
      else if ( iDaxis.eq.9 ) then   ! time axis
         nij(4) = 1
         nijaxs = nd4
         cijaxs = cij(4)
         cij(4) = 'F'
      else if ( iDaxis.eq.3 ) then   ! x axis
         nij(5) = 1
         nijaxs = nd5
         cijaxs = cij(5)
         cij(5) = 'F'
      else if ( iDaxis.eq.4 ) then   ! y axis
         nij(6) = 1
         nijaxs = nd6
         cijaxs = cij(6)
         cij(6) = 'F'
      else if ( iDaxis.eq.5 ) then   ! z axis
         nij(7) = 1
         nijaxs = nd7
         cijaxs = cij(7)
         cij(7) = 'F'
      else if ( iDaxis.eq.6 ) then   ! r axis
         nij(5) = 1
         nijaxs = nd5
         cijaxs = 'r [cm]'
         cij(5) = 'F'
      else if ( iDaxis.eq.7 ) then   ! z (rz) axis
         nij(6) = 1
         nijaxs = nd6
         cijaxs = 'z [cm]'
         cij(6) = 'F'
      else if ( iDaxis.eq.10 ) then   ! angle_r axis
         nij(7) = 1
         nijaxs = nd7
         if ( itaty .gt. 0.0 ) then
            cijaxs = 'Angle [radian]'
         else
            cijaxs = 'Angle [degree]'
         end if
         cij(7) = 'F'
      else if ( iDaxis.eq.12 .or. iDaxis.eq.13 ) then   ! cos or the of angle_p axis
         nij(3) = 1
         nijaxs = nd3
         if ( iDaxis.eq.12 ) then
            cijaxs = 'cos(\theta)'
            if ( itaty .lt. 0.0 ) then ! axis is cos but input format is theta
               do id3 = 1, nd3+ibin(3)
                  db3tmp(id3) = dcos(db3(id3) / 180.d0 * physc(1) )
               end do
               do id3 = 1, nd3+ibin(3)
                  db3(id3) = db3tmp(nd3+ibin(3)-id3+1)
               end do
            end if
         else
            cijaxs = '\theta\,[degree]'
            if ( itaty .gt. 0.0 ) then ! axis is the but input format is cos
               do id3 = 1, nd3+ibin(3)
                  db3tmp(id3) = dacos(db3(id3)) * 180.d0 / physc(1)
               end do
               do id3 = 1, nd3+ibin(3)
                  db3(id3) = db3tmp(nd3+ibin(3)-id3+1)
               end do
            end if
         end if
         cij(3) = 'F'
      else if ( iDaxis.eq.2 ) then    ! reg axis
         nij(5) = 1
         nijaxs = nd5
         cijaxs = 'Serial Num. of Region'
         cij(5) = 'F'
      else if ( iDaxis.eq.8 ) then    ! tet axis
         nij(5) = 1
         nijaxs = nd5
         cijaxs = 'Serial Num. of tetrahedron'
         cij(5) = 'F'
      else if ( iDaxis.eq.15 ) then   ! mass axis
         nij(3) = 1
         nijaxs = nd3
         cijaxs = 'Mass A'
         cij(3) = 'F'
         cij(2) = 'Charge Z'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.16 ) then   ! reg axis of t-yield
         nij(5) = 1
         nijaxs = nd5
         cijaxs = 'Serial Num. of Region'
         cij(5) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.17 ) then   ! charge axis
         nij(2) = 1
         nijaxs = nd2
         cijaxs = 'Charge Z'
         cij(2) = 'F'
         nij(3) = nd3
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.18 ) then   ! chart axis     (matrix)
         nij(2)  = 1
         nij(3)  = 1
         nijaxs  = nd3
         nijaxs2 = nd2
         nijaxs3 = 0
         cijaxs  = 'N: neutron number'
         cijaxs2 = 'Z: proton number'
         cij(2) = 'F'
         cij(3) = 'F'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.19 ) then   ! r axis of t-yield
         nij(5) = 1
         nijaxs = nd5
         cijaxs = 'r [cm]'
         cij(5) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.20 ) then   ! z (rz) axis of t-yield
         nij(6) = 1
         nijaxs = nd6
         cijaxs = 'z [cm]'
         cij(6) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.21 ) then   ! x axis of t-yield
         nij(5) = 1
         nijaxs = nd5
         cijaxs = cij(5)
         cij(5) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.22 ) then   ! y axis of t-yield
         nij(6) = 1
         nijaxs = nd6
         cijaxs = cij(6)
         cij(6) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.23 ) then   ! z axis of t-yield
         nij(7) = 1
         nijaxs = nd7
         cijaxs = cij(7)
         cij(7) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.31 ) then   ! xy (x-y) axis     (matrix)
         nij(5)  = 1
         nij(6)  = 1
         nijaxs  = nd5
         nijaxs2 = nd6
         nijaxs3 = nd7
         cijaxs  = cij(5)
         cijaxs2 = cij(6)
         cij(5) = 'F'
         cij(6) = 'F'
      else if ( iDaxis.eq.32 ) then   ! yz (z-y) axis     (matrix)
         nij(6)  = 1
         nij(7)  = 1
         nijaxs  = nd7
         nijaxs2 = nd6
         nijaxs3 = nd5
         cijaxs  = cij(7)
         cijaxs2 = cij(6)
         cij(6) = 'F'
         cij(7) = 'F'
      else if ( iDaxis.eq.33 ) then   ! xz (z-x) axis     (matrix)
         nij(5)  = 1
         nij(7)  = 1
         nijaxs  = nd7
         nijaxs2 = nd5
         nijaxs3 = nd6
         cijaxs  = cij(7)
         cijaxs2 = cij(5)
         cij(5) = 'F'
         cij(7) = 'F'
      else if ( iDaxis.eq.34 ) then   ! rz axis (z-r)     (matrix)
         nij(5)  = 1
         nij(6)  = 1
         nijaxs  = nd6
         nijaxs2 = nd5
         nijaxs3 = 1
         cijaxs  = cij(6)
         cijaxs2 = cij(5)
         cij(5) = 'F'
         cij(6) = 'F'
      else if ( iDaxis.eq.35 ) then   ! rz axis (z-r) of t-yield (matrix)
         nij(5)  = 1
         nij(6)  = 1
         nijaxs  = nd6
         nijaxs2 = nd5
         nijaxs3 = 1
         cijaxs  = cij(6)
         cijaxs2 = cij(5)
         cij(5) = 'F'
         cij(6) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.36 ) then   ! xy (x-y) axis of t-yield     (matrix)
         nij(5)  = 1
         nij(6)  = 1
         nijaxs  = nd5
         nijaxs2 = nd6
         nijaxs3 = nd7
         cijaxs  = cij(5)
         cijaxs2 = cij(6)
         cij(5) = 'F'
         cij(6) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.37 ) then   ! yz (z-y) axis of t-yield     (matrix)
         nij(6)  = 1
         nij(7)  = 1
         nijaxs  = nd7
         nijaxs2 = nd6
         nijaxs3 = nd5
         cijaxs  = cij(7)
         cijaxs2 = cij(6)
         cij(6) = 'F'
         cij(7) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.38 ) then   ! xz (z-x) axis of t-yield     (matrix)
         nij(5)  = 1
         nij(7)  = 1
         nijaxs  = nd7
         nijaxs2 = nd5
         nijaxs3 = nd6
         cijaxs  = cij(7)
         cijaxs2 = cij(5)
         cij(5) = 'F'
         cij(7) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.51 ) then   ! c1 axis ! S.H. 2021.8.15
         nij(9) = 2
         nijaxs = (nd9-3)/2
         cijaxs = cij(9)
         if ( iduc(ivar) .lt. 10 ) then
            write(cijaxs,'(a1,i1,a6)') 'c',iduc(ivar),'-value'
         else
            write(cijaxs,'(a1,i2,a6)') 'c',iduc(ivar),'-value'
         end if
         cij(9) = 'F'
      end if


C count the number of all meshes
      nfgmax = 0
      do itmp = 2, ntaxis-1
         if ( cij(itmp) .ne. 'F' ) then
            nfgmax = nfgmax + nij(itmp)+ibin(itmp)
         end if
      end do
      allocate( fg(nfgmax) )

C set fg data
      itmp = 0
      if ( cij(2) .ne. 'F' ) then ! energy
         do id2 = 1, nij(2)+ibin(2)
            itmp = itmp + 1
            fg(itmp) = db2(id2)
         end do
      end if
      if ( cij(3) .ne. 'F' ) then ! angle_p / energy2
         do id3 = 1, nij(3)+ibin(3)
            itmp = itmp + 1
            fg(itmp) = db3(id3)
         end do
      end if
      if ( cij(4) .ne. 'F' ) then ! time
         do id4 = 1, nij(4)+ibin(4)
            itmp = itmp + 1
            fg(itmp) = db4(id4)
         end do
      end if
      if ( cij(5) .ne. 'F' ) then ! x,r(z),reg,tetla
         do id5 = 1, nij(5)+ibin(5)
            itmp = itmp + 1
            fg(itmp) = db5(id5)
         end do
      end if
      if ( cij(6) .ne. 'F' ) then ! y,(r)z
         do id6 = 1, nij(6)+ibin(6)
            itmp = itmp + 1
            fg(itmp) = db6(id6)
         end do
      end if
      if (cij(7) .ne. 'F' ) then ! z, angle_r
         do id7 = 1, nij(7)+ibin(7)
            itmp = itmp + 1
            fg(itmp) = db7(id7)
         end do
      end if
      if ( cij(8) .ne. 'F' ) then ! Multiplier/kind/action
         do id8 = 1, nij(8)+ibin(8)
            itmp = itmp + 1
            fg(itmp) = id8
         end do
      end if
      if ( itmp .ne. nfgmax ) write(*,*) 'error itmp .ne. nfgmax' ! for debug

C set fgaxs data
      allocate( fgaxs(nijaxs+1) )
      do itmp = 1, nijaxs+1
         if ( iDaxis.eq.1  .or.         ! axis = eng (e1)
     &        iDaxis.eq.41 ) then       ! axis = e12   (matrix)
            fgaxs(itmp) = db2(itmp)

         else if ( iDaxis.eq.9  .or.    ! axis = t
     &             iDaxis.eq.39 ) then  ! axis = t-e2  (matrix)
            fgaxs(itmp) = db4(itmp)

         else if ( iDaxis.eq.2 .or.     ! axis = reg
     &             iDaxis.eq.8 ) then   ! axis = tet
            if ( itmp .le. nd5 ) then
               fgaxs(itmp) = vl(itmp,1,1)
            else
               fgaxs(itmp) = 1d0
            end if

         else if ( iDaxis.eq.3 .or.     ! axis = x
     &             iDaxis.eq.6 .or.     ! axis = r
     &             iDaxis.eq.19 .or.    ! axis = r of t-yield
     &             iDaxis.eq.21 ) then  ! axis = x of t-yield
            fgaxs(itmp) = db5(itmp)

         else if ( iDaxis.eq.4 .or.     ! axis = y
     &             iDaxis.eq.7 .or.     ! axis = z (rz)
     &             iDaxis.eq.20 .or.    ! axis = z (rz) of t-yield
     &             iDaxis.eq.22 ) then  ! axis = y of t-yield
            fgaxs(itmp) = db6(itmp)

         else if ( iDaxis.eq.5 .or.     ! axis = z (xyz)
     &             iDaxis.eq.23  ) then ! axis = z (xyz) of t-yield
            fgaxs(itmp) = db7(itmp)

         else if ( iDaxis.eq.10 ) then  ! axis = angle_r
            fgaxs(itmp) = db7(itmp)

         else if ( iDaxis.eq.31 .or.    ! axis = xy    (matrix)
     &             iDaxis.eq.36 ) then  ! axis = xy of t-yield    (matrix)
            fgaxs(itmp) = db5(itmp)

         else if ( iDaxis.eq.32 .or.    ! axis = yz    (matrix)
     &             iDaxis.eq.37 ) then  ! axis = yz of t-yield    (matrix)
            fgaxs(itmp) = db7(itmp)

         else if ( iDaxis.eq.33 .or.    ! axis = xz    (matrix)
     &             iDaxis.eq.38 ) then  ! axis = xz of t-yield    (matrix)
            fgaxs(itmp) = db7(itmp)

         else if ( iDaxis.eq.34 .or.    ! axis = rz    (matrix)
     &             iDaxis.eq.35 ) then  ! axis = rz of t-yield
            fgaxs(itmp) = db6(itmp)

         else if ( iDaxis.eq.12 .or.    ! axis = cos (cos theta)
     &             iDaxis.eq.13 ) then  ! axis = the (theta deg)
            fgaxs(itmp) = db3(itmp)

         else if ( iDaxis.eq.14 .or.    ! axis = e2
     &             iDaxis.eq.40 .or.    ! axis = e2-t  (matrix)
     &             iDaxis.eq.42 ) then  ! axis = e21   (matrix)
            fgaxs(itmp) = db3(itmp)

         else if ( iDaxis.eq.15 ) then  ! axis = mass
            fgaxs(itmp) = db3(itmp)

         else if ( iDaxis.eq.16 ) then  ! axis = reg in t-yield
            if ( itmp .le. nd5 ) then
               fgaxs(itmp) = vl(itmp,1,1)
            else
               fgaxs(itmp) = 1d0
            end if

         else if ( iDaxis.eq.17 ) then  ! axis = charge
            fgaxs(itmp) = db2(itmp)

         else if ( iDaxis.eq.18 ) then  ! axis = chart
            fgaxs(itmp) = db3(itmp)

         else if ( iDaxis.eq.51 ) then  !  c1 axis ! S.H. 2021.8.15
            fgaxs(itmp) = cvalue(itmp,ivar)
         end if
      end do

      if ( iDaxis.ge.31 .and. iDaxis.le.50 ) then  ! axis is matrix ! S.H. 2021.8.15
         allocate( fgaxs2(nijaxs2+1) )
         if ( nijaxs3 .gt. 0 ) then
            allocate( fgaxs3(nijaxs3+1) )
         else
         end if

         if ( iDaxis.eq.31 .or.          ! axis = xy    (matrix)
     &        iDaxis.eq.36 ) then        ! axis = xy of t-yield    (matrix)
            do itmp = 1, nijaxs2+1
               fgaxs2(itmp) = db6(itmp)
            end do
            do itmp = 1, nijaxs3+1
               fgaxs3(itmp) = db7(itmp)
            end do

         else if ( iDaxis.eq.32 .or.     ! axis = yz    (matrix)
     &             iDaxis.eq.37 ) then   ! axis = yz of t-yield    (matrix)
            do itmp = 1, nijaxs2+1
               fgaxs2(itmp) = db6(itmp)
            end do
            do itmp = 1, nijaxs3+1
               fgaxs3(itmp) = db5(itmp)
            end do

         else if ( iDaxis.eq.33 .or.     ! axis = xz    (matrix)
     &             iDaxis.eq.38 ) then   ! axis = xz of t-yield    (matrix)
            do itmp = 1, nijaxs2+1
               fgaxs2(itmp) = db5(itmp)
            end do
            do itmp = 1, nijaxs3+1
               fgaxs3(itmp) = db6(itmp)
            end do

         else if ( iDaxis.eq.34 .or.    ! axis = rz    (matrix)
     &             iDaxis.eq.35 ) then  ! axis = rz of t-yield (matrix)
            do itmp = 1, nijaxs2+1
               fgaxs2(itmp) = db5(itmp)
            end do

         else if ( iDaxis.eq.18 ) then  ! axis = chart (matrix)
            do itmp = 1, nijaxs2
               fgaxs2(itmp) = db2(itmp)
            end do



         else if ( iDaxis.eq.39 .or.    ! axis = t-e2  (matrix)
     &             iDaxis.eq.41 ) then  ! axis = e12   (matrix)
            do itmp = 1, nijaxs2+1
               fgaxs2(itmp) = db3(itmp)
            end do

         end if
      end if

!-----------------------------------------------------------------------

C calculate data size of rearranged data
      nanataldata = nij(2)*nij(3)*nij(4)*nij(5)*nij(6)*nij(7)*nij(8)
      np  =nij(1)
      nrst=nij(9)
!      allocate( anataldata(np,nanataldata,nijaxs*nijaxs2,nrst) )
!      allocate( delvol(nanataldata,nijaxs*nijaxs2) )
      allocate( anataldata(np,nanataldata,nijaxs*nijaxs2+1,nrst) )
      allocate( delvol(nanataldata,nijaxs*nijaxs2+1) )

C set rearranged data

      njaxs = nijaxs*nijaxs2 +1

      ianataldata = 0
      do id2 = 1, nij(2)    ! ne
      do id3 = 1, nij(3)    ! angle_p / energy2
      do id4 = 1, nij(4)    ! nt
      do id5 = 1, nij(5)    ! nx / r / reg / tetla
      do id6 = 1, nij(6)    ! ny / z
      do id7 = 1, nij(7)    ! nz / angle_r
      do id8 = 1, nij(8)    ! nm / kind/ action
       ianataldata = ianataldata + 1

       delvol(ianataldata, njaxs) = 0.0d0

       do ijaxs = 1, nijaxs*nijaxs2
       do id1 = 1, nij(1)
       do id9 = 1, nij(9)

        if ( iDaxis.eq.1 ) then                               ! axis = eng (e1)
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,ijaxs,id3,id4,icf(id5,id6,id7),id8,id9)
         delvol(ianataldata,ijaxs) = dw2(ijaxs)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + dw2(ijaxs)


        else if ( iDaxis.eq.14 ) then                         ! axis = eng2(e2)
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,ijaxs,id4,icf(id5,id6,id7),id8,id9)
         delvol(ianataldata,ijaxs) = dw3(ijaxs)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + dw3(ijaxs)

        else if ( iDaxis.eq.9 ) then                          ! axis = t
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,ijaxs,icf(id5,id6,id7),id8,id9)
         delvol(ianataldata,ijaxs) = dw4(ijaxs)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + dw4(ijaxs)

        else if ( iDaxis.eq.2 .or.                            ! axis = reg
     &            iDaxis.eq.3 .or.                            ! axis = x
     &            iDaxis.eq.6 .or.                            ! axis = r (rz)
     &            iDaxis.eq.8 .or.                            ! axis = tet
     &            iDaxis.eq.19 .or.                           ! axis = r (rz) of t-yield
     &            iDaxis.eq.21 ) then                         ! axis = x of t-yield
         anataldata(id1,ianataldata,ijaxs,id9)
     &          = anatalrst(id1,id2,id3,id4,icf(ijaxs,id6,id7),id8,id9)
          delvol(ianataldata,ijaxs) = vl(ijaxs,id6,id7)
          delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                               + vl(ijaxs,id6,id7)

        else if ( iDaxis.eq.4 .or.                            ! axis = y
     &            iDaxis.eq.7 .or.                            ! axis = z (rz)
     &            iDaxis.eq.20 .or.                           ! axis = z (rz) of t-yield
     &            iDaxis.eq.22 ) then                         ! axis = y of t-yield
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(id5,ijaxs,id7),id8,id9)
          delvol(ianataldata,ijaxs) = vl(id5,ijaxs,id7)
          delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                               + vl(id5,ijaxs,id7)

        else if ( iDaxis.eq.5 .or.                            ! axis = z (xyz)
     &            iDaxis.eq.23 ) then                         ! axis = z (xyz) of t-yield
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(id5,id6,ijaxs),id8,id9)
         delvol(ianataldata,ijaxs) = vl(id5,id6,ijaxs)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(id5,id6,ijaxs)

        else if ( iDaxis.eq.10 ) then                         ! axis = angle_r
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(id5,id6,ijaxs),id8,id9)
         delvol(ianataldata,ijaxs) = dw7(ijaxs)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + dw7(ijaxs)

        else if ( iDaxis.eq.11 ) then                         ! axis = angle_p

        else if ( iDaxis.eq.12 .or.                           ! axis = cos (cos theta)
     &            iDaxis.eq.13 ) then                         ! axis = the (theta deg)
         if ( ( iDaxis.eq.12 .and. itaty.gt.0 )
     &          .or. ( iDaxis.eq.13 .and. itaty.lt.0 ) ) then
          anataldata(id1,ianataldata,ijaxs,id9)
     &    = anatalrst(id1,id2,ijaxs,id4,icf(id5,id6,id7),id8,id9)
          delvol(ianataldata,ijaxs) = dw3(ijaxs)
          delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                               + dw3(ijaxs)

         else
          anataldata(id1,ianataldata,ijaxs,id9)
     &   =anatalrst(id1,id2,nijaxs-ijaxs+1,id4,icf(id5,id6,id7),id8,id9)
          delvol(ianataldata,ijaxs) = dw3(nijaxs-ijaxs+1)
          delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                               + dw3(nijaxs-ijaxs+1)

         end if

        else if ( iDaxis.eq.15 ) then                         ! axis = mass
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,ijaxs,id4,icf(id5,id6,id7),id8,id9)
         delvol(ianataldata,ijaxs) = vl(id5,id6,id7)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(id5,id6,id7)

        else if ( iDaxis.eq.16 ) then                         ! axis = reg in t-yield
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(ijaxs,id6,id7),id8,id9)
          delvol(ianataldata,ijaxs) = vl(ijaxs,id6,id7)
          delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                               + vl(ijaxs,id6,id7)

        else if ( iDaxis.eq.17 ) then                         ! axis = charge
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,ijaxs,id3,id4,icf(id5,id6,id7),id8,id9)
         delvol(ianataldata,ijaxs) = vl(id5,id6,id7)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(id5,id6,id7)

        else if ( iDaxis.eq.18 ) then                         ! axis = chart (matrix)
         itmpdata = int(ijaxs-1)
         itmpxx = mod(itmpdata,nd3)+1
         itmpdata = int((ijaxs-itmpxx)/nd3)
         itmpyy = mod(itmpdata,nd2)+1
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,itmpyy,itmpxx,id4,icf(id5,id6,id7),id8,id9)
         delvol(ianataldata,ijaxs) = vl(id5,id6,id7)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(id5,id6,id7)




        else if ( iDaxis.eq.31 .or.                           ! axis = xy    (matrix)
     &            iDaxis.eq.36 ) then                         ! axis = xy of t-yield    (matrix)
         itmpdata = int(ijaxs-1)
         itmpxx = mod(itmpdata,nd5)+1
         itmpdata = int((ijaxs-itmpxx)/nd5)
         itmpyy = mod(itmpdata,nd6)+1
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(itmpxx,itmpyy,id7),id8,id9)
         delvol(ianataldata,ijaxs) = vl(itmpxx,itmpyy,id7)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(itmpxx,itmpyy,id7)

        else if ( iDaxis.eq.32 .or.                           ! axis = yz    (matrix)
     &            iDaxis.eq.37 ) then                         ! axis = yz of t-yield    (matrix)
         itmpdata = int(ijaxs-1)
         itmpzz = mod(itmpdata,nd7)+1
         itmpdata = int((ijaxs-itmpzz)/nd7)
         itmpyy = mod(itmpdata,nd6)+1
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(id5,itmpyy,itmpzz),id8,id9)
         delvol(ianataldata,ijaxs) = vl(id5,itmpyy,itmpzz)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(id5,itmpyy,itmpzz)

        else if ( iDaxis.eq.33 .or.                           ! axis = xz    (matrix)
     &            iDaxis.eq.38 ) then                         ! axis = xz of t-yield    (matrix)
         itmpdata = int(ijaxs-1)
         itmpzz = mod(itmpdata,nd7)+1
         itmpdata = int((ijaxs-itmpzz)/nd7)
         itmpxx = mod(itmpdata,nd5)+1
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(itmpxx,id6,itmpzz),id8,id9)
         delvol(ianataldata,ijaxs) = vl(itmpxx,id6,itmpzz)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(itmpxx,id6,itmpzz)

        else if ( iDaxis.eq.34 .or.                           ! axis = rz    (matrix)
     &            iDaxis.eq.35 ) then                         ! axis = rz of t-yield (matrix)
         itmpdata = int(ijaxs-1)
         itmpyy = mod(itmpdata,nd6)+1
         itmpdata = int((ijaxs-itmpyy)/nd6)
         itmpxx = mod(itmpdata,nd5)+1
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(itmpxx,itmpyy,id7),id8,id9)
         delvol(ianataldata,ijaxs) = vl(itmpxx,itmpyy,id7)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(itmpxx,itmpyy,id7)




        else if ( iDaxis.eq.51 ) then ! c1 axis ! S.H. 2021.8.15
         anataldata(id1,ianataldata,ijaxs,id9)
     &   =anatalrst(id1,id2,id3,id4,icf(id5,id6,id7),id8,iat(id9,ijaxs))
         delvol(ianataldata,ijaxs) = 1d0
         delvol(ianataldata, njaxs) =  1.0d0

        end if

       end do   ! id9
       end do   ! id1
       end do   ! ijaxs
      end do    ! id8
      end do    ! id7
      end do    ! id6
      end do    ! id5
      end do    ! id4
      end do    ! id3
      end do    ! id2

!-----------------------------------------------------------------------
      return
      end subroutine anatal_rearrange

!***********************************************************************

!***********************************************************************
!                                                                      *
      subroutine anatal_rearrange_sum(itaty,iMeVperu,imesh,
     &     iDaxis,ntaxis,anatalrst,
     &     nd1,    nd2,    nd3,    nd4,    nd5,nd6i,nd7i,   nd8,nd9,
     &             db2,dw2,db3,dw3,db4,dw4,db5,db6,db7,dw7,
     &     vl,
     &     nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &     cij,cijaxs,cijaxs2,nanataldata,ierr)
!                                                                      *
!     created by S.Hashimoto on 2021/4/27                              *
!     modified by T.Miura on 2021/8/6                                   *
!                                                                      *
!     anataldata: rearranged data of anatalrst                         *
!     fg: rearranged mesh data of all meshes                           *
!     fgaxs: mesh data of specified axis                               *
!     delvol: delta volume, increments used in calculating sum over    *
!                                                                      *
!***********************************************************************

      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol

      implicit none
!-----------------------------------------------------------------------
      include 'err.inc'
      include 'param-physcnst.inc'
!-----------------------------------------------------------------------

      integer itaty, iMeVperu, imesh
      integer ierr

      integer iDaxis,ntaxis
      integer nd6i,nd7i
      integer nd1,nd2,nd3,nd4,nd5,nd6,nd7,nd8,nd9
      real(8) db2(nd2+1), db3(nd3+1), db4(nd4+1),
     &        db5(nd5+1), db6(nd6i+1), db7(nd7i+1)
      real(8) dw2(nd2), dw3(nd3),   dw4(nd4),   dw7(max(nd7i,1))
      real(8) db3tmp(nd3+1)
      real(8) vl(nd5,max(nd6i,1),max(nd7i,1))
!      real(8) anatalrst(nd1,nd2+1,nd3+1,nd4+1,nd5*nd6*nd7,nd8+1,nd9) ! S.H. 2021.9.1
      real(8) anatalrst(nd1,nd2+1,nd3+1,nd4+1,(nd5+1)*(nd6i+1)*(nd7i+1),
     &                  nd8+1,nd9) ! S.H. 2021.9.1
      integer nij(ntaxis),ibin(ntaxis),nijaxs,nijaxs2,nijaxs3,nfgmax
      character(*) :: cij(ntaxis),cijaxs,cijaxs2
      integer ijaxs,irst,ianataldata,nanataldata

      integer :: mjaxs,njaxs

      integer np, nrst
      integer itmp
      integer id1,id2,id3,id4,id5,id6,id7,id8,id9
      integer itmpdata,itmpxx,itmpyy,itmpzz

*-----------------------------------------------------------------------
      integer maxnvar, maxcvar
      parameter (maxnvar=9, maxcvar=10000)
      integer nvar,iduc(maxnvar), ivar
      real*8 cvalue(maxcvar,maxnvar)
      common /cnvar/ nvar,iduc
      common /ccvar/ cvalue

*-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! input : iDaxis  (axis type)
!       1 energy(energy1) axis
!       2 reg axis
!       3 x axis
!       4 y axis
!       5 z axis (xyz)
!       6 r axis
!       7 z axis (rz)
!       8 tet axis
!       9 time axis
!      10 angle_r axis  (12,13)
!      11 angle_p axis
!      12 cos (cos theta)
!      13 the (theta deg)
!      14 energy2 axis
!      15 mass axis
!      16 reg axis of t-yield
!      17 charge axis
!      18 chart  axis
!      19 kind axis
!      20 let axis
!      21 sed axis
!      22 act axis
!      23 wwg axis
!      31 xy axis     (matrix)
!      32 yz axis     (matrix)
!      33 xz axis     (matrix)
!      34 rz axis     (matrix)
!      35 t-eng axis  (matrix)
!      36 eng-t axis  (matrix)
!      37 t-e1 axis   (matrix)
!      38 e1-t axis   (matrix)
!      39 t-e2 axis   (matrix)
!      40 e2-t axis   (matrix)
!      41 e12 axis    (matrix)
!      42 e21 axis    (matrix)
!      51 c1 axis ! S.H. 2021.8.15
!-----------------------------------------------------------------------
!   dimension
!     nd1:(np): particle / level(t-yield)
!     nd2:(ne): energy / mz(t-yield)
!     nd3:      angle_p(t-cross) / energy2 / mn(t-yield)
!     nd4:(nt): time
!     nd5:(nx): x(xyz) / r(rz) / reg / tetpla
!     nd6:(ny): y(xyz) / z(rz)
!     nd7:(nz): z(xyz) / angle_r(t-track_rz, t-adjoint_rz)
!     nd8:(nm): Multiplier/kind/action
!     nd9:      nrst
!-----------------------------------------------------------------------
!
      integer icf,ix,iy,iz
      integer iat, iad, iaf
!      icf(ix,iy,iz) = ix + ( iy - 1 ) * nd5 + ( iz - 1 ) * nd5 * nd6
      icf(ix,iy,iz) = ix + ( iy - 1 ) * (nd5+1)
     &              + ( iz - 1 ) * (nd5+1) * (nd6i+1)
      iat(iad,iaf) = iad + (iaf-1) * 2

      character yen*1
      yen  = char(92)

      nd6 = max(nd6i,1)
      nd7 = max(nd7i,1)

!-----------------------------------------------------------------------
C default setting of nij, ibin, cij

      nij (1) = nd1 ! np
      ibin(1) = 0
      nij (2) = nd2 ! ne
      if ( nd2 .eq. 0 ) nij (2) = 1 ! ne of t-deposit with output=dose
      ibin(2) = 1
      nij (3) = nd3 ! angle_p / energy2
      if ( nd3 .eq. 0 ) nij (3) = 1
      ibin(3) = 1
      nij (4) = nd4 ! nt
      if ( nd4 .eq. 0 ) nij (4) = 1
      ibin(4) = 1
      if ( imesh .eq. 1 ) then ! mesh = reg
       nij (5) = nd5 ! reg
       ibin(5) = 0
       nij (6) = 1
       ibin(6) = 0
       nij (7) = 1
       ibin(7) = 0
      else if ( imesh .eq. 2 ) then ! mesh = r-z
       nij (5) = nd5 ! r
       ibin(5) = 1
       nij (6) = nd6 ! z
       ibin(6) = 1
       nij (7) = nd7 ! angle_r
       ibin(7) = 1
      else if ( imesh .eq. 3 ) then ! mesh = xyz
       nij (5) = nd5 ! nx
       ibin(5) = 1
       nij (6) = nd6 ! ny
       ibin(6) = 1
       nij (7) = nd7 ! nz
       ibin(7) = 1
      else if ( imesh .eq. 4 ) then ! mesh = tet
       nij (5) = nd5 ! tet
       ibin(5) = 0
       nij (6) = 1
       ibin(6) = 0
       nij (7) = 1
       ibin(7) = 0
      else if ( imesh.eq.5 .or. imesh.eq.6 ) then ! point or ring detector of [t-point]
       nij (5) = nd5 ! number of point or ring detector
       ibin(5) = 0
       nij (6) = 1
       ibin(6) = 0
       nij (7) = 1
       ibin(7) = 0
      end if
      nij (8) = nd8 ! nm / kind/ action
      if ( nd8 .eq. 0 ) nij (8) = 1
      ibin(8) = 0
      nij (9) = nd9 ! nrst
      ibin(9) = 0

      cij(1) = 'Particle'
      if( iMeVperu.eq.0 ) then
         cij(2) = 'Energy [MeV]'
      else
         cij(2) = 'Energy [MeV/n]'
      end if
      if ( nd2 .eq. 0 ) cij(2) = 'F' ! ne of t-deposit with output=dose
      cij(3) = 'F'
      cij(4) = 'Time [nsec]'
      if ( nd4 .eq. 0 ) cij(4) = 'F' ! case that nt is not defined
      if ( imesh .eq. 1 ) then ! mesh = reg
       cij(5) = 'Region'
       cij(6) = 'F'
       cij(7) = 'F'
      else if ( imesh .eq. 2 ) then ! mesh = r-z
       cij(5) = 'r [cm]'
       cij(6) = 'z [cm]'
       cij(7) = 'F'
      else if ( imesh .eq. 3 ) then ! mesh = xyz
       cij(5) = 'x [cm]'
       cij(6) = 'y [cm]'
       cij(7) = 'z [cm]'
      else if ( imesh .eq. 4 ) then ! mesh = tet
       cij(5) = 'Tetla'
       cij(6) = 'F'
       cij(7) = 'F'
      else if ( imesh .eq. 5 ) then ! point detector of [t-point]
       cij(5) = 'Point'
       cij(6) = 'F'
       cij(7) = 'F'
      else if ( imesh .eq. 6 ) then ! ring detector of [t-point]
       cij(5) = 'Ring'
       cij(6) = 'F'
       cij(7) = 'F'
      end if
      cij(8) = 'Multiplier set'
      if ( nd8 .eq. 0 ) cij(8) = 'F' ! tallies without using multiplier
      cij(9) = 'nrst'
      ivar = 1 ! temporary specification in the case of manatally=2

!-----------------------------------------------------------------------

C make data set of axis blank
      nijaxs2 = 1
                                     ! [t-track]
      if ( iDaxis.eq.1 ) then        ! energy(energy1) axis
         nij(2) = 1
         nijaxs = nd2
         cijaxs = cij(2)
         cij(2) = 'F'
      else if ( iDaxis.eq.9 ) then   ! time axis
         nij(4) = 1
         nijaxs = nd4
         cijaxs = cij(4)
         cij(4) = 'F'
      else if ( iDaxis.eq.3 ) then   ! x axis
         nij(5) = 1
         nijaxs = nd5
         cijaxs = cij(5)
         cij(5) = 'F'
      else if ( iDaxis.eq.4 ) then   ! y axis
         nij(6) = 1
         nijaxs = nd6
         cijaxs = cij(6)
         cij(6) = 'F'
      else if ( iDaxis.eq.5 ) then   ! z axis
         nij(7) = 1
         nijaxs = nd7
         cijaxs = cij(7)
         cij(7) = 'F'
      else if ( iDaxis.eq.6 ) then   ! r axis
         nij(5) = 1
         nijaxs = nd5
         cijaxs = 'r [cm]'
         cij(5) = 'F'
      else if ( iDaxis.eq.7 ) then   ! z (rz) axis
         nij(6) = 1
         nijaxs = nd6
         cijaxs = 'z [cm]'
         cij(6) = 'F'
      else if ( iDaxis.eq.10 ) then   ! angle_r axis
         nij(7) = 1
         nijaxs = nd7
         if ( itaty .gt. 0.0 ) then
            cijaxs = 'Angle [radian]'
         else
            cijaxs = 'Angle [degree]'
         end if
         cij(7) = 'F'
      else if ( iDaxis.eq.12 .or. iDaxis.eq.13 ) then   ! cos or the of angle_p axis
         nij(3) = 1
         nijaxs = nd3
         if ( iDaxis.eq.12 ) then
            cijaxs = 'cos(\theta)'
            if ( itaty .lt. 0.0 ) then ! axis is cos but input format is theta
               do id3 = 1, nd3+ibin(3)
                  db3tmp(id3) = dcos(db3(id3) / 180.d0 * physc(1) )
               end do
               do id3 = 1, nd3+ibin(3)
                  db3(id3) = db3tmp(nd3+ibin(3)-id3+1)
               end do
            end if
         else
            cijaxs = '\theta\,[degree]'
            if ( itaty .gt. 0.0 ) then ! axis is the but input format is cos
               do id3 = 1, nd3+ibin(3)
                  db3tmp(id3) = dacos(db3(id3)) * 180.d0 / physc(1)
               end do
               do id3 = 1, nd3+ibin(3)
                  db3(id3) = db3tmp(nd3+ibin(3)-id3+1)
               end do
            end if
         end if
         cij(3) = 'F'
      else if ( iDaxis.eq.2 ) then    ! reg axis
         nij(5) = 1
         nijaxs = nd5
         cijaxs = 'Serial Num. of Region'
         cij(5) = 'F'
      else if ( iDaxis.eq.8 ) then    ! tet axis
         nij(5) = 1
         nijaxs = nd5
         cijaxs = 'Serial Num. of tetrahedron'
         cij(5) = 'F'
      else if ( iDaxis.eq.15 ) then   ! mass axis
         nij(3) = 1
         nijaxs = nd3
         cijaxs = 'Mass A'
         cij(3) = 'F'
         cij(2) = 'Charge Z'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.16 ) then   ! reg axis of t-yield
         nij(5) = 1
         nijaxs = nd5
         cijaxs = 'Serial Num. of Region'
         cij(5) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.17 ) then   ! charge axis
         nij(2) = 1
         nijaxs = nd2
         cijaxs = 'Charge Z'
         cij(2) = 'F'
         nij(3) = nd3
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.18 ) then   ! chart axis     (matrix)
         nij(2)  = 1
         nij(3)  = 1
         nijaxs  = nd3
         nijaxs2 = nd2
         nijaxs3 = 0
         cijaxs  = 'N: neutron number'
         cijaxs2 = 'Z: proton number'
         cij(2) = 'F'
         cij(3) = 'F'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.19 ) then   ! r axis of t-yield
         nij(5) = 1
         nijaxs = nd5
         cijaxs = 'r [cm]'
         cij(5) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.20 ) then   ! z (rz) axis of t-yield
         nij(6) = 1
         nijaxs = nd6
         cijaxs = 'z [cm]'
         cij(6) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.21 ) then   ! x axis of t-yield
         nij(5) = 1
         nijaxs = nd5
         cijaxs = cij(5)
         cij(5) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.22 ) then   ! y axis of t-yield
         nij(6) = 1
         nijaxs = nd6
         cijaxs = cij(6)
         cij(6) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.23 ) then   ! z axis of t-yield
         nij(7) = 1
         nijaxs = nd7
         cijaxs = cij(7)
         cij(7) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.31 ) then   ! xy (x-y) axis     (matrix)
         nij(5)  = 1
         nij(6)  = 1
         nijaxs  = nd5
         nijaxs2 = nd6
         nijaxs3 = nd7
         cijaxs  = cij(5)
         cijaxs2 = cij(6)
         cij(5) = 'F'
         cij(6) = 'F'
      else if ( iDaxis.eq.32 ) then   ! yz (z-y) axis     (matrix)
         nij(6)  = 1
         nij(7)  = 1
         nijaxs  = nd7
         nijaxs2 = nd6
         nijaxs3 = nd5
         cijaxs  = cij(7)
         cijaxs2 = cij(6)
         cij(6) = 'F'
         cij(7) = 'F'
      else if ( iDaxis.eq.33 ) then   ! xz (z-x) axis     (matrix)
         nij(5)  = 1
         nij(7)  = 1
         nijaxs  = nd7
         nijaxs2 = nd5
         nijaxs3 = nd6
         cijaxs  = cij(7)
         cijaxs2 = cij(5)
         cij(5) = 'F'
         cij(7) = 'F'
      else if ( iDaxis.eq.34 ) then   ! rz axis (z-r)     (matrix)
         nij(5)  = 1
         nij(6)  = 1
         nijaxs  = nd6
         nijaxs2 = nd5
         nijaxs3 = 1
         cijaxs  = cij(6)
         cijaxs2 = cij(5)
         cij(5) = 'F'
         cij(6) = 'F'
      else if ( iDaxis.eq.35 ) then   ! rz axis (z-r) of t-yield (matrix)
         nij(5)  = 1
         nij(6)  = 1
         nijaxs  = nd6
         nijaxs2 = nd5
         nijaxs3 = 1
         cijaxs  = cij(6)
         cijaxs2 = cij(5)
         cij(5) = 'F'
         cij(6) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.36 ) then   ! xy (x-y) axis of t-yield     (matrix)
         nij(5)  = 1
         nij(6)  = 1
         nijaxs  = nd5
         nijaxs2 = nd6
         nijaxs3 = nd7
         cijaxs  = cij(5)
         cijaxs2 = cij(6)
         cij(5) = 'F'
         cij(6) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.37 ) then   ! yz (z-y) axis of t-yield     (matrix)
         nij(6)  = 1
         nij(7)  = 1
         nijaxs  = nd7
         nijaxs2 = nd6
         nijaxs3 = nd5
         cijaxs  = cij(7)
         cijaxs2 = cij(6)
         cij(6) = 'F'
         cij(7) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.38 ) then   ! xz (z-x) axis of t-yield     (matrix)
         nij(5)  = 1
         nij(7)  = 1
         nijaxs  = nd7
         nijaxs2 = nd5
         nijaxs3 = nd6
         cijaxs  = cij(7)
         cijaxs2 = cij(5)
         cij(5) = 'F'
         cij(7) = 'F'
         cij(2) = 'Charge Z'
         cij(3) = 'Mass A'
         ibin(2) = 0
         ibin(3) = 0
      else if ( iDaxis.eq.51 ) then   ! c1 axis ! S.H. 2021.8.15
         nij(9) = 2
         nijaxs = (nd9-3)/2
         cijaxs = cij(9)
         if ( iduc(ivar) .lt. 10 ) then
            write(cijaxs,'(a1,i1,a6)') 'c',iduc(ivar),'-value'
         else
            write(cijaxs,'(a1,i2,a6)') 'c',iduc(ivar),'-value'
         end if
         cij(9) = 'F'
      end if


C count the number of all meshes
      nfgmax = 0
      do itmp = 2, ntaxis-1
         if ( cij(itmp) .ne. 'F' ) then
            nfgmax = nfgmax + nij(itmp)+ibin(itmp)
         end if
      end do
      allocate( fg(nfgmax) )

C set fg data
      itmp = 0
      if ( cij(2) .ne. 'F' ) then ! energy
         do id2 = 1, nij(2)+ibin(2)
            itmp = itmp + 1
            fg(itmp) = db2(id2)
         end do
      end if
      if ( cij(3) .ne. 'F' ) then ! angle_p / energy2
         do id3 = 1, nij(3)+ibin(3)
            itmp = itmp + 1
            fg(itmp) = db3(id3)
         end do
      end if
      if ( cij(4) .ne. 'F' ) then ! time
         do id4 = 1, nij(4)+ibin(4)
            itmp = itmp + 1
            fg(itmp) = db4(id4)
         end do
      end if
      if ( cij(5) .ne. 'F' ) then ! x,r(z),reg,tetla
         do id5 = 1, nij(5)+ibin(5)
            itmp = itmp + 1
            fg(itmp) = db5(id5)
         end do
      end if
      if ( cij(6) .ne. 'F' ) then ! y,(r)z
         do id6 = 1, nij(6)+ibin(6)
            itmp = itmp + 1
            fg(itmp) = db6(id6)
         end do
      end if
      if (cij(7) .ne. 'F' ) then ! z, angle_r
         do id7 = 1, nij(7)+ibin(7)
            itmp = itmp + 1
            fg(itmp) = db7(id7)
         end do
      end if
      if ( cij(8) .ne. 'F' ) then ! Multiplier/kind/action
         do id8 = 1, nij(8)+ibin(8)
            itmp = itmp + 1
            fg(itmp) = id8
         end do
      end if
      if ( itmp .ne. nfgmax ) write(*,*) 'error itmp .ne. nfgmax' ! for debug

C set fgaxs data
      allocate( fgaxs(nijaxs+1) )
      do itmp = 1, nijaxs+1
         if ( iDaxis.eq.1  .or.         ! axis = eng (e1)
     &        iDaxis.eq.41 ) then       ! axis = e12   (matrix)
            fgaxs(itmp) = db2(itmp)

         else if ( iDaxis.eq.9  .or.    ! axis = t
     &             iDaxis.eq.39 ) then  ! axis = t-e2  (matrix)
            fgaxs(itmp) = db4(itmp)

         else if ( iDaxis.eq.2 .or.     ! axis = reg
     &             iDaxis.eq.8 ) then   ! axis = tet
            if ( itmp .le. nd5 ) then
               fgaxs(itmp) = vl(itmp,1,1)
            else
               fgaxs(itmp) = 1d0
            end if

         else if ( iDaxis.eq.3 .or.     ! axis = x
     &             iDaxis.eq.6 .or.     ! axis = r
     &             iDaxis.eq.19 .or.    ! axis = r of t-yield
     &             iDaxis.eq.21 ) then  ! axis = x of t-yield
            fgaxs(itmp) = db5(itmp)

         else if ( iDaxis.eq.4 .or.     ! axis = y
     &             iDaxis.eq.7 .or.     ! axis = z (rz)
     &             iDaxis.eq.20 .or.    ! axis = z (rz) of t-yield
     &             iDaxis.eq.22 ) then  ! axis = y of t-yield
            fgaxs(itmp) = db6(itmp)

         else if ( iDaxis.eq.5 .or.     ! axis = z (xyz)
     &             iDaxis.eq.23  ) then ! axis = z (xyz) of t-yield
            fgaxs(itmp) = db7(itmp)

         else if ( iDaxis.eq.10 ) then  ! axis = angle_r
            fgaxs(itmp) = db7(itmp)

         else if ( iDaxis.eq.31 .or.    ! axis = xy    (matrix)
     &             iDaxis.eq.36 ) then  ! axis = xy of t-yield    (matrix)
            fgaxs(itmp) = db5(itmp)

         else if ( iDaxis.eq.32 .or.    ! axis = yz    (matrix)
     &             iDaxis.eq.37 ) then  ! axis = yz of t-yield    (matrix)
            fgaxs(itmp) = db7(itmp)

         else if ( iDaxis.eq.33 .or.    ! axis = xz    (matrix)
     &             iDaxis.eq.38 ) then  ! axis = xz of t-yield    (matrix)
            fgaxs(itmp) = db7(itmp)

         else if ( iDaxis.eq.34 .or.    ! axis = rz    (matrix)
     &             iDaxis.eq.35 ) then  ! axis = rz of t-yield
            fgaxs(itmp) = db6(itmp)

         else if ( iDaxis.eq.12 .or.    ! axis = cos (cos theta)
     &             iDaxis.eq.13 ) then  ! axis = the (theta deg)
            fgaxs(itmp) = db3(itmp)

         else if ( iDaxis.eq.14 .or.    ! axis = e2
     &             iDaxis.eq.40 .or.    ! axis = e2-t  (matrix)
     &             iDaxis.eq.42 ) then  ! axis = e21   (matrix)
            fgaxs(itmp) = db3(itmp)

         else if ( iDaxis.eq.15 ) then  ! axis = mass
            fgaxs(itmp) = db3(itmp)

         else if ( iDaxis.eq.16 ) then  ! axis = reg in t-yield
            if ( itmp .le. nd5 ) then
               fgaxs(itmp) = vl(itmp,1,1)
            else
               fgaxs(itmp) = 1d0
            end if

         else if ( iDaxis.eq.17 ) then  ! axis = charge
            fgaxs(itmp) = db2(itmp)

         else if ( iDaxis.eq.18 ) then  ! axis = chart
            fgaxs(itmp) = db3(itmp)

         else if ( iDaxis.eq.51 ) then  !  c1 axis ! S.H. 2021.8.15
            fgaxs(itmp) = cvalue(itmp,ivar)
         end if
      end do

      if ( iDaxis.ge.31 .and. iDaxis.le.50 ) then  ! axis is matrix ! S.H. 2021.8.15
         allocate( fgaxs2(nijaxs2+1) )
         if ( nijaxs3 .gt. 0 ) then
            allocate( fgaxs3(nijaxs3+1) )
         else
         end if

         if ( iDaxis.eq.31 .or.          ! axis = xy    (matrix)
     &        iDaxis.eq.36 ) then        ! axis = xy of t-yield    (matrix)
            do itmp = 1, nijaxs2+1
               fgaxs2(itmp) = db6(itmp)
            end do
            do itmp = 1, nijaxs3+1
               fgaxs3(itmp) = db7(itmp)
            end do

         else if ( iDaxis.eq.32 .or.     ! axis = yz    (matrix)
     &             iDaxis.eq.37 ) then   ! axis = yz of t-yield    (matrix)
            do itmp = 1, nijaxs2+1
               fgaxs2(itmp) = db6(itmp)
            end do
            do itmp = 1, nijaxs3+1
               fgaxs3(itmp) = db5(itmp)
            end do

         else if ( iDaxis.eq.33 .or.     ! axis = xz    (matrix)
     &             iDaxis.eq.38 ) then   ! axis = xz of t-yield    (matrix)
            do itmp = 1, nijaxs2+1
               fgaxs2(itmp) = db5(itmp)
            end do
            do itmp = 1, nijaxs3+1
               fgaxs3(itmp) = db6(itmp)
            end do

         else if ( iDaxis.eq.34 .or.    ! axis = rz    (matrix)
     &             iDaxis.eq.35 ) then  ! axis = rz of t-yield (matrix)
            do itmp = 1, nijaxs2+1
               fgaxs2(itmp) = db5(itmp)
            end do

         else if ( iDaxis.eq.18 ) then  ! axis = chart (matrix)
            do itmp = 1, nijaxs2
               fgaxs2(itmp) = db2(itmp)
            end do



         else if ( iDaxis.eq.39 .or.    ! axis = t-e2  (matrix)
     &             iDaxis.eq.41 ) then  ! axis = e12   (matrix)
            do itmp = 1, nijaxs2+1
               fgaxs2(itmp) = db3(itmp)
            end do

         end if
      end if

!-----------------------------------------------------------------------

C calculate data size of rearranged data
      nanataldata = nij(2)*nij(3)*nij(4)*nij(5)*nij(6)*nij(7)*nij(8)
      np  =nij(1)
      nrst=nij(9)
!      allocate( anataldata(np,nanataldata,nijaxs*nijaxs2,nrst) )
!      allocate( delvol(nanataldata,nijaxs*nijaxs2) )
      allocate( anataldata(np,nanataldata,nijaxs*nijaxs2+1,nrst) )
      allocate( delvol(nanataldata,nijaxs*nijaxs2+1) )

C set rearranged data

      mjaxs = nijaxs*nijaxs2
      njaxs = nijaxs*nijaxs2 +1

      ianataldata = 0
      do id2 = 1, nij(2)    ! ne
      do id3 = 1, nij(3)    ! angle_p / energy2
      do id4 = 1, nij(4)    ! nt
      do id5 = 1, nij(5)    ! nx / r / reg / tetla
      do id6 = 1, nij(6)    ! ny / z
      do id7 = 1, nij(7)    ! nz / angle_r
      do id8 = 1, nij(8)    ! nm / kind/ action
       ianataldata = ianataldata + 1

       delvol(ianataldata, njaxs) = 0.0d0

       do ijaxs = 1, nijaxs*nijaxs2
       do id1 = 1, nij(1)
       do id9 = 1, nij(9)

        if ( iDaxis.eq.1 ) then                               ! axis = eng (e1)
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,ijaxs,id3,id4,icf(id5,id6,id7),id8,id9)
         if(ijaxs == mjaxs) then
           anataldata(id1,ianataldata,njaxs,id9)
     &     = anatalrst(id1,njaxs,id3,id4,icf(id5,id6,id7),id8,id9)
         endif
         delvol(ianataldata,ijaxs) = dw2(ijaxs)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + dw2(ijaxs)


        else if ( iDaxis.eq.14 ) then                         ! axis = eng2(e2)
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,ijaxs,id4,icf(id5,id6,id7),id8,id9)
         if(ijaxs == mjaxs) then
           anataldata(id1,ianataldata,njaxs,id9)
     &   = anatalrst(id1,id2,njaxs,id4,icf(id5,id6,id7),id8,id9)
         endif
         delvol(ianataldata,ijaxs) = dw3(ijaxs)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + dw3(ijaxs)

        else if ( iDaxis.eq.9 ) then                          ! axis = t
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,ijaxs,icf(id5,id6,id7),id8,id9)
         if(ijaxs == mjaxs) then
           anataldata(id1,ianataldata,njaxs,id9)
     &   = anatalrst(id1,id2,id3,njaxs,icf(id5,id6,id7),id8,id9)
         endif
         delvol(ianataldata,ijaxs) = dw4(ijaxs)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + dw4(ijaxs)

        else if ( iDaxis.eq.2 .or.                            ! axis = reg
     &            iDaxis.eq.3 .or.                            ! axis = x
     &            iDaxis.eq.6 .or.                            ! axis = r (rz)
     &            iDaxis.eq.8 .or.                            ! axis = tet
     &            iDaxis.eq.19 .or.                           ! axis = r (rz) of t-yield
     &            iDaxis.eq.21 ) then                         ! axis = x of t-yield
         anataldata(id1,ianataldata,ijaxs,id9)
     &          = anatalrst(id1,id2,id3,id4,icf(ijaxs,id6,id7),id8,id9)
         if(ijaxs == mjaxs) then
           anataldata(id1,ianataldata,njaxs,id9)
     &         = anatalrst(id1,id2,id3,id4,icf(njaxs,id6,id7),id8,id9)
         endif
          delvol(ianataldata,ijaxs) = vl(ijaxs,id6,id7)
          delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                               + vl(ijaxs,id6,id7)

        else if ( iDaxis.eq.4 .or.                            ! axis = y
     &            iDaxis.eq.7 .or.                            ! axis = z (rz)
     &            iDaxis.eq.20 .or.                           ! axis = z (rz) of t-yield
     &            iDaxis.eq.22 ) then                         ! axis = y of t-yield
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(id5,ijaxs,id7),id8,id9)
         if(ijaxs == mjaxs) then
           anataldata(id1,ianataldata,njaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(id5,njaxs,id7),id8,id9)
         endif
          delvol(ianataldata,ijaxs) = vl(id5,ijaxs,id7)
          delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                               + vl(id5,ijaxs,id7)

        else if ( iDaxis.eq.5 .or.                            ! axis = z (xyz)
     &            iDaxis.eq.23 ) then                         ! axis = z (xyz) of t-yield
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(id5,id6,ijaxs),id8,id9)
         if(ijaxs == mjaxs) then
           anataldata(id1,ianataldata,njaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(id5,id6,njaxs),id8,id9)
         endif
         delvol(ianataldata,ijaxs) = vl(id5,id6,ijaxs)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(id5,id6,ijaxs)

        else if ( iDaxis.eq.10 ) then                         ! axis = angle_r
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(id5,id6,ijaxs),id8,id9)
         if(ijaxs == mjaxs) then
           anataldata(id1,ianataldata,njaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(id5,id6,njaxs),id8,id9)
         endif
         delvol(ianataldata,ijaxs) = dw7(ijaxs)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + dw7(ijaxs)

        else if ( iDaxis.eq.11 ) then                         ! axis = angle_p

        else if ( iDaxis.eq.12 .or.                           ! axis = cos (cos theta)
     &            iDaxis.eq.13 ) then                         ! axis = the (theta deg)
         if ( ( iDaxis.eq.12 .and. itaty.gt.0 )
     &          .or. ( iDaxis.eq.13 .and. itaty.lt.0 ) ) then
          anataldata(id1,ianataldata,ijaxs,id9)
     &    = anatalrst(id1,id2,ijaxs,id4,icf(id5,id6,id7),id8,id9)
         if(ijaxs == mjaxs) then
          anataldata(id1,ianataldata,njaxs,id9)
     &    = anatalrst(id1,id2,njaxs,id4,icf(id5,id6,id7),id8,id9)
         endif
          delvol(ianataldata,ijaxs) = dw3(ijaxs)
          delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                               + dw3(ijaxs)

         else
          anataldata(id1,ianataldata,ijaxs,id9)
     &   =anatalrst(id1,id2,nijaxs-ijaxs+1,id4,icf(id5,id6,id7),id8,id9)
         if(ijaxs == mjaxs) then
           anataldata(id1,ianataldata,njaxs,id9)
     &   =anatalrst(id1,id2,nijaxs-njaxs+1,id4,icf(id5,id6,id7),id8,id9)
         endif
          delvol(ianataldata,ijaxs) = dw3(nijaxs-ijaxs+1)
          delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                               + dw3(nijaxs-ijaxs+1)

         end if

        else if ( iDaxis.eq.15 ) then                         ! axis = mass
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,ijaxs,id4,icf(id5,id6,id7),id8,id9)
         if(ijaxs == mjaxs) then
           anataldata(id1,ianataldata,njaxs,id9)
     &   = anatalrst(id1,id2,njaxs,id4,icf(id5,id6,id7),id8,id9)
         endif
         delvol(ianataldata,ijaxs) = vl(id5,id6,id7)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(id5,id6,id7)

        else if ( iDaxis.eq.16 ) then                         ! axis = reg in t-yield
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(ijaxs,id6,id7),id8,id9)
         if(ijaxs == mjaxs) then
           anataldata(id1,ianataldata,njaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(njaxs,id6,id7),id8,id9)
         endif
          delvol(ianataldata,ijaxs) = vl(ijaxs,id6,id7)
          delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                               + vl(ijaxs,id6,id7)

        else if ( iDaxis.eq.17 ) then                         ! axis = charge
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,ijaxs,id3,id4,icf(id5,id6,id7),id8,id9)
         if(ijaxs == mjaxs) then
           anataldata(id1,ianataldata,njaxs,id9)
     &   = anatalrst(id1,njaxs,id3,id4,icf(id5,id6,id7),id8,id9)
         endif
         delvol(ianataldata,ijaxs) = vl(id5,id6,id7)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(id5,id6,id7)

        else if ( iDaxis.eq.18 ) then                         ! axis = chart (matrix)
         itmpdata = int(ijaxs-1)
         itmpxx = mod(itmpdata,nd3)+1
         itmpdata = int((ijaxs-itmpxx)/nd3)
         itmpyy = mod(itmpdata,nd2)+1
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,itmpyy,itmpxx,id4,icf(id5,id6,id7),id8,id9)
!         if(ijaxs == mjaxs) then
!           anataldata(id1,ianataldata,njaxs,id9)
!     &   = anatalrst(id1,itmpyy,itmpxx,id4,icf(id5,id6,id7),id8,id9)
!         endif
         delvol(ianataldata,ijaxs) = vl(id5,id6,id7)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(id5,id6,id7)




        else if ( iDaxis.eq.31 .or.                           ! axis = xy    (matrix)
     &            iDaxis.eq.36 ) then                         ! axis = xy of t-yield    (matrix)
         itmpdata = int(ijaxs-1)
         itmpxx = mod(itmpdata,nd5)+1
         itmpdata = int((ijaxs-itmpxx)/nd5)
         itmpyy = mod(itmpdata,nd6)+1
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(itmpxx,itmpyy,id7),id8,id9)
!         if(ijaxs == mjaxs) then
!           anataldata(id1,ianataldata,njaxs,id9)
!     &   = anatalrst(id1,id2,id3,id4,icf(itmpxx,itmpyy,id7),id8,id9)
!         endif
         delvol(ianataldata,ijaxs) = vl(itmpxx,itmpyy,id7)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(itmpxx,itmpyy,id7)

        else if ( iDaxis.eq.32 .or.                           ! axis = yz    (matrix)
     &            iDaxis.eq.37 ) then                         ! axis = yz of t-yield    (matrix)
         itmpdata = int(ijaxs-1)
         itmpzz = mod(itmpdata,nd7)+1
         itmpdata = int((ijaxs-itmpzz)/nd7)
         itmpyy = mod(itmpdata,nd6)+1
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(id5,itmpyy,itmpzz),id8,id9)
!         if(ijaxs == mjaxs) then
!           anataldata(id1,ianataldata,njaxs,id9)
!     &   = anatalrst(id1,id2,id3,id4,icf(id5,itmpyy,itmpzz),id8,id9)
!         endif
         delvol(ianataldata,ijaxs) = vl(id5,itmpyy,itmpzz)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(id5,itmpyy,itmpzz)

        else if ( iDaxis.eq.33 .or.                           ! axis = xz    (matrix)
     &            iDaxis.eq.38 ) then                         ! axis = xz of t-yield    (matrix)
         itmpdata = int(ijaxs-1)
         itmpzz = mod(itmpdata,nd7)+1
         itmpdata = int((ijaxs-itmpzz)/nd7)
         itmpxx = mod(itmpdata,nd5)+1
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(itmpxx,id6,itmpzz),id8,id9)
!         if(ijaxs == mjaxs) then
!           anataldata(id1,ianataldata,njaxs,id9)
!     &   = anatalrst(id1,id2,id3,id4,icf(itmpxx,id6,itmpzz),id8,id9)
!         endif
         delvol(ianataldata,ijaxs) = vl(itmpxx,id6,itmpzz)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(itmpxx,id6,itmpzz)

        else if ( iDaxis.eq.34 .or.                           ! axis = rz    (matrix)
     &            iDaxis.eq.35 ) then                         ! axis = rz of t-yield (matrix)
         itmpdata = int(ijaxs-1)
         itmpyy = mod(itmpdata,nd6)+1
         itmpdata = int((ijaxs-itmpyy)/nd6)
         itmpxx = mod(itmpdata,nd5)+1
         anataldata(id1,ianataldata,ijaxs,id9)
     &   = anatalrst(id1,id2,id3,id4,icf(itmpxx,itmpyy,id7),id8,id9)
!         if(ijaxs == mjaxs) then
!           anataldata(id1,ianataldata,njaxs,id9)
!     &   = anatalrst(id1,id2,id3,id4,icf(itmpxx,itmpyy,id7),id8,id9)
!         endif
         delvol(ianataldata,ijaxs) = vl(itmpxx,itmpyy,id7)
         delvol(ianataldata, njaxs) =  delvol(ianataldata, njaxs)
     &                              + vl(itmpxx,itmpyy,id7)




        else if ( iDaxis.eq.51 ) then ! c1 axis ! S.H. 2021.8.15
         anataldata(id1,ianataldata,ijaxs,id9)
     &   =anatalrst(id1,id2,id3,id4,icf(id5,id6,id7),id8,iat(id9,ijaxs))
!         if(ijaxs == mjaxs) then
!           anataldata(id1,ianataldata,njaxs,id9)
!     &   =anatalrst(id1,id2,id3,id4,icf(id5,id6,id7),id8,iat(id9,njaxs))
!         endif
         delvol(ianataldata,ijaxs) = 1d0
         delvol(ianataldata, njaxs) =  1.0d0

        end if

       end do   ! id9
       end do   ! id1
       end do   ! ijaxs
      end do    ! id8
      end do    ! id7
      end do    ! id6
      end do    ! id5
      end do    ! id4
      end do    ! id3
      end do    ! id2

!-----------------------------------------------------------------------
      return
      end subroutine anatal_rearrange_sum

!***********************************************************************
