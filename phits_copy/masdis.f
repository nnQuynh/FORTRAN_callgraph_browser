************************************************************************
*                                                                      *
      subroutine pr_mdis(ida,icd)
*                                                                      *
*                                                                      *
*        Last Revised:     2000 04 25                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write isotope distribution                           *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              ida         : output unit                               *
*                                                                      *
*              icd         : = 30 ; NCASC                              *
*                              31 ; NEVAP                              *
*                              32 ; fission mass distribution          *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      include 'param01.inc'

      common /summas/ sumas(3,0:maxpt,0:maxnt)

*-----------------------------------------------------------------------

      character hzz(103)*2

      data hzz /
     & 2hH , 2hHe, 2hLi, 2hBe, 2hB , 2hC , 2hN , 2hO , 2hF , 2hNe,
     & 2hNa, 2hMg, 2hAl, 2hSi, 2hP , 2hS , 2hCl, 2hAr, 2hK , 2hCa,
     & 2hSc, 2hTi, 2hV , 2hCr, 2hMn, 2hFe, 2hCo, 2hNi, 2hCu, 2hZn,
     & 2hGa, 2hGe, 2hAs, 2hSe, 2hBr, 2hKr, 2hRb, 2hSr, 2hY , 2hZr,
     & 2hNb, 2hMo, 2hTc, 2hRu, 2hRh, 2hPd, 2hAg, 2hCd, 2hIn, 2hSn,
     & 2hSb, 2hTe, 2hI , 2hXe, 2hCs, 2hBa, 2hLa, 2hCe, 2hPr, 2hNd,
     & 2hPm, 2hSm, 2hEu, 2hGd, 2hTb, 2hDy, 2hHo, 2hEr, 2hTm, 2hYb,
     & 2hLu, 2hHf, 2hTa, 2hW , 2hRe, 2hOs, 2hIr, 2hPt, 2hAu, 2hHg,
     & 2hTl, 2hPb, 2hBi, 2hPo, 2hAt, 2hRn, 2hFr, 2hRa, 2hAc, 2hTh,
     & 2hPa, 2hU , 2hNp, 2hPu, 2hAm, 2hCm, 2hBk, 2hCf, 2hEs, 2hFm,
     & 2hMd, 2hNo, 2hLr /

      data ipstep / 12 /

*-----------------------------------------------------------------------
*     Mass Distribuion of NCASC, NEVAP or Fission
*-----------------------------------------------------------------------

         if( icd .eq. 30 ) then

            kk = 1

            write(ida,'(
     &            ''*'',71(''-'') /
     &            ''*     Isotope Distribution of NCASC [mb]''/
     &            ''*'',71(''-''))')


         else if( icd .eq. 31 ) then

            kk= 2

            write(ida,'(
     &            ''*'',71(''-'') /
     &            ''*     Isotope Distribution of NEVAP [mb]''/
     &            ''*'',71(''-''))')


         else if( icd .eq. 32 ) then

            kk= 3

            write(ida,'(
     &            ''*'',71(''-'') /
     &            ''*     Isotope Distribution of fission [mb]''/
     &            ''*'',71(''-''))')

         end if

*-----------------------------------------------------------------------

      do 170 iz = 1, maxpt
      do 110 in = 1, maxnt

         im = in

  110 if( sumas(kk,iz,in) .gt. 0.0 ) goto 120

      go to 170

  120 continue

      do 130 jn = maxnt, im + 1, -1

         jm = jn

  130 if( sumas(kk,iz,jn) .gt. 0.0 ) go to 140

         jm = im

  140 ms = im + iz - 1

      km = jm - im + 1
      lm = ( km - 1 ) / ipstep + 1

      do 160 nm = 1, lm
      n1 = ipstep * ( nm - 1 ) + 1
      n2 = min( ipstep * nm, km )
      n3 = n1 + im - 1
      n4 = n2 + im - 1
      write(ida,62) iz,hzz(iz)
   62 format(/1x,i4,'-',a2,' isotope production')
      write(ida,63) ( ms + n, n = n1, n2 )
   63 format(' reg.',12i10)

      write(ida,65) 0, ( sumas(kk,iz,i), i = n3, n4 )
   65 format(i4,1x,1p12e10.3)

  160 continue
  170 continue


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sm_mdis(ida,icd)
*                                                                      *
*                                                                      *
*        Last Revised:     2000 03 21                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write mass and charge distribution                   *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              ida         : output unit                               *
*                                                                      *
*              icd         : = 30 ; NCASC                              *
*                              31 ; NEVAP                              *
*                              32 ; fission mass distribution          *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      include 'param01.inc'

*-----------------------------------------------------------------------

      common /summas/ sumas(3,0:maxpt,0:maxnt)
      common /summch/ tmas(3,0:maxpt+maxnt), cdis(3,0:maxpt)
!$OMP THREADPRIVATE(/summch/)

*-----------------------------------------------------------------------

      character       squ*1
      data squ       /"'"/

*-----------------------------------------------------------------------
*     Mass Distribuion of NCASC, NEVAP or Fission
*-----------------------------------------------------------------------

         if( icd .eq. 30 ) then

            kk = 1

            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Mass and Charge Distribution of NCASC''/
     &            ''*'',71(''-''))')


         else if( icd .eq. 31 ) then

            kk= 2

            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Mass and Charge Distribution of NEVAP''/
     &            ''*'',71(''-''))')


         else if( icd .eq. 32 ) then

            kk= 3

            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Mass and Charge Distribution of fission''/
     &            ''*'',71(''-''))')

         end if

            write(ida,'(/
     &            ''p: port nofr'')')

*-----------------------------------------------------------------------

            do i = 0, maxpt + maxnt

               tmas(kk,i) = 0.0

            end do

            do i = 0, maxpt

               cdis(kk,i) = 0.0

            end do

*-----------------------------------------------------------------------

               itmax  = 0
               icmax  = 0
               inmax  = 0

            do iz = 0, maxpt
            do in = 0, maxnt

               if( sumas(kk,iz,in) .gt. 0.0 ) then

                  tmas(kk,iz+in) = tmas(kk,iz+in) + sumas(kk,iz,in)
                  cdis(kk,iz)    = cdis(kk,iz) + sumas(kk,iz,in)

                  if( iz + in .gt. itmax ) itmax = iz + in
                  if( iz      .gt. icmax ) icmax = iz
                  if(      in .gt. inmax ) inmax =      in

               end if

            end do
            end do

               cdis(kk,0) = cdis(kk,0) - tmas(kk,0)

*-----------------------------------------------------------------------
*        Write Mass Distribution
*-----------------------------------------------------------------------

         if( icd .eq. 30 ) then

            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Mass Distribution of NCASC''/
     &            ''*'',71(''-''))')

            write(ida,'(/''p: scal(0.7) xorg(0.1) yorg(1.9)'')')
            write(ida,'( ''p: ylog'')')

            write(ida,'( /a1,''Mass Distribution of NCASC''a1)')
     &                    squ, squ

            write(ida,'(/''x: A_{f}'')')
            write(ida,'( ''y: Yield (mb)'')')

            write(ida,'(/''c:  mass       ncasc'')')
            write(ida,'( ''h:  x-0.5    y(NCASC),lh0'')')

         else if( icd .eq. 31 ) then

            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Mass Distribution of NEVAP''/
     &            ''*'',71(''-''))')

            write(ida,'(/''p: scal(0.7) xorg(0.1) yorg(1.9)'')')
            write(ida,'( ''p: ylog'')')

            write(ida,'( /a1,''Mass Distribution of NEVAP''a1)')
     &                    squ, squ

            write(ida,'(/''x: A_{f}'')')
            write(ida,'( ''y: Yield (mb)'')')

            write(ida,'(/''c:  mass       nevap'')')
            write(ida,'( ''h:  x-0.5    y(NEVAP),lh0'')')

         else if( icd .eq. 32 ) then

            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Mass Distribution of fission''/
     &            ''*'',71(''-''))')

            write(ida,'(/''p: scal(0.7) xorg(0.1) yorg(1.9)'')')
            write(ida,'( ''p: ylog'')')

            write(ida,'( /a1,''Mass Distribution of fission''a1)')
     &                    squ, squ

            write(ida,'(/''x: A_{f}'')')
            write(ida,'( ''y: Yield (mb)'')')

            write(ida,'(/''c:  mass       fission'')')
            write(ida,'( ''h:  x-0.5    y(fission),lh0'')')

         end if

            write(ida,'(f7.1,e16.5)') 1.0, 0.0

         do i = 1, itmax + 1

            write(ida,'(f7.1,e16.5)') dble(i), tmas(kk,i)

         end do


*-----------------------------------------------------------------------
*        Write Charge Distribution of NCASC
*-----------------------------------------------------------------------

         if( icd .eq. 30 ) then

            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Charge Distribution of NCASC''/
     &            ''*'',71(''-''))')

            write(ida,'(/''z: yorg(-1.7)'')')
            write(ida,'( ''p: ylog'')')

            write(ida,'(/a1,''Charge Distribution of NCASC''a1)')
     &                    squ, squ

            write(ida,'(/''x: Z_{f}'')')
            write(ida,'( ''y: Yield (mb)'')')

            write(ida,'(/''c:  charge     ncasc'')')
            write(ida,'( ''h:  x-0.5    y(NCASC),lh0'')')

         else if( icd .eq. 31 ) then

            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Charge Distribution of NEVAP''/
     &            ''*'',71(''-''))')

            write(ida,'(/''z: yorg(-1.7)'')')
            write(ida,'( ''p: ylog'')')

            write(ida,'(/a1,''Charge Distribution of NEVAP''a1)')
     &                    squ, squ

            write(ida,'(/''x: Z_{f}'')')
            write(ida,'( ''y: Yield (mb)'')')

            write(ida,'(/''c:  charge     nevap'')')
            write(ida,'( ''h:  x-0.5    y(NEVAP),lh0'')')

         else if( icd .eq. 32 ) then

            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Charge Distribution of fission''/
     &            ''*'',71(''-''))')

            write(ida,'(/''z: yorg(-1.7)'')')
            write(ida,'( ''p: ylog'')')

            write(ida,'(/a1,''Charge Distribution of fission''a1)')
     &                    squ, squ

            write(ida,'(/''x: Z_{f}'')')
            write(ida,'( ''y: Yield (mb)'')')

            write(ida,'(/''c:  charge     fission'')')
            write(ida,'( ''h:  x-0.5    y(fission),lh0'')')

         end if

            write(ida,'(f7.1,e16.5)') 1.0, 0.0

         do i = 1, icmax + 1

            write(ida,'(f7.1,e16.5)') dble(i), cdis(kk,i)

         end do

*-----------------------------------------------------------------------
*        write mass distribution in 2-Dim. of NCASC
*-----------------------------------------------------------------------

         if( icd .eq. 30 ) then

            write(ida,'(/
     &         ''*'',71(''-'') /
     &         ''*     Plot of Mass Distribution in 2-Dim. of NCASC''/
     &         ''*'',71(''-''))')

            write(ida,'(/''newpage:'')')

         else if( icd .eq. 31 ) then

            write(ida,'(/
     &         ''*'',71(''-'') /
     &         ''*     Plot of Mass Distribution in 2-Dim. of NEVAP''/
     &         ''*'',71(''-''))')

            write(ida,'(/''newpage:'')')

         else if( icd .eq. 32 ) then

            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Mass Distribution in 2-Dim. '',
     &            ''of fission''/
     &            ''*'',71(''-''))')

            write(ida,'(/''newpage:'')')

         end if

            write(ida,'(/
     &            ''p: port nofr'')')

*-----------------------------------------------------------------------

            dxmax = dble(inmax+2)
            dymax = dble(icmax+2)
            dform = dymax / dxmax

*-----------------------------------------------------------------------

            write(ida,'( /''p: nosp zlog'')')
            write(ida,'(  ''p: xorg(-0.05) yorg(0.2) afac(0.6)'')')
            write(ida,'(  ''p: xmin(0) xmax('',f5.1,'')'')') dxmax
            write(ida,'(  ''p: ymin(0) ymax('',f5.1,'')'')') dymax
            write(ida,'(  ''p: form('',f7.3,'')'')') dform

            write(ida,'(/''x: N Neutron Number'')')
            write(ida,'( ''y: Z Proton Number'')')

            write(ida,'( ''hc: y = '',i3,'' to 1 by -1 ;'',
     &                      '' x = 1 to '',i3,'' by 1 ;'')')
     &                     icmax+2, inmax+2

         do i = icmax+2, 1, -1

            write(ida,'(10e11.3)')
     &           ( sumas(kk,i,l), l = 1, inmax+2 )

         end do

*-----------------------------------------------------------------------
*        write magic numbers
*        write stable nuclei
*-----------------------------------------------------------------------

            call wmgcstb(ida,icmax,inmax)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sm_mass(iss,iz1,in1,sfact)
*                                                                      *
*                                                                      *
*        Last Revised:     2000 03 21                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to store the mass distribution of the clusters          *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              iss            =1 ; after NCASC                         *
*                             =2 ; after NEVAP                         *
*                             =3 ; before fission                      *
*                                                                      *
*              iz1            : proton number                          *
*              in1            : neutron number                         *
*                                                                      *
*              sfact          : weight facter                          *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      include 'param01.inc'

      common /summas/ sumas(3,0:maxpt,0:maxnt)

*-----------------------------------------------------------------------

               ipro = iz1
               ineu = in1

            if( ipro .le. maxpt .and. ineu .le. maxnt ) then

               sumas(iss,ipro,ineu) = sumas(iss,ipro,ineu) + sfact

            end if


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sm_norm(rmn)
*                                                                      *
*                                                                      *
*        Last Revised:     2000 03 21                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*             to normarize the cross section by the total event number *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      include 'param01.inc'

*-----------------------------------------------------------------------

      common /summas/ sumas(3,0:maxpt,0:maxnt)

*-----------------------------------------------------------------------
*        normalization of total event number : iprun
*-----------------------------------------------------------------------

               do k = 1, 3
               do i = 0, maxpt
               do j = 0, maxnt

                  sumas(k,i,j) = sumas(k,i,j) * rmn

               end do
               end do
               end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sm_mini
*                                                                      *
*                                                                      *
*        Last Revised:     2000 03 21                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*             to initialize the mass distribution cross section        *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      include 'param01.inc'

*-----------------------------------------------------------------------

      common /summas/ sumas(3,0:maxpt,0:maxnt)

*-----------------------------------------------------------------------
*        normalization of total event number : iprun
*-----------------------------------------------------------------------

               do k = 1, 3
               do i = 0, maxpt
               do j = 0, maxnt

                  sumas(k,i,j) = 0.0d0

               end do
               end do
               end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine wmgcstb(ida,icmax,inmax)
*                                                                      *
*       write magic number and stable nuclei in 2D-plot of yield       *
*       last modified by K.Niita on 04/08/2000                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------
*        write magic numbers
*-----------------------------------------------------------------------

            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Magic Numbers in 2-Dim.''/
     &            ''*'',71(''-''))')

         if( 14 .lt. inmax+2 .and. 4 .lt. icmax+2 ) then
            write(ida,'(/''w:2/y(2) x(12) iy(2) ix(1) s(0.7)'')')
         else
            write(ida,'(/''n:2/y(2) x(12) iy(2) ix(1) s(0.7)'')')
         end if

         if( 4 .lt. inmax+2 .and.  14 .lt. icmax+2 ) then
            write(ida,'(/''w:2/x(2) y(12) ix(2) iy(1) s(0.7)'')')
         else
            write(ida,'(/''n:2/x(2) y(12) ix(2) iy(1) s(0.7)'')')
         end if

         if( 18 .lt. inmax+2 .and. 10 .lt. icmax+2 ) then
            write(ida,'(/''w:8/y(8) x(15) iy(2) ix(1) s(0.7)'')')
         else
            write(ida,'(/''n:8/y(8) x(15) iy(2) ix(1) s(0.7)'')')
         end if

         if( 10 .lt. inmax+2 .and.  18 .lt. icmax+2 ) then
            write(ida,'(/''w:8/x(8) y(15) ix(2) iy(1) s(0.7)'')')
         else
            write(ida,'(/''n:8/x(8) y(15) ix(2) iy(1) s(0.7)'')')
         end if

            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''   0  2.5'')')
            write(ida,'( ''  11  2.5'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''   0  1.5'')')
            write(ida,'( ''  11  1.5'')')

            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  2.5   0'')')
            write(ida,'( ''  2.5   11'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  1.5   0'')')
            write(ida,'( ''  1.5   11'')')

            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''   0  8.5'')')
            write(ida,'( ''  14  8.5'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''   0  7.5'')')
            write(ida,'( ''  14  7.5'')')

            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  8.5   0'')')
            write(ida,'( ''  8.5   14'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  7.5   0'')')
            write(ida,'( ''  7.5   14'')')

*-----------------------------------------------------------------------

         if( 11 .lt. inmax+2 .and. 20 .lt. icmax+2 ) then
            write(ida,'(/''w:20/y(20) x(11) iy(2) ix(3) s(0.7)'')')
         else
            write(ida,'(/''n:20/y(20) x(11) iy(2) ix(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  12  20.5'')')
            write(ida,'( ''  52  20.5'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  12  19.5'')')
            write(ida,'( ''  52  19.5'')')

         if( 17 .lt. inmax+2 .and. 28 .lt. icmax+2 ) then
            write(ida,'(/''w:28/y(28) x(17) iy(2) ix(3) s(0.7)'')')
         else
            write(ida,'(/''n:28/y(28) x(17) iy(2) ix(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  18  28.5'')')
            write(ida,'( ''  84  28.5'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  18  27.5'')')
            write(ida,'( ''  84  27.5'')')

         if( 47 .lt. inmax+2 .and. 50 .lt. icmax+2 ) then
            write(ida,'(/''w:50/y(50) x(47) iy(2) ix(3) s(0.7)'')')
         else
            write(ida,'(/''n:50/y(50) x(47) iy(2) ix(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  48   50.5'')')
            write(ida,'( ''  128  50.5'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  48   49.5'')')
            write(ida,'( ''  128  49.5'')')

         if( 79 .lt. inmax+2 .and. 82 .lt. icmax+2 ) then
            write(ida,'(/''w:82/y(82) x(79) iy(2) ix(3) s(0.7)'')')
         else
            write(ida,'(/''n:82/y(82) x(79) iy(2) ix(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  80   82.5'')')
            write(ida,'( ''  150  82.5'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  80   81.5'')')
            write(ida,'( ''  150  81.5'')')

         if( 20 .lt. inmax+2 .and.  7 .lt. icmax+2 ) then
            write(ida,'(/''w:20/x(20) y(7) ix(2) iy(3) s(0.7)'')')
         else
            write(ida,'(/''n:20/x(20) y(7) ix(2) iy(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  20.5   8'')')
            write(ida,'( ''  20.5   30'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  19.5   8'')')
            write(ida,'( ''  19.5   30'')')

         if( 28 .lt. inmax+2 .and.  7 .lt. icmax+2 ) then
            write(ida,'(/''w:28/x(28) y(7) ix(2) iy(3) s(0.7)'')')
         else
            write(ida,'(/''n:28/x(28) y(7) ix(2) iy(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  28.5   8'')')
            write(ida,'( ''  28.5   36'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  27.5   8'')')
            write(ida,'( ''  27.5   36'')')

         if( 50 .lt. inmax+2 .and. 17 .lt. icmax+2 ) then
            write(ida,'(/''w:50/x(50) y(17) ix(2) iy(3) s(0.7)'')')
         else
            write(ida,'(/''n:50/x(50) y(17) ix(2) iy(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  50.5   18'')')
            write(ida,'( ''  50.5   52'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  49.5   18'')')
            write(ida,'( ''  49.5   52'')')

         if( 82 .lt. inmax+2 .and. 25 .lt. icmax+2 ) then
            write(ida,'(/''w:82/x(82) y(25) ix(2) iy(3) s(0.7)'')')
         else
            write(ida,'(/''n:82/x(82) y(25) ix(2) iy(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  82.5   26'')')
            write(ida,'( ''  82.5   84'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  81.5   26'')')
            write(ida,'( ''  81.5   84'')')

         if( 126 .lt. inmax+2 .and. 47 .lt. icmax+2 ) then
            write(ida,'(/''w:126/x(126) y(47) ix(2) iy(3) s(0.7)'')')
         else
            write(ida,'(/''n:126/x(126) y(47) ix(2) iy(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  126.5   48'')')
            write(ida,'( ''  126.5   100'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  125.5   48'')')
            write(ida,'( ''  125.5   100'')')

*-----------------------------------------------------------------------
*        write stable nuclei
*-----------------------------------------------------------------------

            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Stable Nuclei in 2-Dim.''/
     &            ''*'',71(''-''))')

         if( inmax .gt. 100 ) then

            write(ida,'( ''h: x y,n6xxxxx'')')

         else if( inmax .gt. 75 ) then

            write(ida,'( ''h: x y,n6xxxx'')')

         else if( inmax .gt. 50 ) then

            write(ida,'( ''h: x y,n6xxx'')')

         else if( inmax .gt. 30 ) then

            write(ida,'( ''h: x y,n6xx'')')

         else

            write(ida,'( ''h: x y,n6'')')

         end if

            write(ida,'( ''  0   1'')')
            write(ida,'( ''  1   1'')')
            write(ida,'( ''  1   2'')')
            write(ida,'( ''  2   2'')')
            write(ida,'( ''  3   3'')')
            write(ida,'( ''  4   3'')')
            write(ida,'( ''  5   4'')')
            write(ida,'( ''  5   5'')')
            write(ida,'( ''  6   5'')')
            write(ida,'( ''  6   6'')')
            write(ida,'( ''  7   6'')')
            write(ida,'( ''  7   7'')')
            write(ida,'( ''  8   7'')')
            write(ida,'( ''  8   8'')')
            write(ida,'( ''  9   8'')')
            write(ida,'( ''  10  8'')')
            write(ida,'( ''  10  9'')')
            write(ida,'( ''  10  10'')')
            write(ida,'( ''  11  10'')')
            write(ida,'( ''  12  10'')')
            write(ida,'( ''  12  11'')')
            write(ida,'( ''  12  12'')')
            write(ida,'( ''  13  12'')')
            write(ida,'( ''  14  12'')')
            write(ida,'( ''  14  13'')')
            write(ida,'( ''  14  14'')')
            write(ida,'( ''  15  14'')')
            write(ida,'( ''  16  14'')')
            write(ida,'( ''  16  15'')')
            write(ida,'( ''  16  16'')')
            write(ida,'( ''  17  16'')')
            write(ida,'( ''  18  16'')')
            write(ida,'( ''  18  17'')')
            write(ida,'( ''  18  18'')')
            write(ida,'( ''  20  16'')')
            write(ida,'( ''  20  17'')')
            write(ida,'( ''  20  18'')')
            write(ida,'( ''  20  19'')')
            write(ida,'( ''  20  20'')')
            write(ida,'( ''  21  19'')')
            write(ida,'( ''  22  18'')')
            write(ida,'( ''  22  19'')')
            write(ida,'( ''  22  20'')')
            write(ida,'( ''  23  20'')')
            write(ida,'( ''  24  20'')')
            write(ida,'( ''  24  21'')')
            write(ida,'( ''  24  22'')')
            write(ida,'( ''  25  22'')')
            write(ida,'( ''  26  20'')')
            write(ida,'( ''  26  22'')')
            write(ida,'( ''  26  24'')')
            write(ida,'( ''  27  22'')')
            write(ida,'( ''  27  23'')')
            write(ida,'( ''  28  20'')')
            write(ida,'( ''  28  22'')')
            write(ida,'( ''  28  23'')')
            write(ida,'( ''  28  24'')')
            write(ida,'( ''  28  26'')')
            write(ida,'( ''  29  24'')')
            write(ida,'( ''  30  24'')')
            write(ida,'( ''  30  25'')')
            write(ida,'( ''  30  26'')')
            write(ida,'( ''  30  28'')')
            write(ida,'( ''  31  26'')')
            write(ida,'( ''  32  26'')')
            write(ida,'( ''  32  27'')')
            write(ida,'( ''  32  28'')')
            write(ida,'( ''  33  28'')')
            write(ida,'( ''  34  28'')')
            write(ida,'( ''  34  29'')')
            write(ida,'( ''  34  30'')')
            write(ida,'( ''  36  28'')')
            write(ida,'( ''  36  29'')')
            write(ida,'( ''  36  30'')')
            write(ida,'( ''  37  30'')')
            write(ida,'( ''  38  30'')')
            write(ida,'( ''  38  31'')')
            write(ida,'( ''  38  32'')')
            write(ida,'( ''  40  30'')')
            write(ida,'( ''  40  31'')')
            write(ida,'( ''  40  32'')')
            write(ida,'( ''  40  34'')')
            write(ida,'( ''  41  32'')')
            write(ida,'( ''  42  32'')')
            write(ida,'( ''  42  33'')')
            write(ida,'( ''  42  34'')')
            write(ida,'( ''  42  36'')')
            write(ida,'( ''  43  34'')')
            write(ida,'( ''  44  32'')')
            write(ida,'( ''  44  34'')')
            write(ida,'( ''  44  35'')')
            write(ida,'( ''  44  36'')')
            write(ida,'( ''  46  34'')')
            write(ida,'( ''  46  35'')')
            write(ida,'( ''  46  36'')')
            write(ida,'( ''  46  38'')')
            write(ida,'( ''  47  36'')')
            write(ida,'( ''  48  34'')')
            write(ida,'( ''  48  36'')')
            write(ida,'( ''  48  37'')')
            write(ida,'( ''  48  38'')')
            write(ida,'( ''  49  38'')')
            write(ida,'( ''  50  36'')')
            write(ida,'( ''  50  37'')')
            write(ida,'( ''  50  38'')')
            write(ida,'( ''  50  39'')')
            write(ida,'( ''  50  40'')')
            write(ida,'( ''  50  42'')')
            write(ida,'( ''  51  40'')')
            write(ida,'( ''  52  40'')')
            write(ida,'( ''  52  41'')')
            write(ida,'( ''  52  42'')')
            write(ida,'( ''  52  44'')')
            write(ida,'( ''  53  42'')')
            write(ida,'( ''  54  40'')')
            write(ida,'( ''  54  42'')')
            write(ida,'( ''  54  44'')')
            write(ida,'( ''  55  42'')')
            write(ida,'( ''  55  44'')')
            write(ida,'( ''  56  40'')')
            write(ida,'( ''  56  42'')')
            write(ida,'( ''  56  44'')')
            write(ida,'( ''  56  46'')')
            write(ida,'( ''  57  44'')')
            write(ida,'( ''  58  42'')')
            write(ida,'( ''  58  44'')')
            write(ida,'( ''  58  45'')')
            write(ida,'( ''  58  46'')')
            write(ida,'( ''  58  48'')')
            write(ida,'( ''  59  46'')')
            write(ida,'( ''  60  44'')')
            write(ida,'( ''  60  46'')')
            write(ida,'( ''  60  47'')')
            write(ida,'( ''  60  48'')')
            write(ida,'( ''  62  46'')')
            write(ida,'( ''  62  47'')')
            write(ida,'( ''  62  48'')')
            write(ida,'( ''  62  50'')')
            write(ida,'( ''  63  47'')')
            write(ida,'( ''  64  46'')')
            write(ida,'( ''  64  48'')')
            write(ida,'( ''  64  49'')')
            write(ida,'( ''  64  50'')')
            write(ida,'( ''  65  48'')')
            write(ida,'( ''  65  50'')')
            write(ida,'( ''  66  48'')')
            write(ida,'( ''  66  49'')')
            write(ida,'( ''  66  50'')')
            write(ida,'( ''  67  50'')')
            write(ida,'( ''  68  48'')')
            write(ida,'( ''  68  50'')')
            write(ida,'( ''  68  52'')')
            write(ida,'( ''  69  50'')')
            write(ida,'( ''  70  50'')')
            write(ida,'( ''  70  51'')')
            write(ida,'( ''  70  52'')')
            write(ida,'( ''  70  54'')')
            write(ida,'( ''  71  52'')')
            write(ida,'( ''  72  51'')')
            write(ida,'( ''  72  52'')')
            write(ida,'( ''  72  53'')')
            write(ida,'( ''  72  54'')')
            write(ida,'( ''  73  52'')')
            write(ida,'( ''  74  50'')')
            write(ida,'( ''  74  52'')')
            write(ida,'( ''  74  53'')')
            write(ida,'( ''  74  54'')')
            write(ida,'( ''  74  56'')')
            write(ida,'( ''  75  54'')')
            write(ida,'( ''  76  52'')')
            write(ida,'( ''  76  54'')')
            write(ida,'( ''  76  56'')')
            write(ida,'( ''  77  54'')')
            write(ida,'( ''  78  52'')')
            write(ida,'( ''  78  54'')')
            write(ida,'( ''  78  55'')')
            write(ida,'( ''  78  56'')')
            write(ida,'( ''  78  58'')')
            write(ida,'( ''  79  56'')')
            write(ida,'( ''  80  54'')')
            write(ida,'( ''  80  56'')')
            write(ida,'( ''  80  58'')')
            write(ida,'( ''  81  56'')')
            write(ida,'( ''  81  57'')')
            write(ida,'( ''  82  54'')')
            write(ida,'( ''  82  56'')')
            write(ida,'( ''  82  57'')')
            write(ida,'( ''  82  58'')')
            write(ida,'( ''  82  59'')')
            write(ida,'( ''  82  60'')')
            write(ida,'( ''  82  62'')')
            write(ida,'( ''  83  60'')')
            write(ida,'( ''  84  58'')')
            write(ida,'( ''  84  60'')')
            write(ida,'( ''  85  60'')')
            write(ida,'( ''  85  62'')')
            write(ida,'( ''  86  60'')')
            write(ida,'( ''  86  62'')')
            write(ida,'( ''  87  62'')')
            write(ida,'( ''  88  60'')')
            write(ida,'( ''  88  62'')')
            write(ida,'( ''  88  63'')')
            write(ida,'( ''  88  64'')')
            write(ida,'( ''  90  60'')')
            write(ida,'( ''  90  62'')')
            write(ida,'( ''  90  63'')')
            write(ida,'( ''  90  64'')')
            write(ida,'( ''  90  66'')')
            write(ida,'( ''  91  64'')')
            write(ida,'( ''  92  62'')')
            write(ida,'( ''  92  64'')')
            write(ida,'( ''  92  66'')')
            write(ida,'( ''  93  64'')')
            write(ida,'( ''  94  64'')')
            write(ida,'( ''  94  65'')')
            write(ida,'( ''  94  66'')')
            write(ida,'( ''  94  68'')')
            write(ida,'( ''  95  66'')')
            write(ida,'( ''  96  64'')')
            write(ida,'( ''  96  66'')')
            write(ida,'( ''  96  68'')')
            write(ida,'( ''  97  66'')')
            write(ida,'( ''  98  66'')')
            write(ida,'( ''  98  67'')')
            write(ida,'( ''  98  68'')')
            write(ida,'( ''  98  70'')')
            write(ida,'( ''  99  68'')')
            write(ida,'( ''  100 68'')')
            write(ida,'( ''  100 69'')')
            write(ida,'( ''  100 70'')')
            write(ida,'( ''  101 70'')')
            write(ida,'( ''  102 68'')')
            write(ida,'( ''  102 70'')')
            write(ida,'( ''  102 72'')')
            write(ida,'( ''  103 70'')')
            write(ida,'( ''  104 70'')')
            write(ida,'( ''  104 71'')')
            write(ida,'( ''  105 71'')')
            write(ida,'( ''  105 72'')')
            write(ida,'( ''  106 70'')')
            write(ida,'( ''  106 72'')')
            write(ida,'( ''  106 74'')')
            write(ida,'( ''  107 72'')')
            write(ida,'( ''  107 73'')')
            write(ida,'( ''  108 72'')')
            write(ida,'( ''  108 73'')')
            write(ida,'( ''  108 74'')')
            write(ida,'( ''  108 76'')')
            write(ida,'( ''  109 74'')')
            write(ida,'( ''  110 74'')')
            write(ida,'( ''  110 75'')')
            write(ida,'( ''  110 76'')')
            write(ida,'( ''  111 76'')')
            write(ida,'( ''  112 74'')')
            write(ida,'( ''  112 75'')')
            write(ida,'( ''  112 76'')')
            write(ida,'( ''  112 78'')')
            write(ida,'( ''  113 76'')')
            write(ida,'( ''  114 76'')')
            write(ida,'( ''  114 77'')')
            write(ida,'( ''  114 78'')')
            write(ida,'( ''  116 76'')')
            write(ida,'( ''  116 77'')')
            write(ida,'( ''  116 78'')')
            write(ida,'( ''  116 80'')')
            write(ida,'( ''  117 78'')')
            write(ida,'( ''  118 78'')')
            write(ida,'( ''  118 79'')')
            write(ida,'( ''  118 80'')')
            write(ida,'( ''  119 80'')')
            write(ida,'( ''  120 78'')')
            write(ida,'( ''  120 80'')')
            write(ida,'( ''  121 80'')')
            write(ida,'( ''  122 80'')')
            write(ida,'( ''  122 81'')')
            write(ida,'( ''  122 82'')')
            write(ida,'( ''  124 80'')')
            write(ida,'( ''  124 81'')')
            write(ida,'( ''  124 82'')')
            write(ida,'( ''  125 82'')')
            write(ida,'( ''  126 82'')')
            write(ida,'( ''  126 83'')')

*-----------------------------------------------------------------------

      return
      end

