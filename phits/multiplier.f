************************************************************************
*                                                                      *
      subroutine multech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the multiplier section                           *
*       last modified by N.Matuda on 2024/12/09                        *
*                                                                      *
************************************************************************
      use moddas_multiplier
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /multipl/ imltp,iimlt(multmax),inmlt(multmax),
C MATSUDA 2024.11.25 (lagrange: ilmlt)
     &  ilmlt(multmax),
     &  idmlt(multmax),ismlt(multmax),impan(multmax),impat(multmax,6,3),
     &                 jmpat(multmax,6,6,2), imdfl(multmax)

C MATSUDA 2024.11.25 (multplf: file name)
      common /multplf/ imltf(multmax),lmltfile(multmax),mltfile(multmax)
      character mltfile*100
C MATSUDA 2024.12.09 (multpl27,35: x-, y-txt, and epsout)
      common /multpl27/ impxl(multmax),impxt(multmax),
     &                  impyl(multmax),impyt(multmax)
      character impxt*200, impyt*200
      common /multpl35/ impeps(multmax)

*-----------------------------------------------------------------------

      character asfil*100

*-----------------------------------------------------------------------
*     input echo (Multiplier)
*-----------------------------------------------------------------------

            write(iot,'( "[ Multiplier ]")')

               call echprt(1,iot,m,impan,impat,jmpat,multmax,6,6)

*-----------------------------------------------------------------------

            write(iot,'("  number = ",i4)') idmlt(m)

*-----------------------------------------------------------------------
*           iimlt(m) =  1 / -1 : lin-lin / log-log
*                    = -2 /  2 : glow / ghigh (constant value)
*                    =  3 /  4 : lin-log ("xlin" or "ylog")
*                    = -3 / -4 : log-lin ("xlog" or "ylin")
*-----------------------------------------------------------------------

            if( iimlt(m) .eq. 1 ) then
               write(iot,'("   interpolation = lin")')
            else if( iimlt(m) .eq. -1 ) then
               write(iot,'("   interpolation = log")')
            else if( iimlt(m) .eq. -2 ) then
               write(iot,'("   interpolation = glow")')
            else if( iimlt(m) .eq.  2 ) then
               write(iot,'("   interpolation = ghigh")')
            else if( iimlt(m) .eq.  3 ) then
               write(iot,'("   interpolation = xlin")')
            else if( iimlt(m) .eq.  4 ) then
               write(iot,'("   interpolation = ylog")')
            else if( iimlt(m) .eq. -3 ) then
               write(iot,'("   interpolation = xlog")')
            else if( iimlt(m) .eq. -4 ) then
               write(iot,'("   interpolation = ylin")')
            end if

*-----------------------------------------------------------------------

            write(iot,'("  lagrange = ",i4)') ilmlt(m)

*-----------------------------------------------------------------------

            write(iot,'("  ne = ",i4)') inmlt(m)

               igr = inmlt(m)
               igm = ismlt(m) - 1

            do i = 1, igr

               write(iot,'(1p2e13.5)')
     &            gmsh_ismlt(igm+2*i-1), gmsh_ismlt(igm+2*i)

            end do

*-----------------------------------------------------------------------

                  asfil = '  # file name for multiplier'

            if ( imltf(m) .gt. 0 ) then

                  msfile = max( 14, lmltfile(m) )

                  write(iot,'(/"     file = ",100a1)')
     &                  ( mltfile(m)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 28 )

            end if

*-----------------------------------------------------------------------

               if( impxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( impxt(m)(i:i),i = 1, impxl(m) )

               end if

*-----------------------------------------------------------------------

               if( impyl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( impyt(m)(i:i),i = 1, impyl(m) )

               end if

*-----------------------------------------------------------------------

               if( impeps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            impeps(m)

               end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine pmultip(m,np)
*                                                                      *
*       output the multiplier section with interpolation functions     *
*       Last modified by N.Matsuda on 2024/12/09                       *
*                                                                      *
************************************************************************
      use moddas_multiplier
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

*-----------------------------------------------------------------------

      common /multipl/ imltp,iimlt(multmax),inmlt(multmax),
C MATSUDA 2024.11.25 (lagrange: ilmlt)
     &  ilmlt(multmax),
     &  idmlt(multmax),ismlt(multmax),impan(multmax),impat(multmax,6,3),
     &                 jmpat(multmax,6,6,2), imdfl(multmax)
      dimension cdata(8), coefficients(4)
      real*8  cdata, coefficients

C MATSUDA 2024.11.25 (multplf: file name)
      common /multplf/ imltf(multmax),lmltfile(multmax),mltfile(multmax)
      character mltfile*100
C MATSUDA 2024.12.09 (multpl27,35: x-, y-txt, and epsout)
      common /multpl27/ impxl(multmax),impxt(multmax),
     &                  impyl(multmax),impyt(multmax)
      character impxt*200, impyt*200
      common /multpl35/ impeps(multmax)

      character fname*100, fnume*3

      common /ptname/ pname(20), ipln(20)
      character       pname*8

      character kname*9
      character tname*9
      character chau*8

*-----------------------------------------------------------------------

      character yen*1
      yen  = char(92)

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

            do i = 1, impan(m)

                        tname = '         '

                        istyp = impat(m,i,1)
                        inkf0 = impat(m,i,2)
                        isubt = 0

                        ik = 1    ! if( istyp < 0 ) ik > 1
                        kk = 0    !                 kk = 1

               do k = 1, ik

C               NUCLEUS: istyp == 19
                  if( istyp .eq. 19 .and. inkf0 .ne. 0 ) then
                        iz = ichgf(istyp,abs(inkf0))
                        ia = ibryf(istyp,abs(inkf0))
                        call chname(idum,ia,iz,chau)

                        tname = chau(1:8)

C               NORMAL:  istyp != 11, 19
                  else if( istyp .ne. 11 ) then
                        tname = pname(istyp)(1:8)

C               OTHER:   istyp == 11
                  else
                        call kfcname(inkf0,9,tname)

                  end if

               end do

            end do

*-----------------------------------------------------------------------
*        minimum and maximum values
*-----------------------------------------------------------------------

               igr = inmlt(m)
               igm = ismlt(m) - 1

                  xmax = gmsh_ismlt(igm+2*igr-1)
                  xmin = gmsh_ismlt(igm+1)

                  cmax = 0.0
                  cmin = 1.e+33

            do i = 1, igr

                  if( gmsh_ismlt(igm+2*i) .gt. cmax )
     &                        cmax = gmsh_ismlt(igm+2*i)

                  if( gmsh_ismlt(igm+2*i) .lt. cmin )
     &                        cmin = gmsh_ismlt(igm+2*i)

            end do

*-----------------------------------------------------------------------
*     out put unit = 15  : temporary number
*-----------------------------------------------------------------------

      if( imltf(m) .gt. 0 ) then

            fname = mltfile(m)

            iot = 15
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            call multech(iot,m,1,1)

*-----------------------------------------------------------------------

               inum = 0
C              np = 1

C           do ipi = 1, np    ! for part or axis

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

*-----------------------------------------------------------------------

                  form  =  1.0
                  xfac  =  1.0

                  write(iot,'( "set: c1[",f6.3,
     &                            " ] c2[",f6.3," ]")')
     &                    form, xfac
                  write(iot,'( "p: form[c1/c2] xfac[c2]  nosx")')

               if( impxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Energy (MeV)")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (impxt(m)(i:i),i=1,impxl(m))

               end if

               if( impyl(m) .eq. 0 ) then

                  write(iot,'( "y: Coefficients")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (impyt(m)(i:i),i=1,impyl(m))

               end if

               if( iimlt(m) .eq.  1 .or.
     &             iimlt(m) .eq.  3 .or.
     &             iimlt(m) .eq.  4 ) then
                  write(iot,'(/"p: xlin ")', advance='no')
               else
                  write(iot,'(/"p: xlog ")', advance='no')
               end if
               if( iimlt(m) .eq.  1 .or.
     &             iimlt(m) .eq. -3 .or.
     &             iimlt(m) .eq. -4 ) then
                  write(iot,'( "ylin  nosx")')
               else
                  write(iot,'( "ylog  nosx")')
               end if

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') xmin, xmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') cmin, cmax

C              if( cmin .gt. 0.0 .and. cmax .gt. cmin ) then
               if( cmax .gt. cmin ) then

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax
               write(iot,'( "p: ymin[c3] ymax[c4]")')

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"h:  x           y,n3RXX")')

               do i = 1, igr

                  write(iot,'(1p2e13.5)')
     &               gmsh_ismlt(igm+2*i-1), gmsh_ismlt(igm+2*i)

               end do

               write(iot,'(/"msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( iimlt(m) .eq.  1 .or.
     &             iimlt(m) .eq.  3 .or.
     &             iimlt(m) .eq.  4 ) then
                  xfact = xmin + ( xmax - xmin ) * 0.64
               else
                  xfact = log( xmin ) + log( xmax/xmin ) * 0.64
                  xfact = exp( xfact )
               end if
               if( iimlt(m) .eq.  1 .or.
     &             iimlt(m) .eq. -3 .or.
     &             iimlt(m) .eq. -4 ) then
                  cfact = cmin + ( cmax - cmin ) * 0.10
               else
                  cfact = log( cmin ) + log( cmax/cmin ) * 0.10
                  cfact = exp( cfact )
               end if

               write(iot,'(/"w: ",a9,"/ x(",1p1g13.6,") y(",1p1g13.6,
     &                  ") s(2.4)")') tname, xfact, cfact

*-----------------------------------------------------------------------

                  write(iot,'( "")')

C         GROUP-WISE
            if( iimlt(m) .eq. -2 .or. iimlt(m) .eq. 2 ) then

               if( iimlt(m) .eq. -2 ) then
                  write(iot,'( "h:  x           y,lhBZ")')
               else if( iimlt(m) .eq. 2 ) then
                  write(iot,'( "h:  x           y,lhhBZ")')
               end if

               do i = 1, igr

                  write(iot,'(1p2e13.5)')
     &               gmsh_ismlt(igm+2*i-1), gmsh_ismlt(igm+2*i)

               end do

C         POINT-WISE
            else if( iimlt(m) .eq.  1 .or. iimlt(m) .eq. -1 .or.
     &               iimlt(m) .eq.  3 .or. iimlt(m) .eq.  4 .or.
     &               iimlt(m) .eq. -3 .or. iimlt(m) .eq. -4 ) then

               do i = 1, igr-1

C              Lagrange-4, BUT -3 (3 points)
                 if( i .eq. 1 .and. ilmlt(m) .gt. 2 ) then
                    do ii = 1, 3
                      if( iimlt(m) .eq.  1 .or.
     &                    iimlt(m) .eq.  3 .or.
     &                    iimlt(m) .eq.  4 ) then
                        cdata(ii) = gmsh_ismlt(igm+2*(i+ii-1)-1)
                      else
                        cdata(ii) = log(gmsh_ismlt(igm+2*(i+ii-1)-1))
                      end if
                      if( iimlt(m) .eq.  1 .or.
     &                    iimlt(m) .eq. -3 .or.
     &                    iimlt(m) .eq. -4 ) then
                        cdata(4+ii) = gmsh_ismlt(igm+2*(i+ii-1))
                      else
                        cdata(4+ii) = log(gmsh_ismlt(igm+2*(i+ii-1)))
                      end if
                    end do

                    call autofunc_lagrange(3,cdata,coefficients)

C              Lagrange-4, BUT -3 (3 points, right)
                 else if( i .eq. igr-1 .and. ilmlt(m) .gt. 2 ) then
                    do ii = 1, 3
                      if( iimlt(m) .eq.  1 .or.
     &                    iimlt(m) .eq.  3 .or.
     &                    iimlt(m) .eq.  4 ) then
                        cdata(ii) = gmsh_ismlt(igm+2*(i+ii-1)-3)
                      else
                        cdata(ii) = log(gmsh_ismlt(igm+2*(i+ii-1)-3))
                      end if
                      if( iimlt(m) .eq.  1 .or.
     &                    iimlt(m) .eq. -3 .or.
     &                    iimlt(m) .eq. -4 ) then
                        cdata(4+ii) = gmsh_ismlt(igm+2*(i+ii-1)-2)
                      else
                        cdata(4+ii) = log(gmsh_ismlt(igm+2*(i+ii-1)-2))
                      end if
                    end do

                    call autofunc_lagrange(3,cdata,coefficients)

C              Lagrange-4
                 else if( ilmlt(m) .eq. 4 ) then
                    do ii = 1, 4
                      if( iimlt(m) .eq.  1 .or.
     &                    iimlt(m) .eq.  3 .or.
     &                    iimlt(m) .eq.  4 ) then
                        cdata(ii) = gmsh_ismlt(igm+2*(i+ii-1)-3)
                      else
                        cdata(ii) = log(gmsh_ismlt(igm+2*(i+ii-1)-3))
                      end if
                      if( iimlt(m) .eq.  1 .or.
     &                    iimlt(m) .eq. -3 .or.
     &                    iimlt(m) .eq. -4 ) then
                        cdata(4+ii) = gmsh_ismlt(igm+2*(i+ii-1)-2)
                      else
                        cdata(4+ii) = log(gmsh_ismlt(igm+2*(i+ii-1)-2))
                      end if
                    end do

                    call autofunc_lagrange(4,cdata,coefficients)

C              Lagrange-2 or -3 (2 or 3 points)
                 else
                    do ii = 1, ilmlt(m)
                      if( iimlt(m) .eq.  1 .or.
     &                    iimlt(m) .eq.  3 .or.
     &                    iimlt(m) .eq.  4 ) then
                        cdata(ii) = gmsh_ismlt(igm+2*(i+ii-1)-1)
                      else
                        cdata(ii) = log(gmsh_ismlt(igm+2*(i+ii-1)-1))
                      end if
                      if( iimlt(m) .eq.  1 .or.
     &                    iimlt(m) .eq. -3 .or.
     &                    iimlt(m) .eq. -4 ) then
                        cdata(4+ii) = gmsh_ismlt(igm+2*(i+ii-1))
                      else
                        cdata(4+ii) = log(gmsh_ismlt(igm+2*(i+ii-1)))
                      end if
                    end do

                    call autofunc_lagrange(ilmlt(m),cdata,coefficients)

                 end if

C              Lagrange-3 or -4 (center or right)
                 if( ( i .ne. 1 .and. ilmlt(m) .eq. 4 ) .or.
     &               ( i .eq. igr-1 .and. ilmlt(m) .gt. 2 ) ) then

                   write(iot,'( "h: v=[",1pe11.4,",",1pe11.4,",100]")',
     &                          advance='no') cdata(2), cdata(3)
C              Lagrange-2 or -3
                 else

                   write(iot,'( "h: v=[",1pe11.4,",",1pe11.4,",100]")',
     &                          advance='no') cdata(1), cdata(2)

                 end if

                 if( iimlt(m) .eq.  1 .or.
     &               iimlt(m) .eq.  3 .or. iimlt(m) .eq.  4 ) then

                   write(iot,'( " x=[1.*v] ")', advance='no')

                 else

                   write(iot,'( " x=[exp(v)] ")', advance='no')

                 end if

                 if( iimlt(m) .eq.  1 .or.
     &               iimlt(m) .eq. -3 .or. iimlt(m) .eq. -4 ) then

                   if( ilmlt(m) .eq. 2 ) then
                     write(iot,'( "y=[(",1p1e15.8,")*v+(",1p1e15.8,
     &                     ")],l0BZ")')
     &                     coefficients(1), coefficients(2)
                   else if( ilmlt(m) .eq. 3 .or.
     &                     i .eq. 1 .or. i .eq. igr-1 ) then
                     write(iot,'( "y=[(",1p1e15.8,")*v**2+(",1p1e15.8,
     &                     ")*v+(",1p1e15.8,")],l0BZ")')
     &                     coefficients(1), coefficients(2),
     &                     coefficients(3)
                   else if( ilmlt(m) .eq. 4 ) then
                     write(iot,'( "y=[(",1p1e15.8,")*v**3+(",1p1e15.8,
     &                     ")*v**2+(",1p1e15.8,")*v+(",1p1e15.8,
     &                     ")],l0BZ")')
     &                     coefficients(1), coefficients(2),
     &                     coefficients(3), coefficients(4)
                   end if

                else

                   if( ilmlt(m) .eq. 2 ) then
                     write(iot,'( "y=[exp((",1p1e15.8,")*v+(",
     &                     1p1e15.8,"))],l0BZ")') 
     &                     coefficients(1), coefficients(2)
                   else if( ilmlt(m) .eq. 3 .or. 
     &                      i .eq. 1 .or. i .eq. igr-1 ) then
                     write(iot,'( "y=[exp((",1p1e15.8,")*v**2+(",
     &                     1p1e15.8,")*v+(",1p1e15.8,"))],l0BZ")') 
     &                     coefficients(1), coefficients(2),
     &                     coefficients(3)
                   else if( ilmlt(m) .eq. 4 ) then
                     write(iot,'( "y=[exp((",1p1e15.8,")*v**3+(",
     &                     1p1e15.8,")*v**2+(",1p1e15.8,")*v+(",
     &                     1p1e15.8,"))],l0BZ")') 
     &                     coefficients(1), coefficients(2),
     &                     coefficients(3), coefficients(4)
                   end if

                 end if

               end do

            else
            end if

            write(iot,'(/"[END]")')

            close(iot)

         if( impeps(m) .ne. 0 ) then

            iott = 31
            open(iott, file = fname, status = 'unknown' )
            call a_angel(iot,fname)

         end if

C           end do    ! for part or axis

      end if

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine autofunc_lagrange(iclass,cdata,coefficients)
*                                                                      *
*       make a two-four points Lagrange interpolation function         *
*       last modified by N.Matsuda on 2024/12/09                       *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      integer iclass
      dimension cdata(8), coefficients(4)
      real*8  cdata, coefficients
      real*8  coeffi, denomtr, numertr, sub_numertr

*-----------------------------------------------------------------------
*     Lagrange: 
*-----------------------------------------------------------------------

         do 200 k = 1, iclass    ! Linear to Cubic functions

            coeffi = 0.d0
            do 201 m = 1, iclass

               if( ( k .eq. 2 .and. iclass .gt. 2 ) .or.
     &             ( k .eq. 3 .and. iclass .gt. 3 ) ) then
                  numertr = 0.d0
               else
                  numertr = 1.d0
               end if

               denomtr = 1.d0
               do 202 n = 1, iclass

                  if( n .ne. m ) then

                     denomtr = denomtr * ( cdata(m) - cdata(n) )
                     if( k .eq. iclass ) then
                        numertr = numertr * cdata(n)
                     else if( k .eq. 2 ) then
                        numertr = numertr + cdata(n)
                     else if( k .eq. 3 .and. iclass .gt. 3 ) then
                        sub_numertr = 1.d0
                        do 203 i = 1, iclass
                           if( i .ne. m .and. i .ne. n ) then
                              sub_numertr = sub_numertr * cdata(i)
                           end if
  203                   continue
                        numertr = numertr + sub_numertr
                     end if

                  end if

  202          continue

               coeffi = coeffi + cdata(4+m) * ( numertr / denomtr )

  201       continue

            if( mod(k,2) .eq. 1 ) then
               coefficients(k) =        coeffi
            else
               coefficients(k) = -1.0 * coeffi
            end if

  200    continue

*-----------------------------------------------------------------------
      return
      end

