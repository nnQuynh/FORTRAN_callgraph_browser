************************************************************************
*                                                                      *
      subroutine jbook1(id,title,
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to initialize one-dimentional histogram                 *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              id       : histogram id =< 100                          *
*              title    : title of the histgram ( character )          *
*                                                                      *
*              tfac     : total factor                                 *
*                                                                      *
*              ilog     : 0=> x-linear bin, 1=> x-log bin              *
*                                                                      *
*              inum     : =0 ; without event number                    *
*                         =1 ; with N event number                     *
*                         =2 ; with 1/sqrt(N)                          *
*                         =3 ; with Y/sqrt(N)                          *
*                                                                      *
*              ifac     : =0 ; without implicit scaling factor         *
*                         =1 ; with scaling factor 1/dx 1/dw           *
*                                                                      *
*              nx       : number of x-bins                             *
*              xmin     : minimum x-value                              *
*              xmax     : maximum x-value                              *
*              ipx      : x position for print                         *
*                         =1 ; small x position                        *
*                         =2 ; large x position                        *
*                         =3 ; middle x postion                        *
*                                                                      *
*              nw       : number of windows =< 5, or 10                *
*              wi       : edge values of window, wi(10)                *
*              wt       : edge values of window, wt(10)                *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision(a-h, o-z)

      include 'err.inc'

*-----------------------------------------------------------------------

      parameter ( idmax = 90000 )

*-----------------------------------------------------------------------

      common /startf/ iday0,imon0,iyer0,ihor0,imin0,isec0

*-----------------------------------------------------------------------

      dimension vmat(idmax)

      dimension wi(20), wt(20)

      dimension ind(100), ilg(100), inn(100), ifc(100)
      dimension inx(100), inw(100), ixp(100), iqx(100)
      dimension bmi(100), bma(100), bin(100)
      dimension wmi(100,20), wma(100,20)
      dimension gfc(100)

      character title*(*)
      character vtitle*80
      character stitle(100)*80
      dimension ititle(100)

      character br(20)*1

      dimension pxsdd(40)

*-----------------------------------------------------------------------

      data vmat /idmax*0.0/
      save vmat

      data ind / 100*0 /
      save ind, ilg, inn, ifc
      save inw, inx, ixp, iqx
      save bmi, bma, bin
      save wmi, wma

      data gfc / 100*1.0 /
      save gfc

      data indx / 0 /
      save indx

      data ixps / 0 /
      save ixps

      save stitle, ititle

      data br /20*' '/

      common /qmdscm/ qmdscm_h0, qmdscm_d, qmdscm_rcls, iqmdscm

*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------

            if( id .le. 0 .or. id .gt. 100 ) then

           write(ErrCha,*) ' **** Error at jbook1: id is out of range,',
     &                    ' 0 < id < 100'
           ErrID = 'L:114/R:jbook1/F:jbook.f' !E82_001_001
           call ErrWrite(ErrID,ErrCha)
               call parastop( 501 )

            end if

            if( ilog .lt. 0 .or. ilog .gt. 1 ) then

               ilog = 0

            end if

            if( inum .lt. 0 .or. inum .gt. 3 ) then

               inum = 0

            end if

            if( ifac .lt. 0 .or. ifac .gt. 1 ) then

               ifac = 0

            end if

            if( xmin .gt. xmax ) then

               write(ErrCha,*) ' **** Error at jbook1: xmin and xmax ',
     &                    ' are wrong order'
           ErrID = 'L:142/R:jbook1/F:jbook.f' !E82_002_001
           call ErrWrite(ErrID,ErrCha)
               call parastop( 502 )

            end if

            if( ipx .le. 0 .or. ipx  .ge. 4 ) then

               write(ErrCha,*) ' **** Error at jbook1: ipx ',
     &                    ' is wrong.'
           ErrID = 'L:152/R:jbook1/F:jbook.f' !E82_003_001
           call ErrWrite(ErrID,ErrCha)
               call parastop( 503 )

            end if

            if( ilog .eq. 1 .and. xmin .le. 0.0 ) then

           write(ErrCha,*) ' **** Error at jbook1: x-range is wrong,',
     &                    ' as x-bin is log'
           ErrID = 'L:162/R:jbook1/F:jbook.f' !E82_004_001
           call ErrWrite(ErrID,ErrCha)
               call parastop( 504 )

            end if

            if( nw .lt. 0 ) nw = 0

            if( inum .eq. 0 .and. nw .gt. 20 ) then

           write(ErrCha,*) ' **** Error at jbook1: nw should be less ',
     &                    ' than 20, when inum = 0'
           ErrID = 'L:174/R:jbook1/F:jbook.f' !E82_005_001
           call ErrWrite(ErrID,ErrCha)
               call parastop( 505 )

            end if

            if( inum .ne. 0 .and. nw .gt. 10 ) then

           write(ErrCha,*) ' **** Error at jbook1: nw should be less ',
     &                    ' than 10, when inum > 0'
           ErrID = 'L:184/R:jbook1/F:jbook.f' !E82_006_001
           call ErrWrite(ErrID,ErrCha)
               call parastop( 506 )

            end if

            if( ind(id) .ne. 0 ) then

               write(ErrCha,*) ' **** Error at jbook1: double booking'
           ErrID = 'L:193/R:jbook1/F:jbook.f' !E82_007_001
           call ErrWrite(ErrID,ErrCha)
               call parastop( 507 )

            end if

*-----------------------------------------------------------------------
*     initialization of booking
*-----------------------------------------------------------------------
*        store configuration
*-----------------------------------------------------------------------

               indx = indx + 1

               ind(id)   = indx

               ilg(indx) = ilog
               inn(indx) = inum
               ifc(indx) = ifac
               inx(indx) = nx
               inw(indx) = nw
               iqx(indx) = ipx

               ixp(indx) = ixps

               gfc(indx) = gfc(indx) * tfac

*-----------------------------------------------------------------------
*        store bin range
*-----------------------------------------------------------------------

               if( ilog .eq. 0 ) then

                     bmi(indx) = xmin
                     bma(indx) = xmax

               else

                     bmi(indx) = log(xmin)
                     bma(indx) = log(xmax)

               end if

                     bin(indx) = ( bma(indx) - bmi(indx) ) / nx

            do i = 1, nx + 1

               if( ilog .eq. 0 ) then

                     vmat(ixps+i) = xmin + bin(indx) * ( i - 1 )

               else

                     vmat(ixps+i) = xmin * exp( bin(indx) * ( i - 1 ) )

               end if

            end do

*-----------------------------------------------------------------------
*        store window range
*-----------------------------------------------------------------------

            if( nw .gt. 0 ) then

               do i = 1, nw

                  if( wi(i) .lt. wt(i) ) then

                     wmi(indx,i) = wi(i)
                     wma(indx,i) = wt(i)

                  else if( wi(i) .gt. wt(i) ) then

                     wmi(indx,i) = wt(i)
                     wma(indx,i) = wi(i)

                  else if( wi(i) .eq. wt(i) ) then

                     write(ErrCha,*) ' **** Error at jbook1:',
     &                          ' window size is zero'
           ErrID = 'L:274/R:jbook1/F:jbook.f' !E82_008_001
           call ErrWrite(ErrID,ErrCha)
                     call parastop( 508 )

                  end if

               end do

            end if

*-----------------------------------------------------------------------
*        store column position
*-----------------------------------------------------------------------

               mx = nx + 1 + 4

            if( inum .eq. 0 ) then

               if( nw .eq. 0 ) then

                     ixps = ixps + 2 * mx

               else

                     ixps = ixps + ( nw + 1 ) * mx

               end if

            else

               if( nw .eq. 0 ) then

                     ixps = ixps + 3 * mx

               else

                     ixps = ixps + ( 2 * nw + 1 ) * mx

               end if

            end if


            if( ixps .gt. idmax) then

               write(ErrCha,*) ' **** Error at jbook1: over booking.',
     &                    ' Please increase idmax'
           ErrID = 'L:321/R:jbook1/F:jbook.f' !E82_009_001
           call ErrWrite(ErrID,ErrCha)
               call parastop( 509 )

            end if

*-----------------------------------------------------------------------
*        store title of histgram
*-----------------------------------------------------------------------

         vtitle = title//' '

            do i = 80, 1, -1

               if( vtitle(i:i) .ne. ' ' ) goto 100

            end do

               i = 0

  100       continue

               ititle(indx) = i

            do i = 1, ititle(indx)

               stitle(indx)(i:i) = vtitle(i:i)

            end do

            do i = 1, 80

               vtitle(i:i) = ' '

            end do

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------


************************************************************************
*                                                                      *
      entry jbkreset
*                                                                      *
*                                                                      *
*        Last Revised:     2001 03 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              reset the booking                                       *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

         indx = 0
         ixps = 0

      do i = 1, idmax

         vmat(i) = 0.0

      end do

      do i = 1, 100

         ind(i) = 0
         gfc(i) = 1.0

      end do

      do i = 1, 20

         br(i) = ' '

      end do

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------


************************************************************************
*                                                                      *
      entry jfill1(id,x,v,wn,wf)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to fill one dimensional histogram                       *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              id       : histogram id =< 100                          *
*              x        : x-value                                      *
*              v        : window-value                                 *
*              wn       : number weight                                *
*              wf       : factor                                       *
*                                                                      *
*                         total weight = wn * wf                       *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------

            ic = ind(id)

            if( ic .eq. 0 ) then

               write(ErrCha,*) ' **** Error at jfill1: unbooked id'
           ErrID = 'L:440/R:jbook1/F:jbook.f' !E82_010_001
           call ErrWrite(ErrID,ErrCha)
               call parastop( 510 )

            end if

            if( ilg(ic) .eq. 1 .and. x .le. 0.0 ) return

*-----------------------------------------------------------------------
*        fill matrix
*-----------------------------------------------------------------------

                        ix = ixp(ic)
                        mx = inx(ic) + 1
                        lx = mx + 4


                     if( ilg(ic) .eq. 1 ) then

                        xx = log( x )

                     else

                        xx = x

                     end if

*-----------------------------------------------------------------------

                  if( iqmdscm .ne. 1 ) then

                        vmat(ix+mx+3) = vmat(ix+mx+3) + wn * wf
                        vmat(ix+mx+4) = vmat(ix+mx+4) + wn

                     if( xx .lt. bmi(ic) ) then

                        vmat(ix+mx+1) = vmat(ix+mx+1) + wn * wf

                     else if( xx .ge. bma(ic) ) then

                        vmat(ix+mx+2) = vmat(ix+mx+2) + wn * wf

                     end if

                  else

                     if( xx .lt. bmi(ic) ) then

                        vmat(ix+mx+1) = vmat(ix+mx+1) + wn * wf

                     else if( xx .ge. bma(ic) ) then

                        vmat(ix+mx+2) = vmat(ix+mx+2) + wn * wf

                     else

                        vmat(ix+mx+3) = vmat(ix+mx+3) + wn * wf
                        vmat(ix+mx+4) = vmat(ix+mx+4) + wn

                     end if

                  endif

*-----------------------------------------------------------------------

            if( inw(ic) .eq. 0 ) then

                  if( xx .lt. bmi(ic) ) then

                        vmat(ix+lx+mx+1) = vmat(ix+lx+mx+1) + wn

                  else if( xx .ge. bma(ic) ) then

                        vmat(ix+lx+mx+2) = vmat(ix+lx+mx+2) + wn

                  else

                        j = int( ( xx - bmi(ic) ) / bin(ic) ) + 1

                     if( ifc(ic) .eq. 1 ) then

                        fac = 1.0 / ( vmat(ix+j+1) - vmat(ix+j) )

                     else

                        fac = 1.0

                     end if

                        vmat(ix+lx+j) = vmat(ix+lx+j) + wn * wf * fac

                     if( inn(ic) .gt. 0 ) then

                        vmat(ix+2*lx+j) = vmat(ix+2*lx+j) + wn

                     end if

                  end if

*-----------------------------------------------------------------------

            else if( inw(ic) .gt. 0 ) then


               do ii = 1, inw(ic)
               if( v .ge. wmi(ic,ii) .and. v .lt. wma(ic,ii) ) then

                     if( inn(ic) .eq. 0 ) then

                        i = ii

                     else

                        i = 2 * ii - 1

                     end if

                     if( ifc(ic) .eq. 1 ) then

                        fac = 1.0 / ( wma(ic,ii) - wmi(ic,ii) )

                     else

                        fac = 1.0

                     end if

                        vmat(ix+i*lx+mx+3) = vmat(ix+i*lx+mx+3)
     &                                     + wn * wf * fac
                        vmat(ix+i*lx+mx+4) = vmat(ix+i*lx+mx+4)
     &                                     + wn

                  if( xx .lt. bmi(ic) ) then

                        vmat(ix+i*lx+mx+1) = vmat(ix+i*lx+mx+1)
     &                                     + wn * wf * fac

                     if( inn(ic) .gt. 0 ) then

                        vmat(ix+(i+1)*lx+mx+1) = vmat(ix+(i+1)*lx+mx+1)
     &                                         + wn

                     end if

                  else if( xx .ge. bma(ic) ) then

                        vmat(ix+i*lx+mx+2) = vmat(ix+i*lx+mx+2)
     &                                     + wn * wf * fac

                     if( inn(ic) .gt. 0 ) then

                        vmat(ix+(i+1)*lx+mx+2) = vmat(ix+(i+1)*lx+mx+2)
     &                                         + wn

                     end if

                  else

                        j = int( ( xx - bmi(ic) ) / bin(ic) ) + 1

                     if( ifc(ic) .eq. 1 ) then

                        fac = fac / ( vmat(ix+j+1) - vmat(ix+j) )

                     end if

                        vmat(ix+i*lx+j) = vmat(ix+i*lx+j)
     &                                  + wn * wf * fac

                     if( inn(ic) .gt. 0 ) then

                        vmat(ix+(i+1)*lx+j) = vmat(ix+(i+1)*lx+j)
     &                                      + wn

                     end if

                  end if

               end if
               end do

            end if


*-----------------------------------------------------------------------

      return


************************************************************************
*                                                                      *
      entry jscale1(id,scal)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to scal one dimentional histgram data                   *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              id       : histogram id =< 100                          *
*              scal     : global scaling factor                        *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------

            ic = ind(id)

            if( ic .eq. 0 ) then

               write(ErrCha,*) ' **** Error at jscale1: unbooked id'
           ErrID = 'L:657/R:jbook1/F:jbook.f' !E82_011_001
           call ErrWrite(ErrID,ErrCha)
               return

            end if

*-----------------------------------------------------------------------

            gfc(ic) = gfc(ic) * scal

*-----------------------------------------------------------------------

      return


************************************************************************
*                                                                      *
      entry jprint1(id,io)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to print one dimensional histogram                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              id       : histogram id =< 100                          *
*              io       : output unit                                  *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------

            ic = ind(id)

            if( ic .eq. 0 ) then

               write(ErrCha,*) ' **** Error at jprint1: unbooked id'
           ErrID = 'L:701/R:jbook1/F:jbook.f' !E82_012_001
           call ErrWrite(ErrID,ErrCha)
               return

            end if

*-----------------------------------------------------------------------
*        some constants
*-----------------------------------------------------------------------

                  gfac = gfc(ic)

                  ix = ixp(ic)
                  mx = inx(ic) + 1
                  lx = mx + 4

                  ippx = iqx(ic)
                  ilgx = ilg(ic)

*-----------------------------------------------------------------------
*        write title
*-----------------------------------------------------------------------

               write(io,'("#")')

               write(io,'("#",5x,"JBOOK: ",80a1)')
     &              ( stitle(ic)(i:i), i = 1, ititle(ic) )

               write(io,'("#")')

               write(io,'("#",5x,"Histgram ID :",i3)') id

               write(io,'("#",5x,"No. of entries :",1pg15.7)')
     &               vmat(ix+mx+4)

               write(io,'("#",5x,"Date: ",
     &               i4.4,"-",i2.2,"-",i2.2,"   ",
     &               i2.2,":",i2.2,":",i2.2)')
     &               iyer0,imon0,iday0,ihor0,imin0,isec0

               write(io,'("#")')

*-----------------------------------------------------------------------
*        write header
*-----------------------------------------------------------------------

            if( inw(ic) .eq. 0 ) then

               if( inn(ic) .eq. 0 ) then

                     write(io,'("#",5x,"x  ",
     &                     9x,"y01")')

               else

                     write(io,'("#",5x,"x  ",
     &                     9x,"y01",9x,"n01")')

               end if

            else if( inw(ic) .gt. 0 ) then

               if( inn(ic) .eq. 0 ) then

                     write(io,'("#",5x,"x  ",
     &                     20(a1,8x,"y",i2.2))')
     &                     ( br(j), j, j=1,inw(ic) )

               else

                     write(io,'("#",5x,"x  ",
     &                     20(a1,8x,"y"i2.2,9x,"n",i2.2))')
     &                     ( br(j), j, j, j=1,inw(ic) )

               end if

            end if

               write(io,'("#")')

*-----------------------------------------------------------------------
*        print histgrams
*-----------------------------------------------------------------------

         do k = 1, mx - 1

               if( ippx .eq. 1 ) then

                     xpos = vmat(ix+k)

               else if( ippx .eq. 2 ) then

                     xpos = vmat(ix+k+1)

               else if( ippx .eq. 3 ) then

                  if( ilgx .eq. 0 ) then

                     xpos = ( vmat(ix+k) + vmat(ix+k+1) ) / 2.0

                  else

                     xpos = exp( ( log(vmat(ix+k))
     &                           + log(vmat(ix+k+1)) ) / 2.0 )

                  end if

               end if


            if( inw(ic) .eq. 0 ) then

               if( inn(ic) .eq. 0 ) then

                     write(io,'(1x,1p2g12.4)')
     &               xpos, vmat(ix+lx+k) * gfac

               else

                  if( vmat(ix+2*lx+k) .gt. 0.0 ) then

                     if( inn(ic) .eq. 2 ) then

                        vmat(ix+2*lx+k) = 1.0 / sqrt( vmat(ix+2*lx+k) )

                     else if( inn(ic) .eq. 3 ) then

                        vmat(ix+2*lx+k) = vmat(ix+lx+k) * gfac
     &                                  / sqrt( vmat(ix+2*lx+k) )

                     end if

                  end if

                     write(io,'(1x,1p3g12.4)')
     &               xpos, vmat(ix+lx+k) * gfac, vmat(ix+2*lx+k)

               end if

            else if( inw(ic) .gt. 0 ) then

               if( inn(ic) .eq. 0 ) then

                     write(io,'(1x,1p21g12.4)')
     &               xpos, ( vmat(ix+j*lx+k) * gfac, j=1,inw(ic) )

               else

                  if( inn(ic) .eq. 2 .or. inn(ic) .eq. 3 ) then

                     do l = 2, inw(ic)*2, 2

                        if( vmat(ix+l*lx+k) .gt. 0.0 ) then

                           if( inn(ic) .eq. 2 ) then

                              vmat(ix+l*lx+k) = 1.0
     &                                        / sqrt( vmat(ix+l*lx+k) )

                           else if( inn(ic) .eq. 3 ) then

                              vmat(ix+l*lx+k) = vmat(ix+(l-1)*lx+k)
     &                                        * gfac
     &                                        / sqrt( vmat(ix+l*lx+k) )

                           end if

                        end if

                     end do

                  end if

                     write(io,'(1x,1p21g12.4)')
     &               xpos, ( vmat(ix+(2*j-1)*lx+k) * gfac,
     &                             vmat(ix+2*j*lx+k), j=1,inw(ic) )

               end if

            end if

         end do

*-----------------------------------------------------------------------
*        write summary
*-----------------------------------------------------------------------

                     write(io,'("#")')

            if( inw(ic) .eq. 0 ) then

                     write(io,'("#",4x,"total  =",1p1g12.4,
     &                                        "# =",1p1g15.7)')
     &                     vmat(ix+mx+3) * gfac, vmat(ix+mx+4)

                     write(io,'("#",4x,"under  =",1p1g12.4,
     &                                        "# =",1p1g15.7)')
     &                     vmat(ix+mx+1) * gfac, vmat(ix+lx+mx+1)

                     write(io,'("#",4x,"over   =",1p1g12.4,
     &                                        "# =",1p1g15.7)')
     &                     vmat(ix+mx+2) * gfac, vmat(ix+lx+mx+2)

            else if( inw(ic) .gt. 0 ) then

                     write(io,'("#",4x,"total  =",1p1g12.4,
     &                                        "# =",1p1g15.7)')
     &                     vmat(ix+mx+3) * gfac, vmat(ix+mx+4)

                     write(io,'("#",4x,"under  =",1p1g12.4)')
     &                     vmat(ix+mx+1) * gfac

                     write(io,'("#",4x,"over   =",1p1g12.4)')
     &                     vmat(ix+mx+2) * gfac

                     write(io,'("#")')

               if( inn(ic) .eq. 0 ) then

                     write(io,'("#",4x,"win:",
     &                     20(a1,8x,"y",i2.2))')
     &                     ( br(j), j, j=1,inw(ic) )

                     write(io,'("#")')

                     write(io,'("#",4x,"total  =",1p21g12.4)')
     &                     ( vmat(ix+j*lx+mx+3) * gfac, j=1,inw(ic) )

                     write(io,'("#",4x,"    #  =",1p21g12.4)')
     &                     ( vmat(ix+j*lx+mx+4), j=1,inw(ic) )

                     write(io,'("#",4x,"under  =",1p21g12.4)')
     &                     ( vmat(ix+j*lx+mx+1) * gfac, j=1,inw(ic) )

                     write(io,'("#",4x,"over   =",1p21g12.4)')
     &                     ( vmat(ix+j*lx+mx+2) * gfac, j=1,inw(ic) )

               else

                     write(io,'("#",4x,"win:",
     &                     20(a1,8x,"y"i2.2,9x,"#",i2.2))')
     &                     ( br(j), j, j, j=1,inw(ic) )

                     write(io,'("#")')

                     write(io,'("#",4x,"total  =",1p21g12.4)')
     &                     ( vmat(ix+(2*j-1)*lx+mx+3) * gfac,
     &                       vmat(ix+(2*j-1)*lx+mx+4), j=1,inw(ic) )

                     write(io,'("#",4x,"under  =",1p21g12.4)')
     &                     ( vmat(ix+(2*j-1)*lx+mx+1) * gfac,
     &                       vmat(ix+(2*j  )*lx+mx+1), j=1,inw(ic) )

                     write(io,'("#",4x,"over   =",1p21g12.4)')
     &                     ( vmat(ix+(2*j-1)*lx+mx+2) * gfac,
     &                       vmat(ix+(2*j  )*lx+mx+2), j=1,inw(ic) )

               end if

            end if

               write(io,'("#")')


*-----------------------------------------------------------------------

      return


************************************************************************
*                                                                      *
      entry jftot(id,pxst,pxsu,pxso,pfac)
*                                                                      *
*                                                                      *
*        Last Revised:     1999 03 15                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to give total particle production cross section         *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              id       : histogram id =< 100                          *
*              pxst     : total cross section within the range         *
*              pxsu     : total cross section under  the range         *
*              pxso     : total cross section over   the range         *
*              pfac     : normalization factor                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------

            ic = ind(id)

            if( ic .eq. 0 ) then

               write(ErrCha,*) ' **** Error at jftot: unbooked id'
           ErrID = 'L:1000/R:jbook1/F:jbook.f' !E82_013_001
           call ErrWrite(ErrID,ErrCha)
               return

            end if

*-----------------------------------------------------------------------
*        some constants
*-----------------------------------------------------------------------

                  gfac = gfc(ic) * pfac

                  ix = ixp(ic)
                  mx = inx(ic) + 1
                  lx = mx + 4

*-----------------------------------------------------------------------
*        write summary
*-----------------------------------------------------------------------

                  pxst = vmat(ix+mx+3) * gfac
                  pxsu = vmat(ix+mx+1) * gfac
                  pxso = vmat(ix+mx+2) * gfac

                  pxst = pxst - pxsu - pxso

*-----------------------------------------------------------------------

      return


************************************************************************
*                                                                      *
      entry jfddx(id,jd,pxsdd)
*                                                                      *
*                                                                      *
*        Last Revised:     1999 03 16                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to give the double differential cross section           *
*              for one angel                                           *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              id       : histogram id =< 100                          *
*              jd       : angle id                                     *
*              pxsdd    : cross section                                *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------

            ic = ind(id)

         if( ic .ne. 0 ) then

*-----------------------------------------------------------------------
*        some constants
*-----------------------------------------------------------------------

                  gfac = gfc(ic)

                  ix = ixp(ic)
                  mx = inx(ic) + 1
                  lx = mx + 4

*-----------------------------------------------------------------------
*        give ddx
*-----------------------------------------------------------------------

               j = jd

               pxsdd(1)   = vmat(ix+j*lx+mx+1) * gfac
     &                    / vmat(ix+1)

            do k = 1, mx - 1

               pxsdd(k+1) = vmat(ix+j*lx+k) * gfac

            end do

               pxsdd(mx+1)  = 0.0

*-----------------------------------------------------------------------

         else

            do k = 1, 40

               pxsdd(k) = 0.0

            end do

         end if

*-----------------------------------------------------------------------

      return
      end



