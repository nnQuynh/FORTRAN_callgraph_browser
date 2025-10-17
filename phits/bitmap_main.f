
************************************************************************
*                                                                      *
      subroutine bitmap_create_filename(
     &                fname,
     &                bmpfIType, bmpfIndex, numIType,
     &                fileIndex, sfname)
*                                                                      *
*       create bitmap filename.                                        *
*                                                                      *
************************************************************************

        implicit none
        character(1) :: fname(100)
        character(1) :: bmpfIType(*)
        integer :: bmpfIndex(*)
        integer :: numIType
        integer :: fileIndex
        character(len=*) :: sfname

*-----------------------------------------------------------------------

        integer :: isize, ipos, maxdig, newisize, index
        integer :: i
        character(len=20) :: sgid
        character(len=10) :: format
        integer :: tid(numIType), gid(numIType)

        integer :: get_line_length
        integer :: last_char_index
        integer :: max_digits_integers

*-----------------------------------------------------------------------

        call clear_string(sfname)
        isize = get_line_length(fname, 100)
        call char2string(isize, fname, sfname)
        ipos = last_char_index(sfname, '.')
        isize = ipos-1

        maxdig = max_digits_integers(bmpfIndex, numIType)
        write(format,'(''(I'',1I0,''.'',1I0,'')'')') maxdig, maxdig

        index = fileIndex
        do i = numIType, 2, -1
          tid(i) = product(bmpfIndex(1:i-1))
          gid(i) = floor((index-1) / tid(i) * 1.0d0) + 1
          index = index - (gid(i)-1)*tid(i)
        end do
        gid(1) = index


        do i = numIType, 1, -1       ! from outer do-loop to inner
          if ( bmpfIndex(i) .gt. 1 ) then
            write(sgid,format) gid(i)
            newisize = isize + 2 + len_trim(sgid)
            sfname(isize+1:newisize) = '_'//bmpfIType(i)//trim(sgid)
            isize = newisize
          end if
        end do

        sfname(isize+1:isize+4) = '.bmp'
        isize = isize + 4
      end subroutine



************************************************************************
*                                                                      *
      subroutine clear_bitmap_bgrQuad(bgrQuad, isize)
*                                                                      *
*       fill bitmap bgr quadple with white color.                     *
*                                                                      *
************************************************************************

        implicit none
        character(1) :: bgrQuad(*)
        integer :: isize

*-----------------------------------------------------------------------

        integer :: i, ipos

*-----------------------------------------------------------------------

        do i = 1, isize
          ipos = 4*(i-1)+1
          bgrQuad(ipos+0) = achar(255)
          bgrQuad(ipos+1) = achar(255)
          bgrQuad(ipos+2) = achar(255)
          bgrQuad(ipos+3) = achar(0)
        end do
      end subroutine



************************************************************************
*                                                                      *
      subroutine set_bitmap_bgrQuad_hsv(bgrQuad, rcol, ipos)
*                                                                      *
*       set rcol value (hsv color space) to pixel of rgb bitmap.       *
*                                                                      *
************************************************************************

        implicit none
        character(1) :: bgrQuad(*)
        real(8) :: rcol(3)
        integer :: ipos

*-----------------------------------------------------------------------

        integer :: b, g, r
        real(8) :: h, s, v
        real(8) :: rgb(3)

*-----------------------------------------------------------------------

        h = rcol(1)-1.0d0          ! hsv color of angel data has minus 1.
        s = rcol(2)
        v = rcol(3)

        call hsv2rgb(h, s, v, rgb)

        r = floor(rgb(1) * 255.0d0)
        g = floor(rgb(2) * 255.0d0)
        b = floor(rgb(3) * 255.0d0)
        bgrQuad(ipos+0) = achar(b)
        bgrQuad(ipos+1) = achar(g)
        bgrQuad(ipos+2) = achar(r)
        bgrQuad(ipos+3) = achar(0)
      end subroutine



************************************************************************
*                                                                      *
      subroutine set_bitmap_bgrQuad_rgb(bgrQuad, rcol, ipos)
*                                                                      *
*       set rcol value (rgb color space) to pixel of rgb bitmap.       *
*                                                                      *
************************************************************************

        implicit none
        character(1) :: bgrQuad(*)
        real(8) :: rcol(3)
        integer :: ipos

*-----------------------------------------------------------------------

        integer :: b, g, r
        real(8) :: h, s, v
        real(8) :: rgb(3)

*-----------------------------------------------------------------------

        r = rcol(1)-1.0d0          ! hsv color of angel data has minus 1.
        g = rcol(2)
        b = rcol(3)

        r = floor(r * 255.0d0)
        g = floor(g * 255.0d0)
        b = floor(b * 255.0d0)
        bgrQuad(ipos+0) = achar(b)
        bgrQuad(ipos+1) = achar(g)
        bgrQuad(ipos+2) = achar(r)
        bgrQuad(ipos+3) = achar(0)
      end subroutine



************************************************************************
*                                                                      *
      subroutine set_bitmap_bgrQuad_gray(bgrQuad, rcol, ipos, ratio)
*                                                                      *
*       set rcol value (gray scale) to pixel of rgb bitmap.            *
*                                                                      *
*       ratio is not used because it is width ration of point size     *
*       on [hd:] angel cluster plot.                                   *
*                                                                      *
************************************************************************

        implicit none
        character(1) :: bgrQuad(*)
        real(8) :: rcol(3)
        integer :: ipos
        real(8) :: ratio

*-----------------------------------------------------------------------

        integer :: b, g, r
        real(8) :: gray
        real(8) :: rgb(3)

*-----------------------------------------------------------------------

        gray = rcol(1)+2.0d0

        r = floor(gray * 255.0d0)
        g = floor(gray * 255.0d0)
        b = floor(gray * 255.0d0)
        bgrQuad(ipos+0) = achar(b)
        bgrQuad(ipos+1) = achar(g)
        bgrQuad(ipos+2) = achar(r)
        bgrQuad(ipos+3) = achar(0)
      end subroutine



************************************************************************
*                                                                      *
      subroutine set_bitmap_bgrQuad(bgrQuad, rcol, ipos)
*                                                                      *
*       if rcol(1)<0 then gray scale(+2), else hsv(-1)                 *
*                                                                      *
************************************************************************

        implicit none
        character(1) :: bgrQuad(*)
        real(8) :: rcol(3)
        integer :: ipos

*-----------------------------------------------------------------------

        if ( rcol(1).lt.0.0d0 ) then
          call set_bitmap_bgrQuad_gray(bgrQuad, rcol, ipos, 1.0d0)
        else
          call set_bitmap_bgrQuad_hsv(bgrQuad, rcol, ipos)
        end if
      end subroutine



************************************************************************
*                                                                      *
      subroutine bitmap_main(
     &                iubmp_default, filename,
     &                width, height, bgrQuad, rgbaSize)
*                                                                      *
*       create bitmap file and write rgb values.                       *
*                                                                      *
************************************************************************

        implicit none
        integer :: iubmp_default
        character(len=*) :: filename
        integer :: width, height
        integer :: rgbaSize
        character(1) :: bgrQuad(rgbaSize)

*-----------------------------------------------------------------------

        integer :: iubmp
        logical :: isText = .false.
        integer :: ios

*-----------------------------------------------------------------------

        call open_file(iubmp_default, filename, iubmp, ios, isText)

        call write_bitmap_header(iubmp, width, height)
        write(iubmp) bgrQuad

        call close_file(iubmp)
      end subroutine



************************************************************************
*                                                                      *
      subroutine bitmap_fill(
     &                jbt, npoints, rcol,
     &                bmpWidth, bmpHeight, bgrQuad, rgbaSize,
     &                filling)
*                                                                      *
*       fill or draw polygon which points readed from scratch.         *
*                                                                      *
************************************************************************

        implicit none
        integer :: jbt
        integer :: npoints
        real(8) :: rcol(3)
        integer :: bmpWidth
        integer :: bmpHeight
        integer :: rgbaSize
        character(1) :: bgrQuad(rgbaSize)
        logical :: filling

*-----------------------------------------------------------------------

        integer, allocatable :: ix(:), iy(:)
        integer, allocatable :: ipxy(:)
        integer :: i
        integer :: ixo, iyo, ipos

*-----------------------------------------------------------------------

        allocate(ix(npoints))
        allocate(iy(npoints))
        allocate(ipxy(npoints))

        do i = 1, npoints
          read(jbt) ix(i), iy(i)
        end do
        do i = 1, npoints
          ipxy(i) = i
        end do


        if ( filling ) then

          call bitmap_fill_pixel(
     &              bmpWidth, bmpHeight, npoints,
     &              ix, iy, bgrQuad, rgbaSize, rcol)

        else

          call bitmap_draw_pixel(
     &            bmpWidth, bmpHeight, npoints,
     &            ix, iy, bgrQuad, rgbaSize, rcol)

        end if

        deallocate(ipxy)
        deallocate(iy)
        deallocate(ix)
      end subroutine



************************************************************************
*                                                                      *
      subroutine bitmap_fill_pixel(
     &                bmpWidth, bmpHeight, npoints,
     &                ix, iy, bgrQuad, rgbaSize, rcol)
*                                                                      *
*       fill pixels in polygon.                                        *
*                                                                      *
************************************************************************

        implicit none
        integer :: bmpWidth
        integer :: bmpHeight
        integer :: npoints
        integer :: ix(npoints), iy(npoints)
        integer :: rgbaSize
        character(1) :: bgrQuad(rgbaSize)
        real(8) :: rcol(3)

*-----------------------------------------------------------------------

        logical :: inside
        integer :: i, j
        integer :: ixo, iyo, ipos
        integer :: ixw(npoints), iyw(npoints), npw
        integer :: ixw2(npoints), iyw2(npoints), npw2
        integer :: pbuf(bmpWidth*bmpHeight), ip
        integer :: ncall

*-----------------------------------------------------------------------

        call bitmap_fill_remove_same_points(
     &          npoints, ix, iy,
     &          npw, ixw, iyw)
        call bitmap_fill_deduce_points(
     &          npw, ixw, iyw,
     &          npw2, ixw2, iyw2)


        do i = 1, bmpWidth*bmpHeight
          pbuf(i) = 0
        end do


        call bitmap_draw_set_paintbuffer(
     &          bmpWidth, bmpHeight, npw, ixw, iyw, pbuf)

        call bitmap_scanline_set_paintbuffer(
     &          bmpWidth, bmpHeight,
     &          npw, ixw, iyw,
     &          npw2, ixw2, iyw2, pbuf, ncall)


        do j = 1, bmpHeight
          do i = 1, bmpWidth
            ip = (j-1)*bmpWidth + i
            if ( pbuf(ip).gt.0 ) then
              ipos = (j-1)*bmpWidth*4 + (i-1)*4 + 1
              call set_bitmap_bgrQuad(
     &                bgrQuad, rcol, ipos)
            end if
          end do
        end do

      end subroutine



************************************************************************
*                                                                      *
      subroutine bitmap_get_range_polygon(
     &                bmpWidth, bmpHeight,
     &                npoints, ix, iy, ixrange, iyrange)
*                                                                      *
*       return xy range of polygon.                                    *
*                                                                      *
************************************************************************

        implicit none
        integer :: bmpWidth, bmpHeight
        integer :: npoints
        integer :: ix(npoints), iy(npoints)
        integer :: ixrange(2), iyrange(2)

*-----------------------------------------------------------------------

        integer :: ixmin, iymin, ixmax, iymax, i, ix1, iy1

*-----------------------------------------------------------------------

        ixmin = bmpWidth + 1
        iymin = bmpHeight + 1
        ixmax = 0
        iymax = 0
        do i = 1, npoints
          ix1 = ix(i)
          iy1 = iy(i)
          if ( ix1.lt.ixmin ) ixmin = ix1
          if ( iy1.lt.iymin ) iymin = iy1
          if ( ix1.gt.ixmax ) ixmax = ix1
          if ( iy1.gt.iymax ) iymax = iy1
        end do

        ixrange(1) = ixmin
        ixrange(2) = ixmax
        iyrange(1) = iymin
        iyrange(2) = iymax
      end subroutine



************************************************************************
*                                                                      *
      subroutine bitmap_fill_point_is_inside_winding(
     &                npoints, ix, iy,
     &                ixo, iyo,
     &                inside)
*                                                                      *
*       check whether (ixo,iyo) is inside of polygon.                  *
*                                                                      *
************************************************************************

        implicit none
        integer :: npoints
        integer :: ix(npoints), iy(npoints)
        integer :: ixo, iyo
        logical :: inside

*-----------------------------------------------------------------------

        integer :: ix1, iy1, ix2, iy2
        double precision :: a1, a2, b1, b2, wn, cos, theta, dot
        integer :: i

*-----------------------------------------------------------------------

        ix2 = ix(1)
        iy2 = iy(1)
        if ( ixo.eq.ix2 .and. iyo.eq.iy2 ) then
          inside = .true.
          return
        end if

        wn = 0.0d0

        do i = 1, npoints
          if (i.ne.npoints) then
            ix1 = ix(i+1)
            iy1 = iy(i+1)
          else
            ix1 = ix(1)
            iy1 = iy(1)
          end if

          if ( ix1.eq.ix2 .and. iy1.eq.iy2 ) cycle

          if ( ixo.eq.ix1 .and. iyo.eq.iy1 ) then
            inside = .true.
            return
          end if

          a1 = (ix1-ixo)*1.0d0
          a2 = (iy1-iyo)*1.0d0
          b1 = (ix2-ixo)*1.0d0
          b2 = (iy2-iyo)*1.0d0
          cos = (a1*b1 + a2*b2) / sqrt(a1*a1+a2*a2) / sqrt(b1*b1+b2*b2)
          dot = a1*b2 - a2*b1
          if (cos.gt.1.0d0) cos = 1.0d0
          if (cos.lt.-1.0d0) cos = -1.0d0
          theta = acos(cos)

          wn = wn + sign(1.0d0,dot) * theta

          ix2 = ix1
          iy2 = iy1

        end do

        wn = wn/(2.0d0*3.1416d0)

        inside = .true.
        if ( abs(wn) .lt. 0.1d0 ) then
          inside = .false.
        end if
      end subroutine



************************************************************************
*                                                                      *
      subroutine bitmap_fill_remove_same_points(
     &                npoints, ix, iy, npw, ixw, iyw)
*                                                                      *
*       remove points of same pixel coordinates in polygon.            *
*                                                                      *
************************************************************************

        implicit none
        integer :: npoints
        integer :: ipxy(npoints)
        integer :: ix(npoints), iy(npoints)
        integer :: ixw(npoints), iyw(npoints), npw

*-----------------------------------------------------------------------

        integer :: ix0, iy0, ix1, iy1
        integer :: i

*-----------------------------------------------------------------------

        npw = 1
        ixw(1) = ix(1)
        iyw(1) = iy(1)

        do i = 2, npoints
          ix0 = ix(i-1)
          iy0 = iy(i-1)
          ix1 = ix(i)
          iy1 = iy(i)

          if ( .not. ( ix1.eq.ix0 .and. iy1.eq.iy0 ) ) then
            npw = npw + 1
            ixw(npw) = ix1
            iyw(npw) = iy1
          end if
        end do
      end subroutine



************************************************************************
*                                                                      *
      subroutine bitmap_fill_deduce_points(
     &                npoints, ix, iy, npw, ixw, iyw)
*                                                                      *
*       deduce points on linear edges of polygon boundary.             *
*                                                                      *
************************************************************************

        implicit none
        integer :: npoints
        integer :: ix(npoints), iy(npoints)
        integer :: ixw(npoints), iyw(npoints), npw

*-----------------------------------------------------------------------

        integer :: ix0, iy0, ix1, iy1, ix2, iy2
        integer :: i

*-----------------------------------------------------------------------

        npw = 1
        ixw(1) = ix(1)
        iyw(1) = iy(1)

        do i = 2, npoints-1
          ix0 = ix(i-1)
          iy0 = iy(i-1)
          ix1 = ix(i)
          iy1 = iy(i)
          ix2 = ix(i+1)
          iy2 = iy(i+1)

          if ( .not. ( ix1.eq.ix0 .and. ix1.eq.ix2 ) .and.
     &         .not. ( iy1.eq.iy0 .and. iy1.eq.iy2 ) ) then
            npw = npw + 1
            ixw(npw) = ix1
            iyw(npw) = iy1
          end if
        end do

        ix1 = ix(npoints)
        iy1 = iy(npoints)
        npw = npw + 1
        ixw(npw) = ix1
        iyw(npw) = iy1
      end subroutine



************************************************************************
*                                                                      *
      subroutine bitmap_scanline_set_paintbuffer(
     &                bmpWidth, bmpHeight,
     &                np, ix, iy,
     &                npw, ixw, iyw, pbuf, ncall)
*                                                                      *
*       fill paint buffer within xy-range of polygon                   *
*       by pseudo-scanline.                                            *
*                                                                      *
************************************************************************

        implicit none
        integer :: bmpWidth
        integer :: bmpHeight
        integer :: np, ix(np), iy(np)
        integer :: npw, ixw(npw), iyw(npw)
        integer :: pbuf(bmpWidth*bmpHeight)       ! 0:unknown, 1:online, 2:inside, -1:outside
        integer :: ncall

*-----------------------------------------------------------------------

        integer :: ipxy(np)
        integer :: ixrange(2), iyrange(2)
        integer :: ixo, iyo, ip, ip2, pbval
        integer :: i, j, k
        integer :: ixmin, ixmax, iymin, iymax
        logical :: inside

*-----------------------------------------------------------------------

        ncall = 0

        call bitmap_get_range_polygon(
     &          bmpWidth, bmpHeight,
     &          npw, ixw, iyw, ixrange, iyrange)


        ixmin = ixrange(1)
        ixmax = ixrange(2)
        iymin = iyrange(1)
        iymax = iyrange(2)

        do j = iymin, iymax
          do i = ixmin, ixmax
            ixo = i
            iyo = j
            ip = (iyo-1)*bmpWidth + ixo

            if ( pbuf(ip).eq.0 ) then

              call bitmap_fill_point_is_inside_winding(
     &                npw, ixw, iyw, ixo, iyo, inside)
              ncall = ncall + 1

              if ( inside ) then
                pbuf(ip) = 2
              else
                pbuf(ip) = -1
              end if

            end if

            pbval = pbuf(ip)
            if ( pbval.ne.1 ) then
              if ( i.ne.ixmax ) then     ! right
                ixo = i+1
                iyo = j
                ip2 = (iyo-1)*bmpWidth + ixo

                if ( pbuf(ip2).eq.0 ) then
                  pbuf(ip2) = pbval
                end if
              end if

              if ( j.ne.iymax ) then     ! up
                ixo = i
                iyo = j+1
                ip2 = (iyo-1)*bmpWidth + ixo

                if ( pbuf(ip2).eq.0 ) then
                  pbuf(ip2) = pbval

                  if ( i.ne.ixmin ) then    ! up-left
                    ixo = i-1
                    iyo = j+1
                    ip2 = (iyo-1)*bmpWidth + ixo
                    if ( pbuf(ip2).eq.0 ) then
                      pbuf(ip2) = pbval
                    end if
                  end if

                end if
              end if
            end if

          end do
        end do

      end subroutine



************************************************************************
*                                                                      *
      subroutine bitmap_draw_set_paintbuffer(
     &                bmpWidth, bmpHeight,
     &                npoints, ix, iy, pbuf)
*                                                                      *
*       draw boundary of polygon to paint buffer. 1 means boundary.    *
*                                                                      *
************************************************************************

        implicit none
        integer :: bmpWidth
        integer :: bmpHeight
        integer :: npoints
        integer :: ix(npoints), iy(npoints)
        integer :: pbuf(bmpWidth*bmpHeight)

*-----------------------------------------------------------------------

        integer :: i, ip

        do i = 1, npoints
          ip = (iy(i)-1)*bmpWidth + ix(i)
          pbuf(ip) = 1
        end do
      end subroutine



************************************************************************
*                                                                      *
      subroutine bitmap_draw_pixel(
     &      bmpWidth, bmpHeight, npoints,
     &      ix, iy, bgrQuad, rgbaSize, rcol)
*                                                                      *
*       draw outer pixels in polygon.                                  *
*                                                                      *
************************************************************************

        implicit none
        integer :: bmpWidth
        integer :: bmpHeight
        integer :: npoints
        integer :: ipxy(npoints)
        integer :: ix(npoints), iy(npoints)
        integer :: rgbaSize
        character(1) :: bgrQuad(rgbaSize)
        real(8) :: rcol(3)

*-----------------------------------------------------------------------

        integer :: i, j, ip
        integer :: ixo, iyo, ipos
        integer :: ixw(npoints), iyw(npoints), npw
        integer :: pbuf(bmpWidth*bmpHeight)

*-----------------------------------------------------------------------

        call bitmap_fill_remove_same_points(
     &          npoints, ix, iy,
     &          npw, ixw, iyw)

        do i = 1, bmpWidth*bmpHeight
          pbuf(i) = 0
        end do

        call bitmap_draw_set_paintbuffer(
     &          bmpWidth, bmpHeight, npw, ixw, iyw, pbuf)

        do j = 1, bmpHeight
          do i = 1, bmpWidth
            ixo = i
            iyo = j
            ip = (iyo-1)*bmpWidth + ixo
            if ( pbuf(ip).gt.0 ) then
              ipos = (iyo-1)*bmpWidth*4 + (ixo-1)*4 + 1
              call set_bitmap_bgrQuad(
     &                bgrQuad, rcol, ipos)
            end if
          end do
        end do

      end subroutine



************************************************************************
*                                                                      *
      subroutine bitmap_interpolate(jbt, ncount, ixo, iyo, ixo2, iyo2)
*                                                                      *
*       interpolate adjacent points which have skipping pixels         *
*       between them.                                                  *
*                                                                      *
************************************************************************

        implicit none
        integer :: jbt               ! (I)
        integer :: ncount            ! (IO)
        integer :: ixo, iyo          ! (I)
        integer :: ixo2, iyo2        ! (I)

*-----------------------------------------------------------------------

        integer :: iww, ihh, ie
        integer :: i
        integer :: ixot, iyot

*-----------------------------------------------------------------------

        iww = abs(ixo2 - ixo)
        ihh = abs(iyo2 - iyo)

        if ( ihh.le.1 .and. iww.le.1 ) then
          write(jbt) ixo2,iyo2
          ncount = ncount + 1
        else
          if (ixo.eq.ixo2) then
            if ( iyo2.ge.iyo ) then
              do i = iyo+1, iyo2
                write(jbt) ixo2,i
                ncount = ncount + 1
              end do
            else
              do i = iyo-1, iyo2, -1
                write(jbt) ixo2,i
                ncount = ncount + 1
              end do
            end if
          else if (iyo.eq.iyo2) then
            if ( ixo2.ge.ixo ) then
              do i = ixo+1, ixo2
                write(jbt) i,iyo2
                ncount = ncount + 1
              end do
            else
              do i = ixo-1, ixo2, -1
                write(jbt) i,iyo2
                ncount = ncount + 1
              end do
            end if
          else

            ! bresenham
            if (iww.ge.ihh) then
              if (ixo2.gt.ixo) then
                ie = iww
                ixot = ixo
                iyot = iyo
                do while (ixot.lt.ixo2)
                  ixot = ixot + 1
                  ie = ie + ihh*2
                  if (ie.ge.iww*2) then
                    ie = ie - iww*2
                    if (iyo2.gt.iyo) then
                      iyot = iyot + 1
                    else
                      iyot = iyot - 1
                    end if
                  end if
                  write(jbt) ixot, iyot
                  ncount = ncount + 1
                end do
              else
                ie = iww
                ixot = ixo
                iyot = iyo
                do while (ixot.gt.ixo2)
                  ixot = ixot - 1
                  ie = ie + ihh*2
                  if (ie.ge.iww*2) then
                    ie = ie - iww*2
                    if (iyo2.gt.iyo) then
                      iyot = iyot + 1
                    else
                      iyot = iyot - 1
                    end if
                  end if
                  write(jbt) ixot, iyot
                  ncount = ncount + 1
                end do
              end if
            else
              if (iyo2.gt.iyo) then
                ie = ihh
                ixot = ixo
                iyot = iyo
                do while (iyot.lt.iyo2)
                  iyot = iyot + 1
                  ie = ie + iww*2
                  if (ie.ge.ihh*2) then
                    ie = ie - ihh*2
                    if (ixo2.gt.ixo) then
                      ixot = ixot + 1
                    else
                      ixot = ixot - 1
                    end if
                  end if
                  write(jbt) ixot, iyot
                  ncount = ncount + 1
                end do
              else
                ie = ihh
                ixot = ixo
                iyot = iyo
                do while (iyot.gt.iyo2)
                  iyot = iyot - 1
                  ie = ie + iww*2
                  if (ie.ge.ihh*2) then
                    ie = ie - ihh*2
                    if (ixo2.gt.ixo) then
                      ixot = ixot + 1
                    else
                      ixot = ixot - 1
                    end if
                  end if
                  write(jbt) ixot, iyot
                  ncount = ncount + 1
                end do
              end if
            end if
          end if
        end if
      end subroutine



************************************************************************
*                                                                      *
      subroutine write_bitmap_header(iubmp, bmpWidth, bmpHeight)
*                                                                      *
*       write bitmap header                                            *
*                                                                      *
************************************************************************

        implicit none
        integer :: iubmp
        integer :: bmpWidth
        integer :: bmpHeight

*-----------------------------------------------------------------------

        ! Bitmap file header
        !
        type BITMAPFILEHEADER
          integer(kind=2) :: bfType           ! file type = 'BM'
          integer(kind=4) :: bfSize           ! file size (byte)
          integer(kind=2) :: bfReserved1      ! 0
          integer(kind=2) :: bfReserved2      ! 0
          integer(kind=4) :: bfOffBits        ! offset from top of file to image data (byte) = 54
        end type

        ! Bitmap information header
        !
        type BITMAPINFOHEADER
          integer(kind=4) :: biSize            ! size of information header (byte) = 40
          integer(kind=4) :: biWidth           ! width of image (pixel)
          integer(kind=4) :: biHeight          ! height of image (pixel)
          integer(kind=2) :: biPlanes          ! 1
          integer(kind=2) :: biBitCount        ! data-size/pixel (bit) (32=true color)
          integer(kind=4) :: biCompression     ! 0 (no compression)
          integer(kind=4) :: biSizeImage       ! size of image data (byte)
          integer(kind=4) :: biXPixPerMeter    ! horizontal resolution (3780=96dpi)
          integer(kind=4) :: biYPixPerMeter    ! vertical resolution (3780=96dpi)
          integer(kind=4) :: biClrUsed         ! 0
          integer(kind=4) :: biCirImportant    ! 0
        end type

*-----------------------------------------------------------------------

        type(BITMAPFILEHEADER), allocatable :: bmpFile
        type(BITMAPINFOHEADER), allocatable :: bmpInfo

        logical :: is_little_endian

*-----------------------------------------------------------------------

        allocate(bmpFile)
        allocate(bmpInfo)

        bmpFile%bfType = 19778            ! 'BM'
        bmpFile%bfReserved1 = 0
        bmpFile%bfReserved2 = 0
        bmpFile%bfSize = bmpWidth*bmpHeight*4 + 54
        bmpFile%bfOffBits = 54

        bmpInfo%biSize = 40
        bmpInfo%biWidth = bmpWidth
        bmpInfo%biHeight = bmpHeight
        bmpInfo%biPlanes = 1
        bmpInfo%biBitCount = 32
        bmpInfo%biCompression = 0
        bmpInfo%biSizeImage = bmpWidth*bmpHeight*4
        bmpInfo%biXPixPerMeter = 3780
        bmpInfo%biYPixPerMeter = 3780
        bmpInfo%biClrUsed = 0
        bmpInfo%biCirImportant = 0

          call write_int2_le(iubmp, bmpFile%bfType)
          call write_int4_le(iubmp, bmpFile%bfSize)
          call write_int2_le(iubmp, bmpFile%bfReserved1)
          call write_int2_le(iubmp, bmpFile%bfReserved2)
          call write_int4_le(iubmp, bmpFile%bfOffBits)
          call write_int4_le(iubmp, bmpInfo%biSize)
          call write_int4_le(iubmp, bmpInfo%biWidth)
          call write_int4_le(iubmp, bmpInfo%biHeight)
          call write_int2_le(iubmp, bmpInfo%biPlanes)
          call write_int2_le(iubmp, bmpInfo%biBitCount)
          call write_int4_le(iubmp, bmpInfo%biCompression)
          call write_int4_le(iubmp, bmpInfo%biSizeImage)
          call write_int4_le(iubmp, bmpInfo%biXPixPerMeter)
          call write_int4_le(iubmp, bmpInfo%biYPixPerMeter)
          call write_int4_le(iubmp, bmpInfo%biClrUsed)
          call write_int4_le(iubmp, bmpInfo%biCirImportant)

        deallocate(bmpInfo)
        deallocate(bmpFile)

      end subroutine



************************************************************************
*                                                                      *
      subroutine hsv2rgb(h, s, v, rgb)
*                                                                      *
*       convert hsv color value to rgb.                                *
*                                                                      *
************************************************************************

        implicit none
        real(8) :: h, s, v
        real(8),dimension(3) :: rgb

*-----------------------------------------------------------------------

        real(8) :: r, g, b
        real(8) :: h6, f
        integer :: i

*-----------------------------------------------------------------------

        r = v
        g = v
        b = v
        if ( s.gt.0.0d0 ) then
          h6 = h * 6.0d0
          i = floor(h6)
          f = h6 - i
          if ( i.eq.0 ) then
            g = g * (1.0d0-s*(1.0d0-f))
            b = b * (1.0d0-s)
          else if ( i.eq.1 ) then
            r = r * (1.0d0-s*f)
            b = b * (1.0d0-s)
          else if ( i.eq.2 ) then
            r = r * (1.0d0-s)
            b = b * (1.0d0-s*(1.0d0-f))
          else if ( i.eq.3 ) then
            r = r * (1.0d0-s)
            g = g * (1.0d0-s*f)
          else if ( i.eq.4 ) then
            r = r * (1.0d0-s*(1.0d0-f))
            g = g * (1.0d0-s)
          else if ( i.eq.5 ) then
            g = g * (1.0d0-s)
            b = b * (1.0d0-s*f)
          else
            print *, "Illegal hue value in hsv2rgb", h, i
          end if
        end if

        rgb(1) = r
        rgb(2) = g
        rgb(3) = b
      end subroutine

