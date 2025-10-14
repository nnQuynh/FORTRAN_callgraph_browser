************************************************************************
*                                                                      *
      subroutine vtk_create_filename(fname, sfname, is_geom,
     &                   itime, ntime)
*                                                                      *
*       create filename of vtk file                                    *
*                                                                      *
************************************************************************
        implicit none
        character(1) :: fname(100)    ! (I)
        character(len=*) :: sfname    ! (O)
        logical :: is_geom
        integer :: itime
        integer :: ntime

*-----------------------------------------------------------------------

        integer :: isize, ipos, newipos
        integer :: get_line_length
        integer :: last_char_index
        integer :: maxdig, get_digits
        character(len=10) :: format
        character(len=20) :: sit

*-----------------------------------------------------------------------

        call clear_string(sfname)
        isize = get_line_length(fname, 100)
        call char2string(isize, fname, sfname)
        ipos = last_char_index(sfname, '.')
        ipos = ipos - 1

        call clear_string(sfname)
        call char2string(ipos, fname, sfname)
        isize = ipos

        if ( is_geom ) then
          ! geometry

          newipos = ipos + 6
          sfname(ipos+1:newipos) = '_g.vtk'
          isize = newipos

        else
          ! tally

          if ( ntime.le.1 ) then

            newipos = ipos + 4
            sfname(ipos+1:newipos) = '.vtk'
            isize = newipos

          else

            maxdig = get_digits(ntime)
            write(format,'(''(I'',1I0,''.'',1I0,'')'')') maxdig, maxdig
            write(sit,format) itime

            newipos = ipos + 1 + len_trim(sit) + 4
            sfname(ipos+1:newipos) = '_'//trim(sit)//'.vtk'
            isize = newipos

          end if

        end if
      end subroutine



************************************************************************
*                                                                      *
      subroutine vtk_write_polydata(
     &        itfmt,
     &        iuvtk, iuspoly, iusnpts)
*                                                                      *
*       write polydata of geometry to vtk file                         *
*                                                                      *
************************************************************************

        implicit none
        integer :: itfmt
        integer :: iuvtk
        integer :: iuspoly
        integer :: iusnpts

*-----------------------------------------------------------------------

        logical :: isText

*-----------------------------------------------------------------------

        isText = (itfmt.eq.0)

        if ( isText ) then
          call vtk_write_polydata_ascii(
     &            iuvtk, iuspoly, iusnpts)
        else
          call vtk_write_polydata_binary(
     &            iuvtk, iuspoly, iusnpts)
        end if

      end subroutine



************************************************************************
*                                                                      *
      subroutine vtk_write_polydata_ascii(
     &        iuvtk, iuspoly, iusnpts)
*                                                                      *
*       write polydata of geometry to vtk file                         *
*                                                                      *
************************************************************************

        implicit none
        integer :: iuvtk
        integer :: iuspoly
        integer :: iusnpts

*-----------------------------------------------------------------------

        integer :: numPoly
        integer :: npoints, nTotalPoints, nTotalPolygon
        integer :: i, j, icount
        real(8) :: px, py, pz
        integer, allocatable :: np(:)

*-----------------------------------------------------------------------

        numPoly = 0
        do while (.true.)
          read(iusnpts,end=100) npoints
          numPoly = numPoly + 1
        end do
  100   continue

        rewind(iusnpts)

        allocate( np(numPoly) )

        nTotalPoints = 0
        do i = 1, numPoly
          read(iusnpts) npoints
          nTotalPoints = nTotalPoints + npoints
          np(i) = npoints
        end do

*-----------------------------------------------------------------------

        write(iuvtk,'(''# vtk DataFile Version 2.0'')')
        write(iuvtk,'(''vtk output'')')
        write(iuvtk,'(''ASCII'')')
        write(iuvtk,'(''DATASET POLYDATA'')')
        write(iuvtk,'(''POINTS '',I10,'' float'')') nTotalPoints

        nTotalPoints = 0
        nTotalPolygon = 0
        do i = 1, numPoly
          do j = 1, np(i)
            read(iuspoly) px, py, pz
            write(iuvtk,'(1p3e12.4)') px, py, pz
          end do

          ! ignore polygon which has only 1 point.
          if ( np(i).gt.1 ) then
            nTotalPoints = nTotalPoints + np(i)
            nTotalPolygon = nTotalPolygon + 1
          end if
        end do

*-----------------------------------------------------------------------

        write(iuvtk,'(''POLYGONS '',I10,1X,I10)')
     &        nTotalPolygon,(nTotalPoints+nTotalPolygon)

        icount = 0
        do i = 1, numPoly
          if ( np(i).eq.1 ) then
            icount = icount + 1
            cycle
          end if

          write(iuvtk,'(I0)',advance='no') np(i)

          do j = 1, np(i)
            write(iuvtk,'(1X,I0)',advance='no') icount
            icount = icount + 1

            if ( j.eq.np(i) ) then
              write(iuvtk,'()')
            else if ( mod(j,20).eq.0 ) then
              write(iuvtk,'()')
            end if
          end do
        end do

        deallocate( np )

      end subroutine



************************************************************************
*                                                                      *
      subroutine vtk_write_polydata_binary(
     &        iout, iuspoly, iusnpts)
*                                                                      *
*       write polydata of geometry to vtk file                         *
*                                                                      *
************************************************************************

        implicit none
        integer :: iout
        integer :: iuspoly
        integer :: iusnpts

*-----------------------------------------------------------------------

        integer :: numPoly
        integer :: npoints, nTotalPoints, nTotalPolygon
        integer :: i, j, n, idpoint
        real(8) :: pval(3)
        integer, allocatable :: np(:)

        integer :: ncsum, icount
        logical :: isLittleEndian, is_little_endian

*-----------------------------------------------------------------------

        isLittleEndian = is_little_endian(iout)

        numPoly = 0
        do while (.true.)
          read(iusnpts,end=100) npoints
          numPoly = numPoly + 1
        end do
  100   continue

        rewind(iusnpts)

        allocate( np(numPoly) )

        nTotalPoints = 0
        do i = 1, numPoly
          read(iusnpts) npoints
          nTotalPoints = nTotalPoints + npoints
          np(i) = npoints
        end do

*-----------------------------------------------------------------------

        ncsum = 0

        call write_binary_string(iout,
     &          '# vtk DataFile Version 2.0', .true., icount)
        call write_line_break_lf(iout)
        ncsum = ncsum + icount + 1

        call write_binary_string(iout, 'vtk_output', .true., icount)
        call write_line_break_lf(iout)
        ncsum = ncsum + icount + 1

        call write_binary_string(iout, 'BINARY', .true., icount)
        call write_line_break_lf(iout)
        ncsum = ncsum + icount + 1

        call write_binary_string(iout,
     &          'DATASET POLYDATA', .true., icount)
        call write_line_break_lf(iout)
        ncsum = ncsum + icount + 1

        call write_binary_string(iout, 'POINTS ', .false., icount)
        ncsum = ncsum + icount
        call write_integer_as_text_binary(iout, nTotalPoints, icount)
        ncsum = ncsum + icount
        call write_binary_string(iout, ' float', .false., icount)
        ncsum = ncsum + icount
        call write_adjust_4bytes_pre(iout, ncsum, icount)
        ncsum = ncsum + icount
        call write_line_break_lf(iout)
        ncsum = ncsum + 1


        nTotalPoints = 0
        nTotalPolygon = 0
        do i = 1, numPoly
          do j = 1, np(i)
            read(iuspoly) (pval(n), n=1,3)
            call write_real8s_as_real4s_be(
     &              iout, 3, pval, isLittleEndian, icount)
            ncsum = ncsum + icount
          end do

          ! ignore polygon which has only 1 point.
          if ( np(i).gt.1 ) then
            nTotalPoints = nTotalPoints + np(i)
            nTotalPolygon = nTotalPolygon + 1
          end if
        end do

        call write_line_break_lf(iout)
        ncsum = ncsum + 1

*-----------------------------------------------------------------------

        call write_binary_string(iout, 'POLYGONS ', .false., icount)
        ncsum = ncsum + icount
        call write_integer_as_text_binary(iout, nTotalPolygon, icount)
        ncsum = ncsum + icount
        call write_binary_string(iout, ' ', .false., icount)
        ncsum = ncsum + icount
        call write_integer_as_text_binary(
     &          iout, (nTotalPoints+nTotalPolygon), icount)
        ncsum = ncsum + icount
        call write_adjust_4bytes_pre(iout, ncsum, icount)
        ncsum = ncsum + icount
        call write_line_break_lf(iout)
        ncsum = ncsum + 1


        idpoint = 0
        do i = 1, numPoly
          if ( np(i).eq.1 ) then
            icount = icount + 1
            cycle
          end if

          call write_int4_be(iout, np(i))
          ncsum = ncsum + 4

          do j = 1, np(i)
            call write_int4_be(iout, idpoint)
            ncsum = ncsum + 4
            idpoint = idpoint + 1
          end do
        end do

        deallocate( np )

      end subroutine



************************************************************************
*                                                                      *
      subroutine vtk_write_tally(
     &        itfmt,
     &        isunit_meta, isunit, isunit_regmat,
     &        fname)
*                                                                      *
*       write tally values to vtk file                                 *
*                                                                      *
************************************************************************
        implicit none

        integer :: itfmt
        integer :: isunit_meta
        integer :: isunit
        integer :: isunit_regmat
        character(1) :: fname(100)

*-----------------------------------------------------------------------

        character(len=255) :: outfile
        logical :: is_geom = .false.
        integer :: iout, ios
        integer :: iout_default = 94

        integer :: ntime, it
        integer :: ncsum
        logical :: isText
        logical :: isLittleEndian, is_little_endian

*-----------------------------------------------------------------------

        isText = (itfmt.eq.0)

        read(isunit_meta) ntime

        do it = 1, ntime

          call vtk_create_filename(fname, outfile, is_geom, it, ntime)
          call open_file(iout_default, outfile, iout, ios, isText)
          isLittleEndian = .true.
          if ( .not.isText ) then
            isLittleEndian = is_little_endian(iout)
          end if

          ncsum = 0

          call vtk_write_tally_celldata(
     &            isText, isLittleEndian,
     &            isunit_meta, isunit, iout, ncsum)

          rewind(isunit_regmat)
          call vtk_write_cell_regmat(
     &            isText, isLittleEndian,
     &            isunit_regmat, iout, ncsum)

          call close_file(iout)

        end do

      end subroutine



************************************************************************
*                                                                      *
      subroutine vtk_write_tally_celldata(
     &        isText, isLittleEndian,
     &        isunit_meta, isunit, iout, ncsum)
*                                                                      *
*       read tally value on cells from scratch and write to vtk file   *
*                                                                      *
************************************************************************
        implicit none
        logical :: isText
        logical :: isLittleEndian
        integer :: isunit_meta
        integer :: isunit
        integer :: iout
        integer :: ncsum

*-----------------------------------------------------------------------

        integer :: nx, ny, nz
        integer :: ix, iy, iz
        real(8), allocatable :: xval(:), yval(:), zval(:)
        real(8), allocatable :: val(:,:,:)

        integer :: iaxis
        integer :: nparam
        character(len=1), allocatable :: paramName(:)
        integer, allocatable :: paramTotal(:)
        integer, allocatable :: paramCurrent(:)

        integer :: nparamTotal
        integer :: nsum
        integer :: icount, i, n

        character(len=100) :: tallyName

*-----------------------------------------------------------------------

        read(isunit_meta) iaxis
        read(isunit_meta) nparam

        allocate(paramName(nparam))
        allocate(paramTotal(nparam))
        allocate(paramCurrent(nparam))

        do i = 1, nparam
          read(isunit_meta) paramName(i)
        end do
        nparamTotal = 1
        do i = 1, nparam
          read(isunit_meta) paramTotal(i)
          nparamTotal = nparamTotal * paramTotal(i)
        end do

        read(isunit_meta) nx, ny, nz

        nsum = nx*ny*nz

        allocate(xval(nx+1))
        allocate(yval(ny+1))
        allocate(zval(nz+1))
        allocate(val(nx+1,ny+1,nz+1))
        read(isunit_meta) ( xval(ix), ix=1,nx+1 )
        read(isunit_meta) ( yval(iy), iy=1,ny+1 )
        read(isunit_meta) ( zval(iz), iz=1,nz+1 )

*-----------------------------------------------------------------------

        if ( isText ) then
          ! ASCII

          write(iout,'(''# vtk DataFile Version 3.0'')')
          write(iout,'(''vtk_output'')')
          write(iout,'(''ASCII'')')
          write(iout,'(''DATASET RECTILINEAR_GRID'')')
          write(iout,'(''DIMENSIONS '',3(1X,I0))')
     &            (nx+1), (ny+1), (nz+1)

          write(iout,'(''X_COORDINATES '',I0,'' float'')') (nx+1)
          do ix = 1, nx+1
            write(iout,'(1p1e12.4)',advance='no') xval(ix)
            if ( (mod(ix,10).eq.0) .and. (ix.ne.(nx+1)) ) then
              call write_line_break(iout)
            end if
          end do
          call write_line_break(iout)
          write(iout,'(''Y_COORDINATES '',I0,'' float'')') (ny+1)
          do iy = 1, ny+1
            write(iout,'(1p1e12.4)',advance='no') yval(iy)
            if ( (mod(iy,10).eq.0) .and. (iy.ne.(ny+1)) ) then
              call write_line_break(iout)
            end if
          end do
          call write_line_break(iout)
          write(iout,'(''Z_COORDINATES '',I0,'' float'')') (nz+1)
          do iz = 1, nz+1
            write(iout,'(1p1e12.4)',advance='no') zval(iz)
            if ( (mod(iz,10).eq.0) .and. (iz.ne.(nz+1)) ) then
              call write_line_break(iout)
            end if
          end do
          call write_line_break(iout)

          call write_line_break(iout)

          write(iout,'(''CELL_DATA '',I0)') nsum
          write(iout,'(''FIELD FieldData '',I0)') (nparamTotal+2)         ! with region and material

        else
          ! BINARY

          call write_binary_string(iout,
     &            '# vtk DataFile Version 3.0', .true., icount)
          call write_line_break_lf(iout)
          ncsum = ncsum + icount + 1

          call write_binary_string(iout, 'vtk_output', .true., icount)
          call write_line_break_lf(iout)
          ncsum = ncsum + icount + 1

          call write_binary_string(iout, 'BINARY', .true., icount)
          call write_line_break_lf(iout)
          ncsum = ncsum + icount + 1

          call write_binary_string(iout,
     &            'DATASET RECTILINEAR_GRID', .true., icount)
          call write_line_break_lf(iout)
          ncsum = ncsum + icount + 1

          call write_binary_string(iout, 'DIMENSIONS ', .false., icount)
          ncsum = ncsum + icount
          call write_integer_as_text_binary(iout, (nx+1), icount)
          ncsum = ncsum + icount
          call write_binary_string(iout, ' ', .false., icount)
          ncsum = ncsum + icount
          call write_integer_as_text_binary(iout, (ny+1), icount)
          ncsum = ncsum + icount
          call write_binary_string(iout, ' ', .false., icount)
          ncsum = ncsum + icount
          call write_integer_as_text_binary(iout, (nz+1), icount)
          ncsum = ncsum + icount
          call write_line_break_lf(iout)
          ncsum = ncsum + 1

          call write_binary_string(
     &            iout, 'X_COORDINATES ', .false., icount)
          ncsum = ncsum + icount
          call write_integer_as_text_binary(iout, (nx+1), icount)
          ncsum = ncsum + icount
          call write_binary_string(iout, ' float', .false., icount)
          ncsum = ncsum + icount
          call write_adjust_4bytes_pre(iout, ncsum, icount)
          ncsum = ncsum + icount
          call write_line_break_lf(iout)
          ncsum = ncsum + 1

          call write_real8s_as_real4s_be(
     &            iout, (nx+1), xval, isLittleEndian, icount)
          ncsum = ncsum + icount
          call write_line_break_lf(iout)
          ncsum = ncsum + 1

          call write_binary_string(
     &            iout, 'Y_COORDINATES ', .false., icount)
          ncsum = ncsum + icount
          call write_integer_as_text_binary(iout, (ny+1), icount)
          ncsum = ncsum + icount
          call write_binary_string(iout, ' float', .false., icount)
          ncsum = ncsum + icount
          call write_adjust_4bytes_pre(iout, ncsum, icount)
          ncsum = ncsum + icount
          call write_line_break_lf(iout)
          ncsum = ncsum + 1

          call write_real8s_as_real4s_be(
     &            iout, (ny+1), yval, isLittleEndian, icount)
          ncsum = ncsum + icount
          call write_line_break_lf(iout)
          ncsum = ncsum + 1

          call write_binary_string(
     &            iout, 'Z_COORDINATES ', .false., icount)
          ncsum = ncsum + icount
          call write_integer_as_text_binary(iout, (nz+1), icount)
          ncsum = ncsum + icount
          call write_binary_string(iout, ' float', .false., icount)
          ncsum = ncsum + icount
          call write_adjust_4bytes_pre(iout, ncsum, icount)
          ncsum = ncsum + icount
          call write_line_break_lf(iout)
          ncsum = ncsum + 1

          call write_real8s_as_real4s_be(
     &            iout, (nz+1), zval, isLittleEndian, icount)
          ncsum = ncsum + icount
          call write_line_break_lf(iout)
          ncsum = ncsum + 1

          call write_line_break_lf(iout)
          ncsum = ncsum + 1

          call write_binary_string(iout, 'CELL_DATA ', .false., icount)
          ncsum = ncsum + icount
          call write_integer_as_text_binary(iout, nsum, icount)
          ncsum = ncsum + icount
          call write_line_break_lf(iout)
          ncsum = ncsum + 1

          call write_binary_string(
     &            iout, 'FIELD FieldData ', .false., icount)
          ncsum = ncsum + icount
          call write_integer_as_text_binary(
     &            iout, (nparamTotal+2), icount)          ! with region and material
          ncsum = ncsum + icount
          call write_line_break_lf(iout)
          ncsum = ncsum + 1

        end if

*-----------------------------------------------------------------------

        do n = 1, nparamTotal

          do i = 1, nparam
            read(isunit_meta) paramCurrent(i)
          end do

          read(isunit) ( ( ( val(ix,iy,iz), ix=1,nx ),
     &                                      iy=1,ny ),
     &                                      iz=1,nz )

          call vtk_create_tally_name(
     &            nparam, paramName, paramTotal, paramCurrent,
     &            tallyName)


          if ( isText ) then
            ! ASCII

            write(iout,'(A,1X,I0,1X,I0,1X,A)')
     &              trim(tallyName), 1, nsum, 'float'

            icount = 0
            do iz = 1, nz
              do iy = 1, ny
                do ix = 1, nx
                  icount = icount + 1
                  write(iout,'(1p1e12.4)',advance='no') val(ix,iy,iz)
                  if ( (mod(icount,10).eq.0) .and.
     &                 (icount.ne.nsum) ) then
                    call write_line_break(iout)
                  end if
                end do
              end do
            end do
            call write_line_break(iout)

          else
            ! BINARY

            call write_binary_string(
     &              iout, trim(tallyName), .false., icount)
            ncsum = ncsum + icount
            call write_binary_string(iout, ' ', .false., icount)
            ncsum = ncsum + icount
            call write_integer_as_text_binary(iout, 1, icount)
            ncsum = ncsum + icount
            call write_binary_string(iout, ' ', .false., icount)
            ncsum = ncsum + icount
            call write_integer_as_text_binary(iout, nsum, icount)
            ncsum = ncsum + icount
            call write_binary_string(iout, ' float', .false., icount)
            ncsum = ncsum + icount
            call write_adjust_4bytes_pre(iout, ncsum, icount)
            ncsum = ncsum + icount
            call write_line_break_lf(iout)
            ncsum = ncsum + 1

            do iz = 1, nz
              do iy = 1, ny
                call write_real8s_as_real4s_be(
     &                  iout, nx, val(1,iy,iz), isLittleEndian, icount)
                ncsum = ncsum + icount
              end do
            end do
            call write_line_break_lf(iout)
            ncsum = ncsum + 1

          end if

        end do

*-----------------------------------------------------------------------

        deallocate(paramCurrent)
        deallocate(paramTotal)
        deallocate(paramName)

        deallocate(xval)
        deallocate(yval)
        deallocate(zval)
        deallocate(val)

      end subroutine



************************************************************************
*                                                                      *
      subroutine vtk_create_tally_name(
     &        psize, pname, ptotal, pcur, tallyName)
*                                                                      *
*       return tally name for vtk                                      *
*                                                                      *
************************************************************************
        implicit none
        integer :: psize
        character(len=1) :: pname(psize)
        integer :: ptotal(psize)
        integer :: pcur(psize)
        character(len=*) :: tallyName

*-----------------------------------------------------------------------

        logical :: noParam = .true.
        integer :: maxdig, max_digits_integers
        character(len=20) :: sparamid
        character(len=10) :: format
        integer :: ipos, newipos
        integer :: i

*-----------------------------------------------------------------------

        do i = 1, psize
          if ( ptotal(i).gt.1 ) then
            noParam = .false.
          end if
        end do

        if ( noParam ) then

          tallyName = "all"

        else

          maxdig = max_digits_integers(ptotal, psize)
          write(format,'(''(I'',1I0,''.'',1I0,'')'')') maxdig, maxdig

          tallyName = ""
          ipos = 0
          do i = 1, psize
            if ( ptotal(i).gt.1 ) then
              write(sparamid,format) pcur(i)
              newipos = ipos + 1 + len_trim(sparamid)
              tallyName(ipos+1:newipos) = pname(i)//trim(sparamid)
              ipos = newipos
            end if
          end do

        end if

! T.Sato 2024/10/16, add by users suggestion https://meteor.nucl.kyushu-u.ac.jp/phitsforum/t/topic/3116/3
      if (adjustl(trim(tallyName)) .eq. "") then
        tallyName = "d"
      endif

      end subroutine



************************************************************************
*                                                                      *
      subroutine vtk_clip_to_geometry(
     &               isoutg, isoutgmeta,
     &               ioh, np,
     &               iaxis, axisval, ddel)
*                                                                      *
*       convert clip points data to 3d geometry points.                *
*                                                                      *
************************************************************************
        implicit none

        integer :: isoutg, isoutgmeta
        integer :: ioh, np
        integer :: iaxis
        double precision :: axisval, ddel

*-----------------------------------------------------------------------

        integer :: i, k, kk
        integer, parameter :: kc = 3
        double precision :: x(kc), y(kc)
        double precision :: xval, yval, zval
        double precision :: x23, y23, r23, x12, y12, r12
        double precision :: sddel, rdcos
        integer :: numPoints

*-----------------------------------------------------------------------

        sddel = sqrt( ddel )

        k = 1
        x(1) = 0.d0
        y(1) = 0.d0

        numPoints = 0

        do i = 1, np
          read(ioh) x(k), y(k)

          if( k .eq. 3 ) then
            x23 = ( x(3) - x(2) )
            y23 = ( y(3) - y(2) )
            r23 = sqrt( x23**2 + y23**2 )

            if( r23 .lt. sddel ) cycle

            rdcos = ( x12 * x23 + y12 * y23 ) / r12 / r23

            if( rdcos .gt. 0.9999d0 ) then
              x(2) = x(3)
              y(2) = y(3)
              x12 = ( x(2) - x(1) )
              y12 = ( y(2) - y(1) )
              r12 = sqrt( x12**2 + y12**2 )

              if( r12 .lt. sddel ) then
                 k = 2
                 cycle
              end if
            else
              if ( iaxis.eq.1 ) then
                xval = x(1)
                yval = y(1)
                zval = axisval
              else if ( iaxis.eq.2 ) then
                xval = axisval
                yval = y(1)
                zval = x(1)
              else if ( iaxis.eq.3 ) then
                xval = y(1)
                yval = axisval
                zval = x(1)
              end if

              write(isoutg) xval, yval, zval
              numPoints = numPoints + 1

              x(1) = x(2)
              y(1) = y(2)
              x(2) = x(3)
              y(2) = y(3)

              r12  = r23
            end if
          else
            if( k .eq. 2 ) then
              x12 = ( x(2) - x(1) )
              y12 = ( y(2) - y(1) )
              r12 = sqrt( x12**2 + y12**2 )

              if( r12 .lt. sddel ) cycle
            end if

            k = k + 1
          end if

        end do

        if( k.gt.1 ) then
          do kk = 1, k-1
            if ( iaxis.eq.1 ) then
              xval = x(kk)
              yval = y(kk)
              zval = axisval
            else if ( iaxis.eq.2 ) then
              xval = axisval
              yval = y(kk)
              zval = x(kk)
            else if ( iaxis.eq.3 ) then
              xval = y(kk)
              yval = axisval
              zval = x(kk)
            end if

            write(isoutg) xval, yval, zval
            numPoints = numPoints + 1
          end do
        end if

        write(isoutgmeta) numPoints

      end subroutine



************************************************************************
*                                                                      *
      subroutine vtk_set_gshow(
     &                isoutg, isoutgmeta,
     &                nx, ny, nz, xm, ym, zm,
     &                iaxs, iuni, ires, igser,
     &                nr, mr, kr, vl, mtrns)
*                                                                      *
*       get geometry boundaries                                        *
*       and write them into scratch file                               *
*                                                                      *
************************************************************************
        implicit none

        integer :: isoutg, isoutgmeta
        integer :: nx, ny, nz
        double precision :: xm(nx+1), ym(ny+1), zm(nz+1)
        integer :: iaxs, iuni, ires, igser
        integer :: nr, mr
        integer :: kr(mr)
        double precision :: vl(nr)
        integer :: mtrns

*-----------------------------------------------------------------------

        integer :: nxcell, nycell, nzcell
        double precision, allocatable ::xvalues(:),yvalues(:),zvalues(:)
        integer, allocatable :: imat(:)
        integer, allocatable :: idxyreg(:), idxymat(:)

        double precision :: xval, yval, zval
        integer :: ix, iy, iz

*-----------------------------------------------------------------------

        if ( iaxs.eq.1 ) then
          nxcell = nx
          nycell = ny
          nzcell = 1
        else if ( iaxs.eq.2 ) then
          nxcell = 1
          nycell = ny
          nzcell = nz
        else if ( iaxs.eq.3 ) then
          nxcell = nx
          nycell = 1
          nzcell = nz
        end if

        allocate(xvalues(nxcell))
        allocate(yvalues(nycell))
        allocate(zvalues(nzcell))
        allocate(imat(nxcell*nycell*nzcell))
        allocate(idxyreg(nxcell*nycell*nzcell))
        allocate(idxymat(nxcell*nycell*nzcell))

        if ( iaxs .eq. 1 ) then

          do ix = 1, nxcell
            xvalues(ix) = ( xm(ix) + xm(ix+1) ) / 2.0d0
          end do
          do iy = 1, nycell
            yvalues(iy) = ( ym(iy) + ym(iy+1) ) / 2.0d0
          end do

          do iz = 1, nz
            zval = ( zm(iz) + zm(iz+1) ) / 2.0d0
            call vtk_gshow(
     &              isoutg, isoutgmeta,
     &              0,iaxs,iuni,ires,
     &              igser,
     &              nxcell,nycell,nzcell,
     &              xvalues,yvalues,zval,imat(1),
     &              idxyreg, idxymat,
     &              nr,mr,kr,vl,mtrns)
          end do

        else if( iaxs .eq. 2 ) then

          do iz = 1, nzcell
            zvalues(iz) = ( zm(iz) + zm(iz+1) ) / 2.0d0
          end do
          do iy = 1, nycell
            yvalues(iy) = ( ym(iy) + ym(iy+1) ) / 2.0d0
          end do

          do ix = 1, nx
            xval = ( xm(ix) + xm(ix+1) ) / 2.0d0
            call vtk_gshow(
     &              isoutg, isoutgmeta,
     &              0,iaxs,iuni,ires,
     &              igser,
     &              nzcell,nycell,nxcell,
     &              zvalues,yvalues,xval,imat(1),
     &              idxyreg, idxymat,
     &              nr,mr,kr,vl,mtrns)
          end do

        else if( iaxs .eq. 3 ) then

          do iz = 1, nzcell
            zvalues(iz) = ( zm(iz) + zm(iz+1) ) / 2.0d0
          end do
          do ix = 1, nxcell
            xvalues(ix) = ( xm(ix) + xm(ix+1) ) / 2.0d0
          end do

          do iy = 1, ny
            yval = ( ym(iy) + ym(iy+1) ) / 2.0d0
            call vtk_gshow(
     &              isoutg, isoutgmeta,
     &              0,iaxs,iuni,ires,
     &              igser,
     &              nzcell,nxcell,nycell,
     &              zvalues,xvalues,yval,imat(1),
     &              idxyreg, idxymat,
     &              nr,mr,kr,vl,mtrns)
          end do

        end if

*-----------------------------------------------------------------------

        deallocate(idxymat)
        deallocate(idxyreg)
        deallocate(imat)
        deallocate(zvalues)
        deallocate(yvalues)
        deallocate(xvalues)
      end subroutine



************************************************************************
*                                                                      *
      subroutine vtk_set_cell_regmat(
     &                isoutrm, isoutg, isoutgmeta,
     &                nx, ny, nz, xm, ym, zm,
     &                iaxs, iuni, ires, igser,
     &                nr, mr, kr, vl, mtrns)
*                                                                      *
*       get index of region and material on center point of cells      *
*       and write them into scratch file                               *
*                                                                      *
************************************************************************
        implicit none

        integer :: isoutrm
        integer :: isoutg, isoutgmeta
        integer :: nx, ny, nz
        double precision :: xm(nx+1), ym(ny+1), zm(nz+1)
        integer :: iaxs, iuni, ires, igser
        integer :: nr, mr
        integer :: kr(mr)
        double precision :: vl(nr)
        integer :: mtrns

*-----------------------------------------------------------------------

        integer :: nxcell, nycell, nzcell
        double precision, allocatable ::xvalues(:),yvalues(:),zvalues(:)
        integer, allocatable :: imat(:)
        integer, allocatable :: idxyreg(:), idxymat(:)

        double precision :: xval, yval, zval
        integer :: ix, iy, iz

        logical :: isCell = .true.
        logical :: isDebug = .false.
        integer :: i

*-----------------------------------------------------------------------

        if ( isCell ) then

          if ( iaxs.eq.1 ) then
            nxcell = nx
            nycell = ny
            nzcell = 1
          else if ( iaxs.eq.2 ) then
            nxcell = 1
            nycell = ny
            nzcell = nz
          else if ( iaxs.eq.3 ) then
            nxcell = nx
            nycell = 1
            nzcell = nz
          end if

        else

          if ( iaxs.eq.1 ) then
            nxcell = nx+1
            nycell = ny+1
            nzcell = 1
          else if ( iaxs.eq.2 ) then
            nxcell = 1
            nycell = ny+1
            nzcell = nz+1
          else if ( iaxs.eq.3 ) then
            nxcell = nx+1
            nycell = 1
            nzcell = nz+1
          end if

        end if

        allocate(xvalues(nxcell))
        allocate(yvalues(nycell))
        allocate(zvalues(nzcell))
        allocate(imat(nxcell*nycell*nzcell))
        allocate(idxyreg(nxcell*nycell*nzcell))
        allocate(idxymat(nxcell*nycell*nzcell))

*-----------------------------------------------------------------------

        if ( isCell ) then
          ! cell

          if ( iaxs.eq.1 ) then             ! xy asxi

            write(isoutrm) iaxs
            write(isoutrm) nxcell,nycell,nz

            do ix = 1, nxcell
              xvalues(ix) = ( xm(ix) + xm(ix+1) ) / 2.0d0
            end do
            do iy = 1, nycell
              yvalues(iy) = ( ym(iy) + ym(iy+1) ) / 2.0d0
            end do

            do iz = 1, nz
              zval = ( zm(iz) + zm(iz+1) ) / 2.0d0
              call vtk_gshow(
     &                isoutg, isoutgmeta,
     &                0,iaxs,iuni,ires,
     &                igser,
     &                nxcell,nycell,nzcell,
     &                xvalues,yvalues,zval,imat(1),
     &                idxyreg, idxymat,
     &                nr,mr,kr,vl,mtrns)

              write(isoutrm) ( (
     &                idxyreg((iy-1)*nxcell+ix),
     &                        iy=1,nycell ),
     &                        ix=1,nxcell )
              write(isoutrm) ( (
     &                idxymat((iy-1)*nxcell+ix),
     &                        iy=1,nycell ),
     &                        ix=1,nxcell )
            end do

          else if ( iaxs.eq.2 ) then          ! zy axis

            write(isoutrm) iaxs
            write(isoutrm) nx,nycell,nzcell

            do iz = 1, nzcell
              zvalues(iz) = ( zm(iz) + zm(iz+1) ) / 2.0d0
            end do
            do iy = 1, nycell
              yvalues(iy) = ( ym(iy) + ym(iy+1) ) / 2.0d0
            end do

            do ix = 1, nx
              xval = ( xm(ix) + xm(ix+1) ) / 2.0d0
              call vtk_gshow(
     &                isoutg, isoutgmeta,
     &                0,iaxs,iuni,ires,
     &                igser,
     &                nzcell,nycell,nxcell,
     &                zvalues,yvalues,xval,imat(1),
     &                idxyreg, idxymat,
     &                nr,mr,kr,vl,mtrns)

              write(isoutrm) ( (
     &                idxyreg((iy-1)*nzcell+iz),
     &                        iy=1,nycell ),
     &                        iz=1,nzcell )
              write(isoutrm) ( (
     &                idxymat((iy-1)*nzcell+iz),
     &                        iy=1,nycell ),
     &                        iz=1,nzcell )
            end do

          else if ( iaxs.eq.3 ) then          ! zx axis

            write(isoutrm) iaxs
            write(isoutrm) nxcell,ny,nzcell

            do iz = 1, nzcell
              zvalues(iz) = ( zm(iz) + zm(iz+1) ) / 2.0d0
            end do
            do ix = 1, nxcell
              xvalues(ix) = ( xm(ix) + xm(ix+1) ) / 2.0d0
            end do

            do iy = 1, ny
              yval = ( ym(iy) + ym(iy+1) ) / 2.0d0
              call vtk_gshow(
     &                isoutg, isoutgmeta,
     &                0,iaxs,iuni,ires,
     &                igser,
     &                nzcell,nxcell,nycell,
     &                zvalues,xvalues,yval,imat(1),
     &                idxyreg, idxymat,
     &                nr,mr,kr,vl,mtrns)

              if ( isDebug ) then
                write(6,*) 'iy,yval=', iy, yval
                write(6,*) 'nx,nz=', nxcell, nzcell
                do iz = 1, nzcell
                  do ix = 1, nxcell
                   i = (iz-1)*nxcell+ix
                   write(6,'(1X,I0)',advance='no') idxymat(i)
                 end do
                 write(6,*)
                end do
              end if

              write(isoutrm) ( (
     &                idxyreg((ix-1)*nzcell+iz),
     &                        ix=1,nxcell ),
     &                        iz=1,nzcell )
              write(isoutrm) ( (
     &                idxymat((ix-1)*nzcell+iz),
     &                        ix=1,nxcell ),
     &                        iz=1,nzcell )
            end do

          end if

*-----------------------------------------------------------------------

        else
          ! point

          if ( iaxs.eq.1 ) then             ! xy axis

            write(isoutrm) iaxs
            write(isoutrm) nx+1,ny+1,nz

            do iz = 1, nz
              zval = ( zm(iz) + zm(iz+1) ) / 2.0d0
              call vtk_gshow(
     &                isoutg, isoutgmeta,
     &                0,iaxs,iuni,ires,
     &                igser,
     &                nx+1,ny+1,nzcell,
     &                xm,ym,zval,imat(1),
     &                idxyreg, idxymat,
     &                nr,mr,kr,vl,mtrns)

              write(isoutrm) ( (
     &                idxyreg((iy-1)*(nx+1)+ix),
     &                        iy=1,ny+1 ),
     &                        ix=1,nx+1 )
              write(isoutrm) ( (
     &                idxymat((iy-1)*(nx+1)+ix),
     &                        iy=1,ny+1 ),
     &                        ix=1,nx+1 )
            end do

          else if ( iaxs.eq.2 ) then        ! yz axis

            write(isoutrm) iaxs
            write(isoutrm) nx,ny+1,nz+1

            do ix = 1, nx
              xval = ( xm(ix) + xm(ix+1) ) / 2.0d0
              call vtk_gshow(
     &                isoutg, isoutgmeta,
     &                0,iaxs,iuni,ires,
     &                igser,
     &                nz+1,ny+1,nxcell,
     &                zm,ym,xval,imat(1),
     &                idxyreg, idxymat,
     &                nr,mr,kr,vl,mtrns)

              write(isoutrm) ( (
     &                idxyreg((iy-1)*(nz+1)+iz),
     &                        iy=1,ny+1 ),
     &                        iz=1,nz+1 )
              write(isoutrm) ( (
     &                idxymat((iy-1)*(nz+1)+iz),
     &                        iy=1,ny+1 ),
     &                        iz=1,nz+1 )
            end do

          else if ( iaxs.eq.3 ) then        ! xz axis

            write(isoutrm) iaxs
            write(isoutrm) nx+1,ny,nz+1

            do iy = 1, ny
              yval = ( ym(iy) + ym(iy+1) ) / 2.0d0
              call vtk_gshow(
     &                isoutg, isoutgmeta,
     &                0,iaxs,iuni,ires,
     &                igser,
     &                nz+1,nx+1,nycell,
     &                zm,xm,yval,imat(1),
     &                idxyreg, idxymat,
     &                nr,mr,kr,vl,mtrns)
              if ( isDebug ) then
                do iz = 1, nz+1
                  do ix = 1, nx+1
                    i = (iz-1)*(nx+1)+ix
                    write(6,'(1X,I0)',advance='no') idxymat(i)
                  end do
                  write(6,*)
                end do
              end if

              write(isoutrm) ( (
     &                idxyreg((ix-1)*(nz+1)+iz),
     &                        ix=1,nx+1 ),
     &                        iz=1,nz+1 )
              write(isoutrm) ( (
     &                idxymat((ix-1)*(nz+1)+iz),
     &                        ix=1,nx+1 ),
     &                        iz=1,nz+1 )
            end do

          end if

        end if

*-----------------------------------------------------------------------

        deallocate(idxymat)
        deallocate(idxyreg)
        deallocate(imat)
        deallocate(zvalues)
        deallocate(yvalues)
        deallocate(xvalues)
      end subroutine



************************************************************************
*                                                                      *
      subroutine vtk_write_cell_regmat(
     &        isText, isLittleEndian,
     &        isunit, iout, ncsum)
*                                                                      *
*       write index of region and material into vtk file               *
*                                                                      *
************************************************************************

        implicit none
        logical :: isText
        logical :: isLittleEndian
        integer :: isunit
        integer :: iout
        integer :: ncsum

*-----------------------------------------------------------------------

        integer :: iaxis
        integer :: nx, ny, nz
        integer :: ix, iy, iz
        integer, allocatable :: idreg(:,:,:), idmat(:,:,:)
        integer :: nsum
        integer :: icount

*-----------------------------------------------------------------------
*     read
*-----------------------------------------------------------------------

        read(isunit) iaxis
        read(isunit) nx, ny, nz

        allocate(idreg(nx,ny,nz), idmat(nx,ny,nz))

        if ( iaxis.eq.1 ) then

          do iz = 1, nz
            read(isunit) ( ( idreg(ix,iy,iz), iy=1,ny ), ix=1,nx )
            read(isunit) ( ( idmat(ix,iy,iz), iy=1,ny ), ix=1,nx )
          end do

        else if ( iaxis.eq.2 ) then

          do ix = 1, nx
            read(isunit) ( ( idreg(ix,iy,iz), iy=1,ny ), iz=1,nz )
            read(isunit) ( ( idmat(ix,iy,iz), iy=1,ny ), iz=1,nz )
          end do

        else if ( iaxis.eq.3 ) then

          do iy = 1, ny
            read(isunit) ( ( idreg(ix,iy,iz), ix=1,nx ), iz=1,nz )
            read(isunit) ( ( idmat(ix,iy,iz), ix=1,nx ), iz=1,nz )
          end do

        end if

*-----------------------------------------------------------------------
*     write
*-----------------------------------------------------------------------

        nsum = nx * ny * nz

        if ( isText ) then
          ! ASCII

          write(iout,'(A,1X,I0,1X,I0,1X,A)') "region", 1, nsum, "int"
          icount = 0
          do iz = 1, nz
            do iy = 1, ny
              do ix = 1, nx
                icount = icount + 1
                write(iout,'(1X,I0)',advance='no') idreg(ix,iy,iz)
                if ( ( mod(icount,20).eq.0 ) .and. .not.
     &               ( iz.eq.nz .and. iy.eq.ny .and. ix.eq.nx ) ) then
                  write(iout,*)
                end if
              end do
            end do
          end do
          write(iout,*)

          write(iout,'(A,1X,I0,1X,I0,1X,A)') "material", 1, nsum, "int"
          icount = 0
          do iz = 1, nz
            do iy = 1, ny
              do ix = 1, nx
                icount = icount + 1
                write(iout,'(1X,I0)',advance='no') idmat(ix,iy,iz)
                if ( ( mod(icount,20).eq.0 ) .and. .not.
     &               ( iz.eq.nz .and. iy.eq.ny .and. ix.eq.nx ) ) then
                  write(iout,*)
                end if
              end do
            end do
          end do
          write(iout,*)

        else
          ! BINARY

          call write_binary_string(
     &            iout, 'region ', .false., icount)
          ncsum = ncsum + icount
          call write_integer_as_text_binary(iout, 1, icount)
          ncsum = ncsum + icount
          call write_binary_string(iout, ' ', .false., icount)
          ncsum = ncsum + icount
          call write_integer_as_text_binary(iout, nsum, icount)
          ncsum = ncsum + icount
          call write_binary_string(iout, ' int', .false., icount)
          ncsum = ncsum + icount
          call write_adjust_4bytes_pre(iout, ncsum, icount)
          ncsum = ncsum + icount

          call write_line_break_lf(iout)
          ncsum = ncsum + 1

          do iz = 1, nz
            do iy = 1, ny
              call write_int4s_be(
     &                iout, nx, idreg(1,iy,iz), icount)
              ncsum = ncsum + icount
            end do
          end do

          call write_line_break_lf(iout)
          ncsum = ncsum + 1

          call write_binary_string(
     &            iout, 'material ', .false., icount)
          ncsum = ncsum + icount
          call write_integer_as_text_binary(iout, 1, icount)
          ncsum = ncsum + icount
          call write_binary_string(iout, ' ', .false., icount)
          ncsum = ncsum + icount
          call write_integer_as_text_binary(iout, nsum, icount)
          ncsum = ncsum + icount
          call write_binary_string(iout, ' int', .false., icount)
          ncsum = ncsum + icount
          call write_adjust_4bytes_pre(iout, ncsum, icount)
          ncsum = ncsum + icount

          call write_line_break_lf(iout)
          ncsum = ncsum + 1

          do iz = 1, nz
            do iy = 1, ny
              call write_int4s_be(
     &                iout, nx, idmat(1,iy,iz), icount)
              ncsum = ncsum + icount
            end do
          end do

          call write_line_break_lf(iout)
          ncsum = ncsum + 1

        end if

        deallocate(idreg, idmat)

      end subroutine


