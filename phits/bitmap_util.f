
************************************************************************
*                                                                      *
      subroutine char2string(isize, chars, str)
*                                                                      *
*       convert array of character(1) to string of character(len=NN).  *
*                                                                      *
************************************************************************

        implicit none
        integer :: isize
        character :: chars(isize)
        character(len=isize) :: str

*-----------------------------------------------------------------------

        integer :: i

*-----------------------------------------------------------------------

        do i = 1, isize
          str(i:i) = chars(i)
        end do

      end subroutine



************************************************************************
*                                                                      *
      subroutine clear_string(str)
*                                                                      *
*       fill string with space.                                        *
*                                                                      *
************************************************************************

        implicit none
        character(len=*) :: str

*-----------------------------------------------------------------------

        integer :: i, n

*-----------------------------------------------------------------------

        n = len(str)
        do i = 1, n
          str(i:i) = ' '
        end do
      end subroutine



************************************************************************
*                                                                      *
      integer function get_line_length(chr, ichr)
*                                                                      *
*       return printable character length in line.                     *
*                                                                      *
************************************************************************

        implicit none
        character(1) :: chr(ichr)
        integer :: ichr

*-----------------------------------------------------------------------

        integer ic
        integer i

*-----------------------------------------------------------------------

        do i = ichr, 1, -1
          ic = ichar(chr(i))
          if ( ic.ge.33 .and. ic.le.126 ) then         ! within printable ascii code. ic==32 is space.
            get_line_length = i
            return
          end if
        end do

        get_line_length = 0
      end function



************************************************************************
*                                                                      *
      integer function last_char_index(str, ch)
*                                                                      *
*       return last index of character in string.                      *
*                                                                      *
************************************************************************

        implicit none
        character(len=*) :: str
        character(1) :: ch

*-----------------------------------------------------------------------

        integer :: i, nlen

*-----------------------------------------------------------------------

        last_char_index = -1
        nlen = len_trim(str)

        do i = nlen, 1, -1
          if ( str(i:i) .eq. ch ) then
            last_char_index = i
            return
          end if
        end do
      end function



************************************************************************
*                                                                      *
      integer function get_digits(int)
*                                                                      *
*       return digit length of integer.                                *
*                                                                      *
************************************************************************

        implicit none
        integer :: int

*-----------------------------------------------------------------------

        get_digits = floor(log10(int*1.0d0))+1
      end function



************************************************************************
*                                                                      *
      integer function max_digits_integers(ints, num)
*                                                                      *
*       return digit length of maximum value of integers.              *
*                                                                      *
************************************************************************

        implicit none
        integer :: num
        integer :: ints(num)
        integer :: get_digits

*-----------------------------------------------------------------------

        max_digits_integers = get_digits(maxval(ints))
      end function



************************************************************************
*                                                                      *
      subroutine open_file(iunit_default, filename, iunit, ios, isText)
*                                                                      *
*          open file.                                                  *
*                                                                      *
************************************************************************
        implicit none
        integer :: iunit_default       ! (I)
        character(len=*) :: filename   ! (I)
        integer :: iunit               ! (O)
        integer :: ios                 ! (O)
        logical :: isText              ! (I)

*-----------------------------------------------------------------------

        logical :: ok

*-----------------------------------------------------------------------

        ok = .true.
        iunit = iunit_default - 1

        do while ( ok .and. iunit.lt.200 )
          iunit = iunit + 1
          inquire(unit=iunit, opened=ok)
        end do

        if ( len_trim(filename).gt.0 ) then
          if ( isText ) then
            open(unit=iunit, file=filename, iostat=ios)
          else
            open(unit=iunit, file=filename, access='stream',
     &              status='unknown', iostat=ios)
          end if
        else
          open(unit=iunit, access='stream',
     &            status='scratch', iostat=ios)
        end if
      end subroutine



************************************************************************
*                                                                      *
      subroutine close_file(iunit)
*                                                                      *
*          close file.                                                 *
*                                                                      *
************************************************************************

        implicit none
        integer :: iunit

*-----------------------------------------------------------------------

        logical :: ok

*-----------------------------------------------------------------------

        inquire(unit=iunit, opened=ok)

        if ( ok ) then
          close(iunit)
        end if
      end subroutine



************************************************************************
*                                                                      *
      subroutine write_line_break(iunit)
*                                                                      *
*          write line break code.                                      *
*                                                                      *
************************************************************************
        implicit none
        integer :: iunit

*-----------------------------------------------------------------------

        write(iunit,'()')

      end subroutine



************************************************************************
*                                                                      *
      subroutine write_line_break_lf(iunit)
*                                                                      *
*          write line break code.                                      *
*                                                                      *
************************************************************************
        implicit none
        integer :: iunit

*-----------------------------------------------------------------------

        integer(kind=1), parameter :: lf = 10

*-----------------------------------------------------------------------

        write(iunit) lf

      end subroutine



************************************************************************
*                                                                      *
      logical function is_little_endian(iubmp)
*                                                                      *
*          return whether byte order of this machine is little endian  *
*          or big endian.                                              *
*                                                                      *
************************************************************************

        implicit none
        integer :: iubmp

*-----------------------------------------------------------------------

        integer, save :: endian = -1
        integer(kind=4) :: i4
        integer(kind=1) :: i1

*-----------------------------------------------------------------------

        if ( endian.lt.0 ) then
          rewind(iubmp)
          i4 = 1
          write(iubmp) i4
          rewind(iubmp)
          read(iubmp) i1
          endian = i1
          rewind(iubmp)
        end if

        is_little_endian = .true.
        if ( endian.eq.0 ) then
          is_little_endian = .false.
        end if

      end function



************************************************************************
*                                                                      *
      subroutine write_int2_le(iubmp, value)
*                                                                      *
*       write integer(2) value to binary file on little endian.        *
*                                                                      *
************************************************************************

        implicit none
        integer :: iubmp
        integer(kind=2) :: value

*-----------------------------------------------------------------------

        integer(kind=1) :: i1, i2

*-----------------------------------------------------------------------

        i1 = mod(value, 256)
        i2 = value/256

        write(iubmp) i1
        write(iubmp) i2
      end subroutine



************************************************************************
*                                                                      *
      subroutine write_int4_le(iubmp, value)
*                                                                      *
*       write integer(4) value to binary file on little endian.        *
*                                                                      *
************************************************************************

        implicit none
        integer :: iubmp
        integer(kind=4) :: value

*-----------------------------------------------------------------------

        integer(kind=4) :: iv
        integer(kind=1) :: i1, i2, i3, i4

        i1 = mod(value, 256)
        iv = value/256
        i2 = mod(iv, 256)
        iv = iv/256
        i3 = mod(iv, 256)
        iv = iv/256
        i4 = mod(iv, 256)

        write(iubmp) i1
        write(iubmp) i2
        write(iubmp) i3
        write(iubmp) i4
      end subroutine



************************************************************************
*                                                                      *
      subroutine write_int4_be(iubmp, value)
*                                                                      *
*       write integer(4) value to binary file on big endian.           *
*                                                                      *
************************************************************************

        implicit none
        integer :: iubmp
        integer(kind=4) :: value

*-----------------------------------------------------------------------

        integer(kind=4) :: iv
        integer(kind=1) :: i1, i2, i3, i4

*-----------------------------------------------------------------------

        i1 = mod(value, 256)
        iv = value/256
        i2 = mod(iv, 256)
        iv = iv/256
        i3 = mod(iv, 256)
        iv = iv/256
        i4 = mod(iv, 256)

        write(iubmp) i4
        write(iubmp) i3
        write(iubmp) i2
        write(iubmp) i1
      end subroutine


