
************************************************************************
*                                                                      *
      subroutine write_binary_string(iout, str, isTrim, nbytes)
*                                                                      *
*       write string as binary data to file                            *
*                                                                      *
************************************************************************

        implicit none
        integer :: iout
        character(len=*) :: str
        logical :: isTrim
        integer :: nbytes

*-----------------------------------------------------------------------

        integer :: nlen
        integer :: i
        integer(kind=1) :: ia1

*-----------------------------------------------------------------------

        if ( isTrim ) then
          nlen = len_trim(str)
        else
          nlen = len(str)
        end if

        do i = 1, nlen
          ia1 = iachar(str(i:i))
          write(iout) ia1
        end do

        nbytes = nlen

      end subroutine



************************************************************************
*                                                                      *
      subroutine write_integer_as_text_binary(iout, ival, nbytes)
*                                                                      *
*       write string as binary data to file                            *
*                                                                      *
************************************************************************

        implicit none
        integer :: iout
        integer :: ival
        integer :: nbytes

*-----------------------------------------------------------------------

        character(len=100) sval

*-----------------------------------------------------------------------

        call clear_string(sval)
        write(sval,'(I0)') ival

        call write_binary_string(iout, sval, .true., nbytes)

      end subroutine



************************************************************************
*                                                                      *
      subroutine write_real4_be(iout, value, isLittleEndian, nbytes)
*                                                                      *
*       write real(4) value to binary file on big endian.              *
*                                                                      *
************************************************************************

        implicit none
        integer :: iout
        real(4) :: value
        logical :: isLittleEndian
        integer :: nbytes

*-----------------------------------------------------------------------

        integer :: isunit
        integer :: isunit_default = 100
        integer :: ios
        integer(kind=1) :: int1(4)
        integer :: i

*-----------------------------------------------------------------------

        if ( isLittleEndian ) then

          call open_file(isunit_default, "", isunit, ios, .true.)
          write(isunit) value

          read(isunit) (int1(i), i=1,4)
          write(iout)  (int1(i), i=4,1,-1)

          call close_file(isunit)

        else

          write(iout) value

        end if

        nbytes = 4
      end subroutine



************************************************************************
*                                                                      *
      subroutine write_real4s_be(
     &        iout, num, values, isLittleEndian, nbytes)
*                                                                      *
*       write real(4) values to binary file on big endian.             *
*                                                                      *
************************************************************************

        implicit none
        integer :: iout
        integer :: num
        real(4) :: values(num)
        logical :: isLittleEndian
        integer :: nbytes

*-----------------------------------------------------------------------

        integer :: isunit
        integer :: isunit_default = 100
        integer :: ios
        integer(kind=1) :: int1(4)
        integer :: i, n

*-----------------------------------------------------------------------

        if ( isLittleEndian ) then

          call open_file(isunit_default, "", isunit, ios, .true.)
          do n = 1, num
            write(isunit) values(n)
          end do

          rewind(isunit)

          do n = 1, num
            read(isunit) (int1(i), i=1,4)
            write(iout)  (int1(i), i=4,1,-1)
          end do

          call close_file(isunit)

        else

          do n = 1, num
            write(iout) values(n)
          end do

        end if

        nbytes = 4 * num
      end subroutine



************************************************************************
*                                                                      *
      subroutine write_real8s_as_real4s_be(
     &        iout, num, dvalues, isLittleEndian, nbytes)
*                                                                      *
*       convert real(8) values to real(4) and                          *
*       write real(4) values to binary file on big endian.             *
*                                                                      *
************************************************************************

        implicit none
        integer :: iout
        integer :: num
        real(8) :: dvalues(num)
        logical :: isLittleEndian
        integer :: nbytes

*-----------------------------------------------------------------------

        real(4) :: fvalues(num)
        integer :: i

*-----------------------------------------------------------------------

        do i = 1, num
          fvalues(i) = sngl(dvalues(i))
        end do

        call write_real4s_be(iout, num, fvalues, isLittleEndian, nbytes)

      end subroutine



************************************************************************
*                                                                      *
      subroutine write_int4s_be(
     &        iout, num, values, nbytes)
*                                                                      *
*       write integer(4) values to binary file on big endian.          *
*                                                                      *
************************************************************************

        implicit none
        integer :: iout
        integer :: num
        integer(4) :: values(num)
        integer :: nbytes

*-----------------------------------------------------------------------

        integer :: n

*-----------------------------------------------------------------------

        do n = 1, num
          call write_int4_be(iout, values(n))
        end do

        nbytes = 4 * num

      end subroutine



************************************************************************
*                                                                      *
      subroutine write_adjust_4bytes_pre(iout, ncsum, nbytes)
*                                                                      *
*       adjust binay alignment with 4*n-1.                             *
*                                                                      *
************************************************************************

        implicit none
        integer :: iout
        integer :: ncsum
        integer :: nbytes

*-----------------------------------------------------------------------

        integer(kind=1), parameter :: space = 32
        integer :: i, nadd

*-----------------------------------------------------------------------

        nadd = 3 - mod(ncsum, 4)

        do i = 1, nadd
          write(iout) space
        end do

        nbytes = nadd
      end subroutine


