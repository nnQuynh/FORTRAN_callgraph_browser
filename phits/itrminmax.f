
      subroutine setitrmin(itrmin,ib,ie,iarray)

        implicit none

        integer,intent(in) :: ib,ie
        integer,intent(out):: itrmin(1:ie)
        integer,intent(in) :: iarray(ib:ie)

        integer i

        do i = ib,ie
          itrmin(i) = min(itrmin(i), iarray(i))
        enddo

      end subroutine

      subroutine setitrmax(itrmax,ib,ie,iarray)

        implicit none

        integer,intent(in) :: ib,ie
        integer,intent(out):: itrmax(1:ie)
        integer,intent(in) :: iarray(ib:ie)

        integer i

        do i = ib,ie
          itrmax(i) = max(itrmax(i), iarray(i))
        enddo

      end subroutine

      subroutine setitrminmax(itrmin,itrmax,ib,ie,iarray)

        implicit none

        integer,intent(in) :: ib,ie
        integer,intent(out):: itrmin(1:ie)
        integer,intent(out):: itrmax(1:ie)
        integer,intent(in) :: iarray(ib:ie)

        integer i

        do i = ib,ie
          itrmin(i) = min(itrmin(i), iarray(i))
          itrmax(i) = max(itrmax(i), iarray(i))
        enddo

      end subroutine

      subroutine resetitrminmax(itrmin,itrmax,m,iarray)

        implicit none

        integer,intent(in) :: m
        integer,intent(out):: itrmin(m),itrmax(m)
        integer,intent(in) :: iarray(m)

        itrmax(:) = 0
        itrmin(:) = iarray(:)

      end subroutine

      subroutine readitrminmax3(itrmin,itrmax,itrmindef,
     &                          mn1,mn2,mn3,
     &                          mx1,mx2,mx3)

        implicit none

        integer,intent(inout) :: itrmin(3),itrmax(3)
        integer,intent(in) :: itrmindef(3)
        integer,intent(out):: mn1,mn2,mn3
        integer,intent(out):: mx1,mx2,mx3

        if (itrmin(1) .eq. 0) itrmin(:) = itrmindef(:)

        mn1 = itrmin(1)
        mn2 = itrmin(2)
        mn3 = itrmin(3)

        mx1 = itrmax(1)
        mx2 = itrmax(2)
        mx3 = itrmax(3)

      end subroutine

      subroutine readitrminmax4(itrmin,itrmax,itrmindef,
     &                          mn1,mn2,mn3,mn4,
     &                          mx1,mx2,mx3,mx4)

        implicit none

        integer,intent(inout) :: itrmin(4),itrmax(4)
        integer,intent(in) :: itrmindef(4)
        integer,intent(out):: mn1,mn2,mn3,mn4
        integer,intent(out):: mx1,mx2,mx3,mx4

        if (itrmin(1) .eq. 0) itrmin(:) = itrmindef(:)

        mn1 = itrmin(1)
        mn2 = itrmin(2)
        mn3 = itrmin(3)
        mn4 = itrmin(4)

        mx1 = itrmax(1)
        mx2 = itrmax(2)
        mx3 = itrmax(3)
        mx4 = itrmax(4)

      end subroutine

      subroutine readitrminmax5(itrmin,itrmax,itrmindef,
     &                          mn1,mn2,mn3,mn4,mn5,
     &                          mx1,mx2,mx3,mx4,mx5)

        implicit none

        integer,intent(inout) :: itrmin(5),itrmax(5)
        integer,intent(in) :: itrmindef(5)
        integer,intent(out):: mn1,mn2,mn3,mn4,mn5
        integer,intent(out):: mx1,mx2,mx3,mx4,mx5

        if (itrmin(1) .eq. 0) itrmin(:) = itrmindef(:)

        mn1 = itrmin(1)
        mn2 = itrmin(2)
        mn3 = itrmin(3)
        mn4 = itrmin(4)
        mn5 = itrmin(5)

        mx1 = itrmax(1)
        mx2 = itrmax(2)
        mx3 = itrmax(3)
        mx4 = itrmax(4)
        mx5 = itrmax(5)

      end subroutine

      subroutine readitrminmax6(itrmin,itrmax,itrmindef,
     &                          mn1,mn2,mn3,mn4,mn5,mn6,
     &                          mx1,mx2,mx3,mx4,mx5,mx6)

        implicit none

        integer,intent(inout) :: itrmin(6),itrmax(6)
        integer,intent(in) :: itrmindef(6)
        integer,intent(out):: mn1,mn2,mn3,mn4,mn5,mn6
        integer,intent(out):: mx1,mx2,mx3,mx4,mx5,mx6

        if (itrmin(1) .eq. 0) itrmin(:) = itrmindef(:)

        mn1 = itrmin(1)
        mn2 = itrmin(2)
        mn3 = itrmin(3)
        mn4 = itrmin(4)
        mn5 = itrmin(5)
        mn6 = itrmin(6)

        mx1 = itrmax(1)
        mx2 = itrmax(2)
        mx3 = itrmax(3)
        mx4 = itrmax(4)
        mx5 = itrmax(5)
        mx6 = itrmax(6)

      end subroutine

      subroutine readitrminmax7(itrmin,itrmax,itrmindef,
     &                          mn1,mn2,mn3,mn4,mn5,mn6,mn7,
     &                          mx1,mx2,mx3,mx4,mx5,mx6,mx7)

        implicit none

        integer,intent(inout) :: itrmin(7),itrmax(7)
        integer,intent(in) :: itrmindef(7)
        integer,intent(out):: mn1,mn2,mn3,mn4,mn5,mn6,mn7
        integer,intent(out):: mx1,mx2,mx3,mx4,mx5,mx6,mx7

        if (itrmin(1) .eq. 0) itrmin(:) = itrmindef(:)

        mn1 = itrmin(1)
        mn2 = itrmin(2)
        mn3 = itrmin(3)
        mn4 = itrmin(4)
        mn5 = itrmin(5)
        mn6 = itrmin(6)
        mn7 = itrmin(7)

        mx1 = itrmax(1)
        mx2 = itrmax(2)
        mx3 = itrmax(3)
        mx4 = itrmax(4)
        mx5 = itrmax(5)
        mx6 = itrmax(6)
        mx7 = itrmax(7)

      end subroutine

      subroutine readitrminmax8(itrmin,itrmax,itrmindef,
     &                          mn1,mn2,mn3,mn4,mn5,mn6,mn7,mn8,
     &                          mx1,mx2,mx3,mx4,mx5,mx6,mx7,mx8)

        implicit none

        integer,intent(inout) :: itrmin(8),itrmax(8)
        integer,intent(in) :: itrmindef(8)
        integer,intent(out):: mn1,mn2,mn3,mn4,mn5,mn6,mn7,mn8
        integer,intent(out):: mx1,mx2,mx3,mx4,mx5,mx6,mx7,mx8

        if (itrmin(1) .eq. 0) itrmin(:) = itrmindef(:)

        mn1 = itrmin(1)
        mn2 = itrmin(2)
        mn3 = itrmin(3)
        mn4 = itrmin(4)
        mn5 = itrmin(5)
        mn6 = itrmin(6)
        mn7 = itrmin(7)
        mn8 = itrmin(8)

        mx1 = itrmax(1)
        mx2 = itrmax(2)
        mx3 = itrmax(3)
        mx4 = itrmax(4)
        mx5 = itrmax(5)
        mx6 = itrmax(6)
        mx7 = itrmax(7)
        mx8 = itrmax(8)

      end subroutine

