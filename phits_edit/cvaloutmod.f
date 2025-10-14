      module cvaloutmod

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      real(8),allocatable,save,target:: cvalout(:,:)
      character,allocatable,save,target:: cfncvalout(:)*200
      integer imemorycval
      data imemorycval/0/
      save imemorycval
      integer icvalout
      data icvalout/0/
      save icvalout

*-----------------------------------------------------------------------
      contains

      subroutine allocate_cvalout

      implicit real*8 (a-h,o-z)

      common /paraj/ mstz(300), parz(300)
      integer,parameter:: ndatacvalout = 4

      allocate(cvalout(ndatacvalout,iabs(mstz(161))))
      allocate(cfncvalout(iabs(mstz(161))))

      end subroutine


      subroutine deallocate_cvalout

      implicit real*8 (a-h,o-z)

      deallocate(cvalout)
      deallocate(cfncvalout)
      imemorycval = 0

      end subroutine deallocate_cvalout


      end module cvaloutmod
