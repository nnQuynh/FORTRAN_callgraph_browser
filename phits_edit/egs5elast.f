
      subroutine egs5elast(negs)

      implicit none
!      save
      include 'include/egs5_h.f'
      include 'include/egs5_stack.f'

      integer negs

      if(deresid .eq. 0) then
         negs = 0
      else
         negs = 1
      endif

      end
