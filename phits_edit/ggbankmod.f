************************************************************************
      module GGBANKMOD
************************************************************************
      integer,save:: nlaj_bank,nlcaj_bank

      integer,allocatable,save:: laj(:),lcaj(:)
!$OMP THREADPRIVATE(laj,lcaj)

      contains
!------------------------------------------------------------------------
      subroutine ALLOCATE_GGBANK
      implicit none

      allocate(laj(nlaj_bank),lcaj(nlcaj_bank))

      end subroutine ALLOCATE_GGBANK
!------------------------------------------------------------------------
      subroutine DEALLOCATE_GGBANK
      implicit none

      deallocate(laj,lcaj)

      end subroutine DEALLOCATE_GGBANK
!------------------------------------------------------------------------
      subroutine INIT_GGBANK
      use moddas_ggs !frtati20220905

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------
*     mark duplicate mentions of a surface as ambiguity surfaces.
*-----------------------------------------------------------------------
      laj  = 0
      lcaj = 0
      if(nlaj_bank .eq. 0) write(6,'(''** Warning: nlaj_bank = 0. Potent
     &ial memory corruption'')')

      do 400 i = 1, nlja+1
  400    lcaj(i) = 1 !FURUTA llaj+1
      do 440 ic = 1, mxa

      do 430 j = abs(lca(ic))+1,abs(lca(ic+1))-1

      do 410 jb = abs(lca(ic))+1,j

  410    if(abs(lja(jb-1)).eq.abs(lja(j))) goto 420
         goto 430
  420    lcaj(j) = -1 !FURUTA -llaj-1
  430 continue
  440 continue

*-----------------------------------------------------------------------
*     mark logical operators.  replicate lcaj for multitasking.
*-----------------------------------------------------------------------

      do 450 i = 1, nlja
  450    if(lja(i).gt.1000000) lcaj(i) = -1 !FURUTA -llaj-1

*-----------------------------------------------------------------------
*     flag cells that have one or more sq, gq, or tori or are void.
*-----------------------------------------------------------------------
cKN   this part is skipped and move to ggm01


         end subroutine INIT_GGBANK
!-----------------------------------------------------------------------
      end module GGBANKMOD
