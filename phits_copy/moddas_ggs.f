!module moddas_ggs
      module moddas_ggs

      double precision, allocatable :: trf(:,:), scf(:)
      double precision, allocatable :: vcl(:,:,:)
      integer, allocatable :: lsc(:), lca(:), ncl(:), nsfm(:),
     &lja(:), jtr(:), kst(:), ksu(:), ksc(:), nlv(:), idnt(:)
      integer, allocatable :: jun(:), mazp(:,:), mfl(:,:),
     &mazu(:), lat(:,:), ktr(:), idna(:)
      integer, allocatable :: ksm(:)
      integer, allocatable :: idne(:), idns(:)
      integer, allocatable :: lse(:)

!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_deallocate_ggs

      deallocate( trf, scf, vcl, lsc, lca, ncl, nsfm, lja, jtr, kst,
     &ksu, ksc, nlv, idnt, jun, mazp, mfl, mazu, lat, ktr, idna, ksm,
     &idne, idns, lse )

      end subroutine
*-----------------------------------------------------------------------
      end module
