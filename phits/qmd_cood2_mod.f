************************************************************************
*                                                                      *
      module qmd_cood2_mod
*                                                                      *
*                                                                      *
*        Last Revised:     2017 06 14  by T.Ogawa                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to allocate long arrays used in JQMD                    *
*              Length is 800 for JQMD, 13000 for JAMQMD                *
*                                                                      *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)

* Level energy
      logical lfirst
!$OMP THREADPRIVATE(lfirst)

      integer, save :: nqmdm
!$OMP THREADPRIVATE(nqmdm)

      double precision, allocatable :: rr2(:,:), rbij(:,:), pp2(:,:)
      double precision, allocatable :: rha(:,:), rhe(:,:), rhc(:,:)
!$OMP THREADPRIVATE(rr2,rbij,pp2,rha,rhe,rhc)

      double precision, allocatable, save :: ffr(:,:), ffp(:,:), rh3d(:)
!$OMP THREADPRIVATE(ffr, ffp, rh3d)

      double precision, allocatable, save :: qqr(:,:), qqp(:,:),
     & ggr(:,:), ccr(:,:), ggp(:,:), ccp(:,:) !, ffr(:,:), ffp(:,:)
!$OMP THREADPRIVATE(qqr, qqp, ggr, ccr, ggp, ccp) ! , ffr, ffp

      integer, allocatable, save          :: mascl(:), num(:)
      integer, allocatable, save          :: isort(:), isorti(:), it(:)
      double precision, allocatable, save :: rs(:,:), ps(:,:) !, rhoa(:)
!$OMP THREADPRIVATE(mascl, num, isort, isorti, it, rs, ps)

      double precision, allocatable, save :: rhol(:), dpot(:)
      double precision, allocatable, save :: rhoa(:), rhos(:), rhoc(:)
      double precision, allocatable, save :: phase0(:), phase(:)
!$OMP THREADPRIVATE(rhol, dpot, rhoa, rhos, rhoc)
!$OMP THREADPRIVATE(phase0, phase)

      integer, allocatable, save :: inun0(:),iavd0(:),ihis0(:)
      integer, allocatable, save :: ichg0(:),inuc0(:),ibry0(:),inds0(:)
!$OMP THREADPRIVATE(inun0,iavd0,ihis0,ichg0,inuc0,ibry0,inds0)

      double precision, allocatable, save :: f0r(:,:), f0p(:,:),
     & d1r(:,:), d1p(:,:)
!$OMP THREADPRIVATE(d1r, d1p, f0r, f0p)

      double precision, allocatable, save :: r0(:,:), p0(:,:)
!$OMP THREADPRIVATE(r0, p0)

************  Temporary storage for reallocation ***********************
      double precision, allocatable, private :: rr2s(:,:), rbijs(:,:),
     & pp2s(:,:)
      double precision, allocatable, private :: rhas(:,:), rhes(:,:),
     & rhcs(:,:)
!$OMP THREADPRIVATE(rr2s,rbijs,pp2s,rhas,rhes,rhcs)

      double precision, allocatable, private :: ffrs(:,:), ffps(:,:),
     & rh3ds(:)
!$OMP THREADPRIVATE(ffrs, ffps, rh3ds)

      double precision, allocatable, private :: qqrs(:,:), qqps(:,:),
     & ggrs(:,:), ccrs(:,:), ggps(:,:), ccps(:,:)
!$OMP THREADPRIVATE(qqrs, qqps, ggrs, ccrs, ggps, ccps)

      integer, allocatable, private :: mascls(:), nums(:)
      integer, allocatable, private :: isorts(:), isortis(:), its(:)
      double precision, allocatable, private :: rss(:,:), pss(:,:)
!$OMP THREADPRIVATE(mascls, nums, isorts, isortis, its, rss, pss)

      double precision, allocatable, private :: rhols(:), dpots(:)
      double precision, allocatable, private :: rhoas(:), rhoss(:)
     & , rhocs(:)
      double precision, allocatable, private :: phase0s(:), phases(:)
!$OMP THREADPRIVATE(rhols, dpots, rhoas, rhoss, rhocs)
!$OMP THREADPRIVATE(phase0s, phases)

      integer, allocatable, private :: inun0s(:),iavd0s(:), ihis0s(:)
      integer, allocatable, private :: ichg0s(:),inuc0s(:), ibry0s(:),
     & inds0s(:)
!$OMP THREADPRIVATE(inun0s,iavd0s,ihis0s,ichg0s,inuc0s,ibry0s,inds0s)

      double precision, allocatable, private :: f0rs(:,:), f0ps(:,:),
     & d1rs(:,:), d1ps(:,:)
!$OMP THREADPRIVATE(d1rs, d1ps, f0rs, f0ps)

      double precision, allocatable, private :: r0s(:,:), p0s(:,:)
!$OMP THREADPRIVATE(r0s, p0s)
************************************************************************

      data lfirst / .false. /
      data nqmdm  / 800/

      contains

************************************************************************
*                                                                      *
      subroutine coodalloc
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)

       n = 800

       allocate(rr2(n,n),rbij(n,n),pp2(n,n),rha(n,n), rhe(n,n),rhc(n,n),
     &  ffr(3,n), ffp(3,n), rh3d(n),
     &  f0r(3,n), d1r(3,n), f0p(3,n), d1p(3,n),
     &  qqr(3,n), qqp(3,n), ggr(3,n), ccr(3,n), !ffr(3,n), ffp(3,n),
     &  ggp(3,n), ccp(3,n),
     &  mascl(n), num(n), isort(n), isorti(n), it(0:n), ! rhoa(n),
     &  rs(3,n),  ps(3,n),
     &  rhol(n),  dpot(n), rhoa(n), rhos(n), rhoc(n), phase0(n),
     &  phase(n), ! , d1r(3,n), d1p(3,n)
     &  r0(5,n),  p0(6,n), ichg0(n), inuc0(n), ibry0(n),
     &  inds0(n), inun0(n), iavd0(n), ihis0(n))

      return
      end subroutine
************************************************************************

************************************************************************
*                                                                      *
      subroutine coodrealloc
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)

      n = nqmdm
       allocate(rr2s(n,n),rbijs(n,n),pp2s(n,n),rhas(n,n), rhes(n,n),
     &  rhcs(n,n), ffrs(3,n), ffps(3,n), rh3ds(n),
     &  f0rs(3,n), d1rs(3,n), f0ps(3,n), d1ps(3,n),
     &  qqrs(3,n), qqps(3,n), ggrs(3,n), ccrs(3,n),
     &  ggps(3,n), ccps(3,n),
     &  mascls(n), nums(n), isorts(n), isortis(n), its(0:n),
     &  rss(3,n),  pss(3,n),
     &  rhols(n),  dpots(n), rhoas(n), rhoss(n), rhocs(n), phase0s(n),
     &  phases(n),
     &  r0s(5,n),  p0s(6,n), ichg0s(n), inuc0s(n), ibry0s(n),
     &  inds0s(n), inun0s(n), iavd0s(n), ihis0s(n))

       rr2s    = rr2
       rbijs   = rbij
       pp2s    = pp2
       rhas    = rha
       rhes    = rhe
       rhcs    = rhc
       ffrs    = ffr
       ffps    = ffp
       rh3ds   = rh3d
       f0rs    = f0r
       d1rs    = d1r
       f0ps    = f0p
       d1ps    = d1p
       qqrs    = qqr
       qqps    = qqp
       ggrs    = ggr
       ccrs    = ccr
       ggps    = ggp
       ccps    = ccp
       mascls  = mascl
       nums    = num
       isorts  = isort
       isortis = isorti
       its     = it
       rss     = rs
       pss     = ps
       rhols   = rhol
       dpots   = dpot
       rhoas   = rhoa
       rhoss   = rhos
       rhocs   = rhoc
       phase0s = phase0
       phases  = phase
       r0s     = r0
       p0s     = p0
       ichg0s  = ichg0
       inuc0s  = inuc0
       ibry0s  = ibry0
       inds0s  = inds0
       inun0s  = inun0
       iavd0s  = iavd0
       ihis0s  = ihis0

        deallocate(rr2,rbij,pp2,rha, rhe, rhc,
     &  ffr, ffp, rh3d,
     &  f0r, d1r, f0p, d1p,
     &  qqr, qqp, ggr, ccr,
     &  ggp, ccp,
     &  mascl, num, isort, isorti, it,
     &  rs,  ps,
     &  rhol,  dpot, rhoa, rhos, rhoc, phase0, phase,
     &  r0,  p0, ichg0, inuc0, ibry0,
     &  inds0, inun0, iavd0, ihis0)

      m = nqmdm
      nqmdm = nqmdm * 2
      if(nqmdm .gt. 12800) write(*,*) 'Warning. JAMQMD array too big'
      n = nqmdm

       allocate(rr2(n,n),rbij(n,n),pp2(n,n),rha(n,n), rhe(n,n),rhc(n,n),
     &  ffr(3,n), ffp(3,n), rh3d(n),
     &  f0r(3,n), d1r(3,n), f0p(3,n), d1p(3,n),
     &  qqr(3,n), qqp(3,n), ggr(3,n), ccr(3,n), !ffr(3,n), ffp(3,n),
     &  ggp(3,n), ccp(3,n),
     &  mascl(n), num(n), isort(n), isorti(n), it(0:n), ! rhoa(n),
     &  rs(3,n),  ps(3,n),
     &  rhol(n),  dpot(n), rhoa(n), rhos(n), rhoc(n), phase0(n),
     &  phase(n), ! , d1r(3,n), d1p(3,n)
     &  r0(5,n),  p0(6,n), ichg0(n), inuc0(n), ibry0(n),
     &  inds0(n), inun0(n), iavd0(n), ihis0(n))

       rr2    = 0.d0
       rbij   = 0.d0
       pp2    = 0.d0
       rha    = 0.d0
       rhe    = 0.d0
       rhc    = 0.d0
       ffr    = 0.d0
       ffp    = 0.d0
       rh3d   = 0.d0
       f0r    = 0.d0
       d1r    = 0.d0
       f0p    = 0.d0
       d1p    = 0.d0
       qqr    = 0.d0
       qqp    = 0.d0
       ggr    = 0.d0
       ccr    = 0.d0
       ggp    = 0.d0
       ccp    = 0.d0
       mascl  = 0
       num    = 0
       isort  = 0
       isorti = 0
       it     = 0
       rs     = 0.d0
       ps     = 0.d0
       rhol   = 0.d0
       dpot   = 0.d0
       rhoa   = 0.d0
       rhos   = 0.d0
       rhoc   = 0.d0
       phase0 = 0.d0
       phase  = 0.d0
       r0     = 0.d0
       p0     = 0.d0
       ichg0  = 0
       inuc0  = 0
       ibry0  = 0
       inds0  = 0
       inun0  = 0
       iavd0  = 0
       ihis0  = 0

      do iii = 1, m
       do ii = 1, m
       rr2(iii,ii)    = rr2s(iii,ii)
       rbij(iii,ii)   = rbijs(iii,ii)
       pp2(iii,ii)    = pp2s(iii,ii)
       rha(iii,ii)    = rhas(iii,ii)
       rhe(iii,ii)    = rhes(iii,ii)
       rhc(iii,ii)    = rhcs(iii,ii)
       enddo
       do ii = 1, 3
       ffr(ii,iii)    = ffrs(ii,iii)
       ffp(ii,iii)    = ffps(ii,iii)
       f0r(ii,iii)    = f0rs(ii,iii)
       d1r(ii,iii)    = d1rs(ii,iii)
       f0p(ii,iii)    = f0ps(ii,iii)
       d1p(ii,iii)    = d1ps(ii,iii)
       qqr(ii,iii)    = qqrs(ii,iii)
       qqp(ii,iii)    = qqps(ii,iii)
       ggr(ii,iii)    = ggrs(ii,iii)
       ccr(ii,iii)    = ccrs(ii,iii)
       ggp(ii,iii)    = ggps(ii,iii)
       ccp(ii,iii)    = ccps(ii,iii)
       rs(ii,iii)     = rss(ii,iii)
       ps(ii,iii)     = pss(ii,iii)
       enddo
       do ii = 1, 5
       r0(ii,iii)     = r0s(ii,iii)
       enddo
       do ii = 1, 6
       p0(ii,iii)     = p0s(ii,iii)
       enddo
       rh3d(iii)   = rh3ds(iii)
       mascl(iii)  = mascls(iii)
       num(iii)    = nums(iii)
       isort(iii)  = isorts(iii)
       isorti(iii) = isortis(iii)
       it(iii)     = its(iii)
       rhol(iii)   = rhols(iii)
       dpot(iii)   = dpots(iii)
       rhoa(iii)   = rhoas(iii)
       rhos(iii)   = rhoss(iii)
       rhoc(iii)   = rhocs(iii)
       phase0(iii) = phase0s(iii)
       phase(iii)  = phases(iii)
       ichg0(iii)  = ichg0s(iii)
       inuc0(iii)  = inuc0s(iii)
       ibry0(iii)  = ibry0s(iii)
       inds0(iii)  = inds0s(iii)
       inun0(iii)  = inun0s(iii)
       iavd0(iii)  = iavd0s(iii)
       ihis0(iii)  = ihis0s(iii)
      enddo

        deallocate(rr2s,rbijs,pp2s,rhas, rhes, rhcs,
     &  ffrs, ffps, rh3ds,
     &  f0rs, d1rs, f0ps, d1ps,
     &  qqrs, qqps, ggrs, ccrs,
     &  ggps, ccps,
     &  mascls, nums, isorts, isortis, its,
     &  rss,  pss,
     &  rhols,  dpots, rhoas, rhoss, rhocs, phase0s, phases,
     &  r0s,  p0s, ichg0s, inuc0s, ibry0s,
     &  inds0s, inun0s, iavd0s, ihis0s)

      return
      end subroutine
************************************************************************

************************************************************************
*                                                                      *
      subroutine cooddealloc
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)

      if(allocated(rr2)) then
       deallocate(rr2,rbij,pp2,rha, rhe, rhc, ffr, ffp, rh3d, f0r, d1r,
     & f0p, d1p, qqr, qqp, ggr, ccr, ggp, ccp, mascl, num, isort,
     & isorti, it, rs,  ps, rhol, dpot, rhoa, rhos, rhoc, phase0, phase,
     & r0,  p0, ichg0, inuc0, ibry0, inds0, inun0, iavd0, ihis0)
      endif

      return
      end subroutine
************************************************************************

      end module QMD_COOD2_MOD
************************************************************************
