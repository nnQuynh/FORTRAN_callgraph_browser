************************************************************************
      module GGMARRAYMOD
************************************************************************

      implicit none

      integer,save:: lrho,ltmp,lfme,ltth,lawc,lawt,lawn,ltbt,lden,lesa
      integer,save:: lpik,lgwt,ldpt,lrpt,leba,leee,lpru,lpkn,lfim,lrng
      integer,save:: ldrs,lqav,lasp,lear,lqcn,lflc,lx85,legg,ledg,leek
      integer,save:: lwwk,lpxr,lxnm,lpbr,lpbt,lebt,lebd,lfst,lftt,lech
      integer,save:: lrtc,ltgp,lemaxt,lpnt,lqav0,lrng0,lintm,lmat,llxd
      integer,save:: lnmt,ljem,lnsb,ljco,liza,lizn,llmn,lkmm,ljmd,lnpq
      integer,save:: lmbi,lkdr,llme,lkmt,lixl,lnty,lixc,lkxd,lkaw,lipa
      integer,save:: ljmt,ljxs,lnxs,llfc,llmt,lngm,lnht,lipt,lipb,lnpt
      integer,save:: ljpt,lktp,lkxs,lktc,ljss,lixs
      integer,save:: mbmemory ! size of the allocate memory for xxs and exs

      real(8),allocatable,save:: rho(:), tmp(:), fme(:), tth(:), awc(:)
      real(8),allocatable,save:: awt(:), awn(:), tbt(:), den(:), esa(:)
      real(8),allocatable,save:: pik(:), gwt(:), eee(:), pru(:), pkn(:)
      real(8),allocatable,save:: rng(:), drs(:), qav(:), asp(:), ear(:)
      real(8),allocatable,save:: qcn(:), flc(:), edg(:), eek(:), wwk(:)
      real(8),allocatable,save:: pxr(:), xnm(:), pbr(:), pbt(:), fst(:)
      real(8),allocatable,save:: ftt(:), emaxt(:), pnt(:), qav0(:)
      real(8),allocatable,save:: rng0(:), ech(:,:,:), ebd(:,:)
      real(8),allocatable,save:: ebt(:,:), egg(:,:), xse85(:,:)
      real(8),allocatable,save:: fim(:,:), eba(:,:), rptb(:), dptb(:,:)
      real(8),allocatable,save:: tds(:), xss(:), exs(:), tal(:)
      integer,allocatable,save:: mat(:), lxd(:,:), nmt(:), jemi(:)
      integer,allocatable,save:: nsb(:), jcond(:), iza(:), izn(:)
      integer,allocatable,save:: lmn(:), kmm(:), npq(:), kdr(:)
      integer,allocatable,save:: ipan(:), kaw(:), jmt(:), ngmfl(:)
      integer,allocatable,save:: nhtfl(:), ixs(:,:,:), ixl(:,:), kxs(:)
      integer,allocatable,save:: jmd(:), mbi(:), lme(:,:), kmt(:,:)
      integer,allocatable,save:: nty(:), ixc(:,:), kxd(:), jxs(:,:)
      integer,allocatable,save:: nxs(:,:), lfcl(:), lmt(:)
      integer,allocatable,save:: iptal(:,:,:), ktp(:,:), jptal(:,:)
      integer,allocatable,save:: nptb(:), iptb(:,:), itds(:)
      real(8),allocatable,save:: sigmxp(:,:),sigmxd(:,:),sigmxa(:,:)
      real(8),allocatable,save:: sigmxlb(:,:)

      contains
!------------------------------------------------------------------------

! T.Sato 2021/08/10, allocate exact number of dimension
      subroutine ALLOCATE_GGMARRAY(mxa,mtasks,mxt,mix,naw,mxe1,npikmt,
     &  kpt1,kpt2,nee,nmat1,istrg,ke,nbbrem1,mipt,indt,
     &  mpng,mwng,mtop,maxi,ipert,ntal,npert,npkey,mnnm)

      integer,intent(in):: mxa,mtasks,mxt,mix,naw,mxe1,npikmt,
     &  kpt1,kpt2,nee,nmat1,istrg,ke,nbbrem1,mipt,indt,
     &  mpng,mwng,mtop,maxi,ipert,ntal,npert,npkey,mnnm

      allocate(rho(mxa*mtasks), tmp(mxa*mxt), fme(mix*mtasks))
      allocate(tth(mxt), awc(mix), awt(naw))
      allocate(awn(mxe1), tbt(mxe1), den(mxa))
      allocate(esa(mxe1), pik(npikmt*mtasks), gwt(mxa*kpt1*kpt2))
      allocate(eee(nee), pru(nee * nmat1), pkn(nee * nmat1))
      allocate(rng(nee * nmat1), drs(nee * nmat1), qav(nee * nmat1))
      allocate(asp(nee * nmat1 * ( 1 - istrg )))
      allocate(ear(nee * nmat1 * ( 1 - istrg )))
      allocate(qcn(nee * nmat1 * ( 1 - istrg )))
      allocate(flc(nee * nmat1 * ( 1 - istrg )))
      allocate(edg(ke  * nmat1), eek(ke  * nmat1))
      allocate(wwk(ke  * nmat1), pxr(nee * nmat1), xnm(ke  * nmat1))
      allocate(pbr(nee * nmat1), pbt(nee * nmat1))
      allocate(fst(nee * nmat1 * nbbrem1))
      allocate(ftt(nee * nmat1 * nbbrem1))
      allocate(emaxt(mxe1))
      allocate(pnt(nmat1), qav0(nee * nmat1))
      allocate(rng0(nee * nmat1), mat(mxa))
      allocate(lxd(mipt,nmat1), nmt(nmat1 + 1))
      allocate(jemi(nmat1 * max(ke,1)), nsb(nmat1 * ke)) ! 2022/7/21 Ogawa Needed for ITSART
      allocate(jcond(nmat1 * ke), iza(mix))
      allocate(izn(mix), lmn(mix), kmm(mix))
      allocate(npq(0:nmat1+1), kdr(mxe1), kaw(naw)) !20220906frtati fbounds-check nmat1 -> 0:nmat1+1
      allocate(ipan(mxa + 1), jmt(indt))
      allocate(ngmfl(mxe1), nhtfl(mxe1))
      allocate(ech(mpng,mwng,nee/4+1+(nee/4+1)*(nmat1-1)))
      allocate(ebd(mtop, nee*nmat1), fim(mipt + 1, mxa))
      allocate(ebt(mtop, nee*nmat1), lme(mipt, mix))
      allocate(egg(maxi, nee*nmat1))
      allocate(xse85(10, nee*nmat1))
      allocate(eba(mtop, nee*nmat1), rptb(ipert))
      allocate(ixl(3, mxe1), kxs(mxe1))
      allocate(jmd(nmat1 + 2), mbi(nmat1 * ke))
      allocate(kmt(3, indt), nty(0:mxe1)) !20221004frtati fbounds-check 0:mxe1
      allocate(ixc(61,mxe1), kxd(mxe1))
      allocate(jxs(32, 0:mxe1), nxs(16, 0:mxe1), lfcl(mxa + 2)) !20221004frtati fbounds-check 0:mxe1
      allocate(lmt(mix), iptal(8,6,ntal), ktp(mipt,ntal))
      allocate(jptal(18, ntal), nptb(npert + 1))
      allocate(iptb((2+2*npkey), npert), dptb(3, npert*mnnm) )
      allocate(itds(ntal), tds(ntal))  ! T.Sato 2021/08/03

      lme = 0 ! 2021/8/31 Initialization to avoid out-of-range access.
      lmn = 0 ! frtati 2022/03/20
      npq = 0 !20220906frtati fbounds-check
      nty = 0 !20221004frtati fbounds-check
      nxs = 0 !20221004frtati fbounds-check
      jxs = 0 !20221004frtati fbounds-check
      allocate(sigmxp(2,mix),sigmxd(2,mix),sigmxa(2,mix))
      allocate(sigmxlb(2,mxe1))

      end subroutine ALLOCATE_GGMARRAY
!------------------------------------------------------------------------
      subroutine ALLOCATE_IXS(mixs,maxsec,mxe,nt)

      integer,intent(in):: mixs,maxsec,mxe,nt

      allocate(ixs(mixs,maxsec,mxe + nt))

      end subroutine ALLOCATE_IXS
!------------------------------------------------------------------------
      subroutine ALLOCATE_ggmTAL(m)

      integer,intent(in):: m

      if( .not. allocated(tal) ) then
            allocate(tal(m))
      else
            deallocate(tal)
            allocate(tal(m))
      end if

      end subroutine ALLOCATE_ggmTAL
!------------------------------------------------------------------------
      subroutine DEALLOCATE_GGMARRAY

      deallocate(rho, tmp, fme, tth, awc, awt, awn, tbt, den, esa, pik)
      deallocate(gwt, eee, pru, pkn, rng, drs, qav, asp, ear, qcn, flc)
      deallocate(edg, eek, wwk, pxr, xnm, pbr, pbt, fst, ftt, emaxt)
      deallocate(pnt, qav0, rng0, mat, lxd, nmt, jemi, nsb, jcond, iza)
      deallocate(izn, lmn, kmm, npq, kdr, ipan, kaw, jmt, ngmfl, nhtfl)
      deallocate(ech, ebd, ebt, egg, xse85, fim, eba, rptb, ixs, ixl)
      deallocate(kxs, jmd, mbi, lme, kmt, nty, ixc, kxd, jxs, nxs)
      deallocate(lfcl, lmt, iptal, ktp, jptal, iptb, dptb, itds, tds)
      deallocate(tal)
      deallocate(sigmxp,sigmxd,sigmxa,sigmxlb)

      end subroutine DEALLOCATE_GGMARRAY
!------------------------------------------------------------------------
      end module GGMARRAYMOD
