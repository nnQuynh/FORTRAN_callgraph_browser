************************************************************************
      module MMBANKMOD
************************************************************************
      use mod_ompparallel !--- NS 2020.04 del THREADPRIVATE
cKN 2015/03/19 nsos, nsosa, ibksos, iaksos
cFURUTA20160607 itetpos, itetposa, ibtetpos, ibtetposa
      implicit none
      integer,save:: maxbnk_bank
      integer,allocatable,save:: iblz(:,:),name(:,:),nmed(:,:),nty(:,:),
     &     nkf(:,:),nfcs(:,:),ncnt(:,:,:),nzst(:,:),nsos(:,:)
      real(8),allocatable,save:: wt(:,:),u(:,:),v(:,:),w(:,:),e(:,:),
     &          x(:,:),y(:,:),z(:,:),t(:,:),
     &          wtin(:,:),wtnz(:,:),xfcs(:,:)
      real(8),allocatable,save:: ec(:,:),xc(:,:),yc(:,:),zc(:,:),tc(:,:)
      real(8),allocatable,save:: spx(:,:),spy(:,:),spz(:,:)
      integer,allocatable,save:: iblza(:,:),namea(:,:),nmeda(:,:),
     &          ntya(:,:),
     &          nkfa(:,:),nfcsa(:,:),ncnta(:,:,:),nzsta(:,:),nsosa(:,:)
      real(8),allocatable,save:: wta(:,:),ua(:,:),va(:,:),wa(:,:),
     &          ea(:,:),
     &          xa(:,:),ya(:,:),za(:,:),ta(:,:),
     &          wtina(:,:),wtnza(:,:),xfcsa(:,:)
      real(8),allocatable,save:: spxa(:,:),spya(:,:),spza(:,:)
      integer,allocatable,save:: itetpos(:,:),itetposa(:,:)


      integer,save:: ibkblz,ibknam,ibknmd,ibknty,ibknkf,ibknfc,
     &     ibknct
      integer,save:: ibkwt,ibku,ibkv,ibkw,ibke,ibkx,ibky,ibkz,ibkt,
     &     ibkwin,ibkwnz,ibkxfc
      integer,save:: ibkspx,ibkspy,ibkspz,ibkzst,ibksos
      integer,save:: ibkec,ibkxc,ibkyc,ibkzc,ibktc
      integer,save:: iakblz,iaknam,iaknmd,iaknty,iaknkf,iaknfc,
     &     iaknct
      integer,save:: iakwt,iaku,iakv,iakw,iake,iakx,iaky,iakz,iakt,
     &     iakwin,iakwnz,iakxfc
      integer,save:: iakspx,iakspy,iakspz,iakzst,iaksos
      integer,save:: ibtetpos,ibtetposa

      contains
*-----------------------------------------------------------------------
      subroutine ALLOCATE_MMBANK
      implicit none
      integer max00

      max00=maxbnk_bank

!$OMP MASTER
      allocate( iblz(max00,npomp),name(max00,npomp),nmed(max00,npomp),
     &          nty(max00,npomp),
     &          nkf(max00,npomp),nfcs(max00,npomp),ncnt(3,max00,npomp),
     &          nzst(max00,npomp),
     &          nsos(max00,npomp) )
      allocate( wt(max00,npomp),u(max00,npomp),v(max00,npomp),
     &          w(max00,npomp),e(max00,npomp),
     &          x(max00,npomp),y(max00,npomp),z(max00,npomp),
     &          t(max00,npomp),
     &          wtin(max00,npomp),wtnz(max00,npomp),xfcs(max00,npomp) )
      allocate( ec(max00,npomp),xc(max00,npomp),yc(max00,npomp),
     &          zc(max00,npomp),tc(max00,npomp) )
      allocate( spx(max00,npomp),spy(max00,npomp),spz(max00,npomp) )

      allocate( iblza(max00,npomp),namea(max00,npomp),
     &          nmeda(max00,npomp),ntya(max00,npomp),
     &          nkfa(max00,npomp),nfcsa(max00,npomp),
     &          ncnta(3,max00,npomp),nzsta(max00,npomp),
     &          nsosa(max00,npomp) )
      allocate( wta(max00,npomp),ua(max00,npomp),va(max00,npomp),
     &          wa(max00,npomp),ea(max00,npomp),
     &          xa(max00,npomp),ya(max00,npomp),
     &          za(max00,npomp),ta(max00,npomp),
     &          wtina(max00,npomp),wtnza(max00,npomp),
     &          xfcsa(max00,npomp) )
      allocate( spxa(max00,npomp),spya(max00,npomp),spza(max00,npomp) )
      allocate( itetpos(max00,npomp),itetposa(max00,npomp) )
!$OMP END MASTER
!$OMP BARRIER

      end subroutine ALLOCATE_MMBANK
!-------------------------------------------------------------------------
      subroutine DEALLOCATE_MMBANK
!$OMP MASTER  !--- add NS 2020.04 del THREADPRIVATE
      deallocate(iblz,name,nmed,nty,nkf,nfcs,ncnt,nzst,nsos)
      deallocate(wt,u,v,w,e,x,y,z,t,wtin,wtnz,xfcs)
      deallocate(ec,xc,yc,zc,tc)
      deallocate(spx,spy,spz)
      deallocate(iblza,namea,nmeda,ntya,nkfa,nfcsa,ncnta,nzsta,nsosa)
      deallocate(wta,ua,va,wa,ea,xa,ya,za,ta,wtina,wtnza,xfcsa)
      deallocate(spxa,spya,spza)
      deallocate(itetpos,itetposa)
!$OMP END MASTER !--- add NS 2020.04 del THREADPRIVATE
!$OMP BARRIER !--- add NS 2020.04 del THREADPRIVATE
      end subroutine DEALLOCATE_MMBANK
!-------------------------------------------------------------------------

      end module MMBANKMOD
