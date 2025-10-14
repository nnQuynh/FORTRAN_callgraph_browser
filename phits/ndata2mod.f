************************************************************************
      module NDATA2MOD
*
*     Data for ndata=2 for t-yield, read data from yield in XS
*
************************************************************************
      implicit none
      integer,save:: ndatp, ndatn, ndatd, ndata, ndatg, 
     &               mmdatp, mmdatn, mmdatd, mmdata, mmdatg,
     &               iemmap, iemman, iemmad, iemmaa, iemmag,
     &               itallo

      real(8),allocatable,save:: cxproton(:,:,:,:),
     &                          cxneutron(:,:,:,:),
     &                         cxdeuteron(:,:,:,:),
     &                            cxalpha(:,:,:,:),
     &                            cxgamma(:,:,:,:)
      integer,allocatable,save::
     &     jdatp(:), jdatn(:), jdatd(:), jdata(:), jdatg(:),
     &     kndatp(:,:),kndatn(:,:),kndatd(:,:),kndata(:,:),kndatg(:,:),
     &     kldatp(:,:),kldatn(:,:),kldatd(:,:),kldata(:,:),kldatg(:,:),
     &     kedatp(:,:),kedatn(:,:),kedatd(:,:),kedata(:,:),kedatg(:,:),
     &     kidatp(:,:),kidatn(:,:),kidatd(:,:),kidata(:,:),kidatg(:,:)

      integer, save::  indatp(500), indatn(500), indatd(500),
     &     indata(500), indatg(500)

      contains
*-----------------------------------------------------------------------
      subroutine ALLOCATE_NDATA2

      allocate(jdatp(ndatp), jdatn(ndatn), jdatd(ndatd),
     &     jdata(ndata), jdatg(ndatg),
     &     kndatp(ndatp,mmdatp), kldatp(ndatp,mmdatp),
     &     kedatp(ndatp,mmdatp), kidatp(ndatp,mmdatp),
     &     kndatn(ndatn,mmdatn), kldatn(ndatn,mmdatn),
     &     kedatn(ndatn,mmdatn), kidatn(ndatn,mmdatn),
     &     kndatd(ndatd,mmdatd), kldatd(ndatd,mmdatd),
     &     kedatd(ndatd,mmdatd), kidatd(ndatd,mmdatd),
     &     kndata(ndata,mmdata), kldata(ndata,mmdata),
     &     kedata(ndata,mmdata), kidata(ndata,mmdata),
     &     kndatg(ndatg,mmdatg), kldatg(ndatg,mmdatg),
     &     kedatg(ndatg,mmdatg), kidatg(ndatg,mmdatg) )

      allocate(cxproton(ndatp,mmdatp,2,iemmap),
     &        cxneutron(ndatn,mmdatn,2,iemman),
     &       cxdeuteron(ndatd,mmdatd,2,iemmad),
     &          cxalpha(ndata,mmdata,2,iemmaa),
     &          cxgamma(ndatg,mmdatg,2,iemmag) )

      end subroutine ALLOCATE_NDATA2
*-----------------------------------------------------------------------
      subroutine DEALLOCATE_NDATA2

      deallocate( jdatp, jdatn, jdatd, jdata, jdatg,
     &     kndatp, kndatn, kndatd, kndata, kndatg,
     &     kldatp, kldatn, kldatd, kldata, kldatg,
     &     kedatp, kedatn, kedatd, kedata, kedatg,
     &     kidatp, kidatn, kidatd, kidata, kidatg )

      deallocate( cxproton, cxneutron, cxdeuteron, cxalpha, cxgamma )

      end subroutine DEALLOCATE_NDATA2
*-----------------------------------------------------------------------
      end module NDATA2MOD
