************************************************************************
      module COSMICMOD
************************************************************************
      implicit double precision (a-h,o-z)

      real*8,allocatable,save :: angtable(:,:,:),angtable2(:,:,:)
      real*8,allocatable,save :: amidcr(:),ahighcr(:)

      contains
!------------------------------------------------------------------------
      subroutine CosmicAngSetup(j)

      use moddas_source
      implicit double precision (a-h,o-z)
      include 'param.inc'

      parameter(npart=33) ! number of applicable particle
      parameter(nabin=100)
      dimension IangPart(0:npart)
      data IangPart/1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &              0,0,0,0,0,0,0,0,0,4,4,5,5,6/ ! Angular particle ID
      data amin/-1.0d0/  ! minimum angle should be fixed
      data amax/1.0d0/   ! maximum angle should be fixed

      common /isorse/ ngrp(isrc), ngei(isrc), ngea(isrc), ngfe(isrc),
     &                ngft(isrc), ngll(isrc), ngpi(isrc), ngpw(isrc)

      common /cosmicint/ipcosmic(isrc),icenv(isrc)
      common /cosmicreal/solarmod(isrc),rigid(isrc),
     &        depatom(isrc),ground(isrc),environ(isrc)

      data nebin/0/
      save astep


      if(j.eq.1) then ! first time call here
       allocate(ahighcr(0:nabin))
       allocate(amidcr(nabin))
       astep=(amax-amin)/nabin
       do ia=0,nabin
        ahighcr(ia)=amin+astep*ia
        if(ia.ne.0) amidcr(ia)=(ahighcr(ia)+ahighcr(ia-1))*0.5
       enddo
      endif

      if(ngrp(j).gt.nebin) nebin=ngrp(j)

      if(j.eq.1) then
       allocate(angtable(-1:nabin,nebin,j))
       angtable(:,:,:)=0.0 ! initialization
      else
       allocate(angtable2(-1:nabin,nebin,j-1)) ! temporary dimension for storing angtable
       do jj=1,j-1
        do ie=1,nebin
         do ia=-1,nabin
          angtable2(ia,ie,jj)=angtable(ia,ie,jj)
         enddo
        enddo
       enddo
       deallocate(angtable)
       do i=1,100 ! wait for a while to avoid allocation error
       enddo
       allocate(angtable(-1:nabin,nebin,j))
       angtable(:,:,:)=0.0 ! initialization
       do jj=1,j-1
        do ie=1,nebin
         do ia=-1,nabin
          angtable(ia,ie,jj)=angtable2(ia,ie,jj)
         enddo
        enddo
       enddo
       deallocate(angtable2)
      endif

      ip=ipcosmic(j)
      s=solarmod(j)
      r=rigid(j)
      d=depatom(j)
      g=ground(j)
      if(IangPart(ip).ne.0) then
       do ie=1,nebin
        e=(egmin(ngei(j)+ie)+egmax(ngea(j)+ie))*0.5d0
        do ia=1,nabin
         angtable(ia,ie,j)=angtable(ia-1,ie,j)
     &   +getSpecAngFinal(iangpart(ip),s,r,d,e,g,amidcr(ia))  ! angular integrated value
        enddo
        angtable(-1,ie,j)=angtable(nabin,ie,j)*astep*2.0d0*acos(-1.0d0)  ! T.Sato 2021/12/01 for muon and neutron correction case
       enddo
      else ! for heavy ions, consider Earth's shadow effect
       costmp=-getshadowcos(d) ! cos(theta) =< 0
       do ie=1,nebin
        do ia=1,nabin
         if(ahighcr(ia-1).ge.costmp) angtable(ia,ie,j)=
     &   angtable(ia-1,ie,j)+1.0 ! uniform above threshold angle
        enddo
        angtable(-1,ie,j)=1.0d0 ! always 1.0
       enddo
      endif

      do ie=1,nebin ! normalization
       do ia=1,nabin
        angtable(ia,ie,j)=angtable(ia,ie,j)/angtable(nabin,ie,j) ! normalized to 1
       enddo
      enddo

      return
      end subroutine CosmicAngSetup

!------------------------------------------------------------------------
      function getCosmicAng(j,ie)
      implicit double precision (a-h,o-z)
      parameter(nabin=100)

      random = unirn(dummy)

      do ia=1,nabin
       if(random.le.angtable(ia,ie,j)) exit
      enddo

      getCosmicAng=ahighcr(ia-1)
     & +(ahighcr(ia)-ahighcr(ia-1))*unirn(dummy)

      return

      end function getCosmicAng

!------------------------------------------------------------------------
      subroutine DEALLOCATE_angtable
      implicit double precision (a-h,o-z)
      if(allocated(angtable))deallocate(angtable)
      if(allocated(amidcr))deallocate(amidcr)
      if(allocated(ahighcr))deallocate(ahighcr)

      return
      end subroutine DEALLOCATE_angtable
!------------------------------------------------------------------------
      end module COSMICMOD
