************************************************************************
      module TSEDMOD
************************************************************************
      parameter(mparased=9)

      implicit double precision (a-h,o-z)

      real*8,allocatable,save :: eincele(:),cdiamele(:)
      real*8,allocatable,save :: eincstp(:),cdiamstp(:)
      real*8,allocatable,save :: eincion(:),cdiamion(:)
      real*8,allocatable,save :: paraele(:,:,:,:,:)
      real*8,allocatable,save :: parastp(:,:,:,:,:)
      real*8,allocatable,save :: paraion(:,:,:,:,:)

      integer,allocatable,save :: izion(:)
      real*8,allocatable,save :: xmid(:),xwid(:)
      real*8,save :: Apara(mparased)
!$OMP THREADPRIVATE( Apara )
      save npartion,neincion,mcdiamion
      save neincele,mcdiamele
      save neincstp,mcdiamstp
      save nxbinint
      save nauger
      save ip1,ie1,ic1,iaustp,imodel
!$OMP THREADPRIVATE( ip1,ie1,ic1,iaustp,imodel )
      save ratiop,ratioe,ratioc
!$OMP THREADPRIVATE( ratiop,ratioe,ratioc )

      contains
!------------------------------------------------------------------------
      subroutine SetupSED
      implicit double precision (a-h,o-z)

      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn

      data ifirst/0/
      data modelmax/1/ ! maximum number of material

      if(ifirst.eq.0) then ! necessary to read only once
       ifirst=1
       open(100,file=chfn(1)(1:ilfn(1))//'/data/tsed.dat',status='old')
!  read electron parameters
       read(100,*) nauger
       read(100,*) neincele
       read(100,*) mcdiamele
       allocate(eincele(neincele),cdiamele(mcdiamele))
       allocate(paraele(mcdiamele,nauger,neincele,mparased,modelmax))
       read(100,*) (eincele(ie),ie=1,neincele)
       read(100,*) (cdiamele(ic),ic=1,mcdiamele)
       do im=1,modelmax
        do ia=1,nauger
         do ie=1,neincele
          do ic=1,mcdiamele
           read(100,*) (paraele(ic,ia,ie,i,im),i=1,mparased)
          enddo
         enddo
        enddo
       enddo
!  read stopped-electron parameters
       read(100,*) nauger ! should be the same for normal electron
       read(100,*) neincstp
       read(100,*) mcdiamstp
       allocate(eincstp(neincstp),cdiamstp(mcdiamstp))
       allocate(parastp(mcdiamstp,nauger,neincstp,mparased,modelmax))
       read(100,*) (eincstp(ie),ie=1,neincstp)
       read(100,*) (cdiamstp(ic),ic=1,mcdiamstp)
       do im=1,modelmax
        do ia=1,nauger
         do ie=1,neincstp
          do ic=1,mcdiamstp
           read(100,*) (parastp(ic,ia,ie,i,im),i=1,mparased)
          enddo
         enddo
        enddo
       enddo
!  read ion parameters
       read(100,*) neincion
       read(100,*) mcdiamion
       read(100,*) npartion
       allocate(eincion(neincion),cdiamion(mcdiamion))
       allocate(paraion(mcdiamion,npartion,neincion,mparased,modelmax))
       allocate(izion(npartion))
       read(100,*) (eincion(ie),ie=1,neincion)
       read(100,*) (cdiamion(ic),ic=1,mcdiamion)
       read(100,*) (izion(ip),ip=1,npartion)
       do im=1,modelmax
        do ip=1,npartion
         do ie=1,neincion
          do ic=1,mcdiamion
           read(100,*) (paraion(ic,ip,ie,i,im),i=1,mparased)
          enddo
         enddo
        enddo
       enddo
       close(100)
! determine z bin for integration
       nxbinint=120
       allocate(xmid(nxbinint),xwid(nxbinint))
       xwid(1)=1.0d0
       xmid(1)=1.0d0
       do iz=2,nxbinint
        xwid(iz)=xwid(iz-1)*1.1d0
        xmid(iz)=xmid(iz-1)+(xwid(iz-1)+xwid(iz))*0.5d0
       enddo
      endif

      end subroutine SetupSED

!------------------------------------------------------------------------
      subroutine DEALLOCATE_SED
      implicit double precision (a-h,o-z)
      if(allocated(paraele))deallocate(paraele)
      if(allocated(parastp))deallocate(parastp)
      if(allocated(paraion))deallocate(paraion)
      if(allocated(eincele))deallocate(eincele)
      if(allocated(eincstp))deallocate(eincstp)
      if(allocated(eincion))deallocate(eincion)
      if(allocated(cdiamele))deallocate(cdiamele)
      if(allocated(cdiamion))deallocate(cdiamion)
      if(allocated(izion))deallocate(izion)

      return
      end subroutine DEALLOCATE_SED
!------------------------------------------------------------------------

      function sedmean(x,depev)
      implicit double precision (a-h,o-z)
      sedmean=0.0
      A9=0.0 ! mean w-value is needed to be calculated
      if(depev.eq.0) then ! electron mode
       do ie=ie1,ie1+1
        if(ie.eq.ie1) then
         Re=1.0d0-ratioe
        else
         Re=ratioe
        endif
        do ic=ic1,ic1+1
         if(ic.eq.ic1) then
          Rc=1.0d0-ratioc
         else
          Rc=ratioc
         endif
         do i=1,mparased
          if(iaustp.lt.0) then ! stopped electron
           Apara(i)=parastp(ic,abs(iaustp),ie,i,imodel)
          else
           Apara(i)=paraele(ic,abs(iaustp),ie,i,imodel)
          endif
         enddo
         wei=Re*Rc
         sedmean=sedmean+sedfunc(x,depev)*wei
         A9=A9+Apara(9)*wei
        enddo
       enddo
      else ! ion mode
       do ip=ip1,ip1+1
        if(ip.eq.ip1) then
         Rp=1.0d0-ratiop
        else
         Rp=ratiop
        endif
        do ie=ie1,ie1+1
         if(ie.eq.ie1) then
          Re=1.0d0-ratioe
         else
          Re=ratioe
         endif
         do ic=ic1,ic1+1
          if(ic.eq.ic1) then
           Rc=1.0d0-ratioc
          else
           Rc=ratioc
          endif
          do i=1,mparased
           Apara(i)=paraion(ic,ip,ie,i,imodel)
          enddo
          wei=Rp*Re*Rc
          sedmean=sedmean+sedfunc(x,depev)*wei
          A9=A9+Apara(9)*wei
         enddo
        enddo
       enddo
      endif

      Apara(9)=A9
      return

      end function sedmean


      function sedfunc(x,depev)
      implicit double precision (a-h,o-z)

      if(Apara(1).gt.0.0) then
       if(depev.eq.0.0) then ! electron mode
        tmp=min(50.0,ABS(x-apara(2))**apara(3)/(2*apara(2))) ! avoid overflow
        getfirst=Apara(1)*EXP(-tmp)  ! Poisson distribution approximated by Gaussian
       else  ! ion mode
        cst1=depev/Apara(9) ! maxmum number of event when ion pass through the center
        tmp=min(50.0,Apara(2)*(x-cst1*Apara(3)))  ! avoid overflow
        getfirst=apara(1)*x/(exp(tmp)+1)*(2.0/(cst1*Apara(3))**2)  ! core part
       endif
      else
       getfirst=0.0
      endif

      if(Apara(4).gt.0.0) then
       tmp=min(50.0,ABS(x-apara(5))**apara(6)/(2*apara(5)))
       getsecond=apara(4)*EXP(-tmp)
      else
       getsecond=0.0
      endif

      if(Apara(7).gt.0.0) then
       getthird=Apara(7)/(Apara(8)-1.0)*((Apara(8)-1.0)/Apara(8))**x
      else
       getthird=0.0
      endif

      sedfunc=getfirst+getsecond+getthird

      if(sedfunc.lt.1.0e-10) sedfunc=0.0

      return

      end function sedfunc
!------------------------------------------------------------------------
      subroutine getAparaele(CDoriginal,EE,im) ! get electron fitting parameter
!     CDoriginal: Cell diameter in um (negative for excluding Auger peak)
!     EE: Electron energy in MeV
!     im: model name (1:water)
      implicit double precision (a-h,o-z)

      dimension A(2) ! temporary local dimension

      CD=abs(CDoriginal)

      ia=1 ! with Auger peak
      if(CDoriginal.lt.0.0) ia=2 ! without Auger peak

      iaustp=ia ! Auger mode or not, negative for stopped electron
      imodel=im ! model ID

      do ic=1,mcdiamele-1
       if(Cdiamele(ic).ge.CD) exit
      enddo
      if(ic.eq.1) then
       ic1=1
       ratioc=0.0
      else
       ic1=ic-1
       ratioc=min(1.0d0,(log10(CD)-log10(Cdiamele(ic1)))
     &        /(log10(Cdiamele(ic1+1))-log10(Cdiamele(ic1))))
      endif

      do ie=2,neincele-1
       if(eincele(ie).ge.EE) exit
      end do
      ie1=ie-1
      ratioe=max(0.0,min(1.0d0,(log10(EE)-log10(eincele(ie1)))
     &      /(log10(eincele(ie1+1))-log10(eincele(ie1)))))

      end subroutine getAparaele

!------------------------------------------------------------------------
      subroutine getAparastp(CDoriginal,EE,im) ! get stopped electron fitting parameter
!     CDoriginal: Cell diameter in um (negative for excluding Auger peak)
!     EE: Electron energy in MeV
!     im: model name (1:water)
      implicit double precision (a-h,o-z)

      dimension A(2) ! temporary local dimension

      CD=abs(CDoriginal)

      ia=1 ! with Auger peak
      if(CDoriginal.lt.0.0) ia=2 ! without Auger peak

      iaustp=-ia ! Auger mode or not, negative for stopped electron
      imodel=im ! model ID

      do ic=1,mcdiamstp-1
       if(Cdiamstp(ic).ge.CD) exit
      enddo
      if(ic.eq.1) then
       ic1=1
       ratioc=0.0
      else
       ic1=ic-1
       ratioc=min(1.0d0,(log10(CD)-log10(Cdiamstp(ic1)))
     &        /(log10(Cdiamstp(ic1+1))-log10(Cdiamstp(ic1))))
      endif

      do ie=2,neincstp-1
       if(eincstp(ie).ge.EE) exit
      end do
      ie1=ie-1
      ratioe=max(0.0,min(1.0d0,(log10(EE)-log10(eincstp(ie1)))
     &      /(log10(eincstp(ie1+1))-log10(eincstp(ie1)))))

      end subroutine getAparastp


!------------------------------------------------------------------------
      subroutine getAparaion(CD,EE,iz,depev,im) ! get ion fitting parameter
!     CD: Cell diameter in um (always positive)
!     EE: Ion energy in MeV/n
!     iz: Ion charge
!     im: model name (1:water)

      implicit double precision (a-h,o-z)

      dimension A(6) ! temporary local dimension

      imodel=im ! model ID

      do ic=1,mcdiamion-1
       if(Cdiamion(ic).ge.CD) exit
      enddo
      if(ic.eq.1) then
       ic1=1
       ratioc=0.0
      else
       ic1=ic-1
       ratioc=min(1.0d0,(log10(CD)-log10(Cdiamion(ic1)))
     &        /(log10(Cdiamion(ic1+1))-log10(Cdiamion(ic1))))
      endif

      do ie=1,neincion-1
       if(eincion(ie).ge.EE) exit
      end do
      if(ie.eq.1) then
       ie1=1
       ratioe=0.0
      else
       ie1=ie-1
       ratioe=min(1.0d0,(log10(EE)-log10(eincion(ie1)))
     &      /(log10(eincion(ie1+1))-log10(eincion(ie1))))
      endif

      do ip=2,npartion-1
       if(izion(ip).ge.iz) exit
      enddo
      ip1=ip-1
      ratiop=min(1.0d0,(iz-izion(ip1))*1.0d0/(izion(ip1+1)-izion(ip1))) ! linear interpolation

      end subroutine getAparaion

      end module TSEDMOD

