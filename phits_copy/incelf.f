!************************************************************************
!*                                                                      *
!      SUBROUTINE incelf(finput,nopart,kind,eray,aray,bray,gray)
!     SUBROUTINE incelf(finput,nopart,kind,eray,aray,bray,gray,itime)
      SUBROUTINE incelf(finput,nopart,kind,eray,aray,bray,gray)
!*                                                                      *
!*       main routine of iwamoto model                                  *
!*       modified by H.Iwamoto on 04/20/2010                            *
!*       last modified by Y.Sawada on 08/27/2010                        *
!*                                                                      *
!*     input  :                                                         *
!*       finput  : input data                                           *
!*       finput(1)  : target mass                                       *
!*       finput(2)  : target charge                                     *
!*       finput(3)  : projectile energy (MeV)                           *
!*       finput(4)  : 0.0 fixed                                         *
!*       finput(5)  : 1.0 fixed                                         *
!*       finput(6)  : andit ,option of angular distribution of delta    *
!*       finput(7)  : projectile id, ityp - 1                           *
!*       nnn    : dimension size                                       *
!*                                                                      *
!*     output :                                                         *
!*       nopart  : number of out going particles                        *
!*       kind(i) : particle type of i-th particle                       *
!*       eray(i) : energy of i-th particle (MeV)                        *
!*       aray(i),bray(i),gray(i) : unit momentum vector of i-th particle*
!*                                                                      *
!*                                                                      *
!************************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

      dimension  finput(7)
      dimension  kind(nnn),eray(nnn),aray(nnn),bray(nnn),gray(nnn)
      dimension  icoales(nnn),ncoll(5,nnn)
      ntargmas = int(finput(1))
      ntargchg = int(finput(2))
      nproj    = int(finput(7))
      eein     = finput(3)

!   incelfmain(nproj,ntargmas,ntargchg,eein,nopart,kind,eray,
!     &     aray,bray,gray)
!      CALL incelfmain(nproj,ntargmas,ntargchg,eein,nopart,kind,eray,
!    &     aray,bray,gray,itime) ! for debug
       CALL incelfmain(nproj,ntargmas,ntargchg,eein,nopart,kind,
     &     eray,aray,bray,gray,eex,icoales,ncoll,idbreak,bimp)! for debug    eex,icoales add.        

      END SUBROUTINE


!***********************************************************************
!*                                                                     *
!      SUBROUTINE incelfmain(nproj,ntargmas,ntargchg,eein,nemit,kinda,
!     &     energy,alpha,beta,gamm)
!     SUBROUTINE incelfmain(nproj,ntargmas,ntargchg,eein,nemit,kinda,
!     &     energy,alpha,beta,gamm,itime)
      SUBROUTINE incelfmain(nproj,ntargmas,ntargchg,eein,nemit,
     &    kinda,energy,alpha,beta,gamm,eex,icoales,ncoll,idbreak,bimp) ! for debug eex icoales      
!*                                                                     *
!*           Copyright (C) Hiroki IWAMOTO,  2010                       *
!*                                                                     *
!*       Hiroki Iwamoto                                                *
!*       Nuclear Transmutation Technology Group,                       *
!*       Japan Atomic Energy Agency                                    *
!*                                                                     *
!*       Yusuke Sawada                                                 *
!*       School of Enegineering, Kyushu University                     *
!*                                                                     *
!*    -------------------------------------------------------------    *
!*                                                                     *
!*  [ incelfmain ]                                                     *
!*    s  incelfmain    main program of the INC calculation             *
!*    s  relcoll    to calculate the kinematics between two particles  *
!*    s  crosww     to determine collisoin channel                     *
!*    s  incgrnd    to make ground state of the target nucleus         *
!*    s  packinc    to make ground state by random packing method      *
!*    f  eliqq      to calculate liquid drop binding energy (MeV)      *
!*    f  trand      random number generator                            *
!*                                                                     *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'
 
!-----------------------------------------------------------------------

      real*8 inf(5,nnn)

      common /coln01/ iccoll
!$OMP THREADPRIVATE(/coln01/)

      dimension  itimes(8)
      dimension  mstq1(10),r(3,nnn),p(6,nnn),r0(nnn),p0(nnn),p1(nnn)
      dimension  r1(nnn),rrr(nnn),dpdt(3,nnn)
      dimension  ccr(3,nnn),ccp(6,nnn),nmascl(nnn),igroup(nnn)!wata2020/02/19
      dimension  dpcm(3,nnn),rcmclst(3,nnn),pcmclst(5,nnn),pcm0(nnn) ! 2022/10/5 yamaguchi
      dimension  relr(3,nnn) ! 2022/10/5 yamaguchi
      dimension  sp1(3,nnn),sp2(3,nnn),dccp(3,nnn),rrel(3,nnn)
      dimension  ipot(nnn),inio(nnn),ichg(nnn),inds(nnn),iavd(nnn)
      dimension  iavd2(nnn,nnn)   ! 2022/10/5 yamaguchi
      dimension  lcoll(30),nref(nnn),iclst(nnn)
      dimension  icoales(nnn),ncoll(5,nnn),ipknock(nnn)
      dimension  ipotcc(nnn),iniocc(nnn),ininclst(nnn)!watanabe2020/02/27
      dimension  numclst(nnn),nmasclst(nnn) ! 2022/10/5 yamaguchi
      dimension  gamm(nnn),alpha(nnn),beta(nnn),energy(nnn),kinda(nnn)
      dimension  pfr(3,nnn)   ! 2022/10/5 yamaguchi
      character(len=100)  logemit(nnn)

!-------------------------------------------------------------sawada8/27
      integer   inuc(nnn),inun(nnn),ihis(nnn),ibry(nnn)
      integer   massal,massba,ntmax,ncolexb,colexb
      real*8 rr2(nnn,nnn),rt00

      real(8) coulom(6)
!------------------------------------------------------watanabe2018/8/29
      integer   nmasej(nnn),idcas,nbend,gdrflg,nchgpr
!-----------------------------------------------------------------------

!*    Before include inelastic collidion,massal is constant.But,pion   *
!*    is cleated by delta or N* decay.Because of it,massal is variable.*
!*    massba :: number of nmasta + nmaspr                              *
!*    massal :: number of nmasta + nmaspr + pi                         *
      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )

      if( nproj .eq. 1 ) then     ! neutron
         iproj  = 0
	 nchgpr = 0 !watanabe 2018/8/31
         nmaspr = 1
      elseif( nproj .eq. 0 ) then ! proton
         iproj  = 1
	 nchgpr = 1
         nmaspr = 1
      elseif( nproj .eq. 2 ) then ! pion+
         iproj  = 2
	 nchgpr = 1
         nmaspr = 1
      elseif( nproj .eq. 3 ) then ! pion0
         iproj  = 3
	 nchgpr = 0
         nmaspr = 1
      elseif( nproj .eq. 4 ) then ! pion-
         iproj  = 4
	 nchgpr = -1
         nmaspr = 1
      elseif( nproj .eq. 17 ) then ! alpha
         iproj = 5
         nchgpr = 2
         nmaspr = 4
      else
         iproj  = 1
      end if
      nmasta = ntargmas
      nchgta = ntargchg
!      nmaspr = 1
      ein    = eein
      upot   = 45.0d0
      dt     = 1.0d0

      ntmax  = 100 ![default: 100]
!-----------------------------------------------------------------------
      icount = 0
!-----------------------------------------------------------------------

      ipid = 1
      irid = 1
      iaid = 1
      ixid = 1
      
      mstq1(1) = ipid           ! Momentum distribution
      mstq1(2) = irid           ! Spatial distribution
      mstq1(3) = iaid           ! Angular distribution
      mstq1(4) = ixid           ! N-N Elastic cross section
      
      gdrflg = -1 !watanabe 2018/8/30

 112  continue

      ll = 0

!-----------------------------------------------------------------------
!     final time of time evolution
!-----------------------------------------------------------------------

!watanabe2018/8/21      tfin = 70.0 * ( dble(nmasta)/208.0 )**0.16
      tfin = 80.0 * ( dble(nmasta)/208.0 )**0.16
!  cascade no shuuryo hantei wo jikan hantei kara kuukan ni yoru hantei he kaeta.
!-----------------------------------------------------------------------
!     initialization
!-----------------------------------------------------------------------

c$$$      CALL setupelf

      CALL incint(inuc,iavd,inds,inun,ipot,inio,ihis,nref,iclst,
     &                  icoales,r,p,lcoll,ncoll,iclcoll,
     &                  ntag,ichannel,rr2)
      
      do i = 1,nnn
       gamm(i)  = 0.0d0
       alpha(i) = 0.0d0
       beta(i)  = 0.0d0
       energy(i)= 0.0d0
       kinda(i) = 0
       do j = 1,2
        inf(j,i) = 0.0d0
       end do
       do j = 3,5
        inf(j,i) = 0.0d0
       end do
      end do

!-----------------------------------------------------------------------
!     make ground state
!-----------------------------------------------------------------------
      CALL incgrnd(iproj,nmasta,nchgta,nmaspr,mstq1,ein,upot,r,p,
     &     ipot,inio,ichg,bimp,rt00,bmax,tfm,ll,radm,pfm,tt0,ihis,iavd,
     &     inuc,ibry,inds,inun,massal,massba,iclst,ngtgr )
!-----------------------------------------------------------------------

      peex = 0.0d0   ! 2022/10/5 yamaguchi

         ininclst = 0
         numclst = 0
         pfr=0.0d0
      
!-----------------------------------------------------------------------

!---------------------------------------------------------------nogamine
!                1st change momentum to bend polar angle 
!------------------------------------------------------------  2012/4/27
				!arata ni ishoku sita. watanabe 8/27
      nbend = 1
!      CALL deflection( nmasta,nmaspr,ein,upot,ucpot,r,p,iclst,
!     &     inds,nbend,iproj,ichg )
          CALL deflection( nmasta,nmaspr,ein,upot,ucpot,r,p,iclst,
     &	inds,nbend,iproj,ichg,itime,ipot,0,0) ! for debug
        
      coulom(1) = cbar(dble(nchgta),dble(nmasta),1,1) ! p
      coulom(2) = cbar(dble(nchgta),dble(nmasta),1,2) ! d
      coulom(3) = cbar(dble(nchgta),dble(nmasta),1,3) ! t
      coulom(4) = cbar(dble(nchgta),dble(nmasta),2,3) ! h
      coulom(5) = cbar(dble(nchgta),dble(nmasta),2,4) ! a

!-----------------------------------------------------------------yamada
!	 		whether elastic or inelastic
!--------------------------------------------------------------2013/7/22
				!arata ni ishoku sita. 2018/08/21 watanabe
!	CALL collex( rt00,radm,nmasta,nchgta,ein,p,peex,gdrflg )
!      CALL  collex(iproj,rt00,radm,nmasta,nchgta,ein,p,peex,gdrflg)
       CALL collex( iproj,rt00,radm,nmasta,nchgta,ein,p,peex,
     &    ncolexb,colexb,ichg,nmaspr,ngtgr,q,eisq,bmax)!018/11/14katayama
     
       CALL collex_bending( nmasta,nmaspr,ein,upot,ucpot,p,iclst,
     &    iproj,colexb)     
!-----------------------------------------------------------------------

        !initialize

      iavd2 = 0
!	 idcas = 10 !wata2018/08/21
!	 igroup = 0 !wata2020/02/19
!     nmascl = 1
      nmasclst = 0   ! 2022/10/5 yamaguchi
      rcmclst = 0.0d0   ! 2022/10/5 yamaguchi
      pcmclst = 0.0d0   ! 2022/10/5 yamaguchi
      relr = 0.0d0   ! 2022/10/5 yamaguchi
!-----------------------------------------------------------------t.mori
         CALL break_up(iclst,icoales,ipot,inio,ichg,inds,iproj,
     &    nmasta,nchgta,upot,ein,p,r,pfr,idbreak,numclst,ininclst,
     &    relr,pcmclst,rcmclst,nmasclst,nmaspr,idst,idstclst,bimp)!t.mori
!	 		Direct pick up reaction( only pdx )
!--------------------------------------------------------------2016/05/26
!      if( gdrflg .ne. 1 )then

       if( ncolexb .eq. 0) then

        CALL dpickup(iclst,iproj,icoales,ipot,inio,nmaspr,nmasta
     &    ,nchgta,ein,upot,dupot,thupot,aupot,p,r,peex,ichg,inds,tfm
     &    ,radm,numclst,ininclst,nmasclst,pcmclst,rcmclst,relr)

!-----------------------------------------------------------------t.mori
!	\81@\81@\81@\81@\81@ high energy knock out reaction( only pdx )
!--------------------------------------------------------------2016/05/26
!        if( ininclst(1) .eq. 0 )then   !avoid coincidence pickup and knockout
         ipknock(1) = 0  

          CALL knock_clst(iclst,icoales,ipot,nmaspr,nmasta,nchgta,ein,
     &          upot,dupot,thupot,aupot,p,r,ichg,inds,tfm,numclst,
     &          ininclst,rcmclst,pcmclst,relr,nmasclst,ipknock,iproj)
!           CALL knock_clst(iclst,icoales,ipot,nmaspr,nmasta,nchgta,
!     &          upot,dupot,thupot,aupot,p,r,ichg,inds,tfm,numclst,
!     &          ininclst,rcmclst,pcmclst,relr,nmasclst,nknock,nknockp)

        if(ipknock(1).ne.1)then

        CALL dpid(iclst,iproj,icoales,ipot,inio,nmaspr,nmasta
     &    ,nchgta,ein,upot,dupot,thupot,aupot,p,r,peex,ichg,inds,tfm
     &    ,radm,numclst,ininclst,nmasclst,pcmclst,rcmclst,relr,idn1
     &    ,idpid)

        endif

      end if
        
!-----------------------------------------------------------------------
!     start time evolution
!-----------------------------------------------------------------------
      ucpot  = upot   + coul(nchgta,rt00)
      ducpot = dupot  + coul(nchgta,rt00)
      tucpot = thupot + coul(nchgta,rt00)
      hucpot = thupot + 2.0d0* coul(nchgta,rt00)
      aucpot = aupot  + 2.0d0*coul(nchgta,rt00)

      do i1 = 1, nmasta + nmaspr
         do j = 1,3
            dpcm(j,i1) = 0.0d0  !\8FՓˑO\8C\E3\82̃N\83\89\83X\83^\81[\82̉^\93\AE\97ʕω\BB\97\CA
         end do
      end do

      do nt = 1, ntmax

!watanabe2018/8/21         if( dble(nt)*dt .gt. tfin )  exit
         if( dble(nt)*dt .gt. tfin )  exit

!************************************************************watanabe8/21
!*    cascade particle ga sbete kakugai ni arunara exit suru.           *
!*    mada kakunai ni nokotte iru nara cascade ha mada tsudukeru.       *
!	    if( idcas .eq. 0 ) exit
!	    idcas = 0
!************************************************************************

!*************************************************************sawada8/27
!*    to calculate the decay of delta or N*                            *
!!         call pionemm(tfm,ichg,inuc,ibry,inds,inun,iavd,ihis,ipot,
!!     &         iclst,nchgta,massal,massba,upot,ucpot,dt,r,p,rt00,lcoll,
!!     &         igroup,nmascl,ipotcc,ccp)
         call pionemm(tfm,ichg,inuc,ibry,inds,inun,iavd,ihis,ipot,
     &         iclst,nchgta,massal,massba,upot,ucpot,dt,r,p,rt00,lcoll,
     &         numclst,nmasclst,pcmclst) ! 2022/10/5 yamaguchi
!***********************************************************************

!*                                                                     *
!***********************************************************************

!--------------------------------------------------------------------

         do i1 = 1, massal

!----------------------------------------------------------
            if( p(5,i1) .lt. 0.10d0 )cycle     !sawada
            if( iclst(i1) .eq. -1 )cycle     !sawada
!	    if( igroup(i1) .ne. 0 ) cycle  !watanabe20200219
!	    if( nmascl(i1) .gt. 1 ) cycle !watanabe20200227

            e1     = p(4,i1)
            em1    = p(5,i1)
            t1     = e1 - em1
!-----------------------------------------------------------
!            nchg   = ichg(i1)
            r0(i1) = sqrt( r(1,i1)**2 + r(2,i1)**2 + r(3,i1)**2 )
!-----------------------------------------------------------------------
           if( inds(i1) .eq. 1 ) then
            if( iclst(i1) .eq. 4 )then
              if( t1 .gt. aupot ) then
                  ipot(i1) = 1
              else
                  ipot(i1) = 0
              end if
            elseif( iclst(i1) .eq. 3 )then
              if( t1 .gt. thupot ) then
                  ipot(i1) = 1
              else
                  ipot(i1) = 0
              end if
            elseif( iclst(i1) .eq. 2 )then
              if( t1 .gt. thupot ) then
                  ipot(i1) = 1
              else
                  ipot(i1) = 0
              end if
            elseif( iclst(i1) .eq. 1 )then
              if( t1 .gt. dupot ) then
                  ipot(i1) = 1
              else
                  ipot(i1) = 0
              end if
            elseif( iclst(i1) .eq. 0 )then
             if( ichg(i1) .eq. 0 ) then
                if( t1 .gt. upot ) then
                   ipot(i1) = 1
                else
                   ipot(i1) = 0
                end if
             else
                if( t1 .gt. ucpot ) then
                   ipot(i1) = 1
!2012/12/13nogamine                elseif( t1 .lt. upot ) then
                else
                   ipot(i1) = 0
                end if
             end if
            elseif( iclst(i1) .eq. -1 )then
                   ipot(i1) = 0
            end if

           else
                ipot(i1) = 1              ! pion,N* and delta are moving
           end if

          if( numclst(i1).ne.0 )then

              iaa = numclst(i1)
     
	      if( nmasclst(iaa) .eq. 0 )then

              else if( pcmclst(4,iaa) - pcmclst(5,iaa)
     &            .gt. dble(nmasclst(iaa)) * upot )then
                 ipot(i1) = 1
              else
                 ipot(i1) = 0
              end if

          end if

          if( ipot(i1) .eq. 0 ) cycle!hiraoka
          if( inio(i1) .eq. 1 ) cycle
          if( numclst(i1) .ne. 0 
     &           .and. nmasclst(i1) .eq. 0.0d0 ) cycle !wata0428

!***********************************************************watanabe8/21
!         cascade particle ga mada nokotte iruka?		      *
!         nokotte itara cascade ha mada tsuduku.		      *
!---------------------------------------------------------------------*
!	  if( ipot(i1) .eq. 1 .and. inio(i1) .eq. 0 )idcas = 1
!**********************************************************************

!---------------------------------------------------------------wata0916
!               kinematics for composite particles
!-----------------------------------------------------------------------

	  if( numclst(i1) .ne. 0 )then

	     r0(i1) = dsqrt( rcmclst(1,i1)**2 
     &		+ rcmclst(2,i1)**2
     &		+ rcmclst(3,i1)**2 )

	     pcm0(i1) = dsqrt( pcmclst(1,i1)**2
     &  	+ pcmclst(2,i1)**2
     &		+ pcmclst(3,i1)**2 )

	     do j = 1,3
		rcmclst(j,i1) = rcmclst(j,i1)
     &		+ pcmclst(j,i1)
     &		/ dsqrt( pcmclst(5,i1)**2 + pcm0(i1)**2 )
     &          * dt
	     end do

             do ia = 1, nmasta+nmaspr

		if( numclst(ia) .eq. numclst(i1) )then

		   do j = 1,3
                     r(j,ia) = rcmclst(j,i1) + relr(j,ia)
		   end do

		end if

             end do
            
!-----------------------------------------------------------------------
!           relativistic kinematics
!-----------------------------------------------------------------------

	  else if( numclst(i1) .eq. 0 )then
             
            p0(i1) = sqrt( p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )
            r0(i1) = sqrt( r(1,i1)**2 + r(2,i1)**2 + r(3,i1)**2 )

             do j = 1, 3
                r(j,i1) = r(j,i1)
     &              + p(j,i1) / sqrt( em1**2 + p0(i1)**2 ) * dt
             end do

          end if
!-----------------------------------------------------------------------

            r1(i1)  = sqrt( r(1,i1)**2 + r(2,i1)**2 + r(3,i1)**2 )
            p1(i1)  = sqrt( p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )
            p(4,i1) = sqrt( p1(i1)**2 + em1**2 )

!=======================================================================
!           Surface Coalescence
!=======================================================================

            if( r0(i1) .lt. radm .and. r1(i1) .gt. radm .and.
!     &           iclst(i1) .eq. 0 .and. icoales(i1) .eq. 0 .and.
     &           numclst(i1) .eq. 0 .and. icoales(i1) .eq. 0 .and.
     &           inds(i1) .eq. 1 ) then
!               CALL coales(icoales,lcoll,i1,nchgta,massba,rt00,r,p,inds,
!     &              ipot,ichg,iclst,upot,ucpot,ein)
               CALL coales(icoales,lcoll,i1,nchgta,massba,rt00,r,p,inds,
     &           ipot,ichg,iclst,upot,ucpot,ein,numclst,pcmclst,
     &           ininclst,nmasclst,inio)
            end if

!-----------------------------------------------------------------------

!            if( r0(i1) .gt. radm * 1.30 ) then
            if( r0(i1) .gt. bmax + 7.30 ) then
               inio(i1) = 1
            end if

!-----------------------------------------------------------------------

         end do

!-----------------------------------------------------------------------
!            save momenta of nucleons before N-N collision
!-----------------------------------------------------------------------

           do ia = 1, massal
             do j = 1,3
               sp1(j,ia) = p(j,ia)
             end do
           end do

!-----------------------------------------------------------------------
!        collision term
!-----------------------------------------------------------------------

!         CALL relcoll(inuc,iavd,inds,inun,ipot,inio,ichg,nmasta,nchgta,
!     &        nmaspr,massba,massal,rt00,tfm,upot,ucpot,dt,r,p,lcoll,
!     &        ncoll,mstq1,nt,iclst,nnl,ein,ihis,iclcoll,icoales)
!         CALL relcoll(inuc,iavd,inds,inun,ipot,inio,ichg,nmasta,nchgta,
!     &        nmaspr,massba,massal,rt00,tfm,upot,ucpot,dt,r,p,lcoll,
!     &          ncoll,mstq1,nt,iclst,nnl,ein,ihis,iclcoll,icoales,
!     &          igroup,nmascl,ipotcc,ccp)
         CALL relcoll(inuc,iavd,inds,inun,ipot,inio,ichg,nmasta,nchgta,
     &        nmaspr,massba,massal,rt00,tfm,upot,ucpot,dt,r,p,lcoll,
     &          ncoll,mstq1,nt,iclst,nnl,ein,ihis,iclcoll,icoales,
     &          iavd2,ininclst)   ! 2022/10/5 yamaguchi

!------------------------------------------------------------------------
!       To calculate each cluster's momentum after N-N collision
!------------------------------------------------------------------------

	do i1 = 1, nmasta+nmaspr

           if( numclst(i1) .eq. 0 )cycle
	   if( pcmclst(5,i1) .eq. 0 )cycle

	   dpcm(1:3,i1) = 0.0d0

           do ia = 1, nmasta+nmaspr !wata0926

             if( numclst(ia) .eq. 0 )cycle

	     if( numclst(ia) .eq. numclst(i1) )then

		do j =1,3
		  sp2(j,ia) = p(j,ia)
		  dpcm(j,i1) = dpcm(j,i1) + (sp2(j,ia) - sp1(j,ia))
		end do

	     end if

           end do

	   do j = 1,3
		pcmclst(j,i1) = pcmclst(j,i1) + dpcm(j,i1)
	   end do

	   pcmclst(4,i1) = dsqrt( pcmclst(1,i1)**2
     &                          + pcmclst(2,i1)**2
     &                          + pcmclst(3,i1)**2
     &                          + pcmclst(5,i1)**2 )

	end do

        do ia = 1,nmasta+nmaspr

          if( numclst(ia) .eq. 0 ) cycle

          if( numclst(ia) .eq. numclst(i1) ) then
             
             eco = 0

             peco = sqrt(p(4,ia)**2. - p(5,ia)**2.)

             do j = 1,3
               p(j,ia) = pcmclst(j,i1)/nmasclst(i1)
               eco = eco + p(j,ia)**2.
             end do

             p(4,ia) = sqrt( eco + p(5,ia) **2.)

          end if

        end do
        end do

!*****************************************************nogamine*2012/9/05

!*************************************************************sawada8/27
!*    to calculate pion absorption                                     *

!!         CALL pionabb(ein,ichg,inuc,ibry,inds,inun,iavd,ihis,inio,
!!     &                iclst,nmasta,massal,massba,dt,r,p,lcoll)
         CALL pionabb(ein,ichg,inuc,ibry,inds,inun,iavd,ihis,inio,
     &                iclst,nmasta,massal,massba,dt,r,p,lcoll,ininclst)
*                                                                     *
!***********************************************************************
!      end do
!-----------------------------------------------------------------------
!     end time evolution of cascade stage
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!   NN Collision mo Collective Excitation mo okotte inai naraba...
!   ->  mou ikkai keisan yarinaoshi.
!-----------------------------------------------------------------------

!      if( lcoll(7) .eq. 0 .and. gdrflg .eq. 0 ) then!watanabe8/27
!         icount = icount + 1
!!         if( mod(icount,3000) .eq. 0 ) then
!!            write(*,'(27x,'' event ='',3x,i9)') icount
!!         end if
!         if( icount .lt. 20000 )then
!          goto 112
!         else
!!            if( inds(1) .eq. 4 )then
!!              lcoll(7) = lcoll(7) + 1
!!            end if
!         end if
!      else
!c$$$         nkk = nkk + 1 !not given (S.H. on 2012.10.13)
!      end if

*-----------------------------------------------------------------------
*        Final pion decay and Final analysis of clusters
*-----------------------------------------------------------------------

c         CALL fpidecayy(lcoll,upot,ucpot,tfm,ichg,inuc,ibry,inds,inun,
CS.H. on 2012.12.13
c$$$     &        iavd,ihis,massal,massba,r,p,rt00,ipot)
c     &        iavd,ihis,massal,massba,r,p,rt00,ipot,iclst)
c         CALL fpidecayy(lcoll,upot,ucpot,tfm,ichg,inuc,ibry,inds,inun,
CS.H. on 2012.12.13
c$$$     &        iavd,ihis,massal,massba,r,p,rt00,ipot)
!!         CALL fpidecayy(lcoll,upot,ucpot,tfm,ichg,inuc,ibry,inds,inun,
!!    &        iavd,ihis,massal,massba,r,p,rt00,ipot,iclst,itime,nchgta)
!!          CALL fpidecayy(lcoll,upot,ucpot,tfm,ichg,inuc,ibry,inds,inun,
!!     &        iavd,ihis,massal,massba,r,p,rt00,ipot,iclst,nchgta)
!!          CALL fpidecayy(lcoll,upot,ucpot,tfm,ichg,inuc,ibry,inds,inun,
!!     &     iavd,ihis,massal,massba,r,p,rt00,ipot,iclst,nchgta,
!!     &     igroup,nmascl,ipotcc,ccp)
          CALL fpidecayy(lcoll,upot,ucpot,tfm,ichg,inuc,ibry,inds,inun,
     &     iavd,ihis,massal,massba,r,p,rt00,ipot,iclst,nchgta,
     &     numclst,nmasclst,pcmclst) ! 2022/10/5 yamaguchi

!-----------------------------------------------------------------yamada
!                     tunnel effect for heavy target
!-------------------------------------------------------------2013/12/29

! 		CALL tunneff( ipot,upot,ucpot,p,rt00,
!     &         nchgpr,nmaspr,nchgta,nmasta,peex,iclst,ichg )
! 		CALL tunneff( ipot,upot,ucpot,p,rt00,
!     &         nchgpr,nmaspr,nchgta,nmasta,peex,iclst,ichg,ininclst )
!		CALL tunneff( ipot,upot,ucpot,p,rt00,
!    &         nchgpr,nmaspr,nchgta,nmasta,peex,iclst,ichg,itime,inds )

!------------------------------------------------------------------wata
!           put nucleons together which was apart
!--------------------------------------------------------------2016/10/27

!        CALL integration(iclst,ipot,inio,ichg,inds,iproj,ncoll,
!     &  nmasta,nchgta,nmaspr,upot,p,r,ucpot,icoales,ininclst,
!     &  ccp,ccr,nmascl,igroup)!wata
        CALL integration(iclst,ipot,inio,ichg,inds,iproj,ncoll,
     &  nmasta,nchgta,nmaspr,upot,p,r,ucpot,icoales,ininclst,
     &  numclst,pcmclst,rcmclst,nmasclst,nmascl,rt00,peex)!wata
!----------------------------------------------------------------------
!                       give energy of Q-value
!----------------------------------------------------------------------

        CALL qvalue(iclst,nmaspr,nmasta,ein,icoales,
     &         i,upot,p,tfm,nchgta,iproj,ichg,q,ipot,idpid)

!---------------------------------------------------------------nogamine
!                2nd change momentum to bend polar angel 
!------------------------------------------------------------------11/12
							!watanabe 8/27
	   nbend = 2
!           CALL deflection( nmasta,nmaspr,ein,upot,ucpot,r,p,iclst,
!     &	inds,nbend,iproj,ichg )      
       CALL deflection(nmasta,nmaspr,ein,upot,ucpot,r,p,iclst,
     &	inds,nbend,iproj,ichg,itime,ipot,ncolexb,colexb) ! for debug

!        CALL kalbach( nmasta,nchgta,nmaspr,ein,upot,ucpot,p,iclst,
!     & iproj )

!-----------------------------------------------------------------t.mori
!                               Recoil
!------------------------------------------------------------------04/20
        CALL recoil_M(iclst,nmaspr,nmasta,ein,icoales,
     &       i,upot,p,tfm,nchgta,iproj,ichg,q,ipot,idpid)       !t.mori

!!-------------------------------------------------------------wata 8/29
!----------------------------------------------------------------------

!-----------------------------------------------------------------fukuda
!                determination of reaction 
!-------------------------------------------------------------------11/30
							!2018/8/29 wata
!       CALL barrier(p,rt00,nchgta,nmasta,nmaspr,ichg,peex,inds
!     &  ,iclst,ipot,iproj,nmasej,upot,q,tfm,tt0,massal)
!      CALL barrier(p,rt00,nchgta,nmasta,nmaspr,ichg,peex,inds
!    &  ,iclst,ipot,iproj,nmasej,upot,q,tfm,tt0,massal,itime) !for debug

!------------------------------------------------------------------wata
!         calculate the excitation energy of the remnant nucleus
!----------------------------------------------------------------------
							!2018/8/29 wata
!        CALL calext(iclst,nmaspr,nmasta,p,tfm,ipot,tt0,eex,upot)

!-----------------------------------------------------------------fukuda
!                         Recoil Energy
!------------------------------------------------------------------04/26
							!2018/8/29 wata
!       CALL recoil(iclst,nmaspr,nmasta,ein,upot,p,ipot,
!     &     nchgta,iproj,ichg,eex,peex,nmasej,ncoll,icoales,q,
!     &     nreact,mstq1,massal,inds)
!       CALL recoil(iclst,nmaspr,nmasta,ein,upot,p,ipot,
!     &     nchgta,iproj,ichg,eex,peex,nmasej,ncoll,icoales,q,
!     &     nreact,mstq1,massal,inds,itime)   !for debug

!-----------------------------------------------------------------------
!                         Final Output
!-----------------------------------------------------------------------

!      CALL sm_sumo(ein,iproj,nmasta,nchgta,nmaspr,massal,mstq1,rt00,
!     &     tfm,ipot,ichg,inds,r,p,upot,lcoll,ncoll,tt0,iclst,icoales,
!     &     logemit,nemit,inf)
!       CALL sm_sumo(ein,iproj,nmasta,nchgta,massal,rt00,ipot,ichg,inds,
!     &              p,upot,iclst,nemit,inf,eex)! eex for debug fuku
       CALL sm_sumo(ein,iproj,nmasta,nchgta,massal,rt00,ipot,ichg,inds,
     &              p,upot,iclst,nemit,inf,eex,peex)! eex for debug fuku

!-----------------------------------------------------------------------

!     if( nemit .eq. 1 .and. px .eq. 0.0 ) goto 100

      do i = 1, nemit
         nprt = int(inf(1,i))
         nneu = int(inf(2,i))
         px   = dble(inf(3,i))
         py   = dble(inf(4,i))
         pz   = dble(inf(5,i))
         pabs = sqrt( px**2 + py**2 + pz**2 )
         alpha(i)  = px / pabs
         beta(i)   = py / pabs
         gamm(i)   = pz / pabs
         energy(i) = -rmass + sqrt( pabs**2 + rmass**2 )
         if( nneu .eq. 0 .and. nprt .eq. 1 ) then ! proton
            kinda(i) = 0
!            if( energy(i) .gt. ein ) goto 112
         elseif( nneu .eq. 1 .and. nprt .eq. 0 ) then !neutron
            kinda(i) = 1
!            if( energy(i) .gt. ein ) goto 112
         elseif( nneu .eq. 1 .and. nprt .eq. 1 ) then !deuteron
            kinda(i) = 10
            energy(i) = -ddmass + sqrt( pabs**2 + ddmass**2 )
!            if( energy(i) .gt. ein ) goto 112
         elseif( nneu .eq. 2 .and. nprt .eq. 1 ) then !triton
            kinda(i) = 11
            energy(i) = -tmass + sqrt( pabs**2 + tmass**2 )
!            if( energy(i) .gt. ein ) goto 112
         elseif( nneu .eq. 1 .and. nprt .eq. 2 ) then !helion
            kinda(i) = 12
            energy(i) = -hmass + sqrt( pabs**2 + hmass**2 )
!            if( energy(i) .gt. ein ) goto 112
         elseif( nneu .eq. 2 .and. nprt .eq. 2 ) then !alpha
            kinda(i) = 13
            energy(i) = -amass + sqrt( pabs**2 + amass**2 )
!            if( energy(i) .gt. ein ) goto 112
!-------------------------------------------------------------sawada8/27
!           output pion part
         elseif( nneu .eq. 0 .and. nprt .eq. 0 ) then !pion+
              kinda(i) = 2
              energy(i) = -( pmass*1000.0d0 ) 
     &                    + sqrt( pabs**2 + (pmass*1000.0)**2 )
!              if( energy(i) .gt. ein ) goto 112
         elseif( nneu .eq.-1 .and. nprt .eq. 0 ) then !pion0
              kinda(i) = 3
              energy(i) = -( pmass*1000.0d0 ) 
     &                    + sqrt( pabs**2 + (pmass*1000.0)**2 )
!              if( energy(i) .gt. ein ) goto 112
         elseif( nneu .eq. -1 .and. nprt .eq. -1 ) then !pion-
              kinda(i) = 4
              energy(i) = -( pmass*1000.0d0 ) 
     &                    + sqrt( pabs**2 + (pmass*1000.0)**2 )
!              if( energy(i) .gt. ein ) goto 112
!-----------------------------------------------------------------------
         endif

      end do

 5643 continue
      END SUBROUTINE incelfmain

!***********************************************************************
!*                                                                     *
      SUBROUTINE incint(inuc,iavd,inds,inun,ipot,inio,ihis,nref,iclst,
     &                  icoales,r,p,lcoll,ncoll,iclcoll,
     &                  ntag,ichannel,rr2)
!*                                                                     *
!*        Purpose:                                                     *
!*              to initialize the parameters for INC                   *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

      real*8   r(3,nnn),p(6,nnn),rr2(nnn,nnn)
      integer  inio(nnn),ipot(nnn),ichg(nnn),iclst(nnn),icoales(nnn)
!      integer  lcoll(30),ncoll(5,nnn),nref(nnn),iclcoll,iccoll
      integer  lcoll(30),ncoll(5,nnn),nref(nnn),iclcoll
      integer  inds(nnn),ihis(nnn),iavd(nnn),inun(nnn),inuc(nnn)
      integer  ntag,ichannel
!      common /coln01/ iccoll

!      common /cood2r/ rr2(nnn,nnn),rnij(nnn,nnn),pp2(nnn,nnn)
!-----------------------------------------------------------------------
!     set diagonal elements of two-body quantities to be zero
!-----------------------------------------------------------------------

      do i = 1, nnn
       do j = 1, nnn
         rr2(i,j)  = 0.0d0
       end do
      end do

!-----------------------------------------------------------------------
!     initialization
!-----------------------------------------------------------------------

      do i = 1, nnn
         nref(i)    = 0
         iclst(i)   = 0
         icoales(i) = 0
         ichg(i)    = 0
         inuc(i)    = 0
         iavd(i)    = 0
         inds(i)    = 0
         inun(i)    = 0
         ipot(i)    = 0
         inio(i)    = 0
         ihis(i)    = 0
      end do

      do i = 1, 30
         lcoll(i) = 0
      end do

      do j = 1, 5
         do i = 1, nnn
            ncoll(j,i) = 0
            p(j,i)     = 0.0d0
         end do
      end do
      do j = 1, 3
         do i = 1, nnn
            r(j,i) = 0.0d0
         end do
      end do

!      iccoll  = 0
      iclcoll = 0
      ntag    = 0
      ichannel= 0

      END SUBROUTINE incint

!***********************************************************************
!*                                                                     *
!*                         Collision Term                              *
!*                                                                     *
!*                                                                     *
!*                    Last revised : Feb/15/2008                       *
!*                                                                     *
!*                                                                     *
!*  List of subprograms in rough order of relevance with main purpose  *
!*     ( s = SUBROUTINE, f = FUNCTION, b = block data, e = entry )     *
!*                                                                     *
!*                                                                     *
!*  s  relcoll     to calculate the kinematics between two particles   *
!*  s  crosww      to determine collisoin channel                      *
!*                                                                     *
!***********************************************************************
!*                                                                     *
!***********************************************************************
!*                                                                     *
!      SUBROUTINE relcoll(inuc,iavd,inds,inun,ipot,inio,ichg,nmasta,
!     &   nchgta,nmaspr,massba,massal,rt00,tfm,upot,ucpot,dt,r,p,lcoll,
!     &   ncoll,mstq1,nt,iclst,nnl,ein,ihis,iclcoll,icoales)
!      SUBROUTINE relcoll(inuc,iavd,inds,inun,ipot,inio,ichg,nmasta,
!     &   nchgta,nmaspr,massba,massal,rt00,tfm,upot,ucpot,dt,r,p,lcoll,
!     &     ncoll,mstq1,nt,iclst,nnl,ein,ihis,iclcoll,icoales,
!     &     igroup,nmascl,ipotcc,ccp)
      SUBROUTINE relcoll(inuc,iavd,inds,inun,ipot,inio,ichg,nmasta,
     &   nchgta,nmaspr,massba,massal,rt00,tfm,upot,ucpot,dt,r,p,lcoll,
     &     ncoll,mstq1,nt,iclst,nnl,ein,ihis,iclcoll,icoales,
     &     iavd2,ininclst)   ! 2022/10/5 yamaguchi
!*                                                                     *
!*                                                                     *
!*        Purpose:                                                     *
!*              to calculate the kinematics in a collision             *
!*              between two particles                                  *
!*                                                                     *
!*                                                                     *
!*        Variables:                                                   *
!*                                                                     *
!*           [in]                                                      *
!*              inds        : nucleon/pion                             *
!*              ipot        : above/under potential depth              *
!*              inio        : inside/outside the nucleus               *
!*              ichg        : isospin of nucleon                       *
!*              nmasta      : target mass                              *
!*              nchgta      : target charge                            *
!*              rt00        : nuclear radius (fm)                      *
!*              tfm         : Fermi energy (MeV)                       *
!*              upot        : potential depth (MeV)                    *
!*              ucpot       : potential depth + Coulomb barrier (MeV)  *
!*              dt          : time step size (fm/c)                    *
!*              r           : positions of nucleus                     *
!*              p           : momenta of nucleus                       *
!*              lcoll       : counting variables about collision       *
!*              mstq1       : input imformation                        *
!*                                                                     *
!*           [out]                                                     *
!*              inds        : nucleon/pion                             *
!*              ichannel    : channel information back                 *
!*              ipot        : above/under potential depth              *
!*              p           : momenta of nucleus                       *
!*              lcoll       : counting variables about collision       *
!*                                                                     *
!*                                                                     *
!*        Comments : counting variables                                *
!*                                                                     *
!*              lcoll( 1) : all collisions                             *
!*              lcoll( 2) : Pauli-blocked collisions                   *
!*              lcoll( 3) : energetically forbidden collision          *
!*              lcoll( 4) : elastic collision                          *
!*                                                                     *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      integer  inio(nnn),mstq1(10),nmasta,nmaspr,nchgta,nt
      real*8   upot,r(3,nnn),p(6,nnn),dt,tfm,ucpot,rt00,ein
      integer  inds(nnn),ipot(nnn),ichg(nnn),lcoll(30),iclst(nnn)
      integer  ncoll(5,nnn),nnl,ihis(nnn),iavd(nnn),inun(nnn),inuc(nnn)
      integer  icoales(nnn)
      integer  iavd2(nnn,nnn),ininclst(nnn)

      dimension  beta(3),pcm(3)

      common /coln01/ iccoll
!$OMP THREADPRIVATE(/coln01/)
!      common /cood2r/ rr2(nnn,nnn),rnij(nnn,nnn),pp2(nnn,nnn)
      real*8 rr2(nnn,nnn)
!*****************************************************:sawada
      real*8 uuu,accele1,accelem1,accele2,accelem2
      integer iclcoll
!*******************************************************
!      integer, intent(inout) :: igroup(nnn),nmascl(nnn),ipotcc(nnn)
!      real*8, intent(inout) :: ccp(6,nnn)      

      parameter ( sig0 =  55.0, bcmax0 = 1.323142 )
      parameter ( sig1 = 200.0, bcmax1 = 2.523    )
      parameter ( deltar = 4.0 )
!-------------------------------------------------------------sawada8/27
      integer massba,massal
!-----------------------------------------------------------------------
      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )
                  icycle = 0

      iang     = mstq1(3)
      icrs     = mstq1(4)
      ducpot = dupot  + coul(nchgta,rt00)
      tucpot = thupot + coul(nchgta,rt00)
      hucpot = thupot + 2.0d0*coul(nchgta,rt00)
      aucpot = aupot  + 2.0d0*coul(nchgta,rt00)

!      CALL d_ipot(p,rt00,upot,nchgta,massal,ichg,inds,iclst,ipot)
!      CALL d_ipot(p,rt00,upot,nchgta,massal,ichg,inds,iclst,ipot,
!     &     igroup,nmascl,ipotcc,ccp)

      do i1 = 1, massba-1 !do loop for nuleon

!***********************************************************sawada
         if( p(5,i1) .lt. 0.10d0 )cycle
         if( iclst(i1) .eq. -1 )cycle     !sawada
!         if( inds(i1) .eq. 4 )cycle     !sawada
!****************************************************************

!         if( inio(i1) .eq. 1 ) cycle

         do i2 = i1 + 1, massba

!***********************************************************sawada
               if( p(5,i2) .lt. 0.10d0 )cycle
               if( iclst(i2) .eq. -1 )cycle     !sawada
!****************************************************************

!-----------------------------------------------------------------------
!           avoid first collisions within the same nucleus
!-----------------------------------------------------------------------

!            if( iavd(i1) * iavd(i2) .eq. 1 ) cycle

!-----------------------------------------------------------------------
            if( ipot(i1) .eq. 1 .and. ipot(i2) .eq. 1 ) cycle
            if( ipot(i1) .eq. 0 .and. ipot(i2) .eq. 0 ) cycle
!           if( ininclst(i1) .ne. 0 .or. ininclst(i2) .ne. 0 ) cycle
            if( iclst(i1) .ne. 0 .and. iclst(i2) .ne. 0 ) cycle
            if( iclst(i1) .ne. 0 .and. inds(i2) .ne. 1 ) cycle
            if( i1 .le. nmaspr .and. i2 .le. nmaspr ) cycle
            if( inds(i1) .ne. 1 .and. iclst(i2) .ne. 0 ) cycle
            if( iavd2(i1,i2) .eq. 1 )cycle
!            if( il .le. nmaspr .and. iclst(i1) .eq. 0 ) cycle
!            if( i2 .le. nmaspr .and. iclst(i2) .eq. 0 ) cycle

!-----------------------------------------------------------------------
!           the following prescription (deltar) is not covariant,
!           but useful to save cpu time especially below 200 MeV/u
!-----------------------------------------------------------------------

!-------------------------------------------------------------sawada8/27
            if( iclst(i1) .eq. 0 .or. iclst(i2) .eq. 0 ) then
             CALL rr(i1,i2,r,p,rr2)
!             rr2(i1,i2) = deltar**2. - 1.
             if( ein .le. 5000.0 .and. rr2(i1,i2) .gt. deltar**2 ) cycle
            end if
!-----------------------------------------------------------------------
!           now particles are close enough to each other
!-----------------------------------------------------------------------

            px1  = p(1,i1)
            py1  = p(2,i1)
            pz1  = p(3,i1)
            e1   = p(4,i1)
            em1  = p(5,i1)
            iz1  = ichg(i1)
            id1  = inds(i1)
!            inc1 = inuc(i1)
            ipt1 = ipot(i1)
            icl1 = iclst(i1)

            px2  = p(1,i2)
            py2  = p(2,i2)
            pz2  = p(3,i2)
            e2   = p(4,i2)
            em2  = p(5,i2)
            iz2  = ichg(i2)
            id2  = inds(i2)
!            inc2 = inuc(i2)
            ipt2 = ipot(i2)
            icl2 = iclst(i2)
!-------------------------------------------------------------sawada8/27
!     save initial values and use them when pauli blocking happened.
            spx1   = p(1,i1)
            spy1   = p(2,i1)
            spz1   = p(3,i1)
            se1    = p(4,i1)
            sem1   = p(5,i1)
            isz1   = ichg(i1)
            isnd1  = inds(i1)
            isnc1  = inuc(i1)
            ispt1  = ipot(i1)
            iscl1  = iclst(i1)
            spx2   = p(1,i2)
            spy2   = p(2,i2)
            spz2   = p(3,i2)
            se2    = p(4,i2)
            sem2   = p(5,i2)
            isz2   = ichg(i2)
            isnd2  = inds(i2)
            isnc2  = inuc(i2)
            ispt2  = ipot(i2)
            iscl2  = iclst(i2)

!-----------------------------------------------------------------------
!              iclcoll = 0
              uuu = 0.0d0
!******************************************************sawada7/7
!     to accelerate nucleon by 45*xMeV when it collides cluster
           if( inds(i1) .eq. 1 .and. inds(i2) .eq. 1  )then
!            if( iclst(i1) .eq. 0 .and. iclst(i2) .eq. 0 )then
!
!              iclcoll = 0
!
!            elseif( iclst(i1) .ne. 0 .and. iclst(i2) .eq. 0 )then
!
!              iclcoll = 1
!
!             if( iclst(i1) .eq. 1 )then
!              uuu = 45.0d0
!             elseif( iclst(i1) .eq.2 .or. iclst(i1) .eq.3 )then
!              uuu = 90.0d0
!             elseif( iclst(i1) .eq. 4 )then
!              uuu = 135.0d0
!             end if

!              uuu = 0.0d0
!              iclcoll = 0

             if( em1 .lt. em2+1.0d-3 .and. em1 .gt. em2-1.0d-3) then

              px1  = p(1,i1)
              py1  = p(2,i1)
              pz1  = p(3,i1)
              e1   = p(4,i1)
              em1  = p(5,i1)
              px2  = p(1,i2)
              py2  = p(2,i2)
              pz2  = p(3,i2)
              e2   = p(4,i2)
              em2  = p(5,i2)

             elseif( em1 .gt. em2 )then

              uuu = (em1-em2)/rmass*upot
              
              unitx = px2/dsqrt( px2**2 + py2**2 + pz2**2 )
              unity = py2/dsqrt( px2**2 + py2**2 + pz2**2 )
              unitz = pz2/dsqrt( px2**2 + py2**2 + pz2**2 )

              e2 = e2+uuu
              et2 = e2-em2
              pabs = dsqrt( et2*(et2 + 2.0d0*em2) )

              px2 = pabs*unitx
              py2 = pabs*unity
              pz2 = pabs*unitz

!            elseif( iclst(i1) .eq. 0 .and. iclst(i2) .ne. 0 )then
!
!              iclcoll = 2
!
!             if( iclst(i2) .eq. 1 )then
!              uuu = 45.0d0
!             elseif( iclst(i2) .eq. 2 .or. iclst(i2) .eq. 3 )then
!              uuu = 90.0d0
!             elseif( iclst(i2) .eq. 4 )then
!              uuu = 135.0d0
!             end if

             elseif( em1 .lt. em2 )then

              uuu = (em2-em1)/rmass*upot
              
              unitx = px1/dsqrt( px1**2 + py1**2 + pz1**2 )
              unity = py1/dsqrt( px1**2 + py1**2 + pz1**2 )
              unitz = pz1/dsqrt( px1**2 + py1**2 + pz1**2 )

              e1 = e1+uuu
              et1 = e1-em1
              pabs = dsqrt( et1*(et1 + 2.0d0*em1) )

              px1 = pabs*unitx
              py1 = pabs*unity
              pz1 = pabs*unitz

            end if
           end if

             accele1  = e1
             accelem1 = em1
             accele2  = e2
             accelem2 = em2

            s    = (  e1 +  e2 )**2 - ( px1 + px2 )**2
     &           - ( py1 + py2 )**2
     &           - ( pz1 + pz2 )**2

            srt  = sqrt(s)

!-----------------------------------------------------------------------
!     low energy cutoff and max cross section
!     ( 1896.6  =  em1 + em2 + 20.0 MeV ; nucleon case )
!-----------------------------------------------------------------------

           cutoff = 1896.6d0 !not given (S.H. on 2012.10.13)

           if( iclst(i1) .eq. 0 .and. iclst(i2) .eq. 0 )then
            if( em1 .lt. 940.0 .and. em2 .lt. 940.0 ) then
               cutoff =  em1 + em2 + 20.0
               bcmax  = bcmax0
               sig    = sig0
            else
               cutoff =  em1 + em2
               bcmax  = bcmax1
               sig    = sig1
            end if
           end if

            if ( srt .lt. cutoff ) cycle

!-----------------------------------------------------------------------
!           are their impact parameter small enough ?
!-----------------------------------------------------------------------

            dx   = r(1,i1) - r(1,i2)
            dy   = r(2,i1) - r(2,i2)
            dz   = r(3,i1) - r(3,i2)
            rsq  = dx**2 + dy**2 + dz**2

            p12  = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
            p1dr = px1 * dx + py1 * dy + pz1 * dz
            p2dr = px2 * dx + py2 * dy + pz2 * dz
            a12  = 1.0 - ( em1 * em2 / p12 ) ** 2
            b12  = p1dr / em1 - p2dr * em1 / p12
            c12  = rsq + ( p1dr / em1 )**2
            brel = sqrt( abs(c12 - b12**2/a12) )

!-----------------------------------------------------------------------
!           average time-shift of the collision in the fixed frame
!           will particles get closest point in this time interval?
!-----------------------------------------------------------------------

            b21    =   - p2dr / em2 + p1dr * em2 / p12
            t1     = (   p1dr / em1 - b12 / a12 ) * e1 / em1
            t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2

            if( abs( t1 + t2 ) .gt. dt ) cycle

!=======================================================================
!           start knockout process
!           create a cluster
!=======================================================================
!-----------------------------------------------------------------------
!           Lorentz-transformation in i1-i2-c.m. system
!-----------------------------------------------------------------------

               etot12  = e1 + e2

            beta(1) = ( px1 + px2 ) / etot12
            beta(2) = ( py1 + py2 ) / etot12
            beta(3) = ( pz1 + pz2 ) / etot12
            betasq  = beta(1)**2 + beta(2)**2 + beta(3)**2
            gamma   = 1.0d0 / sqrt( 1.0d0 - betasq )

!-----------------------------------------------------------------------
!     transformation of momenta
!-----------------------------------------------------------------------

            p1beta = px1 * beta(1) + py1 * beta(2) + pz1 * beta(3)
            transf = gamma
     &           * ( gamma * p1beta / ( gamma + 1.0 ) - e1 )
            pcm(1) = px1 + beta(1) * transf
            pcm(2) = py1 + beta(2) * transf
            pcm(3) = pz1 + beta(3) * transf
            pcm2   = pcm(1)**2 + pcm(2)**2 + pcm(3)**2
            prcm   = sqrt( pcm2 )

            if( prcm .le. 0.00001d0 ) cycle

!-----------------------------------------------------------------------
!     calculate cross section and collide two particles
!-----------------------------------------------------------------------

!            CALL crosww(i1,i2,inds,ipot,ichg,cutoff,pcm,prcm,srt,
!     &           beta,gamma,accele1,accelem1,accele2,accelem2,iclcoll,
!     &           r,p,upot,ucpot,rt00,tfm,ichannel,sig,iang,
!     &           icrs,mstq1,nmasta,nchgta,iclst,ia2,ia3,ia4,knoc,ein,
!     &           inuc,nt)
            CALL crosww(i1,i2,inds,ipot,ichg,cutoff,pcm,prcm,srt,
     &           beta,gamma,accele1,accelem1,accele2,accelem2,iclcoll,
     &           r,p,upot,ucpot,rt00,tfm,ichannel,sig,iang,
     &           icrs,mstq1,nmasta,nchgta,iclst,ia2,ia3,ia4,knoc,ein,
     &           inuc,nt,ininclst)   ! 2022/10/11 yamaguchi

!-----------------------------------------------------------------------
!     ichannel : channel of this collision
!-----------------------------------------------------------------------
!
!     =  0 ; nothing has happened
!     =  1 ; elastic NN collision
!     =  2 ; N + N -> N + D
!     =  3 ; N + D -> N + D
!     =  4 ; N + N -> N + R
!     =  5 ; N + R -> N + N
!     =  6 ; N + N -> D + D
!     =  7 ; D + D -> N + N
!     =  8 ; N + D -> D + D
!     =  9 ; D + D -> N + D
!     = 10 ; N + R -> D + R
!     = 11 ; D + R -> N + R
!     = 12 ; N + D -> R + D
!     = 13 ; R + D -> N + D
!     = 14 ; N + R -> R + R
!     = 15 ; R + R -> N + R
!     = 99 ; energetically forbidden collision
!
!-----------------------------------------------------------------------

            if( ichannel .eq. 0 ) cycle

            lcoll(1) = lcoll(1) + 1
            ncoll(1,i1) = ncoll(1,i1) + 1
            ncoll(1,i2) = ncoll(1,i2) + 1
            ncoll(1,i1) = max(ncoll(1,i1),ncoll(1,i2))
            ncoll(1,i2) = ncoll(1,i1)

!-----------------------------------------------------------------------
!           a collision has taken place
!-----------------------------------------------------------------------

            ntag  = 0

            if( ichannel .eq. 99 ) then

               ntag  = 1
               lcoll(6) = lcoll(6) + 1
               ncoll(3,i1) = ncoll(3,i1) + 1
               ncoll(3,i2) = ncoll(3,i2) + 1
               ncoll(3,i1) = max(ncoll(3,i1),ncoll(3,i2))
               ncoll(3,i2) = ncoll(3,i1)

               ichannel = 0

            end if

!-----------------------------------------------------------------------
!           Pauli blocked case, reset variables
!-----------------------------------------------------------------------

            if( ntag .eq. 1 ) then

               lcoll(5)    = lcoll(5) + 1
               ncoll(2,i1) = ncoll(2,i1) + 1
               ncoll(2,i2) = ncoll(2,i2) + 1
               ncoll(2,i1) = max(ncoll(2,i1),ncoll(2,i2))
               ncoll(2,i2) = ncoll(2,i1)

               p(1,i1)   = spx1
               p(2,i1)   = spy1
               p(3,i1)   = spz1
               p(4,i1)   = se1
               p(5,i1)   = sem1
               ichg(i1)  = isz1
               inds(i1)  = isnd1
               inuc(i1)  = isnc1
               ipot(i1)  = ispt1
               iclst(i1) = iscl1

               p(1,i2)   = spx2
               p(2,i2)   = spy2
               p(3,i2)   = spz2
               p(4,i2)   = se2
               p(5,i2)   = sem2
               ichg(i2)  = isz2
               inds(i2)  = isnd2
               inuc(i2)  = isnc2
               ipot(i2)  = ispt2
               iclst(i2) = iscl2

!-----------------------------------------------------------------------
!           collision is really happened
!-----------------------------------------------------------------------

            else

               ihis1 = 1
               ihis2 = 1

               if( ihis(i1) .lt. 0 ) ihis1 = -1
               if( ihis(i2) .lt. 0 ) ihis2 = -1

               if( ichannel .gt. 1 ) then
                   ihis1 = -1
                   ihis2 = -1
               end if

               ihis0 = abs(ihis(i1)) + abs(ihis(i2)) + 1
               ihis(i1) = ihis2*ihis0
               ihis(i2) = ihis2*ihis0

               inun(i1) = i2
               inun(i2) = i1

               iavd(i1) = iavd(i1) + iavd(i1) / abs(iavd(i1))
               iavd(i2) = iavd(i2) + iavd(i2) / abs(iavd(i2))

!-----------------------------------------------------------------------
!              counting of collision
!-----------------------------------------------------------------------

c$$$               iccoll = iccoll + 1 !not given (S.H. on 2012.10.13)

               if( ichannel .eq. 1 ) then
                  lcoll(7) = lcoll(7) + 1
                  ncoll(4,i1) = ncoll(4,i1) + 1
                  ncoll(4,i2) = ncoll(4,i2) + 1
                  ncoll(4,i1) = max(ncoll(4,i1),ncoll(4,i2))
                  ncoll(4,i2) = ncoll(4,i1)
                  nnl = nnl + 1
                  iavd2(i1,i2) = 1

                 if( iclst(i1) + iclst(i2) .eq. 0 )then
                  if( ( em1 .lt. 940.0 .and. em2 .lt. 940.0 ) .and.
     &                ichg(i1) + ichg(i2) .eq. 1 ) then
                     lcoll(3) = lcoll(3) + 1
                  end if

                  if( ( em1 .gt. 940.0 .and. em2 .lt. 940.0 ) .and.
     &                ( em1 .lt. 940.0 .and. em2 .gt. 940.0 ) .and.
     &                ichg(i1) + ichg(i2) .eq. 1 ) then
                     lcoll(2) = lcoll(2) + 1
                  end if

                  if( em1 .gt. 940.0 .and. em2 .gt. 940.0 )
!-----------------------------------------------------------------------
     &               lcoll(4) = lcoll(4) + 1
                  end if
                 end if

               if( ichannel .eq.  2 ) lcoll(8)  = lcoll(8)  + 1
               if( ichannel .eq.  3 ) lcoll(9)  = lcoll(9)  + 1
               if( ichannel .eq.  4 ) lcoll(10) = lcoll(10) + 1
               if( ichannel .eq.  5 ) lcoll(11) = lcoll(11) + 1
               if( ichannel .eq.  6 ) lcoll(12) = lcoll(12) + 1
               if( ichannel .eq.  7 ) lcoll(13) = lcoll(13) + 1
               if( ichannel .eq.  8 ) lcoll(14) = lcoll(14) + 1
               if( ichannel .eq.  9 ) lcoll(15) = lcoll(15) + 1
               if( ichannel .eq. 10 ) lcoll(16) = lcoll(16) + 1
               if( ichannel .eq. 11 ) lcoll(17) = lcoll(17) + 1
               if( ichannel .eq. 12 ) lcoll(18) = lcoll(18) + 1
               if( ichannel .eq. 13 ) lcoll(19) = lcoll(19) + 1
               if( ichannel .eq. 14 ) lcoll(20) = lcoll(20) + 1
               if( ichannel .eq. 15 ) lcoll(21) = lcoll(21) + 1
            end if

!-----------------------------------------------------------------sawada
               if( ichannel .ge. 2 .and. ichannel .le. 15 ) then
                  lcoll(7) = lcoll(7) + 1
               end if
!-----------------------------------------------------------------------
!***********************************************************************
!*                                                                     *
!*    start nucleon knock_out 

!          if( ichannel .eq. 1 )then
!           if( inds(i1) .eq. 1 .and. inds(i2) .eq.1 )then
!            if( iclst(i2) .eq. 0 )then
!             CALL knockout(i1,i2,ichannel,inds,ipot,ichg,r,p,rt00,upot,
!     &                    nchgta,massba,iclst,icoales)
!            end if
!           end if
!          end if
!*                                                                     *
!***********************************************************************

         end do
      end do

      END SUBROUTINE relcoll

!***********************************************************************
!*                                                                     *
!      SUBROUTINE crosww(i1,i2,inds,ipot,ichg,cutoff,pcm,prcm,srt,
!     &           beta,gamma,accele1,accelem1,accele2,accelem2,iclcoll,
!     &      r,p,upot,ucpot,rt00,tfm,ichannel,sig,
!     &     iang,icrs,mstq1,nmasta,nchgta,
!     &     iclst,ia2,ia3,ia4,knoc,ein,inuc,nt)
      SUBROUTINE crosww(i1,i2,inds,ipot,ichg,cutoff,pcm,prcm,srt,
     &           beta,gamma,accele1,accelem1,accele2,accelem2,iclcoll,
     &      r,p,upot,ucpot,rt00,tfm,ichannel,sig,
     &     iang,icrs,mstq1,nmasta,nchgta,
     &     iclst,ia2,ia3,ia4,knoc,ein,inuc,nt,ininclst)   ! 2022/10/11 yamaguchi
!*                                                                     *
!*                                                                     *
!*        Purpose:                                                     *
!*              to determine collisoin channel                         *
!*                                                                     *
!*        Variables:                                                   *
!*           [in]                                                      *
!*              i1,i2       : identificator of particle 1 and 2        *
!*              inds        : nucleon/pion                             *
!*              ipot        : above/under potential depth              *
!*              ichg        : isospin of nucleon                       *
!*              cutoff      : cutoff energy                            *
!*                            em1 + em2 + 0.02 GeV for N + N           *
!*                            em1 + em2            for the others      *
!*              pcm(3)      : momentum coordinates of one particle     *
!*                            in cm frame                              *
!*              prcm        : sqrt(pcm)                                *
!*              srt         : sqrt of s                                *
!*              beta(3)     : beta of cm frame to reference frame      *
!*              gamma       : gamma of above beta(3)                   *
!*              ucpot       : potential depth + Coulomb barrier (MeV)  *
!*              r           : positions of nucleus                     *
!*              p           : momenta of nucleus                       *
!*              tfm         : Fermi energy (MeV)                       *
!*              sig         : max cross section at cutoff energy       *
!*              iang        : angular distribution, old or new         *
!*              icrs        : elastic cross section, Cugnon or Niita   *
!*              mstq1       : input information flag                   *
!*              nmasta      : target mass                              *
!*              nchgta      : target charge                            *
!*              ein         : incident energy in MeV                   *
!*                                                                     *
!*           [out]                                                     *
!*              p           : momenta of nucleus                       *
!*              pcm(3)      : momentum coordinates of one particle     *
!*                            in cm frame                              *
!*              inds        : nucleon/pion                             *
!*              ipot        : above/under potential depth              *
!*              ichannel    : channel information back                 *
!*                                                                     *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      integer  i1,i2,ichg(nnn),iang,icrs,mstq1(10),nchgta,nmasta
      integer  iclst(nnn),ia2,ia3,ia4,knoc
      integer  inds(nnn),ipot(nnn),ichannel
      real*8   cutoff,prcm,srt,beta(3),gamma,r(3,nnn),upot,ucpot,tfm,sig
      real*8   p(6,nnn),pcm(3),ein,rt00

      integer inuc(nnn) !sawada8/27
!*********************:sawada
      real*8 tfmc(nnn),accele1,accelem1,accele2,accelem2,uu,u,tfm1,tfm2
      real*8 dtfm,thtfm,atfm,dxdst,thxdst,axdst
      real*8 dcross,thcross,across
      integer iclcoll,nt
      integer ininclst(nnn)  ! 2022/10/11 yamaguchi
!*****************************

      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )

!-----------------------------------------------------------------------
      ducpot = dupot  + coul(nchgta,rt00)
      tucpot = thupot + coul(nchgta,rt00)
      hucpot = thupot + 2.0d0*coul(nchgta,rt00)
      aucpot = aupot  + 2.0d0*coul(nchgta,rt00)
      u    = 0.0d0
      uu   = 0.0d0
      tfm1 = tfm
      tfm2 = tfm
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     generate seeds
!-----------------------------------------------------------------------

      CALL genseed(iseed)


!-----------------------------------------------------------------------
!     choice of channels
!-----------------------------------------------------------------------
!     (1) ratio of NN -> DD to NN-> NR
         fcdd = 0.0d0
         fnnr = 1.5d0*(1.0d0-fcdd)
         fndd = 1.5d0*fcdd
!     (2) ND -> DD (2-2), (4-2)
         inddd = 1
!     (3) ND -> RD (2-3), (5-2)
         indrd = 1
!     (4) NR -> DR  (3-2), (5-1)
         inrdr = 1
!     (5) NR -> RR (3-3), (6-1)
         inrrr = 1
!-----------------------------------------------------------------------
      e1   = accele1
      em1  = accelem1
      iz1  = ichg(i1)
      id1  = inds(i1)

      e2   = accele2
      em2  = accelem2
      iz2  = ichg(i2)
      id2  = inds(i2)

!-----------------------------------------------------------------------

      rm2  = rmass**2
      pr   = prcm
      c2   = pcm(3) / pr

      csrt = srt - cutoff

      pri = prcm
      prf = sqrt( 0.25*(srt)**2 - rm2 )

      asrt = srt - em1 - em2
      pra  = prcm

!-------------------------------------------------------------sawada8/27
!           MeV->GeV
      srt  = srt / 1000.0
      asrt = asrt/ 1000.0
      pra  = pra / 1000.0
      prf  = prf / 1000.0
      pri  = pri / 1000.0
!      csrt = csrt / 1000.0      ! (MeV) -> (GeV)
!-----------------------------------------------------------------------

!-------------------------------------------------------------sawada8/27
!      plsq = (srt)**4 / ( 4.0 * (rmass/1000.0)**2 )
!     &     - (srt)**2
!-----------------------------------------------------------------------
!      pl    = 0.0d0
!      pcc   = 0.0d0

!      if( plsq .gt. 0.0 ) pl = sqrt( plsq )

!      if( iclcoll .ne. 0 )then
!       emm1 = em1/1000.0d0
!       emm2 = em2/1000.0d0
!
!       pcsq =  (srt - emm1 - emm2) *
!     &        (srt - emm1 + emm2) *
!     &        (srt + emm1 - emm2) *
!     &        (srt + emm1 + emm2)
!
!       if( pcsq .gt. 0.0 ) pcc = sqrt( pcsq )
!
!        pc = ( 1.0d0/( 2.0d0*srt ) ) * pcc
!        pl = pc*srt/emm2
!      end if

      em1pc = em1*1.0d-3
      em2pc = em2*1.0d-3
       pc = 1.0d0/( 2.0d0*srt ) * 
     &      sqrt(  (srt - em1pc - em2pc)  *
     &             (srt - em1pc + em2pc)  *
     &             (srt + em1pc - em2pc)  *
     &             (srt + em1pc + em2pc)  )

       pl = pc*srt*1.0d3/em2          !\8E\C0\8C\B1\8E\BA\8Cn\89^\93\AE\97\CA
      
      el   = -rmass + sqrt( (pl*1000.0)**2 + rmass**2 ) ! (MeV)

      csrt = csrt / 1000.0      ! (MeV) -> (GeV)      
!-----------------------------------------------------------------------
!     random number
!-----------------------------------------------------------------------

      x1 = rn(0)
      x2 = rn(0)

!-----------------------------------------------------------------------
!     Elastic cross section (icrs) :
!       J.Cugnon et al. Nucl. Phys. A620, 475 (1997)
!-----------------------------------------------------------------------

      if( iz1 .eq. iz2 ) then
         if( pl .lt. 0.8d0 ) then
            sigel = 23.5d0 + 1000.0d0 * ( pl - 0.7d0 )**4
         elseif( pl .lt. 2.0 ) then
            sigel = 1250.0d0 / ( pl + 50.0d0 )
     &           - 4.0d0 * ( pl - 1.3d0 ) **2
         else
            sigel = 77.0d0 / ( pl + 1.5d0 )
         end if
      else
         if( pl .lt. 0.8d0 ) then
            sigel = 33.0d0 + 196.0d0 * abs( pl - 0.95d0 )**2.5
         elseif( pl .lt. 2.0d0 ) then
            sigel = 31.0d0 / sqrt( pl )
         else
            sigel = 77.0d0 / ( pl + 1.5d0 )
         end if
      end if

*-----------------------------------------------------------------------
*           elastic cross section;       modified Cugnon
*-----------------------------------------------------------------------
!            if( iz1 .eq. iz2 ) then
!               if( csrt .lt. 0.4286 ) then
!                  sigel = 35.0 / ( 1. + csrt * 100.0 )  +  20.0
!               else
!                  sigel = ( - atan( ( csrt - 0.4286 ) * 1.5 - 0.8 )
!     &                  *   2. / pi + 1.0 ) * 9.65 + 7.0
!               end if
!            else
!               if( csrt .lt. 0.4286 ) then
!                  sigel = 28.0 / ( 1. + csrt * 100.0 )  +  27.0
!               else
!                  sigel = ( - atan( ( csrt - 0.4286 ) * 1.5 - 0.8 )
!     &                  *   2. / pi + 1.0 ) * 12.34 + 10.0
!               end if
!            end if
*-----------------------------------------------------------------------

      sigtt = sigel

!-----------------------------------------------------------------------
!   When minimum distance "rdst" is smallter than "xdst",
!   the pair nucleons are allowed to scatter.
!-----------------------------------------------------------------------

      cross = sigel / 10.0      ! (mb) -> (fm^2)
      xdst = dsqrt( cross / pi )
!-----------------------------------------------------------8/17 sawada
!      if( iclst(i1) .eq. 1 .or. iclst(i2) .eq. 1 ) then
!       xdst = 2.0d0*xdst
!      elseif( iclst(i1) .eq. 2 .or. iclst(i2) .eq. 2 .or.
!     &           iclst(i1) .eq. 3 .or. iclst(i2) .eq. 3 ) then
!        xdst = 3.0d0*xdst
!      elseif( iclst(i1) .eq. 4 .or. iclst(i2) .eq. 4 ) then
!        xdst = 3.5d0*xdst
!      end if
!----------------------------------------------------------------------
      
!-----------------------------------------------------------8/17 sawada
!      dcross  = 2.00d0*cross
!      thcross = 3.00d0*cross
!      across  = 3.50d0*cross
!!      if( iclcoll .ne. 0 )then
!       if( iclst(i1) .eq. 1 .or. iclst(i2) .eq. 1 ) then
!        cross = dcross
!       elseif( iclst(i1) .eq. 2 .or. iclst(i2) .eq. 2 ) then
!        cross = thcross
!       elseif( iclst(i1) .eq. 3 .or. iclst(i2) .eq. 3 ) then
!        cross = thcross
!       elseif( iclst(i1) .eq. 4 .or. iclst(i2) .eq. 4 ) then
!        cross = across
!       end if
!      end if
!-----------------------------------------------------------8/17 sawada

!      xdst = sqrt( cross / pi )

      rdstsq = ( r(1,i1) - r(1,i2) )**2
     &       + ( r(2,i1) - r(2,i2) )**2
     &       + ( r(3,i1) - r(3,i2) )**2

      rdst = dsqrt( rdstsq )

      if( rdst .lt. xdst) then
         ichannel = 1
!======================================================================!

      if( iclst(i1) .ne. 0 .or. iclst(i2) .ne. 0 )goto 199
      if( ininclst(i1) .ne. 0 .or. ininclst(i2) .ne. 0 )goto 199   ! 2022/10/11 yamaguchi
!      if( iclcoll .ne. 0 )goto 199

!=======================================================================
!     Start of Inelastic Collisions, at least one pion production

      if( x1 .gt. sigel/sig ) then

          ichannel = 0
!-------------------------------------------------------------sawada8/27
          if( srt .le. 2.0d0 * rmass/1000.0 + pmass) return
!-----------------------------------------------------------------------

!=======================================================================
!     [1] Nucleon-Nucleon Inelastic Collisions
!             (1-1) N + N -> N + D
!             (1-2) N + N -> N + R
!             (1-3) N + N -> D + D
         if( id1 .eq. 1 .and. id2 .eq. 1 ) then
             signd = 0.0d0
             signr = 0.0d0
             sigdd = 0.0d0
!-----------------------------------------------------------------------
!             (1-1) N + N -> N + D
!-----------------------------------------------------------------------
             z11 = soo(1,srt)
             z10 = soo(2,srt)

             if( iz1 .eq. iz2 ) then
                signd = 2.0*z11+z10
             else
                signd = z11 + 0.5*z10
             end if
!-----------------------------------------------------------------------
             if( x1 .lt. ( sigel + signd ) / sig ) then
                ichannel = 2
                idd1 = 2
                idd2 = 1
                izz1 = iz1
                izz2 = iz2
                if( x2 .gt. 0.5d0 ) then
                   izz1 = iz2
                   izz2 = iz1
                end if

!-------------------------------------------------------------sawada8/27
!     in INC.f,rmass=938.93[MeV],but in inelastic part,rmass=0.938[GeV].
                dem1 = rmass/1000.0 + pmass
                dem2 = rmass/1000.0
!-----------------------------------------------------------------------

                call resmass(1,idd1,idd2,srt,dem1,dem2,gam0,fqrq)
                if( ( iz1 + iz2 .ne. 1 ) .and. 
     &              ( x1 .lt. ( sigel + 0.5d0*z11+z10)/sig)) then
                    izz1 = 3 * iz1 - 1
                    izz2 = iz2 - 2 * iz1 + 1
                end if
                goto 100
             end if
!-----------------------------------------------------------------------
!             (1-2) N + N -> N + R
!-----------------------------------------------------------------------
             z01 = soo(3,srt)
             signr = fnnr * z01
!-----------------------------------------------------------------------
             if( x1 .lt. ( sigel + signd + signr ) / sig ) then
                ichannel = 4
                idd1 = 3
                idd2 = 1
                izz1 = iz1
                izz2 = iz2
                if( x2 .gt. 0.5 ) then
                   izz1 = iz2
                   izz2 = iz1
                end if

!-------------------------------------------------------------sawada8/27
                dem1 = rmass/1000.0 + pmass
                dem2 = rmass/1000.0
!-----------------------------------------------------------------------

                call resmass(1,idd1,idd2,srt,dem1,dem2,gam0,fqrq)
                goto 100
             end if
!-------------------------------------------------------------sawada8/27
             if( srt .le. 2.0*rmass/1000.0+2.0*pmass ) return
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!             (1-3) N + N -> D + D
!-----------------------------------------------------------------------
             sigdd = fndd * z01
             if( x1 .lt. ( sigel + signd + signr + sigdd ) / sig ) then
                ichannel = 6
                idd1 = 2
                idd2 = 2
                izz1 = iz1
                izz2 = iz2
!-------------------------------------------------------------sawada8/27
                dem1 = rmass/1000.0 + pmass
                dem2 = rmass/1000.0 + pmass
!-----------------------------------------------------------------------
                call resmass(1,idd1,idd1,srt,dem1,dem2,gam0,fqrq)
                if( iz1+iz2 .eq. 1 ) then
                   if( x1 .lt. ( sigel + signd + signr
     &                         + 7.0d0/10.0d0*sigdd ) /sig ) then
                      izz1 =  2
                      izz2 = -1
                   end if
                else
                   if( x1 .lt. ( sigel + signd + signr
     &                         + 3.0d0/7.0d0*sigdd ) / sig ) then
                      izz1 = ( iz1 + iz2 ) / 2 + 1
                      izz2 = ( iz1 + iz2 ) / 2 - 1
                   end if
                end if
                goto 100
             end if
!-----------------------------------------------------------------------
             return
         end if
!=======================================================================
!     [2] Nucleon-Delta; Inelastic Collisions
!             (2-1) N + D -> N + N
!             (2-2) N + D -> D + D
!             (2-3) N + D -> R + D
         if( id1 + id2 .eq. 3 ) then
            sigdn = 0.0d0
            sigdd = 0.0d0
            sigdr = 0.0d0
!-----------------------------------------------------------------------
!           (2-1) N + D -> N + N
            if( ( iz1 + iz2 + 3 ) / 3 .eq. 1 ) then
!-----------------------------------------------------------------------
               signd = 0.0d0
               z11 = soo(1,srt)
               z10 = soo(2,srt)
               if( iz1 + iz2 .eq. 1 ) then
                  signd = 0.5d0*z11 + 0.25*z10
                  sigiv = 2.0d0
               end if
               if( iz1 .eq. iz2 ) then
                  signd = 1.5d0*z11
                  sigiv = 4.0d0
               end if
               if( ( iz1 .ne. iz2 ) .and. ( iz1 + iz2 .ne. 1 ) ) then
                  signd = 0.5d0*z11 + z10
                  sigiv = 4.0d0
               end if
!-------------------------------------------------------------sawada8/27
               pidn  = pideno(2,rmass/1000.0+pmass,rmass/1000.0,srt)
!-----------------------------------------------------------------------
               sigdn = signd*prf**2/pri**2/pidn/sigiv
!-----------------------------------------------------------------------
               if( x1 .lt. ( sigel + sigdn ) / sig ) then
                  ichannel = 3
                  idd1 = 1
                  idd2 = 1
                  izz1 = iz1
                  izz2 = iz2
!-------------------------------------------------------------sawada8/27
                  dem1 = rmass/1000.0
                  dem2 = rmass/1000.0
!-----------------------------------------------------------------------
                  if( ( iz1+iz2 .eq. 0 ) .or.
     &                ( iz1+iz2 .eq. 2 ) ) then
                     izz1 = ( iz1 + iz2 ) / 2
                     izz2 = ( iz1 + iz2 ) / 2
                  end if
                  goto 100
               end if
            end if
!-------------------------------------------------------------sawada8/27
            if( srt .lt. 2.0*rmass/1000.0+2.0*pmass ) return
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!           (2-2) N + D -> D + D
            if( inddd .ne. 0 ) then
!-----------------------------------------------------------------------
!-------------------------------------------------------------sawada8/27
               srtd = srt - ( em1 + em2 - 2.0d0 * rmass )/1000.0
!-----------------------------------------------------------------------
               z11 = soo(1,srtd)
               z10 = soo(2,srtd)
               if( ( iz1+iz2 .eq. 3 ) .or. ( iz1+iz2 .eq. -1 ) ) then
                  sigdd = 1.5*z11
               else
                  sigdd = 2.0*z11+z10
               end if
!-----------------------------------------------------------------------
               if( x1 .lt. ( sigel + sigdn + sigdd ) / sig ) then
                     ichannel = 8
                     idd1 = 2
                     idd2 = 2
                  if( id1 .eq. 1 ) then
                     izz1 = iz1
                     izz2 = iz2
                  else
                     izz1 = iz2
                     izz2 = iz1
                  end if
!-------------------------------------------------------------sawada8/27
                     dem1 = rmass/1000.0 + pmass
                     dem2 = rmass/1000.0 + pmass
!-----------------------------------------------------------------------
                  call resmass(1,idd1,idd1,srt,dem1,dem2,gam0,fqrq)
               if( x1 .gt. ( sigel + sigdn + 1.5 * z11 ) / sig ) then
                  if( iz1 + iz2 .eq. 1 ) then
                     izz1 =  3 * iz1 - 1
                     izz2 =  iz2 - 2 * iz1 +1
                  else if( iz1 .eq. iz2 ) then
                     izz1 = ( iz1 + iz2 ) / 2 + 1
                     izz2 = ( iz1 + iz2 ) / 2 - 1
                  else
                     izz1 = ( iz1 + iz2 ) / 2
                     izz2 = ( iz1 + iz2 ) / 2
                  end if
               end if
                  goto 100
               end if
            end if
!-----------------------------------------------------------------------
!           (2-3) N + D -> R + D
            if( indrd .ne. 0 ) then
!-----------------------------------------------------------------------
               z01   = soo(3,srtd)
               sigdr = 1.5 * z01
!-----------------------------------------------------------------------
               if( x1 .lt.
     &              ( sigel+sigdn+sigdd+sigdr ) / sig ) then
                  ichannel = 12
                  idd1 = 3
                  idd2 = 2
                  if( id1 .eq. 1 ) then
                     izz1 = iz1
                     izz2 = iz2
                  else
                     izz1 = iz2
                     izz2 = iz1
                  end if
!-------------------------------------------------------------sawada8/27
                     dem1 = rmass/1000.0 + pmass
                     dem2 = rmass/1000.0 + pmass
!-----------------------------------------------------------------------
                  call resmass(1,idd1,idd2,srt,dem1,dem2,gam0,fqrq)
                  goto 100
               end if
            end if
!-----------------------------------------------------------------------
            return
         end if
!-----------------------------------------------------------------------
!=======================================================================
!        [3] Nucleon - Nstar  ;  Inelastic Collisions
!                        (3-1) N + R -> N + N
!                        (3-2) N + R -> D + R
!                        (3-3) N + R -> R + R
         if( ( id1 + id2 .eq. 4 ) .and. ( id1 .ne. 2 ) ) then
!=======================================================================
            sigrn = 0.0
            sigdr = 0.0
            sigrr = 0.0
!-----------------------------------------------------------------------
!           (3-1) N + R -> N + N
!-----------------------------------------------------------------------
            z01   = soo(3,srt)
            signr = fnnr * z01
!-------------------------------------------------------------sawada8/27
            pidn  = pideno(3,rmass/1000.0+pmass,rmass/1000.0,srt)
!-----------------------------------------------------------------------
            if( iz1 .eq. iz2 ) then
               sigiv = 2.0
            else
               sigiv = 1.0
            end if
            sigrn = signr * prf**2 / pri**2 / pidn / sigiv
!-----------------------------------------------------------------------
            if( x1 .lt. ( sigel + sigrn ) / sig ) then
               ichannel = 5
               idd1 = 1
               idd2 = 1
               izz1 = iz1
               izz2 = iz2
!-------------------------------------------------------------sawada8/27
               dem1 = rmass/1000.0
               dem2 = rmass/1000.0
!-----------------------------------------------------------------------
               goto 100
            end if
!-------------------------------------------------------------sawada8/27
!            if( srt .le. 2.0 * rmass + 2.0 * pmass ) return
            if( srt .le. 2.0 * rmass/1000.0 + 2.0 * pmass ) return
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!           (3-2) N + R -> D + R
            if( inrdr .ne. 0 ) then
!-----------------------------------------------------------------------
!-------------------------------------------------------------sawada8/27
               srtd  = srt - ( em1 + em2 - 2.0 * rmass )/1000.0
!-----------------------------------------------------------------------
               z11   = soo(1,srtd)
               z10   = soo(2,srtd)
               sigdr = 2.0 * z11 + z10
!-----------------------------------------------------------------------
               if( x1 .lt. ( sigel + sigrn + sigdr ) / sig ) then
                  ichannel = 10
                  idd1 = 2
                  idd2 = 3
                  if( id1 .eq. 1 ) then
                     izz1 = iz1
                     izz2 = iz2
                  else
                     izz1 = iz2
                     izz2 = iz1
                  end if
!-------------------------------------------------------------sawada8/27
                  dem1 = rmass/1000.0 + pmass
                  dem2 = rmass/1000.0 + pmass
!-----------------------------------------------------------------------
                  call resmass(1,idd2,idd1,srt,dem2,dem1,gam0,fqrq)
                  if( x1 .gt. ( sigel+sigrn+1.5*z11 ) / sig ) then
                     if( iz1 + iz2 .eq. 1 ) then
                        izzd = izz1
                        izz1 = izz2
                        izz2 = izzd
                     else 
                        izz1 = 3 * ( iz1 + iz2 ) / 2 -1
                        izz2 = iz1 + iz2 - izz1
                     end if
                  end if
                  goto 100
               end if
            end if
!-----------------------------------------------------------------------
!           (3-3) N + R -> R + R
            if( inrrr .ne. 0 ) then
!-----------------------------------------------------------------------
               z01   = soo(3,srtd)
               sigrr = 1.5 * z01
!-----------------------------------------------------------------------
               if( x1 .lt.
     &              ( sigel+sigrn+sigdr+sigrr ) / sig ) then
                  ichannel = 14
                  idd1 = 3
                  idd2 = 3
                  if( id1 .eq. 1 ) then
                     izz1 = iz1
                     izz2 = iz2
                  else
                     izz1 = iz2
                     izz2 = iz1
                  end if
!-------------------------------------------------------------sawada8/27
                  dem1 = rmass/1000.0 + pmass
                  dem2 = rmass/1000.0 + pmass
!-----------------------------------------------------------------------
                  call resmass(1,idd1,idd2,srt,dem1,dem2,gam0,fqrq)
                  goto 100
               end if
            end if
!-----------------------------------------------------------------------
            return
         end if
!-----------------------------------------------------------------------
!=======================================================================
!        [4] Delta - Delta  ;  Inelastic Collisions
!                        (4-1) D + D -> N + N
!                        (4-2) D + D -> D + N
         if( ( id1 .eq. 2 ) .and. ( id2 .eq. 2 ) ) then
!======================================================================
            sigdnn = 0.0
            sigdnd = 0.0
!-----------------------------------------------------------------------
!           (4-1) D + D -> N + N
            if( ( iz1 + iz2 + 3 ) / 3 .eq. 1 .and.
     &           fndd .gt. 0.0 ) then
!-----------------------------------------------------------------------
               z01    = soo(3,srt)
               sigdd  = fndd * z01
               if( ( iz1 + iz2 .eq. 1 ) .or. ( iz1 .eq. iz2 ) ) then
                  sigiv = 4.0
               else 
                  sigiv = 8.0
               end if
               sigdnn = sigdd * prf**2 / pri**2 / sigiv
!-----------------------------------------------------------------------
               if( x1 .lt. ( sigel + sigdnn ) / sig ) then
                  ichannel = 7
                  idd1 = 1
                  idd2 = 1
!-------------------------------------------------------------sawada8/27
                  dem1 = rmass/1000.0
                  dem2 = rmass/1000.0
!-----------------------------------------------------------------------
                  if( iz1 + iz2 .eq. 1 ) then
                     izz1 = 1
                     izz2 = 0
                  else
                     izz1 = ( iz1 + iz2 ) / 2
                     izz2 = ( iz1 + iz2 ) / 2
                  end if
                  goto 100
               end if
            end if
!-----------------------------------------------------------------------
!           (4-2) D + D -> D + N
            if( inddd .ne. 0 ) then
!-----------------------------------------------------------------------
               if( ( iz1 + iz2 .ge. -1 ) .and.
     &              ( iz1 + iz2 .le. 3 ) ) then
!-------------------------------------------------------------sawada8/27
                  dem1 = rmass/1000.0
                  if( x2 .gt. 0.5 ) then
                     dem2 = em2/1000.0
                  else
                     dem2 = em1/1000.0
                  end if
                  srtd  = srt - ( dem2 - rmass/1000.0 )
!-----------------------------------------------------------------------
                  z11   = soo(1,srtd)
                  z10   = soo(2,srtd)
                  prfd  = sqrt( ( srt**2 - dem1**2 - dem2**2 )**2
     &                 - 4.0 * ( dem1 * dem2 )**2 ) / ( 2.0 * srt )
                  if( ( ( iz1 + iz2 .eq. 3 ) .or.
     &                 ( iz1 + iz2 .eq. -1 ) ) .or.
     &                 ( ( iz1 + iz2 .eq. 1 ) .and.
     &                 ( iz1 .eq. 1 .or. iz2 .eq. 1 ) ) ) then
                     signdd = 3.0 * z11
                     icc    = 1
                  elseif( iz1 + iz2 .eq. 1 ) then
                     signdd = z11 + 2.0 * z10
                     icc    = 2
                  else
                     signdd = 2.0 * z11 + z10
                     icc    = 3
                  end if
!-------------------------------------------------------------sawada8/27

                  pidn   = pideno(2,rmass/1000.0+pmass,dem2,srt)
!-----------------------------------------------------------------------
                  if( iz1 .eq. iz2 ) then
                     sigiv = 1.0
                  else
                     sigiv = 2.0
                  end if
                  sigdnd = signdd * prfd**2 / pri**2 / pidn /sigiv
!-----------------------------------------------------------------------
                  if( x1 .lt. ( sigel + sigdnn + sigdnd ) / sig ) then
                     ichannel = 9
                     idd1 = 1
                     idd2 = 2
                     izz1 =  iz1
                     izz2 =  iz2
                     if( ( icc .eq. 2 ) .or.
     &                    ( ( icc .eq. 3 ) .and. 
     &                    ( x1 .gt.
     &                    ( sigel+sigdnn+1.5*z11 ) / sig ) ) ) then
                        if( iz1 + iz2 .eq. 1 ) then
                           if( x2 .gt. 0.5 ) then
                              izz1 =  1
                              izz2 =  0
                           else
                              izz1 =  0
                              izz2 =  1
                           end if
                        else if( iz1 .eq. iz2 ) then
                           izz1 = ( iz1 + iz2 ) / 2 + 1
                           izz2 = ( iz1 + iz2 ) / 2 - 1
                        else
                           izz1 = ( iz1 + iz2 ) / 2
                           izz2 = ( iz1 + iz2 ) / 2
                        end if
                     end if
                     if( ( izz1 + 2 ) / 2 .ne. 1 ) then
                        izzd = izz1
                        izz1 = izz2
                        izz2 = izzd
                     end if
                     goto 100
                  end if
               end if
            end if
!-----------------------------------------------------------------------
            return
         end if
!-----------------------------------------------------------------------
!=======================================================================
!        [5] Delta - Nstar  ;  Inelastic Collisions
!                        (5-1) D + R -> N + R
!                        (5-2) R + D -> N + D
         if( ( id1 .eq. 2 ) .and. ( id2 .eq. 3 ) .or.
     &       ( id1 .eq. 3 ) .and. ( id2 .eq. 2 ) ) then
!=======================================================================
            sigdnr = 0.0
            sigrnd = 0.0
!-----------------------------------------------------------------------
!           (5-1) D + R -> N + R
            if( inrdr .ne. 0 ) then
!-----------------------------------------------------------------------
!-------------------------------------------------------------sawada8/27
               if( ( iz1 + iz2 + 3 ) / 3 .eq. 1 ) then
                  dem1 = rmass/1000.0
                  if( id1 .eq. 2 ) then
                     dem2 = em2/1000.0
                     izz1 = iz1
                     izz2 = iz2
                  else
                     dem2 = em1/1000.0
                     izz1 = iz2
                     izz2 = iz1
                  end if
                  srtd  = srt - ( dem2 - rmass/1000.0 )
!-----------------------------------------------------------------------
                  z11   = soo(1,srtd)
                  z10   = soo(2,srtd)
                  prfd  = sqrt( ( srt**2 - dem1**2 - dem2**2 )**2
     &                 - 4.0 * ( dem1 * dem2 )**2 ) / ( 2.0 * srt )
                  if( iz1 + iz2 .eq. 1 ) then
                     sigrdr = 2.0 * z11 + z10
                  else if( iz1 .eq. iz2 ) then
                     sigrdr = 1.5 * z11
                  else
                     sigrdr = 0.5 * z11 + z10
                  end if
!-------------------------------------------------------------sawada8/27
                  pidn   = pideno(2,rmass/1000.0+pmass,dem2,srt)
!-----------------------------------------------------------------------
                  sigiv  = 2.0
                  sigdnr = sigrdr * prfd**2 / pri**2 / pidn / sigiv
!-----------------------------------------------------------------------
                  if( x1 .lt. ( sigel + sigdnr ) / sig ) then
                     ichannel = 11
                     idd1 = 1
                     idd2 = 3
                     if( ( iz1 + iz2 .eq. 1 ) .and.
     &                    ( x1 .gt. ( sigel + 1.5 * z11 ) / sig ) ) then
                        izzd = izz1
                        izz1 = izz2
                        izz2 = izzd
                     end if
                     if( ( iz1 + iz2 .ne. 1 ) .and.
     &                    ( iz1 .ne. iz2 ) ) then
                        izz1 = ( iz1 + iz2 ) / 2
                        izz2 = ( iz1 + iz2 ) / 2
                     end if
                     goto 100
                  end if
               end if
            end if
!-----------------------------------------------------------------------
!           (5-2) R + D -> N + D
            if( indrd .ne. 0 ) then
!-----------------------------------------------------------------------
!-------------------------------------------------------------sawada/8/27
               if( id1 .eq. 3 ) then
                  dem1 = rmass/1000.0
                  dem2 = em2/1000.0
                  izz1 = iz1
                  izz2 = iz2
               else
                  dem1 = rmass/1000.0
                  dem2 = em1/1000.0
                  izz1 = iz2
                  izz2 = iz1
               end if
               srtd  = srt - ( dem2 - rmass/1000.0 )
!-----------------------------------------------------------------------

               z01   = soo(3,srtd)
               prfd  = sqrt( ( srt**2 - dem1**2 - dem2**2 )**2
     &              - 4.0 * ( dem1 * dem2 )**2 ) / ( 2.0 * srt )
               sigdrd = 1.5 * z01
!-------------------------------------------------------------sawada8/27
               pidn   = pideno(3,rmass/1000.0+pmass,dem2,srt)
!-----------------------------------------------------------------------
               sigiv  = 1.0
               sigrnd = sigdrd * prfd**2 / pri**2 / pidn / sigiv
!-----------------------------------------------------------------------
               if( x1 .lt. ( sigel + sigdnr + sigrnd ) / sig ) then
                  ichannel = 13
                  idd1 = 1
                  idd2 = 2
                  goto 100
               end if
            end if
!-----------------------------------------------------------------------
            return
         end if
!-----------------------------------------------------------------------
!=======================================================================
!        [6] Nstar - Nstar  ;  Inelastic Collisions
!                        (6-1) R + R -> N + R
         if( ( id1 .eq. 3 ) .and. ( id2 .eq. 3 ) ) then
!=======================================================================
            sigrnr = 0.0
!-----------------------------------------------------------------------
!           (6-1) R + R -> N + R
            if( inrrr .ne. 0 ) then
!-----------------------------------------------------------------------
!-------------------------------------------------------------sawada8/27
               dem1 = rmass/1000.0
               dem2 = em2/1000.0
               srtd  = srt - ( dem2 - rmass/1000.0 )
!-----------------------------------------------------------------------
               z01   = soo(3,srtd)
               prfd  = sqrt( ( srt**2 - dem1**2 - dem2**2 )**2
     &              - 4.0 * ( dem1 * dem2 )**2 ) / ( 2.0 * srt )
               sigrrr = 1.5 * z01
!-------------------------------------------------------------sawada8/27
               pidn   = pideno(3,rmass/1000.0+pmass,dem2,srt)
!-----------------------------------------------------------------------
               if( iz1 .eq. iz2 ) then
                  sigiv = 0.5
               else
                  sigiv = 1.0
               end if
               sigrnr = sigrrr * prfd**2 / pri**2 / pidn / sigiv
!-----------------------------------------------------------------------
               if( x1 .lt. ( sigel + sigrnr ) / sig ) then
                  ichannel = 15
                  izz1 = iz1
                  izz2 = iz2
                  idd1 = 1
                  idd2 = 3
                  goto 100
               end if
            end if
!-----------------------------------------------------------------------
            return
         end if
!-----------------------------------------------------------------------
!=======================================================================
         return
!=======================================================================
  100    continue

         if( x2 .gt. 0.5 ) then
            em1 = dem1
            em2 = dem2
            iz1 = izz1
            iz2 = izz2
            id1 = idd1
            id2 = idd2
         else
            em1 = dem2
            em2 = dem1
            iz1 = izz2
            iz2 = izz1
            id1 = idd2
            id2 = idd1
         end if

         pr  = sqrt( ( srt**2 - em1**2 - em2**2 )**2
     &                - 4.0 * ( em1 * em2 )**2 ) / ( 2.0 * srt )

      pr = pr * 1000.0  
      em1 = em1*1000.0
      em2 = em2*1000.0

!----------------------------------------------------------------------


!======================================================================!

      end if

199   continue
!======================================================================!

!-----------------------------------------------------------------------
!     angular distribution
!-----------------------------------------------------------------------

         if( ichannel .eq. 1 ) then

!-------------------------------------------------------------sawada8/27
            ta  = -2.0 * (pra)**2
!-----------------------------------------------------------------------
            x   = rn(0)

!-----------------------------------------------------------------------
!           New Cugnon
!-----------------------------------------------------------------------

            if( iang .eq. 1 ) then

!-----------------------------------------------------------------------
!              p + n collision
!-----------------------------------------------------------------------

               if( iz1 .ne. iz2 .or. id1 .ne. id2 ) then

                  if( pl .le. 0.225 ) then

                     c1 = 1.0 - 2.0 * x

                  else

                     if( pl .le. 0.6d0 ) then
                        a = 6.2d0 * ( pl - 0.225d0 ) / 0.375d0
                     elseif( pl .le. 1.6 ) then
                        a =  - 1.63d0 * pl + 7.16d0
                     elseif( pl .le. 2.0d0 ) then
                        a = 5.5d0 * pl**8 / ( 7.7d0 + pl**8 )
                     else
                        a = 5.34d0 + 0.67d0 * ( pl - 2.0d0 )
                     end if

                     tt1 = log( ( 1.0d0 - x )
     &                    * exp( 2.0d0 * a * ta ) + x ) / a

                     bprob = 1.0d0
                     if( pl .gt. 0.8d0 ) bprob = 0.64d0 / pl / pl

                     bprob = bprob / ( 1.0d0 + bprob )

                     if( rn(0) .gt. bprob ) then
                        c1 = 1.0d0 - tt1 / ta
                     else
                        c1 = - 1.0d0 + tt1 / ta
                     end if

                     if( abs(c1) .gt. 1.0d0 ) c1 = 1.0 - 2.0 * x
                          
                  end if

!-----------------------------------------------------------------------
!              p + p collision
!-----------------------------------------------------------------------

               else

                  if( pl .le. 0.0d0 ) then
                     c1 = 1.0d0 - 2.0d0 * x
                  else                      
                     if( pl .le. 2.0d0 ) then
                        a = 5.5d0 * pl**8 / ( 7.7d0 + pl**8 )
                     else
                        a = 5.334d0 + 0.67d0 * ( pl - 2.0d0 )
                     end if
                     tt1 = log( (1.0d0-x) * exp(2.0d0*a*ta) + x )
     &                    / a
                     c1  = 1.0d0 - tt1 / ta
                     if( abs(c1) .gt. 1.0d0 ) c1 = 2.0 * x - 1.0
                  end if

               end if

!-----------------------------------------------------------------------
!     Old Cugnon
!-----------------------------------------------------------------------

            elseif( iang .eq. 0 ) then

!-------------------------------------------------------------sawada8/27
               as  = ( 3.65d0 * asrt )**6
!-----------------------------------------------------------------------
               a   = 6.0d0 * as / (1.0d0 + as)
!-------------------------------------------------------------sawada8/27
               ta  = -2.0d0 * (pra)**2
!-----------------------------------------------------------------------
               x   = rn(0)
               t1  = dlog( (1.0d0-x) * exp(2.0d0*a*ta) + x )  /  a
               c1  = 1.0d0 - t1/ta
               if( abs(c1) .gt. 1.0d0 ) c1 = 2.0d0 * x - 1.0d0
               
            end if

!-----------------------------------------------------------------------

         else

!-----------------------------------------------------------------------

            srtd  = srt
            asrtd = asrt
            asrtn = asrt
            prad  = pra
!-----------------------------------------------------------------------

            if( ichannel .eq.  4 .or. ichannel .eq.  5 .or.
     &           ichannel .eq. 12 .or. ichannel .eq. 13 .or.
     &           ichannel .eq. 14 .or. ichannel .eq. 15 ) then

               srtd  = srt  - ( smass - dmass )
               asrtn = asrt - ( smass - dmass )

               if( asrtn .gt. 0.0 ) then

                  asrtd = asrtn
!-------------------------------------------------------------sawada8/27
                  prad  = sqrt(( srtd**2 - ( p(5,i1)/1000.0 )**2
     &                       - ( p(5,i2)/1000.0 )**2 )**2 
     &                       - 4.0*( p(5,i1)/1000.0*p(5,i2)/1000.0 )**2)
     &                    / ( 2.0 * srtd )                        ![GeV]
!-----------------------------------------------------------------------
               end if

            end if

!-----------------------------------------------------------------------

            if( rn(0) .lt. 0.5 ) then

               if( asrtn .gt. 0.0 ) then
                  as  = ( 3.65 * asrtd )**6
                  a   = as / (1.0 + as) * srtd**4 * 0.14
                  ta  = -2.0 * prad**2
                  x   = rn(0)
                  t1  = log( (1-x)*exp( max(-50.d0,2.*a*ta) ) + x ) / a
                  c1  = 1.0 - t1/ta
                  if(abs(c1).gt.1.0) c1 = 2.0 * x - 1.0
               else
                  c1 = 2.0 * rn(0) - 1.0
               end if

            else

               if( srtd .lt. 2.14 ) then

                  c1 = 2.0 * rn(0) - 1.0

               else

                  if( srtd .gt. 2.4 ) then

                     b1 = 0.06
                     b3 = 0.4

                  else

                     b1 = 29.0286 - 23.749  * srtd + 4.86549*srtd**2
                     b3 =-30.3283 + 25.5257 * srtd - 5.30129*srtd**2

                  end if

                  pp3 = b1 / (3. * b3)
                  qq3 = 0.5 * (0.5 - rn(0)) / b3
                  pq3 = sqrt(qq3**2 + pp3**3)
                  uu  = (-qq3 + pq3 )**(1./3.)
                  vv  = ( qq3 + pq3 )**(1./3.)
                  c1  = uu - vv
                  if( abs(c1) .gt. 1. ) c1 = c1 / abs(c1)

               end if

            end if

         end if

!-----------------------------------------------------------------------
!     set the new momentum coordinates
!-----------------------------------------------------------------------

            if( pcm(1) .eq. 0.0 .and. pcm(2) .eq. 0.0 ) then
               t2 = 0.0
            else
               t2 = atan2(pcm(2),pcm(1))
            end if

            abst2 = sqrt(pcm(2)**2.+pcm(1)**2.)

!            iball = 0

!            if(iclst(1) .ne. 0) then
!              if( i1 .le. 4) iball = i1
!              if( i2 .le. 4) iball = i2

 !             if( iball .eq. 0 .or. iclst(iball) .lt.3) then
 !                 t1 = 2.0 * pi *rn(0)
 !             else
 !                 !t1 = pi * (0.5-rn(0))
 !                 t1 = 2.0 * pi *rn(0)
 !             end if

 !           else

            t1  = 2.0 * pi * rn(0)

!            end if

!            if( pcm(1) .eq. 0.0 .and. pcm(2) .eq. 0.0 ) then
!               t2 = 0.0
!            else
!               t2 = atan2(pcm(2),pcm(1))
!            end if

            s1   = dsqrt( 1.0 - c1**2 )
            s2  =  dsqrt( 1.0 - c2**2 )
            ct1  = dcos(t1)
            st1  = dsin(t1)
            ct2  = dcos(t2)
            st2  = dsin(t2)
            ss   = c2 * s1 * ct1  +  s2 * c1

            pcm(1) = pr * ( ss*ct2 - s1*st1*st2 )
            pcm(2) = pr * ( ss*st2 + s1*st1*ct2 )
            pcm(3) = pr * ( c1*c2  - s1*s2 *ct1 )

         end if

!-------------------------------------------------------------sawada8/27
*-----------------------------------------------------------------------
*     charge and state identification
*-----------------------------------------------------------------------
         if( ichannel .gt. 1 ) then
               inds(i1) = id1
               inds(i2) = id2
            if( inds(i1).eq.1 ) then
               inuc(i1) = 1
            else
               inuc(i1) = 0
            end if
            if( inds(i2).eq.1 ) then
               inuc(i2) = 1
            else
               inuc(i2) = 0
            end if
               ichg(i1) = iz1
               ichg(i2) = iz2
               p(5,i1)  = em1
               p(5,i2)  = em2
         end if

*-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Lorentz-transformation into reference frame
!-----------------------------------------------------------------------

      e1cm    = dsqrt( em1**2 + pcm(1)**2 + pcm(2)**2 + pcm(3)**2 )
      p1beta  = pcm(1)*beta(1) + pcm(2)*beta(2) + pcm(3)*beta(3)
      transf  = gamma * ( gamma * p1beta / (gamma + 1.0d0) + e1cm )
      
      p(1,i1) = beta(1) * transf + pcm(1)
      p(2,i1) = beta(2) * transf + pcm(2)
      p(3,i1) = beta(3) * transf + pcm(3)
      
      p(4,i1) = dsqrt( p(5,i1)**2 + p(1,i1)**2
     &                           + p(2,i1)**2 + p(3,i1)**2 )

         e2cm    = dsqrt(em2**2 + pcm(1)**2 + pcm(2)**2 + pcm(3)**2)

         transf  = gamma * (-gamma * p1beta / (gamma + 1.0) + e2cm)

         p(1,i2) = beta(1) * transf - pcm(1)
         p(2,i2) = beta(2) * transf - pcm(2)
         p(3,i2) = beta(3) * transf - pcm(3)

         p(4,i2) = dsqrt( p(5,i2)**2 + p(1,i2)**2
     &        + p(2,i2)**2
     &        + p(3,i2)**2 )
!******************************************************sawada7/7
!     to deaccelerate nucleon by 45*xMeV when it collides cluster
              px1  = p(1,i1)
              py1  = p(2,i1)
              pz1  = p(3,i1)
              e1   = p(4,i1)
              em1  = p(5,i1)
              px2  = p(1,i2)
              py2  = p(2,i2)
              pz2  = p(3,i2)
              e2   = p(4,i2)
              em2  = p(5,i2)
!              uu   = 0.0d0
!              u    = 0.0d0

             if( em1 .lt. em2+1.0d-3 .and. em1 .gt. em2-1.0d-3 ) then

              px1  = p(1,i1)
              py1  = p(2,i1)
              pz1  = p(3,i1)
              e1   = p(4,i1)
              px2  = p(1,i2)
              py2  = p(2,i2)
              pz2  = p(3,i2)
              e2   = p(4,i2)
              
!             if( iclcoll .eq. 1 )then
!
!              if( iclst(i1) .eq. 1 )then
!               uu = 45.0d0
!              elseif( iclst(i1) .eq.2 .or. iclst(i1) .eq.3 )then
!               uu = 90.0d0
!              elseif( iclst(i1) .eq. 4 )then
!               uu = 135.0d0
!              end if

             elseif( em1 .gt. em2 )then

              uuu = (em1-em2)/rmass*upot
              
              unitx = px2/dsqrt( px2**2 + py2**2 + pz2**2 )
              unity = py2/dsqrt( px2**2 + py2**2 + pz2**2 )
              unitz = pz2/dsqrt( px2**2 + py2**2 + pz2**2 )

              e2 = p(4,i2)-uuu
              et2 = e2-em2
                if( et2 .le. 0.0d0)then
                 ichannel = 99
                 return
                end if
              pabs = dsqrt( et2*(et2 + 2.0d0*em2) )

              p(1,i2) = pabs*unitx
              p(2,i2) = pabs*unity
              p(3,i2) = pabs*unitz
              p(4,i2) = dsqrt( p(5,i2)**2 + p(1,i2)**2
     &                          + p(2,i2)**2
     &                          + p(3,i2)**2 )

!             elseif( iclcoll .eq. 2 )then
!
!              if( iclst(i2) .eq. 1 )then
!               u = 45.0d0
!              elseif( iclst(i2) .eq.2 .or. iclst(i2) .eq.3 )then
!               u = 90.0d0
!              elseif( iclst(i2) .eq. 4 )then
!               u = 135.0d0
!              end if

             elseif( em1 .lt. em2 )then

              uuu = (em2-em1)/rmass*upot
              
              unitx = px1/dsqrt( px1**2 + py1**2 + pz1**2 )
              unity = py1/dsqrt( px1**2 + py1**2 + pz1**2 )
              unitz = pz1/dsqrt( px1**2 + py1**2 + pz1**2 )

              e1 = p(4,i1)-uuu
              et1 = e1-em1
                if( et1 .le. 0.0d0)then
                 ichannel = 99
                 return
                end if
              pabs = dsqrt( et1*(et1 + 2.0d0*em1) )

              p(1,i1) = pabs*unitx
              p(2,i1) = pabs*unity
              p(3,i1) = pabs*unitz
              p(4,i1) = dsqrt( p(5,i1)**2 + p(1,i1)**2
     &                          + p(2,i1)**2
     &                          + p(3,i1)**2 )
             end if

!**************************************
            t1 = p(4,i1) - p(5,i1)
            t2 = p(4,i2) - p(5,i2)
!-----------------------------------------------------------------------
!     deuteron    sawada8/27
!-----------------------------------------------------------------------
!           if( inds(i2) .eq. 1 )then
!             if( iclst(i2) .eq. 1 ) then
!                 if(  t2  .gt. ducpot ) then
!                 if(  t2  .gt. dupot ) then
!                    ipot(i2) = 1
!                 else
!                    ipot(i2) = 0
!                 end if
!-----------------------------------------------------------------------
!     triton    sawada8/27
!-----------------------------------------------------------------------
!             elseif( iclst(i2) .eq. 2 )then
!                 if(  t2  .gt. tucpot ) then
!                 if(  t2  .gt. thupot ) then
!                    ipot(i2) = 1
!                 else
!                    ipot(i2) = 0
!                 end if
!-----------------------------------------------------------------------
!    He-3    sawada8/27
!-----------------------------------------------------------------------
!             elseif( iclst(i2) .eq. 3 )then
!                 if(  t2  .gt. hucpot ) then
!                 if(  t2  .gt. thupot ) then
!                    ipot(i2) = 1
!                 else
!                    ipot(i2) = 0
!                 end if
!-----------------------------------------------------------------------
!     alpha    sawada8/27
!-----------------------------------------------------------------------
!             elseif( iclst(i2) .eq. 4 ) then
!                 if(  t2  .gt. aucpot ) then
!                 if(  t2  .gt. aupot ) then
!                    ipot(i2) = 1
!                 else
!                    ipot(i2) = 0
!                 end if
!-----------------------------------------------------------------------
!     neutron or proton
!-----------------------------------------------------------------------
!             elseif( iclst(i2) .eq. 0 ) then
!                 if( iz2 .eq. 1 ) then
!                    if( t2 .gt. ucpot ) then
!                       ipot(i2) = 1
!                    else
!                       ipot(i2) = 0
!                    end if
!                 elseif( iz2 .eq. 0 ) then
!                    if( t2 .gt. upot ) then
!                       ipot(i2) = 1
!                    else
!                       ipot(i2) = 0
!                    end if
!                 end if
!             elseif( iclst(i2) .eq. -1 ) then
!                       ipot(i2) = 0
!             end if
!           else
!               ipot(i2) = 1              ! pion,N* and delta are moving
!           end if
!-----------------------------------------------------------------------
!     save whether nucleon is above Potential or not
!-----------------------------------------------------------------------
!           if( inds(i1) .eq. 1 ) then
!             if( iclst(i1) .eq. 4 )then    !sawada
!               if( t1 .gt. aucpot ) then
!               if( t1 .gt. aupot ) then
!                  ipot(i1) = 1
!               else
!                  ipot(i1) = 0
!               end if
!             elseif( iclst(i1) .eq. 3 )then
!               if( t1 .gt. hucpot ) then
!               if( t1 .gt. thupot ) then
!                  ipot(i1) = 1
!               else
!                  ipot(i1) = 0
!               end if
!             elseif( iclst(i1) .eq. 2 )then
!               if( t1 .gt. tucpot ) then
!               if( t1 .gt. thupot ) then
!                  ipot(i1) = 1
!               else
!                  ipot(i1) = 0
!               end if
!             elseif( iclst(i1) .eq. 1 )then
!               if( t1 .gt. ducpot ) then
!               if( t1 .gt. dupot ) then
!                  ipot(i1) = 1
!               else
!                  ipot(i1) = 0
!               end if
!             elseif( iclst(i1) .eq. 0 )then
!              if( ichg(i1) .eq. 0 ) then
!                 if( t1 .gt. upot ) then
!                   ipot(i1) = 1
!                 else
!                   ipot(i1) = 0
!                 end if
!              else
!                 if( t1 .gt. ucpot ) then
!                   ipot(i1) = 1
!                elseif( t1 .lt. upot ) then
!                 else
!                   ipot(i1) = 0
!                 end if
!              end if
!             elseif( iclst(i1) .eq. -1 ) then
!                      ipot(i1) = 0
!             end if
!-----------------------------------------------------------------sawada
!           else
!                ipot(i1) = 1              ! pion,N* and delta are moving
!           end if

!-----------------------------------------------------------------------
!     Energetically forbidden collision ( Pauli Blocking )
!-----------------------------------------------------------------------
!-------------------------------------------------sawada7/23
            tfm1  = tfm
            tfm2  = tfm
            dtfm  = 2.0d0*tfm
            thtfm = 3.0d0*tfm
            atfm  = 4.0d0*tfm

            if( inds(i1) .eq. 1 )then
             if( iclst(i1) .eq. 0 )then
              tfm1 = tfm
             elseif( iclst(i1) .eq. 1 )then
              tfm1 = dtfm
             elseif( iclst(i1) .eq. 2 .or. iclst(i1) .eq. 3 )then
              tfm1 = thtfm
             elseif( iclst(i1) .eq. 4 )then
              tfm1 = atfm
             end if
            else
              tfm1 = 0.0d0
            end if
            if( inds(i2) .eq. 1 )then
             if( iclst(i2) .eq. 0 )then
              tfm2 = tfm
             elseif( iclst(i2) .eq. 1 )then
              tfm2 = dtfm
             elseif( iclst(i2) .eq. 2 .or. iclst(i2) .eq. 3 )then
              tfm2 = thtfm
             elseif( iclst(i2) .eq. 4 )then
              tfm2 = atfm
             end if
            else
              tfm2 = 0.0d0
            end if

            tfm1 = tfm
            tfm2 = tfm

         if( t1 .le. tfm1 .or. t2 .le. tfm2 ) then
            ichannel = 99
         end if
!         if( t1 .le. 0.0d0 .or. 
!     &        t2 .le. 0.0d0 ) then
!            ichannel = 99
!         end if

!-----------------------------------------------------------------------

      END SUBROUTINE crosww



!***********************************************************************
!*                                                                     *
      real*8 FUNCTION coul(nchgta,rt00)
!*                                                                     *
!*      Coulomb potential                                              *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------
!     proton emission from cascade process
!-----------------------------------------------------------------------

      coul = hbar / ee * dble(nchgta) / rt00 / 2.0d0

!-----------------------------------------------------------------------

      END FUNCTION coul

!***********************************************************************
!*                                                                     *
!      SUBROUTINE coales(icoales,lcoll,i1,nchgta,massba,rt00,r,p,inds,
!     &                  ipot,ichg,iclst,upot,ucpot,ein)
      SUBROUTINE coales(icoales,lcoll,i1,nchgta,massba,rt00,r,p,inds,
     &     ipot,ichg,iclst,upot,ucpot,ein,numclst,pcmclst,
     &     ininclst,nmasclst,inio)
!*                                                                     *
!*        Purpose:                                                     *
!*                                                                     *
!*           to calculate surface-coalescence process                  *
!*                                                                     *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      integer  i1,nchgta,massba,lcoll(30),ichg(nnn)
      integer  ipot(nnn),icoales(nnn),iclst(nnn),inds(nnn)
      integer  inio(nnn) ! 2022/10/5 yamaguchi
      real*8   r(3,nnn),p(6,nnn),upot,ucpot,rt00
      real*8   pcmclst(5,nnn)

      dimension  pclst(3)
      integer numclst(nnn),ininclst(nnn),nmasclst(nnn) ! 2022/10/5 yamaguchi
      
      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )
!-----------------------------------------------------------------------
      ducpot = dupot  + coul(nchgta,rt00)
      tucpot = thupot + coul(nchgta,rt00)
      hucpot = thupot + 2.0d0*coul(nchgta,rt00)
      aucpot = aupot  + 2.0d0*coul(nchgta,rt00)

!-----------------------------------------------------------------------
!     coalescence radius (MeV/c fm)
!-----------------------------------------------------------------------
! Au: 1000 determined (11.20)
! Ta: 1000 determined (11.19)
! Nb: 1000 determined (11.18)
! V:  1000 determined (11.17)
! Al: 1200 determined (11.18)
! C:  1500 determined (11.18)

      if( inds(i1) .ne. 1 )return
!      if( iclst(i1) .ne. 0 )return
      if( iclst(i1) .eq. -1 )return   ! 2022/10/5 yamaguchi

!-----------------------------------------------------------------------
!      ia2 = 0
!      ia3 = 0
!      ia4 = 0

!      do j = 1, nmaspr + nmasta

!         if( j .eq. i1 ) cycle

!         rsq = ( r(1,i1) - r(1,j) )**2 + ( r(2,i1) - r(2,j) )**2
!     &        + ( r(3,i1) - r(3,j) )**2

!         rij = sqrt( rsq )

!         psq = ( p(1,i1) - p(1,j) )**2 + ( p(2,i1) - p(2,j) )**2
!     &        + ( p(3,i1) - p(3,j) )**2

!         pij = sqrt( psq )

!         h = rij * pij

!--------------------------------------------------
!      cluster radius cut-off (2009.11.16)
!            NOTE! very sensitive parameter!
!--------------------------------------------------
!         if( ein .ge. 200.0d0 ) then
!            if( amas .le. 12.0d0 ) then
!               if( rij .gt. 2.1d0 ) h = h0*1.2d0
!            else
!               if( rij .gt. 1.5d0 ) h = h0*1.2d0
!            end if
!         else
!            if( rij .gt. 2.3d0 ) h = h0*1.2d0
!         end if
!--------------------------------------------------


!-----------------------------------------------------------------------
!         If elastic collision was occured at least one time and
!         distance of j-th nucleon is near the i1-th nucleon,
!         surface coalescence process is started
!-----------------------------------------------------------------------

!         if( h .lt. h0 .and. ipot(j) .eq. 0 .and. lcoll(4) .ne. 0 ) then

!            if( ia2 .eq. 0 ) then
!               ia2 = j
!            end if
!            if( ia2 .gt. 0 .and. j .ne. ia2 .and. ia3 .eq. 0 ) then
!               ia3 = j
!            end if
!            if( ia2 .gt. 0 .and. ia3 .gt. 0 .and.
!     &           j .ne. ia2 .and. j .ne. ia3 .and. ia4 .eq. 0 ) then
!               ia4 = j
!            end if
!         end if

!      end do

!***********************************************************************
!change coalescence parameter
!2010/11/8 sawada
!       hd  ::  coalescence parameter of deuteron
!       hth ::  coalescence parameter of triton or 3He
!       ha  ::  coalescence parameter of alpha
!      hd  = 1300.0d0
!      hd  = 1600.0d0
!      hth = 1800.0d0
!      ha  = 2000.0d0
!2012/10/22  nogamine
      hd  =2000.0d0
      hth =2250.0d0
!      ha  = 2500.0d0
      ha  =3500.0d0
!***********************************************************************
!     alpha

      h0 = ha
!---2022/10/12 yamaguchi
      h2 = 1.0d6 
      h3 = 1.0d6
      h4 = 1.0d6
!---      
      ia2 = 0
      ia3 = 0
      ia4 = 0

      do j = 1, massba
!         if( iclst(j) .ne. 0 )cycle
         if( numclst(j) .ne. 0 )cycle
         if( j .eq. i1 ) cycle
         if( p(5,j) .lt. 0.10d0 ) cycle
         if( inds(j) .ne. 1 )cycle
         if( iclst(j) .eq. -1 )cycle     !sawada

         rsq = ( r(1,i1) - r(1,j) )**2 + ( r(2,i1) - r(2,j) )**2
     &        + ( r(3,i1) - r(3,j) )**2
         rij = sqrt( rsq )

         psq = ( p(1,i1) - p(1,j) )**2 + ( p(2,i1) - p(2,j) )**2
     &        + ( p(3,i1) - p(3,j) )**2
         pij = sqrt( psq )

         h = rij * pij

!         if( h .lt. ha .and. ipot(j) .eq. 0 .and. lcoll(7) .ne. 0 ) then
         if( h .lt. h0 .and. ipot(j) .eq. 0 .and. lcoll(7) .ne. 0 ) then
!            if( ia2 .eq. 0 ) then
!               ia2 = j
!            end if
!            if( ia2 .gt. 0 .and. j .ne. ia2 .and. ia3 .eq. 0 ) then
!               ia3 = j
!            end if
!            if( ia2 .gt. 0 .and. ia3 .gt. 0 .and.
!     &           j .ne. ia2 .and. j .ne. ia3 .and. ia4 .eq. 0 ) then
!               ia4 = j
!            end if
!---2022/10/12 yamaguchi
           if( h .lt. h2)then
              h4  = h3
              h3  = h2
              h2  = h
              ia4 = ia3
              ia3 = ia2
              ia2 = j
           else if(h.lt.h3)then
              h4  = h3
              h3  = h
              ia4 = ia3
              ia3 = j
           else if(h.lt.h4)then
              h4  = h
              ia4 = j
           end if
!---2022/10/12 yamaguchi
         end if
      end do

!3He and triton
      if( ia4 .eq. 0 )then
         h0 = hth

!---2022/10/12 yamaguchi
      h2 = 1.0d6 
      h3 = 1.0d6
      h4 = 1.0d6
!---      
         
       ia2 = 0
       ia3 = 0
       ia4 = 0

       do j = 1, massba
!          if( iclst(j) .ne. 0 )cycle
          if( numclst(j) .ne. 0 )cycle
          if( j .eq. i1 ) cycle
          if( p(5,j) .lt. 0.10d0 ) cycle
          if( inds(j) .ne. 1 )cycle
          if( iclst(j) .eq. -1 )cycle     !sawada

          rsq = ( r(1,i1) - r(1,j) )**2 + ( r(2,i1) - r(2,j) )**2
     &         + ( r(3,i1) - r(3,j) )**2
          rij = sqrt( rsq )

          psq = ( p(1,i1) - p(1,j) )**2 + ( p(2,i1) - p(2,j) )**2
     &         + ( p(3,i1) - p(3,j) )**2
          pij = sqrt( psq )

          h = rij * pij

!          if( h .lt. hth .and. ipot(j) .eq. 0 .and. lcoll(7) .ne. 0)then
          if( h .lt. h0 .and. ipot(j) .eq. 0 .and. lcoll(7) .ne. 0)then
!             if( ia2 .eq. 0 ) then
!                ia2 = j
!             end if
!             if( ia2 .gt. 0 .and. j .ne. ia2 .and. ia3 .eq. 0 ) then
!                ia3 = j
!             end if
!             if( ia2 .gt. 0 .and. ia3 .gt. 0 .and.
!     &            j .ne. ia2 .and. j .ne. ia3 .and. ia4 .eq. 0 ) then
!                ia4 = j
!             end if
           if( h .lt. h2)then
              h4  = h3
              h3  = h2
              h2  = h
              ia4 = ia3
              ia3 = ia2
              ia2 = j
           else if(h.lt.h3)then
              h4  = h3
              h3  = h
              ia4 = ia3
              ia3 = j
           else if(h.lt.h4)then
              h4  = h
              ia4 = j
           end if
!---2022/10/12 yamaguchi
        end if
       end do
      end if

!deuteron
      if( ia4 .eq. 0 .and. ia3 .eq. 0 )then
         h0 = hd

!---2022/10/12 yamaguchi
      h2 = 1.0d6 
      h3 = 1.0d6
      h4 = 1.0d6
!---      
         
       ia2 = 0
       ia3 = 0
       ia4 = 0

       do j = 1, massba
!          if( iclst(j) .ne. 0 )cycle
          if( numclst(j) .ne. 0 )cycle
          if( j .eq. i1 ) cycle
          if( p(5,j) .lt. 0.10d0 ) cycle
          if( inds(j) .ne. 1 )cycle
          if( iclst(j) .eq. -1 )cycle     !sawada

          rsq = ( r(1,i1) - r(1,j) )**2 + ( r(2,i1) - r(2,j) )**2
     &         + ( r(3,i1) - r(3,j) )**2
          rij = sqrt( rsq )

          psq = ( p(1,i1) - p(1,j) )**2 + ( p(2,i1) - p(2,j) )**2
     &         + ( p(3,i1) - p(3,j) )**2
          pij = sqrt( psq )

          h = rij * pij

!          if( h .lt. hd .and. ipot(j) .eq. 0 .and. lcoll(7) .ne. 0 )then
          if( h .lt. h0 .and. ipot(j) .eq. 0 .and. lcoll(7) .ne. 0 )then
!             if( ia2 .eq. 0 ) then
!                ia2 = j
!             end if
!             if( ia2 .gt. 0 .and. j .ne. ia2 .and. ia3 .eq. 0 ) then
!                ia3 = j
!             end if
!             if( ia2 .gt. 0 .and. ia3 .gt. 0 .and.
!     &            j .ne. ia2 .and. j .ne. ia3 .and. ia4 .eq. 0 ) then
!                ia4 = j
!             end if
           if( h .lt. h2)then
              h4  = h3
              h3  = h2
              h2  = h
              ia4 = ia3
              ia3 = ia2
              ia2 = j
           else if(h.lt.h3)then
              h4  = h3
              h3  = h
              ia4 = ia3
              ia3 = j
           else if(h.lt.h4)then
              h4  = h
              ia4 = j
           end if
!---2022/10/12 yamaguchi
          end if
       end do
      end if

!-----------------------------------------------------------------------
!     Helium 4 production
!-----------------------------------------------------------------------

      if( ia2 .gt. 0 .and. ia3 .gt. 0 .and. ia4 .gt. 0 ) then 
 
         iclchg = ichg(i1) + ichg(ia2) + ichg(ia3) + ichg(ia4)

         if( iclchg .eq. 2 ) then

            do k = 1, 3
               pclst(k) = p(k,i1) + p(k,ia2) + p(k,ia3) + p(k,ia4)
               pcmclst(k,i1) = pclst(k)               
            end do
!     pc1 = sqrt( pclst(1)**2 + pclst(2)**2 + pclst(3)**2 )
            pc1 = pclst(1)**2 + pclst(2)**2 + pclst(3)**2            
            tc1 = -amass + sqrt(amass**2 + pc1)
!            tc0 = tc1 + 28.29d0
!            pc0 = sqrt( tc0**2 + 2.0d0*amass*tc0 )
!            do k = 1, 3
!!               p(k,i1) = pclst(k)/pc1*pc0
!               p(k,i1) = pclst(k)*0.25d0
!               p(k,ia2) = p(k,i1)
!               p(k,ia3) = p(k,i1)
!               p(k,ia4) = p(k,i1)               
!            end do

            !**************************************************2012/09/13
            !to preserve moment after reproduction
!            p4ia2 = p(4,ia2)
!            p4ia3 = p(4,ia3)
!            p4ia4 = p(4,ia4)
            !****************************************************nogamine

            pcmclst(5,i1) = amass            
!            p(5,i1)  = amass
!            p(5,ia2) = 0.0d0
!            p(5,ia3) = 0.0d0
!            p(5,ia4) = 0.0d0
!            p(4,ia2) = 0.0d0
!            p(4,ia3) = 0.0d0
!            p(4,ia4) = 0.0d0
!-----------------------------------------------------------------------

!            p(4,i1) = sqrt( amass**2
!     &           + p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )
            pcmclst(4,i1) = dsqrt( pcmclst(5,i1)**2 + pc1 )
            
!            iclst(i1)  = 4
!            iclst(ia2) = 4
!            iclst(ia3) = 4
!            iclst(ia4) = 4

            ininclst(i1)  = 4
            ininclst(ia2) = 4
            ininclst(ia3) = 4
            ininclst(ia4) = 4
            nmasclst(i1)  = 4

            numclst(i1) = i1
            numclst(ia2) = i1
            numclst(ia3) = i1
            numclst(ia4) = i1
            
            icoales(i1)  = 2
            icoales(ia2) = 2
            icoales(ia3) = 2
            icoales(ia4) = 2

!            if( ichg(ia2) .eq. 1 ) then


!               if( p(4,i1) - amass .gt. aucpot ) then
            if( pcmclst(4,i1) - amass .gt. aucpot ) then

                  ipot(ia2) = 1
                  ipot(ia3) = 1
                  ipot(ia4) = 1
                  ipot(i1)  = 1
                  inio(ia2) = 1
                  inio(ia3) = 1
                  inio(ia4) = 1
                  inio(i1)  = 1
               
!                  ipot(ia2) = 0
!                  ipot(ia3) = 0
!                  ipot(ia4) = 0
!                  ipot(i1)  = 1

               else

                  ipot(ia2) = 0
                  ipot(ia3) = 0
                  ipot(ia4) = 0
                  ipot(i1)  = 0

               end if

!            else
!
!               if( p(4,i1) - amass .gt. aucpot ) then
!
!                  ipot(ia2) = 0
!                  ipot(ia3) = 0
!                  ipot(ia4) = 0
!                  ipot(i1)  = 1
!
!               else
!
!                  ipot(ia2) = 0
!                  ipot(ia3) = 0
!                  ipot(ia4) = 0
!                  ipot(i1)  = 0
!
!               end if
!            end if

!-----------------------------------------------------------------------
!     reproduction from Helium 4
           if( ipot(i1) .eq. 0 )then
!-----------------------------------------------------------------------

!*     Helium 3 or Triton production from Helium 4                     *
               iclchg = ichg(i1) + ichg(ia2) + ichg(ia4)
               do k = 1, 3
                  pclst(k) = pclst(k) - p(k,ia3)
                  pcmclst(k,i1) = pclst(k)                  
               end do
!     pc1 = sqrt( pclst(1)**2 + pclst(2)**2 + pclst(3)**2 )
               pc1 = pclst(1)**2 + pclst(2)**2 + pclst(3)**2               
!               tc1 = -tmass + sqrt( tmass**2 + pc1**2)

               if( iclchg .eq. 2 )then
!     tc0 = tc1 + 7.7181d0 !bindimg e = 7.7181[MeV]3He
                  tc1 = -hmass + dsqrt( ( hmass )**2 + pc1 )                  
               elseif( iclchg .eq. 1)then
!     tc0 = tc1 + 8.4820d0 !bindimg e = 8.482[MeV]triton
                  tc1 = -tmass + dsqrt( ( tmass )**2 + pc1 )                   
               end if

!               pc0 = sqrt( tc0**2 + 2.0d0 * tmass * tc0 )
!               do k = 1, 3
!!     p(k,i1)  = pclst(k) / pc1 * pc0
!                  p(k,i1)  = pclst(k)/3.0d0                  
!                  p(k,ia2)  = p(k,i1)                  
!                  p(k,ia4)  = p(k,i1)                  
!               end do

               if( iclchg .eq. 2 )then
                  pcmclst(5,i1) = hmass
               elseif( iclchg .eq. 1 )then
                  pcmclst(5,i1) = tmass
               end if
               
!               p(5,i1)  = tmass
!               p(5,ia2) = 0.0d0
!               p(4,ia2) = 0.0d0
               p(5,ia3) = rmass
!               p(4,ia3) = p4ia3
!               p(5,ia4) = 0.0d0
!               p(4,ia4) = 0.0d0

!               p(4,i1) = dsqrt( tmass**2 
!     &            + p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )
               pcmclst(4,i1) = dsqrt( pcmclst(5,i1)**2 + pc1 )
               if( iclchg .eq. 1 .and. 
!     &              p(4,i1) - tmass .gt. tucpot ) then
     &              pcmclst(4,i1) - tmass .gt. tucpot ) then
               
                  icoales(i1)  = 2
                  icoales(ia2) = 2
                  icoales(ia3) = 0
                  icoales(ia4) = 2

!                  iclst(i1)  = 2
!                  iclst(ia2) = 2
!                  iclst(ia3) = 0
!                  iclst(ia4) = 2
                  ininclst(i1)  = 2
                  ininclst(ia2) = 2
                  ininclst(ia3) = 0
                  ininclst(ia4) = 2

                  nmasclst(i1)  = 3

                  numclst(i1) = i1
                  numclst(ia2) = i1
                  numclst(ia3) = 0
                  numclst(ia4) = i1
                  
                  ipot(i1)  = 1
                  ipot(ia2) = 1
                  ipot(ia3) = 0
                  ipot(ia4) = 1
                  inio(ia2) = 1
                  inio(ia4) = 1
                  inio(i1)  = 1

               else if( iclchg .eq. 2 .and.
!     &              p(4,i1) - hmass .gt. hucpot ) then
     &                 pcmclst(4,i1) - hmass .gt. hucpot ) then
                  
                  icoales(i1)  = 2
                  icoales(ia2) = 2
                  icoales(ia3) = 0
                  icoales(ia4) = 2

!                  iclst(i1)  = 3
!                  iclst(ia2) = 3
!                  iclst(ia3) = 0
!                  iclst(ia4) = 3

                  ininclst(i1)  = 3
                  ininclst(ia2) = 3
                  ininclst(ia3) = 0
                  ininclst(ia4) = 3

                  nmasclst(i1)  = 3

                  numclst(i1) = i1
                  numclst(ia2) = i1
                  numclst(ia3) = 0
                  numclst(ia4) = i1
                  
                  ipot(i1)  = 1
                  ipot(ia2) = 1
                  ipot(ia3) = 0
                  ipot(ia4) = 1
                  inio(ia2) = 1
                  inio(ia4) = 1
                  inio(i1)  = 1                  
                  
               else

!*     Deuteron production from Helium 4                               *

                  iclchg = ichg(i1) + ichg(ia4)
                  do k = 1, 3
                     pclst(k) = pclst(k) - p(k,ia2)
                     pcmclst(k,i1) = pclst(k)                     
                  end do
!     pc1 = sqrt( pclst(1)**2 + pclst(2)**2 + pclst(3)**2 )
                  pc1 = pclst(1)**2 + pclst(2)**2 + pclst(3)**2                  
                  tc1 = -ddmass + 
     &                   sqrt( ddmass**2 + pc1)
!                  tc0 = tc1 + 2.2245d0 !bindimg e = 2.2245[MeV]
!                  pc0 = sqrt( tc0**2 + 2.0d0 * ddmass * tc0 )
!                  do k = 1, 3
!!                     p(k,i1)  = pclst(k) / pc1 * pc0
!                     p(k,i1)  = pclst(k)*0.5d0
!                     p(k,ia4) = p(k,i1)
!                  end do
!     p(5,i1)  = ddmass
                  pcmclst(5,i1) = ddmass                  
                  p(5,ia2) = rmass
!                  p(4,ia2) = p4ia2
                  p(5,ia3) = rmass
!                  p(4,ia3) = p4ia3
!                  p(5,ia4) = 0.0d0
!                  p(4,ia4) = 0.0d0

!                  p(4,i1) = dsqrt( ddmass**2 
!     &               + p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )
                  pcmclst(4,i1) = dsqrt( pcmclst(5,i1)**2 + pc1 )
                  if( iclchg .eq.1 .and. 
!     &               p(4,i1) - ddmass .gt. ducpot ) then
     &                 pcmclst(4,i1) - ddmass .gt. ducpot ) then
                  
                    icoales(i1)  = 2
                    icoales(ia2) = 0
                    icoales(ia3) = 0
                    icoales(ia4) = 2

!                    iclst(i1)  = 1
!                    iclst(ia2) = 0
!                    iclst(ia3) = 0
!                    iclst(ia4) = 1
                    ininclst(i1)  = 1
                    ininclst(ia2) = 0
                    ininclst(ia3) = 0
                    ininclst(ia4) = 1

                    nmasclst(i1)  = 2

                    numclst(i1) = i1
                    numclst(ia2) = 0
                    numclst(ia3) = 0
                    numclst(ia4) = i1
                    
                    ipot(i1)  = 1
                    ipot(ia2) = 0
                    ipot(ia3) = 0
                    ipot(ia4) = 1
                    inio(ia4) = 1
                    inio(i1)  = 1                  
                    
                  else
!-----------------------------------------------------------------------
!     Proton or Neutron production from Helium 4
!-----------------------------------------------------------------------
                     do k = 1, 3
!                        p(k,i1) = pclst(k) - p(k,ia4)
                        pcmclst(k,i1) = 0.0d0                        
                     end do

                     pcmclst(5,i1) = 0.0d0                     
                     p(5,i1)  = rmass
                     p(5,ia2) = rmass
!                     p(4,ia2) = p4ia2
                     p(5,ia3) = rmass
!                     p(4,ia3) = p4ia3
                     p(5,ia4) = rmass
!                     p(4,ia4) = p4ia4

                     pcmclst(4,i1) = 0.0d0                     
                     p(4,i1) = dsqrt( p(5,i1)**2 
     &                 + p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )

                     icoales(i1)  = 0
                     icoales(ia2) = 0
                     icoales(ia3) = 0
                     icoales(ia4) = 0

!                     iclst(i1)    = 0
!                     iclst(ia2)   = 0
!                     iclst(ia3)   = 0
!                     iclst(ia4)   = 0
                     ininclst(i1)    = 0
                     ininclst(ia2)   = 0
                     ininclst(ia3)   = 0
                     ininclst(ia4)   = 0

                     nmasclst(i1)  = 0

                     numclst(i1) = 0
                     numclst(ia2) = 0
                     numclst(ia3) = 0
                     numclst(ia4) = 0
                     
                     ipot(i1)     =  0
                     ipot(ia2)    =  0
                     ipot(ia3)    =  0
                     ipot(ia4)    =  0

                     if( ichg(i1) .eq.1 .and. 
     &                  p(4,i1) - p(5,i1) .gt. ucpot ) then

                        ipot(i1)  = 1
                        ipot(ia2) = 0
                        ipot(ia3) = 0
                        ipot(ia4) = 0
                        inio(i1)  = 1
                        
                     else if( ichg(i1) .eq. 0 .and. 
     &                  p(4,i1) - p(5,i1) .gt. upot ) then

                        ipot(i1)  = 1
                        ipot(ia2) = 0
                        ipot(ia3) = 0
                        ipot(ia4) = 0
                        inio(i1)  = 1                        
                     end if

                  end if
               end if

!-----------------------------------------------------------------------
!     end reproduction from Helium 4
           end if
!-----------------------------------------------------------------------

         end if

!-----------------------------------------------------------------------
!     Helium 3 production
!-----------------------------------------------------------------------

      elseif( ia2 .gt. 0 .and. ia3 .gt. 0 .and. ia4 .eq. 0 ) then

         iclchg = ichg(i1) + ichg(ia2) + ichg(ia3)

         if( iclchg .eq. 2 ) then

            do k = 1, 3
               pclst(k) = p(k,i1) + p(k,ia2) + p(k,ia3)
               pcmclst(k,i1) = pclst(k)               
            end do
!     pc1 = sqrt( pclst(1)**2 + pclst(2)**2 + pclst(3)**2 )
            pc1 = pclst(1)**2 + pclst(2)**2 + pclst(3)**2            
            tc1 = -hmass + sqrt( hmass**2 + pc1)
!            tc0 = tc1 + 7.7181d0 !bindimg e = 7.7181[MeV]
!            pc0 = sqrt( tc0**2 + 2.0d0 * hmass * tc0 )
!            do k = 1, 3
!!     p(k,i1)  = pclst(k) / pc1 * pc0
!               p(k,i1)  = pclst(k)/3.0d0
!               p(k,ia2) = p(k,i1)
!               p(k,ia3) = p(k,i1)               
!            end do

            !**************************************************2012/09/13
!            p4ia2 = p(4,ia2)
!            p4ia3 = p(4,ia3)
            !****************************************************nogamine
            pcmclst(5,i1)  = hmass
!            p(5,i1) = hmass
!            p(5,ia2) = 0.0d0
!            p(5,ia3) = 0.0d0
!            p(4,ia2) = 0.0d0
!            p(4,ia3) = 0.0d0
!            p(4,i1) = sqrt( p(5,i1)**2
!     &           + p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )
            pcmclst(4,i1) = dsqrt( pcmclst(5,i1)**2 + pc1 )

!            iclst(i1)  = 3
!            iclst(ia2) = 3
!            iclst(ia3) = 3
            ininclst(i1)  = 3
            ininclst(ia2) = 3
            ininclst(ia3) = 3

            nmasclst(i1)  = 3

            numclst(i1) = i1
            numclst(ia2) = i1
            numclst(ia3) = i1
            
            icoales(i1)  = 2
            icoales(ia2) = 2
            icoales(ia3) = 2

!            if( ichg(ia2) .eq. 1 ) then

               if( p(4,i1) - hmass .gt. hucpot ) then

                  ipot(ia2) = 1
                  ipot(ia3) = 1
                  ipot(i1)  = 1
                  inio(ia2) = 1
                  inio(ia3) = 1
                  inio(i1)  = 1

               else

                  ipot(ia2) = 0
                  ipot(ia3) = 0
                  ipot(i1)  = 0

               end if

!            else
!
!               if( p(4,i1) - hmass .gt. hucpot ) then
!
!
!                  ipot(ia2) = 0
!                  ipot(ia3) = 0
!                  ipot(i1)  = 1
!
!               else
!
!                  ipot(ia2) = 0
!                  ipot(ia3) = 0
!                  ipot(i1)  = 0
!
!               end if
!            endif

!-----------------------------------------------------------------------
!     reproduction from Helium 3
           if( ipot(i1) .eq. 0 )then
!-----------------------------------------------------------------------
!*     Deuteron production from Helium 3                               *
               iclchg = ichg(i1) + ichg(ia3)
               do k = 1, 3
                  pclst(k) = pclst(k) -p(k,ia2)
                  pcmclst(k,i1) = pclst(k)
               end do

!               p(5,i1)  = ddmass
!               p(5,ia2) = rmass
!               p(4,ia2) = p4ia2
!               p(5,ia3) = 0.0d0
!               p(4,ia3) = 0.0d0

!     pc1 = sqrt( pclst(1)**2 + pclst(2)**2 + pclst(3)**2 )
               pc1 = pclst(1)**2 + pclst(2)**2 + pclst(3)**2               
               tc1 = -ddmass + 
     &                      sqrt( ddmass**2 + pc1)
!               tc0 = tc1 + 2.2245d0 !bindimg e = 2.2245[MeV]
!               pc0 = sqrt( tc0**2 + 2.0d0 * ddmass * tc0 )
!               do k = 1, 3
!!     p(k,i1)  = pclst(k) / pc1 * pc0
!                  p(k,i1)  = pclst(k)*0.5d0
!                  p(k,ia3) = p(k,i1)                  
!               end do

               pcmclst(5,i1)  = ddmass
               p(5,ia2) = rmass
               
!               p(4,i1) = dsqrt( ddmass**2 
!     &            + p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )

               pcmclst(4,i1) = dsqrt( pcmclst(5,i1)**2 + pc1 )              
               if( iclchg .eq. 1 .and.
!     &               p(4,i1) - ddmass .gt. ducpot ) then
     &              pcmclst(4,i1) - ddmass .gt. ducpot ) then
               
                 icoales(i1)  = 2
                 icoales(ia2) = 0
                 icoales(ia3) = 2

!                 iclst(i1)  = 1
!                 iclst(ia2) = 0
!                 iclst(ia3) = 1
                 ininclst(i1)  = 1
                 ininclst(ia2) = 0
                 ininclst(ia3) = 1

                 nmasclst(i1)  = 2

                 numclst(i1) = i1
                 numclst(ia2) = 0
                 numclst(ia3) = i1
                 
                 ipot(i1)  = 1
                 ipot(ia2) = 0
                 ipot(ia3) = 1
                 inio(ia3) = 1
                 inio(i1)  = 1
                 
               else
!-----------------------------------------------------------------------
!     Proton or Neutron production from Helium 3
!-----------------------------------------------------------------------

                  do k = 1, 3
!                     p(k,i1) = pclst(k) - p(k,ia3)
                     pcmclst(k,i1) = 0.0d0                     
                  end do

                  pcmclst(5,i1) = 0.0d0                  
                  p(5,i1)  = rmass
                  p(5,ia2) = rmass
!                  p(4,ia2) = p4ia2
                  p(5,ia3) = rmass
!                  p(4,ia3) = p4ia3

                  pcmclst(4,i1) = 0.0d0                  
                  p(4,i1) = dsqrt( p(5,i1)**2 
     &               + p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )
                  icoales(i1)  = 0
                  icoales(ia2) = 0
                  icoales(ia3) = 0

!                  iclst(i1)    = 0
!                  iclst(ia2)   = 0
!                  iclst(ia3)   = 0
                  ininclst(i1)    = 0
                  ininclst(ia2)   = 0
                  ininclst(ia3)   = 0

                  nmasclst(i1)  = 0

                  numclst(i1) = 0
                  numclst(ia2) = 0
                  numclst(ia3) = 0
                  
                  ipot(i1)     = 0
                  ipot(ia2)    = 0
                  ipot(ia3)    = 0

                  if( ichg(i1) .eq.1 .and. 
     &               p(4,i1) - p(5,i1) .gt. ucpot ) then

                     ipot(i1)  = 1
                     ipot(ia2) = 0
                     ipot(ia3) = 0
                     inio(i1)  = 1
                     
                  else if( ichg(i1) .eq. 0 .and. 
     &               p(4,i1) - p(5,i1) .gt. upot ) then

                     ipot(i1)  = 1
                     ipot(ia2) = 0
                     ipot(ia3) = 0
                     inio(i1)  = 1
                     
                  end if
               end if

!-----------------------------------------------------------------------
!     end reproduction from Helium 3
           end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     triton production
!-----------------------------------------------------------------------

         elseif( iclchg .eq. 1 ) then

            do k = 1, 3
               pclst(k) = p(k,i1) + p(k,ia2) + p(k,ia3)
               pcmclst(k,i1) = pclst(k)               
            end do
!     pc1 = sqrt( pclst(1)**2 + pclst(2)**2 + pclst(3)**2 )
            pc1 = pclst(1)**2 + pclst(2)**2 + pclst(3)**2            
            tc1 = -tmass + sqrt(tmass**2 + pc1)
!            tc0 = tc1 + 8.4820d0
!            pc0 = sqrt( tc0**2 + 2.0d0*tmass*tc0 )
!            do k = 1, 3
!!     p(k,i1) = pclst(k)/pc1*pc0
!               p(k,i1)  = pclst(k)/3.0d0
!               p(k,ia2) = p(k,i1)
!               p(k,ia3) = p(k,i1)               
!            end do

            !**************************************************2012/09/13
!            p4ia2 = p(4,ia2)
!            p4ia3 = p(4,ia3)
            !****************************************************nogamine

            pcmclst(5,i1)  = tmass            
!            p(5,i1) = tmass
!            p(5,ia2) = 0.0d0
!            p(5,ia3) = 0.0d0
!            p(4,ia2) = 0.0d0
!            p(4,ia3) = 0.0d0

!            p(4,i1) = sqrt( tmass**2
!     &           + p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )
            pcmclst(4,i1) = dsqrt( pcmclst(5,i1)**2 + pc1 )
            
!            iclst(i1)  = 2
!            iclst(ia2) = 2
!            iclst(ia3) = 2
            ininclst(i1)  = 2
            ininclst(ia2) = 2
            ininclst(ia3) = 2

            nmasclst(i1)  = 3

            numclst(i1) = i1
            numclst(ia2) = i1
            numclst(ia3) = i1
            
            icoales(i1)  = 2
            icoales(ia2) = 2
            icoales(ia3) = 2

!            if( ichg(ia2) .eq. 1 ) then
!               if( p(4,i1) - tmass .gt. tucpot ) then
               if( pcmclst(4,i1) - tmass .gt. tucpot ) then
                  ipot(ia2) = 1
                  ipot(ia3) = 1
                  ipot(i1)  = 1
                  inio(ia2) = 1
                  inio(ia3) = 1
                  inio(i1)  = 1                  
               else
                  ipot(ia2) = 0
                  ipot(ia3) = 0
                  ipot(i1)  = 0
               end if
!            else
!               if( p(4,i1) - tmass .gt. tucpot ) then
!                  ipot(ia2) = 0
!                  ipot(ia3) = 0
!                  ipot(i1)  = 1
!               else
!                  ipot(ia2) = 0
!                  ipot(ia3) = 0
!                  ipot(i1)  = 0
!               end if
!            endif

!-----------------------------------------------------------------------
!     reproduction from triton
           if( ipot(i1) .eq. 0 )then
!-----------------------------------------------------------------------
!*     Deuteron production from triton                                 *
               iclchg = ichg(i1) + ichg(ia2)
               do k = 1, 3
                  pclst(k) = pclst(k) - p(k,ia3)
                  pcmclst(k,i1) = pclst(k)                  
               end do

!               p(5,i1)  = ddmass
!               p(5,ia2) = 0.0d0
!               p(4,ia2) = 0.0d0
!               p(5,ia3) = rmass
!               p(4,ia3) = p4ia3

!     pc1 = sqrt( pclst(1)**2 + pclst(2)**2 + pclst(3)**2 )
               pc1 = pclst(1)**2 + pclst(2)**2 + pclst(3)**2               
               tc1 = -ddmass + 
     &                      sqrt( ddmass**2 + pc1)
!               tc0 = tc1 + 2.2245d0 !bindimg e = 2.2245[MeV]
!               pc0 = sqrt( tc0**2 + 2.0d0 * ddmass * tc0 )
!               do k = 1, 3
!!     p(k,i1)  = pclst(k) / pc1 * pc0
!                  p(k,i1)  = pclst(k)*0.5d0
!                  p(k,ia2) = p(k,i1)                  
!               end do

               pcmclst(5,i1)  = ddmass               
               p(5,ia3) = rmass
!               p(4,i1) = dsqrt( ddmass**2 
!     &            + p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )
               pcmclst(4,i1) = dsqrt( pcmclst(5,i1)**2 + pc1 )               
               if( iclchg .eq. 1 .and.
!     &               p(4,i1) - ddmass .gt. ducpot ) then
     &               pcmclst(4,i1) - ddmass .gt. ducpot ) then               
                 icoales(i1)  = 2
                 icoales(ia2) = 2
                 icoales(ia3) = 0

!                 iclst(i1)  = 1
!                 iclst(ia2) = 1
!                 iclst(ia3) = 0
                 ininclst(i1)  = 1
                 ininclst(ia2) = 1
                 ininclst(ia3) = 0

                 nmasclst(i1)  = 2

                 numclst(i1) = i1
                 numclst(ia2) = i1
                 numclst(ia3) = 0
                 
                 ipot(i1)  = 1
                 ipot(ia2) = 1
                 ipot(ia3) = 0
                 inio(ia2) = 1
                 inio(i1)  = 1                 
               else
!-----------------------------------------------------------------------
!     Proton or Neutron production from triton
!-----------------------------------------------------------------------

                  do k = 1, 3
!                     p(k,i1) = pclst(k) - p(k,ia2)
                     pcmclst(k,i1) = 0.0d0                     
                  end do

                  pcmclst(5,i1) = 0.0d0                  
                  p(5,i1)  = rmass
                  p(5,ia2) = rmass
!                  p(4,ia2) = p4ia2
                  p(5,ia3) = rmass
!                  p(4,ia3) = p4ia3

                  pcmclst(4,i1) = 0.0d0                  
                  p(4,i1) = dsqrt( p(5,i1)**2 
     &               + p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )

                  icoales(i1)  = 0
                  icoales(ia2) = 0
                  icoales(ia3) = 0

!                  iclst(i1)    = 0
!                  iclst(ia2)   = 0
!                  iclst(ia3)   = 0
                  ininclst(i1)    = 0
                  ininclst(ia2)   = 0
                  ininclst(ia3)   = 0

                  nmasclst(i1)  = 0

                  numclst(i1) = 0
                  numclst(ia2) = 0
                  numclst(ia3) = 0
                  
                  ipot(i1)     = 0
                  ipot(ia2)    = 0
                  ipot(ia3)    = 0

                  if( ichg(i1) .eq.1 .and. 
     &               p(4,i1) - p(5,i1) .gt. ucpot ) then

                     ipot(i1)  = 1
                     ipot(ia2) = 0
                     ipot(ia3) = 0
                     inio(i1)  = 1
                     
                  else if( ichg(i1) .eq. 0 .and. 
     &               p(4,i1) - p(5,i1) .gt. upot ) then

                     ipot(i1)  = 1
                     ipot(ia2) = 0
                     ipot(ia3) = 0
                     inio(i1)  = 1                     
                  end if
               end if

!-----------------------------------------------------------------------
!     end reproduction from triton
           end if
!-----------------------------------------------------------------------

         end if

!-----------------------------------------------------------------------
!     deuteron production
!-----------------------------------------------------------------------

      elseif( ia2 .gt. 0 .and. ia3 .eq. 0 .and. ia4 .eq. 0 ) then

         iclchg = ichg(i1) + ichg(ia2) !not given (S.H. on 2012.10.13)

         if( iclchg .eq. 1 ) then

            do k = 1, 3
               pclst(k) = p(k,i1) + p(k,ia2)
               pcmclst(k,i1) = pclst(k)                
            end do
!     pc1 = sqrt( pclst(1)**2 + pclst(2)**2 + pclst(3)**2 )
            pc1 = pclst(1)**2 + pclst(2)**2 + pclst(3)**2            
            tc1 = -ddmass + sqrt( ddmass**2 + pc1)
!            tc0 = tc1 + 2.2245d0 !bindimg e = 2.2245[MeV]
!            pc0 = sqrt( tc0**2 + 2.0d0 * ddmass * tc0 )
!            do k = 1, 3
!!     p(k,i1)  = pclst(k) / pc1 * pc0
!               p(k,i1)  = pclst(k)*0.5d0
!               p(k,ia2) = p(k,i1)               
!            end do
            !**************************************************2012/09/13
!            p4ia2 = p(4,ia2)
            !****************************************************nogamine
            pcmclst(5,i1)  = ddmass
!            p(5,i1) = ddmass
!            p(5,ia2) = 0.0d0
!            p(4,ia2) = 0.0d0

!            p(4,i1) = sqrt( ddmass**2
!     &           + p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )
            pcmclst(4,i1) = dsqrt( pcmclst(5,i1)**2 + pc1 )

!            iclst(i1)  = 1
!            iclst(ia2) = 1
            ininclst(i1)  = 1
            ininclst(ia2) = 1

            nmasclst(i1)  = 2

            numclst(i1) = i1
            numclst(ia2) = i1
            
            icoales(i1)  = 2
            icoales(ia2) = 2

!            if( ichg(ia2) .eq. 1 ) then
               if( pcmclst(4,i1) - ddmass .gt. ducpot ) then
                  ipot(ia2) = 1
                  ipot(i1)  = 1
                  inio(ia2) = 1
                  inio(i1)  = 1
               else
                  ipot(ia2) = 0
                  ipot(i1)  = 0
               end if
!            else
!               if( p(4,i1) - ddmass .gt. ducpot ) then
!                  ipot(ia2) = 0
!                  ipot(i1)  = 1
!               else
!                  ipot(ia2) = 0
!                  ipot(i1)  = 0
!               end if
!            endif

!-----------------------------------------------------------------------
!     reproduction from deuteron
           if( ipot(i1) .eq. 0 )then
!-----------------------------------------------------------------------
!     Proton or Neutron production from deuteron
!-----------------------------------------------------------------------

               do k = 1, 3
!                  p(k,i1) = pclst(k) - p(k,ia2)
                  pcmclst(k,i1) = 0.0d0                  
               end do

               pcmclst(5,i1) = 0.0d0               
               p(5,i1)  = rmass
               p(5,ia2) = rmass
!               p(4,ia2) = p4ia2

               pcmclst(4,i1) = 0.0d0               
               p(4,i1) = dsqrt( p(5,i1)**2 
     &            + p(1,i1)**2 + p(2,i1)**2 + p(3,i1)**2 )

               p(5,i1)  = rmass
               p(5,ia2) = rmass

               icoales(i1)  = 0
               icoales(ia2) = 0

!               iclst(i1)    = 0
!               iclst(ia2)   = 0
               ininclst(i1)    = 0
               ininclst(ia2)   = 0

               nmasclst(i1)  = 0

               numclst(i1) = 0
               numclst(ia2) = 0
               
               ipot(i1)     = 0
               ipot(ia2)    = 0

               if( ichg(i1) .eq.1 .and. 
     &           p(4,i1) - p(5,i1) .gt. ucpot ) then

                  ipot(i1)  = 1
                  ipot(ia2) = 0
                  inio(i1)  = 1
                  
               else if( ichg(i1) .eq. 0 .and. 
     &            p(4,i1) - p(5,i1) .gt. upot ) then

                  ipot(i1)  = 1
                  ipot(ia2) = 0
                  inio(i1)  = 1                  
               end if

!-----------------------------------------------------------------------
!     end reproduction from deuteron
           end if
!-----------------------------------------------------------------------

         end if

      end if

!      return
      END SUBROUTINE coales

!***********************************************************************
!*                                                                     *
!*                          ground State                               *
!*                                                                     *
!*                                                                     *
!*                    Last revised : Feb/14/2008                       *
!*                                                                     *
!*                                                                     *
!*  List of subprograms in rough order of relevance with main purpose  *
!*     ( s = SUBROUTINE, f = FUNCTION, b = block data, e = entry )     *
!*                                                                     *
!*                                                                     *
!*  s  incgrnd    to make ground state of the target nucleus           *
!*  s  packinc   to make ground state by random packing method         *
!*  f  eliqq      to calculate liquid drop binding energy (MeV)        *
!*  f  ws        to calculate the Woods-Saxon potential                *
!*  f  dws       to calculate derivative of the Woods-Saxon potential  *
!*  f  trand     random number generator                               *
!*                                                                     *
!*                                                                     *
!***********************************************************************
!*                                                                     *
!***********************************************************************
!*                                                                     *
!      SUBROUTINE incgrnd(iproj,nmasta,nchgta,nmaspr,mstq1,ein,upot,r,p,
!     &     ipot,inio,ichg,bimp,rt00,bmax,tfm,ll,radm,pfm,tt0,ihis,iavd)
      SUBROUTINE incgrnd(iproj,nmasta,nchgta,nmaspr,mstq1,ein,upot,r,p,
     &     ipot,inio,ichg,bimp,rt00,bmax,tfm,ll,radm,pfm,tt0,ihis,iavd,
     &     inuc,ibry,inds,inun,massal,massba,iclst,ngtgr )
!*                                                                     *
!*                                                                     *
!*       Purpose:                                                      *
!*            make the ground state of target and projectile           *
!*            for INC mode using the random packing method             *
!*                                                                     *
!*                                                                     *
!*       Variables:                                                    *
!*          [in]                                                       *
!*             iproj       : projectile information                    *
!*             nmasta      : target mass                               *
!*             nchgta      : proton mass in the target nucleus         *
!*             nmaspr      : projectile mass                           *
!*             mstq1       : input information                         *
!*             ein         : incident energy (MeV)                     *
!*             upot        : potential depth (MeV)                     *
!*                                                                     *
!*          [out]                                                      *
!*             r           : positions of nucleus                      *
!*             p           : momenta of nucleus                        *
!*             ipot        : above/under the Fermi surface             *
!*             inio        : inside/outside the nucleus                *
!*             ichg        : isospin of nucleon                        *
!*             bimp        : impact parameter(fm)                      *
!*             rt00        : Woods-Saxon radius(fm)                    *
!*             bmax        : Maximum of impact parameter (fm)          *
!*             tfm         : Fermi energy (MeV)                        *
!*             radm        : cutoff radius (fm)                        *
!*             pfm         : Fermi momentum (MeV/c)                    *
!*             tt0         : Total kinetic energy of ground state      *
!*                                                                     *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      integer   ipot(nnn),inio(nnn),ichg(nnn),ihis(nnn),iavd(nnn)
      integer   nmasta,nchgta,nmaspr,mstq1(10),ll,iproj,iclst(nnn)
      integer   ngtgr
!*************************************************************sawada8/27
      integer   inuc(nnn),ibry(nnn),inun(nnn),massal,massba
!***********************************************************************

      real*8    ein,upot
      real*8    r(3,nnn),p(6,nnn),bimp,bmax,tfm,radm,pfm,tt0,rt00,saa
      dimension inds(nnn)


!-----------------------------------------------------------------------
!     generate seeds
!-----------------------------------------------------------------------

      CALL genseed(iseed)

!-------------------------------------------------------------sawada8/27
      massba = nmasta +  nmaspr
      massal = massba
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Nucleon or Pion
!-----------------------------------------------------------------------

      idntp = 0  ! Now nucleon only

      if( idntp .eq. 0 ) then
         inutp = 1
         ibrtp = 1
         indtp = 1
      elseif( idntp .ne. 0 ) then
         inutp = 0
         ibrtp = 0
         indtp = 4
      end if

!-----------------------------------------------------------------------

      ipid = mstq1(1)           ! Momentum distribution
      irid = mstq1(2)           ! Spatial distribution

!-----------------------------------------------------------------------
!     set coordinate of each particle
!-----------------------------------------------------------------------

      do i = 1, massal

!-----------------------------------------------------------------------
!       set variables
!          indtp: nucleon/pion              [ 4: pion     1: nucleon ]
!          ipot: above/under potential      [ 1: above,   0: under   ]
!          inio: inside/outside the nucleus [ 1: outside, 0: inside  ]
!          ichg: isospin of nucleon         [ 1: up(p)    0: down(n) ]
!-----------------------------------------------------------------------

         if( i .eq. 1 ) iavd(i) = -1  !projectile
         if( i .ne. 1 ) iavd(i) =  1  !target

         ihis(i) = 0
         inio(i) = 0

!*************************************************************sawada8/27
         inuc(i) = inutp
         ibry(i) = ibrtp
         inun(i) = 0
!***********************************************************************

         inds(i) = indtp

         do j =1, 3
            r(j,i)  = 0.0d0
            p(j,i)  = 0.0d0
         end do

         p(4,i)  = 0.0d0
         p(5,i)  = rmass

         if( i .eq. 1 ) then
            ipot(i) = 1
            if( iproj .eq. 1 ) then!proton
               ichg(i) = 1
            elseif( iproj .eq. 0 ) then!neuteron
               ichg(i) = 0
            elseif( iproj .eq. 2 ) then!pi+
               ichg(i) = 1
               inds(i) = 4
               p(5,i)  = pmass*1000.0d0
            elseif( iproj .eq. 3 ) then!pi0
               ichg(i) = 0
               inds(i) = 4
               p(5,i)  = pmass*1000.0d0
            elseif( iproj .eq. 4 ) then!pi-
               ichg(i) = -1
               inds(i) = 4
               p(5,i)  = pmass*1000.0d0
            elseif( iproj .eq. 5) then !alpha
               ichg(i) = 2
               p(5,i) = amass
               iclst(i) = 4
            end if
         else
            ipot(i) = 0
            if( i .le. nchgta + 1 ) then
               ichg(i) = 1
            else
               ichg(i) = 0
            end if
         end if

      end do

!-----------------------------------------------------------------------
!     set positions and momenta of nucleons
!-----------------------------------------------------------------------

      CALL packinc(massal,nmasta,nmaspr,nchgta,ichg,inds,ein,upot,
     &     r,p,ipid,irid,bimp,rt00,bmax,tfm,ll,radm,pfm,tt0,iproj,ngtgr)

!-----------------------------------------------------------------------
      END SUBROUTINE incgrnd



!***********************************************************************
!*                                                                     *
!      SUBROUTINE packinc(massal,nmasta,nchgta,ichg,ein,upot,r,p,ipid,
!     &     irid,bimp,rt00,bmax,tfm,ll,radm,pfm,tt0)
      SUBROUTINE packinc(massal,nmasta,nmaspr,nchgta,ichg,inds,ein,upot,
     &     r,p,ipid,irid,bimp,rt00,bmax,tfm,ll,radm,pfm,tt0,iproj,ngtgr)
!*                                                                     *
!*                                                                     *
!*        Purpose:                                                     *
!*                                                                     *
!*              to sample positions and momenta of nucleus             *
!*              according to the Woods-Saxon type distribution.        *
!*                                                                     *
!*                                                                     *
!*        Variables:                                                   *
!*                                                                     *
!*           [in]                                                      *
!*              massal      : total mass                               *
!*              nmasta      : target mass                              *
!*              nchgta      : proton mass in the target nucleus        *
!*              ichg        : isospin of nucleon                       *
!*              ein         : incident energy (MeV)                    *
!*              upot        : potential depth (MeV)                    *
!*              r           : positions of nucleus                     *
!*              p           : momenta   of nucleus                     *
!*              ipid        : momentum distribution                    *
!*              irid        : spartial distribution                    *
!*                                                                     *
!*           [out]                                                     *
!*              r           : positions of nucleus                     *
!*              p           : momenta   of nucleus                     *
!*              bimp        : impact parameter(fm)                     *
!*              rt00        : Woods-Saxon radius(fm)                   *
!*              bmax        : maximum impact parameter (fm)            *
!*              tfm         : Fermi energy (MeV)                       *
!*              radm        : cutoff radius (fm)                       *
!*              pfm         : Fermi momentum (MeV/c)                   *
!*                                                                     *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      integer  massal,nmaspr,nmasta,nchgta,ichg(nnn),inds(nnn)
     &         ,ipid,irid,ll,ngtgr
      real*8   ein,upot
      real*8   r(3,nnn),p(6,nnn)
      real*8   bimp,rt00,bmax,tfm,radm,tt0,saa,pfm

!-----------------------------------------------------------------------

      dimension  cmangp(nnn)

!-----------------------------------------------------------------------
!     generate seeds
!-----------------------------------------------------------------------

      CALL genseed(iseed)

!-----------------------------------------------------------------------
!     determine position of nucleons
!               Woods-Saxon radius (rt00)
!               cut-off radius     (radm)
!               maximum value      (rrmax)
!-----------------------------------------------------------------------

      irid = 1
      if( irid .eq. 1 ) then

      !---- Woods-Saxon parameter( Ref. Phys.Rev.C vol.1, 4 (1970) )----
         saa   = 0.54d0         ![fm] ( default: 0.54  )
         ddif2 = 1.00d0         ![fm] ( default: 1.00  )
         dsam2 = 2.25d0         ![fm] ( default: 2.25  )
      !-----------------------------------------------------------------

            rt00 = ( 0.978 + 0.0206*dble( nmasta )**(1./3.) )
     &           * dble( nmasta )**(1./3.)
!-----------------------------------------------------------------------

!         if( nmasta .le. 12 ) then
!            radm = rt00 + 8.0 * saa
!         else
            radm = rt00 + 5.0 * saa
!         end if

         rrmax = 1.0 / ( 1.0 + exp( - rt00 / saa ) )

         do i = nmaspr+1, massal
            rwod = -1.0
            do while( rn(0) * rrmax .gt. rwod )
               rsqr = 10.0
               do while( rsqr .gt. 1.0 )
                  rx = 1.0 - 2.0 * rn(0)
                  ry = 1.0 - 2.0 * rn(0)
                  rz = 1.0 - 2.0 * rn(0)

                  rsqr = rx*rx + ry*ry + rz*rz
               end do
               rrr  = radm * sqrt(rsqr)
               rwod = 1.0
     &              / ( 1.0 + exp( ( rrr - rt00 ) / saa ) )
            end do
            r(1,i) = radm * rx
            r(2,i) = radm * ry
            r(3,i) = radm * rz
         end do

      end if

!-----------------------------------------------------------------------
!     Determine momentum of nucelons stochastically
!            Fermi momentum is determined by the square well potential
!            and binding energy calculated by the liquid drop model.
!            Momentum distribution follows the uniform Fermi gas model.
!-----------------------------------------------------------------------

      nz  = nchgta
      nn  = nmasta - nchgta

      tt0  = 0.0
      tfm = upot - eliqq(nz,nn) / dble(nmasta)
      efm = tfm + rmass
      pfm = sqrt( efm**2 - rmass**2 ) ! [MeV/c]
 9989 continue

!-----------------------------------------------------------------------
!   case 1 :  Uniform Fermi Gas distribution
!-----------------------------------------------------------------------
cKN 2016/07/21 for test

      ipid = 1
c     if( rn(0) .gt. 0.5 ) ipid = 1


      if( ipid .eq. 1 ) then

         do i = nmaspr+1, massal
            psqr = 10.0
            do while( psqr .gt. 1.0 )
               px = 1.0 - 2.0 * rn(0)
               py = 1.0 - 2.0 * rn(0)
               pz = 1.0 - 2.0 * rn(0)
               psqr = px*px + py*py + pz*pz
            end do
            p(1,i) = pfm* px
            p(2,i) = pfm* py
            p(3,i) = pfm* pz
            p(4,i) = sqrt( psqr * pfm**2 + p(5,i)**2 )
            tt0 = tt0 - rmass
     &           + sqrt( rmass**2 + p(1,i)**2 + p(2,i)**2
     &           + p(3,i)**2 )
         end do

!-----------------------------------------------------------------------
!     case 2 :  Fermi-Dirac distribution
!               Fermi radius       (pt00)
!               cut-off radius     (cutp)
!               Maxmum value       (ppmax)
!-----------------------------------------------------------------------

      elseif( ipid .eq. 0 ) then

         pt00 = pfm/2.0d0
         cutp = pfm
         spp  = 50.0
         ppmax = 1.0 / ( 1.0 + exp( - pt00 / spp ) )
         do i = nmaspr+1, massal
            pfrm = -1.0
            do while( rn(0) * ppmax .gt. pfrm )
               psqr = 10.0
               do while( psqr .gt. 1.0 )
                  px = 1.0 - 2.0 * rn(0)
                  py = 1.0 - 2.0 * rn(0)
                  pz = 1.0 - 2.0 * rn(0)
                  psqr = px*px + py*py + pz*pz
               end do
               ppp  = cutp * sqrt(psqr)
               pfrm = 1.0
     &              / ( 1.0 + exp( ( ppp - pt00 ) / spp ) )
            end do
            p(1,i) = cutp * px
            p(2,i) = cutp * py
            p(3,i) = cutp * pz
            p(4,i) = sqrt( psqr * cutp**2 + p(5,i)**2 )
            tt0 = tt0 + ( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )
     &           / 2.0 / rmass
         end do

      end if

!-----------------------------------------------------------------------
!     initial incident energy and momentum control
!-----------------------------------------------------------------------
      em1    = p(5,1)
      pin    = ein + em1

      if( inds(1) .eq. 1 )then
       p(4,1) = ( ein + upot*(em1/rmass) ) + em1 ! without refraction
       if( iproj .eq. 5 ) then
          p(4,1) = ( ein ) + em1
       end if
      elseif( inds(1) .eq. 4 )then
       p(4,1) = ein + em1            ! without refraction
      end if

      e1     = p(4,1)
      p(1,1) = 0.0
      p(2,1) = 0.0
      p(3,1) = sqrt( e1**2 - em1**2 )
!-----------------------------------------------------------------------
!     impact parameter control
!-----------------------------------------------------------------------

      rrand = 10.0
      do while( rrand .gt. 1.0 )
         xrand = 1.0 - 2.0 * rn(0)
         yrand = 1.0 - 2.0 * rn(0)
         rrand = xrand * xrand + yrand * yrand
      end do
      bmax = radm
      if( nmaspr .eq. 4) then
        bmax = radm + 3.2
      end if

      r(1,1) = xrand * bmax
      r(2,1) = yrand * bmax
!      zmax = radm + 1.0
      zmax = bmax
      rrand = dsqrt( rrand ) !watanabe8/27
      bimp = rrand * bmax
      alpha = dasin( bimp / zmax )
      r(3,1) = - zmax * dcos( alpha )

      if( bimp .ge. bmax * 0.90) then
        ngtgr = 1
      else
        ngtgr = 0
      end if


      if( iproj .eq. 5) then
              r(3,1) = r(3,1) -3.40d0
      end if
 
!-----------------------------------------------------------------------

      END SUBROUTINE



!***********************************************************************
!*                                                                     *
      FUNCTION eliqq(nz,nn)
!*                                                                     *
!*                                                                     *
!*        Purpose:                                                     *
!*              to calculate liquid drop binding energy ( MeV )        *
!*                                                                     *
!*                                                                     *
!*        Variables:                                                   *
!*          [in]                                                       *
!*              nz    :   proton number                                *
!*              nn    :   neutron number                               *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      parameter( bvol=15.56, bsur=17.23, bsym=46.57, coul=1.2 )
      
      eliqq = 0.0

      if( nn + nz .le. 1 ) return

      mm = nn + nz
      a  = 1.0 * nn + nz
      rc = 1.24 * a**(1./3.)

      eliqq = ( bvol * a
     &     - bsur * a**(2./3.)
     &     - bsym / 2. * ( nn - nz ) * ( nn - nz ) / a
     &     - 3./5. * nz * coul * nz * coul / rc )

      if( mm.eq. 2 .and. nz.eq. 1 ) eliqq = 1.11 * a
      if( mm.eq. 3 .and. nz.eq. 1 ) eliqq = 2.83 * a
      if( mm.eq. 3 .and. nz.eq. 2 ) eliqq = 2.57 * a
      if( mm.eq. 4 .and. nz.eq. 1 ) eliqq = 1.40 * a
      if( mm.eq. 4 .and. nz.eq. 2 ) eliqq = 7.07 * a
      if( mm.eq. 4 .and. nz.eq. 3 ) eliqq = 1.20 * a
      if( mm.eq. 5 .and. nz.eq. 2 ) eliqq = 5.48 * a
      if( mm.eq. 5 .and. nz.eq. 3 ) eliqq = 5.27 * a
      if( mm.eq. 6 .and. nz.eq. 4 ) eliqq = 4.49 * a
      if( mm.eq. 6 .and. nz.eq. 2 ) eliqq = 4.88 * a
      if( mm.eq. 6 .and. nz.eq. 3 ) eliqq = 5.33 * a
      if( mm.eq. 8 .and. nz.eq. 4 ) eliqq = 7.06 * a
      if( mm.eq.11 .and. nz.eq. 3 ) eliqq = 4.14 * a
      if( mm.eq.12 .and. nz.eq. 6 ) eliqq = 7.68 * a
      if( mm.eq.16 .and. nz.eq. 8 ) eliqq = 7.98 * a

      END FUNCTION eliqq



!***********************************************************************
!*                                                                     *
!*                          Output on file                             *
!*                                                                     *
!*                                                                     *
!*                    Last revised : Feb/15/2008                       *
!*                                                                     *
!*                                                                     *
!*  List of subprograms in rough order of relevance with main purpose  *
!*     ( s = SUBROUTINE, f = FUNCTION, b = block data, e = entry )     *
!*                                                                     *
!*                                                                     *
!*  s  sm_sumo   to write final output on file                         *
!*                                                                     *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!*                                                                     *
!      SUBROUTINE sm_sumo(ein,iproj,nmasta,nchgta,nmaspr,massal,mstq1,
!     &   rt00,tfm,ipot,ichg,inds,r,p,upot,lcoll,ncoll,tt0,iclst,icoales,
!     &     logemit,nemit,inf)

!      SUBROUTINE sm_sumo(ein,iproj,nmasta,nchgta,massal,
!     &   rt00,ipot,ichg,inds,p,upot,iclst,nemit,inf,eex)!eex for debug   fuku
      SUBROUTINE sm_sumo(ein,iproj,nmasta,nchgta,massal,
     &   rt00,ipot,ichg,inds,p,upot,iclst,nemit,inf,eex,peex)!eex for debug   fuku
!*                                                                     *
!*                                                                     *
!*       Purpose:                                                      *
!*                                                                     *
!*             to write final output on file                           *
!*                                                                     *
!*                                                                     *
!*        Variables:                                                   *
!*                                                                     *
!*           [in]                                                      *
!*              ein         : incident energy                          *
!*              nmasta      : target mass                              *
!*              nchgta      : target charge                            *
!*              nmaspr      : projectile mass                          *
!*              mstq1       : input information                        *
!*              tfm         : Fermi energy                             *
!*              ipot        :
!*                                                                     *
!*           [out]                                                     *
!*              p           : momenta of nucleus                       *
!*              pcm(3)      : momentum coordinates of one particle     *
!*                            in cm frame                              *
!*              inds        : nucleon/pion                             *
!*              ipot        : above/under potential depth              *
!*              ichannel    : channel information back                 *
!*                                                                     *
!*                                                                     *
!*        Output Variables:                                            *
!*                                                                     *
!*              iz          : proton number of mother                  *
!*              in          : neutron number of mother                 *
!*              p(i,j)      : momentum vector of mother                *
!*              et          : energy of mother                         *
!*              rm          : rest mass of mother                      *
!*                                                                     *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'
!-----------------------------------------------------------------------
      integer, intent(in)    :: iproj,nmasta,nchgta,massal,ipot(nnn)
     &                         ,inds(nnn),iclst(nnn)
      integer, intent(out)   :: nemit
      integer, intent(inout) :: ichg(nnn)
      real(8), intent(in)    :: ein,rt00,upot,peex
      real(8), intent(out)   :: inf(5,nnn)
      real(8), intent(inout) :: p(6,nnn)

      real(8),  intent(out)   :: eex !for debug fuku
!-----------------------------------------------------------------------

!      integer  nmasta,nchgta,nmaspr,iproj,ipot(nnn),lcoll(30)
!      integer  ncoll(5,nnn),mstq1(10)
!      integer  ichg(nnn),iclst(nnn),icoales(nnn),nemit
      integer  ineu(nnn)
!     real*8   upot,ein,r(3,nnn),tfm,tt0,p(6,nnn),inf(5,nnn),coulom(6)
      real(8) coulom(6)      
      real(8) resp(3),pin(3) !fuku for debug
!      character(len = 100 ) logemit(nnn) 

      parameter( eps = 1.0d-10 )

!      integer inds(nnn)
!      real*8  rt00
!-----------------------------------------------------------------------
      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )
!-----------------------------------------------------------------------
      ducpot = dupot  + coul(nchgta,rt00)
      tucpot = thupot + coul(nchgta,rt00)
      hucpot = thupot + 2.0d0*coul(nchgta,rt00)
      aucpot = aupot  + 2.0d0*coul(nchgta,rt00)

!watanabe8/30      CALL d_ipot(p,rt00,upot,nchgta,massal,ichg,inds,iclst,ipot)

!-----------------------------------------------------------------------
!     Coulomb barrier
!-----------------------------------------------------------------------

!     coulom(1) = cbar(dble(nchgta),dble(nmasta),1,1) ! p
      
      coulom(1) = cgbar(dble(nchgta),dble(nmasta),1,1,peex) ! p
      coulom(2) = cbar(dble(nchgta),dble(nmasta),1,2) ! d
      coulom(3) = cbar(dble(nchgta),dble(nmasta),1,3) ! t
      coulom(4) = cbar(dble(nchgta),dble(nmasta),2,3) ! h
      coulom(5) = cbar(dble(nchgta),dble(nmasta),2,4) ! a
!-----------------------------------------------------------------------

      px_tmp_max = 0.0d0
      py_tmp_max = 0.0d0
      pz_tmp_max = 0.0d0
      px_tmp_min = 0.0d0
      py_tmp_min = 0.0d0
      pz_tmp_min = 0.0d0

      do i = 1,nnn
       do j = 1,5
        inf(j,i) = 0.0d0
       end do
      end do

      nemit = 0
!-----------------------------------------------------------------------
!     calculate residual nuclei ( for debug.) fuku
!-----------------------------------------------------------------------
      emitt = 0.0d0
      ttp  = 0.0d0
      ttn  = 0.0d0
      ttdeu = 0.0d0
      tttri = 0.0d0
      tthel = 0.0d0
      ttalp = 0.0d0
      if( iproj .eq. 0 .or. iproj .eq. 1 )then
      pin(3) = dsqrt( ein**2 + 2.0d0*ein*rmass )
      end if

      resp(1) = 0.0d0
      resp(2) = 0.0d0
      resp(3) = pin(3)

      if( iproj .eq. 1 ) then   !t.mori0601        
         nchgres = nchgta + 1
         nneures = nmasta - nchgta
      elseif( iproj .eq. 0 ) then
         nchgres = nchgta
         nneures = nmasta - nchgta + 1
      elseif ( iproj .eq. 5) then
         nchgres = nchgta + 2
         nneures = nmasta - nchgta + 2
      endif
!--------------------------------------------------end of debug processing     

      do i = 1, massal

!**************************************************sawada
       if( p(5,i) .lt. 0.10d0 )cycle
       if( iclst(i) .eq. -1 )cycle     !sawada
!******************************************************

       if( inds(i) .eq. 1 ) then    !sawada8/27

!-----------------------------------------------------------------------
!        Deuteron case
!-----------------------------------------------------------------------
         if( iclst(i) .eq. 1 .and. ipot(i) .eq. 1 ) then

            if( abs(p(1,i)) .gt. px_tmp_min
     &           .and. abs(p(1,i)) .lt. px_tmp_max
     &           .and.
     &           abs(p(2,i)) .gt. py_tmp_min
     &           .and. abs(p(2,i)) .lt. py_tmp_max
     &           .and.
     &           abs(p(3,i)) .gt. pz_tmp_min
     &           .and. abs(p(3,i)) .lt. pz_tmp_max ) then

               cycle

            else

               tinn1 = p(4,i) - p(5,i)
               pinn1 = sqrt( tinn1*(tinn1 + 2.0d0*ddmass) )
               tout1 = p(4,i) - p(5,i) - upot*2.!dupot
               pout1 = sqrt( tout1*(tout1 + 2.0d0*ddmass) )
               px_tmp = abs( p(1,i) )
               py_tmp = abs( p(2,i) )
               pz_tmp = abs( p(3,i) )

               do k = 1, 3
                  p(k,i) = p(k,i) * pout1 / pinn1
               end do

               pabs   = sqrt( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )
               
              !--------------------------------------------------
              ! Coulomb barrier
              !--------------------------------------------------

               tdeuteron = -ddmass + sqrt( pabs**2 + ddmass**2 )

              !--------------------------------------------------

               p(4,i) = dsqrt( pabs**2 + ddmass**2 )
               et     = p(4,i)
               ichg(i) = 1
               ineu(i) = 1

               nemit = nemit + 1

               inf(1,nemit) = dble(ichg(i))
               inf(2,nemit) = dble(ineu(i))
               inf(3,nemit) = p(1,i)
               inf(4,nemit) = p(2,i)
               inf(5,nemit) = p(3,i)
!               inf(3,nemit) = p(1,i)*2.0d0
!               inf(4,nemit) = p(2,i)*2.0d0
!               inf(5,nemit) = p(3,i)*2.0d0

               px_tmp_max = px_tmp + eps
               px_tmp_min = px_tmp - eps
               py_tmp_max = py_tmp + eps
               py_tmp_min = py_tmp - eps
               pz_tmp_max = pz_tmp + eps
               pz_tmp_min = pz_tmp - eps

!-------------------------------------------
! for debug eex fuku
!-------------------------------------------
               ttdeu = ttdeu - ddmass
     &              + dsqrt( ddmass**2 + 
     &               ( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 ) )               
               nchgres = nchgres - ichg(i)
               nneures = nneures - ineu(i)
               resp(1) = resp(1) - p(1,i)
               resp(2) = resp(2) - p(2,i)
               resp(3) = resp(3) - p(3,i)               
!--------------------------------------------
               
               cycle

            end if

!-----------------------------------------------------------------------
!         Triton case
!-----------------------------------------------------------------------

         elseif( iclst(i) .eq. 2 .and. ipot(i) .eq. 1 ) then

            if( abs(p(1,i)) .gt. px_tmp_min
     &           .and. abs(p(1,i)) .lt. px_tmp_max
     &           .and.
     &           abs(p(2,i)) .gt. py_tmp_min
     &           .and. abs(p(2,i)) .lt. py_tmp_max
     &           .and.
     &           abs(p(3,i)) .gt. pz_tmp_min
     &           .and. abs(p(3,i)) .lt. pz_tmp_max ) then

               cycle

            else

                  tinn1 = p(4,i) - p(5,i)
                  pinn1 = sqrt( tinn1*(tinn1 + 2.0d0 * tmass) )
                  tout1 = p(4,i) - p(5,i) - upot*3.!thupot
                  pout1 = sqrt( tout1*(tout1 + 2.0d0 * tmass) )

               px_tmp = abs(p(1,i))
               py_tmp = abs(p(2,i))
               pz_tmp = abs(p(3,i))

               do k = 1, 3
                  p(k,i) = p(k,i) * pout1 / pinn1
               end do

               pabs   = sqrt( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )

              !--------------------------------------------------
              ! Coulomb barrier
              !--------------------------------------------------

               ttriton = -tmass + sqrt( pabs**2 + tmass**2 )

!               if( ttriton .le. coulom(3) ) cycle

              !--------------------------------------------------

               p(4,i) = sqrt( pabs**2 + tmass**2 )
               et     = p(4,i)

               ichg(i) = 1
               ineu(i) = 2

               nemit = nemit + 1

               inf(1,nemit) = dble(ichg(i))
               inf(2,nemit) = dble(ineu(i))
               inf(3,nemit) = p(1,i)
               inf(4,nemit) = p(2,i)
               inf(5,nemit) = p(3,i)

               px_tmp_max = px_tmp + eps
               px_tmp_min = px_tmp - eps
               py_tmp_max = py_tmp + eps
               py_tmp_min = py_tmp - eps
               pz_tmp_max = pz_tmp + eps
               pz_tmp_min = pz_tmp - eps
               
!-------------------------------------------
! for debug eex fuku
!-------------------------------------------
              tttri = tttri - tmass
     &              + dsqrt( tmass**2 + 
     &               ( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 ) )               
               nchgres = nchgres - ichg(i)
               nneures = nneures - ineu(i)
               resp(1) = resp(1) - p(1,i)
               resp(2) = resp(2) - p(2,i)
               resp(3) = resp(3) - p(3,i)               
!--------------------------------------------

               cycle

            end if

!-----------------------------------------------------------------------
!         3He case
!-----------------------------------------------------------------------

         elseif( iclst(i) .eq. 3 .and. ipot(i) .eq. 1 ) then

            if( abs(p(1,i)) .gt. px_tmp_min
     &           .and. abs(p(1,i)) .lt. px_tmp_max
     &           .and.
     &           abs(p(2,i)) .gt. py_tmp_min
     &           .and. abs(p(2,i)) .lt. py_tmp_max
     &           .and.
     &           abs(p(3,i)) .gt. pz_tmp_min
     &           .and. abs(p(3,i)) .lt. pz_tmp_max ) then

               cycle

            else

                  tinn1 = p(4,i) - p(5,i)
                  pinn1 = sqrt( tinn1*(tinn1 + 2.0d0*hmass) )
                  tout1 = p(4,i) - p(5,i) - upot*3.!thupot
                  pout1 = sqrt( tout1*(tout1 + 2.0d0*hmass) )

               px_tmp = abs(p(1,i))
               py_tmp = abs(p(2,i))
               pz_tmp = abs(p(3,i))

               do k = 1, 3
                  p(k,i) = p(k,i) * pout1 / pinn1
               end do

               pabs   = sqrt( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )

              !--------------------------------------------------
              ! Coulomb barrier
              !--------------------------------------------------

               thelion = -hmass + sqrt( pabs**2 + hmass**2 )

!               if( thelion .le. coulom(4) ) cycle

              !--------------------------------------------------

               p(4,i) = sqrt( pabs**2 + hmass**2 )
               et     = p(4,i)

               ichg(i) = 2
               ineu(i) = 1

               nemit = nemit + 1

               inf(1,nemit) = dble(ichg(i))
               inf(2,nemit) = dble(ineu(i))
               inf(3,nemit) = p(1,i)
               inf(4,nemit) = p(2,i)
               inf(5,nemit) = p(3,i)

               px_tmp_max = px_tmp + eps
               px_tmp_min = px_tmp - eps
               py_tmp_max = py_tmp + eps
               py_tmp_min = py_tmp - eps
               pz_tmp_max = pz_tmp + eps
               pz_tmp_min = pz_tmp - eps

!-------------------------------------------
! for debug eex fuku 
!-------------------------------------------
               tthel = tthel - hmass
     &              + dsqrt( hmass**2 + 
     &               ( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 ) )               
               nchgres = nchgres - ichg(i)
               nneures = nneures - ineu(i)
               resp(1) = resp(1) - p(1,i)
               resp(2) = resp(2) - p(2,i)
               resp(3) = resp(3) - p(3,i)               
!--------------------------------------------

               cycle

            end if

!-----------------------------------------------------------------------
!         4He case
!-----------------------------------------------------------------------

         elseif( iclst(i) .eq. 4 .and. ipot(i) .eq. 1 ) then

            if( abs(p(1,i)) .gt. px_tmp_min
     &           .and. abs(p(1,i)) .lt. px_tmp_max
     &           .and.
     &           abs(p(2,i)) .gt. py_tmp_min
     &           .and. abs(p(2,i)) .lt. py_tmp_max
     &           .and.
     &           abs(p(3,i)) .gt. pz_tmp_min
     &           .and. abs(p(3,i)) .lt. pz_tmp_max ) then

               cycle

            else

                  tinn1 = p(4,i) - p(5,i)
                  pinn1 = dsqrt( tinn1*(tinn1 + 2.0d0*amass) )
                  tout1 = p(4,i) - p(5,i) - upot*4.
                  pout1 = dsqrt( tout1*(tout1 + 2.0d0*amass) )

               px_tmp = abs(p(1,i))
               py_tmp = abs(p(2,i))
               pz_tmp = abs(p(3,i))

               do k = 1, 3
                  p(k,i) = p(k,i) * pout1 / pinn1
               end do

               pabs   = dsqrt( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )

              !--------------------------------------------------
              ! Coulomb barrier
              !--------------------------------------------------

               talpha = -amass + sqrt( pabs**2 + amass**2 )

!               if( talpha .le. coulom(5) ) cycle

               !--------------------------------------------------

               p(4,i) = sqrt( pabs**2 + amass**2 )
               et     = p(4,i)

               ichg(i) = 2
               ineu(i) = 2

               nemit = nemit + 1

               inf(1,nemit) = dble(ichg(i))
               inf(2,nemit) = dble(ineu(i))
               inf(3,nemit) = p(1,i)
               inf(4,nemit) = p(2,i)
               inf(5,nemit) = p(3,i)

               px_tmp_max = px_tmp + eps
               px_tmp_min = px_tmp - eps
               py_tmp_max = py_tmp + eps
               py_tmp_min = py_tmp - eps
               pz_tmp_max = pz_tmp + eps
               pz_tmp_min = pz_tmp - eps

!-------------------------------------------
!     for debug eex fuku
!-------------------------------------------
              ttalp = ttalp - amass
     &              + dsqrt( amass**2 + 
     &               ( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 ) )               
               nchgres = nchgres - ichg(i)
               nneures = nneures - ineu(i)
               resp(1) = resp(1) - p(1,i)
               resp(2) = resp(2) - p(2,i)
               resp(3) = resp(3) - p(3,i)               
!--------------------------------------------
               cycle

            end if

!-----------------------------------------------------------------------
!         proton or neutron case
!-----------------------------------------------------------------------

         elseif( iclst(i) .eq. 0 .and. ipot(i) .eq. 1 ) then

            if( ichg(i) .eq. 0 ) ineu(i) = 1
            if( ichg(i) .eq. 1 ) ineu(i) = 0

            tinn = p(4,i) - p(5,i)
            pinn = sqrt( tinn*(tinn + 2.0d0*p(5,i)) )
            tout = p(4,i) - p(5,i) - upot
            pout = sqrt( tout*(tout + 2.0d0*p(5,i)) )

            do k = 1, 3
               p(k,i) = p(k,i) * pout / pinn
            end do

            pabs = sqrt( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )
               p(4,i) = sqrt( pabs**2 + rmass**2 )
               et     = p(4,i)

!-----------------------------------------------------------------------
!             proton
!-----------------------------------------------------------------------

!            if( ichg(i) .eq. 1 .and. et - rmass .le. ein ) then !hiroki
            if( ichg(i) .eq. 1 ) then !hiroki

              !--------------------------------------------------
              ! Coulomb barrier
              !--------------------------------------------------

               tproton = -rmass + et

               !if( tproton .le. coulom(1) ) cycle

              !--------------------------------------------------


               nemit = nemit + 1

               inf(1,nemit) = dble(ichg(i))
               inf(2,nemit) = dble(ineu(i))
               inf(3,nemit) = p(1,i)
               inf(4,nemit) = p(2,i)
               inf(5,nemit) = p(3,i)
!-------------------------------------------
!     for debug eex fuku
!-------------------------------------------
               ttp = ttp - rmass
     &              + dsqrt( rmass**2
     &              + ( p(1,i)**2 + p(2,i)**2 + p(3,i)**2) )
               nchgres = nchgres - ichg(i)
               nneures = nneures - ineu(i)
               resp(1) = resp(1) - p(1,i)
               resp(2) = resp(2) - p(2,i)
               resp(3) = resp(3) - p(3,i)               
!--------------------------------------------
!-----------------------------------------------------------------------
!             neutron
!-----------------------------------------------------------------------

!            elseif( ineu(i) .eq. 1 .and. et - rmass .le. ein ) then
            elseif( ineu(i) .eq. 1 ) then

               nemit = nemit + 1

               inf(1,nemit) = dble(ichg(i))
               inf(2,nemit) = dble(ineu(i))
               inf(3,nemit) = p(1,i)
               inf(4,nemit) = p(2,i)
               inf(5,nemit) = p(3,i)
!-------------------------------------------
!     for debug eex fuku
!-------------------------------------------
              ttn = ttn - rmass
     &              + dsqrt( rmass**2
     &              + ( p(1,i)**2 + p(2,i)**2 + p(3,i)**2) )               
               nchgres = nchgres - ichg(i)
               nneures = nneures - ineu(i)
               resp(1) = resp(1) - p(1,i)
               resp(2) = resp(2) - p(2,i)
               resp(3) = resp(3) - p(3,i)               
!--------------------------------------------
            end if

!-----------------------------------------------------------------------

         end if

       elseif( inds(i) .eq. 4 .and. ipot(i) .eq. 1 ) then  !!!sawada8/27
!-----------------------------------------------------------------------
!     output pion part
!-----------------------------------------------------------------------
!     ichg =  1 :: pion+ ( kf = 211 )
!     ichg = -1 :: pion- ( kf =-211 )
!     ichg =  0 :: pion0 ( kf = 111 )
!-----------------------------------------------------------------------

            tinn = p(4,i) - p(5,i)
            pinn = sqrt( tinn*(tinn + 2.0d0*pmass*1000.0) )
            tout = p(4,i) - p(5,i)
            pout = sqrt( tout*(tout + 2.0d0*pmass*1000.0) )

            do k = 1, 3
               p(k,i) = p(k,i) * pout / pinn
            end do

            pabs = sqrt( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )
            et   = sqrt( pabs**2 + ( pmass*1000.0 )**2 )

            if( ichg(i) .eq. 1 ) then

             nemit = nemit + 1

             inf(1,nemit) = 0
             inf(2,nemit) = 0
             inf(3,nemit) = p(1,i)
             inf(4,nemit) = p(2,i)
             inf(5,nemit) = p(3,i)

            elseif( ichg(i) .eq. 0 ) then

             nemit = nemit + 1

             inf(1,nemit) = 0
             inf(2,nemit) = -1
             inf(3,nemit) = p(1,i)
             inf(4,nemit) = p(2,i)
             inf(5,nemit) = p(3,i)

            elseif( ichg(i) .eq. -1 ) then

             nemit = nemit + 1

             inf(1,nemit) = -1
             inf(2,nemit) = -1
             inf(3,nemit) = p(1,i)
             inf(4,nemit) = p(2,i)
             inf(5,nemit) = p(3,i)
            end if

       end if                             !sawada8/27

      end do

!---eex calculate for debug fuku-----      
!-----------------------------------------------------------------------
!     mass of remnants from cascade process (for debug)
!-----------------------------------------------------------------------

      nmasres  = nchgres + nneures

!-----------------------------------------------------------------------
!     mass energy of remnants from cascade process (for debug)
!-----------------------------------------------------------------------

      rmasres  = rmass * dble( nchgres + nneures )

!-----------------------------------------------------------------------
!     momentum of residual nuclei from cascade process (for debug)
!-----------------------------------------------------------------------

      if( dabs(resp(1)) .lt. eps ) resp(1) = 0.0d0
      if( dabs(resp(2)) .lt. eps ) resp(2) = 0.0d0
      if( dabs(resp(3)) .lt. eps ) resp(3) = 0.0d0

      resppp  = dsqrt( resp(1)**2 + resp(2)**2 + resp(3)**2 )

!-----------------------------------------------------------------------
!     recoil energy of remnant from cascade process (for debug)
!-----------------------------------------------------------------------

      resrec  = - rmasres  + dsqrt( rmasres**2 + resppp**2 )

      emitt = ttp + ttn + ttdeu + tttri + tthel + ttalp
      
!---------------------------------------------------------------------------
!     excitation energy of residual nuclei from cascade process (for debug)
!     (2)Bertini
!---------------------------------------------------------------------------    
               resext = ein - resrec - emitt
!     &        - (eliq(nchgta,nmasta-nchgta) - eliq(nchgres,nneures))

               if( resext .lt. 1.0d-10 )resext = 0.0d0
               eex = resext
!-----------------------------------------------------end of debug processing
               
      END SUBROUTINE sm_sumo


!***********************************************************************
!*                                                                     *
      real*8  FUNCTION cbar(z,a,iz,ia)
!*                                                                     *
!*      Calculate Coulomb potential                                    *
!*                                                                     *
!*     input :                                                         *
!*       a   :   mass of nucleus #1                    (IN)            *
!*       z   :   charge  of nucleus #1                 (IN)            *
!*      ia   :   mass of nucleus #2 (emitted)          (IN)            *
!*      iz   :   charge  of nucleus #2 (emitted)       (IN)            *
!*                                                                     *
!*     output :                                                        *
!*     cbar  :   Coulomb potential  [MeV]              (OUT)           *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      parameter( rc=1.70d0, ee=137.0359895d0, hbarc=197.327053d0 )

      cbar=hbarc/ee
      rcal = 10.0

      if( rcal .eq. 10.0 ) then
         rr = 1.5d0
      else
         rr = rcal
      endif

      r1 = rr*a**.333333d0
      r2 = rr*dble(ia)**.333333d0
      r0 = r1+r2

      cbar = cbar * dble(iz) * z / r0

      END FUNCTION cbar

!***********************************************************************
!*                                                                     *
      SUBROUTINE genseed(iseed)
!*                                                                     *
!*                                                                     *
!*       Purpose:                                                      *
!*                                                                     *
!*             to generate seed from data_and_time                     *
!*                                                                     *
!*                                                                     *
!***********************************************************************

!      g77 can't use CALL date_and_time

!      dimension  dt(8)
!      CALL date_and_time(values=idt)
!      iseed = idt(6) * 100 + idt(7)
      iseed = 100

!      iseed = 10.0*rn(0) !idt(6) * 100 + idt(7)
!      iseed = int( rn(iseed) * 10.0 * rn(iseed)*10.0 )

      END SUBROUTINE genseed


!***********************************************************************
! to calculate inelastic collidion,some subroutines are added  by sawada
!*                                                                     *
!***********************************************************************


************************************************************************
*                                                                      *
!      subroutine pionemm(tfm,ichg,inuc,ibry,inds,inun,iavd,ihis,ipot,
!     &      iclst,nchgta,massal,massba,upot,ucpot,dt,r,p,rt00,lcoll)
!      subroutine pionemm(tfm,ichg,inuc,ibry,inds,inun,iavd,ihis,ipot,
!     &     iclst,nchgta,massal,massba,upot,ucpot,dt,r,p,rt00,lcoll,
!     &     igroup,nmascl,ipotcc,ccp)            
      subroutine pionemm(tfm,ichg,inuc,ibry,inds,inun,iavd,ihis,ipot,
     &     iclst,nchgta,massal,massba,upot,ucpot,dt,r,p,rt00,lcoll,
     &     numclst,nmasclst,pcmclst) ! 2022/10/5 yamaguchi
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate the decay of delta or N*                   *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              dt          : time interval for decay                   *
*                                                                      *
*        Comments : counting variables                                 *
*                                                                      *
*              lcoll(22) : D + N -> pi                                 *
*              lcoll(23) : R + N -> pi                                 *
*              lcoll(24) : R + D -> pi                                 *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

!      include 'param00.inc'
!      include 'param01.inc'
!      include 'param02.inc'
      include 'param-incelf.inc'
*-----------------------------------------------------------------------

      parameter     ( epse = 0.0001 )
      parameter     ( prmin = 0.0001 )
      parameter     ( hbc = 0.1973   )

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

      dimension       beta(3), pcm(3)

!*************************************************************sawada8/27
      integer   ichg(nnn),inuc(nnn),ibry(nnn),inds(nnn),inun(nnn)
      integer   ihis(nnn),iavd(nnn),ipot(nnn)
      integer   iclst(nnn) !(S.H. on 2012.12.13)
      integer   lcoll(30)

      integer     massal,massba
      real*8      r(3,nnn),p(6,nnn),dt,tfm,upot,ucpot,rt00

!      integer, intent(inout) :: igroup(nnn),nmascl(nnn),ipotcc(nnn)
      integer, intent(inout) :: numclst(nnn),nmasclst(nnn)
      real*8, intent(inout) :: pcmclst(5,nnn)      
!***********************************************************************

      CALL genseed(iseed)
!-----------------------------------------------------------------------
!     final time of time evolution
!-----------------------------------------------------------------------

      tfin = 70.0 * ( dble(massba)/208.0 )**0.16
      dtmax = tfin

!-------------------------------------------------------------sawada8/27
!     here, p = [MeV]. [MeV]->[GeV].
      do k = 1,massal
        do j = 1,5
         p(j,k) = p(j,k)/1000.0d0
        end do
      end do
!-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*        output unit
*-----------------------------------------------------------------------

            ieo = 6

*-----------------------------------------------------------------------
*        zero set of counter
*-----------------------------------------------------------------------

                  lcoll(22)= 0
                  lcoll(23)= 0
                  lcoll(24)= 0

*-----------------------------------------------------------------------
*        factor of N* -> D + pi to N* -> N + pi
*-----------------------------------------------------------------------

                  frdp = 0.4

*-----------------------------------------------------------------------

      do 800 i  = 1, massba

            if( inds(i) .ne. 2 .and. inds(i) .ne. 3 ) goto 800

*-----------------------------------------------------------------------
                  px1 = p(1,i) 
                  py1 = p(2,i) 
                  pz1 = p(3,i) 
                  e1  = p(4,i) 
                  em1 = p(5,i) 

                  iz1 = ichg(i)
                  in1 = inuc(i)
                  id1 = inds(i)
                  iu1 = inun(i)

                  beta(1)  = - px1 / e1
                  beta(2)  = - py1 / e1
                  beta(3)  = - pz1 / e1
                  betasq   =   beta(1)**2 + beta(2)**2 + beta(3)**2
                  gamma    =   1.0 / sqrt( 1.0 - betasq )

*-----------------------------------------------------------------------
*           can the Delta or N* decay in this time step ?
*-----------------------------------------------------------------------

                  dt0  = dt / gamma
!-------------------------------------------------------------sawada8/27
!                 [GeV]
                  dem1 = em1
!-----------------------------------------------------------------------
                  dem2 = 0.0

                  call resmass(0,id1,1,em1,dem1,dem2,gam,fqrq)

                  w    = exp( - dt0 * gam / hbc )


            if( dt .le. dtmax .and. rn(0) .lt. w )  goto 800

*-----------------------------------------------------------------------

               if( id1 .eq. 2 ) then

*                    ---------  Delta -> Nucleon + pi ---------
!-------------------------------------------------------------sawada8/27
                     dnmass = rmass/1000.0
!-----------------------------------------------------------------------

                     inid   = 1
                     idid   = 1

               else if( id1 .eq. 3 ) then

!-------------------------------------------------------------sawada8/27
                  if( ( em1 .gt. rmass/1000.0+2.0 * pmass ) .and.
     &                ( rn(0) .lt. frdp )) then
!-----------------------------------------------------------------------

*                    ---------  N* -> Delta + pi ---------
!-------------------------------------------------------------sawada8/27
                     dem1 = rmass/1000.0 + pmass
                     dem2 = pmass
!-----------------------------------------------------------------------


                     call resmass(1,2,1,em1,dem1,dem2,gam0,fqrq)

                     dnmass = dem1
                     inid   = 0
                     idid   = 2

                  else
*                    ---------  N* -> Nucleon + pi ---------
!-------------------------------------------------------------sawada8/27
                     dnmass = rmass/1000.0
!-----------------------------------------------------------------------

                     inid   = 1
                     idid   = 1

                  end if

               end if

*-----------------------------------------------------------------------
*           now the Delta or N* may decay
*-----------------------------------------------------------------------

                  j = massal + 1

                  if( j .gt. nnn ) then

                     write(ieo,'('' Error: too many pions,'',
     &                           '' use larger nnn value.'')')
                     write(ieo,'('' ====='')')

                     call parastop( 222 )

                  end if


*-----------------------------------------------------------------------
*           check potential energy
*-----------------------------------------------------------------------

                  eini =        e1
                  etwo =        e1

*-----------------------------------------------------------------------
*           charge and state identification
*-----------------------------------------------------------------------

                        ichg(j) = 0

                        xx = rn(0)

               if( inds(i) .eq. 2 ) then  !delta

                  if( ichg(i) .eq. 2 ) then

                        ichg(i) =  1
                        ichg(j) =  1

                  else if( ichg(i) .eq. -1 ) then

                        ichg(i) =  0
                        ichg(j) = -1

                  else

                     if( xx .gt. 0.66666667 ) then

                        ichg(j) =  2 * ichg(i) - 1
                        ichg(i) =  1 - ichg(i)

                     end if

                  end if

               else if( inid .eq. 1 ) then

                     if( xx .gt. 0.33333 ) then

                        ichg(j) =  2 * ichg(i) - 1
                        ichg(i) =  1 - ichg(i)

                     end if

               else

                  if( ichg(i) .eq. 1 ) then

                     if( xx .lt. 0.5 ) then

                        ichg(i) =  2
                        ichg(j) = -1

                     else if( xx .lt. 0.66666667 ) then

                        ichg(i) =  0
                        ichg(j) =  1

                     else

                        ichg(i) =  1
                        ichg(j) =  0

                     end if

                  else

                     if( xx .lt. 0.5 ) then

                        ichg(i) = -1
                        ichg(j) =  1

                     else if( xx .lt. 0.66666667 ) then

                        ichg(i) =  1
                        ichg(j) = -1

                     else

                        ichg(i) =  0
                        ichg(j) =  0

                     end if

                  end if

               end if

*-----------------------------------------------------------------------

                        inuc(i) = inid
                        inds(i) = idid
                        inun(i) = j
                        p(5,i)  = dnmass

                        inuc(j) = 0
                        inun(j) = i
                        inds(j) = 4
                        p(5,j)  = pmass

                        massal  = massal + 1

*-----------------------------------------------------------------------
*           the outgoing nucleon momentum : pcm(i) in the Delta cms
*           and determine the pion position
*-----------------------------------------------------------------------

                  ntag = 1

                  pr = max( prmin,
     &                     .25 * ( em1**2 - dnmass**2 - pmass**2 )**2
     &                                    - dnmass**2 * pmass**2 )
                  pr = sqrt( pr ) / em1

                  rr = 10.0

               do while( ( rr .lt. 0.001 ) .or. ( rr .gt. 1.0 ) )

                  xx = 1. - 2. * rn(0)
                  yy = 1. - 2. * rn(0)
                  zz = 1. - 2. * rn(0)
                  rr = sqrt( xx**2 + yy**2 + zz**2 )

               end do

                  pcm(1)  = pr * xx / rr
                  pcm(2)  = pr * yy / rr
                  pcm(3)  = pr * zz / rr

                  r(1,j) = r(1,i)
                  r(2,j) = r(2,i)
                  r(3,j) = r(3,i)

*-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     p(5,i)=dnmass,p(5,j)=pmass [GeV] therefore p(1,i)-p(4,i) = [GeV]
!-----------------------------------------------------------------------

         do 7000 kk = 1, 4

*-----------------------------------------------------------------------
*           Lorentz-transformation into reference frame
*-----------------------------------------------------------------------

                  e1cm   = sqrt( dnmass**2 + pcm(1)**2
     &                                     + pcm(2)**2
     &                                     + pcm(3)**2 )

                  p1beta = pcm(1) * beta(1)
     &                   + pcm(2) * beta(2)
     &                   + pcm(3) * beta(3)

                  transf = gamma * ( gamma * p1beta
     &                   / ( gamma + 1 ) - e1cm )

                  p(1,i) = beta(1) * transf + pcm(1)
                  p(2,i) = beta(2) * transf + pcm(2)
                  p(3,i) = beta(3) * transf + pcm(3)

                  p(4,i) = sqrt( p(5,i)**2 + p(1,i)**2
     &                         + p(2,i)**2 + p(3,i)**2 )


                  e2cm   = sqrt( pmass**2 + pcm(1)**2
     &                                    + pcm(2)**2
     &                                    + pcm(3)**2 )

                  transf = gamma * ( -gamma * p1beta 
     &                   / ( gamma + 1 ) - e2cm )

                  p(1,j) = beta(1) * transf - pcm(1)
                  p(2,j) = beta(2) * transf - pcm(2)
                  p(3,j) = beta(3) * transf - pcm(3)

                  p(4,j) = sqrt( p(5,j)**2 + p(1,j)**2
     &                         + p(2,j)**2 + p(3,j)**2 )

*-----------------------------------------------------------------------
*   CHECK ENERGY OF TWO PATICLES
*-----------------------------------------------------------------------

                  efin =  p(4,i) + p(4,j)
*-----------------------------------------------------------------------

            if( abs( eini - efin ) .gt. epse ) then

                  cona = ( eini - efin + etwo ) / gamma

                  fac2 = 1. / ( 4.0 * cona**2 * pr**2 )
     &                 * ( ( cona**2 - ( p(5,i)**2 + p(5,j)**2 ) )**2
     &                     - 4.0 * p(5,i)**2 * p(5,j)**2 )

               if( fac2 .gt. 0.0 ) then

                  fact = sqrt( fac2 )

                  pcm(1) = pcm(1) * fact
                  pcm(2) = pcm(2) * fact
                  pcm(3) = pcm(3) * fact

               else

                  if( dt .gt. dtmax ) ntag = 0
                  goto 9000

               end if

            else

                  ntag = 0
                  goto 9000

            end if

 7000    continue

*-----------------------------------------------------------------------
*           check final pauli-blocking of out going nucleon
*-----------------------------------------------------------------------

 9000          if( ntag .eq. 0 .and. dt .lt. dtmax .and. inds(i) .eq. 1
     &                 .and. p(4,i)-p(5,i) .lt. tfm/1000.0 )then
                    ntag = 1
               end if

!-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*              it does not decay, reset variables
*-----------------------------------------------------------------------

               if( ntag .eq. 1 ) then

!-------------------------------------------------------------sawada8/27
                  p(1,i)  = px1 
                  p(2,i)  = py1
                  p(3,i)  = pz1
                  p(4,i)  = e1 
                  p(5,i)  = em1
!-----------------------------------------------------------------------
                  ichg(i) = iz1
                  inuc(i) = in1
                  inun(i) = iu1
                  inds(i) = id1

                  ichg(j) = 0
                  inuc(j) = 0
                  inun(j) = 0
                  inds(j) = 0

                  massal  = massal - 1
*-----------------------------------------------------------------------
*              Decay, count the summary variables
*-----------------------------------------------------------------------

               else

                     iavd(i) = iavd(i) + iavd(i)/abs(iavd(i))

                     ihis(i) = - abs(ihis(i)) - 1

                  if( id1 .eq. 2 ) ihis(j) =   abs(ihis(i))
                  if( id1 .eq. 3 ) ihis(j) = - abs(ihis(i))

                  if( id1 .eq. 2 ) then

                     lcoll(22) = lcoll(22) + 1
                  else if( inid .eq. 1 ) then

                     lcoll(23) = lcoll(23) + 1

                  else if( inid .eq. 0 ) then

                     lcoll(24) = lcoll(24) + 1

                  end if

               end if
*-----------------------------------------------------------------------

  800  continue

*-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     here, p = [GeV]. [GeV]->[MeV].
       do k = 1,massal
        do j = 1,5
         p(j,k) = p(j,k)*1000.0d0
        end do
       end do

!      CALL d_ipot(p,rt00,upot,nchgta,massal,ichg,inds,iclst,ipot)
!       CALL d_ipot(p,rt00,upot,nchgta,massal,ichg,inds,iclst,ipot,
!     &      igroup,nmascl,ipotcc,ccp)
       CALL d_ipot(p,rt00,upot,nchgta,massal,ichg,inds,iclst,ipot,
     &      numclst,nmasclst,pcmclst)

      return
      end


************************************************************************
*                                                                      *
!      subroutine pionabb(ein,ichg,inuc,ibry,inds,inun,iavd,ihis,inio,
!     &                   iclst,nmasta,massal,massba,dt,r,p,lcoll)
      subroutine pionabb(ein,ichg,inuc,ibry,inds,inun,iavd,ihis,inio,
     &     iclst,nmasta,massal,massba,dt,r,p,lcoll,ininclst)      
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate pion absorption                            *
*                                                                      *
*        Comments : counting variables                                 *
*                                                                      *
*              lcoll(25) : N + pi -> D                                 *
*              lcoll(26) : N + pi -> R                                 *
*              lcoll(27) : D + pi -> R                                 *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------
      include 'param-incelf.inc'

*-----------------------------------------------------------------------
      parameter     ( epse = 0.0001 )
*-----------------------------------------------------------------------
      parameter     ( delpi  = 4.0 )
      parameter     ( pirr   = 0.977, pir0  = 2.07)
      parameter     ( pir1   = 2.52,  pir2  = 1.49)
*-----------------------------------------------------------------------
*
*     delpi         : maximum spatial distance for which a absorption
*                     still can occur only for 200 MeV
*
*     bcmax         : maximum impact parameter for
*                     pirr = 0.977 fm corresponds  to  30 mb
*                     pir0 = 2.07  fm corresponds  to 135 mb
*                     pir1 = 2.52  fm corresponds  to 200 mb
*                     pir2 = 1.49  fm corresponds  to  70 mb
*
*-----------------------------------------------------------------------

      integer massal,massba,nmasta
      integer ichg(nnn),inuc(nnn),ibry(nnn),inds(nnn),
     &     inun(nnn),iavd(nnn),ihis(nnn),inio(nnn),iclst(nnn),
     &     ininclst(nnn)      
      integer lcoll(30)

      real*8 r(3,nnn),p(6,nnn),ein,dt
      real*8 rr2(nnn,nnn)
      save srt,px1,py1 !(S.H. on 2012.12.16)
!$OMP THREADPRIVATE(srt,px1,py1)
      elab = ein/1000.0

      CALL genseed(iseed)

*-----------------------------------------------------------------------
*        zero set of counter
*-----------------------------------------------------------------------

               lcoll(25) = 0
               lcoll(26) = 0
               lcoll(27) = 0

      if( massal .eq. massba ) return

!-------------------------------------------------------------sawada8/27
!     [MeV]->[GeV] 
       do i = 1,massal
         do j = 1,5
          p(j,i) = p(j,i) / 1000.0
         end do
       end do
!-----------------------------------------------------------------------

*-----------------------------------------------------------------------

      do 800 i1  = massba + 1, massal

                  px1 = p(1,i1)
                  py1 = p(2,i1)
                  pz1 = p(3,i1)

                  e1  = p(4,i1)
                  em1 = p(5,i1)

                  iz1 = ichg(i1)
                  in1 = inuc(i1)
                  id1 = inds(i1)
                  iu1 = inun(i1)
                  ia1 = iavd(i1)
                  ih1 = ihis(i1)
*-----------------------------------------------------------------------
*        look for an absorbent nucleon
*-----------------------------------------------------------------------

         do 600 i2 = massba, 1, -1


               if( inio(i2) .eq. 1 ) cycle
               if( iclst(i2) .ne. 0 )cycle
               if( ininclst(i2) .ne. 0 )cycle
               
*-----------------------------------------------------------------------
*           avoid second collisions for the same pairs
*-----------------------------------------------------------------------

               if( i2 .eq. inun(i1) .and.
     &             i1 .eq. inun(i2) ) goto 600

*-----------------------------------------------------------------------
*           look only for nucleons or Delta + pion with proper charge
*-----------------------------------------------------------------------

               if( inds(i2) .eq. 3 ) goto 600

               if( inds(i2) .eq. 2 .and.
     &           ( iz1 + ichg(i2) + 2 ) / 2 .ne. 1 ) goto 600

*-----------------------------------------------------------------------
*           delpi is appled only for 5000 MeV
*-----------------------------------------------------------------------

                  CALL rr(i1,i2,r,p,rr2)

               if( elab .le. 5.0 .and.
     &             rr2(i1,i2) .gt. delpi**2 ) goto 600

*-----------------------------------------------------------------------
*           check impact parameter and collision time
*-----------------------------------------------------------------------

                  px2 = p(1,i2)
                  py2 = p(2,i2)
                  pz2 = p(3,i2)

                  e2  = p(4,i2)
                  em2 = p(5,i2)

                  iz2 = ichg(i2)
                  in2 = inuc(i2)
                  id2 = inds(i2)
                  iu2 = inun(i2)
                  ia2 = iavd(i2)
                  ih2 = ihis(i2)

                  s   = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2
     &                                             - (pz1+pz2)**2
                  srt = sqrt(s)
*-----------------------------------------------------------------------
*           impact parameter
*-----------------------------------------------------------------------

                  dx   = r(1,i1) - r(1,i2)
                  dy   = r(2,i1) - r(2,i2)
                  dz   = r(3,i1) - r(3,i2)
                  rsq  = dx**2 + dy**2 + dz**2

                  p12  = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
                  p1dr = px1 * dx + py1 * dy + pz1 * dz
                  p2dr = px2 * dx + py2 * dy + pz2 * dz
                  a12  = 1.0 - ( em1 * em2 / p12 ) ** 2
                  b12  = p1dr / em1 - p2dr * em1 / p12
                  c12  = rsq + ( p1dr / em1 )**2

                  brel = sqrt( abs(c12 - b12**2/a12) )

*-----------------------------------------------------------------------
*           bcmax of maximum cross section
*-----------------------------------------------------------------------

                  izzz = ( iz2 + iz1 + 2 ) / 2

                  pind = 0.0
                  pird = 0.0

               if( inds(i2) .eq. 1 ) then

                  pind = pir0**2

                  if( izzz .ne. 1 )   pind = pir1**2

                  if( ( izzz .eq. 1 ) .and.
     &                ( iz1  .ne. 0 ) ) pind = pir2**2

                  if( izzz .eq. 1 ) pird = pirr**2

               else

                  pird = pirr**2

               end if

                  pirn = sqrt( pind + pird )

*-----------------------------------------------------------------------
*           is their impact parameter small enough ?
*-----------------------------------------------------------------------

               if( brel .gt. pirn ) goto 600

*-----------------------------------------------------------------------
*           average time-shift of the collision in the fixed frame
*           will particles get closest point in this time interval ?
*-----------------------------------------------------------------------

                  b21 =   - p2dr / em2 + p1dr * em2 / p12
                  t1  = (   p1dr / em1 - b12 / a12 ) * e1 / em1
                  t2  = ( - p2dr / em2 - b21 / a12 ) * e2 / em2

               if ( abs( t1 + t2 ) .gt. dt )  goto 600

*-----------------------------------------------------------------------
*           now the pion may be absorbed in this time step
*           Check the cross section with phase space factor
*-----------------------------------------------------------------------

                  xx = rn(0)
                  yy = 0.0
                  zz = 0.0

              if( inds(i2) .eq. 1 ) then

                  dem1 = srt
                  dem2 = 0.0

                  call resmass(0,2,1,srt,dem1,dem2,gamdl,fqrq)

                  yy = fqrq 
     &               / ( 4.0 * ( srt - dmass )**2 / gamdl**2 + 1.0 )
                  yy = ( pind * yy ) / ( pind + pird )

              end if

              if( izzz .eq. 1 ) then

                  dem1 = srt
                  dem2 = 0.0

                  call resmass(0,3,1,srt,dem1,dem2,gamrs,fqrq)

                  zz  = fqrq
     &                / ( 4. * ( srt - smass )**2 / gamrs**2 + 1.0 )
                  zz = ( pird * zz ) / ( pind + pird )

              end if

               if( inds(i2) .eq. 1 .and. xx .gt. yy + zz ) goto 600
               if( inds(i2) .eq. 2 .and. xx .gt.      zz ) goto 600

*-----------------------------------------------------------------------
*           check potential energy
*-----------------------------------------------------------------------

                  eini = e1 + e2

*-----------------------------------------------------------------------
*           set variables
*-----------------------------------------------------------------------

               if( izzz .eq. 1 .and. xx .le. zz ) then

                  inds(i2) = 3

               else

                  inds(i2) = 2

               end if

                  p(1,i2)  = px1 + px2
                  p(2,i2)  = py1 + py2
                  p(3,i2)  = pz1 + pz2
                  p(5,i2)  = srt

                  pps      = p(1,i2)**2 + p(2,i2)**2 + p(3,i2)**2
                  p(4,i2)  = sqrt( p(5,i2)**2 + pps )

                  ichg(i2) = iz2 + iz1
                  inuc(i2) = 0

                  inun(i2) = 0
                  iavd(i2) = iavd(i2) + iavd(i2)/abs(iavd(i2))
                  ihis(i2) = - abs(ihis(i2)) - abs(ihis(i1)) - 1

                  ichg(i1) = 0
                  inds(i1) = 0
                  inun(i1) = 0
                  iavd(i1) = 0
                  ihis(i1) = 0

*-----------------------------------------------------------------------

         do 7000 kk = 1, 4

*-----------------------------------------------------------------------
*           check final energy
*-----------------------------------------------------------------------

                  efin = p(4,i2)

*-----------------------------------------------------------------------
*           pion is absorbed, counting varables
*-----------------------------------------------------------------------

            if( abs( eini - efin ) .lt. epse ) then

                  if( inds(i2) .eq. 2 ) then

                        lcoll(25) = lcoll(25) + 1
                        lcoll(7) = lcoll(7) + 1

                  else if( inds(i2) .eq. 3 ) then

                     if( id2 .eq. 1 ) then

                        lcoll(26) = lcoll(26) + 1
                        lcoll(7) = lcoll(7) + 1

                     else if( id2 .eq. 2 ) then

                        lcoll(27) = lcoll(27) + 1
                        lcoll(7) = lcoll(7) + 1

                     end if

                  end if

                  goto 800

*-----------------------------------------------------------------------
            else

                  enew = eini

                  srt  = sqrt( enew**2 - pps )

!-------------------------------------------------------------sawada8/27
               if( srt .gt. rmass/1000.0 + pmass ) then
!-----------------------------------------------------------------------
                  p(5,i2)  = srt
                  p(4,i2) = enew

               else

                  goto 8000

               end if

            end if

 7000    continue

*-----------------------------------------------------------------------
*              reset variables
*-----------------------------------------------------------------------

 8000             p(1,i2)  = px2
                  p(2,i2)  = py2
                  p(3,i2)  = pz2

                  p(5,i2)  = em2
                  p(4,i2)  = e2 

                  ichg(i2) = iz2
                  inuc(i2) = in2
                  inds(i2) = id2
                  inun(i2) = iu2
                  iavd(i2) = ia2
                  ihis(i2) = ih2

                  ichg(i1) = iz1
                  inuc(i1) = in1
                  inds(i1) = id1
                  inun(i1) = iu1
                  iavd(i1) = ia1
                  ihis(i1) = ih1

*-----------------------------------------------------------------------

  600   continue

  800 continue

*-----------------------------------------------------------------------
*        renumbering of pion part
*-----------------------------------------------------------------------

                     inpion    = 0

         do 900 i = massba + 1, massal

            if( inds(i) .eq. 0 ) then

               do 950 j = massal, i + 1, -1

                  if( inds(j) .eq. 4 ) then

                     r(1,i)    =   r(1,j)
                     r(2,i)    =   r(2,j)
                     r(3,i)    =   r(3,j)
                     p(1,i)    =   p(1,j)
                     p(2,i)    =   p(2,j)
                     p(3,i)    =   p(3,j)
                     p(4,i)    =   p(4,j)
                     p(5,i)    =   p(5,j)
                     ichg(i)   =   ichg(j)
                     inds(i)   =   inds(j)
                     ihis(i)   =   ihis(j)
                     inds(j)   =   0
                     ihis(j)   =   0

                  do 100 k = 1, massal

                     rr2(i,k)  =   rr2 (j,k)


                    rr2 (k,i)  =   rr2 (i,k)

  100             continue

                     inpion    = inpion + 1

                     goto 900

                  end if

  950          continue

            else

                     inpion = inpion + 1

            end if

  900    continue


                     massal = massba + inpion

*-----------------------------------------------------------------------

!-------------------------------------------------------------sawada8/27
!     [GeV]->[MeV] 
      do i = 1,massal
        do j = 1,5
         p(j,i) = p(j,i) * 1000.0d0
        end do

      end do

!-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
c     subroutine fpidecayy(lcoll,upot,ucpot,tfm,ichg,inuc,ibry,inds,
CS.H. on 2012.12.13
c$$$     &           inun,iavd,ihis,massal,massba,r,p,rt00,ipot)
!     subroutine fpidecayy(lcoll,upot,ucpot,tfm,ichg,inuc,ibry,inds,
!    &           inun,iavd,ihis,massal,massba,r,p,rt00,ipot,iclst,itime,
!    &           nchgta)
!      subroutine fpidecayy(lcoll,upot,ucpot,tfm,ichg,inuc,ibry,inds,
!     &           inun,iavd,ihis,massal,massba,r,p,rt00,ipot,iclst,
!     &           nchgta)
!      subroutine fpidecayy(lcoll,upot,ucpot,tfm,ichg,inuc,ibry,inds,
!     &           inun,iavd,ihis,massal,massba,r,p,rt00,ipot,iclst,
!     &           nchgta,igroup,nmascl,ipotcc,ccp)      
      subroutine fpidecayy(lcoll,upot,ucpot,tfm,ichg,inuc,ibry,inds,
     &           inun,iavd,ihis,massal,massba,r,p,rt00,ipot,iclst,
     &           nchgta,numclst,nmasclst,pcmclst) ! 2022/10/5 yamaguchi      
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate final decay of the resonances              *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param-incelf.inc'

      integer massal,massba
      integer ichg(nnn),inuc(nnn),ibry(nnn),inds(nnn),
     &        inun(nnn),iavd(nnn),ihis(nnn)
      integer lcoll(30)
      integer, intent(inout) :: ipot(nnn)
      integer iclst(nnn) !(S.H. on 2012.12.13)

      real*8 r(3,nnn),p(6,nnn),upot,ucpot,tfm,timepi,rt00

!      integer, intent(inout) :: igroup(nnn),nmascl(nnn),ipotcc(nnn)
      integer, intent(inout) :: numclst(nnn),nmasclst(nnn)
!      real*8, intent(inout) :: ccp(6,nnn)      
      real*8, intent(inout) :: pcmclst(5,nnn)      
      
      icolt = 1

*-----------------------------------------------------------------------
*        final pion decay
*-----------------------------------------------------------------------
            if( icolt .eq. 1 ) then

                     timepi = 1000.

                     icol22 = lcoll(22)
                     icol23 = lcoll(23)
                     icol24 = lcoll(24)


!           call pionemm(tfm,ichg,inuc,ibry,inds,inun,iavd,ihis,ipot,
!     &     iclst,nchgta,massal,massba,upot,ucpot,timepi,r,p,rt00,lcoll)
!           call pionemm(tfm,ichg,inuc,ibry,inds,inun,iavd,ihis,ipot,
!     &      iclst,nchgta,massal,massba,upot,ucpot,timepi,r,p,rt00,lcoll,
!     &      igroup,nmascl,ipotcc,ccp)
           call pionemm(tfm,ichg,inuc,ibry,inds,inun,iavd,ihis,ipot,
     &      iclst,nchgta,massal,massba,upot,ucpot,timepi,r,p,rt00,lcoll,
     &      numclst,nmasclst,pcmclst) ! 2022/10/5 yamaguchi

                     ncol22 = lcoll(22)
                     ncol23 = lcoll(23)
                     ncol24 = lcoll(24)

!           call pionemm(tfm,ichg,inuc,ibry,inds,inun,iavd,ihis,ipot,
!     &     iclst,nchgta,massal,massba,upot,ucpot,timepi,r,p,rt00,lcoll)
!           call pionemm(tfm,ichg,inuc,ibry,inds,inun,iavd,ihis,ipot,
!     &      iclst,nchgta,massal,massba,upot,ucpot,timepi,r,p,rt00,lcoll,
!     &      igroup,nmascl,ipotcc,ccp)
           call pionemm(tfm,ichg,inuc,ibry,inds,inun,iavd,ihis,ipot,
     &      iclst,nchgta,massal,massba,upot,ucpot,timepi,r,p,rt00,lcoll,
     &      numclst,nmasclst,pcmclst) ! 2022/10/5 yamaguchi

                     lcoll(22) = ncol22 + lcoll(22)
                     lcoll(23) = ncol23 + lcoll(23)
                     lcoll(24) = ncol24 + lcoll(24)

            end if

*-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      do i = massba + 1,massal
            ipot(i) = 1
      end do
!-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine resmass(ipw,isd1,isd2,srt,dem1,dem2,gam0,fqrq)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calcurate the mass and width of Delta and N*         *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              ipw = 0     : only width                                *
*              ipw = 1     : determine the mass                        *
*                                                                      *
*              isd1 = 2    : for Delta                                 *
*              isd1 = 3    : for N*(1440)                              *
*                                                                      *
*              isd2 = 1    : mass2 = dem2                              *
*              isd2 = 2    : for Delta                                 *
*              isd2 = 3    : for N*(1440)                              *
*                                                                      *
*              srt         : sqrt(s)                                   *
*                                                                      *
*              dem1        : minimum mass of resonance 1               *
*                            mass of resonance                         *
*                                                                      *
*              dem2        : minimum mass of resonance 2               *
*                            mass of resonance                         *
*                                                                      *
*              gam0        : width                                     *
*                                                                      *
*              fqrq        : phase space factor for absorption         *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param-incelf.inc'
*-----------------------------------------------------------------------

      dimension      rmss(3), gamr(3),
     &               bet2(3), qqr2(3)

      dimension      isdm(3)

*-----------------------------------------------------------------------

      data           rmss / 0.938, 1.232, 1.440 /
      data           gamr / 0.000, 0.110, 0.200 /
      data           bet2 / 1.000, 0.090, 0.274 /
      data           qqr2 / 1.000, 0.051936, 0.158125 /
      data           isdm / 0, 1, 1 /

*-----------------------------------------------------------------------

      CALL genseed(iseed)

            gam0   = 0.0
            fqrq   = 0.0

            rm2    = ( rmass/1000.0d0 )**2
            pm2    = pmass**2

*-----------------------------------------------------------------------
*     calculate only the width
*-----------------------------------------------------------------------

         if( ipw .eq. 0 ) then

            qq2  = 0.25 * ( ( dem1**2 - rm2 + pm2 ) / dem1 )**2 - pm2

            if( qq2 .le. 0.0 ) return

            form  = ( 1. + qqr2(isd1) / bet2(isd1) )
     &            / ( 1. + qq2        / bet2(isd1) )

            gam0  =   sqrt( qq2 / qqr2(isd1) )**3
     &              * rmss(isd1) / dem1 * gamr(isd1) * form**2

            fqrq  = qqr2(isd1) / qq2

            return

         end if


*-----------------------------------------------------------------------
*     maximum phase space factor
*-----------------------------------------------------------------------

            dem10 = dem1
            dem20 = dem2

            p02 = ( srt - dem10 - dem20 )
     &          * ( srt + dem10 + dem20 )
     &          * ( srt - dem10 + dem20 )
     &          * ( srt + dem10 - dem20 )

            if( p02 .le. 0.0 ) p02 = 0.0001

*-----------------------------------------------------------------------
*     determine the mass
*-----------------------------------------------------------------------

            idem = 0

  200       idem = idem + 1

            if( idem .gt. 500 )  return


            dem1 = dem10 + ( srt - dem10 - dem20 ) * rn(0)
            dem2 = dem20 + ( srt - dem1  - dem20 ) * rn(0)
     &                                             * isdm(isd2)


*-----------------------------------------------------------------------
*     check of phase space factor
*-----------------------------------------------------------------------

            pr2 = ( srt - dem1 - dem2 )
     &          * ( srt + dem1 + dem2 )
     &          * ( srt - dem1 + dem2 )
     &          * ( srt + dem1 - dem2 )

            gmcn =  pr2 / p02

            if( rn(0) .gt. gmcn ) goto 200

*-----------------------------------------------------------------------
*     calculate the width of Delta or N* by Moniz
*     check whether this mass satisfies the bright wigner distribution
*-----------------------------------------------------------------------

               demax1 = srt - dem20

               qq2  = 0.25 * ( ( dem1**2 - rm2 + pm2 ) / dem1 )**2 - pm2

               if( qq2 .le. 0.0 ) return

               form  = ( 1. + qqr2(isd1) / bet2(isd1) )
     &               / ( 1. + qq2        / bet2(isd1) )

               gam2  = ( sqrt( qq2 / qqr2(isd1) )**3
     &                 * rmss(isd1) / dem1 * gamr(isd1) * form**2 )**2


               bwtop = 1.0

            if( demax1 .lt. rmss(isd1) ) then

               bwtop = 0.25 * gam2 
     &               / ( ( rmss(isd1) - demax1 )**2 + 0.25 * gam2 )

            end if

               fmcn  = 0.25 * gam2 
     &               / ( ( rmss(isd1) - dem1   )**2 + 0.25 * gam2 )

            if( rn(0) * bwtop .gt. fmcn )  goto 200

*-----------------------------------------------------------------------
*     for second resonances
*-----------------------------------------------------------------------

         if( isd2 .ge. 1 ) then

               demax2 = srt - dem1

               qq2  = 0.25 * ( ( dem2**2 - rm2 + pm2 ) / dem2 )**2 - pm2

               if( qq2 .le. 0.0 ) return

               form  = ( 1. + qqr2(isd2) / bet2(isd2) )
     &               / ( 1. + qq2        / bet2(isd2) )

               gam2  = ( sqrt( qq2 / qqr2(isd2) )**3
     &                 * rmss(isd2) / dem2 * gamr(isd2) * form**2 )**2

               bwtop = 1.0

            if( demax2 .lt. rmss(isd2) ) then

               bwtop = 0.25 * gam2 
     &               / ( ( rmss(isd2) - demax2 )**2 + 0.25 * gam2 )

            end if

               fmcn  = 0.25 * gam2 
     &               / ( ( rmss(isd2) - dem2   )**2 + 0.25 * gam2 )

            if( rn(0) * bwtop .gt. fmcn ) goto 200

         end if

*-----------------------------------------------------------------------

      return
      end

!***********************************************************************
!*                                                                     *
      subroutine rr(i1,i2,r,p,rr2)
!*                                                                     *
!***********************************************************************

      implicit real*8(a-h,o-z)
      INCLUDE 'param-incelf.inc'

      real*8 r(3,nnn),p(6,nnn),rr2(nnn,nnn)

            rbrb = 0.0
            bij2 = 0.0
            rij2 = 0.0
            pij2 = 0.0

            eij  = p(4,i1) + p(4,i2)

            do l = 1, 3
               rij  =   r(l,i1) - r(l,i2)
               pij  =   p(l,i1) - p(l,i2)
               bij  = ( p(l,i1) + p(l,i2) ) /  eij
               rbrb = rbrb + rij * bij
               bij2 = bij2 + bij * bij
               rij2 = rij2 + rij * rij
               pij2 = pij2 + pij * pij
            end do

            rbrb      = 1.0 * rbrb
            gij2      = 1. / ( 1. - bij2 )

            rr2(i1,i2)  = rij2 + gij2 * rbrb * rbrb
            rr2(i2,i1)  = rr2(i1,i2)

      end

!***********************************************************************
!*                                                                     *
      subroutine knockout(i1,i2,ichannel,inds,ipot,ichg,r,p,rt00,upot,
     &                    nchgta,massba,iclst,icoales)
!*                                                                     *
!*    start knock out process to make cluster                          *
!*                                                                     *
!***********************************************************************
      implicit real*8(a-h,o-z)
      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      integer  i1,i2,ichg(nnn),inds(nnn),ipot(nnn),iclst(nnn)
      integer  icoales(nnn)
      integer  ia2,ia3,ia4,knoc,ichannel,nchgta,massba
      real*8   r(3,nnn),p(6,nnn),upot,rt00
      real*8   ea2,ea3,ea4

      integer inuc(nnn) !sawada8/27
      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )
      ucpot  = upot   + coul(nchgta,rt00)
      ducpot = dupot  + coul(nchgta,rt00)
      tucpot = thupot + coul(nchgta,rt00)
      hucpot = thupot + 2.0d0*coul(nchgta,rt00)
      aucpot = aupot  + 2.0d0*coul(nchgta,rt00)

!-----------------------------------------------------------------------


!=======================================================================
!           start knockout process
!           create a cluster
!=======================================================================

!------------------------------------------------------------12/9 sawada
             ia2 = 0
             ia3 = 0
             ia4 = 0
            knoc = 0
             ea2 = 0.0d0
             ea3 = 0.0d0
             ea4 = 0.0d0
!-----------------------------------------------------------------------
!*2010/11/4 change knockout parameter                                  *
!*       rd  :: knockout parameter of deuteron                         *
!*       rth :: knockout parameter of triton or He-3                   *
!*       ra  :: knockout parameter of alpha                            *

               rd  = 0.8d0
               rth = 1.0d0
               ra  = 1.0d0

!               rd  = 1.2d0
!               rth = 1.3d0
!               ra  = 1.4d0!*********sawada

!*2012/10/22 change knockout parameter   !nogamine                     *
!               rd  = 1.2d0
!               rth = 1.5d0
!               ra  = 1.7d0

!----------------------------------------------------------------------

             if( iclst(i2) .ne. 0 )return
             if( inds(i2) .ne. 1 )return
             if( ichannel .eq. 0 )return
             if( ichannel .eq. 99 )return
             if( icoales(i2) .ne. 0 )return
!alpha
             do j = 1, massba
               if( j .eq. i2 ) cycle
               if( j .eq. i1 ) cycle
               if( ipot(j) .eq. 1 ) cycle
               if( iclst(j) .ne. 0 ) cycle
               if( p(5,j) .lt. 0.10d0 ) cycle
               if( inds(j) .ne. 1 ) cycle
!-------------------------------------------------------------def2010/11
               rsq = ( r(1,i2) - r(1,j) )**2
     &              + ( r(2,i2) - r(2,j) )**2
     &              + ( r(3,i2) - r(3,j) )**2
!-----------------------------------------------------------------------
               rij = sqrt( rsq )

               if( rij .lt. ra ) then
                  if( ia2 .eq. 0 ) ia2 = j
                  if( ia2 .gt. 0 .and. j .ne. ia2 .and.
     &                 ia3 .eq. 0 ) ia3 = j
                  if( ia2 .gt. 0 .and. ia3 .gt. 0 .and.
     &                 j .ne. ia2 .and. j .ne. ia3 .and.
     &                 ia4 .eq. 0 ) ia4 = j
               end if
              if( ia4 .gt. 0 ) knoc = 3
             end do

!3He & triton
            if( ia4 .eq. 0 )then
             ia2 = 0
             ia3 = 0
             ia4 = 0

             do j = 1, massba
               if( j .eq. i2 ) cycle
               if( j .eq. i1 ) cycle
               if( ipot(j) .eq. 1 ) cycle
               if( iclst(j) .ne. 0 ) cycle
               if( inds(j) .ne. 1 ) cycle
               if( p(5,j) .lt. 0.10d0 ) cycle
!-------------------------------------------------------------def2010/11
               rsq = ( r(1,i2) - r(1,j) )**2
     &              + ( r(2,i2) - r(2,j) )**2
     &              + ( r(3,i2) - r(3,j) )**2
!-----------------------------------------------------------------------
               rij = sqrt( rsq )

               if( rij .lt. rth ) then
                  if( ia2 .eq. 0 ) ia2 = j
                  if( ia2 .gt. 0 .and. j .ne. ia2 .and.
     &                 ia3 .eq. 0 ) ia3 = j
                  if( ia2 .gt. 0 .and. ia3 .gt. 0 .and.
     &                 j .ne. ia2 .and. j .ne. ia3 .and.
     &                 ia4 .eq. 0 ) ia4 = j
               end if
             end do
              if( ia3 .gt. 0 ) knoc = 2
            end if

!deuteron
            if( ia4 .eq. 0 .and. ia3 .eq. 0 )then
             ia2 = 0
             ia3 = 0
             ia4 = 0

             do j = 1, massba
               if( j .eq. i2 ) cycle
               if( j .eq. i1 ) cycle
               if( ipot(j) .eq. 1 ) cycle
               if( iclst(j) .ne. 0 ) cycle
               if( inds(j) .ne. 1 ) cycle
               if( p(5,j) .lt. 0.10d0 ) cycle
!-------------------------------------------------------------def2010/11
               rsq = ( r(1,i2) - r(1,j) )**2
     &              + ( r(2,i2) - r(2,j) )**2
     &              + ( r(3,i2) - r(3,j) )**2
!-----------------------------------------------------------------------
               rij = sqrt( rsq )

               if( rij .lt. rd ) then
                  knoc = 1
                  if( ia2 .eq. 0 ) ia2 = j
                  if( ia2 .gt. 0 .and. j .ne. ia2 .and.
     &                 ia3 .eq. 0 ) ia3 = j
                  if( ia2 .gt. 0 .and. ia3 .gt. 0 .and.
     &                 j .ne. ia2 .and. j .ne. ia3 .and.
     &                 ia4 .eq. 0 ) ia4 = j
               end if
             end do
              if( ia2 .gt. 0 ) knoc = 1
            end if

!-----------------------------------------------------------------------
!     Helium 4 production
!-----------------------------------------------------------------------

             if( ia2 .gt. 0 .and. ia3 .gt. 0 .and. ia4 .gt. 0 ) then

               if( iclst(i2) .eq. 0 .and. iclst(ia2) .eq. 0 .and.
     &              iclst(ia3) .eq. 0 .and. iclst(ia4) .eq. 0 ) then
                  
                  iclchg = ichg(i2) + ichg(ia2)
     &                 + ichg(ia3) + ichg(ia4)

                  if( iclchg .eq. 2 ) then
                     iclst(i2)  = 4
                     iclst(ia2) = -1
                     iclst(ia3) = -1
                     iclst(ia4) = -1
                     inds(i2)    = 1
                     inds(ia2)   = 1
                     inds(ia3)   = 1
                     inds(ia4)   = 1
                     ipot(ia2)   = 0
                     ipot(ia3)   = 0
                     ipot(ia4)   = 0
                     ea2 = sqrt( p(1,ia2)**2 + p(2,ia2)**2
     &                    + p(3,ia2)**2 + rmass**2 )
                     ea3 = sqrt( p(1,ia3)**2 + p(2,ia3)**2
     &                    + p(3,ia3)**2 + rmass**2 )
                     ea4 = sqrt( p(1,ia4)**2 + p(2,ia4)**2
     &                    + p(3,ia4)**2 + rmass**2 )
                     p(4,i2) = dsqrt(p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+
     &                                    rmass**2) +ea2 +ea3 +ea4
!     &                                    + 28.29 !Binding E = 28.29[MeV]
                     p(5,i2) = amass

                     et2 = p(4,i2) - amass
                     pabs = dsqrt( et2*( et2 + 2.0d0*amass ) )

                     xa = p(1,i2) + p(1,ia2) + p(1,ia3) + p(1,ia4)
                     ya = p(2,i2) + p(2,ia2) + p(2,ia3) + p(2,ia4)
                     za = p(3,i2) + p(3,ia2) + p(3,ia3) + p(3,ia4)

                     pclst0 = sqrt( xa**2 + ya**2 + za**2 )

                     p(1,i2)  = xa / pclst0 * pabs
                     p(2,i2)  = ya / pclst0 * pabs
                     p(3,i2)  = za / pclst0 * pabs

                     do k = 1,3
                      p(k,ia2) = 0.0d0
                      p(k,ia3) = 0.0d0
                      p(k,ia4) = 0.0d0
                     end do

                     p(4,ia2) = 0.0d0
                     p(4,ia3) = 0.0d0
                     p(4,ia4) = 0.0d0
                     p(5,ia2) = 0.0d0
                     p(5,ia3) = 0.0d0
                     p(5,ia4) = 0.0d0
                     if( et2 .gt. aucpot )then
                      ipot(i2) = 1
                     else
                      ipot(i2) = 0
                     end if

                  else
                     iclst(i2)  = 0
                     iclst(ia2) = 0
                     iclst(ia3) = 0
                     iclst(ia4) = 0
                  end if
               end if

!-----------------------------------------------------------------------
!     Helium 3 production
!-----------------------------------------------------------------------

             elseif( ia2 .gt. 0 .and. ia3 .gt. 0 .and. ia4 .eq. 0 ) then

               if( iclst(i2) .eq. 0 .and. iclst(ia2) .eq. 0 .and.
     &              iclst(ia3) .eq. 0 ) then

                  iclchg = ichg(i2) + ichg(ia2) + ichg(ia3)

                  if( iclchg .eq. 2 ) then
                     iclst(i2)  = 3
                     iclst(ia2) = -1
                     iclst(ia3) = -1
                     inds(i2)    = 1
                     inds(ia2)   = 1
                     inds(ia3)   = 1
                     ipot(ia2)   = 0
                     ipot(ia3)   = 0

                     ea2 = sqrt( p(1,ia2)**2 + p(2,ia2)**2
     &                    + p(3,ia2)**2 + rmass**2 )
                     ea3 = sqrt( p(1,ia3)**2 + p(2,ia3)**2
     &                    + p(3,ia3)**2 + rmass**2 )
                     p(4,i2) = dsqrt(p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+
     &                                    rmass**2) + ea2 + ea3
!     &                                    + 7.7181d0 !Binding E
                     p(5,i2) = hmass

                     et2 = p(4,i2) - hmass
                     pabs = dsqrt( et2*(et2 + 2.0d0*hmass ) )

                     xh = p(1,i2) + p(1,ia2) + p(1,ia3)
                     yh = p(2,i2) + p(2,ia2) + p(2,ia3)
                     zh = p(3,i2) + p(3,ia2) + p(3,ia3)

                     pclst0 = sqrt( xh**2 + yh**2 + zh**2 )

                     p(1,i2)  = xh / pclst0 * pabs
                     p(2,i2)  = yh / pclst0 * pabs
                     p(3,i2)  = zh / pclst0 * pabs

                     do k = 1,3
                      p(k,ia2) = 0.0d0
                      p(k,ia3) = 0.0d0
                     end do

                     p(4,ia2) = 0.0d0
                     p(4,ia3) = 0.0d0
                     p(5,ia2) = 0.0d0
                     p(5,ia3) = 0.0d0
                     if( et2 .gt. hucpot )then
                      ipot(i2) = 1
                     else
                      ipot(i2) = 0
                     end if

!-----------------------------------------------------------------------
!     Triton production
!-----------------------------------------------------------------------

                  elseif( iclchg .eq. 1 ) then
                     iclst(i2)  = 2
                     iclst(ia2) = -1
                     iclst(ia3) = -1
                     inds(i2)    = 1
                     inds(ia2)   = 1
                     inds(ia3)   = 1
                     ipot(ia2)   = 0
                     ipot(ia3)   = 0
                     ea2 = sqrt( p(1,ia2)**2 + p(2,ia2)**2
     &                    + p(3,ia2)**2 + rmass**2 )
                     ea3 = sqrt( p(1,ia3)**2 + p(2,ia3)**2
     &                    + p(3,ia3)**2 + rmass**2 )
                     p(4,i2) = dsqrt(p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+
     &                                    rmass**2) +ea2 +ea3
!     &                                    + 8.4820d0 !Binding E
                     p(5,i2) = tmass

                     et2 = p(4,i2) - tmass
                     pabs = dsqrt( et2*(et2 + 2.0d0*tmass ) )

                     xt = p(1,i2) + p(1,ia2) + p(1,ia3)
                     yt = p(2,i2) + p(2,ia2) + p(2,ia3)
                     zt = p(3,i2) + p(3,ia2) + p(3,ia3)

                     pclst0 = sqrt( xt**2 + yt**2 + zt**2 )

                     p(1,i2)  = xt / pclst0 * pabs
                     p(2,i2)  = yt / pclst0 * pabs
                     p(3,i2)  = zt / pclst0 * pabs

                     do k = 1,3
                      p(k,ia2) = 0.0d0
                      p(k,ia3) = 0.0d0
                     end do
                     p(4,ia2) = 0.0d0
                     p(4,ia3) = 0.0d0
                     p(5,ia2) = 0.0d0
                     p(5,ia3) = 0.0d0
                     if( et2 .gt. tucpot )then
                      ipot(i2) = 1
                     else
                      ipot(i2) = 0
                     end if

                  else
                     iclst(i2)  = 0
                     iclst(ia2) = 0
                     iclst(ia3) = 0
                  end if

               end if

!-----------------------------------------------------------------------
!     deuteron production
!-----------------------------------------------------------------------

             elseif( ia2 .gt. 0 .and. ia3 .eq. 0 .and. ia4 .eq. 0 ) then

               if( iclst(i2) .eq. 0 .and. iclst(ia2) .eq. 0 ) then

                  iclchg = ichg(i2) + ichg(ia2)

                  if( iclchg .eq. 1 ) then

                     iclst(i2)  = 1
                     iclst(ia2) = -1
                     inds(i2)    = 1
                     inds(ia2)   = 1
                     ipot(ia2)   = 0
                     ea2 = sqrt( p(1,ia2)**2 + p(2,ia2)**2
     &                    + p(3,ia2)**2 + rmass**2 )

                     p(4,i2) = dsqrt(p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+
     &                                    rmass**2) + ea2
!     &                                    + 2.2245d0 !Binding E
                     p(5,i2) = ddmass

                     et2 = p(4,i2) - ddmass
                     pabs = dsqrt( et2*(et2 + 2.0d0*ddmass ) )

                     xd = p(1,i2) + p(1,ia2)
                     yd = p(2,i2) + p(2,ia2)
                     zd = p(3,i2) + p(3,ia2)

                     pclst0 = sqrt( xd**2 + yd**2 + zd**2 )

                     p(1,i2)  = xd / pclst0 * pabs
                     p(2,i2)  = yd / pclst0 * pabs
                     p(3,i2)  = zd / pclst0 * pabs

                     do k = 1,3
                      p(k,ia2) = 0.0d0
                     end do

                     p(4,ia2) =0.0d0
                     p(5,ia2) =0.0d0
                     if( et2 .gt. ducpot )then
                      ipot(i2) = 1
                     else
                      ipot(i2) = 0
                     end if

                  else
                     iclst(i2)  = 0
                     iclst(ia2) = 0
                  end if

               end if

             else

               iclst(i2) = 0

             end if
!-------------------------------------------------------------end sawada

      return
      end


!***********************************************************************
!*                                                                     *
!     SUBROUTINE d_ipot(p,rt00,upot,nchgta,massal,ichg,inds,iclst,ipot)
!      SUBROUTINE d_ipot(p,rt00,upot,nchgta,massal,ichg,inds,iclst,ipot,
!     &     igroup,nmascl,ipotcc,ccp)      
      SUBROUTINE d_ipot(p,rt00,upot,nchgta,massal,ichg,inds,iclst,ipot,
     &     numclst,nmasclst,pcmclst) ! 2022/10/5 yamaguchi     
!*                                                                     *
!*                                                                     *
!***********************************************************************

      implicit real*8(a-h,o-z)
      INCLUDE 'param-incelf.inc'

      integer  massal,nchgta,i1
      integer  ipot(nnn),ichg(nnn),iclst(nnn),inds(nnn)
      real*8   p(6,nnn),upot,rt00
!      integer, intent(in) :: igroup(nnn),nmascl(nnn)
      integer, intent(in) :: numclst(nnn),nmasclst(nnn)   ! 2022/10/5 yamaguchi
!      integer, intent(inout) :: ipotcc(nnn)
!      real*8, intent(in) :: ccp(6,nnn)
      real*8, intent(in) :: pcmclst(5,nnn)
!-----------------------------------------------------------------------
      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )
!-----------------------------------------------------------------------
      ucpot  = upot   + coul(nchgta,rt00)
      ducpot = dupot  + coul(nchgta,rt00)
      tucpot = thupot + coul(nchgta,rt00)
      hucpot = thupot + 2.0d0*coul(nchgta,rt00)
      aucpot = aupot  + 2.0d0*coul(nchgta,rt00)

      do i1 = 1, massal
!----------------------------------------------------------For cluster
         if( numclst(i1) .ne. 0 )then
            
            if( nmasclst(i1) .ne. 0 )then
               
               ccp02 = pcmclst(1,i1)**2 + pcmclst(2,i1)**2
     &              + pcmclst(3,i1)**2

               cct0 = dsqrt( ccp02 + pcmclst(5,i1)**2 ) - pcmclst(5,i1)
               
               if( cct0 .lt. upot*dble(nmasclst(i1)) )then
                  
                  ipot(i1) = 0

               else

                  ipot(i1) = 1

               end if

            end if
!---------------------------------------------------------------------
         else
            
       e1     = p(4,i1)
       em1    = p(5,i1)
       t1     = e1 - em1
!----------------------------------------------------------
!       if( p(5,i1) .lt. 0.10d0 )cycle     !sawada
!       if( iclst(i1) .eq. -1 )cycle     !sawada
!-----------------------------------------------------------
            if( inds(i1) .eq. 1 ) then
                  if( iclst(i1) .eq. 4 )then    !sawada
!                    if( t1 .gt. aucpot ) then
                    if( t1 .gt. aupot ) then
                        ipot(i1) = 1
                    else
                        ipot(i1) = 0
                    end if
                  elseif( iclst(i1) .eq. 3 )then
!                    if( t1 .gt. hucpot ) then
                    if( t1 .gt. thupot ) then
                        ipot(i1) = 1
                    else
                        ipot(i1) = 0
                    end if
                  elseif( iclst(i1) .eq. 2 )then
!                    if( t1 .gt. tucpot ) then
                    if( t1 .gt. thupot ) then
                        ipot(i1) = 1
                    else
                        ipot(i1) = 0
                    end if
                  elseif( iclst(i1) .eq. 1 )then
!                    if( t1 .gt. ducpot ) then
                    if( t1 .gt. dupot ) then
                        ipot(i1) = 1
                    else
                        ipot(i1) = 0
                    end if
                  elseif( iclst(i1) .eq. 0 )then
                    if( ichg(i1) .eq. 0 ) then
                      if( t1 .gt. upot ) then
                         ipot(i1) = 1
                      else
                         ipot(i1) = 0
                      end if
                    else
                      if( t1 .gt. ucpot ) then
                         ipot(i1) = 1
!2012/12/13nogamine                      elseif( t1 .lt. upot ) then
                      else
                         ipot(i1) = 0
                      end if
                    end if
                  elseif( iclst(i1) .eq. -1 )then
                        ipot(i1) = 0
                  end if
!-----------------------------------------------------------------sawada
            else
                ipot(i1) = 1              ! pion,N* and delta are moving
            end if
         end if

      end do
!-----------------------------------------------------------------------

      end SUBROUTINE d_ipot


!******************************************************2018/8/21watanabe
!*                                                                     *
!      SUBROUTINE collex( rt00,radm,nmasta,nchgta,ein,p,peex,gdrflg )
      SUBROUTINE collex(iproj,rt00,radm,nmasta,nchgta,ein,p,peex,
     &    ncolexb,colexb,ichg,nmaspr,ngtgr,q,eisq,bmax ) !fuku
!*                                                                     *
!*                                                                     *
!*        Purpose:                                                     *
!*                                                                     *
!*              to give collective excitation correction	       *
!*                                                                     *
!*                                            2013/7/17yamada	       *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      integer, intent(in)    :: nmasta, nchgta,iproj,nmaspr,ngtgr
      integer, intent(out)   :: ncolexb,colexb
      integer, intent(inout) :: ichg(nnn)
      real*8, intent(in)    :: ein,rt00,radm,q
      real*8, intent(out)   :: peex,eisq
      real*8, intent(inout) :: p(6,nnn)

!-----------------------------------------------------------------------

!      if( gdrflg .eq. -1 )then !watanabe2018/8/30

        ncolexb = 0
        colexb  = 0
	gdrflg = 0 !watanabe 2018/8/30
!         dein = ein - coul(nchgta,rt00)
        if(iproj .eq. 1)then ! p fuku
           dein = ein - coul(nchgta,rt00)
        elseif(iproj .eq. 0)then ! n 
           dein = ein
        else
           return
        endif
        exr = dein
        peex = 0.0d0
        
!--------------------------------2020/4/27ktym_pn_probability
        pgtgrr = rn(0)
!hisatomi
        if(iproj .eq. 1)then
           pgtgri = 6.0d-4      !IAS
           pgtgrj = 9.0d-5*(nmasta - ( 2.0d0 * nchgta)) + 4.1d-3 !L=0
           pgtgrk = 9.0d-4*(nmasta - ( 2.0d0 * nchgta)) + 4.12d-3 !L=1
           pgtgrl = 9.0d-4*(nmasta - ( 2.0d0 * nchgta)) + 3.12d-3 !L=2
!hisatomi
!           pgtgri = 0!6.0d-4      !IAS
!           pgtgrj = 0!9.0d-5*(nmasta - ( 2.0d0 * nchgta)) + 4.1d-3 !L=0
!           pgtgrk = 0!9.0d-4*(nmasta - ( 2.0d0 * nchgta)) + 4.12d-2 !L=1
!           pgtgrl = 0!9.0d-4*(nmasta - ( 2.0d0 * nchgta)) + 3.12d-2 !L=2


        else if(iproj .eq. 0)then
           pgtgri = 5.0d-2      !L=3
           pgtgrj = 7.0d-3      !L=0
           pgtgrk = 5.0d-2      !L=1
           pgtgrl = 5.0d-2      !L=2
        end if
        
!------------------------------------------------------------
        peliel = rn(0)

        if((ngtgr .eq. 0) .or. (ngtgr .eq. 1 .and.  
     &     pgtgrr .gt. pgtgri + pgtgrj + pgtgrk + pgtgrl))then

!------------------------------------------------------------
!-----------------------------------------------------------------------
!	Low excitation energy[MeV] at Breit-Wigner peak
!-----------------------------------------------------------------------
        expeak = 0.0433d0 * dble(nmasta)**0.8069d0

!-----------------------------------------------------------------------
!	Energy-dependent Low excitation cross-section[mb]
!-----------------------------------------------------------------------

        peakcs = -0.3693 * nmasta + 78.773
        wp = 100.			!width parameter[MeV]

        if( nmasta .le. 63 )then
		egr1 = 0.8273 * nmasta - 23.849
        elseif( nmasta .le. 200 )then
		egr1 = 30.
        else
		egr1 = 2.2727 * nmasta - 420.
        end if
        csmax = wp**2.
     &         / ( wp**2. / 4. )
        csstd = wp**2.
     &         / ( ( egr1 - ein )**2. + wp**2 / 4. ) / csmax

        csein = csstd * peakcs

!-----------------------------------------------------------------------
!	Calculation of Breit-Wigner area[mb] 
!	centering on energy-dependent peak.
!-----------------------------------------------------------------------

	rwp = 3.				!width parameter [MeV]
	unibw = 0.001
	bw = 0.
	crax = 0.

        do i = 1, 30000
          bw = bw + unibw
          bwmax = rwp**2.
     &             / ( rwp**2. / 4. )! / 4.
	  bwstd = rwp**2.
     &             / ( ( expeak - bw)**2. + rwp**2 / 4. ) / bwmax! / 4.

	  crax = crax + bwstd * csein * unibw

	end do

!-----------------------------------------------------------------------
!	Magnification of the area[mb] 
!	based on ratio of collision and passing through.
!-----------------------------------------------------------------------

	    crsmg = -( 0.0008d0 * exp( 0.0091d0 * dble(nmasta) ) ) * ein
     &            + ( 0.0067d0 * dble(nmasta) + 1.5724d0 )

	     dpcr = crax * crsmg

!-----------------------------------------------------------------------

	     tcrs = pi * radm**2 *10. !total cross section[mb]

             if( iproj .eq. 1 )then ! fuku
                dppa = dpcr / tcrs !rate of low vibration particle
             else
                dppa = dpcr / tcrs + 0.1d0 !rate of low vibration particle
             end if
             
	     grpa = 0.065d0		!rate of giant resonance particle
	     gdrr = ((60. * nchgta * ( nmasta - nchgta ) / nmasta)
     &             / ( pi * 4.  )) / (10. * pi * bmax**2.)     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!hisatomi  
             if (ein .ge. 80 .and. ein .le.150) then
             grpa = 0.015d0
             else if (ein .gt. 150) then 
             grpa = 0.0001
             end if 
           colexb = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!hisatomi

	   peliel = rn(0)

	if( peliel .gt. dppa + grpa + gdrr ) then

	   peex = 0.0

	else if( peliel .gt. grpa + gdrr ) then

!-----------------------------------------------------------------------
!	Determination of low excitation energy[MeV]
!-----------------------------------------------------------------------

        colexb = 1          !for collex_bending     
        rwp = 3.			!width parameter [MeV]
        eex = rn(0) * exr
        pex0 = rn(0)
        pex1 = -10.0

        do while ( pex0 .gt. pex1 )
          eex = rn(0) * exr
          pex0 = rn(0)
          pex1 = rwp**2.
     &            / ( ( expeak - eex )**2. + rwp**2. / 4. ) / 4.

        end do
      
        peex = eex

	 else if( peliel .gt. gdrr ) then

!-----------------------------------------------------------------------
!	Determination of high excitation energy[MeV]
!	arise from giant resonance
!-----------------------------------------------------------------------

           colexb = 1          !for collex_bending       
	   rwp = 8.			!width parameter [MeV]
	   eisq = 65.0d0 * nmasta**( -1. / 3. )		!iso scalar quadrupole

	   eex = rn(0) * exr
	   pgr = rn(0)
	   pex0 = rn(0)
           pex1 = -10.0

	   do while ( pex0 .gt. pex1 )

	      eex = rn(0) * exr
	      pex0 = rn(0)
	      pex1 = rwp**2.
     &             / ( ( eisq - eex )**2. + rwp**2. / 4. ) / 4.
	   end do

	   peex = eex

	else
!---------------------------------------------------------------------
!GDR hisatomi
!---------------------------------------------------------------------
           colexb = 0          !for collex_bending
           rwp = 8.             !width parameter [MeV] INC-LOW
           rgdr = 1.18d0 * nmasta**( 1. / 3. ) !

           eex = 0
           pgr = rn(0)
           pex0 = 10.
           pex1 = -10.0

           do while ( pex0 .gt. pex1 )

              eex = rn(0) * exr
!              sigm = (60 * nchgta * ( nmasta - nchgta ) / nmasta)
!     &             / ( pi * rwp /2 )
              ugdr = (( 3 * 36.8 ) / 17 ) * nmasta **( -1. / 3. )  
              egdre = 0.0768
              egdra = 0.7 * ( rgdr**2. ) / ( 8 * 36.8 )
              egdrb = 1. + ugdr - (1. + egdre + 3.*ugdr)*egdre
     &              / (1. + egdre + ugdr)

              egdr = 6.58212 / (sqrt(egdra * egdrb))
        
!              egdr = 11.0d0
              pex0 = rn(0)
              pex1 = 1
     &             / ( 1 + (( eex**2 - egdr**2 ) 
     &             / (eex * rwp))**2 )

!                 pex1 = rwp**2
!    &             / ( ( eisq - eex )**2 + rwp**2 / 4. ) / 4.
           end do
           peex = eex
        end if
!hisatomi

!---------------------------------------------------------------------
!kyodaikyoumei happen 2021/1 katayama
!---------------------------------------------------------------------
! IAS  & L=3 no reiki
!---------------------------------
        else if( ngtgr.eq.1 .and. pgtgrr.gt.pgtgrj+pgtgrk+pgtgrl) then

           if(iproj .eq. 1)then
              ncolexb = 4       !IAS no reiki
              rwp  = 0.01       !width parameter [MeV]
              eisq = -6.8d-4*ein*(dble(nmasta)-(2.0*dble(nchgta)))
     &             + 0.3826*(dble(nmasta)-(2.0*dble(nchgta))) + 4.0283
           else if(iproj .eq. 0)then
              ncolexb = 3       !IAS no reiki
              rwp  = 20.0       !width parameter [MeV]
              eisq = 30.0
           end if

           eex = 0
           pgr = rn(0)
           pex0 = 10.
           pex1 = -10.0
           do while ( pex0 .gt. pex1 )
              eex = rn(0) * exr
              pex0 = rn(0)
              pex1 = rwp**2
     &             / ( ( eisq - eex )**2 + rwp**2 / 4. ) / 4.
           end do
           peex = eex
!------------------------------
!          L=0 no reiki
!------------------------------
        else if( ngtgr.eq.1 .and. pgtgrr .gt. pgtgrk + pgtgrl) then
           ncolexb = 8          !L=0 no reiki
           if(iproj .eq. 1)then
              rwp  = 4.2        !width parameter [MeV]
              eisq = -6.0d-4*ein*(dble(nmasta)-(2.0*dble(nchgta)))
     &             + 0.2213*(dble(nmasta)-(2.0*dble(nchgta))) + 11.128
           else if(iproj .eq. 0)then
              rwp  = 20.0       !width parameter [MeV]
              eisq = 18.0
           end if

              eex = 0
              pgr = rn(0)
              pex0 = 10.
              pex1 = -10.0
           do while ( pex0 .gt. pex1 )
              eex = rn(0) * exr
              pex0 = rn(0)
              pex1 = rwp**2
     &                / ( ( eisq - eex )**2 + rwp**2 / 4. ) / 4.
           end do
           peex = eex
!--------------------------------
!     L=1 no reiki
!--------------------------------
        else if( ngtgr .eq. 1 .and. pgtgrr .gt. pgtgrl ) then
           ncolexb = 1          !L=1 no reiki
           if(iproj .eq. 1)then
              rwp  = 10.0       !width parameter [MeV]
              eisq = -6.0d-4*ein*(dble(nmasta)-(2.0*dble(nchgta)))
     &            + 0.1038*(dble(nmasta)-(2.0*dble(nchgta))) + 23.858
           else if(iproj .eq. 0)then
              rwp  = 20.0       !width parameter [MeV]
              eisq = 18.0
           end if

           eex = 0
           pgr = rn(0)
           pex0 = 10.
           pex1 = -10.0
           do while ( pex0 .gt. pex1 )
              eex = rn(0) * exr
              pex0 = rn(0)
              if(eex .gt. eisq )then !2019/5/10/katayama/renzoku
                 pex1 = (eisq + 20.0)**2
     &                / (( eisq - eex )**2 + (eisq + 20.0 )**2)
              else
                 pex1 = rwp**2
     &                / ( ( eisq - eex )**2 + rwp**2 / 4. ) / 4.
              end if
           end do
           peex = eex
!--------------------------------
!          L=2 no reiki
!--------------------------------
       else
          ncolexb = 2           !L=2 no reiki
          if(iproj .eq. 1)then
             rwp  = 14.0        !width parameter [MeV]
             eisq = -3.0d-4*ein*(dble(nmasta)-(2.0*dble(nchgta)))
     &           + 0.074*(dble(nmasta)-(2.0*dble(nchgta))) + 29.56
          else if(iproj .eq. 0)then
             rwp  = 20.0        !width parameter [MeV]
             eisq = 22.4
          end if
           eex = 0
           pgr = rn(0)
           pex0 = 10.
           pex1 = -10.0
           do while ( pex0 .gt. pex1 )
              eex = rn(0) * exr
              pex0 = rn(0)
              if(eex .gt. eisq )then !2019/5/10/katayama/renzoku
                 pex1 = (eisq + 20.0)**2
     &                / (( eisq - eex )**2 + (eisq + 20.0 )**2)
              else
                 pex1 = rwp**2
     &                / ( ( eisq - eex )**2 + rwp**2 / 4. ) / 4.
              end if
           end do
           peex = eex
        end if
!     print'(I3)',ncolexb
!-----------------------------------------------------------------------
        if( ngtgr .eq. 1 .and. ncolexb .ne. 0 )then
           if(iproj .eq. 1)then
              ichg(1) = 0       !iproj ga proton
              ichg(nmaspr + nchgta + 1) = 1 !iproj ga proton
           else if(iproj .eq. 0)then
              ichg(1) = 1       !iproj ga neutron
              ichg(nmaspr + nchgta + 1) = 0 !iproj ga neutron
           end if
        end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

      tout = p(4,1) - p(5,1)
      pout = dsqrt( tout * ( tout + 2.0d0 * p(5,1) ) )

      tinn = p(4,1) - p(5,1) - peex
      pinn = dsqrt( tinn * ( tinn + 2.0d0 * p(5,1) ) )

      do k = 1, 3
         p(k,1) = p(k,1) * pinn / pout
      end do

      p1tot = dsqrt( p(1,1)**2 + p(2,1)**2 + p(3,1)**2 )
      p(4,1) = dsqrt( p1tot**2 + p(5,1)**2 )


!      end if

	END SUBROUTINE collex


!***********************************************************************
!*                                                                     *
      SUBROUTINE deflection(nmasta,nmaspr,ein,upot,ucpot,r,p,
     &	iclst,inds,nbend,iproj,ichg,itime,ipot,ncolexb,colexb)
!     SUBROUTINE deflection( nmasta,nmaspr,ein,upot,ucpot,r,p,iclst,
!    &	inds,nbend,iproj,ichg,itime,ipot )
!*                                                                     *
!*                                                                     *
!*        Purpose:                                                     *
!*                                                                     *
!*              to calculate the angle of orientation bent by 	     *
!*		    nuclear potential             				     *
!*                                                                     *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      integer, intent(in)    :: nmasta,nmaspr,iclst(nnn),inds(nnn),
     &  nbend,iproj,ichg(nnn),ncolexb,colexb
      real*8, intent(in)    :: upot,ucpot,ein,r(3,nnn)
      real*8, intent(inout) :: p(6,nnn)
      real*8 ke
      integer,intent(inout) :: ipot(nnn)
!-----------------------------------------------------------------------
      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )
!-----------------------------------------------------------------------
! factor for ke cal.
! p(4,i) is total E inside nucleus when 1st deflection and outside nucleus
! when 2nd deflection.
!      fac = dble(abs(nbend - 2))
      !fac = 1.0d0  ! 2022/10/5 yamaguchi
      fac = 1
!-----------------------------------------------------------------------

	do i = 1 , nmasta+nmaspr
      
          pp0 = 0.0d0
!          psi = 0.0d0
!          chi = 0.0d0
          phi = 0.0d0
          the = 0.0d0

         ome1 = 0.0d0
         ome2 = 0.0d0
         ome3 = 0.0d0

        dome1 = 0.0d0
        dome2 = 0.0d0
        dome3 = 0.0d0

          opt = 0.0d0
         fthe = 0.0d0
           ke = 0.0d0
           h = 0.0d0
!           if( iproj.eq.1 .and. i.eq.1 .and. ichg(1).eq.0 ) cycle !wata

!*****************************************watanabe 8/27
          if( p(5,i) .lt. 0.1d0 ) cycle
	  if( inds(i) .eq. 4 ) cycle
!          if( idpid(i) .eq. 1 ) cycle
!	  if( iclst(i) .ne. 0 ) cycle
! pion ha Deflection nasi.
!*******************************************************

!          if( iproj .eq. 1 .and. nbend .eq. 2 )then! for kalbach
!             if( iclst(i) .eq. 1 .and. i .eq. 1 )cycle !t.mori
!          end if
      
           ke = p(4,i) - p(5,i) - fac*upot

          if( iclst(i) .eq. 1 )then
           ke = p(4,i) - fac*dupot
	  else if( iclst(i) .eq. 2 .or. iclst(i) .eq. 3 )then
           ke = p(4,i) - hmass - fac*thupot
	  else if( iclst(i) .eq. 4 )then
           ke = p(4,i) - amass - fac*aupot
	  end if

         if( ke .le. 0.0d0 ) cycle

          pp0 = dsqrt( p(1,i)**2. + p(2,i)**2. + p(3,i)**2. )

         ome1 = p(1,i) / pp0
         ome2 = p(2,i) / pp0
         ome3 = p(3,i) / pp0

          opt = rn(0)

        if( iproj .eq. 5) then

            h = -0.006d0 * dlog(dble(nmasta)) + 0.005d0
            g = -0.0013 * ke + h
          if( iclst(i) .eq. 1 )then
            h =  0.001d0 * dlog(dble(nmasta)) + 0.02! - 0.007d0
            g = -0.00005 * ke + h
          else if( iclst(i) .eq. 2 )then
            h = -0.006d0 * dlog(dble(nmasta)) + 0.1!0.005d0
            g = -0.0006 * ke + h
          else if( iclst(i) .eq. 3 )then
            h = -0.006d0 * dlog(dble(nmasta)) + 0.05!0.005d0
            g = -0.0006 * ke + h
          else if( iclst(i) .eq. 4 )then
            h = 0.01d0 * dlog(dble(nmasta)) - 0.04d0
            g = -0.0012 * ke *(ke/400.)+ h
          end if
          
        else
            h = -0.006d0 * dlog(dble(nmasta)) + 0.005d0
            g = -0.0013 * ke + h

          if( iclst(i) .eq. 1 )then
           h =  0.001d0 * dlog(dble(nmasta)) - 0.007d0
           g = -0.001 * ke + h
          else if( iclst(i) .eq. 2 )then
            h = -0.006d0 * dlog(dble(nmasta)) + 0.005d0
            g = -0.0012 * ke + h
          else if( iclst(i) .eq. 3 )then
            h = -0.006d0 * dlog(dble(nmasta)) + 0.005d0
            g = -0.0012 * ke + h
          else if( iclst(i) .eq. 4 )then
            h = 0.01d0 * dlog(dble(nmasta)) - 0.04d0
            g = -0.0012 * ke + h
          end if
          
          if(nbend.eq.2)then
          if(iclst(i).eq.1)then 
          
               theta=0
             if( sqrt( p(1,i)**2. + p(2,i)**2. + p(3,i)**2. ).eq.0 )then
                    theta = 0.0
             else
                     theta = 57.3*acos( p(3,i)      
     &                   / sqrt( p(1,i)**2. + p(2,i)**2. + p(3,i)**2. ))
             endif
           endif
           endif
           endif

	   if( nbend .eq. 1 )then

!-----------------------------------------------------------------------
!	Determination of incident energy-dependent
!	angle[deg] range of deflection.
!-----------------------------------------------------------------------

              eintot = ein!dble(nmaspr)*ein !t.mori0601
              
	      prdef = 0.1818d0 * exp( 0.056d0 * eintot )
!	      prdef = 0.2214d0 * exp( 0.0553d0 * ein )

	      agpe = datan( -pi / 180.0d0 / ( prdef * g ) ) /pi * 180.0d0	![deg]

           fthemax = dsin( agpe / 180.0d0 * pi)
     &    		 * exp( prdef * g * agpe )
!     &                    +exp(0.1*agpe))!202440708

	       agr = log( 0.05d0 * fthemax ) / prdef / g

	       if( agr .gt. 180.0d0 )then
		    agr = 180.0d0
	       end if

!-----------------------------------------------------------------------
!	Determination of angle[deg] of deflection
!	to nucleous inward direction
!	based on elastic scattering angular distribution.
!-----------------------------------------------------------------------

            fthe = -10.0

            do while ( opt .ge. fthe )

               the =  agr * rn(0)

               opt = rn(0)

!               fthe = dsin(the / 180.0d0 * pi) *				![rad]
!     &                exp( 0.1818d0 * exp( 0.056d0 * ein ) *
!     &              ( -0.0013d0 * ke + h ) * the )

               fthe = dsin( the / 180.0d0 * pi) *
     &                exp( prdef * g * the ) / fthemax

            end do

            the = the / 180.0d0 * pi

!-----------------------------------------------------------------------
!	              determination of parameter 'phi'
!-----------------------------------------------------------------------

	      if( r(1,1) .gt. 0. .and. r(2,1) .gt. 0. )then		!\91\E61\8Fی\C0

               phi = pi / 2. * rn(0)
!	         num = 1

	      elseif( r(1,1) .gt. 0. .and. r(2,1) .lt. 0. )then		!\91\E62\8Fی\C0

               phi = -pi / 2. * rn(0)
!	         num = 2

	      elseif( r(1,1) .lt. 0. .and. r(2,1) .lt. 0. )then		!\91\E63\8Fی\C0

               phi = pi + pi / 2. * rn(0)
!	         num = 3

	      elseif( r(1,1) .lt. 0. .and. r(2,1) .gt. 0.)then		!\91\E64\8Fی\C0

               phi = -pi - pi / 2. * rn(0)
!	         num = 4

	      elseif( r(1,1) .ge. 0. .and. r(2,1) .eq. 0. )then		!x\8E\B2\8F\E30<x<1

               phi = 0.
!	         num = 5

	      elseif( r(1,1) .lt. 0. .and. r(2,1) .eq. 0. )then		!x\8E\B2\8F\E3-1<x<0

               phi = pi
!	         num = 6

	      elseif( r(1,1) .eq. 0. .and. r(2,1) .gt. 0. )then		!y\8E\B2\8F\E30<y<1

               phi = 3. / 2. * pi
!	         num = 7

	      elseif( r(1,1) .eq. 0. .and. r(2,1) .lt. 0. )then		!y\8E\B2\8F\E3-1<y<0

               phi = 1. / 2. * pi
!	         num = 8

	      end if

!	write(*,*)'location number is',num

!-----------------------------------------------------------------------

	else

          if( ncolexb .eq. 0) then  
!-----------------------------------------------------------------------
!	Determination of angle[deg] of deflection
!	to nucleous outward direction
!	based on elastic scattering angular distribution.
!-----------------------------------------------------------------------

            fthe = -10.0

!            themax = 180./pi*abs(acos(p(3,i)/pp0))
            themax = 180.

            do while ( opt .ge. fthe )

              if(g .lt. -2.5) then
                  the = 0
                  exit
              end if

               the = 180.0d0 * rn(0)

               if( iclst(i) .eq. 0) then
                   the = themax * rn(0)
               end if

               opt = rn(0)
              fthe = dsin(the / 180.0d0 * pi) *
     &                exp( g * the )

            end do
!-----------------------------------------------------------------
! for GR_deflection katayama add.
!-----------------------------------------------------------------

         else if( ncolexb .eq. 8)then  ! L = 0
          agr = 25.0d0
          bkein = dsqrt(2.0d0 * 938.27231 * ein)
     &         /(6.58211 * 2.99792 * 10.0)
!          bkout = dsqrt(2.0d0 * 939.56541 * (ein - eisq))
!     &         /(6.58211 * 2.99792 * 10.0)
          bkout = dsqrt(2.0d0 * 939.56541 * ke)
     &         /(6.58211 * 2.99792 * 10.0)
!          qqrr = dsqrt(4.0**2 * (bkein**2 + bkout**2
!     &         - (2.0d0 * (bkein * bkout))
!     &         * dcos(the / 180.0d0 * pi)))
          fthe = -10.0
          do while ( opt .ge. fthe )
             the = agr * rn(0)
             qqrr = dsqrt(4.0**2 * (bkein**2 + bkout**2
     &         - (2.0d0 * (bkein * bkout))
     &         * dcos(the / 180.0d0 * pi)))
             opt = rn(0)
             fthe = 50.0d0
     &            *dsin(the / 180.0 * pi)*(dsin(qqrr)/qqrr)**2
          end do
          
         else if( ncolexb .eq. 1)then !L = 1
          agr = 25.0d0
          bkein = dsqrt(2.0d0 * 938.27231 * ein)
     &         /(6.58211 * 2.99792 * 10.0)
!          bkout = dsqrt(2.0d0 * 939.56541 * (ein - eisq))
!     &         /(6.58211 * 2.99792 * 10.0)
          bkout = dsqrt(2.0d0 * 939.56541 * ke)
     &         /(6.58211 * 2.99792 * 10.0)
!          qqrr = dsqrt(5.5**2 * (bkein**2 + bkout**2
!     &         - (2.0d0 * (bkein * bkout))
!     &         * dcos(the / 180.0d0 * pi)))
          fthe = -10.0
          do while ( opt .ge. fthe )
             the = agr * rn(0)
             qqrr = dsqrt(5.5**2 * (bkein**2 + bkout**2
     &         - (2.0d0 * (bkein * bkout))
     &         * dcos(the / 180.0d0 * pi)))
             opt = rn(0)
             fthe = dsin(the / 180.0 * pi ) * 5.0d0
     &            * ((dsin(qqrr)/qqrr**2) - (dcos(qqrr)/qqrr))**2
          end do
  
 
         else if( ncolexb .eq. 2)then  ! L = 2
          agr = 25.0d0
          bkein = dsqrt(2.0d0 * 938.27231 * ein)
     &         /(6.58211 * 2.99792 * 10.0)
!          bkout = dsqrt(2.0d0 * 939.56541 * (ein - eisq))
!     &         /(6.58211 * 2.99792 * 10.0)
          bkout = dsqrt(2.0d0 * 939.56541 * ke)
     &         /(6.58211 * 2.99792 * 10.0)
!          qqrr = dsqrt(7.0**2 * (bkein**2 + bkout**2
!     &         - (2.0d0 * (bkein * bkout))
!     &         * dcos(the / 180.0d0 * pi)))
          fthe = -10.0
          do while ( opt .ge. fthe )
             the = agr * rn(0)
             qqrr = dsqrt(7.0**2 * (bkein**2 + bkout**2
     &         - (2.0d0 * (bkein * bkout))
     &         * dcos(the / 180.0d0 * pi)))
             opt = rn(0)
             fthe = dsin(the / 180.0 * pi) *
     &            10.0d0 * ((((3.0d0 / qqrr** 3) - (1.0d0 / qqrr))
     &            * dsin(qqrr)) - ((3.0d0 / qqrr**2) * dcos(qqrr)))**2
          end do
          
         else if( ncolexb .eq. 4)then  !  IAS              
          agr = 25.0d0
          bkein = dsqrt(2.0d0 * 938.27231 * ein)
     &         /(6.58211 * 2.99792 * 10.0)
!          bkout = dsqrt(2.0d0 * 939.56541 * (ein - eisq))
!     &         /(6.58211 * 2.99792 * 10.0)
          bkout = dsqrt(2.0d0 * 939.56541 * ke)
     &         /(6.58211 * 2.99792 * 10.0)
!          qqrr = dsqrt(7.5**2 * (bkein**2 + bkout**2
!     &         - (2.0d0 * (bkein * bkout))
!     &         * dcos(the / 180.0d0 * pi)))
          fthe = -10.0
          do while ( opt .ge. fthe )
             the = agr * rn(0)
             qqrr = dsqrt(7.5**2 * (bkein**2 + bkout**2
     &         - (2.0d0 * (bkein * bkout))
     &         * dcos(the / 180.0d0 * pi)))
             opt = rn(0)
             fthe = 8.0d0 *
     &            dsin(the / 180.0 * pi)*(dsin(qqrr)/qqrr)**2
          end do

         else if( ncolexb .eq. 3)then ! L=3
          agr = 25.0d0
          bkein = dsqrt(2.0d0 * 938.27231 * ein)
     &         /(6.58211 * 2.99792 * 10.0)
!          bkout = dsqrt(2.0d0 * 939.56541 * (ein - eisq))
!     &         /(6.58211 * 2.99792 * 10.0)
          bkout = dsqrt(2.0d0 * 939.56541 * ke)
     &         /(6.58211 * 2.99792 * 10.0)
!          qqrr = dsqrt(5.5d0**2 * (bkein**2 + bkout**2
!     &         - (2.0d0 * (bkein * bkout))
!     &         * dcos(the / 180.0d0 * pi)))
          fthe = -10.0
          do while ( opt .ge. fthe )
             the = agr * rn(0)
             qqrr = dsqrt(5.5d0**2 * (bkein**2 + bkout**2
     &         - (2.0d0 * (bkein * bkout))
     &         * dcos(the / 180.0d0 * pi)))
             opt = rn(0)
             fthe = dsin(the / 180.0 * pi) *
     &            10.0d0 * ((((3.0d0 / qqrr** 3) - (1.0d0 / qqrr))
     &            * dsin(qqrr)) - ((3.0d0 / qqrr**2) * dcos(qqrr)))**2
          end do
          
        end if
          
          the = the / 180.0d0 * pi

          phi = 2.0d0 * pi * rn(0)
          
          if (colexb .eq. 1) then
              phi = 0
          end if
       end if

!-----------------------------------------------------------------------
!           set the new momentum coordinates after deflection
!-----------------------------------------------------------------------

         if( dabs(ome3) .ge. 1.0d0 ) then

            dome1 = dsin(the) * dcos(phi)

            dome2 = dsin(the) * dsin(phi)

            dome3 = dcos(the)

         else

            dome1 = ( ome1 * ome3 * dsin(the) * dcos(phi)
     &            - ome2 * dsin(the) * dsin(phi))
     &            / dsqrt( 1.0d0 - ome3**2.0d0) + ome1 * dcos(the)

            dome2 = ( ome2 * ome3 * dsin(the) * dcos(phi)
     &            + ome1 * dsin(the) * dsin(phi) ) 
     &            / dsqrt( 1.0d0 - ome3**2.0d0) + ome2 * dcos(the)

            dome3 = ome3 * dcos(the) - dsin(the) * dcos(phi)
     &            * dsqrt( 1.0d0 - ome3**2.0d0)

         end if

         p(1,i) = pp0 * dome1

         p(2,i) = pp0 * dome2

         p(3,i) = pp0 * dome3

!-------------------------------------------8----------------------------

	if( ( nbend .eq. 1 ) .and. ( i .eq. 1 ) )return
        
        enddo

      END SUBROUTINE deflection

!***********************************************************************
!*                                                                     *
      SUBROUTINE kalbach( nmasta,nchgta,nmaspr,ein,upot,ucpot,p,iclst,
     & iproj )
!*                                                                     *
!*                                                                     *
!*        Purpose:                                                     *
!*                                                                     *
!*              to calculate the angle of orientation bent by 	       *
!*		    kalbach                                            *
!*                                                                     *
!*                                                                     *
!***********************************************************************

      implicit real(8) (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      integer, intent(in)    :: nmasta, nmaspr,nchgta,iclst(nnn),iproj
      real(8), intent(in)    :: upot,ucpot,ein
      real(8), intent(inout) :: p(6,nnn)
      real :: ke,ibach,nbbach


      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )

!-----------------------------------------------------------------------
	if( iproj .eq. 1 )then 
        do i = 1 , nmasta + nmaspr
           if( i .ne. 1 ) cycle!t.mori
           if( iclst(i) .ne. 1 ) cycle

           pp0   = 0.0d0
           psi = 0.0d0
           phi = 0.0d0
           the = 0.0d0

           ome1 = 0.0d0
           ome2 = 0.0d0
           ome3 = 0.0d0

           dome1 = 0.0d0
           dome2 = 0.0d0
           dome3 = 0.0d0

           opt = 0.0d0
           fthe = 0.0d0
           ke = 0.0d0


           if( iclst(i) .eq. 0 ) cycle!t.mori
           if( p(5,i) .le. 0.1d0 ) cycle

           ke = p(4,i) - p(5,i) - upot

!------------------------------sonoda11/12

	  if( iclst(i) .eq. 1 )then
           ke = p(4,i) - p(5,i) - dupot
	  else if( iclst(i) .eq. 2 .or. iclst(i) .eq. 3)then
           ke = p(4,i) - p(5,i) - thupot
	  else if( iclst(i) .eq. 4 )then
           ke = p(4,i) - p(5,i) - aupot
	  else if( iclst(i) .eq. 5 )then
           ke = p(4,i) - p(5,i) - 6.0d0*upot
	  else if( iclst(i) .eq. 6 .or. iclst(i) .eq. 7 )then
           ke = p(4,i) - p(5,i) - 7.0d0*upot
	  else if( iclst(i) .eq. 8 )then
           ke = p(4,i) - p(5,i) -8.0d0*upot
	  else if( iclst(i) .eq. 9 )then
           ke = p(4,i) - p(5,i) - 9.0d0*upot
	  else if( iclst(i) .eq. 10 .or. iclst(i) .eq. 11 )then
           ke = p(4,i) - p(5,i) - 10.0d0*upot
	  else if( iclst(i) .eq. 12 .or. iclst(i) .eq. 13 )then
           ke = p(4,i) - p(5,i) - 11.0d0*upot
	  else if( iclst(i) .eq. 14 )then
           ke = p(4,i) - p(5,i) - 12.0d0*upot
	  end if


!------------------------------sonoda11/12

           if( ke .le. 0.0d0 ) cycle

!================================================sonoda1220

!========================================================sonoda0217
!               parameter of Kalbach equation                    ==
!==================================================================

      	acbach = dble(nmasta) !target

      	if(iclst(i) .eq. 1) then !d
      	abbach = acbach - 1.0d0
      	else if(iclst(i) .eq. 2) then !t
      	abbach = acbach - 2.0d0
      	else if(iclst(i) .eq. 3) then !h
      	abbach = acbach - 2.0d0
      	else if(iclst(i) .eq. 4) then !alpha
      	abbach = acbach - 3.0d0
      	else if(iclst(i) .eq. 5) then !6Li
      	abbach = acbach - 5.0d0
      	else if(iclst(i) .eq. 6 .or. iclst(i) .eq. 7) then !7Li,7Be
      	abbach = acbach - 6.0d0
      	else if(iclst(i) .eq. 8) then !8Be
      	abbach = acbach - 7.0d0
      	else if(iclst(i) .eq. 9) then !9Be
      	abbach = acbach - 8.0d0
      	else if(iclst(i) .eq. 10 .or. iclst(i) .eq. 11) then !10Be,10B
      	abbach = acbach - 9.0d0
      	else if(iclst(i) .eq. 12 .or. iclst(i) .eq. 13) then !11B,11C
      	abbach = acbach - 10.0d0
      	else if(iclst(i) .eq. 4) then !12C
      	abbach = acbach - 11.0d0
      	end if

        ncbach = dble(nmasta - nchgta)

        if(iclst(i) .eq. 1) then !d
        nbbach = ncbach -1.0d0
        else if(iclst(i) .eq. 2) then !t
        nbbach = ncbach -2.0d0
        else if(iclst(i) .eq. 3) then !h
        nbbach = ncbach -1.0d0
        else if(iclst(i) .eq. 4) then !alpha
        nbbach = ncbach -2.0d0
        else if(iclst(i) .eq. 5) then !6Li
        nbbach = ncbach -3.0d0
        else if(iclst(i) .eq. 6) then !7Li
        nbbach = ncbach -4.0d0
        else if(iclst(i) .eq. 7) then !7Be
        nbbach = ncbach -3.0d0
        else if(iclst(i) .eq. 8) then !8Be
        nbbach = ncbach -4.0d0
        else if(iclst(i) .eq. 9) then !9Be
        nbbach = ncbach -5.0d0
        else if(iclst(i) .eq. 10) then !10Be
        nbbach = ncbach -6.0d0
        else if(iclst(i) .eq. 11) then !10B
        nbbach = ncbach -5.0d0
        else if(iclst(i) .eq. 12) then !11B
        nbbach = ncbach -6.0d0
        else if(iclst(i) .eq. 13) then !11C
        nbbach = ncbach -5.0d0
        else if(iclst(i) .eq. 14) then !12C
        nbbach = ncbach -6.0d0
        end if

      	zcbach = dble(nchgta)
      	zbbach = zcbach
      
      	if(iclst(i) .eq. 1) then !d
      	ibach = 2.224628d0
      	else if(iclst(i) .eq. 2) then !t
      	ibach = 8.481960d0
      	else if(iclst(i) .eq. 3) then !h
      	ibach = 7.718183d0
      	else if(iclst(i) .eq. 4) then !alpha
      	ibach = 28.29599d0
      	else if(iclst(i) .eq. 5) then !6Li
      	ibach = 31.980d0
      	else if(iclst(i) .eq. 6) then !7Li
      	ibach = 39.265d0
      	else if(iclst(i) .eq. 7) then !7Be
      	ibach = 36.716d0
      	else if(iclst(i) .eq. 8) then !8Be
      	ibach = 56.480d0
      	else if(iclst(i) .eq. 9) then !9Be
      	ibach = 58.14d0
      	else if(iclst(i) .eq. 10) then !10Be
      	ibach = 61.137d0
      	else if(iclst(i) .eq. 11) then !10B
      	ibach = 64.8d0
      	else if(iclst(i) .eq. 12) then !11B
      	ibach = 76.23d0
      	else if(iclst(i) .eq. 13) then !11C
      	ibach = 72.543d0
      	else if(iclst(i) .eq. 14) then !12C
      	ibach = 92.160d0
      	end if
 

      	sabach = 15.68d0 * (acbach - abbach)

      	sabach = sabach - 28.07d0 * ((ncbach - zcbach)**2
     &        /acbach - (nbbach - zbbach)**2/abbach)

      	sabach = sabach - 18.56d0 * (acbach**(2./3.) - abbach**(2./3.))

      	sabach = sabach + 33.22d0 * ((ncbach - zcbach)**2
     &         /acbach**(4./3.) - (nbbach - zbbach)**2
     &         /abbach**(4./3.))

      	sabach = sabach - 0.717d0 * (zcbach**2/acbach**(1./3.)
     &         - zbbach**2/abbach**(1./3.))

      	sabach = sabach + 1.211d0 * (zcbach**2/acbach
     &         - zbbach**2/abbach)

      	sbbach = sabach - ibach

     	eabach = ein + sabach

      	ebbach = ke  + sbbach

      	e1bach = min(eabach,130.0d0) !130\81}10MeV

      	e2bach = min(eabach,41.0d0) !41\81}5MeV

        x1bach = e1bach * ebbach / eabach

	x2bach = e2bach * ebbach / eabach

	abach  = 4.0d-2 * x1bach + 1.8d-6 * x1bach**3
     &         + 6.7d-7 * x2bach**4
        pp0 = dsqrt( p(1,i)**2.0d0 + p(2,i)**2.0d0 + p(3,i)**2.0d0 )

        ome1 = p(1,i) / pp0
        ome2 = p(2,i) / pp0
        ome3 = p(3,i) / pp0

        the = 180.0d0 * rn(0)

        opt = rn(0)

	fthe = dsin(the / 180.0d0 * pi) * 
     &      exp(abach * dcos(the / 180.0d0 * pi))
     &       / exp(abach)

          do while ( opt .ge. fthe )

            the = 180.0d0 * rn(0)

            opt = rn(0)

	   fthe = dsin(the / 180.0d0 * pi) * 
     &       exp(abach * dcos(the / 180.0d0 * pi))
     &       / exp(abach)

          end do

        the = the / 180.0d0 * pi

        phi = 2.0d0 * pi * rn(0)

          if( dabs(ome3) .ge. 1.0d0 ) then

          dome1 = dsin(the) * dcos(phi)

          dome2 = dsin(the) * dsin(phi)

          dome3 = dcos(the)

          else

          dome1 = ( ome1 * ome3 * dsin(the) * dcos(phi)
     &            - ome2 * dsin(the) * dsin(phi))
     &            / dsqrt( 1.0d0 - ome3**2.0d0) + ome1 * dcos (the)

          dome2 = ( ome2 * ome3 * dsin(the) * dcos(phi)
     &            + ome1 * dsin(the) * dsin(phi) ) 
     &            / dsqrt( 1.0d0 - ome3**2.0d0) + ome2 * dcos(the)

          dome3 = ome3 * dcos(the) - dsin(the) * dcos(phi)
     &            * dsqrt( 1.0d0 - ome3**2.0d0)

          end if

         p(1,i) = pp0 * dome1

         p(2,i) = pp0 * dome2

         p(3,i) = pp0 * dome3

           if( p(1,i) .le. 0.0d0 .and. p(2,i) .ge. 0.0d0) then

             phi = datan( p(2,i) / p(1,i))

           else if(p(1,i) .le. 0.0d0 .and. p(2,i) .le. 0.0d0) then

             phi = 180.0d0 + ( datan( (p(2,i)) / p(1,i) )
     &                             / pi * 180.0d0)

             phi = phi / 180.0d0 * pi

           else if(p(1,i) .ge. 0.0d0 .and. p(2,i) .le. 0.0d0) then

             phi = datan( p(2,i) / p(1,i))

             phi =  360.0d0 + phi / pi * 180.0d0

             phi = phi / 180.0d0 * pi

           else

             phi = datan( p(2,i) / p(1,i) )

           end if

!             psi = pasi / 180.0d0 * pi
             
         end do

	end if !sono1

      END SUBROUTINE kalbach

!**********************************************************************
!*
      SUBROUTINE dpickup(iclst,iproj,icoales,ipot,inio,nmaspr,nmasta
     &    ,nchgta,ein,upot,dupot,thupot,aupot,p,r,peex,ichg,inds,tfm
     &    ,radm,numclst,ininclst,nmasclst,pcmclst,rcmclst,relr)

!*                                                                    *
!*                                                                    *
!*      Purpose:                                                      *
!*                                                                    *
!*            incident particle pickup nucleon                        *
!*                                                                    *
!*                                                                    *
!**********************************************************************

!      use NGSDATAMOD, only : energm

      implicit real(8) (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      real(8), intent(inout) :: p(6,nnn),r(3,nnn),peex,
     &		pcmclst(5,nnn),rcmclst(3,nnn),relr(3,nnn)
      real(8), intent(in) :: upot,dupot,thupot,aupot,ein,tfm,radm

!----------------------------sonoda11/13
      integer, intent(in) :: nmaspr,nmasta,iproj,inds(nnn)
!----------------------------sonoda11/13
      integer, intent(inout) :: iclst(nnn),icoales(nnn),ipot(nnn),
     &     inio(nnn),ichg(nnn),numclst(nnn),ininclst(nnn)
      integer isel(8),dbpass,idn2,idn1,idp1,nmasclst(nnn)

!******************Initialization*********************************

        i = 0
        j = 0
        k = 0

        ! Pick up reaction number
        ! 1 : deuteron / 2 : triron / 3 : 3He / 4 : alpha / 0 : no
        dbpass = 0

        ! Picked up Particle number
        idn1 = 0
        idn2 = 0
        idp1 = 0

        rsq1 = 10.0d0
        rsq2 = 20.0d0

        ! Highest grand energy
!       tn1 = nmasta + nmaspr     ! 1st neutron 
!       tn2 = nmasta + nmaspr - 1 ! 2rd neutron
!       tp1 = nchgta + nmaspr     ! 1st proton

        npdp = nchgta !Z
        nndp = nmasta - nchgta !N

!*****************************************************************

        ! This subroutine is for p incident.
        !0603/tmori for low energy p in

!        if( iproj .eq. 1 .and. ein .lt. 400.0d0  )then
!         if( ichg(1) .eq. 1 .and. iclst(1) .eq. 0 )then

        if( ein .lt. 400.0d0  )then
         if( iproj .eq. 1 .or. iproj .eq. 0 )then
!++++++++++++++++++Determine reaction+++++++++++++++++++++++++++++

	  peliel = rn(0)
	  frand = rn(0)

          brannp = 0.0
          rrate0 = 0.0d0
          einp = p(4,1) - rmass - upot

	  pcmclst = 0.0
	  rcmclst = 0.0
	  relr = 0.0

	  numclst = 0
	  nmasclst = 0.0

!-calcurate-Q-value---------------------------------------------------

!          CALL setup ! call gem's data

          amcpr  = energmelf(1,1,2)
          amcta  = energmelf(nchgta,nmasta,2)

          amcem  = energmelf(1,2,2)
          amcrs  = energmelf(nchgta,nmasta - 1,2)
 
          qd = amcem + amcrs - amcpr - amcta

          amcem  = energmelf(1,3,2)
          amcrs  = energmelf(nchgta,nmasta - 2,2 )

          qt = amcem + amcrs - amcpr - amcta

          amcem  = energmelf(2,3,2)
          amcrs  = energmelf(nchgta - 1,nmasta - 2,2 )

          qh = amcem + amcrs - amcpr - amcta

          amcem  = energmelf(2,4,2)
          amcrs  = energmelf(nchgta - 1,nmasta - 3,2 )

          qa = amcem + amcrs - amcpr - amcta

!---------------------------------------------------------------------
!--------------------parameter--each--pickup--reaction----------------
!---------------------------------------------------------------------
!     Probability for direct pickup(deuteron) depending on incident E
!     yamaguchi
!---------------------------------------------------------------------

!          if( ein .lt. 79.0 )then
!             dcrs = 0.1        ! p,d
!          elseif( ein .ge. 79.0 .and. ein .lt. 195.0 )then
!             dcrs = 0.0
!          elseif( ein .ge. 195.0 .and. ein .lt. 346.0 )then
!             dcrs = 0.0
!          else
!             dcrs = 0.1
!          end if 
!             dcr1=-0.00003*ein+0.11
!             dcr2=1.001669449**ein
!             dcr3=exp(-0.003*(ein**dcr2))
!             dcrs=dcr1*dcr3

!---for deuteron parameter
          dcrs = 0.11*exp(-0.00015*(ein)**1.9)
          
!---for triton parameter
          if( ein .lt. 70.0 )then ! for Al (hh0 seems to be low)
          if( nchgta .eq. 6 )then
             th0 = 1.0d-2
             th1 = 0.0d0
             th2 = 0.0d0
          else if( nchgta .eq. 14 )then
             th0 = 1.0d-2
             th1 = 0.0d0
             th2 = 1.4d-2
          else if( nchgta .eq. 28 )then
             th0 = 1.0d-2
             th1 = 0.0d0
             th2 = 1.4d-2
          else if( nchgta .eq. 40 )then
             th0 = 1.4d-2
             th1 = 0.0d0
             th2 = 1.4d-2
          else if( nchgta .eq. 79 )then
             th0 = 1.4d-2
             th1 = 0.0d0
             th2 = 1.4d-2
          else
             th0 = 1.0d-2
             th1 = 0.0d0
             th2 = 1.4d-2
          end if
          
          else                  
             if( nchgta .eq. 27 )then ! Co
                th0 = 1.00d-4
                th1 = 3.0d-3
                th2 = 1.0d-2
             else if( nchgta .eq. 28 )then ! Ni temporary
                th0 = 1.00d-4
                th1 = 3.0d-3
                th2 = 1.0d-2
             else ! Au
                th0 = 1.00d-4
                th1 = 3.0d-3
                th2 = 5.0d-3
             end if
          end if
          
          tcrs = th0 + th1 + th2
          th0 = th0/tcrs
          th1 = th1/tcrs
          th2 = th2/tcrs
          
!---for 3He parameter
          if( ein .lt. 70.0 )then ! for Al (hh0 seems to be low)
          if( nchgta .eq. 6 )then
             hh0 = 1.0d-2
             hh1 = 0.0d0
             hh2 = 0.0d0
          else if( nchgta .eq. 14 )then
             hh0 = 1.0d-2
             hh1 = 0.0d0
             hh2 = 1.4d-2
          else if( nchgta .eq. 28 )then
             hh0 = 1.0d-2
             hh1 = 0.0d0
             hh2 = 1.4d-2
          else if( nchgta .eq. 40 )then
             hh0 = 1.4d-2
             hh1 = 0.0d0
             hh2 = 1.4d-2
          else if( nchgta .eq. 79 )then
             hh0 = 1.4d-2
             hh1 = 0.0d0
             hh2 = 1.4d-2
          else
             hh0 = 1.0d-2
             hh1 = 0.0d0
             hh2 = 1.4d-2
          end if
          else                  
             if( nchgta .eq. 27 )then ! Co
                hh0 = 1.00d-4
                hh1 = 3.0d-3
                hh2 = 1.0d-2
             else if( nchgta .eq. 28 )then ! Ni temporary
                hh0 = 1.00d-4
                hh1 = 3.0d-3
                hh2 = 1.0d-2
             else ! Au
                hh0 = 1.00d-4
                hh1 = 3.0d-3
                hh2 = 5.0d-3
             end if
          end if
          hcrs = hh0 + hh1 + hh2
          hh0 = hh0/hcrs
          hh1 = hh1/hcrs
          hh2 = hh2/hcrs
!---for alpha parameter serch
          if( ein .lt. 65.0 )then
             if( nchgta .eq. 6 )then ! for C (ah0 seems to be low)
                ah0 = 2.0d-3
                ah1 = 2.7d-3
                ah2 = 0.0d0          
                ah3 = 2.85d-2
             else if( nchgta .eq. 39 )then ! Y
                ah0 = 1.9d-3
                ah1 = 1.75d-3
                ah2 = 0.0d0          
                ah3 = 5.5d-2
             else if( nchgta .eq. 83 )then ! Bi
                ah0 = 5.0d-4
                ah1 = 1.0d-3
                ah2 = 0.0d0          
                ah3 = 4.0d-2
             else ! Ni
                ah0 = 2.0d-3
                ah1 = 5.0d-3
                ah2 = 0.0d0          
                ah3 = 5.0d-2
             end if
          else if( ein .lt. 70.0 )then ! Ni
             ah0 = 9.5d-4
             ah1 = 2.5d-3
             ah2 = 0.0d0          
             ah3 = 4.05d-2             
          else if( ein .lt. 100.0 )then ! Ni
             ah0 = 4.5d-4
             ah1 = 1.0d-3
             ah2 = 0.0d0          
             ah3 = 3.15d-2             
          else if( ein .lt. 150.0 )then ! Co
             ah0 = 1.0d-4
             ah1 = 2.5d-4
             ah2 = 0.0d0          
             ah3 = 1.5d-2             
          else if( ein .lt. 190.0 )then ! Co
             ah0 = 7.0d-5
             ah1 = 2.0d-4
             ah2 = 0.0d0          
             ah3 = 4.0d-3
          else if( ein .lt. 250.0 )then ! Co
             ah0 = 1.0d-5
             ah1 = 3.0d-5
             ah2 = 0.0d0          
             ah3 = 1.0d-3
          else
             ah0 = 1.0d-6
             ah1 = 1.0d-5
             ah2 = 0.0d0          
             ah3 = 1.0d-4
          end if
          acrs = ah0 + ah1 + ah2 +ah3
          ah0 = ah0/acrs
          ah1 = ah1/acrs
          ah2 = ah2/acrs
          ah3 = ah3/acrs
          
!          acrs = 0.05
 
          if( einp .lt. qd )dcrs = 0.0
          if( einp .lt. qt )tcrs = 0.0
          if( einp .lt. qh )hcrs = 0.0
          if( einp .lt. qa )acrs = 0.0

          drang = dcrs
          trang = drang + tcrs
          hrang = trang + hcrs
          arang = hrang + acrs

!++start++deuteron+++++++++++++++++++++++++++++++++++++++++++++++++++++
! 

	   if( peliel .lt. drang ) then ! deuteron

             dbpass = 1

             if( qd .gt. 0.0 )then
               pexmax0 = einp - qd
             else
               pexmax0 = einp
             end if

             pexmax = min(50.0d0,pexmax0)!Limit excitation energy 50MeV

             rsq1 = 1000.0d0
             rsq2 = 2000.0d0

             do j = nmaspr+1, nmasta + nmaspr

               if( ipot(j) .eq. 1 ) cycle
               if( iclst(j) .ne. 0 ) cycle
               if( inds(j) .ne. 1 ) cycle
               if( ichg(j) .eq. ichg(1) ) cycle
               if( p(5,j) .lt. 0.10d0 ) cycle

                rsq = ( r(1,1) - r(1,j) )**2
     &              + ( r(2,1) - r(2,j) )**2
     &              + ( r(3,1) - r(3,j) )**2

               rsq0 = sqrt(rsq)

                if( rsq0 .lt. rsq1 )then
                  rsq1 = rsq0 
                  idn1 = j
                end if

             end do
!williams
	      eex = pexmax * rn(0)
	      apex = 10.0d0
	      bpex = -10.0d0
              rrate0 = - 10.0d0

!              do j = 1,100
!
!               eex = pexmax*dble(j)/100.0d0
!
!               wisd = (eex)/(pexmax)! Williams state density
!
!               gw = (48.6)/dble(nmasta)
!               fadp =  gw/( eex**2.0 + gw )!Factor direct pick up
!
!               if( eex .ge. 35. )then
!                exlim = exp( - (eex - 35.)/10.0 )
!               else
!                exlim = 1.0d0
!               end if
!
!               rrate1 = wisd * fadp * exlim
!
!               if( rrate1 .gt. rrate0 )rrate0 = rrate1 !\8Dő\E5\92l
!
!              end do

!	      eex = pexmax * rn(0)
	      apex = 10.0d0
	      bpex = -10.0d0

	      do while( apex .gt. bpex )

	       apex = rn(0)!Random number

	       eex = pexmax * rn(0)!Excitation energy
!              wisd = (eex)/(pexmax)
               wisd = eex*1.0d-2    !for 27Al
!              wisd = eex*2.5d-2    !for 58Ni
!               wisd = 0.0
!              gw = (48.6)/dble(nmasta)
!              fadp =  gw/( eex**2.0 + gw ) !Factor direct pick up
!              fadp =  25.0d0/( eex**2 + 25.0d0 ) !Factor direct pick up for 27Al width 10 MeV
               fadp =  100.0d0/( eex**2 + 100.0d0 ) !Factor direct pick up for 27Al width 20 MeV
!              fadp =  64.0d0/( eex**2 + 64.0d0 ) !Factor direct pick up for 58Ni

               if( eex .ge. 35. )then
!               exlim = exp( - (eex - 35.)/10.0 )
                exlim = (55.0d0 - eex)*0.05d0
               else
                exlim = 1.0d0
               end if

!              bpex = wisd * fadp * exlim / rrate0 !Height of distribution func. at eex
!              bpex = wisd + fadp !Height of distribution func. at eex
               bpex = wisd*exlim + fadp !Height of distribution func. at eex
               
	      end do

	      dpex = eex
              epick = dpex + rmass
              ppick = dsqrt( epick**2.0d0 - rmass**2.0d0 )

              ux2 = p(1,1)
              uy2 = p(2,1)
              uz2 = p(3,1)
              ptot0 = dsqrt( ux2**2.0 + uy2**2.0 + uz2**2.0  )

	      tc2 = einp + dupot - dpex
             
              pc2 = sqrt( tc2**2 + 2.0d0 * ddmass * tc2 )/2.0d0

              p(1,1) = pc2 *  ux2 / ptot0
              p(2,1) = pc2 *  uy2 / ptot0
              p(3,1) = pc2 *  uz2 / ptot0

              do k = 1, 3
               p(k,idn1)  = p(k,1) !pickup\82\B3\82\EA\82钆\90\AB\8Eq\82̉^\93\AE\97ʁA\97z\8Eq\82Ɠ\AF\82\B6
              end do

              p(4,1) = sqrt( rmass **2
     &         + p(1,1)**2 + p(2,1)**2 + p(3,1)**2 )

              p(4,idn1) = p(4,1)

	      p(5,1) = rmass
	      p(5,idn1) = rmass

	      do j = 1,3
		pcmclst(j,1) = p(j,1) + p(j,idn1) !\8Fd\97z\8Eq\82̏d\90S\82̉^\93\AE\97\CA
	      end do

	      pcmclst(5,1) = ddmass
	      pcmclst(4,1) = dsqrt( pcmclst(1,1)**2
     &				+ pcmclst(2,1)**2
     &				+ pcmclst(3,1)**2
     &				+ pcmclst(5,1)**2 )


!sonoda density 0125----------------------------------------------

            rcmx = r(1,1)
            rcmy = r(2,1)
            rcmz = r(3,1)

	    rcmclst(1,1) = rcmx !\8Fd\97z\8Eq\82̏d\90S\82̈ʒu
	    rcmclst(2,1) = rcmy
	    rcmclst(3,1) = rcmz


            rmax = 1.0
            dr = 10.0

            do while( dr .gt. 1.0d0 )
             dx = rn(0)
             dy = rn(0)
             dz = rn(0)
             dr = dsqrt( dx**2 + dy**2 + dz**2 )
            end do

            r(1,1)  = rmax*dx + rcmx 
            r(2,1)  = rmax*dy + rcmy 
            r(3,1)  = rmax*dz + rcmz 

            relr(1,1)  = rmax*dx
            relr(2,1)  = rmax*dy
            relr(3,1)  = rmax*dz

            dr = 10.0

            do while( dr .gt. 1.0d0 )
             dx = rn(0)
             dy = rn(0)
             dz = rn(0)
             dr = dsqrt( dx**2.0 + dy**2.0d0 + dz**2.0 )
            end do

            r(1,idn1)  = rmax*dx + rcmx 
            r(2,idn1)  = rmax*dy + rcmy 
            r(3,idn1)  = rmax*dz + rcmz 

            relr(1,idn1)  = rmax*dx
            relr(2,idn1)  = rmax*dy
            relr(3,idn1)  = rmax*dz


!            ichg(1) = 1
!            ichg(idn1) = 0

            iclst(1) = 0
            iclst(idn1) = 0
             
 
            
	    ininclst(1) = 1
	    ininclst(idn1) = 1

            numclst(1) = 1
            numclst(idn1) = 1

	    nmasclst(1) = 2

            icoales(1) = 1
            icoales(idn1) = 1

            inio(1) = 0
            inio(idn1) = 0

            !ipot(1) = 1
            !ipot(idn1) = 1 

                
                ! check moving particle
!               theta=0
!             if( sqrt( p(1,1)**2. + p(2,1)**2. + p(3,1)**2. ).eq.0 )then
!                    theta = 0.0
!             else
!                     theta = 57.3*acos( p(3,1)      
!     &                   / sqrt( p(1,1)**2. + p(2,1)**2. + p(3,1)**2. ))
!             endif
           
             
          

!++++++++++++++++++++++++++++++++++++++++++++++++++++++end++deuteron+++

!++start++triton+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 

        else if( peliel .lt. trang ) then ! tmori 10/15

             dbpass = 2

             if( qt .gt. 0.0 )then
               pexmax0 = einp - qt
             else
               pexmax0 = einp
             end if

             pexmax = min(90.0d0,pexmax0)

             rsq1 = 1000.0d0
             rsq2 = 2000.0d0

             do j = nmaspr+1, nmasta + nmaspr
               if( ipot(j) .eq. 1 ) cycle
               if( iclst(j) .ne. 0 ) cycle
               if( inds(j) .ne. 1 ) cycle
               if( ichg(j) .eq. ichg(1) ) cycle
               if( p(5,j) .lt. 0.10d0 ) cycle

                rsq = ( r(1,1) - r(1,j) )**2
     &              + ( r(2,1) - r(2,j) )**2
     &              + ( r(3,1) - r(3,j) )**2

               rsq0 = sqrt(rsq)

                if( rsq0 .lt. rsq1 )then
                   rsq2 = rsq1
                   rsq1 = rsq0

                   idn2 = idn1
                   idn1 = j
                else if( rsq0 .lt. rsq2  )then
                  rsq2 = rsq0
                  idn2 = j
                end if

             end do

!************ come here if pauli blocked
1520         continue
!**************************************

!--- Breit Wigner
           
           !do i = 1, 2
	      eex = pexmax * rn(0)
	      apex = 10.0d0
	      bpex = -10.0d0

              rrate0 = - 10.0d0

	      do while( apex .gt. bpex )

	       apex = rn(0)

	       eex = pexmax * rn(0)

               fadp =  400.0d0/( eex**2 + 400.0d0 ) !Factor direct pick up for 27Al width 2 MeV
               bpex = fadp !Height of distribution func. at eex only breitwigner
            end do

	     ! if(i .eq. 1)
               dpex1 = eex/2.
             ! if(i .eq. 2)
               dpex2 = eex/2.
               
!williams
      
	      apex = 10.0d0
	      bpex = -10.0d0
	      pexmax =  min(upot*2.,pexmax0) 

	      do while( apex .gt. bpex )

	       apex = rn(0)
	       eex = pexmax * rn(0)
               wisd = ((eex + 1.)**2./(pexmax + 1.0)**2.)

               if( eex .ge. tfm*2. )then
!                exlim = exp( - (eex - 35.)/10.0 )
                exlim = (upot *2. - eex )/(upot*2. - tfm*2.)
               else
                exlim = 1.0d0
               end if

               bpex = wisd * exlim

	      end do

	      wdpex1 = eex/2.
	      wdpex2 = eex/2.

          if(  frand .gt. th2+th1 )then            
            tc11 = upot - dpex1  !181202 de_dephth add. fukuda
            tc12 = upot - dpex2
            
            if(dpex1+dpex2 .gt. einp)goto 1520

            idad = 10
         elseif(  frand .gt. th2 )then
          
            tc11 = upot - dpex1
            tc12 = upot - wdpex2 
            if(dpex1+wdpex2 .gt. einp)goto 1520
            idad = 11
         else
              
            tc11 = upot - wdpex1 
            tc12 = upot - wdpex2 
            
            if(wdpex1+wdpex2 .gt. einp)goto 1520
            idad = 12
         end if

          pc01=sqrt(p(1,idn1)**2.+p(2,idn1)**2.+p(3,idn1)**2.) 
          pc02=sqrt(p(1,idn2)**2.+p(2,idn2)**2.+p(3,idn2)**2.)
          pc11 = sqrt(tc11**2.0d0 + 2.0d0 * rmass * tc11 )
          pc12 = sqrt(tc12**2.0d0 + 2.0d0 * rmass * tc12 )
          do k = 1,3
             p(k,idn1) = p(k,idn1)*pc11/pc01
             p(k,idn2) = p(k,idn2)*pc12/pc02
          enddo
           
                      p(4,1)=sqrt(rmass**2.0+
     &           p(1,1)**2.0+p(2,1)**2.0+p(3,1)**2.0)
            p(4,idn1)=sqrt(rmass**2.0+
     &           p(1,idn1)**2.0+p(2,idn1)**2.0+p(3,idn1)**2.0)
            p(4,idn2)=sqrt(rmass**2.0+
     &           p(1,idn2)**2.0+p(2,idn2)**2.0+p(3,idn2)**2.0)
             tc1=p(4,1)-p(5,1)
             tc2=p(4,idn1)-p(5,idn1)
             tc3=p(4,idn2)-p(5,idn2)

	     tca = tc1 + tc2 + tc3  ! + 2.*(upot-tfm)
             pca = sqrt( tca**2.0 + 2.0d0 * hmass * tca )/3.0d0

             pp0 = dsqrt( p(1,1)**2. + p(2,1)**2. + p(3,1)**2. )

             ome1 = p(1,1) / pp0
             ome2 = p(2,1) / pp0
             ome3 = p(3,1) / pp0

!	     tc2 = einp + thupot - dpex
!             pc2 = sqrt( tc2**2 + 2.0d0 * tmass * tc2 )/3.0d0

             p(1,1) = pca * ome1
             p(2,1) = pca * ome2
             p(3,1) = pca * ome3

             do k = 1, 3
              p(k,idn1)  = p(k,1)
              p(k,idn2)  = p(k,1)
             end do

            p(4,1) = sqrt( rmass **2
     &         + p(1,1)**2 + p(2,1)**2 + p(3,1)**2 )

            p(4,idn1) = p(4,1)
            p(4,idn2) = p(4,1)

	    p(5,1) = rmass
	    p(5,idn1) = rmass
	    p(5,idn2) = rmass

	      do j = 1,3
		pcmclst(j,1) = p(j,1) + p(j,idn1) + p(j,idn2) !\8Fd\90S\82̉^\93\AE\97\CA
	      end do

	      pcmclst(5,1) = tmass
	      pcmclst(4,1) = dsqrt( pcmclst(1,1)**2
     &				+ pcmclst(2,1)**2
     &				+ pcmclst(3,1)**2
     &				+ pcmclst(5,1)**2 )

!-----------------------------------------------------------

            rcmx = (r(1,1) + r(1,idn1) + r(1,idn2))/3.0
            rcmy = (r(2,1) + r(2,idn1) + r(2,idn2))/3.0
            rcmz = (r(3,1) + r(3,idn1) + r(3,idn2))/3.0

	    rcmclst(1,1) = rcmx
	    rcmclst(2,1) = rcmy
	    rcmclst(3,1) = rcmz

            rmax = 1.0
            dr = 10.0

            do while( dr .gt. 1.0d0 )
             dx = rn(0)
             dy = rn(0)
             dz = rn(0)
             dr = dsqrt( dx**2.0 + dy**2.0d0 + dz**2.0 )
            end do

            r(1,1)  = rmax*dx + rcmx 
            r(2,1)  = rmax*dy + rcmy 
            r(3,1)  = rmax*dz + rcmz 

            relr(1,1)  = rmax*dx
            relr(2,1)  = rmax*dy
            relr(3,1)  = rmax*dz


            dr = 10.0

            do while( dr .gt. 1.0d0 )
             dx = rn(0)
             dy = rn(0)
             dz = rn(0)
             dr = dsqrt( dx**2.0 + dy**2.0d0 + dz**2.0 )
            end do


            r(1,idn1)  = rmax*dx + rcmx 
            r(2,idn1)  = rmax*dy + rcmy 
            r(3,idn1)  = rmax*dz + rcmz 

            relr(1,idn1)  = rmax*dx
            relr(2,idn1)  = rmax*dy
            relr(3,idn1)  = rmax*dz


            dr = 10.0

            do while( dr .gt. 1.0d0 )
             dx = rn(0)
             dy = rn(0)
             dz = rn(0)
             dr = dsqrt( dx**2.0 + dy**2.0d0 + dz**2.0 )
            end do

            r(1,idn2)  = rmax*dx + rcmx 
            r(2,idn2)  = rmax*dy + rcmy 
            r(3,idn2)  = rmax*dz + rcmz 

            relr(1,idn2)  = rmax*dx
            relr(2,idn2)  = rmax*dy
            relr(3,idn2)  = rmax*dz

!-----------------------------------------------------------

            ichg(1) = 1
            ichg(idn1) = 0
            ichg(idn2) = 0

            iclst(1) = 0
            iclst(idn1) = 0
            iclst(idn2) = 0

            ininclst(1) = 2
            ininclst(idn1) = 2
            ininclst(idn2) = 2

            numclst(1) = 1
            numclst(idn1) = 1
            numclst(idn2) = 1

	    nmasclst(1) = 3

            icoales(1) = 1
            icoales(idn1) = 1
            icoales(idn2) = 1

            inio(1) = 0
            inio(idn1) = 0
            inio(idn2) = 0
!++++++++++++++++++++++++++++++++++++++++++++++++++++++end++triton+++++

!++start++helium3++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
        else if( peliel .lt. hrang ) then ! tmori 10/20

             dbpass = 3

             if( qh .gt. 0.0 )then
               pexmax0 = einp - qh
             else
               pexmax0 = einp
             end if

             pexmax = min(90.0d0,pexmax0)

             rsq1 = 1000.0d0
             rsq2 = 2000.0d0

             do j = nmaspr+1, nmasta + nmaspr
               if( ipot(j) .eq. 1 ) cycle
               if( iclst(j) .ne. 0 ) cycle
               if( inds(j) .ne. 1 ) cycle
               if( ichg(j) .eq. ichg(1) ) cycle
               if( p(5,j) .lt. 0.10d0 ) cycle

                rsq = ( r(1,1) - r(1,j) )**2
     &              + ( r(2,1) - r(2,j) )**2
     &              + ( r(3,1) - r(3,j) )**2

               rsq0 = sqrt(rsq)

                if( rsq0 .lt. rsq1 )then
                  rsq1 = rsq0
                  idn1 = j
                end if

             end do

             rsq1 = 1000.0d0        
             rsq2 = 20.0d0

             do j = nmaspr + 1, nmasta + nmaspr
               if( ipot(j) .eq. 1 ) cycle
               if( iclst(j) .ne. 0 ) cycle
               if( inds(j) .ne. 1 ) cycle
               if( ichg(j) .ne. ichg(1) ) cycle
               if( p(5,j) .lt. 0.10d0 ) cycle

                rsq = ( r(1,1) - r(1,j) )**2
     &              + ( r(2,1) - r(2,j) )**2
     &              + ( r(3,1) - r(3,j) )**2

               rsq0 = sqrt(rsq)

                if( rsq0 .lt. rsq1 )then
                  rsq1 = rsq0
                  idp1 = j
                end if

             end do
!*********** come here if pauli blocked
1530       continue
!************************************             
!--- Breit Wigner
            ! do i = 1, 2
	      eex = pexmax * rn(0)
	      apex = 10.0d0
	      bpex = -10.0d0

              rrate0 = - 10.0d0

	      do while( apex .gt. bpex )

	       apex = rn(0)

	       eex = pexmax * rn(0)

               wisd = 0.0d0

               fadp =  400.0d0/( eex**2. + 400.0d0 ) !Factor direct pick up for 27Al width 2 MeV
               bpex = fadp !Height of distribution func. at eex only breitwigner
            end do

	      !if(i .eq. 1)
              dpex1 = eex/2.
              !if(i .eq. 2)
              dpex2 = eex/2.
!williams

	      apex = 10.0d0
	      bpex = -10.0d0
	      pexmax =  min(upot*2.,pexmax0) 

	      do while( apex .gt. bpex )

	       apex = rn(0)
	       eex = pexmax * rn(0)

               wisd = ((eex + 1.0d0)/(pexmax + 1.0d0))
          
               if( eex .ge. tfm*2. )then
!                exlim = exp( - (eex - 35.)/10.0 )
                exlim = (upot *2. - eex )/(upot*2. - tfm*2.)
               else
                exlim = 1.0d0
               end if

               bpex = wisd *exlim

	      end do

	      wdpex1 = eex/2.
	      wdpex2 = eex/2.

!---------------------------------------------------------fukuda 18.11.13
!     nucleus energy  = upot-eex
!------------------------------------------------------------------------

           if(  frand .gt. hh2+hh1 )then            
            tc11 = upot - dpex1  !181202 de_dephth add. fukuda
            tc12 = upot - dpex2
            
            if(dpex1+dpex2 .gt. einp)goto 1530

            idad = 10
         elseif(  frand .gt. hh2 )then
          
            tc11 = upot - dpex1
            tc12 = upot - wdpex2 
            if(dpex1+wdpex2 .gt. einp)goto 1530
            idad = 11
         else
              
            tc11 = upot - wdpex1 
            tc12 = upot - wdpex2 
            
            if(wdpex1+wdpex2 .gt. einp)goto 1530
            idad = 12
         end if

          pc01=sqrt(p(1,idn1)**2.+p(2,idn1)**2.+p(3,idn1)**2.) 
          pc02=sqrt(p(1,idp1)**2.+p(2,idp1)**2.+p(3,idp1)**2.)
          pc11 = sqrt(tc11**2.0d0 + 2.0d0 * rmass * tc11 )
          pc12 = sqrt(tc12**2.0d0 + 2.0d0 * rmass * tc12 )
          do k = 1,3
             p(k,idn1) = p(k,idn1)*pc11/pc01
             p(k,idp1) = p(k,idp1)*pc12/pc02
          enddo
               
            p(4,1)=sqrt(rmass**2.0+
     &           p(1,1)**2.0+p(2,1)**2.0+p(3,1)**2.0)
            p(4,idn1)=sqrt(rmass**2.0+
     &           p(1,idn1)**2.0+p(2,idn1)**2.0+p(3,idn1)**2.0)
            p(4,idp1)=sqrt(rmass**2.0+
     &           p(1,idp1)**2.0+p(2,idp1)**2.0+p(3,idp1)**2.0)
             tc1=p(4,1)-p(5,1)
             tc2=p(4,idn1)-p(5,idn1)
             tc3=p(4,idp1)-p(5,idp1)

	     tca = tc1 + tc2 + tc3  ! + 2.*(upot-tfm)
             pca = sqrt( tca**2.0 + 2.0d0 * hmass * tca )/3.0d0
          
             pp0 = dsqrt( p(1,1)**2. + p(2,1)**2. + p(3,1)**2. )

             ome1 = p(1,1) / pp0
             ome2 = p(2,1) / pp0
             ome3 = p(3,1) / pp0

!	     tc2 = einp + thupot - dpex
!             pc2 = sqrt( tc2**2 + 2.0d0 * hmass * tc2 )/3.0d0

             p(1,1) = pca * ome1
             p(2,1) = pca * ome2
             p(3,1) = pca * ome3

             do k = 1, 3
              p(k,idn1)  = p(k,1)
              p(k,idp1)  = p(k,1)
             end do

            p(4,1) = sqrt( rmass **2
     &         + p(1,1)**2 + p(2,1)**2 + p(3,1)**2 )

            p(4,idn1) = p(4,1)
            p(4,idp1) = p(4,1)

!sonoda density 0125----------------------------------------------

	    p(5,1) = rmass
	    p(5,idn1) = rmass
	    p(5,idp1) = rmass

	      do j = 1,3
		pcmclst(j,1) = p(j,1) + p(j,idn1) + p(j,idp1) !\8Fd\90S\82̉^\93\AE\97\CA
	      end do

	      pcmclst(5,1) = hmass
	      pcmclst(4,1) = dsqrt( pcmclst(1,1)**2
     &				+ pcmclst(2,1)**2
     &				+ pcmclst(3,1)**2
     &				+ pcmclst(5,1)**2 )

!-----------------------------------------------------------

            rcmx = (r(1,1) + r(1,idn1) + r(1,idp1))/3.0
            rcmy = (r(2,1) + r(2,idn1) + r(2,idp1))/3.0
            rcmz = (r(3,1) + r(3,idn1) + r(3,idp1))/3.0

	    rcmclst(1,1) = rcmx
	    rcmclst(2,1) = rcmy
	    rcmclst(3,1) = rcmz

            rmax = 1.0

            dr = 10.0

            do while( dr .gt. 1.0d0 )
             dx = rn(0)
             dy = rn(0)
             dz = rn(0)
             dr = dsqrt( dx**2.0 + dy**2.0d0 + dz**2.0 )
            end do

            r(1,1)  = rmax*dx + rcmx
            r(2,1)  = rmax*dy + rcmy
            r(3,1)  = rmax*dz + rcmz

            relr(1,1)  = rmax*dx
            relr(2,1)  = rmax*dy
            relr(3,1)  = rmax*dz

            dr = 10.0

            do while( dr .gt. 1.0d0 )
             dx = rn(0)
             dy = rn(0)
             dz = rn(0)
             dr = dsqrt( dx**2.0 + dy**2.0d0 + dz**2.0 )
            end do

            r(1,idn1)  = rmax*dx + rcmx
            r(2,idn1)  = rmax*dy + rcmy
            r(3,idn1)  = rmax*dz + rcmz

            relr(1,idn1)  = rmax*dx
            relr(2,idn1)  = rmax*dy
            relr(3,idn1)  = rmax*dz

            dr = 10.0

            do while( dr .gt. 1.0d0 )
             dx = rn(0)
             dy = rn(0)
             dz = rn(0)
             dr = dsqrt( dx**2.0 + dy**2.0d0 + dz**2.0 )
            end do

            r(1,idp1)  = rmax*dx + rcmx
            r(2,idp1)  = rmax*dy + rcmy
            r(3,idp1)  = rmax*dz + rcmz

            relr(1,idp1)  = rmax*dx
            relr(2,idp1)  = rmax*dy
            relr(3,idp1)  = rmax*dz

!-----------------------------------------------------------

            ichg(1) = 1
            ichg(idn1) = 0
            ichg(idp1) = 1

            iclst(1) = 0
            iclst(idn1) = 0
            iclst(idp1) = 0

            ininclst(1) = 3
            ininclst(idn1) = 3
            ininclst(idp1) = 3

            numclst(1) = 1
            numclst(idn1) = 1
            numclst(idp1) = 1

	    nmasclst(1) = 3

            icoales(1) = 1
            icoales(idn1) = 1
            icoales(idp1) = 1

            inio(1) = 0
            inio(idn1) = 0
            inio(idp1) = 0

!++++++++++++++++++++++++++++++++++++++++++++++++++++++end++helium3++++

!++start++alpha++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

        else if( peliel .lt. arang ) then

             dbpass = 4

             if( qa .gt. 0.0 )then
               pexmax0 = einp - qa
             else
               pexmax0 = einp
             end if

             pexmax = min(45.0d0,pexmax0)

             rsq1 = 1000.0d0
             rsq2 = 2000.0d0

             do j = nmaspr+1, nmasta + nmaspr
               if( ipot(j) .eq. 1 ) cycle
               if( iclst(j) .ne. 0 ) cycle
               if( inds(j) .ne. 1 ) cycle
               if( ichg(j) .eq. ichg(1) ) cycle
               if( p(5,j) .lt. 0.10d0 ) cycle

                rsq = ( r(1,1) - r(1,j) )**2
     &              + ( r(2,1) - r(2,j) )**2
     &              + ( r(3,1) - r(3,j) )**2

               rsq0 = sqrt(rsq)

                if( rsq0 .lt. rsq1 )then
                  rsq2 = rsq1
                  rsq1 = rsq0
                  idn2 = idn1
                  idn1 = j
                else if( rsq0 .lt. rsq2  )then
                  rsq2 = rsq0
                  idn2 = j
                end if

             end do

             rsq1 = 1000.0d0
             rsq2 = 2000.0d0

             do j = nmaspr + 1, nmasta + nmaspr
               if( ipot(j) .eq. 1 ) cycle
               if( iclst(j) .ne. 0 ) cycle
               if( inds(j) .ne. 1 ) cycle
               if( ichg(j) .ne. ichg(1) ) cycle
               if( p(5,j) .lt. 0.10d0 ) cycle

                rsq = ( r(1,1) - r(1,j) )**2
     &              + ( r(2,1) - r(2,j) )**2
     &              + ( r(3,1) - r(3,j) )**2

               rsq0 = sqrt(rsq)

                if( rsq0 .lt. rsq1 )then
                  rsq1 = rsq0
                  idp1 = j
                end if

             end do

!****** come here if pauli blocked
 1540      continue
!*********************************     
!--- Breit Wigner
             do i = 1, 3
	      eex = pexmax * rn(0)
	      apex = 10.0d0
	      bpex = -10.0d0

              rrate0 = - 10.0d0

	      do while( apex .gt. bpex )

	       apex = rn(0)

	       eex = pexmax * rn(0)

               wisd = 0.0d0

               fadp =  4.0d0/( eex**2 + 4.0d0 ) !Factor direct pick up for 27Al width 2 MeV
               bpex = fadp !Height of distribution func. at eex only breitwigner
            end do

	      if(i .eq. 1)dpex1 = eex
              if(i .eq. 2)dpex2 = eex
              if(i .eq. 3)dpex3 = eex
              
           end do
!---------------------------------------------------------fukuda 18.12.28
!     determin depth deep (using williams state density)
!------------------------------------------------------------------------
!---Williams          
              apex = 10.0d0
              bpex = -10.0d0
              pexmax =  min(upot*3.,pexmax0) !Limit excitation energy 45MeV*3
              do while( apex .gt. bpex )
                 apex = rn(0) !Random number
                 eex = pexmax * rn(0)!*3.0d0 !Excitation energy

                 wisd = (eex/pexmax)**2.0d0
                 
!--------------------------------------------- 20.01.15
            if( eex .ge. tfm*3. )then
!               exlim = exp( - (eex - 35.)/10.0 )
                exlim = (upot *3. - eex )/(upot*3. - tfm*3.)
               else
                exlim = 1.0d0
             end if
!------------------------------------------------------
             !bpex = wisd
             bpex = wisd * exlim
          enddo

          wdpex1 = eex/3.
          wdpex2 = eex/3.
          wdpex3 = eex/3.

!---------------------------------------------------------fukuda 18.11.13
!     nucleus energy  = upot-eex
!------------------------------------------------------------------------
              idad = 999
           if(  frand .gt. ah3+ah2+ah1 )then
            
            tc11 = upot - dpex1  !181202 de_dephth add. fukuda
            tc12 = upot - dpex2
            tc13 = upot - dpex3 
            
            if(dpex1+dpex2+dpex3 .gt. einp)goto 1540

            idad = 10
         elseif(  frand .gt. ah3+ah2 )then
          
            if( nmasta * frand .lt. nchgta)then     
              tc11 = upot - dpex1
              tc12 = upot - dpex2 
              tc13 = upot - wdpex3 
              ! if(dpex1+dpex2+tfm-tc13 .gt. einp)goto 1540           
              if(dpex1+dpex2+wdpex3 .gt. einp)goto 1540 !200115
           else                 !nchg        
              tc11 = upot - dpex1
              tc12 = upot - wdpex2 
              tc13 = upot - dpex3 
!              if(dpex1+dpex2+tfm-tc13 .gt. einp)goto 1540
              if(dpex1+wdpex2+dpex3 .gt. einp)goto 1540
               
               endif !nchgta
               idad = 11
            else if(  frand .gt. ah3 )then
               if( nmasta * frand .lt. nchgta)then       
              tc11 = upot - dpex1 
              tc12 = upot - wdpex2 
              tc13 = upot - wdpex3
!              if(dpex1-tc12-tc13+2.*tfm .gt. einp)goto 1540
              if(dpex1+wdpex2+wdpex3 .gt. einp)goto 1540
               else             !nchg         
              tc11 = upot - wdpex1  
              tc12 = upot - wdpex2 
              tc13 = upot - dpex3 
             
!              if(dpex1-tc12-tc13+2.*tfm .gt. einp)goto 1540
               if(wdpex1+wdpex2+dpex3 .gt. einp)goto 1540                    
               endif!endif!nchgta
               idad = 12
            else
              
            tc11 = upot - wdpex1 
            tc12 = upot - wdpex2 
            tc13 = upot - wdpex3 
            
!            if(tc11+tc12+tc13+einp - 3.*tfm .lt. 0.0d0)goto 1540
           if(wdpex1+wdpex2+wdpex3 .gt. einp)goto 1540
            idad = 13
         end if

          pc01=sqrt(p(1,idn1)**2.+p(2,idn1)**2.+p(3,idn1)**2.) 
          pc02=sqrt(p(1,idn2)**2.+p(2,idn2)**2.+p(3,idn2)**2.)
          pc03=sqrt(p(1,idp1)**2.+p(2,idp1)**2.+p(3,idp1)**2.)
          pc11 = sqrt(tc11**2.0d0 + 2.0d0 * rmass * tc11 )
          pc12 = sqrt(tc12**2.0d0 + 2.0d0 * rmass * tc12 )
          pc13 = sqrt(tc13**2.0d0 + 2.0d0 * rmass * tc13 )
          do k = 1,3
             p(k,idn1) = p(k,idn1)*pc11/pc01
             p(k,idn2) = p(k,idn2)*pc12/pc02
             p(k,idp1) = p(k,idp1)*pc13/pc03
          enddo

!              do j = 1,100
!
!               eex = pexmax*dble(j)/100.0d0
!               gw = (48.6)/dble(nmasta)
!               wisd = ((eex + 1.0)/(pexmax + 1.0))**5.0
!
!               fadp =  gw/( eex**2.0 + gw )!Factor direct pick up
!
!               if( eex .ge. 35. )then
!                exlim = exp( - (eex - 35.)/10.0 )
!               else
!                exlim = 1.0d0
!               end if
!
!               rrate1 = wisd *( fadp**2.0d0 )* exlim
!
!               if( rrate1 .gt. rrate0 )rrate0 = rrate1
!
!              end do
!
!	      eex = pexmax * rn(0)
!	      apex = 10.0d0
!	      bpex = -10.0d0
!
!	      do while( apex .gt. bpex )
!
!	       apex = rn(0)
!	       eex = pexmax * rn(0)
!               gw = (48.6)/dble(nmasta)
!               wisd = ((eex + 1.0)/(pexmax + 1.0))**5.0
!
!               fadp =  gw/( eex**2.0 + gw )!Factor direct pick up
!
!               if( eex .ge. 35. )then
!                exlim = exp( - (eex - 35.)/10.0 )
!               else
!                exlim = 1.0d0
!               end if
!
!               bpex = wisd *( fadp**2.0d0 )* exlim/rrate0
!
!	      end do
!
!	     dpex = eex

!------------------------------------------------------------
!     cal K.E. and momentum
!------------------------------------------------------------
            p(4,1)=sqrt(rmass**2.0+
     &           p(1,1)**2.0+p(2,1)**2.0+p(3,1)**2.0)
            p(4,idn1)=sqrt(rmass**2.0+
     &           p(1,idn1)**2.0+p(2,idn1)**2.0+p(3,idn1)**2.0)
            p(4,idn2)=sqrt(rmass**2.0+
     &           p(1,idn2)**2.0+p(2,idn2)**2.0+p(3,idn2)**2.0)
            p(4,idp1)=sqrt(rmass**2.0+
     &           p(1,idp1)**2.0+p(2,idp1)**2.0+p(3,idp1)**2.0)
             tc1=p(4,1)-p(5,1)
             tc2=p(4,idn1)-p(5,idn1)
             tc3=p(4,idn2)-p(5,idn2)
             tc4=p(4,idp1)-p(5,idp1)

	     tca = tc1 + tc2 + tc3 +tc4 ! + 3.*(upot-tfm)
             pca = sqrt( tca**2.0 + 2.0d0 * amass * tca )/4.0d0
             pp0 = dsqrt( p(1,1)**2. + p(2,1)**2. + p(3,1)**2. )

             ome1 = p(1,1) / pp0
             ome2 = p(2,1) / pp0
             ome3 = p(3,1) / pp0

!	     tc2 = einp + aupot - dpex
!             pc2 = sqrt( tc2**2 + 2.0d0 * amass * tc2 )/4.0d0

             p(1,1) = pca * ome1
             p(2,1) = pca * ome2
             p(3,1) = pca * ome3

             do k = 1, 3
              p(k,idn1)  = p(k,1)
              p(k,idn2)  = p(k,1)
              p(k,idp1)  = p(k,1)
             end do

            p(4,1) = sqrt( rmass **2
     &         + p(1,1)**2 + p(2,1)**2 + p(3,1)**2 )

            p(4,idn1) = p(4,1)
            p(4,idn2) = p(4,1)
            p(4,idp1) = p(4,1)

	    p(5,1) = rmass
	    p(5,idn1) = rmass
	    p(5,idn2) = rmass
	    p(5,idp1) = rmass

	      do j = 1,3
		pcmclst(j,1) = p(j,1) + p(j,idn1) + p(j,idn2)
     &			+ p(j,idp1) !\8Fd\90S\82̉^\93\AE\97\CA
	      end do

	      pcmclst(5,1) = amass
	      pcmclst(4,1) = dsqrt( pcmclst(1,1)**2
     &				+ pcmclst(2,1)**2
     &				+ pcmclst(3,1)**2
     &				+ pcmclst(5,1)**2 )

!-----------------------------------------------------------


            rcmx = (r(1,1) + r(1,idn1) + r(1,idn2)
     &  + r(1,idp1) )/4.0

            rcmy = (r(2,1) + r(2,idn1) + r(2,idn2)
     &  + r(2,idp1) )/4.0

            rcmz = (r(3,1) + r(3,idn1) + r(3,idn2)
     &  + r(3,idp1) )/4.0

	    rcmclst(1,1) = rcmx
	    rcmclst(2,1) = rcmy
	    rcmclst(3,1) = rcmz


            rmax = 1.0

            dr = 10.0


            do while( dr .gt. 1.0d0 )
             dx = rn(0)
             dy = rn(0)
             dz = rn(0)
             dr = dsqrt( dx**2.0 + dy**2.0d0 + dz**2.0 )
            end do

            r(1,1)  = rmax*dx + rcmx
            r(2,1)  = rmax*dy + rcmy
            r(3,1)  = rmax*dz + rcmz

            relr(1,1)  = rmax*dx
            relr(2,1)  = rmax*dy
            relr(3,1)  = rmax*dz


            dr = 10.0

            do while( dr .gt. 1.0d0 )
             dx = rn(0)
             dy = rn(0)
             dz = rn(0)
             dr = dsqrt( dx**2.0 + dy**2.0d0 + dz**2.0 )
            end do

            r(1,idn1)  = rmax*dx + rcmx
            r(2,idn1)  = rmax*dy + rcmy
            r(3,idn1)  = rmax*dz + rcmz

            relr(1,idn1)  = rmax*dx
            relr(2,idn1)  = rmax*dy
            relr(3,idn1)  = rmax*dz

                                  
            dr = 10.0

            do while( dr .gt. 1.0d0 )
             dx = rn(0)
             dy = rn(0)
             dz = rn(0)
             dr = dsqrt( dx**2.0 + dy**2.0d0 + dz**2.0 )
            end do

            r(1,idn2)  = rmax*dx + rcmx
            r(2,idn2)  = rmax*dy + rcmy
            r(3,idn2)  = rmax*dz + rcmz

            relr(1,idn2)  = rmax*dx
            relr(2,idn2)  = rmax*dy
            relr(3,idn2)  = rmax*dz


            dr = 10.0

            do while( dr .gt. 1.0d0 )
             dx = rn(0)
             dy = rn(0)
             dz = rn(0)
             dr = dsqrt( dx**2.0 + dy**2.0d0 + dz**2.0 )
            end do

            r(1,idp1)  = rmax*dx + rcmx
            r(2,idp1)  = rmax*dy + rcmy
            r(3,idp1)  = rmax*dz + rcmz

            relr(1,idp1)  = rmax*dx
            relr(2,idp1)  = rmax*dy
            relr(3,idp1)  = rmax*dz

!-----------------------------------------------------------

!            ichg(1) = 1
!            ichg(idn1) = 0
!            ichg(idn2) = 0
!            ichg(idp1) = 1

            iclst(1) = 0
            iclst(idn1) = 0
            iclst(idn2) = 0
            iclst(idp1) = 0

            ininclst(1) = 4
            ininclst(idn1) = 4
            ininclst(idn2) = 4
            ininclst(idp1) = 4

            numclst(1) = 1
            numclst(idn1) = 1
            numclst(idn2) = 1
            numclst(idp1) = 1

	    nmasclst(1) = 4

            icoales(1) = 1
            icoales(idn1) = 1
            icoales(idn2) = 1
            icoales(idp1) = 1

            inio(1) = 0
            inio(idn1) = 0
            inio(idn2) = 0
            inio(idp1) = 0


            !ipot(1) = 1
            !ipot(idn1) = 1 
            !ipot(idn2) = 1 
            !ipot(idp1) = 1 
        end if
            
!++++++++++++++++++++++++++++++++++++++++++++++++++++++end++alpha++++++

        end if
       end if

      END SUBROUTINE dpickup
      
!**********************************************************************
!*
      SUBROUTINE dpid(iclst,iproj,icoales,ipot,inio,nmaspr,nmasta
     &    ,nchgta,ein,upot,dupot,thupot,aupot,p,r,peex,ichg,inds,tfm
     &    ,radm,numclst,ininclst,nmasclst,pcmclst,rcmclst,relr,idn1
     &    ,idpid)

!*                                                                    *
!*                                                                    *
!*      Purpose:                                                      *
!*                                                                    *
!*            incident particle pickup nucleon                        *
!*                                                                    *
!*                                                                    *
!**********************************************************************

      implicit real(8) (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      real(8), intent(inout) :: p(6,nnn),r(3,nnn),peex,
     &		pcmclst(5,nnn),rcmclst(3,nnn),relr(3,nnn)
      real(8), intent(in) :: upot,dupot,thupot,aupot,ein,tfm,radm

!----------------------------sonoda11/13
      integer, intent(in) :: nmaspr,nmasta,iproj,inds(nnn)
!----------------------------sonoda11/13
      integer, intent(inout) :: iclst(nnn),icoales(nnn),ipot(nnn),
     &     inio(nnn),ichg(nnn),numclst(nnn),ininclst(nnn)
      integer isel(8),dbpass,idn2,idn1,idp1,nmasclst(nnn)

      idpid = 0

      if (iproj .ne. 1 ) return
      if (ein .lt. 450.) return

      dpidp = rn(0)
      sas = 0.01*nmasta + 0.38
      asa = (ein - 700.) / 20.
      fdpid = 0.0137 * exp((ein-580.)**2./2./(110.**2.))
      gdpid = 0.310  * exp( -0.005 * ein )
      pdpid = fdpid / (1.+exp(asa)) + gdpid / (1.+exp(-asa))
      if (dpidp .gt. pdpid*sas) return

        idpid = 1

        i = 0
        j = 0
        k = 0

        dbpass = 0

        idn1 = 0
        idn2 = 0
        idp1 = 0

        rsq1 = 10.0d0
        rsq2 = 20.0d0

        npdp = nchgta !Z
        nndp = nmasta - nchgta !N

        peliel = rn(0)
        frand = rn(0)

        brannp = 0.0
        rrate0 = 0.0d0
        einp = p(4,1) - rmass - upot

        pcmclst = 0.0
        rcmclst = 0.0  
        relr = 0.0

        numclst = 0
        nmasclst = 0.0

 1232  continue  

        rsq1 = 1000.0d0
        rsq2 = 2000.0d0
   
        do j = 2, nmasta+nmaspr
          if( ipot(j) .eq. 1 ) cycle 
          if( iclst(j) .ne. 0 ) cycle
          if( ininclst(j) .ne. 0 ) cycle
          if( inds(j) .ne. 1 ) cycle
          if( ichg(j) .ne. 1 ) cycle
          if( p(5,j) .lt. 0.10d0 ) cycle
                 
    
          rsq = ( r(1,1) - r(1,j) )**2.
     &       + ( r(2,1) - r(2,j) )**2.
     &       + ( r(3,1) - r(3,j) )**2.

          rsq0 = sqrt(rsq)

          if( rsq0 .lt. rsq1 )then
            rsq1 = rsq0 
            idn1 = j
          end if

        end do

        delmass = ddmass+290 
        pxall = p(1,1)! + p(1,idn1)
        pyall = p(2,1)! + p(2,idn1)
        pzall = p(3,1)! + p(3,idn1)
        eall = dsqrt(pxall**2+pyall**2+pzall**2+delmass**2)
        ed0 = ( delmass**2 + (ddmass)**2 - 139.6**2)/2/delmass

        pd0 = dsqrt(ed0**2-ddmass**2)
           
        pall = dsqrt(eall**2-delmass**2)
        gammm = eall / delmass

        vamm = pall / eall
        vx = pxall/eall
        vy = pyall/eall
        vz = pzall/eall

        dithe =0.0
        aithe =0.0
        opp =0.0
        aaaithe = -10.0


        do while(opp .gt. aaaithe)
          aithe  = pi*rn(0)
          opp    = rn(0) 
          aaaithe = dsin(aithe) 
            
        end do

        dithe  =  aithe  
      
        px1 = pd0*dsin(dithe)
        pz1 = pd0*dcos(dithe)
        pv = px1*vx+pz1*vz
        lj = (gammm**2/(gammm+1))*(pv)+gammm*ed0
        px = px1+lj*vx
        pz = pz1+lj*vz           

        pcmclst(1,1) = px
        pcmclst(2,1) = 0.0
        pcmclst(3,1) = pz

        theaa = datan(px/pz)
        theaaa = 180*theaa/pi
               
        do k = 1, 3
          p(k,1)  = pcmclst(k,1)/2.0d0 
          p(k,idn1)  = p(k,1) 
        end do

        pcmclst(5,1) = ddmass

        pabs =dsqrt(pcmclst(1,1)**2.+pcmclst(2,1)**2.+pcmclst(3,1)**2.)
    
        pcmclst(4,1) = dsqrt(ddmass**2.+pabs**2.)+dupot

        edpid = pcmclst(4,1)-pcmclst(5,1)-dupot

        ddithe = dithe/pi*180     

        if(ddithe.gt.0.0.and.ddithe.lt.90.0)then
          
          sigg1  = 40.
          gaws    = -10.0
          gpex   = 10.
          
          do while(gpex.gt.gaws)
            gpex = rn(0)
            xxx = ein*rn(0)
            gaws = exp(-((xxx-edpid)**2.)/2./(sigg1**2.))
          enddo

        else
          
          sigg2  = 100.
          gaws    = -10.0
          gpex   = 10.

          do while(gpex.gt.gaws)

            gpex = rn(0)
            xxx = 670*rn(0) 
            gaws = exp(-((xxx-edpid)**2.)/2./(sigg2**2.))

          enddo

        endif

        if(gaws .gt. gpex) then
          edpid2 = xxx+dupot
        end if

        pdpid2 = (edpid2*(edpid2+2.*ddmass))**0.5
        pdpid = (edpid*(edpid+2.*ddmass))**0.5
        
        pcmclst(1,1) = (pdpid2/pdpid)*px
        pcmclst(2,1) = 0.0
        pcmclst(3,1) = (pdpid2/pdpid)*pz
        pcmclst(4,1) = edpid2 + ddmass
        pcmclst(5,1) = ddmass
       
        do i = 1,4
          p(i,1) = pcmclst(i,1)/2
          p(i,idn1) = p(i,1)
        enddo     
         
        rcmx = r(1,idn1)
        rcmy = r(2,idn1)
        rcmz = r(3,idn1)

        rcmclst(1,1) = rcmx !\8Fd\97z\8Eq\82̏d\90S\82̈ʒu
        rcmclst(2,1) = rcmy
        rcmclst(3,1) = rcmz

        rmax = 1.0
        dr = 10.0

        do while( dr .gt. 1.0d0 )
          dx = rn(0)
          dy = rn(0)
          dz = rn(0)
          dr = dsqrt( dx**2 + dy**2 + dz**2 )
        end do

        r(1,1)  = rmax*dx + rcmx 
        r(2,1)  = rmax*dy + rcmy 
        r(3,1)  = rmax*dz + rcmz 

        relr(1,1)  = rmax*dx
        relr(2,1)  = rmax*dy
        relr(3,1)  = rmax*dz

        dr = 10.0

        do while( dr .gt. 1.0d0 )
          dx = rn(0)
          dy = rn(0)
          dz = rn(0)
          dr = dsqrt( dx**2.0 + dy**2.0d0 + dz**2.0 )
        end do

        r(1,idn1)  = rmax*dx + rcmx 
        r(2,idn1)  = rmax*dy + rcmy 
        r(3,idn1)  = rmax*dz + rcmz 

        relr(1,idn1)  = rmax*dx
        relr(2,idn1)  = rmax*dy
        relr(3,idn1)  = rmax*dz

        ichg(1) = 1
        ichg(idn1) = 0

        iclst(1) = 1
        iclst(idn1) = 1
 
        ininclst(1) = 1
        ininclst(idn1) = 1

        numclst(1) = 1
        numclst(idn1) = 1

        nmasclst(1) = 2

        icoales(1) = 1
        icoales(idn1) = 1

        inio(1) = 0
        inio(idn1) = 0

      END SUBROUTINE dpid

!**********************************************************************
!*
      SUBROUTINE knock_clst(iclst,icoales,ipot,nmaspr,nmasta,nchgta,ein,
     &  upot,dupot,thupot,aupot,p,r,ichg,inds,tfm,numclst,
     &  ininclst,rcmclst,pcmclst,relr,nmasclst,ipknock,iproj)
!      SUBROUTINE knock_clst(iclst,icoales,ipot,nmaspr,nmasta,nchgta,
!     &  upot,dupot,thupot,aupot,p,r,ichg,inds,tfm,numclst,
!     &  ininclst,rcmclst,pcmclst,relr,nmasclst,nknock,nknockp)
!sonoda1204
!*                                                                    *
!*                                                                    *
!*      Purpose:                                                      *
!*                                                                    *
!*               Knock out a real cluster around fermi surface        *
!*                                                                    *
!*                                                                    *
!**********************************************************************

      use NGSDATAMOD, only : energm
      implicit real(8) (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      real(8), intent(inout) :: p(6,nnn),r(3,nnn),pcmclst(5,nnn),
     &		rcmclst(3,nnn),relr(3,nnn)
      real(8), intent(in) :: upot,dupot,thupot,aupot,tfm,ein
      real(8) pcm(3),beta(3),kcpd,kcpt,kcph,kcpa
!----------------------------sonoda11/13
      integer, intent(in) :: nmaspr,nmasta,nchgta
!----------------------------sonoda11/13
      integer, intent(inout) :: iclst(nnn),icoales(nnn),ipot(nnn),
     &     ichg(nnn),inds(nnn),numclst(nnn),ininclst(nnn),
     &	   nmasclst(nnn),ipknock(nnn)

      iknock = 0
      i = 0
      j = 0
      ikn1 = 0
      ikn2 = 0
      ikp1 = 0
      ikp2 = 0
!---------------------------------------------------------------------
!     Probability for knockout(deuteron) depending on incident E
!     yamaguchi
!---------------------------------------------------------------------
      !  if( ein .lt. 55.0 )then
      !     dcrs = 0.0          ! p,d
      !  elseif( ein .ge. 55.0 .and. ein .lt. 79.0 )then
      !     dcrs = 0.0
       ! elseif( ein .ge. 79.0 .and. ein .lt. 195.0 )then
       !    dcrs = 0.25
       ! elseif( ein .ge. 195.0 .and. ein .lt. 346.0 )then
       !    dcrs = 0.15
       ! else
       !    dcrs = 0.0
       ! end if
      

      if (iproj .eq. 5) then
        dcrs = 0.2d0
        tcrs = 0.1d0
        hcrs = 0.05d0
        acrs = 0.01

        tcrs = tcrs * 4.!nmasclst(1)
        hcrs = hcrs * 4.!nmasclst(1)
        acrs = acrs * 4.!nmasclst(1)

      else
        dcrs = 0.016+0.17/(1+exp((ein-200)/15))
        if( ein .lt. 70.0 ) then
          tcrs = 0.0
        else 
          tcrs = 0.01
        end if
           
        hcrs = tcrs  

        if( ein .lt. 70.0 )then
           acrs = 0.0d0
        else if( ein .lt. 100.0 )then
           acrs = 0.181d0
        else if( ein .lt. 170.0 )then
           acrs = 0.3d0
        else if( ein .lt. 250.0 )then
           acrs = 0.05d0
        else
           acrs = 2.5d-3
        end if
      end if

        npdp = nchgta
        nndp = nmasta - nchgta
	kcpass = 0
	pbcz = 0.0d0

        pra = nmaspr

!+++++++++++++sort grand energy+++++++++++++++++++++++++++++


        pferm1 = 1000.0d0
        pferm2 = 2000.0d0

        einp = p(4,1) - rmass - upot
        pexmax = einp

           do j = nmaspr + 1 , nchgta + nmaspr

           ! if these condition is not reached ,
           ! the particle is removed from this sort
           if( ipot(j) .eq. 1 )cycle
           if( iclst(j) .ne. 0 )cycle
           if( inds(j) .ne. 1 )cycle
           if( p(5,j) .lt. 0.10d0 )cycle
           if( ichg(j) .ne. 1 )cycle

            rsq = ( r(1,j) - r(1,1) )**2
     &  + ( r(2,j) - r(2,1) )**2 + ( r(3,j) - r(3,1) )**2

            rsq = dsqrt(rsq)

              if( rsq .lt. pferm1  )then

               ikp2 = ikp1
               ikp1 = j
               pferm2 = pferm1
               pferm1 = rsq

              else if( rsq .lt. pferm2 )then

               ikp2 = j
               pferm2 = rsq

              end if

            end do

        pferm1 = 1000.0d0
        pferm2 = 2000.0d0

           !neutron
!           do j = nchgta + nmaspr + 1, nmasta + nmaspr - 1
           do j = nchgta + nmaspr + 1, nmasta + nmaspr

           if( ipot(j) .eq. 1 )cycle
           if( iclst(j) .ne. 0 )cycle
           if( inds(j) .ne. 1 )cycle
           if( p(5,j) .lt. 0.10d0 )cycle
           if( ichg(j) .ne. 0 )cycle

            rsq = ( r(1,j) - r(1,1) )**2
     &  + ( r(2,j) - r(2,1) )**2 + ( r(3,j) - r(3,1) )**2

            rsq = dsqrt(rsq)

              if( rsq .lt. pferm1  )then

               ikn2 = ikn1
               ikn1 = j
               pferm2 = pferm1
               pferm1 = rsq

              else if( rsq .lt. pferm2 )then

               ikn2 = j
               pferm2 = rsq

              end if

           end do


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!         CALL setup ! call gem's data

          amcta  = energmelf(nchgta,nmasta)

          amcem  = energmelf(1,2)
          amcrs  = energmelf(nchgta - 1,nmasta - 2)
          qd = amcem + amcrs - amcta

          amcem  = energmelf(1,3)
          amcrs  = energmelf(nchgta - 1,nmasta - 3)
          qt = amcem + amcrs - amcta

          amcem  = energmelf(2,3)
          amcrs  = energmelf(nchgta - 2,nmasta - 3)
          qh = amcem + amcrs - amcta

          amcem  = energmelf(2,4)
          amcrs  = energmelf(nchgta - 2,nmasta - 4)
          qa = amcem + amcrs - amcta


          if( einp .lt. qd )dcrs = 0.0
          if( einp .lt. qt )tcrs = 0.0
          if( einp .lt. qh )hcrs = 0.0
          if( einp .lt. qa )acrs = 0.0

          drang = dcrs
          trang = drang + tcrs
          hrang = trang + hcrs
          arang = hrang + acrs

!##############product a cluster near fermi surface#########

	frand = rn(0)

        if( frand .lt. drang )then
!          nknock = nknock + 1   !yamaguchi used to count total knockout event
          iknock = 1
          pexmax0 = einp - qd

          pexmax = min(pexmax0,50.0)

	  eex = pexmax * rn(0)
	  apex = 10.0d0
	  bpex = -10.0d0

          do iex = 1, 2
             
             do while( apex .gt. bpex )

                apex = rn(0) !Random number

                eex = pexmax * rn(0) !Excitation energy
!         wisd = (eex/pexmax)**3.0d0
!         wisd = eex/pexmax
                wisd = eex*1.0d-2
                fadp = 100.0d0/(eex**2 + 100.0d0)
!         bpex = wisd
                bpex = wisd + fadp

             end do

             apex = 10.0
             bpex = -10.0
             if( iex .eq. 1 ) dpex = eex
             if( iex .eq. 2 ) dpex2 = eex

!            dpex = eex
             
          end do
          pnc = 2.0d0

        else if( frand .lt. trang )then

         iknock = 2

         pexmax0 = einp - qt

         pexmax = min(pexmax0,50.)

	 eex = pexmax * rn(0)
	 apex = 10.0d0
	 bpex = -10.0d0

	 do while( apex .gt. bpex )

	  apex = rn(0)!Random number

	  eex = pexmax * rn(0)!Excitation energy
          wisd = ((eex - 1.)/(pexmax - 1.))**5.0d0
          bpex = wisd

	 end do

         pnc = 3.0d0
	 dpex = eex

        else if( frand .lt. hrang  )then

         iknock = 3

         pexmax0 = einp - qh

         pexmax = min(pexmax0,50.)

	 eex = pexmax * rn(0)
	 apex = 10.0d0
	 bpex = -10.0d0

	 do while( apex .gt. bpex )

	  apex = rn(0)!Random number

	  eex = pexmax * rn(0)!Excitation energy
          wisd = ((eex - 1.)/(pexmax - 1.))**5.0d0
          bpex = wisd

	 end do

         pnc = 3.0d0
	 dpex = eex

        else if( frand .lt. arang  )then
         iknock = 4

         pexmax0 = einp - qa

         coul =  cbar(dble(nchgta),dble(nmasta),2,4)

         pexmax1 = min(pexmax0,tfm*4.)
!         pexmax = min(pexmax0,50.)
         pexmax = min(pexmax1,einp-qa-coul)

!	 eex = pexmax * rn(0)

         do i = 1, 4
            apex = 10.0d0
            bpex = -10.0d0

            do while( apex .gt. bpex )

               apex = rn(0)     !Random number

               eex = pexmax * rn(0) !Excitation energy
               wisd = ((eex + 1.)/(pexmax + 1.))**3.0d0
               bpex = wisd

            end do

            if(i .eq. 1)dpex  = eex/4.0d0
            if(i .eq. 1)dpex2 = eex/4.0d0
            if(i .eq. 1)dpex3 = eex/4.0d0
            if(i .eq. 1)dpex4 = eex/4.0d0
            exit
         end do
         pnc = 1.0d0
!	  dpex = eex

        end if

!###########################################################

!***********************************************************

      if( iknock .ne. 0 )then
!***********************************************************

!------------------------start collision-------------------!

            e1   = p(4,1)
            em1  = p(5,1)

            px1  = p(1,1)
            py1  = p(2,1)
            pz1  = p(3,1)
            iz1  = ichg(1)

!           e2   = tfm + rmass - bpex/pnc
            e2   = tfm + rmass - dpex/pnc
            em2  = rmass
            p2 = dsqrt( e2**2.0d0 - rmass**2.0d0 )

            p02 = dsqrt( p(1,ikn1)**2.0 + p(2,ikn1)**2.0
     &    + p(3,ikn1) **2.0d0 )

            px2  = p2*p(1,ikn1)/p02
            py2  = p2*p(2,ikn1)/p02
            pz2  = p2*p(3,ikn1)/p02
            iz2  = ichg(ikn1)

            etot12  = e1 + e2

            beta(1) = ( px1 + px2 ) / etot12
            beta(2) = ( py1 + py2 ) / etot12
            beta(3) = ( pz1 + pz2 ) / etot12

            betasq  = beta(1)**2 + beta(2)**2 + beta(3)**2
            gamma   = 1.0d0 / dsqrt( 1.0d0 - betasq )

!-----------------------------------------------------------------------
!             transformation of momenta
!-----------------------------------------------------------------------

            p1beta = px1 * beta(1) + py1 * beta(2) + pz1 * beta(3)

            transf = gamma
     &           * ( gamma * p1beta / ( gamma + 1.0 ) - e1 )

            pcm(1) = px1 + beta(1) * transf
            pcm(2) = py1 + beta(2) * transf
            pcm(3) = pz1 + beta(3) * transf
            pcm2   = pcm(1)**2 + pcm(2)**2 + pcm(3)**2
            prcm   = sqrt( pcm2 )

            s    = (  e1 +  e2 )**2 - ( px1 + px2 )**2
     &           - ( py1 + py2 )**2
     &           - ( pz1 + pz2 )**2

            srt  = sqrt(s)

            cutoff = em1 + em2 + 20.0
            bcmax  = 1.323142
            sig    = 55.0

!-----------------------------------------------------------------------
            rm2  = rmass**2
            pr   = prcm
            c2   = pcm(3) / pr

            csrt = srt - cutoff
            pri = prcm
            prf = sqrt( 0.25*(srt)**2 - rm2 )

            asrt = srt - em1 - em2
            pra  = prcm
!-------------------------------------------------------------sawada8/27
!           MeV/GeV
            srt  = srt / 1000.0
            asrt = asrt / 1000.0
            pra  = pra / 1000.0
            prf  = prf / 1000.0
            pri  = pri / 1000.0
            csrt = csrt / 1000.0
!-----------------------------------------------------------------------
            plsq = (srt)**4 / ( 4.0 * (rmass/1000.0)**2 )
     &     - (srt)**2
!-----------------------------------------------------------------------
            pl    = 0.0d0
            pcc   = 0.0d0

            if( plsq .gt. 0.0 ) pl = sqrt( plsq )
!-----------------------------------------------------------------------
!     2019.7.10 Comment out because difm is used only here
!              if( difm .ne. 0 )then
!
!               emm1 = em1/1000.0d0
!               emm2 = em2/1000.0d0
!
!               pcsq =  (srt - emm1 - emm2) *
!     &          (srt - emm1 + emm2) *
!     &          (srt + emm1 - emm2) *
!     &          (srt + emm1 + emm2)
!
!                if( pcsq .gt. 0.0 ) pcc = sqrt( pcsq )
!
!                pc = ( 1.0d0/( 2.0d0*srt ) ) * pcc
!                pl = pc*srt/emm2
!              end if
!-----------------------------------------------------------------------
            el   = -rmass + dsqrt( (pl*1000.0)**2 + rmass**2 ) 

!-----------------------------------------------------------------------
!           New Cugnon
!-----------------------------------------------------------------------
            ta  = -2.0 * (pra)**2
            x   = rn(0)
!-----------------------------------------------------------------------
!           collision
!-----------------------------------------------------------------------

                  if( pl .le. 0.225 ) then
                     c1 = 1.0 - 2.0 * x
                  else
                     if( pl .le. 0.6d0 ) then
                        a = 6.2d0 * ( pl - 0.225d0 ) / 0.375d0
                     else if( pl .le. 1.6 ) then
                        a =  - 1.63d0 * pl + 7.16d0
                     else if( pl .le. 2.0d0 ) then
                        a = 5.5d0 * pl**8 / ( 7.7d0 + pl**8 )
                     else
                        a = 5.34d0 + 0.67d0 * ( pl - 2.0d0 )
                     end if
                     tt1 = log( ( 1.0d0 - x )
     &                    * exp( 2.0d0 * a * ta ) + x ) / a
                     bprob = 1.0d0
                     if( pl .gt. 0.8d0 ) bprob = 0.64d0 / pl / pl
                     bprob = bprob / ( 1.0d0 + bprob )
                     if( rn(0) .gt. bprob ) then
                        c1 = 1.0d0 - tt1 / ta
                     else
                        c1 = - 1.0d0 + tt1 / ta
                     end if
                     if( abs(c1) .gt. 1.0d0 ) c1 = 1.0 - 2.0 * x
                  end if

            t1  = 2.0 * pi * rn(0)

            if( pcm(1) .eq. 0.0 .and. pcm(2) .eq. 0.0 ) then
               t2 = 0.0
            else
               t2 = atan2(pcm(2),pcm(1))
            end if

            s1   = dsqrt( 1.0 - c1**2 )
            s2  =  dsqrt( 1.0 - c2**2 )
            ct1  = dcos(t1)
            st1  = dsin(t1)
            ct2  = dcos(t2)
            st2  = dsin(t2)
            ss   = c2 * s1 * ct1  +  s2 * c1

            pcm(1) = pr * ( ss*ct2 - s1*st1*st2 )
            pcm(2) = pr * ( ss*st2 + s1*st1*ct2 )
            pcm(3) = pr * ( c1*c2  - s1*s2 *ct1 )

            testpcm = pcm(1) ** 2 + pcm(2) **2 + pcm(3) ** 2

            e1cm    = dsqrt( em1**2 + pcm(1)**2 + pcm(2)**2 
     &  + pcm(3)**2 )
            p1beta  = pcm(1)*beta(1) + pcm(2)*beta(2) 
     &  + pcm(3)*beta(3)
         transf  = gamma * ( gamma * p1beta / (gamma 
     &  + 1.0d0) + e1cm )

         px1 = beta(1) * transf + pcm(1)
         py1 = beta(2) * transf + pcm(2)
         pz1 = beta(3) * transf + pcm(3)
         e1 = dsqrt( rmass**2. + px1**2. + py1**2. + pz1**2. )

         e2cm    = dsqrt(em2**2 + pcm(1)**2 + pcm(2)**2 
     &  + pcm(3)**2)

         transf  = gamma * (-gamma * p1beta / (gamma 
     &  + 1.0) + e2cm)

         px2 = beta(1) * transf - pcm(1)
         py2 = beta(2) * transf - pcm(2)
         pz2 = beta(3) * transf - pcm(3)
         e2 = dsqrt( rmass**2. + px2**2. + py2**2. + pz2**2. )

!        write(*,*)dble(difm) * 45.0d0

         if( iknock .eq. 1 )then

!             ptot = e2 + rmass + tfm 
!    &                - (1.0 - 1.0/pnc)*dpex
            ptot = e2 + rmass + tfm - dpex2/pnc
            trdel = 0!1/(1+exp((ptot-650)/70))

         else if( iknock .eq. 2 )then

              ptot = e2 + 2.0*(rmass + tfm)
     &                 - (1.0 - 1.0/pnc)*dpex

         else if( iknock .eq. 3 )then

              ptot = e2 + 2.0*(rmass + tfm)
     &                 - (1.0 - 1.0/pnc)*dpex

         else if( iknock .eq. 4 )then

              ptot = e2 + 3.0*(rmass + tfm)
!     &                 - (1.0 - 1.0/pnc)*dpex
     &                 - (dpex2+dpex3+dpex4)
        end if

        if( iknock .eq. 1 )then !d
         cge2 = ptot - ddmass
         fe2 = tfm * 2.0
        else if( iknock .eq. 2 )then !t
         cge2 = ptot - tmass
         fe2 = tfm * 3.0
        else if( iknock .eq. 3 )then !3He
         cge2 = ptot - hmass
         fe2 = tfm * 3.0
        else if( iknock .eq. 4 )then !a
         cge2 = ptot - amass
!         fe2 = tfm * 4.0
         fe2 = upot *4.0 + cbar(dble(nchgta),dble(nmasta),2,4) ! a !tfm ->  coulomb 190627        
         end if

         cge1 = e1 - rmass
         fe1 =  tfm * 1.0 
       ! Pauli Blocking
       if( cge1 .gt. fe1 .and. cge2 .gt. fe2  )then
!         nknockp = nknockp + 1 !yamaguchi used to count actual knockout event
           p(1,1)     = px1
           p(2,1)     = py1
           p(3,1)     = pz1
           p(4,1)     = e1
           p(5,1)     = rmass

         if( iknock .eq. 1 )then ! deuteron

                     icoales(ikn1) = 1
                     icoales(ikp1) = 1

                     p(4,ikn1) = ptot/2.0
                     p(5,ikn1) = rmass

                     et2 = p(4,ikn1) - rmass
                     pabs = dsqrt( et2*(et2 + 2.0d0*rmass) )
          if(rn(0).lt.trdel)then
            p(1,ikn1) = px2
            p(2,ikn1) = py2
            p(3,ikn1) = pz2
          else
!                     xa = px2
!                     ya = py2
!                     za = pz2

                     xa = px2 + p(1,ikp1)
                     ya = py2 + p(2,ikp1)
                     za = pz2 + p(3,ikp1)

                     pclst0 = sqrt( xa**2 + ya**2 + za**2 )

                     p(1,ikn1)  = xa / pclst0 * pabs
                     p(2,ikn1)  = ya / pclst0 * pabs
                     p(3,ikn1)  = za / pclst0 * pabs

                     do k = 1,4
                      p(k,ikp1) = p(k,ikn1)
                     end do

		     do j = 1,3
			pcmclst(j,ikn1) = p(j,ikn1) + p(j,ikp1)
		     end do

		     pcmclst(5,ikn1) = ddmass
		     pcmclst(4,ikn1) = dsqrt( pcmclst(1,ikn1)**2
     &				+ pcmclst(2,ikn1)**2
     &				+ pcmclst(3,ikn1)**2
     &				+ pcmclst(5,ikn1)**2 )

                     ipknock(1) = 1
		     ininclst(ikn1) = 1
		     ininclst(ikp1) = 1

                     numclst(ikn1) = ikn1
                     numclst(ikp1) = ikn1

		     nmasclst(ikn1) = 2

                     ipot(ikn1) = 1
                     ipot(ikp1) = 1
          endif
 
         else if( iknock .eq. 2 )then ! Triton

                     icoales(ikn1) = 1
                     icoales(ikn2) = 1
                     icoales(ikp1) = 1

                     p(4,ikn1) = ptot/3.0d0
                     p(5,ikn1) = rmass

                     et2 = p(4,ikn1) - rmass
                     pabs = dsqrt( et2*(et2 + 2.0d0*rmass) )

                     xa = px2 + p(1,ikn2) + p(1,ikp1)
                     ya = py2 + p(2,ikn2) + p(2,ikp1)
                     za = pz2 + p(3,ikn2) + p(3,ikp1)

                     pclst0 = sqrt( xa**2 + ya**2 + za**2 )

                     p(1,ikn1)  = xa / pclst0 * pabs
                     p(2,ikn1)  = ya / pclst0 * pabs
                     p(3,ikn1)  = za / pclst0 * pabs

                     do k = 1,4
                      p(k,ikn2) = p(k,ikn1)
                      p(k,ikp1) = p(k,ikn1)
                     end do

		     do j = 1,3
			pcmclst(j,ikn1) = p(j,ikn1) + p(j,ikn2)
     &				+ p(j,ikp1)
		     end do

		     pcmclst(5,ikn1) = tmass
		     pcmclst(4,ikn1) = dsqrt( pcmclst(1,ikn1)**2
     &				+ pcmclst(2,ikn1)**2
     &				+ pcmclst(3,ikn1)**2
     &				+ pcmclst(5,ikn1)**2 )


                     ipknock(1) = 1
                     ininclst(ikn1) = 2
                     ininclst(ikn2) = 2
                     ininclst(ikp1) = 2

                     numclst(ikn1) = ikn1
                     numclst(ikn2) = ikn1
                     numclst(ikp1) = ikn1

		     nmasclst(ikn1) = 3

                     ipot(ikn1) = 1
                     ipot(ikn2) = 1
                     ipot(ikp1) = 1

         else if( iknock .eq. 3 )then ! Helium3

                     icoales(ikn1) = 1
                     icoales(ikp1) = 1
                     icoales(ikp2) = 1

                     p(4,ikn1) = ptot/3.0d0
                     p(5,ikn1) = rmass

                     et2 = p(4,ikn1) - rmass
                     pabs = dsqrt( et2*(et2 + 2.0d0*rmass) )

                     xa = px2 + p(1,ikp1) + p(1,ikp2)
                     ya = py2 + p(2,ikp1) + p(2,ikp2)
                     za = pz2 + p(3,ikp1) + p(3,ikp2)

                     pclst0 = sqrt( xa**2 + ya**2 + za**2 )

                     p(1,ikn1)  = xa / pclst0 * pabs
                     p(2,ikn1)  = ya / pclst0 * pabs
                     p(3,ikn1)  = za / pclst0 * pabs

                     do k = 1,4
                      p(k,ikp1) = p(k,ikn1)
                      p(k,ikp2) = p(k,ikn1)
                     end do

		     do j = 1,3
			pcmclst(j,ikn1) = p(j,ikn1) + p(j,ikp1)
     &				+ p(j,ikp2)
		     end do

		     pcmclst(5,ikn1) = hmass
		     pcmclst(4,ikn1) = dsqrt( pcmclst(1,ikn1)**2
     &				+ pcmclst(2,ikn1)**2
     &				+ pcmclst(3,ikn1)**2
     &				+ pcmclst(5,ikn1)**2 )


                     ipknock(1) = 1
                     ininclst(ikn1) = 3
                     ininclst(ikn2) = 3
                     ininclst(ikp1) = 3

                     numclst(ikn1) = ikn1
                     numclst(ikp1) = ikn1
                     numclst(ikp2) = ikn1

		     nmasclst(ikn1) = 3

                     ipot(ikn1) = 1
                     ipot(ikp1) = 1
                     ipot(ikp2) = 1

         else if( iknock .eq. 4 )then ! Alpha

                     icoales(ikn1) = 1
                     icoales(ikn2) = 1
                     icoales(ikp1) = 1
                     icoales(ikp2) = 1

                     p(4,ikn1) = ptot/4.0d0
!                     p(5,ikn1) = amass
                     p(5,ikn1) = rmass

                     et2 = p(4,ikn1) - rmass
                     pabs = dsqrt( et2*(et2 + 2.0d0*rmass) )

                     xa = px2! + p(1,ikn2) + p(1,ikp1) + p(1,ikp2)
                     ya = py2! + p(2,ikn2) + p(2,ikp1) + p(2,ikp2)
                     za = pz2! + p(3,ikn2) + p(3,ikp1) + p(3,ikp2)

                     pclst0 = sqrt( xa**2 + ya**2 + za**2 )

                     p(1,ikn1)  = xa / pclst0 * pabs
                     p(2,ikn1)  = ya / pclst0 * pabs
                     p(3,ikn1)  = za / pclst0 * pabs

                     do k = 1,4
                      p(k,ikn2) = p(k,ikn1)
                      p(k,ikp1) = p(k,ikn1)
                      p(k,ikp2) = p(k,ikn1)
                     end do

                     ppp0 = 0

		     do j = 1,3
			pcmclst(j,ikn1) = p(j,ikn1) + p(j,ikn2)
     &				+ p(j,ikp1) + p(j,ikp2)
                        ppp0 = ppp0 + pcmclst(j,ikn1)**2.
		     end do

		     pcmclst(5,ikn1) = amass
		     pcmclst(4,ikn1) = dsqrt( pcmclst(1,ikn1)**2
     &				+ pcmclst(2,ikn1)**2
     &				+ pcmclst(3,ikn1)**2
     &				+ pcmclst(5,ikn1)**2 )
!     &                          + upot*4.
 
                     ppp0 = sqrt(ppp0)
                     ppp1 = sqrt(pcmclst(4,ikn1)**2. -
     &                                    pcmclst(5,ikn1)**2.)

                     do j = 1,3
                       pcmclst(j,ikn1) = pcmclst(j,ikn1) * ppp1/ppp0
                     end do
                     

                     ininclst(ikn1) = 4
                     ininclst(ikn2) = 4
                     ininclst(ikp1) = 4
                     ininclst(ikp2) = 4

                     numclst(ikn1) = ikn1
                     numclst(ikn2) = ikn1
                     numclst(ikp1) = ikn1
                     numclst(ikp2) = ikn1

		     nmasclst(ikn1) = 4

                     ipknock(1) = 1
                     ipot(ikn1) = 1
                     ipot(ikn2) = 1
                     ipot(ikp1) = 1
                     ipot(ikp2) = 1


         end if

		     rcmclst(1,ikn1) = r(1,ikn1)
		     rcmclst(2,ikn1) = r(2,ikn1)
		     rcmclst(3,ikn1) = r(3,ikn1)

         if( ikn1 .ne. 0 )then ! deuteron


                     rdis = 3.0d0
                     psqr = 10.0

                      do while( psqr .gt. 1.0 )

                      rx = 1.0 - 2.0 * rn(0)
                      ry = 1.0 - 2.0 * rn(0)
                      rz = 1.0 - 2.0 * rn(0)
                      psqr = rx*rx + ry*ry + rz*rz

                      end do

                      r(1,ikn1) = r(1,ikn1) + rdis * rx
                      r(2,ikn1) = r(2,ikn1) + rdis * ry
                      r(3,ikn1) = r(3,ikn1) + rdis * rz

                      relr(1,ikn1) = rdis * rx
                      relr(2,ikn1) = rdis * ry
                      relr(3,ikn1) = rdis * rz

         end if

         if( ikn2 .ne. 0 )then

                     rdis = 2.0d0
                     psqr = 10.0

                      do while( psqr .gt. 1.0 )

                      rx = 1.0 - 2.0 * rn(0)
                      ry = 1.0 - 2.0 * rn(0)
                      rz = 1.0 - 2.0 * rn(0)
                      psqr = rx*rx + ry*ry + rz*rz

                      end do

                      r(1,ikn2) = rcmclst(1,ikn1) + rdis * rx
                      r(2,ikn2) = rcmclst(2,ikn1) + rdis * ry
                      r(3,ikn2) = rcmclst(3,ikn1) + rdis * rz

                      relr(1,ikn2) = rdis * rx
                      relr(2,ikn2) = rdis * ry
                      relr(3,ikn2) = rdis * rz

         end if

         if( ikp1 .ne. 0 )then

                     rdis = 2.0d0
                     psqr = 10.0

                      do while( psqr .gt. 1.0 )

                      rx = 1.0 - 2.0 * rn(0)
                      ry = 1.0 - 2.0 * rn(0)
                      rz = 1.0 - 2.0 * rn(0)
                      psqr = rx*rx + ry*ry + rz*rz

                      end do

!                      r(1,ikp1) = r(1,ikn1) + rdis * rx
!                      r(2,ikp1) = r(2,ikn1) + rdis * ry
!                      r(3,ikp1) = r(3,ikn1) + rdis * rz
                      r(1,ikp1) = rcmclst(1,ikn1) + rdis * rx
                      r(2,ikp1) = rcmclst(2,ikn1) + rdis * ry
                      r(3,ikp1) = rcmclst(3,ikn1) + rdis * rz

                      relr(1,ikp1) = rdis * rx
                      relr(2,ikp1) = rdis * ry
                      relr(3,ikp1) = rdis * rz

         end if

         if( ikp2 .ne. 0 )then

                     rdis = 2.0d0
                     psqr = 10.0

                      do while( psqr .gt. 1.0 )
                      rx = 1.0 - 2.0 * rn(0)
                      ry = 1.0 - 2.0 * rn(0)
                      rz = 1.0 - 2.0 * rn(0)
                      psqr = rx*rx + ry*ry + rz*rz
                      end do

                      r(1,ikp2) = rcmclst(1,ikn1) + rdis * rx
                      r(2,ikp2) = rcmclst(2,ikn1) + rdis * ry
                      r(3,ikp2) = rcmclst(3,ikn1) + rdis * rz

                      relr(1,ikp2) = rdis * rx
                      relr(2,ikp2) = rdis * ry
                      relr(3,ikp2) = rdis * rz

         end if

        else

        end if

       end if
!------------------------end collision----------------------!

      END SUBROUTINE knock_clst

!**********************************************2018/08/28ishoku watanabe
!*                                                                     *
!      SUBROUTINE tunneff( ipot,upot,ucpot,p,rt00,
!     &  	nchgpr,nmaspr,nchgta,nmasta,peex,iclst,ichg )
      SUBROUTINE tunneff( ipot,upot,ucpot,p,rt00,
     &  	nchgpr,nmaspr,nchgta,nmasta,peex,iclst,ichg,ininclst )
!     SUBROUTINE tunneff( ipot,upot,ucpot,p,rt00,
!    &  	nchgpr,nmaspr,nchgta,nmasta,peex,iclst,ichg,itime,inds )
!*                                                                     *
!*                                                                     *
!*        Purpose:                                                     *
!*                                                                     *
!*              to determinate tunnel-effect-adjusted ipot             *
!                   : the particle is above the potential or not       *
!*                                                                     *
!*	     					2013/12/29yamada       *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      integer, intent(in)    :: nchgpr,nmaspr,nchgta,nmasta,iclst(nnn),
     &     ininclst(nnn)
      real*8, intent(in)    :: rt00,upot,ucpot,peex
      integer, intent(inout) :: ipot(nnn),ichg(nnn)
      real*8, intent(inout) :: p(6,nnn)
      integer inds(nnn)
!-----------------------------------------------------------------------

	do i = 1, nmaspr+nmasta

!	 if( iclst(i).ne.0 .or. ichg(i).ne.1 ) cycle
	 if( iclst(i).ne.0 .or. ininclst(i).ne.0 .or. ichg(i).ne.1 ) cycle
         
	 if( p(4,i) - p(5,i) .gt. upot )then

!	   if( nchgta .lt. 50 )then

!   		if( p(4,i) - p(5,i) .gt. ucpot )then
!	           ipot(i) = 1
!		else
!		   ipot(i) = 0
!	        end if

!	   else

! 	      barcentr = upot + coul(nchgta,rt00) * 3.0d0 !centrifugal force barrier. baria no takasa.
!   		if( p(4,i) - p(5,i) .gt. barcentr )then
!-----------------------------------------------------------------------
!                  above the modified Coulomb potential
!-----------------------------------------------------------------------

!	         ipot(i) = 1

!		else
!-----------------------------------------------------------------------
!                  under the modified Coulomb potential
!                            tunnel effect
!-----------------------------------------------------------------------

!		   coulba = barcentr - upot	!modified Coulomb barrier
		   eexit = p(4,i) - p(5,i)	!intranuclear kinetic energy
		   ueexit = eexit - upot	!extranuclear kinetic energy
!		   coulcut = 7.0
!		   coulcut = 0.0166 * nmasta + 3.7252
!               coulcut = cgbar(dble(nchgta),dble(nmasta),1,1,peex)
!	write(*,*)coulcut,peex
!		   if( ueexit .gt. coulcut )then

	            pcoul = rn(0)
                      nmascl=1
                      nchgcl=1
                  CALL gamow( rt00,nchgpr,nmaspr,nchgta,nmasta,
     &		            upot,nmascl,nchgcl,barcentr,eexit,ueexit,transm )

	            if( pcoul .le. transm )then
		         ipot(i) = 1

	            else
	               ipot(i) = 0

		    end if

!                   else

!		       ipot(i) = 0

!	           end if

!		end if

!          end if

	 else

             ipot(i) = 0

	 end if

	end do

      END SUBROUTINE tunneff

!******************************************************2018/8/28watanabe
!*                                                                     *
      SUBROUTINE gamow(rt00,nchgpr,nmaspr,nchgta,nmasta,upot,
     &     		nmascl,nchgcl,barcentr,eexit,ueexit,transm )
!*                                                                     *
!*                                                                     *
!*        Purpose:                                                     *
!*                                                                     *
!*              to calculate the Gamow transmission factor             *
!*                                                                     *
!*						2013/7/9yamada         *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      integer, intent(in)   :: nchgpr,nmaspr,nchgta,nmasta,nmascl
      real*8, intent(in)    :: rt00,upot,eexit,ueexit,barcentr
      real*8, intent(out)   :: transm

      real(8) coulom(6)      
      coulom(1) = cbar(dble(nchgta),dble(nmasta),1,1) ! p
      coulom(2) = cbar(dble(nchgta),dble(nmasta),1,2) ! d
      coulom(3) = cbar(dble(nchgta),dble(nmasta),1,3) ! t
      coulom(4) = cbar(dble(nchgta),dble(nmasta),2,3) ! h
      coulom(5) = cbar(dble(nchgta),dble(nmasta),2,4) ! a


      if(nmascl.eq.1.and.nchgcl.eq.1)then
      coulomm = coulom(2)
      elseif(nmascl.eq.2.and.nchgcl.eq.1)then
      coulomm = coulom(2)
      elseif(nmascl.eq.3.and.nchgcl.eq.1)then
      coulomm = coulom(3)
      elseif(nmascl.eq.3.and.nchgcl.eq.2)then
      coulomm = coulom(4)
      elseif(nmascl.eq.4.and.nchgcl.eq.2)then
      coulomm = coulom(5)
      endif
!-----------------------------------------------------------------------
      if( eexit .gt. upot*nmascl)then !.and. eexit .lt. barcentr ) then
                
!         ueexit = 17.806

!------------------------------------WS---------------------------------


!	     transm = 1. / ( coul(nchgta, rt00) * 3.0d0 ) * ueexit

!	   	trmax = 1. / ( exp( coul(nchgta, rt00)
!     &		* 2. - coul(nchgta, rt00) * 3.0d0 ) + 1. )				!WS
!	     transm = 1. / ( exp( coul(nchgta, rt00)
!     &		* 2. - ueexit ) + 1. ) / trmax		!WS

!	   	trmax = 1. / ( exp( ( coul(nchgta, rt00)
!     &		* 2.55 - coul(nchgta, rt00) * 3.0d0 ) / 3. )  + 1. )			!WS2
!	     transm = 1. /
!     &		( exp( ( coul(nchgta, rt00)
!     &		* 2.55 - ueexit )
!     &		/ 3. )  + 1. ) / trmax				!WS2

!	      trmax = 1. / ( exp( ( coul(nchgta, rt00)
!     &		- coul(nchgta, rt00) * 3.0d0 ) / 6. )  + 1. )				!WS3
!	     transm = 1. /
!     &		( exp( ( coul(nchgta, rt00)
!     &		- ueexit )
!     &		/ 6. )  + 1. ) / trmax				!WS3

!	  	trmax = 1. / ( exp( ( coul(nchgta, rt00)
!    &		* 2. - coul(nchgta, rt00) * 3.0d0 ) / 6. )  + 1. )			!WS4
!	     transm = 1. / ( exp( ( coul(nchgta, rt00)
!     &		* 2. - ueexit )
!     &		/ 6. )  + 1. ) / trmax				!WS4


!-----------------------------------Gamow-------------------------------

         if( nchgta .ge. 50 )then

            gamop = ( 0.391 - 0.0025 * nchgta ) * 2.0
            vgamow = hbar / ee * nchgcl * dble(nchgta) / ueexit
     &           / rt00 / 2.0 * 3.0d0

         else
            
	     gamop = ( 1.27 - 0.0201 * nchgta ) * 2.0
             vgamow = hbar / ee * nchgcl * dble(nchgta) / ueexit
     &          / rt00

         end if
!	write(*,*)gamop
!	stop

         pgamow = dsqrt( 2.0 * nchgcl *nmascl* rmass * ueexit )

         egamow = dsqrt( 1. / vgamow )

!          transm = dexp( -2.0 * pgamow * rt00 * vgamow / hbar
!     &	     * ( ( dacos( egamow ) )
!     &  	     - egamow * dsqrt( abs( 1.0 - egamow**2.0 ) ) ) )		!gamow
!          transm = dexp( -1.5 * pgamow * rt00 * vgamow / hbar
!     &	     * ( ( dacos( egamow ) )
!     &  	     - egamow * dsqrt( abs( 1.0 - egamow**2.0 ) ) ) )		!gamow2
!          transm = dexp( -1.0 * pgamow * rt00 * vgamow / hbar
!     &	     * ( ( dacos( egamow ) )
!     &  	     - egamow * dsqrt( abs( 1.0 - egamow**2.0 ) ) ) )		!gamow3
!          transm = dexp( -0.5 * pgamow * rt00 * vgamow / hbar
!     &	     * ( ( dacos( egamow ) )
!     &  	     - egamow * dsqrt( abs( 1.0 - egamow**2.0 ) ) ) )		!gamow4
!          transm = dexp( -0.4 * pgamow * rt00 * vgamow / hbar
!     &	     * ( ( dacos( egamow ) )
!     &  	     - egamow * dsqrt( abs( 1.0 - egamow**2.0 ) ) ) )		!gamow5
!          transm = dexp( -0.3 * pgamow * rt00 * vgamow / hbar
!     &	     * ( ( dacos( egamow ) )
!     &  	     - egamow * dsqrt( abs( 1.0 - egamow**2.0 ) ) ) )		!gamow6
!          transm = dexp( -0.2 * pgamow * rt00 * vgamow / hbar
!     &	     * ( ( dacos( egamow ) )
!     &  	     - egamow * dsqrt( abs( 1.0 - egamow**2.0 ) ) ) )		!gamow7
!          transm = dexp( -0.39 * pgamow * rt00 * vgamow / hbar
!     &	     * ( ( dacos( egamow ) )
!     &  	     - egamow * dsqrt( abs( 1.0 - egamow**2.0 ) ) ) )		!gamow8 
!       transm = dexp( -gamop * pgamow * rt00 * vgamow / hbar
!     &     *( ( dacos( egamow ) )
!     &  	     - egamow * dsqrt( abs( 1.0 - egamow**2.0 ) ) ) )

!       if(egamow .gt. 1.0d0)then
!          transm = 1.0
!       endif
!        tansm =1
!-------------------------------------GEM-------------------------------

!		trmax = 1. - 7. / ( coul(nchgta, rt00) * 3.0d0 )
!	     transm = ( 1. - 7. / ueexit) / trmax					!GEM

            
!	write(*,*)vgamow
!	write(*,*)dacos( egamow )
!	write(*,*)egamow * dsqrt( abs( 1.0 - egamow**2.0 ) )
 
!	stop
   
      cmas = dble(nchgcl)/dble(nmascl)
 
!      coula = cmas + (1-cmas) / ( 1 + (exp(coulomm - ueexit )/3.0))
      coulx = 1+exp(-(ueexit/coulomm-1)/0.4)
       
      couly = 1+exp(-(ueexit-1.4*coulomm)
     &/(12+nmasta/200))
     
      coulz = 1/coulx/couly
 
!      dlc1 =1+exp(-(eexit-upot*nmascl-(200/nmasta)**0.7)
!     &/abs(nmasta-70)*25)
!
!      dlc2 =1-0.12*exp(-((eexit-upot*nmascl-6)/15)**2)
!
!      dlc3 =1+exp(-(eexit-upot*nmascl-1.5)) 
!  
!           dlc =dlc2/dlc1/dlc3

       transm =coulz
!       transm = 1
        
!        transm =1.0

      elseif( eexit .le. upot*nmascl ) then
!       endif
         transm = 0.0d0

!      elseif( eexit.ge. barcentr ) then
        
!         transm = 1.0d0

      end if

      END SUBROUTINE gamow

!**********************************************2018/8/28 ishoku watanabe
!*                                                                     *
      real*8  FUNCTION cgbar( z,a,iz,ia,ex )
!*                                                                     *
!*      Calculate Coulomb potential                                    *
!*                                  2014/1/1yamada		       *
!*                                                                     *
!*     input :                                                         *
!*       a   :   mass of nucleus #1                    (IN)            *
!*       z   :   charge  of nucleus #1                 (IN)            *
!*      ia   :   mass of nucleus #2 (emitted)          (IN)            *
!*      iz   :   charge  of nucleus #2 (emitted)       (IN)            *
!*      ex   :   excited energy of parent nuclei(MeV)  (IN)            *
!*                                                                     *
!*    output :                                                         *
!*     cgbar  :   modified Coulomb potential  [MeV]    (OUT)           *
!*                                                                     *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      parameter( rc=1.70d0, ee=137.0359895d0, hbarc=197.327053d0 )

!-----------------------------------------------------------------------
!     transmission probability set
!-----------------------------------------------------------------------

      ck = dostg(1,z-1)     !proton

      if( ia .eq. 2 .and. iz .eq. 1 )ck = ck + 0.06   !deuteron
      if( ia .eq. 3 .and. iz .eq. 1 )ck = ck + 0.12   !triton

      if( ia .eq. 4 .and. iz .eq. 2 )ck = dostg(2,z-2)!alpha
      if( ia .eq. 3 .and. iz .eq. 2 )ck = ck - 0.06   !helium-3

      if( ia .le. 4 .and. iz .le. 2 )then
!-----------------------------------------------------------------------
!     Dostrovsky's parameter set
!-----------------------------------------------------------------------

         cgbar  = hbarc / rc / ee * ck
         r2    = dble(ia)**.333333d0

         if( ia .le. 4 .and. iz .le. 2 ) r2 = 1.2d0 / rc
         if( ia .eq. 1 .and. iz .eq. 1 ) r2 = 0.d0

         r1 = a**.333333d0

         if( a .le. 4 .and. z .le. 2) r1 = 1.2d0/rc

         r0 = r1 + r2

      else
!-----------------------------------------------------------------------
!     Matsuse's parameter set ...PRC26(1982)2338
!-----------------------------------------------------------------------

         cgbar = hbarc / ee

         r1 = 1.12d0 * a**0.333333d0 - 0.86d0 / a**0.333333d0
         r2 = 1.12d0 * dble(ia)**.333333d0 - 0.86d0/dble(ia)**.333333d0

         r0 = r1 + r2 + 3.75d0

      end if
     
      cgbar = cgbar * dble(iz) * z / r0 / ( 1.d0 + 0.005d0 * ex / ia )

      END FUNCTION cgbar
      

!******************************************************2018/8/29watanabe
!*                                                                     *
       SUBROUTINE barrier(p,rt00,nchgta,nmasta,nmaspr,ichg,peex,inds
     &  ,iclst,ipot,iproj,nmasej,upot,q,tfm,tt0,massal)
!      SUBROUTINE barrier(p,rt00,nchgta,nmasta,nmaspr,ichg,peex,inds
!    &  ,iclst,ipot,iproj,nmasej,upot,q,tfm,tt0,massal,itime)
!*							11/27 fukuda   *
!*                                                                     *
!*                                                                     *
!*      [in out]                                                       *
!*       p     : momentum(1-3), energy(4) and mass(5)                  *
!*       nmaspr: mass number of projectile                             *
!*       nmasta: mass number of target                                 *
!*       nchgta: atomic number of target                               *
!*       iclst : identical number of clster                            *
!*       nmasej : mass number of ejectile                              *
!*       ncoll : number of collisions                                  *
!*       ichg  : charge of particle                                    *
!*        q    : Q-value                                               *
!*      [internal variable]                                            *
!*       keje(*):ejectile energy                                       *
!*       t1  :ejectile energy                                          *
!*       reactche : id of reaction chenge or not                       *
!*                                                                     *
!*                                                                     *
!***********************************************************************

      implicit real*8(a-h,o-z)
      INCLUDE 'param-incelf.inc'

      integer  massal,nchgta,i1,nmasta,nmaspr,iproj,i,ioutm,emimasstot
      integer  ipot(nnn),ichg(nnn),iclst(nnn),inds(nnn)
      integer, intent(in)::  nmasej(nnn)
      real*8, intent(out) :: q
      real*8, intent(inout) :: p(6,nnn)
      real*8  rt00,coulom(15),peex,t1(nnn),
     &  keje(nnn),k1,upot
      logical  reactche
!---------------------------------------------------------------------
      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!     if the reaction change, return here 
!---------------------------------------------------------------------

 1001  continue 

!----------------------------------------------------------------------
!     initialize
!----------------------------------------------------------------------     

      emimasstot = 0
      emiketot   = 0.0d0

      do i = 1, massal
         keje(i) = 0.0d0
      end do
!----------------------------------------------------------------------   
!----------------------------------------------------------------------
!     ejectile go out
!----------------------------------------------------------------------
            nej   = 0

            do i = 1, massal    ! check moving particle
             
               ioutm = 0

               if( p(5,i) .lt. 0.1d0 ) cycle
	       if( ipot(i) .eq. 0 )cycle

		 ioutm = nmasej(i)

                 tinn1 = p(4,i) - p(5,i)
                 tout1 = p(4,i) - p(5,i) - dble(ioutm)*upot !tfm                
		 k1    = tout1
                 keje(i) = k1
                 if(  k1 .gt. 0.0d0 )then
                    emimasstot = emimasstot + ioutm
                    nej   = nej + 1
                 end if
!                if( k1 .le. 0.0d0)ipot(i) = 0 
           end do
!-----------------------------------------------------------------------
!     call Q-value
!-----------------------------------------------------------------------

      CALL Q_value(p,nchgta,ipot,nmaspr,nmasta,nmasej
     &        ,ichg,q,iproj,upot)

!----------------------------------------------------------------------

!----------------------------------------------------------------------
!     Coulomb barrier
!----------------------------------------------------------------------

!      coulom(1) = cbar(dble(nchgta),dble(nmasta),1,1) ! p
      coulom(1) = cgbar(dble(nchgta),dble(nmasta),1,1,peex)! p   !2014/1/1yamada
      coulom(2) = cbar(dble(nchgta),dble(nmasta),1,2) ! d
      coulom(3) = cbar(dble(nchgta),dble(nmasta),1,3) ! t
      coulom(4) = cbar(dble(nchgta),dble(nmasta),2,3) ! h
      coulom(5) = cbar(dble(nchgta),dble(nmasta),2,4) ! a

!--------------------------------------------------------------------
       reactche = .false.      !chenge reaction
      do i1 = 1,  massal
!--------------------------------------------------------------------
       if( p(5,i1) .lt. 0.10d0 )cycle     !sawada
       if( ipot(i1) .eq. 0 )cycle
       if( emimasstot .eq. 0)cycle
!--------------------------------------------------------------------
!      K.E. of ejectiles
!--------------------------------------------------------------------
       t1(i1)     =  keje(i1) - q * dble(nmasej(i1))/dble(emimasstot)
     &  +(emimasstot-1)*(upot-tfm)* dble(nmasej(i1))/dble(emimasstot)

	!fukugou ryuusi nyuusya nara.. -1 -> nmaspr

!      write(*,*)(nej -1)*(upot-tfm)* dble(nmasej(i1))/dble(emimasstot)
!     &,nej
!       stop

       if( keje(i1) .lt. 0.0d0 .or. t1(i1) .lt. 0.0d0 )then
          ipot(i1) = 0
          reactche = .true. !reaction chenge
          cycle
       end if

!--------------------------------------------------------------------
!              Coulomb barrier
!--------------------------------------------------------------------
            if( inds(i1) .eq. 1 ) then
                  if( iclst(i1) .eq. 0 )then
                    if( ichg(i1) .eq. 0 ) then
                      if( t1(i1) .gt. 0.0 ) then
                         ipot(i1) = 1
                      else
                         ipot(i1) = 0
                         reactche = .true.
                      end if
                    else
                      if( t1(i1) .gt. coulom(1) ) then
                         ipot(i1) = 1
                      else 
                         ipot(i1) = 0
                         reactche = .true.
                      end if
                    end if
                  else if( iclst(i1) .eq. 1 )then
                    if( t1(i1) .gt. coulom(2) ) then
                        ipot(i1) = 1
                    else
                       ipot(i1) = 0
                       reactche = .true.
                    end if
                  else if( iclst(i1) .eq. 2 )then
                    if( t1(i1) .gt. coulom(3) ) then
                        ipot(i1) = 1
                    else
                       ipot(i1) = 0
                       reactche = .true.
                    end if
                  else if( iclst(i1) .eq. 3 )then
                    if( t1(i1) .gt. coulom(4) ) then
                        ipot(i1) = 1
                    else
                       ipot(i1) = 0
                       reactche = .true.
                    end if
                  else if( iclst(i1) .eq. 4 )then
                    if( t1(i1) .gt. coulom(5) ) then
                        ipot(i1) = 1
                    else
                       ipot(i1) = 0
                       reactche = .true.
                    end if
                  end if
!-----------------------------------------------------------------sawada
            else
                ipot(i1) = 1              ! pion,N* and delta are moving
            end if
         end do       
!-----------------------------------------------------------------------
         if(reactche) goto 1001

!------------------------------------------------------------------
!             Scale Ejectile Energy
!-------------------------------------------------------------------
         etot = 0.0d0
         pp    = 0.0d0

         do i1 = 1 ,  massal
            
            if( p(5,i1) .lt. 0.10d0 )cycle     !sawada
            if( ipot(i1) .eq. 0 )cycle

            if( emimasstot .ne. 0 )then

               p(4,i1) = t1(i1) + p(5,i1)

            else

               t1(i1) = p(4,i1) - p(5,i1)

            end if

                pp     = dsqrt(t1(i1)**2 + 2.0d0*t1(i1)*p(5,i1))
                ptot   = dsqrt(p(1,i1)**2 + p(2,i1)**2
     &               +p(3,i1)**2)
                
                p(1,i1) = pp * p(1,i1)/ptot
                p(2,i1) = pp * p(2,i1)/ptot
                p(3,i1) = pp * p(3,i1)/ptot

         end do

             
         
      END SUBROUTINE barrier

!**********************************************************************
!*
      SUBROUTINE Q_value(p,nchgta,ipot,nmaspr,nmasta,nmasej
     &                   ,ichg,q,iproj,upot)
!*
!*      Purpose: caliculate Q value
!*    
!***********************************************************************
      use NGSDATAMOD, only : energm
      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'
!---------------------------------------------------------------------

      real*8, intent(in) :: p(6,nnn)
      integer, intent(in) :: ipot(nnn),nmaspr,nmasta,ichg(nnn)
     &                       ,nmasej(nnn),iproj
      real*8, intent(out) :: q
!      integer, intent(out) :: iresc,iresm   

      real*8 onke(nnn),qemi
      integer nmastot,iprc,ioutc,ioutm

!      CALL setup ! call gem's data
      
!---------------------------------------------------------------------
      q  = 0.0d0
      qemi = 0.0d0
      iprc = 0
       if( iproj .eq. 1) then! p
               iprc = 1
         else if( iproj .eq. 0 ) then! n
               iprc = 0
	 else
	       iprc = 0
         end if
!---------------------------------------------------------------------
	     nmastot = nmaspr + nmasta
             iresm = nmastot
             iresc = nchgta + iprc

             do i = 1, nmastot  ! check moving particle

                ioutm = 0
                ioutc = 0
                               
                if( p(5,i) .lt. 0.1d0 ) cycle
		if( ipot(i) .eq. 0 ) cycle

                ioutm = nmasej(i)
		ioutc = ichg(i)
                
                onke(i) = p(4,i) - p(5,i)
     &                  - dble(ioutm)*upot
		 
                 
                 if( onke(i) .gt. 0.0d0 )then

                  iresm = iresm - ioutm
                  iresc = iresc - ioutc

                  qemi = energmelf(ioutc,ioutm) + qemi
                  
               end if

              end do            !check moving particle
!--------------------------------------------------------------
            
             qpr  = energmelf(iprc,nmaspr)
             qta  = energmelf(nchgta,nmasta)
             qrs  = energmelf(iresc,iresm)

             q   = qemi + qrs - ( qpr + qta )

      END SUBROUTINE Q_value

      
!**********************************************************************
!								      *
      SUBROUTINE calext(iclst,nmaspr,nmasta,p,tfm,ipot,tt0,eex,upot)
!*                                                                    *
!*                                                                    *
!*      Purpose:                                                      *
!*                                                                    *
!*            to caluculate the excitation energy (E*)                *
!*                                                                    *
!*                                                                    *
!*      Variables:                                                    *
!*                                                                    *
!*         [in]                                                       *
!*            iclst(nnn)  : cluster ID                                *
!*            nmaspr                                                  *
!*            nmasta                                                  *
!*            rmass                                                   *
!*            i           : i-th particle                             *
!*            upot        : potential depth (MeV)                     *
!*            p           : momenta of nucleus                        *
!*	      tfm         : fermi energy			      *
!*            tt0         : sum of inner energy of target nucleons    *
!*                                                                    *
!**********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      real*8, intent(inout) :: p(6,nnn)
      real*8, intent(in) :: tfm,tt0
      real*8, intent(out) :: eex
      integer, intent(in) :: nmaspr,nmasta
      integer, intent(inout) :: iclst(nnn),ipot(nnn)
      integer i,nmastot,nmasrem
      real*8 tbartot,tbar

      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )

!*********************************************************************

	tbartot = 0.0
	nmasrem = 0

	nmastot = nmaspr + nmasta

	do i = 1, nmastot
           
		if( p(5,i) .lt. 0.1d0 ) cycle
		if( ipot(i) .eq. 1 ) cycle

		if( iclst(i) .eq. 0 )then!p,n
			tbar = p(4,i) - rmass
			iremm = 1
		else if( iclst(i) .eq. 1 )then!d
			tbar = p(4,i) - ddmass
			iremm = 2
		else if( iclst(i) .eq. 2 .or. iclst(i) .eq. 3 )then!t,h
			tbar = p(4,i) - tmass
			iremm = 3
		else if( iclst(i) .eq. 4 )then!a
			tbar = p(4,i) - amass
			iremm = 4
		end if

		tbartot = tbartot + tbar !remnant inner energy
		nmasrem = nmasrem + iremm !remnant mass number

	end do

	eex = tbartot - ( tt0 - (nmasta - nmasrem)*tfm )

        
	if( eex .lt. 1.0d-10 )eex = 0.0d0


!========================================================================
!
!	i.e.
!
!	eex = \87\94[i\81\B8Arem,p>pf](Ti-Tf)
!	    + { \87\94[i\81\B8Arem,p<=pf]Ti - [ \87\94[i\81\B8At]Ti0 - (At-Aremf)*Tf ] }
!	
!	Ti : remnant i-th nucleon
!	Ti0 : target i-th nucleon
!	Tf : fermi enegy
!	At : target mass number
!	Arem : remnat mass number
!	Aremf : the number of nucleons under fermi surface
!
!========================================================================

      END SUBROUTINE calext


!*******************************************************************
!*								   *
      SUBROUTINE writenmasej(iclst,nmasej,massal,p,inds,ichg)
!     SUBROUTINE writenmasej(iclst,nmasej,massal,p,inds,ichg,itime,ipot)
!*								   *
!*      purpose: write mass number of ejectile			   *
!*      							   *
!*******************************************************************

      implicit integer(a-z)

      INCLUDE 'param-incelf.inc'

      INTEGER, intent(in) :: iclst(nnn),massal,inds(nnn)
      INTEGER, intent(out) :: nmasej(nnn)
      INTEGER, intent(inout) :: ichg(nnn)
      REAL*8, intent(in) :: p(6,nnn)
      integer ipot(nnn)
      do i = 1, massal

         if(iclst(i).gt.4)cycle
         if(p(5,i).lt.0.1d0)cycle

	 if(inds(i).eq.4)then !pion
	    nmasej(i) = 0 !watanabe 8/29
	    cycle
         elseif(iclst(i).eq.0)then !nucleon
            nmasej(i) = 1
            cycle
         elseif(iclst(i).eq.1)then !deuteron
            nmasej(i) = 2
            cycle
         elseif(iclst(i).eq.2 .or. iclst(i).eq.3)then !triton/3He
            nmasej(i) = 3
            if(iclst(i).eq.3)ichg(i)=2
            cycle
         elseif(iclst(i).eq.4)then !alpha
            nmasej(i) = 4
            ichg(i) = 2
            cycle
         end if
         
      end do
      
      
      END SUBROUTINE writenmasej

!**********************************************************************
!
      SUBROUTINE recoil(iclst,nmaspr,nmasta,ein,upot,p,ipot,
     &     nchgta,iproj,ichg,eex,peex,nmasej,ncoll,icoales,q,
     &     nreact,mstq1,massal,inds)
!      SUBROUTINE recoil(iclst,nmaspr,nmasta,ein,upot,p,ipot,
!     &     nchgta,iproj,ichg,eex,peex,nmasej,ncoll,icoales,q,
!     &     nreact,mstq1,massal,inds,itime)   !for debug
!*
!*      
!*    Purpose: caliculate recoil and output
!*
!*      [in]
!*       mstq1 : method of cal. distribution no keijou wo kimeru.
!*       upot  : potential depth
!*       ein   : K.E. of projectile
!*       peex  : collective excitation 
!*       eex   : excitation of residual nucleus
!*       q     : Q-value
!*       iproj : identical number of projectile
!*
!*      [in out]
!*       p     : momentum(1-3), energy(4) and mass(5)
!*       nmaspr: mass number of projectile
!*       nmasta: mass number of target
!*       nchgta: atomic number of target
!*       iclst : identical number of clster
!*       nmasej : mass number of ejectile
!*       ncoll : number of collisions
!*       ichg  : charge of particle
!*
!*      [internal variable]
!*       massal : number of nmasta + nmaspr + pi
!*       iprc  : charge of projectile
!*       iresc : atomic number of residual nucleus
!*       iresm : mass number of residual nucleus
!*       emike : K.E. of ejectile
!*       px(,y,z)tot : sum of ejectile momentum
!*       resp  : momentum (,energy,mass) of residual nucleus
!*       peje  : sum of ejectile momentum (or mass)
!*       
!***********************************************************************

      implicit real*8 (a-h,o-z)

      INCLUDE 'param-incelf.inc'
      INCLUDE 'param-physcnst.inc'

!---------------------------------------------------------------------
      integer,intent(in)::mstq1(10)
      integer, intent(inout) :: nmaspr,nmasta,nchgta,iclst(nnn)
     &                    ,nmasej(nnn),ncoll(5,nnn),ichg(nnn)
      real*8, intent(in) :: upot,ein,peex,eex,q
      real*8, intent(inout) :: p(6,nnn)
      
      integer i,j,ipot(nnn),massal,icoales(nnn),iputm,emimasstot
      integer ineu(nnn),inds(nnn)
      real*8 onke(nnn),resp(5),peje(5)
      real*8 k1,ke,emike,masres
      
      parameter ( eps    = 1.0d-10 )
      parameter ( dupot  =  90.0d0 )
      parameter ( thupot = 135.0d0 )
      parameter ( aupot  = 180.0d0 )

!--------------------------RECOIL---------------------------------------
!      call setup
!-------------------------select cal. method --------------------------
      idScal = 1 ! 0:ns  1:se 
!      idrec  = 1 ! 0:off(if Ein<eex+peex+q) 1:on     
!------------------------initialize-------------------------------------

      iprc = 0
      emike = 0.0d0
      emiketot = 0.0d0
      pxtot    = 0.0d0
      pytot    = 0.0d0
      pztot    = 0.0d0

!------------------------------watanabe 8/29
      npion = 0 !pion no kazu wo count suru.
!-------------------------------------------

         do i = 1,massal
            onke(i)  = 0.0d0
         end do

         do i = 1,5
            peje(i) = 0.0d0
         end do

!----------------------------------------------------------------------
!                projectile information
!----------------------------------------------------------------------
        
         if( iproj .eq. 0 ) then! neutron
            iprc = 0
            eintot = ein
	    pp00 = dsqrt( eintot * ( eintot + 2.0d0*rmass ) )
         else if( iproj .eq. 1 ) then! proton
            iprc = 1
            eintot = ein
	    pp00 = dsqrt( eintot * ( eintot + 2.0d0*rmass ) )
         else if( iproj .eq. 2 ) then! pion+
            iprc = 1
            eintot = ein
	    pp00 = dsqrt( eintot * ( eintot + 2.0d0*pmass*1000. ) )
         else if( iproj .eq. 3 ) then! pion0
            iprc = 0
            eintot = ein
	    pp00 = dsqrt( eintot * ( eintot + 2.0d0*pmass*1000. ) )
         else if( iproj .eq. 4 ) then! pion-
            iprc = -1
            eintot = ein
	    pp00 = dsqrt( eintot * ( eintot + 2.0d0*pmass*1000. ) )
         end if


         iresm = massal
         iresc = nchgta + iprc

!-------------------------------------------------------------------
!               Combine ejectile particles
!--------------------------------------------------------------------

              do i = 1 , massal

                 if( p(5,i) .le. 0.1d0 ) cycle
                 if( ipot(i) .eq. 0  ) cycle
                 
                 ineu(i) = nmasej(i) - ichg(i)

                 peje(1) = peje(1) + p(1,i)
                 peje(2) = peje(2) + p(2,i)
                 peje(3) = peje(3) + p(3,i)
                 peje(5) = peje(5) + p(5,i)

                 iresm = iresm - nmasej(i)
                 iresc = iresc - ichg(i)
                 nipot = i

		 if( inds(i) .eq. 4 )npion = npion + 1

              end do

              pabsej = sqrt( peje(1)**2 + peje(2)**2 + peje(3)**2 )
              emike = sqrt( pabsej**2 + peje(5)**2 ) - peje(5)

              resp(1) = -peje(1)
              resp(2) = -peje(2)
              resp(3) = pp00 - peje(3)

!              if( idrec .eq. 0)then
!                 if(ein-Q-eex-peex .lt. 0.0)then
!                    write(*,*)ein-Q-eex-peex,Q,eex,peex,ll
!                 goto 1000
!                 endif
!              endif

!-------------------------------------------------------------------
!            virtual 2-body collision
!-------------------------------------------------------------------         
 
              CALL recoil_2b(peje,iresm,massal,nmaspr,nmasta
     &            ,nchgta,p,ipot,eex,peex,masres,emike,pp00,iresc,ein
     &            ,pabsej,resp)

!------------------------------------------------------------------
!             energy conservation
!------------------------------------------------------------------
          if(idScal .ne. 0)then
             if(peje(5) .lt. 0.1d0)goto 999
             
              do i = 1, massal
                 if ( p(5,i) .le. 0.1d0)cycle
                 if (ipot(i) .eq. 0) cycle
                 pabs   = dsqrt( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )
                 em1    = p(5,i)
                 k1     = dsqrt( pabs**2 + em1**2 ) - em1
                 onke(i)= k1
                 p(4,i) = k1 + em1   
                 emiketot = emiketot + k1
              end do

             respp = dsqrt( resp(1)**2 + resp(2)**2 + resp(3)**2 )
             reske = dsqrt( respp**2 + resp(5)**2 ) - resp(5)

             enecm = emiketot + reske

!--------------------------------------------------------------------
!   energy syuusi ni aratani pion no kouryo wo kuwaeta. watanabe 8/29
!--------------------------------------------------------------------
!            totke = eintot - q - eex - peex - pmass*1000.* npion
             totke = eintot - q - eex - peex   !revised by yamaguchi
!--------------------------------------------------------------------

             emiketot = emiketot
             
              do i = 1, massal

                 if ( p(5,i) .le. 0.1d0)cycle
                 if (ipot(i) .eq. 0) cycle
                 if( inds(i) .eq. 4 ) cycle
                 eke    = totke *(onke(i)/enecm)
                 p(4,i) = eke + p(5,i)

                 pp     = sqrt( p(4,i)**2 - p(5,i)**2 )
                 ptot   = dsqrt( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )

                 p(1,i) = pp * p(1,i) / ptot 
                 p(2,i) = pp * p(2,i) / ptot 
                 p(3,i) = pp * p(3,i) / ptot

!                 pxtot  = pxtot + p(1,i)
!                 pytot  = pytot + p(2,i)
!                 pztot  = pztot + p(3,i)

              end do

!              resp(1) = -pxtot
!              resp(2) = -pytot
!              resp(3) = pp00 - pztot
              

 999          continue

           end if

! 1000      continue

      END SUBROUTINE recoil

      
!**********************************************************************
!*
      SUBROUTINE recoil_2b(peje,iresm,massal,nmaspr,nmasta
     &            ,nchgta,p,ipot,eex,peex,masres,emike,pp00,iresc,ein
     &            ,pabsej,resp)
      
!*
!*
!*    purpose:caliculate 2body-collision
!*
!*    [in out]
!*    peje   : ejectile momentum
!*    iresm  : mass number of residual nucleus
!*    p      : momentum(1-3), energy(4) and mass(5)
!*
!*    [in]
!*    emike   : K.E. of ejectile
!*    eex     : excitation of residual nucleus
!*    peex    : collective excitation
!*    massal : mass number of nmasta + nmaspr + pi
!*    nmaspr  : mass number of projectile      
!*    nmasta  : mass number of target
!*    nchg... : charge
!*    pp00    : momentum of projectile
!*    iresm   : mass number of residual nucleus
!*    iresc   : charge of residual nucleus 
!*      
!*    [internal variable]
!*    masres  : mass of residual nucleus [MeV]
!*    i       : counter
!*    Etot    : Total Energy [MeV]    
!*
!*********************************************************************
      use NGSDATAMOD, only : energm
      implicit real*8 (a-h,o-z)
      
      INCLUDE 'param-incelf.inc'
      
!--------------------------------------------------------------------
      real*8,intent(in) :: emike,eex,peex,pp00,pabsej
      real*8,intent(inout)::peje(5),p(6,nnn),resp(5)
      real*8 masres,massej
      integer,intent(inout):: iresm,iresc,nmasta,nchgta,nmaspr
      integer,intent(in)   :: massal
      integer i,ipot(nnn)
!----------------------------------------------------------------------
!      CALL setup                ! call gem's data
      rmasres   = 0.0d0
      rhedefeje = 0.0d0
      eke       = 0.0d0
      ekerel    = 0.0d0
      Etot      = 0.0d0
      A1        = 0.0d0
      A2        = 0.0d0
      B1        = 0.0d0

!----------------------------------------------------------------
      if( peje(5) .lt. 0.1d0)then
         resp(1) = 0.0d0
         resp(2) = 0.0d0
         resp(3) =  pp00
         return
      end if
!----------------------------------------------------------------
!     Total Energy (K.E.(proj) + Mass(proj) + Mass(Target))
!----------------------------------------------------------------

      Etot   = emike + peje(5)
     &     + rmass * dble(iresm) + energmelf(iresc,iresm) + eex +peex

!----------------------------      
!     mass after reaction
!----------------------------
      
      masres = rmass * dble(iresm) + energmelf(iresc,iresm) + eex +peex
      massej = peje(5)
      peje(5) =massej           ! rmass
      
!----------------------------
!     angle of ejectile
!----------------------------
      
      thedefeje = acos(peje(3)/sqrt(peje(1)**2+peje(2)**2+peje(3)**2))
      costhe = cos(thedefeje)
      
!----------------------------
!     Cal. ejectile energy
!----------------------------
      
      A1   = Etot**2.d0 - pabsej**2.d0 + massej**2.d0 - masres**2.d0
      B1   = Etot**2.d0-(pabsej**2.d0)*(costhe**2.d0)
      squ  = A1**2.d0-4.d0*B1*massej**2.d0
      if(squ .lt. 0.0d0)then
         eke = peje(5)
         else
      eke  =(A1*Etot+pabsej*costhe*dsqrt(squ))
     &     /(2.d0*B1) !
      endif
      
      peje(4)  = eke
      eke      = eke - peje(5)
      pp       = dsqrt( eke**2.d0 + 2.d0*eke*peje(5))

!-----normalization ----------------------------------------
      pabs  = dsqrt( peje(1)**2 + peje(2)**2 + peje(3)**2)
      emi_ke = dsqrt( pabs**2 + peje(5)**2 ) - peje(5)
!-----------------------------------------------------------
      
      if( massej .gt. 0.1d0 ) then
      peje(1) = peje(1)*pp/pabs
      peje(2) = peje(2)*pp/pabs 
      peje(3) = peje(3)*pp/pabs
      end if
      
      do i=1,massal
          if( p(5,i) .le. 0.1d0 ) cycle
          if( ipot(i) .eq. 0  ) cycle
         p(1,i) = p(1,i)*pp/pabs
         p(2,i) = p(2,i)*pp/pabs
         p(3,i) = p(3,i)*pp/pabs
      end do
      
      resp(1) = -peje(1)
      resp(2) = -peje(2)
      resp(3) = pp00 - peje(3)
      
      resp(5) = masres

      END SUBROUTINE recoil_2b


!**********************************************************************
!*
!      SUBROUTINE integration(iclst,ipot,inio,ichg,inds,iproj,ncoll,
!     &  nmasta,nchgta,nmaspr,upot,p,r,ucpot,icoales,ininclst,
!     &  ccp,ccr,nmascl,igroup)!watanabe1027
      SUBROUTINE integration(iclst,ipot,inio,ichg,inds,iproj,ncoll,
     &  nmasta,nchgta,nmaspr,upot,p,r,ucpot,icoales,ininclst,
     &  numclst,pcmclst,rcmclst,nmasclst,nmascl,rt00,peex)!wata1027
!*                                                                    *
!*                                                                    *
!*      purpose: to integrate nucleons in cluster                     *
!*                                                                    *
!*                                                                    *
!*                                                                    *
!**********************************************************************

      implicit real(8) (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      real(8), intent(inout) :: p(6,nnn),r(3,nnn)

      real(8), intent(in) :: upot,ucpot,rt00,peex
!      real(8) :: ccp(6,nnn),ccr(3,nnn)
      real(8) :: pcmclst(5,nnn),rcmclst(3,nnn)

      integer, intent(in) :: nmasta,nchgta,nmaspr,iproj

      integer, intent(inout) :: iclst(nnn),icoales(nnn),
     &                       ipot(nnn),inio(nnn),ininclst(nnn),
     &                       ichg(nnn),inds(nnn),ncoll(5,nnn),
     &                       nmascl(nnn)

!      integer :: i,j,numclst(nnn),nmascl(nnn),igroup(nnn)
      integer :: i,j,numclst(nnn),nmasclst(nnn)
      
      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )

         !initialize

        i = 0
        j = 0

!###########################################################

!	i = 0
	nmasal = nmaspr + nmasta
        nmascl = 0

!!----------------------------------------------------------------
!!         MINUS UPOT     NUCLEON
!!-----------------------------------------------------------------
!
!	do i = 1, nmasal
!
!	  if( ipot(i) .eq. 0 )cycle
!
!	  p(4,i) = p(4,i) - upot * dble(nmascl(i))
!
!          pp = sqrt( p(4,i)**2 - p(5,i)**2 )
!
!          pabs = dsqrt( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )
!
!          p(1,i) = pp * p(1,i) / pabs
!          p(2,i) = pp * p(2,i) / pabs
!          p(3,i) = pp * p(3,i) / pabs
!
!	end do
!
!!----------------------------------------------------------------
!!         INTEGRATIPON PART
!!----------------------------------------------------------------
!
!	do i = 1, nmasal
!
!	   if( nmascl(i) .gt. 1 )then
!
!	        do j = 1, 3
!	           p(j,i) = ccp(j,i)
!	           r(j,i) = ccr(j,i)
!	        end do
!
!	        p(4,i) = ccp(4,i)
!	        p(5,i) = ccp(5,i)
!
!	        if( p(4,i) - p(5,i) .gt. 0.0d0 )then
!		    ipot(i) = 1
!	        else
!		    ipot(i) = 0
!	        end if
!
!	        iclst(i) = ininclst(i)
!
!	   else if( igroup(i) .ne. 0 .and. nmascl(i) .le. 1 )then
!
!		do j = 1, 3
!		  p(j,i) = 0.0d0
!		  r(j,i) = -1000.0
!		end do
!
!		p(4:5,i) = 0.0d0
!		ipot(i) = 0
!
!	   end if
!
!	end do

	do i = 1, nmaspr+nmasta
           
		if( numclst(i) .eq. 0 ) cycle

!		do j = i+1, nmaspr+nmasta
		do j = 1, nmaspr+nmasta
                   
!		    if( pcmclst(4,j) .eq. 0.0d0 ) cycle
		    if( pcmclst(4,i) .lt. 1.0d-5 ) cycle
                    
!		    if( numclst(j) .eq. numclst(i) ) then
		    if( numclst(j) .eq. numclst(i) .and.
     &                   j .ne. numclst(j)) then
                       ichg(i) = ichg(i) + ichg(j)
                       ncoll(4,i) = ncoll(4,i) + ncoll(4,j)
		    end if

		end do

		do j = 1, 3
			p(j,i) = pcmclst(j,i)
			r(j,i) = rcmclst(j,i)
		end do

		p(4,i) = pcmclst(4,i)

		p(5,i) = pcmclst(5,i)

		if( p(4,i) - p(5,i) .gt. 
     &                        dble(nmasclst(i))*upot )then

!	         if( nchgta .lt. 50 )then

!   		   if( p(4,i) - p(5,i) .gt. upot*nmasclst(i)+
!     &                  coul(nchgta,rt00)*ichg(i))then
                    
!	                ipot(i) = 1
!		   else
!		        ipot(i) = 0
!	           end if
!!
!	         else

! 	          barcentr = upot*nmasclst(i) + 
!     &                              coul(nchgta,rt00) * 3.0d0 !centrifugal force barrier. baria no takasa.
!   		  if( p(4,i) - p(5,i) .gt. barcentr )then
!-----------------------------------------------------------------------
!                  above the modified Coulomb potential
!-----------------------------------------------------------------------

!	              ipot(i) = 1

!		  else
!-----------------------------------------------------------------------
!                  under the modified Coulomb potential
!                            tunnel effect
!-----------------------------------------------------------------------
!		   coulba = barcentr - upot*nmasclst(i)	!modified Coulomb barrier
		   eexit = p(4,i) - p(5,i)	!intranuclear kinetic energy
		   ueexit = eexit - nmasclst(i)*upot	!extranuclear kinetic energy
!		   coulcut = 7.0
!		   coulcut = 0.0166 * nmasta + 3.7252
!               coulcut = cgbar(dble(nchgta),dble(nmasta),1,2,peex)
!		   if( ueexit .gt. coulcut )then
                      
	            pcoul = rn(0)

                  nmascl(i) = nmasclst(i)
                  nchgcl = ichg(i)
!                  CALL gamow( rt00,nchgpr,nmaspr,nchgta,nmasta,
!     &		      upot,nmasclst(i),nchgcl,barcentr,eexit,ueexit,transm )

		  ipot(i) = 1

	 else

             ipot(i) = 0

	 end if
		

		iclst(i) = ininclst(i)

		if( pcmclst(4,i) .lt. 1.0d-5 )then
			do j = 1, 3
				!p(j,i) = 0.0d0
				!r(j,i) = -1000.0
			end do

			!p(4:5,i) = 0.0d0
			!ipot(i) = 0

		end if

	end do

        
      END SUBROUTINE integration

!***********************************************************************
!
      SUBROUTINE  qvalue(iclst,nmaspr,nmasta,ein,icoales,
     &            i,upot,p,tfm,nchgta,iproj,ichg,q,ipot,idpid)
!*                                                                    *
!*                                                                    *
!*      Purpose:                                                      *
!*                                                                    *
!*            to give energy of Q value                               *
!*                                                                    *
!*                                                                    *
!*      Variables:                                                    *
!*                                                                    *
!*         [in]                                                       *
!*            iclst(nnn)     :                                        *
!*            icoales(nnn)   :                                        *
!*            nmaspr                                                  *
!*            nmasta                                                  *
!*            rmass                                                   *
!*            eall                                                    *
!*            i           : i-th particle                             *
!*            upot        : potential depth (MeV)                     *
!*            p           : momenta   of nucleus                      *
!*            q\81@         : Q value                                   *
!*                                                                    *
!*                                                                    *
!**********************************************************************
      use NGSDATAMOD, only : energm

      implicit real(8) (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      real(8), intent(inout) :: p(6,nnn)
      real(8), intent(in) :: upot,ein
      integer, intent(in) :: nmaspr,nmasta,nchgta,idpid!idst!idst 190730add.
      integer, intent(inout) :: iclst(nnn),ichg(nnn),ipot(nnn)
      integer i,ipotr(nnn),nmastot,icoales(nnn)
      real(8) onke(nnn),q

      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )

!--------------------------RECOIL---------------------------------------
!                       CHECK EMITTED PARTICLE
!    IF ONLY INCIDENT PARTICLE IS EMITTED , ADOPT RECOIL FUNCTION
!-----------------------------------------------8/19-------t.mori-------


        iprc = 0
        iemm = 0
        iemc = 0
        ipotr = 0

        q = 0.0d0
        qemi = 0.0d0

        totke = 0.0d0

        nmastot = nmaspr + nmasta

        onke = 0.0d0

!        CALL setup ! call gem's data

         if( iproj .eq. 1) then! p
               iprc = 1
               eintot = ein
         else if( iproj .eq. 0 ) then! n
               iprc = 0
               eintot = ein
         else if (iproj .eq. 2 ) then! d
               iprc = 1
               eintot = 2.0d0 * ein
         else if (iproj .eq. 3 ) then!\83\BF
               iprc = 2
               eintot = 4.0d0 * ein
         else if (iproj .eq. 5) then
               iprc = 2
               eintot = ein
         end if

         iresm = nmastot
         iresc = nchgta + iprc

             do i = 1, nmastot    ! check moving particle

               ioutm = 0
               ioutc = 0

               if( p(5,i) .lt. 0.1d0 ) cycle
               if( ipot(i) .eq. 0)cycle !191212 fukuda

	       if( iclst(i) .eq. 0 )then
!                  onke(i) = p(4,i) - rmass
                  onke(i) = p(4,i) - rmass - upot
                 ioutm = 1
                if( ichg(i) .eq. 1 )then
                 ioutc = 1
                else
                 ioutc = 0
                end if

	       else if( iclst(i) .eq. 1 )then
!                 onke(i) = p(4,i) - ddmass
               if (idpid .eq. 1)then 
                 onke(i) = p(4,i) - ddmass - dupot
                 ioutm = 2
                 ioutc = 2
               else
                 onke(i) = p(4,i) - ddmass - dupot
                 ioutm = 2
                 ioutc = 1
               end if
	       else if( iclst(i) .eq. 2 )then
!                 onke(i) = p(4,i) - tmass
                 onke(i) = p(4,i) - tmass - thupot
                 ioutm = 3
                 ioutc = 1
	       else if( iclst(i) .eq. 3 )then
!                 onke(i) = p(4,i) - hmass
                 onke(i) = p(4,i) - hmass - thupot
                 ioutm = 3
                 ioutc = 2
	       else if( iclst(i) .eq. 4 )then
!                 onke(i) = p(4,i) - amass
                 onke(i) = p(4,i) - amass - aupot
                 ioutm = 4
                 ioutc = 2
	       end if

               if( onke(i) .gt. 0.0d0 )then

                 iresm = iresm - ioutm
                 iresc = iresc - ioutc 

                 qemi = energmelf(ioutc,ioutm) + qemi

               end if

             end do

             qpr  = energmelf(iprc,nmaspr)
             qta  = energmelf(nchgta,nmasta)
             qrs  = energmelf(iresc,iresm)

             q =  qrs + qemi - qpr - qta


      END SUBROUTINE qvalue


!***********************************************************************
!
      SUBROUTINE  recoil_M(iclst,nmaspr,nmasta,ein,icoales,
     &     i,upot,p,tfm,nchgta,iproj,ichg,q,ipot,idpid)    !t.mori
!*                                                                    *
!*                                                                    *
!*      Purpose:                                                      *
!*                                                                    *
!*            to give energy of recoil                                *
!*                                                                    *
!*                                                                    *
!*      Variables:                                                    *
!*                                                                    *
!*         [in]                                                       *
!*            iclst(nnn)     :                                        *
!*            icoales(nnn)   :                                        *
!*            nmaspr                                                  *
!*            nmasta                                                  *
!*            rmass                                                   *
!*            dome1                                                   *
!*            dome2                                                   *
!*            dome3                                                   *
!*            eall                                                    *
!*            i           : i-th particle                             *
!*            upot        : potential depth (MeV)                     *
!*            p           : momenta   of nucleus                      *
!*                                                                    *
!*                                                                    *
!**********************************************************************

      implicit real(8) (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      real(8), intent(inout) :: p(6,nnn)
      real(8), intent(in) :: upot,ein,q
      integer, intent(in) :: nmaspr,nmasta,nchgta,iclst(nnn),ichg(nnn)
      integer, intent(inout) :: ipot(nnn),idpid!191212 fukuda
      integer i,ipotr(nnn),nmastot,icoales(nnn)
      real(8) onke(nnn),coulom(6)

!      parameter( dupot  =  90.0d0 )
!      parameter( thupot = 135.0d0 )
!      parameter( aupot  = 180.0d0 )
      dupot = upot*2.
      thupot = upot*3.
      aupot = upot *4.

      coulom(1) = cbar(dble(nchgta),dble(nmasta),1,1) ! p
      coulom(2) = cbar(dble(nchgta),dble(nmasta),1,2) ! d
      coulom(3) = cbar(dble(nchgta),dble(nmasta),1,3) ! t
      coulom(4) = cbar(dble(nchgta),dble(nmasta),2,3) ! h
      coulom(5) = cbar(dble(nchgta),dble(nmasta),2,4) ! a

!--------------------------RECOIL---------------------------------------
!                       CHECK EMITTED PARTICLE
!    IF ONLY INCIDENT PARTICLE IS EMITTED , ADOPT RECOIL FUNCTION
!-----------------------------------------------8/19-------t.mori-------

        rmastot = 0.0d0

        pxtot = 0.0d0
        pytot = 0.0d0
        pztot = 0.0d0

        emike = 0.0d0

        nmastot = nmaspr + nmasta

        
         do i = 1,nmastot
              ipotr(i) = 0
              onke(i) = 0.0d0
         end do

!        CALL setup ! call gem's data

         if( iproj .eq. 1) then! p
               eintot = ein
         else if( iproj .eq. 0 ) then! n
               eintot = ein
         else if (iproj .eq. 2 ) then! d
               eintot = 2.0d0 * ein
         else if (iproj .eq. 3 ) then!\83\BF
               eintot = 4.0d0 * ein
         else if( iproj .eq. 5 ) then
               eintot = ein
         end if


             do i = 1, nmastot    ! check moving particle

               if( p(5,i) .lt. 0.1d0 ) cycle
               if( ipot(i) .eq. 0)cycle !191212 fukuda

	       if( iclst(i) .eq. 0 )then
                  onke(i) = p(4,i) - rmass - upot
	       else if( iclst(i) .eq. 1 )then
                 onke(i) = p(4,i) - ddmass - 2.0d0 * upot
	       else if( iclst(i) .eq. 2 )then
                 onke(i) = p(4,i) - tmass - 3.0d0 * upot
	       else if( iclst(i) .eq. 3 )then
                 onke(i) = p(4,i) - hmass - 3.0d0 * upot
                   
	       else if( iclst(i) .eq. 4 )then
                 onke(i) = p(4,i) - amass - 4.0d0 * upot
	       end if

                 if( onke(i) .gt. 0.0d0 )then

                 ipotr(i) = 1

                 tinn1 =  p(4,i) - p(5,i)
                 pinn1 =  dsqrt( tinn1*(tinn1 + 2.0d0 * p(5,i) ) )

                 tout1 =  onke(i)
                 pout1 =  dsqrt( tout1*(tout1 + 2.0d0 * p(5,i) ) )

                 px1 =  p(1,i) * pout1/pinn1
                 py1 =  p(2,i) * pout1/pinn1
                 pz1 =  p(3,i) * pout1/pinn1
                 pionpx=-px1
                 pionpy=0.0
                 pionpz=dsqrt(ein*(ein+2*938.9))-pz1
                 pionpp=dsqrt(pionpx**2+pionpz**2)
                 pionke=dsqrt(pionpp**2+139.6**2)-139.6
                 if(idpid .eq. 1)then
                 if(iclst(i) .eq. 1)then
                 pxtot = pxtot + px1 + pionpx
                 pytot = pytot + py1
                 pztot = pztot + pz1 + pionpz
                 end if
                 else
                 pxtot = pxtot + px1
                 pytot = pytot + py1
                 pztot = pztot + pz1
                 end if
                 !if(dpidp .lt. 0.02)then
                 !if(iclst(i) .eq. 1)then
                 !rmastot = rmastot + (p(5,i)+139.6)/rmass
                 !end if
                 !else
                 rmastot = rmastot + p(5,i)/rmass
                 !end if
                 
                 if(idpid .eq. 1)then
                 if(iclst(i) .eq. 1)then
                 emike = emike + onke(i)+pionke!+139.6
                 end if
                 else
                 emike = emike + onke(i) 
                 end if
                 end if

              end do
              

             emimass = rmass * rmastot
             resmass = rmass * dble( nmastot ) - emimass

	     pp00 = dsqrt( ( eintot + dble(nmaspr) * rmass )**2
     &                    -  ( dble(nmaspr) * rmass )**2 )

             pxtot = - pxtot
             pytot = - pytot
             pztot = pp00 - pztot

             respp = dsqrt( pxtot**2 + pytot**2 + pztot**2 )

             reske = dsqrt( respp**2 + resmass**2 ) - resmass

             enecm = emike + reske !
!             emike = emike !- q -peexs !191109 -q com. out fukuda
             emike = emike - q !-peexs !191109 -q com. out fukuda
             
             if( emike .lt. 0.0 )then
                emike = 0.0d0
                ipot = 0                
             end if
             if( emike .gt. eintot - q)then
             emike = eintot - q !-peexs!200124add. c
             end if

             do i = 1, nmaspr + nmasta

               if( p(5,i) .lt. 0.1d0 ) cycle
               if( ipotr(i) .eq. 0 )cycle

      if(iclst(i).eq.0)then
      coulomm = coulom(2)
      elseif(iclst(i).eq.1)then
      coulomm = coulom(2)
      elseif(iclst(i).eq.2)then
      coulomm = coulom(3)
      elseif(iclst(i).eq.3)then
      coulomm = coulom(4)
      elseif(iclst(i).eq.4)then
      coulomm = coulom(5)
      endif

                 eke = emike * (onke(i)/enecm)
             if(iclst(i) .eq. 1)then
             if(eke .gt. 350)then 
             end if
             end if
             coul1 = 1.+exp(-(eke/coulomm-1.)/0.4)

             coul2 = 1.+exp(-(eke-1.4*coulomm)
     &                          /(8.+nmasta/200.))

             coul3 = (1/coul1/coul2)!-(exp(-eke/10)/10)
             ep=rn(0)
                 if( coul3 .lt. ep)then
                 ipot(i) = 0
                 p(4,i)=p(5,i)
                 else
               
                 p(4,i) = eke + p(5,i)

  	         if( iclst(i) .eq. 0 )then

                  p(4,i) =  p(4,i) + upot

	         else if( iclst(i) .eq. 1 )then

                  p(4,i) =  p(4,i) + 2.0d0*upot
                  
	         else if( iclst(i) .eq. 2 .or. iclst(i) .eq. 3 )then

                  p(4,i) =  p(4,i) + 3.0d0*upot

	         else if( iclst(i) .eq. 4 )then

                  p(4,i) =  p(4,i) + 4.0d0*upot

               end if
               end if
                 pp = sqrt( p(4,i) ** 2 - p(5,i) ** 2 )

                 ptot = dsqrt( p(1,i)**2.0d0 +  p(2,i)**2.0d0 +
     &                p(3,i)**2.0d0 )

                 p(1,i) = pp * p(1,i) / ptot 
                 p(2,i) = pp * p(2,i) / ptot 
                 p(3,i) = pp * p(3,i) / ptot
               
                theaa = datan(dsqrt(p(1,i)**2.+p(2,i)**2.)/p(3,i))
                theaaa = 180.*theaa/pi

              end do 
!============================================================12/18

      END SUBROUTINE recoil_M
      
!**********************************************************************
!*
      SUBROUTINE break_up(iclst,icoales,ipot,inio,ichg,inds,iproj,
     &    nmasta,nchgta,upot,ein,p,r,pfr,idbreak,numclst,ininclst,
     &    relr,pcmclst,rcmclst,nmasclst,nmaspr,idst,idstclst,binp)!t.mori idst add.fuku 190730
!*                                                                    *
!*                                                                    *
!*      incident alpha is devided to nucleons                         *
!*                                                                    *
!*                                                                    *
!*                                                                    *
!**********************************************************************

      implicit real(8) (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      real(8) :: p(6,nnn),r(3,nnn),pfr(3,nnn),pcmclst(5,nnn),
     &           rcmclst(3,nnn),relr(3,nnn),bpara_clst(nnn),inppn
      real(8) :: Efr(4)
      real(8), intent(in) :: upot,ein
!----------------------------sonoda11/13
      integer, intent(in) :: nmasta,nchgta,iproj,nmaspr
      integer, intent(out) :: idst
!----------------------------sonoda11/13
      integer                iclst(nnn),icoales(nnn),
     &                       ipot(nnn),inio(nnn),ininclst(nnn),
     &                       ichg(nnn),inds(nnn),i,j,k,io_clst(nnn)
      integer :: idbreak,nmasclst(nnn),numclst(nnn),ia1,ia2,
     &		skipia(nnn),nmasal

      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )

	  nmasal = nmaspr + nmasta

	  do j = 1,3
		relr(j,1:nmasal) = 0.0d0
          end do

          do i = 1, nmasal
             pcmclst(1:5,i) = 0.0d0
             rcmclst(1:3,i) = 0.0d0
	  end do

          nmasclst(1:nmasal) = 0
	  numclst(1:nmasal) = 0

      if( iclst(1) .eq. 4 .and. iproj .eq. 5 )then

        i = 0
        j = 0
        idbreak = 0
	ininclst(1:nmasal) = 0

!*******************************************************************

        crx = r(1,1) !position of cm of alpha
        cry = r(2,1)
        crz = r(3,1)

!-------------------------------------------------------------

         saa   = 0.54d0         ![fm] ( default: 0.54  )
         ddif2 = 1.00d0         ![fm] ( default: 1.00  )
         dsam2 = 2.25d0         ![fm] ( default: 2.25  )

         rgl = 1.348d0          ![fm] (Glauber model)

         radma = rgl + 5.0 * saa 

         cgl = ( rgl * sqrt(pi) ) ** ( -3.0d0 )

         rrmax = cgl * exp(-(1.0/rgl**2.0)) !Glauber model

         rwod = -1.0

         do while( rn(0) * rrmax .gt. rwod )

           rsqr = 10.0

           do while( rsqr .gt. 1.0 )
             rx = 1.0 - 2.0 * rn(0)
             ry = 1.0 - 2.0 * rn(0)
             rz = 1.0 - 2.0 * rn(0)
             rsqr = rx*rx + ry*ry + rz*rz
           end do

           rrr  = radma * sqrt(rsqr)

           rwod = cgl * exp(-(rrr**2.0d0/rgl**2.0d0)) !Glauber model

         end do

         r(1,1) = radma * rx + crx
         r(2,1) = radma * ry + cry
         r(3,1) = radma * rz + crz

         r(1,2) = -radma * rx + crx
         r(2,2) = +radma * ry + cry
         r(3,2) = -radma * rz + crz

         r(1,3) = -radma * rx + crx
         r(2,3) = -radma * ry + cry
         r(3,3) = radma * rz + crz

         r(1,4) = +radma * rx + crx
         r(2,4) = -radma * ry + cry
         r(3,4) = -radma * rz + crz

!-------------------------------------------------------------

        oein = p(4,1)! - aupot !WIjOÅÌGlM[

        expp = sqrt( oein**2 - amass**2 ) !momentum extra Nucleus

        inpp = sqrt( p(1,1)**2 + p(2,1)**2 + p(3,1)**2 )

        pnax = expp * p(1,1)/inpp
        pnay = expp * p(2,1)/inpp
        pnaz = expp * p(3,1)/inpp

        upot2 = 14.

        do i = 1,4
          p(5,i) = rmass
        end do

        padma = 250. 
       
        RR = 1.348d0
        rho0 = (RR**3.0/pi/hbar) ** 3.0d0        
        
        pxsum = 0
        pysum = 0
        pzsum = 0
        pfr = 0
        Efr = 0

        do i = 1,4

          pgla = -1.0
!
          do while( rn(0) * rho0 .gt. pgla)

            psqr = 10.0

            do while( psqr .gt. 1.0 )

              px = 1.0 - 2.0 * rn(0)
              py = 1.0 - 2.0 * rn(0)
              pz = 1.0 - 2.0 * rn(0)

              psqr = px*px + py*py + pz*pz

            end do

            pfr(1,i) =  padma*px
            pfr(2,i) =  padma*py
            pfr(3,i) =  padma*pz
            
            ppp = padma * sqrt(psqr)
            pgla = rho0 * exp(-1.0 * (RR*ppp/hbar)**2.0d0)
           
          end do
                        
          pxsum = pxsum + pfr(1,i)
          pysum = pysum + pfr(2,i)
          pzsum = pzsum + pfr(3,i)
            
          Efr(i) = sqrt(p(5,i)**2. + ppp**2.)
            
        end do
        
        deltau = 0!upot-upot2
        expp2 = expp+upot*4.
        betaa = expp2/oein

        betaax = -expp2*(p(1,1)/expp)/oein
        betaay = -expp2*(p(2,1)/expp)/oein
        betaaz = -expp2*(p(3,1)/expp)/oein
        gam = 1./sqrt(1.-betaa**2.)

        do i = 1,4            

          pps = pfr(1,i)**2. + pfr(2,i)**2. + pfr(3,i)**2.

          ei = sqrt( pps + p(5,i) **2.)

          pdb = pfr(1,i)*betaax + pfr(2,i)*betaay + pfr(3,i)*betaaz

          p(1,i) = pfr(1,i) + gam*(gam/(gam+1.)*pdb-ei)*betaax
          p(2,i) = pfr(2,i) + gam*(gam/(gam+1.)*pdb-ei)*betaay
          p(3,i) = pfr(3,i) + gam*(gam/(gam+1.)*pdb-ei)*betaaz
	  exppn = sqrt( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )

          p(5,i) = rmass

          p(4,i) = sqrt( p(1,i)**2 + p(2,i)**2 +
     &             p(3,i)**2 + p(5,i)**2 )! + upot

          exppn2 = sqrt(p(4,i)**2. - p(5,i)**2.)

          do jjj = 1,3
            p(jjj,i) = p(jjj,i) * exppn2 / exppn
          end do

        end do         
 
!-----------------------------------------------------------Uèª¯ 190724 fuku.
!        np_clst = 0
!        nn_clst = 0
!        bpara_th = rt00 + 0.54*5.!5.        !!rt00 habetsu no parameter nisuru

!        bpara = 0   

!        do i = 1, 4 !judge bparameter

!          bpara_pre = 
!     &     (betaax * r(1,i) + betaay * r(2,i) + betaaz * r(3,i))/betaa
 
!          bpara_pre2 = r(1,i) **2. + r(2,i) **2. +r(3,i)**2.

!          bpara_clst(i) = dsqrt(bpara_pre2 - bpara_pre **2.)

!          bpara = bpara+bpara_clst(i)

!          if(bpara_clst(i).ge.bpara_th)then
!             io_clst(i) = 1
!             if(i .eq. 1 .or. i.eq. 2)np_clst = np_clst + 1
!             if(i .eq. 3 .or. i.eq. 4)nn_clst = nn_clst + 1
!          else
!             io_clst(i) = 0
!          endif
!       enddo
    
       apara = nmasta**(1./3.)/(27.**(1./3.))
       bpara = sqrt(crx**2. + cry**2.)
       bwdd = 0.08/(1.+exp((binp-1.*apara)/0.1))
       bwtp = bwdd + 0.08/(1.+exp((binp-1.*apara)/0.1))
       bwhn = bwtp + 0.005/(1.+exp((binp-3.*apara)/0.5))
       bwpn = bwhn + 0!0.01/(1.+exp((bpara-3.*apara)/0.5))
       bwal = bwpn+0.001!/(1.+exp((bpara-6.*apara)/0.5))

       frand2 = rn(0)

       if(frand2 .lt. bwdd/bwal) then
          idbreak = 1
       else if(frand2 .lt. bwtp/bwal) then
          idbreak = 2
       else if(frand2 .lt. bwhn/bwal) then
          idbreak = 3
       else if(frand2 .lt. bwpn/bwal) then
          idbreak = 4             
       else
          idbreak = 0
       end if
 
       if( idbreak .eq. 4 )then!2p2n
          ichg(1:2) = 1
          ichg(3:4) = 0
          iclst(1:4) = 0
          numclst(1:4) = 0
          icoales(1:4) = 0
          inio(1:4) = 0

       else if( idbreak .eq. 3 )then!3He and n
          ichg(1:2) = 1
          ichg(3:4) = 0
          iclst(1:4) = 3!
	  ininclst(1:2) = 3
          numclst(1:2) = 1
          ininclst(3) = 3
          ininclst(4) = 0
          iclst(4) = 0
          numclst(3) = 1
          numclst(4) = 0
          idst = 4         !190730 add.
          idstclst = 1
          icoales(1:4) = 0
          inio(1:4) = 0
	  pcmclst(5,1) = hmass
          nmasclst(1) = 3

       else if( idbreak .eq. 2 )then!t and p

          rsqr0 = 1000000.  
          ip = 0

          do i = 1,4
            rsqr = sqrt(r(1,i)**2. + r(2,i)**2.)
            if ( rsqr .lt. rsqr0) then
              ip = i
              rsqr0 = rsqr
            end if
          end do

          rx4 = r(1,4)
          ry4 = r(2,4)
          rz4 = r(3,4)

          ip = 4
 
          do j = 1,3
            r(j,4) = r(j,ip)
          end do

          r(1,ip) = rx4
          r(2,ip) = ry4
          r(3,ip) = rz4

          ichg(1) = 1
          ichg(2:3) = 0
          ichg(4) = 1 
          ininclst(4) = 0
          iclst(1) = 2 !
          iclst(2) = 2 !
          numclst(4) = 0
          numclst(1:3) = 1
          pcmclst(5,1) = tmass
          nmasclst(1) = 3
          idst = 2         !190730 add.
          idstclst = 1
          ininclst(1:3) = 2
          iclst(1:4) = 0
          iclst(1:3) = 2!
          icoales(1:4) = 0
          inio(1:4) = 0

       else if( idbreak .eq. 1 )then!d+d

          ichg(1:2) = 1
          ichg(3:4) = 0
          iclst(1:4) = 1!
	  ininclst(1:4) = 0
          numclst(1) = 1
          numclst(2) = 2
          pcmclst(5,1:2) = ddmass
          nmasclst(1:2) = 2
          ininclst(1:4) = 1
               
          if(io_clst(1) .eq. 1 )then
               
            idstclst = 1
            if(io_clst(3) .eq.1)then !13-24
              numclst(3) = 1
              numclst(4) = 2
              idst = 4
            elseif(io_clst(4) .eq.1)then !14-23
              numclst(3) = 2
              numclst(4) = 1
              idst = 3
            endif

          elseif(io_clst(2) .eq. 1)then
             idstclst = 2
             if(io_clst(3) .eq.1)then !23 -14
               numclst(3) = 2
               numclst(4) = 1
               idst = 4
             elseif(io_clst(4) .eq.1)then !24-13
               numclst(3) = 1
               numclst(4) = 2
               idst = 3
             endif
          endif
            
          numclst(3) = 1
          numclst(4) = 2
          inio(1:4) = 0
          icoales(1:4) = 0
       
       else if(idbreak .eq. 0 ) then !a
          ichg(1:2) = 1
          ichg(3:4) = 0
          iclst(1:4) = 4!
          ininclst(1:4) = 4
	  numclst(1:4) = 1
          icoales(1:4) = 0
          inio(1:4) = 0
	  nmasclst(1) = 4
          pcmclst(5,1) = amass
            ! for test 191030
          pcmclst(4,1) = oein + upot*4.
          expp2 = sqrt(pcmclst(4,1)**2. - pcmclst(5,1)**2.)

          pcmclst(1,1) = expp2*pnax/expp !* 4.
          pcmclst(2,1) = expp2*pnay/expp !* 4.
          pcmclst(3,1) = expp2*pnaz/expp !* 4.
       endif                  !idbreak
       
!--------------------calculate rcm of cluster---------------------
        
	skipia(1:nmaspr) = 0

	do ia1 = 1, nmaspr-1

	  if( numclst(ia1) .eq. 0 ) cycle
	  if( skipia(ia1) .eq. 1 ) cycle

          if(idbreak .eq. 0)exit !191030 youkentou
          
	  do ia2 = ia1, nmaspr

	    if( skipia(ia2) .eq. 1 ) cycle

	    if( numclst(ia2) .eq. numclst(ia1) )then

		skipia(ia2) = 1 !skip i in the following roop

		do j = 1,3
		  pcmclst(j,ia1) = pcmclst(j,ia1) + p(j,ia2)
		  rcmclst(j,ia1) = rcmclst(j,ia1)
     &				+ r(j,ia2)/dble(nmasclst(ia1))
		end do

	    end if

	  end do

!	  pcsum=dsqrt(pcmclst(1,ia1)**2+pcmclst(2,ia1)**2+pcmclst(3,ia1)**2.)
	  pcmclst(4,ia1) = dsqrt( pcmclst(1,ia1)**2 
     &		+ pcmclst(2,ia1)**2 + pcmclst(3,ia1)**2 
     &			+ pcmclst(5,ia1)**2 )
 
          pppd = sqrt(pcmclst(4,ia1)**2. - pcmclst(5,ia1)**2.)
                    
	end do

	skipia(1:nmaspr) = 0

	do ia1 = 1, nmaspr-1

	  if( numclst(ia1) .eq. 0 ) cycle
	  if( skipia(ia1) .eq. 1 ) cycle

	  do ia2 = ia1, nmaspr

	    if( skipia(ia2) .eq. 1 ) cycle

	    if( numclst(ia2) .eq. numclst(ia1) )then

		skipia(ia2) = 1 !skip i in the following roop

		do j = 1,3
		  relr(j,ia2) = r(j,ia2) - rcmclst(j,ia1)
		end do

	    end if

	  end do

	end do

!----------------------------------------------------------------------


      else if( iclst(1) .eq. 1 .and. iproj .eq. 2 )then

        i = 0
        j = 0
        idbreak = 0!not breakup
        nmasclst(:) = 0

        !reaction factor[per]

        bdpn = 0.30d0  ! d to p and n


	frand = rn(0)

        if( frand .lt. bdpn )then

          idbreak = 1

        end if

!###########################################################

!***********************************************************

        if( idbreak .eq. 1 )then

!*******************************************************************

          crx = r(1,1)
          cry = r(2,1)
          crz = r(3,1)

!-------------------------------------------------------------

         saa   = 0.54d0         ![fm] ( default: 0.54  )
         rt00 = 2.15d0

         radma = rt00 + 5.0 * saa

         rrmax = ( exp( - 0.2317 * 1.7 ) - exp( - 1.202 * 1.7 ) )
     &       / 1.7

            rhult = -1.0

            do while( rn(0) * rrmax .gt. rhult 
     &   .or. rrr .le. 0.5 )

               rsqr = 10.0

               do while( rsqr .gt. 1.0 )
                  rx = 1.0 - 2.0 * rn(0)
                  ry = 1.0 - 2.0 * rn(0)
                  rz = 1.0 - 2.0 * rn(0)
                  rsqr = rx*rx + ry*ry + rz*rz
               end do

               rrr  = radma * sqrt(rsqr)
               rhult = ( exp( -0.2317 * rrr ) - exp( -1.202 * rrr ) )
     &              / rrr

            end do

            r(1,1) = radma * rx + crx
            r(2,1) = radma * ry + cry
            r(3,1) = radma * rz + crz

            r(1,2) = - radma * rx + crx
            r(2,2) = - radma * ry + cry
            r(3,2) = - radma * rz + crz


!***********decide a distribution rate of incident energy*********0927

	  oein = p(4,1) - dupot !WIjOÅÌGlM[

	  expp = sqrt( oein**2 - ddmass**2 ) !momentum extra Nucleus

	  inpp = sqrt( p(1,1)**2 + p(2,1)**2 + p(3,1)**2 )

          pndx = expp * p(1,1)/inpp /2.
          pndy = expp * p(2,1)/inpp /2.
          pndz = expp * p(3,1)/inpp /2.

!          pndx = p(1,1)/2.
!          pndy = p(2,1)/2.
!          pndz = p(3,1)/2.

!	  tind = p(4,1) - ddmass - dupot

!        dpotr = - ( 1438.481 * exp( - 3.110 * rrr ) / rrr
!     &      - 570.316 * exp( - 1.550 * rrr ) / rrr )


        tfmd = 15.0 - eliqq(1,1) / dble(2)

        efmd = tfmd + rmass

        pfmd = sqrt( efmd**2. - rmass**2. ) ! [MeV/c]

         do i = 1,2

!           tpn1 = -1.0
!           tpn2 = -1.0

!	   do while( tpn1 .lt. 0.0 .or. tpn1 .gt. tind 
!     &          .or. tpn2 .lt. 0.0 .or. tpn2 .gt. tind )

            psqr = 10.0

            do while( psqr .gt. 1.0 )

               px = 1.0 - 2.0 * rn(0)
               py = 1.0 - 2.0 * rn(0)
               pz = 1.0 - 2.0 * rn(0)

               psqr = px*px + py*py + pz*pz

            end do

            pfr(1,i) =  pfmd*px
            pfr(2,i) =  pfmd*py
            pfr(3,i) =  pfmd*pz

            p(1,i) = pndx + pfr(1,i)
            p(2,i) = pndy + pfr(2,i)
            p(3,i) = pndz + pfr(3,i)
	    exppn = sqrt( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )

            p(5,i) = rmass

            p(4,i) = sqrt( p(1,i)**2 + p(2,i)**2 +
     &               p(3,i)**2 + rmass**2 ) + upot

	    inppn = dsqrt( p(4,i)**2 - rmass**2 )

	    p(1,i) = inppn * p(1,i)/exppn
	    p(2,i) = inppn * p(2,i)/exppn
	    p(3,i) = inppn * p(3,i)/exppn

!            pfr(1,i+1) =  - pfmd*px
!            pfr(2,i+1) =  - pfmd*py
!            pfr(3,i+1) =  - pfmd*pz

!            p(1,i+1) = pndx + pfr(1,i+1)
!            p(2,i+1) = pndy + pfr(2,i+1)
!            p(3,i+1) = pndz + pfr(3,i+1)

!	    tpn1 = p(4,i) - rmass - upot
!	    tpn2 = p(4,i+1) - rmass - upot

!           end do

	 end do

            ichg(1) = 1
            ichg(2) = 0

            iclst(1:2) = 0
            numclst(1:2) = 0

            inio(1:2) = 0

            icoales(1:2) = 0

        end if

      end if

      END SUBROUTINE break_up

!*****************************************************2018/10/19katayama
!*                                                                     *
      SUBROUTINE collex_bending( nmasta,nmaspr,ein,upot,ucpot,p,
     &      iclst,iproj,colexb)
!*                                                                     *
!*                                                                     *
!*        Purpose:                                                     *
!*                                                                     *
!*              to calculate the angle of orientation bent by          *
!*              gamow teller giant resonance                           *
!*                                                                     *
!*                                                                     *
!***********************************************************************

      implicit real(8) (a-h,o-z)

      INCLUDE 'param-incelf.inc'

!-----------------------------------------------------------------------

      integer, intent(in)    :: nmasta, nmaspr,iproj,iclst(nnn)
      integer                :: colexb
      real(8), intent(in)    :: upot,ucpot,ein
      real(8), intent(inout) :: p(6,nnn)
      real(8) :: ke


      parameter( dupot  =  90.0d0 )
      parameter( thupot = 135.0d0 )
      parameter( aupot  = 180.0d0 )
      parameter( cupot  = 540.0d0 )

!-----------------------------------------------------------------------

        if( colexb .eq. 1 )then
          pp0 = 0.0d0
          psi = 0.0d0
          chi = 0.0d0
          phi = 0.0d0
          the = 0.0d0

         ome1 = 0.0d0
         ome2 = 0.0d0
         ome3 = 0.0d0

        dome1 = 0.0d0
        dome2 = 0.0d0
        dome3 = 0.0d0

          opt = 0.0d0
         fthe = 0.0d0                
           ke = 0.0d0
            h = 0.0d0


           ke = p(4,1) - p(5,1) - upot ! t.mori0425

          if( iclst(1) .eq. 1 )then
           ke = p(4,1) - ddmass - dupot
          else if( iclst(1) .eq. 2 .or. iclst(1) .eq. 3)then
           ke = p(4,1) - hmass - thupot
          else if( iclst(1) .eq. 4 )then
           ke = p(4,1) - amass - aupot
!          else if( iclst(1) .eq. 14 )then
!           ke = p(4,1) - cmass12 - cupot
          end if

          pp0 = dsqrt( p(1,1)**2. + p(2,1)**2. + p(3,1)**2. )

         ome1 = p(1,1) / pp0
         ome2 = p(2,1) / pp0
         ome3 = p(3,1) / pp0

          opt = rn(0)

            h = -0.006d0 * dlog(dble(nmasta)) + 0.005d0  !wata
            g = -0.0013 * ke + h
!           g = -0.00497   !2018/10/19katayama

!-----------------------------------------------------------------------
!       Determination of angle[deg] of coll-ex bend
!       based on elastic scattering angular distribution.
!-----------------------------------------------------------------------
              fthe = -10.0

              do while ( opt .ge. fthe )

                 the = 180.0d0 * rn(0)

                 opt = rn(0)

                 fthe = dsin(the / 180.0d0 * pi) *
     &                ( 1.0d0 / (sqrt(2*pi)) ) *
     &                exp( -((the - 5.0d0)**2)/2 )

              end do
!hisatomi-----------------------------------------------------
!             the = 5
!hisatomi------------------------------------------              
              the = the / 180.0d0 * pi

              phi = 2.0d0 * pi * rn(0)

!-----------------------------------------------------------------------
!           set the new momentum coordinates after deflection
!-----------------------------------------------------------------------

         if( dabs(ome3) .ge. 1.0d0 ) then

            dome1 = dsin(the) * dcos(phi)

            dome2 = dsin(the) * dsin(phi)

            dome3 = dcos(the)

         else

            dome1 = ( ome1 * ome3 * dsin(the) * dcos(phi)
     &            - ome2 * dsin(the) * dsin(phi))
     &            / dsqrt( 1.0d0 - ome3**2.0d0) + ome1 * dcos (the)

            dome2 = ( ome2 * ome3 * dsin(the) * dcos(phi)
     &            + ome1 * dsin(the) * dsin(phi) )
     &            / dsqrt( 1.0d0 - ome3**2.0d0) + ome2 * dcos(the)

            dome3 = ome3 * dcos(the) - dsin(the) * dcos(phi)
     &            * dsqrt( 1.0d0 - ome3**2.0d0)

         end if
         p(1,1) = pp0 * dome1

         p(2,1) = pp0 * dome2

         p(3,1) = pp0 * dome3


        end if

      END SUBROUTINE collex_bending


