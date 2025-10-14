
************************************************************************
*                                                                      *
      subroutine egs5init2
*                                                                      *
*       reset edep (dE)                                                *
*                                                                      *
************************************************************************

      implicit none
!      save
      include 'include/egs5_h.f'
      include 'include/egs5_epcont.f'
      integer iegsxray
      common /egs5cmn8/iegsxray
!$OMP THREADPRIVATE(/egs5cmn8/)
      edep = 0.0d0
      iegsxray = 0

      end

!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
      subroutine copyin_egs5 !<-  copyin_tst
      USE egs5_brempr_mod !<- 2015.08xx allocatable
      USE egs5_eiicom_mod !<- 2015.08xx allocatable
      USE egs5_edge_mod   !<- 2015.08xx allocatable
      implicit none
c
      include 'include/egs5_h.f'
!      save
cc-------------------------
c> cmn3
      real*8 thard, tinel, tmscat, hardstep, sig, scpow, dedx, sig0, ams
      common /egs5cmn3/thard,tinel,tmscat,hardstep,sig,scpow,dedx,sig0,
     $                 ams
!$OMP THREADPRIVATE(/egs5cmn3/)
c> cmn4
      real*8 k1i,k1r,k1s
      common /egs5cmn4/k1i,k1r,k1s
!$OMP THREADPRIVATE(/egs5cmn4/)
c> cmn5
      integer ircode
      common /egs5cmn5/ircode
!$OMP THREADPRIVATE(/egs5cmn5/)
c> cmn6
      real*8 pstep,dpmfp,gmfp
      common /egs5cmn6/pstep,dpmfp,gmfp
!$OMP THREADPRIVATE(/egs5cmn6/)
c> cmn8
       integer iegsxray
       common /egs5cmn8/iegsxray
!$OMP THREADPRIVATE(/egs5cmn8/)
      common /egs5cmn10/denstepold,denstepnew,deinit
      real*8            denstepold,denstepnew,deinit
!$OMP THREADPRIVATE(/egs5cmn10/)

cc-------------------------
c.ada_2
!! include 'include/egs5_brempr.f' !<- allocatable 2015.08xx
c.1    (ibrdst,iprdst, nbrspl,ibrslp,fbrspl)
!! include 'include/egs5_edge.f'
c.2    (exray(),eauger(), nxray,nauger , iz, iextp )
!! include 'include/egs5_eiicom.f'
c.3    ( feispl, ieispl, neispl  )
      include 'include/egs5_userpr.f'
c.4    ( b, blc,omega0,ems,bms,gms, iskpms )
      include 'include/egs5_uservr.f'
c.5    ( cexptr  <- use future private )
      include 'include/randomm.f'
c.6
c
cc-------------------------
c.user
      include 'include/egs5_epcont.f'
      include 'include/egs5_uphiot.f' ! ( theta, sinthe, costhe  )  !<- need copyin ( set sub.hatch )
      include 'include/egs5_useful.f' ! ( RM, medium, medold , iblold ) !<- no need copyin.
c
cc-------------------------
c
! T.Sato 2015/4/17 : add /egs5cmn10/
!>ada.2015.08xx : exchange threadprivate list , from COMMON to variable.
! /BREMPR_priv/ -> nbrspl, ibrspl, iprdst, ibrdst, fbrspl,
! /EIICOM_priv/ -> feispl,ieispl,neispl ,
! /EDGE_priv/   -> exray,eauger,ebind,nauger,nxray,iz,iextp,
!!!!!! exray,eauger,ebind,nauger,nxray,iz,iextp,
c
cc-------------------------
c
!$OMP parallel
!$OMPz copyin(
!$OMPz /egs5cmn3/ ,
!$OMPz /egs5cmn4/ ,
!$OMPz /egs5cmn5/ ,
!$OMPz /egs5cmn6/ ,
!$OMPz /egs5cmn8/ ,
!$OMPz /egs5cmn10/ ,
cccccc /BREMPR_priv/,
!$OMPz fbrspl, ibrdst, iprdst, ibrspl, nbrspl,
cccccc /EDGE_priv/,
!$OMPz exray,eauger,ebind,nauger,nxray,iz,iextp,
cccccc /EIICOM_priv/,
!$OMPz feispl,ieispl,neispl ,
!$OMPz /USERPR/ ,
!$OMPz /USERVR/ ,
!$OMPz /RLUXCOM/,/RLUXDAT/  ,
!$OMPz /EPCONT/ , /EPCONTedep/ ,
!$OMPz /UPHIOT/ ,
!$OMPz /USEFUL/ )
!$OMP end parallel
cc-------------------------
c
      return
      end subroutine copyin_egs5
c     end subroutine copyin_tst
!-----------------------------------------------------------------------
