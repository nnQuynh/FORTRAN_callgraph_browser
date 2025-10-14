************************************************************************
*                                                                      *
      subroutine egs5esteps0
*                                                                      *
*       set electron steps of tinel, tmscat, thart zero
*       last modified by H.Iwase on 2012/3/14
*                                                                      *
************************************************************************
      implicit none
!      save
      include 'include/egs5_h.f'
      include 'include/egs5_epcont.f'

      real*8 thard,tinel,tmscat,hardstep,sig,scpow,dedx,sig0,ams
      real*8 k1i,k1r,k1s
      common /egs5cmn3/thard,tinel,tmscat,hardstep,sig,scpow,dedx,sig0,
     $                 ams
!$OMP THREADPRIVATE(/egs5cmn3/)
      common /egs5cmn4/k1i,k1r,k1s
!$OMP THREADPRIVATE(/egs5cmn4/)

      thard  = 0d0
      tinel  = 0d0
      tmscat = 0d0
      k1i    = 0d0
      k1r    = 0d0
      k1s    = 0d0

      return
      end
