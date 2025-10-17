cnais 2023/01/31 samepage
c allocatable_pro
      character*11, allocatable:: chp_n(:)
      character*40, allocatable:: chpchl_n(:,:)
      character*40, allocatable:: chpchm_n(:,:)
      real(8),allocatable:: tott_n(:,:)
      real(8),allocatable:: samewtt(:)
      real(8),allocatable:: totnorm_n(:)
      real(8),allocatable:: tr_n(:,:)
c use at nr*nm > 0
      character*10000, allocatable:: dum1_n(:)
      character*10000, allocatable:: dum2_n(:)
      integer, allocatable :: lng1_n(:)
      integer, allocatable :: lng2_n(:)

      logical :: nmnr_flag
      logical :: samepage_flag
      common /samepage/np_n,nm_n,na_n,nt_n,nr_n,nx_n,ny_n,nz_n,
     &                 ne_n,ne1_n,ne2_n,nact_n,
     &                 maxtott_n,nmnr_flag
      common /samestep/iloopmode,nsame,
     &                npstep,nmstep,nastep,ntstep,nrstep,
     &                nxstep,nystep,nzstep,
     &                nestep,ne1step,ne2step,nactstep

      parameter ( ndatamax = 36 )

      character*3 :: chld(ndatamax)
      character*5 :: chmd(ndatamax)
      character*3 :: chno
      character*200 :: angeltitle, changelsub(10)


