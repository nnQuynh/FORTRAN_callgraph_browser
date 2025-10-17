c***********************************************************************
c...Include file of JAM
c***********************************************************************

C...Double precision.
      implicit double precision(a-h, o-z)

c...Maximum number of all particles.

      include 'param00.inc'


c...The event record.
      common/coodrp/r(5,nnn),p(6,nnn)
!$OMP THREADPRIVATE(/coodrp/)
      common/jamevnt1/v(5,nnn),k(11,nnn)

      common /jamevnt2/nv,nbary,nmeson
      save /jamevnt1/,/jamevnt2/
!$OMP THREADPRIVATE(/jamevnt1/)
!$OMP THREADPRIVATE(/jamevnt2/)

c...Deuteron coalesence.
      common/jamevnt3/weid(nnn)
      save /jamevnt3/
!$OMP THREADPRIVATE(/jamevnt3/)

c...Parton properites or decay orientation.
      common/jamjet1/vq(10,nnn),kq(2,nnn)
      save /jamjet1/
!$OMP THREADPRIVATE(/jamjet1/)

c...Save two-body collisision information.
      common/jamsave1/rcp(5,6),pcp(6,6),vcp(5,6),kcp(11,6) ! T.Sato 2022/03/25, change pcp(5,6) to pcp(6,6)
      common/jamsave2/vqcp(10,6),kqcp(2,6)
      save /jamsave1/,/jamsave2/
!$OMP THREADPRIVATE(/jamsave1/)
!$OMP THREADPRIVATE(/jamsave2/)

c...Max. collision number in each time step.
      parameter (mxcoll=500000)

c...Collision matrix.
      common/jamcol/coll(6,mxcoll),icoll(3,mxcoll),mentry
      save /jamcol/
!$OMP THREADPRIVATE(/jamcol/)

c...Box by maru  cell.
      common/jamcell/rcell(3,27)
      save /jamcell/
!$OMP THREADPRIVATE(/jamcell/)

c=======================================================================
