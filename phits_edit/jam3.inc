c...common block for buu

c...Number of spatial grid.
      parameter(maxx=1,maxz=1)

c...Force acting on hadrons.
      common/jamevf/force(3,nnn)
!$OMP THREADPRIVATE(/jamevf/)

c...Baryon density.
      common/jambuu1/rhob(-maxx:maxx,-maxx:maxx,-maxz:maxz)
!$OMP THREADPRIVATE(/jambuu1/)

c...Charge density.
      common/jambuu2/rhoz(-maxx:maxx,-maxx:maxx,-maxz:maxz)
!$OMP THREADPRIVATE(/jambuu2/)

c....Coulomb potential.
      common/jambuu3/cop(-maxx:maxx,-maxx:maxx,-maxz:maxz)
!$OMP THREADPRIVATE(/jambuu3/)

c...4-momentum of the grid.
      common/jambuu4/avp(4,-maxx:maxx,-maxx:maxx,-maxz:maxz)
!$OMP THREADPRIVATE(/jambuu4/)


c....end of jam3.inc

