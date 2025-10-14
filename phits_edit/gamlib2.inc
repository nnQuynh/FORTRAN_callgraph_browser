      common/ bigtp  /wtnone,wtfin(len1)
!$OMP THREADPRIVATE(/ bigtp  /)
      common/ data1  /c1,c2,fj(21),cra(250)
!$OMP THREADPRIVATE(/ data1  /)
      common/ dataps /pz(100),pn(150),sz(100),sn(150),con(2),amean(240)
!$OMP THREADPRIVATE(/ dataps /)
      common/ ga     /p2ux,ax,de,ux,exx,tn,c2ssq
!$OMP THREADPRIVATE(/ ga     /)
      common/ new    /nkk,lvls,pmax,addon,ps(500)
!$OMP THREADPRIVATE(/ new    /)
      common/ pss    /psjsum,psj(lms),nstat
!$OMP THREADPRIVATE(/ pss    /)
      common/ repgam /finv,iz,ia,e1,eminxx,isp0,ipar0,icc,irep,ratio
!$OMP THREADPRIVATE(/ repgam /)
      common/ sumry  /sumr(16),tempr(16),ncount(5)
!$OMP THREADPRIVATE(/ sumry  /)
