************************************************************************
*                                                                      *
*     'ggsparam.inc'                                                   *
*                                                                      *
*   definition of parameter value for GG                               *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      parameter ( zero = 0.0d0, one = 1.0d0, two = 2.0d0, huge = 1d37 )
      parameter ( half = 0.5d0, three = 3.0d0 )
      parameter ( mlgc = 20000, mxlv = 10 )
      parameter ( pi  = 3.1415926535898d0 )
      parameter ( pie = 3.1415926535898d0 )
      parameter ( ndp2 = 2 )
      parameter ( epss = 1.0d-10 )

      common /gg002/ mlja, lscf, ltrf, lfir, lkst, lksu,
     &               llca, llja, llsc, lncl, lnsf, ljtr,
     &               mlaj, llaj, llcj, lksc, lnlv, lidt
      common /gg003/ tpp(64), ang(3), lgc(mlgc+1)
!$OMP THREADPRIVATE(/gg003/)
      common /gg004/ mxa, mxafs, mxj, mxtr, nlja, nljc
      common /gg005/ xxx, yyy, zzz, uuu, vvv, www, tme, erg,
     &               dls, wgt, vel, dtc,
     &               icl, iii, jjj, kkk, jsu, iap, jgp, ipt,
     &               mtp, iexp, iex, idx,
     &               npa, ncp
!$OMP THREADPRIVATE(/gg005/)
      common /gg006/ dti(mlgc+1), iti(mlgc+1), nlt, jap, kdb, nlaj
!$OMP THREADPRIVATE(/gg006/)
      common /gg007/ krflg
      common /gg008/ junf,nmzu,mlaf,ljun,lmfl,lmzp,llaf,nmaz,
     &               lmzu,lmaz,nlat,llat,lktr,lida
      common /gg009/ lksm
      common /gg010/ kfq
!$OMP THREADPRIVATE(/gg010/)
      common /gg011/ lide, lids
      common /gg012/ lvcl
      common /gg013/ coincd
      common /gg014/ lev, levp
!$OMP THREADPRIVATE(/gg014/)

      common /gg015/ udt(10,0:mxlv), ttn, totm, siga, totgp1,
     &               pfp, stp, pmf, tpd(7), ergace, eg0, ycn,
     &               vtr(3), cmult, ssr, uold(3), deb,
     &               wg0, totnp
!$OMP THREADPRIVATE(/gg015/)

      common /gg016/ llse, klse, nlse !FURUTA



*-----------------------------------------------------------------------
