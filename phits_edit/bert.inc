csib********************************************************************
csib*                                                                  *
csib*       unlabeled common   =>   include file 'bert'                *
csib*                              ('99.02.07)                         *
csib********************************************************************

      common /abcom/ amasno,andit,absec,abx,a,ans,any,awd(3),begru
!$OMP THREADPRIVATE(/abcom/)
      common /ccom/  crsc(2500),ctofe,casesn,c(3),cc(12),ce(21),
     1               cfepn(6),clcfe,clsm,curr(11),col(24),com,com1,
     2               com2,com3,com4,crdt(25),cst
!$OMP THREADPRIVATE(/ccom/)
      common /dcom/  dncms,d(6),dco(3)
!$OMP THREADPRIVATE(/dcom/)
      common /ecom/  einc,e(4),eco(2),efrn,efrp,erand,esps(481),ex
!$OMP THREADPRIVATE(/ecom/)
      common /fcom/  fcn,fcp,ffptfn(3),fmax(7),fmpn(6),frand,ftr,fvnp(3)
!$OMP THREADPRIVATE(/fcom/)
      common /hcom/  hvn(3),hvp(3)
!$OMP THREADPRIVATE(/hcom/)
      common /icom/  i1,i2,i3,i4,i5,i6,ifc,ifca,ifcc,ik,in,inc,inpt,
     1               iout(6),ip,ipec(13),isw(13),
     2               it,itote,itoti,itot2,itot3,iv
!$OMP THREADPRIVATE(/icom/)
      common /klmn/  ln,nrt,ka,knot,med,nmas,no,nor,not,nwds
!$OMP THREADPRIVATE(/klmn/)
      common /opcom/ part(6),pnms,poms,prtin,out(40),p1oe1,p2,pacnt,
     1               pecnt,pgcnt,pgvc(441),plvc(961),pm(4),pmac(3),
     2               pnbc(5),pnidk(23),polc,pols,ppan(3),ppmda,ppnda,
     3               ppnna,pppda,pt(48),pxyz(12)
!$OMP THREADPRIVATE(/opcom/)
      common /rcom/  rcpmv,rlke
!$OMP THREADPRIVATE(/rcom/)
      common /scom/  sf,sqnm,s(37),sigq,sn,snt,sopc,sops,space(178),
     1               strkp
!$OMP THREADPRIVATE(/scom/)
      common /tcom/  trsym,tffn(3),tffp(3),thpn(3)
!$OMP THREADPRIVATE(/tcom/)
      common /uvwxz/ zee,univ,unive,univer,value1,value2,value3,
     1               vnvp(3),xi(3),wkrpn(6)
!$OMP THREADPRIVATE(/uvwxz/)

c-----------------------------------------------------------------------
