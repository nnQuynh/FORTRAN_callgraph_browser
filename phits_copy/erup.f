************************************************************************
*                                                                      *
      subroutine erup(iqstep,lvlopt,
     &                apr,zpr,ex,erec,alp00,bet00,gam00,
     &                atar,ztar,npcle,nhole,efermi,
     &                npart,epart,hepart,cosq2)
*                                                                      *
*       main control routine of dres and qdres                         *
*       evaporation part of PHITS                                      *
*       last modified by K.Niita on 07/02/2000                         *
*                                                                      *
*     input :                                                          *
*                                                                      *
*        iqstep     : =2, evaporation,                                 *
*                     =3, preequ. + evap.                              *
*        lvlopt     : level density options                            *
*                     =1, 8/A                                          *
*                     =2, Baba's parameters                            *
*                     =3, Ignatyuk's parameters                        *
*        apr, zpr   : mass and proton number of mather                 *
*        ex         : excitation energy of mather (MeV)                *
*        erec       : recoil energy of mather (MeV)                    *
*        alp00,bet00,gam00  : unit momentum vector of mather           *
*        atar,ztar  : mass and proton number of grand mather nucleus   *
*                     before cascade,                                  *
*                     which are used in pre-equilibrium                *
*        npcle,nhole : particle and hole number for pre-equ.           *
*        efermi      : fermi energy for pre-equ.                       *
*                                                                      *
*     output:                                                          *
*                                                                      *
*        npart(i)    : number of out going particles                   *
*                i = 1:proton, 2:neutron, 3:deuteron,                  *
*                    4:triton, 5:He3, 6:alpha                          *
*        epart(j,i)  : energy of j-th particles                        *
*                i = 1:proton, 2:neutron                               *
*        hepart(j,i) : energy of j-th particles                        *
*                i = 1:deuteron, 2:triton, 3:He3, 4:alpha              *
*        cosq2(k,j,i): unit momentum vectors of j-th particle          *
*                i = 1:proton, 2:neutron, 3:deuteron,                  *
*                    4:triton, 5:He3, 6:alpha                          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension nparq2(6),eparq2(100,6)

      dimension coslbr(3)
      dimension npart(6),epart(100,2),hepart(100,4)

      dimension npartq(19),epartq(100,2),heparq(100,17)
      dimension cosq(3,100,6)
      dimension cosq2(3,100,6)

*-----------------------------------------------------------------------
*        input values of angle
*-----------------------------------------------------------------------

               coslbr(1) = alp00
               coslbr(2) = bet00
               coslbr(3) = gam00

*-----------------------------------------------------------------------
*        initialization and zero set
*-----------------------------------------------------------------------

 200        do 802 i=1,6
               npart(i)=0
               npartq(i)=0
               nparq2(i)=0
 802        continue

            do 803 i=1,100
            do 804 j=1,2
               epart(i,j)  = 0.d+0
               epartq(i,j) = 0.d+0
               eparq2(i,j) = 0.d+0
 804        continue
 803        continue

            do 805 i=1,3
            do 805 j=1,100
            do 805 k=1,6
               cosq2(i,j,k) =0.d+0
 805        continue

*-----------------------------------------------------------------------
*        check initial values
*-----------------------------------------------------------------------

             uu = ex

         if( uu .lt. -9000.d+0 ) then

            apr = 0.d+0
            zpr = 0.d+0
            return

         end if

            m2 = nint(apr)
            m3 = nint(zpr)

*-----------------------------------------------------------------------
*     Preequilibrium
*-----------------------------------------------------------------------

      if( iqstep .eq. 3 ) then

            m2 = atar
            m3 = ztar

            call qdres(m2,m3,ex,npartq,epartq,hevsuq,
     &                 uu,erecq,heparq,cosq,
     &                 npcle,nhole,efermi,lvlopt)

            zpr = zpr - npartq(2) - npartq(3) - npartq(4)
     &          - 2.0 * ( npartq(5) + npartq(6) )

            apr = apr - npartq(1) - npartq(2)
     &          - 2.0 * npartq(3)
     &          - 3.0 * ( npartq(4) + npartq(5) )
     &          - 4.0 * npartq(6)

            m2 = apr
            m3 = zpr
            ex = uu

*-----------------------------------------------------------------------

         do i = 1, 6
            nparq2(i) = npartq(i)
         end do

         do 8021 j = 1, 100

            do 8022 k = 1, 2
               eparq2(j,k) = epartq(j,k)
 8022       continue

            do 8023 k = 1, 4
               eparq2(j,k+2) = heparq(j,k)
 8023       continue

 8021    continue

      end if

*-----------------------------------------------------------------------
*     call DRES : evaporation code
*-----------------------------------------------------------------------

            iang = 1

            call dres(m2,m3,ex,npart,epart,hevsum,
     &                uu,erec,hepart,iang,
     &                coslbr,cosq2,lvlopt)

*-----------------------------------------------------------------------

            zpr = m3 - npart(2) - npart(3) - npart(4)
     &          - 2.0 * ( npart(5) + npart(6) )

            apr = m2 - npart(1) - npart(2)
     &          - 2.0 * npart(3)
     &          - 3.0 * ( npart(4) + npart(5) )
     &          - 4.0 * npart(6)

            ex  = uu

            zpr = max( 0.0d0, zpr )
            apr = max( 0.0d0, apr )

*-----------------------------------------------------------------------
*     Preequilibrium + evaporation
*-----------------------------------------------------------------------

      if( iqstep .eq. 3 ) then

         do 8888 i = 1, 6

               nkind = i
               nadd  = npartq(i)

            if( nadd .le. 0 ) goto 8888

               nevp = npart(i)

            do 8889 j = 1, nadd

               if( nkind .le. 2 )  epart( nevp + j,i)   = eparq2(j,i)
               if( nkind .gt. 2 ) hepart( nevp + j,i-2) = eparq2(j,i-2)

               cosq2(1, nevp + j,i) = cosq(1,j,i)
               cosq2(2, nevp + j,i) = cosq(2,j,i)
               cosq2(3, nevp + j,i) = cosq(3,j,i)

 8889       continue

               npart(i) = npart(i) + npartq(i)

 8888    continue

      end if

*-----------------------------------------------------------------------

      return
      end
************************************************************************
*                                                                      *
      subroutine dres(m2,m3,t1,npart,epart,sochpe,
     &                u,erec,hepart,iang,
     &                coslbr,cosq2,lvlopt)
*                                                                      *
*       main routine of dres: evaporation part of PHITS                *
*       modified by K.Niita on 07/02/2000                              *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      parameter( y0 = 1.5, b0 = 8.0 )

*-----------------------------------------------------------------------

      common/prab/pp0(1001),pp1(1001),pp2(1001),
     &             rmass(300),alph(300),bet(300)

      common/lvel/aprime(240)

*-----------------------------------------------------------------------

      dimension cosq2(3,100,6)
      dimension ia(6),iz(6),fla(6),flz(6),exmass(6)
      dimension npart(6),epart(100,2),hepart(100,4)
      dimension coslbp(3), coslbr(3)
      dimension c(3)
      dimension zmass(6),
     &    rho(6),r(6),q(6),s(6),sos(6),strun(6),eye0(6),eye1(6),
     &    smom1(6),flkcou(6),ccoul(6),thresh(6),smalla(6),omega(6)

      data um    / 931.504d+0 /
     &     ia    / 1, 1, 2, 3, 3, 4 /
     &     iz    / 0, 1, 1, 1, 2, 2 /
     &     rho   / 2*0.d+0, 4*0.70588d+0 /
     &     omega / 1.d+0, 1.d+0, 3.d+0, 3.d+0, 3.d+0, 2.d+0 /
     &     eps9  / 1.0d-9 /

      data nkey  / 0 /
      save nkey !FURUTA
!$OMP THREADPRIVATE(nkey)
c ----------------------------------------------------------------------
c
cABE 2019/05/28, commentout following two lines
      fkey=0.d+0
      nbe8=0
      nrneep=0
      dmino1=1.d+0
      dmino0=0.d+0
      emnum=energy(dmino1,dmino0)
      emn=emnum+um
      emh=um+energy(dmino1,dmino1)
      emhn=emh-emn
      do 5 k=1,6
      fla(k)=ia(k)
      flz(k)=iz(k)
      exmass(k)=energy(fla(k),flz(k))
      zmass(k)=exmass(k)+fla(k)*um
    5 continue
      flkcou(1)=0.d+0
      ccoul(1)=1.d+0
c
   10 continue
      do 15 i=1,6
      npart(i)=0
   15 smom1(i)=0.d+0
      ja=m2
      jz=m3
      u=t1
      if(ja.eq.2 .and. jz.eq.2) go to 110
      if(ja.ge.jz.and.jz.gt.0)go to 30
      write(6,25)
   25 format(' cascade residual nucleus has mass no. less than z')
      return
   30 a=ja
      z=jz
      eneraz=energy(a,z)
      rnmass=a*um+eneraz
      if (erec) 35,35,40
   35 vrnsq=0.d+0
      vcm=0.d+0
      go to 45
   40 vrnsq=2.0d+0*erec/rnmass
      if (vrnsq.le.0.) then
         vcm=0.d+0
      else
         vcm=sqrt(vrnsq)
         vcm=dsqrt(vrnsq)
      endif
   45 if(ja.eq.8.and.jz.eq.4) go to 210
      flkcou(2)=dost(1,z-flz(2))
      flkcou(3)=flkcou(2)+.06d+0
      flkcou(4)=flkcou(2)+.12d+0
      flkcou(6)=dost(2,z-flz(6))
      flkcou(5)=flkcou(6)-.06d+0
      ccou2=dost(3,z-flz(2))
      ccoul(2)=ccou2+1.d+0
      ccoul(3)=ccou2*1.5d+0+3.0d+0
      ccoul(4)=ccou2+3.0d+0
      ccoul(6)=dost(4,z-flz(6))*2.0d+0+2.0d+0
      ccoul(5)=2.0*ccoul(6)-1.0d+0
      sigma=0.d+0
      paire=0.d+0
      if(a.ge.10.)paire=10.5d+0/sqrt(a-1.d+0)
      corr=0.d+0
      do 50 j=1,6
      q(j)=0.d+0
      r(j)=0.d+0
      s(j)=0.d+0
      sos(j)=0.d+0
      zz=z-flz(j)
      aa=a-fla(j)
      nn=aa-zz
      if(zz.lt.flz(j) .or. nn.lt.(ia(j)-iz(j)))go to 50
      q(j)=energy(aa,zz)+exmass(j)-eneraz
      iaaa=nint(aa)

      if (lvlopt.lt.1 .or. lvlopt.gt.3) then
         write(6,61) lvlopt
   61    format(/' *** error message from s.dres ***'
     &   /' invalid level density parameter option was found.'
     &   /' lvlopt =',i5)
         call parastop( 835 )
      endif

      if (lvlopt.eq.1) then
         aprim = aa/b0
      elseif (lvlopt.eq.2) then
         if (aa.le.240.d+0) then
            aprim=aprime(iaaa)
         else
            aprim=aa/b0
         endif
      elseif (lvlopt.eq.3) then
         aprim = aignt(aa,zz,u)
         if (aprim.le.0.) go to 50
      endif

      smalla(j)=(1.d+0+y0*(1.d+0-2.d+0*zz/aa)**2)*aprim
      mm=ja-ia(j)

      coulef=0.846927d+0*flkcou(j)*flz(j)*zz/(rmass(mm)+rho(j))
      if(j.gt.1) coulef=coulef/(1.0d+0+0.005d+0*u/flz(j))
      thresh(j)=q(j)+coulef
      izz=zz
      dmino2=2.d+0
      if(fkey.ne.1.d+0)
     &   corr=(2.d+0-dmod(zz,dmino2)-dmod(aa-zz,dmino2))*paire
      arg=u-thresh(j)-corr
      if (arg.le.0.) go to 50
      s(j)=dsqrt(smalla(j)*arg)*2.0d+0
      sos(j)=10.0d+0*s(j)
   50 continue
      n1=1
      ses=dmax1(s(1),s(2),s(3),s(4),s(5),s(6))

      if (ses.eq.0.) go to 110
      if (ses.ge.125.) go to 65
      n1=2
      ses=50.d+0
65    do 100 j=1,6
      if(s(j).eq.0.) go to 100
      js=sos(j)+1.d+0
      mm=ja-ia(j)
      if(n1.ne.1.and.js.lt.1000) go to 70
      sas=dexp(s(j)-ses)
      eye1(j)=(s(j)*(s(j)-3.d+0)+3.d+0)*sas/smalla(j)/smalla(j)*.25d+0
      go to 75
70    fjs=js
      strun(j)=fjs-1.d+0
      eye1(j)=(pp1(js)+(pp1(js+1)-pp1(js))*(sos(j)-strun(j)))/
     &        smalla(j)**2
75    if (j-1) 85,85,80
80    r(j)=ccoul(j)*rmass(mm)**2*eye1(j)
      go to 100
85    if (n1.ne.1.and.js.lt.1000) then
         eye0(j)=(pp0(js)+(pp0(js+1)-pp0(js))*(sos(j)-strun(j)))/
     &           smalla(j)
      else
         eye0(j)=(s(j)-1.0d+0)*0.5d+0*sas/smalla(j)
      endif
      r(j)=dmax1(0.0d+0,
     &      (rmass(mm)**2*alph(mm)*(eye1(j)+bet(mm)*eye0(j))))
  100 sigma=sigma+r(j)
  105 ncount=0
      if (sigma.gt.0.) go to 135
  110 continue
      do 115 j=1,6
      if(ja.eq.ia(j).and.jz.eq.iz(j)) go to 120
  115 continue
      if(ja.eq.8.and.jz.eq.4) go to 210
      go to 235
120   jemiss=j
      eps=u+erec
      npart(jemiss)=npart(jemiss)+1
      nrneep=nrneep+1
      smom1(jemiss)=smom1(jemiss)+eps
      k67=npart(jemiss)
      if (iang.gt.0) then
        cosq2(1,k67,jemiss) = coslbr(1)
        cosq2(2,k67,jemiss) = coslbr(2)
        cosq2(3,k67,jemiss) = coslbr(3)
      endif
      if (jemiss-2) 125,125,130
125   epart(k67,jemiss)=eps
      go to 230
130   kemiss=jemiss-2
      hepart(k67,kemiss)=eps
      go to 230
  135 uran=sigma*unirn(dummy)
      sum=0.d+0
      do 140 j=1,6
      jemiss=j
      sum=r(j)+sum
      if(sum.ge.uran) go to 145
  140 continue
      go to 110
145   aa=a-fla(jemiss)
      zz=z-flz(jemiss)
      if (aa.lt.zz) go to 110
      eneraz=energy(aa,zz)
      rnmass=aa*um+eneraz
      js=sos(jemiss)+1.0
      if (js.ge.1000) then
         ratio2=(-15.d+0+s(jemiss)*(15.d+0+s(jemiss)*(s(jemiss)-6.d+0)))
     &         /(s(jemiss)*(s(jemiss)-3.d+0)+3.d+0)/2.d+0/smalla(jemiss)
      else
         ratio2=(pp2(js)+(pp2(js+1)-pp2(js))*
     &          (sos(jemiss)-strun(jemiss)))/smalla(jemiss)
      endif
      epsav=ratio2*2.0d+0
      if (jemiss-1) 160,160,165
160   mm=ja-ia(jemiss)
      epsav=(epsav+bet(mm))/(1.0+bet(mm)*eye0(jemiss)/eye1(jemiss))
165   e1=exprnf(v)/2.0d+0
      e2=exprnf(v)/2.0d+0
      eps=(e1+e2)*epsav+thresh(jemiss)-q(jemiss)
      if (iang.gt.0) go to 901
        coscm = unirn(dummy)
        plormi= unirn(dummy)
        if (plormi - 0.5) 900, 900, 901
  900 coscm = -coscm
  901 ar = a - dble(ia(jemiss))
      zr = z - dble(iz(jemiss))
      be = zr * emhn + ar * emnum - energy(ar, zr)
      rnmass = zr * emh + (ar - zr)*emn -be
      vcmeps=2.0d+0*eps/(zmass(jemiss)*(1.0d+0 + zmass(jemiss)/rnmass))
      if(vcmeps.le.0.) go to 3330
      vcmevp=dsqrt(vcmeps)
      go to 3332
 3330 vcmevp=0.d+0
      write(6,3334)
 3334 format(15x,'vcmevp is negative and change to zero, in s.dres.')
 3332 continue
c
      if (iang.eq.0) go to 311
      c1 = vcmevp
      c2 = zmass(jemiss) *c1/rnmass
      call gtiso(c(1), c(2), c(3))
      coscm = coslbr(1)*c(1) + coslbr(2)*c(2) + coslbr(3)*c(3)
      temp = 2.0d+0 * vcm * coscm
             slp = vrnsq + c1 * (c1 - temp)
             slr = vrnsq + c2 * (c2 + temp)
             slp0= dsqrt(slp)
             slr0= dsqrt(slr)
      do 310 i = 1, 3
         temp = vcm * coslbr(i)
         coslbp(i) = (temp - c1 * c(i))/slp0
         coslbr(i) = (temp + c2 * c(i))/slr0
  310 continue
      vlbeps = slp
      vcm = slr0
      vrnsq = vcm * vcm
      unew  = u -eps - q(jemiss)
      go to 312
  311 continue
      coscm=2.d+0*(unirn(dummy)-0.5d+0)
      vrnsq=vcmeps*zmass(jemiss)*zmass(jemiss)/(rnmass*rnmass)+vcm*vcm+
     1      2.0d+0*zmass(jemiss)/rnmass*vcmevp*vcm*coscm
      vlbeps=vcmeps+vcm*vcm-2.d+0*vcmevp*vcm*coscm
      unew=u-0.5d+0*vcmeps*(zmass(jemiss)*zmass(jemiss)/rnmass+
     1 zmass(jemiss))-q(jemiss)
  312 continue
      eps=0.5d+0*zmass(jemiss)*vlbeps
      if (unew) 175,170,170
170   u=unew
c
      if (iang.eq.0) vcm = dsqrt(vrnsq)
      if (vrnsq.le.0.) then
         vcm=0.d+0
         erec=0.d+0
      else
         vcm=dsqrt(vrnsq)
         erec=0.5d+0*rnmass*vrnsq
      endif
      go to 185
175   ncount=ncount+1
      if (ncount-10) 165,165,180
180   sigma=sigma-r(jemiss)
      r(jemiss)=0.0
      go to 105
185   npart(jemiss)=npart(jemiss)+1
      ja=aa
      jz=zz
      a=aa
      z=zz
      smom1(jemiss)=smom1(jemiss)+eps
      if (npart(jemiss).le.0) go to 240
      k67=npart(jemiss)
      if (iang.gt.0) then
        cosq2(1,k67,jemiss) = coslbp(1)
        cosq2(2,k67,jemiss) = coslbp(2)
        cosq2(3,k67,jemiss) = coslbp(3)
      endif
      if (jemiss-2) 195,195,200
195   epart(k67,jemiss)=eps
      go to 45
200   kemiss=jemiss-2
      hepart(k67,kemiss)=eps
      go to 45
210   eps=0.d+0
      if (u.gt.0.) eps=0.5d+0*(u+0.093d+0)
      eps=0.5d+0*(u+0.093d+0)
      nbe8=nbe8+1
      if (iang.eq.0) go to 1231
          call gtiso(c(1), c(2), c(3))
          coscm = coslbr(1)*c(1) + coslbr(2)*c(2) + coslbr(3)*c(3)
      go to 1232
 1231 continue
      coscm=unirn(dummy)
 1232 continue
      vcmeps=2.0d+0*eps/zmass(6)
      if(vcmeps.le.0.) go to 5550
      vcmevp=dsqrt(vcmeps)
      go to 5552
 5550 vcmevp=0.d+0
      write(6,3334)
 5552 continue
      if (iang.gt.0) c1 = vcmevp /vcm
      vlbeps=vcmeps+vrnsq+2.0d+0*vcmevp*vcm*coscm
      nop=0
  215 eps=0.5d+0*zmass(6)*vlbeps
      smom1(6)=smom1(6)+eps
      npart(6)=npart(6)+1
      k67=npart(6)
      hepart(k67,4)=eps
      if (iang.eq.0) go to 361
        c2 = vcm/dsqrt(vlbeps)
        cosq2(1,k67,6) = c2 * (coslbr(1) + c1 *c(1))
        cosq2(2,k67,6) = c2 * (coslbr(2) + c1 *c(2))
        cosq2(3,k67,6) = c2 * (coslbr(3) + c1 *c(3))
  361 continue
      if (nop) 220,220,225
  220 vlbeps=vcmeps+vrnsq-2.0d+0*vcmevp*vcm*coscm
      if (iang.gt.0) c1 = -c1
      nop=1
      go to 215
  225 continue
  230 erec=0.d+0
      u=0.d+0
  235 sochpe=smom1(3)+smom1(5)+smom1(6)+smom1(4)
      return
c
  240 write(6,245)
  245 format(1h ,'npart(jemiss) negative in dres')
      call exit
      end


************************************************************************
*                                                                      *
      subroutine qdres(m2,m3,t1,npart,epart,sochpe,
     &                 u,erec,hepart,cosq,
     &                 npcle,nhole,efermi,lvlopt)
*                                                                      *
*       main routine of qdres: pre-equilibrium part of PHITS           *
*       modified by K.Niita on 07/02/2000                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h, o-z)

*-----------------------------------------------------------------------

      parameter( y0 = 1.5, b0 = 8.0 )

*-----------------------------------------------------------------------

      common/lvel/aprime(240)

*-----------------------------------------------------------------------

      common/cnst/flkcou(6),ccoul(6),g,pai
!$OMP THREADPRIVATE(/cnst/)
      common/zzzz/zex,zbe,zbeta,qint(19),npz,nhz,npjz
!$OMP THREADPRIVATE(/zzzz/)
      common/qqhvy/flah(19),flzh(19),spin(19),bind(19)

      common/ymo/eincq
!$OMP THREADPRIVATE(/ymo/)

*-----------------------------------------------------------------------

      dimension prbab(20)
      dimension omega2(6),rho2(6),exmass(6)
      dimension cosq(3,100,6)
      dimension q(6),npart(19),epart(100,2),smom1(19),
     1 hepart(100,17),ramdhv(19),zmass(19)

      data zmass/  939.5124,  938.7298, 1876.0177, 2809.2727,
     1 2809.2539, 3728.1883, 3728.1883, 3739.5800, 3739.5800, 3736.26,
     2 3736.2600, 4667.1150, 4667.1150, 4667.4050, 4667.4050, 5600.96,
     3 5600.9600, 6532.9250, 6532.9250/

       data emh,emn,eum,emhn,emnum/ 938.7298, 939.5124, 931.145,
     1 -0.7826, 8.3674/

      data rho2,omega2,exmass/0.,0.,
     1 4*0.70588,1.,1.,3.,3.,3.,2.,
     2 8.3675489, 7.5851116, 13.727994, 15.838178, 15.819549,
     4 3.6092443/

      data flah/1., 1., 2., 3., 3., 4., 4., 4., 4., 4., 4., 5., 5., 5.,
     &          5., 6., 6., 7., 7./

      data flzh/0., 1., 1., 1., 2., 2., 2., 1., 1., 3., 3., 2., 2., 3.,
     &          3., 3., 3., 3., 3./

      data bind/8.367, 7.585, 13.73, 15.84, 15.82, 3.60, 31.13, 16.18,
     &21.48,12.86,16.07,12.87,28.76,13.16,29.91,15.86,29.25,16.98,24.97/

      data spin/2., 2., 3., 2., 2., 1., 28., 5., 6., 5., 7., 4.,
     &         11., 4.,13., 3.,37., 4., 40./

c ----------------------------------------------------------------------

         qa = dble( m2 )
c
cKN test
      npcle = min( npcle, 5 )
      nhole = min( nhole, 5 )
cKN test

c*****  flag set to true if we are penetrating the coulomb barrier
c********** start of initialisations
      ncount=0
      do 100 kk=1,19
      smom1(kk)=0.0
      npart(kk)=0
 100  continue
      qlimit=0.0
      jend1=0
      jend2=0
c********** nucleus ja(mass number),jz(atomic number) excited to u .
 3333 continue
      u=t1
      ja=m2
      jz=m3
      einc=eincq
      jflq=0
      np=npcle
      nh=nhole
      pai=acos(-1.)
      if (np.lt.1) then
c---mri
         return
c-----
      endif
      if (u.lt.0) then
         return
      endif
 2    continue
      rmdmns=0.
      rmdpls=0.
      nex=np+nh
      nexold=nex
c
c***********************************************************
c
c***********************************************************
c
      if (jz.gt.ja) then
         write(6,61) jz,ja
   61    format(/' *** error message from s.qdres ***'
     &   /' atomic and mass number produced nucleus were illegal.'
     &   /' z =',i5,'  a =',i5)
         call parastop( 821 )
      endif
      if (jend2.gt.0) go to 777
c
         a = dble(ja)
         z = dble(jz)

      if (lvlopt.lt.1 .or. lvlopt.gt.3) then
        write(6,* ) ' error !!! --- invarid level density parm option '
         write(6,62) lvlopt
   62    format(/' *** error message from s.qdres ***'
     &   /' invalid level density parameter option was used.'
     &   /' lvlopt =',i5)
         call parastop( 822 )
      endif

       if (lvlopt.eq.1) aprim = a/b0
       if (lvlopt.eq.2)  then
       if (ja .le.240) then
          aprim=aprime(ja)
       else
          aprim=a/b0
       end if
       end if
       if (lvlopt.eq.3)  aprim=aignt(a,z,u )
         g = 6.*aprim/pai/pai
       neq = aint(sqrt(2.*g*u))
       if ((u/nex).le.qlimit)  then
         jend1 = 1
         go to 777
       endif
       if (nex.ge.15)then
         jend1 = 2
         go to 777
       endif
       if (nex.ge.neq)then
         jend1 = 2
         go to 777
       endif

 778   continue

c*preq ----- added next lines
      if(z.ge.a) then
      write(6,6401)
 6401 format('  nucleus has mass no. less than or',
     1       ' equal to z')
      return
      end if
      if(z.lt.1.0)then
      write(6,6602)
 6602 format('  charge number is less than 0 ')
      return
      end if
       bmassd=energy(a,z)
       be=z*emhn + a* emnum -bmassd
       rnmass=z*emh + (a-z)*emn -be
       if(erec.le.0.0) go to 21
       vrnsq=2.0*erec/rnmass
       vcm=sqrt(vrnsq)
       go to 22
  21     continue
        vrnsq=0.0
        vcm=0.0
  22     continue
c****** coulomb barrier constants. ***********************
       flkcou(2) = dost(1,(z-flzh(2)))
       flkcou(3) = flkcou(2)+0.06
       flkcou(4) = flkcou(2)+0.12
       flkcou(6) = dost(2,(z-flzh(6)))
       flkcou(5) = flkcou(6)-0.06
       ccou2     = dost(3,(z-flzh(2)))
       ccoul(2)  = 1.0 + ccou2
       ccoul(3)  = 3.0 + ccou2*1.5
       ccoul(4)  = 3.0 + ccou2
       ccoul(3)=ccoul(3)/3.
       ccoul(4)=ccoul(4)/3.
       ccoul(5)=1.
       ccoul(6)=1.
       do 910 j = 1, 6
          zz1     = z - flzh(j)
          aa1     = a - flah(j)
c*preq ----- added next lines
      if(zz1.ge.aa1) then
      write(6,6402)
 6402 format('  nucleus has mass no. less than or',
     1       ' equal to z')
      return
      end if
      if(zz1.lt.1.0)then
      write(6,6603)
 6603 format('  charge number is less than 0 ')
      return
      end if
          bind(j) = energy(aa1,zz1) -bmassd  + exmass(j)
  910  continue
       do 911 jjj=1,6
       q(jjj)=bind(jjj)
  911  continue
c
c*****************************************************************
c
c       for the 6 clusters.1=n,2=p,3=d,4=t,5=he3,6=he4
c
c*****************************************************************
       if ((np+nh-1).le.0) then
           jend1 = 3
           go to 777
       endif
       call qramda(z,a,u,np,nh,-1,rmdmns,efermi)
       call qramda(z,a,u,np,nh,0,rmdzro,efermi)
       call qramda(z,a,u,np,nh,1,rmdpls,efermi)
       call qneutr(z,a,u,q(1),np,nh,rmdneu)
       call qprton(z,a,u,q(2),np,nh,rmdpro)
c--- type 19
       f=3.5-13*qa**(-0.5)
       f=dmax1(0.2d+0,f)
c----------
       rmdmns=f*rmdmns
       rmdzro=f*rmdzro
       rmdpls=f*rmdpls
c------------
       do 304 i=3,6
       call qheavy(i,z,a,u,np,nh,heavyp)
       ramdhv(i) = heavyp
  304  continue
       if (rmdmns.ge.rmdpls) then
           jend1 = 2
           go to 777
       endif
       subtot=0.
       do 355 i=3,6
 355   subtot=subtot+ramdhv(i)
       totsum=rmdmns+rmdzro+rmdpls+rmdneu+rmdpro+subtot
       prbab(1)=0.
       prbab(2)=rmdmns/totsum
       prbab(3)=rmdzro/totsum+prbab(2)
       prbab(4)=rmdpls/totsum+prbab(3)
       prbab(5)=rmdneu/totsum+prbab(4)
       prbab(6)=rmdpro/totsum+prbab(5)
       do 306 i=3,5
  306  prbab(4+i)=ramdhv(i)/totsum+prbab(3+i)
       prbab(10)=1.
       randm =unirn(dummy)
       do 305 jj=1,9
        if ((prbab(jj).le.randm).and.(randm.lt.prbab(jj+1))) go to 310
  305  continue
  310  continue
       go to (320,330,340,350,360,370,370,370,370,370,370,370,370,370,
     &       370,370) jj
c ********* ramda minus *****
  320  np=np-1
       nh=nh-1
       go to 2
c ********* ramda zero  *****
  330  go to 2
c ********* ramda plus  *****
  340  np=np+1
       nh=nh+1
       go to 2
c*********  neutron emission. *******
  350  k=jj-3
       aaa=a-flah(k)
       zzz=z-flzh(k)
       if (aaa.lt.zzz) go to 777
c----------
       randm = unirn(dummy)
       call emte(k,z,a,np,nh,u,randm,eps,jend2)
       unew = u - eps
       if(unew.lt.bind(k)) then
       go to 2
       endif
       if(eps.lt.0.0) then
         goto 2
       end if
       np=np-1
       if (np.le.0) then
          jend2 = 1
       endif
       nh=nh
       go to 63
c   *****  proton emission.  *******
  360  k=jj-3
       aaa=a-flah(k)
       zzz=z-flzh(k)
       if (aaa.lt.zzz) go to 777
c----------
       randm = unirn(dummy)
       call emte(k,z,a,np,nh,u,randm,eps,jend2)
       unew = u - eps
       if(unew.lt.bind(k)) then
       go to 2
       endif
       if(eps.lt.0.0.or.k.lt.0)    then
         goto 2
       end if
       np=np-1
       if (np.le.0) then
          jend2 = 2
       endif
       go to 63
c********   heavy particle emission   *****************
  370  k=jj-3
       aaa=a-flah(k)
       zzz=z-flzh(k)
       if (aaa.lt.zzz) go to 777
c----------
       randm = unirn(dummy)
      call emte(k,z,a,np,nh,u,randm,eps,jend2)
       unew = u - eps
       if(unew.lt.0.0) then
       go to 2
       endif
       if(eps.lt.0.0.or.k.lt.0)    then
       go to 2
       end if
       if (k.gt.7) then
       end if
       np=np-int(flah(k))
       if (np.le.0) then
          jend2 = 3
       endif
       nh=nh
       go to 63
c********* pick random emission direction and hence new recoil etc.
c
c********* estimated final excitation energy
63     continue
       unew=u-eps
c********* was the process energetically possible ?
       if(unew.lt.0.0) then
       unew=0.
       go to 2
       endif
c********* calculate recoil of nucleus
       rndy = unirn(dummy)
       coscm= (rndy+rndy-1.)
       plormi=(rndy+rndy-1.)
       if(plormi.lt.0.5) coscm=-coscm
       ar=a-flah(k)
       zr=z-flzh(k)
c********* mass of new nucleus.
c*preq ----- added next lines
      if(zr.ge.ar) then
      write(6,6403)
 6403 format('  nucleus has mass no. less than or',
     1       ' equal to z')
      return
      end if
      if(zr.lt.1.0)then
      write(6,6604)
 6604 format('  charge number is less than 0 ')
      return
      end if
       be= zr*emhn + ar* emnum - energy(ar,zr)
       rnmass= zr*emh + (ar-zr)*emn - be
c********* cms velocity**2 (total)
       temp =zmass(k)/rnmass
       vcmeps= 2.0*eps/((1.0+temp)*zmass(k))
       vcmevp=sqrt(vcmeps)
c********* recoil of residual nucleus
       vrnsq2=vcmeps*temp*temp+vrnsq+2.0*vcmevp*vcm*temp*coscm
c********* recoil of cluster
       vlbeps=vcmeps + vrnsq -2.0*vcmevp*vcm*coscm
      eps=0.5*zmass(k)*vlbeps
c
c********************************************************************
c
c*********************************************************************
c
  620  continue
       uold=u
       u=unew
       erec=0.5*vrnsq2*rnmass
       erec=0.0
       ja=ja-int(flah(k))
       jz=jz-int(flzh(k))
       smom1(k)=smom1(k)+eps
       npart(k)=npart(k)+1
cq    if uold/nexold >  7 eps>0 else eps<0
cq       uold: excitation energy before transition
cq       nexold: number of excitons before transition
cq       uold/nexold >  7 (binding energy)---> direct process
cq       uold/nexold <= 7 (binding energy)---> compound pr0ocess
cq     see gudimer's refference
       if (uold/nexold.le. 7.) then
         feps=-1.
       else
         feps= 1.
       endif
        eeps=eps*feps
       if(k.gt.2) go to 51
c
       randm = unirn(dummy)
       epart(npart(k),k)=eps
      sign=2.0
c
 49   continue
cq    kalbach anguler distribution
           call qklbch(k,eeps,exang,ja,jz,einc)
           cosine=cos(exang)
c*preq --- added & changed next lines
      if(k.le.6)then

         theta=unirn(dummy)*6.28318
         sinet=sqrt(1.-cosine**2)

        cosq(1,npart(k),k)=sinet*sin(theta)
        cosq(2,npart(k),k)=sinet*cos(theta)
        cosq(3,npart(k),k)=cosine

      end if
      sign=0.0
       go to 2
51     continue
       randm = unirn(dummy)
           call qklbch(k,eeps,exang,ja,jz,einc)
           cosine=cos(exang)
       hepart(npart(k),(k-2))=eps
      sign=2.0
       go to 2
  105 format(f5.1,i6,f12.4,f12.4,f10.4)
c *** output state at the end****
 777   continue
       call endcde(jend1,jend2)
       npout=np
       nhout=nh
       jend1=0
       jend2=0
       go to 999
  999 return
       end


************************************************************************
*                                                                      *
      subroutine qramda(z,a,ex,npart,nhole,ksign,ramda,efermi)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/cnst/flkcou(6),ccoul(6),g,pai
!$OMP THREADPRIVATE(/cnst/)

c ----------------------------------------------------------------------
      nex=npart+nhole
      an=a-z
      r0=1.2
      tf=efermi
      tcol=tf+ex/nex+3./5.*tf
      xx=tf/tcol
      if (xx.ge.0.5) then
         eta=1.-7.*xx/5.+2./5.*xx*(2.-1./xx)**2.5
      else
         eta=1-7.*xx/5.
      endif
      if (eta.le.0.) then
         write(6,61) eta,xx
   61    format(/' *** error message from s.qramda ***'
     &   /' negative eta was detected.'
     &   /' eta =',1pe13.6,'  xx =',e13.6)
         call parastop( 818 )
      endif
      call qsgm(z,a,tcol,sgm)
      sigma=eta*sgm
      if (sigma.le.0.) then
         write(6,62) sigma,eta,sgm
   62    format(/' *** error message from s.qramda ***'
     &   /' negative sigma was detected.'
     &   /' sigma =',1pe13.6,'  eta =',e13.6,'  sgm =',e13.6)
         call parastop( 819 )
      endif
      vrel=1.383549d7*sqrt(tcol)
c ******** 1.383549e7 <=  v=sqrt(2t/m) (m) *****
      rc=0.6
      rmdhar=4.553638d+0/sqrt(tcol)
c ******** 4.553638 <=  ramda-har = h-har/sqrt(2mt) (fm) ***
      const3=1.d-45
      vint=const3*4.d+0/3.d+0*pai*(2.d+0*rc+rmdhar)**3.
c *******  vint (cubic meter) ************
      al1=alarge(npart,nhole)
      al2=alarge(npart+1,nhole+1)
      if (ksign) 10,20,30
   10 ramdam=sigma*vrel/vint*term3(g,ex,npart,nhole,al1,al2)
      ramda=ramdam
      return
   20 ramda0=sigma*vrel/vint*(nex+1.)/nex
     &      *term2(g,ex,npart,nhole,al1,al2)
      ramda=ramda0
      return
   30 ramdap=sigma*vrel/vint
      ramda=ramdap
      return
      end


************************************************************************
*                                                                      *
      subroutine qneutr(z,a,ex,be,np,nh,ramda)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c
      common/cnst/flkcou(6),ccoul(6),g,pai
!$OMP THREADPRIVATE(/cnst/)

c ----------------------------------------------------------------------
      aa    = a - 1.
      r0    = 1.5
      rz    = aa**(1./3.)
      r     = r0 * rz
      alpha = 0.76 + 1.93 * a**(-1./3.)
      beta  =(1.66 * aa**(-2./3.) - 0.050)/alpha
c ******* r = (cm) ********
      const1=2.268e32
      smlrmd = const1   * 2. * (a-1.) * np * (np+nh-1) * pai * r**2
     &     * alpha / (pai**2) / a / g
c
      vn=0.
      ramda = -smlrmd * qq(np,nh,ex,be,beta,vn)
      if (ramda.lt.0) then
        ramda=0.
      endif
      if (be.ge.ex)then
        iflg1=1
      else
        iflg1=0
      endif
      if(iflg1.eq.1) then
        ramda=0.
      endif
      return
      end


************************************************************************
*                                                                      *
      subroutine qprton(z,a,ex,be,np,nh,ramda)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c
      common/cnst/flkcou(6),ccoul(6),g,pai
!$OMP THREADPRIVATE(/cnst/)
      common/qqhvy/flah(19),flzh(19),spin(19),bind(19)

      dimension omega2(6),rho2(6)
      data rho2,omega2/0.,0.,
     1 4*0.70588,1.,1.,3.,3.,3.,2./
c ----------------------------------------------------------------------
      aa    = a - 1.
      r0    = 1.5
      rz    = aa**(1./3.)
      r     = r0*rz
      alpha = ccoul(2)
      const2=2.268e32
      smlrmd=const2  *2.* (a-1) * np * (np+nh-1) * pai * r**2 * alpha
     &     /(pai**2)/a/g
      flzx  = 1.
      zz   = z - flzx
      vn   = 0.846927 * zz * flzx/(rz+rho2(2))
      no=2
      vn   = vn /(1.0 + 0.005 * ex / flzh(no))
      beta = -flkcou(2)*vn
      tini = abs(beta)
      ramda=-smlrmd*qq(np,nh,ex,be,beta,tini)
      if (ramda.lt.0) then
      ramda=0.
      endif
      if (be.ge.ex)then
        iflg1=1
      else
        iflg1=0
      endif
      if (tini .ge.ex)then
        iflg2=1
      else
        iflg2=0
      endif
      if(iflg1.eq.1 .or. iflg2.eq.1) then
        ramda=0.
      endif
      return
      end


************************************************************************
*                                                                      *
      subroutine qheavy(no,z,a,ex,np,nh,hermd)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c

      common/zzzz/zex,zbe,zbeta,qint(19),npz,nhz,npjz
!$OMP THREADPRIVATE(/zzzz/)
      common/cnst/flkcou(6),ccoul(6),g,pai
!$OMP THREADPRIVATE(/cnst/)
      common/qqhvy/flah(19),flzh(19),spin(19),bind(19)

      dimension ramda(19)
      dimension omega2(6),rho2(6)
      external qhe
      data rho2,omega2/0.,0.,
     1 4*0.70588,1.,1.,3.,3.,3.,2./
c ----------------------------------------------------------------------
c
c  *****  probability of heavy particle emission  ****
c  ****  no = 3- d, 4- t, 5- he3, 6-he4, 7,,,,,,
c
      if (no.le.2) then
         write(6,61) no
   61    format(/' *** error message from s.qheavy ***'
     &   /' illegal particle type was detected, and it must be ',
     &    'greater than 3.'
     &   /' no =',i3,'  (1/2=proton/neutron)')
         call parastop( 820 )
      endif
      zz  = z - flzh(no)
      aa  = a - flah(no)
      npj = flah(no)
      r0  = 1.5e-15
c ****** unit = (m) ***
      rz  = aa**(1./3.)
      r   = r0 * rz
      vol = 4./3.*pai*r**3
      um  = flah(no) * aa/a
      um  = 1.6606e-27*um
c ****** unit = (kg) *******
      if (no.gt.6) then
      vn  = 0.846927 * flzh(no) * zz / (rz + flah(no)**(1./3.))
      else
      vn  = 0.846927 * flzh(no) * zz / (rz + rho2(no))
      endif
      vn = vn / (1.0 + 0.005 * ex / flzh(no))
      if (no.le.6) then
                   alpha = ccoul(no)
                   ba    = flkcou(no)
      else
                   alpha = 1.
                   ba    = 1.
      endif
      beta = -ba * vn
        be = energy(aa,zz) + bind(no) - energy(a,z)
      if (((np+nh-npj) .le. 1) .or. (np-npj).le.0) then
        ramda(no) = 0.
        hermd     = 0.
        return
      endif
      nl  = np + nh - 1
      gj  = npj**3 * (npj/a)**(npj-1)
      bu1 = fact2(np,npj)
      bu2 = fact2(nl,npj)
      bd1 = fact(npj)
      bd2 = fact(npj-1)
      smlrmd = log10(bu1) + log10(bu2) + log10(r**2)
     &    -(log10(vol) + log10((2.*um)**0.5) + log10(bd1) + log10(bd2))
      smlrmd = spin(no) * gj* 10**smlrmd * pai * alpha/2
      x1  = abs(beta)
      x2  = ex - be
            if (x1+be .lt.0.) then
            continue
            else
            continue
            endif
      smlrmd = smlrmd * sqrt(1.602e-13)
c****************** note ********************************************
c*********************************************************************
cccccccccc
      zex   = ex
      zbe   = be
      zbeta = beta
      npz   = np
      nhz   = nh
      npjz  = npj
cccccccccc
      enei  = energy(a,z)
      eneo  = energy(aa,zz)
      if(x1+be.lt.0.)then
          hermd     = 0.
          return
      endif
      zqhe      = qqhev2( npj,np,nh,ex,be,beta,x1,x2)
      ramda(no) = smlrmd * qqhev2(npj,np,nh,ex,be,beta,x1,x2)
      epshe     = abs(zqhe) * 0.01
      call simp (x1,x2,qhe,epshe,hermd,icon)
           epshe  = abs(hermd)*0.01
      call simp (x1,x2,qhe,epshe,hermd,icon)
           qint(no) = hermd
           hermd    = smlrmd * hermd
  10  format(1h ,a,5x,3i10)
  11  format(1h ,a,2e12.5,i6)
      if  (be.ge.ex) then
        iflg1 = 1
      else
        iflg1 = 0
      endif
      if (x1.ge.ex) then
        iflg2 = 1
      else
        iflg2 = 0
      endif
      if(iflg1.eq.1 .or. iflg2.eq.1) then
        hermd = 0.
      endif
      return
      end


************************************************************************
*                                                                      *
      subroutine emte(npsgn,z,a,np,nh,ex,randm,eemt,jend2)
c ----------------------------------------------------------------------
c     if coul.potential>= excitation ene.-binding ene. (proton) eemt=-1.
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c
      common/zzzz/zex,zbe,zbeta,qint(19),npz,nhz,npjz
!$OMP THREADPRIVATE(/zzzz/)
      common/cnst/flkcou(6),ccoul(6),g,pai
!$OMP THREADPRIVATE(/cnst/)
      common/qqhvy/flah(19),flzh(19),spin(19),bind(19)

      dimension rho(6)
      dimension omega2(6),rho2(6)
      external qhe
      data rho /0.,0., 4*0.70588 /
      data rho2,omega2 /0.,0.,4*0.70588,1.,1.,3.,3.,3.,2./
      data eps/1.e-2/
      data zex,zbe,zbeta,npz,nhz,npjz/0.,0.,0.,0,0,0/
c ----------------------------------------------------------------------
c
      eemt=0.
      jend2 = 0
      go to (10,20,30,30,30,30,30,30,30,30,30,30,30,30,
     & 30,30,30,30,30) npsgn
c
c neutron
   10 zz=z
      aa=a-1.
      vn=0.
      r0=1.5e-13
      r=r0*a**(1./3.)
      alpha=0.76+1.93*a**(-1./3.)
      beta=(1.66*a**(-2./3.)-0.050)/alpha
      be=bind(npsgn)
      lpflg=1
      tini=0.
      go to 100
c proton
   20 zz=z-1.
      aa=a-1.
      flzx=1.
      r0=1.5e-13
      rz=aa**(1./3.)
      r=r0*rz
      vn=0.846927*flzx*zz/(rz+rho(2))
      beta=-flkcou(2)*vn
      alpha=ccoul(2)
      be=bind(npsgn)
      tini=abs(beta)
      lpflg=1
      if(tini.ge.ex-be) then
csk --------------------------------- nov. 28 1995 -- blok 1 -- start --
      eemt=-1.
csk --------------------------------- nov. 28 1995 -- blok 1 --  end  --
      return
      endif
      goto 100
c
   30 zz=z-flzh(npsgn)
      aa=a-flah(npsgn)
      rz=aa**(1./3.)
      vn=0.846927*flzh(npsgn)*zz/(flah(npsgn)**(1./3.)+rz)
      be=bind(npsgn)
      mep=int(flah(npsgn))
      npj=mep
      if (npsgn.le.6) then
         alpha=ccoul(npsgn)
         beta=-flkcou(npsgn)*vn
      else
         alpha=1.
         beta=-vn
      endif
cccccc
      zbeta=beta
      zbe=be
      npz=np
      nhz=nh
      npjz=npj
      zex=ex
cccccc
      lpflg=2
      tini = abs(beta)
      if(tini.ge.ex-be) then
         npsgn=-npsgn
         return
      endif
      x1 = tini
      x2 = ex-be
      epshe=qint(npsgn)*1.e-2
      call simp(x1,x2,qhe,epshe,con2,icon)
      con  = con2
      go to 102
c
  100 continue
      zbeta=beta
      zbe=be
      npz=np
      nhz=nh
      npjz=npj
      zex=ex
102   continue
      tnew=0.
      ceff = 0.1
177   t0=tini+((ex-be)-tini)*0.1
      if (t0.le.0..or.t0.gt.(ex-be)) then
         kec = kec +1
c ============< csk >================
c ===================================
         if (kec.gt.20) then
            write(6,61)
   61       format(/' *** error message from s.emte ***'
     &      /' number of invalid initial energy was exceeded 20.')
            call parastop( 834 )
         endif
         eemt = -1.
         return
      endif
      kec = 0
c
      tmax=tini
      tmin=ex-be
      ll=1
  200 continue
      if (ll.gt.10) go to 400
      if (lpflg.eq.1) then
          val = qq(np,nh,ex,be,beta,t0)
          con = qq(np,nh,ex,be,beta,tini)
            f = (val-con)/(-con)-randm
         if (f.ge.0.) tmin = dmin1(t0,tmin)
         if (f.lt.0.) tmax = dmax1(t0,tmax)
      fprime=(np+nh)*(np+nh-1)*(t0+beta)/((ex-be-tini)*(t0+beta)
     &      +(ex-be-tini)**2)*((ex-be-t0)/(ex-be-tini))**(np+nh-2)
      else
cq  goto new method instead of newton method
      call solran(randm,tini,eemt,lpflg,con)
      return
      endif
      const5=1.d-50
      if (fprime.eq.0..or. fprime.lt.const5  ) go to 150
      go to 300
c
  150 rand2 = unirn(dummy)
      ram=rand2
      t0=tini+((ex-be)-tini)*ram
      ll=ll+1
      go to 200
c
 300  continue
      tnew=t0-f/fprime
      if (tnew.lt.0..or.tnew.gt.(ex-be)) go to 150
      deltat=abs(tnew-t0)
      if (deltat.le.eps) then
         eemt=tnew
         return
      else
         t0=tnew
         ll=ll+1
         go to 200
      endif
 400  continue
      call solran(randm,tini,eemt,lpflg,con)
      if (eemt.lt.0..or.eemt.gt.(ex-be)) eemt=-1.
      return
      end


************************************************************************
*                                                                      *
      subroutine qklbch(k,eb,ang,ja,jz,einc)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c
      external qsddcs

      common/aze/aba,zba,ea
!$OMP THREADPRIVATE(/aze/)
      common/sddcs/smla,fmsd,sddcs(0:36),angdeg(0:36)
!$OMP THREADPRIVATE(/sddcs/)

      dimension amass(6),smass(6),fla(6),flz(6)
      data amass/1.0, 1.0, 1.0, 1.0, 0.0, 0.0/
      data smass/0.5, 1.0, 1.0, 1.0, 1.0, 2.0/
      data  fla /1.0, 1.0, 2.0, 3.0, 3.0, 4.0/
      data  flz /0.0, 1.0, 1.0, 1.0, 2.0, 2.0/
      data n2, n3,  c1,     c2,     c3,  et1, et3/
     &      3,  4, .04, 1.8e-6, 6.7e-7, 130., 41./
      data riba, ribb / 0., 0./
      data angdeg / 0.,   5.,  10.,  15.,  20.,  25.,  30.,  35.,
     &             40.,  45.,  50.,  55.,  60.,  65.,  70.,  75.,
     &             80.,  85.,  90.,  95., 100., 105., 110., 115.,
     &            120., 125., 130., 135., 140., 145., 150., 155.,
     &            160., 165., 170., 175., 180./
c ----------------------------------------------------------------------
c
      s(ab,zb,rnb,ac,zc,rnc,rib) =
     & 15.68 * (ac-ab) - 28.07 * ((rnc-zc)**2 /ac-(rnb-zb)**2 /ab)
     &-18.56 * (ac**0.66666 - ab**0.66666)
     &+33.22 * ((rnc-zc)**2 /ac**1.33333 - (rnb-zb)**2 /ab**1.33333)
     &-0.717 * (zc**2 /ac**.33333 - zb**2 /ab**.33333)
     &+1.211 * (zc**2 /ac-zb**2 /ab) - rib
c
      call qaze(ja,jz,einc)
c
      rma = 1.
       ac = aba + 1.
       zc = zba+1.
      rma=amass(k)
      smlb=smass(k)
      abb=ac-fla(k)
      zbb=zc-flz(k)
c
      rnba = aba - zba
      rnbb = abb - zbb
       rnc =  ac - zc
c
      if (eb.gt.0) then
c---type 17
          fmsd = 0.
c-----------------
          iflg=0
      else
c--- type 17
          fmsd = 1.
c-----------------
          iflg=1
            eb = abs(eb)
      endif
c
      pai = acos(-1.)
       sa = s(aba,zba,rnba,ac,zc,rnc,riba)
       sb = s(abb,zbb,rnbb,ac,zc,rnc,ribb)
      eadsh = ea + sa
      ebdsh = eb + sb
         e1 = min(eadsh,et1)
         e3 = min(eadsh,et3)
         x1 = e1 * ebdsh/eadsh
         x3 = e3 * ebdsh/eadsh
c----yoshizawa
       smla = c1 * x1 + c2 * x1**n2 + c3 * rma * smlb * x3**n3
c--- type 17
c
c-------------
c
      rmax =  pai
      rmin = 0.
       eps = 0.0001
             call simp (rmin,rmax,qsddcs,eps,y,ill)
      total = y
c
      do 100 i = 0, 36
            th = angdeg(i)/180. * pai
c
          rmax = th
          rmin = 0.
           eps = 0.0001
                 call simp (rmin,rmax,qsddcs,eps,y,ill)
      sddcs(i) = y
      sddcs(i) = sddcs(i)/total
 100  continue
      call qsolag(ang)
      return
      end


************************************************************************
*                                                                      *
      subroutine endcde(i,j)
c ----------------------------------------------------------------------
      common/code/icdcnt(10)
      data icdcnt/10*0/ !FURUTA
c ----------------------------------------------------------------------
c
      n = 10*j + i
      if (n.ge.10) go to 2
      if (n.eq.1)  icdcnt(1) = icdcnt(1) + 1
      if (n.eq.2)  icdcnt(2) = icdcnt(2) + 1
      if (n.eq.3)  icdcnt(3) = icdcnt(3) + 1
      go to 5
 2        continue
      if (n.lt.20)               icdcnt(4) = icdcnt(4) + 1
      if (n.ge.20 .and. n.lt.30) icdcnt(5) = icdcnt(5) + 1
      if (n.ge.30)               icdcnt(6) = icdcnt(6) + 1
 5    return
      end


************************************************************************
*                                                                      *
      subroutine qsgm(z,a,e,sgm)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c ----------------------------------------------------------------------
c
      c=2.9979e8
      v=1.383549e7*sqrt(e)
      betav=v/c
      sgmpp=(10.63/betav**2-29.92/betav+42.9)*1.e-31
      sgmnp=(34.10/betav**2-82.2/betav+82.2)*1.e-31
      sgm=sgmpp*((a-z)**2+z**2)/a**2+sgmnp*2.*(a-z)*z/a**2
      return
      end


************************************************************************
*                                                                      *
      subroutine simp(a,b,fx,eps,x,is)
c ---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      external fx
c ---------------------------------------------------------------------
c
c*preq --- added next lines
      if(a.ge.b) then
      x=0.
      return
      end if
      n=6
      maxx=100
c
      s2=sss(a,b,fx,n)
  100 continue
      s1=s2
      n=n*3
      s2=sss(a,b,fx,n)
      dif=abs(s1-s2)
      if (n.gt.maxx.or.dif.lt.eps) then
         go to 200
      else
         go to 100
      endif
c
  200 continue
      x=s2
      if (dif.lt.eps) then
         is=0
      else
         is=1
      endif
      return
      end


************************************************************************
*                                                                      *
       subroutine solran(rn,xini,x1,lpflg,con)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      external qhe
c
      common/zzzz/zex,zbe,zbeta,qint(19),npz,nhz,npjz
!$OMP THREADPRIVATE(/zzzz/)
      dimension x(0:100),y(0:100)
c ----------------------------------------------------------------------
c
901   data nx/20/
      xmax = zex - zbe
      x(0) = xini
      y(0) = 0.
        dx = (xmax-xini)/nx
      do 100 i = 0, nx
           x(i)= x(0) + dx * i
           xx  = x(i)
        if (lpflg.eq.1) then
           yy  =  qq(npz,nhz,zex,zbe,zbeta,xini)
           y(i)= (qq(npz,nhz,zex,zbe,zbeta,xx)-yy)/(-yy)
        else
           val = qqhev2(npjz,npz,nhz,zex,zbe,zbeta,xini,xx)
           val = abs(val)
         epshe = val * 1.e-2
                 call simp (xini,xx,qhe,epshe,val2,icon)
         epshe = val2*1.e-2
                 call simp (xini,xx,qhe,epshe,vall,icon)
          y(i) = vall/con
        endif
 100  continue
c
        y1 = rn
      nreg = -1
      do 200 i = 0, nx
             if (y1.ge.y(i)) nreg=i
 200  continue
       ya = y(nreg+1)
       yb = y(nreg)
c
        a = x(nreg+1)
        b = x(nreg)
       y1 = rn
            call m3(x1,ya,yb,a,b,y1)
  11  continue
      end


************************************************************************
*                                                                      *
      subroutine qaze(ja,jz,einc)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/aze/aba,zba,ea
!$OMP THREADPRIVATE(/aze/)

      data ncount/0/ !FURUTA
      save ncount    !FURUTA
!$OMP THREADPRIVATE(ncount)
c ----------------------------------------------------------------------
c
      if (ncount.ne.1010) then
          aba = ja
          zba = jz
           ea = einc
       ncount = 1010
      endif
      return
      end


************************************************************************
*                                                                      *
      subroutine qsolag(ang)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/sddcs/smla,fmsd,sddcs(0:36),angdeg(0:36)
!$OMP THREADPRIVATE(/sddcs/)
      dimension a(1)
c ----------------------------------------------------------------------
c
      a(1)=unirn(dummy)
      pai = acos(-1.)
      ran = a(1)
c
      if(sddcs(0).lt.-.0001) then
      endif
c
      if (ran. le. sddcs(0)) then
          iangrg = 0
      else
        do 10 i = 1, 36
             if (ran. gt. sddcs(i-1)) iangrg=i
 10   continue
      endif
             if (ran. gt. sddcs(36)) then
             if (ran. le. 1.) iangrg=36
      endif
c
      if (iangrg. eq. 0) then
           ang = angdeg(0)
      else
           ang = angdeg(iangrg-1) + ( ran - sddcs(iangrg-1) )
     &         / (  sddcs(iangrg) -  sddcs(iangrg-1) )
     &         * ( angdeg(iangrg) - angdeg(iangrg-1) )
      endif
      ang = ang/180. * pai
      return
      end


************************************************************************
*                                                                      *
      subroutine m3(x1,ya,yb,a,b,y1)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c ----------------------------------------------------------------------
c
      if (ya.eq.yb) then
        x1 = a
      else
        x1 = (y1-ya)*(a-b)/(ya-yb) + a
      endif
      return
      end


************************************************************************
*                                                                      *
      function aignt(a,z,u)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c ----------------------------------------------------------------------
c
      alpha=0.058025d+0
      beta=-5.09059d+0
      gamma=0.40d+0/a**(1./3.)
c
      aa=alpha*a*(1.d+0-beta*a**(-1./3.))
      if(u.eq.0.) then
        aignt=aa
      else
        aignt=aa*(1.d+0+tuyysh(a,z)/u*(1.d+0-exp(-gamma*u)))
      endif
      return
      end


************************************************************************
*                                                                      *
      function alarge(np,nh)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c
      alarge=(np**2+nh**2+np-nh)/4.-nh/2.
      return
      end


************************************************************************
*                                                                      *
      function qqhevy(npj,np,nh,ex,be,beta,x1,x2)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c
      c0=(-1.)**(np+nh-npj-1)
      sum1=0.
      sum2=0.
      nmax=np+nh-2
      do 100 k=0,nmax
      sum1=sum1+uuu(k,x1,beta,be)*ddd(k,ex,beta,be,np,nh,npj)
      sum2=sum2+uuu(k,x2,beta,be)*ddd(k,ex,beta,be,np,nh,npj)
  100 continue
      qqhevy=c0*(sum2-sum1)
      return
      end


************************************************************************
*                                                                      *
      function qqhev2(npj,np,nh,ex,be,beta,x1,x2)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c
      sum = 0.
      mmm = 50
      dx=(x2-x1)/mmm
      do 100 i=1,mmm+1
      x = x1 + (i-1)*dx
      if(x+be.lt.0.)then
      jjj=jjj+1
      end if
      y = ((ex-be-x)/ex)**(np+nh-npj-1)*(x+be)**(npj-1)*(x+be)**(-0.5)
     &   *(x+beta)/ex**npj
      sum = sum + dx*y
  100 continue
      qqhev2=sum
      return
      end


************************************************************************
*                                                                      *
      function ddd(k,ex,beta,bind,np,nh,npj)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c
      sum = 0.
       nl = np + nh - 1
      minim = nl - npj
      maxmn =  k - npj + 1
       nmax = min0(k,minim)
       nmin = max0(0,maxmn)
           do 100 i = nmin, nmax
            c1 = fact(nl-npj) / ( fact(nl-npj-i) * fact(i) )
            c2 = fact(npj-1) /  ( fact(npj-1-k+i)* fact(k-i))
            s1 = c1 * c2 * (bind-ex)**(nl-npj-i) * bind**(npj-1-(k-i))
           sum = sum + s1
  100 continue
      ddd = sum
      return
      end


************************************************************************
*                                                                      *
      function uuu(k,t,beta,bind)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c
      sum = 0
       do 100 l = 0, k
          combi = fact(k)/( fact(k-l) * fact(l) )
             s1 = ( t + bind )**(l+0.5) * combi * (-bind)**(k-l)/(2*l+1)
     &          * ( beta - bind * (k+1)/(k+1-l) )
            sum = sum + s1
 100  continue
            uuu = 2. * ( ( t + bind )**(k+1.5) /( 2 * k + 3 ) + sum )
      return
      end


************************************************************************
*                                                                      *
      function qhe(x)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c
      common/zzzz/zex,zbe,zbeta,qint(19),npz,nhz,npjz
!$OMP THREADPRIVATE(/zzzz/)
c ----------------------------------------------------------------------
c
      qhe = ( (zex - zbe - x)/zex )**(npz + nhz - npjz - 1)
     &    * ( x + zbe )**( npjz - 1 ) * ( x + zbe )**(-0.5)
     &    * ( x + zbeta )/zex**npjz
      return
      end


************************************************************************
*                                                                      *
      function qq(np,nh,ex,be,beta,t)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c
      qq=-1./(np+nh-1)*((ex-be-t)/ex)**(np+nh-1)*(t+beta)
     &   -1./(np+nh-1)/(np+nh)*((ex-be-t)/ex)**(np+nh-1)*(ex-be-t)
      qq=1.602e-13*qq
      return
      end


************************************************************************
*                                                                      *
      function qsddcs(x)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c
      common/sddcs/smla,fmsd,sddcs(0:36),angdeg(0:36)
!$OMP THREADPRIVATE(/sddcs/)
c ----------------------------------------------------------------------
      pai = acos(-1.)
c
      qsddcs=2*pai*sin(x)*1.0
c-----------
      return
      end


************************************************************************
*                                                                      *
      function fact(ix)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
c ----------------------------------------------------------------------
c
      if (ix.lt.0) then
         write(6,61) ix
   61    format(/' *** error message from s.fact ***'
     &   /' illegal index was detected.'
     &   /' ix =',i5)
         call parastop( 829 )
      endif
      fact=1.
      if (ix.eq.0) go to 5
      do 100 i=1,ix
 100   fact=fact*i
  5   return
      end


************************************************************************
*                                                                      *
      function fact2(ix,ia)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
c ----------------------------------------------------------------------
c
      if (ix.lt.0) then
         write(6,61) ix
   61    format(/' *** error message from s.fact2 ***'
     &   /' negative index was detected.'
     &   /' ix =',i5)
         call parastop( 828 )
      endif
      fact2=1.
      if (ix.eq.0) go to 5
      do 100 i=0,ia-1
 100   fact2=fact2*(ix-i)
  5   return
      end


************************************************************************
*                                                                      *
      function sss(a,b,f,n)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      external f
c ----------------------------------------------------------------------
c
      istep=n/3-2
      rh=(b-a)/dble(n)
      s=0.0
c
      do 10 i=0,istep
         x0=a+rh*3*i
         x1=x0+rh
         x2=x1+rh
         x3=x2+rh
         s=s+3*(f(x1)+f(x2))+2*f(x3)
   10 continue
c
      s=s+3*(f(a+(n-2)*rh)+f(a+(n-1)*rh))
      sss=(f(a)+f(b)+s)*rh*3.0/8.0
c
      return
      end


************************************************************************
*                                                                      *
      function term2(g,ex,np,nh,alrg1,alrg2)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c
      ne=np+nh
      term2=((g*ex-alrg1)/(g*ex-alrg2))**(ne+1)
     &     *(np*(np-1)+4.*np*nh+nh*(nh-1))/(g*ex-alrg1)
      return
      end


************************************************************************
*                                                                      *
      function term3(g,ex,np,nh,alrg1,alrg2)
c ----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c
      ne=np+nh
      term3=((g*ex-alrg1)/(g*ex-alrg2))**(ne+1)
     &     *(np*nh*(ne+1)*(ne-2))/(g*ex-alrg1)**2
      return
      end


************************************************************************
*                                                                      *
      function dost(i,z)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension tt(4,7)

      data tt/0.36, 0.77,  0.08, 0.0,
     &       0.51, 0.81,  0.0,  0.0,
     &       0.60, 0.85, -0.06, 0.0,
     &       0.66, 0.89, -0.10, 0.0,
     &       0.68, 0.93, -0.10, 0.0,
     &       0.69, 0.97, -0.10, 0.0,
     &       0.69, 1.00, -0.10, 0.0/

c ----------------------------------------------------------------------
c
      if (z-70.0) 15,5,5
5     dost=tt(i,7)
10    return
15    if (z-10.0) 20,20,25
20    dost=tt(i,1)
      go to 10
25    n=.1*z+1.0
      x=10*n
      x=(x-z)*.1
      dost=x*tt(i,n-1)+(1.0-x)*tt(i,n)
      go to 10
      end


************************************************************************
*                                                                      *
      function tuyysh(a,z)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/inout/in,io
      dimension pz0(112),qn0(160),rh(3)
c
      data rh/0.0086,3.,4./
      data pz0/
     &  3.3764,-2.0706, 0.3421, 0.4083, 0.1452,-1.2844,-1.2340,-1.0438,
     &  0.0437, 0.4501, 0.4202,-0.0212,-0.3736,-1.1844,-1.4664,-0.8469,
     & -0.5253,-0.1130,-0.0585,-0.2347, 0.5049, 0.7244, 0.3841, 0.1588,
     & -0.5640,-0.9816,-1.5457,-2.1470,-1.9357,-0.9058,-0.5359, 0.1827,
     &  0.7647, 1.2352, 1.6049, 1.5653, 1.6677, 1.3279, 1.2463, 1.1141,
     &  1.3291, 0.9777, 0.9554, 0.2360,-0.1586,-0.9014,-1.2590,-2.1312,
     & -2.6974,-3.5037,-2.9306,-2.1833,-1.6754,-1.1022,-0.4575, 0.0086,
     &  0.6777, 1.1335, 1.7750, 1.9301, 2.4789, 2.4892, 2.7213, 2.4900,
     &  2.7562, 2.5067, 2.6205, 2.3627, 2.4455, 2.2247, 2.0710, 1.7852,
     &  1.3680, 0.5275,-0.2235,-0.6404,-1.1082,-1.4288,-1.8189,-2.1098,
     & -2.4707,-2.5869,-2.2005,-1.6725,-1.2519,-0.6844,-0.0149, 0.4435,
     &  1.1010, 1.4378, 1.8098, 1.9595, 2.1862, 2.2157, 2.5622, 2.5422,
     &  2.8144, 2.8661, 2.9730, 2.9568, 2.9825, 2.9703, 3.0822, 3.1122,
     &  3.2130, 3.1392, 3.4331, 2.9866, 2.4663, 1.5000, 1.0000, 0.6000/
c
      data (qn0(k),k=1,80)/
     &  3.3051,-2.1766, 0.2606, 0.3153, 0.3039,-1.2476,-1.3224,-0.7845,
     &  0.2746, 0.6473, 0.5921, 0.1284,-0.3480,-1.3306,-1.4650,-0.9938,
     & -0.7800,-0.3090,-0.2201,-0.2014, 0.6951, 1.0069, 0.9542, 0.6077,
     & -0.0837,-0.5971,-1.4749,-2.2756,-2.3044,-1.5800,-1.2231,-0.4875,
     &  0.0292, 0.4112, 0.9253, 1.0433, 1.5009, 1.3770, 1.5770, 1.5156,
     &  1.7003, 1.3946, 1.4593, 0.7689, 0.4636,-0.2965,-0.7784,-1.7119,
     & -2.4656,-3.5084,-3.2470,-2.3216,-1.8101,-1.1622,-0.6612,-0.0522,
     &  0.6513, 1.1860, 1.8995, 2.1728, 2.6235, 2.6382, 2.9261, 2.7763,
     &  3.0569, 2.8629, 3.0762, 2.7049, 2.7712, 2.2735, 2.2122, 1.6328,
     &  1.4487, 0.6925, 0.2889,-0.4018,-0.8149,-1.4666,-1.8890,-2.5454/
      data (qn0(k),k=81,160)/
     & -3.1376,-3.7687,-3.1836,-2.4571,-1.8545,-1.1582,-0.3941, 0.1315,
     &  0.8144, 1.1536, 1.4597, 1.8282, 2.2986, 2.4462, 2.6887, 2.8585,
     &  3.1322, 3.0510, 3.2439, 3.1636, 3.3765, 3.2778, 3.4899, 3.2868,
     &  3.4559, 2.9406, 3.3515, 2.4736, 1.9056, 1.1040, 0.9672, 0.6274,
     &  0.4052, 0.0463,-0.2541,-0.7419,-1.1454,-1.5727,-1.9297,-2.2692,
     & -2.6569,-2.9763,-3.3309,-3.6385,-4.0750,-4.1558,-3.5375,-2.9367,
     & -2.3591,-1.7801,-1.1584,-0.6248,-0.0099, 0.5738, 0.8592, 1.8652,
     &  2.1544, 2.3168, 2.3602, 2.4541, 2.5025, 2.5072, 2.5165, 2.4991,
     &  2.5762, 2.4918, 2.5224, 2.4905, 2.6218, 2.5551, 2.6352, 2.6116,
     &  2.6372, 2.5132, 2.4140, 2.4063, 2.5635, 2.4000, 2.5000, 2.4000/
c ----------------------------------------------------------------------
c
      iz=int(z)
      n=int(a-z)
c
      if(iz.gt.112.or.n.gt.160) then
         write(io,61) iz,n
   61    format(/' *** error message from s.tuyysh ***'
     &   /' requested atomic or mass number was not allowed in mass ',
     &    'formula.'
     &   /' atomic number (z)=',i5,'   neutron number (n)=',i5)
         call parastop( 814 )
      end if
      pz=(2*z/a)**(2./3.)*pz0(iz)
      qn=(2*n/a)**(2./3.)*qn0(n)
      pq=pz+qn
      x=pz*qn
      if(pq.gt.0) then
        alpha=0.003
      else
        alpha=0.
      end if
c
      tuyysh=pq
     &  -rh(1)*a**(2./3.)*(x-rh(3)+(x**2-2*rh(2)*x+rh(3)**2)**.5)
     &  -alpha*pq**2*z*z/a
      return
      end


************************************************************************
*                                                                      *
      block data qdres1
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/lvel/aprime(240)
c ----------------------------------------------------------------------
c
      data aprime/
     &  0.13, 0.25, 0.38, 0.50, 0.63, 0.75, 0.88, 1.00, 1.13, 1.25,
     &  1.38, 1.50, 1.63, 1.75, 1.88, 2.00, 2.13, 2.25, 2.38, 3.94,
     &  2.63, 2.75, 2.88, 3.55, 4.35, 3.25, 3.38, 3.96, 3.63, 3.75,
     &  3.88, 4.82, 4.44, 4.43, 4.43, 4.42, 4.63, 5.66, 5.81, 5.95,
     &  5.49, 6.18, 7.11, 6.96, 7.20, 7.73, 6.41, 6.85, 6.77, 6.91,
     &  7.26, 7.20, 6.86, 8.06, 7.81, 7.82, 8.41, 8.13, 7.19, 8.35,
     &  8.13, 8.02, 8.93, 8.90, 9.69, 9.65,10.55, 9.38, 9.72,10.66,
     & 11.98,12.76,12.10,12.86,13.03,12.81,12.54,12.65,12.00,12.69,
     & 14.05,13.33,13.28,13.23,13.17, 8.66,11.09,10.40,13.47,10.17,
     & 12.22,11.62,12.95,13.15,13.57,12.87,16.16,14.71,15.69,14.09,
     & 18.56,16.22,16.67,17.13,17.00,16.86,15.33,15.61,16.77,17.93,
     & 17.45,16.97,17.88,17.58,15.78,16.83,17.49,16.03,15.08,16.74,
     & 17.74,17.43,18.14,17.06,19.01,17.02,17.02,17.02,18.51,17.20,
     & 16.75,16.97,16.94,16.91,17.69,15.55,14.56,14.35,16.55,18.29,
     & 17.80,17.05,21.31,19.15,19.51,19.87,20.39,20.90,21.85,22.89,
     & 25.68,24.64,24.91,23.24,22.85,22.46,21.98,21.64,21.75,21.85,
     & 21.77,21.69,23.74,21.35,23.03,20.66,21.81,20.77,22.18,22.58,
     & 22.55,21.45,21.16,21.02,20.87,22.09,22.00,21.28,23.05,21.70,
     & 21.45,22.28,23.00,22.11,23.56,22.83,24.88,22.64,23.27,23.89,
     & 23.92,23.94,21.16,22.30,21.75,21.19,20.72,20.24,21.34,19.00,
     & 17.93,17.85,15.70,13.54,11.78,10.02,10.98,10.28,11.72,13.81,
     & 14.46,15.30,16.16,16.99,17.84,18.68,19.53,20.37,21.22,22.06,
     & 22.91,23.75,24.60,25.44,26.29,27.13,27.98,28.82,29.67,30.71,
     & 30.53,31.45,29.63,30.15,30.65,30.27,29.52,30.08,29.80,29.87/

      end

************************************************************************
*                                                                      *
      block data dres1

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common/prab/pp0(1001),pp1(1001),pp2(1001),
     &             rmass(300),alph(300),bet(300)

*-----------------------------------------------------------------------

      data ( bet(i), i = 1, 150 ) /
     &  5.9851295D-01,  4.3446904D-01,  3.5651958D-01,
     &  3.0811018D-01,  2.7411449D-01,  2.4846697D-01,
     &  2.2818285D-01,  2.1159416D-01,  1.9768369D-01,
     &  1.8579006D-01,  1.7546159D-01,  1.6637766D-01,
     &  1.5830314D-01,  1.5106136D-01,  1.4451629D-01,
     &  1.3856125D-01,  1.3311148D-01,  1.2809837D-01,
     &  1.2346560D-01,  1.1916685D-01,  1.1516327D-01,
     &  1.1142206D-01,  1.0791546D-01,  1.0461956D-01,
     &  1.0151374D-01,  9.8580360D-02,  9.5803559D-02,
     &  9.3169987D-02,  9.0667367D-02,  8.8285267D-02,
     &  8.6013973D-02,  8.3845198D-02,  8.1771255D-02,
     &  7.9785466D-02,  7.7881694D-02,  7.6054275D-02,
     &  7.4298143D-02,  7.2608948D-02,  7.0982277D-02,
     &  6.9414377D-02,  6.7901850D-02,  6.6441357D-02,
     &  6.5029979D-02,  6.3665092D-02,  6.2344000D-02,
     &  6.1064497D-02,  5.9824400D-02,  5.8621697D-02,
     &  5.7454597D-02,  5.6321200D-02,  5.5220097D-02,
     &  5.4149598D-02,  5.3108297D-02,  5.2095097D-02,
     &  5.1108498D-02,  5.0147396D-02,  4.9210798D-02,
     &  4.8297599D-02,  4.7406897D-02,  4.6537697D-02,
     &  4.5689099D-02,  4.4860397D-02,  4.4050798D-02,
     &  4.3259598D-02,  4.2485997D-02,  4.1729398D-02,
     &  4.0989198D-02,  4.0264796D-02,  3.9555699D-02,
     &  3.8861297D-02,  3.8181096D-02,  3.7514597D-02,
     &  3.6861397D-02,  3.6220998D-02,  3.5592996D-02,
     &  3.4976996D-02,  3.4372598D-02,  3.3779498D-02,
     &  3.3197299D-02,  3.2625698D-02,  3.2064397D-02,
     &  3.1512998D-02,  3.0971300D-02,  3.0438900D-02,
     &  2.9915698D-02,  2.9401299D-02,  2.8895497D-02,
     &  2.8398097D-02,  2.7908798D-02,  2.7427498D-02,
     &  2.6953798D-02,  2.6487596D-02,  2.6028797D-02,
     &  2.5576998D-02,  2.5132097D-02,  2.4693996D-02,
     &  2.4262499D-02,  2.3837399D-02,  2.3418598D-02,
     &  2.3005798D-02,  2.2598997D-02,  2.2197999D-02,
     &  2.1802697D-02,  2.1412998D-02,  2.1028597D-02,
     &  2.0649496D-02,  2.0275597D-02,  1.9906797D-02,
     &  1.9542798D-02,  1.9183699D-02,  1.8829297D-02,
     &  1.8479500D-02,  1.8134300D-02,  1.7793398D-02,
     &  1.7456897D-02,  1.7124597D-02,  1.6796496D-02,
     &  1.6472396D-02,  1.6152300D-02,  1.5836097D-02,
     &  1.5523598D-02,  1.5214998D-02,  1.4909998D-02,
     &  1.4608499D-02,  1.4310699D-02,  1.4016200D-02,
     &  1.3725098D-02,  1.3437398D-02,  1.3152897D-02,
     &  1.2871597D-02,  1.2593497D-02,  1.2318399D-02,
     &  1.2046300D-02,  1.1777200D-02,  1.1510998D-02,
     &  1.1247698D-02,  1.0987200D-02,  1.0729399D-02,
     &  1.0474399D-02,  1.0221999D-02,  9.9721998D-03,
     &  9.7249970D-03,  9.4803981D-03,  9.2381984D-03,
     &  8.9983977D-03,  8.7610967D-03,  8.5260980D-03,
     &  8.2933977D-03,  8.0629997D-03,  7.8348964D-03/
      data ( bet(i), i = 151, 300 ) /
     &  7.6088980D-03,  7.3851980D-03,  7.1635991D-03,
     &  6.9439970D-03,  6.7265965D-03,  6.5111965D-03,
     &  6.2977970D-03,  6.0863979D-03,  5.8768988D-03,
     &  5.6692995D-03,  5.4635964D-03,  5.2597970D-03,
     &  5.0577968D-03,  4.8575997D-03,  4.6591982D-03,
     &  4.4624992D-03,  4.2675994D-03,  4.0742978D-03,
     &  3.8828000D-03,  3.6927999D-03,  3.5044998D-03,
     &  3.3177999D-03,  3.1327000D-03,  2.9491999D-03,
     &  2.7671000D-03,  2.5865999D-03,  2.4075999D-03,
     &  2.2300999D-03,  2.0539998D-03,  1.8793000D-03,
     &  1.7060998D-03,  1.5341998D-03,  1.3636998D-03,
     &  1.1945998D-03,  1.0267999D-03,  8.6029992D-04,
     &  6.9519994D-04,  5.3129997D-04,  3.6869990D-04,
     &  2.0739998D-04,  4.7299996D-05, -1.1159999D-04,
     & -2.6929984D-04, -4.2579998D-04, -5.8109988D-04,
     & -7.3529989D-04, -8.8829990D-04, -1.0401998D-03,
     & -1.1908999D-03, -1.3406000D-03, -1.4891997D-03,
     & -1.6366998D-03, -1.7830997D-03, -1.9284999D-03,
     & -2.0728000D-03, -2.2160998D-03, -2.3583998D-03,
     & -2.4996998D-03, -2.6399998D-03, -2.7792999D-03,
     & -2.9175999D-03, -3.0549997D-03, -3.1914997D-03,
     & -3.3269997D-03, -3.4615998D-03, -3.5952998D-03,
     & -3.7279997D-03, -3.8598997D-03, -3.9908998D-03,
     & -4.1209981D-03, -4.2502992D-03, -4.3786988D-03,
     & -4.5062974D-03, -4.6329983D-03, -4.7588982D-03,
     & -4.8839971D-03, -5.0082989D-03, -5.1316991D-03,
     & -5.2543990D-03, -5.3762980D-03, -5.4974966D-03,
     & -5.6177974D-03, -5.7373978D-03, -5.8562979D-03,
     & -5.9743971D-03, -6.0917996D-03, -6.2083974D-03,
     & -6.3243993D-03, -6.4395964D-03, -6.5540969D-03,
     & -6.6678971D-03, -6.7810975D-03, -6.8934970D-03,
     & -7.0052966D-03, -7.1163997D-03, -7.2267987D-03,
     & -7.3365979D-03, -7.4456967D-03, -7.5541995D-03,
     & -7.6620989D-03, -7.7692978D-03, -7.8758970D-03,
     & -7.9818964D-03, -8.0872998D-03, -8.1919990D-03,
     & -8.2961991D-03, -8.3997995D-03, -8.5026994D-03,
     & -8.6050965D-03, -8.7068975D-03, -8.8081993D-03,
     & -8.9088976D-03, -9.0089962D-03, -9.1084987D-03,
     & -9.2074983D-03, -9.3059987D-03, -9.4038993D-03,
     & -9.5012970D-03, -9.5980987D-03, -9.6944980D-03,
     & -9.7902976D-03, -9.8854974D-03, -9.9802986D-03,
     & -1.0074597D-02, -1.0168299D-02, -1.0261599D-02,
     & -1.0354400D-02, -1.0446597D-02, -1.0538395D-02,
     & -1.0629699D-02, -1.0720596D-02, -1.0810897D-02,
     & -1.0900799D-02, -1.0990199D-02, -1.1079200D-02,
     & -1.1167698D-02, -1.1255797D-02, -1.1343397D-02,
     & -1.1430498D-02, -1.1517297D-02, -1.1603497D-02,
     & -1.1689398D-02, -1.1774797D-02, -1.1859797D-02,
     & -1.1944398D-02,  5.0781054D+17, -9.8159995D-54,
     &  0.0d0, 0.0d0, 0.0d0/

      data ( alph(i), i = 1, 150 )/
     &  2.6899996D+00,  2.2918415D+00,  2.0981865D+00,
     &  1.9758234D+00,  1.8886700D+00,  1.8221197D+00,
     &  1.7689228D+00,  1.7249994D+00,  1.6878471D+00,
     &  1.6558266D+00,  1.6278133D+00,  1.6030045D+00,
     &  1.5808096D+00,  1.5607824D+00,  1.5425758D+00,
     &  1.5259209D+00,  1.5105982D+00,  1.4964323D+00,
     &  1.4832792D+00,  1.4710178D+00,  1.4595480D+00,
     &  1.4487839D+00,  1.4386530D+00,  1.4290934D+00,
     &  1.4200497D+00,  1.4114771D+00,  1.4033327D+00,
     &  1.3955812D+00,  1.3881903D+00,  1.3811312D+00,
     &  1.3743792D+00,  1.3679113D+00,  1.3617077D+00,
     &  1.3557501D+00,  1.3500214D+00,  1.3445063D+00,
     &  1.3391924D+00,  1.3340664D+00,  1.3291178D+00,
     &  1.3243351D+00,  1.3197088D+00,  1.3152313D+00,
     &  1.3108931D+00,  1.3066874D+00,  1.3026075D+00,
     &  1.2986469D+00,  1.2947998D+00,  1.2910595D+00,
     &  1.2874222D+00,  1.2838821D+00,  1.2804356D+00,
     &  1.2770777D+00,  1.2738047D+00,  1.2706137D+00,
     &  1.2674999D+00,  1.2644606D+00,  1.2614937D+00,
     &  1.2585945D+00,  1.2557611D+00,  1.2529917D+00,
     &  1.2502832D+00,  1.2476330D+00,  1.2450390D+00,
     &  1.2424994D+00,  1.2400122D+00,  1.2375755D+00,
     &  1.2351875D+00,  1.2328472D+00,  1.2305517D+00,
     &  1.2283001D+00,  1.2260914D+00,  1.2239227D+00,
     &  1.2217951D+00,  1.2197056D+00,  1.2176533D+00,
     &  1.2156372D+00,  1.2136555D+00,  1.2117090D+00,
     &  1.2097950D+00,  1.2079124D+00,  1.2060623D+00,
     &  1.2042408D+00,  1.2024498D+00,  1.2006874D+00,
     &  1.1989517D+00,  1.1972437D+00,  1.1955624D+00,
     &  1.1939058D+00,  1.1922750D+00,  1.1906681D+00,
     &  1.1890850D+00,  1.1875248D+00,  1.1859865D+00,
     &  1.1844702D+00,  1.1829758D+00,  1.1815023D+00,
     &  1.1800489D+00,  1.1786146D+00,  1.1772003D+00,
     &  1.1758051D+00,  1.1744289D+00,  1.1730700D+00,
     &  1.1717281D+00,  1.1704044D+00,  1.1690979D+00,
     &  1.1678076D+00,  1.1665325D+00,  1.1652737D+00,
     &  1.1640310D+00,  1.1628027D+00,  1.1615896D+00,
     &  1.1603909D+00,  1.1592064D+00,  1.1580353D+00,
     &  1.1568785D+00,  1.1557350D+00,  1.1546040D+00,
     &  1.1534863D+00,  1.1523809D+00,  1.1512880D+00,
     &  1.1502066D+00,  1.1491375D+00,  1.1480799D+00,
     &  1.1470346D+00,  1.1459999D+00,  1.1449757D+00,
     &  1.1439629D+00,  1.1429596D+00,  1.1419678D+00,
     &  1.1409864D+00,  1.1400137D+00,  1.1390524D+00,
     &  1.1380997D+00,  1.1371565D+00,  1.1362228D+00,
     &  1.1352987D+00,  1.1343832D+00,  1.1334772D+00,
     &  1.1325788D+00,  1.1316900D+00,  1.1308088D+00,
     &  1.1299362D+00,  1.1290722D+00,  1.1282158D+00,
     &  1.1273670D+00,  1.1265268D+00,  1.1256933D+00,
     &  1.1248684D+00,  1.1240501D+00,  1.1232395D+00/
      data ( alph(i), i = 151, 300 )/
     &  1.1224356D+00,  1.1216393D+00,  1.1208496D+00,
     &  1.1200666D+00,  1.1192904D+00,  1.1185217D+00,
     &  1.1177588D+00,  1.1170025D+00,  1.1162519D+00,
     &  1.1155081D+00,  1.1147709D+00,  1.1140394D+00,
     &  1.1133137D+00,  1.1125946D+00,  1.1118803D+00,
     &  1.1111727D+00,  1.1104698D+00,  1.1097736D+00,
     &  1.1090822D+00,  1.1083965D+00,  1.1077156D+00,
     &  1.1070404D+00,  1.1063709D+00,  1.1057062D+00,
     &  1.1050463D+00,  1.1043911D+00,  1.1037416D+00,
     &  1.1030970D+00,  1.1024570D+00,  1.1018219D+00,
     &  1.1011906D+00,  1.1005650D+00,  1.0999432D+00,
     &  1.0993261D+00,  1.0987139D+00,  1.0981054D+00,
     &  1.0975018D+00,  1.0969028D+00,  1.0963068D+00,
     &  1.0957165D+00,  1.0951290D+00,  1.0945463D+00,
     &  1.0939674D+00,  1.0933924D+00,  1.0928221D+00,
     &  1.0922546D+00,  1.0916920D+00,  1.0911322D+00,
     &  1.0905771D+00,  1.0900249D+00,  1.0894766D+00,
     &  1.0889320D+00,  1.0883913D+00,  1.0878534D+00,
     &  1.0873194D+00,  1.0867891D+00,  1.0862617D+00,
     &  1.0857382D+00,  1.0852175D+00,  1.0847006D+00,
     &  1.0841875D+00,  1.0836763D+00,  1.0831690D+00,
     &  1.0826654D+00,  1.0821638D+00,  1.0816660D+00,
     &  1.0811710D+00,  1.0806799D+00,  1.0801907D+00,
     &  1.0797043D+00,  1.0792217D+00,  1.0787420D+00,
     &  1.0782642D+00,  1.0777903D+00,  1.0773191D+00,
     &  1.0768499D+00,  1.0763845D+00,  1.0759211D+00,
     &  1.0754604D+00,  1.0750027D+00,  1.0745478D+00,
     &  1.0740948D+00,  1.0736446D+00,  1.0731974D+00,
     &  1.0727520D+00,  1.0723104D+00,  1.0718699D+00,
     &  1.0714331D+00,  1.0709982D+00,  1.0705652D+00,
     &  1.0701351D+00,  1.0697069D+00,  1.0692816D+00,
     &  1.0688591D+00,  1.0684376D+00,  1.0680199D+00,
     &  1.0676031D+00,  1.0671892D+00,  1.0667772D+00,
     &  1.0663681D+00,  1.0659609D+00,  1.0655556D+00,
     &  1.0651522D+00,  1.0647516D+00,  1.0643520D+00,
     &  1.0639553D+00,  1.0635605D+00,  1.0631676D+00,
     &  1.0627775D+00,  1.0623884D+00,  1.0620022D+00,
     &  1.0616169D+00,  1.0612345D+00,  1.0608540D+00,
     &  1.0604744D+00,  1.0600977D+00,  1.0597229D+00,
     &  1.0593491D+00,  1.0589781D+00,  1.0586081D+00,
     &  1.0582409D+00,  1.0578747D+00,  1.0575104D+00,
     &  1.0571480D+00,  1.0567875D+00,  1.0564289D+00,
     &  1.0560713D+00,  1.0557165D+00,  1.0553627D+00,
     &  1.0550108D+00,  1.0546598D+00,  1.0543118D+00,
     &  1.0539637D+00,  1.0536184D+00,  1.0532751D+00,
     &  1.0529327D+00,  1.0525923D+00,  1.0522528D+00,
     &  1.0519152D+00,  1.0515795D+00,  1.0512457D+00,
     &  1.0509119D+00,  1.0505810D+00,  1.0502510D+00,
     &  1.0499229D+00,  3.9846233D+28,  2.1710545D-56,
     &  0.0d0, 0.0d0, 0.0d0/

      data ( rmass(i), i = 1, 150 )/
     &  9.9999988D-01,  1.2599201D+00,  1.4422493D+00,
     &  1.5874004D+00,  1.7099752D+00,  1.8171196D+00,
     &  1.9129305D+00,  2.0000000D+00,  2.0800829D+00,
     &  2.1544342D+00,  2.2239790D+00,  2.2894278D+00,
     &  2.3513336D+00,  2.4101419D+00,  2.4662113D+00,
     &  2.5198412D+00,  2.5712805D+00,  2.6207409D+00,
     &  2.6684008D+00,  2.7144165D+00,  2.7589235D+00,
     &  2.8020391D+00,  2.8438663D+00,  2.8844986D+00,
     &  2.9240170D+00,  2.9624958D+00,  2.9999990D+00,
     &  3.0365887D+00,  3.0723162D+00,  3.1072321D+00,
     &  3.1413803D+00,  3.1748018D+00,  3.2075338D+00,
     &  3.2396107D+00,  3.2710657D+00,  3.3019266D+00,
     &  3.3322210D+00,  3.3619747D+00,  3.3912106D+00,
     &  3.4199514D+00,  3.4482164D+00,  3.4760256D+00,
     &  3.5033970D+00,  3.5303478D+00,  3.5568924D+00,
     &  3.5830469D+00,  3.6088257D+00,  3.6342402D+00,
     &  3.6593046D+00,  3.6840305D+00,  3.7084293D+00,
     &  3.7325106D+00,  3.7562847D+00,  3.7797623D+00,
     &  3.8029518D+00,  3.8258619D+00,  3.8485003D+00,
     &  3.8708763D+00,  3.8929958D+00,  3.9148664D+00,
     &  3.9364967D+00,  3.9578905D+00,  3.9790564D+00,
     &  3.9999990D+00,  4.0207253D+00,  4.0412397D+00,
     &  4.0615473D+00,  4.0816545D+00,  4.1015654D+00,
     &  4.1212845D+00,  4.1408167D+00,  4.1601667D+00,
     &  4.1793385D+00,  4.1983356D+00,  4.2171621D+00,
     &  4.2358227D+00,  4.2543201D+00,  4.2726574D+00,
     &  4.2908392D+00,  4.3088684D+00,  4.3267479D+00,
     &  4.3444805D+00,  4.3620701D+00,  4.3795185D+00,
     &  4.3968287D+00,  4.4140043D+00,  4.4310465D+00,
     &  4.4479589D+00,  4.4647446D+00,  4.4814034D+00,
     &  4.4979410D+00,  4.5143566D+00,  4.5306540D+00,
     &  4.5468349D+00,  4.5629015D+00,  4.5788565D+00,
     &  4.5946999D+00,  4.6104355D+00,  4.6260643D+00,
     &  4.6415882D+00,  4.6570082D+00,  4.6723280D+00,
     &  4.6875477D+00,  4.7026682D+00,  4.7176933D+00,
     &  4.7326221D+00,  4.7474585D+00,  4.7622023D+00,
     &  4.7768555D+00,  4.7914190D+00,  4.8058949D+00,
     &  4.8202839D+00,  4.8345871D+00,  4.8488064D+00,
     &  4.8629436D+00,  4.8769979D+00,  4.8909721D+00,
     &  4.9048672D+00,  4.9186840D+00,  4.9324236D+00,
     &  4.9460869D+00,  4.9596748D+00,  4.9731894D+00,
     &  4.9866304D+00,  4.9999990D+00,  5.0132971D+00,
     &  5.0265245D+00,  5.0396833D+00,  5.0527735D+00,
     &  5.0657959D+00,  5.0787525D+00,  5.0916424D+00,
     &  5.1044674D+00,  5.1172295D+00,  5.1299267D+00,
     &  5.1425619D+00,  5.1551361D+00,  5.1676483D+00,
     &  5.1801004D+00,  5.1924934D+00,  5.2048273D+00,
     &  5.2171030D+00,  5.2293205D+00,  5.2414818D+00,
     &  5.2535868D+00,  5.2656364D+00,  5.2776308D+00,
     &  5.2895718D+00,  5.3014584D+00,  5.3132915D+00/
      data ( rmass(i), i = 151, 300 )/
     &  5.3250732D+00,  5.3368025D+00,  5.3484802D+00,
     &  5.3601074D+00,  5.3716841D+00,  5.3832121D+00,
     &  5.3946896D+00,  5.4061193D+00,  5.4175005D+00,
     &  5.4288340D+00,  5.4401207D+00,  5.4513607D+00,
     &  5.4625549D+00,  5.4737024D+00,  5.4848061D+00,
     &  5.4958639D+00,  5.5068779D+00,  5.5178471D+00,
     &  5.5287743D+00,  5.5396576D+00,  5.5504980D+00,
     &  5.5612965D+00,  5.5720539D+00,  5.5827694D+00,
     &  5.5934439D+00,  5.6040773D+00,  5.6146717D+00,
     &  5.6252251D+00,  5.6357393D+00,  5.6462154D+00,
     &  5.6566515D+00,  5.6670504D+00,  5.6774101D+00,
     &  5.6877327D+00,  5.6980181D+00,  5.7082663D+00,
     &  5.7184782D+00,  5.7286530D+00,  5.7387924D+00,
     &  5.7488966D+00,  5.7589645D+00,  5.7689972D+00,
     &  5.7789955D+00,  5.7889595D+00,  5.7988892D+00,
     &  5.8087845D+00,  5.8186464D+00,  5.8284760D+00,
     &  5.8382721D+00,  5.8480349D+00,  5.8577652D+00,
     &  5.8674631D+00,  5.8771296D+00,  5.8867645D+00,
     &  5.8963680D+00,  5.9059401D+00,  5.9154806D+00,
     &  5.9249916D+00,  5.9344711D+00,  5.9439211D+00,
     &  5.9533405D+00,  5.9627314D+00,  5.9720917D+00,
     &  5.9814234D+00,  5.9907255D+00,  5.9999990D+00,
     &  6.0092440D+00,  6.0184603D+00,  6.0276489D+00,
     &  6.0368099D+00,  6.0459423D+00,  6.0550480D+00,
     &  6.0641260D+00,  6.0731773D+00,  6.0822010D+00,
     &  6.0911980D+00,  6.1001692D+00,  6.1091137D+00,
     &  6.1180325D+00,  6.1269245D+00,  6.1357918D+00,
     &  6.1446323D+00,  6.1534481D+00,  6.1622391D+00,
     &  6.1710043D+00,  6.1797457D+00,  6.1884623D+00,
     &  6.1971531D+00,  6.2058210D+00,  6.2144642D+00,
     &  6.2230835D+00,  6.2316790D+00,  6.2402506D+00,
     &  6.2487984D+00,  6.2573233D+00,  6.2658253D+00,
     &  6.2743044D+00,  6.2827606D+00,  6.2911940D+00,
     &  6.2996044D+00,  6.3079920D+00,  6.3163586D+00,
     &  6.3247023D+00,  6.3330250D+00,  6.3413248D+00,
     &  6.3496037D+00,  6.3578606D+00,  6.3660955D+00,
     &  6.3743105D+00,  6.3825035D+00,  6.3906755D+00,
     &  6.3988266D+00,  6.4069576D+00,  6.4150677D+00,
     &  6.4231577D+00,  6.4312267D+00,  6.4392757D+00,
     &  6.4473047D+00,  6.4553137D+00,  6.4633036D+00,
     &  6.4712725D+00,  6.4792223D+00,  6.4871531D+00,
     &  6.4950638D+00,  6.5029564D+00,  6.5108290D+00,
     &  6.5186834D+00,  6.5265179D+00,  6.5343342D+00,
     &  6.5421314D+00,  6.5499105D+00,  6.5576715D+00,
     &  6.5654135D+00,  6.5731373D+00,  6.5808430D+00,
     &  6.5885315D+00,  6.5962009D+00,  6.6038532D+00,
     &  6.6114883D+00,  6.6191053D+00,  6.6267042D+00,
     &  6.6342859D+00,  6.6418514D+00,  6.6493988D+00,
     &  6.6569290D+00, -5.5638099D-67,  3.4426720D+07,
     &  0.0d0, 0.0d0, 0.0d0/


      data ( pp0(i),pp1(i),pp2(i),i=1,50 )/
     &  0.0000000D+00,  0.0000000D+00,  0.0000000D+00,
     &  5.1557212D-25,  6.3526568D-28,  0.0000000D+00,
     &  2.2062775D-24,  1.0740024D-26,  3.2685078D-03,
     &  5.3135927D-24,  5.7417734D-26,  7.3648617D-03,
     &  1.0116792D-23,  1.9173910D-25,  1.2902647D-02,
     &  1.6938207D-23,  4.9483238D-25,  2.0033035D-02,
     &  2.6149250D-23,  1.0851791D-24,  2.8619513D-02,
     &  3.8177107D-23,  2.1272313D-24,  3.8642701D-02,
     &  5.3512361D-23,  3.8416224D-24,  5.0061021D-02,
     &  7.2717678D-23,  6.5172178D-24,  6.2842131D-02,
     &  9.6437476D-23,  1.0525275D-23,  7.6943576D-02,
     &  1.2540889D-22,  1.6336072D-23,  9.2329681D-02,
     &  1.6047421D-22,  2.4538421D-23,  1.0896289D-01,
     &  2.0259478D-22,  3.5862412D-23,  1.2680721D-01,
     &  2.5286665D-22,  5.1205984D-23,  1.4582402D-01,
     &  3.1253876D-22,  7.1665893D-23,  1.6597921D-01,
     &  3.8303208D-22,  9.8573787D-23,  1.8723571D-01,
     &  4.6596300D-22,  1.3353783D-22,  2.0955956D-01,
     &  5.6316762D-22,  1.7849147D-22,  2.3291522D-01,
     &  6.7673108D-22,  2.3574910D-22,  2.5726932D-01,
     &  8.0901919D-22,  3.0807207D-22,  2.8258818D-01,
     &  9.6271464D-22,  3.9874225D-22,  3.0883926D-01,
     &  1.1408568D-21,  5.1164958D-22,  3.3599013D-01,
     &  1.3468877D-21,  6.5139221D-22,  3.6400956D-01,
     &  1.5847034D-21,  8.2339045D-22,  3.9286703D-01,
     &  1.8587106D-21,  1.0340200D-21,  4.2253220D-01,
     &  2.1738914D-21,  1.2907634D-21,  4.5297587D-01,
     &  2.5358756D-21,  1.6023862D-21,  4.8416930D-01,
     &  2.9510216D-21,  1.9791346D-21,  5.1608509D-01,
     &  3.4265052D-21,  2.4329690D-21,  5.4869568D-01,
     &  3.9704314D-21,  2.9778252D-21,  5.8197492D-01,
     &  4.5919348D-21,  3.6299141D-21,  6.1589724D-01,
     &  5.3013165D-21,  4.4080778D-21,  6.5043777D-01,
     &  6.1101861D-21,  5.3341679D-21,  6.8557239D-01,
     &  7.0316233D-21,  6.4335071D-21,  7.2127777D-01,
     &  8.0803623D-21,  7.7354086D-21,  7.5753140D-01,
     &  9.2729822D-21,  9.2737545D-21,  7.9431129D-01,
     &  1.0628153D-20,  1.1087674D-20,  8.3159631D-01,
     &  1.2166867D-20,  1.3222304D-20,  8.6936599D-01,
     &  1.3912751D-20,  1.5729667D-20,  9.0760070D-01,
     &  1.5892356D-20,  1.8669664D-20,  9.4628143D-01,
     &  1.8135534D-20,  2.2111203D-20,  9.8538983D-01,
     &  2.0675825D-20,  2.6133489D-20,  1.0249081D+00,
     &  2.3550921D-20,  3.0827490D-20,  1.0648193D+00,
     &  2.6803140D-20,  3.6297596D-20,  1.1051064D+00,
     &  3.0480013D-20,  4.2663530D-20,  1.1457558D+00,
     &  3.4634879D-20,  5.0062473D-20,  1.1867495D+00,
     &  3.9327620D-20,  5.8651500D-20,  1.2280750D+00,
     &  4.4625412D-20,  6.8610444D-20,  1.2697182D+00,
     &  5.0603621D-20,  8.0144857D-20,  1.3116646D+00/
      data ( pp0(i),pp1(i),pp2(i),i=51,100 )/
     &  5.7346726D-20,  9.3489814D-20,  1.3539019D+00,
     &  6.4949601D-20,  1.0891373D-19,  1.3964176D+00,
     &  7.3518417D-20,  1.2672309D-19,  1.4392004D+00,
     &  8.3172235D-20,  1.4726762D-19,  1.4822388D+00,
     &  9.4044594D-20,  1.7094618D-19,  1.5255213D+00,
     &  1.0628495D-19,  1.9821326D-19,  1.5690384D+00,
     &  1.2006082D-19,  2.2958687D-19,  1.6127796D+00,
     &  1.3555975D-19,  2.6565681D-19,  1.6567354D+00,
     &  1.5299191D-19,  3.0709432D-19,  1.7008972D+00,
     &  1.7259237D-19,  3.5466325D-19,  1.7452564D+00,
     &  1.9462448D-19,  4.0923210D-19,  1.7898035D+00,
     &  2.1938293D-19,  4.7178830D-19,  1.8345318D+00,
     &  2.4719739D-19,  5.4345379D-19,  1.8794327D+00,
     &  2.7843696D-19,  6.2550276D-19,  1.9244995D+00,
     &  3.1351433D-19,  7.1938223D-19,  1.9697256D+00,
     &  3.5289117D-19,  8.2673425D-19,  2.0151043D+00,
     &  3.9708397D-19,  9.4942136D-19,  2.0606279D+00,
     &  4.4667014D-19,  1.0895582D-18,  2.1062927D+00,
     &  5.0229529D-19,  1.2495399D-18,  2.1520910D+00,
     &  5.6468120D-19,  1.4320813D-18,  2.1980181D+00,
     &  6.3463473D-19,  1.6402603D-18,  2.2440691D+00,
     &  7.1305761D-19,  1.8775602D-18,  2.2902384D+00,
     &  8.0095744D-19,  2.1479283D-18,  2.3365211D+00,
     &  8.9945965D-19,  2.4558314D-18,  2.3829126D+00,
     &  1.0098237D-18,  2.8063244D-18,  2.4294090D+00,
     &  1.1334542D-18,  3.2051264D-18,  2.4760065D+00,
     &  1.2719201D-18,  3.6587050D-18,  2.5227003D+00,
     &  1.4269759D-18,  4.1743744D-18,  2.5694866D+00,
     &  1.6005788D-18,  4.7604029D-18,  2.6163626D+00,
     &  1.7949149D-18,  5.4261328D-18,  2.6633244D+00,
     &  2.0124254D-18,  6.1821204D-18,  2.7103691D+00,
     &  2.2558340D-18,  7.0402905D-18,  2.7574921D+00,
     &  2.5281841D-18,  8.0141088D-18,  2.8046923D+00,
     &  2.8328703D-18,  9.1187779D-18,  2.8519659D+00,
     &  3.1736819D-18,  1.0371459D-17,  2.8993101D+00,
     &  3.5548467D-18,  1.1791518D-17,  2.9467230D+00,
     &  3.9810846D-18,  1.3400803D-17,  2.9942007D+00,
     &  4.4576590D-18,  1.5223961D-17,  3.0417414D+00,
     &  4.9904427D-18,  1.7288789D-17,  3.0893440D+00,
     &  5.5859896D-18,  1.9626627D-17,  3.1370049D+00,
     &  6.2516061D-18,  2.2272798D-17,  3.1847229D+00,
     &  6.9954449D-18,  2.5267138D-17,  3.2324953D+00,
     &  7.8265968D-18,  2.8654503D-17,  3.2803202D+00,
     &  8.7551989D-18,  3.2485461D-17,  3.3281965D+00,
     &  9.7925570D-18,  3.6816949D-17,  3.3761225D+00,
     &  1.0951275D-17,  4.1713117D-17,  3.4240952D+00,
     &  1.2245407D-17,  4.7246174D-17,  3.4721146D+00,
     &  1.3690618D-17,  5.3497423D-17,  3.5201778D+00,
     &  1.5304363D-17,  6.0558382D-17,  3.5682850D+00,
     &  1.7106135D-17,  6.8532073D-17,  3.6164322D+00/
      data ( pp0(i),pp1(i),pp2(i),i=101,150 )/
     &  1.9117613D-17,  7.7534378D-17,  3.6646204D+00,
     &  2.1362872D-17,  8.7693427D-17,  3.7126265D+00,
     &  2.3869071D-17,  9.9160427D-17,  3.7609081D+00,
     &  2.6666146D-17,  1.1209811D-16,  3.8092251D+00,
     &  2.9787528D-17,  1.2669203D-16,  3.8575764D+00,
     &  3.3270528D-17,  1.4315077D-16,  3.9059610D+00,
     &  3.7156662D-17,  1.6170883D-16,  3.9543772D+00,
     &  4.1492213D-17,  1.8262985D-16,  4.0028257D+00,
     &  4.6328732D-17,  2.0621003D-16,  4.0513039D+00,
     &  5.1723617D-17,  2.3278223D-16,  4.0998125D+00,
     &  5.7740833D-17,  2.6272061D-16,  4.1483488D+00,
     &  6.4451624D-17,  2.9644544D-16,  4.1969137D+00,
     &  7.1935305D-17,  3.3442830D-16,  4.2455053D+00,
     &  8.0280221D-17,  3.7719987D-16,  4.2941236D+00,
     &  8.9584746D-17,  4.2535496D-16,  4.3427677D+00,
     &  9.9958438D-17,  4.7956210D-16,  4.3914366D+00,
     &  1.1152326D-16,  5.4057177D-16,  4.4401293D+00,
     &  1.2441499D-16,  6.0922612D-16,  4.4888458D+00,
     &  1.3878488D-16,  6.8647065D-16,  4.5375853D+00,
     &  1.5480119D-16,  7.7336628D-16,  4.5863476D+00,
     &  1.7265131D-16,  8.7110372D-16,  4.6351318D+00,
     &  1.9254381D-16,  9.8101874D-16,  4.6839371D+00,
     &  2.1471086D-16,  1.1046099D-15,  4.7327633D+00,
     &  2.3941069D-16,  1.2435585D-15,  4.7816095D+00,
     &  2.6693121D-16,  1.3997499D-15,  4.8304749D+00,
     &  2.9759253D-16,  1.5752987D-15,  4.8793602D+00,
     &  3.3175041D-16,  1.7725761D-15,  4.9282637D+00,
     &  3.6980167D-16,  1.9942412D-15,  4.9771852D+00,
     &  4.1218698D-16,  2.2432746D-15,  5.0261250D+00,
     &  4.5939764D-16,  2.5230183D-15,  5.0750828D+00,
     &  5.1197932D-16,  2.8372175D-15,  5.1240568D+00,
     &  5.7053979D-16,  3.1900705D-15,  5.1730471D+00,
     &  6.3575498D-16,  3.5862799D-15,  5.2220535D+00,
     &  7.0837704D-16,  4.0311213D-15,  5.2710762D+00,
     &  7.8924264D-16,  4.5305048D-15,  5.3201141D+00,
     &  8.7928203D-16,  5.0910407D-15,  5.3691673D+00,
     &  9.7953092D-16,  5.7201456D-15,  5.4182348D+00,
     &  1.0911406D-15,  6.4261273D-15,  5.4673166D+00,
     &  1.2153921D-15,  7.2182860D-15,  5.5164127D+00,
     &  1.3537096D-15,  8.1070370D-15,  5.5655222D+00,
     &  1.5076780D-15,  9.1040490D-15,  5.6146450D+00,
     &  1.6790590D-15,  1.0222387D-14,  5.6637821D+00,
     &  1.8698123D-15,  1.1476673D-14,  5.7129307D+00,
     &  2.0821171D-15,  1.2883286D-14,  5.7620926D+00,
     &  2.3183965D-15,  1.4460557D-14,  5.8112659D+00,
     &  2.5813452D-15,  1.6229002D-14,  5.8604517D+00,
     &  2.8739595D-15,  1.8211591D-14,  5.9096498D+00,
     &  3.1995706D-15,  2.0434033D-14,  5.9588585D+00,
     &  3.5618820D-15,  2.2925086D-14,  6.0080795D+00,
     &  3.9650120D-15,  2.5716950D-14,  6.0573101D+00/
      data ( pp0(i),pp1(i),pp2(i),i=151,200 )/
     &  4.4135431D-15,  2.8845639D-14,  6.1065531D+00,
     &  4.9125573D-15,  3.2351451D-14,  6.1558056D+00,
     &  5.4677230D-15,  3.6279471D-14,  6.2050686D+00,
     &  6.0853219D-15,  4.0680140D-14,  6.2543421D+00,
     &  6.7723503D-15,  4.5609871D-14,  6.3036251D+00,
     &  7.5365807D-15,  5.1131771D-14,  6.3529177D+00,
     &  8.3866494D-15,  5.7316402D-14,  6.4022198D+00,
     &  9.3321652D-15,  6.4242719D-14,  6.4515314D+00,
     &  1.0383797D-14,  7.1998939D-14,  6.5008526D+00,
     &  1.1553411D-14,  8.0683832D-14,  6.5501823D+00,
     &  1.2854186D-14,  9.0407716D-14,  6.5995207D+00,
     &  1.4300779D-14,  1.0129392D-13,  6.6488676D+00,
     &  1.5909471D-14,  1.1348051D-13,  6.6982241D+00,
     &  1.7698360D-14,  1.2712157D-13,  6.7475882D+00,
     &  1.9687553D-14,  1.4238946D-13,  6.7969599D+00,
     &  2.1899393D-14,  1.5947692D-13,  6.8463402D+00,
     &  2.4358716D-14,  1.7859921D-13,  6.8957281D+00,
     &  2.7093111D-14,  1.9999703D-13,  6.9451237D+00,
     &  3.0133231D-14,  2.2393936D-13,  6.9945269D+00,
     &  3.3513143D-14,  2.5072658D-13,  7.0439377D+00,
     &  3.7270690D-14,  2.8069468D-13,  7.0933552D+00,
     &  4.1447917D-14,  3.1421876D-13,  7.1427803D+00,
     &  4.6091549D-14,  3.5171806D-13,  7.1922121D+00,
     &  5.1253469D-14,  3.9366096D-13,  7.2416515D+00,
     &  5.6991304D-14,  4.4057059D-13,  7.2910967D+00,
     &  6.3369232D-14,  4.9303145D-13,  7.3405495D+00,
     &  7.0458233D-14,  5.5169633D-13,  7.3900080D+00,
     &  7.8337510D-14,  6.1729446D-13,  7.4394732D+00,
     &  8.7094719D-14,  6.9064014D-13,  7.4889441D+00,
     &  9.6827494D-14,  7.7264296D-13,  7.5384226D+00,
     &  1.0764415D-13,  8.6431865D-13,  7.5879059D+00,
     &  1.1966491D-13,  9.6680042D-13,  7.6373959D+00,
     &  1.3302358D-13,  1.0813572D-12,  7.6868916D+00,
     &  1.4786848D-13,  1.2094006D-12,  7.7363930D+00,
     &  1.6436456D-13,  1.3525119D-12,  7.7858992D+00,
     &  1.8269490D-13,  1.5124516D-12,  7.8354120D+00,
     &  2.0306283D-13,  1.6911889D-12,  7.8849297D+00,
     &  2.2569425D-13,  1.8909206D-12,  7.9344540D+00,
     &  2.5083988D-13,  2.1140997D-12,  7.9839821D+00,
     &  2.7877836D-13,  2.3634627D-12,  8.0335159D+00,
     &  3.0981896D-13,  2.6420654D-12,  8.0830555D+00,
     &  3.4430510D-13,  2.9533173D-12,  8.1325989D+00,
     &  3.8261825D-13,  3.3010261D-12,  8.1821480D+00,
     &  4.2518192D-13,  3.6894376D-12,  8.2317019D+00,
     &  4.7246636D-13,  4.1232929D-12,  8.2812614D+00,
     &  5.2499384D-13,  4.6078818D-12,  8.3308249D+00,
     &  5.8334414D-13,  5.1491060D-12,  8.3803921D+00,
     &  6.4816105D-13,  5.7535513D-12,  8.4299650D+00,
     &  7.2015928D-13,  6.4285660D-12,  8.4795427D+00,
     &  8.0013253D-13,  7.1823485D-12,  8.5291243D+00/
      data ( pp0(i),pp1(i),pp2(i),i=201,250 )/
     &  8.8896186D-13,  8.0240441D-12,  8.5787106D+00,
     &  9.8762491D-13,  8.9638583D-12,  8.6283007D+00,
     &  1.0972091D-12,  1.0013173D-11,  8.6778946D+00,
     &  1.2189191D-12,  1.1184683D-11,  8.7274933D+00,
     &  1.3540939D-12,  1.2492556D-11,  8.7770958D+00,
     &  1.5042186D-12,  1.3952588D-11,  8.8267031D+00,
     &  1.6709438D-12,  1.5582396D-11,  8.8763142D+00,
     &  1.8561003D-12,  1.7401636D-11,  8.9259281D+00,
     &  2.0617206D-12,  1.9432220D-11,  8.9755468D+00,
     &  2.2900614D-12,  2.1698601D-11,  9.0251694D+00,
     &  2.5436276D-12,  2.4228036D-11,  9.0747957D+00,
     &  2.8251984D-12,  2.7050903D-11,  9.1244249D+00,
     &  3.1378615D-12,  3.0201106D-11,  9.1740589D+00,
     &  3.4850404D-12,  3.3716446D-11,  9.2236958D+00,
     &  3.8705384D-12,  3.7639059D-11,  9.2733364D+00,
     &  4.2985745D-12,  4.2015905D-11,  9.3229799D+00,
     &  4.7738332D-12,  4.6899401D-11,  9.3726273D+00,
     &  5.3015118D-12,  5.2347904D-11,  9.4222784D+00,
     &  5.8873809D-12,  5.8426555D-11,  9.4719324D+00,
     &  6.5378441D-12,  6.5207895D-11,  9.5215902D+00,
     &  7.2600051D-12,  7.2772857D-11,  9.5712509D+00,
     &  8.0617535D-12,  8.1211593D-11,  9.6209145D+00,
     &  8.9518401D-12,  9.0624647D-11,  9.6705818D+00,
     &  9.9399794D-12,  1.0112404D-10,  9.7202520D+00,
     &  1.1036950D-11,  1.1283467D-10,  9.7699251D+00,
     &  1.2254713D-11,  1.2589570D-10,  9.8196011D+00,
     &  1.3606545D-11,  1.4046227D-10,  9.8692799D+00,
     &  1.5107166D-11,  1.5670722D-10,  9.9189625D+00,
     &  1.6772944D-11,  1.7482323D-10,  9.9686470D+00,
     &  1.8621993D-11,  1.9502497D-10,  1.0018334D+01,
     &  2.0674462D-11,  2.1755166D-10,  1.0068026D+01,
     &  2.2952668D-11,  2.4266966D-10,  1.0117719D+01,
     &  2.5481409D-11,  2.7067659D-10,  1.0167415D+01,
     &  2.8288163D-11,  3.0190295D-10,  1.0217114D+01,
     &  3.1403449D-11,  3.3671754D-10,  1.0266816D+01,
     &  3.4861114D-11,  3.7553138D-10,  1.0316520D+01,
     &  3.8698711D-11,  4.1880188D-10,  1.0366226D+01,
     &  4.2957929D-11,  4.6703930D-10,  1.0415936D+01,
     &  4.7684995D-11,  5.2081184D-10,  1.0465649D+01,
     &  5.2931215D-11,  5.8075189D-10,  1.0515363D+01,
     &  5.8753474D-11,  6.4756489D-10,  1.0565080D+01,
     &  6.5214945D-11,  7.2203621D-10,  1.0614799D+01,
     &  7.2385667D-11,  8.0504026D-10,  1.0664520D+01,
     &  8.0343343D-11,  8.9755159D-10,  1.0714245D+01,
     &  8.9174210D-11,  1.0006558D-09,  1.0763971D+01,
     &  9.8973899D-11,  1.1155612D-09,  1.0813700D+01,
     &  1.0984852D-10,  1.2436145D-09,  1.0863432D+01,
     &  1.2191578D-10,  1.3863153D-09,  1.0913163D+01,
     &  1.3530629D-10,  1.5453332D-09,  1.0962899D+01,
     &  1.5016487D-10,  1.7225279D-09,  1.1012637D+01/
      data ( pp0(i),pp1(i),pp2(i),i=251,300 )/
     &  1.6665222D-10,  1.9199708D-09,  1.1062378D+01,
     &  1.8494657D-10,  2.1399682D-09,  1.1112120D+01,
     &  2.0524568D-10,  2.3850886D-09,  1.1161863D+01,
     &  2.2776886D-10,  2.6581917D-09,  1.1211610D+01,
     &  2.5275937D-10,  2.9624621D-09,  1.1261357D+01,
     &  2.8048697D-10,  3.3014456D-09,  1.1311109D+01,
     &  3.1125147D-10,  3.6790906D-09,  1.1360861D+01,
     &  3.4538417D-10,  4.0997925D-09,  1.1410616D+01,
     &  3.8325387D-10,  4.5684452D-09,  1.1460372D+01,
     &  4.2526893D-10,  5.0904987D-09,  1.1510130D+01,
     &  4.7188253D-10,  5.6720211D-09,  1.1559890D+01,
     &  5.2359672D-10,  6.3197660D-09,  1.1609653D+01,
     &  5.8096927D-10,  7.0412511D-09,  1.1659416D+01,
     &  6.4461836D-10,  7.8448430D-09,  1.1709182D+01,
     &  7.1522921D-10,  8.7398675D-09,  1.1758949D+01,
     &  7.9356233D-10,  9.7366915D-09,  1.1808719D+01,
     &  8.8046126D-10,  1.0846861D-08,  1.1858491D+01,
     &  9.7686126D-10,  1.2083230D-08,  1.1908263D+01,
     &  1.0837993D-09,  1.3460102D-08,  1.1958038D+01,
     &  1.2024259D-09,  1.4993404D-08,  1.2007814D+01,
     &  1.3340169D-09,  1.6700849D-08,  1.2057592D+01,
     &  1.4799868D-09,  1.8602172D-08,  1.2107372D+01,
     &  1.6419053D-09,  2.0719323D-08,  1.2157154D+01,
     &  1.8215116D-09,  2.3076733D-08,  1.2206938D+01,
     &  2.0207356D-09,  2.5701596D-08,  1.2256722D+01,
     &  2.2417175D-09,  2.8624171D-08,  1.2306510D+01,
     &  2.4868296D-09,  3.1878137D-08,  1.2356297D+01,
     &  2.7587039D-09,  3.5500978D-08,  1.2406088D+01,
     &  3.0602578D-09,  3.9534388D-08,  1.2455877D+01,
     &  3.3947274D-09,  4.4024784D-08,  1.2505671D+01,
     &  3.7656989D-09,  4.9023807D-08,  1.2555466D+01,
     &  4.1771564D-09,  5.4588916D-08,  1.2605262D+01,
     &  4.6335060D-09,  6.0784032D-08,  1.2655059D+01,
     &  5.1396398D-09,  6.7680332D-08,  1.2704857D+01,
     &  5.7009863D-09,  7.5356979D-08,  1.2754658D+01,
     &  6.3235603D-09,  8.3902080D-08,  1.2804460D+01,
     &  7.0140267D-09,  9.3413576D-08,  1.2854261D+01,
     &  7.7797857D-09,  1.0400043D-07,  1.2904066D+01,
     &  8.6290299D-09,  1.1578413D-07,  1.2953873D+01,
     &  9.5708579D-09,  1.2889944D-07,  1.3003680D+01,
     &  1.0615345D-08,  1.4349661D-07,  1.3053489D+01,
     &  1.1773668D-08,  1.5974257D-07,  1.3103299D+01,
     &  1.3058219D-08,  1.7782315D-07,  1.3153111D+01,
     &  1.4482740D-08,  1.9794510D-07,  1.3202924D+01,
     &  1.6062458D-08,  2.2033828D-07,  1.3252736D+01,
     &  1.7814266D-08,  2.4525849D-07,  1.3302552D+01,
     &  1.9756889D-08,  2.7299023D-07,  1.3352369D+01,
     &  2.1911081D-08,  3.0384996D-07,  1.3402187D+01,
     &  2.4299862D-08,  3.3818969D-07,  1.3452006D+01,
     &  2.6948747D-08,  3.7640098D-07,  1.3501826D+01/
      data ( pp0(i),pp1(i),pp2(i),i=301,350 )/
     &  2.9886024D-08,  4.1891929D-07,  1.3551648D+01,
     &  3.3143056D-08,  4.6622904D-07,  1.3601471D+01,
     &  3.6754610D-08,  5.1886894D-07,  1.3651296D+01,
     &  4.0759232D-08,  5.7743824D-07,  1.3701120D+01,
     &  4.5199656D-08,  6.4260331D-07,  1.3750947D+01,
     &  5.0123248D-08,  7.1510527D-07,  1.3800775D+01,
     &  5.5582532D-08,  7.9576847D-07,  1.3850603D+01,
     &  6.1635717D-08,  8.8550951D-07,  1.3900433D+01,
     &  6.8347333D-08,  9.8534747D-07,  1.3950264D+01,
     &  7.5788932D-08,  1.0964168D-06,  1.4000096D+01,
     &  8.4039868D-08,  1.2199780D-06,  1.4049930D+01,
     &  9.3188021D-08,  1.3574327D-06,  1.4099765D+01,
     &  1.0333082D-07,  1.5103396D-06,  1.4149599D+01,
     &  1.1457632D-07,  1.6804324D-06,  1.4199435D+01,
     &  1.2704436D-07,  1.8696401D-06,  1.4249274D+01,
     &  1.4086760D-07,  2.0801044D-06,  1.4299112D+01,
     &  1.5619315D-07,  2.3142093D-06,  1.4348952D+01,
     &  1.7318428D-07,  2.5746049D-06,  1.4398793D+01,
     &  1.9202162D-07,  2.8642371D-06,  1.4448635D+01,
     &  2.1290572D-07,  3.1863829D-06,  1.4498478D+01,
     &  2.3605867D-07,  3.5446847D-06,  1.4548323D+01,
     &  2.6172665D-07,  3.9431907D-06,  1.4598166D+01,
     &  2.9018275D-07,  4.3864056D-06,  1.4648012D+01,
     &  3.2172937D-07,  4.8793354D-06,  1.4697859D+01,
     &  3.5670195D-07,  5.4275451D-06,  1.4747707D+01,
     &  3.9547200D-07,  6.0372204D-06,  1.4797556D+01,
     &  4.3845165D-07,  6.7152423D-06,  1.4847405D+01,
     &  4.8609741D-07,  7.4692562D-06,  1.4897257D+01,
     &  5.3891534D-07,  8.3077630D-06,  1.4947107D+01,
     &  5.9746645D-07,  9.2402142D-06,  1.4996961D+01,
     &  6.6237243D-07,  1.0277113D-05,  1.5046815D+01,
     &  7.3432227D-07,  1.1430139D-05,  1.5096668D+01,
     &  8.1407973D-07,  1.2712274D-05,  1.5146523D+01,
     &  9.0249125D-07,  1.4137947D-05,  1.5196380D+01,
     &  1.0004942D-06,  1.5723184D-05,  1.5246237D+01,
     &  1.1091297D-06,  1.7485843D-05,  1.5296096D+01,
     &  1.2295495D-06,  1.9445724D-05,  1.5345954D+01,
     &  1.3630306D-06,  2.1624859D-05,  1.5395814D+01,
     &  1.5109890D-06,  2.4047738D-05,  1.5445674D+01,
     &  1.6749918D-06,  2.6741560D-05,  1.5495536D+01,
     &  1.8567789D-06,  2.9736562D-05,  1.5545398D+01,
     &  2.0582766D-06,  3.3066390D-05,  1.5595261D+01,
     &  2.2816193D-06,  3.6768397D-05,  1.5645124D+01,
     &  2.5291738D-06,  4.0884101D-05,  1.5694988D+01,
     &  2.8035629D-06,  4.5459659D-05,  1.5744855D+01,
     &  3.1076934D-06,  5.0546354D-05,  1.5794721D+01,
     &  3.4447839D-06,  5.6201199D-05,  1.5844587D+01,
     &  3.8184053D-06,  6.2487525D-05,  1.5894456D+01,
     &  4.2325128D-06,  6.9475762D-05,  1.5944324D+01,
     &  4.6914884D-06,  7.7244113D-05,  1.5994193D+01/
      data ( pp0(i),pp1(i),pp2(i),i=351,400 )/
     &  5.2001906D-06,  8.5879554D-05,  1.6044052D+01,
     &  5.7640027D-06,  9.5478681D-05,  1.6093933D+01,
     &  6.3888883D-06,  1.0614886D-04,  1.6143799D+01,
     &  7.0814585D-06,  1.1800941D-04,  1.6193665D+01,
     &  7.8490384D-06,  1.3119291D-04,  1.6243546D+01,
     &  8.6997452D-06,  1.4584669D-04,  1.6293411D+01,
     &  9.6425729D-06,  1.6213443D-04,  1.6343292D+01,
     &  1.0687490D-05,  1.8023804D-04,  1.6393173D+01,
     &  1.1845541D-05,  2.0035966D-04,  1.6443039D+01,
     &  1.3128964D-05,  2.2272386D-04,  1.6492920D+01,
     &  1.4551323D-05,  2.4758000D-04,  1.6542801D+01,
     &  1.6127640D-05,  2.7520582D-04,  1.6592682D+01,
     &  1.7874583D-05,  3.0590896D-04,  1.6642563D+01,
     &  1.9810584D-05,  3.4003169D-04,  1.6692429D+01,
     &  2.1956104D-05,  3.7795468D-04,  1.6742310D+01,
     &  2.4333785D-05,  4.2010029D-04,  1.6792191D+01,
     &  2.6968744D-05,  4.6693766D-04,  1.6842072D+01,
     &  2.9888790D-05,  5.1898882D-04,  1.6891968D+01,
     &  3.3124757D-05,  5.7683280D-04,  1.6941849D+01,
     &  3.6710771D-05,  6.4111361D-04,  1.6991730D+01,
     &  4.0684681D-05,  7.1254652D-04,  1.7041611D+01,
     &  4.5088425D-05,  7.9192570D-04,  1.7091492D+01,
     &  4.9968454D-05,  8.8013406D-04,  1.7141388D+01,
     &  5.5376222D-05,  9.7815227D-04,  1.7191269D+01,
     &  6.1368788D-05,  1.0870695D-03,  1.7241150D+01,
     &  6.8009307D-05,  1.2080960D-03,  1.7291046D+01,
     &  7.5367832D-05,  1.3425760D-03,  1.7340927D+01,
     &  8.3521911D-05,  1.4920027D-03,  1.7390823D+01,
     &  9.2557486D-05,  1.6580352D-03,  1.7440704D+01,
     &  1.0256980D-04,  1.8425162D-03,  1.7490601D+01,
     &  1.1366437D-04,  2.0474929D-03,  1.7540482D+01,
     &  1.2595803D-04,  2.2752383D-03,  1.7590378D+01,
     &  1.3958036D-04,  2.5282786D-03,  1.7640274D+01,
     &  1.5467482D-04,  2.8094193D-03,  1.7690155D+01,
     &  1.7140039D-04,  3.1217760D-03,  1.7740051D+01,
     &  1.8993320D-04,  3.4688103D-03,  1.7789948D+01,
     &  2.1046840D-04,  3.8543672D-03,  1.7839844D+01,
     &  2.3322215D-04,  4.2827129D-03,  1.7889740D+01,
     &  2.5843387D-04,  4.7585964D-03,  1.7939636D+01,
     &  2.8636912D-04,  5.2872859D-03,  1.7989532D+01,
     &  3.1732186D-04,  5.8746263D-03,  1.8039429D+01,
     &  3.5161781D-04,  6.5271184D-03,  1.8089325D+01,
     &  3.8961763D-04,  7.2519854D-03,  1.8139221D+01,
     &  4.3172133D-04,  8.0572329D-03,  1.8189117D+01,
     &  4.7837174D-04,  8.9517720D-03,  1.8239014D+01,
     &  5.3005922D-04,  9.9454857D-03,  1.8288910D+01,
     &  5.8732741D-04,  1.1049356D-02,  1.8338806D+01,
     &  6.5077888D-04,  1.2275577D-02,  1.8388702D+01,
     &  7.2108023D-04,  1.3637695D-02,  1.8438614D+01,
     &  7.9897069D-04,  1.5150752D-02,  1.8488510D+01/
      data ( pp0(i),pp1(i),pp2(i),i=401,450 )/
     &  8.8526914D-04,  1.6831446D-02,  1.8538406D+01,
     &  9.8088221D-04,  1.8698335D-02,  1.8588318D+01,
     &  1.0868148D-03,  2.0772010D-02,  1.8638214D+01,
     &  1.2041803D-03,  2.3075350D-02,  1.8688126D+01,
     &  1.3342111D-03,  2.5633764D-02,  1.8738022D+01,
     &  1.4782737D-03,  2.8475456D-02,  1.8787918D+01,
     &  1.6378809D-03,  3.1631757D-02,  1.8837830D+01,
     &  1.8147093D-03,  3.5137448D-02,  1.8887726D+01,
     &  2.0106155D-03,  3.9031167D-02,  1.8937637D+01,
     &  2.2276568D-03,  4.3355804D-02,  1.8987549D+01,
     &  2.4681115D-03,  4.8158985D-02,  1.9037445D+01,
     &  2.7345039D-03,  5.3493604D-02,  1.9087357D+01,
     &  3.0296305D-03,  5.9418388D-02,  1.9137268D+01,
     &  3.3565881D-03,  6.5998495D-02,  1.9187180D+01,
     &  3.7188083D-03,  7.3306441D-02,  1.9237076D+01,
     &  4.1200891D-03,  8.1422567D-02,  1.9286987D+01,
     &  4.5646466D-03,  9.0436161D-02,  1.9336899D+01,
     &  5.0571375D-03,  1.0044622D-01,  1.9386810D+01,
     &  5.6027360D-03,  1.1156303D-01,  1.9436722D+01,
     &  6.2071569D-03,  1.2390852D-01,  1.9486633D+01,
     &  6.8767406D-03,  1.3761854D-01,  1.9536545D+01,
     &  7.6185092D-03,  1.5284365D-01,  1.9586456D+01,
     &  8.4402412D-03,  1.6975111D-01,  1.9636368D+01,
     &  9.3505494D-03,  1.8852657D-01,  1.9686279D+01,
     &  1.0358974D-02,  2.0937622D-01,  1.9736191D+01,
     &  1.1476088D-02,  2.3252887D-01,  1.9786102D+01,
     &  1.2713600D-02,  2.5823867D-01,  1.9836014D+01,
     &  1.4084477D-02,  2.8678775D-01,  1.9885925D+01,
     &  1.5603080D-02,  3.1848919D-01,  1.9935837D+01,
     &  1.7285325D-02,  3.5369086D-01,  1.9985748D+01,
     &  1.9148827D-02,  3.9277858D-01,  2.0035660D+01,
     &  2.1213111D-02,  4.3618107D-01,  2.0085587D+01,
     &  2.3499802D-02,  4.8437393D-01,  2.0135498D+01,
     &  2.6032839D-02,  5.3788537D-01,  2.0185410D+01,
     &  2.8838746D-02,  5.9730166D-01,  2.0235336D+01,
     &  3.1946912D-02,  6.6327369D-01,  2.0285248D+01,
     &  3.5389870D-02,  7.3652405D-01,  2.0335159D+01,
     &  3.9203662D-02,  8.1785470D-01,  2.0385086D+01,
     &  4.3428212D-02,  9.0815622D-01,  2.0434998D+01,
     &  4.8107732D-02,  1.0084162D+00,  2.0484924D+01,
     &  5.3291194D-02,  1.1197329D+00,  2.0534836D+01,
     &  5.9032839D-02,  1.2433243D+00,  2.0584763D+01,
     &  6.5392733D-02,  1.3805418D+00,  2.0634674D+01,
     &  7.2437406D-02,  1.5328865D+00,  2.0684601D+01,
     &  8.0240607D-02,  1.7020235D+00,  2.0734512D+01,
     &  8.8883936D-02,  1.8898029D+00,  2.0784439D+01,
     &  9.8457754D-02,  2.0982771D+00,  2.0834366D+01,
     &  1.0906219D-01,  2.3297234D+00,  2.0884277D+01,
     &  1.2080812D-01,  2.5866718D+00,  2.0934204D+01,
     &  1.3381851D-01,  2.8719282D+00,  2.0984131D+01/
      data ( pp0(i),pp1(i),pp2(i),i=451,500 )/
     &  1.4822918D-01,  3.1886091D+00,  2.1034042D+01,
     &  1.6419089D-01,  3.5401716D+00,  2.1083969D+01,
     &  1.8187040D-01,  3.9304552D+00,  2.1133896D+01,
     &  2.0145261D-01,  4.3637199D+00,  2.1183823D+01,
     &  2.2314215D-01,  4.8446932D+00,  2.1233734D+01,
     &  2.4716562D-01,  5.3786249D+00,  2.1283661D+01,
     &  2.7377409D-01,  5.9713392D+00,  2.1333588D+01,
     &  3.0324554D-01,  6.6293011D+00,  2.1383514D+01,
     &  3.3588785D-01,  7.3596869D+00,  2.1433441D+01,
     &  3.7204206D-01,  8.1704597D+00,  2.1483368D+01,
     &  4.1208577D-01,  9.0704584D+00,  2.1533295D+01,
     &  4.5643723D-01,  1.0069491D+01,  2.1583221D+01,
     &  5.0555962D-01,  1.1178450D+01,  2.1633148D+01,
     &  5.5996585D-01,  1.2409413D+01,  2.1683075D+01,
     &  6.2022406D-01,  1.3775792D+01,  2.1733002D+01,
     &  6.8696332D-01,  1.5292470D+01,  2.1782928D+01,
     &  7.6088041D-01,  1.6975952D+01,  2.1832855D+01,
     &  8.4274691D-01,  1.8844589D+01,  2.1882782D+01,
     &  9.3341732D-01,  2.0918716D+01,  2.1932709D+01,
     &  1.0338373D+00,  2.3220901D+01,  2.1982635D+01,
     &  1.1450567D+00,  2.5776199D+01,  2.2032562D+01,
     &  1.2682343D+00,  2.8612411D+01,  2.2082489D+01,
     &  1.4046555D+00,  3.1760391D+01,  2.2132431D+01,
     &  1.5557442D+00,  3.5254379D+01,  2.2182358D+01,
     &  1.7230768D+00,  3.9132385D+01,  2.2232285D+01,
     &  1.9083986D+00,  4.3436554D+01,  2.2282211D+01,
     &  2.1136417D+00,  4.8213684D+01,  2.2332138D+01,
     &  2.3409481D+00,  5.3515686D+01,  2.2382080D+01,
     &  2.5926876D+00,  5.9400208D+01,  2.2432007D+01,
     &  2.8714857D+00,  6.5931152D+01,  2.2481934D+01,
     &  3.1802483D+00,  7.3179489D+01,  2.2531876D+01,
     &  3.5221958D+00,  8.1223938D+01,  2.2581802D+01,
     &  3.9008932D+00,  9.0151871D+01,  2.2631729D+01,
     &  4.3202868D+00,  1.0006023D+02,  2.2681671D+01,
     &  4.7847490D+00,  1.1105656D+02,  2.2731598D+01,
     &  5.2991219D+00,  1.2326025D+02,  2.2781540D+01,
     &  5.8687639D+00,  1.3680373D+02,  2.2831467D+01,
     &  6.4996128D+00,  1.5183395D+02,  2.2881409D+01,
     &  7.1982412D+00,  1.6851399D+02,  2.2931335D+01,
     &  7.9719286D+00,  1.8702478D+02,  2.2981277D+01,
     &  8.8287363D+00,  2.0756709D+02,  2.3031204D+01,
     &  9.7775898D+00,  2.3036368D+02,  2.3081146D+01,
     &  1.0828371D+01,  2.5566170D+02,  2.3131073D+01,
     &  1.1992028D+01,  2.8373535D+02,  2.3181015D+01,
     &  1.3280681D+01,  3.1488892D+02,  2.3230942D+01,
     &  1.4707746D+01,  3.4946021D+02,  2.3280884D+01,
     &  1.6288071D+01,  3.8782349D+02,  2.3330826D+01,
     &  1.8038147D+01,  4.3039478D+02,  2.3380753D+01,
     &  1.9976166D+01,  4.7763477D+02,  2.3430695D+01,
     &  2.2122314D+01,  5.3005542D+02,  2.3480637D+01/
      data ( pp0(i),pp1(i),pp2(i),i=501,550 )/
     &  2.4498947D+01,  5.8822412D+02,  2.3530563D+01,
     &  2.7130768D+01,  6.5277100D+02,  2.3580505D+01,
     &  3.0045212D+01,  7.2439429D+02,  2.3630447D+01,
     &  3.3272583D+01,  8.0386963D+02,  2.3680389D+01,
     &  3.6846466D+01,  8.9205688D+02,  2.3730316D+01,
     &  4.0804077D+01,  9.8991016D+02,  2.3780258D+01,
     &  4.5186569D+01,  1.0984883D+03,  2.3830200D+01,
     &  5.0039566D+01,  1.2189658D+03,  2.3880142D+01,
     &  5.5413544D+01,  1.3526458D+03,  2.3930084D+01,
     &  6.1364410D+01,  1.5009734D+03,  2.3980011D+01,
     &  6.7954056D+01,  1.6655527D+03,  2.4029953D+01,
     &  7.5251053D+01,  1.8481626D+03,  2.4079895D+01,
     &  8.3331268D+01,  2.0507773D+03,  2.4129837D+01,
     &  9.2278732D+01,  2.2755862D+03,  2.4179779D+01,
     &  1.0218651D+02,  2.5250186D+03,  2.4229721D+01,
     &  1.1315764D+02,  2.8017693D+03,  2.4279663D+01,
     &  1.2530615D+02,  3.1088281D+03,  2.4329605D+01,
     &  1.3875838D+02,  3.4495115D+03,  2.4379547D+01,
     &  1.5365419D+02,  3.8274985D+03,  2.4429489D+01,
     &  1.7014841D+02,  4.2468672D+03,  2.4479431D+01,
     &  1.8841248D+02,  4.7121523D+03,  2.4529373D+01,
     &  2.0863628D+02,  5.2283750D+03,  2.4579315D+01,
     &  2.3102995D+02,  5.8011016D+03,  2.4629257D+01,
     &  2.5582625D+02,  6.4365156D+03,  2.4679199D+01,
     &  2.8328271D+02,  7.1414766D+03,  2.4729141D+01,
     &  3.1368481D+02,  7.9235820D+03,  2.4779083D+01,
     &  3.4734863D+02,  8.7912773D+03,  2.4829025D+01,
     &  3.8462354D+02,  9.7539180D+03,  2.4878967D+01,
     &  4.2589673D+02,  1.0821883D+04,  2.4928909D+01,
     &  4.7159741D+02,  1.2006687D+04,  2.4978867D+01,
     &  5.2219995D+02,  1.3321109D+04,  2.5028809D+01,
     &  5.7822998D+02,  1.4779312D+04,  2.5078751D+01,
     &  6.4026953D+02,  1.6397020D+04,  2.5128693D+01,
     &  7.0896289D+02,  1.8191656D+04,  2.5178635D+01,
     &  7.8502319D+02,  2.0182566D+04,  2.5228592D+01,
     &  8.6924048D+02,  2.2391199D+04,  2.5278534D+01,
     &  9.6248901D+02,  2.4841344D+04,  2.5328476D+01,
     &  1.0657368D+03,  2.7559387D+04,  2.5378418D+01,
     &  1.1800564D+03,  3.0574605D+04,  2.5428375D+01,
     &  1.3066338D+03,  3.3919469D+04,  2.5478317D+01,
     &  1.4467834D+03,  3.7629984D+04,  2.5528259D+01,
     &  1.6019597D+03,  4.1746098D+04,  2.5578201D+01,
     &  1.7737732D+03,  4.6312109D+04,  2.5628159D+01,
     &  1.9640071D+03,  5.1377168D+04,  2.5678101D+01,
     &  2.1746357D+03,  5.6995770D+04,  2.5728043D+01,
     &  2.4078445D+03,  6.3228371D+04,  2.5778000D+01,
     &  2.6660535D+03,  7.0142000D+04,  2.5827942D+01,
     &  2.9519419D+03,  7.7811062D+04,  2.5877884D+01,
     &  3.2684751D+03,  8.6318000D+04,  2.5927841D+01,
     &  3.6189375D+03,  9.5754375D+04,  2.5977783D+01/
      data ( pp0(i),pp1(i),pp2(i),i=551,600 )/
     &  4.0069644D+03,  1.0622150D+05,  2.6027740D+01,
     &  4.4365781D+03,  1.1783206D+05,  2.6077682D+01,
     &  4.9122422D+03,  1.3071081D+05,  2.6127640D+01,
     &  5.4388828D+03,  1.4499619D+05,  2.6177582D+01,
     &  6.0219648D+03,  1.6084162D+05,  2.6227539D+01,
     &  6.6675312D+03,  1.7841750D+05,  2.6277481D+01,
     &  7.3822852D+03,  1.9791262D+05,  2.6327438D+01,
     &  8.1736250D+03,  2.1953644D+05,  2.6377380D+01,
     &  9.0497695D+03,  2.4352112D+05,  2.6427338D+01,
     &  1.0019793D+04,  2.7012437D+05,  2.6477280D+01,
     &  1.1093750D+04,  2.9963187D+05,  2.6527237D+01,
     &  1.2282781D+04,  3.3236044D+05,  2.6577179D+01,
     &  1.3599207D+04,  3.6866144D+05,  2.6627136D+01,
     &  1.5056676D+04,  4.0892456D+05,  2.6677078D+01,
     &  1.6670289D+04,  4.5358194D+05,  2.6727036D+01,
     &  1.8456773D+04,  5.0311294D+05,  2.6776993D+01,
     &  2.0434641D+04,  5.5804894D+05,  2.6826935D+01,
     &  2.2624387D+04,  6.1897950D+05,  2.6876892D+01,
     &  2.5048703D+04,  6.8655831D+05,  2.6926834D+01,
     &  2.7732707D+04,  7.6151019D+05,  2.6976791D+01,
     &  3.0704207D+04,  8.4463906D+05,  2.7026749D+01,
     &  3.3993988D+04,  9.3683650D+05,  2.7076691D+01,
     &  3.7636133D+04,  1.0390912D+06,  2.7126648D+01,
     &  4.1668363D+04,  1.1524990D+06,  2.7176605D+01,
     &  4.6132457D+04,  1.2782760D+06,  2.7226547D+01,
     &  5.1074641D+04,  1.4177720D+06,  2.7276505D+01,
     &  5.6546109D+04,  1.5724800D+06,  2.7326462D+01,
     &  6.2603523D+04,  1.7440580D+06,  2.7376419D+01,
     &  6.9309562D+04,  1.9343460D+06,  2.7426361D+01,
     &  7.6733812D+04,  2.1453820D+06,  2.7476318D+01,
     &  8.4953000D+04,  2.3794270D+06,  2.7526276D+01,
     &  9.4052250D+04,  2.6389880D+06,  2.7576233D+01,
     &  1.0412587D+05,  2.9268450D+06,  2.7626175D+01,
     &  1.1527806D+05,  3.2460810D+06,  2.7676132D+01,
     &  1.2762431D+05,  3.6001140D+06,  2.7726089D+01,
     &  1.4129237D+05,  3.9927350D+06,  2.7776047D+01,
     &  1.5642375D+05,  4.4281470D+06,  2.7826004D+01,
     &  1.7317512D+05,  4.9110120D+06,  2.7875961D+01,
     &  1.9171981D+05,  5.4464960D+06,  2.7925903D+01,
     &  2.1224969D+05,  6.0403320D+06,  2.7975861D+01,
     &  2.3497731D+05,  6.6988740D+06,  2.8025818D+01,
     &  2.6013781D+05,  7.4291670D+06,  2.8075775D+01,
     &  2.8799156D+05,  8.2390260D+06,  2.8125732D+01,
     &  3.1882675D+05,  9.1371130D+06,  2.8175690D+01,
     &  3.5296237D+05,  1.0133034D+07,  2.8225647D+01,
     &  3.9075169D+05,  1.1237442D+07,  2.8275604D+01,
     &  4.3258556D+05,  1.2462145D+07,  2.8325562D+01,
     &  4.7889681D+05,  1.3820240D+07,  2.8375519D+01,
     &  5.3016437D+05,  1.5326247D+07,  2.8425476D+01,
     &  5.8691869D+05,  1.6996256D+07,  2.8475418D+01/
      data ( pp0(i),pp1(i),pp2(i),i=601,650 )/
     &  6.4974669D+05,  1.8848144D+07,  2.8525375D+01,
     &  7.1929812D+05,  2.0901664D+07,  2.8575333D+01,
     &  7.9629244D+05,  2.3178800D+07,  2.8625290D+01,
     &  8.8152569D+05,  2.5703872D+07,  2.8675247D+01,
     &  9.7587937D+05,  2.8503856D+07,  2.8725204D+01,
     &  1.0803290D+06,  3.1608672D+07,  2.8775162D+01,
     &  1.1959540D+06,  3.5051472D+07,  2.8825119D+01,
     &  1.3239520D+06,  3.8869040D+07,  2.8875092D+01,
     &  1.4656440D+06,  4.3102144D+07,  2.8925049D+01,
     &  1.6224950D+06,  4.7795984D+07,  2.8975006D+01,
     &  1.7961280D+06,  5.3000704D+07,  2.9024963D+01,
     &  1.9883370D+06,  5.8771840D+07,  2.9074921D+01,
     &  2.2011080D+06,  6.5171024D+07,  2.9124878D+01,
     &  2.4366410D+06,  7.2266560D+07,  2.9174835D+01,
     &  2.6973710D+06,  8.0134176D+07,  2.9224792D+01,
     &  2.9859910D+06,  8.8857840D+07,  2.9274750D+01,
     &  3.3054850D+06,  9.8530640D+07,  2.9324707D+01,
     &  3.6591530D+06,  1.0925578D+08,  2.9374664D+01,
     &  4.0506520D+06,  1.2114770D+08,  2.9424637D+01,
     &  4.4840250D+06,  1.3433325D+08,  2.9474594D+01,
     &  4.9637510D+06,  1.4895309D+08,  2.9524551D+01,
     &  5.4947860D+06,  1.6516314D+08,  2.9574509D+01,
     &  6.0826160D+06,  1.8313629D+08,  2.9624466D+01,
     &  6.7333140D+06,  2.0306416D+08,  2.9674423D+01,
     &  7.4536020D+06,  2.2515926D+08,  2.9724396D+01,
     &  8.2509190D+06,  2.4965717D+08,  2.9774353D+01,
     &  9.1335020D+06,  2.7681894D+08,  2.9824310D+01,
     &  1.0110467D+07,  3.0693427D+08,  2.9874268D+01,
     &  1.1191903D+07,  3.4032384D+08,  2.9924225D+01,
     &  1.2388979D+07,  3.7734400D+08,  2.9974197D+01,
     &  1.3714057D+07,  4.1838874D+08,  3.0024155D+01,
     &  1.5180822D+07,  4.6389581D+08,  3.0074112D+01,
     &  1.6804416D+07,  5.1434957D+08,  3.0124069D+01,
     &  1.8601600D+07,  5.7028813D+08,  3.0174042D+01,
     &  2.0590944D+07,  6.3230669D+08,  3.0223999D+01,
     &  2.2792976D+07,  7.0106598D+08,  3.0273956D+01,
     &  2.5230448D+07,  7.7729869D+08,  3.0323914D+01,
     &  2.7928496D+07,  8.6181606D+08,  3.0373886D+01,
     &  3.0914992D+07,  9.5551846D+08,  3.0423843D+01,
     &  3.4220752D+07,  1.0594033D+09,  3.0473801D+01,
     &  3.7879904D+07,  1.1745766D+09,  3.0523773D+01,
     &  4.1930192D+07,  1.3022636D+09,  3.0573730D+01,
     &  4.6413440D+07,  1.4438239D+09,  3.0623688D+01,
     &  5.1375920D+07,  1.6007642D+09,  3.0673660D+01,
     &  5.6868832D+07,  1.7747546D+09,  3.0723618D+01,
     &  6.2948880D+07,  1.9676465D+09,  3.0773575D+01,
     &  6.9678784D+07,  2.1814920D+09,  3.0823547D+01,
     &  7.7128000D+07,  2.4185661D+09,  3.0873505D+01,
     &  8.5373392D+07,  2.6813911D+09,  3.0923462D+01,
     &  9.4500032D+07,  2.9727624D+09,  3.0973434D+01/
      data ( pp0(i),pp1(i),pp2(i),i=651,700 )/
     &  1.0460206D+08,  3.2957791D+09,  3.1023392D+01,
     &  1.1578373D+08,  3.6538760D+09,  3.1073349D+01,
     &  1.2816035D+08,  4.0508616D+09,  3.1123322D+01,
     &  1.4185963D+08,  4.4909527D+09,  3.1173279D+01,
     &  1.5702288D+08,  4.9788396D+09,  3.1223251D+01,
     &  1.7380648D+08,  5.5196959D+09,  3.1273209D+01,
     &  1.9238357D+08,  6.1192765D+09,  3.1323181D+01,
     &  2.1294573D+08,  6.7839590D+09,  3.1373138D+01,
     &  2.3570504D+08,  7.5208008D+09,  3.1423096D+01,
     &  2.6089621D+08,  8.3376333D+09,  3.1473068D+01,
     &  2.8877901D+08,  9.2431360D+09,  3.1523026D+01,
     &  3.1964083D+08,  1.0246935D+10,  3.1572998D+01,
     &  3.5380019D+08,  1.1359691D+10,  3.1622955D+01,
     &  3.9160909D+08,  1.2593222D+10,  3.1672928D+01,
     &  4.3345766D+08,  1.3960638D+10,  3.1722885D+01,
     &  4.7977702D+08,  1.5476457D+10,  3.1772858D+01,
     &  5.3104486D+08,  1.7156780D+10,  3.1822815D+01,
     &  5.8778957D+08,  1.9019448D+10,  3.1872787D+01,
     &  6.5059635D+08,  2.1084246D+10,  3.1922745D+01,
     &  7.2011264D+08,  2.3373091D+10,  3.1972717D+01,
     &  7.9705472D+08,  2.5910288D+10,  3.2022675D+01,
     &  8.8221568D+08,  2.8722766D+10,  3.2072647D+01,
     &  9.7647360D+08,  3.1840383D+10,  3.2122604D+01,
     &  1.0808000D+09,  3.5296231D+10,  3.2172577D+01,
     &  1.1962696D+09,  3.9126979D+10,  3.2222534D+01,
     &  1.3240727D+09,  4.3373289D+10,  3.2272507D+01,
     &  1.4655263D+09,  4.8080212D+10,  3.2322464D+01,
     &  1.6220879D+09,  5.3297689D+10,  3.2372437D+01,
     &  1.7953710D+09,  5.9081081D+10,  3.2422409D+01,
     &  1.9871611D+09,  6.5491739D+10,  3.2472366D+01,
     &  2.1994342D+09,  7.2597635D+10,  3.2522339D+01,
     &  2.4343772D+09,  8.0474210D+10,  3.2572296D+01,
     &  2.6944110D+09,  8.9204916D+10,  3.2622269D+01,
     &  2.9822141D+09,  9.8882421D+10,  3.2672241D+01,
     &  3.3007519D+09,  1.0960929D+11,  3.2722198D+01,
     &  3.6533051D+09,  1.2149929D+11,  3.2772171D+01,
     &  4.0435057D+09,  1.3467845D+11,  3.2822128D+01,
     &  4.4753715D+09,  1.4928655D+11,  3.2872101D+01,
     &  4.9533542D+09,  1.6547833D+11,  3.2922073D+01,
     &  5.4823731D+09,  1.8342556D+11,  3.2972031D+01,
     &  6.0678799D+09,  2.0331836D+11,  3.3022003D+01,
     &  6.7159040D+09,  2.2536756D+11,  3.3071976D+01,
     &  7.4331136D+09,  2.4980685D+11,  3.3121933D+01,
     &  8.2269020D+09,  2.7689517D+11,  3.3171906D+01,
     &  9.1054408D+09,  3.0691957D+11,  3.3221878D+01,
     &  1.0077774D+10,  3.4019803D+11,  3.3271835D+01,
     &  1.1153916D+10,  3.7708320D+11,  3.3321808D+01,
     &  1.2344947D+10,  4.1796580D+11,  3.3371780D+01,
     &  1.3663130D+10,  4.6327877D+11,  3.3421753D+01,
     &  1.5122031D+10,  5.1350202D+11,  3.3471710D+01/
      data ( pp0(i),pp1(i),pp2(i),i=701,750 )/
     &  1.6736674D+10,  5.6916751D+11,  3.3521683D+01,
     &  1.8523681D+10,  6.3086461D+11,  3.3571655D+01,
     &  2.0501451D+10,  6.9924671D+11,  3.3621613D+01,
     &  2.2690337D+10,  7.7503765D+11,  3.3671585D+01,
     &  2.5112871D+10,  8.5903999D+11,  3.3721558D+01,
     &  2.7793990D+10,  9.5214286D+11,  3.3771530D+01,
     &  3.0761288D+10,  1.0553318D+12,  3.3821487D+01,
     &  3.4045309D+10,  1.1696981D+12,  3.3871460D+01,
     &  3.7679845D+10,  1.2964541D+12,  3.3921432D+01,
     &  4.1702306D+10,  1.4369402D+12,  3.3971405D+01,
     &  4.6154084D+10,  1.5926433D+12,  3.4021362D+01,
     &  5.1080991D+10,  1.7652095D+12,  3.4071335D+01,
     &  5.6533725D+10,  1.9564667D+12,  3.4121307D+01,
     &  6.2568395D+10,  2.1684363D+12,  3.4171280D+01,
     &  6.9247042D+10,  2.4033624D+12,  3.4221252D+01,
     &  7.6638519D+10,  2.6637280D+12,  3.4271210D+01,
     &  8.4818723D+10,  2.9522888D+12,  3.4321182D+01,
     &  9.3871931D+10,  3.2720961D+12,  3.4371155D+01,
     &  1.0389121D+11,  3.6265326D+12,  3.4421127D+01,
     &  1.1497964D+11,  4.0193449D+12,  3.4471100D+01,
     &  1.2725132D+11,  4.4546864D+12,  3.4521057D+01,
     &  1.4083247D+11,  4.9371614D+12,  3.4571030D+01,
     &  1.5586276D+11,  5.4718701D+12,  3.4621002D+01,
     &  1.7249685D+11,  6.0644645D+12,  3.4670975D+01,
     &  1.9090578D+11,  6.7212096D+12,  3.4720947D+01,
     &  2.1127889D+11,  7.4490472D+12,  3.4770920D+01,
     &  2.3382570D+11,  8.2556685D+12,  3.4820892D+01,
     &  2.5877814D+11,  9.1495984D+12,  3.4870850D+01,
     &  2.8639278D+11,  1.0140284D+13,  3.4920822D+01,
     &  3.1695359D+11,  1.1238193D+13,  3.4970795D+01,
     &  3.5077489D+11,  1.2454926D+13,  3.5020767D+01,
     &  3.8820440D+11,  1.3803338D+13,  3.5070740D+01,
     &  4.2962702D+11,  1.5297675D+13,  3.5120712D+01,
     &  4.7546866D+11,  1.6953722D+13,  3.5170685D+01,
     &  5.2620067D+11,  1.8788955D+13,  3.5220657D+01,
     &  5.8234457D+11,  2.0822790D+13,  3.5270630D+01,
     &  6.4447762D+11,  2.3076692D+13,  3.5320587D+01,
     &  7.1323858D+11,  2.5574467D+13,  3.5370560D+01,
     &  7.8933433D+11,  2.8342472D+13,  3.5420532D+01,
     &  8.7354710D+11,  3.1409951D+13,  3.5470505D+01,
     &  9.6674264D+11,  3.4809301D+13,  3.5520477D+01,
     &  1.0698787D+12,  3.8576373D+13,  3.5570450D+01,
     &  1.1840153D+12,  4.2750980D+13,  3.5620422D+01,
     &  1.3103257D+12,  4.7377163D+13,  3.5670395D+01,
     &  1.4501093D+12,  5.2503761D+13,  3.5720367D+01,
     &  1.6048015D+12,  5.8184878D+13,  3.5770340D+01,
     &  1.7759920D+12,  6.4480461D+13,  3.5820313D+01,
     &  1.9654404D+12,  7.1456965D+13,  3.5870285D+01,
     &  2.1750927D+12,  7.9188007D+13,  3.5920258D+01,
     &  2.4071058D+12,  8.7755158D+13,  3.5970230D+01/
      data ( pp0(i),pp1(i),pp2(i),i=751,800 )/
     &  2.6638612D+12,  9.7248797D+13,  3.6020203D+01,
     &  2.9479991D+12,  1.0776908D+14,  3.6070175D+01,
     &  3.2624377D+12,  1.1942702D+14,  3.6120148D+01,
     &  3.6104086D+12,  1.3234556D+14,  3.6170120D+01,
     &  3.9954867D+12,  1.4666097D+14,  3.6220093D+01,
     &  4.4216279D+12,  1.6252425D+14,  3.6270065D+01,
     &  4.8932114D+12,  1.8010268D+14,  3.6320038D+01,
     &  5.4150803D+12,  1.9958164D+14,  3.6370010D+01,
     &  5.9925982D+12,  2.2116656D+14,  3.6419983D+01,
     &  6.6316948D+12,  2.4508499D+14,  3.6469955D+01,
     &  7.3389373D+12,  2.7158917D+14,  3.6519928D+01,
     &  8.1215902D+12,  3.0095829D+14,  3.6569901D+01,
     &  8.9876920D+12,  3.3350233D+14,  3.6619873D+01,
     &  9.9461387D+12,  3.6956449D+14,  3.6669846D+01,
     &  1.1006774D+13,  4.0952433D+14,  3.6719818D+01,
     &  1.2180494D+13,  4.5380329D+14,  3.6769791D+01,
     &  1.3479349D+13,  5.0286792D+14,  3.6819763D+01,
     &  1.4916680D+13,  5.5723550D+14,  3.6869736D+01,
     &  1.6507251D+13,  6.1747886D+14,  3.6919708D+01,
     &  1.8267385D+13,  6.8423285D+14,  3.6969681D+01,
     &  2.0215170D+13,  7.5820078D+14,  3.7019653D+01,
     &  2.2370606D+13,  8.4016164D+14,  3.7069626D+01,
     &  2.4755806D+13,  9.3097953D+14,  3.7119614D+01,
     &  2.7395298D+13,  1.0316106D+15,  3.7169586D+01,
     &  3.0316144D+13,  1.1431152D+15,  3.7219559D+01,
     &  3.3548359D+13,  1.2666677D+15,  3.7269531D+01,
     &  3.7125110D+13,  1.4035693D+15,  3.7319504D+01,
     &  4.1083141D+13,  1.5552619D+15,  3.7369476D+01,
     &  4.5463051D+13,  1.7233430D+15,  3.7419449D+01,
     &  5.0309838D+13,  1.9095825D+15,  3.7469421D+01,
     &  5.5673245D+13,  2.1159414D+15,  3.7519394D+01,
     &  6.1608336D+13,  2.3445923D+15,  3.7569382D+01,
     &  6.8176012D+13,  2.5979428D+15,  3.7619354D+01,
     &  7.5443701D+13,  2.8786600D+15,  3.7669327D+01,
     &  8.3485994D+13,  3.1896985D+15,  3.7719299D+01,
     &  9.2385451D+13,  3.5343332D+15,  3.7769272D+01,
     &  1.0223339D+14,  3.9161906D+15,  3.7819244D+01,
     &  1.1313090D+14,  4.3392906D+15,  3.7869217D+01,
     &  1.2518982D+14,  4.8080827D+15,  3.7919205D+01,
     &  1.3853392D+14,  5.3275075D+15,  3.7969177D+01,
     &  1.5330012D+14,  5.9030245D+15,  3.8019150D+01,
     &  1.6963993D+14,  6.5406897D+15,  3.8069122D+01,
     &  1.8772106D+14,  7.2472204D+15,  3.8119095D+01,
     &  2.0772905D+14,  8.0300384D+15,  3.8169067D+01,
     &  2.2986917D+14,  8.8973855D+15,  3.8219055D+01,
     &  2.5436862D+14,  9.8583845D+15,  3.8269028D+01,
     &  2.8147873D+14,  1.0923146D+16,  3.8319000D+01,
     &  3.1147747D+14,  1.2102866D+16,  3.8368973D+01,
     &  3.4467327D+14,  1.3409953D+16,  3.8418945D+01,
     &  3.8140598D+14,  1.4858156D+16,  3.8468933D+01/
      data ( pp0(i),pp1(i),pp2(i),i=801,850 )/
     &  4.2205275D+14,  1.6462704D+16,  3.8518906D+01,
     &  4.6703072D+14,  1.8240468D+16,  3.8568878D+01,
     &  5.1680107D+14,  2.0210145D+16,  3.8618851D+01,
     &  5.7187436D+14,  2.2392443D+16,  3.8668823D+01,
     &  6.3281538D+14,  2.4810304D+16,  3.8718811D+01,
     &  7.0024959D+14,  2.7489152D+16,  3.8768784D+01,
     &  7.7486847D+14,  3.0457146D+16,  3.8818756D+01,
     &  8.5743761D+14,  3.3745485D+16,  3.8868729D+01,
     &  9.4880337D+14,  3.7388734D+16,  3.8918716D+01,
     &  1.0499034D+15,  4.1425183D+16,  3.8968689D+01,
     &  1.1617744D+15,  4.5897262D+16,  3.9018661D+01,
     &  1.2855634D+15,  5.0851962D+16,  3.9068634D+01,
     &  1.4225401D+15,  5.6341355D+16,  3.9118622D+01,
     &  1.5741095D+15,  6.2423123D+16,  3.9168594D+01,
     &  1.7418256D+15,  6.9161171D+16,  3.9218567D+01,
     &  1.9274082D+15,  7.6626271D+16,  3.9268539D+01,
     &  2.1327602D+15,  8.4896866D+16,  3.9318527D+01,
     &  2.3599876D+15,  9.4059921D+16,  3.9368500D+01,
     &  2.6114201D+15,  1.0421157D+17,  3.9418472D+01,
     &  2.8896357D+15,  1.1545855D+17,  3.9468445D+01,
     &  3.1974872D+15,  1.2791897D+17,  3.9518433D+01,
     &  3.5381304D+15,  1.4172361D+17,  3.9568405D+01,
     &  3.9150584D+15,  1.5701761D+17,  3.9618378D+01,
     &  4.3321346D+15,  1.7396150D+17,  3.9668365D+01,
     &  4.7936345D+15,  1.9273319D+17,  3.9718338D+01,
     &  5.3042889D+15,  2.1352983D+17,  3.9768311D+01,
     &  5.8693348D+15,  2.3656982D+17,  3.9818283D+01,
     &  6.4945661D+15,  2.6209498D+17,  3.9868271D+01,
     &  7.1863865D+15,  2.9037339D+17,  3.9918243D+01,
     &  7.9518914D+15,  3.2170192D+17,  3.9968216D+01,
     &  8.7989234D+15,  3.5640937D+17,  4.0018204D+01,
     &  9.7361669D+15,  3.9486012D+17,  4.0068176D+01,
     &  1.0773230D+16,  4.3745781D+17,  4.0118149D+01,
     &  1.1920738D+16,  4.8464947D+17,  4.0168137D+01,
     &  1.3190454D+16,  5.3693042D+17,  4.0218109D+01,
     &  1.4595390D+16,  5.9484940D+17,  4.0268082D+01,
     &  1.6149940D+16,  6.5901408D+17,  4.0318069D+01,
     &  1.7870045D+16,  7.3009792D+17,  4.0368042D+01,
     &  1.9773321D+16,  8.0884673D+17,  4.0418015D+01,
     &  2.1879281D+16,  8.9608672D+17,  4.0468002D+01,
     &  2.4209498D+16,  9.9273331D+17,  4.0517975D+01,
     &  2.6787857D+16,  1.0998005D+18,  4.0567947D+01,
     &  2.9640768D+16,  1.2184106D+18,  4.0617935D+01,
     &  3.2797470D+16,  1.3498100D+18,  4.0667908D+01,
     &  3.6290300D+16,  1.4953743D+18,  4.0717896D+01,
     &  4.0155054D+16,  1.6566320D+18,  4.0767868D+01,
     &  4.4431321D+16,  1.8352751D+18,  4.0817841D+01,
     &  4.9162910D+16,  2.0331751D+18,  4.0867828D+01,
     &  5.4398303D+16,  2.2524089D+18,  4.0917801D+01,
     &  6.0191132D+16,  2.4952757D+18,  4.0967773D+01/
      data ( pp0(i),pp1(i),pp2(i),i=851,900 )/
     &  6.6600735D+16,  2.7643207D+18,  4.1017761D+01,
     &  7.3692774D+16,  3.0623664D+18,  4.1067734D+01,
     &  8.1539851D+16,  3.3925365D+18,  4.1117722D+01,
     &  9.0222488D+16,  3.7582946D+18,  4.1167694D+01,
     &  9.9829471D+16,  4.1634734D+18,  4.1217667D+01,
     &  1.1045927D+17,  4.6123204D+18,  4.1267654D+01,
     &  1.2222082D+17,  5.1095427D+18,  4.1317627D+01,
     &  1.3523450D+17,  5.6603496D+18,  4.1367615D+01,
     &  1.4963357D+17,  6.2705170D+18,  4.1417587D+01,
     &  1.6556563D+17,  6.9464385D+18,  4.1467560D+01,
     &  1.8319369D+17,  7.6951982D+18,  4.1517548D+01,
     &  2.0269840D+17,  8.5246434D+18,  4.1567520D+01,
     &  2.2427955D+17,  9.4434646D+18,  4.1617508D+01,
     &  2.4815799D+17,  1.0461292D+19,  4.1667480D+01,
     &  2.7457836D+17,  1.1588791D+19,  4.1717468D+01,
     &  3.0381121D+17,  1.2837772D+19,  4.1767441D+01,
     &  3.3615582D+17,  1.4221324D+19,  4.1817413D+01,
     &  3.7194341D+17,  1.5753940D+19,  4.1867401D+01,
     &  4.1154047D+17,  1.7451676D+19,  4.1917374D+01,
     &  4.5535243D+17,  1.9332317D+19,  4.1967361D+01,
     &  5.0382777D+17,  2.1415549D+19,  4.2017334D+01,
     &  5.5746298D+17,  2.3723221D+19,  4.2067322D+01,
     &  6.1680706D+17,  2.6279489D+19,  4.2117294D+01,
     &  6.8246762D+17,  2.9111127D+19,  4.2167282D+01,
     &  7.5511696D+17,  3.2247779D+19,  4.2217255D+01,
     &  8.3549868D+17,  3.5722306D+19,  4.2267242D+01,
     &  9.2443584D+17,  3.9571089D+19,  4.2317215D+01,
     &  1.0228387D+18,  4.3834432D+19,  4.2367203D+01,
     &  1.1317147D+18,  4.8556967D+19,  4.2417175D+01,
     &  1.2521777D+18,  5.3788144D+19,  4.2467148D+01,
     &  1.3854627D+18,  5.9582729D+19,  4.2517136D+01,
     &  1.5329314D+18,  6.6001378D+19,  4.2567108D+01,
     &  1.6960956D+18,  7.3111313D+19,  4.2617096D+01,
     &  1.8766245D+18,  8.0986930D+19,  4.2667068D+01,
     &  2.0763650D+18,  8.9710684D+19,  4.2717056D+01,
     &  2.2973625D+18,  9.9373861D+19,  4.2767029D+01,
     &  2.5418796D+18,  1.1007763D+20,  4.2817017D+01,
     &  2.8124155D+18,  1.2193401D+20,  4.2866989D+01,
     &  3.1117433D+18,  1.3506707D+20,  4.2916977D+01,
     &  3.4429228D+18,  1.4961426D+20,  4.2966949D+01,
     &  3.8093449D+18,  1.6572777D+20,  4.3016937D+01,
     &  4.2147590D+18,  1.8357626D+20,  4.3066910D+01,
     &  4.6633136D+18,  2.0334644D+20,  4.3116898D+01,
     &  5.1595991D+18,  2.2524520D+20,  4.3166870D+01,
     &  5.7086930D+18,  2.4950162D+20,  4.3216858D+01,
     &  6.3162149D+18,  2.7636946D+20,  4.3266830D+01,
     &  6.9883805D+18,  3.0612965D+20,  4.3316818D+01,
     &  7.7320670D+18,  3.3909375D+20,  4.3366806D+01,
     &  8.5548855D+18,  3.7560668D+20,  4.3416779D+01,
     &  9.4652514D+18,  4.1604985D+20,  4.3466766D+01/
      data ( pp0(i),pp1(i),pp2(i),i=901,950 )/
     &  1.0472481D+19,  4.6084688D+20,  4.3516739D+01,
     &  1.1586880D+19,  5.1046585D+20,  4.3566727D+01,
     &  1.2819847D+19,  5.6542581D+20,  4.3616699D+01,
     &  1.4183998D+19,  6.2630124D+20,  4.3666687D+01,
     &  1.5693288D+19,  6.9372942D+20,  4.3716660D+01,
     &  1.7363156D+19,  7.6841458D+20,  4.3766647D+01,
     &  1.9210667D+19,  8.5113839D+20,  4.3816620D+01,
     &  2.1254756D+19,  9.4276525D+20,  4.3866608D+01,
     &  2.3516320D+19,  1.0442536D+21,  4.3916595D+01,
     &  2.6018474D+19,  1.1566637D+21,  4.3966568D+01,
     &  2.8786833D+19,  1.2811716D+21,  4.4016556D+01,
     &  3.1849703D+19,  1.4190781D+21,  4.4066528D+01,
     &  3.5238415D+19,  1.5718255D+21,  4.4116516D+01,
     &  3.8987627D+19,  1.7410100D+21,  4.4166489D+01,
     &  4.3135688D+19,  1.9283994D+21,  4.4216476D+01,
     &  4.7725015D+19,  2.1359534D+21,  4.4266464D+01,
     &  5.2802559D+19,  2.3658405D+21,  4.4316437D+01,
     &  5.8420219D+19,  2.6204631D+21,  4.4366425D+01,
     &  6.4635486D+19,  2.9024824D+21,  4.4416397D+01,
     &  7.1511884D+19,  3.2148453D+21,  4.4466385D+01,
     &  7.9119766D+19,  3.5608155D+21,  4.4516373D+01,
     &  8.7536888D+19,  3.9440082D+21,  4.4566345D+01,
     &  9.6849365D+19,  4.3684272D+21,  4.4616333D+01,
     &  1.0715240D+20,  4.8385053D+21,  4.4666306D+01,
     &  1.1855137D+20,  5.3591575D+21,  4.4716293D+01,
     &  1.3116278D+20,  5.9358164D+21,  4.4766281D+01,
     &  1.4511564D+20,  6.5745124D+21,  4.4816254D+01,
     &  1.6055259D+20,  7.2819153D+21,  4.4866241D+01,
     &  1.7763143D+20,  8.0654155D+21,  4.4916214D+01,
     &  1.9652685D+20,  8.9331916D+21,  4.4966202D+01,
     &  2.1743196D+20,  9.8943093D+21,  4.5016190D+01,
     &  2.4056055D+20,  1.0958812D+22,  4.5066162D+01,
     &  2.6614903D+20,  1.2137809D+22,  4.5116150D+01,
     &  2.9445901D+20,  1.3443614D+22,  4.5166122D+01,
     &  3.2577970D+20,  1.4889869D+22,  4.5216110D+01,
     &  3.6043180D+20,  1.6491668D+22,  4.5266098D+01,
     &  3.9876926D+20,  1.8265740D+22,  4.5316071D+01,
     &  4.4118388D+20,  2.0230611D+22,  4.5366058D+01,
     &  4.8810913D+20,  2.2406791D+22,  4.5416046D+01,
     &  5.4002522D+20,  2.4817005D+22,  4.5466019D+01,
     &  5.9746216D+20,  2.7486410D+22,  4.5516006D+01,
     &  6.6100739D+20,  3.0442870D+22,  4.5565994D+01,
     &  7.3131055D+20,  3.3717257D+22,  4.5615967D+01,
     &  8.0908969D+20,  3.7343740D+22,  4.5665955D+01,
     &  8.9514053D+20,  4.1360181D+22,  4.5715927D+01,
     &  9.9034184D+20,  4.5808494D+22,  4.5765915D+01,
     &  1.0956670D+21,  5.0735108D+22,  4.5815903D+01,
     &  1.2121923D+21,  5.6191440D+22,  4.5865875D+01,
     &  1.3411086D+21,  6.2234437D+22,  4.5915863D+01,
     &  1.4837334D+21,  6.8927151D+22,  4.5965851D+01/
      data ( pp0(i),pp1(i),pp2(i),i=951,1001 )/
     &  1.6415243D+21,  7.6339400D+22,  4.6015839D+01,
     &  1.8160940D+21,  8.4548562D+22,  4.6065811D+01,
     &  2.0092261D+21,  9.3640357D+22,  4.6115799D+01,
     &  2.2228943D+21,  1.0370954D+23,  4.6165787D+01,
     &  2.4592818D+21,  1.1486117D+23,  4.6215759D+01,
     &  2.7208044D+21,  1.2721163D+23,  4.6265747D+01,
     &  3.0101342D+21,  1.4088982D+23,  4.6315735D+01,
     &  3.3302276D+21,  1.5603827D+23,  4.6365707D+01,
     &  3.6843549D+21,  1.7281522D+23,  4.6415695D+01,
     &  4.0761348D+21,  1.9139549D+23,  4.6465683D+01,
     &  4.5095700D+21,  2.1197298D+23,  4.6515656D+01,
     &  4.9890877D+21,  2.3476220D+23,  4.6565643D+01,
     &  5.5195892D+21,  2.6000102D+23,  4.6615631D+01,
     &  6.1064938D+21,  2.8795252D+23,  4.6665604D+01,
     &  6.7557958D+21,  3.1890825D+23,  4.6715591D+01,
     &  7.4741334D+21,  3.5319102D+23,  4.6765579D+01,
     &  8.2688386D+21,  3.9115838D+23,  4.6815552D+01,
     &  9.1480358D+21,  4.3320615D+23,  4.6865540D+01,
     &  1.0120701D+22,  4.7977279D+23,  4.6915527D+01,
     &  1.1196773D+22,  5.3134391D+23,  4.6965515D+01,
     &  1.2387241D+22,  5.8845712D+23,  4.7015488D+01,
     &  1.3704274D+22,  6.5170790D+23,  4.7065475D+01,
     &  1.5161314D+22,  7.2175567D+23,  4.7115463D+01,
     &  1.6773256D+22,  7.9933064D+23,  4.7165436D+01,
     &  1.8556551D+22,  8.8524159D+23,  4.7215424D+01,
     &  2.0529425D+22,  9.8038392D+23,  4.7265411D+01,
     &  2.2712022D+22,  1.0857493D+24,  4.7315399D+01,
     &  2.5126636D+22,  1.2024364D+24,  4.7365372D+01,
     &  2.7797928D+22,  1.3316601D+24,  4.7415359D+01,
     &  3.0753181D+22,  1.4747688D+24,  4.7465347D+01,
     &  3.4022578D+22,  1.6332540D+24,  4.7515335D+01,
     &  3.7639505D+22,  1.8087667D+24,  4.7565308D+01,
     &  4.1640904D+22,  2.0031354D+24,  4.7615295D+01,
     &  4.6067636D+22,  2.2183870D+24,  4.7665283D+01,
     &  5.0964904D+22,  2.4567627D+24,  4.7715271D+01,
     &  5.6382730D+22,  2.7207483D+24,  4.7765244D+01,
     &  6.2376427D+22,  3.0130923D+24,  4.7815231D+01,
     &  6.9007207D+22,  3.3368419D+24,  4.7865219D+01,
     &  7.6342715D+22,  3.6953705D+24,  4.7915207D+01,
     &  8.4457985D+22,  4.0924113D+24,  4.7965179D+01,
     &  9.3435785D+22,  4.5321033D+24,  4.8015167D+01,
     &  1.0336784D+23,  5.0190247D+24,  4.8065155D+01,
     &  1.1435555D+23,  5.5582484D+24,  4.8115143D+01,
     &  1.2651101D+23,  6.1553914D+24,  4.8165115D+01,
     &  1.3995847D+23,  6.8166726D+24,  4.8215103D+01,
     &  1.5483512D+23,  7.5489818D+24,  4.8265091D+01,
     &  1.7129293D+23,  8.3599434D+24,  4.8315079D+01,
     &  1.8949987D+23,  9.2580058D+24,  4.8365051D+01,
     &  2.0964184D+23,  1.0252519D+25,  4.8415039D+01,
     &  2.3192443D+23,  1.1353843D+25,  4.8465027D+01,
     &  2.5657511D+23,  1.2573445D+25,  4.8515015D+01/



      end

