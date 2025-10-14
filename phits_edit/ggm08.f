************************************************************************
*                                                                      *
      subroutine cp_coll(eein,wgti,ireg,imat,ffac,rhi,icmm,elrt,itypin)
*                                                                      *
*        calculate the collision of p,d,a with a nucleus.           *
*                                                                      *
*     This routine calculates everything about a collision             *
*     involving an incident charged particle in the                    *
*     tabular energy range                                             *
*     use global_data                                                  *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for  REDUCTION_COUNTER
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2
      use GGMARRAYMOD !2020ASTOM
      use GGMBANKMOD !FURUTA

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clusttp/ nixcos, ixcoss(nnn)
!$OMP THREADPRIVATE(/clusttp/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

C for USE_MOD_COUNTER

      common /rcomon/ rcasc

*-----------------------------------------------------------------------

      dimension usave(3)
      character ht*10

*-----------------------------------------------------------------------

           nclst  = -1
           nclsts = 0
           nixcos = 0  ! for T-Point

*-----------------------------------------------------------------------

! T.Sato 2021/08/14, add deuteron and alpha
           if(itypin.eq.1) then ! proton
            iptin = 9
            kfin = 2212
            izin = 1
            inin = 0
           elseif(itypin.eq.15) then ! deuteron
            iptin = 31
            kfin = 1000002
            izin = 1
            inin = 1
           elseif(itypin.eq.18) then ! alpha
            iptin = 34
            kfin = 2000004
            izin = 2
            inin = 2
           endif

            nter = 0
            ipt  = iptin
            iels = 0
            ifis = 0

            erg  = eein
            eg0  = eein
            icl  = ireg
            wgt  = wgti
            wg0  = wgti
            wgf  = wgti
            tme  = 0.0

            usave(1) = uuu
            usave(2) = vvv
            usave(3) = www

            uuu  = 0.0
            vvv  = 0.0
            www  = 1.0

            do i = 0, 20

               numpal(i) = 0
               rumpal(i) = 0.0d0

            end do

            kcoll = 0

*-----------------------------------------------------------------------
*     save the incoming direction.
*-----------------------------------------------------------------------

            uold(1) = 0.0
            uold(2) = 0.0
            uold(3) = 1.0

            ntyn = 0
            ns   = 0
            mk   = imat
            m    = jmd(1+mk)
            m1   = m

*-----------------------------------------------------------------------
*     sample the nuclide, using the cumulative total cross section.
*-----------------------------------------------------------------------

         if(npq(mk).eq.1) goto 20

*-----------------------------------------------------------------------

         if( icmm .ne. 0 ) then

               m = m1 + icmm - 1

         else

            c = rang()*totm

            do 10 m = m1, jmd(1+mk+1) - 2
              if(lme(ipt,m).ne.0) c = c-rtc(5,lme(ipt,m))*fme(m) ! T.Sato 2023/08/17, avoid access violation by checking lme.ne.0
   10          if(c.lt.0.) goto 20

         end if

*-----------------------------------------------------------------------

   20    iex = lme(ipt,m)
         iexp = m
         mtp = 2
         mpan = ipan(icl)+m-m1

*-----------------------------------------------------------------------
*           mother nucleus
*-----------------------------------------------------------------------

            mathz = iza(iexp) / 1000
            mathn = iza(iexp) - 1000 * mathz - mathz

*-----------------------------------------------------------------------
*     generate photons if appropriate
*-----------------------------------------------------------------------

         if( kpt(2) .ne. 0 .and. gwt(icl) .ne. -1.e6 .and.
     &     ( ( jxs(12,iex) .ne. 0 .and. npikmt .eq. 0 ) .or.
     &         totgp1 .ne. 0. ) ) then

                  call cp_gam(ireg,imat,ffac,iptin)

         end if

*-----------------------------------------------------------------------
*     generate other secondary particles if appropriate
*-----------------------------------------------------------------------

         if(nxs(7,iex).ne.0) then

                  call acecp(ireg,imat)

         end if

*-----------------------------------------------------------------------
*     for now, only allow implicit capture by weight reduction
*-----------------------------------------------------------------------

      if(rtc(3,iex).gt.0.) then

        t1 = wgt*rtc(3,iex)/rtc(5,iex)
        wgt = wgt-t1

                     wga = wgf - wgt
                     aevts(iaevt+59,ireg) =
     &               aevts(iaevt+59,ireg) + wga
                     bevts(ibevt+59,imat) =
     &               bevts(ibevt+59,imat) + wga
      end if

                     wgf = wgt

*-----------------------------------------------------------------------
*     decide on elastic vs nonelastic reaction
*-----------------------------------------------------------------------

         if(nxs(5,iex).eq.0) goto 60
         ic = ktc(1,iex)+jxs(1,iex)-1

*-----------------------------------------------------------------------
*     el is the elastic cross section
*-----------------------------------------------------------------------

         n = nxs(3,iex)
         el = xss(ic+3*n)+rtc(1,iex)*(xss(ic+3*n+1)-xss(ic+3*n))

         elrt = el / rtc(5,iex)

         r = rang()*rtc(5,iex)-el
         if(r.ge.0.) goto 40

*-----------------------------------------------------------------------
*     elastic scattering sampled
*-----------------------------------------------------------------------

   60    ixre = 0
         mtp = 2
         ntyn = -99
         cmult = 1.

               kcoll = 3

*-----------------------------------------------------------------------
*     sample the center-of-mass scattering cosine (c)
*-----------------------------------------------------------------------

         c = acecos(jxs(9,iex),nint(xss(jxs(8,iex))))

         ac = -c

*        have not yet derived the kinematic formulae.
*        for initial implementation, simply use the neutron elastic
*        scattering formulae, with awr being the ratio to the ipt mass.

         yy = awn(iex)*gpt(1)/gpt(ipt)
         t1 = 1.+yy*(yy+2.*c)
         colout(1,1) = erg*t1/(1.+yy)**2
         colout(2,1) = (1.+yy*c)/sqrt(t1)

         ns = 1

*-----------------------------------------------------------------------
*        booking of proton and target proton if exists
*           need to consider light recoil particles here.  for example,
*           code needs to account for p+p elastic scattering.
*        if target is proton ns = 2
*-----------------------------------------------------------------------
                     aevts(iaevt+60,ireg) =
     &               aevts(iaevt+60,ireg) + wgt
                     bevts(ibevt+60,imat) =
     &               bevts(ibevt+60,imat) + wgt
               if( iza(iexp) .eq. 1001 ) then

                     yy = awn(iex)*gpt(1)/gpt(ipt)
                     t1 = 1.+yy*(yy+2.*ac)
                     colout(1,2) = eein*t1/(1.+yy)**2
                     colout(2,2) = (1.+yy*ac)/sqrt(t1)

                     ns = ns + 1

               end if

         goto 500

*-----------------------------------------------------------------------
*     inelastic scattering sampled
*-----------------------------------------------------------------------

   40 continue
         tt = 0.
         kx = 0
         if(nxs(5,iex).eq.1) then
            ixre = 1
            goto 150
         end if

*-----------------------------------------------------------------------
*     sample the inelastic reaction.
*-----------------------------------------------------------------------

      do 120 ixre = 1, nxs(5,iex)

*-----------------------------------------------------------------------
*     calculate the cross section for reaction whose index is ixre.
*     data:  ie,ne,(xs(i),i=1,ne)   xs(1) corresponds to es(ie).
*-----------------------------------------------------------------------

         is = jxs(7,iex)+nint(xss(jxs(6,iex)+ixre-1))
         ic = ktc(2,iex)+1-nint(xss(is-1))
         if(ic.lt.1.or.ic.gt.nint(xss(is)).or.ic.eq.nint(xss(is)).and.
     &    rtc(2,iex).ne.0.) goto 120
         t = xss(ic+is)
         if(rtc(2,iex).ne.0.) t = t+rtc(2,iex)*
     &    (xss(ic+is+1)-xss(ic+is))
         if(t.eq.0.) goto 120
         if(t.gt.tt) kx = ixre
         tt = max(t,tt)
         r = r-t
         if(r.lt.0.) goto 150
  120    continue

*-----------------------------------------------------------------------
*     handle failure of the reaction cross sections to add up.
*-----------------------------------------------------------------------

         if(kx.eq.0) goto 60

         call zaid(2,ht,ixl(1,iex))

         if( rcasc .le. 50. ) then

            ErrCha = ''
            ErrID = 'L:325/R:cp_coll/F:ggm08.f' !W04_006_001
            call ErrWrite(ErrID,ErrCha)

            write(*,*) 'Warning : no reaction mt found in co_coll.'//
     &                 ' collision resampled. zaid = '//ht
            write(*,'(''  ncasc = '',i3,'' erg = '',e13.5)')
     &      int(rcasc), erg
         end if

         ixre = kx

*-----------------------------------------------------------------------
*     get the reaction type number.
*-----------------------------------------------------------------------

  150    ntyn = nint(xss(jxs(5,iex)+ixre-1))
         mtp = nint(xss(jxs(3,iex)+ixre-1))
         q = xss(jxs(4,iex)+ixre-1)
         ia = jxs(9,iex)
         ka = nint(xss(jxs(8,iex)+ixre))
         id = jxs(11,iex)
         kd = nint(xss(jxs(10,iex)+ixre-1))

         if(abs(ntyn).ge.100) goto 190
         cmult = abs(ntyn)
         ns = cmult

      do 170 i = 1, ns
         call xstcas(i,iptin,q,ia,ka,id,kd)

  170    if(kdb.ne.0) return

         nter = 17

      goto 500

*-----------------------------------------------------------------------
*     high energy reaction other than fission with energy-dependent
*     multiplicity.  kalbach-87 (law 44) endf/b-vi.

*     the location of the multiplicity table relative to dlw is
*     abs(ntyn)-100.  acefcn gets non-fission reaction multiplicity.
*-----------------------------------------------------------------------

  190    l = jxs(11,iex)+abs(ntyn)-101

         cmult = acefcn(l,erg,ln)
         cmult = aint(cmult+rang())

         if(cmult.eq.0.) then

            nter = 16

         else

            do 200 i = 1, int(cmult)
               call xstcas(i,iptin,q,ia,ka,id,kd)

  200       if(kdb.ne.0) return

            nter = 17

         end if

         ns = cmult

         goto 500

*-----------------------------------------------------------------------
*     bank all particles same as primary from the collision.
*-----------------------------------------------------------------------

  500 continue

      if( ns .gt. 0 ) then

         if( kcoll .ne. 3 ) kcoll = 4

         do k = 1, ns

                  call dtcos(colout(2,k),uold,uuu,0,irdm)
                  erg = colout(1,k)

*-----------------------------------------------------------------------
*           bank all particles same as primary
*-----------------------------------------------------------------------

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 1

! T.Sato 2021/08/14, deuteron and alpha
                     jclusts(0,nclsts) = 0
                     jclusts(1,nclsts) = izin
                     jclusts(2,nclsts) = inin
                     jclusts(3,nclsts) = itypin
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = izin
                     jclusts(6,nclsts) = izin+inin
                     jclusts(7,nclsts) = kfin
                     jclusts(8,nclsts) = 0
                     rms = rmtyp(itypin,kfin)


                     rmg = rms / 1000.
                     pr  = sqrt( erg * ( erg + 2.0 * rms ) ) / 1000.
                     ett = sqrt( pr**2 + rmg**2 )

                     pxl = pr * uuu
                     pyl = pr * vvv
                     pzl = pr * www

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = pxl
                     qclusts(2,nclsts)  = pyl
                     qclusts(3,nclsts)  = pzl
                     qclusts(4,nclsts)  = ett
                     qclusts(5,nclsts)  = rmg
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = erg
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = tme
                     qclusts(10,nclsts) = 0.0d0
                     qclusts(11,nclsts) = 0.0d0
                     qclusts(12,nclsts) = 0.0d0
                     aevts(iaevt+53,ireg) =
     &               aevts(iaevt+53,ireg) + wgt
                     bevts(ibevt+53,imat) =
     &               bevts(ibevt+53,imat) + wgt

                     aevts(iaevt+61,ireg) =
     &               aevts(iaevt+61,ireg) + wgt
                     bevts(ibevt+61,imat) =
     &               bevts(ibevt+61,imat) + wgt
c 2024/4/10 Do not count because they are primary
c                     numpal(1) = numpal(1) + 1
c                     rumpal(1) = rumpal(1) + wgt

         end do

      end if

*-----------------------------------------------------------------------

  900 continue

             uuu = usave(1)
             vvv = usave(2)
             www = usave(3)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine cp_gam(ireg,imat,ffac,iptin)
*                                                                      *
*        generate and bank photons from a proton collision.            *
*        last modified by K.Niita on 2009/10/06                        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for  REDUCTION_COUNTER
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2
      use GGMARRAYMOD !2020ASTOM
      use GGMBANKMOD !FURUTA

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      include 'err.inc'

      common /rcomon/ rcasc

*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)


*-----------------------------------------------------------------------

      character ht*10
      dimension yd(2)
      data f/1.00000000001d0/

      data iwarning/0/ ! T.Sato 2022/12/07

*-----------------------------------------------------------------------
*     initialize.
*-----------------------------------------------------------------------

         kp = mpan
         kx = -100 ! T.Sato 2021/05/12
*-----------------------------------------------------------------------
*     calculate the probability of production of a photon.
*-----------------------------------------------------------------------

         l = jxs(12,iex)+ktc(1,iex)
         gp = (xss(l-1)+rtc(1,iex)*(xss(l)-xss(l-1)))
     &      / rtc(5,iex)

         if(gp.eq.0.) return

*-----------------------------------------------------------------------
*     Decide how many photons to make.
*-----------------------------------------------------------------------
!
!     Next setmat commented out -- don't know if need to keep
!         from nearest two integers to gp.  weight of each particle is
!         simply the incident particle weight.


*-----------------------------------------------------------------------
*     same as prdgam
*-----------------------------------------------------------------------

         gp = wgt*gp
         ni = 1
         t1 = max(gwt(icl),-gwt(icl)*wg0) * ffac

*-----------------------------------------------------------------------
*        determine split number, ni, or roulette.
*-----------------------------------------------------------------------

  130    if(gp.lt.t1) goto 140

         if(t1.ne.0.) ni = min(10,int(gp/(5.*t1)+1.))
         goto 150

  140    if(gp.le.t1*rang()) return
         gp = t1

*-----------------------------------------------------------------------
*     generate photons and write them to the bank.
*-----------------------------------------------------------------------

  150 continue

         ipt = 2
         gp = gp/ni
         es = erg
         st = totm
         vel = slite

*-----------------------------------------------------------------------

      do 450 k = 1, ni

         wgt = gp

*-----------------------------------------------------------------------
*     sample from full photon-production data, if any.
*-----------------------------------------------------------------------

         ipsc = 102

*-----------------------------------------------------------------------
*     sample the reaction that produced the photon.
*-----------------------------------------------------------------------

         ixre = 1
         if(nxs(6,iex).eq.1) goto 350
         r = 0.
         l = jxs(12,iex)+ktc(1,iex)
         r = rang()*(xss(l-1)+rtc(1,iex)*(xss(l)-xss(l-1)))
         s = 0.
         tt = 0.
         md = 0
      do 290 ixre = 1, nxs(6,iex)
         l = jxs(15,iex)+nint(xss(jxs(14,iex)+ixre-1))+1
         if(nint(xss(l-2)).eq.13) goto 250

*-----------------------------------------------------------------------
*     mf=12 -- partial cross section = yield * neutron cross section.
*-----------------------------------------------------------------------

      do 230 jy = 1, 2
  230    yd(jy) = acefcn(l,xss(jxs(1,iex)
     &            + ktc(1,iex)+jy-2)*f,ln)

         if(yd(1).eq.0..and.yd(2).eq.0.) goto 290

         ix = nint(xss(l-1))

         if(ix.eq.md) goto 240
         md = ix
         x1 = 0.
         x2 = 0.
         is = jxs(21,iex)+1
         if(ix.gt.0)
     &       is = jxs(7,iex)+nint(xss(jxs(6,iex)+ix-1))
         ic = min(ktc(1,iex)-nint(xss(is-1))+1,nint(xss(is)))
         if(ic.gt.0) x1 = xss(is+ic)
         ic = min(ic+1,nint(xss(is)))
         if(ic.gt.0) x2 = xss(is+ic)
  240    t = x1*yd(1)+(x2*yd(2)-x1*yd(1))*rtc(1,iex)
         goto 270

*-----------------------------------------------------------------------
*     mf=13 -- partial cross section is direct from table.
*-----------------------------------------------------------------------

  250    ic = ktc(1,iex)-nint(xss(l-1))+1
         if(ic.lt.1) goto 290
         if(rtc(1,iex).ne.0.) goto 260
         if(ic.gt.nint(xss(l))) goto 290
         t = xss(ic+l)
         goto 270

  260    if(ic.ge.nint(xss(l))) goto 290
         t = xss(ic+l)+rtc(1,iex)*(xss(ic+l+1)-xss(ic+l))
  270    if(t.eq.0.) goto 290
         if(t.gt.tt) kx = ixre
         tt = max(t,tt)
         s = s + t
         if(s.gt.r) goto 350
  290    continue

*-----------------------------------------------------------------------
*     handle failure of the reaction cross sections to add up.
*-----------------------------------------------------------------------

      call zaid(2,ht,ixl(1,iex))

         if( iwarning.eq.0 ) then
            write(ErrCha,*) 'Warning : in cp_gam, '//
     &      'no photon-production xsec was found in nuclear data '//
     &      'library for '//ht
            ErrID = 'L:669/R:cp_gam/F:ggm08.f' !W04_007_001
            call ErrWrite(ErrID,ErrCha)
            iwarning=1
         end if

         if(kx.eq.-100) kx=1 ! T.Sato 2021/05/12, to avoid uninitilization occurred in 7Li data

         ixre = kx

*-----------------------------------------------------------------------
*     sample the energy and direction of the photon.
*-----------------------------------------------------------------------

  350    erg = es
         mtp = nint(xss(jxs(13,iex)+ixre-1))
         ia = jxs(17,iex)

         if(jxs(16,iex).ne.0) then
            ka = nint(xss(jxs(16,iex)+ixre-1))
         else
            ka = 0
         end if

         id = jxs(19,iex)
         kd = nint(xss(jxs(18,iex)+ixre-1))
         call xstcas(1,iptin,zero,ia,ka,id,kd)

*-----------------------------------------------------------------------

         if(colout(1,1).eq.-huge) goto 450

         call dtcos(colout(2,1),uold,uuu,lev,irt)

         erg = colout(1,1)

*-----------------------------------------------------------------------
*     booking
*-----------------------------------------------------------------------

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 4

                     jclusts(0,nclsts) = 0
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(3,nclsts) = 14
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = 0
                     jclusts(6,nclsts) = 0
                     jclusts(7,nclsts) = 22
                     jclusts(8,nclsts) = 0

                     pr   = erg / 1000.0

                     pxl = pr * uuu
                     pyl = pr * vvv
                     pzl = pr * www

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = pxl
                     qclusts(2,nclsts)  = pyl
                     qclusts(3,nclsts)  = pzl
                     qclusts(4,nclsts)  = pr
                     qclusts(5,nclsts)  = 0.0
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = erg
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = 0.0
                     qclusts(11,nclsts) = 0.0
                     qclusts(12,nclsts) = 0.0
                     aevts(iaevt+58,ireg) =
     &               aevts(iaevt+58,ireg) + wgt
                     bevts(ibevt+58,imat) =
     &               bevts(ibevt+58,imat) + wgt
                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

  450 continue

*-----------------------------------------------------------------------

         wgt = wg0
         erg = eg0
         ipt = iptin
         totm = st

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine acecp(ireg,imat)
*                                                                      *
*        generate and bank other particles from nuclear reactions      *
*        (very similar to prdgam for neutron-induced photons).         *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for  REDUCTION_COUNTER
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2
      use GGMARRAYMOD !2020ASTOM
      use GGMBANKMOD !FURUTA

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      include 'err.inc'

      common /rcomon/ rcasc

*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
C MATSUDA 2022.12.12  for T-Point
      common /clusttp/ nixcos, ixcoss(nnn)
!$OMP THREADPRIVATE(/clusttp/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)

*-----------------------------------------------------------------------

      character ht*10
      dimension yd(2)

      data f/1.00000000001d0/

*-----------------------------------------------------------------------

         es = erg
         in = ipt

*-----------------------------------------------------------------------
*     do loop over number of other secondary particles
*-----------------------------------------------------------------------

      do 10 i = 1, nxs(7,iex)

*-----------------------------------------------------------------------

         ip = nint(xss(jxs(30,iex)+i-1))

*-----------------------------------------------------------------------
*        TEMPORARY STATEMENT to avoid trying to produce particles that
*        we are not transporting.  will go away after RCL writes an
*        appropriate expung for charged particle tables (also need to
*        account for heating in this new expung).
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*     determine total production cross section for this particle (t)
*-----------------------------------------------------------------------

         is = ixs(1,i,iex)+1
         ic = ktc(1,iex)+1-nint(xss(is-1))

         if(ic.lt.1.or.ic.gt.nint(xss(is)).or.ic.eq.nint(xss(is)).and.
     &    rtc(1,iex).ne.0.) goto 10

         t = xss(ic+is)
         if(rtc(1,iex).ne.0.) t = t+rtc(1,iex)*
     &    (xss(ic+is+1)-xss(ic+is))

         if(t.eq.0.) goto 10

*-----------------------------------------------------------------------
*         r is the ratio of total particle production cross section
*         to total interaction cross section and represents the average
*         number of particles produced per collision.
*-----------------------------------------------------------------------

         r = t/rtc(5,iex)

*-----------------------------------------------------------------------
*         produce 'n' of this type of secondary particles.  sample n
*         from nearest two integers to r.  weight of each particle is
*         simply the incident weight.
*-----------------------------------------------------------------------

            n   = int(r+rang())
            wgf = wgt






*-----------------------------------------------------------------------
*  TEMPORARY WRITE
*     if(n.gt.0.and.in.eq.9)write(iuo,11)
*    & n,int(xss(jxs(30,iex)+i-1))
*  11 format(' PRODUCING ',i3,' sec particles with IPT = ',i2)
*-----------------------------------------------------------------------

         if(n.eq.0) goto 10

         ipt = ip

*-----------------------------------------------------------------------

      do 40 j = 1, n

*-----------------------------------------------------------------------
*         sample reaction that produced the particle
*-----------------------------------------------------------------------

         ixre = 1
         if(nint(xss(jxs(31,iex)+i-1)).eq.1) goto 70

         a = 0.
         tt = 0.
         b = rang()*t

      do 80 ixre = 1, nint(xss(jxs(31,iex)+i-1))
         l = ixs(5,i,iex)+
     &    nint(xss(ixs(4,i,iex)+ixre-1))+1

            if(nint(xss(l-2)).ne.12) then

               write(ErrCha,*) 'Error : in acecp, '//
     &                    'mftype=12 expected but not found'
              ErrID = 'L:907/R:acecp/F:ggm08.f' !E04_026_001
              call ErrWrite(ErrID,ErrCha)

              call parastop( 416 )

            end if

*-----------------------------------------------------------------------
*         mf=12; partial cross section = yield * reaction cross section
*-----------------------------------------------------------------------

      do 90 jy = 1, 2

*-----------------------------------------------------------------------
*        This call to ACEFCN incorporates the same f "safety factor" as
*        is used in prdgam. This is probably there to avoid problems at
*        discontinuities in yield tables. It is unknown whether it is
*        strictly needed or not.
*-----------------------------------------------------------------------

   90    yd(jy) = acefcn(l,xss(jxs(1,iex)
     &          + ktc(1,iex)+jy-2)*f,ln)

         if(yd(1).eq.0..and.yd(2).eq.0.) goto 80

         ix = nint(xss(l-1))
         is = jxs(7,iex)+nint(xss(jxs(6,iex)+ix-1))
         x1 = 0.
         x2 = 0.

*-----------------------------------------------------------------------
*        RCL writes:
*        For the logic of the next 4 lines of the code, I originally
*        copied 4 lines from prdgam (ag4xq.9 through ag.133).  Upon
*        looking at the logic contained therein, I don't agree with
*        the treatment if we are looking up a cross section at an
*        energy beyond the last tabulated value for the particular
*        reaction.  prdgam, in such cases, always uses the last value
*        in the table.  I think it is better to use the logic found
*        in xstcol (ac.76 through ac.78).  That is what I have tried
*        to implement here.
*-----------------------------------------------------------------------

         ic = ktc(1,iex)-nint(xss(is-1))+1
         if(ic.gt.0.and.ic.le.nint(xss(is))) x1 = xss(is+ic)
         ic = ic+1
         if(ic.gt.0.and.ic.le.nint(xss(is))) x2 = xss(is+ic)
         s = x1*yd(1)+rtc(1,iex)*(x2*yd(2)-x1*yd(1))
         if(s.eq.0.) goto 80
         if(s.gt.tt) kx = ixre
         tt = max(s,tt)
         a = a+s
         if(a.gt.b) goto 70

   80 continue

*-----------------------------------------------------------------------
*        handle failure of reaction cross sections to sum
*        to the total.  treatment is same as in prdgam,
*        whereby a warning is issued (1st time) and the reaction
*        having the largest cross section is used.
*-----------------------------------------------------------------------

         call zaid(2,ht,ixl(1,iex))

         write(*,'(''Warning : in acecp, ipt='',i3,'' erg='',e13.4,//
     &   ''no particle-production mt found in acecp. zaid = '',a)')
     &   ipt, erg, ht

         ixre = kx

   70 continue

*-----------------------------------------------------------------------
*         call xstcas to sample secondary energy
*-----------------------------------------------------------------------

         erg = es
         ntyn = nint(xss(ixs(3,i,iex)+ixre-1))
         mtp = nint(xss(ixs(2,i,iex)+ixre-1))
         ia = ixs(7,i,iex)
         ka = nint(xss(ixs(6,i,iex)+ixre-1))
         id = ixs(9,i,iex)
         kd = nint(xss(ixs(8,i,iex)+ixre-1))
         call xstcas(1,in,zero,ia,ka,id,kd)
         call dtcos(colout(2,1),uold,uuu,0,irdm)

*-----------------------------------------------------------------------
*  TEMPORARY WRITE STATEMENT
*        if(in.eq.9) write(iuo,281) colout(1,1),colout(2,1)
* 281    format(5x,"sampled energy=",1pe12.5," and cosine=",e12.5,
*    &    " lab")
*-----------------------------------------------------------------------

         erg = colout(1,1)

*-----------------------------------------------------------------------
*     booking
*-----------------------------------------------------------------------

               call iptch(ipt,ityp,ktyp,ipat,ichg,ibry,iprt,inut,rmss)

               if( ityp .gt. 0 ) then

                     nclsts = nclsts + 1
                     iclusts(nclsts) = ipat
                     nixcos = nixcos + 1  ! MATSUDA 2022.12.12  for T-Point

                     jclusts(0,nclsts) = 0
                     jclusts(1,nclsts) = iprt
                     jclusts(2,nclsts) = inut
                     jclusts(3,nclsts) = ityp
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = ichg
                     jclusts(6,nclsts) = ibry
                     jclusts(7,nclsts) = ktyp
                     jclusts(8,nclsts) = 0

                     rms = rmss
                     rmg = rms / 1000.
                     pr  = sqrt( erg * ( erg + 2.0 * rms ) ) / 1000.
                     ett = sqrt( pr**2 + rmg**2 )

                     pxl = pr * uuu
                     pyl = pr * vvv
                     pzl = pr * www

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = pxl
                     qclusts(2,nclsts)  = pyl
                     qclusts(3,nclsts)  = pzl
                     qclusts(4,nclsts)  = ett
                     qclusts(5,nclsts)  = rmg
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = erg
                     qclusts(8,nclsts)  = wgf / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = 0.0
                     qclusts(11,nclsts) = 0.0
                     qclusts(12,nclsts) = 0.0
                     ixcoss(nixcos) = 0  ! MATSUDA 2022.12.12  for T-Point

                  if( in .eq. 1 ) then

                     if( ityp .eq. 1 ) then
                           icnm = 50
                     else if( ityp .ge. 3 .and. ityp .le. 5 ) then
                           icnm = 51
                     else
                           icnm = 52
                     end if

                  else if( in.eq.9.or.in.eq.31.or.in.eq.34 ) then ! T.Sato 2021/08/14, deuteron and alpha

                     if( ityp .eq. 2 ) then
                           icnm = 54
                     else if( ityp .ge. 3 .and. ityp .le. 5 ) then
                           icnm = 55
                     else
                           icnm = 56
                     end if
                     aevts(iaevt+61,ireg) =
     &               aevts(iaevt+61,ireg) + wgf
                     bevts(ibevt+61,imat) =
     &               bevts(ibevt+61,imat) + wgf
                  end if
                     aevts(iaevt+icnm,ireg) =
     &               aevts(iaevt+icnm,ireg) + wgf
                     bevts(ibevt+icnm,imat) =
     &               bevts(ibevt+icnm,imat) + wgf
                     numpal(ityp) = numpal(ityp) + 1
                     rumpal(ityp) = rumpal(ityp) + wgf

               end if

*-----------------------------------------------------------------------

   40 continue
   10 continue

         ipt = in
         erg = es

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sig_tot(ityp,sigt,sigaa,ein,mk)
*                                                                      *
*        this routine loads various cross sections in the rtc array    *
*        at energy erg for material mk.  it also loads totm and siga.  *
*                                                                      *
*        ktc(1,iex)=index in cross-section table (also ktc(2,iex))     *
*        rtc(1,iex)=table interpolation fraction (also rtc(2,iex))     *
*        rtc(3,iex)=absorption cross section                           *
*        rtc(4,iex)=total cross section [also rtc(5,iex)]              *
*        rtc(6,iex)=energy                                             *
*                                                                      *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************

      use MEMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      use GGMBANKMOD !FURUTA
      use moddas_material

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      common /eparm/  esmax, esmin, emin(20)
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)

      common /ndemax/ dnmax(20) ! T.Sato 2021/08/14
      common /emode/  emodem, ge1, ge2, iemode

*-----------------------------------------------------------------------
*     initialization of this material and particle
*-----------------------------------------------------------------------

         erg  = ein
*-----------------------------------------------------------------------

         totm = 0.
         siga = 0.

*-----------------------------------------------------------------------
*     now proton, deuteron, and alpha only  T.Sato 2021/08/14
*-----------------------------------------------------------------------
         if( ityp .eq. 1 .and. kpt(9) .ne. 0 ) then
            ipt = 9
         elseif( ityp .eq. 15 .and. kpt(31) .ne. 0 ) then
            ipt = 31
         elseif( ityp .eq. 18 .and. kpt(34) .ne. 0 ) then
            ipt = 34
         else
            return
         end if

*-----------------------------------------------------------------------
*     loop over number of nuclides in this material
*-----------------------------------------------------------------------

         isigc(0) = jmd(1+mk+1) - jmd(1+mk)

         mii = 0

      do m = jmd(1+mk), jmd(1+mk+1) - 1

*-----------------------------------------------------------------------

         mii = mii + 1
         siggc(mii) = 0.0d0
         siggcn(mii) = 0.d0
         siggce(mii) = 0.d0

         lem = m-jmd(1+mk)+1

cfrtati 2021/12/17 [data max] for ityp = 15, 18
         if(ityp.eq.1) then ! [data max] type
          dsmax = das_kmate(kmate(mk)+(lem-1)*5+11)
         else if(ityp.eq.15) then
          dsmax = das_kmate(kmate(mk)+(lem-1)*5+14)
         else if(ityp.eq.18) then
          dsmax = das_kmate(kmate(mk)+(lem-1)*5+15)
         endif

      if( dsmax .gt. emin(ityp) ) then

*-----------------------------------------------------------------------

         iex = lme(ipt,m)
         if(rtc(6,iex).ne.erg) then
            n = nxs(3,iex)
            ic = jxs(1,iex)
            ib = ic+n-1
            if(erg.lt.rtc(6,iex)) ib = ic+ktc(1,iex)
            if(erg.gt.rtc(6,iex)) ic = ic+ktc(1,iex)-1
   60       if(ib-ic.eq.1) goto 80
            ih = (ic+ib)/2
            if(erg.lt.xss(ih)) goto 70
            ic = ih
            goto 60
   70       ib = ih
            goto 60

*-----------------------------------------------------------------------
*           get the bin index, the interpolation fraction,
*           and the absorption and total cross sections.
*-----------------------------------------------------------------------

   80       ktc(1,iex) = ic-jxs(1,iex)+1
            ex = (erg-xss(ic))/(xss(ib)-xss(ic))
            rtc(1,iex) = max(zero,min(one,ex))
            ktc(2,iex) = ktc(1,iex)
            rtc(2,iex) = rtc(1,iex)
            rtc(3,iex) = xss(ic+2*n)+rtc(1,iex)*
     &         (xss(ic+2*n+1)-xss(ic+2*n))
            rtc(4,iex) = xss(ic+n)+rtc(1,iex)*
     &         (xss(ib+n)-xss(ic+n))
            rtc(5,iex) = rtc(4,iex)
            rtc(6,iex) = erg
cfrtati 2021/12/17 elastic cross-section
            rtcel(iex) = xss(ic+3*n)+rtc(1,iex)*
     &         (xss(ic+3*n+1)-xss(ic+3*n))
         end if
*-----------------------------------------------------------------------

         totm = totm+rtc(5,iex)*fme(m)
         siga = siga+rtc(3,iex)*fme(m)

         siggc(mii) = rtc(5,iex)*fme(m)
         siggcn(mii) = (rtc(5,iex)-rtcel(iex))*fme(m) ! frtati 2021/12/17

         isigza(mii) = iza(m)
         if( iemode.gt.0 ) then
           siggcmx(1,mii) = sigmxlb(1,iex)*fme(m)
           siggcmx(2,mii) = sigmxlb(2,iex)*fme(m)
           itz = iza(m) / 1000
           ita = iza(m) - 1000 * itz
           if ( ita.eq.0 .and. itz.eq.6 ) ita = 12 ! frtati 2022/03/11
           call sigrc(1,ein,ita,itz,sigt,signe,sigel)
           siggce(mii) = sigel*fme(m)
         end if

      else

        if( iza(m).ne.1001 ) then ! frtati 2022/03/11
          if(ityp.eq.1) then
            siggcmx(1,mii) = sigmxp(1,m)*fme(m)
            siggcmx(2,mii) = sigmxp(2,m)*fme(m)
          else if(ityp.eq.15) then
            siggcmx(1,mii) = sigmxd(1,m)*fme(m)
          else if(ityp.eq.18) then
            siggcmx(1,mii) = sigmxa(1,m)*fme(m)
          end if
        else
          siggcmx(1,mii) = 0.d0
          siggcmx(2,mii) = 0.d0
        end if

      end if

*-----------------------------------------------------------------------

      end do

*-----------------------------------------------------------------------

         sigt  = totm
         sigaa = siga

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine pnctot(mk,totpn,ergi)
*                                                                      *
*        calculate the photonuclear cross section in material mk.      *
*                                                                      *
*        local variables: ii - isotope index, it - table index,        *
*          ib - upper energy index, ic - lower energy index,           *
*                                                                      *
*        global variables:                                             *
*          jmd - material isotope indices                              *
*          lmn - photonuclear isotope table indices                    *
*          erg - current photon energy                                 *
*          pnt(mk) - the lowest photonuclear threshold in material mk  *
*          ktc(1,table) - table index above current energy             *
*          rtc(1,table) - current linear table interpolation factor    *
*          rtc(2,table) - current total cross section                  *
*          rtc(6,table) - current energy being interpolated            *
*                                                                      *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      use MEMBANKMOD ! frtati 2021/12/17
      use moddas_material
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------
      common /eparm/  esmax, esmin, emin(20)
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)

cABE 2024/05/07
      common /pnixs/ totmpni,pnixsm(pnlmax),ipnikind(pnlmax),egypni,
     &               dlibmax
!$OMP THREADPRIVATE(/pnixs/)
cABE end

*-----------------------------------------------------------------------

         erg = ergi
         totpn = 0.
cABE 2024/05/07
         totmpni = 0.0d0
cABE end

*-----------------------------------------------------------------------
*        return if below photonuclear threshold for material.
*-----------------------------------------------------------------------

         if(erg.le.pnt(mk)) then ! T.Sato 2022/03/27
          isigc(0)=0
          return
         endif

*-----------------------------------------------------------------------
*        accumulate the cross sections by nuclide in material.
*-----------------------------------------------------------------------

         mii = 0
         isigc(0) = jmd(1+mk+1) - jmd(1+mk)

      do 60 ii = jmd(1+mk), jmd(1+mk+1) - 1

         mii = mii + 1
         siggc(mii) = 0.0d0
         lem = ii-jmd(1+mk)+1
         dsmax = das_kmate(kmate(mk)+(lem-1)*5+13)
         if( dsmax.gt.emin(14) ) then

         it = lmn(ii)
         if(it.eq.0) goto 60
         if(erg.eq.rtc(6,it)) goto 50
         rtc(6,it) = erg

*-----------------------------------------------------------------------
*        find index of current energy.
*-----------------------------------------------------------------------

         ic = jxs(1,it)
         ib = ic+nxs(3,it)-1

*-----------------------------------------------------------------------
*        if below photonuclear threshold, use zero values.
*-----------------------------------------------------------------------

         if(erg.ge.xss(ic)) goto 10
         ktc(1,it) = 0
         rtc(1,it) = 0.
         rtc(2,it) = 0.
         goto 50

*-----------------------------------------------------------------------
*        if above last energy value, use boundary values.
*-----------------------------------------------------------------------

   10    if(erg.lt.xss(ib)) goto 20
         ktc(1,it) = ic-jxs(1,it)
         rtc(1,it) = 0.
         rtc(2,it) = xss(jxs(2,it)+ktc(1,it))
         goto 50

*-----------------------------------------------------------------------
*        else handle normal value within table.
*-----------------------------------------------------------------------

   20    if(ib-ic.eq.1) goto 40
         ih = (ic+ib)/2
         if(erg.lt.xss(ih)) goto 30
         ic = ih
         goto 20

   30    ib = ih
         goto 20

   40    ktc(1,it) = ib-jxs(1,it)
         rtc(1,it) = (erg-xss(ic))/(xss(ib)-xss(ic))
         k = ktc(1,it)+jxs(2,it)
         rtc(2,it) = (xss(k)-xss(k-1))*rtc(1,it)+xss(k-1)

*-----------------------------------------------------------------------
*        accumulate the photonuclear microscopic xs for the material.
*-----------------------------------------------------------------------

   50    totpn = totpn+rtc(2,it)*fme(ii)
         siggc(mii) = rtc(2,it)*fme(ii) ! frtati 2021/12/17
cABE 2024/05/07
         pnixsm(mii) = rtc(2,it)
         totmpni = totmpni + pnixsm(mii) * fme(ii)
cABE end
         isigza(mii) = iza(ii) ! S.H. 2024.10.15

         end if ! frtati 2021/12/17
   60    continue

*-----------------------------------------------------------------------

      return
      end

************************************************************************                                                                      *
      subroutine SET_LIBCSMAX() ! frtati 2021/12/17
*     set non-elastic cross-section maximum for high-energy event mode *
************************************************************************
      use GGMARRAYMOD
      use MEMBANKMOD
      use moddas_material
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------
      include 'param.inc'
      include 'param01.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------
      common /eparm/  esmax, esmin, emin(20)
      common /ndemax/ dnmax(20)
      common /emode/  emodem, ge1, ge2, iemode
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)
      common /natnuc/ natnn(maxpt), natnm(maxpt,10), patnn(maxpt,10)
      parameter( nemd = 300 )
      parameter( smxmargin = 1.1d0 )

      integer,dimension(1000) :: ia_za, ia_nn

*-----------------------------------------------------------------------

      sigmxlb = 0.d0
      sigmxp = 0.d0
      sigmxd = 0.d0
      sigmxa = 0.d0

      if( dnmax(1).gt.emin(1) .or. dnmax(15).gt.emin(15) .or.
     &    dnmax(18).gt.emin(18) ) then

        if( iemode.gt.0 ) then
          do iex = 1, mxe
            do km = 1, mix
              if( lme(9,km).eq.iex .or. lme(31,km).eq.iex .or.
     &            lme(34,km).eq.iex ) then
                i1_nesz = jxs(1,iex)
                i2_nnes = nxs(3,iex)
                t0_non = 0.d0
                t1_tot = 0.d0
                do k = 1, i2_nnes
                  i4_tot = i1_nesz + i2_nnes + k - 1
                  i6_ela = i1_nesz + 3*i2_nnes + k - 1
                  tt_non = xss(i4_tot)-xss(i6_ela)
                  if( tt_non.gt.t0_non ) t0_non = tt_non
                  if( lme(9,km).eq.iex ) then
                    ityp = 1
                    i3_ein = i1_nesz + k - 1
                    ein = xss(i3_ein)
                    itz = iza(km) / 1000
                    ita = iza(km) - 1000 * itz
                    if( ita.ne.0 ) then
                      call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)
                      tt_ela = sigel
                    else
                      tt_ela = 0.d0
                      do in = 1, natnn(itz)
                       i0_ita = natnm(itz,in)
                       call sigrc(ityp,ein,i0_ita,itz,sigt,signe,sigel)
                       tt_ela = tt_ela+0.01*patnn(itz,in)*sigel
                      end do
                    end if
                    tt_tot = tt_non + tt_ela
                    if( tt_tot.gt.t1_tot ) t1_tot = tt_tot
                  end if
                end do
                sigmxlb(1,iex) = t0_non
                sigmxlb(2,iex) = t1_tot
                exit
              end if
            end do
          end do
        end if

        ia_za = 0
        ia_nn = 0
        i0_em = 1
        mn = 1
        do km = 1, mix
          if ( km.ge.jmd(1+mn+1) ) mn = mn + 1
          i0_lem = km-jmd(1+mn)+1
          im = 0
          do i = 1, i0_em
            if( iza(km).eq.ia_za(i) ) im = i
          end do
          if( im.ne.0 ) then
            sigmxp(1,km) = sigmxp(1,ia_nn(im))
            sigmxp(2,km) = sigmxp(2,ia_nn(im))
            sigmxd(1,km) = sigmxd(1,ia_nn(im))
            sigmxa(1,km) = sigmxa(1,ia_nn(im))
          else
            ia_za(i0_em) = iza(km)
            ia_nn(i0_em) = km
            i0_em = i0_em + 1
            itz = iza(km) / 1000
            ita = iza(km) - 1000 * itz
            dsmax = das_kmate(kmate(mn)+(i0_lem-1)*5+11)
            if( dsmax.lt.dnmax(1) ) then
              ityp = 1
              edel = log(dnmax(1)) / dble(nemd - 1)
              t0_non = 0.d0
              t1_tot = 0.d0
              do ie = 1, nemd
                ein = dexp( (ie-1)*edel )
                if( ita.ne.0 ) then
                  call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)
                  tt_non = signe
                  tt_tot = sigt
                else
                  tt_non = 0.d0
                  tt_tot = 0.d0
                  do in = 1, natnn(itz)
                    i0_ita = natnm(itz,in)
                    call sigrc(ityp,ein,i0_ita,itz,sigt,signe,sigel)
                    tt_non = tt_non+0.01*patnn(itz,in)*signe
                    tt_tot = tt_tot+0.01*patnn(itz,in)*sigt
                  end do
                end if
                if( tt_non.gt.t0_non ) t0_non = tt_non
                if( tt_tot.gt.t1_tot ) t1_tot = tt_tot
              end do
              sigmxp(1,km) = t0_non
              sigmxp(2,km) = t1_tot
            end if
            dsmax = das_kmate(kmate(mn)+(i0_lem-1)*5+14)
            if( dsmax.lt.dnmax(15) ) then
              ap = 2.d0
              zp = 1.d0
              edel = log(dnmax(15)/esmin)/dble(nemd-1)
              t0_non = 0.d0
              do ie = 1, nemd
                ein = dexp((ie-1)*edel)*esmin*ap
                if( ita.ne.0 ) then
                  call sighi(ap,zp,ein,dble(ita),dble(itz)
     &                      ,signe,sigel,bmax)
                  tt_non = signe
                else
                  tt_non = 0.d0
                  do in = 1, natnn(itz)
                    i0_ita = natnm(itz,in)
                    call sighi(ap,zp,ein,dble(i0_ita),dble(itz)
     &                        ,signe,sigel,bmax)
                    tt_non = tt_non+0.01*patnn(itz,in)*signe
                  end do
                end if
                if( tt_non.gt.t0_non ) t0_non = tt_non
              end do
              sigmxd(1,km) = t0_non
            end if
            dsmax = das_kmate(kmate(mn)+(i0_lem-1)*5+15)
            if( dsmax.lt.dnmax(18) ) then
              ap = 4.d0
              zp = 2.d0
              edel = log(dnmax(18)/esmin)/dble(nemd-1)
              t0_non = 0.d0
              do ie = 1, nemd
                ein = dexp((ie-1)*edel)*esmin*ap
                if( ita.ne.0 ) then
                  call sighi(ap,zp,ein,dble(ita),dble(itz)
     &                      ,signe,sigel,bmax)
                  tt_non = signe
                else
                  tt_non = 0.d0
                  do in = 1, natnn(itz)
                    i0_ita = natnm(itz,in)
                    call sighi(ap,zp,ein,dble(i0_ita),dble(itz)
     &                        ,signe,sigel,bmax)
                    tt_non = tt_non+0.01*patnn(itz,in)*signe
                  end do
                end if
                if( tt_non.gt.t0_non ) t0_non = tt_non
              end do
              sigmxa(1,km) = t0_non
            end if
          end if
        end do

        sigmxlb(:,:) = smxmargin*sigmxlb(:,:)
        sigmxp(:,:) = smxmargin*sigmxp(:,:)
        sigmxd(:,:) = smxmargin*sigmxd(:,:)
        sigmxa(:,:) = smxmargin*sigmxa(:,:)
      end if

      return
      end

************************************************************************                                                                      *
      subroutine datamaxsummary() ! frtati 2022/03/20
*     set maximum energy of library for data max                       *
************************************************************************
      use GGMARRAYMOD
      use moddas_material
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------
      common /ndemax/ dnmax(20)
      common /dpnmaxcom/ dpnmax
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1g/ kmat(kvlmax)
      common /kmat1ka/ kmathd(kvlmax), kmathe(kvlmax)
      common /kmat1o/ iom1, iom2, iom3

      real*8, allocatable :: a0_elibmx(:)
      character ht*10
      character pjname(5)*12
      data pjname/'proton','neutron','photonuclear','deuteron','alpha'/

*-----------------------------------------------------------------------
      i1_flag = 0
      i2_flag = 0
      allocate(a0_elibmx(1:mxe))
      a0_elibmx = 0.d0
      do i = 1, mxe
        a0_elibmx(i) = xss(jxs(1,i)+nxs(3,i)-1)
      end do
      do m = 1, mipt
        if(m.eq.1 .or. m.eq.9 .or. m.eq.31.or. m.eq.34 ) then
          mn = 1
          do km = 1, mix
            if ( km.ge.jmd(1+mn+1) ) mn = mn + 1
            lem_tem = km-jmd(1+mn)+1
            iwarning = 0
            if( lme(m,km).ne.0 ) then
              dlibmax = 0.d0
              t0_elibmx = 0.d0
              if( m.eq.1 ) then
                if( das_kmate(kmate(mn)+(lem_tem-1)*5+12) .gt.
     &              a0_elibmx(lme(1,km)) ) then
                  iwarning = 1
                  dlibmax = das_kmate(kmate(mn)+(lem_tem-1)*5+12)
                  t0_elibmx = a0_elibmx(lme(1,km))
                end if
              else if( m.eq.9 ) then
                if( das_kmate(kmate(mn)+(lem_tem-1)*5+11) .gt.
     &              a0_elibmx(lme(9,km)) ) then
                  iwarning = 1
                  dlibmax = das_kmate(kmate(mn)+(lem_tem-1)*5+11)
                  t0_elibmx = a0_elibmx(lme(9,km))
                end if
              else if( m.eq.31 ) then
                if( das_kmate(kmate(mn)+(lem_tem-1)*5+14) .gt.
     &              0.5d0*a0_elibmx(lme(31,km)) ) then
                  iwarning = 1
                  dlibmax = das_kmate(kmate(mn)+(lem_tem-1)*5+14)
                  t0_elibmx = 0.5d0*a0_elibmx(lme(31,km))
                end if
              else if( m.eq.34 ) then
                if( das_kmate(kmate(mn)+(lem_tem-1)*5+15) .gt.
     &              0.25d0*a0_elibmx(lme(34,km)) ) then
                  iwarning = 1
                  dlibmax = das_kmate(kmate(mn)+(lem_tem-1)*5+15)
                  t0_elibmx = 0.25d0*a0_elibmx(lme(34,km))
                end if
              end if
              if( iwarning.eq.1 ) then
                i1_flag = 1
                call zaid(2,ht,ixl(1,lme(m,km)))
                ErrID = 'L:1688/R:datamaxsummary/F:ggm08.f'
                write(ErrCha,'(" matID =",i5,": Maximum energy of ",
     &          a10," =",f8.2," is smaller than dmax =",f8.2)')
     &          nmt(mn),ht,t0_elibmx,dlibmax
                call ErrWrite(ErrID,ErrCha)
              end if
            end if
          end do
        end if
      end do
      mn = 1
      do km = 1, mix
        if ( km.ge.jmd(1+mn+1) ) mn = mn + 1
        lem_tem = km-jmd(1+mn)+1
        if( lmn(km).ne.0 ) then
          if( das_kmate(kmate(mn)+(lem_tem-1)*5+13) .gt.
     &        a0_elibmx(lmn(km)) )then
            i2_flag = 1
            call zaid(2,ht,ixl(1,lmn(km)))
            ErrID = 'L:1707/R:datamaxsummary/F:ggm08.f'
            write(ErrCha,'(" matID =",i5,": Maximum energy of ",
     &      a10," =",f8.2," is smaller than dmax =",f8.2)') nmt(mn),
     &      ht,a0_elibmx(lmn(km)),
     &      das_kmate(kmate(mn)+(lem_tem-1)*5+13)
            call ErrWrite(ErrID,ErrCha)
          end if
        end if
      end do

      if ( i1_flag.eq.1 ) then ! set data max

      do m = 1, mipt
        if ( m.eq.1 .or. m.eq.9 .or. m.eq.31.or. m.eq.34 ) then
          mn = 1
          do km = 1, mix
            if ( km.ge.jmd(1+mn+1) ) mn = mn + 1
            lem_tem = km-jmd(1+mn)+1
            if( lme(m,km).ne.0 ) then
              if( m.eq.1 ) then
                if( das_kmate(kmate(mn)+(lem_tem-1)*5+12) .gt.
     &              a0_elibmx(lme(1,km)) ) then
                 das_kmate(kmate(mn)+(lem_tem-1)*5+12) =
     &           a0_elibmx(lme(1,km))
                end if
              else if( m.eq.9 ) then
                if( das_kmate(kmate(mn)+(lem_tem-1)*5+11) .gt.
     &              a0_elibmx(lme(9,km)) ) then
                 das_kmate(kmate(mn)+(lem_tem-1)*5+11) =
     &           a0_elibmx(lme(9,km))
                end if
              else if( m.eq.31 ) then
                if( das_kmate(kmate(mn)+(lem_tem-1)*5+14) .gt.
     &              0.5d0*a0_elibmx(lme(31,km)) ) then
                  das_kmate(kmate(mn)+(lem_tem-1)*5+14) =
     &            0.5d0*a0_elibmx(lme(31,km))
                end if
              else if( m.eq.34 ) then
                if( das_kmate(kmate(mn)+(lem_tem-1)*5+15) .gt.
     &              0.25d0*a0_elibmx(lme(34,km)) ) then
                  das_kmate(kmate(mn)+(lem_tem-1)*5+15) =
     &            0.25d0*a0_elibmx(lme(34,km))
                end if
              end if
            end if
          end do
        end if
      end do
      do m = 1, mxmat
        t1_dmax = 0.d0
        t2_dmax = 1.d10
        t3_dmax = 0.d0
        t4_dmax = 1.d10
        t5_dmax = 0.d0
        t6_dmax = 1.d10
        t7_dmax = 0.d0
        t8_dmax = 1.d10
        nel = nint( dnel_das(kmat0+m) )
        do l = 1, nel
          iz = nint( zz_das(kmat(m)+l) )
          ia = nint( a_das(kmat(m)+l) )
          mn = 1
          do km = 1, mix
            if( km .ge. jmd(1+mn+1) ) mn = mn + 1
            lem_tem = km-jmd(1+mn)+1
            itz = iza(km) / 1000
            ita = iza(km) - 1000 * itz
            if( m.eq.mn .and. iza(km).eq.1001 ) then
              das_kmatd(kmatd(m)+6) =
     &        das_kmate(kmate(m)+(lem_tem-1)*5+11)
              das_kmatd(kmatd(m)+7) =
     &        das_kmate(kmate(m)+(lem_tem-1)*5+12)
              das_kmatd(kmatd(m)+9) =
     &        das_kmate(kmate(m)+(lem_tem-1)*5+14)
              das_kmatd(kmatd(m)+10) =
     &        das_kmate(kmate(m)+(lem_tem-1)*5+15)
            end if
            if( (m.eq.mn .and. iz.eq.itz) .and.
     &          (ia.eq.ita .or. ita.eq.0) ) then
              das_kmatd(kmatd(m)+(l-1)*5+11) =
     &        das_kmate(kmate(m)+(lem_tem-1)*5+11)
              das_kmatd(kmatd(m)+(l-1)*5+12) =
     &        das_kmate(kmate(m)+(lem_tem-1)*5+12)
              das_kmatd(kmatd(m)+(l-1)*5+14) =
     &        das_kmate(kmate(m)+(lem_tem-1)*5+14)
              das_kmatd(kmatd(m)+(l-1)*5+15) =
     &        das_kmate(kmate(m)+(lem_tem-1)*5+15)
              exit
            end if
          end do
          if ( das_kmatd(kmatd(m)+(l-1)*5+12).gt.t1_dmax ) then
            t1_dmax = das_kmatd(kmatd(m)+(l-1)*5+11)
          end if
          if ( das_kmatd(kmatd(m)+(l-1)*5+12).lt.t2_dmax ) then
            t2_dmax = das_kmatd(kmatd(m)+(l-1)*5+11)
          end if
          if ( das_kmatd(kmatd(m)+(l-1)*5+12).gt.t3_dmax ) then
            t3_dmax = das_kmatd(kmatd(m)+(l-1)*5+12)
          end if
          if ( das_kmatd(kmatd(m)+(l-1)*5+12).lt.t4_dmax ) then
            t4_dmax = das_kmatd(kmatd(m)+(l-1)*5+12)
          end if
          if ( das_kmatd(kmatd(m)+(l-1)*5+14).gt.t5_dmax ) then
            t5_dmax = das_kmatd(kmatd(m)+(l-1)*5+14)
          end if
          if ( das_kmatd(kmatd(m)+(l-1)*5+14).lt.t6_dmax ) then
            t6_dmax = das_kmatd(kmatd(m)+(l-1)*5+14)
          end if
          if ( das_kmatd(kmatd(m)+(l-1)*5+15).gt.t7_dmax ) then
            t7_dmax = das_kmatd(kmatd(m)+(l-1)*5+15)
          end if
          if ( das_kmatd(kmatd(m)+(l-1)*5+15).lt.t8_dmax ) then
            t8_dmax = das_kmatd(kmatd(m)+(l-1)*5+15)
          end if
        end do
        das_kmatd(kmatd(m)+1) = t1_dmax
        das_kmatd(kmatd(m)+2) = t2_dmax
        das_kmatd(kmatd(m)+3) = t3_dmax
        das_kmatd(kmatd(m)+4) = t4_dmax
        das_kmathd(kmathd(m)+3) = t5_dmax
        das_kmathd(kmathd(m)+4) = t6_dmax
        das_kmathd(kmathd(m)+5) = t7_dmax
        das_kmathd(kmathd(m)+6) = t8_dmax
      end do

      end if

      if( i2_flag.eq.1 ) then

      mn = 1
      do km = 1, mix
        if ( km .ge. jmd(1+mn+1) ) mn = mn + 1
        lem_tem = km-jmd(1+mn)+1
        if( lmn(km).ne.0 ) then
          if( das_kmate(kmate(mn)+(lem_tem-1)*5+13) .gt.
     &        a0_elibmx(lmn(km)) ) then
            das_kmate(kmate(mn)+(lem_tem-1)*5+13) = a0_elibmx(lmn(km))
          end if
        end if
      end do
      do m = 1, mxmat
        t1_dmax = 0.d0
        t2_dmax = 1.d10
        nel = nint( dnel_das(kmat0+m) )
        do l = 1, nel
          iz = nint( zz_das(kmat(m)+l) )
          ia = nint( a_das(kmat(m)+l) )
          mn = 1
          do km = 1, mix
            if( km .ge. jmd(1+mn+1) ) mn = mn + 1
            lem_tem = km-jmd(1+mn)+1
            itz = iza(km) / 1000
            ita = iza(km) - 1000 * itz
            if( (m.eq.mn .and. iz.eq.itz) .and.
     &          (ia.eq.ita .or. ita.eq.0) ) then
              das_kmatd(kmatd(m)+(l-1)*5+13) =
     &        das_kmate(kmate(m)+(lem_tem-1)*5+13)
              exit
            end if
          end do
          if ( das_kmatd(kmatd(m)+(l-1)*5+13).gt.t1_dmax ) then
            t1_dmax = das_kmatd(kmatd(m)+(l-1)*5+13)
          end if
          if ( das_kmatd(kmatd(m)+(l-1)*5+13).lt.t2_dmax ) then
            t2_dmax = das_kmatd(kmatd(m)+(l-1)*5+13)
          end if
        end do
        das_kmathd(kmathd(m)+1) = t1_dmax
        das_kmathd(kmathd(m)+2) = t2_dmax
      end do

      end if

      deallocate(a0_elibmx)

*-----------------------------------------------------------------------
*     Output datamax for each nuclei, T.Sato 2022/03/21
*-----------------------------------------------------------------------
      iuo = iom1 ! frtati 2022/03/22
      write(iuo,1010)
      write(iuo,1020)
 1010 format(/'*** Summary of maximum library energy for each nuclide')
 1020 format('   matID iZ  iA    dmax projectile')
      do i = 1, mxa ! frtati 2022/03/25 display only used material
       do m = 1, mxmat ! material loop
        if( mat(i).eq.m ) then
         nel = nint( dnel_das(kmat0+m) ) ! number of nuclide in each material
         do l = 1, nel
          iztmp = nint( zz_das(kmat(m)+l) )
          iatmp = nint( a_das(kmat(m)+l) )
          do itmp=1,5 ! proton, neutron, photonuclear, deuteron, alpha
           if(das_kmatd(kmatd(m)+(l-1)*5+itmp+10).gt.1.0e-3)
     &     write(iuo,'(i8,i3,i4,f8.1,1x,a12)')
     &     nmt(m),iztmp,iatmp,das_kmatd(kmatd(m)+(l-1)*5+itmp+10),
     &     pjname(itmp)
          enddo
         enddo
         exit
        endif
       enddo
      enddo

      return
      end

************************************************************************
      subroutine displaymtinfo(iot,m) ! frtati 2023/12/07
************************************************************************
      use GGMBANKMOD
      use MEMBANKMOD
      use GGMARRAYMOD, only:jmd, lme
      use partmod
      use moddas_mesh
      use moddas_tally
      use moddas_material
      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------
      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------
      common /multipl/ imltp,iimlt(multmax),inmlt(multmax),
C MATSUDA 2024.11.25 (lagrange: ilmlt)
     &  ilmlt(multmax),
     &  idmlt(multmax),ismlt(multmax),impan(multmax),impat(multmax,6,3),
     &                 jmpat(multmax,6,6,2), imdfl(multmax)
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall84/ itmto(itlmax,6,4)
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /eparm/  esmax, esmin, emin(20)
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)
      common /ptname/ pname(20), ipln(20)
      common /cMeVperu/ iMeVperu

      real*8, allocatable :: fac_mto(:,:)
      character :: pname*8, elmnt(104)*3
      character :: cht*200, chn*3, chp*8

*-----------------------------------------------------------------------
      data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------
cfrtati 2024/01/05 check library needed for multiplier
      do k = 1, itmlp(m)
        do j = 1, itmpn(m,k)
          do jj = 1, itmln(m,k)
            do it = 1, itmst(m)
              nn = it * 2 - 1
              iiin = mltp(1+nn,itmli(m,k)/13+jj)
              iifn = mltp(1+nn+1,itmli(m,k)/13+jj)
              do ii = iiin, iifn, 2
                if( nint(slib(ii)).eq.5 ) then
                  imat = nint(slib(ii+1))
                else if( nint(slib(ii)).eq.6 ) then
                  if( itmpt(m,k,j,1).ne.2 .and.
     &                itmpt(m,k,j,1).ne.14 ) then
                    kmat = idnm(imat)
                    do mm = jmd(1+kmat), jmd(1+kmat+1) - 1
                      dsmax = 0.d0
                      lem_tem = mm-jmd(1+kmat)+1
                      if( itmpt(m,k,j,1).eq.1 ) then
                        dsmax = das_kmate(kmate(mm)+(lem_tem-1)*5+11)
                      else if( itmpt(m,k,j,1).eq.15 ) then
                        dsmax = das_kmate(kmate(mm)+(lem_tem-1)*5+14)
                      else if( itmpt(m,k,j,1).eq.18 ) then
                        dsmax = das_kmate(kmate(mm)+(lem_tem-1)*5+15)
                      end if
                      if( dsmax.le.emin(itmpt(m,k,j,1)) ) then
                        write(ErrCha,'("Library needed for multiplier",
     &                                 " is not loaded")')
                        ErrID = 'L:1995/R:displaymtinfo/F:ggm08.f'
                        call ErrWrite(ErrID,ErrCha)
                        stop
                      end if
                    end do
                  end if
                end if
              end do
            end do
          end do
        end do
      end do

      if( maxval(itmto(m,:,1)).ne.1 ) return

      ialloc_ggmbank = 0
      ialloc_membank = 0
      if( .not. allocated(rtc) ) then
        call ALLOCATE_GGMBANK
        call INIT_GGMBANK
        ialloc_ggmbank = 1
      end if
      if( .not. allocated(siggc) ) then
        call ALLOCATE_MEMBANK
        ialloc_membank = 1
      end if
      cmvpu = 1.d0
      cinv = 1.d0
      mulnum = -250
      do k = 1, itmlp(m)
        if( itmto(m,k,1).eq.1 ) then
          allocate(fac_mto(0:itmto(m,k,3),itmst(m)))
          do j = 1, itmpn(m,k)
            ityp = itmpt(m,k,j,1)
            ktyp = itmpt(m,k,j,2)
            ipat(m,1,1) = ityp
            ipat(m,1,2) = ktyp
            do jj = 1, itmln(m,k)
              mat = mltp(1,itmli(m,k)/13+jj)
              if( mat.eq.-1 ) then ! when mat is all
               do matinput = 1, mxmat
                  ntotlme = 0
                do melement = jmd(1+matinput), jmd(1+matinput+1) - 1
                   ntotlme = ntotlme + lme(1,melement) ! only for neutron
                end do
                if ( ntotlme .ne. 0 ) exit
               end do
               mat = matinput
              end if
              do ie = 0, itmto(m,k,3)
                if( iMeVperu.eq.1 ) then
                  if( ityp.lt.15 ) then
                    cmvpu = 1.d0
                  else if( ityp.eq.15 ) then
                    cmvpu = 2.d0
                  else if( ityp.eq.16 .or. ityp.eq.17 ) then
                    cmvpu = 3.d0
                  else if( ityp.eq.18 ) then
                    cmvpu = 4.d0
                  end if
                end if
                erg_mto = cmvpu * gmsh(itmto(m,k,4)+ie)
                if( ityp.ne.20 ) then
                  call fmfac(m,1,erg_mto,1.d0,fac_mto(ie,:))
                end if
              end do
              do it = 1, itmst(m)
                nn = it * 2 - 1
                iiin = mltp(1+nn,itmli(m,k)/13+jj)
                iifn = mltp(1+nn+1,itmli(m,k)/13+jj)
                itabuse = 0
                ia = 0
                cht = " "
                do ii = iiin, iifn, 2
                  if( nint(slib(ii)).eq.4 ) then
                    cinv = 1.d0/slib(ii+1)
                  else if( nint(slib(ii)).eq.6 ) then
                    write(chn,'(I3)') nint(slib(ii+1))
                    chn = adjustL(chn)
                    cht = trim(cht) // " " // trim(chn)
                    ia = ia + 1
                  else if( nint(slib(ii)).eq.3 ) then
                    cht = trim(cht) // ":"
                  else if( nint(slib(ii)).eq.2 ) then
                    if( nint(slib(ii-2)).eq.5 .and.
     &                  nint(slib(ii-1)).lt.-100 ) then
                      itabuse = 1
                      if( ityp.eq.20 ) itabuse = 2
                      idtab =  nint(slib(ii-1))
                    end if
                  else if( nint(slib(ii)).eq.5 ) then
                     imat = nint(slib(ii+1))
                     if ( imat .ne. 0 ) mat = idnm( imat )
                  end if
                end do
                if( mltp(1,itmli(m,k)/13+jj).eq.-1 ) then
                  write(iot,'("# material ID =",I6,"(all), mset",I1)')
     &            idmn(mat), it
                else
                  write(iot,'("# material ID =",I6,", mset",I1)')
     &            idmn(mat), it
                end if
                if( itabuse.eq.0 ) then
                  if( ia.gt.1 ) then
                    write(iot,'(A)') '# MT = (' // trim(cht) // ' )'
                  else
                    write(iot,'(A)') "# MT =" // trim(cht)
                  end if
                else
                  write(iot,'("# ID = ",I4)') idtab
                end if
                if( itabuse.ne.2 ) then
                  mulnum = mulnum - 1
                  write(iot,'("[ Multiplier ]")')
                  write(iot,'(" number = ",I4)') mulnum
                  write(iot,'(" interpolation = lin")')
                  write(iot,'(" part = ",a8)') pname(ityp)
                  write(iot,'(" ne = ",I6)') itmto(m,k,3) + 1
                  do ie = 0, itmto(m,k,3)
                    erg_mto = cmvpu * gmsh(itmto(m,k,4)+ie)
                    write(iot,'(2es10.3)') erg_mto, fac_mto(ie,it)*cinv
                  end do
                  write(iot,*)
                else
                  ityp_save = ityp
                  ktyp_save = ktyp
                  do kk = 1, imltp
                    if( idmlt(kk).eq.idtab ) then
!                      mulnum = mulnum - 1   ! not necessary to change multiplier ID
                      ityp = impat(kk,1,1)
                      ktyp = impat(kk,1,2)
                      chp = pname(ityp)
                      if( ityp.eq.19 ) then
                        jtyp = abs(ktyp)/1000000
                        if( ktyp.lt.0 ) ktyp = -ktyp + 1
                        chp = elmnt(jtyp)
                      end if
                      write(iot,'("[ Multiplier ]")')
                      write(iot,'(" number = ",I4)') mulnum
                      write(iot,'(" interpolation = lin")')
                      write(iot,'(" part = ",a8)') chp
                      write(iot,'(" ne = ",I6)') itmto(m,k,3) + 1
                      do ie = 0, itmto(m,k,3)
                        erg_mto =  gmsh(itmto(m,k,4)+ie)
                        write(iot,'(2es10.3)')
     &                  erg_mto, dosf(erg_mto,idtab)
!  for output reciprocal of dose conversion coefficient used for [t-wwg]
!        dcctmp=0.0
!        if(dosf(erg_mto,idtab).ne.0) dcctmp=1.0/dosf(erg_mto,idtab)
!        write(iot,'(2es10.3)') erg_mto, dcctmp
                      end do
                      write(iot,*)
                    end if
                  end do
                  ityp = ityp_save
                  ktyp = ktyp_save
                end if
              end do
            end do
          end do
          deallocate(fac_mto)
        end if
      end do

      if( ialloc_ggmbank.eq.1 ) call DEALLOCATE_GGMBANK
      if( ialloc_membank.eq.1 ) call DEALLOCATE_MEMBANK

      return
      end

************************************************************************
*                                                                      *
      function f_ionp(iz,im,jj,iom,ierr)
*                                                                      *
*        Computes mean ionization potentials in ev.                    *
*        See SPAR Manual, p. 5.                                        *
*                                                                      *
*        im = 1 Original SPAR model.                                   *
*        im = 2 From Berger and Selzer
*                                                                      *
*        jj=1 gaseous element                                          *
*        jj=2 condensed element                                        *
*        jj=3 gaseous mixture                                          *
*        jj=4 condensed mixture                                        *
*                                                                      *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter (p0=0.d0, p2=2.0d0, p3=3.0d0)
      parameter (a1=9.76d0, a2=58.8d0, a3=-0.19d0)
      parameter (b1=10.3d0, b2=1.0d0, b3=0.793d0, b4=p2/p3)
      dimension pv(13), pm(10:100), ph(9,4), ff(4)

*-----------------------------------------------------------------------
*     mean ionization potentials.
*-----------------------------------------------------------------------

      data pv /18.7d0,42.d0,39.d0,60.d0,68.d0,78.d0,99.5d0,98.5d0,
     & 117.d0,140.d0,150.d0,157.d0,163.d0/
      data (pm(i),i=10,100) /1.37d+02,1.49d+02,1.56d+02,1.66d+02,
     & 1.73d+02,1.73d+02,1.80d+02,1.5929d+02,1.88d+02,1.90d+02,1.91d+02,
     & 2.16d+02,2.33d+02,2.45d+02,2.57d+02,2.72d+02,2.86d+02,2.97d+02,
     & 3.11d+02,3.22d+02,3.30d+02,3.34d+02,3.50d+02,3.47d+02,3.48d+02,
     & 3.57d+02,3.52d+02,3.63d+02,3.66d+02,3.79d+02,3.93d+02,4.17d+02,
     & 4.24d+02,4.28d+02,4.41d+02,4.49d+02,4.70d+02,4.70d+02,4.69d+02,
     & 4.88d+02,4.88d+02,4.87d+02,4.85d+02,4.91d+02,4.82d+02,4.88d+02,
     & 4.91d+02,5.01d+02,5.23d+02,5.35d+02,5.46d+02,5.60d+02,5.74d+02,
     & 5.80d+02,5.91d+02,6.14d+02,6.28d+02,6.50d+02,6.58d+02,6.74d+02,
     & 6.84d+02,6.94d+02,7.05d+02,7.18d+02,7.27d+02,7.36d+02,7.46d+02,
     & 7.57d+02,7.90d+02,7.90d+02,8.00d+02,8.10d+02,8.23d+02,8.23d+02,
     & 8.30d+02,8.25d+02,7.94d+02,8.27d+02,8.26d+02,8.41d+02,8.47d+02,
     & 8.78d+02,8.90d+02,9.02d+02,9.21d+02,9.34d+02,9.39d+02,9.52d+02,
     & 9.66d+02,9.80d+02,9.94d+02/
      data ph /1.92d+01,4.18d+01,3.40d+01,3.86d+01,4.90d+01,7.00d+01,
     & 8.20d+01,9.50d+01,1.15d+02,2.18d+01,4.18d+01,4.00d+01,6.37d+01,
     & 7.60d+01,7.80d+01,9.09d+01,1.05d+02,1.12d+02,1.92d+01,4.18d+01,
     & 3.40d+01,3.86d+01,4.90d+01,7.00d+01,8.20d+01,9.70d+01,1.15d+02,
     & 1.92d+01,4.18d+01,4.52d+01,7.20d+01,8.59d+01,8.10d+01,8.20d+01,
     & 1.06d+02,1.12d+02/

*-----------------------------------------------------------------------

         ff(1) = 1.
         ff(2) = 1.
         ff(3) = 1.13
         ff(4) = 1.13

         z = iz

      if(im.eq.1) then
         if(iz.gt.13) then
            p = a1*z+a2*z**a3
         else
            p = pv(iz)
         end if
      else if(im.eq.2) then
         if(iz.le.9) then
            p = ph(iz,jj)
         else if(iz.le.100) then
            p = ff(jj)*pm(iz)
         else
            p = b1*z*(b2-b3/z**b4)
         endif
      else
         p = p0

         write(iom,'(/''Error : illegal ionization energy model'',
     &                '' designator: '',i10)') im
         ierr = 1
         return

      end if

         f_ionp = p

*-----------------------------------------------------------------------

      return
      end
