************************************************************************
*                                                                      *
      subroutine sctgam(eein,wgti,ireg,imat)
*                                                                      *
*        calculate a collision of a photon with an atom.               *
*        last modified by K.Niita on 2012/12/27                        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
     &                      ,iaevt, ibevt, jaevt, jbevt
     &                      ,aevts,aevtr,bevts,bevtr
C for REDUCTION_COUNTER
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      parameter ( rpmass = 938.27d0, rnmass = 939.58d0 )

*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn),  qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
C MATSUDA 2022.12.12  for T-Point
      common /clusttp/ nixcos, ixcoss(nnn)
!$OMP THREADPRIVATE(/clusttp/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)



      common /kmat1g/ kmat(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

*-----------------------------------------------------------------------

      common /wparm/  swtm(20), wc1(20), wc2(20) !FURUTA
      common /wparm0/ wc01(20), wc02(20)         !FURUTA
!$OMP THREADPRIVATE(/wparm0/)


      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)

*-----------------------------------------------------------------------
      common /adjoint/ iadjnt

*-----------------------------------------------------------------------

      dimension usave(3), suu(3,2)
      dimension xv(3)

*-----------------------------------------------------------------------
*     Adjoint mode :  cKN 2015/08/06
*-----------------------------------------------------------------------

         if( iadjnt .eq. 1 ) then  ! T.Sato 2024/01/07, change .eq.0 to .ne.1 because iadjnt = 2 is introduced

            call sctgam_a(eein,wgti,ireg,imat)

            return

         end if

*-----------------------------------------------------------------------

            nclst = -1
            nclsts = 0
            nixcos = 0  ! MATSUDA 2022.12.12  for T-Point

*-----------------------------------------------------------------------

            if( mcal .ne. 0 ) return

            nter = 0
            ipt  = 2
            idx  = 0

            erg  = eein
            eg0  = eein
            icl  = ireg
            wgt  = wgti
            wg0  = wgti
            wga  = 0.0d0
            tme  = 0.0d0
            ecmp = 0.0d0

            usave(1) = uuu
            usave(2) = vvv
            usave(3) = www

            uuu  = 0.0
            vvv  = 0.0
            www  = 1.0

            do i = 0, 20
               numpal(i) = 0
               rumpal(i) = 0.0d0
            enddo

            dxd = 0.0
            xv(1) = 0.0
            xv(2) = 0.0
            xv(3) = 0.0

*-----------------------------------------------------------------------
*     save the incoming direction.
*-----------------------------------------------------------------------

            uold(1) = 0.0
            uold(2) = 0.0
            uold(3) = 1.0

            ntyn = 0
            jsu  = 0

            mk   = imat
            m    = jmd(1+mk)
            m1   = m

*-----------------------------------------------------------------------
*     photonuclear reaction
*-----------------------------------------------------------------------




*-----------------------------------------------------------------------
*     sample cumulative cross section to find the collision nuclide.
*-----------------------------------------------------------------------

         if( npq(mk) .eq. 1 ) goto 20

         t1 = rang() * totm
         c = 0.

         do m = m1, jmd(1+mk+1) - 2

            c = c + rtc(4,lme(2,m)) * fme(m)
            if( t1 .lt. c ) goto 20

         end do

   20    iex = lme(2,m)
         iexp = m
         mtp = -2
         mpan = ipan(icl) + m - m1

         do i = 1, 4

            tpd(i) = rtc(i,iex)

         end do

         tpd(7) = rtc(7,iex)

*-----------------------------------------------------------------------
*     if incident energy greater than emcf, use simple physics,
*     calling elcgam for electrons and bremsstrahlung for implicit
*     capture and skipping coherent scatter and fluorescence.
*-----------------------------------------------------------------------

      call cputime(10)

      if( erg .gt. emcf(2) ) then

*-----------------------------------------------------------------------

               t2 = tpd(3) - tpd(2)
               t1 = wgt * t2 / tpd(4)
               e  = erg

               sw  = wgt
               wgt = t1

               if( ides .eq. 0 )
     &             call elcgam(1,e,ireg,imat,xv,dxd)

               wgt = sw

*-----------------------------------------------------------------------
*        sample analog capture if required.
*-----------------------------------------------------------------------

         if( wc01(14) .eq. 0. ) then

            if( t2 .ge. rang()*tpd(4) ) then
                     aevts(iaevt+28,ireg) =
     &               aevts(iaevt+28,ireg) + wgt
                     bevts(ibevt+28,imat) =
     &               bevts(ibevt+28,imat) + wgt

                     aevts(iaevt+29,ireg) =
     &               aevts(iaevt+29,ireg) + wgt * erg
                     bevts(ibevt+29,imat) =
     &               bevts(ibevt+29,imat) + wgt * erg
                  wgt = 0.0d0

                  goto 910

            end if

*-----------------------------------------------------------------------
*        otherwise simulate capture by weight reduction.
*-----------------------------------------------------------------------

         else
                     aevts(iaevt+28,ireg) =
     &               aevts(iaevt+28,ireg) + t1
                     bevts(ibevt+28,imat) =
     &               bevts(ibevt+28,imat) + t1

                     aevts(iaevt+29,ireg) =
     &               aevts(iaevt+29,ireg) + t1 * erg
                     bevts(ibevt+29,imat) =
     &               bevts(ibevt+29,imat) + t1 * erg
                  wgt = wgt - t1

         end if

*-----------------------------------------------------------------------
*        terminate or continue the photon according to weight cutoffs.
*        sample the type of collision.
*-----------------------------------------------------------------------

            t2 = ( tpd(4) - t2 ) * rang()

            if( t2 .ge. tpd(2) ) goto 330

            s = erg / gpt(3)

*-----------------------------------------------------------------------
*     detailed treatment as in mcp.
*        erg < emcf(2)
*-----------------------------------------------------------------------

      else

*-----------------------------------------------------------------------

            t2 = tpd(4) * rang()

            if( t2 .ge. tpd(2) ) goto 250

            s = erg / gpt(3)

            if( t2 .ge. tpd(1) ) goto 180

      end if

*-----------------------------------------------------------------------
*     incoherent scattering
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        sample scattering angle, cs, from the klein-nishina formula.
*-----------------------------------------------------------------------

            it = 0

  100       call klnishi(s,t4,cs)

*-----------------------------------------------------------------------
*        if energy .gt. emcf treat as mcg and bypass form factors.
*-----------------------------------------------------------------------

         if( erg .gt. emcf(2) ) goto 170

*-----------------------------------------------------------------------
*        calculate form factors for rejection.
*-----------------------------------------------------------------------

            t3 = s*sqrt(1.-cs)
            if(nxs(2,iex).gt.1) goto 110
            t3 = 96.9014*t3
            t6 = 1.-1./(1.+t3**2)**4
            if(it.eq.0) t8 = 1.-1./(1.+(fscon*s)**2)**4
            goto 160

  110       t3 = 29.1445*t3
            if(t3.ge.vic(minc)) goto 170

         do 120 i1 = 2, minc
  120       if(t3.lt.vic(i1)) goto 130

  130       ib = jxs(2,iex)+i1-1
            t6 = xss(ib)+(xss(ib-1)-xss(ib))
     &         *(t3-vic(i1))/(vic(i1-1)-vic(i1))
            if(it.ne.0) goto 160
            t3 = 41.2166*s
            t8 = nxs(2,iex)
            if(t3.ge.vic(minc)) goto 160
         do 140 i = i1, minc
  140       if(t3.lt.vic(i)) goto 150
  150       ib = jxs(2,iex)+i-1
            t8 = xss(ib)+(xss(ib-1)-xss(ib))
     &         *(t3-vic(i))/(vic(i-1)-vic(i))
  160       t3 = t8*rang()
            it = 1
            if(t3.ge.t6) goto 100

*-----------------------------------------------------------------------
*           compton energy loss
*-----------------------------------------------------------------------

  170    continue

            erg = gpt(3)*t4

                     aevts(iaevt+62,ireg) =
     &               aevts(iaevt+62,ireg) + wgt * ( eg0 - erg )
                     bevts(ibevt+62,imat) =
     &               bevts(ibevt+62,imat) + wgt * ( eg0 - erg )
                     ecmp = ecmp + wgt * ( eg0 - erg )

            cs = min(max(-one,cs),one)

            call dtcos(cs,uold,uuu,0,irdm)

            atmrc(1,4) = atmrc(1,4) + wgt / wg0   ! S.Abe 2017/11/21

            if( ides .eq. 0 )
     &            call elcgam(2,cs,ireg,imat,xv,dxd)

            ntyn = 1
            mtp = -1
            ipsc = 7

            goto 900

*-----------------------------------------------------------------------
*     coherent scattering
*-----------------------------------------------------------------------

  180       t5 = 1698.8038 * s**2

*-----------------------------------------------------------------------
*        sample scattering angle.
*-----------------------------------------------------------------------

            if(t5.lt.wco(mcoh)) goto 190
            t7 = xss(jxs(3,iex)+mcoh-1)
            goto 220
  190    do 200 i = 2, mcoh
  200       if(t5.lt.wco(i)) goto 210

  210       t3 = (t5-wco(i))/(wco(i-1)-wco(i))
            ib = jxs(3,iex)+i-1
            t7 = xss(ib)+t3*(xss(ib-1)-xss(ib))
  220       t3 = t7*rang()
            ib = jxs(3,iex)-1
         do 230 i = 2, mcoh
  230       if(t3.lt.xss(ib+i)) goto 240

  240       t3 = (t3-xss(ib+i))/(xss(ib+i-1)-xss(ib+i))
            cs = 1.-2.*(wco(i)+t3*(wco(i-1)-wco(i)))/t5
            t3 = 1.+cs**2
            if(t3.le.2.*rang()) goto 220
            ipsc = 6
            call dtcos(cs,uold,uuu,0,irdm)
            ntyn = 2

            atmrc(1,6) = atmrc(1,6) + wgt / wg0   ! S.Abe 2017/11/21

            goto 900

*-----------------------------------------------------------------------
*     capture with fluorescence
*-----------------------------------------------------------------------

  250    if( t2 .ge. tpd(3) ) goto 330
                     aevts(iaevt+28,ireg) =
     &               aevts(iaevt+28,ireg) + wgt
                     bevts(ibevt+28,imat) =
     &               bevts(ibevt+28,imat) + wgt

                     aevts(iaevt+29,ireg) =
     &               aevts(iaevt+29,ireg) + wgt * erg
                     bevts(ibevt+29,imat) =
     &               bevts(ibevt+29,imat) + wgt * erg
            er = erg
            ie = iex
            ik = ipt

            call flugam(er,ie,mk,ik,ireg,imat,xv,dxd)

*-----------------------------------------------------------------------

         if( nclsts .gt. 0 ) then

               eout = 0.0

            do k = 1, nclsts

c T.Sato 2014/2/21 bug for photon transport using weight control
               eout = eout + qclusts(7,k) * qclusts(8,k)

            end do

            if( eout .lt. eg0 ) then

                     ers = eg0 - eout

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 7
                     nixcos = nixcos + 1  ! MATSUDA 2022.12.12  for T-Point

                     jclusts(0,nclsts) = 0
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(3,nclsts) = 12

cKN 2012/12/06 : dead electron !!!
                     jclusts(4,nclsts) = -11

                     jclusts(5,nclsts) = -1
                     jclusts(6,nclsts) = 0
                     jclusts(7,nclsts) = 11
                     jclusts(8,nclsts) = 0

                     rms = gpt(3)
                     rmg = rms / 1000.
                     pr  = sqrt( ers * ( ers + 2.0 * rms ) ) / 1000.
                     ett = sqrt( pr**2 + rmg**2 )

                     pxl = 0.
                     pyl = 0.
                     pzl = pr

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = pxl
                     qclusts(2,nclsts)  = pyl
                     qclusts(3,nclsts)  = pzl
                     qclusts(4,nclsts)  = ett
                     qclusts(5,nclsts)  = rmg
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = ers
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd
                     ixcoss(nixcos) = 0  ! MATSUDA 2022.12.12  for T-Point

            end if

         end if

*-----------------------------------------------------------------------

            goto 910

*-----------------------------------------------------------------------
*     pair production
*-----------------------------------------------------------------------

  330    continue
                     aevts(iaevt+36,ireg) =
     &               aevts(iaevt+36,ireg) + wgt
                     bevts(ibevt+36,imat) =
     &               bevts(ibevt+36,imat) + wgt

                     aevts(iaevt+37,ireg) =
     &               aevts(iaevt+37,ireg) + wgt * eein
                     bevts(ibevt+37,imat) =
     &               bevts(ibevt+37,imat) + wgt * eein
            ntyn =  5
            mtp  = -4

            if( ides .eq. 0 ) then

               cs = 0.
               call elcgam(3,cs,ireg,imat,xv,dxd)

               if( cs .eq. 0. ) goto 910

            end if

*-----------------------------------------------------------------------
*     p-annihilation product two photons
*-----------------------------------------------------------------------

            if( gpt(3) .ge. elc(2) ) then

                  erg  = gpt(3)
                  ipsc = 105

                  call isos(uuu,0)

                     suu(1,1) = uuu
                     suu(2,1) = vvv
                     suu(3,1) = www

                     suu(1,2) = -uuu
                     suu(2,2) = -vvv
                     suu(3,2) = -www

*-----------------------------------------------------------------------
*              photon heating flag : jclusts(4,nclsts) = 1
*-----------------------------------------------------------------------

                  do k = 1, 2

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 4
                     nixcos = nixcos + 1  ! MATSUDA 2022.12.12  for T-Point

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(3,nclsts) = 14
                     jclusts(4,nclsts) = 1
                     jclusts(5,nclsts) = 0
                     jclusts(6,nclsts) = 0
                     jclusts(7,nclsts) = 22
                     jclusts(8,nclsts) = 0

                     pr   = erg / 1000.0

                     pxl = pr * suu(1,k)
                     pyl = pr * suu(2,k)
                     pzl = pr * suu(3,k)

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
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd
                     ixcoss(nixcos) = 0  ! MATSUDA 2022.12.12  for T-Point
                     aevts(iaevt+30,ireg) =
     &               aevts(iaevt+30,ireg) + wgt
                     bevts(ibevt+30,imat) =
     &               bevts(ibevt+30,imat) + wgt

                     aevts(iaevt+31,ireg) =
     &               aevts(iaevt+31,ireg) + wgt * erg
                     bevts(ibevt+31,imat) =
     &               bevts(ibevt+31,imat) + wgt * erg
                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

                  end do

            end if

                     goto 910

*-----------------------------------------------------------------------
*     out going photon from coherent and incoherent
*     photon heating flag : jclusts(4,nclsts) = 1
*-----------------------------------------------------------------------

  900 continue

                     nclsts = nclsts + 1

                     iclusts(nclsts) = 4
                     nixcos = nixcos + 1  ! MATSUDA 2022.12.12  for T-Point

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(3,nclsts) = 14
                     jclusts(4,nclsts) = 1
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
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd
                     ixcoss(nixcos) = 0  ! MATSUDA 2022.12.12  for T-Point

*-----------------------------------------------------------------------
*     finish
*-----------------------------------------------------------------------

  910 continue

            kdecay(4) = 0

            uuu = usave(1)
            vvv = usave(2)
            www = usave(3)

            rncnt(10) = rncnt(10) + 1.0
            call cputime(10)

*-----------------------------------------------------------------------

  999 continue

      return
      end

************************************************************************
*                                                                      *
      subroutine sctgam_a(eein,wgti,ireg,imat)
*                                                                      *
*        Adjoint mode :                                                *
*        calculate a collision of a photon with an atom.               *
*        last modified by K.Niita on 2015/08/05                        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
     &                      ,iaevt, ibevt, jaevt, jbevt
     &                      ,aevts,aevtr,bevts,bevtr
C for REDUCTION_COUNTER
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2

      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      parameter ( rpmass = 938.27, rnmass = 939.58 )

*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn),  qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

C for USE_MOD_COUNTER


      common /kmat1g/ kmat(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

*-----------------------------------------------------------------------

      common /wparm/  swtm(20), wc1(20), wc2(20) !FURUTA
      common /wparm0/ wc01(20), wc02(20)         !FURUTA
!$OMP THREADPRIVATE(/wparm0/)


*-----------------------------------------------------------------------
      common /adjoint/ iadjnt
      common /adjenrg/ adjemax

*-----------------------------------------------------------------------

      dimension usave(3), suu(3,2)
      dimension xv(3)

*-----------------------------------------------------------------------

            nclst = -1
            nclsts = 0

*-----------------------------------------------------------------------

            if( mcal .ne. 0 ) return

            nter = 0
            ipt  = 2
            idx  = 0

            erg  = eein
            eg0  = eein
            icl  = ireg
            wgt  = wgti
            wg0  = wgti
            wga  = 0.0d0
            tme  = 0.0d0
            ecmp = 0.0d0

            usave(1) = uuu
            usave(2) = vvv
            usave(3) = www

            uuu  = 0.0
            vvv  = 0.0
            www  = 1.0

            do i = 0, 20
               numpal(i) = 0
               rumpal(i) = 0.0d0
            enddo

            dxd = 0.0
            xv(1) = 0.0
            xv(2) = 0.0
            xv(3) = 0.0

*-----------------------------------------------------------------------
*     save the incoming direction.
*-----------------------------------------------------------------------

            uold(1) = 0.0
            uold(2) = 0.0
            uold(3) = 1.0

            ntyn = 0
            jsu  = 0

            mk   = imat
            m    = jmd(1+mk)
            m1   = m
	
*-----------------------------------------------------------------------
*        implicit treatment for photoelectric
cAlex 2020/06/30  this needs to be before selection of collision nuclide
*-----------------------------------------------------------------------

            t1 = wgt * ( siga ) / ( totm )
                     aevts(iaevt+28,ireg) =
     &               aevts(iaevt+28,ireg) + t1
                     bevts(ibevt+28,imat) =
     &               bevts(ibevt+28,imat) + t1

                     aevts(iaevt+29,ireg) =
     &               aevts(iaevt+29,ireg) + t1 * erg
                     bevts(ibevt+29,imat) =
     &               bevts(ibevt+29,imat) + t1 * erg

                  wgt = t1

*-----------------------------------------------------------------------
*     sample cumulative cross section to find the collision nuclide.
*-----------------------------------------------------------------------

         if( npq(mk) .eq. 1 ) goto 20

         t1 = rang() * siga
         c = 0.

         do m = m1, jmd(1+mk+1) - 2

            c = c + (rtc(8,lme(2,m)) + rtc(9,lme(2,m))
     &            + rtc(2,lme(2,m))
     &            - rtc(1,lme(2,m))) * fme(m)
            if( t1 .lt. c ) goto 20

         end do

   20    iex = lme(2,m)
         iexp = m
         mtp = -2
         mpan = ipan(icl) + m - m1

*-----------------------------------------------------------------------
*     if incident energy greater than emcf, use simple physics,
*     calling elcgam for electrons and bremsstrahlung for implicit
*     capture and skipping coherent scatter and fluorescence.
*-----------------------------------------------------------------------

      call cputime(10)

*-----------------------------------------------------------------------
cAlex 2020/06/30  only need nuclide adjoint XS from now on
*        tpd8   : adjoint incoherent
*        tpd9   : adjoint incoherent + adjoint pp
*        tpd10  : adjoint incoherent + adjoint pp + adjoint (forward) coherent
*-----------------------------------------------------------------------

            tpd8  = rtc(8,iex)
            tpd9  = tpd8 + rtc(9,iex)
            tpd10 = tpd9 + ( rtc(2,iex) - rtc(1,iex) )

*-----------------------------------------------------------------------
*        sample the type of collision.
*-----------------------------------------------------------------------
cAlex 2020/06/30  sample from three possible scattering collisions
            t2 = tpd10 * rang()
            s = erg / gpt(3)


            if( t2 .lt. tpd8 )  goto 90   ! incoherent
            if( t2 .lt. tpd9 )  goto 330  ! pair production
            if( t2 .lt. tpd10 ) goto 180  ! coherent

*-----------------------------------------------------------------------

            goto 900

*-----------------------------------------------------------------------
*     incoherent scattering
*-----------------------------------------------------------------------

   90    continue

*-----------------------------------------------------------------------
*        sample scattering angle, cs, from the klein-nishina formula.
*-----------------------------------------------------------------------

        emax = adjemax
        rkp  = erg / gpt(3)
        efn  = erg

cAlex 2020/06/30 set zmin depending on whether E' < E_m / (1+2E_m/E_e)'
        alx = emax / ( 1 + 2.0 * emax / gpt(3) )
        if( efn .gt. alx ) then
            zmin = efn / emax
        else
            zmin = 1 - 2.0 * rkp
        end if

*-----------------------------------------------------------------------

            it = 0
            itt = 0

  100       continue

            itt = itt + 1
            if( itt .gt. 10000 ) goto 111

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

cAlex 2020/06/30 corrected adjoint KN rejection sampling
            zv = rang() * ( 1.0 - zmin ) + zmin
            rkn = ( 1.0 / zv + zv - 1.0
     &          + ( 1.0 + zv / rkp - 1.0 / rkp )**2 ) /
     &           (1.0/zmin+zmin-1.0+(1.0+zmin/rkp-1.0/rkp)**2)

            if( rkn .lt. rang() ) goto 100

  111    continue

            ein = efn / zv
            s   = ein / gpt(3)
            t4  = rkp
            cs  = 1. + zv / rkp - 1. / rkp

*-----------------------------------------------------------------------
*        if energy .gt. emcf treat as mcg and bypass form factors.
*-----------------------------------------------------------------------

         if( ein .gt. emcf(2) ) goto 170
         if( ein .gt. emax ) goto 910

*-----------------------------------------------------------------------
*        calculate form factors for rejection.
*-----------------------------------------------------------------------

            t3 = s*sqrt(1.-cs)
            if(nxs(2,iex).gt.1) goto 110
            t3 = 96.9014*t3
            t6 = 1.-1./(1.+t3**2)**4
            if(it.eq.0) t8 = 1.-1./(1.+(fscon*s)**2)**4
            goto 160

  110       t3 = 29.1445*t3
            if(t3.ge.vic(minc)) goto 170

         do 120 i1 = 2, minc
  120       if(t3.lt.vic(i1)) goto 130

  130       ib = jxs(2,iex)+i1-1
            t6 = xss(ib)+(xss(ib-1)-xss(ib))
     &         *(t3-vic(i1))/(vic(i1-1)-vic(i1))
            if(it.ne.0) goto 160
cAlex 2020/06/30 corrected adjoint KN rejection sampling
            if( efn .gt. alx ) then
                t3 = 29.1445*s*sqrt((1.-zmin)/rkp)
            else
                t3 = 41.2166*s
            end if
            t8 = nxs(2,iex)
            if(t3.ge.vic(minc)) goto 160
         do 140 i = i1, minc
  140       if(t3.lt.vic(i)) goto 150
  150       ib = jxs(2,iex)+i-1
            t8 = xss(ib)+(xss(ib-1)-xss(ib))
     &         *(t3-vic(i))/(vic(i-1)-vic(i))
  160       t3 = t8*rang()
            it = 1

            if( itt .gt. 10000 ) goto 170
            if(t3.ge.t6) goto 100

*-----------------------------------------------------------------------
*           compton energy loss
*-----------------------------------------------------------------------

  170    continue

            erg = ein

                     aevts(iaevt+62,ireg) =
     &               aevts(iaevt+62,ireg) + wgt * ( eg0 - erg )
                     bevts(ibevt+62,imat) =
     &               bevts(ibevt+62,imat) + wgt * ( eg0 - erg )
                     ecmp = ecmp + wgt * ( eg0 - erg )
cAlex 2020/06/30 corrected adjoint KN rejection sampling
            if( efn .gt. alx ) then
                cs = min(max(1.+zmin/rkp-1/rkp,cs),one)
            else
                cs = min(max(-one,cs),one)
            end if

            call dtcos(cs,uold,uuu,0,irdm)

cAlex 2020/06/30 remove legacy code relating to relaxation electrons

            ntyn = 1
            mtp = -1
            ipsc = 17

            goto 900

*-----------------------------------------------------------------------
*     coherent scattering
*-----------------------------------------------------------------------

  180       t5 = 1698.8038 * s**2

*-----------------------------------------------------------------------
*        sample scattering angle.
*-----------------------------------------------------------------------

            if(t5.lt.wco(mcoh)) goto 190
            t7 = xss(jxs(3,iex)+mcoh-1)
            goto 220
  190    do 200 i = 2, mcoh
  200       if(t5.lt.wco(i)) goto 210

  210       t3 = (t5-wco(i))/(wco(i-1)-wco(i))
            ib = jxs(3,iex)+i-1
            t7 = xss(ib)+t3*(xss(ib-1)-xss(ib))
  220       t3 = t7*rang()
            ib = jxs(3,iex)-1
         do 230 i = 2, mcoh
  230       if(t3.lt.xss(ib+i)) goto 240

  240       t3 = (t3-xss(ib+i))/(xss(ib+i-1)-xss(ib+i))
            cs = 1.-2.*(wco(i)+t3*(wco(i-1)-wco(i)))/t5
            t3 = 1.+cs**2
            if(t3.le.2.*rang()) goto 220
            ipsc = 6
            call dtcos(cs,uold,uuu,0,irdm)
            ntyn = 2

            goto 900

*-----------------------------------------------------------------------
*     pair production
*-----------------------------------------------------------------------

  330    continue
         goto 910

                     aevts(iaevt+36,ireg) =
     &               aevts(iaevt+36,ireg) + wgt
                     bevts(ibevt+36,imat) =
     &               bevts(ibevt+36,imat) + wgt

                     aevts(iaevt+37,ireg) =
     &               aevts(iaevt+37,ireg) + wgt * eein
                     bevts(ibevt+37,imat) =
     &               bevts(ibevt+37,imat) + wgt * eein
            ntyn =  5
            mtp  = -4

            if( ides .eq. 0 ) then

               cs = 0.
               call elcgam(3,cs,ireg,imat,xv,dxd)

               if( cs .eq. 0. ) goto 910

            end if

*-----------------------------------------------------------------------
*     p-annihilation product two photons
*-----------------------------------------------------------------------

            if( gpt(3) .ge. elc(2) ) then

                  erg  = gpt(3)
                  ipsc = 105

                  call isos(uuu,0)

                     suu(1,1) = uuu
                     suu(2,1) = vvv
                     suu(3,1) = www

                     suu(1,2) = -uuu
                     suu(2,2) = -vvv
                     suu(3,2) = -www

*-----------------------------------------------------------------------
*              photon heating flag : jclusts(4,nclsts) = 1
*-----------------------------------------------------------------------

                  do k = 1, 2

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 4

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(3,nclsts) = 14
                     jclusts(4,nclsts) = 1
                     jclusts(5,nclsts) = 0
                     jclusts(6,nclsts) = 0
                     jclusts(7,nclsts) = 22
                     jclusts(8,nclsts) = 0

                     pr   = erg / 1000.0

                     pxl = pr * suu(1,k)
                     pyl = pr * suu(2,k)
                     pzl = pr * suu(3,k)

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
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd
                     aevts(iaevt+30,ireg) =
     &               aevts(iaevt+30,ireg) + wgt
                     bevts(ibevt+30,imat) =
     &               bevts(ibevt+30,imat) + wgt

                     aevts(iaevt+31,ireg) =
     &               aevts(iaevt+31,ireg) + wgt * erg
                     bevts(ibevt+31,imat) =
     &               bevts(ibevt+31,imat) + wgt * erg
                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

                  end do

            end if

                     goto 910

*-----------------------------------------------------------------------
*     out going photon from coherent and incoherent
*     photon heating flag : jclusts(4,nclsts) = 1
*-----------------------------------------------------------------------

  900 continue

                     nclsts = nclsts + 1

                     iclusts(nclsts) = 4

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(3,nclsts) = 14
                     jclusts(4,nclsts) = 1
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
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd

*-----------------------------------------------------------------------
*     finish
*-----------------------------------------------------------------------

  910 continue

            kdecay(4) = 0

            uuu = usave(1)
            vvv = usave(2)
            www = usave(3)

            rncnt(10) = rncnt(10) + 1.0
            call cputime(10)

*-----------------------------------------------------------------------

  999 continue

      return
      end

************************************************************************
cfrtati 2021/12/17 added argument icmm                                 *
*                                                                      *
      subroutine sctpni(eein,wgti,ireg,imat,icmm)
*                                                                      *
*        calculate a collision of a photon with an atom.               *
*        last modified by K.Niita on 2012/12/27                        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
     &                      ,iaevt, ibevt, jaevt, jbevt
     &                      ,aevts,aevtr,bevts,bevtr
C for REDUCTION_COUNTER
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2
      use GGMBANKMOD !FURUTA
      use NGSDATAMOD, only : bindeg
      use GGMARRAYMOD !2020ASTOM
      use moddas_material

      use udm_Parameter
      use udm_Utility
      use udm_Manager

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      parameter ( rpmass = 938.27, rnmass = 939.58 )

*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn),  qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /pnint/  ipnint
      common /pionflag/ iqdflag, ipionflag, istrflag, iphreturn
!$OMP THREADPRIVATE(/pionflag/)

C for USE_MOD_COUNTER


      common /kmat1g/ kmat(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

*-----------------------------------------------------------------------

      common /wparm/  swtm(20), wc1(20), wc2(20) !FURUTA
      common /wparm0/ wc01(20), wc02(20)         !FURUTA
!$OMP THREADPRIVATE(/wparm0/)

*-----------------------------------------------------------------------

      common /qmdflg/ irqmd
      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)
      common /tarzmn/ itarz, itarm, itarn
!$OMP THREADPRIVATE(/tarzmn/)

C for USE_MOD_COUNTER

      common /pnixs/ totmpni,pnixsm(pnlmax),ipnikind(pnlmax),egypni,
     &               dlibmax
!$OMP THREADPRIVATE(/pnixs/)
      common /nrfmem/ spis, spgr, levabs, lflgnrf, mpole
!$OMP THREADPRIVATE(/nrfmem/)
      common /nrfdir/ usave(3)
!$OMP THREADPRIVATE(/nrfdir/)

      common /cascid/ jcasc
!$OMP THREADPRIVATE(/cascid/)

*-----------------------------------------------------------------------

      dimension suu(3,2)
      dimension xv(3)

*-----------------------------------------------------------------------
      integer igdrqd ! receiver of junction "jgdrqd"
      integer jgdrqd ! function, judging GDR or QD reactions

*-----------------------------------------------------------------------

            nclst = -1
            if( imuinthit .ne. 1 ) then
               nclsts = 0
            endif


*-----------------------------------------------------------------------

            if( mcal .ne. 0 ) return

            nter = 0
            ipt  = 2
            idx  = 0

            erg  = eein
            eg0  = eein
            icl  = ireg
            wgt  = wgti
            wg0  = wgti
            wga  = 0.0d0
            tme  = 0.0d0
            ecmp = 0.0d0

            usave(1) = uuu
            usave(2) = vvv
            usave(3) = www

            uuu  = 0.0
            vvv  = 0.0
            www  = 1.0

            if( imuinthit .ne. 1 ) then
               do i = 0, 20
                  numpal(i) = 0
                  rumpal(i) = 0.0d0
               enddo
            endif
            dxd = 0.0
            xv(1) = 0.0
            xv(2) = 0.0
            xv(3) = 0.0

*-----------------------------------------------------------------------
*     save the incoming direction.
*-----------------------------------------------------------------------

            uold(1) = 0.0
            uold(2) = 0.0
            uold(3) = 1.0

            ntyn = 0
            jsu  = 0

            mk   = imat
            m    = jmd(1+mk)
            m1   = m


            icntm = 1   ! S.Abe 2015/12/08, corrected

*-----------------------------------------------------------------------
*     photonuclear reaction
*-----------------------------------------------------------------------




*-----------------------------------------------------------------------
*     sample cumulative cross section to find the collision nuclide.
*-----------------------------------------------------------------------

      if( imuinthit .ne. 1 ) then
         if( npq(mk) .eq. 1 ) goto 20
       if( icmm.ne.0 ) then
          m = m1 + icmm - 1
       else

         t1 = rang() * totmpni
         c = 0.

         do m = m1, jmd(1+mk+1) - 2
ccABE bugfix @2014/11/11
            icntm = m - m1 + 1
            c = c + pnixsm(icntm) * fme(m)

            if( t1 .lt. c ) goto 20

         end do

       end if ! frtati 2021/12/17

   20    iex = lme(2,m)
         iexp = m
         mtp = -2
         mpan = ipan(icl) + m - m1

         do i = 1, 4

            tpd(i) = rtc(i,iex)

         end do

ccABE bugfix @2014/11/11
         tpd(7) = pnixsm(icntm)

      endif

*-----------------------------------------------------------------------
*     photonuclear GDR part
*-----------------------------------------------------------------------

      if( ipnint .ne. 0 .or. imuinthit .eq. 1 ) then

*-----------------------------------------------------------------------
*           mother nmucleus
*-----------------------------------------------------------------------
         if( imuinthit .eq. 1 ) then
            izm = itarz
            ims = itarm
            inm = itarn
         else
            izm = iza(iexp) / 1000
            ims = iza(iexp) - 1000 * izm

            if( ims .eq. 0 ) then

               lemm  = nint( dnel_das(kmat0+imat) )

               sek = 0.0
               do lem = 1, lemm
                  izt = zz_das(kmat(imat)+lem)
                  if( izt .eq. izm )
     &               sek = sek + den_das(kmat(imat)+lem)
               end do

               sekd = sek * rang()

               sek = 0.0
               do lem = 1, lemm
                  izt = zz_das(kmat(imat)+lem)
                  if( izt .eq. izm ) then

                     sek = sek + den_das(kmat(imat)+lem)
                     if( sek .ge. sekd ) goto 55

                  end if

               end do

   55          continue

               ims = nint( a_das(kmat(imat)+lem) )

            end if

            inm = ims - izm

         endif

cABE add @2024/04/26
         mathz = izm
         mathn = inm
cabe end

*-----------------------------------------------------------------------

            call cputime(22)
            if( imuinthit .eq. 1 ) then
             call jgdrqd_vpi(izm,inm,erg,igdrqd) ! judging function, GDR or QD or PD?
            else
             igdrqd = jgdrqd(izm,inm,erg) ! judging function, GDR or QD or PD?
            endif

            ipim = 0    ! S.Abe 2015/08/10

*-----------------------------------------------------------------------
*     user defined photonuclear reaction (call user_defined_interaction(-22,i))
*-----------------------------------------------------------------------
            if(udm_int_num .gt. 0) then
              mat_Z(1)=izm
              mat_A(1)=izm+inm
              do i = 1, udm_int_num
                ! generate final states ------------
                udm_Kin=eein
                udm_kf_incident=22
                udm_logical=.false.
                call user_defined_interaction(-22,i) ! special action: -22
                ! ----------------------------------
                if(udm_logical) goto 999 ! final state is generated
              enddo
            endif

*-----------------------------------------------------------------------
*     photonuclear GDR reaction or NRF
*-----------------------------------------------------------------------

            if ( igdrqd .ge. 1 ) then ! GDR or NRF reaction occurs  2014/8/25 ogawa modified

             if(igdrqd .eq. 2) lflgnrf = 1 ! NRF reaction

cABE relocate 2024/04/26 
c               mathz = izm
c               mathn = inm
cABE end
               be = bindeg(mathz,mathn)
               tams = ( mathz * rpmass + mathn * rnmass - be )

               tapm = sqrt(tams**2 + 2.0d0 * tams * erg)

               plmx = 0.0
               plmy = 0.0
               plmz = erg

               pabs = plmx**2 + plmy**2 + plmz**2
               etot = sqrt( pabs + tapm**2 )
               elab = etot - tapm

               nclst = 1
               iclust(nclst) = 0

               jclust(0,nclst) = 0
               jclust(1,nclst) = mathz
               jclust(2,nclst) = mathn
               jclust(3,nclst) = 19
               jclust(4,nclst) = 0
               jclust(5,nclst) = mathz
               jclust(6,nclst) = mathz + mathn
               jclust(7,nclst) = mathz * 1000000 + mathz + mathn
               jclust(8,nclst) = 0

               qclust(0,nclst)  = 0.0
               qclust(1,nclst)  = plmx / 1000.
               qclust(2,nclst)  = plmy / 1000.
               qclust(3,nclst)  = plmz / 1000.
               qclust(4,nclst)  = etot / 1000.
               qclust(5,nclst)  = tapm / 1000.
               qclust(6,nclst)  = tapm - tams
               qclust(7,nclst)  = elab
               qclust(8,nclst)  = wgt / wg0
               qclust(9,nclst)  = tme
               qclust(10,nclst) = 0.0d0
               qclust(11,nclst) = 0.0d0
               qclust(12,nclst) = 0.0d0

               call nevap(0)

               if( imuinthit .eq. 1 ) then
                  rnpnt(5) = rnpnt(5) + 1.0
                  rnpntr(5) = rnpntr(5) + 1.0
               else
                  rnpnt(1) = rnpnt(1) + 1.0
                  rnpntr(1) = rnpntr(1) + 1.0
               endif

*-----------------------------------------------------------------------
*     Quasideuteron process
*-----------------------------------------------------------------------

            else if ( igdrqd == 0 ) then
               iqdflag = 1 ! flag of quasi-deuteron process

  102          continue

               jcasc = 4      ! S.Abe 2019/03/13, JQMD flag for therdec

               if (irqmd .eq. 0) then
                  call jqmdin(14,22,erg,ims,izm,bmax0)
               else
                  call jqmdinR(14,22,erg,ims,izm,bmax0)
               endif
               call jqmdchk(erg,ims,izm)
               if( nclst .lt. 0 ) then
                ipim = ipim + 1
                if( ipim .le. 20 .and. nclst .gt. -2 ) goto 102
               endif

               call nevap(0)

               if( imuinthit .eq. 1 ) then
                  rnpnt(6) = rnpnt(6) + 1.0
                  rnpntr(6) = rnpntr(6) + 1.0
               else
                  rnpnt(2) = rnpnt(2) + 1.0
                  rnpntr(2) = rnpntr(2) + 1.0
               endif
               iqdflag = 0

*-----------------------------------------------------------------------
*     Pion produciton process
*-----------------------------------------------------------------------

            else if ( igdrqd == -1 ) then

               ipionflag = 1 ! flag of pion production process

  103          continue

               jcasc = 4      ! S.Abe 2019/03/13, JQMD flag for therdec

               if (irqmd .eq. 0) then
                  call jqmdin(14,22,erg,ims,izm,bmax0)
               else
                  call jqmdinR(14,22,erg,ims,izm,bmax0)
               endif
               call jqmdchk(erg,ims,izm)
               if( nclst .lt. 0 ) then
                ipim = ipim + 1
                if( ipim .le. 20 .and. nclst .gt. -2 ) goto 103
               endif

               call nevap(0)

               if( imuinthit .eq. 1 ) then
                  rnpnt(7) = rnpnt(7) + 1.0
                  rnpntr(7) = rnpntr(7) + 1.0
               else
                  rnpnt(3) = rnpnt(3) + 1.0
                  rnpntr(3) = rnpntr(3) + 1.0
               endif
               ipionflag = 0

*-----------------------------------------------------------------------
*     String produciton process
*-----------------------------------------------------------------------

            elseif (igdrqd .eq. -2) then

               istrflag = 1 ! flag of string production process
               jcasc = 3      ! S.Abe 2019/03/13, JAM flag for therdec

               call jamin(0,14,22,erg,ims,izm,bmax0)

               if( imuinthit .eq. 1 ) then
                  rnpnt(8) = rnpnt(8) + 1.0
               else
                  rnpnt(4) = rnpnt(4) + 1.0
               endif

               if (iphreturn .ne. 1) then

                  call nevap(0)

                  if( imuinthit .eq. 1 ) then
                     rnpntr(8) = rnpntr(8) + 1.0
                  else
                     rnpntr(4) = rnpntr(4) + 1.0
                  endif

               else

                  jcasc = 4      ! S.Abe 2019/03/13, JQMD flag for therdec
                  ipionflag = 1

  104             continue
                  if(irqmd .eq. 0) then
                     call jqmdin( 14, 22, erg, ims, izm, bmax0)
                  else
                     call jqmdinR( 14, 22, erg, ims, izm, bmax0)
                  endif

               call jqmdchk(erg,ims,izm)
                  if( nclst .lt. 0 ) then
                   ipim = ipim + 1
                   if( ipim .le. 20 .and. nclst .gt. -2 ) goto 104
                  endif

                  call nevap(0)

                  if( imuinthit .eq. 1 ) then
                     rnpnt(7) = rnpnt(7) + 1.0
                     rnpntr(7) = rnpntr(7) + 1.0
                  else
                     rnpnt(3) = rnpnt(3) + 1.0
                     rnpntr(3) = rnpntr(3) + 1.0
                  endif

               endif

               istrflag = 0
               iphreturn = 0

*-----------------------------------------------------------------------
            end if


            kdecay(4) = 0

            uuu = usave(1)
            vvv = usave(2)
            www = usave(3)

            if( imuinthit .ne. 1 ) then
               rncnt(22) = rncnt(22) + 1.0
               call cputime(22)
            endif

      end if

*-----------------------------------------------------------------------

  999 continue

      return
      end


************************************************************************
*                                                                      *
      subroutine elcgam(m,cs,ireg,imat,xv,dxd)
*                                                                      *
*        make one or two electrons in a photon collision of type m.    *
*        m=1  for photoelectric effect.                                *
*        m=2  for incoherent (compton) scattering.                     *
*        m=3  for pair production.                                     *
*        m=4  for auger emission from positrons, electrons or photons. *
*        cs = electron energy for m=1,4, cosine for m=2, flag for m=3. *
*                                                                      *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for REDUCTION_COUNTER
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
C MATSUDA 2022.12.12  for T-Point
      common /clusttp/ nixcos, ixcoss(nnn)
!$OMP THREADPRIVATE(/clusttp/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

C for USE_MOD_COUNTER

*-----------------------------------------------------------------------
      common /elpos/  ielpos
!$OMP THREADPRIVATE(/elpos/)
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)

*-----------------------------------------------------------------------

      dimension ee(13),ce(21,13),ep(9),dp(21,9)
      data ee/100.,50.,20.,10.,5.,2.,1.,.5,.2,.1,.05,.02,.01/
      data (ce(i, 1),i=1,21)/    5.0000000e-1,1.8391734e+2,3.7051644e+2,
     1 5.6409171e+2,8.5433494e+2,8.6549644e+2,1.1960445e+3,1.8133120e+3,
     2 3.3154915e+3,1.1560527e+4,7.7358920e+4,7.7358920e+4,7.7358920e+4,
     3 7.7358920e+4,7.7358920e+4,7.7358920e+4,7.7358920e+4,7.7358920e+4,
     4 7.7358920e+4,7.7358920e+4,7.7358920e+4/
      data (ce(i, 2),i=1,21)/    5.0001000e-1,1.5734214e+2,3.1179562e+2,
     1 4.6600627e+2,6.8708044e+2,7.2630997e+2,9.3864458e+2,1.2726402e+3,
     2 1.8585851e+3,3.1036588e+3,7.2354569e+3,1.9539631e+4,1.9539631e+4,
     3 1.9539631e+4,1.9539631e+4,1.9539631e+4,1.9539631e+4,1.9539631e+4,
     4 1.9539631e+4,1.9539631e+4,1.9539631e+4/
      data (ce(i, 3),i=1,21)/    5.0008000e-1,7.3991679e+1,1.4697244e+2,
     1 2.2059103e+2,2.9650433e+2,3.6916714e+2,4.5806222e+2,5.6159689e+2,
     2 5.9186199e+2,6.8870573e+2,8.0858436e+2,9.5808897e+2,1.1453783e+3,
     3 1.3799180e+3,1.6705798e+3,2.0208947e+3,2.4192528e+3,2.8236151e+3,
     4 3.1519086e+3,3.2220501e+3,3.2220501e+3/
      data (ce(i, 4),i=1,21)/    5.0030000e-1,2.5031580e+1,4.8446990e+1,
     1 7.2297589e+1,9.6665639e+1,1.2168910e+2,1.4755735e+2,1.7420693e+2,
     2 2.0162045e+2,2.2944535e+2,2.6029701e+2,2.8916076e+2,3.2332973e+2,
     3 3.6023191e+2,3.8754624e+2,4.2971674e+2,4.8031309e+2,5.4202665e+2,
     4 6.1881714e+2,7.1677975e+2,8.4578033e+2/
      data (ce(i, 5),i=1,21)/    5.0108000e-1,8.5127000e+0,1.5526810e+1,
     1 2.2556960e+1,2.9680480e+1,3.6933560e+1,4.4343450e+1,5.1940279e+1,
     2 5.9746320e+1,6.7800979e+1,7.6153979e+1,8.4843859e+1,9.3893880e+1,
     3 1.0344828e+2,1.1366273e+2,1.2457488e+2,1.3645768e+2,1.4961923e+2,
     4 1.6563791e+2,1.8436704e+2,2.3213951e+2/
      data (ce(i, 6),i=1,21)/    5.0528999e-1,2.9538700e+0,4.7418300e+0,
     1 6.4446200e+0,8.1125700e+0,9.7666099e+0,1.1419390e+1,1.3080480e+1,
     2 1.4758300e+1,1.6461250e+1,1.8198290e+1,1.9979720e+1,2.1817860e+1,
     3 2.3728780e+1,2.5733030e+1,2.7860760e+1,3.0156810e+1,3.2697300e+1,
     4 3.5626950e+1,3.9331470e+1,4.7791030e+1/
      data (ce(i, 7),i=1,21)/    5.1518000e-1,1.7938600e+0,2.5771400e+0,
     1 3.2745600e+0,3.9273000e+0,4.5526800e+0,5.1606300e+0,5.7578800e+0,
     2 6.3496200e+0,6.9402800e+0,7.5340199e+0,8.1350700e+0,8.7480799e+0,
     3 9.3785600e+0,1.0033510e+1,1.0722500e+1,1.1459730e+1,1.2268640e+1,
     4 1.3194470e+1,1.4353420e+1,1.6972980e+1/
      data (ce(i, 8),i=1,21)/    5.3681000e-1,1.2748100e+0,1.6511800e+0,
     1 1.9657000e+0,2.2487300e+0,2.5125400e+0,2.7637900e+0,3.0066900e+0,
     2 3.2442600e+0,3.4788700e+0,3.7125900e+0,3.9473700e+0,4.1852200e+0,
     3 4.4284100e+0,4.6797000e+0,4.9428000e+0,5.2230800e+0,5.5293500e+0,
     4 5.8785000e+0,6.3138799e+0,7.2922699e+0/
      data (ce(i, 9),i=1,21)/    5.8986000e-1,9.4238999e-1,1.1027300e+0,
     1 1.2324200e+0,1.3469100e+0,1.4522500e+0,1.5516100e+0,1.6469600e+0,
     2 1.7396700e+0,1.8307700e+0,1.9211600e+0,2.0116400e+0,2.1030300e+0,
     3 2.1962100e+0,2.2922600e+0,2.3926000e+0,2.4992800e+0,2.6156400e+0,
     4 2.7480400e+0,2.9128300e+0,3.2821700e+0/
      data (ce(i,10),i=1,21)/    6.4590000e-1,8.5676000e-1,9.5119999e-1,
     1 1.0272100e+0,1.0941100e+0,1.1555200e+0,1.2133500e+0,1.2687700e+0,
     2 1.3225900e+0,1.3754400e+0,1.4278200e+0,1.4802100e+0,1.5330900e+0,
     3 1.5869800e+0,1.6425100e+0,1.7004800e+0,1.7620900e+0,1.8292500e+0,
     4 1.9056400e+0,2.0006700e+0,2.2135200e+0/
      data (ce(i,11),i=1,21)/    7.0787000e-1,8.3725999e-1,8.9556000e-1,
     1 9.4266000e-1,9.8425999e-1,1.0225600e+0,1.0587400e+0,1.0935100e+0,
     2 1.1273700e+0,1.1607100e+0,1.1938500e+0,1.2270900e+0,1.2607400e+0,
     3 1.2951300e+0,1.3306700e+0,1.3679000e+0,1.4075900e+0,1.4510100e+0,
     4 1.5006000e+0,1.5626000e+0,1.7026900e+0/
      data (ce(i,12),i=1,21)/    7.8624000e-1,8.6494000e-1,9.0017000e-1,
     1 9.2852999e-1,9.5349000e-1,9.7640999e-1,9.9800000e-1,1.0187000e+0,
     2 1.0388100e+0,1.0585700e+0,1.0781600e+0,1.0977700e+0,1.1175700e+0,
     3 1.1377600e+0,1.1585700e+0,1.1803100e+0,1.2034400e+0,1.2286600e+0,
     4 1.2573800e+0,1.2931400e+0,1.3733800e+0/
      data (ce(i,13),i=1,21)/    8.3683000e-1,8.9150999e-1,9.1594999e-1,
     1 9.3561000e-1,9.5291000e-1,9.6878000e-1,9.8372000e-1,9.9804000e-1,
     2 1.0119400e+0,1.0255900e+0,1.0391300e+0,1.0526600e+0,1.0663200e+0,
     3 1.0802400e+0,1.0945900e+0,1.1095700e+0,1.1254900e+0,1.1428500e+0,
     4 1.1625900e+0,1.1871600e+0,1.2422200e+0/
      data ep/200.,80.,40.,20.,10.,6.,4.,3.,2./
      data (dp(i,1),i=1,21)/       .50000,.47443,.44894,.42358,.39836,
     1 .37338,.34863,.32410,.29982,.27581,.25210,.22872,.20568,.18287,
     2 .16023,.13762,.11495,.09196,.06816,.04245,.00000/
      data (dp(i,2),i=1,21)/       .50000,.47577,.45159,.42753,.40359,
     1 .37982,.35621,.33274,.30944,.28629,.26323,.24027,.21735,.19441,
     2 .17135,.14807,.12443,.10011,.07446,.04614,.00000/
      data (dp(i,3),i=1,21)/       .50000,.47705,.45413,.43123,.40838,
     1 .38560,.36288,.34020,.31753,.29485,.27218,.24948,.22670,.20379,
     2 .18069,.15727,.13336,.10857,.08217,.05237,.00000/
      data (dp(i,4),i=1,21)/       .50000,.47861,.45720,.43578,.41431,
     1 .39280,.37125,.34962,.32793,.30614,.28417,.26204,.23967,.21699,
     2 .19381,.16996,.14517,.11897,.09064,.05895,.00000/
      data (dp(i,5),i=1,21)/       .50000,.47973,.45943,.43907,.41864,
     1 .39813,.37751,.35676,.33583,.31465,.29318,.27137,.24920,.22657,
     2 .20330,.17918,.15390,.12696,.09752,.06325,.00000/
      data (dp(i,6),i=1,21)/       .50000,.47973,.45946,.43919,.41892,
     1 .39864,.37830,.35772,.33692,.31587,.29454,.27275,.25049,.22764,
     2 .20390,.17930,.15363,.12622,.09608,.06173,.00000/
      data (dp(i,7),i=1,21)/       .50000,.47962,.45925,.43887,.41850,
     1 .39812,.37777,.35718,.33635,.31545,.29451,.27329,.25155,.22932,
     2 .20625,.18221,.15717,.13011,.10020,.06496,.00000/
      data (dp(i,8),i=1,21)/       .50000,.47915,.45830,.43744,.41659,
     1 .39574,.37489,.35403,.33318,.31248,.29120,.26904,.24600,.22301,
     2 .19878,.17326,.14770,.12050,.09134,.05931,.00000/
      data (dp(i,9),i=1,21)/       .50000,.47500,.45000,.42500,.40000,
     1 .37500,.35000,.32500,.30000,.27500,.25000,.22500,.20000,.17500,
     2 .15000,.12500,.10000,.07500,.05000,.02500,.00000/

*-----------------------------------------------------------------------

      dimension usave(3)
      dimension xv(3)

*-----------------------------------------------------------------------
*     save photon parameters while electrons are produced.
*-----------------------------------------------------------------------

      ipt9 = ipt
      erg9 = erg
      wgt9 = wgt
      npa9 = npa

      usave(1) = uuu
      usave(2) = vvv
      usave(3) = www

*-----------------------------------------------------------------------

         ipt  = 3
         npa = 1
         mm  = m

*-----------------------------------------------------------------------
*     roulette or split according to enum bias.
*-----------------------------------------------------------------------

         wga = wgt

         if(enum.eq.1.)go to 40
         if(ipt9.ne.2)go to 40

         npa = enum + rang()
         if( npa .eq. 0 ) goto 200
         wgt = wgt / enum

         wga = npa * wgt

*-----------------------------------------------------------------------
*     branch according to photon collision type.
*-----------------------------------------------------------------------

   40 goto (50,90,100,160) m

*-----------------------------------------------------------------------
*     m=1 -- make a photoelectron.
*-----------------------------------------------------------------------

   50 continue

         erg = cs

      do 60 k = 2, 13
   60    if(cs.ge.ee(k)) goto 80

   70    am = 1.-2.*rang()
         if(rang().lt.am**2) goto 70

         call dtcos(am,uold,uuu,0,irdm)

         goto 170

   80    ek = ee(k-1)-cs
         e1 = cs-ee(k)
         r = rang()*20.+1.
         l = r
         pl = r-l
         am = min(one,max(-one,(1.-(ee(k-1)-ee(k))/((1.-pl)*e1*ce(l,k-1)
     &      + pl*e1*ce(l+1,k-1)+(1.-pl)*ek*ce(l,k)+pl*ek*ce(l+1,k)))
     &      * (eg0+gpt(3))/sqrt(eg0*(eg0+2.*gpt(3)))))

         call dtcos(am,uold,uuu,0,irdm)

         goto 170

*-----------------------------------------------------------------------
*     m=2 -- make a compton recoil electron.
*-----------------------------------------------------------------------

   90    bt = eg0/erg
         erg = eg0-erg
         ec = min((bt-cs)/sqrt(bt*(bt-2.*cs)+1.),one)
         b = -sqrt((1.-ec**2)/(1.-cs**2))
         a = ec-b*cs
         uuu = a*uold(1)+b*uuu
         vvv = a*uold(2)+b*vvv
         www = a*uold(3)+b*www
         t1 = 1./sqrt(uuu**2+vvv**2+www**2)
         uuu = t1*uuu
         vvv = t1*vvv
         www = t1*www

         goto 170

*-----------------------------------------------------------------------
*     m=3 -- make an electron-positron pair.
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*     first sample the electron.
*-----------------------------------------------------------------------

  100    s = min(erg/gpt(3),ep(1))

      do 110 k = 2, 8
  110    if(s.ge.ep(k)) goto 120

  120    ek = ep(k-1)-s
         e1 = s-ep(k)
         r = rang()*20.+1.
         l = r
         pl = r-l

         em = erg-2.*gpt(3)
         erg = em*((1.-pl)*e1*dp(l,k-1)+pl*e1*dp(l+1,k-1)
     &       + (1.-pl)*ek*dp(l,k)+pl*ek*dp(l+1,k))/(ep(k-1)-ep(k))
         if(rang().gt..5) erg = em-erg
         b = sqrt(eg0*(eg0+2.*gpt(3)))/(eg0+gpt(3))
         r = 2.*rang()
         e1 = min(max(-one,(r-1.-b)/(r*b-1.-b)),one)

         call dtcos(e1,uold,uuu,0,irdm)

*-----------------------------------------------------------------------

         if( kpt(3) .ne. 0 ) then

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 7
                     nixcos = nixcos + 1  ! MATSUDA 2022.12.12  for T-Point

                     jclusts(0,nclsts) = 0
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(3,nclsts) = 12
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = -1
                     jclusts(6,nclsts) = 0
                     jclusts(7,nclsts) = 11
                     jclusts(8,nclsts) = 0

                     rms = gpt(3)
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
                     qclusts(8,nclsts)  = wga / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd
                     ixcoss(nixcos) = 0  ! MATSUDA 2022.12.12  for T-Point
                     aevts(iaevt+38,ireg) =
     &               aevts(iaevt+38,ireg) + wga
                     bevts(ibevt+38,imat) =
     &               bevts(ibevt+38,imat) + wga

                     aevts(iaevt+39,ireg) =
     &               aevts(iaevt+39,ireg) + wga * erg
                     bevts(ibevt+39,imat) =
     &               bevts(ibevt+39,imat) + wga * erg
                     numpal(12) = numpal(12) + 1
                     rumpal(12) = rumpal(12) + wga

         else

            if( erg .gt. elc(2) ) then

               do j = 1, npa

                  call brmgam(zero,ireg,imat,xv,dxd)

               end do

            end if

         end if

*-----------------------------------------------------------------------
*     then sample the positron.
*-----------------------------------------------------------------------

         erg = em - erg

         if( kpt(3) .eq. 0 ) cs = 1.

         r = 2.*rang()
         e2 = min(max(-one,(r-1.-b)/(r*b-1.-b)),one)
         b = -sqrt((1.-e2**2)/(1.-e1**2))
         a = e2-b*e1
         uuu = a*uold(1)+b*uuu
         vvv = a*uold(2)+b*vvv
         www = a*uold(3)+b*www
         t1 = 1./sqrt(uuu**2+vvv**2+www**2)
         uuu = t1*uuu
         vvv = t1*vvv
         www = t1*www

         mm = 5

         goto 170

*-----------------------------------------------------------------------
*     m=4 -- make an auger electron.
*-----------------------------------------------------------------------

  160    erg = cs

         call isos(uuu,0)

         goto 170

*-----------------------------------------------------------------------
*     bank the new electron
*-----------------------------------------------------------------------

  170  continue

         if( kpt(3) .ne. 0 ) then

                  if( mm .eq. 3 .or. mm .eq. 5 ) then
                     iv = 38
                     iw = 39
                  else if( mm .eq. 2 ) then
                     iv = 40
                     iw = 41
                  else if( mm .eq. 1 ) then
                     iv = 42
                     iw = 43
                  else if( mm .eq. 4 ) then
                     iv = 44
                     iw = 45
                  end if

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 7
                     nixcos = nixcos + 1  ! MATSUDA 2022.12.12  for T-Point

                     jclusts(0,nclsts) = 0
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(4,nclsts) = 0
                     jclusts(6,nclsts) = 0
                     jclusts(8,nclsts) = 0

                  if( mm .ne. 5 ) then

                     jclusts(3,nclsts) =  12
                     jclusts(7,nclsts) =  11
                     jclusts(5,nclsts) =  -1

                     numpal(12) = numpal(12) + 1
                     rumpal(12) = rumpal(12) + wga

                  else

                     jclusts(3,nclsts) =  13
                     jclusts(7,nclsts) = -11
                     jclusts(5,nclsts) =   1

                     numpal(13) = numpal(13) + 1
                     rumpal(13) = rumpal(13) + wga

                  end if

                     ers = max( 0.0d0, erg )
                     rms = gpt(3)
                     rmg = rms / 1000.
                     pr  = sqrt( ers * ( ers + 2.0 * rms ) ) / 1000.
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
                     qclusts(7,nclsts)  = ers
                     qclusts(8,nclsts)  = wga / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd
                     ixcoss(nixcos) = 0  ! MATSUDA 2022.12.12  for T-Point
                     aevts(iaevt+iv,ireg) =
     &               aevts(iaevt+iv,ireg) + wga
                     bevts(ibevt+iv,imat) =
     &               bevts(ibevt+iv,imat) + wga

                     aevts(iaevt+iw,ireg) =
     &               aevts(iaevt+iw,ireg) + wga * ers
                     bevts(ibevt+iw,imat) =
     &               bevts(ibevt+iw,imat) + wga * ers

*-----------------------------------------------------------------------
*     make ttb photons.
*-----------------------------------------------------------------------

         else

            if( erg .gt. elc(2) ) then

               do j = 1, npa

                  call brmgam(zero,ireg,imat,xv,dxd)

               end do

            end if

         end if

*-----------------------------------------------------------------------
         if( m .eq. 1 ) then
            atmrc(1,3) = atmrc(1,3) + wga / wg0
         elseif( m .eq. 3 ) then
            atmrc(1,5) = atmrc(1,5) + wga / wg0
         elseif( m .eq. 4 ) then
            if( ielpos .eq.  0 ) atmrc(1,2) = atmrc(1,2) + wga / wg0
            if( ielpos .eq. -1 ) atmrc(2,2) = atmrc(2,2) + wga / wg0
            if( ielpos .eq.  1 ) atmrc(3,2) = atmrc(3,2) + wga / wg0
         endif

*-----------------------------------------------------------------------

  200 continue

         ipt = ipt9
         erg = erg9
         wgt = wgt9
         npa = npa9

         uuu = usave(1)
         vvv = usave(2)
         www = usave(3)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine brmgam(d,ireg,imat,xv,dxd)
*                                                                      *
*        generate bremsstrahlung photons with the thick-target         *
*        bremsstrahlung approximation.                                 *
*        d = length of path segment.                                   *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for REDUCTION_COUNTER
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
C MATSUDA 2022.12.12  for T-Point
      common /clusttp/ nixcos, ixcoss(nnn)
!$OMP THREADPRIVATE(/clusttp/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

C for USE_MOD_COUNTER

*-----------------------------------------------------------------------
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)

*-----------------------------------------------------------------------

      dimension xv(3)
      dimension usave(3)

*-----------------------------------------------------------------------

         erg9 = erg
         wgt9 = wgt
         ipt9 = ipt
         npa9 = npa

         usave(1) = uuu
         usave(2) = vvv
         usave(3) = www

*-----------------------------------------------------------------------
*     find the electron energy group.
*-----------------------------------------------------------------------

         m = 1
         ib = nee
   10    if(ib-m.eq.1) goto 30
         ih = (m+ib)/2
         if(erg.ge.eee(ih)) goto 20
         m = ih
         goto 10
   20    ib = ih
         goto 10
   30    f = (eee(m)-erg)/(eee(m)-eee(m+1))

*-----------------------------------------------------------------------
*     decide how many photons to make.
*-----------------------------------------------------------------------

         mkc = imat
         pr = (1.-f)*pbt(m+nee*(mkc-1))+f*pbt(m+1+nee*(mkc-1))

         pr = pr * abs(bnum)
         ni = pr + rang()

*-----------------------------------------------------------------------
*     no photon
*-----------------------------------------------------------------------

      if( ni .eq. 0 ) return

*-----------------------------------------------------------------------
*     prepare to make photons; protect variables from tallyd, dxtran.
*-----------------------------------------------------------------------

         u = uold(1)
         v = uold(2)
         w = uold(3)

         uold(1) = uuu
         uold(2) = vvv
         uold(3) = www

         if( abs(bnum) .gt. 0. .and. bnum .ne. 1. )
     &                                 wgt = wgt / abs(bnum)

         wt = wgt

         ipt = 2
         vel = slite
         npa = 1
         ncp = 0
         jsu = 0
         idx = 0
         ipsc = 13

         deb = sqrt(erg*(erg+2.*gpt(3)))/(erg+gpt(3))

*-----------------------------------------------------------------------
*     make photons and bank them.
*     approximation:  photon direction equals electron direction.
*-----------------------------------------------------------------------

      do 160 n = 1, ni

         rn = rang()
         pk = qpol(rn,ebt(1,m+nee*(mkc-1)),rkt,ntop)
         pa = qpol(rn,ebt(1,m+1+nee*(mkc-1)),rkt,ntop)

         erg = ((1.-f)*pk+f*pa)*erg9

         if( erg .lt .elc(2) ) goto 160

         wgt = wt

         if(mbi(mkc).eq.0) goto 150

         do 110 j1 = 2, ntop - 1
  110    if(pk.le.rkt(j1)) goto 120

  120    do 130 j = 2, ntop - 1
  130    if(pa.le.rkt(j)) goto 140

  140    wgt = wt*((1.-f)*ftt(m+nee*(mkc-1))*bbrem(j1)+
     &    f*ftt(m+1+nee*(mkc-1))*bbrem(j))

  150    continue

*-----------------------------------------------------------------------
*           bank all photons
*-----------------------------------------------------------------------

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 4
                     nixcos = nixcos + 1  ! MATSUDA 2022.12.12  for T-Point

                     jclusts(0,nclsts) = ipsc
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
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd
                     ixcoss(nixcos) = 0  ! MATSUDA 2022.12.12  for T-Point
                     aevts(iaevt+26,ireg) =
     &               aevts(iaevt+26,ireg) + wgt
                     bevts(ibevt+26,imat) =
     &               bevts(ibevt+26,imat) + wgt

                     aevts(iaevt+27,ireg) =
     &               aevts(iaevt+27,ireg) + wgt * erg
                     bevts(ibevt+27,imat) =
     &               bevts(ibevt+27,imat) + wgt * erg
                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt
         atmrc(2,5) = atmrc(2,5) + wgt / wg0

*-----------------------------------------------------------------------

  160 continue

*-----------------------------------------------------------------------

         uold(1) = u
         uold(2) = v
         uold(3) = w

         uuu = usave(1)
         vvv = usave(2)
         www = usave(3)

         erg = erg9
         wgt = wgt9
         ipt = ipt9
         npa = npa9

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine flugam(er,ie,mk,ik,ireg,imat,xv,dxd)
*                                                                      *
*        generate fluourescent photons or auger electrons to fill      *
*        vacancies.  called by sctgam and kxrelc.                      *
*                                                                      *
*        input arguments (not changed by subroutine):                  *
*           er=energy of absorbed photon or edge energy for kxrelc.    *
*           ie=element that is undergoing the interaction; for kxrelc  *
*              it is the highest z element.                            *
*          mk=material.                                                *
*          ik=ipt particle type of parent particle used for bookkeeping*
*              and variance reduction; =2 for sctgam; =3 for kxrelc.   *
*                                                                      *
*        last modified by K.Niita on 2012/12/27                        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for REDUCTION_COUNTER
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

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

C for USE_MOD_COUNTER

*-----------------------------------------------------------------------
      common /elpos/  ielpos
!$OMP THREADPRIVATE(/elpos/)
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)

*-----------------------------------------------------------------------

      dimension usave(3)
      dimension xv(3)

*-----------------------------------------------------------------------

         ntyn = 3
         mtp = -3
         rp = 0.

         usave(1) = uuu
         usave(2) = vvv
         usave(3) = www

         erg9 = erg
         wgt9 = wgt
         ipt9 = ipt
         npa9 = npa

*-----------------------------------------------------------------------
*     if incident particle is electron from kxrelc, setup photon
*     parameters and sample auger or photoelectric electrons.
*-----------------------------------------------------------------------

      if( ik .eq. 2 ) goto 10

         ipt = 2
         vel = slite
         ncp = 0
         idx = 0
         if(rang().gt.wwk(mk)) goto 130

*-----------------------------------------------------------------------
*     integrated tiger series (its) method before its3.0.
*-----------------------------------------------------------------------

      if(nxs(16,lme(3,jmd(1+mk))).eq.3) goto 10
      if(iphot.ne.0.or.eek(mk).le.elc(2)) goto 900

               erg = eek(mk)
               call isos(uuu,0)
               ipsc = 102

*-----------------------------------------------------------------------
*        bank first photon from electron, electron x-ray
*-----------------------------------------------------------------------

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 4
                     nixcos = nixcos + 1  ! MATSUDA 2022.12.12  for T-Point

                     jclusts(0,nclsts) = ipsc
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
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd
                     ixcoss(nixcos) = 0  ! MATSUDA 2022.12.12  for T-Point
                     aevts(iaevt+32,ireg) =
     &               aevts(iaevt+32,ireg) + wgt
                     bevts(ibevt+32,imat) =
     &               bevts(ibevt+32,imat) + wgt

                     aevts(iaevt+33,ireg) =
     &               aevts(iaevt+33,ireg) + wgt * erg
                     bevts(ibevt+33,imat) =
     &               bevts(ibevt+33,imat) + wgt * erg
                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt
            if( ielpos .eq. -1 ) atmrc(2,1) = atmrc(2,1) + wgt / wg0
            if( ielpos .eq.  1 ) atmrc(3,1) = atmrc(3,1) + wgt / wg0

            goto 900

*-----------------------------------------------------------------------
*     check for electron emission from impact ionzation.
*-----------------------------------------------------------------------

   10    if(nxs(2,ie).le.11) goto 140
         i1 = jxs(4,ie)+2
         if(er.le.xss(i1-1)) goto 140

      do 20 ka = i1, i1-3+nxs(4,ie)
   20    if(er.le.xss(ka)) goto 30

   30    k0 = ka+nxs(4,ie)-1
   40    ka = k0
         t3 = xss(ka)*rang()
         ka = ka+nxs(4,ie)
         if(t3.lt.xss(ka)) goto 50
         if(ik.eq.2) goto 140
         goto 40

*-----------------------------------------------------------------------
*     at least single fluorescence or auger electron.
*-----------------------------------------------------------------------

   50    ka = ka-1
         if(t3.lt.xss(ka)) goto 50

         rp = 1.
         eg = xss(ka-2*nxs(4,ie)+1)
         if(kpt(3).eq.0) goto 60

         if(dbcn(20).ne.0..and.ik.ne.2) r = rang()
         if(ik.ne.2) goto 60
         if(rang().gt.wwk(mk)) goto 130

   60    e = xss(nxs(4,ie)+ka+1)

         if(e.lt.elc(2)) goto 140

*-----------------------------------------------------------------------
*           photon first fluorescence
*-----------------------------------------------------------------------

            if( ik .eq. 2 ) then
                     aevts(iaevt+34,ireg) =
     &               aevts(iaevt+34,ireg) + wgt
                     bevts(ibevt+34,imat) =
     &               bevts(ibevt+34,imat) + wgt

                     aevts(iaevt+35,ireg) =
     &               aevts(iaevt+35,ireg) + wgt * e
                     bevts(ibevt+35,imat) =
     &               bevts(ibevt+35,imat) + wgt * e
                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

            end if

*-----------------------------------------------------------------------

         ed = 0.
         if(nxs(2,ie).le.30) goto 100
         i = ka+1-2*nxs(4,ie)
         if(i.ne.i1.and.i.ne.i1+1) goto 100
         i = i1-2+nxs(4,ie)
         ka = nxs(4,ie)+i+1
         if(rang().ge.xss(ka)/(1.-xss(i))) goto 100

*-----------------------------------------------------------------------
*     double fluorescence.  bank the second photon.
*-----------------------------------------------------------------------

         erg = xss(nxs(4,ie)+ka)

         if(erg.lt.elc(2)) goto 100

         ntyn = 4

         call isos(uuu,0)

         sw = wgt

*-----------------------------------------------------------------------

            if( ik .eq. 2 ) then
                     aevts(iaevt+34,ireg) =
     &               aevts(iaevt+34,ireg) + wgt
                     bevts(ibevt+34,imat) =
     &               bevts(ibevt+34,imat) + wgt

                     aevts(iaevt+35,ireg) =
     &               aevts(iaevt+35,ireg) + wgt * erg
                     bevts(ibevt+35,imat) =
     &               bevts(ibevt+35,imat) + wgt * erg
                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

                  if( nter .ne. 0 ) then

                     nter = 0
                     goto 100

                  end if

            else if( ik .eq. 3 ) then
                     aevts(iaevt+32,ireg) =
     &               aevts(iaevt+32,ireg) + wgt
                     bevts(ibevt+32,imat) =
     &               bevts(ibevt+32,imat) + wgt

                     aevts(iaevt+33,ireg) =
     &               aevts(iaevt+33,ireg) + wgt * erg
                     bevts(ibevt+33,imat) =
     &               bevts(ibevt+33,imat) + wgt * erg
                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

            end if

*-----------------------------------------------------------------------
*     bank second photon from both
*-----------------------------------------------------------------------
*              photon heating flag : jclusts(4,nclsts) = 1
*-----------------------------------------------------------------------

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 4
                     nixcos = nixcos + 1  ! MATSUDA 2022.12.12  for T-Point

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(3,nclsts) = 14

                  if( ik .eq. 2 ) then

                     jclusts(4,nclsts) = 1

                  else

                     jclusts(4,nclsts) = 0

                  end if

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
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd
                     ixcoss(nixcos) = 0  ! MATSUDA 2022.12.12  for T-Point

*-----------------------------------------------------------------------

         npa = 1
         ed  = erg
         wgt = sw

*-----------------------------------------------------------------------
*     continue with first fluorescent photon.
*-----------------------------------------------------------------------

  100    if(ides.eq.0.and.ik.eq.2)
     &   call elcgam(1,eg0-eg,ireg,imat,xv,dxd)


         erg  = e
         ipsc = 99+ntyn

         call isos(uuu,0)

*-----------------------------------------------------------------------
*     bank first photon from both
*-----------------------------------------------------------------------
*              photon heating flag : jclusts(4,nclsts) = 1
*-----------------------------------------------------------------------

      if( ik .eq. 2 .or. ik .eq. 3 ) then

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 4
                     nixcos = nixcos + 1  ! MATSUDA 2022.12.12  for T-Point

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(3,nclsts) = 14

                  if( ik .eq. 2 ) then

                     jclusts(4,nclsts) = 1

                  else

                     jclusts(4,nclsts) = 0

                  end if

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
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd
                     ixcoss(nixcos) = 0  ! MATSUDA 2022.12.12  for T-Point
         if( ik .eq. 2 ) atmrc(1,1) = atmrc(1,1) + wgt / wg0
         if( ik .eq. 3 .and. ielpos .eq. -1 )
     &                   atmrc(2,1) = atmrc(2,1) + wgt / wg0
         if( ik .eq. 3 .and. ielpos .eq.  1 )
     &                   atmrc(3,1) = atmrc(3,1) + wgt / wg0


         goto 900

      end if

*-----------------------------------------------------------------------
*     make auger electrons.
*     set ipt=ik for bookkeeping.
*-----------------------------------------------------------------------

  130    ipt = ik

         if(ides.eq.0.and.er.ge.edg(mk))
     &      call elcgam(4,eek(mk),ireg,imat,xv,dxd)

*-----------------------------------------------------------------------
*     make photoelectric electrons.
*-----------------------------------------------------------------------

  140 if(ik.ne.2) goto 900

         if( ides .eq. 0 )
     &    call elcgam(1,eg0-rp*edg(mk),ireg,imat,xv,dxd)

         nter = 12

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------

         uuu = usave(1)
         vvv = usave(2)
         www = usave(3)

         erg = erg9
         wgt = wgt9
         ipt = ipt9
         npa = npa9

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine jqmdchk(eein,mmas,mchg)
*                                                                      *
*        check conservation after call JQMD for photonuclear reaction  *
*        this routine is made with reference lines 357-601, ncasc.f    *
*                                                                      *
*        made by S.Abe on 2015/08/10                                   *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg

      implicit real*8 (a-h,o-z)

      include 'param00.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

      common /qmdflg/ irqmd

*-----------------------------------------------------------------------

      data rmnuc / 938.95d0 /

*-----------------------------------------------------------------------
*        check of the total energy conservation 2008/06/20
*-----------------------------------------------------------------------

      if( nclst .le. 0 .and. irqmd .eq. 0) return

*-----------------------------------------------------------------------

      eout = 0.d0
      ebin = 0.d0
      exct = 0.d0

      do i = 1, nclst

       exct = exct + qclust(6,i)
       eout = eout + qclust(7,i)

       if( iclust(i) .eq. 0 ) then
        ebin = ebin + bindeg(jclust(1,i),jclust(2,i))
       else if( iclust(i) .gt. 2 ) then
        eout = eout + qclust(5,i) * 1000.d0
       end if

      end do

      eres = eein - eout + ebin - bindeg(mchg,mmas-mchg)

*-----------------------------------------------------------------------

      if( eres .lt. 0.d0 ) then
       nclst = -1
       return
      else if( exct .eq. 0.d0 ) then
       erat = 0.d0
      else
       erat = eres / exct
      end if

      do i = 1, nclst
       qclust(6,i) = qclust(6,i) * erat
      end do

*-----------------------------------------------------------------------
      return
      end

