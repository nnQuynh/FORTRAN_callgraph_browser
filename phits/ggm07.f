************************************************************************
*                                                                      *
      subroutine sctelc(eein,wgti,ireg,imat,xv,rhi,eezo)
*                                                                      *
*        calculate a collision of a electron with an atom.             *
*        last modified by K.Niita on 2011/08/30                        *
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

      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, no, mtel
!$OMP THREADPRIVATE(/elect/)
      common /elpos/  ielpos
!$OMP THREADPRIVATE(/elpos/)

*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
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


*-----------------------------------------------------------------------
      common /dedxfac/ dedxfd
!$OMP THREADPRIVATE(/dedxfac/)

*-----------------------------------------------------------------------
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)

*-----------------------------------------------------------------------

      dimension xv(3)
      dimension usave(3), suu(3,2)

*-----------------------------------------------------------------------

           nclst  = -1
           nclsts = 0

*-----------------------------------------------------------------------

            nter = 0
            ipt  = 3

            erg  = eein
            eg0  = eint

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

            end do

            d = delc
            mtel = imat

*-----------------------------------------------------------------------
*        scattering energy loss
*-----------------------------------------------------------------------
                  aevts(iaevt+48,ireg) =
     &            aevts(iaevt+48,ireg) + wgt * d * qs
                  bevts(ibevt+48,imat) =
     &            bevts(ibevt+48,imat) + wgt * d * qs
*-----------------------------------------------------------------------
*     save the incoming direction respect to the direction u,v,w.
*     rev. 2010/01/04 by K.N
*-----------------------------------------------------------------------

            uold(1) = uuu
            uold(2) = vvv
            uold(3) = www

*-----------------------------------------------------------------------
*        lose energy to ionization along the track.
*-----------------------------------------------------------------------

         if( erg .lt. elc(3) ) goto 170

         if( erg .lt. eee(ngp+1) ) goto 100

         goto 120

  100    ngp = ngp + 1

         if( erg .lt. eee(ngp+1) ) goto 100

  120    continue

*-----------------------------------------------------------------------
*     produce secondary particles and lose energy to bremsstrahlung.
*-----------------------------------------------------------------------

         npa = 1

         f = (eee(ngp)-erg) / (eee(ngp)-eee(ngp+1))

*-----------------------------------------------------------------------
*        generate knock-on electrons along path segment d.
*-----------------------------------------------------------------------

         if( rnok .ne. 0. ) then

                  nclstsb = nclsts

            call knoelc(d,f,no,am,xv,ireg,imat,rhi)

            if( nclsts .gt. nclstsb ) then
                  eout = 0.d0
               do k = nclstsb + 1, nclsts
c T.Sato 2014/2/21 bug for photon transport using weight control
                  eout = eout + qclusts(7,k) * qclusts(8,k)
               end do
                  erg = max( 0.d0, erg - eout )
            end if

         end if

*-----------------------------------------------------------------------
*        generate x-rays and/or auger electrons along path segment d.
*-----------------------------------------------------------------------

         if( erg .ge. edg(mkc) .and.
     &       eek(mkc) .gt. min(elc(2),elc(3)) ) then

                  nclstsb = nclsts

            call kxrelc(d,f,xv,ireg,imat,rhi)

            if( nclsts .gt. nclstsb ) then
                  eout = 0.d0
               do k = nclstsb + 1, nclsts
c T.Sato 2014/2/21 bug for photon transport using weight control
                  eout = eout + qclusts(7,k) * qclusts(8,k)
               end do
                  erg = max( 0.d0, erg - eout )
            end if

         end if

*-----------------------------------------------------------------------
*        generate a bremsstrahlung photon along path segment d.
*-----------------------------------------------------------------------

                  nclstsb = nclsts

            call brmelc(d,qo,f,no,am,xv,ireg,imat,rhi)

            if( nclsts .gt. nclstsb ) then
                  eout = 0.d0
               do k = nclstsb + 1, nclsts
c T.Sato 2014/2/21 bug for photon transport using weight control
                  eout = eout + qclusts(7,k) * qclusts(8,k)
               end do
                  erg = max( 0.d0, erg - eout )
            end if

*-----------------------------------------------------------------------
*        check energy bin
*-----------------------------------------------------------------------

         if( erg .lt. elc(3) ) goto 170

         if( erg .ge. eee(ngp+1) ) goto 140

  130    ngp = ngp + 1

         if( erg .lt. eee(ngp+1) ) goto 130

  140    continue

*-----------------------------------------------------------------------
*     process the end of the energy substep.
*-----------------------------------------------------------------------

         ns = ns - 1

*-----------------------------------------------------------------------
*     dbcn(18).ne.0:  its3.0 nearest group boundary data treatment.
*-----------------------------------------------------------------------

         if( dbcn(18) .ne. 0. ) then

            if( ns .eq. 0 ) then

               nq = ngp

               if( erg .lt. 0.5*(eee(ngp)+eee(ngp+1) ).and.
     &             ngp .lt. nee-1 ) nq = ngp + 1

               ns = nsb(mkc)
               qs = esloss(nq,rhi,qsexx)

            end if

*-----------------------------------------------------------------------
*     dbcn(18)=0:  bin-centered data treatment.
*-----------------------------------------------------------------------

         else

            if( n1 .ne. ngp ) then

               nq = ngp
               ns = nsb(mkc)
               qs = esloss(nq,rhi,qsexx)

            end if

         end if

*-----------------------------------------------------------------------
*     sequential electron or positron trnasfer
*-----------------------------------------------------------------------

         goto 500

*-----------------------------------------------------------------------
*     energy cut-off
*-----------------------------------------------------------------------

  170 continue

            ns = 0

*-----------------------------------------------------------------------
*     make ttb photons.
*-----------------------------------------------------------------------

            if( iphot .eq. 0 .and. erg .gt. elc(2) ) then

                        nclstsb = nclsts

               call brmgam(d,ireg,imat,xv,0.5*d)

                  if( nclsts .gt. nclstsb ) then
                        eout = 0.d0
                     do k = nclstsb + 1, nclsts
c T.Sato 2014/2/21 bug for photon transport using weight control
                        eout = eout + qclusts(7,k) * qclusts(8,k)
                     end do
                        erg = max( 0.d0, erg - eout )
                  end if

            end if

*-----------------------------------------------------------------------
*     generate annihilation photons for positron below energy cutoff.
*-----------------------------------------------------------------------

            if( ielpos .eq. 1 .and. gpt(3) .ge. elc(2) ) then

                  ergp = gpt(3)
                  call isos(uuu,0)

                     suu(1,1) = uuu
                     suu(2,1) = vvv
                     suu(3,1) = www

                     suu(1,2) = -uuu
                     suu(2,2) = -vvv
                     suu(3,2) = -www

*-----------------------------------------------------------------------

                  do k = 1, 2

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 4

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(3,nclsts) = 14
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = 0
                     jclusts(6,nclsts) = 0
                     jclusts(7,nclsts) = 22
                     jclusts(8,nclsts) = 0

                     pr   = ergp / 1000.0

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
                     qclusts(7,nclsts)  = ergp
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = xv(1) * d * 0.5
                     qclusts(11,nclsts) = xv(2) * d * 0.5
                     qclusts(12,nclsts) = xv(3) * d * 0.5
                     aevts(iaevt+30,ireg) =
     &               aevts(iaevt+30,ireg) + wgt
                     bevts(ibevt+30,imat) =
     &               bevts(ibevt+30,imat) + wgt

                     aevts(iaevt+31,ireg) =
     &               aevts(iaevt+31,ireg) + wgt * ergp
                     bevts(ibevt+31,imat) =
     &               bevts(ibevt+31,imat) + wgt * ergp
                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

                  end do
                  atmrc(3,5) = atmrc(3,5) + wgt / wg0    ! S.Abe 2017/02/08

! T.Sato 2015/2/26 include original positron as a dead particle
                     nclsts = nclsts + 1
                     iclusts(nclsts) = 7      ! 7:e-,e+, 4:photon

                     jclusts(0,nclsts) = 0
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(3,nclsts) = 13   ! 12:e-, 13:e+, 14:photon
                     jclusts(4,nclsts) = -1   ! dead particle
                     jclusts(5,nclsts) = 1    ! charge
                     jclusts(6,nclsts) = 0
                     jclusts(7,nclsts) = -11  ! kf code
                     jclusts(8,nclsts) = 0

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = 0.0
                     qclusts(2,nclsts)  = 0.0
                     qclusts(3,nclsts)  = 0.0
                     qclusts(4,nclsts)  = 0.0
                     qclusts(5,nclsts)  = 0.511008000000000*1.0d-3
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = erg
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = xv(1) * d * 0.5
                     qclusts(11,nclsts) = xv(2) * d * 0.5
                     qclusts(12,nclsts) = xv(3) * d * 0.5

                  return

            end if


*-----------------------------------------------------------------------
*     bank the final electron or positron
*-----------------------------------------------------------------------

  500 continue

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 7

                     jclusts(0,nclsts) = 0
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(4,nclsts) = 0
                     jclusts(6,nclsts) = 0

                  if( ielpos .eq. -1 ) then

                     jclusts(3,nclsts) = 12
                     jclusts(5,nclsts) = -1
                     jclusts(7,nclsts) = 11

                  else if( ielpos .eq. 1 ) then

                     jclusts(3,nclsts) = 13
                     jclusts(5,nclsts) =  1
                     jclusts(7,nclsts) = -11

                  end if

                     ers = max( 0.0d0, erg )

                     rms = gpt(3)
                     rmg = rms / 1000.
                     pr  = sqrt( ers * ( ers + 2.0 * rms ) ) / 1000.
                     ett = sqrt( pr**2 + rmg**2 )

                     uuu = 0.0
                     vvv = 0.0
                     www = 1.0

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
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = 0.0d0
                     qclusts(11,nclsts) = 0.0d0
                     qclusts(12,nclsts) = 0.0d0


*-----------------------------------------------------------------------
*     finish
*-----------------------------------------------------------------------

            kdecay(4) = 0

            uuu = usave(1)
            vvv = usave(2)
            www = usave(3)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine eletot(pmff,icli,ein,mki,rhi)
*                                                                      *
*        calculate a collision of a electron with an atom.             *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, no, mtel
!$OMP THREADPRIVATE(/elect/)

      data ns / 0 /

*-----------------------------------------------------------------------
*     initialization of this material and particle
*-----------------------------------------------------------------------

               erg  = ein
               icl  = icli
               mkc  = mki
               pmf  = huge
               pmff = huge
               mtel = mki

               if( mkc .le. 0 ) return

*-----------------------------------------------------------------------
*     ngp : find the electron energy group.
*-----------------------------------------------------------------------

      if( ns .eq. 0 ) then

               ngp = 1
               ib = nee
   10          if(ib-ngp.eq.1) goto 25
               ih = (ngp+ib)/2
               if(erg.ge.eee(ih)) goto 20
               ngp = ih
               goto 10
   20          ib = ih
               goto 10
   25          continue

*-----------------------------------------------------------------------
*     set up the material index
*-----------------------------------------------------------------------

               nq = ngp

               if( dbcn(18) .ne. 0. .and.
     &             erg .lt. 0.5*(eee(ngp)+eee(ngp+1)) .and.
     &             ngp .lt. nee-1 ) nq = ngp + 1

               qs = esloss(nq,rhi,qsexx)
               ns = nsb(mkc)

      else

               qs = esloss(nq,rhi,qsexx)

      end if



*-----------------------------------------------------------------------
*     calculate the distances to energy substep.
*-----------------------------------------------------------------------

               n1   = ngp
               pmf  = drs(nq+nee*(mkc-1))/rhi
               pmff = pmf

               qsex = qsexx

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine defelc(aq,dk,fd,ae,ik,
     &                 uu1,uu2,uu3,ua1,ua2,ua3)
*                                                                      *
*        electron scatter and deflection correction.                   *
*        return new cosine, aq, and direction, uuu, for particle.      *
*        return new uold for bremsstrahlung photons for tallyd, dxtran.*
*                                                                      *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, no, mtel
!$OMP THREADPRIVATE(/elect/)

*-----------------------------------------------------------------------

      dimension ue(3), uu(3), uo(3)

*-----------------------------------------------------------------------

         ue(1) = ua1
         ue(2) = ua2
         ue(3) = ua3

         uo(1) = ua1
         uo(2) = ua2
         uo(3) = ua3

         if(ik.ne.1.and.dbcn(17).ne.0.) goto 40

*-----------------------------------------------------------------------
*     dbcn(17)=0: proper delta scattering. xtm:95-236,319
*-----------------------------------------------------------------------

         if(dk.ge.pmf) goto 10

         ep = 1.d-10
         if(rang().gt.max(ep,1.-egg(maxi,nq+nee*(mkc-1)))**
     1    (dk/pmf)) goto 10

         if(ik.eq.1) return

         goto 60

   10    rn = rang()

      do 20 m = 2, maxi
   20    if(egg(m,nq+nee*(mkc-1)).gt.rn) goto 30

         if(dk.lt.pmf) goto 10
         return

   30    a = calph(m)+(calph(m-1)-calph(m))
     &     * (egg(m,nq+nee*(mkc-1))-rn)
     &     / (egg(m,nq+nee*(mkc-1))-egg(m-1,nq+nee*(mkc-1)))
         if(dk.lt.pmf) a = 1.-(1.-a)*dk/pmf
         if(ik.eq.1) aq = a
         if(ik.gt.1) call dtcos(a,uo,ue,0,irdm)
         goto 60

*-----------------------------------------------------------------------
*     dbcn(17)>0: schlumberger (sdr) scattering method. xtm:95-267
*-----------------------------------------------------------------------

   40    if(dbcn(17).lt.0.) goto 50

         tpp(1) = fd*(1.-ae)
         tpp(2) = sqrt((2.*fd-fd*tpp(1))/(1.+ae))
         tpp(3) = 1.-tpp(1)-ae*tpp(2)

         ue(1) = tpp(2)*uu1+tpp(3)*ua1
         ue(2) = tpp(2)*uu2+tpp(3)*ua2
         ue(3) = tpp(2)*uu3+tpp(3)*ua3

         goto 60

*-----------------------------------------------------------------------
*     dbcn(17)<0: its3.0 methodology.
*-----------------------------------------------------------------------

   50 if(fd.le..5) goto 60

         ue(1) = uu1
         ue(2) = uu2
         ue(3) = uu3

*-----------------------------------------------------------------------
*     get new direction, uuu, for cosine aq.
*-----------------------------------------------------------------------

   60    aq = min(max(aq,-one),one)
         call dtcos(aq,ue,uu,0,irdm)

         uu1 = uu(1)
         uu2 = uu(2)
         uu3 = uu(3)

*-----------------------------------------------------------------------
*     change uold for bremsstrahlung photons for tallyd, dxtran.
*-----------------------------------------------------------------------

         if( ik .ne. 2 ) return

         ua1 = ue(1)
         ua2 = ue(2)
         ua3 = ue(3)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function esloss(nq,rhi,qsex)
*                                                                      *
*        sample the electron energy loss rate, with optional straggling*
*                                                                      *
*        istrg : 0; sampled straggling for electron energy loss        *
*                1; expaected value straggling for electron energy loss*
*                                                                      *
*          nq, mkc, erg                                                *
*                                                                      *
*        modified by K.Niita on 2009/09/30                             *
*        last modified by T.Sato on 2011/09/05                         *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

      parameter (mlam=5001)
      common /lancom/ eqlm(mlam)


*-----------------------------------------------------------------------

      parameter (bs=euler-1.5d+0)
      parameter (bc=.989571788d+0,h2=2.08643410365460141d-2)

*-----------------------------------------------------------------------
*     Modified by T.Sato on 2011/9/5 *****
*     define expected dEdx for electron
*-----------------------------------------------------------------------

         qsex = qav0(nq+nee*(mkc-1))*rhi

*-----------------------------------------------------------------------
*     bc = probability that lambda.lt.100., sampling landau function.
*     h2 =  2/f, with f = integral of borsch-supan asymptotic form.
*-----------------------------------------------------------------------

         if(istrg.eq.0) goto 10
         esloss = qav(nq+nee*(mkc-1))*rhi
         return

   10    if(rang().le.bc) goto 20
         w = 1./(sqrt(1.+h2*rang())-1.)
         fl = w+log(w)+bs
         goto 30

   20    f = (mlam-1)*rang()
         j = f
         fl = eqlm(j+1)+(f-j)*(eqlm(j+2)-eqlm(j+1))

   30    if(fl.gt.flc(nq+nee*(mkc-1))) goto 10

   40    t1 = 2.*rang()-1.
         t2 = t1**2+rang()**2
         if(t2.gt.1.) goto 40

         q = qav(nq+nee*(mkc-1))
     &     + asp(nq+nee*(mkc-1))*(fl-ear(nq+nee*(mkc-1)))
     &     + qcn(nq+nee*(mkc-1))*sqrt(-2.*log(t2)/t2)*t1

         if(q*drs(nq+nee*(mkc-1))*nsb(mkc).ge.erg) goto 10
         if(q.lt.0.) q = 0.
         esloss = q*rhi

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine knoelc(d,f,no,ae,xv,ireg,imat,rhi)
*                                                                      *
*        generate knock-on electrons along path segment d.             *
*        f = energy interpolation factor.                              *
*        ae = deflection cosine of electron over path segment d.       *
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

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)

C for USE_MOD_COUNTER

*-----------------------------------------------------------------------
      common /elpos/  ielpos
!$OMP THREADPRIVATE(/elpos/)
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)

*-----------------------------------------------------------------------

      dimension xv(3)
      dimension usave(3)

*-----------------------------------------------------------------------

         erg9 = erg
         wgt9 = wgt

         usave(1) = uuu
         usave(2) = vvv
         usave(3) = www

*-----------------------------------------------------------------------

            ni = rang() + d * rnok * rhi
     &                  * ( (1.-f) * pkn(ngp+nee*(mkc-1) )
     &                         + f * pkn(ngp+1+nee*(mkc-1) ) )

            if( ni .eq. 0 ) return

            wgt = wgt / rnok

      do 140 n = 1, ni
            rn = rang()
            p = (1.-f)*pru(ngp)+f*pru(ngp+1)

            es = ecf(3)/eee(ngp)
            s = eee(ngp)/gpt(3)

            if(rn.gt.p) goto 70
            ag = (s/(s+1.))**2
            bs = 1.-ag
            e1 = 1./es-2.
            if(rn.gt.(e1+bs*log(2.*es))/(p+ag*(.5-es))) goto 50
            b = 1./(1.-es*bs)
   40       as = 1./(2.+rang()*e1)
            if(rang().gt.(1.-as*bs)*b) goto 40
            goto 80

   50       e2 = 2.-1./(1.-es)
            b = 1./(1.-.5*bs)
   60       as = 1./(2.-rang()*e2)
            if(rang().gt.(1.-as*bs)*b) goto 60
            as = 1.-as
            goto 80
   70       as = .5-rang()*(.5-es)

   80       erg = as * eee(ngp)

            if( erg .lt. elc(3) ) goto 140

*-----------------------------------------------------------------------

            fd = rang()
            dk = fd * d
            dxd = d - dk
            aq = sqrt(as*(s+2.)/(as*s+2.))

            call defelc(aq,dk,fd,ae,3,
     &                 uuu,vvv,www,uold(1),uold(2),uold(3))

*-----------------------------------------------------------------------

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 7

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
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd
                     aevts(iaevt+46,ireg) =
     &               aevts(iaevt+46,ireg) + wgt
                     bevts(ibevt+46,imat) =
     &               bevts(ibevt+46,imat) + wgt

                     aevts(iaevt+47,ireg) =
     &               aevts(iaevt+47,ireg) + wgt * erg
                     bevts(ibevt+47,imat) =
     &               bevts(ibevt+47,imat) + wgt * erg
                     numpal(12) = numpal(12) + 1
                     rumpal(12) = rumpal(12) + wgt
         if( ielpos .eq. -1 ) atmrc(2,4) = atmrc(2,4) + wgt / wg0
         if( ielpos .eq.  1 ) atmrc(3,4) = atmrc(3,4) + wgt / wg0

*-----------------------------------------------------------------------

  140 continue

*-----------------------------------------------------------------------

         erg = erg9
         wgt = wgt9

         uuu = usave(1)
         vvv = usave(2)
         www = usave(3)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine kxrelc(d,f,xv,ireg,imat,rhi)
*                                                                      *
*        generate x-rays and/or auger electrons along path segment d.  *
*        f = energy interpolation factor.                              *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

      dimension xv(3)

*-----------------------------------------------------------------------

         eg09 = eg0
         erg9 = erg
         wgt9 = wgt

         eg0 = erg

*-----------------------------------------------------------------------

            ni = rang() + d * rhi * xnm(mkc)
     &                  * ( (1.-f) * pxr(ngp+nee*(mkc-1) )
     &                         + f * pxr(ngp+1+nee*(mkc-1) ) )

            if( ni . eq. 0 ) return

*-----------------------------------------------------------------------

      do 110 n = 1, ni

            dk  = rang() * d
            dxd = d - dk

            wgt = wgt / xnm(mkc)

*-----------------------------------------------------------------------
*        er=energy edge energy for kxrelc.
*        er=xss(jxs(4,ix)+nxs(4,ix)-1)*1.00001
*        1.0001 ensures round-off does not give spurious results.
*-----------------------------------------------------------------------

            er = edg(mkc) * 1.0001
            ix = iex
            mk = mkc
            ik = ipt

*-----------------------------------------------------------------------
*        integrated tiger series (its) method before its3.0.
*-----------------------------------------------------------------------

         if( nxs(16,lme(3,jmd(1+mkc))) .ne. 3 ) then

               call flugam(er,ix,mk,ik,ireg,imat,xv,dxd)

*-----------------------------------------------------------------------
*        assume the electron has enough energy to excite the atom
*        above the threshold of the highest k shell energy of the
*        highest z atom of the material.
*-----------------------------------------------------------------------

         else if( iphot .eq. 0 .and. eek(mkc) .gt. elc(2) ) then

               jm = jmd(1+mkc)

            do j = jmd(1+mkc), jmd(1+mkc+1)-1

               if( iza(j) .gt. iza(jm) ) jm = j

            end do

               ix  = lme(2,jm)
               iex = ix

               if( er .gt. eg0 ) er = eg0

               call flugam(er,ix,mk,ik,ireg,imat,xv,dxd)

         end if

*-----------------------------------------------------------------------

  110 continue

*-----------------------------------------------------------------------

         erg = erg9
         eg0 = eg09
         wgt = wgt9

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine brmelc(d,qo,f,no,ae,xv,ireg,imat,rhi)
*                                                                      *
*        generate a bremsstrahlung photon along path segment d.        *
*        qo = energy loss per unit path length at ngp for eg0.         *
*        f = energy interpolation factor.                              *
*        no = value of ngp appropriate for eg0.                        *
*        ae = deflection cosine of electron over path segment d.       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        two kinds of bremsstrahlung bias are allowed (not together):  *
*           bnum:  produce bnum times the analog number of photons;    *
*           numb:  produce 1 photon per electron substep.              *
*        el is dimensioned to keep track of first, possibly real event *
*        energy loss, and second and higher order events to keep track *
*        of energy loss only for production purposes.                  *
*        ib = 0/1 / no/yes physical bremsstrahlung photons produced.   *
*        im = 0 for pre-its3.0 physics; = 3 for its3.0 physics.        *
*                                                                      *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for  REDUCTION_COUNTER
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

      dimension usave(3)
      dimension xv(3)
      dimension el(2)

*-----------------------------------------------------------------------

         wgt9 = wgt
         erg9 = erg

         usave(1) = uuu
         usave(2) = vvv
         usave(3) = www

*-----------------------------------------------------------------------
*     decide whether to make a photon.
*-----------------------------------------------------------------------

            pr = d * rhi
     &         * ( (1.-f) * pbr(ngp+nee*(mkc-1))
     &               + f  * pbr(ngp+1+nee*(mkc-1)) )

            im = nxs(16,lme(3,jmd(1+mkc)))
            r1 = rang()
            p1 = exp(-pr)
            ib = 1
            if((im.eq.0.and.pr+r1.lt.1.).or.(im.eq.3.and.r1.lt.p1))ib=0
            if(numb+ib.eq.0) return

*-----------------------------------------------------------------------
*     set up to produce none, one, or more photons.
*-----------------------------------------------------------------------

            na = 1
            nb = 1
            wb = wgt

*-----------------------------------------------------------------------
*     do splitting or roulette if bnum biasing is used.
*-----------------------------------------------------------------------

            if(bnum.eq.1.) goto 10
            if(dbcn(20).ne.0.) goto 10
            nb = abs(bnum)+rang()
            if(nb.eq.0) na = 0
            nb = max(nb,1)
            if(bnum.ne.0.) wb = wgt / abs(bnum)
   10       rb = nb
            w0 = wb
            if(numb.eq.1) wb = wb-wb*p1
            ps = p1

*-----------------------------------------------------------------------
*     nb = number of bremsstrahlung bias photons. nb>1 if bnum>1.
*     note nb loop done outside of gpblcm(i) save because position
*     is sampled for each photon.
*-----------------------------------------------------------------------

      do 240 n1 = 1, nb

*-----------------------------------------------------------------------

            ie = min(n1,2)
            el(ie) = 0.
            if(im.eq.0) goto 16
            p1 = ps
            ap = ps
            fa = 0.

*-----------------------------------------------------------------------
*     prepare to make a photon.
*-----------------------------------------------------------------------

   16 continue

            npa = na
            wgt = wb

*-----------------------------------------------------------------------
*     sample the site and get the electron energy.
*-----------------------------------------------------------------------

            fd  = rang()
            dk  = fd * d
            dxd = d - dk

            es = eg0 - qo * dk - el(ie)

            if( es .lt. elc(3) ) goto 200

            m = no

   40       if( es .ge. eee(m+1) ) goto 50
            m = m + 1
            goto 40

*-----------------------------------------------------------------------
*     sample an analog photon for correct electron track energy loss.
*-----------------------------------------------------------------------

   50       rn = rang()
            pk = qpol(rn,eba(1,m+nee*(mkc-1)),rkt,ntop)
            pa = qpol(rn,eba(1,m+1+nee*(mkc-1)),rkt,ntop)
            pj = pa
     &         + (es-eee(m+1))*(pk-pa)/(eee(m)-eee(m+1))

            erg = pj * es

*-----------------------------------------------------------------------
*     always accumulate electron track energy loss.
*-----------------------------------------------------------------------
*        accumulate radiative energy loss for first photon or all
*        of them for bnum<0.
*-----------------------------------------------------------------------

            el(ie) = el(ie) + erg

            if( n1 .gt.1 .and. bnum .lt. 0. ) el(1) = el(1) + erg

*-----------------------------------------------------------------------
*     return if no photon can be followed.
*-----------------------------------------------------------------------

            if( iphot .ne. 0 .or.
     &          es .lt. elc(2) .or.
     &          npa .eq. 0 ) goto 200

*-----------------------------------------------------------------------
*     sample a biased photon energy, if bbrem is used.
*-----------------------------------------------------------------------

            if( mbi(mkc) .eq. 0 ) goto 100

            rn = rang()
            pk = qpol(rn,ebd(1,m+nee*(mkc-1)),rkt,ntop)
            pa = qpol(rn,ebd(1,m+1+nee*(mkc-1)),rkt,ntop)
            pj = pa
     &         + (es-eee(m+1))*(pk-pa)/(eee(m)-eee(m+1))

         do 60 j1 = 2, ntop - 1
   60       if(pk.le.rkt(j1)) goto 70

   70    do 80 j = 2, ntop - 1
   80       if(pa.le.rkt(j)) goto 90

   90       pb = fst(m+1+nee*(mkc-1)) * bbrem(j)
     &         + (es-eee(m+1)) * (fst(m+nee*(mkc-1))
     *         * bbrem(j1)
     *         - fst(m+1+nee*(mkc-1))*bbrem(j))
     &         / (eee(m)-eee(m+1))

            wgt = wgt * pb
            erg = pj * es

*-----------------------------------------------------------------------
*     return if this photon cannot be followed.
*-----------------------------------------------------------------------

  100       if( erg .lt. elc(2) ) goto 200

*-----------------------------------------------------------------------
*     do splitting or roulette if bnum biasing is used.
*-----------------------------------------------------------------------

            if( bnum .eq. 1. ) goto 110
            if( dbcn(20) .eq. 0. ) goto 110

            npa = bnum + rang()

            if( npa .eq. 0 ) goto 200

            wgt = wgt / bnum

  110       ipt = 2

*-----------------------------------------------------------------------
*     sample the photon direction.
*-----------------------------------------------------------------------

            deb = sqrt(es*(es+2.*gpt(3))) / (es+gpt(3))
            rn  = rang()

*-----------------------------------------------------------------------
*     use the detailed angular distribution model ...
*-----------------------------------------------------------------------

      if( ibad .ne. 0 ) goto 160

            k = (m-1)/nstp+1
            ek = eee((k-1)*nstp+1)-es
            e1 = es-eee(k*nstp+1)

         do 140 j=2,nwng
  140       if(pj.ge.rka(j)) goto 150

  150       rj = rka(j-1)-pj
            sj = pj-rka(j)
            l = rn*(mpng-1)
            pl = rn*(mpng-1)-l
            ql = 1.-pl
            i = l+1
            mc = k+(nee/nstp+1)*(mkc-1)
            t1 = ql*sj*e1*ech(i,j-1,mc)+pl*rj*ek*ech(i+1,j,mc+1)+
     &       pl*sj*e1*ech(i+1,j-1,mc)+ql*rj*e1*ech(i,j,mc)+
     &       ql*sj*ek*ech(i,j-1,mc+1)+ql*rj*ek*ech(i,j,mc+1)+
     &       pl*sj*ek*ech(i+1,j-1,mc+1)+pl*rj*e1*ech(i+1,j,mc)
            if(im.eq.0)
     &      cx = 4.*t1/((rka(j-1)-rka(j))*eee((k-1)*nstp+1))

            if(im.eq.3)
     &      cx = 2.*t1/((rka(j-1)-rka(j))
     &         * (eee((k-1)*nstp+1)-eee(k*nstp+1)))

            am = (cx-1.+deb) / (cx*deb+1.-deb)

            goto 170

*-----------------------------------------------------------------------
*     ... or the simple model.
*-----------------------------------------------------------------------

  160       am = (2.*rn-1.-deb) / (2.*rn*deb-1.-deb)

*-----------------------------------------------------------------------
*     sample the photon direction relative to an intermediate
*     electron direction.  change uold for tallyd, dxtran.
*-----------------------------------------------------------------------

  170       u = uold(1)
            v = uold(2)
            w = uold(3)

            call defelc(am,dk,fd,ae,2,
     &                 uuu,vvv,www,uold(1),uold(2),uold(3))

*-----------------------------------------------------------------------
*     bank the photon.
*-----------------------------------------------------------------------

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 4

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 0
                     jclusts(3,nclsts) = 14
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = 0
                     jclusts(6,nclsts) = 0
                     jclusts(7,nclsts) = 22
                     jclusts(8,nclsts) = 0

                     pra  = erg / 1000.0

                     pxl = pra * uuu
                     pyl = pra * vvv
                     pzl = pra * www

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = pxl
                     qclusts(2,nclsts)  = pyl
                     qclusts(3,nclsts)  = pzl
                     qclusts(4,nclsts)  = pra
                     qclusts(5,nclsts)  = 0.0
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = erg
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = xv(1) * dxd
                     qclusts(11,nclsts) = xv(2) * dxd
                     qclusts(12,nclsts) = xv(3) * dxd
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
         if( ielpos .eq. -1 ) atmrc(2,3) = atmrc(2,3) + wgt / wg0
         if( ielpos .eq.  1 ) atmrc(3,3) = atmrc(3,3) + wgt / wg0


*-----------------------------------------------------------------------

            if( dbcn(20) .ne. 0. .and. bnum .gt. 1. ) wgt = wgt * bnum

            uold(1) = u
            uold(2) = v
            uold(3) = w

*-----------------------------------------------------------------------

  200 continue

            erg = erg9
            wgt = wgt9

            uuu = usave(1)
            vvv = usave(2)
            www = usave(3)

*-----------------------------------------------------------------------
*     reduce the electron energy by the amount of the radiation loss.
*-----------------------------------------------------------------------

            if( im .eq. 0 .or. ib .eq. 0 ) goto 240

            if( numb .eq. 1 ) wb = w0
            fa = fa + 1.
            p1 = p1 * pr / fa
            ap = ap + p1

            if( r1 .gt. ap ) goto 16

*-----------------------------------------------------------------------

  240 continue

*-----------------------------------------------------------------------

            if( ib .eq. 0 ) return

            if( dbcn(20) .eq. 0. .and. bnum .lt. 0. )
     &      el(1) = el(1) / rb
                     aevts(iaevt+49,ireg) =
     &               aevts(iaevt+49,ireg) + wgt * el(1)
                     bevts(ibevt+49,imat) =
     &               bevts(ibevt+49,imat) + wgt * el(1)
*-----------------------------------------------------------------------

            erg = erg9
            wgt = wgt9

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      block data landau
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      parameter (mlam=5001)
      common /lancom/ eqlm(mlam)

*-----------------------------------------------------------------------
*     the array eqlm(5001) gives the boundaries of 5000 equiprobable
*     lambda bins.  lambda is the dimensionless independent variable
*     in the landau theory of electron energy loss straggling.
*     the tabular data cover the range labmda = (-4,100).
*     the range (lambda.gt.100.) is sampled analytically.
*-----------------------------------------------------------------------

      data (eqlm(i),i=   1,  48)/-4.00000d+00,-2.88786d+00,-2.78578d+00,
     1 -2.72011d+00,-2.67110d+00,-2.63151d+00,-2.59805d+00,-2.56895d+00,
     2 -2.54309d+00,-2.51977d+00,-2.49847d+00,-2.47884d+00,-2.46060d+00,
     3 -2.44356d+00,-2.42753d+00,-2.41240d+00,-2.39805d+00,-2.38440d+00,
     4 -2.37137d+00,-2.35891d+00,-2.34695d+00,-2.33545d+00,-2.32437d+00,
     5 -2.31368d+00,-2.30335d+00,-2.29335d+00,-2.28365d+00,-2.27424d+00,
     6 -2.26509d+00,-2.25619d+00,-2.24753d+00,-2.23908d+00,-2.23084d+00,
     7 -2.22279d+00,-2.21493d+00,-2.20725d+00,-2.19972d+00,-2.19236d+00,
     8 -2.18514d+00,-2.17807d+00,-2.17113d+00,-2.16432d+00,-2.15763d+00,
     9 -2.15106d+00,-2.14461d+00,-2.13826d+00,-2.13202d+00,-2.12588d+00/
      data (eqlm(i),i=  49,  96)/-2.11983d+00,-2.11388d+00,-2.10802d+00,
     1 -2.10224d+00,-2.09654d+00,-2.09093d+00,-2.08539d+00,-2.07993d+00,
     2 -2.07454d+00,-2.06922d+00,-2.06397d+00,-2.05878d+00,-2.05366d+00,
     3 -2.04860d+00,-2.04360d+00,-2.03865d+00,-2.03377d+00,-2.02893d+00,
     4 -2.02416d+00,-2.01943d+00,-2.01475d+00,-2.01013d+00,-2.00555d+00,
     5 -2.00101d+00,-1.99653d+00,-1.99208d+00,-1.98768d+00,-1.98332d+00,
     6 -1.97901d+00,-1.97473d+00,-1.97049d+00,-1.96630d+00,-1.96213d+00,
     7 -1.95801d+00,-1.95392d+00,-1.94986d+00,-1.94584d+00,-1.94186d+00,
     8 -1.93790d+00,-1.93398d+00,-1.93009d+00,-1.92623d+00,-1.92240d+00,
     9 -1.91860d+00,-1.91483d+00,-1.91109d+00,-1.90737d+00,-1.90368d+00/
      data (eqlm(i),i=  97, 144)/-1.90002d+00,-1.89639d+00,-1.89278d+00,
     1 -1.88920d+00,-1.88564d+00,-1.88210d+00,-1.87859d+00,-1.87510d+00,
     2 -1.87164d+00,-1.86820d+00,-1.86478d+00,-1.86138d+00,-1.85801d+00,
     3 -1.85465d+00,-1.85132d+00,-1.84801d+00,-1.84471d+00,-1.84144d+00,
     4 -1.83819d+00,-1.83495d+00,-1.83174d+00,-1.82854d+00,-1.82536d+00,
     5 -1.82220d+00,-1.81906d+00,-1.81593d+00,-1.81283d+00,-1.80974d+00,
     6 -1.80666d+00,-1.80360d+00,-1.80056d+00,-1.79754d+00,-1.79453d+00,
     7 -1.79153d+00,-1.78855d+00,-1.78559d+00,-1.78264d+00,-1.77971d+00,
     8 -1.77679d+00,-1.77388d+00,-1.77099d+00,-1.76812d+00,-1.76525d+00,
     9 -1.76240d+00,-1.75957d+00,-1.75674d+00,-1.75393d+00,-1.75114d+00/
      data (eqlm(i),i= 145, 192)/-1.74835d+00,-1.74558d+00,-1.74282d+00,
     1 -1.74007d+00,-1.73734d+00,-1.73462d+00,-1.73190d+00,-1.72920d+00,
     2 -1.72652d+00,-1.72384d+00,-1.72117d+00,-1.71852d+00,-1.71588d+00,
     3 -1.71324d+00,-1.71062d+00,-1.70801d+00,-1.70541d+00,-1.70282d+00,
     4 -1.70024d+00,-1.69767d+00,-1.69511d+00,-1.69256d+00,-1.69002d+00,
     5 -1.68749d+00,-1.68497d+00,-1.68246d+00,-1.67996d+00,-1.67746d+00,
     6 -1.67498d+00,-1.67250d+00,-1.67004d+00,-1.66758d+00,-1.66513d+00,
     7 -1.66270d+00,-1.66027d+00,-1.65784d+00,-1.65543d+00,-1.65303d+00,
     8 -1.65063d+00,-1.64824d+00,-1.64586d+00,-1.64349d+00,-1.64112d+00,
     9 -1.63877d+00,-1.63642d+00,-1.63408d+00,-1.63174d+00,-1.62942d+00/
      data (eqlm(i),i= 193, 240)/-1.62710d+00,-1.62479d+00,-1.62248d+00,
     1 -1.62019d+00,-1.61790d+00,-1.61562d+00,-1.61334d+00,-1.61108d+00,
     2 -1.60881d+00,-1.60656d+00,-1.60431d+00,-1.60207d+00,-1.59984d+00,
     3 -1.59761d+00,-1.59539d+00,-1.59318d+00,-1.59097d+00,-1.58877d+00,
     4 -1.58658d+00,-1.58439d+00,-1.58221d+00,-1.58003d+00,-1.57786d+00,
     5 -1.57570d+00,-1.57354d+00,-1.57139d+00,-1.56924d+00,-1.56710d+00,
     6 -1.56497d+00,-1.56284d+00,-1.56072d+00,-1.55860d+00,-1.55649d+00,
     7 -1.55439d+00,-1.55229d+00,-1.55019d+00,-1.54811d+00,-1.54602d+00,
     8 -1.54394d+00,-1.54187d+00,-1.53981d+00,-1.53774d+00,-1.53569d+00,
     9 -1.53363d+00,-1.53159d+00,-1.52955d+00,-1.52751d+00,-1.52548d+00/
      data (eqlm(i),i= 241, 288)/-1.52345d+00,-1.52143d+00,-1.51942d+00,
     1 -1.51740d+00,-1.51540d+00,-1.51339d+00,-1.51140d+00,-1.50940d+00,
     2 -1.50742d+00,-1.50543d+00,-1.50346d+00,-1.50148d+00,-1.49951d+00,
     3 -1.49755d+00,-1.49559d+00,-1.49363d+00,-1.49168d+00,-1.48973d+00,
     4 -1.48779d+00,-1.48585d+00,-1.48392d+00,-1.48199d+00,-1.48006d+00,
     5 -1.47814d+00,-1.47622d+00,-1.47431d+00,-1.47240d+00,-1.47049d+00,
     6 -1.46859d+00,-1.46669d+00,-1.46480d+00,-1.46291d+00,-1.46103d+00,
     7 -1.45914d+00,-1.45727d+00,-1.45539d+00,-1.45352d+00,-1.45166d+00,
     8 -1.44980d+00,-1.44794d+00,-1.44608d+00,-1.44423d+00,-1.44238d+00,
     9 -1.44054d+00,-1.43870d+00,-1.43687d+00,-1.43503d+00,-1.43320d+00/
      data (eqlm(i),i= 289, 336)/-1.43138d+00,-1.42956d+00,-1.42774d+00,
     1 -1.42592d+00,-1.42411d+00,-1.42230d+00,-1.42050d+00,-1.41870d+00,
     2 -1.41690d+00,-1.41510d+00,-1.41331d+00,-1.41152d+00,-1.40974d+00,
     3 -1.40796d+00,-1.40618d+00,-1.40440d+00,-1.40263d+00,-1.40086d+00,
     4 -1.39910d+00,-1.39734d+00,-1.39558d+00,-1.39382d+00,-1.39207d+00,
     5 -1.39032d+00,-1.38857d+00,-1.38683d+00,-1.38508d+00,-1.38335d+00,
     6 -1.38161d+00,-1.37988d+00,-1.37815d+00,-1.37642d+00,-1.37470d+00,
     7 -1.37298d+00,-1.37126d+00,-1.36955d+00,-1.36784d+00,-1.36613d+00,
     8 -1.36442d+00,-1.36272d+00,-1.36102d+00,-1.35932d+00,-1.35762d+00,
     9 -1.35593d+00,-1.35424d+00,-1.35255d+00,-1.35087d+00,-1.34919d+00/
      data (eqlm(i),i= 337, 384)/-1.34751d+00,-1.34583d+00,-1.34416d+00,
     1 -1.34249d+00,-1.34082d+00,-1.33915d+00,-1.33749d+00,-1.33583d+00,
     2 -1.33417d+00,-1.33251d+00,-1.33086d+00,-1.32921d+00,-1.32756d+00,
     3 -1.32592d+00,-1.32427d+00,-1.32263d+00,-1.32099d+00,-1.31936d+00,
     4 -1.31772d+00,-1.31609d+00,-1.31446d+00,-1.31283d+00,-1.31121d+00,
     5 -1.30959d+00,-1.30797d+00,-1.30635d+00,-1.30473d+00,-1.30312d+00,
     6 -1.30151d+00,-1.29990d+00,-1.29830d+00,-1.29669d+00,-1.29509d+00,
     7 -1.29349d+00,-1.29189d+00,-1.29030d+00,-1.28871d+00,-1.28712d+00,
     8 -1.28553d+00,-1.28394d+00,-1.28236d+00,-1.28077d+00,-1.27919d+00,
     9 -1.27762d+00,-1.27604d+00,-1.27447d+00,-1.27290d+00,-1.27133d+00/
      data (eqlm(i),i= 385, 432)/-1.26976d+00,-1.26819d+00,-1.26663d+00,
     1 -1.26507d+00,-1.26351d+00,-1.26195d+00,-1.26040d+00,-1.25884d+00,
     2 -1.25729d+00,-1.25574d+00,-1.25419d+00,-1.25265d+00,-1.25110d+00,
     3 -1.24956d+00,-1.24802d+00,-1.24648d+00,-1.24495d+00,-1.24341d+00,
     4 -1.24188d+00,-1.24035d+00,-1.23882d+00,-1.23729d+00,-1.23577d+00,
     5 -1.23425d+00,-1.23273d+00,-1.23121d+00,-1.22969d+00,-1.22817d+00,
     6 -1.22666d+00,-1.22515d+00,-1.22364d+00,-1.22213d+00,-1.22062d+00,
     7 -1.21911d+00,-1.21761d+00,-1.21611d+00,-1.21461d+00,-1.21311d+00,
     8 -1.21161d+00,-1.21012d+00,-1.20862d+00,-1.20713d+00,-1.20564d+00,
     9 -1.20415d+00,-1.20267d+00,-1.20118d+00,-1.19970d+00,-1.19822d+00/
      data (eqlm(i),i= 433, 480)/-1.19674d+00,-1.19526d+00,-1.19378d+00,
     1 -1.19231d+00,-1.19083d+00,-1.18936d+00,-1.18789d+00,-1.18642d+00,
     2 -1.18495d+00,-1.18349d+00,-1.18202d+00,-1.18056d+00,-1.17910d+00,
     3 -1.17764d+00,-1.17618d+00,-1.17472d+00,-1.17327d+00,-1.17181d+00,
     4 -1.17036d+00,-1.16891d+00,-1.16746d+00,-1.16601d+00,-1.16457d+00,
     5 -1.16312d+00,-1.16168d+00,-1.16024d+00,-1.15880d+00,-1.15736d+00,
     6 -1.15592d+00,-1.15449d+00,-1.15305d+00,-1.15162d+00,-1.15019d+00,
     7 -1.14876d+00,-1.14733d+00,-1.14590d+00,-1.14447d+00,-1.14305d+00,
     8 -1.14162d+00,-1.14020d+00,-1.13878d+00,-1.13736d+00,-1.13594d+00,
     9 -1.13453d+00,-1.13311d+00,-1.13170d+00,-1.13028d+00,-1.12887d+00/
      data (eqlm(i),i= 481, 528)/-1.12746d+00,-1.12605d+00,-1.12465d+00,
     1 -1.12324d+00,-1.12184d+00,-1.12043d+00,-1.11903d+00,-1.11763d+00,
     2 -1.11623d+00,-1.11483d+00,-1.11343d+00,-1.11204d+00,-1.11064d+00,
     3 -1.10925d+00,-1.10786d+00,-1.10646d+00,-1.10507d+00,-1.10369d+00,
     4 -1.10230d+00,-1.10091d+00,-1.09953d+00,-1.09814d+00,-1.09676d+00,
     5 -1.09538d+00,-1.09400d+00,-1.09262d+00,-1.09124d+00,-1.08987d+00,
     6 -1.08849d+00,-1.08712d+00,-1.08574d+00,-1.08437d+00,-1.08300d+00,
     7 -1.08163d+00,-1.08026d+00,-1.07889d+00,-1.07753d+00,-1.07616d+00,
     8 -1.07480d+00,-1.07344d+00,-1.07207d+00,-1.07071d+00,-1.06935d+00,
     9 -1.06800d+00,-1.06664d+00,-1.06528d+00,-1.06393d+00,-1.06257d+00/
      data (eqlm(i),i= 529, 576)/-1.06122d+00,-1.05987d+00,-1.05852d+00,
     1 -1.05717d+00,-1.05582d+00,-1.05447d+00,-1.05312d+00,-1.05178d+00,
     2 -1.05043d+00,-1.04909d+00,-1.04775d+00,-1.04641d+00,-1.04507d+00,
     3 -1.04373d+00,-1.04239d+00,-1.04105d+00,-1.03971d+00,-1.03838d+00,
     4 -1.03704d+00,-1.03571d+00,-1.03438d+00,-1.03305d+00,-1.03172d+00,
     5 -1.03039d+00,-1.02906d+00,-1.02773d+00,-1.02640d+00,-1.02508d+00,
     6 -1.02375d+00,-1.02243d+00,-1.02111d+00,-1.01979d+00,-1.01846d+00,
     7 -1.01715d+00,-1.01583d+00,-1.01451d+00,-1.01319d+00,-1.01188d+00,
     8 -1.01056d+00,-1.00925d+00,-1.00793d+00,-1.00662d+00,-1.00531d+00,
     9 -1.00400d+00,-1.00269d+00,-1.00138d+00,-1.00007d+00,-9.98766d-01/
      data (eqlm(i),i= 577, 624)/-9.97460d-01,-9.96155d-01,-9.94850d-01,
     1 -9.93547d-01,-9.92245d-01,-9.90943d-01,-9.89642d-01,-9.88343d-01,
     2 -9.87044d-01,-9.85746d-01,-9.84448d-01,-9.83152d-01,-9.81857d-01,
     3 -9.80562d-01,-9.79268d-01,-9.77975d-01,-9.76683d-01,-9.75392d-01,
     4 -9.74101d-01,-9.72812d-01,-9.71523d-01,-9.70235d-01,-9.68948d-01,
     5 -9.67661d-01,-9.66376d-01,-9.65091d-01,-9.63807d-01,-9.62524d-01,
     6 -9.61242d-01,-9.59960d-01,-9.58679d-01,-9.57399d-01,-9.56120d-01,
     7 -9.54842d-01,-9.53564d-01,-9.52287d-01,-9.51011d-01,-9.49736d-01,
     8 -9.48462d-01,-9.47188d-01,-9.45915d-01,-9.44643d-01,-9.43372d-01,
     9 -9.42101d-01,-9.40831d-01,-9.39562d-01,-9.38294d-01,-9.37026d-01/
      data (eqlm(i),i= 625, 672)/-9.35759d-01,-9.34493d-01,-9.33228d-01,
     1 -9.31963d-01,-9.30699d-01,-9.29436d-01,-9.28174d-01,-9.26912d-01,
     2 -9.25651d-01,-9.24391d-01,-9.23131d-01,-9.21872d-01,-9.20614d-01,
     3 -9.19357d-01,-9.18100d-01,-9.16844d-01,-9.15589d-01,-9.14335d-01,
     4 -9.13081d-01,-9.11828d-01,-9.10575d-01,-9.09323d-01,-9.08072d-01,
     5 -9.06822d-01,-9.05572d-01,-9.04323d-01,-9.03075d-01,-9.01827d-01,
     6 -9.00580d-01,-8.99334d-01,-8.98088d-01,-8.96843d-01,-8.95599d-01,
     7 -8.94355d-01,-8.93112d-01,-8.91870d-01,-8.90628d-01,-8.89387d-01,
     8 -8.88147d-01,-8.86907d-01,-8.85668d-01,-8.84430d-01,-8.83192d-01,
     9 -8.81955d-01,-8.80718d-01,-8.79482d-01,-8.78247d-01,-8.77012d-01/
      data (eqlm(i),i= 673, 720)/-8.75778d-01,-8.74545d-01,-8.73312d-01,
     1 -8.72080d-01,-8.70849d-01,-8.69618d-01,-8.68387d-01,-8.67158d-01,
     2 -8.65929d-01,-8.64700d-01,-8.63472d-01,-8.62245d-01,-8.61018d-01,
     3 -8.59792d-01,-8.58567d-01,-8.57342d-01,-8.56117d-01,-8.54894d-01,
     4 -8.53670d-01,-8.52448d-01,-8.51226d-01,-8.50005d-01,-8.48784d-01,
     5 -8.47563d-01,-8.46344d-01,-8.45125d-01,-8.43906d-01,-8.42688d-01,
     6 -8.41471d-01,-8.40254d-01,-8.39038d-01,-8.37822d-01,-8.36607d-01,
     7 -8.35392d-01,-8.34178d-01,-8.32965d-01,-8.31752d-01,-8.30539d-01,
     8 -8.29327d-01,-8.28116d-01,-8.26905d-01,-8.25695d-01,-8.24486d-01,
     9 -8.23276d-01,-8.22068d-01,-8.20860d-01,-8.19652d-01,-8.18445d-01/
      data (eqlm(i),i= 721, 768)/-8.17239d-01,-8.16033d-01,-8.14827d-01,
     1 -8.13622d-01,-8.12418d-01,-8.11214d-01,-8.10011d-01,-8.08808d-01,
     2 -8.07606d-01,-8.06404d-01,-8.05203d-01,-8.04002d-01,-8.02802d-01,
     3 -8.01602d-01,-8.00403d-01,-7.99204d-01,-7.98005d-01,-7.96808d-01,
     4 -7.95610d-01,-7.94414d-01,-7.93217d-01,-7.92022d-01,-7.90826d-01,
     5 -7.89631d-01,-7.88437d-01,-7.87243d-01,-7.86050d-01,-7.84857d-01,
     6 -7.83665d-01,-7.82473d-01,-7.81281d-01,-7.80090d-01,-7.78900d-01,
     7 -7.77709d-01,-7.76520d-01,-7.75331d-01,-7.74142d-01,-7.72954d-01,
     8 -7.71766d-01,-7.70579d-01,-7.69392d-01,-7.68205d-01,-7.67020d-01,
     9 -7.65834d-01,-7.64649d-01,-7.63464d-01,-7.62280d-01,-7.61097d-01/
      data (eqlm(i),i= 769, 816)/-7.59913d-01,-7.58730d-01,-7.57548d-01,
     1 -7.56366d-01,-7.55185d-01,-7.54003d-01,-7.52823d-01,-7.51643d-01,
     2 -7.50463d-01,-7.49283d-01,-7.48104d-01,-7.46926d-01,-7.45748d-01,
     3 -7.44570d-01,-7.43393d-01,-7.42216d-01,-7.41040d-01,-7.39864d-01,
     4 -7.38688d-01,-7.37513d-01,-7.36338d-01,-7.35164d-01,-7.33990d-01,
     5 -7.32816d-01,-7.31643d-01,-7.30470d-01,-7.29298d-01,-7.28126d-01,
     6 -7.26955d-01,-7.25784d-01,-7.24613d-01,-7.23443d-01,-7.22273d-01,
     7 -7.21103d-01,-7.19934d-01,-7.18765d-01,-7.17597d-01,-7.16429d-01,
     8 -7.15261d-01,-7.14094d-01,-7.12927d-01,-7.11761d-01,-7.10595d-01,
     9 -7.09429d-01,-7.08264d-01,-7.07099d-01,-7.05935d-01,-7.04771d-01/
      data (eqlm(i),i= 817, 864)/-7.03607d-01,-7.02443d-01,-7.01280d-01,
     1 -7.00118d-01,-6.98955d-01,-6.97793d-01,-6.96632d-01,-6.95471d-01,
     2 -6.94310d-01,-6.93149d-01,-6.91989d-01,-6.90829d-01,-6.89670d-01,
     3 -6.88511d-01,-6.87352d-01,-6.86194d-01,-6.85036d-01,-6.83878d-01,
     4 -6.82721d-01,-6.81564d-01,-6.80407d-01,-6.79251d-01,-6.78095d-01,
     5 -6.76939d-01,-6.75784d-01,-6.74629d-01,-6.73475d-01,-6.72320d-01,
     6 -6.71166d-01,-6.70013d-01,-6.68860d-01,-6.67707d-01,-6.66554d-01,
     7 -6.65402d-01,-6.64250d-01,-6.63098d-01,-6.61947d-01,-6.60796d-01,
     8 -6.59646d-01,-6.58495d-01,-6.57345d-01,-6.56196d-01,-6.55046d-01,
     9 -6.53897d-01,-6.52748d-01,-6.51600d-01,-6.50452d-01,-6.49304d-01/
      data (eqlm(i),i= 865, 912)/-6.48157d-01,-6.47009d-01,-6.45863d-01,
     1 -6.44716d-01,-6.43570d-01,-6.42424d-01,-6.41278d-01,-6.40133d-01,
     2 -6.38988d-01,-6.37843d-01,-6.36699d-01,-6.35554d-01,-6.34411d-01,
     3 -6.33267d-01,-6.32124d-01,-6.30981d-01,-6.29838d-01,-6.28696d-01,
     4 -6.27553d-01,-6.26412d-01,-6.25270d-01,-6.24129d-01,-6.22988d-01,
     5 -6.21847d-01,-6.20707d-01,-6.19567d-01,-6.18427d-01,-6.17287d-01,
     6 -6.16148d-01,-6.15009d-01,-6.13870d-01,-6.12732d-01,-6.11594d-01,
     7 -6.10456d-01,-6.09318d-01,-6.08181d-01,-6.07044d-01,-6.05907d-01,
     8 -6.04770d-01,-6.03634d-01,-6.02498d-01,-6.01362d-01,-6.00227d-01,
     9 -5.99091d-01,-5.97956d-01,-5.96822d-01,-5.95687d-01,-5.94553d-01/
      data (eqlm(i),i= 913, 960)/-5.93419d-01,-5.92285d-01,-5.91152d-01,
     1 -5.90019d-01,-5.88886d-01,-5.87753d-01,-5.86621d-01,-5.85489d-01,
     2 -5.84357d-01,-5.83225d-01,-5.82093d-01,-5.80962d-01,-5.79831d-01,
     3 -5.78701d-01,-5.77570d-01,-5.76440d-01,-5.75310d-01,-5.74180d-01,
     4 -5.73051d-01,-5.71921d-01,-5.70792d-01,-5.69664d-01,-5.68535d-01,
     5 -5.67407d-01,-5.66279d-01,-5.65151d-01,-5.64023d-01,-5.62896d-01,
     6 -5.61768d-01,-5.60641d-01,-5.59515d-01,-5.58388d-01,-5.57262d-01,
     7 -5.56136d-01,-5.55010d-01,-5.53884d-01,-5.52759d-01,-5.51633d-01,
     8 -5.50508d-01,-5.49384d-01,-5.48259d-01,-5.47135d-01,-5.46011d-01,
     9 -5.44887d-01,-5.43763d-01,-5.42639d-01,-5.41516d-01,-5.40393d-01/
      data (eqlm(i),i= 961,1008)/-5.39270d-01,-5.38147d-01,-5.37025d-01,
     1 -5.35903d-01,-5.34781d-01,-5.33659d-01,-5.32537d-01,-5.31416d-01,
     2 -5.30294d-01,-5.29173d-01,-5.28053d-01,-5.26932d-01,-5.25811d-01,
     3 -5.24691d-01,-5.23571d-01,-5.22451d-01,-5.21332d-01,-5.20212d-01,
     4 -5.19093d-01,-5.17974d-01,-5.16855d-01,-5.15736d-01,-5.14618d-01,
     5 -5.13499d-01,-5.12381d-01,-5.11263d-01,-5.10145d-01,-5.09028d-01,
     6 -5.07910d-01,-5.06793d-01,-5.05676d-01,-5.04559d-01,-5.03443d-01,
     7 -5.02326d-01,-5.01210d-01,-5.00094d-01,-4.98978d-01,-4.97862d-01,
     8 -4.96746d-01,-4.95631d-01,-4.94516d-01,-4.93400d-01,-4.92285d-01,
     9 -4.91171d-01,-4.90056d-01,-4.88942d-01,-4.87827d-01,-4.86713d-01/
      data (eqlm(i),i=1009,1056)/-4.85599d-01,-4.84486d-01,-4.83372d-01,
     1 -4.82258d-01,-4.81145d-01,-4.80032d-01,-4.78919d-01,-4.77806d-01,
     2 -4.76694d-01,-4.75581d-01,-4.74469d-01,-4.73357d-01,-4.72245d-01,
     3 -4.71133d-01,-4.70021d-01,-4.68909d-01,-4.67798d-01,-4.66687d-01,
     4 -4.65576d-01,-4.64465d-01,-4.63354d-01,-4.62243d-01,-4.61133d-01,
     5 -4.60022d-01,-4.58912d-01,-4.57802d-01,-4.56692d-01,-4.55582d-01,
     6 -4.54472d-01,-4.53363d-01,-4.52253d-01,-4.51144d-01,-4.50035d-01,
     7 -4.48926d-01,-4.47817d-01,-4.46709d-01,-4.45600d-01,-4.44492d-01,
     8 -4.43383d-01,-4.42275d-01,-4.41167d-01,-4.40059d-01,-4.38951d-01,
     9 -4.37844d-01,-4.36736d-01,-4.35629d-01,-4.34522d-01,-4.33414d-01/
      data (eqlm(i),i=1057,1104)/-4.32307d-01,-4.31201d-01,-4.30094d-01,
     1 -4.28987d-01,-4.27881d-01,-4.26774d-01,-4.25668d-01,-4.24562d-01,
     2 -4.23456d-01,-4.22350d-01,-4.21244d-01,-4.20139d-01,-4.19033d-01,
     3 -4.17928d-01,-4.16822d-01,-4.15717d-01,-4.14612d-01,-4.13507d-01,
     4 -4.12402d-01,-4.11297d-01,-4.10193d-01,-4.09088d-01,-4.07984d-01,
     5 -4.06879d-01,-4.05775d-01,-4.04671d-01,-4.03567d-01,-4.02463d-01,
     6 -4.01359d-01,-4.00256d-01,-3.99152d-01,-3.98049d-01,-3.96945d-01,
     7 -3.95842d-01,-3.94739d-01,-3.93636d-01,-3.92533d-01,-3.91430d-01,
     8 -3.90327d-01,-3.89225d-01,-3.88122d-01,-3.87020d-01,-3.85917d-01,
     9 -3.84815d-01,-3.83713d-01,-3.82611d-01,-3.81509d-01,-3.80407d-01/
      data (eqlm(i),i=1105,1152)/-3.79305d-01,-3.78204d-01,-3.77102d-01,
     1 -3.76000d-01,-3.74899d-01,-3.73798d-01,-3.72696d-01,-3.71595d-01,
     2 -3.70494d-01,-3.69393d-01,-3.68292d-01,-3.67191d-01,-3.66091d-01,
     3 -3.64990d-01,-3.63889d-01,-3.62789d-01,-3.61689d-01,-3.60588d-01,
     4 -3.59488d-01,-3.58388d-01,-3.57288d-01,-3.56188d-01,-3.55088d-01,
     5 -3.53988d-01,-3.52888d-01,-3.51788d-01,-3.50688d-01,-3.49589d-01,
     6 -3.48489d-01,-3.47390d-01,-3.46290d-01,-3.45191d-01,-3.44092d-01,
     7 -3.42992d-01,-3.41893d-01,-3.40794d-01,-3.39695d-01,-3.38596d-01,
     8 -3.37497d-01,-3.36399d-01,-3.35300d-01,-3.34201d-01,-3.33103d-01,
     9 -3.32004d-01,-3.30906d-01,-3.29807d-01,-3.28709d-01,-3.27610d-01/
      data (eqlm(i),i=1153,1200)/-3.26512d-01,-3.25414d-01,-3.24316d-01,
     1 -3.23218d-01,-3.22120d-01,-3.21022d-01,-3.19924d-01,-3.18826d-01,
     2 -3.17728d-01,-3.16630d-01,-3.15533d-01,-3.14435d-01,-3.13337d-01,
     3 -3.12240d-01,-3.11142d-01,-3.10045d-01,-3.08948d-01,-3.07850d-01,
     4 -3.06753d-01,-3.05656d-01,-3.04559d-01,-3.03461d-01,-3.02364d-01,
     5 -3.01267d-01,-3.00170d-01,-2.99073d-01,-2.97976d-01,-2.96879d-01,
     6 -2.95782d-01,-2.94686d-01,-2.93589d-01,-2.92492d-01,-2.91395d-01,
     7 -2.90299d-01,-2.89202d-01,-2.88106d-01,-2.87009d-01,-2.85913d-01,
     8 -2.84816d-01,-2.83720d-01,-2.82623d-01,-2.81527d-01,-2.80430d-01,
     9 -2.79334d-01,-2.78238d-01,-2.77142d-01,-2.76045d-01,-2.74949d-01/
      data (eqlm(i),i=1201,1248)/-2.73853d-01,-2.72757d-01,-2.71661d-01,
     1 -2.70564d-01,-2.69468d-01,-2.68372d-01,-2.67276d-01,-2.66180d-01,
     2 -2.65084d-01,-2.63988d-01,-2.62892d-01,-2.61797d-01,-2.60701d-01,
     3 -2.59605d-01,-2.58509d-01,-2.57413d-01,-2.56317d-01,-2.55221d-01,
     4 -2.54126d-01,-2.53030d-01,-2.51934d-01,-2.50838d-01,-2.49743d-01,
     5 -2.48647d-01,-2.47551d-01,-2.46456d-01,-2.45360d-01,-2.44264d-01,
     6 -2.43169d-01,-2.42073d-01,-2.40977d-01,-2.39882d-01,-2.38786d-01,
     7 -2.37691d-01,-2.36595d-01,-2.35499d-01,-2.34404d-01,-2.33308d-01,
     8 -2.32213d-01,-2.31117d-01,-2.30022d-01,-2.28926d-01,-2.27830d-01,
     9 -2.26735d-01,-2.25639d-01,-2.24544d-01,-2.23448d-01,-2.22353d-01/
      data (eqlm(i),i=1249,1296)/-2.21257d-01,-2.20162d-01,-2.19066d-01,
     1 -2.17971d-01,-2.16875d-01,-2.15780d-01,-2.14684d-01,-2.13589d-01,
     2 -2.12493d-01,-2.11397d-01,-2.10302d-01,-2.09206d-01,-2.08111d-01,
     3 -2.07015d-01,-2.05920d-01,-2.04824d-01,-2.03728d-01,-2.02633d-01,
     4 -2.01537d-01,-2.00441d-01,-1.99346d-01,-1.98250d-01,-1.97154d-01,
     5 -1.96059d-01,-1.94963d-01,-1.93867d-01,-1.92771d-01,-1.91676d-01,
     6 -1.90580d-01,-1.89484d-01,-1.88388d-01,-1.87293d-01,-1.86197d-01,
     7 -1.85101d-01,-1.84005d-01,-1.82909d-01,-1.81813d-01,-1.80717d-01,
     8 -1.79621d-01,-1.78525d-01,-1.77429d-01,-1.76333d-01,-1.75237d-01,
     9 -1.74141d-01,-1.73045d-01,-1.71949d-01,-1.70853d-01,-1.69756d-01/
      data (eqlm(i),i=1297,1344)/-1.68660d-01,-1.67564d-01,-1.66468d-01,
     1 -1.65371d-01,-1.64275d-01,-1.63179d-01,-1.62082d-01,-1.60986d-01,
     2 -1.59890d-01,-1.58793d-01,-1.57697d-01,-1.56600d-01,-1.55503d-01,
     3 -1.54407d-01,-1.53310d-01,-1.52214d-01,-1.51117d-01,-1.50020d-01,
     4 -1.48923d-01,-1.47827d-01,-1.46730d-01,-1.45633d-01,-1.44536d-01,
     5 -1.43439d-01,-1.42342d-01,-1.41245d-01,-1.40148d-01,-1.39051d-01,
     6 -1.37954d-01,-1.36856d-01,-1.35759d-01,-1.34662d-01,-1.33565d-01,
     7 -1.32467d-01,-1.31370d-01,-1.30272d-01,-1.29175d-01,-1.28077d-01,
     8 -1.26980d-01,-1.25882d-01,-1.24784d-01,-1.23687d-01,-1.22589d-01,
     9 -1.21491d-01,-1.20393d-01,-1.19295d-01,-1.18197d-01,-1.17099d-01/
      data (eqlm(i),i=1345,1392)/-1.16001d-01,-1.14903d-01,-1.13804d-01,
     1 -1.12706d-01,-1.11608d-01,-1.10509d-01,-1.09411d-01,-1.08312d-01,
     2 -1.07214d-01,-1.06115d-01,-1.05016d-01,-1.03918d-01,-1.02819d-01,
     3 -1.01720d-01,-1.00621d-01,-9.95221d-02,-9.84230d-02,-9.73239d-02,
     4 -9.62247d-02,-9.51254d-02,-9.40261d-02,-9.29267d-02,-9.18273d-02,
     5 -9.07278d-02,-8.96282d-02,-8.85286d-02,-8.74288d-02,-8.63291d-02,
     6 -8.52292d-02,-8.41293d-02,-8.30293d-02,-8.19293d-02,-8.08291d-02,
     7 -7.97290d-02,-7.86287d-02,-7.75284d-02,-7.64280d-02,-7.53275d-02,
     8 -7.42270d-02,-7.31264d-02,-7.20257d-02,-7.09249d-02,-6.98241d-02,
     9 -6.87232d-02,-6.76222d-02,-6.65212d-02,-6.54201d-02,-6.43189d-02/
      data (eqlm(i),i=1393,1440)/-6.32176d-02,-6.21163d-02,-6.10149d-02,
     1 -5.99134d-02,-5.88118d-02,-5.77101d-02,-5.66084d-02,-5.55066d-02,
     2 -5.44047d-02,-5.33027d-02,-5.22007d-02,-5.10986d-02,-4.99963d-02,
     3 -4.88940d-02,-4.77917d-02,-4.66892d-02,-4.55867d-02,-4.44840d-02,
     4 -4.33813d-02,-4.22785d-02,-4.11756d-02,-4.00727d-02,-3.89696d-02,
     5 -3.78665d-02,-3.67632d-02,-3.56599d-02,-3.45565d-02,-3.34530d-02,
     6 -3.23494d-02,-3.12457d-02,-3.01419d-02,-2.90380d-02,-2.79341d-02,
     7 -2.68300d-02,-2.57258d-02,-2.46216d-02,-2.35172d-02,-2.24128d-02,
     8 -2.13082d-02,-2.02036d-02,-1.90989d-02,-1.79940d-02,-1.68891d-02,
     9 -1.57840d-02,-1.46789d-02,-1.35737d-02,-1.24683d-02,-1.13629d-02/
      data (eqlm(i),i=1441,1488)/-1.02573d-02,-9.15170d-03,-8.04596d-03,
     1 -6.94011d-03,-5.83415d-03,-4.72809d-03,-3.62193d-03,-2.51566d-03,
     2 -1.40929d-03,-3.02806d-04, 8.03781d-04, 1.91048d-03, 3.01728d-03,
     3  4.12419d-03, 5.23121d-03, 6.33834d-03, 7.44557d-03, 8.55292d-03,
     4  9.66038d-03, 1.07680d-02, 1.18756d-02, 1.29834d-02, 1.40913d-02,
     5  1.51994d-02, 1.63075d-02, 1.74157d-02, 1.85241d-02, 1.96326d-02,
     6  2.07412d-02, 2.18499d-02, 2.29587d-02, 2.40676d-02, 2.51767d-02,
     7  2.62859d-02, 2.73952d-02, 2.85046d-02, 2.96141d-02, 3.07238d-02,
     8  3.18335d-02, 3.29434d-02, 3.40534d-02, 3.51636d-02, 3.62738d-02,
     9  3.73842d-02, 3.84947d-02, 3.96053d-02, 4.07161d-02, 4.18269d-02/
      data (eqlm(i),i=1489,1536)/ 4.29379d-02, 4.40490d-02, 4.51603d-02,
     1  4.62716d-02, 4.73831d-02, 4.84947d-02, 4.96065d-02, 5.07183d-02,
     2  5.18303d-02, 5.29424d-02, 5.40547d-02, 5.51670d-02, 5.62795d-02,
     3  5.73922d-02, 5.85049d-02, 5.96178d-02, 6.07308d-02, 6.18440d-02,
     4  6.29572d-02, 6.40706d-02, 6.51842d-02, 6.62978d-02, 6.74117d-02,
     5  6.85256d-02, 6.96397d-02, 7.07539d-02, 7.18682d-02, 7.29827d-02,
     6  7.40973d-02, 7.52120d-02, 7.63269d-02, 7.74419d-02, 7.85571d-02,
     7  7.96724d-02, 8.07878d-02, 8.19034d-02, 8.30191d-02, 8.41349d-02,
     8  8.52509d-02, 8.63671d-02, 8.74833d-02, 8.85998d-02, 8.97163d-02,
     9  9.08330d-02, 9.19499d-02, 9.30669d-02, 9.41840d-02, 9.53013d-02/
      data (eqlm(i),i=1537,1584)/ 9.64188d-02, 9.75364d-02, 9.86541d-02,
     1  9.97720d-02, 1.00890d-01, 1.02008d-01, 1.03127d-01, 1.04245d-01,
     2  1.05364d-01, 1.06483d-01, 1.07602d-01, 1.08721d-01, 1.09840d-01,
     3  1.10959d-01, 1.12079d-01, 1.13199d-01, 1.14319d-01, 1.15439d-01,
     4  1.16559d-01, 1.17679d-01, 1.18800d-01, 1.19920d-01, 1.21041d-01,
     5  1.22162d-01, 1.23283d-01, 1.24405d-01, 1.25526d-01, 1.26648d-01,
     6  1.27770d-01, 1.28892d-01, 1.30014d-01, 1.31136d-01, 1.32258d-01,
     7  1.33381d-01, 1.34504d-01, 1.35627d-01, 1.36750d-01, 1.37873d-01,
     8  1.38996d-01, 1.40120d-01, 1.41244d-01, 1.42368d-01, 1.43492d-01,
     9  1.44616d-01, 1.45740d-01, 1.46865d-01, 1.47990d-01, 1.49114d-01/
      data (eqlm(i),i=1585,1632)/ 1.50240d-01, 1.51365d-01, 1.52490d-01,
     1  1.53616d-01, 1.54742d-01, 1.55868d-01, 1.56994d-01, 1.58120d-01,
     2  1.59246d-01, 1.60373d-01, 1.61500d-01, 1.62627d-01, 1.63754d-01,
     3  1.64881d-01, 1.66009d-01, 1.67136d-01, 1.68264d-01, 1.69392d-01,
     4  1.70521d-01, 1.71649d-01, 1.72778d-01, 1.73906d-01, 1.75035d-01,
     5  1.76164d-01, 1.77294d-01, 1.78423d-01, 1.79553d-01, 1.80683d-01,
     6  1.81813d-01, 1.82943d-01, 1.84073d-01, 1.85204d-01, 1.86335d-01,
     7  1.87466d-01, 1.88597d-01, 1.89728d-01, 1.90860d-01, 1.91991d-01,
     8  1.93123d-01, 1.94255d-01, 1.95387d-01, 1.96520d-01, 1.97652d-01,
     9  1.98785d-01, 1.99918d-01, 2.01051d-01, 2.02185d-01, 2.03318d-01/
      data (eqlm(i),i=1633,1680)/ 2.04452d-01, 2.05586d-01, 2.06720d-01,
     1  2.07855d-01, 2.08989d-01, 2.10124d-01, 2.11259d-01, 2.12394d-01,
     2  2.13529d-01, 2.14665d-01, 2.15800d-01, 2.16936d-01, 2.18072d-01,
     3  2.19209d-01, 2.20345d-01, 2.21482d-01, 2.22619d-01, 2.23756d-01,
     4  2.24893d-01, 2.26030d-01, 2.27168d-01, 2.28306d-01, 2.29444d-01,
     5  2.30582d-01, 2.31721d-01, 2.32860d-01, 2.33999d-01, 2.35138d-01,
     6  2.36277d-01, 2.37416d-01, 2.38556d-01, 2.39696d-01, 2.40836d-01,
     7  2.41977d-01, 2.43117d-01, 2.44258d-01, 2.45399d-01, 2.46540d-01,
     8  2.47682d-01, 2.48823d-01, 2.49965d-01, 2.51107d-01, 2.52249d-01,
     9  2.53392d-01, 2.54534d-01, 2.55677d-01, 2.56820d-01, 2.57964d-01/
      data (eqlm(i),i=1681,1728)/ 2.59107d-01, 2.60251d-01, 2.61395d-01,
     1  2.62539d-01, 2.63683d-01, 2.64828d-01, 2.65973d-01, 2.67118d-01,
     2  2.68263d-01, 2.69408d-01, 2.70554d-01, 2.71700d-01, 2.72846d-01,
     3  2.73992d-01, 2.75139d-01, 2.76286d-01, 2.77433d-01, 2.78580d-01,
     4  2.79727d-01, 2.80875d-01, 2.82023d-01, 2.83171d-01, 2.84319d-01,
     5  2.85468d-01, 2.86616d-01, 2.87765d-01, 2.88915d-01, 2.90064d-01,
     6  2.91214d-01, 2.92364d-01, 2.93514d-01, 2.94664d-01, 2.95815d-01,
     7  2.96965d-01, 2.98116d-01, 2.99268d-01, 3.00419d-01, 3.01571d-01,
     8  3.02723d-01, 3.03875d-01, 3.05027d-01, 3.06180d-01, 3.07333d-01,
     9  3.08486d-01, 3.09639d-01, 3.10793d-01, 3.11946d-01, 3.13100d-01/
      data (eqlm(i),i=1729,1776)/ 3.14255d-01, 3.15409d-01, 3.16564d-01,
     1  3.17719d-01, 3.18874d-01, 3.20030d-01, 3.21185d-01, 3.22341d-01,
     2  3.23497d-01, 3.24654d-01, 3.25810d-01, 3.26967d-01, 3.28125d-01,
     3  3.29282d-01, 3.30440d-01, 3.31597d-01, 3.32756d-01, 3.33914d-01,
     4  3.35073d-01, 3.36231d-01, 3.37391d-01, 3.38550d-01, 3.39710d-01,
     5  3.40869d-01, 3.42030d-01, 3.43190d-01, 3.44351d-01, 3.45511d-01,
     6  3.46673d-01, 3.47834d-01, 3.48996d-01, 3.50157d-01, 3.51320d-01,
     7  3.52482d-01, 3.53645d-01, 3.54807d-01, 3.55971d-01, 3.57134d-01,
     8  3.58298d-01, 3.59462d-01, 3.60626d-01, 3.61790d-01, 3.62955d-01,
     9  3.64120d-01, 3.65285d-01, 3.66450d-01, 3.67616d-01, 3.68782d-01/
      data (eqlm(i),i=1777,1824)/ 3.69948d-01, 3.71115d-01, 3.72281d-01,
     1  3.73448d-01, 3.74616d-01, 3.75783d-01, 3.76951d-01, 3.78119d-01,
     2  3.79287d-01, 3.80456d-01, 3.81625d-01, 3.82794d-01, 3.83963d-01,
     3  3.85133d-01, 3.86303d-01, 3.87473d-01, 3.88643d-01, 3.89814d-01,
     4  3.90985d-01, 3.92156d-01, 3.93327d-01, 3.94499d-01, 3.95671d-01,
     5  3.96843d-01, 3.98016d-01, 3.99189d-01, 4.00362d-01, 4.01535d-01,
     6  4.02709d-01, 4.03883d-01, 4.05057d-01, 4.06231d-01, 4.07406d-01,
     7  4.08581d-01, 4.09756d-01, 4.10932d-01, 4.12108d-01, 4.13284d-01,
     8  4.14460d-01, 4.15637d-01, 4.16814d-01, 4.17991d-01, 4.19168d-01,
     9  4.20346d-01, 4.21524d-01, 4.22702d-01, 4.23881d-01, 4.25060d-01/
      data (eqlm(i),i=1825,1872)/ 4.26239d-01, 4.27418d-01, 4.28598d-01,
     1  4.29778d-01, 4.30958d-01, 4.32139d-01, 4.33320d-01, 4.34501d-01,
     2  4.35682d-01, 4.36864d-01, 4.38046d-01, 4.39228d-01, 4.40411d-01,
     3  4.41594d-01, 4.42777d-01, 4.43960d-01, 4.45144d-01, 4.46328d-01,
     4  4.47512d-01, 4.48697d-01, 4.49882d-01, 4.51067d-01, 4.52253d-01,
     5  4.53438d-01, 4.54625d-01, 4.55811d-01, 4.56998d-01, 4.58185d-01,
     6  4.59372d-01, 4.60559d-01, 4.61747d-01, 4.62935d-01, 4.64124d-01,
     7  4.65313d-01, 4.66502d-01, 4.67691d-01, 4.68881d-01, 4.70071d-01,
     8  4.71261d-01, 4.72451d-01, 4.73642d-01, 4.74833d-01, 4.76025d-01,
     9  4.77217d-01, 4.78409d-01, 4.79601d-01, 4.80794d-01, 4.81987d-01/
      data (eqlm(i),i=1873,1920)/ 4.83180d-01, 4.84373d-01, 4.85567d-01,
     1  4.86762d-01, 4.87956d-01, 4.89151d-01, 4.90346d-01, 4.91541d-01,
     2  4.92737d-01, 4.93933d-01, 4.95130d-01, 4.96326d-01, 4.97523d-01,
     3  4.98720d-01, 4.99918d-01, 5.01116d-01, 5.02314d-01, 5.03513d-01,
     4  5.04712d-01, 5.05911d-01, 5.07110d-01, 5.08310d-01, 5.09510d-01,
     5  5.10711d-01, 5.11911d-01, 5.13112d-01, 5.14314d-01, 5.15515d-01,
     6  5.16718d-01, 5.17920d-01, 5.19123d-01, 5.20325d-01, 5.21529d-01,
     7  5.22732d-01, 5.23936d-01, 5.25141d-01, 5.26345d-01, 5.27550d-01,
     8  5.28755d-01, 5.29961d-01, 5.31167d-01, 5.32373d-01, 5.33580d-01,
     9  5.34786d-01, 5.35994d-01, 5.37201d-01, 5.38409d-01, 5.39617d-01/
      data (eqlm(i),i=1921,1968)/ 5.40826d-01, 5.42035d-01, 5.43244d-01,
     1  5.44453d-01, 5.45663d-01, 5.46873d-01, 5.48084d-01, 5.49295d-01,
     2  5.50506d-01, 5.51717d-01, 5.52929d-01, 5.54141d-01, 5.55354d-01,
     3  5.56567d-01, 5.57780d-01, 5.58993d-01, 5.60207d-01, 5.61421d-01,
     4  5.62636d-01, 5.63851d-01, 5.65066d-01, 5.66281d-01, 5.67497d-01,
     5  5.68714d-01, 5.69930d-01, 5.71147d-01, 5.72364d-01, 5.73582d-01,
     6  5.74800d-01, 5.76018d-01, 5.77237d-01, 5.78456d-01, 5.79675d-01,
     7  5.80895d-01, 5.82115d-01, 5.83335d-01, 5.84556d-01, 5.85777d-01,
     8  5.86999d-01, 5.88220d-01, 5.89442d-01, 5.90665d-01, 5.91888d-01,
     9  5.93111d-01, 5.94334d-01, 5.95558d-01, 5.96783d-01, 5.98007d-01/
      data (eqlm(i),i=1969,2016)/ 5.99232d-01, 6.00458d-01, 6.01683d-01,
     1  6.02909d-01, 6.04136d-01, 6.05362d-01, 6.06589d-01, 6.07817d-01,
     2  6.09045d-01, 6.10273d-01, 6.11501d-01, 6.12730d-01, 6.13959d-01,
     3  6.15189d-01, 6.16419d-01, 6.17649d-01, 6.18880d-01, 6.20111d-01,
     4  6.21343d-01, 6.22574d-01, 6.23806d-01, 6.25039d-01, 6.26272d-01,
     5  6.27505d-01, 6.28739d-01, 6.29973d-01, 6.31207d-01, 6.32442d-01,
     6  6.33677d-01, 6.34912d-01, 6.36148d-01, 6.37384d-01, 6.38621d-01,
     7  6.39858d-01, 6.41095d-01, 6.42333d-01, 6.43571d-01, 6.44809d-01,
     8  6.46048d-01, 6.47287d-01, 6.48527d-01, 6.49767d-01, 6.51007d-01,
     9  6.52248d-01, 6.53489d-01, 6.54730d-01, 6.55972d-01, 6.57214d-01/
      data (eqlm(i),i=2017,2064)/ 6.58457d-01, 6.59699d-01, 6.60943d-01,
     1  6.62186d-01, 6.63431d-01, 6.64675d-01, 6.65920d-01, 6.67165d-01,
     2  6.68411d-01, 6.69657d-01, 6.70903d-01, 6.72150d-01, 6.73397d-01,
     3  6.74644d-01, 6.75892d-01, 6.77141d-01, 6.78389d-01, 6.79638d-01,
     4  6.80888d-01, 6.82138d-01, 6.83388d-01, 6.84639d-01, 6.85890d-01,
     5  6.87141d-01, 6.88393d-01, 6.89645d-01, 6.90898d-01, 6.92151d-01,
     6  6.93404d-01, 6.94658d-01, 6.95912d-01, 6.97167d-01, 6.98422d-01,
     7  6.99677d-01, 7.00933d-01, 7.02189d-01, 7.03445d-01, 7.04702d-01,
     8  7.05960d-01, 7.07218d-01, 7.08476d-01, 7.09735d-01, 7.10994d-01,
     9  7.12253d-01, 7.13513d-01, 7.14773d-01, 7.16034d-01, 7.17295d-01/
      data (eqlm(i),i=2065,2112)/ 7.18556d-01, 7.19818d-01, 7.21080d-01,
     1  7.22343d-01, 7.23606d-01, 7.24870d-01, 7.26134d-01, 7.27398d-01,
     2  7.28663d-01, 7.29928d-01, 7.31193d-01, 7.32459d-01, 7.33726d-01,
     3  7.34992d-01, 7.36260d-01, 7.37527d-01, 7.38795d-01, 7.40064d-01,
     4  7.41333d-01, 7.42602d-01, 7.43872d-01, 7.45142d-01, 7.46412d-01,
     5  7.47683d-01, 7.48954d-01, 7.50226d-01, 7.51498d-01, 7.52771d-01,
     6  7.54044d-01, 7.55318d-01, 7.56591d-01, 7.57866d-01, 7.59140d-01,
     7  7.60416d-01, 7.61691d-01, 7.62967d-01, 7.64244d-01, 7.65520d-01,
     8  7.66798d-01, 7.68075d-01, 7.69354d-01, 7.70632d-01, 7.71911d-01,
     9  7.73191d-01, 7.74471d-01, 7.75751d-01, 7.77032d-01, 7.78313d-01/
      data (eqlm(i),i=2113,2160)/ 7.79594d-01, 7.80876d-01, 7.82159d-01,
     1  7.83442d-01, 7.84725d-01, 7.86009d-01, 7.87293d-01, 7.88578d-01,
     2  7.89863d-01, 7.91149d-01, 7.92435d-01, 7.93721d-01, 7.95008d-01,
     3  7.96295d-01, 7.97583d-01, 7.98871d-01, 8.00160d-01, 8.01449d-01,
     4  8.02738d-01, 8.04028d-01, 8.05319d-01, 8.06609d-01, 8.07901d-01,
     5  8.09193d-01, 8.10485d-01, 8.11777d-01, 8.13070d-01, 8.14364d-01,
     6  8.15658d-01, 8.16952d-01, 8.18247d-01, 8.19543d-01, 8.20838d-01,
     7  8.22135d-01, 8.23431d-01, 8.24729d-01, 8.26026d-01, 8.27324d-01,
     8  8.28623d-01, 8.29922d-01, 8.31221d-01, 8.32521d-01, 8.33822d-01,
     9  8.35122d-01, 8.36424d-01, 8.37725d-01, 8.39028d-01, 8.40330d-01/
      data (eqlm(i),i=2161,2208)/ 8.41634d-01, 8.42937d-01, 8.44241d-01,
     1  8.45546d-01, 8.46851d-01, 8.48156d-01, 8.49462d-01, 8.50769d-01,
     2  8.52076d-01, 8.53383d-01, 8.54691d-01, 8.55999d-01, 8.57308d-01,
     3  8.58617d-01, 8.59927d-01, 8.61237d-01, 8.62548d-01, 8.63859d-01,
     4  8.65171d-01, 8.66483d-01, 8.67795d-01, 8.69108d-01, 8.70422d-01,
     5  8.71736d-01, 8.73050d-01, 8.74365d-01, 8.75681d-01, 8.76997d-01,
     6  8.78313d-01, 8.79630d-01, 8.80947d-01, 8.82265d-01, 8.83583d-01,
     7  8.84902d-01, 8.86222d-01, 8.87541d-01, 8.88862d-01, 8.90182d-01,
     8  8.91503d-01, 8.92825d-01, 8.94147d-01, 8.95470d-01, 8.96793d-01,
     9  8.98117d-01, 8.99441d-01, 9.00766d-01, 9.02091d-01, 9.03416d-01/
      data (eqlm(i),i=2209,2256)/ 9.04742d-01, 9.06069d-01, 9.07396d-01,
     1  9.08724d-01, 9.10052d-01, 9.11380d-01, 9.12709d-01, 9.14039d-01,
     2  9.15369d-01, 9.16700d-01, 9.18031d-01, 9.19362d-01, 9.20694d-01,
     3  9.22027d-01, 9.23360d-01, 9.24694d-01, 9.26028d-01, 9.27362d-01,
     4  9.28697d-01, 9.30033d-01, 9.31369d-01, 9.32706d-01, 9.34043d-01,
     5  9.35381d-01, 9.36719d-01, 9.38057d-01, 9.39396d-01, 9.40736d-01,
     6  9.42076d-01, 9.43417d-01, 9.44758d-01, 9.46100d-01, 9.47442d-01,
     7  9.48785d-01, 9.50128d-01, 9.51472d-01, 9.52816d-01, 9.54161d-01,
     8  9.55506d-01, 9.56852d-01, 9.58199d-01, 9.59545d-01, 9.60893d-01,
     9  9.62241d-01, 9.63589d-01, 9.64938d-01, 9.66287d-01, 9.67637d-01/
      data (eqlm(i),i=2257,2304)/ 9.68988d-01, 9.70339d-01, 9.71691d-01,
     1  9.73043d-01, 9.74395d-01, 9.75748d-01, 9.77102d-01, 9.78456d-01,
     2  9.79811d-01, 9.81166d-01, 9.82522d-01, 9.83878d-01, 9.85235d-01,
     3  9.86593d-01, 9.87951d-01, 9.89309d-01, 9.90668d-01, 9.92028d-01,
     4  9.93388d-01, 9.94748d-01, 9.96109d-01, 9.97471d-01, 9.98833d-01,
     5  1.00020d+00, 1.00156d+00, 1.00292d+00, 1.00429d+00, 1.00565d+00,
     6  1.00702d+00, 1.00838d+00, 1.00975d+00, 1.01112d+00, 1.01249d+00,
     7  1.01385d+00, 1.01522d+00, 1.01659d+00, 1.01796d+00, 1.01933d+00,
     8  1.02070d+00, 1.02207d+00, 1.02345d+00, 1.02482d+00, 1.02619d+00,
     9  1.02757d+00, 1.02894d+00, 1.03032d+00, 1.03169d+00, 1.03307d+00/
      data (eqlm(i),i=2305,2352)/ 1.03444d+00, 1.03582d+00, 1.03720d+00,
     1  1.03858d+00, 1.03996d+00, 1.04133d+00, 1.04271d+00, 1.04410d+00,
     2  1.04548d+00, 1.04686d+00, 1.04824d+00, 1.04962d+00, 1.05101d+00,
     3  1.05239d+00, 1.05377d+00, 1.05516d+00, 1.05655d+00, 1.05793d+00,
     4  1.05932d+00, 1.06071d+00, 1.06209d+00, 1.06348d+00, 1.06487d+00,
     5  1.06626d+00, 1.06765d+00, 1.06904d+00, 1.07043d+00, 1.07182d+00,
     6  1.07322d+00, 1.07461d+00, 1.07600d+00, 1.07740d+00, 1.07879d+00,
     7  1.08019d+00, 1.08158d+00, 1.08298d+00, 1.08438d+00, 1.08578d+00,
     8  1.08717d+00, 1.08857d+00, 1.08997d+00, 1.09137d+00, 1.09277d+00,
     9  1.09417d+00, 1.09558d+00, 1.09698d+00, 1.09838d+00, 1.09979d+00/
      data (eqlm(i),i=2353,2400)/ 1.10119d+00, 1.10259d+00, 1.10400d+00,
     1  1.10541d+00, 1.10681d+00, 1.10822d+00, 1.10963d+00, 1.11104d+00,
     2  1.11244d+00, 1.11385d+00, 1.11526d+00, 1.11668d+00, 1.11809d+00,
     3  1.11950d+00, 1.12091d+00, 1.12232d+00, 1.12374d+00, 1.12515d+00,
     4  1.12657d+00, 1.12798d+00, 1.12940d+00, 1.13082d+00, 1.13223d+00,
     5  1.13365d+00, 1.13507d+00, 1.13649d+00, 1.13791d+00, 1.13933d+00,
     6  1.14075d+00, 1.14217d+00, 1.14360d+00, 1.14502d+00, 1.14644d+00,
     7  1.14787d+00, 1.14929d+00, 1.15072d+00, 1.15214d+00, 1.15357d+00,
     8  1.15500d+00, 1.15642d+00, 1.15785d+00, 1.15928d+00, 1.16071d+00,
     9  1.16214d+00, 1.16357d+00, 1.16500d+00, 1.16644d+00, 1.16787d+00/
      data (eqlm(i),i=2401,2448)/ 1.16930d+00, 1.17074d+00, 1.17217d+00,
     1  1.17361d+00, 1.17504d+00, 1.17648d+00, 1.17792d+00, 1.17936d+00,
     2  1.18079d+00, 1.18223d+00, 1.18367d+00, 1.18511d+00, 1.18656d+00,
     3  1.18800d+00, 1.18944d+00, 1.19088d+00, 1.19233d+00, 1.19377d+00,
     4  1.19522d+00, 1.19666d+00, 1.19811d+00, 1.19956d+00, 1.20100d+00,
     5  1.20245d+00, 1.20390d+00, 1.20535d+00, 1.20680d+00, 1.20825d+00,
     6  1.20970d+00, 1.21115d+00, 1.21261d+00, 1.21406d+00, 1.21551d+00,
     7  1.21697d+00, 1.21843d+00, 1.21988d+00, 1.22134d+00, 1.22280d+00,
     8  1.22425d+00, 1.22571d+00, 1.22717d+00, 1.22863d+00, 1.23009d+00,
     9  1.23155d+00, 1.23302d+00, 1.23448d+00, 1.23594d+00, 1.23741d+00/
      data (eqlm(i),i=2449,2496)/ 1.23887d+00, 1.24034d+00, 1.24180d+00,
     1  1.24327d+00, 1.24474d+00, 1.24620d+00, 1.24767d+00, 1.24914d+00,
     2  1.25061d+00, 1.25208d+00, 1.25356d+00, 1.25503d+00, 1.25650d+00,
     3  1.25797d+00, 1.25945d+00, 1.26092d+00, 1.26240d+00, 1.26387d+00,
     4  1.26535d+00, 1.26683d+00, 1.26831d+00, 1.26979d+00, 1.27127d+00,
     5  1.27275d+00, 1.27423d+00, 1.27571d+00, 1.27719d+00, 1.27868d+00,
     6  1.28016d+00, 1.28164d+00, 1.28313d+00, 1.28461d+00, 1.28610d+00,
     7  1.28759d+00, 1.28908d+00, 1.29057d+00, 1.29206d+00, 1.29355d+00,
     8  1.29504d+00, 1.29653d+00, 1.29802d+00, 1.29951d+00, 1.30101d+00,
     9  1.30250d+00, 1.30400d+00, 1.30549d+00, 1.30699d+00, 1.30849d+00/
      data (eqlm(i),i=2497,2544)/ 1.30998d+00, 1.31148d+00, 1.31298d+00,
     1  1.31448d+00, 1.31598d+00, 1.31748d+00, 1.31899d+00, 1.32049d+00,
     2  1.32199d+00, 1.32350d+00, 1.32500d+00, 1.32651d+00, 1.32801d+00,
     3  1.32952d+00, 1.33103d+00, 1.33254d+00, 1.33405d+00, 1.33556d+00,
     4  1.33707d+00, 1.33858d+00, 1.34009d+00, 1.34161d+00, 1.34312d+00,
     5  1.34463d+00, 1.34615d+00, 1.34767d+00, 1.34918d+00, 1.35070d+00,
     6  1.35222d+00, 1.35374d+00, 1.35526d+00, 1.35678d+00, 1.35830d+00,
     7  1.35982d+00, 1.36134d+00, 1.36287d+00, 1.36439d+00, 1.36591d+00,
     8  1.36744d+00, 1.36897d+00, 1.37049d+00, 1.37202d+00, 1.37355d+00,
     9  1.37508d+00, 1.37661d+00, 1.37814d+00, 1.37967d+00, 1.38121d+00/
      data (eqlm(i),i=2545,2592)/ 1.38274d+00, 1.38427d+00, 1.38581d+00,
     1  1.38734d+00, 1.38888d+00, 1.39042d+00, 1.39195d+00, 1.39349d+00,
     2  1.39503d+00, 1.39657d+00, 1.39811d+00, 1.39965d+00, 1.40120d+00,
     3  1.40274d+00, 1.40428d+00, 1.40583d+00, 1.40737d+00, 1.40892d+00,
     4  1.41047d+00, 1.41201d+00, 1.41356d+00, 1.41511d+00, 1.41666d+00,
     5  1.41821d+00, 1.41976d+00, 1.42132d+00, 1.42287d+00, 1.42442d+00,
     6  1.42598d+00, 1.42753d+00, 1.42909d+00, 1.43065d+00, 1.43221d+00,
     7  1.43376d+00, 1.43532d+00, 1.43688d+00, 1.43844d+00, 1.44001d+00,
     8  1.44157d+00, 1.44313d+00, 1.44470d+00, 1.44626d+00, 1.44783d+00,
     9  1.44939d+00, 1.45096d+00, 1.45253d+00, 1.45410d+00, 1.45567d+00/
      data (eqlm(i),i=2593,2640)/ 1.45724d+00, 1.45881d+00, 1.46038d+00,
     1  1.46196d+00, 1.46353d+00, 1.46510d+00, 1.46668d+00, 1.46826d+00,
     2  1.46983d+00, 1.47141d+00, 1.47299d+00, 1.47457d+00, 1.47615d+00,
     3  1.47773d+00, 1.47931d+00, 1.48090d+00, 1.48248d+00, 1.48406d+00,
     4  1.48565d+00, 1.48724d+00, 1.48882d+00, 1.49041d+00, 1.49200d+00,
     5  1.49359d+00, 1.49518d+00, 1.49677d+00, 1.49836d+00, 1.49995d+00,
     6  1.50155d+00, 1.50314d+00, 1.50474d+00, 1.50633d+00, 1.50793d+00,
     7  1.50953d+00, 1.51113d+00, 1.51273d+00, 1.51433d+00, 1.51593d+00,
     8  1.51753d+00, 1.51913d+00, 1.52074d+00, 1.52234d+00, 1.52395d+00,
     9  1.52555d+00, 1.52716d+00, 1.52877d+00, 1.53038d+00, 1.53199d+00/
      data (eqlm(i),i=2641,2688)/ 1.53360d+00, 1.53521d+00, 1.53682d+00,
     1  1.53843d+00, 1.54005d+00, 1.54166d+00, 1.54328d+00, 1.54489d+00,
     2  1.54651d+00, 1.54813d+00, 1.54975d+00, 1.55137d+00, 1.55299d+00,
     3  1.55461d+00, 1.55623d+00, 1.55786d+00, 1.55948d+00, 1.56111d+00,
     4  1.56273d+00, 1.56436d+00, 1.56599d+00, 1.56762d+00, 1.56925d+00,
     5  1.57088d+00, 1.57251d+00, 1.57414d+00, 1.57577d+00, 1.57741d+00,
     6  1.57904d+00, 1.58068d+00, 1.58232d+00, 1.58395d+00, 1.58559d+00,
     7  1.58723d+00, 1.58887d+00, 1.59051d+00, 1.59215d+00, 1.59380d+00,
     8  1.59544d+00, 1.59709d+00, 1.59873d+00, 1.60038d+00, 1.60203d+00,
     9  1.60367d+00, 1.60532d+00, 1.60697d+00, 1.60863d+00, 1.61028d+00/
      data (eqlm(i),i=2689,2736)/ 1.61193d+00, 1.61358d+00, 1.61524d+00,
     1  1.61689d+00, 1.61855d+00, 1.62021d+00, 1.62187d+00, 1.62353d+00,
     2  1.62519d+00, 1.62685d+00, 1.62851d+00, 1.63017d+00, 1.63184d+00,
     3  1.63350d+00, 1.63517d+00, 1.63683d+00, 1.63850d+00, 1.64017d+00,
     4  1.64184d+00, 1.64351d+00, 1.64518d+00, 1.64685d+00, 1.64853d+00,
     5  1.65020d+00, 1.65188d+00, 1.65355d+00, 1.65523d+00, 1.65691d+00,
     6  1.65859d+00, 1.66027d+00, 1.66195d+00, 1.66363d+00, 1.66531d+00,
     7  1.66700d+00, 1.66868d+00, 1.67037d+00, 1.67205d+00, 1.67374d+00,
     8  1.67543d+00, 1.67712d+00, 1.67881d+00, 1.68050d+00, 1.68219d+00,
     9  1.68389d+00, 1.68558d+00, 1.68728d+00, 1.68897d+00, 1.69067d+00/
      data (eqlm(i),i=2737,2784)/ 1.69237d+00, 1.69407d+00, 1.69577d+00,
     1  1.69747d+00, 1.69917d+00, 1.70087d+00, 1.70258d+00, 1.70428d+00,
     2  1.70599d+00, 1.70769d+00, 1.70940d+00, 1.71111d+00, 1.71282d+00,
     3  1.71453d+00, 1.71624d+00, 1.71796d+00, 1.71967d+00, 1.72139d+00,
     4  1.72310d+00, 1.72482d+00, 1.72654d+00, 1.72825d+00, 1.72997d+00,
     5  1.73170d+00, 1.73342d+00, 1.73514d+00, 1.73686d+00, 1.73859d+00,
     6  1.74031d+00, 1.74204d+00, 1.74377d+00, 1.74550d+00, 1.74723d+00,
     7  1.74896d+00, 1.75069d+00, 1.75242d+00, 1.75416d+00, 1.75589d+00,
     8  1.75763d+00, 1.75936d+00, 1.76110d+00, 1.76284d+00, 1.76458d+00,
     9  1.76632d+00, 1.76806d+00, 1.76981d+00, 1.77155d+00, 1.77330d+00/
      data (eqlm(i),i=2785,2832)/ 1.77504d+00, 1.77679d+00, 1.77854d+00,
     1  1.78029d+00, 1.78204d+00, 1.78379d+00, 1.78554d+00, 1.78730d+00,
     2  1.78905d+00, 1.79081d+00, 1.79256d+00, 1.79432d+00, 1.79608d+00,
     3  1.79784d+00, 1.79960d+00, 1.80136d+00, 1.80313d+00, 1.80489d+00,
     4  1.80665d+00, 1.80842d+00, 1.81019d+00, 1.81196d+00, 1.81373d+00,
     5  1.81550d+00, 1.81727d+00, 1.81904d+00, 1.82081d+00, 1.82259d+00,
     6  1.82436d+00, 1.82614d+00, 1.82792d+00, 1.82970d+00, 1.83148d+00,
     7  1.83326d+00, 1.83504d+00, 1.83682d+00, 1.83861d+00, 1.84039d+00,
     8  1.84218d+00, 1.84397d+00, 1.84576d+00, 1.84755d+00, 1.84934d+00,
     9  1.85113d+00, 1.85292d+00, 1.85472d+00, 1.85651d+00, 1.85831d+00/
      data (eqlm(i),i=2833,2880)/ 1.86011d+00, 1.86191d+00, 1.86371d+00,
     1  1.86551d+00, 1.86731d+00, 1.86911d+00, 1.87092d+00, 1.87272d+00,
     2  1.87453d+00, 1.87634d+00, 1.87814d+00, 1.87995d+00, 1.88176d+00,
     3  1.88358d+00, 1.88539d+00, 1.88720d+00, 1.88902d+00, 1.89084d+00,
     4  1.89265d+00, 1.89447d+00, 1.89629d+00, 1.89811d+00, 1.89994d+00,
     5  1.90176d+00, 1.90358d+00, 1.90541d+00, 1.90724d+00, 1.90906d+00,
     6  1.91089d+00, 1.91272d+00, 1.91455d+00, 1.91639d+00, 1.91822d+00,
     7  1.92005d+00, 1.92189d+00, 1.92373d+00, 1.92557d+00, 1.92741d+00,
     8  1.92925d+00, 1.93109d+00, 1.93293d+00, 1.93477d+00, 1.93662d+00,
     9  1.93847d+00, 1.94031d+00, 1.94216d+00, 1.94401d+00, 1.94586d+00/
      data (eqlm(i),i=2881,2928)/ 1.94772d+00, 1.94957d+00, 1.95142d+00,
     1  1.95328d+00, 1.95514d+00, 1.95699d+00, 1.95885d+00, 1.96071d+00,
     2  1.96258d+00, 1.96444d+00, 1.96630d+00, 1.96817d+00, 1.97003d+00,
     3  1.97190d+00, 1.97377d+00, 1.97564d+00, 1.97751d+00, 1.97938d+00,
     4  1.98126d+00, 1.98313d+00, 1.98501d+00, 1.98689d+00, 1.98877d+00,
     5  1.99064d+00, 1.99253d+00, 1.99441d+00, 1.99629d+00, 1.99818d+00,
     6  2.00006d+00, 2.00195d+00, 2.00384d+00, 2.00573d+00, 2.00762d+00,
     7  2.00951d+00, 2.01140d+00, 2.01330d+00, 2.01519d+00, 2.01709d+00,
     8  2.01899d+00, 2.02089d+00, 2.02279d+00, 2.02469d+00, 2.02659d+00,
     9  2.02850d+00, 2.03040d+00, 2.03231d+00, 2.03422d+00, 2.03613d+00/
      data (eqlm(i),i=2929,2976)/ 2.03804d+00, 2.03995d+00, 2.04186d+00,
     1  2.04378d+00, 2.04569d+00, 2.04761d+00, 2.04953d+00, 2.05145d+00,
     2  2.05337d+00, 2.05529d+00, 2.05721d+00, 2.05914d+00, 2.06106d+00,
     3  2.06299d+00, 2.06492d+00, 2.06685d+00, 2.06878d+00, 2.07071d+00,
     4  2.07265d+00, 2.07458d+00, 2.07652d+00, 2.07845d+00, 2.08039d+00,
     5  2.08233d+00, 2.08428d+00, 2.08622d+00, 2.08816d+00, 2.09011d+00,
     6  2.09205d+00, 2.09400d+00, 2.09595d+00, 2.09790d+00, 2.09985d+00,
     7  2.10181d+00, 2.10376d+00, 2.10572d+00, 2.10767d+00, 2.10963d+00,
     8  2.11159d+00, 2.11355d+00, 2.11552d+00, 2.11748d+00, 2.11944d+00,
     9  2.12141d+00, 2.12338d+00, 2.12535d+00, 2.12732d+00, 2.12929d+00/
      data (eqlm(i),i=2977,3024)/ 2.13126d+00, 2.13324d+00, 2.13521d+00,
     1  2.13719d+00, 2.13917d+00, 2.14115d+00, 2.14313d+00, 2.14511d+00,
     2  2.14710d+00, 2.14908d+00, 2.15107d+00, 2.15306d+00, 2.15505d+00,
     3  2.15704d+00, 2.15903d+00, 2.16102d+00, 2.16302d+00, 2.16502d+00,
     4  2.16701d+00, 2.16901d+00, 2.17101d+00, 2.17302d+00, 2.17502d+00,
     5  2.17702d+00, 2.17903d+00, 2.18104d+00, 2.18305d+00, 2.18506d+00,
     6  2.18707d+00, 2.18908d+00, 2.19110d+00, 2.19311d+00, 2.19513d+00,
     7  2.19715d+00, 2.19917d+00, 2.20119d+00, 2.20321d+00, 2.20524d+00,
     8  2.20726d+00, 2.20929d+00, 2.21132d+00, 2.21335d+00, 2.21538d+00,
     9  2.21741d+00, 2.21945d+00, 2.22148d+00, 2.22352d+00, 2.22556d+00/
      data (eqlm(i),i=3025,3072)/ 2.22760d+00, 2.22964d+00, 2.23168d+00,
     1  2.23373d+00, 2.23577d+00, 2.23782d+00, 2.23987d+00, 2.24192d+00,
     2  2.24397d+00, 2.24603d+00, 2.24808d+00, 2.25014d+00, 2.25219d+00,
     3  2.25425d+00, 2.25631d+00, 2.25838d+00, 2.26044d+00, 2.26250d+00,
     4  2.26457d+00, 2.26664d+00, 2.26871d+00, 2.27078d+00, 2.27285d+00,
     5  2.27492d+00, 2.27700d+00, 2.27908d+00, 2.28116d+00, 2.28324d+00,
     6  2.28532d+00, 2.28740d+00, 2.28948d+00, 2.29157d+00, 2.29366d+00,
     7  2.29575d+00, 2.29784d+00, 2.29993d+00, 2.30202d+00, 2.30412d+00,
     8  2.30621d+00, 2.30831d+00, 2.31041d+00, 2.31251d+00, 2.31462d+00,
     9  2.31672d+00, 2.31882d+00, 2.32093d+00, 2.32304d+00, 2.32515d+00/
      data (eqlm(i),i=3073,3120)/ 2.32726d+00, 2.32938d+00, 2.33149d+00,
     1  2.33361d+00, 2.33573d+00, 2.33785d+00, 2.33997d+00, 2.34209d+00,
     2  2.34421d+00, 2.34634d+00, 2.34847d+00, 2.35060d+00, 2.35273d+00,
     3  2.35486d+00, 2.35699d+00, 2.35913d+00, 2.36127d+00, 2.36340d+00,
     4  2.36554d+00, 2.36769d+00, 2.36983d+00, 2.37197d+00, 2.37412d+00,
     5  2.37627d+00, 2.37842d+00, 2.38057d+00, 2.38272d+00, 2.38488d+00,
     6  2.38703d+00, 2.38919d+00, 2.39135d+00, 2.39351d+00, 2.39567d+00,
     7  2.39784d+00, 2.40000d+00, 2.40217d+00, 2.40434d+00, 2.40651d+00,
     8  2.40868d+00, 2.41086d+00, 2.41303d+00, 2.41521d+00, 2.41739d+00,
     9  2.41957d+00, 2.42175d+00, 2.42394d+00, 2.42612d+00, 2.42831d+00/
      data (eqlm(i),i=3121,3168)/ 2.43050d+00, 2.43269d+00, 2.43488d+00,
     1  2.43708d+00, 2.43927d+00, 2.44147d+00, 2.44367d+00, 2.44587d+00,
     2  2.44807d+00, 2.45027d+00, 2.45248d+00, 2.45469d+00, 2.45690d+00,
     3  2.45911d+00, 2.46132d+00, 2.46353d+00, 2.46575d+00, 2.46797d+00,
     4  2.47019d+00, 2.47241d+00, 2.47463d+00, 2.47686d+00, 2.47908d+00,
     5  2.48131d+00, 2.48354d+00, 2.48577d+00, 2.48800d+00, 2.49024d+00,
     6  2.49247d+00, 2.49471d+00, 2.49695d+00, 2.49919d+00, 2.50144d+00,
     7  2.50368d+00, 2.50593d+00, 2.50818d+00, 2.51043d+00, 2.51268d+00,
     8  2.51494d+00, 2.51719d+00, 2.51945d+00, 2.52171d+00, 2.52397d+00,
     9  2.52623d+00, 2.52850d+00, 2.53076d+00, 2.53303d+00, 2.53530d+00/
      data (eqlm(i),i=3169,3216)/ 2.53757d+00, 2.53985d+00, 2.54212d+00,
     1  2.54440d+00, 2.54668d+00, 2.54896d+00, 2.55124d+00, 2.55353d+00,
     2  2.55581d+00, 2.55810d+00, 2.56039d+00, 2.56268d+00, 2.56497d+00,
     3  2.56727d+00, 2.56957d+00, 2.57187d+00, 2.57417d+00, 2.57647d+00,
     4  2.57877d+00, 2.58108d+00, 2.58339d+00, 2.58570d+00, 2.58801d+00,
     5  2.59032d+00, 2.59264d+00, 2.59496d+00, 2.59728d+00, 2.59960d+00,
     6  2.60192d+00, 2.60424d+00, 2.60657d+00, 2.60890d+00, 2.61123d+00,
     7  2.61356d+00, 2.61590d+00, 2.61823d+00, 2.62057d+00, 2.62291d+00,
     8  2.62525d+00, 2.62760d+00, 2.62994d+00, 2.63229d+00, 2.63464d+00,
     9  2.63699d+00, 2.63934d+00, 2.64170d+00, 2.64406d+00, 2.64641d+00/
      data (eqlm(i),i=3217,3264)/ 2.64878d+00, 2.65114d+00, 2.65350d+00,
     1  2.65587d+00, 2.65824d+00, 2.66061d+00, 2.66298d+00, 2.66536d+00,
     2  2.66773d+00, 2.67011d+00, 2.67249d+00, 2.67487d+00, 2.67726d+00,
     3  2.67965d+00, 2.68203d+00, 2.68442d+00, 2.68682d+00, 2.68921d+00,
     4  2.69161d+00, 2.69400d+00, 2.69640d+00, 2.69881d+00, 2.70121d+00,
     5  2.70362d+00, 2.70603d+00, 2.70844d+00, 2.71085d+00, 2.71326d+00,
     6  2.71568d+00, 2.71810d+00, 2.72052d+00, 2.72294d+00, 2.72536d+00,
     7  2.72779d+00, 2.73022d+00, 2.73265d+00, 2.73508d+00, 2.73752d+00,
     8  2.73995d+00, 2.74239d+00, 2.74483d+00, 2.74727d+00, 2.74972d+00,
     9  2.75217d+00, 2.75461d+00, 2.75707d+00, 2.75952d+00, 2.76197d+00/
      data (eqlm(i),i=3265,3312)/ 2.76443d+00, 2.76689d+00, 2.76935d+00,
     1  2.77181d+00, 2.77428d+00, 2.77675d+00, 2.77922d+00, 2.78169d+00,
     2  2.78416d+00, 2.78664d+00, 2.78912d+00, 2.79160d+00, 2.79408d+00,
     3  2.79657d+00, 2.79905d+00, 2.80154d+00, 2.80403d+00, 2.80652d+00,
     4  2.80902d+00, 2.81152d+00, 2.81402d+00, 2.81652d+00, 2.81902d+00,
     5  2.82153d+00, 2.82404d+00, 2.82655d+00, 2.82906d+00, 2.83157d+00,
     6  2.83409d+00, 2.83661d+00, 2.83913d+00, 2.84165d+00, 2.84418d+00,
     7  2.84671d+00, 2.84924d+00, 2.85177d+00, 2.85430d+00, 2.85684d+00,
     8  2.85938d+00, 2.86192d+00, 2.86446d+00, 2.86701d+00, 2.86956d+00,
     9  2.87211d+00, 2.87466d+00, 2.87721d+00, 2.87977d+00, 2.88233d+00/
      data (eqlm(i),i=3313,3360)/ 2.88489d+00, 2.88746d+00, 2.89002d+00,
     1  2.89259d+00, 2.89516d+00, 2.89773d+00, 2.90031d+00, 2.90288d+00,
     2  2.90546d+00, 2.90804d+00, 2.91063d+00, 2.91321d+00, 2.91580d+00,
     3  2.91839d+00, 2.92099d+00, 2.92358d+00, 2.92618d+00, 2.92878d+00,
     4  2.93138d+00, 2.93399d+00, 2.93659d+00, 2.93920d+00, 2.94182d+00,
     5  2.94443d+00, 2.94705d+00, 2.94966d+00, 2.95229d+00, 2.95491d+00,
     6  2.95754d+00, 2.96016d+00, 2.96279d+00, 2.96543d+00, 2.96806d+00,
     7  2.97070d+00, 2.97334d+00, 2.97598d+00, 2.97863d+00, 2.98127d+00,
     8  2.98392d+00, 2.98658d+00, 2.98923d+00, 2.99189d+00, 2.99455d+00,
     9  2.99721d+00, 2.99987d+00, 3.00254d+00, 3.00521d+00, 3.00788d+00/
      data (eqlm(i),i=3361,3408)/ 3.01055d+00, 3.01323d+00, 3.01591d+00,
     1  3.01859d+00, 3.02127d+00, 3.02396d+00, 3.02665d+00, 3.02934d+00,
     2  3.03203d+00, 3.03473d+00, 3.03743d+00, 3.04013d+00, 3.04283d+00,
     3  3.04554d+00, 3.04825d+00, 3.05096d+00, 3.05367d+00, 3.05639d+00,
     4  3.05911d+00, 3.06183d+00, 3.06455d+00, 3.06728d+00, 3.07001d+00,
     5  3.07274d+00, 3.07547d+00, 3.07821d+00, 3.08095d+00, 3.08369d+00,
     6  3.08643d+00, 3.08918d+00, 3.09193d+00, 3.09468d+00, 3.09743d+00,
     7  3.10019d+00, 3.10295d+00, 3.10571d+00, 3.10848d+00, 3.11124d+00,
     8  3.11401d+00, 3.11678d+00, 3.11956d+00, 3.12234d+00, 3.12512d+00,
     9  3.12790d+00, 3.13069d+00, 3.13347d+00, 3.13627d+00, 3.13906d+00/
      data (eqlm(i),i=3409,3456)/ 3.14185d+00, 3.14465d+00, 3.14746d+00,
     1  3.15026d+00, 3.15307d+00, 3.15588d+00, 3.15869d+00, 3.16150d+00,
     2  3.16432d+00, 3.16714d+00, 3.16996d+00, 3.17279d+00, 3.17562d+00,
     3  3.17845d+00, 3.18128d+00, 3.18412d+00, 3.18696d+00, 3.18980d+00,
     4  3.19264d+00, 3.19549d+00, 3.19834d+00, 3.20120d+00, 3.20405d+00,
     5  3.20691d+00, 3.20977d+00, 3.21263d+00, 3.21550d+00, 3.21837d+00,
     6  3.22124d+00, 3.22412d+00, 3.22700d+00, 3.22988d+00, 3.23276d+00,
     7  3.23565d+00, 3.23854d+00, 3.24143d+00, 3.24432d+00, 3.24722d+00,
     8  3.25012d+00, 3.25302d+00, 3.25593d+00, 3.25884d+00, 3.26175d+00,
     9  3.26467d+00, 3.26758d+00, 3.27050d+00, 3.27343d+00, 3.27635d+00/
      data (eqlm(i),i=3457,3504)/ 3.27928d+00, 3.28222d+00, 3.28515d+00,
     1  3.28809d+00, 3.29103d+00, 3.29397d+00, 3.29692d+00, 3.29987d+00,
     2  3.30282d+00, 3.30578d+00, 3.30873d+00, 3.31170d+00, 3.31466d+00,
     3  3.31763d+00, 3.32060d+00, 3.32357d+00, 3.32655d+00, 3.32953d+00,
     4  3.33251d+00, 3.33549d+00, 3.33848d+00, 3.34147d+00, 3.34447d+00,
     5  3.34746d+00, 3.35046d+00, 3.35347d+00, 3.35647d+00, 3.35948d+00,
     6  3.36249d+00, 3.36551d+00, 3.36853d+00, 3.37155d+00, 3.37457d+00,
     7  3.37760d+00, 3.38063d+00, 3.38366d+00, 3.38670d+00, 3.38974d+00,
     8  3.39278d+00, 3.39583d+00, 3.39888d+00, 3.40193d+00, 3.40499d+00,
     9  3.40805d+00, 3.41111d+00, 3.41417d+00, 3.41724d+00, 3.42031d+00/
      data (eqlm(i),i=3505,3552)/ 3.42338d+00, 3.42646d+00, 3.42954d+00,
     1  3.43263d+00, 3.43571d+00, 3.43880d+00, 3.44190d+00, 3.44499d+00,
     2  3.44809d+00, 3.45120d+00, 3.45430d+00, 3.45741d+00, 3.46052d+00,
     3  3.46364d+00, 3.46676d+00, 3.46988d+00, 3.47301d+00, 3.47613d+00,
     4  3.47927d+00, 3.48240d+00, 3.48554d+00, 3.48868d+00, 3.49183d+00,
     5  3.49498d+00, 3.49813d+00, 3.50128d+00, 3.50444d+00, 3.50760d+00,
     6  3.51077d+00, 3.51394d+00, 3.51711d+00, 3.52028d+00, 3.52346d+00,
     7  3.52664d+00, 3.52983d+00, 3.53302d+00, 3.53621d+00, 3.53940d+00,
     8  3.54260d+00, 3.54580d+00, 3.54901d+00, 3.55222d+00, 3.55543d+00,
     9  3.55865d+00, 3.56187d+00, 3.56509d+00, 3.56831d+00, 3.57154d+00/
      data (eqlm(i),i=3553,3600)/ 3.57478d+00, 3.57801d+00, 3.58125d+00,
     1  3.58449d+00, 3.58774d+00, 3.59099d+00, 3.59425d+00, 3.59750d+00,
     2  3.60076d+00, 3.60403d+00, 3.60730d+00, 3.61057d+00, 3.61384d+00,
     3  3.61712d+00, 3.62040d+00, 3.62369d+00, 3.62698d+00, 3.63027d+00,
     4  3.63356d+00, 3.63686d+00, 3.64017d+00, 3.64347d+00, 3.64679d+00,
     5  3.65010d+00, 3.65342d+00, 3.65674d+00, 3.66006d+00, 3.66339d+00,
     6  3.66672d+00, 3.67006d+00, 3.67340d+00, 3.67674d+00, 3.68009d+00,
     7  3.68344d+00, 3.68680d+00, 3.69015d+00, 3.69352d+00, 3.69688d+00,
     8  3.70025d+00, 3.70362d+00, 3.70700d+00, 3.71038d+00, 3.71376d+00,
     9  3.71715d+00, 3.72054d+00, 3.72394d+00, 3.72734d+00, 3.73074d+00/
      data (eqlm(i),i=3601,3648)/ 3.73415d+00, 3.73756d+00, 3.74097d+00,
     1  3.74439d+00, 3.74781d+00, 3.75124d+00, 3.75467d+00, 3.75810d+00,
     2  3.76154d+00, 3.76498d+00, 3.76842d+00, 3.77187d+00, 3.77533d+00,
     3  3.77878d+00, 3.78224d+00, 3.78571d+00, 3.78918d+00, 3.79265d+00,
     4  3.79613d+00, 3.79961d+00, 3.80309d+00, 3.80658d+00, 3.81007d+00,
     5  3.81357d+00, 3.81707d+00, 3.82057d+00, 3.82408d+00, 3.82759d+00,
     6  3.83111d+00, 3.83463d+00, 3.83815d+00, 3.84168d+00, 3.84521d+00,
     7  3.84875d+00, 3.85229d+00, 3.85583d+00, 3.85938d+00, 3.86293d+00,
     8  3.86649d+00, 3.87005d+00, 3.87362d+00, 3.87718d+00, 3.88076d+00,
     9  3.88433d+00, 3.88792d+00, 3.89150d+00, 3.89509d+00, 3.89868d+00/
      data (eqlm(i),i=3649,3696)/ 3.90228d+00, 3.90589d+00, 3.90949d+00,
     1  3.91310d+00, 3.91672d+00, 3.92034d+00, 3.92396d+00, 3.92759d+00,
     2  3.93122d+00, 3.93485d+00, 3.93849d+00, 3.94214d+00, 3.94579d+00,
     3  3.94944d+00, 3.95310d+00, 3.95676d+00, 3.96043d+00, 3.96410d+00,
     4  3.96777d+00, 3.97145d+00, 3.97513d+00, 3.97882d+00, 3.98251d+00,
     5  3.98621d+00, 3.98991d+00, 3.99362d+00, 3.99733d+00, 4.00104d+00,
     6  4.00476d+00, 4.00848d+00, 4.01221d+00, 4.01594d+00, 4.01968d+00,
     7  4.02342d+00, 4.02716d+00, 4.03091d+00, 4.03467d+00, 4.03843d+00,
     8  4.04219d+00, 4.04596d+00, 4.04973d+00, 4.05351d+00, 4.05729d+00,
     9  4.06108d+00, 4.06487d+00, 4.06866d+00, 4.07246d+00, 4.07627d+00/
      data (eqlm(i),i=3697,3744)/ 4.08008d+00, 4.08389d+00, 4.08771d+00,
     1  4.09153d+00, 4.09536d+00, 4.09919d+00, 4.10303d+00, 4.10687d+00,
     2  4.11072d+00, 4.11457d+00, 4.11842d+00, 4.12228d+00, 4.12615d+00,
     3  4.13002d+00, 4.13389d+00, 4.13777d+00, 4.14166d+00, 4.14555d+00,
     4  4.14944d+00, 4.15334d+00, 4.15724d+00, 4.16115d+00, 4.16507d+00,
     5  4.16898d+00, 4.17291d+00, 4.17683d+00, 4.18077d+00, 4.18470d+00,
     6  4.18865d+00, 4.19260d+00, 4.19655d+00, 4.20050d+00, 4.20447d+00,
     7  4.20843d+00, 4.21241d+00, 4.21638d+00, 4.22037d+00, 4.22435d+00,
     8  4.22835d+00, 4.23234d+00, 4.23634d+00, 4.24035d+00, 4.24436d+00,
     9  4.24838d+00, 4.25240d+00, 4.25643d+00, 4.26046d+00, 4.26450d+00/
      data (eqlm(i),i=3745,3792)/ 4.26855d+00, 4.27259d+00, 4.27665d+00,
     1  4.28070d+00, 4.28477d+00, 4.28884d+00, 4.29291d+00, 4.29699d+00,
     2  4.30107d+00, 4.30516d+00, 4.30926d+00, 4.31336d+00, 4.31746d+00,
     3  4.32157d+00, 4.32569d+00, 4.32981d+00, 4.33394d+00, 4.33807d+00,
     4  4.34221d+00, 4.34635d+00, 4.35050d+00, 4.35465d+00, 4.35881d+00,
     5  4.36297d+00, 4.36714d+00, 4.37132d+00, 4.37550d+00, 4.37968d+00,
     6  4.38387d+00, 4.38807d+00, 4.39227d+00, 4.39648d+00, 4.40069d+00,
     7  4.40491d+00, 4.40913d+00, 4.41336d+00, 4.41760d+00, 4.42184d+00,
     8  4.42608d+00, 4.43034d+00, 4.43459d+00, 4.43886d+00, 4.44312d+00,
     9  4.44740d+00, 4.45168d+00, 4.45596d+00, 4.46026d+00, 4.46455d+00/
      data (eqlm(i),i=3793,3840)/ 4.46885d+00, 4.47316d+00, 4.47748d+00,
     1  4.48180d+00, 4.48612d+00, 4.49045d+00, 4.49479d+00, 4.49913d+00,
     2  4.50348d+00, 4.50784d+00, 4.51220d+00, 4.51656d+00, 4.52094d+00,
     3  4.52531d+00, 4.52970d+00, 4.53409d+00, 4.53848d+00, 4.54289d+00,
     4  4.54729d+00, 4.55171d+00, 4.55613d+00, 4.56055d+00, 4.56498d+00,
     5  4.56942d+00, 4.57386d+00, 4.57831d+00, 4.58277d+00, 4.58723d+00,
     6  4.59170d+00, 4.59617d+00, 4.60065d+00, 4.60514d+00, 4.60963d+00,
     7  4.61413d+00, 4.61864d+00, 4.62315d+00, 4.62766d+00, 4.63219d+00,
     8  4.63672d+00, 4.64125d+00, 4.64579d+00, 4.65034d+00, 4.65490d+00,
     9  4.65946d+00, 4.66403d+00, 4.66860d+00, 4.67318d+00, 4.67777d+00/
      data (eqlm(i),i=3841,3888)/ 4.68236d+00, 4.68696d+00, 4.69156d+00,
     1  4.69617d+00, 4.70079d+00, 4.70542d+00, 4.71005d+00, 4.71469d+00,
     2  4.71933d+00, 4.72398d+00, 4.72864d+00, 4.73330d+00, 4.73797d+00,
     3  4.74265d+00, 4.74733d+00, 4.75202d+00, 4.75672d+00, 4.76142d+00,
     4  4.76613d+00, 4.77085d+00, 4.77557d+00, 4.78030d+00, 4.78504d+00,
     5  4.78978d+00, 4.79453d+00, 4.79929d+00, 4.80406d+00, 4.80883d+00,
     6  4.81360d+00, 4.81839d+00, 4.82318d+00, 4.82798d+00, 4.83278d+00,
     7  4.83759d+00, 4.84241d+00, 4.84724d+00, 4.85207d+00, 4.85691d+00,
     8  4.86176d+00, 4.86661d+00, 4.87147d+00, 4.87634d+00, 4.88122d+00,
     9  4.88610d+00, 4.89099d+00, 4.89588d+00, 4.90079d+00, 4.90570d+00/
      data (eqlm(i),i=3889,3936)/ 4.91061d+00, 4.91554d+00, 4.92047d+00,
     1  4.92541d+00, 4.93036d+00, 4.93531d+00, 4.94027d+00, 4.94524d+00,
     2  4.95022d+00, 4.95520d+00, 4.96019d+00, 4.96519d+00, 4.97019d+00,
     3  4.97521d+00, 4.98023d+00, 4.98525d+00, 4.99029d+00, 4.99533d+00,
     4  5.00038d+00, 5.00544d+00, 5.01050d+00, 5.01558d+00, 5.02066d+00,
     5  5.02574d+00, 5.03084d+00, 5.03594d+00, 5.04105d+00, 5.04617d+00,
     6  5.05130d+00, 5.05643d+00, 5.06157d+00, 5.06672d+00, 5.07188d+00,
     7  5.07704d+00, 5.08221d+00, 5.08739d+00, 5.09258d+00, 5.09778d+00,
     8  5.10298d+00, 5.10819d+00, 5.11341d+00, 5.11864d+00, 5.12388d+00,
     9  5.12912d+00, 5.13437d+00, 5.13963d+00, 5.14490d+00, 5.15018d+00/
      data (eqlm(i),i=3937,3984)/ 5.15546d+00, 5.16075d+00, 5.16605d+00,
     1  5.17136d+00, 5.17668d+00, 5.18200d+00, 5.18734d+00, 5.19268d+00,
     2  5.19803d+00, 5.20338d+00, 5.20875d+00, 5.21412d+00, 5.21951d+00,
     3  5.22490d+00, 5.23030d+00, 5.23571d+00, 5.24112d+00, 5.24655d+00,
     4  5.25198d+00, 5.25742d+00, 5.26287d+00, 5.26833d+00, 5.27380d+00,
     5  5.27927d+00, 5.28476d+00, 5.29025d+00, 5.29576d+00, 5.30127d+00,
     6  5.30679d+00, 5.31231d+00, 5.31785d+00, 5.32340d+00, 5.32895d+00,
     7  5.33451d+00, 5.34009d+00, 5.34567d+00, 5.35126d+00, 5.35686d+00,
     8  5.36246d+00, 5.36808d+00, 5.37371d+00, 5.37934d+00, 5.38498d+00,
     9  5.39064d+00, 5.39630d+00, 5.40197d+00, 5.40765d+00, 5.41334d+00/
      data (eqlm(i),i=3985,4032)/ 5.41904d+00, 5.42475d+00, 5.43046d+00,
     1  5.43619d+00, 5.44193d+00, 5.44767d+00, 5.45343d+00, 5.45919d+00,
     2  5.46496d+00, 5.47074d+00, 5.47654d+00, 5.48234d+00, 5.48815d+00,
     3  5.49397d+00, 5.49980d+00, 5.50564d+00, 5.51149d+00, 5.51735d+00,
     4  5.52321d+00, 5.52909d+00, 5.53498d+00, 5.54088d+00, 5.54679d+00,
     5  5.55270d+00, 5.55863d+00, 5.56457d+00, 5.57051d+00, 5.57647d+00,
     6  5.58243d+00, 5.58841d+00, 5.59440d+00, 5.60039d+00, 5.60640d+00,
     7  5.61242d+00, 5.61844d+00, 5.62448d+00, 5.63053d+00, 5.63658d+00,
     8  5.64265d+00, 5.64873d+00, 5.65482d+00, 5.66091d+00, 5.66702d+00,
     9  5.67314d+00, 5.67927d+00, 5.68541d+00, 5.69156d+00, 5.69772d+00/
      data (eqlm(i),i=4033,4080)/ 5.70389d+00, 5.71007d+00, 5.71626d+00,
     1  5.72247d+00, 5.72868d+00, 5.73490d+00, 5.74114d+00, 5.74738d+00,
     2  5.75364d+00, 5.75991d+00, 5.76618d+00, 5.77247d+00, 5.77877d+00,
     3  5.78508d+00, 5.79140d+00, 5.79774d+00, 5.80408d+00, 5.81043d+00,
     4  5.81680d+00, 5.82317d+00, 5.82956d+00, 5.83596d+00, 5.84237d+00,
     5  5.84879d+00, 5.85522d+00, 5.86167d+00, 5.86812d+00, 5.87459d+00,
     6  5.88106d+00, 5.88755d+00, 5.89405d+00, 5.90057d+00, 5.90709d+00,
     7  5.91362d+00, 5.92017d+00, 5.92673d+00, 5.93330d+00, 5.93988d+00,
     8  5.94647d+00, 5.95308d+00, 5.95969d+00, 5.96632d+00, 5.97296d+00,
     9  5.97961d+00, 5.98628d+00, 5.99295d+00, 5.99964d+00, 6.00634d+00/
      data (eqlm(i),i=4081,4128)/ 6.01306d+00, 6.01978d+00, 6.02652d+00,
     1  6.03326d+00, 6.04003d+00, 6.04680d+00, 6.05358d+00, 6.06038d+00,
     2  6.06719d+00, 6.07402d+00, 6.08085d+00, 6.08770d+00, 6.09456d+00,
     3  6.10143d+00, 6.10832d+00, 6.11521d+00, 6.12212d+00, 6.12905d+00,
     4  6.13598d+00, 6.14293d+00, 6.14989d+00, 6.15687d+00, 6.16385d+00,
     5  6.17085d+00, 6.17787d+00, 6.18489d+00, 6.19193d+00, 6.19899d+00,
     6  6.20605d+00, 6.21313d+00, 6.22022d+00, 6.22733d+00, 6.23444d+00,
     7  6.24158d+00, 6.24872d+00, 6.25588d+00, 6.26305d+00, 6.27024d+00,
     8  6.27744d+00, 6.28465d+00, 6.29188d+00, 6.29912d+00, 6.30637d+00,
     9  6.31364d+00, 6.32092d+00, 6.32821d+00, 6.33552d+00, 6.34284d+00/
      data (eqlm(i),i=4129,4176)/ 6.35018d+00, 6.35753d+00, 6.36490d+00,
     1  6.37228d+00, 6.37967d+00, 6.38708d+00, 6.39450d+00, 6.40193d+00,
     2  6.40939d+00, 6.41685d+00, 6.42433d+00, 6.43182d+00, 6.43933d+00,
     3  6.44685d+00, 6.45439d+00, 6.46194d+00, 6.46951d+00, 6.47709d+00,
     4  6.48469d+00, 6.49230d+00, 6.49993d+00, 6.50757d+00, 6.51522d+00,
     5  6.52289d+00, 6.53058d+00, 6.53828d+00, 6.54600d+00, 6.55373d+00,
     6  6.56148d+00, 6.56924d+00, 6.57702d+00, 6.58481d+00, 6.59262d+00,
     7  6.60044d+00, 6.60828d+00, 6.61614d+00, 6.62401d+00, 6.63190d+00,
     8  6.63980d+00, 6.64772d+00, 6.65565d+00, 6.66360d+00, 6.67157d+00,
     9  6.67955d+00, 6.68755d+00, 6.69557d+00, 6.70360d+00, 6.71164d+00/
      data (eqlm(i),i=4177,4224)/ 6.71971d+00, 6.72779d+00, 6.73588d+00,
     1  6.74400d+00, 6.75213d+00, 6.76027d+00, 6.76844d+00, 6.77661d+00,
     2  6.78481d+00, 6.79302d+00, 6.80125d+00, 6.80950d+00, 6.81777d+00,
     3  6.82605d+00, 6.83434d+00, 6.84266d+00, 6.85099d+00, 6.85934d+00,
     4  6.86771d+00, 6.87609d+00, 6.88450d+00, 6.89291d+00, 6.90135d+00,
     5  6.90981d+00, 6.91828d+00, 6.92677d+00, 6.93528d+00, 6.94380d+00,
     6  6.95235d+00, 6.96091d+00, 6.96949d+00, 6.97808d+00, 6.98670d+00,
     7  6.99534d+00, 7.00399d+00, 7.01266d+00, 7.02135d+00, 7.03006d+00,
     8  7.03878d+00, 7.04753d+00, 7.05629d+00, 7.06508d+00, 7.07388d+00,
     9  7.08270d+00, 7.09154d+00, 7.10040d+00, 7.10927d+00, 7.11817d+00/
      data (eqlm(i),i=4225,4272)/ 7.12709d+00, 7.13602d+00, 7.14498d+00,
     1  7.15395d+00, 7.16294d+00, 7.17196d+00, 7.18099d+00, 7.19004d+00,
     2  7.19912d+00, 7.20821d+00, 7.21732d+00, 7.22645d+00, 7.23561d+00,
     3  7.24478d+00, 7.25397d+00, 7.26318d+00, 7.27242d+00, 7.28167d+00,
     4  7.29094d+00, 7.30024d+00, 7.30956d+00, 7.31889d+00, 7.32825d+00,
     5  7.33763d+00, 7.34703d+00, 7.35645d+00, 7.36589d+00, 7.37535d+00,
     6  7.38483d+00, 7.39434d+00, 7.40386d+00, 7.41341d+00, 7.42298d+00,
     7  7.43257d+00, 7.44218d+00, 7.45182d+00, 7.46147d+00, 7.47115d+00,
     8  7.48085d+00, 7.49057d+00, 7.50032d+00, 7.51009d+00, 7.51987d+00,
     9  7.52969d+00, 7.53952d+00, 7.54938d+00, 7.55926d+00, 7.56916d+00/
      data (eqlm(i),i=4273,4320)/ 7.57909d+00, 7.58903d+00, 7.59901d+00,
     1  7.60900d+00, 7.61902d+00, 7.62906d+00, 7.63912d+00, 7.64921d+00,
     2  7.65932d+00, 7.66946d+00, 7.67962d+00, 7.68980d+00, 7.70001d+00,
     3  7.71024d+00, 7.72049d+00, 7.73077d+00, 7.74108d+00, 7.75140d+00,
     4  7.76176d+00, 7.77213d+00, 7.78253d+00, 7.79296d+00, 7.80341d+00,
     5  7.81389d+00, 7.82439d+00, 7.83492d+00, 7.84547d+00, 7.85604d+00,
     6  7.86665d+00, 7.87727d+00, 7.88793d+00, 7.89861d+00, 7.90931d+00,
     7  7.92004d+00, 7.93080d+00, 7.94159d+00, 7.95240d+00, 7.96323d+00,
     8  7.97410d+00, 7.98499d+00, 7.99590d+00, 8.00684d+00, 8.01781d+00,
     9  8.02881d+00, 8.03984d+00, 8.05089d+00, 8.06196d+00, 8.07307d+00/
      data (eqlm(i),i=4321,4368)/ 8.08420d+00, 8.09536d+00, 8.10655d+00,
     1  8.11777d+00, 8.12901d+00, 8.14029d+00, 8.15159d+00, 8.16292d+00,
     2  8.17427d+00, 8.18566d+00, 8.19708d+00, 8.20852d+00, 8.21999d+00,
     3  8.23149d+00, 8.24302d+00, 8.25458d+00, 8.26617d+00, 8.27779d+00,
     4  8.28944d+00, 8.30112d+00, 8.31283d+00, 8.32456d+00, 8.33633d+00,
     5  8.34813d+00, 8.35996d+00, 8.37182d+00, 8.38371d+00, 8.39563d+00,
     6  8.40758d+00, 8.41956d+00, 8.43158d+00, 8.44362d+00, 8.45570d+00,
     7  8.46781d+00, 8.47995d+00, 8.49212d+00, 8.50432d+00, 8.51655d+00,
     8  8.52882d+00, 8.54112d+00, 8.55345d+00, 8.56582d+00, 8.57821d+00,
     9  8.59064d+00, 8.60311d+00, 8.61560d+00, 8.62813d+00, 8.64069d+00/
      data (eqlm(i),i=4369,4416)/ 8.65329d+00, 8.66592d+00, 8.67859d+00,
     1  8.69128d+00, 8.70402d+00, 8.71678d+00, 8.72958d+00, 8.74242d+00,
     2  8.75529d+00, 8.76819d+00, 8.78113d+00, 8.79411d+00, 8.80712d+00,
     3  8.82017d+00, 8.83325d+00, 8.84637d+00, 8.85952d+00, 8.87271d+00,
     4  8.88594d+00, 8.89920d+00, 8.91250d+00, 8.92583d+00, 8.93921d+00,
     5  8.95262d+00, 8.96606d+00, 8.97955d+00, 8.99307d+00, 9.00663d+00,
     6  9.02023d+00, 9.03387d+00, 9.04754d+00, 9.06125d+00, 9.07501d+00,
     7  9.08880d+00, 9.10263d+00, 9.11649d+00, 9.13040d+00, 9.14435d+00,
     8  9.15834d+00, 9.17236d+00, 9.18643d+00, 9.20054d+00, 9.21469d+00,
     9  9.22887d+00, 9.24310d+00, 9.25737d+00, 9.27169d+00, 9.28604d+00/
      data (eqlm(i),i=4417,4464)/ 9.30043d+00, 9.31487d+00, 9.32935d+00,
     1  9.34387d+00, 9.35843d+00, 9.37304d+00, 9.38768d+00, 9.40237d+00,
     2  9.41711d+00, 9.43189d+00, 9.44671d+00, 9.46157d+00, 9.47648d+00,
     3  9.49143d+00, 9.50643d+00, 9.52147d+00, 9.53655d+00, 9.55168d+00,
     4  9.56686d+00, 9.58208d+00, 9.59735d+00, 9.61266d+00, 9.62802d+00,
     5  9.64342d+00, 9.65887d+00, 9.67437d+00, 9.68992d+00, 9.70551d+00,
     6  9.72115d+00, 9.73684d+00, 9.75257d+00, 9.76835d+00, 9.78418d+00,
     7  9.80006d+00, 9.81599d+00, 9.83197d+00, 9.84799d+00, 9.86407d+00,
     8  9.88020d+00, 9.89637d+00, 9.91260d+00, 9.92887d+00, 9.94520d+00,
     9  9.96158d+00, 9.97801d+00, 9.99449d+00, 1.00110d+01, 1.00276d+01/
      data (eqlm(i),i=4465,4512)/ 1.00443d+01, 1.00609d+01, 1.00777d+01,
     1  1.00945d+01, 1.01113d+01, 1.01282d+01, 1.01452d+01, 1.01622d+01,
     2  1.01793d+01, 1.01964d+01, 1.02136d+01, 1.02308d+01, 1.02481d+01,
     3  1.02654d+01, 1.02828d+01, 1.03003d+01, 1.03178d+01, 1.03353d+01,
     4  1.03530d+01, 1.03706d+01, 1.03884d+01, 1.04062d+01, 1.04240d+01,
     5  1.04420d+01, 1.04599d+01, 1.04780d+01, 1.04961d+01, 1.05142d+01,
     6  1.05324d+01, 1.05507d+01, 1.05690d+01, 1.05874d+01, 1.06059d+01,
     7  1.06244d+01, 1.06430d+01, 1.06616d+01, 1.06803d+01, 1.06991d+01,
     8  1.07179d+01, 1.07368d+01, 1.07558d+01, 1.07748d+01, 1.07939d+01,
     9  1.08131d+01, 1.08323d+01, 1.08516d+01, 1.08710d+01, 1.08904d+01/
      data (eqlm(i),i=4513,4560)/ 1.09099d+01, 1.09294d+01, 1.09491d+01,
     1  1.09688d+01, 1.09885d+01, 1.10083d+01, 1.10282d+01, 1.10482d+01,
     2  1.10683d+01, 1.10884d+01, 1.11085d+01, 1.11288d+01, 1.11491d+01,
     3  1.11695d+01, 1.11900d+01, 1.12105d+01, 1.12311d+01, 1.12518d+01,
     4  1.12726d+01, 1.12934d+01, 1.13143d+01, 1.13353d+01, 1.13564d+01,
     5  1.13775d+01, 1.13988d+01, 1.14200d+01, 1.14414d+01, 1.14629d+01,
     6  1.14844d+01, 1.15060d+01, 1.15277d+01, 1.15495d+01, 1.15713d+01,
     7  1.15932d+01, 1.16152d+01, 1.16373d+01, 1.16595d+01, 1.16817d+01,
     8  1.17041d+01, 1.17265d+01, 1.17490d+01, 1.17716d+01, 1.17943d+01,
     9  1.18170d+01, 1.18399d+01, 1.18628d+01, 1.18858d+01, 1.19089d+01/
      data (eqlm(i),i=4561,4608)/ 1.19321d+01, 1.19554d+01, 1.19788d+01,
     1  1.20023d+01, 1.20258d+01, 1.20495d+01, 1.20732d+01, 1.20970d+01,
     2  1.21210d+01, 1.21450d+01, 1.21691d+01, 1.21933d+01, 1.22176d+01,
     3  1.22420d+01, 1.22665d+01, 1.22910d+01, 1.23157d+01, 1.23405d+01,
     4  1.23654d+01, 1.23904d+01, 1.24154d+01, 1.24406d+01, 1.24659d+01,
     5  1.24913d+01, 1.25168d+01, 1.25423d+01, 1.25680d+01, 1.25938d+01,
     6  1.26197d+01, 1.26457d+01, 1.26718d+01, 1.26981d+01, 1.27244d+01,
     7  1.27508d+01, 1.27773d+01, 1.28040d+01, 1.28308d+01, 1.28576d+01,
     8  1.28846d+01, 1.29117d+01, 1.29389d+01, 1.29662d+01, 1.29937d+01,
     9  1.30212d+01, 1.30489d+01, 1.30767d+01, 1.31046d+01, 1.31326d+01/
      data (eqlm(i),i=4609,4656)/ 1.31608d+01, 1.31890d+01, 1.32174d+01,
     1  1.32459d+01, 1.32745d+01, 1.33033d+01, 1.33322d+01, 1.33612d+01,
     2  1.33903d+01, 1.34195d+01, 1.34489d+01, 1.34784d+01, 1.35081d+01,
     3  1.35378d+01, 1.35677d+01, 1.35978d+01, 1.36279d+01, 1.36582d+01,
     4  1.36886d+01, 1.37192d+01, 1.37499d+01, 1.37808d+01, 1.38117d+01,
     5  1.38429d+01, 1.38741d+01, 1.39055d+01, 1.39371d+01, 1.39688d+01,
     6  1.40006d+01, 1.40326d+01, 1.40647d+01, 1.40970d+01, 1.41294d+01,
     7  1.41620d+01, 1.41947d+01, 1.42275d+01, 1.42606d+01, 1.42938d+01,
     8  1.43271d+01, 1.43606d+01, 1.43942d+01, 1.44280d+01, 1.44620d+01,
     9  1.44961d+01, 1.45304d+01, 1.45649d+01, 1.45995d+01, 1.46343d+01/
      data (eqlm(i),i=4657,4704)/ 1.46692d+01, 1.47044d+01, 1.47397d+01,
     1  1.47751d+01, 1.48108d+01, 1.48466d+01, 1.48825d+01, 1.49187d+01,
     2  1.49550d+01, 1.49915d+01, 1.50282d+01, 1.50651d+01, 1.51022d+01,
     3  1.51394d+01, 1.51769d+01, 1.52145d+01, 1.52523d+01, 1.52903d+01,
     4  1.53285d+01, 1.53668d+01, 1.54054d+01, 1.54442d+01, 1.54832d+01,
     5  1.55223d+01, 1.55617d+01, 1.56013d+01, 1.56410d+01, 1.56810d+01,
     6  1.57212d+01, 1.57616d+01, 1.58022d+01, 1.58430d+01, 1.58841d+01,
     7  1.59253d+01, 1.59668d+01, 1.60085d+01, 1.60504d+01, 1.60925d+01,
     8  1.61349d+01, 1.61775d+01, 1.62203d+01, 1.62633d+01, 1.63066d+01,
     9  1.63501d+01, 1.63939d+01, 1.64379d+01, 1.64821d+01, 1.65265d+01/
      data (eqlm(i),i=4705,4752)/ 1.65713d+01, 1.66162d+01, 1.66614d+01,
     1  1.67069d+01, 1.67526d+01, 1.67986d+01, 1.68448d+01, 1.68913d+01,
     2  1.69380d+01, 1.69851d+01, 1.70323d+01, 1.70799d+01, 1.71277d+01,
     3  1.71758d+01, 1.72242d+01, 1.72728d+01, 1.73217d+01, 1.73710d+01,
     4  1.74205d+01, 1.74702d+01, 1.75203d+01, 1.75707d+01, 1.76214d+01,
     5  1.76723d+01, 1.77236d+01, 1.77752d+01, 1.78271d+01, 1.78793d+01,
     6  1.79318d+01, 1.79846d+01, 1.80377d+01, 1.80912d+01, 1.81450d+01,
     7  1.81991d+01, 1.82535d+01, 1.83083d+01, 1.83634d+01, 1.84189d+01,
     8  1.84747d+01, 1.85309d+01, 1.85874d+01, 1.86442d+01, 1.87014d+01,
     9  1.87590d+01, 1.88170d+01, 1.88753d+01, 1.89340d+01, 1.89930d+01/
      data (eqlm(i),i=4753,4800)/ 1.90524d+01, 1.91123d+01, 1.91725d+01,
     1  1.92331d+01, 1.92941d+01, 1.93554d+01, 1.94172d+01, 1.94794d+01,
     2  1.95420d+01, 1.96051d+01, 1.96685d+01, 1.97324d+01, 1.97967d+01,
     3  1.98614d+01, 1.99266d+01, 1.99922d+01, 2.00582d+01, 2.01247d+01,
     4  2.01917d+01, 2.02591d+01, 2.03269d+01, 2.03953d+01, 2.04641d+01,
     5  2.05334d+01, 2.06032d+01, 2.06735d+01, 2.07443d+01, 2.08155d+01,
     6  2.08873d+01, 2.09596d+01, 2.10324d+01, 2.11058d+01, 2.11796d+01,
     7  2.12540d+01, 2.13290d+01, 2.14044d+01, 2.14805d+01, 2.15571d+01,
     8  2.16342d+01, 2.17120d+01, 2.17903d+01, 2.18692d+01, 2.19487d+01,
     9  2.20288d+01, 2.21095d+01, 2.21908d+01, 2.22727d+01, 2.23553d+01/
      data (eqlm(i),i=4801,4848)/ 2.24385d+01, 2.25223d+01, 2.26068d+01,
     1  2.26919d+01, 2.27777d+01, 2.28642d+01, 2.29513d+01, 2.30392d+01,
     2  2.31277d+01, 2.32170d+01, 2.33070d+01, 2.33977d+01, 2.34891d+01,
     3  2.35812d+01, 2.36742d+01, 2.37678d+01, 2.38623d+01, 2.39575d+01,
     4  2.40535d+01, 2.41504d+01, 2.42480d+01, 2.43464d+01, 2.44457d+01,
     5  2.45459d+01, 2.46468d+01, 2.47487d+01, 2.48514d+01, 2.49550d+01,
     6  2.50595d+01, 2.51649d+01, 2.52712d+01, 2.53785d+01, 2.54867d+01,
     7  2.55959d+01, 2.57061d+01, 2.58172d+01, 2.59293d+01, 2.60425d+01,
     8  2.61567d+01, 2.62719d+01, 2.63882d+01, 2.65056d+01, 2.66240d+01,
     9  2.67435d+01, 2.68642d+01, 2.69860d+01, 2.71090d+01, 2.72331d+01/
      data (eqlm(i),i=4849,4896)/ 2.73584d+01, 2.74849d+01, 2.76126d+01,
     1  2.77416d+01, 2.78718d+01, 2.80033d+01, 2.81361d+01, 2.82702d+01,
     2  2.84056d+01, 2.85424d+01, 2.86806d+01, 2.88202d+01, 2.89611d+01,
     3  2.91036d+01, 2.92474d+01, 2.93928d+01, 2.95397d+01, 2.96881d+01,
     4  2.98380d+01, 2.99896d+01, 3.01427d+01, 3.02975d+01, 3.04540d+01,
     5  3.06121d+01, 3.07720d+01, 3.09336d+01, 3.10969d+01, 3.12621d+01,
     6  3.14291d+01, 3.15980d+01, 3.17687d+01, 3.19414d+01, 3.21161d+01,
     7  3.22927d+01, 3.24714d+01, 3.26522d+01, 3.28350d+01, 3.30200d+01,
     8  3.32072d+01, 3.33965d+01, 3.35882d+01, 3.37821d+01, 3.39784d+01,
     9  3.41771d+01, 3.43782d+01, 3.45817d+01, 3.47878d+01, 3.49965d+01/
      data (eqlm(i),i=4897,4944)/ 3.52078d+01, 3.54217d+01, 3.56384d+01,
     1  3.58578d+01, 3.60801d+01, 3.63052d+01, 3.65333d+01, 3.67644d+01,
     2  3.69985d+01, 3.72358d+01, 3.74762d+01, 3.77199d+01, 3.79669d+01,
     3  3.82173d+01, 3.84711d+01, 3.87285d+01, 3.89895d+01, 3.92541d+01,
     4  3.95226d+01, 3.97948d+01, 4.00710d+01, 4.03512d+01, 4.06355d+01,
     5  4.09239d+01, 4.12167d+01, 4.15139d+01, 4.18155d+01, 4.21217d+01,
     6  4.24326d+01, 4.27483d+01, 4.30690d+01, 4.33946d+01, 4.37254d+01,
     7  4.40615d+01, 4.44030d+01, 4.47501d+01, 4.51028d+01, 4.54613d+01,
     8  4.58258d+01, 4.61964d+01, 4.65733d+01, 4.69566d+01, 4.73465d+01,
     9  4.77432d+01, 4.81469d+01, 4.85577d+01, 4.89758d+01, 4.94015d+01/
      data (eqlm(i),i=4945,4992)/ 4.98349d+01, 5.02762d+01, 5.07258d+01,
     1  5.11837d+01, 5.16503d+01, 5.21258d+01, 5.26105d+01, 5.31046d+01,
     2  5.36084d+01, 5.41222d+01, 5.46463d+01, 5.51810d+01, 5.57267d+01,
     3  5.62837d+01, 5.68523d+01, 5.74329d+01, 5.80260d+01, 5.86318d+01,
     4  5.92509d+01, 5.98836d+01, 6.05305d+01, 6.11921d+01, 6.18687d+01,
     5  6.25609d+01, 6.32694d+01, 6.39947d+01, 6.47373d+01, 6.54980d+01,
     6  6.62774d+01, 6.70762d+01, 6.78951d+01, 6.87349d+01, 6.95965d+01,
     7  7.04806d+01, 7.13882d+01, 7.23203d+01, 7.32778d+01, 7.42619d+01,
     8  7.52736d+01, 7.63141d+01, 7.73847d+01, 7.84867d+01, 7.96215d+01,
     9  8.07907d+01, 8.19957d+01, 8.32384d+01, 8.45205d+01, 8.58438d+01/
      data (eqlm(i),i=4993,5001)/ 8.72105d+01, 8.86227d+01, 9.00828d+01,
     1  9.15931d+01, 9.31565d+01, 9.47757d+01, 9.64538d+01, 9.81940d+01,
     2  1.00000d+02/
      end
