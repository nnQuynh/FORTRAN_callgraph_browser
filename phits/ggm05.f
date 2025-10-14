************************************************************************
*                                                                      *
      subroutine sctneu(eein,wgti,ireg,imat,ffac,icmm,elrt)
*                                                                      *
*        calculate the collision of a neutron with a nucleus.          *
*        last modified by K.Niita on 2016/08/08                        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for #ifdef REDUCTION_COUNTER
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2

      use GGMBANKMOD !FURUTA
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
      common /eparm/  esmax, esmin, emin(20)


      common /fsmul/  ridfs, iidfs
      common /elreg/  elarg(kvlmax,4), mlrgn, melrg(kvlmax,2), ielusr
      common /regdc/  idrg(kvlmax), idgr(kvmmax)


*-----------------------------------------------------------------------

      dimension vr(3)
      dimension usave(3)

*-----------------------------------------------------------------------

           nclst  = -1
           nclsts = 0
           nixcos = 0  ! MATSUDA 2022.12.12  for T-Point

*-----------------------------------------------------------------------

            nter = 0
            ipt  = 1
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
            mk   = imat
            m    = jmd(1+mk)
            m1   = m

*-----------------------------------------------------------------------
*     sample the nuclide, using the cumulative total cross section.
*-----------------------------------------------------------------------

         if( npq(mk) .eq. 1 ) goto 20

*-----------------------------------------------------------------------

         if( icmm .ne. 0 ) then

               m = m1 + icmm - 1

         else

               c = rang() * totm

            do m = m1, jmd(1+mk+1) - 2

               c = c - rtc(5,lme(1,m)) * fme(m)

               if( c .lt. 0. ) goto 20

            end do

        end if

*-----------------------------------------------------------------------

   20       iex  = lme(1,m)
            iexp = m
            mtp  = 2
            mpan = ipan(icl) + m - m1

*-----------------------------------------------------------------------
*           iet ; S(a,b)
*-----------------------------------------------------------------------

            iet  = lmt(m)
               if( iet.ne.0 ) then
                 if( erg .gt. esa(iet) )   iet = 0
               end if
               if( rtc(11,iex) .ge. 0. ) iet = 0

*-----------------------------------------------------------------------
*           mother nucleus
*-----------------------------------------------------------------------

            mathz = iza(iexp) / 1000
            mathn = iza(iexp) - 1000 * mathz - mathz
            izaid = mathz * 1000 + mathn + mathz

*-----------------------------------------------------------------------
*     user elastic and inelastic option if exist
*-----------------------------------------------------------------------

            ioel = 0
            icels = 1

      if( mlrgn .gt. 0 ) then

            do jj = 1, mlrgn

              do kk = melrg(jj,1), melrg(jj,2)

                 if( ireg .eq. idgr(kk) ) then

                    ioel = jj
                    goto 300

                 end if

              end do

            end do

  300       continue

      end if

*-----------------------------------------------------------------------
*     capture cross sections by user defined function
*-----------------------------------------------------------------------

      if( ioel .ne. 0 ) then

         if( ielusr .eq. 1 ) then

            call usrelst1(3,icels,ioel,erg,sigca,sigal)

         else

            call usrelst2(3,icels,ioel,erg,sigca,sigal)

         end if

      end if

*-----------------------------------------------------------------------
*     make photons
*        in prdgam, data is stored
*-----------------------------------------------------------------------

      if( ( jxs(12,iex) .ne. 0 .and. npikmt .eq. 0 ) .or.
     &     totgp1 .ne. 0. )then
       if( gwt(icl) .ne. -1.e6 ) then

        call prdgam(ireg,imat,ffac)

       end if
      endif

*-----------------------------------------------------------------------
*     make charged particles
*-----------------------------------------------------------------------

         if(nxs(7,iex).ne.0) then

                  call acecp(ireg,imat)

         end if

*-----------------------------------------------------------------------
*     capture
*-----------------------------------------------------------------------
*     sample analog capture if required.
*                    rtc(5,iex) : total cross section
*        nter = 12 ; rtc(3,iex) : capture
*        nter = 14 ; rtc(8,iex) : capture to fission
*-----------------------------------------------------------------------

            if(erg.gt.emcf(1))go to 120

         if( icels .ne. 0 ) then

            t1=rtc(3,iex)+rtc(8,iex)
            r=rang()*rtc(5,iex)

         else

            t1=sigca
            r=rang()*sigal

         end if

            if(t1.lt.r) goto 160

            i=3
            if(r.ge.rtc(3,iex))i=4
            nter=2*i+6

            wgt = 0.0d0
            goto 125

*-----------------------------------------------------------------------
*     otherwise simulate capture by weight reduction.
*     implicit capture
*-----------------------------------------------------------------------

  120       continue

         if( icels .ne. 0 ) then

            t8=wgt/rtc(5,iex)
            t1=t8*(rtc(3,iex)+rtc(8,iex))

         else

            t8=wgt/sigal
            t1=t8*sigca

         end if

            wgt=wgt-t1
            wgt = max( 0.0d0, wgt )

*-----------------------------------------------------------------------
*        summary of capture
*        wga is the capture loss probability
*        wgf is the weight after absorption
*-----------------------------------------------------------------------

  125 continue

                     wga = wgf - wgt
                     aevts(iaevt+15,ireg) =
     &               aevts(iaevt+15,ireg) + wga
                     bevts(ibevt+15,imat) =
     &               bevts(ibevt+15,imat) + wga

                     aevts(iaevt+16,ireg) =
     &               aevts(iaevt+16,ireg) + wga * eein
                     bevts(ibevt+16,imat) =
     &               bevts(ibevt+16,imat) + wga * eein
                     wgf = wgt

*-----------------------------------------------------------------------
*              analog capture
*-----------------------------------------------------------------------

               if( wgf .le. 0.0d0 ) then

                     kcoll = 6

                     goto 900

               end if

*-----------------------------------------------------------------------
*     make neutrons
*-----------------------------------------------------------------------
*     s(a,b) case :
*        in xstsab, data is stored
*-----------------------------------------------------------------------

  160    if( iet .ne. 0 ) then

               call xstsab(ireg,imat)

               goto 900

         end if

*-----------------------------------------------------------------------
*     sample the velocity of the target nucleus. (x-6:jsh:rep-92-163)
*-----------------------------------------------------------------------

               ssr=0.

               rtc(2,iex)=rtc(1,iex)
               ktc(2,iex)=ktc(1,iex)

         if( icels .ne. 0 ) then

            if(awn(iex).le.1.5.or.erg.le.400.*ttn)
     &         call samvel(erg,rtc(2,iex),ar,vr,ktc(2,iex))

            if( ssr .ne. 0. ) then

               s=1./ssr
               vr(1)=s*vr(1)
               vr(2)=s*vr(2)
               vr(3)=s*vr(3)

            end if

         end if

*-----------------------------------------------------------------------
*     sample elastic or inelastic scattering.
*-----------------------------------------------------------------------

            call xstcol(0,ac,ireg,elrt)

               if( ntyn .eq. -99 ) iels = 1
               if( ntyn .eq.  19 ) ifis = 1

               if(kdb.ne.0) goto 900
               if(nter.ne.0.and.nter.ne.12) goto 900

               ni = cmult

*-----------------------------------------------------------------------
*           induced fission option : number of neutrons
*-----------------------------------------------------------------------

            if( ifis .eq.1 .and. iidfs .eq. 1 ) then

                  tnbar = acenu(jxs(2,iex))
                  nj    = 0

               if( izaid .eq. 94239 .or. izaid .eq. 94241 ) then

                  ni = nSmpNuDistDataPu239_241_MC(tnbar)
                  nj = 1

               else if( izaid .eq. 92232 .or. izaid .eq. 92234 .or.
     &                  izaid .eq. 92236 .or. izaid .eq. 92238 ) then

                  ni = nSmpNuDistDataU232_234_236_238_MC(tnbar)
                  nj = 1

               else if( izaid .eq. 92233 .or. izaid .eq. 92235 ) then

                  ni = nSmpNuDistDataU233_235_MC(tnbar)
                  nj = 1

               end if

                  if( ni .le. 0 ) goto 900

               if( nj .eq. 1 ) then

                     ridfs = 1.0d0

                  do k = 1, ni

                     colout(1,k) = SmpWatt(eein, izaid)
                     colout(2,k) = 2.*rang()-1.
                     colout(3,k) = 0.0

                  end do

               end if

            end if

*-----------------------------------------------------------------------
*           nter=12 case: capture
*-----------------------------------------------------------------------

            if( nter .eq. 12 ) then

                  wgt = 0.0d0
                  goto 125

            end if

*-----------------------------------------------------------------------

               if( iels .eq. 1 ) then

                     kcoll = 3
                     aevts(iaevt+17,ireg) =
     &               aevts(iaevt+17,ireg) + wgt
                     bevts(ibevt+17,imat) =
     &               bevts(ibevt+17,imat) + wgt
*-----------------------------------------------------------------------

               else if( ifis .eq. 1 ) then

                     kcoll = 5
                     aevts(iaevt+19,ireg) =
     &               aevts(iaevt+19,ireg) + wgt
                     bevts(ibevt+19,imat) =
     &               bevts(ibevt+19,imat) + wgt
                     rumpal(2) = rumpal(2) - wgt

*-----------------------------------------------------------------------
cKN 2014/12/02; fission turn off, treated as capture

                  if( itfxs .ne. 0 ) then

                     kcoll = 6
                     nter = 12
                     wgt = 0.0d0
                     ifis = 0
                     goto 125

                  end if

*-----------------------------------------------------------------------

               else if( ni .gt. 0 ) then

                     kcoll = 4
                     aevts(iaevt+18,ireg) =
     &               aevts(iaevt+18,ireg) + wgt
                     bevts(ibevt+18,imat) =
     &               bevts(ibevt+18,imat) + wgt

                     aevts(iaevt+20,ireg) =
     &               aevts(iaevt+20,ireg) + wgt * ni
                     bevts(ibevt+20,imat) =
     &               bevts(ibevt+20,imat) + wgt * ni
                     rumpal(2) = rumpal(2) + wgt * ( ni - 1 )

               end if

*-----------------------------------------------------------------------

               if( ni .gt. 1 ) then

                     numpal(2) = numpal(2) + ( ni - 1 )

               end if

*-----------------------------------------------------------------------
*     bank all the neutrons from the collision.
*     number of neutrons
*-----------------------------------------------------------------------

      if( ni .gt. 0 ) then

         do k = 1, ni

*-----------------------------------------------------------------------
*           prepare to bank non-thermal neutrons.
*-----------------------------------------------------------------------

               if( ssr .eq. 0. ) then

                  call dtcos(colout(2,k),uold,uuu,0,irdm)
                  erg = colout(1,k)

                  ipsc = 5

*-----------------------------------------------------------------------
*           prepare to bank thermal free gas collision neutrons.
*-----------------------------------------------------------------------

               else

                  call dtcos(colout(2,k),vr,uuu,0,irdm)

                  s=sqrt(colout(1,k)*ar)
                  uuu=s*uuu+vtr(1)
                  vvv=s*vvv+vtr(2)
                  www=s*www+vtr(3)
                  t=uuu**2+vvv**2+www**2

                  s=1./sqrt(t)
                  uuu=uuu*s
                  vvv=vvv*s
                  www=www*s
                  erg=t/ar

                  ipsc = 4

               end if

*-----------------------------------------------------------------------
*        prepare to bank delayed neutrons.  colout(3,k) = time delay.
*        if biased, colout(3,1) = -weight adjusted for bias.
*        tme is delayed time ( shakes is 10 * ns )
*-----------------------------------------------------------------------

               if( colout(3,k) .le. 0. ) then

                  tme = 0.0

               else

                  ifis = 2

                  tme = colout(3,k) * 10.
                  if(colout(3,1).lt.0.) wgt = -colout(3,1)

               end if

*-----------------------------------------------------------------------
*           bank all neutron
*-----------------------------------------------------------------------

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 2
                     nixcos = nixcos + 1  ! MATSUDA 2022.12.12  for T-Point

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 1
                     jclusts(3,nclsts) = 2
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = 0
                     jclusts(6,nclsts) = 1
                     jclusts(7,nclsts) = 2112
                     jclusts(8,nclsts) = 0

                     rms = gpt(1)
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
                     ixcoss(nixcos) = ixcos  ! MATSUDA 2022.12.12  for T-Point

               if( ifis .ne. 0 ) then
                     aevts(iaevt+21,ireg) =
     &               aevts(iaevt+21,ireg) + wgt
                     bevts(ibevt+21,imat) =
     &               bevts(ibevt+21,imat) + wgt
                     rumpal(2) = rumpal(2) + wgt

               end if

               if( ni .eq. 1 ) then

                  if( erg .gt. eg0 ) then
                     aevts(iaevt+22,ireg) =
     &               aevts(iaevt+22,ireg) + wgt * ( erg - eg0 )
                     bevts(ibevt+22,imat) =
     &               bevts(ibevt+22,imat) + wgt * ( erg - eg0 )
                  else
                     aevts(iaevt+23,ireg) =
     &               aevts(iaevt+23,ireg) + wgt * ( eg0 - erg )
                     bevts(ibevt+23,imat) =
     &               bevts(ibevt+23,imat) + wgt * ( eg0 - erg )
                  end if

                  if( iels .eq. 0 .and. ifis .eq. 0 ) then
                     aevts(iaevt+24,ireg) =
     &               aevts(iaevt+24,ireg) + wgt
                     bevts(ibevt+24,imat) =
     &               bevts(ibevt+24,imat) + wgt
                  end if

               end if

*-----------------------------------------------------------------------

               if( nter .ne. 0 ) nter = 0
               wgt = wgf

         end do

      end if

*-----------------------------------------------------------------------

  900 continue

            kdecay(4) = ifis

             uuu = usave(1)
             vvv = usave(2)
             www = usave(3)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sctneut(eein,wgti,ireg,imat,ffac1,ffac2,icmm,elrt)
*                                                                      *
*        calculate the collision of a neutron with a nucleus           *
*        in event generrator mode.                                     *
*        modified by K.Niita on   2005/11/01                           *
*        modified by Y.Iwamoto on 2007/01/04                           *
*        modified by K.Niita on   2008/12/10                           *
*        modified by K.Niita on   2009/06/02                           *
*        modified by K.Niita on   2009/08/21                           *
*        modified by K.Niita on   2010/01/27                           *
*        modified by T.Ogawa on   2014/03/23                           *
*        modified by K.Niita on   2016/08/08                           *
*                                                                      *
*        parameters:                                                   *
*           ge1 (1MeV)   : cutoff energy for evapolation model         *
*           ge2 (1MeV)   : cutoff excitation energy for gamma decay of *
*                          residual nuclei                             *
*                                                                      *
*        Note: switching by ge1 and ge2 is meaningful only when        *
*              igamma = 1, where gamma emission is sometimes incomplete*
*              therefore excitation energy is left after nevap.        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for REDUCTION_COUNTER
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2

      use GGMBANKMOD !FURUTA
      use NGSDATAMOD, only : bindeg
      use GGMARRAYMOD !2020ASTOM
      use moddas_material
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( rpmass = 938.27, rnmass = 939.58 )

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

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
C MATSUDA 2023.03.15  for T-Point
      common /clusttp/ nixcos, ixcoss(nnn)
!$OMP THREADPRIVATE(/clusttp/)

      common /emode/  emodem, ge1, ge2, iemode
      common /emode1/ nevhin, nlowlv, nefiss, ntwidt, mtprec,ipcnt(70)
!$OMP THREADPRIVATE(/emode1/)

      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /eparm/  esmax, esmin, emin(20)
      common /swich3/ ielst, jelst, kelst
!$OMP THREADPRIVATE(/swich3/)
      common /ismtfl/ ismtflg
!$OMP THREADPRIVATE(/ismtfl/)
      common /tcntl/  icntl, inucr
      common /elreg/  elarg(kvlmax,4), mlrgn, melrg(kvlmax,2), ielusr
      common /regdc/  idrg(kvlmax), idgr(kvmmax)



      common /kmat1g/ kmat(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)
      common /usaves/  usave(3)
!$OMP THREADPRIVATE(/usaves/)
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma

*-----------------------------------------------------------------------

      dimension vr(3)
      dimension dirc(3)

      data nevhin / 0 /
      data ntwidt / 0 /
      data nlowlv / 0 /
      data nefiss / 0 /
      data mtprec / 0 /

      mtprec = 0

*-----------------------------------------------------------------------

           nclst  = 0
           nclsts = 0
           nixcos = 0  ! MATSUDA 2023.03.15  for T-Point

*-----------------------------------------------------------------------

            kelst = 0
            kcoll = 0

            nter = 0
            ipt  = 1
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

            ffac1 = 0.d0
            ffac2 = 0.d0

*-----------------------------------------------------------------------
*     save the incoming direction.
*-----------------------------------------------------------------------

            uold(1) = 0.0
            uold(2) = 0.0
            uold(3) = 1.0

            ntyn = 0
            mk   = mat(icl)
            m    = jmd(1+mk)
            m1   = m

*-----------------------------------------------------------------------
*     sample the nuclide, using the cumulative total cross section.
*-----------------------------------------------------------------------

         if( npq(mk) .eq. 1 ) goto 20

*-----------------------------------------------------------------------

         if( icmm .ne. 0 ) then

               m = icmm

         else

               c = rang() * totm

            do m = m1, jmd(1+mk+1) - 2

               c = c - rtc(5,lme(1,m)) * fme(m)

               if( c .lt. 0. ) goto 20

            end do

        end if

*-----------------------------------------------------------------------

   20       iex  = lme(1,m)
            iexp = m
            mtp  = 2
            mpan = ipan(icl) + m - m1

*-----------------------------------------------------------------------
*           iet = 0 ; S(a,b)
*-----------------------------------------------------------------------

            iet  = lmt(m)
               if( iet.ne.0 ) then
                 if( erg .gt. esa(iet) )   iet = 0
               end if
               if( rtc(11,iex) .ge. 0. ) iet = 0

*-----------------------------------------------------------------------
*        only for inucr=13 : gamma production by nuclear data
*-----------------------------------------------------------------------

         if( icntl .eq. 1 .and. inucr .eq. 13 ) then

                     nclst  = 0
                     nclsts = 0

            if(   gwt(icl) .ne. -1.e6 .and.
     &        ( ( jxs(12,iex) .ne. 0 .and. npikmt .eq. 0 ) .or.
     &            totgp1 .ne. 0. ) ) then

                     ffag = 1.0

                  call prdgam(ireg,imat,ffag)

            end if

                  sek = 0.0d0

               do k = 1, nclsts

                  if( jclusts(3,k) .eq. 14 ) then

                     sek = sek + qclusts(7,k) * qclusts(8,k)

                  end if

               end do

                  ffac2 = sek

                  nclst  = 0
                  nclsts = 0

         end if

*-----------------------------------------------------------------------
*        target mass (MeV)
*        expand natural nucleus and choose one nucleus
*-----------------------------------------------------------------------

                     izm = iza(iexp) / 1000
                     ims = iza(iexp) - 1000 * izm

                  if( ims .eq. 0 ) then

                     lemm  = nint( dnel_das(kmat0+imat) )

                        sek = 0.0
                     do lem = 1, lemm
                            izt = zz_das(kmat(imat)+lem)
                        if( izt .eq. izm )
     &                     sek = sek + den_das(kmat(imat)+lem)
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

   55                   continue

                        ims = nint( a_das(kmat(imat)+lem) )

                  end if

                     inm = ims - izm

                     be = bindeg(izm,inm)

                     tams = ( izm * rpmass + inm * rnmass - be )
                     rnms = rnmass

                     rpm = rnms / ( tams + rnms )
                     rtm = tams / ( tams + rnms )

                     pin0 = sqrt( 2. * rnms * erg )

                     pipc = pin0 - rpm * pin0
                     pitc = - rtm * pin0

                     egpc = pipc**2 / 2. / rnms
                     egtc = pitc**2 / 2. / tams

                     if( izm .eq. 1 .and. inm .eq. 0 ) then
                        ipad = 1
                     else if( izm .eq. 1 .and. inm .eq. 1 ) then
                        ipad = 15
                     else if( izm .eq. 1 .and. inm .eq. 2 ) then
                        ipad = 16
                     else if( izm .eq. 2 .and. inm .eq. 1 ) then
                        ipad = 17
                     else if( izm .eq. 2 .and. inm .eq. 2 ) then
                        ipad = 18
                     else
                        ipad = 19
                     end if

*-----------------------------------------------------------------------
*           mother nucleus
*-----------------------------------------------------------------------

                  mathz = izm
                  mathn = inm

*-----------------------------------------------------------------------
*     user elastic and inelastic option if exist
*-----------------------------------------------------------------------

            ioel = 0
            icels = 1

      if( mlrgn .gt. 0 ) then

            do jj = 1, mlrgn

              do kk = melrg(jj,1), melrg(jj,2)

                 if( ireg .eq. idgr(kk) ) then

                    ioel = jj
                    goto 300

                 end if

              end do

            end do

  300       continue

      end if

*-----------------------------------------------------------------------
*     capture cross sections by user defined function
*-----------------------------------------------------------------------

      if( ioel .ne. 0 ) then

         if( ielusr .eq. 1 ) then

            call usrelst1(3,icels,ioel,erg,sigca,sigal)

         else

            call usrelst2(3,icels,ioel,erg,sigca,sigal)

         end if

      end if

*-----------------------------------------------------------------------
*     capture
*-----------------------------------------------------------------------
*     only analog capture is considered.
*                    rtc(5,iex) : total cross section
*        nter = 12 ; rtc(3,iex) : capture
*        nter = 14 ; rtc(8,iex) : capture to fission
*-----------------------------------------------------------------------

         if( icels .ne. 0 ) then

            t1 = rtc(3,iex) + rtc(8,iex)
            r = rang() * rtc(5,iex)

         else

            t1 = sigca
            r = rang() * sigal

         end if

            if( t1 .lt. r ) goto 160

*-----------------------------------------------------------------------

            i = 3
            if( icels .ne. 0 .and. r .ge. rtc(3,iex) ) i = 4
            nter = 2 * i + 6

*-----------------------------------------------------------------------
*     fission or capture
*-----------------------------------------------------------------------
            if( nter .eq. 12 ) goto 125
            if( nter .eq. 14 ) then
               ni = 0
               goto 126
            end if

*-----------------------------------------------------------------------
*     fission case
*-----------------------------------------------------------------------
	
  126 continue
*-----------------------------------------------------------------------

                     tapm = egpc + egtc + rnms + tams

                     be = bindeg(izm,inm+1)
                     tapl = ( izm * rpmass + ( inm + 1 ) * rnmass - be )

*-----------------------------------------------------------------------

                     plmx = 0.0
                     plmy = 0.0
                     plmz = tapm / ( tams + rnms ) * pin0

                     pabs = plmx**2 + plmy**2 + plmz**2
                     etot = sqrt( pabs + tapm**2 )
                     elab = etot - tapm

                     nclst = nclst + 1
                     iclust(nclst) = 0

                     jclust(0,nclst) = 0
                     jclust(1,nclst) = izm
                     jclust(2,nclst) = inm + 1
                     jclust(3,nclst) = 19
                     jclust(4,nclst) = 0
                     jclust(5,nclst) = izm
                     jclust(6,nclst) = izm + inm + 1
                     jclust(7,nclst) = izm * 1000000 + izm + inm + 1
                     jclust(8,nclst) = 0

                     qclust(0,nclst)  = 0.0
                     qclust(1,nclst)  = plmx / 1000.
                     qclust(2,nclst)  = plmy / 1000.
                     qclust(3,nclst)  = plmz / 1000.
                     qclust(4,nclst)  = etot / 1000.
                     qclust(5,nclst)  = tapm / 1000.
                     qclust(6,nclst)  = tapm - tapl
                     qclust(7,nclst)  = elab
                     qclust(8,nclst)  = wgt / wg0
                     qclust(9,nclst)  = tme
                     qclust(10,nclst) = 0.0d0
                     qclust(11,nclst) = 0.0d0
                     qclust(12,nclst) = 0.0d0

*-----------------------------------------------------------------------

                           nevhin = 0
                           nlowlv = 1
                           nefiss = 1

                     call nevap(2)

                           nevhin = 0
                           nlowlv = 0
                           nefiss = 0

*-----------------------------------------------------------------------

               do k = 1, nclsts

*-----------------------------------------------------------------------
*                 nucleus and exc. energy > 1.0 MeV
*                 decay to gamma
*-----------------------------------------------------------------------

                  if( iclusts(k) .eq. 0 .and.
     &              qclusts(6,k) .gt. ge2 .and. abs(igamma) .ne. 3) then ! prohibit deexcitation of isomers

                     pxc = qclusts(1,k) * 1000.
                     pyc = qclusts(2,k) * 1000.
                     pzc = qclusts(3,k) * 1000.
                     exc = qclusts(6,k)
                     enk = qclusts(7,k)

                     tamr = qclusts(5,k) * 1000.
                     tamp = tamr + exc
                     pal  = sqrt( pxc**2 + pyc**2 + pzc**2 )

                     css = 2.d0 * rang() - 1.0d0
                     pgam = ( tamp**2 - tamr**2 ) / 2.
     &                    / ( enk + tamp - pal * css )

                     dirc(1) = pxc / pal
                     dirc(2) = pyc / pal
                     dirc(3) = pzc / pal

                     call dtcos(css,dirc,uuu,0,irdm)

                     pgax = uuu * pgam
                     pgay = vvv * pgam
                     pgaz = www * pgam

                     plmx = pxc - pgax
                     plmy = pyc - pgay
                     plmz = pzc - pgaz

                     pabs = plmx**2 + plmy**2 + plmz**2
                     etot = sqrt( pabs + tamr**2 )
                     elab = etot - tamr

                     qclusts(1,k)  = plmx / 1000.
                     qclusts(2,k)  = plmy / 1000.
                     qclusts(3,k)  = plmz / 1000.
                     qclusts(4,k)  = etot / 1000.
                     qclusts(5,k)  = tamr / 1000.
                     qclusts(6,k)  = 0.0
                     qclusts(7,k)  = elab

                     ippad = jclusts(3,k)

                     numpal(ippad) = numpal(ippad) + 1
                     rumpal(ippad) = rumpal(ippad) + wgt

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

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = pgax / 1000.
                     qclusts(2,nclsts)  = pgay / 1000.
                     qclusts(3,nclsts)  = pgaz / 1000.
                     qclusts(4,nclsts)  = pgam / 1000.
                     qclusts(5,nclsts)  = 0.0
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = pgam
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = 0.0
                     qclusts(11,nclsts) = 0.0
                     qclusts(12,nclsts) = 0.0

                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

*-----------------------------------------------------------------------

                  else

                     ippad = jclusts(3,k)

                     numpal(ippad) = numpal(ippad) + 1
                     rumpal(ippad) = rumpal(ippad) + wgt

                  end if

               end do

*-----------------------------------------------------------------------
*           fission case
*-----------------------------------------------------------------------

                     kcoll = 5

                     aevts(iaevt+19,ireg) =
     &               aevts(iaevt+19,ireg) + wgt
                     bevts(ibevt+19,imat) =
     &               bevts(ibevt+19,imat) + wgt

                  goto 900

*-----------------------------------------------------------------------
*        summary of capture
*-----------------------------------------------------------------------

  125 continue
                     aevts(iaevt+15,ireg) =
     &               aevts(iaevt+15,ireg) + wgt
                     bevts(ibevt+15,imat) =
     &               bevts(ibevt+15,imat) + wgt

                     aevts(iaevt+16,ireg) =
     &               aevts(iaevt+16,ireg) + wgt * eein
                     bevts(ibevt+16,imat) =
     &               bevts(ibevt+16,imat) + wgt * eein
*-----------------------------------------------------------------------

                     tapm = egpc + egtc + rnms + tams

                     be = bindeg(izm,inm+1)
                     tapl = ( izm * rpmass + ( inm + 1 ) * rnmass - be )

*-----------------------------------------------------------------------
*           low excitation, only gamma decay
*-----------------------------------------------------------------------

            if( tapm - tapl .lt. ge1 .and. abs(igamma) .ne. 3 .and.
     &       iemode .eq. 1 .and. eein .le. emodem ) then

                     ecmp = tapm - sqrt( tapl * ( 2. * tapm - tapl ) )
                     pgam = sqrt( 2.* tapl * ecmp )

                     css = 2.d0 * rang() - 1.0d0

                     call dtcos(css,uold,uuu,0,irdm)

                     pgax = uuu * pgam
                     pgay = vvv * pgam
                     pgaz = www * pgam

                     plmx = -pgax
                     plmy = -pgay
                     plmz = -pgaz + tapl / ( tams + rnms ) * pin0

                     pabs = plmx**2 + plmy**2 + plmz**2
                     etot = sqrt( pabs + tapl**2 )
                     elab = etot - tapl

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 0

                     jclusts(0,nclsts) = 0
                     jclusts(1,nclsts) = izm
                     jclusts(2,nclsts) = inm + 1
                     jclusts(3,nclsts) = 19
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = izm
                     jclusts(6,nclsts) = izm + inm + 1
                     jclusts(7,nclsts) = izm * 1000000 + izm + inm + 1
                     jclusts(8,nclsts) = 0

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = plmx / 1000.
                     qclusts(2,nclsts)  = plmy / 1000.
                     qclusts(3,nclsts)  = plmz / 1000.
                     qclusts(4,nclsts)  = etot / 1000.
                     qclusts(5,nclsts)  = tapl / 1000.
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = elab
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = tme
                     qclusts(10,nclsts) = 0.0d0
                     qclusts(11,nclsts) = 0.0d0
                     qclusts(12,nclsts) = 0.0d0

                     numpal(19) = numpal(19) + 1
                     rumpal(19) = rumpal(19) + wgt

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

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = pgax / 1000.
                     qclusts(2,nclsts)  = pgay / 1000.
                     qclusts(3,nclsts)  = pgaz / 1000.
                     qclusts(4,nclsts)  = pgam / 1000.
                     qclusts(5,nclsts)  = 0.0
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = pgam
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = 0.0
                     qclusts(11,nclsts) = 0.0
                     qclusts(12,nclsts) = 0.0

                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

*-----------------------------------------------------------------------
*           high excitation, gamma and/or charge particle decay
*-----------------------------------------------------------------------

            else

                     plmx = 0.0
                     plmy = 0.0
                     plmz = tapm / ( tams + rnms ) * pin0

                     pabs = plmx**2 + plmy**2 + plmz**2
                     etot = sqrt( pabs + tapm**2 )
                     elab = etot - tapm

                     nclst = nclst + 1
                     iclust(nclst) = 0

                     jclust(0,nclst) = 0
                     jclust(1,nclst) = izm
                     jclust(2,nclst) = inm + 1
                     jclust(3,nclst) = 19
                     jclust(4,nclst) = 0
                     jclust(5,nclst) = izm
                     jclust(6,nclst) = izm + inm + 1
                     jclust(7,nclst) = izm * 1000000 + izm + inm + 1
                     jclust(8,nclst) = 0

                     qclust(0,nclst)  = 0.0
                     qclust(1,nclst)  = plmx / 1000.
                     qclust(2,nclst)  = plmy / 1000.
                     qclust(3,nclst)  = plmz / 1000.
                     qclust(4,nclst)  = etot / 1000.
                     qclust(5,nclst)  = tapm / 1000.
                     qclust(6,nclst)  = tapm - tapl
                     qclust(7,nclst)  = elab
                     qclust(8,nclst)  = wgt / wg0
                     qclust(9,nclst)  = tme
                     qclust(10,nclst) = 0.0d0
                     qclust(11,nclst) = 0.0d0
                     qclust(12,nclst) = 0.0d0

*-----------------------------------------------------------------------

                           nevhin = 1
                           nlowlv = 1
                           nefiss = 0

                     if(iemode .ge. 2) call xsielsch

                     call nevap(2)

                           nevhin = 0
                           nlowlv = 0
                           nefiss = 0

*-----------------------------------------------------------------------

               do k = 1, nclsts

*-----------------------------------------------------------------------
*                 nucleus and exc. energy > 1.0 MeV
*                 decay to gamma
*-----------------------------------------------------------------------

                  if( iclusts(k) .eq. 0 .and.
     &              qclusts(6,k) .gt. ge2 .and. abs(igamma) .ne. 3) then ! prohibit deexcitation of isomers

                     pxc = qclusts(1,k) * 1000.
                     pyc = qclusts(2,k) * 1000.
                     pzc = qclusts(3,k) * 1000.
                     exc = qclusts(6,k)
                     enk = qclusts(7,k)

                     tamr = qclusts(5,k) * 1000.
                     tamp = tamr + exc
                     pal  = sqrt( pxc**2 + pyc**2 + pzc**2 )

                     css = 2.d0 * rang() - 1.0d0
                     pgam = ( tamp**2 - tamr**2 ) / 2.
     &                    / ( enk + tamp - pal * css )

                     dirc(1) = pxc / pal
                     dirc(2) = pyc / pal
                     dirc(3) = pzc / pal

                     call dtcos(css,dirc,uuu,0,irdm)

                     pgax = uuu * pgam
                     pgay = vvv * pgam
                     pgaz = www * pgam

                     plmx = pxc - pgax
                     plmy = pyc - pgay
                     plmz = pzc - pgaz

                     pabs = plmx**2 + plmy**2 + plmz**2
                     etot = sqrt( pabs + tamr**2 )
                     elab = etot - tamr

                     qclusts(1,k)  = plmx / 1000.
                     qclusts(2,k)  = plmy / 1000.
                     qclusts(3,k)  = plmz / 1000.
                     qclusts(4,k)  = etot / 1000.
                     qclusts(5,k)  = tamr / 1000.
                     qclusts(6,k)  = 0.0
                     qclusts(7,k)  = elab

                     ippad = jclusts(3,k)

                     numpal(ippad) = numpal(ippad) + 1
                     rumpal(ippad) = rumpal(ippad) + wgt

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

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = pgax / 1000.
                     qclusts(2,nclsts)  = pgay / 1000.
                     qclusts(3,nclsts)  = pgaz / 1000.
                     qclusts(4,nclsts)  = pgam / 1000.
                     qclusts(5,nclsts)  = 0.0
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = pgam
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = 0.0
                     qclusts(11,nclsts) = 0.0
                     qclusts(12,nclsts) = 0.0

                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

*-----------------------------------------------------------------------

                  else

                     ippad = jclusts(3,k)

                     numpal(ippad) = numpal(ippad) + 1
                     rumpal(ippad) = rumpal(ippad) + wgt

                  end if

               end do

*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------
*           analog capture
*-----------------------------------------------------------------------

                  kcoll = 6

                  goto 900

*-----------------------------------------------------------------------
*     make neutrons
*-----------------------------------------------------------------------
*     s(a,b) case :
*        in xstsab, data is stored
*-----------------------------------------------------------------------

  160    if( iet .ne. 0 ) then

               call xstsab(ireg,imat)

               iels = 1

               goto 910

         end if

*-----------------------------------------------------------------------
*     sample the velocity of the target nucleus. (x-6:jsh:rep-92-163)
*-----------------------------------------------------------------------

               ssr=0.

               rtc(2,iex)=rtc(1,iex)
               ktc(2,iex)=ktc(1,iex)

         if( icels .ne. 0 ) then

            if(awn(iex).le.1.5.or.erg.le.400.*ttn)
     &         call samvel(erg,rtc(2,iex),ar,vr,ktc(2,iex))

            if( ssr .ne. 0. ) then

               s=1./ssr
               vr(1)=s*vr(1)
               vr(2)=s*vr(2)
               vr(3)=s*vr(3)

            end if

         end if

*-----------------------------------------------------------------------
*     sample elastic or inelastic scattering.
*-----------------------------------------------------------------------

            ismtflg = 0
            call xstcol(0,ac,ireg,elrt)
               if( ismtflg .eq. 1 ) goto 900 ! advanced event generator mode

               if( ntyn .eq. -99 ) iels = 1 ! elastic
               if( ntyn .eq.  19 ) ifis = 1 ! fission

               if( iels .eq. 1 ) kelst = 1

               if( kdb .ne. 0 ) goto 900
               if( nter .ne. 0 .and. nter .ne. 12 ) goto 900

               ni = cmult ! # of neutrons

*-----------------------------------------------------------------------
*           nter=12 case: capture
*           ifis = 1  case: fission
*-----------------------------------------------------------------------

            if( ifis .eq. 1 ) then
             if(itfxs.eq.1) then  ! T.Sato 2024/09/22 consider nonu option (nonu = 0 -> itfxs = 1 and no fission)
              goto 125  ! simple capture
             else
              goto 126  ! cause fission
             endif
            endif
            if( nter .eq. 12 ) goto 125

*-----------------------------------------------------------------------
*     bank all the neutrons from the collision.
*     number of neutrons
*-----------------------------------------------------------------------

      if( ni .le. 0 ) goto 900

*-----------------------------------------------------------------------

                  j = 1

! here the larger energy is taken as the first neutron energy while the smaller energy is discarded.
! this process is justified because subsequent neutrons sampled by evaporation model is low-energy
             if( ni .eq. 2 ) then
                   enim = 0.0d+0
                do k = 1, ni
                   if( colout(1,k) .gt. enim ) then
                      enim = colout(1,k)
                      j = k
                   end if
                end do
             else if( ni .ge. 3 ) then ! here the smallest energy is taken as the first neutron energy while the others are discarded.
                   enim = 1.0d+10
                do k = 1, ni
                   if( colout(1,k) .lt. enim ) then
                      enim = colout(1,k)
                      j = k
                   end if
                end do
             end if

*-----------------------------------------------------------------------
*           prepare to bank non-thermal neutrons.
*-----------------------------------------------------------------------

               if( ssr .eq. 0. ) then

                  call dtcos(colout(2,j),uold,uuu,0,irdm)
                  erg = colout(1,j)

*-----------------------------------------------------------------------
*           prepare to bank thermal free gas collision neutrons.
*-----------------------------------------------------------------------

               else

                  call dtcos(colout(2,j),vr,uuu,0,irdm)

                  s=sqrt(colout(1,j)*ar)
                  uuu=s*uuu+vtr(1)
                  vvv=s*vvv+vtr(2)
                  www=s*www+vtr(3)
                  t=uuu**2+vvv**2+www**2

                  s=1./sqrt(t)
                  uuu=uuu*s
                  vvv=vvv*s
                  www=www*s
                  erg=t/ar

               end if

*-----------------------------------------------------------------------
*           bank first neutron
*-----------------------------------------------------------------------

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 2

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 1
                     jclusts(3,nclsts) = 2
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = 0
                     jclusts(6,nclsts) = 1
                     jclusts(7,nclsts) = 2112
                     jclusts(8,nclsts) = 0

                     rms = rnmass
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

                     numpal(2) = numpal(2) + 1
                     rumpal(2) = rumpal(2) + wgt

*-----------------------------------------------------------------------
*              elastic : booking proton if target is proton
*-----------------------------------------------------------------------

               if( iels .eq. 1 ) then

                     kcoll = 3
                     aevts(iaevt+17,ireg) =
     &               aevts(iaevt+17,ireg) + wgt
                     bevts(ibevt+17,imat) =
     &               bevts(ibevt+17,imat) + wgt
*-----------------------------------------------------------------------

               else if( ni .gt. 0 ) then

                     kcoll = 4
                     aevts(iaevt+18,ireg) =
     &               aevts(iaevt+18,ireg) + wgt
                     bevts(ibevt+18,imat) =
     &               bevts(ibevt+18,imat) + wgt

                     aevts(iaevt+20,ireg) =
     &               aevts(iaevt+20,ireg) + wgt * ni
                     bevts(ibevt+20,imat) =
     &               bevts(ibevt+20,imat) + wgt * ni
               end if

*-----------------------------------------------------------------------

  910 continue

*-----------------------------------------------------------------------
*        elastic
*-----------------------------------------------------------------------

         if( iels .eq. 1 ) then

*-----------------------------------------------------------------------
*              elastic : booking of target, proton and the others
*-----------------------------------------------------------------------

                     pxl = qclusts(1,nclsts) * 1000.
                     pyl = qclusts(2,nclsts) * 1000.
                     pzl = qclusts(3,nclsts) * 1000.

                     pxc = - pxl
                     pyc = - pyl
                     pzc = - pzl + pin0

*-----------------------------------------------------------------------

                     pabs = pxc**2 + pyc**2 + pzc**2
                     etot = sqrt( pabs + tams**2 )
                     elab = etot - tams

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 0

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = izm
                     jclusts(2,nclsts) = inm
                     jclusts(3,nclsts) = ipad
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = izm
                     jclusts(6,nclsts) = izm + inm
                     jclusts(8,nclsts) = 0

                  if( ipad .eq. 1 ) then
                     jclusts(7,nclsts) = 2212
                  else
                     jclusts(7,nclsts) = izm * 1000000 + izm + inm
                  end if

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = pxc / 1000.
                     qclusts(2,nclsts)  = pyc / 1000.
                     qclusts(3,nclsts)  = pzc / 1000.
                     qclusts(4,nclsts)  = etot / 1000.
                     qclusts(5,nclsts)  = tams / 1000.
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = elab
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = tme
                     qclusts(10,nclsts) = 0.0d0
                     qclusts(11,nclsts) = 0.0d0
                     qclusts(12,nclsts) = 0.0d0

                     numpal(ipad) = numpal(ipad) + 1
                     rumpal(ipad) = rumpal(ipad) + wgt

                     goto 900

*-----------------------------------------------------------------------
*        inelastic
*-----------------------------------------------------------------------

         else   ! case inelastic

                     pxl = qclusts(1,nclsts) * 1000. ! x momentum. conversion from GeV/c to MeV/c
                     pyl = qclusts(2,nclsts) * 1000.
                     pzl = qclusts(3,nclsts) * 1000.

                     pxc = - pxl
                     pyc = - pyl
                     pzc = - pzl + rpm * pin0 ! pin0: sqrt(2 mass E_kin) ref: L1707,

                     pal = sqrt( pxc**2 + pyc**2 + pzc**2 )
                     egpp = pal**2 / 2. / rnms ! classical kinetic energy

                     tapm = sqrt( ( egpc - egpp + egtc + tams )**2
     &                              - pal**2 ) !

*-----------------------------------------------------------------------
*              low excitation inelastic, only gamma decay
*-----------------------------------------------------------------------

! from here to "goto 900", one gamma emission and target recoil is determined.
            if( tapm - tams .le. ge1 .and. abs(igamma) .eq. 1) then  ! Energy balance, gives excitation energy

                     egtp = pal**2 / 2. / tapm

                     css = 2.d0 * rang() - 1.0d0
                     pgam = ( tapm**2 - tams**2 ) / 2.
     &                    / ( egtp + tapm - pal * css )  ! gamma ray momentum (only one gamma is considered)

                     dirc(1) = pxc / pal
                     dirc(2) = pyc / pal
                     dirc(3) = pzc / pal

                     call dtcos(css,dirc,uuu,0,irdm)

                     pgax = uuu * pgam  ! gamma momentum
                     pgay = vvv * pgam
                     pgaz = www * pgam

                     plmx = pxc - pgax  ! product momentum
                     plmy = pyc - pgay
                     plmz = pzc - pgaz + rtm * pin0

                     pabs = plmx**2 + plmy**2 + plmz**2
                     etot = sqrt( pabs + tams**2 )  ! system total energy
                     elab = etot - tams ! emission particle energy

                     nclsts = nclsts + 1  ! produce one particle
                     iclusts(nclsts) = 0

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = izm
                     jclusts(2,nclsts) = inm
                     jclusts(3,nclsts) = ipad
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = izm
                     jclusts(6,nclsts) = izm + inm
                     jclusts(7,nclsts) = izm * 1000000 + izm + inm
                     jclusts(8,nclsts) = 0

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = plmx / 1000.
                     qclusts(2,nclsts)  = plmy / 1000.
                     qclusts(3,nclsts)  = plmz / 1000.
                     qclusts(4,nclsts)  = etot / 1000.
                     qclusts(5,nclsts)  = tams / 1000.
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = elab
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = tme
                     qclusts(10,nclsts) = 0.0d0
                     qclusts(11,nclsts) = 0.0d0
                     qclusts(12,nclsts) = 0.0d0

                     numpal(ipad) = numpal(ipad) + 1
                     rumpal(ipad) = rumpal(ipad) + wgt

*-----------------------------------------------------------------------

                     nclsts = nclsts + 1  ! only one gamma is considered.
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

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = pgax / 1000.
                     qclusts(2,nclsts)  = pgay / 1000.
                     qclusts(3,nclsts)  = pgaz / 1000.
                     qclusts(4,nclsts)  = pgam / 1000.
                     qclusts(5,nclsts)  = 0.0
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = pgam
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = 0.0
                     qclusts(11,nclsts) = 0.0
                     qclusts(12,nclsts) = 0.0

                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

                     goto 900

*-----------------------------------------------------------------------
*           high excitation inelastic, gamma and/or charge paticle decay
*-----------------------------------------------------------------------

            else  !  (n,Xn) reaction

                     plmx = pxc
                     plmy = pyc
                     plmz = pzc + tapm / ( tams + rnms ) * pin0

                     pabs = plmx**2 + plmy**2 + plmz**2
                     etot = sqrt( pabs + tapm**2 )
                     elab = etot - tapm

                     nclst = nclst + 1
                     iclust(nclst) = 0

                     jclust(0,nclst) = 0
                     jclust(1,nclst) = izm
                     jclust(2,nclst) = inm
                     jclust(3,nclst) = ipad
                     jclust(4,nclst) = 0
                     jclust(5,nclst) = izm
                     jclust(6,nclst) = izm + inm
                     jclust(7,nclst) = izm * 1000000 + izm + inm
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

*-----------------------------------------------------------------------

                           nevhin = ni
                           nlowlv = 1
                           nefiss = ifis

                     call nevap(2)

                           nevhin = 0
                           nlowlv = 0
                           nefiss = 0

*-----------------------------------------------------------------------

               do k = 1, nclsts

                  if( iclusts(k) .eq. 0 .and.
     &              qclusts(6,k) .gt. ge2 .and. abs(igamma) .ne. 3) then ! prohibit deexcitation of isomers

                     pxc = qclusts(1,k) * 1000.
                     pyc = qclusts(2,k) * 1000.
                     pzc = qclusts(3,k) * 1000.
                     exc = qclusts(6,k)
                     enk = qclusts(7,k)

                     tamr = qclusts(5,k) * 1000.
                     tamp = tamr + exc
                     pal  = sqrt( pxc**2 + pyc**2 + pzc**2 )

                     css = 2.d0 * rang() - 1.0d0
                     pgam = ( tamp**2 - tamr**2 ) / 2.
     &                    / ( enk + tamp - pal * css )

                     dirc(1) = pxc / pal
                     dirc(2) = pyc / pal
                     dirc(3) = pzc / pal

                     call dtcos(css,dirc,uuu,0,irdm)

                     pgax = uuu * pgam
                     pgay = vvv * pgam
                     pgaz = www * pgam

                     plmx = pxc - pgax
                     plmy = pyc - pgay
                     plmz = pzc - pgaz

                     pabs = plmx**2 + plmy**2 + plmz**2
                     etot = sqrt( pabs + tamr**2 )
                     elab = etot - tamr

                     qclusts(1,k)  = plmx / 1000.
                     qclusts(2,k)  = plmy / 1000.
                     qclusts(3,k)  = plmz / 1000.
                     qclusts(4,k)  = etot / 1000.
                     qclusts(5,k)  = tamr / 1000.
                     qclusts(6,k)  = 0.0
                     qclusts(7,k)  = elab

                     ippad = jclusts(3,k)

                     numpal(ippad) = numpal(ippad) + 1
                     rumpal(ippad) = rumpal(ippad) + wgt

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

                     qclusts(0,nclsts)  = 0.0
                     qclusts(1,nclsts)  = pgax / 1000.
                     qclusts(2,nclsts)  = pgay / 1000.
                     qclusts(3,nclsts)  = pgaz / 1000.
                     qclusts(4,nclsts)  = pgam / 1000.
                     qclusts(5,nclsts)  = 0.0
                     qclusts(6,nclsts)  = 0.0
                     qclusts(7,nclsts)  = pgam
                     qclusts(8,nclsts)  = wgt / wg0
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = 0.0
                     qclusts(11,nclsts) = 0.0
                     qclusts(12,nclsts) = 0.0

                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

*-----------------------------------------------------------------------

                  else

                     ippad = jclusts(3,k)

                     numpal(ippad) = numpal(ippad) + 1
                     rumpal(ippad) = rumpal(ippad) + wgt

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------

            if( icntl .eq. 1 .and.
     &        ( inucr .eq. 12 .or. inucr .eq. 13 ) ) then

                  sek = 0.0d0

               do k = 1, nclsts

                  if( inucr .eq. 12 .and. jclusts(5,k) .ne. 0 ) then

                     sek = sek + qclusts(7,k) * qclusts(8,k)

                  end if

                  if( inucr .eq. 13 .and. jclusts(3,k) .eq. 14 ) then

                     sek = sek + qclusts(7,k) * qclusts(8,k)

                  end if

               end do

                  ffac1 = sek

            end if

*-----------------------------------------------------------------------

             uuu = usave(1)
             vvv = usave(2)
             www = usave(3)

c-----------------------------------------------------------------------


      return
      end


************************************************************************
*                                                                      *
      subroutine prdgam(ireg,imat,ffac)
*                                                                      *
*        generate and bank photons from a neutron collision.           *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for  REDUCTION_COUNTER
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2

      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      use moddas_material
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


      common /eparm/  esmax, esmin, emin(20)
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)

*-----------------------------------------------------------------------

      character ht*10
      dimension en(30),yd(2)

      data f/1.00000000001d0/

      data en / 1.39e-10, 1.52e-7, 4.14e-7, 1.13e-6, 3.06e-6,
     1 8.32e-6, 2.26e-5, 6.14e-5, 1.67e-4, 4.54e-4, 1.235e-3, 3.35e-3,
     2 9.12e-3, 2.48e-2, 6.76e-2, 0.184, 0.303, 0.500, 0.823, 1.353,
     3 1.738, 2.232, 2.865, 3.68, 6.07, 7.79, 10.0, 12.0, 13.5, 15.0/

      data iwarning/0/ ! T.Sato 2022/12/07

*-----------------------------------------------------------------------
*     initialize.
*-----------------------------------------------------------------------

      kp=mpan
      nt=0
      rfact = 0d0

*-----------------------------------------------------------------------
*     with pikmt, select the nuclide and generate exactly one photon.
*-----------------------------------------------------------------------

         if(npikmt.ne.0) then

            io=iex
            t=rang()*totgp1
            mk=imat

            do m=jmd(1+mk),jmd(1+mk+1)-1

*-----------------------------------------------------------------------

         lem = m-jmd(1+mk)+1
         dsmax = das_kmate(kmate(mk)+(lem-1)*5+12)

      if( dsmax .gt. emin(2) ) then

*-----------------------------------------------------------------------

               t=t-tgp(lme(1,m))*fme(m)
               if(t.lt.0.) goto 20

            end if
*-----------------------------------------------------------------------

            end do

                write(ErrCha,*) 'Error : in prdgam, no isotope found'
                ErrID = 'L:2181/R:prdgam/F:ggm05.f' !E04_010_001
                call ErrWrite(ErrID,ErrCha)

                  call parastop( 416 )

   20       iex=lme(1,m)
            gp=totgp1/totm
            kp=ipan(icl)+m-jmd(1+mk)
            ni=1
            lg=1
            nt=nxs(15,iex)

            go to 150

         end if

*-----------------------------------------------------------------------
*     calculate the probability of production of a photon.
*-----------------------------------------------------------------------

         l = jxs(12,iex)+ktc(1,iex)
         gp = (xss(l-1)+rtc(1,iex)*(xss(l)-xss(l-1)))
     &      / rtc(5,iex)

         if(gp.eq.0.) return

*-----------------------------------------------------------------------
*     find the neutron energy group in the en table.
*     below emin, no photon
*-----------------------------------------------------------------------

         if(jxs(13,iex).ne.0) goto 90

         ib = 31
         ic = 1
   40    if(ib-ic.eq.1) goto 60
         ih = (ic+ib)/2
         if(erg.lt.en(ih)) goto 50
         ic = ih
         goto 40

   50    ib = ih
         goto 40

   60    lg = jxs(12,iex)+nxs(3,iex)+(ic-1)*20-1

      do 70 it = 1, 20
   70    if(xss(lg+it).gt.elc(2)) goto 80

         return

   80    gp = gp*(1.-(it-1)*.05)

*-----------------------------------------------------------------------
*     modify gamma production by ratio of unresolved total xsec to
*     smooth (average) total xsec.
*-----------------------------------------------------------------------

   90    if(rtc(11,iex).lt.0.) goto 95

         ic = ktc(1,iex)+jxs(1,iex)-1
         n = nxs(3,iex)
         tz = rtc(5,iex)/
     &    (xss(ic+n)+rtc(1,iex)*(xss(ic+n+1)-xss(ic+n)))
         gp = gp*tz

*-----------------------------------------------------------------------
*     decide how many photons to make.
*     gp : photon production probability
*-----------------------------------------------------------------------

   95    gp = wgt*gp
         ni = 1
!           produce analog gammas for pulse height tally.
!           get target weight, t1, for dxtran, importances, or windows.

         t1 = max(gwt(icl),-gwt(icl)*wg0) * ffac

*-----------------------------------------------------------------------
*     determine split number, ni, or roulette.
*-----------------------------------------------------------------------

  130    if(gp.lt.t1) goto 140

         if(t1.ne.0.) ni = min(10,int(gp/(5.*t1)+1.))
         if(rtc(11,iex).ge.0.) gp = gp/tz
         goto 150

  140    if(gp.le.t1*rang()) return
         gp = t1
         if(rtc(11,iex).ge.0.) gp = gp/tz

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
         ts = -1.

*-----------------------------------------------------------------------
*     sample from full photon-production data, if any.
*-----------------------------------------------------------------------

         ipsc = 102
         if(jxs(13,iex).eq.0) goto 380

*-----------------------------------------------------------------------
*     sample the reaction that produced the photon.
*-----------------------------------------------------------------------

         ixre = 1
      if(nxs(6,iex).eq.1.and.rtc(11,iex).lt.0.) goto 350
         r = 0.
      if(nt.ne.0) goto 190
         l = jxs(12,iex)+ktc(1,iex)
         rfact = rang()
         r = rfact*(xss(l-1)+rtc(1,iex)*(xss(l)-xss(l-1)))
         s = 0.
         tt = 0.
  190    md = 0

      do 290 ixre = 1, nxs(6,iex)

         ts = -1.
         if(nt.eq.0) goto 220

      do 200 ik = 1, nt
  200    if(ixre.eq.nint(xss(jxs(32,iex)+2*ik-2))) goto 210
         goto 290

  210    pik(kpik+ik) = 0.
  220    l = jxs(15,iex)+nint(xss(jxs(14,iex)+ixre-1))+1
         if(nint(xss(l-2)).eq.13) goto 250

*-----------------------------------------------------------------------
*     mf=12 -- partial cross section = yield * neutron cross section.
*-----------------------------------------------------------------------

      do 230 jy=1,2
  230    yd(jy) = acefcn(l,xss(jxs(1,iex)
     &          + ktc(1,iex)+jy-2)*f,ln)
         if(yd(1).eq.0..and.yd(2).eq.0.) goto 290
         ix = nint(xss(l-1))
         if(ix.eq.md) goto 240
         md = ix
         x1 = 0.
         x2 = 0.
         is = jxs(21,iex)+1
         if(ix.gt.0) is = jxs(7,iex)
     &                  + nint(xss(jxs(6,iex)+ix-1))
         ic = min(ktc(1,iex)-nint(xss(is-1))+1,nint(xss(is)))
         if(ic.gt.0) x1 = xss(is+ic)
         ic = min(ic+1,nint(xss(is)))
         if(ic.gt.0) x2 = xss(is+ic)
  240    t = x1*yd(1)+(x2*yd(2)-x1*yd(1))*rtc(1,iex)
         if(t.eq.0.) goto 290

*-----------------------------------------------------------------------
*     compute change in photon production due to change in capture
*     or fission when using unresolved tables rather than smooth.
*     multiply photon weight by {partial-prob-table/partial-smooth}*
*     {total-smooth/total-prob-table} for unbiased.
*-----------------------------------------------------------------------

         if(rtc(11,iex).lt.0.) goto 270
         mt = nint(xss(jxs(13,iex)+ixre-1))/1000
         if(mt.eq.18.or.mt.eq.19) tf = rtc(12,iex)/
     &    (x1+(x2-x1)*rtc(1,iex))
         if(mt.eq.18.or.mt.eq.19) ts = tf
         if(mt.eq.102) tc = rtc(14,iex)
     &                    / (x1+(x2-x1)*rtc(1,iex))
         if(mt.eq.102) ts = tc
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

         if(rtc(11,iex).ge.0.) then
            write(ErrCha,*) 'Error : in prdgam, '//
     &    'mf=13 photon prod not allowed with unresolved prob tables.'
           ErrID = 'L:2382/R:prdgam/F:ggm05.f' !E04_011_001
           call ErrWrite(ErrID,ErrCha)

            goto 450
            call parastop( 417 )

         end if

         t = xss(ic+l)+rtc(1,iex)*(xss(ic+l+1)-xss(ic+l))
  270    if(t.eq.0.) goto 290
         if(nt.eq.0) goto 280
         pik(kpik+ik) = t*xss(jxs(32,iex)+ik*2-1)
         r = r+pik(kpik+ik)
         goto 290
  280    if(t.gt.tt) kx = ixre
         tt = max(t,tt)
         s = s+t
         if(s.gt.r) goto 350
  290    continue

*-----------------------------------------------------------------------

         if(nt.eq.0) goto 320
         if(r.le.0.) goto 450
         ik = 1
         if(nt.eq.1) goto 310
         t = rang()*r

      do 300 ik = 1, nt
         t = t-pik(kpik+ik)
  300    if(t.lt.0.) goto 310

            write(ErrCha,*) 'Error : in prdgam, '//
     &                 'no photon mt selected.'
            ErrID = 'L:2416/R:prdgam/F:ggm05.f' !E04_012_001
            call ErrWrite(ErrID,ErrCha)

            call parastop( 418 )

  310    ixre = nint(xss(jxs(32,iex)+ik*2-2))
         wgt = wgt*r/(tgp(iex)*xss(jxs(32,iex)+ik*2-1))

         if(rtc(11,iex).lt.0.) goto 350
         mt = nint(xss(jxs(13,iex)+ixre-1))/1000
         if(mt.eq.18.or.mt.eq.19) ts = tf
         if(mt.eq.102) ts = tc
         goto 350

*-----------------------------------------------------------------------
*     handle failure of the reaction cross sections to add up.
*-----------------------------------------------------------------------

  320    call zaid(2,ht,ixl(1,iex))

         if( iwarning.eq.0 ) then

          if ( s.le.r .and. (r-s).lt.r*1d-6
     &           .and. (1d0-rfact).lt.1d-6 ) then
            write(ErrCha,*) 'Warning : in prdgam, '//
     &      'a photon production in a neutron reaction has failed '//
     &      'due to random number.'
            ErrID = 'L:2443/R:prdgam/F:ggm05.f'
            call ErrWrite(ErrID,ErrCha)

          else
            write(ErrCha,*) 'Warning : in prdgam, '//
     &      'no photon-production xsec was found in nuclear data '//
     &      'library for '//ht
            ErrID = 'L:2450/R:prdgam/F:ggm05.f' !W04_002_001
            call ErrWrite(ErrID,ErrCha)

          end if

            iwarning=1
         end if

         ixre = kx

*-----------------------------------------------------------------------
*     sample the energy and direction of the photon.
*-----------------------------------------------------------------------

  350    erg = es

         if(ts.ge.0.) wgt = wgt*ts
         mtp = nint(xss(jxs(13,iex)+ixre-1))
         ia = jxs(17,iex)
         ka = 0
         if(jxs(16,iex).ne.0)
     &          ka = nint(xss(jxs(16,iex)+ixre-1))
         id = jxs(19,iex)
         kd = nint(xss(jxs(18,iex)+ixre-1))

         call xstcas(1,1,zero,ia,ka,id,kd)

         if(ixcos.ne.0) ipsc = 8

*-----------------------------------------------------------------------
*     energy and angle, if energy is smaller then ecut, neglect
*-----------------------------------------------------------------------

         if(colout(1,1).lt.elc(2)) goto 450

            call dtcos(colout(2,1),uold,uuu,0,irdm)

            erg = colout(1,1)

         goto 390

*-----------------------------------------------------------------------
*     otherwise sample from equi-probable energies table.
*-----------------------------------------------------------------------

  380 continue

         if(it.eq.20) erg = xss(lg+20)
         if(it.ne.20) erg = xss(lg+it+int(rang()*(21-it)))

      if(jxs(23,iex).ne.0.and.rtc(11,iex).ge.0.) then

         write(ErrCha,*) 'Error : in prdgam, '//
     &   'old photon prod not allowed with unresolved prob tables.'
         ErrID = 'L:2504/R:prdgam/F:ggm05.f' !E04_013_001
         call ErrWrite(ErrID,ErrCha)

         call parastop( 419 )

      end if

         call isos(uuu,0)

*-----------------------------------------------------------------------
*     booking
*-----------------------------------------------------------------------

  390 continue

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
                     qclusts(10,nclsts) = 0.0
                     qclusts(11,nclsts) = 0.0
                     qclusts(12,nclsts) = 0.0
                     ixcoss(nixcos) = ixcos  ! MATSUDA 2022.12.12  for T-Point
                     aevts(iaevt+13,ireg) =
     &               aevts(iaevt+13,ireg) + wgt
                     bevts(ibevt+13,imat) =
     &               bevts(ibevt+13,imat) + wgt

                     aevts(iaevt+14,ireg) =
     &               aevts(iaevt+14,ireg) + wgt * erg
                     bevts(ibevt+14,imat) =
     &               bevts(ibevt+14,imat) + wgt * erg
                     numpal(14) = numpal(14) + 1
                     rumpal(14) = rumpal(14) + wgt

  450 continue

*-----------------------------------------------------------------------

         wgt = wg0
         erg = eg0
         ipt = 1
         totm = st

         if(npikmt.ne.0) iex = io

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine xstcas(ls,ip,q,ia,ka,id,kd)
*                                                                      *
*        sample the emission energy and scattering angle from the      *
*        appropriate law data in the laboratory coordinate system.     *
*                                                                      *
*        ls - the current particle index in the colout array           *
*        ip - the ipt particle type for the incident particle          *
*         q - the q value for the reaction being sampled               *
*        ia - the first word of the relevant AND block in the xss array*
*        ka - the offset to the first word of the table in AND block   *
*        id - the first word of the relevant DLW block in the xss array*
*        kd - the offset to the first word of the table in DLW block   *
*                                                                      *
*        returns the sample emission parameters in the lab system:     *
*          colout(1,ls) - the emission energy                          *
*          colout(2,ls) - the emission scattering angle                *
*                                                                      *
*        Law 4/44/61 makes use of a biased distribution which          *
*        can affect the outgoing particle weight, wgt.                 *
*                                                                      *
*        ipsc, kdb, and tpd also modified.                             *
*                                                                      *
*        Last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      parameter (ep=0.000001)
      character ht*10

*-----------------------------------------------------------------------
*     select the law.
*-----------------------------------------------------------------------

         nx = 0
         colout(3,ls) = 0.
   10    t1 = rang()
         n = id+kd
         goto 25

   20    n = id+nint(xss(n-1))
   25    if(nint(xss(n-1)).eq.0) goto 30
         t1 = t1-acefcn(n+2,erg,ln)
         if(t1.ge.0.) goto 20

*-----------------------------------------------------------------------
*     use the selected law to sample the energy (and possibly angle).
*     if law samples without error, go to sample angle or coordinate
*     transform as appropriate.
*-----------------------------------------------------------------------

   30    colout(1,ls) = -huge
         lw = nint(xss(n))
         iw = id-1+nint(xss(n+1))

         if(lw.eq.1) goto 40
         if(lw.eq.2) goto 60
         if(lw.eq.3) goto 70
         if(lw.eq.4) goto 80
         if(lw.eq.5) goto 160
         if(lw.eq.7) goto 170
         if(lw.eq.9) goto 190
         if(lw.eq.11) goto 210
         if(lw.eq.22) goto 230
         if(lw.eq.24) goto 242
         if(lw.eq.33) goto 70
         if(lw.eq.44) goto 80
         if(lw.eq.61) goto 80
         if(lw.eq.66) goto 245
         if(lw.eq.67) goto 255
         goto 300

*-----------------------------------------------------------------------
* >>>>>  law 1 (from endf Law 1) -- tabular equiprobable energy bins.
*-----------------------------------------------------------------------

   40    call xsrtbl(iw,ic,r,ln)
         nt = nint(xss(iw+ln))
         iw = iw+ln+nt*(ic-1)
         k = int(rang()*(nt-1)+1)
         fr = rang()
         if(r.ne.0.) goto 50

*-----------------------------------------------------------------------
*        sample from single table.
*-----------------------------------------------------------------------

         colout(1,ls) = xss(iw+k)+fr*(xss(iw+k+1)-xss(iw+k))
         goto 260

*-----------------------------------------------------------------------
*        sample by scaled interpolation between tables.
*-----------------------------------------------------------------------

   50    t1 = xss(iw+1)+r*(xss(iw+1+nt)-xss(iw+1))
         i = iw
         if(rang().le.r) i = i+nt
         colout(1,ls) = t1
     &                + (xss(i+k)+fr*(xss(i+k+1)-xss(i+k))-xss(i+1))
     &                * (xss(iw+nt)+r*(xss(iw+2*nt)-xss(iw+nt))-t1)
     &                / (xss(i+nt)-xss(i+1))
         goto 260

*-----------------------------------------------------------------------
* >>>>>  law 2 -- discrete photon lines.
*-----------------------------------------------------------------------
*        photon production from neutrons only.
*        (see endf-102 rev. 2/97 manual page 12.3).
*-----------------------------------------------------------------------

   60    colout(1,ls) = xss(iw+1)
         if(nint(xss(iw)).eq.2) colout(1,ls) = colout(1,ls)
     &                         + erg*awn(iex)/(awn(iex)+1.)
         goto 260

*-----------------------------------------------------------------------
* >>>>>  law 3 & 33 -- level scattering.
*        law 3 applies only for neutron in-neutron out scattering.
*        law 33 for a more general combination of particle types.
*-----------------------------------------------------------------------

   70    colout(1,ls) = xss(iw+1)*(erg-xss(iw))
         goto 260

*-----------------------------------------------------------------------
* >>>>>  law 4 (from endf law 1) -- continuous erg tabular distribution.
* >>>>>  law 44 (from endf law 1) -- kalbach-87 correlated formalism.
* >>>>>  law 61 (from endf law 1) -- correlated tab energy-angle dist.
*-----------------------------------------------------------------------

   80    call xsrtbl(iw,ic,r,ln)
         nr = nint(xss(iw))
         lb = iw-1+ln+ic
         lc = id+nint(xss(lb))
         ld = lc
         lf = lc

*-----------------------------------------------------------------------
*        xss(lc,ld,lf) is an overloaded variable.  it contains the
*        number of points in the emission distribution and if part
*        of the continuum has been expunged, it contains 0.5 times
*        the cumulative probability of the portion expunged.
*-----------------------------------------------------------------------

         np = int(xss(lc)+ep)
         mp = np
         jj = nint(xss(lc-1))
         nd = 0
         if(jj.lt.9999) nd = jj/10
         if(r.eq.0.) goto 90
         ld = id+nint(xss(lb+1))
         mp = int(xss(ld)+ep)
         t1 = xss(lc+nd+1)+r*(xss(ld+nd+1)-xss(lc+nd+1))
         t2 = xss(lc+np)+r*(xss(ld+mp)-xss(lc+np))
         ra = rang()
         if(ra.ge.r) goto 90
         lf = ld
         jj = nint(xss(ld-1))

         if(jj.lt.9999.and.jj/10.ne.nd) then
            write(ErrCha,*) 'Error : in xstcas, '//
     &      'wrong number of discrete lines for law 4/44'
            ErrID = 'L:2756/R:xstcas/F:ggm05.f' !E04_014_001
            call ErrWrite(ErrID,ErrCha)

            call parastop( 420 )
         end if

         if(jj.gt.10000) nd = 0
   90    if(jj.eq.9999) return
         jj = jj-nd*10
         r1 = rang()

*-----------------------------------------------------------------------
*        check for a hit in the discrete part.
*        use histogram or corresponding-point interpolation.
*-----------------------------------------------------------------------

         if(nd.eq.0) goto 97

      do 95 ih = 1, nd
         cc = xss(lc+2*np+ih)+r*(xss(ld+2*mp+ih)-xss(lc+2*np+ih))
   95    if(cc.ge.r1) goto 96
         if((np.eq.nd).or.(mp.eq.nd)) goto 98
         goto 97

   96    t = xss(lc+ih)+r*(xss(ld+ih)-xss(lc+ih))
         if(t.lt.0.) t = erg*awn(iex)/(awn(iex)+1.)-t
         if(lw.ne.44) goto 152
         tpd(1) = xss(lc+ih+3*np)+r*(xss(ld+ih+3*mp)-xss(lc+ih+3*np))
         tpd(2) = xss(lc+ih+4*np)+r*(xss(ld+ih+4*mp)-xss(lc+ih+4*np))
         goto 152

*-----------------------------------------------------------------------
*        handle a hit in the continuous part.
*        use histogram or unit-base interpolation inside table.
*        use scaled interpolation between tables.
*-----------------------------------------------------------------------
   98    if (int(xss(lc)+ep).lt.int(xss(ld)+ep)) then
         lf=ld
         else
         lf=lc
         endif
   97    np = int(xss(lf)+ep)
         ns = nd+1
         if(nd.ne.0) r1 = 1.-(1.-r1)*(1.-xss(lf+2*np+nd))/(1.-cc)

*-----------------------------------------------------------------------
*        adjust for the energy cutoff if necessary.
*        see subroutine expung for law 4/44 distribution. 0 < wc < 1
*-----------------------------------------------------------------------

         if(jj.lt.9999) goto 100
         jj = jj/10000
         ns = nint(xss(lf-1))-jj*10000
         wc = 2.*(xss(lf)-np)
         r1 = 1.-wc*(1.-r1)
         wgt = wgt*wc

*-----------------------------------------------------------------------
*        find the energy bin in the continuous part.
*-----------------------------------------------------------------------

  100    ic = lf+2*np+ns
         ib = lf+3*np
  110    if(ib-ic.eq.1) goto 130
         ih = (ic+ib)/2
         if(r1.lt.xss(ih)) goto 120
         ic = ih
         goto 110
  120    ib = ih
         goto 110
  130    ln = ic-2*np
         fa = xss(ln+np)
         ea = xss(ln)

*-----------------------------------------------------------------------
*        sample from linear-linear interpolated bin.
*-----------------------------------------------------------------------

         if(jj.eq.1) goto 140
         bb = (xss(ln+np+1)-fa)/(xss(ln+1)-ea)
         if(bb.eq.0) goto 140
         t = ea+(sqrt(max(zero,fa**2+2.*bb*(r1-xss(ic))))-fa)/bb
         if(lw.ne.44) goto 150
         fb = (t-xss(ln))/(xss(ln+1)-xss(ln))
         tpd(1) = xss(ic+np)+fb*(xss(ic+1+np)-xss(ic+np))
         tpd(2) = xss(ic+2*np)+fb*(xss(ic+1+2*np)-xss(ic+2*np))
         goto 150

*-----------------------------------------------------------------------
*        sample from histogram bin.
*-----------------------------------------------------------------------
  140    continue
         if( fa .eq. 0.d0) then
            t = 0.d0
         else
            t = ea+(r1-xss(ic))/fa
         end if
         if(lw.ne.44) goto 150
         tpd(1) = xss(ic+np)
         tpd(2) = xss(ic+2*np)

*-----------------------------------------------------------------------
*        use scaled interpolation between energies.
*-----------------------------------------------------------------------

  150    if(r.ne.0.)
     &   t = t1+(t-xss(lf+nd+1))*(t2-t1)/(xss(lf+np)-xss(lf+nd+1))
  152    colout(1,ls) = t
         if(lw.eq.4) goto 260
         if(lw.eq.44) goto 153
         if(lw.eq.61) goto 156
         goto 295

*-----------------------------------------------------------------------
*        sample law 44 -- kalbach-87 angular systematics.
*        tpd(1)=r, tpd(2)=a
*-----------------------------------------------------------------------

  153    if(ka.ne.-1.or.ntyn.ge.0) goto 295
         ipsc = 14
         ixcos = 0  ! MATSUDA 2022.12.12  for T-Point
         if(rang().le.tpd(1)) goto 155
         t1 = (2.*rang()-1.)*sinh(tpd(2))
         colout(2,ls) = log(t1+sqrt(t1**2+1.))/tpd(2)
         goto 280

  155    r2 = rang()
         colout(2,ls) = log(r2*exp(tpd(2))+(1.-r2)*exp(-tpd(2)))/tpd(2)
         goto 280

*-----------------------------------------------------------------------
*        sample law 61 -- tabulated angular distribution.
*-----------------------------------------------------------------------

  156    ipsc = 16

*-----------------------------------------------------------------------
*        if jj=1 (i.e., histogram on e-primes) always use "ic."
*        if jj=2 (lin-lin on e-primes), use the distribution for
*        the e-prime closest to "r1" (in cdf space).
*-----------------------------------------------------------------------

         lb = ic+np
         if(jj.ne.1.and.(xss(ib)-r1.lt.r1-xss(ic))) lb = lb+1

*-----------------------------------------------------------------------
*        sample from appropriate distribution.  unlike the and block, in
*        law 61 only isotropic or tabular angular information is passed.
*-----------------------------------------------------------------------

         ixcos = 0
         lm = nint(xss(lb))
         if(lm.eq.0) colout(2,ls) = 2.*rang()-1.
         if(lm.gt.0) colout(2,ls) = acecos(id,-lm)

         if(lm.lt.0) then
            write(ErrCha,*) 'Error : in xstcas, '//
     &      'data: law 61 contained a bad table pointer'
            ErrID = 'L:2914/R:xstcas/F:ggm05.f' !E04_015_001
            call ErrWrite(ErrID,ErrCha)

            call parastop( 421 )

         end if

         goto 280

*-----------------------------------------------------------------------
* >>>>>  law 5 (from endf law 5) -- general evaporation spectrum.
*-----------------------------------------------------------------------

  160    t1 = acefcn(iw,erg,ln)
         i = iw+ln+1+int(rang()*(nint(xss(iw+ln))-1))
         colout(1,ls) = t1*(xss(i)+rang()*(xss(i+1)-xss(i)))
         goto 260

*-----------------------------------------------------------------------
* >>>>>  law 7 (from endf law 7) -- simple maxwell fission spectrum.
*-----------------------------------------------------------------------

  170    t1 = acefcn(iw,erg,ln)
         t3 = erg-xss(iw+ln)
         if(t3.le.0.) goto 260
         ilaw = 0
  180    t4 = rang()**2
         ilaw = ilaw + 1
         t2 = t4+rang()**2
         if(t2.gt.1.) goto 180
         t2 = log(rang())*t4/t2
         colout(1,ls) = -t1*(t2+log(rang()))

*-----------------------------------------------------------------------
*        reject if outside range 0 ... e-u
*-----------------------------------------------------------------------

         if( colout(1,ls).gt.t3 .and. ilaw .gt. 1000 ) then
             colout(1,ls) = t3
             write(ErrCha,*) 'Warning at law 7, t3 =', t3
             ErrID = 'L:2954/R:xstcas/F:ggm05.f' !W04_003_001
             call ErrWrite(ErrID,ErrCha)
         end if
         if(colout(1,ls).gt.t3) goto 180
         goto 260

*-----------------------------------------------------------------------
* >>>>>  law 9 (from endf law 9) -- evaporation spectrum.
*-----------------------------------------------------------------------

  190    t1 = acefcn(iw,erg,ln)
         t2 = erg-xss(iw+ln)
         if(t2.le.0.) goto 260
         ilaw = 0
  200    fr = rang()
         ilaw = ilaw + 1
         colout(1,ls) = -t1*log(fr*rang())

*-----------------------------------------------------------------------
*       reject if outside range 0 ... e-u
*-----------------------------------------------------------------------

         if( colout(1,ls).gt.t2 .and. ilaw .gt. 1000 ) then
             colout(1,ls) = t2
             write(ErrCha,*) 'Warning at law 9, t2 =', t2
             ErrID = 'L:2979/R:xstcas/F:ggm05.f' !W04_003_002
             call ErrWrite(ErrID,ErrCha)
         end if
         if(colout(1,ls).gt.t2) goto 200
         goto 260

*-----------------------------------------------------------------------
* >>>>>  law 11 (from endf law 11) -- energy dependent watt spectrum.
*-----------------------------------------------------------------------

  210    t1 = acefcn(iw,erg,ln)
         t2 = acefcn(iw+ln,erg,lb)
         if(erg.le.xss(iw+ln+lb)) goto 260
         t5 = sqrt((1.+.125*t1*t2)**2-1.)+1.+.125*t1*t2
         ilaw = 0
  220    t = -log(rang())
         ilaw = ilaw + 1
         colout(1,ls) = t1*t5*t
         if(((1.-t5)*(1.+t)-log(rang()))**2.gt.t2*colout(1,ls).and.
     &       ilaw .lt. 1000 ) goto 220
         if( ilaw .ge. 1000 )
     &  write(ErrCha,*) 'Warning at law 11, colout(1,) =', colout(1,ls)
        ErrID = 'L:3001/R:xstcas/F:ggm05.f' !W04_003_003
        call ErrWrite(ErrID,ErrCha)

         goto 260

*-----------------------------------------------------------------------
* >>>>>  law 22 (from uk law 2) -- tabular linear functions.
*-----------------------------------------------------------------------

  230    call xsrtbl(iw,ic,r,ln)
         ie = id-1+nint(xss(iw+ln+ic-1))
         nf = nint(xss(ie))
  235    iw = ie
         fr = rang()
  240    iw = iw+1
         fr = fr-xss(iw)
         if(fr.ge.0.) goto 240
         if(iw.gt.ie+nf) goto 235
         colout(1,ls) = xss(iw+2*nf)*(erg-xss(iw+nf))
         goto 260

*-----------------------------------------------------------------------
* >>>>>  law 24 (from uk law 6) -- tabular energy multipliers.
*-----------------------------------------------------------------------

  242    call xsrtbl(iw,ic,r,ln)
         i = iw+ln+1+nint(xss(iw+ln))*(ic-1)
     &     + int(rang()*(nint(xss(iw+ln))-1))
         colout(1,ls) = erg*(xss(i)+rang()*(xss(i+1)-xss(i)))
         goto 260

*-----------------------------------------------------------------------
* >>>>>  law 66 (from endf law 6) -- n-body phase space distribution.
*-----------------------------------------------------------------------

  245    nb = nint(xss(iw))
         ap = xss(iw+1)
         if(ipt.gt.2) ap = ap*gpt(1)/gpt(ipt)

  246    r = rang()**2
         s = r+rang()**2
         if(s.gt.1.) goto 246
         x = -r*log(s)/s-log(rang())

  248    r = rang()**2
         s = r+rang()**2
         if(s.gt.1.) goto 248
         p = rang()
         goto (252,251,250) nb - 2

  250    p = p*rang()*rang()
  251    p = p*rang()
  252    y = -r*log(s)/s-log(p)
         t = x/(x+y)
         aw = awn(iex)
         colout(1,ls) = t*((ap-1.)/ap)*(erg*aw/(aw+1.)+q)
         colout(2,ls) = 2.*rang()-1.
         ixcos = 0
         goto 280

*-----------------------------------------------------------------------
* >>>>>  law 67 (endf/b-vi law 7) -- correleted energy-angle scatter.
*-----------------------------------------------------------------------

  255    call xsrtbl(iw,ic,r,ln)
         cs = acecos(ia,ka)
         ipsc = 15
         colout(2,ls) = cs
         colout(1,ls) = acecs6(0,id,iw,ic,r,cs)
         if(ntyn.le.0) goto 295
         if(colout(1,ls).lt.0.) goto 300
         goto 290

*-----------------------------------------------------------------------
*     if not correlated energy-angle, calculate the cosine.
*-----------------------------------------------------------------------

  260    colout(2,ls) = acecos(ia,ka)

*-----------------------------------------------------------------------
*     adjust if energy and cosine are given in center-of-mass system.
*-----------------------------------------------------------------------

  280 if(colout(1,ls).lt.0.) goto 300
      ergace = colout(1,ls)
      if(ntyn.ge.0) goto 290

*-----------------------------------------------------------------------
*        formulas below are from p. 2 of x-6:res-93-68.
*        these formulas assume two-body kinematics.
*        seamon's formulas specify atomic weight ratios (to neutron)
*        for incident particle (a), awr=gpt(ipt_incident)/gpt(1)
*        for exiting particle (b), awr=gpt(ipt)/gpt(1)
*        for target (a), awr=awn(iex)
*-----------------------------------------------------------------------

         a1 = gpt(ip)/gpt(1)
         a2 = gpt(ipt)/gpt(1)
         a3 = awn(iex)
         t1 = colout(1,ls)
         t2 = erg*a1*a2/(a3+a1)**2
         t3 = 2.*sqrt(a2*a1*erg*colout(1,ls))*colout(2,ls)/(a3+a1)
         t4 = t1+t2+t3
         s1 = colout(2,ls)*sqrt(colout(1,ls)/t4)
         s2 = sqrt(a1*a2*erg/t4)/(a3+a1)
         colout(1,ls) = t4
         colout(2,ls) = s1+s2

*-----------------------------------------------------------------------
*     resample energy up to 100 times if > emx(1).
*-----------------------------------------------------------------------

  290    if(colout(1,ls).le.emx(ipt)) return

         if( ipt .gt. 3 ) return

         nx=nx+1
         if(nx.lt.100)go to 10

         ErrCha =""
         ErrID = 'L:3121/R:xstcas/F:ggm05.f' !W04_004_001
         call ErrWrite(ErrID,ErrCha)

         write(*,'(''Warning : in xstcas, '',
     &    ''erg>emx 100 times in one collision.''/''ipt='',i4,3e13.4)')
     &     ipt, emx(ipt),colout(1,ls), erg

         return

*-----------------------------------------------------------------------
*      print debug information for cross-section table errors.
*-----------------------------------------------------------------------

  295    colout(1,ls) = huge
  300    call zaid(2,ht,ixl(1,iex))

         ErrCha = ''
         ErrID = 'L:3138/R:xstcas/F:ggm05.f' !E04_016_001
         call ErrWrite(ErrID,ErrCha)
         write(*,310)ht,erg,ixre,mtp,ntyn,lw,colout(1,ls)

  310    format(/30h error in cross-section table ,a10/12h energy in =,
     1    1pe12.4,5x,16hreaction index =,i3,5x,4hmt =,i4,5x,4hty =,i4,
     2    5x,5hlaw =,i3,5x,12henergy out =,1pe12.4)

         if(colout(1,ls).eq.-huge) then
            write(ErrCha,*) 'Error in xstcas, '//
     &      'an inappropriate or non-existent law was selected.'
          ErrID = 'L:3149/R:xstcas/F:ggm05.f' !E04_017_001
          call ErrWrite(ErrID,ErrCha)

            return
            call parastop( 423 )

         end if
         if(colout(1,ls).lt.0.) then
            write(ErrCha,*) 'Error in xstcas, '//
     &      'emission energy was negative.'
            ErrID = 'L:3159/R:xstcas/F:ggm05.f' !E04_018_001
            call ErrWrite(ErrID,ErrCha)

            call parastop( 423 )
         end if

         if(colout(1,ls).eq.huge) then
            write(ErrCha,*) 'Error in xstcas, '//
     &      'faulty cross-section data.'
            ErrID = 'L:3168/R:xstcas/F:ggm05.f' !E04_019_001
            call ErrWrite(ErrID,ErrCha)

            call parastop( 423 )
         end if

         if(colout(1,ls).gt.erg) then
            write(ErrCha,*) 'Error in xstcas, '//
     &      'emission energy exceeds incident energy.'
            ErrID = 'L:3177/R:xstcas/F:ggm05.f' !E04_020_001
            call ErrWrite(ErrID,ErrCha)

            call parastop( 423 )
         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine xsrtbl(l,ii,r,ln)
*                                                                      *
*        get interpolation parameters r and ii of the energy in the    *
*        table at xss(l).                                              *
*        data:  nr,nbttc(i=1,nr),int(i=1,nr),nf,e(i=1,nf)              *
*        Last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

         nr = nint(xss(l))
         ie = l+2*nr+1
         nf = nint(xss(ie))
         ln = 2*(nr+1)+nf
         r = 0.

*-----------------------------------------------------------------------
*     use extreme value if energy is off either end of the table.
*-----------------------------------------------------------------------

         if(erg.ge.xss(ie+nf)) goto 70
         if(erg.le.xss(ie+1)) goto 80

*-----------------------------------------------------------------------
*     binary search for the location of the energy in the table.
*-----------------------------------------------------------------------

         ic = ie+1
         ib = ie+nf
   10    if(ib-ic.eq.1) goto 30
         ih = (ic+ib)/2
         if(erg.lt.xss(ih)) goto 20
         ic = ih
         goto 10
   20    ib = ih
         goto 10
   30    ii = ic-ie

*-----------------------------------------------------------------------
*     calculate interpolation fraction r unless int=1 (histogram).
*-----------------------------------------------------------------------

         if(nr.eq.0) goto 60

      do 40 n = 1, nr
   40    if(ib-ie.le.nint(xss(l+n))) goto 50
         n = nr
   50    if(nint(xss(l+nr+n)).eq.1) return
   60    if(erg-xss(ic).lt.1.e-6*(xss(ib)-xss(ic))) return
         r = (erg-xss(ic))/(xss(ib)-xss(ic))
         return

   70    ii = nf
         return
   80    ii = 1

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function acecos(ia,ka)
*                                                                      *
*        sample a cosine from the distribution at ia,ka.  return the   *
*        cosine (acecos) and the index in xss of the first word of the *
*        cosine table used (ixcos). ka=0: isotropic. ka<0: tabular only*
*        32-equiprobable bin cosine data:                              *
*           ne,(e(i),i=1,ne),(lmu(i),i=1,ne),(mu(lmu+j),j=1,33)        *
*        tabular probability angular distribution data:                *
*           ne,(e(i),i=1,ne),(lmu(i),i=1,ne),jj,np,                    *
*           (mu(lmu+j),j=1,np),(pdf(lmu+j),j=1,np),(cdf(lmu+j),j=1,np) *
*           where jj is the interpolation flag at lmu(i).              *
*                                                                      *
*        Last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------
*     find the cosine table by binary search on the energy table.
*-----------------------------------------------------------------------

         if(ka.eq.0) goto 40
         lm = ka
         if(ka.lt.0) goto 50
         ic = ia+ka
         n = nint(xss(ic-1))
         ib = ic-1+n
   10    if(ib-ic.eq.1) goto 30
         ih = (ic+ib)/2
         if(erg.lt.xss(ih)) goto 20
         ic = ih
         goto 10
   20    ib = ih
         goto 10

*-----------------------------------------------------------------------
*     sample between adjoining tables by interpolation fraction.
*-----------------------------------------------------------------------

   30    if(rang()*(xss(ib)-xss(ic)).lt.erg-xss(ic)) ic = ib
         lm = nint(xss(ic+n))
         if(lm.eq.0) goto 40
         if(lm.lt.0) goto 50

*-----------------------------------------------------------------------
*     sample from table of 32 equiprobable cosine groups.
*-----------------------------------------------------------------------

         t1 = rang()*32.
         kr = t1
         ixcos = ia-1+lm
         acecos = xss(ixcos+kr)+(t1-kr)*(xss(ixcos+kr+1)-xss(ixcos+kr))
         return

*-----------------------------------------------------------------------
*     isotropic case
*-----------------------------------------------------------------------

   40    acecos = 2.*rang()-1.
         ixcos = 0
         return

*-----------------------------------------------------------------------
*     tabular probability angular distribution.
*-----------------------------------------------------------------------

   50    k = ia-1-lm
         ixcos = -k
         jj = nint(xss(k))
         np = nint(xss(k+1))

*-----------------------------------------------------------------------
cKN p+d -> 3He
*-----------------------------------------------------------------------

         if( np .le. 1 ) then
           acecos = 1.0
           return
         end if

         rn = rang()

*-----------------------------------------------------------------------
*     binary search of cumulative density function.
*-----------------------------------------------------------------------

         ic = k+2*np+2
         ib = k+3*np+1

   60    if(ib-ic.eq.1) goto 80
         ih = (ic+ib)/2
         if(rn.lt.xss(ih)) goto 70
         ic = ih
         goto 60

   70    ib = ih
         goto 60

   80    fa = xss(ic-np)
         ca = xss(ic-2*np)
cABE change @2014/08/14 to avoid break
         if (fa .ne. 0) acecos = ca+(rn-xss(ic))/fa
         if(jj.eq.1) return

         bb = (xss(ic-np+1)-fa)/(xss(ic-2*np+1)-ca)
         if(bb.ne.0.)
     &    acecos = ca+(sqrt(max(zero,fa**2+2.*bb*(rn-xss(ic))))-fa)/bb

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine samvel(e2,r2,ar,vr,k2)
*                                                                      *
*        sample the velocity of the target nucleus for the neutron free*
*        gas thermal treatment.                                        *
*        should be used for elastic scattering only.                   *
*        e2=relative energy in rest frame of target motion.            *
*        r2, k2 are interpolation fraction, index of energy in xss.    *
*                                                                      *
*        Last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

      dimension vr(3)

*-----------------------------------------------------------------------
c        return if in the unresolved prob table treatment.
*-----------------------------------------------------------------------

         if(rtc(11,iex).ge.0.) return

*-----------------------------------------------------------------------
*    sample the velocity of the target nucleus.
*-----------------------------------------------------------------------

         ar = awn(iex)/ttn
         ycn = sqrt(erg*ar)

   10    if(rang()*(ycn+1.12837917).le.ycn) goto 20
         r1 = rang()
         z2 = -log(r1*rang())
         goto 30

   20    r1 = rang()**2
         s = r1+rang()**2
         if(s.gt.1.) goto 20
         z2 = -r1*log(s)/s-log(rang())

   30    z = sqrt(z2)
         c = 2.*rang()-1.
         x2 = ycn**2+z2-2.*ycn*z*c
         if((rang()*(ycn+z))**2.gt.x2) goto 10

         call dtcos(c,uuu,vtr,0,irdm)

         vtr(1) = z*vtr(1)
         vtr(2) = z*vtr(2)
         vtr(3) = z*vtr(3)

*-----------------------------------------------------------------------
*     calculate functions of the target velocity.
*-----------------------------------------------------------------------

         vr(1) = ycn*uuu-vtr(1)
         vr(2) = ycn*vvv-vtr(2)
         vr(3) = ycn*www-vtr(3)
         ssr = sqrt(vr(1)**2+vr(2)**2+vr(3)**2)
         e2 = x2/ar
         ic = ktc(1,iex)+jxs(1,iex)-1
         if(e2.le.xss(ic+1).and.e2.ge.xss(ic)) goto 70
         if(e2.gt.xss(ic)) goto 50

   40    if(ic.eq.jxs(1,iex)) goto 60
         ic = ic-1
         if(e2.lt.xss(ic)) goto 40
         goto 60

   50    if(ic.eq.jxs(1,iex)+nxs(3,iex)-2) goto 60
         ic = ic+1
         if(e2.gt.xss(ic+1)) goto 50

   60    k2 = ic-jxs(1,iex)+1

   70    if(nty(iex).ne.2)
     &       r2 = max(zero,min(one,(e2-xss(ic))/(xss(ic+1)-xss(ic))))

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine xsielsch
*                                                                      *
*        Determine reaction channel whose neutron multiplicity is 0    *
*        Reaction channle number is sent as mtprec in common block     *
*        Activated when emode = 2                                      *
*                                                                      *
*        Last modified by T.Ogawa on 2014/07/29                        *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      common /emode1/ nevhin, nlowlv, nefiss, ntwidt, mtprec,ipcnt(70)
!$OMP THREADPRIVATE(/emode1/)

      double precision sigch(10)

*-----------------------------------------------------------------------
*     calculate the cross section for reactions MT = 102 ~ 110.
*     >200 is already included as (n,na), (n,np), etc.
*-----------------------------------------------------------------------

      i      = 1
      sigsum = 0.d0
      sigch  = 0.d0
      do while
     & (nint(xss(jxs(3,iex)+nxs(5,iex)-1 + i)) .ge. 102 .and.
     &  nint(xss(jxs(3,iex)+nxs(5,iex)-1 + i)) .le. 110)
         is = jxs(7,iex)+nint(xss(jxs(6,iex)+nxs(5,iex)-1
     &    +i))
         ic = ktc(1,iex)+1-nint(xss(is-1))
         if(ic.lt.1.or.ic.gt.nint(xss(is)).or.ic.eq.nint(xss(is)).and.
     &    rtc(1,iex).ne.0.) then
          sigch(i) = 0.d0
         elseif(rtc(1,iex).ne.0.) then
          sigch(i) = xss(ic+is) + rtc(1,iex) * (xss(ic+is+1)-xss(ic+is))
         else
          sigch(i) = xss(ic+is)
         endif
         i = i + 1
      enddo

      ransig = sum(sigch(1:10)) * unirn(dummy)

      do j = 1, i - 1
       if(ransig .le. sigch(j)) then
       mtprec = nint(xss(jxs(3,iex)+nxs(5,iex)-1 + j))
       exit
       endif
       ransig = ransig - sigch(j)
      enddo

      return
      end

************************************************************************
*                                                                      *
*                                                                      *
      subroutine xstcol(jq,ac,ireg,elrt)
*                                                                      *
*        sample an elastic or inelastic neutron collision.             *
*        jq=0 normally.                                                *
*        jq=2 for kcode sampling of fission energy only.               *
*                                                                      *
*        Last modified by K.Niita on 2009/09/30                        *
*        Last modified by T.Ogawa on 2014/04/16                        *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use SEGMENT
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      include 'err.inc'

      common /rcomon/ rcasc
      common /elreg/  elarg(kvlmax,4), mlrgn, melrg(kvlmax,2), ielusr
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /emode/  emodem, ge1, ge2, iemode
      common /ismtfl/ ismtflg
!$OMP THREADPRIVATE(/ismtfl/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /emode1/ nevhin, nlowlv, nefiss, ntwidt, mtprec,ipcnt(70)
!$OMP THREADPRIVATE(/emode1/)

*-----------------------------------------------------------------------

      character ht*10

*-----------------------------------------------------------------------

         nt = 0
         colout(3,1) = 0.
   10    if(jq.eq.2) goto 80

*-----------------------------------------------------------------------
*     user elastic and inelastic option if exist
*-----------------------------------------------------------------------

            ioel = 0
            icels = 1

      if( mlrgn .gt. 0 ) then

            do jj = 1, mlrgn

              do kk = melrg(jj,1), melrg(jj,2)

                 if( ireg .eq. idgr(kk) ) then

                    ioel = jj
                    goto 300

                 end if

              end do

            end do

  300       continue

      end if

*-----------------------------------------------------------------------
*     branching for user define cross sections
*-----------------------------------------------------------------------

      if( ioel .ne. 0 ) then

         if( ielusr .eq. 1 ) then

            call usrelst1(4,icels,ioel,erg,outdd,sigdd)

         else

            call usrelst2(4,icels,ioel,erg,outdd,sigdd)

         end if

          if( icels .eq. 0 ) then

             if( outdd .lt. 0.0 ) then

                ielch = 0
                goto 60

             else if( outdd .gt. 0.0 .and. outdd .lt. 10.0 ) then

                ielch = 1
                goto 60

             else

                ielch = 1
                r = rang() * sigdd
                goto 100

             end if

          end if

      end if

*-----------------------------------------------------------------------
*     sample elastic vs. inelastic collision.
*-----------------------------------------------------------------------

      if(nxs(5,iex).eq.0) goto 60

*-----------------------------------------------------------------------
*     inelastic collision.
*-----------------------------------------------------------------------

      if(rtc(11,iex).lt.0.) goto 15

*-----------------------------------------------------------------------
*     probability tables for unresolved resonances.
*-----------------------------------------------------------------------

         if(jxs(23,iex).eq.0) then
            write(ErrCha,*) 'Error : in xstcol, '//
     &      'flag conflict regarding unresolved prob tables.'
            ErrID = 'L:3673/R:xstcol/F:ggm05.f' !E04_021_001
            call ErrWrite(ErrID,ErrCha)

            call parastop( 424 )
         end if

         if(erg.ne.eg0) then
            write(ErrCha,*) 'Error : in xstcol, '//
     &      'erg/eg0 conflict for unresolved resonances'
            ErrID = 'L:3682/R:xstcol/F:ggm05.f' !E04_022_001
            call ErrWrite(ErrID,ErrCha)

            call parastop( 425 )
         end if

         el = rtc(11,iex)
         sr = rtc(4,iex)-rtc(3,iex)-rtc(8,iex)-el
         goto 50

*-----------------------------------------------------------------------
*     no probability tables (usual case).
*-----------------------------------------------------------------------

   15    ic = ktc(2,iex)+jxs(1,iex)-1
         n = nxs(3,iex)
         el = xss(ic+3*n)+rtc(2,iex)*(xss(ic+3*n+1)-xss(ic+3*n))
         if(erg.ne.eg0) goto 20
         sr = rtc(4,iex)-rtc(3,iex)-rtc(8,iex)-el
         goto 30

   20    sr = xss(ic+n)-xss(ic+2*n)+rtc(2,iex)
     &      * (xss(ic+n+1)-xss(ic+n)-xss(ic+2*n+1)+xss(ic+2*n))-el
         if(lfcl(icl).eq.0) goto 30
         l = jxs(21,iex)
         if(l.eq.0) goto 30
         j = l+2+ktc(2,iex)-nint(xss(l))
         sr = sr-xss(j)-rtc(2,iex)*(xss(j+1)-xss(j))

   30    if(awn(iex)*erg.gt.500.*tbt(iex)) goto 50
         a2 = awn(iex)*erg/tbt(iex)
         if(a2.lt.4.) goto 40
         el = el*a2/(a2+.5)
         goto 50

   40    a = sqrt(a2)
         b = 25.*a
         i = b
         el = el*a/(thgf(i)+(b-i)*(thgf(i+1)-thgf(i)))

   50    r = rang()*(sr+el)-el

cKN 2018/02/19 elastic rate
         elrt = el / ( sr + el )

         if(r.ge.0.) goto 100

*-----------------------------------------------------------------------
*     elastic case
*-----------------------------------------------------------------------

   60    ixre = 0
         mtp = 2
         ntyn = -99
         cmult = 1.

*-----------------------------------------------------------------------
*     sample the neutron output direction and calculate its energy.
*-----------------------------------------------------------------------

         if( icels .ne. 0 ) then

               c = acecos(jxs(9,iex),nint(xss(jxs(8,iex))))

         else

            if( ielch .eq. 0 ) then

               if( ielusr .eq. 1 ) then

                  call usrelst1(2,icels,ioel,erg,c,dum)

               else

                  call usrelst2(2,icels,ioel,erg,c,dum)

               end if

               if( icels .ne. 0 ) then

                  c = acecos(jxs(9,iex),nint(xss(jxs(8,iex))))

               else

                  colout(1,1) = erg
                  colout(2,1) = c
                  ac = -c

                  return

               end if

            else if( ielch .eq. 1 ) then

                  c = acecos(jxs(9,iex),nint(xss(jxs(8,iex))))

            end if

         end if

*-----------------------------------------------------------------------
*        user defined elastic option usrelst1
*-----------------------------------------------------------------------

cKN for target proton

         ac = -c
         if(awn(iex).lt.1.) goto 70

         t1 = 1.+awn(iex)*(awn(iex)+2.*c)
         colout(1,1) = erg*t1/(1.+awn(iex))**2
         colout(2,1) = (1.+awn(iex)*c)/sqrt(t1)

      return

*-----------------------------------------------------------------------
*     special calculation for hydrogen.
*-----------------------------------------------------------------------

   70    colout(1,1) = .5*erg*(1.+c)
         colout(2,1) = sqrt(.5+.5*c)

         return

*-----------------------------------------------------------------------
*     inelastic case.
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*     sample the reaction, using cumulative partial cross sections.
*-----------------------------------------------------------------------

   80    if(erg.ne.eg0) goto 90
         r = rang()*rtc(8,iex)
         goto 100

   90    if(rtc(11,iex).ge.0.) then
            write(ErrCha,*) 'Error : in xstcol, '//
     &      'fission conflict for unresolved resonances'
         ErrID = 'L:3821/R:xstcol/F:ggm05.f' !E04_023_001
         call ErrWrite(ErrID,ErrCha)

            call parastop( 426 )
         end if

         l = jxs(21,iex)
         j = l+2+ktc(2,iex)-nint(xss(l))
         r = rang()*(xss(j)+rtc(2,iex)*(xss(j+1)-xss(j)))

  100    tt = 0.
         kx = 0

*-----------------------------------------------------------------------

      do 120 ixre = 1, nxs(5,iex)

         if(lfcl(icl).eq.0) goto 110
         if(jq.ne.2.eqv.nint(xss(jxs(5,iex)+ixre-1)).eq.19)
     &         goto 120

*-----------------------------------------------------------------------
*     calculate the cross section for reaction whose index is ixre.
*     data:  ie,ne,(xs(i),i=1,ne)   xs(1) corresponds to es(ie).
*-----------------------------------------------------------------------

  110    is = jxs(7,iex)+nint(xss(jxs(6,iex)+ixre-1))
         ic = ktc(2,iex)+1-nint(xss(is-1))
         if(ic.lt.1.or.ic.gt.nint(xss(is)).or.ic.eq.nint(xss(is)).and.
     &    rtc(2,iex).ne.0.) goto 120
         t = xss(ic+is)
         if(rtc(2,iex).ne.0.) t = t+rtc(2,iex)
     &                               * (xss(ic+is+1)-xss(ic+is))
         if(rtc(11,iex).lt.0.) goto 115
         mn = nint(xss(jxs(3,iex)+ixre-1))
         if(mn.lt.18.or.mn.gt.19) goto 115
         if(jq.eq.2) goto 150
         t = rtc(12,iex)

  115    if(t.eq.0.) goto 120
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
          ErrID = 'L:3879/R:xstcol/F:ggm05.f' !W04_004_002
          call ErrWrite(ErrID,ErrCha)

            write(*,*) 'Warning : no reaction mt found in xstcol.'//
     &                 ' collision resampled. zaid = '//ht
            write(*,'(''  ncasc = '',i3,'' erg = '',e13.5)')
     &      int(rcasc), erg
         end if

         ixre = kx

*-----------------------------------------------------------------------
*     get the reaction type number and parameters for xstcas.
*-----------------------------------------------------------------------

  150    ntyn = nint(xss(jxs(5,iex)+ixre-1))
         mtp = nint(xss(jxs(3,iex)+ixre-1))
         q = xss(jxs(4,iex)+ixre-1)
         ia = jxs(9,iex)
         ka = nint(xss(jxs(8,iex)+ixre))
         id = jxs(11,iex)
         kd = nint(xss(jxs(10,iex)+ixre-1))

         mtprec = mtp ! 2014/7/25 ogawa

*-----------------------------------------------------------------------
*     check that the energy is within the band for this reaction.
*-----------------------------------------------------------------------

         if(nty(iex).ne.2) goto 160
         l = jxs(11,iex)+2+nint(xss(jxs(10,iex)+ixre-1))
         ie = 2*nint(xss(l))+l
         if(erg.le.xss(ie+2).or.erg.ge.xss(ie+1+nint(xss(ie+1))))
     &      goto 180

*-----------------------------------------------------------------------
*     sample the energy and scattering angle of emerging neutrons.
*-----------------------------------------------------------------------

  160    ns = 1
         if(abs(ntyn).ge.100) goto 190
         if(ntyn.ne.19) then
!           (n,xn) reactions.
            ns = abs(ntyn)
            goto 163
         else if(nxs(8,iex).ne.0.and.dnb.ne.zero) then
!           sample prompt and delayed fission.
            call xstdly(jq)
            return
         else if(jq.eq.2) then
!           jq=2: call from colidk; get prompt fission energy only.
            cmult = aint(rtc(10,iex)+rang())
            call xstcas(1,1,q,ia,ka,id,kd)
            return
         else
!           sample prompt fission.
            fn = acenu(jxs(2,iex))
            ns=fn+rang()
         endif
  163    cmult = ns

         IF(ntyn. ne. 19 .and. ns .ge. 2 .and.
     &    iemode .ge. 2 .and. erg .le. emodem .and. ! cKN 2018/02/12
     &    mathz  .ge. 3) then ! 2014 4 T.Ogawa

          ismtflg = 0
          call mainsmt(erg, ns, ntyn, q, ia, ka, id, kd) ! advanced event generator mode. Return to "sctneut" just after this mode
          IF(ismtflg .eq. 1) return
         ENDIF

      do 170 i = 1, ns
         call xstcas(i,1,q,ia,ka,id,kd) ! sample neutron energy, angle etc.
  170    if(kdb.ne.0) return

         return

*-----------------------------------------------------------------------
*     try again, up to 100 times, if the energy is out of the
*     reaction band.
*-----------------------------------------------------------------------

  180    if(erg.le.xss(jxs(1,iex)).or.erg.gt.xss(jxs(1,iex)+
     &    nxs(3,iex)-1)) goto 160

         nt = nt+1
         if( nt .le. 100 ) goto 10

         write(ErrCha,*) 'Error : in xstcol, '//
     &   'energy was not within the band for the reaction 100 times.'
         ErrID = 'L:3968/R:xstcol/F:ggm05.f' !W04_005_001
         call ErrWrite(ErrID,ErrCha)

          call parastop( 427 )

*-----------------------------------------------------------------------
*     high energy reaction other than fission with energy-dependent
*     multiplicity.  kalbach-87 (law 44) endf/b-vi.
*-----------------------------------------------------------------------
*     the location of the multiplicity table relative to dlw is
*     abs(ntyn)-100.  acefcn gets non-fission reaction multiplicity.
*-----------------------------------------------------------------------

  190    l = jxs(11,iex)+abs(ntyn)-101
         cmult = acefcn(l,erg,ln)
         if(cmult.lt.1.) goto 210
         cmult = aint(cmult+rang())

      do 200 i = 1, int(cmult)
         call xstcas(i,1,q,ia,ka,id,kd)
  200    if(kdb.ne.0) return
         return

*-----------------------------------------------------------------------
*     multiplicity < 1.0 ; reduce weight by "absorption."
*-----------------------------------------------------------------------
*     analog capture case.
*-----------------------------------------------------------------------

  210    if(erg.gt.emcf(1)) goto 220
         if(cmult+rang().le.1.) goto 280
         goto 270

*-----------------------------------------------------------------------
*     otherwise simulate capture by weight reduction.
*-----------------------------------------------------------------------

  220    if(cmult.eq.0.) goto 280
         wgt = wgt*cmult

*-----------------------------------------------------------------------
*     sample the energy and scattering angle of emerging neutrons.
*-----------------------------------------------------------------------

  270    cmult = 1.
         call xstcas(1,1,q,ia,ka,id,kd)

         return

*-----------------------------------------------------------------------
*     kill particles of zero multiplicity or by analog capture.
*     neutron capture loss
*-----------------------------------------------------------------------

  280    nter = 12

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine xstdly(jq)
*                                                                      *
*        sample prompt and delayed neutrons.                           *
*        jq=0/2=normal/kcode sampling of fission energy only.          *
*        colout(1,i)=energy; colout(2,i)=cosine scattering angle;      *
*        colout(3,i)=delay time of delayed neutrons.                   *
*        if delayed neutron biasing (4th phys:n entry) then prompt     *
*           neutrons are followed by delayed neutrons and the biased   *
*           weight is -colout(3,1).                                    *
*                                                                      *
*        Last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

         q = xss(jxs(4,iex)+ixre-1)
         ia = jxs(9,iex)
         ka = nint(xss(jxs(8,iex)+ixre))

*-----------------------------------------------------------------------
*     get total nu.
*-----------------------------------------------------------------------

         tn = acenu(jxs(2,iex))

*-----------------------------------------------------------------------
*     get delayed nu.
*-----------------------------------------------------------------------

         dn = acenu(jxs(24,iex))

*-----------------------------------------------------------------------
*     biased sampling of delayed neutrons.
*-----------------------------------------------------------------------
         if(dnb.le.zero) goto 10
         pn = tn-dn
         np = pn+rang()
         db = min(zero+10.999999-np,dnb)
         nt = np+db+rang()
         cmult = nt
         goto 20

*-----------------------------------------------------------------------
*     analog sampling of delayed neutrons.
*-----------------------------------------------------------------------

   10    nt = tn+rang()
         cmult = nt
         if(jq.ne.0) nt = 1
         np = 100

*-----------------------------------------------------------------------
*     sample prompt and delayed neutrons.
*-----------------------------------------------------------------------

   20 do 70 i = 1, nt

         if(i.gt.np) goto 40
         if(np.ne.100) goto 30
         if(rang()*tn.le.dn) goto 40

*-----------------------------------------------------------------------
*     get prompt neutrons.
*-----------------------------------------------------------------------

   30    id = jxs(11,iex)
         kd = nint(xss(jxs(10,iex)+ixre-1))
         call xstcas(i,1,q,ia,ka,id,kd)

         colout(3,i) = 0.

         goto 70

*-----------------------------------------------------------------------
*     sample delayed neutrons.
*-----------------------------------------------------------------------

   40    ns = jxs(25,iex)
         rn = rang()
         j = 0
         pg = 0.

      do 50 ix = 1, nxs(8,iex)
         pg = pg+acefcn(ns+1,erg,j)
         if(rn.le.pg) goto 60
   50    ns = ns+1+j
         ix = nxs(8,iex)
         ns = ns-1-j
   60    id = jxs(27,iex)
         kd = nint(xss(jxs(26,iex)+ix-1))
         call xstcas(i,1,q,0,0,id,kd)

         colout(3,i) = -log(rang())/xss(ns)

   70    continue

*-----------------------------------------------------------------------
*     store biased delayed neutron weight in 1st prompt neutron time.
*-----------------------------------------------------------------------

         if(np.ne.100) colout(3,1) = -wgt*dn/db

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine xstsab(ireg,imat)
*                                                                      *
*        calculate exit energy and direction for an s(a,b) collision.  *
*        Last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for  REDUCTION_COUNTER
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2

      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
C MATSUDA 2022.12.12  for T-Point
      common /clusttp/ nixcos, ixcoss(nnn)
!$OMP THREADPRIVATE(/clusttp/)

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)


      common /sabad/  isaba
      common /sansd/  isans

*-----------------------------------------------------------------------
*     select elastic or inelastic scattering.
*-----------------------------------------------------------------------

         ipsc = 9
         ntyn = 0
         ixcos = 0
         mixed_incoherent = 0
         if (isans .eq. 0) then
           ! Turn off SANS treatment (default)
           if(jxs(4,iet).eq.0) goto 10
           r = rang()
           if(r*rtc(7,iet).gt.rtc(8,iet)) then
             if(r*rtc(7,iet).lt.(rtc(8,iet)+rtc(9,iet)))
     &           mixed_incoherent = 1
             goto 60
           endif
         else
           ITSANS = jxs(10,iet)
           if( (jxs(4,iet).ne.0).or.(ITSANS.ne.0) ) then
             ! from xstneu:  rtc(7,iet,)  = inelast + elast + SANS + elast2 xsec for Sab
             !               rtc(8,iet,)  = inelast                         xsec for Sab
             !               rtc(9,iet,)  = elast2                          xsec for Sab
             !               rtc(10,iet,) = SANS                            xsec for Sab
             r = rang()
             if( r*rtc(7,iet).gt.(rtc(8,iet) + rtc(10,iet))) then
               !    -> elast
               if(r*rtc(7,iet).lt.(rtc(8,iet)+rtc(9,iet)+rtc(10,iet)))
     &           mixed_incoherent = 1
               goto 60
             else if (r*rtc(7,iet).gt.rtc(8,iet))then
               !    -> SANS
               SANS_MODEL = xss(ITSANS)
               if (SANS_MODEL .eq. 1) then
                 ! Analytical model
                 ntyn = 3
                 SANS_A1 = xss(ITSANS+2)
                 SANS_b1 = xss(ITSANS+3)
                 SANS_A2 = xss(ITSANS+4)
                 SANS_b2 = xss(ITSANS+5)
                 SANS_Q0 = xss(ITSANS+6)
                 SANS_k0 = sqrt(erg/2.072e-9) !Angst^-1
                 !
                 SANS_s1 = SANS_A1*      SANS_Q0**(SANS_b1+2.0)
     &                                   /(SANS_b1+2.0)
                 SANS_s2 = SANS_A2*(2.0*SANS_k0)**(SANS_b2+2.0)
     &                                   /(SANS_b2+2.0)
                 SANS_s3 = SANS_A2*      SANS_Q0**(SANS_b2+2.0)
     &                                   /(SANS_b2+2.0)
                 !
                 ! Sample change in wavevector
                 r = rang()
                 if (r.lt.SANS_s1/(SANS_s1+SANS_s2-SANS_s3)) then
                  SANS_q =  (r*(SANS_s1+SANS_s2-SANS_s3)
     &             *(SANS_b1+2.0)/SANS_A1)
     &             **(1.0/(SANS_b1+2.0))
                 else
                   SANS_q =
     &             ((r*(SANS_s1+SANS_s2-SANS_s3)-SANS_s1+SANS_s3)
     &             *(SANS_b2+2.0)/SANS_A2)
     &             **(1.0/(SANS_b2+2.0))
                 endif
                 ! Compute cosine
                 cs = 1.0 - 0.5*(SANS_q/SANS_k0)**2.0
                 goto 130
               else
                 ! -> error, other SANS models not implemented
                 goto 140
               endif
             endif
           endif
           ! -> inelast
      endif

*-----------------------------------------------------------------------
*     get the sampling parameters and energy for the inelastic case.
*-----------------------------------------------------------------------

   10    if(nxs(7,iet).ne.0) goto 20
*-----------------------------------------------------------------------
*     IFENG = 0: Discrete distribution
*-----------------------------------------------------------------------
         kr = rang()*nxs(4,iet)
         goto 50

   20    if(nxs(7,iet).gt.1) goto 55
*-----------------------------------------------------------------------
*     IFENG = 1: Skewed discrete distribution
*     weights: 1, 4, 10, 10, ..., 10, 10, 4, 1
*-----------------------------------------------------------------------
         r = rang()*(nxs(4,iet)-3)
         if(r.lt.1.) goto 30
         kr = r+1.
         goto 50

   30    if(r.lt..5) goto 40
         kr = nxs(4,iet)-2
         if(r.lt..6) kr = kr+1
         goto 50

   40    kr = 1
         if(r.lt..1) kr = 0

   50    ll = nxs(3,iet)
         ni = nxs(4,iet)*(ll+2)
         kx = (ktc(1,iet)-1)*ni+jxs(3,iet)+kr*(ll+2)
         erg = xss(kx)+rtc(1,iet)*(xss(kx+ni)-xss(kx))

*-----------------------------------------------------------------------

      if( isaba .ne. 0 ) then
         ! Smear outgoing energy in logarithmic scale.
         ! Improves the spectrum when discrete files are used,
         ! but this should not be applied to
         ! continuous energy files.

         erg0 = erg
         kr2 = max(0,kr-1)
         kx2 = (ktc(1,iet)-1)*ni+jxs(3,iet)+kr2*(ll+2)
         erg2 = xss(kx2)+rtc(1,iet)*(xss(kx2+ni)-xss(kx2))
         kr3 = min(kr+1,nxs(4,iet))
         kx3 = (ktc(1,iet)-1)*ni+jxs(3,iet)+kr3*(ll+2)
         erg3 = xss(kx3)+rtc(1,iet)*(xss(kx3+ni)-xss(kx3))

         erg4 = ( log(erg0) + log(erg2) ) / 2.0
         erg5 = ( log(erg3) + log(erg0) ) / 2.0
         erg = exp( erg4 + rang() * ( erg5 - erg4 ) )

      end if

*-----------------------------------------------------------------------

         if(erg.lt.1.e-11) erg = 1.e-11
         kx = kx+1
         ri = rtc(1,iet)
         ntyn = 1
         goto (140,140,70,80) nxs(2,iet) + 1

*-----------------------------------------------------------------------
*     IFENG = 2: Continuous distribution
*     JIMD 2021/10/11
*     Documentation:
*       - NJOY2016 Manual
*         https://github.com/njoy/NJOY2016-manual
*       - ACE Format Manual
*         https://github.com/NuclearData/ACEFormat
*       - OpenMC S(a,b) physics manual and implementation
*         https://docs.openmc.org/en/stable/methods/neutron_physics.html#sab-tables
*-----------------------------------------------------------------------
55       if(nxs(7,iet).gt.2) goto 140 ! IFENG > 2 not implemented
         ri   = rtc(1,iet)       ! interpolation factor for energy table
         ie   = ktc(1,iet)       ! index of XS table
         ITIE = jxs(1,iet)       ! Index of inelastic energy table
         ITXE = jxs(3,iet)       ! Index of inelastic energy/angle distributions table
         ne   = nint(xss(ITIE))  ! number of outgoing energies
         NIEB = nxs(4,iet)       ! number of incident energies
         NIL  = nxs(3,iet)       ! inelastic dimensioning parameter
         ergn = erg              ! incident energy
         erg1 = xss(ITIE+ie)     ! previous point in the energy table
         erg2 = xss(ITIE+ie+1)   ! next point in the energy table

         loc1 = ITIE + nint(xss(ITXE-1+ie))
         nn1  =        nint(xss(ITXE-1+ie+ne))   ! number of outgoing energy bins

         loc2 = ITIE + nint(xss(ITXE-1+ie+1))
         nn2  =        nint(xss(ITXE-1+ie+ne+1)) ! number of outgoing energy bins

         !
         ! Interpolate the maximum and minimum energies of the distribution
         !
         Emin1 = xss(loc1)
         Emin2 = xss(loc2)
         Emax1 = xss(loc1 +(nn1-1)*(NIL+2))
         Emax2 = xss(loc2 +(nn2-1)*(NIL+2))
         Emin  = Emin1 + ri*(Emin2-Emin1)
         Emax  = Emax1 + ri*(Emax2-Emax1)
         !
         ! Choose outgoing energy distribution
         !
         if( rang().ge.ri) then
            locn  = loc1
            nn    = nn1
         else
            locn  = loc2
            nn    = nn2
         endif
         !
         ! Binary search of outgoing energy bin
         !
         r = rang()
         ic = 1
         ib = nn
         do
            if ((ib-ic).eq.1) exit
            kx = (ic+ib)/2
            loc = locn + (kx-1)*(NIL+2)
            if ( r.lt.xss(loc+2)) then
               ib = kx
            else
               ic = kx
            endif
         enddo
         kx = ib
         loc = locn + (kx-1)*(NIL+2)
         !
         ! Sample inside of bin using method in
         ! Edward, Rathkopf and Smidt (UCRL-JC-104791)
         !
         x  = xss(loc)
         y  = xss(loc+1)
         xl = xss(loc  -NIL-2)
         yl = xss(loc+1-NIL-2)
         r1 = rang()
         r2 = rang()
         if ( r2*(yl+y) < r1*(y-yl) + yl ) then
            erg = xl + (x-xl)*r1
         else
            erg = xl + (x-xl)*(1.0-r1)
         endif
         kx = loc+3-(NIL+2)
         ri = (erg - xl)/(x-xl)
         rtc(2,iet) = ri ! save outgoing energy interpolation factor
         !
         ! Scale energy exchange and apply to incident energy
         !
         if ( locn.eq.loc1 ) then
            if ( erg.lt.erg1) then
              erg = ergn + (erg-erg1)/(erg1 - Emin1)*(ergn - Emin)
            else
              erg = ergn + (erg-erg1)/(Emax1 - erg1)*(Emax - ergn)
            end if
         else
            if ( erg.lt.erg2) then
              erg = ergn + (erg-erg2)/(erg2 - Emin2)*(ergn - Emin)
            else
              erg = ergn + (erg-erg2)/(Emax2 - erg2)*(Emax - ergn)
            end if
         endif
         !
         ! Indices for angular distribution
         !
         ntyn = 1
         ll = NIL - 2
         ni = (NIL+2)
         if(erg.lt.1.e-11) erg = 1.e-11
         goto 80

*-----------------------------------------------------------------------
*     set up the sampling parameters for the elastic case.
*-----------------------------------------------------------------------

   60    ntyn = 2

         if (nxs(5,iet).eq.5.and.mixed_incoherent.eq.1) then
           kx = jxs(9,iet)+(ktc(3,iet)-1)
     &        * (abs(nxs(8,iet))+1)
           ri = rtc(5,iet)
           ni = nxs(8,iet)+1
           ll = nxs(8,iet)
         else
           kx = jxs(6,iet)+(ktc(2,iet)-1)
     &        * (abs(nxs(6,iet))+1)
           ri = rtc(4,iet)
           ni = nxs(6,iet)+1
           ll = nxs(6,iet)
         endif
         if (nxs(5,iet).eq.2) goto 70
         if (nxs(5,iet).eq.3) goto 80
         if (nxs(5,iet).eq.4) goto 90
         if (nxs(5,iet).eq.5.and.mixed_incoherent.eq.1) goto 80
         if (nxs(5,iet).eq.5) goto 90
         goto 140

*-----------------------------------------------------------------------
* >>>>>  equally-probable angle bins.
*-----------------------------------------------------------------------
         ! IDPNI/IDPNC = 2
   70    if(rang().le.ri) kx = kx+ni
         f = ll*rang()+1.
         j = kx+f
         cs = xss(j)+(xss(j-1)-xss(j))*(f-aint(f))
         ixcos = -kx
         goto 130

*-----------------------------------------------------------------------
* >>>>>  equally-probable discrete angles.
*-----------------------------------------------------------------------
         ! IDPNI/IDPNC = 3
   80    j = kx+(ll+1)*rang()
         cs = xss(j)+ri*(xss(j+ni)-xss(j))
         ixcos = kx

         ! JIMD 2021/10/11
         ! Smear scattering cosine using the Hendricks-Prael
         ! algorithm (LA-11952). Always apply if IFENG=2.
         ! Modified to preserve the average cosine in the bin.
         ! For IFENG=2 the interpolation is done between
         ! outgoing energies

         if((nxs(7,iet).gt.1).or.(isaba.eq.1)) then
           cs0 = cs
           ! left width
           if( j .le. kx ) then
              ! first panel
              cs1 = -1.0
           else
              j11 = j-1
              cs1=xss(j11)+ri*(xss(j11+ni)-xss(j11))
           end if
           csm1 = ( cs0 - cs1 ) / 2.0

           ! right width
           if( j .ge. kx + ll ) then
              ! last panel
              cs2 = 1.0
           else
              j22 = j+1
              cs2=xss(j22)+ri*(xss(j22+ni)-xss(j22))
           end if
           csm2 = ( cs2 - cs0 ) / 2.0

           csm = min(csm1, csm2)
           cs = cs0 + 2.0 * csm * (rang() - 0.5)
         end if

      go to 130

*-----------------------------------------------------------------------
* >>>>>  exact treatment of coherent elastic scattering
*-----------------------------------------------------------------------
         ! IDPNC = 4
         ! Sample Bragg edge between the first and the incident energy

   90    ic = jxs(5,iet)-1
         ib = jxs(5,iet)+ktc(2,iet)-1
         if(ib.eq.ic+1) goto 120
         pr = xss(ib)*rang()

  100    if(ib-ic.eq.1) goto 120
         ih = (ic+ib)/2
         if(pr.lt.xss(ih)) goto 110
         ic = ih
         goto 100

  110    ib = ih
         goto 100

  120    cs = 1.-2.*xss(ib-nint(xss(jxs(4,iet))))/erg

*-----------------------------------------------------------------------
*     sample the azimuthal angle.
*     and booking
*-----------------------------------------------------------------------

  130 continue

                     call dtcos(cs,uold,uuu,0,irdm)

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 2
                     nixcos = nixcos + 1  ! MATSUDA 2022.12.12  for T-Point

                     jclusts(0,nclsts) = ipsc
                     jclusts(1,nclsts) = 0
                     jclusts(2,nclsts) = 1
                     jclusts(3,nclsts) = 2
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = 0
                     jclusts(6,nclsts) = 1
                     jclusts(7,nclsts) = 2112
                     jclusts(8,nclsts) = 0

                     rms = gpt(1)
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
                     qclusts(10,nclsts) = 0.0
                     qclusts(11,nclsts) = 0.0
                     qclusts(12,nclsts) = 0.0
                     ixcoss(nixcos) = ixcos  ! MATSUDA 2022.12.12  for T-Point

               if( ntyn .eq. 2 ) then

                     kcoll = 3
                     aevts(iaevt+17,ireg) =
     &               aevts(iaevt+17,ireg) + wgt
                     bevts(ibevt+17,imat) =
     &               bevts(ibevt+17,imat) + wgt
               else

                     kcoll = 4

                  if( erg .gt. eg0 ) then
                     aevts(iaevt+22,ireg) =
     &               aevts(iaevt+22,ireg) + wgt * ( erg - eg0 )
                     bevts(ibevt+22,imat) =
     &               bevts(ibevt+22,imat) + wgt * ( erg - eg0 )
                  else
                     aevts(iaevt+23,ireg) =
     &               aevts(iaevt+23,ireg) + wgt * ( eg0 - erg )
                     bevts(ibevt+23,imat) =
     &               bevts(ibevt+23,imat) + wgt * ( eg0 - erg )
                  end if
                     aevts(iaevt+24,ireg) =
     &               aevts(iaevt+24,ireg) + wgt
                     bevts(ibevt+24,imat) =
     &               bevts(ibevt+24,imat) + wgt
               end if

      return

*-----------------------------------------------------------------------

  140    write(ErrCha,*) 'Error : in xstsab, '//
     &       'inappropriate distribution required.'
           ErrID = 'L:4623/R:xstsab/F:ggm05.f' !E04_024_001
           call ErrWrite(ErrID,ErrCha)

            call parastop( 428 )

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine isos(d,l)
*                                                                      *
*        sample a direction d isotropically.                           *
*        Last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

      dimension d(3)

*-----------------------------------------------------------------------

   10    t1 = 2.*rang()-1.
         t2 = 2.*rang()-1.
         r = t1**2+t2**2
         if(r.gt.1.) goto 10

         d(1) = 2.*r-1.
         t3 = sqrt((1.-d(1)**2)/r)
         d(2) = t1*t3
         d(3) = t2*t3


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function acecs6(ii,id,iw,jc,r,cs)
*                                                                      *
*        sample energy acecs6 given angle cs, law 67 (endf/b-vi law 7).*
*        ii=0 acecs6 called from xstcas (transport): save id,iw,jc,r   *
*             and random numbers generated in sampling energy;         *
*        ii=1 acecs6 called from calcps (next-event estimator):        *
*             use id,iw,jc,r, and random numbers saved from transport. *
*        Last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      dimension el(2),eh(2),lx(2)

*-----------------------------------------------------------------------

         t = 0.
         if(ii.ne.0) goto 10
         id0 = id
         iw0 = iw
         ic0 = jc
         rr0 = r

c%%<2010/03/19:kosako: for poor data at threshold energy of law67
         ixloop = 0
c%%>
*-----------------------------------------------------------------------
*     go to appropriate table for incident energy
*-----------------------------------------------------------------------

   10    jw = iw0+2*nint(xss(iw0))+1
         ne = nint(xss(jw))
         lx(1) = id0+nint(xss(jw+ne+ic0))-1
         lx(2) = 0
         ix = 1
         ir = 0
         if(rr0.eq.0.) goto 20
         lx(2) = id0+nint(xss(jw+ne+ic0+1))-1
         if(ii.eq.0) rnb(1) = rang()
         ir = 1
         if(rnb(1).le.rr0) ix = 2

   20 do 110 m = 1, 2

         if(lx(m).eq.0) goto 110
         le = lx(m)
         mu = nint(xss(le))
         nm = nint(xss(le+1))

*-----------------------------------------------------------------------
*     find appropriate table for sampled cosine
*-----------------------------------------------------------------------

      do 30 iq = 1, nm - 1
   30    if(xss(le+iq+2).ge.cs) goto 40

   40    if(mu.eq.1) goto 50
         ir = ir+1
         if(ii.eq.0) rnb(ir) = rang()
         if(rnb(ir).le.(cs-xss(le+iq+1))/(xss(le+iq+2)-
     &    xss(le+iq+1))) iq = iq+1

*-----------------------------------------------------------------------
*     sample from tabulated energy distribution
*-----------------------------------------------------------------------

   50    lb = id0+nint(xss(le+nm+iq+1))
         jj = nint(xss(lb-1))
         np = nint(xss(lb))
         el(m) = xss(lb+1)
         eh(m) = xss(lb+np)
         if(m.ne.ix) goto 110
         ir = ir+1
         if(ii.eq.0) rnb(ir) = rang()
         ic = lb+2*np+1
         ib = lb+3*np
         call bnsrch(rnb(ir),ic,ib,ig)
   80    l = ic-2*np
         fa = xss(l+np)
         ea = xss(l)
         if(jj.eq.1) goto 90
         bb = (xss(l+np+1)-fa)/(xss(l+1)-ea)
         if(bb.eq.0.) goto 90
         t = ea+(sqrt(max(zero,fa**2+2.*bb*(rnb(ir)-xss(ic))))-fa)/bb
         goto 100

*-----------------------------------------------------------------------
c%%<2010/03/19:kosako: for poor data at threshold energy of law67
c%%90 t=ea+(rnb(ir)-xss(ic))/fa
*-----------------------------------------------------------------------
c
c%%>
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
c%%<2010/03/19:kosako: for poor data at threshold energy of law67
c%%<2011/11/14:kosako: revised
c%%90 t=ea+(rnb(ir)-xss(ic))/fa
*-----------------------------------------------------------------------

   90 continue
      if (fa.eq.0.) then
         ixloop=ixloop+1
         if (ixloop.eq.1) then
            if (cs.lt.xss(le+iq+1)) then
               if (iq.gt.1) then
                  iq=iq-1
                  go to 50
               endif
            else
               if (iq.lt.nm) then
                  iq=iq+1
                  go to 50
               endif
            endif
            if (m.eq.1) then
               go to 110
            endif
         elseif (ixloop.eq.2) then
            if (m.eq.1) then
               go to 110
            elseif (m.eq.2.and.t.gt.0.) then
               go to 110
            endif
         endif
         if (ixloop.gt.0) then

           write(*,*) '# Warning: in acecs6, '//
     &        'energy sampling of law67 was illegal.'

            le=nint(xss(le))
            nm=nint(xss(le+1))
            lb=id0+nint(xss(le+nm+iq+1))
            np=nint(xss(lb))
            write(*,*) ' le,mu,nm,iq,lb,np=',le,mu,nm,iq,lb,np
            do i=1,np
               write(*,'(i4,1p,3e13.5)') i,xss(lb+i),xss(lb+np+i),
     &                                   xss(lb+2*np+i)
            enddo
         endif
      else
         t=ea+(rnb(ir)-xss(ic))/fa
      endif
c%%>
*-----------------------------------------------------------------------

  100    acecs6 = t

  110 continue

*-----------------------------------------------------------------------
*     use scaled interpolation between energies
*-----------------------------------------------------------------------

         if(rr0.eq.0.) return
         t1 = el(1)+rr0*(el(2)-el(1))
         t2 = eh(1)+rr0*(eh(2)-eh(1))
         acecs6 = t1+(t-el(ix))*(t2-t1)/(eh(ix)-el(ix))

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine bnsrch(tv,il,ih,ig)
*                                                                      *
*        find the value tv in the array slice xss(il..ih)              *
*                                                                      *
*     Preconditions:                                                   *
*        tv is a real value within the array                           *
*        xss(il..ih) is an numerically ascending ordered array list    *
*        il is the first word of the list (index low)                  *
*        ih is the last word of the list (index high)                  *
*                                                                      *
*        tv and all xss(i) are read only variables                     *
*                                                                      *
*     Postconditions:                                                  *
*        il/ih bracket the test value tv such that ih-il.eq.1 and      *
*          xss(il) < tv < xss(ih)                                      *
*                                                                      *
*        Occasionally a value outside the array will be called for     *
*        and a reasonable value, the extreme edge of the array,        *
*        will be returned along with an non-zero error flag ig         *
*          ig =  0 is normal state - value found in array slice        *
*          ig = -1 is warning state - value below first array value    *
*          ig =  1 is warning state - value above last array value     *
*          ig =  2 is warning state - incoming indices equal tv=xss(il)*
*          ig = -3 is warning state - incoming indices equal tv<xss(il)*
*          ig =  3 is warning state - incoming indices equal tv>xss(il)*
*        also il = ih = edge value for exit warning condition          *
*                                                                      *
*        il, ih and ig are modified return values                      *
*                                                                      *
*        Last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------
*     warning state - value below table
*-----------------------------------------------------------------------

      if(tv.lt.xss(il)) then
         ig = -1
         ih = il

*-----------------------------------------------------------------------
*     warning state - value above table
*-----------------------------------------------------------------------

      else if(tv.gt.xss(ih)) then
         ig = 1
         il = ih

*-----------------------------------------------------------------------
*     warning state - incoming indices are equal
*-----------------------------------------------------------------------

      else if(il.eq.ih) then

         if(tv.eq.xss(il)) then
            ig = 2
         else if(tv.lt.xss(il)) then
            ig = -3
         else
            ig = 3
         endif

*-----------------------------------------------------------------------
*     else handle normal value within table
*-----------------------------------------------------------------------

      else
         ig = 0
 10      if(ih-il.eq.1) return
         im = (il+ih)/2
         if(tv.lt.xss(im)) goto 20
         il = im
         goto 10
 20      ih = im
         goto 10

*-----------------------------------------------------------------------
*     end of cases
*-----------------------------------------------------------------------
      endif

*-----------------------------------------------------------------------

      return
      end


