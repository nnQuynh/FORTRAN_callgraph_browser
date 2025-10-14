************************************************************************
*                                                                      *
      subroutine ncasc(ipos,ityp,ktyp,eein,mmas,mchg,bmax0)
*                                                                      *
*                                                                      *
*       control routine of intra-nuclear cascade                       *
*            modified by K.Niita and S. Hashimoto on 2011/08/02        *
*       last modified by T.Ogawa                  on 2015/02/25        *
*                                                                      *
*        call subroutine : bertin, jamin, jqmd                         *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       ipos    : =0, call from nreac.f ,=1 from tally                 *
*                                                                      *
*       ityp    : particle type of projectile                          *
*       ktyp    : kf code of projectile                                *
*       eein    : energy of projectile (MeV)                           *
*       mmas    : mass of target                                       *
*       mchg    : charge of target                                     *
*       bmax0    : max impact parameter                                *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclst   : total number of out going particles and nuclei      *
*        iclust(nclst)                                                 *
*        jclust(i,nclst)                                               *
*        qclust(i,nclst)                                               *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg
C for USE_MOD_COUNTER
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
C for REDUCTION_COUNTER
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'

      parameter ( rpmass = 938.27, rnmass = 939.58 )
      parameter ( pi = 3.1415926535898d0 )

*-----------------------------------------------------------------------

      common /cascid/ jcasc
!$OMP THREADPRIVATE(/cascid/)
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /cincl/  inclg, inclv
      common /cincelf/incelf
      common /ceinc/  einclmin, einclmax, eielfmin, eielfmax
      common /cidwba/ idwba
      common /cdwba/  bdwba, idwbahit
!$OMP THREADPRIVATE(/cdwba/)

*-----------------------------------------------------------------------


      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

      common /othjmp/ nojmp(2000,2), nojmn

      common /engch/  ejamnu, ejampi, eisobar, eqmdnu, eqmdmn, ejamqmd
      common /qmdflg/ irqmd

      common /vriab1/ b, llnow, ntnow
!$OMP THREADPRIVATE(/vriab1/)
      common /rqmdaux/ ibsv
!$OMP THREADPRIVATE(/rqmdaux/)
      common /coln01/ iccoll
!$OMP THREADPRIVATE(/coln01/)
      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
!$OMP THREADPRIVATE(/const1/)
      common /qmdswt/ ijmq, imosc, prat
!$OMP THREADPRIVATE(/qmdswt/)

! T.Sato 2015/03/09, to avoid forced stop due to parastop
      common /errorcom/ncascerr
!$OMP THREADPRIVATE(/errorcom/)
      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
!$OMP THREADPRIVATE(/const0/)
      common /framtr/ betafr(0:2), gammfr(0:2)
!$OMP THREADPRIVATE(/framtr/)

      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

*-----------------------------------------------------------------------
      data rmnuc / 938.95 /

! T.Sato 2015/03/09, initialization for error index
      ncascerr = 0

*-----------------------------------------------------------------------
*        zero set
*-----------------------------------------------------------------------
                  kcoll = 0

               do i = 1, 4
                  kdecay(i) = 0
               end do

*-----------------------------------------------------------------------
*        choice of changing energy from Bertini to JAM
*-----------------------------------------------------------------------

                  cenucl = ejamnu
                  cepion = ejampi
                  ceqmdn = eqmdnu

*-----------------------------------------------------------------------
*        additional option of discrete levels using DWBA calculation
*-----------------------------------------------------------------------

         if ( idwba .ne. 0 .and. idwbahit .eq. 1 ) then

            call dwbain(ityp,eein,mmas,mchg)

         else if ( idwba .eq. 0 ) then

            idwbahit = 0

         end if

*-----------------------------------------------------------------------
*        jcasc
*-----------------------------------------------------------------------
*           nucleus -> QMD,  Kaon and A < 6 -> JAM
*-----------------------------------------------------------------------

         if( ( incelf .ge. 1 .and. ityp .le. 2 .and.
     &        eein .ge. eielfmin .and. eein .lt. eielfmax ) .or.
     &        ( incelf .eq. 1 .and. ityp .eq. 18 .and.
     &        eein / dble(ibryf(ityp,ktyp)) .ge. eielfmin .and.
     &        eein / dble(ibryf(ityp,ktyp)) .lt. eielfmax ) ) then

            jcasc = 6

         else if( ( inclg .ge. 1 .and. ityp .le. 5 .and.
     &           eein .ge. einclmin .and. eein .lt. einclmax ) .or.
     &           ( inclg .eq. 1 .and.
     &           ityp .ge. 15 .and. ityp .le. 18 .and.
     &           eein / dble(ibryf(ityp,ktyp)) .ge. einclmin .and.
     &           eein / dble(ibryf(ityp,ktyp)) .lt. einclmax ) ) then

            jcasc = 5

C S.H. When A < 5 in pion-induced reactions, INCL -> JAM. (2014.8.14)
            if ( ityp .ge. 3 .and. ityp .le. 5 .and. mmas .lt. 5 ) then
               jcasc = 3
            end if

         else if( ityp .ge. 15 .or.
     &           ( ityp .le. 2 .and. eein .ge. eqmdnu .and.
     &           eein .lt. ejamnu ) ) then

            jcasc = 4

         else if( ( ityp .le. 2 .and. eein .gt. ejamnu ) .or.
     &           ( ( ityp .eq. 3 .or. ityp .eq. 5 ) .and.
     &           eein .ge. ejampi ) .or.
     &           ityp .eq. 4 .or.
     &           ityp .ge. 8 .or. mmas .lt. 6 ) then

            jcasc = 3

         else if( icasc .eq. 2 .and.
     &           ityp .le. 2 .and. eein .le. eisobar ) then

            jcasc = 2

         else

            jcasc = 1

         end if

         if ( idwbahit .eq. 1  ) jcasc = 0

*-----------------------------------------------------------------------
*     jcasc = 6 : INC-ELF
*-----------------------------------------------------------------------

      if( jcasc .eq. 6 ) then
            if( ipos .eq. 0 ) call cputime(21)
*-----------------------------------------------------------------------

           call incelfin(ityp,eein,mmas,mchg)

*-----------------------------------------------------------------------

            if( ipos .eq. 0 ) then

                  rnint(120+ityp) = rnint(120+ityp) + 1.0

               if( nclst .ge. 0 ) then
                  rnintr(120+ityp) = rnintr(120+ityp) + 1.0
               end if

                  rncnt(21) = rncnt(21) + 1.0

                  call cputime(21)

            end if

*-----------------------------------------------------------------------
*     jcasc = 5 : INCL
*-----------------------------------------------------------------------

      else if( jcasc .eq. 5 ) then

            if( ipos .eq. 0 ) call cputime(20)

*-----------------------------------------------------------------------

            call inclin(ityp,eein,mmas,mchg)

*-----------------------------------------------------------------------

            if( ipos .eq. 0 ) then

                  rnint(100+ityp) = rnint(100+ityp) + 1.0

               if( nclst .ge. 0 ) then
                  rnintr(100+ityp) = rnintr(100+ityp) + 1.0
               end if

                  rncnt(20) = rncnt(20) + 1.0
                  call cputime(20)

            end if

*-----------------------------------------------------------------------
*     jcasc = 4 : QMD
*-----------------------------------------------------------------------

      else if( jcasc .eq. 4 ) then


         ibsv = 0 ! >1: retry run
         ijmq = 0 ! ijmq = 0 for JQMD, ijmq = 1 for JAMQMD

*-------- target - projectile fusion mode --------------------------------




*-----------------------------------------------------------------------

 100     if( eein / dble(ibryf(ityp,ktyp)) .gt. eqmdmn ) then

            if( eein / dble(ibryf(ityp,ktyp)) .le. ejamqmd ) then

               if( ipos .eq. 0 ) call cputime(18)
               if(irqmd .eq. 0) then

                  call jqmdin(ityp,ktyp,eein,mmas,mchg,bmax0) ! legacy JQMD mode

               else

                  call jqmdinR(ityp,ktyp,eein,mmas,mchg,bmax0) ! JQMD Ver.2 mode

               endif

               if( ipos .eq. 0 ) then

                  rnint(60+ityp) = rnint(60+ityp) + 1.0

                  if( nclst .ge. 0 ) then
                     rnintr(60+ityp) = rnintr(60+ityp) + 1.0
                  end if

                  rncnt(18) = rncnt(18) + 1.0
                  call cputime(18)

               endif

            else

               ijmq  = 1
               imosc = 0
               prat  = 0.d0

               if( ipos .eq. 0 ) call cputime(19)
               call jamin(1,ityp,ktyp,eein,mmas,mchg,bmax0)

! T.Sato 2015/03/09, to avoid forced stop due to parastop
               if(ncascerr.ne.0) then
                nclst = -2
                ncascerr = 0
               endif

               if( ipos .eq. 0 ) then

                  rnint(80+ityp) = rnint(80+ityp) + 1.0

                  if( nclst .ge. 0 ) then
                     rnintr(80+ityp) = rnintr(80+ityp) + 1.0
                  end if

                  rncnt(19) = rncnt(19) + 1.0
                  call cputime(19)

               endif

            end if


         else

                  nclst = 0

         end if

*-----------------------------------------------------------------------
*        check of the total energy conservation 2008/06/20
*-----------------------------------------------------------------------

  200    if( nclst .le. 0 .and. irqmd .eq. 0) return

            ibrym = ibryf(ityp,ktyp)
            projm = rmtyp(ityp,ktyp)

            if( ityp .lt. 3 ) then

                  pima = 0.d0

            else if( ityp .lt. 15 ) then

               if( ibrym .eq. 0 ) then

                  pima = projm

               else if( ibrym .lt. 0 ) then

                  pima = projm + rmnuc

               else if( ibrym .gt. 0 ) then

                  pima = projm - rmnuc

               end if

           else

                  mpms = ibryf(ityp,ktyp)
                  mpch = ichgf(ityp,ktyp)
                  pima = -bindeg(mpch,mpms-mpch)

           end if

               eout = 0.d0
               ebin = 0.d0
               exct = 0.d0

            do i = 1, nclst

                  exct = exct + qclust(6,i)
                  eout = eout + qclust(7,i)

               if( iclust(i) .eq. 0 ) then

                  ebin = ebin + bindeg(jclust(1,i),jclust(2,i))

               else if( iclust(i) .gt. 2 ) then

                  eout = eout + qclust(5,i) * 1000.0

               end if

            end do

               eres = eein + pima - eout
     &              + ebin - bindeg(mchg,mmas-mchg)

           if(irqmd .eq. 0 .and. ijmq .eq. 0) then  ! 2014/6/13 ogawa   Use the scheme same as ever.

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

           else  ! 2014/6/13 ogawa  In JQMD Ver2, Retry with same b unless reaction is elastic.

            if( nclst .le. -2) then ! strange nucleus. Give up

             return

            elseif( nclst .lt. 0) then ! Too much negative excitation. Retry with same b

              ibsv = 1    ! retry with same b
              goto 100

            elseif( eres .lt. 0.d0 .and. nclst .gt. 2) then ! nclst ==2 is exempted in case (eout > eein)

           if(ijmq .eq. 1 .and. abs(eres) .lt. eein * 1.d-1 .and. imosc
     &      .lt. 10) then ! 2016/11/04 recover central collision events (energy conservation is often unconserved)

            rmta = dble(mstapr) * 938.27 + dble(massta-mstapr) * 939.58
     &       - bindeg(mstapr,massta-mstapr)
            rmpr = dble(msprpr) * 938.27 + dble(masspr-msprpr) * 939.58
     &       - bindeg(msprpr,masspr-msprpr)

            if( sqrt(rmpr**2 + rmta**2 + 2.d0*rmta*(rmpr+eein)) - rmta
     &       - rmpr .gt. 0.2d6 .or. nclst .gt. 100) then ! Apply momentum scaling when energy is high ( sqrt(s)/(m1+m2) > 7 )

            imosc = imosc + 1
            if(prat .eq. 0.d0) then
             prat  = abs(eres) / (sum(qclust(1:3,1:nclst)) * 1.d3)
            endif
            call cldistR

               ifrm = 0

               do i = 1, nclst

                     px = qclust(1,i)
                     py = qclust(2,i)
                     pz = qclust(3,i)
                     et = qclust(4,i)
                     rm = qclust(5,i)

                     gamm = gammfr(ifrm)
                     beta = betafr(ifrm)

                     pz = pz * gamm - beta * gamm * et
                     et = sqrt( px**2 + py**2 + pz**2 + rm**2 )

                     qclust(3,i)  = pz
                     qclust(4,i)  = et
                     qclust(7,i)  = ( et - rm ) * 1000.

               end do


*-----------------------------------------------------------------------
*           x-y rotation by random
*-----------------------------------------------------------------------

                  theta = 2.0d0 * pi * rn(0)

               do i = 1, nclst

                  px = qclust(1,i)
                  py = qclust(2,i)

                  qclust(1,i)  = px * cos( theta ) - py * sin( theta )
                  qclust(2,i)  = px * sin( theta ) + py * cos( theta )

               end do

            goto 200
           endif
           endif

              ibsv = 1    ! retry with same b
              goto 100

            elseif( eres .lt. 0.d0 .and. nclst .le. 2) then !
               if(unirn(dummy) .le. 0.85d0) then
                 ibsv = 1
                 goto 100
               endif

              nclst = -1
              return

            elseif( product(qclust(6,1:nclst),
     &       mask = jclust(6,1:nclst).gt.4) .eq. 0.d0) then

              if(nclst .eq. 2 ) then
               if(unirn(dummy) .le. 0.85d0) then
                 ibsv = 1
                 goto 100
               endif

                nclst = -1
                return
              else
              goto 100
                ibsv = 1
                goto 100
              endif

            elseif( qclust(6,1) .lt. dble(jclust(6,1))*0.5d0 .and.
     &              qclust(6,2) .lt. dble(jclust(6,2))*0.5d0 .and.
     &              nclst .le. 2 .and. ijmq .eq. 1) then ! 2016/11/01 reject if projectile and target pass by with small excitation

                nclst = -1
                return

            else

             if(exct .eq. 0.d0) then
               erat = 0.d0
             else
               erat = eres / exct
             endif

             do i = 1, nclst

               qclust(6,i) = qclust(6,i) * erat

             end do

            end if

           endif

            if(((jclust(6,2) .eq. mmas               .and.
     &           jclust(5,2) .eq. mchg               .and.
     &           jclust(6,1) .eq. ibryf(ityp,ktyp)   .and.
     &           jclust(5,1) .eq. ichgf(ityp,ktyp))  .or.
     &          (jclust(6,1) .eq. mmas               .and.
     &           jclust(5,1) .eq. mchg               .and.
     &           jclust(6,2) .eq. ibryf(ityp,ktyp)   .and.
     &           jclust(5,2) .eq. ichgf(ityp,ktyp))) .and.
     &           irqmd .eq. 1 ) then ! 2014/5/26 Ogawa. Judge whether excitation energy is high. In case of no nucleon removal

             if(jclust(6,1) .gt. 100 ) then
              remin1n = 0.3d0 * dble(jclust(6,1)) ! large nuclei often get excited. Do not check, otherwise spuriously disintegrated
              remin1p = 0.3d0 * dble(jclust(6,1))
              remin1a = 0.3d0 * dble(jclust(6,1))
             elseif(jclust(6,1) .gt. 1 ) then
              call esepa(jclust(1,1), jclust(2,1), remin1n, remin1p,
     &         remin1a)
             else
              remin1n = 0.d0
              remin1p = 0.d0
              remin1a = 0.d0
             endif

             if(jclust(6,2) .gt. 100 ) then
              remin2n = 0.3d0 * dble(jclust(6,2)) ! large nuclei often get excited. Do not check, otherwise spuriously disintegrated
              remin2p = 0.3d0 * dble(jclust(6,2))
              remin2a = 0.3d0 * dble(jclust(6,2))
             elseif(jclust(6,2) .gt. 1 ) then
              call esepa(jclust(1,2), jclust(2,2), remin2n, remin2p,
     &         remin2a)
             else
              remin2n = 0.d0
              remin2p = 0.d0
              remin2a = 0.d0
             endif

             if( qclust(6,1) .le. min(remin1n, remin1p, remin1a) .and.
     &        qclust(6,2) .le. min(remin2n, remin2p, remin2a) ) then

               if(unirn(dummy) .le. 0.85d0) then
                 ibsv = 1
                 goto 100
               endif

               nclst = -1
               return

             elseif( (qclust(6,1) .lt. min(remin1n, remin1p, remin1a)
     &      .and. qclust(6,2) .ge. min(remin2n, remin2p, remin2a)) .or.
     &           (qclust(6,1) .ge. min(remin1n, remin1p, remin1a)
     &      .and. qclust(6,2) .lt. min(remin2n, remin2p, remin2a))) then

              if(unirn(dummy) .le. 0.5d0) then
               nclst = -1
               return
              endif

             endif
            endif


 210        ijmq = 0



*-----------------------------------------------------------------------
*     jcasc = 3 : JAM
*-----------------------------------------------------------------------

      else if( jcasc .eq. 3 ) then

            if( ipos .eq. 0 ) call cputime(17)

*-----------------------------------------------------------------------

            call jamin(0,ityp,ktyp,eein,mmas,mchg,bmax0)

*-----------------------------------------------------------------------

! T.Sato 2015/03/09, to avoid forced stop due to parastop
            if(ncascerr.ne.0) then
             nclst = -2
             ncascerr = 0
            endif
! End of revision by T.Sato

            if( ipos .eq. 0 ) then

               if( ityp .eq. 11 ) then
!$OMP CRITICAL (jmn_crit)
                  if( nojmn .eq. 0 ) then

                        nojmn = nojmn + 1
                        nojmp(nojmn,1) = 1
                        nojmp(nojmn,2) = ktyp

                  else

                     do m = 1, nojmn

                        if( ktyp .eq. nojmp(m,2) ) then

                           nojmp(m,1) = nojmp(m,1) + 1
                           goto 300

                        end if

                     end do

                        nojmn = nojmn + 1
                        nojmp(nojmn,1) = 1
                        nojmp(nojmn,2) = ktyp

  300                continue

                  end if
!$OMP END CRITICAL (jmn_crit)
               end if
               rnint(40+ityp) = rnint(40+ityp) + 1.0

               if( nclst .ge. 0 ) then
                  rnintr(40+ityp) = rnintr(40+ityp) + 1.0
               end if

                  rncnt(17) = rncnt(17) + 1.0

                  call cputime(17)

            end if

*-----------------------------------------------------------------------
*     jcasc = 1 : Bertini
*     jcasc = 2 : Isobar
*-----------------------------------------------------------------------

      else if( jcasc .eq. 1 .or. jcasc .eq. 2 ) then

            if( ipos .eq. 0 .and. jcasc .eq. 1 ) call cputime(15)
            if( ipos .eq. 0 .and. jcasc .eq. 2 ) call cputime(16)

*-----------------------------------------------------------------------

            call bertin(jcasc,ityp,eein,mmas,mchg)

*-----------------------------------------------------------------------

            if( ipos .eq. 0 ) then

               if( jcasc .eq. 1 ) then
                  rnint(0+ityp) = rnint(0+ityp) + 1.0

                  if( nclst .ge. 0 ) then
                     rnintr(0+ityp) = rnintr(0+ityp) + 1.0
                  end if

                  rncnt(15) = rncnt(15) + 1.0

                  call cputime(15)

               else if( jcasc .eq. 2 ) then

                  rnint(20+ityp) = rnint(20+ityp) + 1.0

                  if( nclst .ge. 0 ) then
                     rnintr(20+ityp) = rnintr(20+ityp) + 1.0
                  end if

                  rncnt(16) = rncnt(16) + 1.0

                  call cputime(16)

               end if

            end if

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine nfrgdat(ifg,ein)
*                                                                      *
*                                                                      *
*       fragment production by Frag Data                               *
*       modified by K.Niita on 2003/11/17                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       ifg     : data id                                              *
*       ein     : energy of projectile / nucleon (MeV)                 *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclst   : total number of out going particles and nuclei      *
*        iclust(nclst)                                                 *
*        jclust(i,nclst)                                               *
*        qclust(i,nclst)                                               *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg

      use fragdatamod
      use moddas_fragdata

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'
      include 'param-physcnst.inc'

      parameter ( pi  = 3.1415926535898d0 )

*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)

      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)
      character frgfl*200

!--*-----------------------------------------------------------------------

               ne   = ifgs01_nei(ifg)
               neo  = ifgs04_neo(ifg)
               nag  = ifgs06_nag(ifg)
               nfrg = ifgs08_nfrg(ifg)
               kne  = ifgs02_kne(ifg)
               kxs  = ifgs03_kxs(ifg)
               kef  = ifgs05_kef(ifg)
               kaf  = ifgs07_kaf(ifg)
               kim  = ifgs09_kim(ifg)
               ks0  = ifgs10_ks0(ifg)

               if ( neo .gt. 0 ) then
                  if ( ifgdf(ifg,1) .eq. 5 ) then
                     nef = neo
                  else
                     nef = neo + 1
                  end if
               else
                  nef = neo
               end if
               if ( iabs(nag) .gt. 0 ) then
                  if ( ifgdf(ifg,1) .eq. 5 ) then
                     naf = iabs(nag)
                  else
                     naf = iabs(nag) + 1
                  end if
               else
                  naf = iabs(nag)
               end if

*-----------------------------------------------------------------------
*        inicident total energy
*-----------------------------------------------------------------------

               ipch = ifgdf(ifg,2) / 1000
               ipms = ifgdf(ifg,2) - ifgdf(ifg,2) / 1000 * 1000
               ipkf = ipch * 1000000 + ipms

               bein = bindeg(ipch,ipms-ipch)
               eint = ein * ipms
               etin = ein * ipms - bein

               esum = 0.0d0

*-----------------------------------------------------------------------
*        inicident energy if ein is out of e range, return
*-----------------------------------------------------------------------

            if( ein .lt. frgne(kne+1) ) return

         do i = 1, ne

            if( ein .lt. frgne(kne+i+1) ) goto 100

         end do

            return

  100    ie = i

*-----------------------------------------------------------------------

               ee1 = frgne(kne+ie)
               ee2 = frgne(kne+ie+1)
               eed = ee2 - ee1


               mm1 = ifrge_ksf(ks0+ie)
               mm2 = ifrge_ksf(ks0+ie+1)
               k11 = ifrge_ks1(ks0+ie)
               k12 = ifrge_ks2(ks0+ie)
               k13 = ifrge_ks3(ks0+ie)
               k21 = ifrge_ks1(ks0+ie+1)
               k22 = ifrge_ks2(ks0+ie+1)
               k23 = ifrge_ks3(ks0+ie+1)

*-----------------------------------------------------------------------
*        fragments
*-----------------------------------------------------------------------

         if( ifgdf(ifg,1) .ne. 2 ) then

            nclst = 0

         end if

         do k = 1, nfrg

                  wgt = ( ( ee2 - ein ) * frgsf(mm1+k)
     &                  - ( ee1 - ein ) * frgsf(mm2+k) )
     &                  / eed

            if( wgt .gt. 0.0d0 ) then

*-----------------------------------------------------------------------
C case 1
            if( neo .gt. 0 .and. nag .ne. 0 ) then

                  rnd = unirn(dummy)

               do l = 1, nef-1
               do m = 1, naf-1

                  crs = ( ( ee2 - ein )
     &                * frgdd(k12+(k-1)*neo*iabs(nag)+(l-1)*iabs(nag)+m)
     &                - ( ee1 - ein )
     &                * frgdd(k22+(k-1)*neo*iabs(nag)+(l-1)*iabs(nag)+m)
     &                ) / eed

                     if( rnd .le. crs ) goto 300

               end do
               end do

                  l = nef-1
                  m = naf-1

  300             ien = l
                  iag = m

                  eout = frgef(kef+ien)
     &                 + unirn(dummy) * ( frgef(kef+ien+1)
     &                                  - frgef(kef+ien) )

                  csth = dcos( frgaf(kaf+iag+1) )
     &                 + unirn(dummy) * ( dcos( frgaf(kaf+iag) )
     &                                  - dcos( frgaf(kaf+iag+1) ) )

                  if ( ifgdf(ifg,1) .eq. 5 ) then

                     call FragData_engang(1,ien,iag,k,ie,ifg
     &                    ,eout,csth,wxs)
                     wgt = wgt*wxs

                  end if

                  phih = 2.d0 * pi * unirn(dummy)

*-----------------------------------------------------------------------
C case 2
            else if ( neo .eq. 0 .and. nag .ne. 0 ) then

                  ems = ( ( ee2 - ein ) * frgee(k11+(k-1)*2+1)
     &                  - ( ee1 - ein ) * frgee(k21+(k-1)*2+1) )
     &                  / eed

                  eds = ( ( ee2 - ein ) * frgee(k11+(k-1)*2+2)
     &                  - ( ee1 - ein ) * frgee(k21+(k-1)*2+2) )
     &                  / eed

  400             eout =  ems + eds * gaurn(dummy)
                  if( eout .le. 0 ) goto 400

                  rnd = unirn(dummy)

               do m = 1, naf-1

                  ang = ( ( ee2 - ein ) * frgdx(k13+(k-1)*iabs(nag)+m)
     &                  - ( ee1 - ein ) * frgdx(k23+(k-1)*iabs(nag)+m) )
     &                  / eed

                  if( rnd .le. ang ) goto 200

               end do

                  m = naf-1

  200             iag = m

                  csth = dcos( frgaf(kaf+iag+1) )
     &                 + unirn(dummy) * ( dcos( frgaf(kaf+iag) )
     &                                  - dcos( frgaf(kaf+iag+1) ) )

                  if ( ifgdf(ifg,1) .eq. 5 ) then

                     ien = 1
                     call FragData_engang(2,ien,iag,k,ie,ifg
     &                    ,eout,csth,wxs)
                     wgt = wgt*wxs

                  end if

                  phih = 2.d0 * pi * unirn(dummy)

*-----------------------------------------------------------------------
C case 3
            else if ( neo .eq. 0 .and. nag .eq. 0 ) then

                  ems = ( ( ee2 - ein ) * frgee(k11+(k-1)*2+1)
     &                  - ( ee1 - ein ) * frgee(k21+(k-1)*2+1) )
     &                  / eed

                  eds = ( ( ee2 - ein ) * frgee(k11+(k-1)*2+2)
     &                  - ( ee1 - ein ) * frgee(k21+(k-1)*2+2) )
     &                  / eed

  500             eout =  ems + eds * gaurn(dummy)
                  if( eout .le. 0 ) goto 500

                  csth = 1d0 + unirn(dummy) * ( -2d0 )

                  phih = 2.d0 * pi * unirn(dummy)

*-----------------------------------------------------------------------
C case 4
            else if ( neo .gt. 0 .and. nag .eq. 0 ) then

                  rnd = unirn(dummy)

               do l = 1, nef-1

                  crs = ( ( ee2 - ein )
     &                  * frgdx(k13+(k-1)*neo+l)
     &                  - ( ee1 - ein )
     &                  * frgdx(k23+(k-1)*neo+l) )
     &                  / eed

                     if( rnd .le. crs ) goto 600

               end do

                  l = nef-1

  600             ien = l

                  eout = frgef(kef+ien)
     &                 + unirn(dummy) * ( frgef(kef+ien+1)
     &                                  - frgef(kef+ien) )

                  if ( ifgdf(ifg,1) .eq. 5 ) then

                     iag = 1
                     call FragData_engang(4,ien,iag,k,ie,ifg
     &                    ,eout,csth,wxs)
                     wgt = wgt*wxs

                  end if

                  csth = 1d0 + unirn(dummy) * ( -2d0 )

                  phih = 2.d0 * pi * unirn(dummy)

*-----------------------------------------------------------------------
C case 5
             else if ( neo .lt. 0 .and. nag .ne. 0 ) then

                  rnd = unirn(dummy)

               do l = 1, iabs(neo)
               do m = 1, iabs(nag)

                  crs = ( ( ee2 - ein )
     &                 * frgdd(k12+(k-1)*iabs(neo)*iabs(nag)
     &                 +(l-1)*iabs(nag)+m)
     &                 - ( ee1 - ein )
     &                 * frgdd(k22+(k-1)*iabs(neo)*iabs(nag)
     &                 +(l-1)*iabs(nag)+m) )
     &                 / eed

                     if( rnd .le. crs ) goto 700

               end do
               end do

                  l = iabs(neo)
                  m = iabs(nag)

  700             ien = l
                  iag = m

                  eout = frgef(kef+ien)

                  csth = dcos( frgaf(kaf+iag+1) )
     &                 + unirn(dummy) * ( dcos( frgaf(kaf+iag) )
     &                                  - dcos( frgaf(kaf+iag+1) ) )

                  phih = 2.d0 * pi * unirn(dummy)


*-----------------------------------------------------------------------
C case 6
             else if ( neo .lt. 0 .and. nag .eq. 0 ) then

                  rnd = unirn(dummy)

               do l = 1, iabs(neo)

                  crs = ( ( ee2 - ein )
     &                  * frgdx(k13+(k-1)*iabs(neo)+l)
     &                  - ( ee1 - ein )
     &                  * frgdx(k23+(k-1)*iabs(neo)+l) )
     &                  / eed

                  if( rnd .le. crs ) goto 800

               end do

                  l = iabs(neo)

  800             ien = l

                  eout = frgef(kef+ien)

                  csth = 1d0 + unirn(dummy) * ( -2d0 )

                  phih = 2.d0 * pi * unirn(dummy)

            end if

*-----------------------------------------------------------------------
*           booking
*-----------------------------------------------------------------------

            if( eout .gt. 0.0d0 ) then

CS.Hashimoto revised for extension of kind of particles. (2014.9.24) (start)
                  kf = ifrgm(kim+k)

               if( kf .eq. 2212 ) then ! proton

                  ibary = 1
                  ipid  = 1
                  ippad = 1
                  ipprt = 1
                  ipneu = 0
                  ipchg = 1

                  bene = 0.0d0
                  rms = rstms(1)

               else if( kf .eq. 2112 ) then ! neutron

                  ibary = 1
                  ipid  = 2
                  ippad = 2
                  ipprt = 0
                  ipneu = 1
                  ipchg = 0

                  bene = 0.0d0
                  rms = rstms(2)

               else if( kf .eq. 211 ) then ! pion+

                  ibary = 0
                  ipid  = 3
                  ippad = 3
                  ipprt = 0
                  ipneu = 0
                  ipchg = 1

                  bene = 0.0d0
                  rms = rstms(3)

               else if( kf .eq. 111 ) then ! pion0

                  ibary = 0
                  ipid  = 3
                  ippad = 4
                  ipprt = 0
                  ipneu = 0
                  ipchg = 0

                  bene = 0.0d0
                  rms = rstms(4)

               else if( kf .eq. -211 ) then ! pion-

                  ibary = 0
                  ipid  = 3
                  ippad = 5
                  ipprt = 0
                  ipneu = 0
                  ipchg = -1

                  bene = 0.0d0
                  rms = rstms(5)

               else if( kf .eq. -13 ) then ! muon+

                  ibary = 0
                  ipid  = 6
                  ippad = 6
                  ipprt = 0
                  ipneu = 0
                  ipchg = 1

                  bene = 0.0d0
                  rms = rstms(6)

               else if( kf .eq. 13 ) then ! muon-

                  ibary = 0
                  ipid  = 6
                  ippad = 7
                  ipprt = 0
                  ipneu = 0
                  ipchg = -1

                  bene = 0.0d0
                  rms = rstms(7)

               else if( kf .eq. 321 ) then ! kaon+

                  ibary = 0
                  ipid  = 5
                  ippad = 8
                  ipprt = 0
                  ipneu = 0
                  ipchg = 1

                  bene = 0.0d0
                  rms = rstms(8)

               else if( kf .eq. 311 ) then ! kaon0

                  ibary = 0
                  ipid  = 5
                  ippad = 9
                  ipprt = 0
                  ipneu = 0
                  ipchg = 0

                  bene = 0.0d0
                  rms = rstms(9)

               else if( kf .eq. -321 ) then ! kaon-

                  ibary = 0
                  ipid  = 5
                  ippad = 10
                  ipprt = 0
                  ipneu = 0
                  ipchg = -1

                  bene = 0.0d0
                  rms = rstms(10)

               else if( kf .eq. 11 ) then ! electron

                  ibary = 0
                  ipid  = 7
                  ippad = 12
                  ipprt = 0
                  ipneu = 0
                  ipchg = -1

                  bene = 0d0
                  rms = rstms(12)

               else if( kf .eq. -11 ) then ! positron

                  ibary = 0
                  ipid  = 7
                  ippad = 13
                  ipprt = 0
                  ipneu = 0
                  ipchg = 1

                  bene = 0d0
                  rms = rstms(13)

               else if( kf .eq. 22 ) then ! photon

                  ibary = 0
                  ipid  = 4
                  ippad = 14
                  ipprt = 0
                  ipneu = 0
                  ipchg = 0

                  bene = 0d0
                  rms = rstms(14)

               else if( kf .eq. 1000002 ) then ! deuteron

                  ibary = 2
                  ipid  = 0
                  ippad = 15
                  ipprt = 1
                  ipneu = 1
                  ipchg = 1

                  bene = bindeg(ipprt,ipneu)
                  rms = rstms(15)

                  eout = eout * ibary

               else if( kf .eq. 1000003 ) then ! triton

                  ibary = 3
                  ipid  = 0
                  ippad = 16
                  ipprt = 1
                  ipneu = 2
                  ipchg = 1

                  bene = bindeg(ipprt,ipneu)
                  rms = rstms(16)

                  eout = eout * ibary

               else if( kf .eq. 2000003 ) then ! 3He

                  ibary = 3
                  ipid  = 0
                  ippad = 17
                  ipprt = 2
                  ipneu = 1
                  ipchg = 2

                  bene = bindeg(ipprt,ipneu)
                  rms = rstms(17)

                  eout = eout * ibary

               else if( kf .eq. 2000004 ) then ! alpha

                  ibary = 4
                  ipid  = 0
                  ippad = 18
                  ipprt = 2
                  ipneu = 2
                  ipchg = 2

                  bene = bindeg(ipprt,ipneu)
                  rms = rstms(18)

                  eout = eout * ibary

               else

                  ipid  = 0
                  ibary = mod(kf,1000000)
                  ipprt = ( kf - ibary ) / 1000000 + 0.1
                  ipneu = ibary - ipprt
                  ippad = 19
                  ipchg = ipprt

                  bene = bindeg(ipprt,ipneu)
                  rms = 0.93827d0 * ipprt + 0.93958d0 * ipneu
     &                - bene / 1000.0d0

                  eout = eout * ibary

               end if

                  tmas = rms / 1000.d0

                  etot = tmas + eout / 1000.d0
                  psqr = dsqrt( etot**2 - tmas**2 )

                  snth = dsqrt( 1.d0 - csth**2 )

                  pclx = psqr * snth * dcos( phih )
                  pcly = psqr * snth * dsin( phih )
                  pclz = psqr * csth

                  esum = esum + ( eout - bene ) * wgt

*-----------------------------------------------------------------------

                  nclst = nclst + 1
                  ii    = nclst

                  iclust(ii)    = ipid

                  jclust(0,ii)  = 0
                  jclust(1,ii)  = ipprt
                  jclust(2,ii)  = ipneu
                  jclust(3,ii)  = ippad
                  jclust(4,ii)  = 0
                  jclust(5,ii)  = ipchg
                  jclust(6,ii)  = ibary
                  jclust(7,ii)  = kf
                  jclust(8,ii)  = 0

                  qclust(0,ii)  = 1.0
                  qclust(1,ii)  = pclx
                  qclust(2,ii)  = pcly
                  qclust(3,ii)  = pclz
                  qclust(4,ii)  = etot
                  qclust(5,ii)  = tmas
                  qclust(6,ii)  = 0.0d0
                  qclust(7,ii)  = eout
                  qclust(8,ii)  = wgt
                  qclust(9,ii)  = 0.0d0
                  qclust(10,ii) = 0.0d0
                  qclust(11,ii) = 0.0d0
                  qclust(12,ii) = 0.0d0

                  numpat(ippad) = numpat(ippad) + 1
                  rumpat(ippad) = rumpat(ippad) + wgt

            end if

*-----------------------------------------------------------------------

            end if

         end do

*-----------------------------------------------------------------------
*        flux correction for opt = 3
*-----------------------------------------------------------------------

            if( etin - esum .gt. 0.0d0 .and.
     &          ifgdf(ifg,1) .eq. 3 ) then

                  wgt  = ( etin - esum ) / eint
                  tmas = rmtyp(19,ipkf) / 1000.d0
                  etot = tmas + eint / 1000.d0
                  psqr = dsqrt( etot**2 - tmas**2 )

                  nclst = nclst + 1
                  ii    = nclst

                  iclust(ii)    = 0

                  jclust(0,ii)  = 0
                  jclust(1,ii)  = ipch
                  jclust(2,ii)  = ipms - ipch
                  jclust(3,ii)  = 19
                  jclust(4,ii)  = 0
                  jclust(5,ii)  = ipch
                  jclust(6,ii)  = ipms
                  jclust(7,ii)  = ipkf
                  jclust(8,ii)  = 0

                  qclust(0,ii)  = 1.0
                  qclust(1,ii)  = 0.0d0
                  qclust(2,ii)  = 0.0d0
                  qclust(3,ii)  = psqr
                  qclust(4,ii)  = etot
                  qclust(5,ii)  = tmas
                  qclust(6,ii)  = 0.0d0
                  qclust(7,ii)  = eint
                  qclust(8,ii)  = wgt
                  qclust(9,ii)  = 0.0d0
                  qclust(10,ii) = 0.0d0
                  qclust(11,ii) = 0.0d0
                  qclust(12,ii) = 0.0d0

                  numpat(19) = numpat(19) + 1
                  rumpat(19) = rumpat(19) + wgt

            end if
CS.Hashimoto revised for extension of kind of particles. (2014.9.24) (end)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine esepa(iz,in,esn,esp,esa)
*                                                                      *
*                                                                      *
*       created by T.Ogawa on 2016/09/26                               *
*                                                                      *
*       Purpose:                                                       *
*                                                                      *
*              to give particle separation energy of nucleus iz and    *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       iz      : nucleus proton number                                *
*       in      : nucleus neutron number                               *
*                                                                      *
*     output  :                                                        *
*                                                                      *
*       esn     : neutron separation energy (MeV)                      *
*       esp     : proton separation energy  (MeV)                      *
*       esa     : alpha separation energy   (MeV)                      *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------
      if(iz + in .eq. 1) then ! nucleon
        esn = 1.d20 ! dummy large number
        esp = 1.d20 ! dummy large number
        esa = 1.d20 ! dummy large number
        return
      else
        remin0 = bindeg( iz, in ) ! define initial binding E
        remin1 = bindeg( iz, in - 1)
        remin2 = bindeg( iz -1, in)
        esn = remin0 - remin1
        esp = remin0 - remin2
     &   + vcoul(dble(iz), dble(iz + in)
     &   , 1, 1, dostg(1,dble(iz - 1)), 2) ! proton separation energy
      endif

      if(iz .ge. 2 .and. in .ge. 2 .and. iz + in .gt. 4) then
        remin3 = bindeg( iz -2, in -2)
        esa = remin0 - remin3 - bindeg(2,2)
     &   +vcoul(dble(iz), dble(iz + in)
     &   , 2, 4, dostg(2,dble(iz - 2)), 6) ! alpha separation energy
      else
        esa = 1.d20 ! dummy large number
      endif

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------

      return
      end

