************************************************************************
*                                                                      *
      subroutine splitng(ncol)
*                                                                      *
*                                                                      *
*     play the splitting of particle for special crossing with counter *
*                                                                      *
*     2004/12/05 : last revised by K.Niita                             *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use moddas_region
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /inout/  in,io

*-----------------------------------------------------------------------

      common /bnkmem/ maxbnk, maxbn2, rtrckflp
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /mulchck/ mulover
!$OMP THREADPRIVATE(/mulchck/)

*-----------------------------------------------------------------------

      common /regdc/ idrg(kvlmax), idgr(kvmmax)
      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      common /splreg/ isptn, npreg(6), mnspt(6,0:20),
     &                ipgrc(6), ipgrt(6), ksplt(6), isplt(6),
     &                ispct(6,9), ispem(6,2)
      common /splrge/ espem(6,2)

*-----------------------------------------------------------------------

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      data eps/1.d-4/

*-----------------------------------------------------------------------
*        for after crossing and isptn is not zero
*        not for neutrino
*-----------------------------------------------------------------------

         if( ncol .ne. 10 ) return
         if( isptn .eq. 0 ) return

         if( abs( ktyp ) .eq. 12 .or. abs( ktyp ) .eq. 14 ) return

*-----------------------------------------------------------------------
*        no process : zero weight
*-----------------------------------------------------------------------

         if( wt(ibkwt+no,ipomp+1) .le. 0.0 ) return

*-----------------------------------------------------------------------
*     do loop of isptn
*-----------------------------------------------------------------------

      do 1000 ll = 1, isptn

*-----------------------------------------------------------------------
*        check particle
*-----------------------------------------------------------------------

         if( mnspt(ll,20) .lt. 19 ) then

            do j = 1, mnspt(ll,20)

               if( ityp .eq. mnspt(ll,j) ) goto 20

            end do

               goto 1000

   20       continue

         end if

*-----------------------------------------------------------------------
*        check energy
*-----------------------------------------------------------------------

                enrg = ec(ibkec+no,ipomp+1)

            if( ityp .ge. 15 ) enrg = enrg / dble( mtyp )

            if( ispem(ll,1) .ne. 0 ) then

               if( enrg .lt. espem(ll,1) ) goto 1000

            end if

            if( ispem(ll,2) .ne. 0 ) then

               if( enrg .gt. espem(ll,2) ) goto 1000

            end if

*-----------------------------------------------------------------------
*        check counter
*-----------------------------------------------------------------------

         do m = 1, 3

            if( ispct(ll,m) .ne. 0 ) then

               if( ncnt(ibknct+m,no,ipomp+1) .lt. ispct(ll,2*m+2) .or.
     &             ncnt(ibknct+m,no,ipomp+1) .gt. ispct(ll,2*m+3) )
     &            goto 1000

            end if

         end do

*-----------------------------------------------------------------------
*        check crossing regions
*-----------------------------------------------------------------------

                     idsm = ipgrt(ll)
                     jdsm = -1

                     kdsm = ksplt(ll)
                     ldsm = 0

            do ir = 1, npreg(ll)

                     ldsm = ldsm + 1
                     sfac = das_ksplt(kdsm+ldsm)

                     jdsm = jdsm + 1
                     ntrn = idas_ipgrt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_ipgrt(idsm+jdsm)

                     jj = 0

                  do k = 1, ntrn

                     call tregck(iblz1,ilev1,ilat1,
     &                           mtrn,idas_ipgrt(idsm+jdsm+1),jj,ic1)

                  end do

                     jdsm = jdsm + mtrn

                     jdsm = jdsm + 1
                     ntrn = idas_ipgrt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_ipgrt(idsm+jdsm)

               if( ic1 .ne. 0 ) then

                     jj  = 0

                  do k = 1, ntrn

                     call tregck(iblz2,ilev2,ilat2,
     &                           mtrn,idas_ipgrt(idsm+jdsm+1),jj,ic2)

                     if( ic2 .ne. 0 ) goto 400

                  end do

               end if

                     jdsm = jdsm + mtrn

            end do

               goto 1000

*-----------------------------------------------------------------------
*        set up initial data and condition
*-----------------------------------------------------------------------

  400    continue

            nabov = 0
            weit  = wt(ibkwt+no,ipomp+1)

            rj = sfac

*-----------------------------------------------------------------------
*     splitting
*-----------------------------------------------------------------------

      if( rj .gt. 1.0d0 ) then

*-----------------------------------------------------------------------

               weight = weit / rj
               wt(ibkwt+no,ipomp+1) = weight

*-----------------------------------------------------------------------

c S.H. revised to avoid an error of a large rj value (2022.12.16)
              if ( rj .ge. 2147483647d0 ) then
               goto 9010
              else
               mmm = int(rj) - 1
              end if

               rp = rj - dble( mmm + 1 )

               if( rp .ge. eps ) then

                  r = unirn(dummy)
                  if( r .le. rp ) mmm = mmm + 1

               end if

*-----------------------------------------------------------------------

         if( mmm .gt. 0 ) then

               nabov = mmm

               if( nabov .gt. 30 .and. mulover .le. 5 ) then
                 mulover = mulover + 1
                 write(*,*) "Warning : Particle is split by", nabov,
     &"Check splitting."
                 if(mulover .gt. 5) write(*,*) "This message is not show
     &n anymore"
               endif

               if( nabov .gt. maxbn2 ) goto 9010

            do m = 1, nabov

               ntya(iaknty+m,ipomp+1) = nty(ibknty+no,ipomp+1)
               nkfa(iaknkf+m,ipomp+1) = nkf(ibknkf+no,ipomp+1)

               ea(iake+m,ipomp+1)   = ec(ibkec+no,ipomp+1)
               ta(iakt+m,ipomp+1)   = tc(ibktc+no,ipomp+1)
               wta(iakwt+m,ipomp+1) = wt(ibkwt+no,ipomp+1)

               ua(iaku+m,ipomp+1) = u(ibku+no,ipomp+1)
               va(iakv+m,ipomp+1) = v(ibkv+no,ipomp+1)
               wa(iakw+m,ipomp+1) = w(ibkw+no,ipomp+1)

               namea(iaknam+m,ipomp+1) = name(ibknam+no,ipomp+1)

               xa(iakx+m,ipomp+1) = xc(ibkxc+no,ipomp+1) +
     &                              u(ibku+no,ipomp+1) * parz(28)
               ya(iaky+m,ipomp+1) = yc(ibkyc+no,ipomp+1) +
     &                              v(ibkv+no,ipomp+1) * parz(28)
               za(iakz+m,ipomp+1) = zc(ibkzc+no,ipomp+1) +
     &                              w(ibkw+no,ipomp+1) * parz(28)

               iblza(iakblz+m,ipomp+1) = iblz(ibkblz+no,ipomp+1)
               nmeda(iaknmd+m,ipomp+1) = nmed(ibknmd+no,ipomp+1)
               wtina(iakwin+m,ipomp+1) = wtin(ibkwin+no,ipomp+1)
               wtnza(iakwnz+m,ipomp+1) = wtnz(ibkwnz+no,ipomp+1)

               spxa(iakspx+m,ipomp+1) = spx(ibkspx+no,ipomp+1)
               spya(iakspy+m,ipomp+1) = spy(ibkspy+no,ipomp+1)
               spza(iakspz+m,ipomp+1) = spz(ibkspz+no,ipomp+1)

               nzsta(iakzst+m,ipomp+1) = nzst(ibkzst+no,ipomp+1)
               nsosa(iaksos+m,ipomp+1) = nsos(ibksos+no,ipomp+1)

               ncnta(iaknct+1,m,ipomp+1) = ncnt(ibknct+1,no,ipomp+1)
               ncnta(iaknct+2,m,ipomp+1) = ncnt(ibknct+2,no,ipomp+1)
               ncnta(iaknct+3,m,ipomp+1) = ncnt(ibknct+3,no,ipomp+1)

               nfcsa(iaknfc+m,ipomp+1) = nfcs(ibknfc+no,ipomp+1)
               xfcsa(iakxfc+m,ipomp+1) = xfcs(ibkxfc+no,ipomp+1)
              itetposa(ibtetposa+m,ipomp+1)=itetpos(ibtetpos+no,ipomp+1)

            end do

         end if

*-----------------------------------------------------------------------
*     russian roulette
*-----------------------------------------------------------------------

      else if( rj .lt. 1.0d0 ) then

*-----------------------------------------------------------------------

            r = unirn(dummy)

         if( r .lt. rj ) then

            wt(ibkwt+no,ipomp+1) = weit / rj

         else

            wt(ibkwt+no,ipomp+1) =  0.0d0

         end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

 1000 continue

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*     error message
*-----------------------------------------------------------------------

 9010 continue

         write(io,12) nabov,maxbn2,nocas
   12    format(/' *** error messge from splitng ***'
     &   /' storage size (nabov) of particles is greater ',
     &    'than maxbn2.'
     &   /' nabov  =',i10
     &   /' maxbn2 =',i10
     &   /' nocas  =',i10)
         call parastop( 837 )

         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine wwindw(ncol)
*                                                                      *
*                                                                      *
*     play the splitting or russian roulette of particle with the cell *
*     weight window when particle reach the boundary surface           *
*     or after collision                                               *
*                                                                      *
*        iwit =  1   : (default) zero weight -> weight cutoff          *
*        iwit =  0   : outside the window                              *
*        iwit = -1   : inside the window                               *
*        iwit = -10  : kill the particle                               *
*                                                                      *
*     2005/12/20 : last revised by K.Niita                             *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use moddas_mesh
      use moddas_variance_reduction
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'err.inc'
*-----------------------------------------------------------------------

      common /inout/  in,io

*-----------------------------------------------------------------------

      common /bnkmem/ maxbnk, maxbn2, rtrckflp
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)

      dimension jlat(5,10)

*-----------------------------------------------------------------------

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

      common /regdc/ idrg(kvlmax), idgr(kvmmax)
      common /paraj/ mstz(300), parz(300)

      common /mulchck/ mulover
!$OMP THREADPRIVATE(/mulchck/)
*-----------------------------------------------------------------------

      common /impmsg/ iimpn, isimp, iswct, isstr, maxip,
     &                mnimp(7,0:20),
     &                kfimp(7), inimc(7), inimt(7)
      common /wparm/  swtm(20), wc1(20), wc2(20) !FURUTA
      common /wparm0/ wc01(20), wc02(20)         !FURUTA
!$OMP THREADPRIVATE(/wparm0/)

      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww

*-----------------------------------------------------------------------

      common /wwxyzp/ iwxpp, iwypp, iwzpp
!$OMP THREADPRIVATE(/wwxyzp/)
      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)

*-----------------------------------------------------------------------


      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------
cKN 2024/03/26

      common /bnkmew/ iwwbnk

*-----------------------------------------------------------------------
cKN 2024/03/26

      common /wwindp/ wupn, wsurvn, mxspln, mwhere, mvoww
      common /wwindp0/ wupn0, wsurvn0, mxspln0, mwhere0, mvoww0
!$OMP THREADPRIVATE(/wwindp0/)

*-----------------------------------------------------------------------

      data eps/1.d-4/

      data imessage1/0/  ! T.Sato 2024/03/26
      data imessage2/0/  ! T.Sato 2024/03/26
      data imessage3/0/  ! T.Sato 2024/05/26

*-----------------------------------------------------------------------
cKN 2024/03/26

         if( ncol .eq. 4 ) then

               wupn0   = wupn
               wsurvn0 = wsurvn
               mxspln0 = mxspln
               mwhere0 = mwhere
               mvoww0  = mvoww

               return

         end if

cKN 2024/03/26
*-----------------------------------------------------------------------
*        for after crossing or collision
*        not for neutrino
*        no process : zero weight or equal of cell importance
*-----------------------------------------------------------------------

         if( ncol .eq. 10 .and. idmn(max(0,mat)) .eq. 0 .and.
     &       mvoww0 .ne. 0 ) then

            if( idmn( nmed(ibknmd+no,ipomp+1) ) .eq. 0 ) return

         end if

*-----------------------------------------------------------------------

         if( abs( ktyp ) .eq. 12 .or. abs( ktyp ) .eq. 14 ) return
         if( ncol .eq. 10 .and. nfcs(ibknfc+no,ipomp+1) .gt. 0 ) return
         if( wt(ibkwt+no,ipomp+1) .le. 0.0 ) return

         if( ncol .ne. 10 .and. ncol .ne. 13 .and.
     &       ncol .ne. 14 .and. ncol .ne. 16 ) return

         if( ncol .eq. 10 .and.
     &        .not.(iwmsh .eq. 1 .or. iwmsh .eq. 4) ) return
         if( ncol .eq. 16 .and. iwmsh .ne. 3 ) return

         if( mwhere0 .eq. -1 .and.
     &     ( ncol .eq. 10 .or. ncol .eq. 16 ) ) return
         if( mwhere0 .eq.  1 .and.
     &     ( ncol .ne. 10 .and. ncol .ne. 16 ) ) return

*-----------------------------------------------------------------------
*     Special for WWG
*-----------------------------------------------------------------------
cKN 2024/03/26

         if( wupn0.lt.1.0e7 ) then ! HDD has been already used. no more split

         if( iwwbnk.ne.0.and.nomax-no.gt.maxbnk/3 ) then

            mxspln0 = 2
            mwhere0 = 1
            wupn0   = 1000.0d0
            wsurvn0 = wupn * 0.6

            if(imessage3.eq.0) then ! first time
             imessage3=1
      write(ErrCha,'("Particle multiplication warning: level 3.",
     &" Increase in maxbnk is encouraged if you use [weight window]")')
             ErrID = 'L:521/R:wwindw/F:celimp.f'
             call ErrWrite(ErrID,ErrCha)
            endif

         else if( iwwbnk.ne.0.and.nomax-no.gt.maxbnk/5 ) then

            mxspln0 = 2
            mwhere0 = 1
            wupn0   = 200.d0
            wsurvn0 = wupn * 0.6

!            if(imessage2.eq.0) then ! first time
!             imessage2=1
!      write(ErrCha,'("Particle multiplication warning: level 2.",
!     &" Increase in iwwbnk is encouraged if you use [weight window]")')
!             ErrID = 'L:576/R:wwindw/F:celimp.f'
!             call ErrWrite(ErrID,ErrCha)
!            endif

         else if( iwwbnk.ne.0.and.nomax-no.gt.maxbnk/10 ) then

            mxspln0 = 3
!            mwhere0 = mwhere0
            wupn0   = 50.d0
            wsurvn0 = wupn * 0.6

!            if(imessage1.eq.0) then ! first time
!             imessage1=1
!      write(ErrCha,'("Particle multiplication warning: level 1.",
!     &" Increase in maxbnk is encouraged if you use [weight window]")')
!             ErrID = 'L:591/R:wwindw/F:celimp.f'
!             call ErrWrite(ErrID,ErrCha)
!            endif

         end if
         end if

cKN 2024/03/26
*-----------------------------------------------------------------------
*        initial nabov
*-----------------------------------------------------------------------

            nabob = nabov

*-----------------------------------------------------------------------
*        set up initial data and condition
*-----------------------------------------------------------------------

                  enrg = ec(ibkec+no,ipomp+1)
                  if( ityp .ge. 15 ) enrg = enrg / dble( mtyp )

                  timn = abs(tc(ibktc+no,ipomp+1))
                  weit = wt(ibkwt+no,ipomp+1)

            if( ncol .eq. 10 ) then

               if( iwmsh .eq. 1 ) then

                  jblz = iblz2
                  jlev = ilev2

                  do i = 1, 5
                  do j = 1, jlev
                     jlat(i,j) = ilat2(i,j)
                  end do
                  end do

               else if( iwmsh .eq. 4 ) then

                jblz = iblz2
                iwxp = iii
                iwzp = kkk

               end if

            else if( ncol .eq. 16 ) then

                  iwxp = iwxpp
                  iwyp = iwypp
                  iwzp = iwzpp

            else

               if( iwmsh .eq. 1 ) then

                  jblz = iblz1
                  jlev = ilev1

                  do i = 1, 5
                  do j = 1, jlev
                     jlat(i,j) = ilat1(i,j)
                  end do
                  end do

               else if( iwmsh .eq. 3 ) then

                  iwctl = 0

                  call wwxdis(iwctl,dpx,dwwt,
     &                        xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                        zc(ibkzc+no,ipomp+1),
     &                        u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                        w(ibkw+no,ipomp+1),
     &                        iwxnm(1),iwynm(1),iwznm(1),
     &                        das_iwxrg(iwxrg(1)),das_iwyrg(iwyrg(1)),
     &                        das_iwzrg(iwzrg(1)))

                  iwxp = iwxpp
                  iwyp = iwypp
                  iwzp = iwzpp

               else if( iwmsh .eq. 4 ) then

                jblz = iblz1
                iwxp = iii
                iwzp = kkk

               end if

            end if

*-----------------------------------------------------------------------
*        not weight window particles, produced particles
*-----------------------------------------------------------------------
cKN 2018/01/08  !!!! previous, wrong position of goto 500

         if( iwwin(ityp) .eq. 0 ) goto 500
         if( ncol .eq. 13 ) goto 500

*-----------------------------------------------------------------------
*        ncol = 14, fcl < 0 ; goto 500
*-----------------------------------------------------------------------

         if( ncol .eq. 14 .and. nfcs(ibknfc+no,ipomp+1) .lt. 0 ) then

            nfcs(ibknfc+no,ipomp+1) = 0

            goto 500

         end if

*-----------------------------------------------------------------------
****  ncol = 10, 16, 14 (fcl>0)  ; cross surface
*-----------------------------------------------------------------------
*        iwit =  1   : (default) zero weight -> weight cutoff
*        iwit =  0   : outside the window
*        iwit = -1   : inside the window
*        iwit = -10  : kill the particle
*-----------------------------------------------------------------------

               call awwd(ityp,jblz,jlev,jlat,ii2,
     &                   iwxp,iwyp,iwzp,
     &                   enrg,timn,weit,cbup,cbsv,cbdn,iwit)

*-----------------------------------------------------------------------
*     iwit = 0 : outside the window
*         for boundary or after collisions, ( ncol = 10, 13, 14, 16 )
*-----------------------------------------------------------------------

      if( iwit .eq. 0 ) then

*-----------------------------------------------------------------------
*        splitting
*-----------------------------------------------------------------------

         if( weit .gt. cbup ) then

*-----------------------------------------------------------------------

                  rj = weit / cbup

c S.H. revised to avoid an error of a large rj value (2022.12.13)
               if ( rj + 1d0 .le. dble(mxspln0) ) then
                  mmm = int(rj) + 1
               else
                  mmm = mxspln0
               end if

*-----------------------------------------------------------------------

                        weight = weit / dble(mmm)
                        wt(ibkwt+no,ipomp+1) = weight

                  if( mmm - 1 .gt. 0 ) then

                        nabov = nabov + mmm - 1

                        if( nabov .gt. 30 .and. mulover .le. 5 ) then
                          mulover = mulover + 1
                          write(*,*) "Warning : Particle is split by",
     &                                nabov, "Check [Weight Window]."
                          if(mulover .gt. 5) write(*,*) "This message is
     & not shown anymore"
                        endif

                        if( nabov .gt. maxbn2 ) goto 9010

                     do m = nabob + 1, nabob + mmm - 1

                        ntya(iaknty+m,ipomp+1) = nty(ibknty+no,ipomp+1)
                        nkfa(iaknkf+m,ipomp+1) = nkf(ibknkf+no,ipomp+1)

                        ea(iake+m,ipomp+1)   = ec(ibkec+no,ipomp+1)
                        ta(iakt+m,ipomp+1)   = tc(ibktc+no,ipomp+1)
                        wta(iakwt+m,ipomp+1) = wt(ibkwt+no,ipomp+1)

                        ua(iaku+m,ipomp+1) = u(ibku+no,ipomp+1)
                        va(iakv+m,ipomp+1) = v(ibkv+no,ipomp+1)
                        wa(iakw+m,ipomp+1) = w(ibkw+no,ipomp+1)

                        namea(iaknam+m,ipomp+1) =name(ibknam+no,ipomp+1)

                        xa(iakx+m,ipomp+1) = xc(ibkxc+no,ipomp+1) +
     &                                     u(ibku+no,ipomp+1) * parz(28)
                        ya(iaky+m,ipomp+1) = yc(ibkyc+no,ipomp+1) +
     &                                     v(ibkv+no,ipomp+1) * parz(28)
                        za(iakz+m,ipomp+1) = zc(ibkzc+no,ipomp+1) +
     &                                     w(ibkw+no,ipomp+1) * parz(28)

                        iblza(iakblz+m,ipomp+1) =iblz(ibkblz+no,ipomp+1)
                        nmeda(iaknmd+m,ipomp+1) =nmed(ibknmd+no,ipomp+1)
                        wtina(iakwin+m,ipomp+1) =wtin(ibkwin+no,ipomp+1)
                        wtnza(iakwnz+m,ipomp+1) =wtnz(ibkwnz+no,ipomp+1)

                        spxa(iakspx+m,ipomp+1) = spx(ibkspx+no,ipomp+1)
                        spya(iakspy+m,ipomp+1) = spy(ibkspy+no,ipomp+1)
                        spza(iakspz+m,ipomp+1) = spz(ibkspz+no,ipomp+1)

                        nzsta(iakzst+m,ipomp+1) =nzst(ibkzst+no,ipomp+1)
                        nsosa(iaksos+m,ipomp+1) =nsos(ibksos+no,ipomp+1)

                        ncnta(iaknct+1,m,ipomp+1) =
     &                                         ncnt(ibknct+1,no,ipomp+1)
                        ncnta(iaknct+2,m,ipomp+1) =
     &                                         ncnt(ibknct+2,no,ipomp+1)
                        ncnta(iaknct+3,m,ipomp+1) =
     &                                         ncnt(ibknct+3,no,ipomp+1)

                        nfcsa(iaknfc+m,ipomp+1) =nfcs(ibknfc+no,ipomp+1)
                        xfcsa(iakxfc+m,ipomp+1) =xfcs(ibkxfc+no,ipomp+1)
                        itetposa(ibtetposa+m,ipomp+1)=
     &                                   itetpos(ibtetpos+no,ipomp+1)

                     end do

                        wtxwwp(iswwp+ityp,1,ii2) =
     &                  wtxwwp(iswwp+ityp,1,ii2) + weit
                        wtxwwp(iswwp+ityp,2,ii2) =
     &                  wtxwwp(iswwp+ityp,2,ii2) + weight
                        wtxwwp(iswwp+ityp,3,ii2) =
     &                  wtxwwp(iswwp+ityp,3,ii2) + 1.d+0

                  end if

*-----------------------------------------------------------------------
*        russian roulette
*-----------------------------------------------------------------------

         else if( weit .lt. cbdn ) then

*-----------------------------------------------------------------------

                  if( cbsv .gt. 0.0d0 ) then

                        ri = max( weit / cbsv, 1.0d0 / dble( mxspln0 ) )
                        ccsv = weit / ri

                        r  = unirn(dummy)

                        k = 3
                        if( ncol .ne. 10 .and. ncol .ne. 16 ) k = 7

                     if( r .lt. ri ) then

                        wtxwwp(iswwp+ityp,k+1,ii2) =
     &                  wtxwwp(iswwp+ityp,k+1,ii2) + ccsv - weit
                        wtxwwp(iswwp+ityp,k+2,ii2) =
     &                  wtxwwp(iswwp+ityp,k+2,ii2) + 1.d+0

                        wt(ibkwt+no,ipomp+1) = ccsv

                     else

                        wtxwwp(iswwp+ityp,k+3,ii2) =
     &                  wtxwwp(iswwp+ityp,k+3,ii2) + weit
                        wtxwwp(iswwp+ityp,k+4,ii2) =
     &                  wtxwwp(iswwp+ityp,k+4,ii2) + 1.d+0

                        wt(ibkwt+no,ipomp+1) =  0.0d0

                     end if

                  end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*     iwit = 1 : zero factor, weight cutoff with 1-for-2 roullet limit
*-----------------------------------------------------------------------

      else if( iwit .eq. 1 ) then

                        cinimp = wtin(ibkwin+no,ipomp+1)

                     if( ncol .eq. 10 ) then

                        cotimp = aimp(ityp,iblz2,ilev2,ilat2,ii4)

                     else

                        cotimp = aimp(ityp,iblz1,ilev1,ilat1,ii4)

                     end if

                        rjsuv = wc01(ityp) / cotimp * cinimp
                        rjcut = wc02(ityp) / cotimp * cinimp

               if( rjsuv .gt. 0.0d0 ) then

                  if( weit .lt. rjcut ) then

                        ri = max( weit / rjsuv, 0.5d0 )
                        rrsuv = weit / ri

                        r  = unirn(dummy)

                     if( r .lt. ri ) then

                        wtxcut(iswct+ityp,1,ii4) =
     &                  wtxcut(iswct+ityp,1,ii4) + rrsuv - weit
                        wtxcut(iswct+ityp,2,ii4) =
     &                  wtxcut(iswct+ityp,2,ii4) + 1.d+0

                        wt(ibkwt+no,ipomp+1) = rrsuv

                     else

                        wtxcut(iswct+ityp,3,ii4) =
     &                  wtxcut(iswct+ityp,3,ii4) + weit
                        wtxcut(iswct+ityp,4,ii4) =
     &                  wtxcut(iswct+ityp,4,ii4) + 1.d+0

                        wt(ibkwt+no,ipomp+1) =  0.0d0

                     end if

                  end if

               end if

*-----------------------------------------------------------------------
*     iwit = -10 : kill this particles
*-----------------------------------------------------------------------

      else if( iwit .eq. -10 ) then

                        k = 3
                        if( ncol .ne. 10 ) k = 7

                        wtxwwp(iswwp+ityp,k+3,ii2) =
     &                  wtxwwp(iswwp+ityp,k+3,ii2) + weit
                        wtxwwp(iswwp+ityp,k+4,ii2) =
     &                  wtxwwp(iswwp+ityp,k+4,ii2) + 1.d+0

                        wt(ibkwt+no,ipomp+1) = 0.0d0

      end if

*-----------------------------------------------------------------------
*     check above bank only for collision
*-----------------------------------------------------------------------

  500 continue

      if( nabob .gt. 0 ) then

*-----------------------------------------------------------------------
*        splitting
*-----------------------------------------------------------------------

         do 300 ll = 1, nabob

                        weit = wta(iakwt+ll,ipomp+1)
                        itys = ntya(iaknty+ll,ipomp+1)

                  if( wta(iakwt+ll,ipomp+1) .le. 0.0d0 ) goto 300
                  if( iwwin(itys) .eq. 0 ) goto 300
                  if( nfcsa(iaknfc+ll,ipomp+1) .ne. 0 ) goto 300

*-----------------------------------------------------------------------
*                 weight window for above particle
*-----------------------------------------------------------------------

                     enrg = ea(iake+ll,ipomp+1)
                     timn = ta(iakt+ll,ipomp+1)

                     call awwd(itys,jblz,jlev,jlat,ii3,
     &                         iwxp,iwyp,iwzp,
     &                         enrg,timn,weit,cbup,cbsv,cbdn,jwit)

*-----------------------------------------------------------------------

               if( jwit .eq. 0 .and. weit .gt. cbup ) then

                  rj = weit / cbup

c S.H. revised to avoid an error of a large rj value (2022.12.16)
               if ( rj + 1d0 .le. dble(mxspln0) ) then
                  mmm = int(rj) + 1
               else
                  mmm = mxspln0
               end if

*-----------------------------------------------------------------------

                        weight = weit / dble(mmm)
                        wta(iakwt+ll,ipomp+1) = weight

                  if( mmm - 1 .gt. 0 ) then

                        naboa = nabov
                        nabov = nabov + mmm - 1

                        if( nabov .gt. 30 .and. mulover .le. 5 ) then
                          mulover = mulover + 1
                          write(*,*) "Warning : Particle is split by",
     &                                nabov, "Check [Weight Window]."
                          if(mulover .gt. 5) write(*,*) "This message is
     & not shown anymore"
                        endif

                        if( nabov .gt. maxbn2 ) goto 9010

                     do m = naboa + 1, naboa + mmm - 1

                        ea(iake+m,ipomp+1) = ea(iake+ll,ipomp+1)
                        namea(iaknam+m,ipomp+1)=namea(iaknam+ll,ipomp+1)
                        ntya(iaknty+m,ipomp+1) = ntya(iaknty+ll,ipomp+1)
                        nkfa(iaknkf+m,ipomp+1) = nkfa(iaknkf+ll,ipomp+1)

                        ua(iaku+m,ipomp+1) = ua(iaku+ll,ipomp+1)
                        va(iakv+m,ipomp+1) = va(iakv+ll,ipomp+1)
                        wa(iakw+m,ipomp+1) = wa(iakw+ll,ipomp+1)
                        wta(iakwt+m,ipomp+1) = wta(iakwt+ll,ipomp+1)

                        xa(iakx+m,ipomp+1) = xa(iakx+ll,ipomp+1)
                        ya(iaky+m,ipomp+1) = ya(iaky+ll,ipomp+1)
                        za(iakz+m,ipomp+1) = za(iakz+ll,ipomp+1)
                        ta(iakt+m,ipomp+1) = ta(iakt+ll,ipomp+1)

                        iblza(iakblz+m,ipomp+1)=iblza(iakblz+ll,ipomp+1)
                        nmeda(iaknmd+m,ipomp+1)=nmeda(iaknmd+ll,ipomp+1)
                        wtina(iakwin+m,ipomp+1)=wtina(iakwin+ll,ipomp+1)
                        wtnza(iakwnz+m,ipomp+1)=wtnza(iakwnz+ll,ipomp+1)

                        spxa(iakspx+m,ipomp+1) = spxa(iakspx+ll,ipomp+1)
                        spya(iakspy+m,ipomp+1) = spya(iakspy+ll,ipomp+1)
                        spza(iakspz+m,ipomp+1) = spza(iakspz+ll,ipomp+1)

                        nzsta(iakzst+m,ipomp+1)=nzsta(iakzst+ll,ipomp+1)
                        nsosa(iaksos+m,ipomp+1)=nsosa(iaksos+ll,ipomp+1)

                        ncnta(iaknct+1,m,ipomp+1) =
     &                                        ncnta(iaknct+1,ll,ipomp+1)
                        ncnta(iaknct+2,m,ipomp+1) =
     &                                        ncnta(iaknct+2,ll,ipomp+1)
                        ncnta(iaknct+3,m,ipomp+1) =
     &                                        ncnta(iaknct+3,ll,ipomp+1)

                        nfcsa(iaknfc+m,ipomp+1)=nfcsa(iaknfc+ll,ipomp+1)
                        xfcsa(iakxfc+m,ipomp+1)=xfcsa(iakxfc+ll,ipomp+1)
                        itetposa(ibtetposa+m,ipomp+1)=
     &                                    itetposa(ibtetposa+ll,ipomp+1)

                     end do

                        wtxwwp(iswwp+itys,12,ii3) =
     &                  wtxwwp(iswwp+itys,12,ii3) + weit
                        wtxwwp(iswwp+itys,13,ii3) =
     &                  wtxwwp(iswwp+itys,13,ii3) + weight
                        wtxwwp(iswwp+itys,14,ii3) =
     &                  wtxwwp(iswwp+itys,14,ii3) + 1.d+0

                  end if

               end if

*-----------------------------------------------------------------------

  300    continue

*-----------------------------------------------------------------------
*        russian roulette
*-----------------------------------------------------------------------

               m = 0

         do 100 ll = 1, nabov

                        weit = wta(iakwt+ll,ipomp+1)
                        itys = ntya(iaknty+ll,ipomp+1)

                  if( wta(iakwt+ll,ipomp+1) .le. 0.0d0 ) goto 100

                        iupnm = 1

                  if( iwwin(itys) .eq. 0 ) goto 200

                  if( nfcsa(iaknfc+ll,ipomp+1) .lt. 0 .and.
     &              ( ncol .eq. 13 .or. ncol .eq. 14 ) ) then

                        nfcsa(iaknfc+ll,ipomp+1) = 0

                        goto 200

                  else if( nfcsa(iaknfc+ll,ipomp+1) .gt. 0 .and.
     &                     ncol .eq. 10 ) then

                        goto 200

                  end if

*-----------------------------------------------------------------------
*                 weight window for above particle
*-----------------------------------------------------------------------

                     enrg = ea(iake+ll,ipomp+1)
                     timn = ta(iakt+ll,ipomp+1)

                     call awwd(itys,jblz,jlev,jlat,ii3,
     &                         iwxp,iwyp,iwzp,
     &                         enrg,timn,weit,cbup,cbsv,cbdn,jwit)

                     k = 3
                     if( ncol .ne. 10 ) k = 7

*-----------------------------------------------------------------------

                     if( jwit .eq. -10 ) then

                           goto 220

*-----------------------------------------------------------------------

                     else if( jwit .eq. 0 .and. weit .lt. cbdn ) then

                           iupnm = 1
                           if( cbsv .le. 0.0 ) goto 200

                           ri = max( weit / cbsv, 1.d0 / dble( mxspln0))
                           ccsv = weit / ri

                           r  = unirn(dummy)

                        if( r .lt. ri ) then

                           goto 210

                        else

                           goto 220

                        end if

*-----------------------------------------------------------------------

                     else if( jwit .eq. 1 ) then

                           cinimp = wtina(iakwin+ll,ipomp+1)

                        if( ncol .eq. 10 ) then

                           cotimp = aimp(itys,iblz2,ilev2,ilat2,ii4)

                        else

                           cotimp = aimp(itys,iblz1,ilev1,ilat1,ii4)

                        end if

                           rjsuv = wc01(itys) / cotimp * cinimp
                           rjcut = wc02(itys) / cotimp * cinimp

                        if( rjsuv .gt. 0.0d0 ) then

                           if( weit .lt. rjcut ) then

                                 ri = max( weit / rjsuv, 0.5d0 )
                                 ccsv = weit / ri

                                 r  = unirn(dummy)

                              if( r .lt. ri ) then

                                 iupnm = 1

                                 wtxcut(iswct+itys,1,ii4) =
     &                           wtxcut(iswct+itys,1,ii4) + ccsv - weit
                                 wtxcut(iswct+itys,2,ii4) =
     &                           wtxcut(iswct+itys,2,ii4) + 1.d+0

                                 wta(iakwt+ll,ipomp+1) = ccsv

                              else

                                 iupnm = 0

                                 wtxcut(iswct+itys,3,ii4) =
     &                           wtxcut(iswct+itys,3,ii4) + weit
                                 wtxcut(iswct+itys,4,ii4) =
     &                           wtxcut(iswct+itys,4,ii4) + 1.d+0

                                 wta(iakwt+ll,ipomp+1) = 0.0d0

                              end if

                           end if

                        end if

*-----------------------------------------------------------------------

                     end if

                           goto 200

*-----------------------------------------------------------------------

  210                iupnm = 1

                     wtxwwp(iswwp+itys,k+1,ii3) =
     &               wtxwwp(iswwp+itys,k+1,ii3) + ccsv - weit
                     wtxwwp(iswwp+itys,k+2,ii3) =
     &               wtxwwp(iswwp+itys,k+2,ii3) + 1.d+0

                     wta(iakwt+ll,ipomp+1) = ccsv

                     goto 200

*-----------------------------------------------------------------------

  220                iupnm = 0

                     wtxwwp(iswwp+itys,k+3,ii3) =
     &               wtxwwp(iswwp+itys,k+3,ii3) + weit
                     wtxwwp(iswwp+itys,k+4,ii3) =
     &               wtxwwp(iswwp+itys,k+4,ii3) + 1.d+0

                     wta(iakwt+ll,ipomp+1) = 0.0d0

                     goto 200

*-----------------------------------------------------------------------

  200             continue

                  if( iupnm .eq. 1 ) then

                     m = m + 1

                     ea(iake+m,ipomp+1) = ea(iake+ll,ipomp+1)
                     namea(iaknam+m,ipomp+1) = namea(iaknam+ll,ipomp+1)
                     ntya(iaknty+m,ipomp+1) = ntya(iaknty+ll,ipomp+1)
                     nkfa(iaknkf+m,ipomp+1) = nkfa(iaknkf+ll,ipomp+1)

                     ua(iaku+m,ipomp+1) = ua(iaku+ll,ipomp+1)
                     va(iakv+m,ipomp+1) = va(iakv+ll,ipomp+1)
                     wa(iakw+m,ipomp+1) = wa(iakw+ll,ipomp+1)
                     wta(iakwt+m,ipomp+1) = wta(iakwt+ll,ipomp+1)

                     xa(iakx+m,ipomp+1) = xa(iakx+ll,ipomp+1)
                     ya(iaky+m,ipomp+1) = ya(iaky+ll,ipomp+1)
                     za(iakz+m,ipomp+1) = za(iakz+ll,ipomp+1)
                     ta(iakt+m,ipomp+1) = ta(iakt+ll,ipomp+1)

                     iblza(iakblz+m,ipomp+1) = iblza(iakblz+ll,ipomp+1)
                     nmeda(iaknmd+m,ipomp+1) = nmeda(iaknmd+ll,ipomp+1)
                     wtina(iakwin+m,ipomp+1) = wtina(iakwin+ll,ipomp+1)
                     wtnza(iakwnz+m,ipomp+1) = wtnza(iakwnz+ll,ipomp+1)

                     spxa(iakspx+m,ipomp+1) = spxa(iakspx+ll,ipomp+1)
                     spya(iakspy+m,ipomp+1) = spya(iakspy+ll,ipomp+1)
                     spza(iakspz+m,ipomp+1) = spza(iakspz+ll,ipomp+1)

                     nzsta(iakzst+m,ipomp+1) = nzsta(iakzst+ll,ipomp+1)
                     nsosa(iaksos+m,ipomp+1) = nsosa(iaksos+ll,ipomp+1)

                     ncnta(iaknct+1,m,ipomp+1) =
     &                                      ncnta(iaknct+1,ll,ipomp+1)
                     ncnta(iaknct+2,m,ipomp+1) =
     &                                      ncnta(iaknct+2,ll,ipomp+1)
                     ncnta(iaknct+3,m,ipomp+1) =
     &                                      ncnta(iaknct+3,ll,ipomp+1)

                     nfcsa(iaknfc+m,ipomp+1) = nfcsa(iaknfc+ll,ipomp+1)
                     xfcsa(iakxfc+m,ipomp+1) = xfcsa(iakxfc+ll,ipomp+1)
                     itetposa(ibtetposa+m,ipomp+1)=
     &                                   itetposa(ibtetposa+ll,ipomp+1)

                  end if

  100    continue

            nabov = m

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*     error message
*-----------------------------------------------------------------------

 9010 continue

         write(io,12) nabov,maxbn2,nocas
   12    format(/' *** error messge from wwindw ***'
     &   /' storage size (nabov) of particles is greater ',
     &    'than maxbn2.'
     &   /' nabov  =',i10
     &   /' maxbn2 =',i10
     &   /' nocas  =',i10)
         call parastop( 837 )

         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine celimp(ncol)
*                                                                      *
*                                                                      *
*     play the splitting or russian roulette of particle with the cell *
*     importances when particle reach the boundary surface.            *
*                                                                      *
*     2002/04/25 : last revised by K.Niita                             *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use moddas_variance_reduction
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /inout/  in,io

*-----------------------------------------------------------------------

      common /bnkmem/ maxbnk, maxbn2, rtrckflp
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /mulchck/ mulover
!$OMP THREADPRIVATE(/mulchck/)

*-----------------------------------------------------------------------

      common /regdc/ idrg(kvlmax), idgr(kvmmax)
      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      common /impmsg/ iimpn, isimp, iswct, isstr, maxip,
     &                mnimp(7,0:20),
     &                kfimp(7), inimc(7), inimt(7)

      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)


*-----------------------------------------------------------------------

      data eps/1.d-4/

*-----------------------------------------------------------------------
*        for after crossing and icimp is not zero
*        not for weight window
*        not for neutrino
*-----------------------------------------------------------------------

         if( ncol .ne. 10 ) return
         if( icimp(ityp) .eq. 0 ) return

         if( iwwin(ityp) .ne. 0 ) return

         if( abs( ktyp ) .eq. 12 .or. abs( ktyp ) .eq. 14 ) return

*-----------------------------------------------------------------------
*        no process : zero weight or equal of cell importance
*-----------------------------------------------------------------------

         if( wt(ibkwt+no,ipomp+1) .le. 0.0 ) return

*-----------------------------------------------------------------------
*        set up initial data and condition
*-----------------------------------------------------------------------

            nabov  = 0

            cinimp = aimp(ityp,iblz1,ilev1,ilat1,ii1)
            cotimp = aimp(ityp,iblz2,ilev2,ilat2,ii2)
            weit   = wt(ibkwt+no,ipomp+1)

*-----------------------------------------------------------------------
*        error : previous cell has zero importance
*-----------------------------------------------------------------------

            if( cinimp .le. 0.0 ) goto 9001

*-----------------------------------------------------------------------
*     splitting ( for no-void )
*-----------------------------------------------------------------------

            rj = cotimp / wtnz(ibkwnz+no,ipomp+1)

      if( rj .gt. 1.0d0 .and. nmed(ibknmd+no,ipomp+1) .gt. 0 ) then

*-----------------------------------------------------------------------

               weight = weit / rj
               wt(ibkwt+no,ipomp+1) = weight

*-----------------------------------------------------------------------

c S.H. revised to avoid an error of a large rj value (2022.12.16)
              if ( rj .ge. 2147483647d0 ) then
               goto 9010
              else
               mmm = int(rj) - 1
              end if

               rp = rj - dble( mmm + 1 )

               if( rp .ge. eps ) then

                  r = unirn(dummy)
                  if( r .le. rp ) mmm = mmm + 1

               end if

*-----------------------------------------------------------------------

         if( mmm .gt. 0 ) then

               nabov = mmm

               if( nabov .gt. 30 .and. mulover .le. 5) then
                 mulover = mulover + 1
                 write(*,*) "Warning : Particle is split by", nabov,
     &"Check [Importance]."
                 if(mulover .gt. 5) write(*,*) "This message is not show
     &n anymore"
               endif

               if( nabov .gt. maxbn2 ) goto 9010

            do m = 1, nabov

               ntya(iaknty+m,ipomp+1) = nty(ibknty+no,ipomp+1)
               nkfa(iaknkf+m,ipomp+1) = nkf(ibknkf+no,ipomp+1)

               ea(iake+m,ipomp+1)   = ec(ibkec+no,ipomp+1)
               ta(iakt+m,ipomp+1)   = tc(ibktc+no,ipomp+1)
               wta(iakwt+m,ipomp+1) = wt(ibkwt+no,ipomp+1)

               ua(iaku+m,ipomp+1) = u(ibku+no,ipomp+1)
               va(iakv+m,ipomp+1) = v(ibkv+no,ipomp+1)
               wa(iakw+m,ipomp+1) = w(ibkw+no,ipomp+1)

               namea(iaknam+m,ipomp+1) = name(ibknam+no,ipomp+1)

               xa(iakx+m,ipomp+1) = xc(ibkxc+no,ipomp+1) +
     &                              u(ibku+no,ipomp+1) * parz(28)
               ya(iaky+m,ipomp+1) = yc(ibkyc+no,ipomp+1) +
     &                              v(ibkv+no,ipomp+1) * parz(28)
               za(iakz+m,ipomp+1) = zc(ibkzc+no,ipomp+1) +
     &                              w(ibkw+no,ipomp+1) * parz(28)

               iblza(iakblz+m,ipomp+1) = iblz(ibkblz+no,ipomp+1)
               nmeda(iaknmd+m,ipomp+1) = nmed(ibknmd+no,ipomp+1)
               wtina(iakwin+m,ipomp+1) = wtin(ibkwin+no,ipomp+1)
               wtnza(iakwnz+m,ipomp+1) = wtnz(ibkwnz+no,ipomp+1)

               spxa(iakspx+m,ipomp+1) = spx(ibkspx+no,ipomp+1)
               spya(iakspy+m,ipomp+1) = spy(ibkspy+no,ipomp+1)
               spza(iakspz+m,ipomp+1) = spz(ibkspz+no,ipomp+1)

               nzsta(iakzst+m,ipomp+1) = nzst(ibkzst+no,ipomp+1)
               nsosa(iaksos+m,ipomp+1) = nsos(ibksos+no,ipomp+1)

               ncnta(iaknct+1,m,ipomp+1) = ncnt(ibknct+1,no,ipomp+1)
               ncnta(iaknct+2,m,ipomp+1) = ncnt(ibknct+2,no,ipomp+1)
               ncnta(iaknct+3,m,ipomp+1) = ncnt(ibknct+3,no,ipomp+1)

               nfcsa(iaknfc+m,ipomp+1) = nfcs(ibknfc+no,ipomp+1)
               xfcsa(iakxfc+m,ipomp+1) = xfcs(ibkxfc+no,ipomp+1)
              itetposa(ibtetposa+m,ipomp+1)=itetpos(ibtetpos+no,ipomp+1)

            end do

         end if

*-----------------------------------------------------------------------

            wtximp(isimp+ityp,1,ii2) = wtximp(isimp+ityp,1,ii2) + weit
            wtximp(isimp+ityp,2,ii2) = wtximp(isimp+ityp,2,ii2) + weight
            wtximp(isimp+ityp,3,ii2) = wtximp(isimp+ityp,3,ii2) + 1.d+0

      end if

*-----------------------------------------------------------------------
*     russian roulette
*-----------------------------------------------------------------------

            ri = cotimp / cinimp

            if( abs( ri - 1.0d0 ) .le. eps ) return

      if( ri .lt. 1.0d0 ) then

*-----------------------------------------------------------------------

            r = unirn(dummy)

         if( r .lt. ri ) then

            wtximp(isimp+ityp,4,ii2) =
     &      wtximp(isimp+ityp,4,ii2) + weit / ri - wt(ibkwt+no,ipomp+1)
            wtximp(isimp+ityp,5,ii2) =
     &      wtximp(isimp+ityp,5,ii2) + 1.d+0

            wt(ibkwt+no,ipomp+1) = weit / ri

         else

            wt(ibkwt+no,ipomp+1) =  0.0d0
            wtximp(isimp+ityp,6,ii2) = wtximp(isimp+ityp,6,ii2) + weit
            wtximp(isimp+ityp,7,ii2) = wtximp(isimp+ityp,7,ii2) + 1.d+0

         end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*     error message
*-----------------------------------------------------------------------

 9001 continue

         write(io,'(/''Error in cell importance,''/
     &               '' particle is propagated from zero importance'',
     &               '' region ='',i6/
     &               '' particle ='',i3/
     &               '' ncol     ='',i3)') iblz1, ityp, ncol

         call parastop( 838 )

         return

*-----------------------------------------------------------------------

 9010 continue

         write(io,12) nabov,maxbn2,nocas
   12    format(/' *** error messge from celimp ***'
     &   /' storage size (nabov) of particles is greater ',
     &    'than maxbn2.'
     &   /' nabov  =',i10
     &   /' maxbn2 =',i10
     &   /' nocas  =',i10)
         call parastop( 837 )

         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine wtcutof(ncol)
*                                                                      *
*     play the russian roulette of particle by the weight              *
*                                                                      *
*     2002/04/25 : last revised by K.Niita                             *
*                                                                      *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use moddas_variance_reduction
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /inout/  in,io

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)

      common /wparm/  swtm(20), wc1(20), wc2(20) !FURUTA
      common /wparm0/ wc01(20), wc02(20)         !FURUTA
!$OMP THREADPRIVATE(/wparm0/)
      common /impmsg/ iimpn, isimp, iswct, isstr, maxip,
     &                mnimp(7,0:20),
     &                kfimp(7), inimc(7), inimt(7)
      common /fclmsg/ ifcln, isfcl, maxrg,
     &                mnfcl(6,0:20), mrfcl(kvlmax),
     &                kfcls(6), inflc(6), inflt(6)

cKN 2024/03/26
      common /wwindp/ wupn, wsurvn, mxspln, mwhere, mvoww
      common /wwindp0/ wupn0, wsurvn0, mxspln0, mwhere0, mvoww0
!$OMP THREADPRIVATE(/wwindp0/)

      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww


*-----------------------------------------------------------------------
*        for after crossing or reactions
*        no weight cutoff after crossing with forced collision
*-----------------------------------------------------------------------

            if( iwt .eq. 0 .or.
     &        ( ncol .ne. 10 .and. ncol .ne. 13 .and.
     &          ncol .ne. 14 ) ) return

            if( ncol .eq. 10 .and. nfcs(ibknfc+no,ipomp+1) .gt. 0 )
     &                                     return

            if( wt(ibkwt+no,ipomp+1) .le. 0.0d0 ) return

*-----------------------------------------------------------------------
*        weight window particle
*-----------------------------------------------------------------------

            if( iwwin(ityp) .ne. 0 ) goto 500

*-----------------------------------------------------------------------
*        ncol = 14, fcl < 0
*-----------------------------------------------------------------------

         if( ncol .eq. 14 .and. nfcs(ibknfc+no,ipomp+1) .lt. 0 ) then

            nfcs(ibknfc+no,ipomp+1) = 0

            goto 500

         end if

*-----------------------------------------------------------------------
*        set up initial data and condition
*-----------------------------------------------------------------------

               cinimp = wtin(ibkwin+no,ipomp+1)

            if( ncol .eq. 10 ) then

               cotimp = aimp(ityp,iblz2,ilev2,ilat2,ii2)

            else

               cotimp = aimp(ityp,iblz1,ilev1,ilat1,ii2)

            end if

               if( ncol .ne. 10 .and. cotimp .le. 0.0d0 ) goto 9001
               if( cotimp .le. 0.0d0 ) wt(ibkwt+no,ipomp+1) = 0.0d0

               weit = wt(ibkwt+no,ipomp+1)

*-----------------------------------------------------------------------
*        ncol = 10 : russian roulette after boundary crossing
*-----------------------------------------------------------------------

         if( ncol .eq. 10 .and. cotimp .gt. 0.0d0 ) then

                        rjsuv = wc01(ityp) / cotimp * cinimp
                        rjcut = wc02(ityp) / cotimp * cinimp

               if( rjsuv .gt. 0.0d0 ) then

                  if( weit .lt. rjcut ) then

                        ri = weit / rjsuv
                        r  = unirn(dummy)

                     if( r .lt. ri ) then

                        wtxcut(iswct+ityp,1,ii2) =
     &                  wtxcut(iswct+ityp,1,ii2) + rjsuv - weit
                        wtxcut(iswct+ityp,2,ii2) =
     &                  wtxcut(iswct+ityp,2,ii2) + 1.d+0

                        wt(ibkwt+no,ipomp+1) = rjsuv

                     else

                        wtxcut(iswct+ityp,3,ii2) =
     &                  wtxcut(iswct+ityp,3,ii2) + weit
                        wtxcut(iswct+ityp,4,ii2) =
     &                  wtxcut(iswct+ityp,4,ii2) + 1.d+0

                        wt(ibkwt+no,ipomp+1) =  0.0d0

                     end if

                  end if

               end if

         end if

*-----------------------------------------------------------------------
*     ncol = 14
*     with sequential particle
*-----------------------------------------------------------------------

         if( ncol .eq. 14 ) then

                        rjsuv = wc01(ityp) / cotimp * cinimp
                        rjcut = wc02(ityp) / cotimp * cinimp

               if( rjsuv .gt. 0.0 ) then

                  if( weit .lt. rjcut ) then

                        ri = weit / rjsuv
                        r  = unirn(dummy)

                     if( r .lt. ri ) then

                        wtxcut(iswct+ityp,5,ii2) =
     &                  wtxcut(iswct+ityp,5,ii2) + rjsuv - weit
                        wtxcut(iswct+ityp,6,ii2) =
     &                  wtxcut(iswct+ityp,6,ii2) + 1.d+0

                        wt(ibkwt+no,ipomp+1) = rjsuv

                     else

                        wtxcut(iswct+ityp,7,ii2) =
     &                  wtxcut(iswct+ityp,7,ii2) + weit
                        wtxcut(iswct+ityp,8,ii2) =
     &                  wtxcut(iswct+ityp,8,ii2) + 1.d+0

                        wt(ibkwt+no,ipomp+1) =  0.0d0

                     end if

                  end if

               end if

         end if

*-----------------------------------------------------------------------
*     check above bank
*-----------------------------------------------------------------------

  500 continue

      if( nabov .gt. 0 ) then

               m = 0

            do 100 ll = 1, nabov

                     itys  = ntya(iaknty+ll,ipomp+1)

               if( wta(iakwt+ll,ipomp+1) .le. 0.0d0 ) goto 100

                     iupnm = 1

               if( nfcsa(iaknfc+ll,ipomp+1) .lt. 0 .and.
     &           ( ncol .eq. 13 .or. ncol .eq. 14 ) ) then

                     nfcsa(iaknfc+ll,ipomp+1) = 0

                     goto 200

               else if( nfcsa(iaknfc+ll,ipomp+1) .gt. 0 .and.
     &                  ncol .eq. 10 ) then

                     goto 200

               end if

               if( iwwin(ntya(iaknty+ll,ipomp+1)) .ne. 0 ) then

                  if( mwhere0 .eq. 0 .or.
     &              ( mwhere0 .eq. -1 .and. ncol .ne. 10 ) .or.
     &              ( mwhere0 .eq.  1 .and. ncol .eq. 10 ) ) goto 200

               end if

*-----------------------------------------------------------------------

                     cinimp = wtina(iakwin+ll,ipomp+1)

                  if( ncol .eq. 10 ) then

                     cotimp = aimp(itys,iblz2,ilev2,ilat2,ii2)

                  else

                     cotimp = aimp(itys,iblz1,ilev1,ilat1,ii2)

                  end if

                     rjsuv  = wc01(itys) / cotimp * cinimp
                     rjcut  = wc02(itys) / cotimp * cinimp

                  if( rjsuv .le. 0.0 ) goto 200

                     weit = wta(iakwt+ll,ipomp+1)

               if( weit .lt. rjcut ) then

                     ri = weit / rjsuv
                     r  = unirn(dummy)

                  if( r .lt. ri ) then

                     wtxcut(iswct+itys,5,ii2) =
     &               wtxcut(iswct+itys,5,ii2) + rjsuv - weit
                     wtxcut(iswct+itys,6,ii2) =
     &               wtxcut(iswct+itys,6,ii2) + 1.d+0

                     wta(iakwt+ll,ipomp+1) = rjsuv

                     iupnm = 1

                  else

                     wtxcut(iswct+itys,7,ii2) =
     &               wtxcut(iswct+itys,7,ii2) + wta(iakwt+ll,ipomp+1)
                     wtxcut(iswct+itys,8,ii2) =
     &               wtxcut(iswct+itys,8,ii2) + 1.d+0

                     iupnm = 0

                  end if

               end if

*-----------------------------------------------------------------------

  200          continue

               if( iupnm .eq. 1 ) then

                     m = m + 1

                     ea(iake+m,ipomp+1) = ea(iake+ll,ipomp+1)
                     namea(iaknam+m,ipomp+1) = namea(iaknam+ll,ipomp+1)
                     ntya(iaknty+m,ipomp+1) = ntya(iaknty+ll,ipomp+1)
                     nkfa(iaknkf+m,ipomp+1) = nkfa(iaknkf+ll,ipomp+1)

                     ua(iaku+m,ipomp+1) = ua(iaku+ll,ipomp+1)
                     va(iakv+m,ipomp+1) = va(iakv+ll,ipomp+1)
                     wa(iakw+m,ipomp+1) = wa(iakw+ll,ipomp+1)
                     wta(iakwt+m,ipomp+1) = wta(iakwt+ll,ipomp+1)

                     xa(iakx+m,ipomp+1) = xa(iakx+ll,ipomp+1)
                     ya(iaky+m,ipomp+1) = ya(iaky+ll,ipomp+1)
                     za(iakz+m,ipomp+1) = za(iakz+ll,ipomp+1)
                     ta(iakt+m,ipomp+1) = ta(iakt+ll,ipomp+1)

                     iblza(iakblz+m,ipomp+1) = iblza(iakblz+ll,ipomp+1)
                     nmeda(iaknmd+m,ipomp+1) = nmeda(iaknmd+ll,ipomp+1)
                     wtina(iakwin+m,ipomp+1) = wtina(iakwin+ll,ipomp+1)
                     wtnza(iakwnz+m,ipomp+1) = wtnza(iakwnz+ll,ipomp+1)

                     spxa(iakspx+m,ipomp+1) = spxa(iakspx+ll,ipomp+1)
                     spya(iakspy+m,ipomp+1) = spya(iakspy+ll,ipomp+1)
                     spza(iakspz+m,ipomp+1) = spza(iakspz+ll,ipomp+1)

                     nzsta(iakzst+m,ipomp+1) = nzsta(iakzst+ll,ipomp+1)
                     nsosa(iaksos+m,ipomp+1) = nsosa(iaksos+ll,ipomp+1)

                     ncnta(iaknct+1,m,ipomp+1) =
     &                                       ncnta(iaknct+1,ll,ipomp+1)
                     ncnta(iaknct+2,m,ipomp+1) =
     &                                       ncnta(iaknct+2,ll,ipomp+1)
                     ncnta(iaknct+3,m,ipomp+1) =
     &                                       ncnta(iaknct+3,ll,ipomp+1)

                     nfcsa(iaknfc+m,ipomp+1) = nfcsa(iaknfc+ll,ipomp+1)
                     xfcsa(iakxfc+m,ipomp+1) = xfcsa(iakxfc+ll,ipomp+1)
                     itetposa(ibtetposa+m,ipomp+1)=
     &                                  itetposa(ibtetposa+ll,ipomp+1)

               end if

  100       continue

               nabov = m

      end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

 9001 continue

         write(io,'(/''Error in weight cutoff,''/
     &               '' particle is generated in zero importance'',
     &               '' region ='',i6,i6/
     &               '' particle ='',i3/
     &               '' ncol     ='',i3)') iblz1, iblz2, ityp, ncol

         call parastop( 838 )

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine fclsta(ncol,mark,markp)
*                                                                      *
*     control of forced collisions                                     *
*                                                                      *
*     2002/04/25 : last revised by K.Niita                             *
*                                                                      *
*        nfcs(m) =  0 : without forced collisions                      *
*                =  1 : forced collisions for pass through             *
*                =  2 : forced collisions                              *
*                = -1 : normal collision without weight cutoff         *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use moddas_region
      use moddas_variance_reduction
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /inout/  in,io

*-----------------------------------------------------------------------

      common /bnkmem/ maxbnk, maxbn2, rtrckflp
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /fclmsg/ ifcln, isfcl, maxrg,
     &                mnfcl(6,0:20), mrfcl(kvlmax),
     &                kfcls(6), inflc(6), inflt(6)

      dimension     idas(1)
      equivalence ( das, idas )


      dimension dfcl(20)
      dimension ifcl(20)

*-----------------------------------------------------------------------

      data eps/1.d-4/

*-----------------------------------------------------------------------
*     for after crossing or reactions
*-----------------------------------------------------------------------

            if( ifcln .eq. 0 ) return
            if( wt(ibkwt+no,ipomp+1) .le. 0.0d0 ) return

            if( ncol .ne. 10 .and. ncol .ne. 13 .and.
     &          ncol .ne. 14 ) return

*-----------------------------------------------------------------------
*        check particle and region after boundary crossing
*-----------------------------------------------------------------------

            if( ncol .eq. 10 ) then

               if( ifcls(ityp) .eq. 0 .or.
     &             mrfcl(idgr(iblz2)) .eq. 0 ) then

                     nfcs(ibknfc+no,ipomp+1) = 0
                     return

               end if

            end if

*-----------------------------------------------------------------------
*        check particle and region and store the fcl factor
*-----------------------------------------------------------------------

                  do i = 1, 20

                     dfcl(i) = 0.0d0
                     ifcl(i) = 0

                  end do

                     icd = 0

            do 100 j = 1, ifcln

                        icp = 0

                        do l = 1, mnfcl(j,20)

                           if( mnfcl(j,l) .eq. ityp ) icp = 1

                        end do

                        if( ncol .eq. 10 .and. icp .eq. 0 ) goto 100

                        ii = 0

                        idsm = inflt(j)
                        jdsm = 0

                        kdsm = kfcls(j)

                  do i = 1, mnfcl(j,0)

                        jj = 0

                        jdsm = jdsm + 1
                        ntrn = idas_inflt(idsm+jdsm)
                        jdsm = jdsm + 1
                        mtrn = idas_inflt(idsm+jdsm)

                        rfcls = das_kfcls(kdsm-1+i)

                     do 200 k = 1, ntrn

                           ii = ii + 1

                        call tregck(iblz2,ilev2,ilat2,
     &                              mtrn,idas_inflt(idsm+jdsm+1),jj,icc)

                        if( icc .ne. 0 .and.
     &                      abs( rfcls ) .gt. 1.0d-8 ) then

                              icd = icd + 1

                           do l = 1, mnfcl(j,20)

                              dfcl(mnfcl(j,l)) = rfcls
                              ifcl(mnfcl(j,l)) = ii

                              if( ncol .eq. 10 .and.
     &                            mnfcl(j,l) .eq. ityp ) goto 300

                           end do

                              goto 100

                        end if

  200                continue

                        jdsm = jdsm + mtrn

                  end do

  100       continue

                  if( ncol .eq. 10 .or. icd .eq. 0 ) then

                     nfcs(ibknfc+no,ipomp+1) = 0
                     return

                  end if

  300          continue

*-----------------------------------------------------------------------
*        after boundary crossing or ncol = 14
*-----------------------------------------------------------------------

         m = nabov

         if( ncol .eq. 10 .or. ncol .eq. 14 ) then

                  ita = nty(ibknty+no,ipomp+1)

            if( ncol .eq. 14 .and. dfcl(ita) .lt. -1.0d-8 ) then

                     nfcs(ibknfc+no,ipomp+1) = -1

            else if( ncol .eq. 10 .or. dfcl(ita) .gt. 1.0d-8 ) then

                  if( ncol .eq. 10 ) then

                     xc(ibkxc+no,ipomp+1) = xc(ibkxc+no,ipomp+1) +
     &                                    u(ibku+no,ipomp+1) * parz(28)
                     yc(ibkyc+no,ipomp+1) = yc(ibkyc+no,ipomp+1) +
     &                                    v(ibkv+no,ipomp+1) * parz(28)
                     zc(ibkzc+no,ipomp+1) = zc(ibkzc+no,ipomp+1) +
     &                                    w(ibkw+no,ipomp+1) * parz(28)

                     mark  = 1
                     markp = 0

                  end if

                  call fcldis(mark,markp,ipas,prob,probinv,disx,
     &                        nty(ibknty+no,ipomp+1),
     &                        nkf(ibknkf+no,ipomp+1),e(ibke+no,ipomp+1),
     &                        xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                        zc(ibkzc+no,ipomp+1),
     &                        u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                        w(ibkw+no,ipomp+1),
     &                        iblz(ibkblz+no,ipomp+1),
     &                        nmed(ibknmd+no,ipomp+1),
     &                        abs(t(ibkt+no,ipomp+1)),
     &                            wt(ibkwt+no,ipomp+1))

               if( ipas .eq. 1 ) then

                     nfcs(ibknfc+no,ipomp+1) = 1
                     xfcs(ibkxfc+no,ipomp+1) = 1.0d+10

                     r = unirn(dummy)
                     ap = abs( dfcl(ityp) )

                  if( r .le. ap ) then

                     m = m + 1

                     if( m .gt. maxbn2 ) goto 9010

                     nfcsa(iaknfc+m,ipomp+1) = 2
                     xfcsa(iakxfc+m,ipomp+1) = disx

                     ntya(iaknty+m,ipomp+1) = nty(ibknty+no,ipomp+1)
                     nkfa(iaknkf+m,ipomp+1) = nkf(ibknkf+no,ipomp+1)

                     ea(iake+m,ipomp+1)   = ec(ibkec+no,ipomp+1)
                     ta(iakt+m,ipomp+1)   = tc(ibktc+no,ipomp+1)
                     wta(iakwt+m,ipomp+1) = wt(ibkwt+no,ipomp+1)
     &                            * ( 1.0d0 - prob ) / ap
                     if(probinv .gt. 0.d0)
     &                wta(iakwt+m,ipomp+1) =
     &                             wt(ibkwt+no,ipomp+1) * probinv / ap

                     ua(iaku+m,ipomp+1) = u(ibku+no,ipomp+1)
                     va(iakv+m,ipomp+1) = v(ibkv+no,ipomp+1)
                     wa(iakw+m,ipomp+1) = w(ibkw+no,ipomp+1)

                     namea(iaknam+m,ipomp+1) = name(ibknam+no,ipomp+1)

                     xa(iakx+m,ipomp+1) = xc(ibkxc+no,ipomp+1)
                     ya(iaky+m,ipomp+1) = yc(ibkyc+no,ipomp+1)
                     za(iakz+m,ipomp+1) = zc(ibkzc+no,ipomp+1)

                     iblza(iakblz+m,ipomp+1) = iblz(ibkblz+no,ipomp+1)
                     nmeda(iaknmd+m,ipomp+1) = nmed(ibknmd+no,ipomp+1)
                     wtina(iakwin+m,ipomp+1) = wtin(ibkwin+no,ipomp+1)
                     wtnza(iakwnz+m,ipomp+1) = wtnz(ibkwnz+no,ipomp+1)

                     spxa(iakspx+m,ipomp+1) = spx(ibkspx+no,ipomp+1)
                     spya(iakspy+m,ipomp+1) = spy(ibkspy+no,ipomp+1)
                     spza(iakspz+m,ipomp+1) = spz(ibkspz+no,ipomp+1)

                     nzsta(iakzst+m,ipomp+1) = nzst(ibkzst+no,ipomp+1)
                     nsosa(iaksos+m,ipomp+1) = nsos(ibksos+no,ipomp+1)

                     ncnta(iaknct+1,m,ipomp+1)=ncnt(ibknct+1,no,ipomp+1)
                     ncnta(iaknct+2,m,ipomp+1)=ncnt(ibknct+2,no,ipomp+1)
                     ncnta(iaknct+3,m,ipomp+1)=ncnt(ibknct+3,no,ipomp+1)
                     itetposa(ibtetposa+m,ipomp+1)=
     &                                 itetpos(ibtetpos+no,ipomp+1)
                     wtxfcl(isfcl+ityp,1,ifcl(ityp)) =
     &               wtxfcl(isfcl+ityp,1,ifcl(ityp)) +
     &                                           wta(iakwt+m,ipomp+1)
                     wtxfcl(isfcl+ityp,2,ifcl(ityp)) =
     &               wtxfcl(isfcl+ityp,2,ifcl(ityp)) + 1.0d0

                  end if

                     wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * prob

               end if

            end if

         end if

*-----------------------------------------------------------------------
*     nabov > 0 and
*     ncol=10  or  ncol=13,14 and dfcl > 0
*     ncol=13,14 and dfcl < 0 => nfcsa = -1 for no weight cutoff
*-----------------------------------------------------------------------

      if( nabov .gt. 0 ) then

         do i = 1, nabov

                  ita = ntya(iaknty+i,ipomp+1)

            if( ( ncol .eq. 13 .or. ncol .eq. 14 ) .and.
     &            dfcl(ita) .lt. -1.0d-8 ) then

                     nfcsa(iaknfc+i,ipomp+1) = -1

            else if( ncol .eq. 10 .or. dfcl(ita) .gt. 1.0d-8 ) then

                     mark  = 1
                     markp = 0

                  call fcldis(mark,markp,ipas,prob,probinv,disx,
     &                        ntya(iaknty+i,ipomp+1),
     &                        nkfa(iaknkf+i,ipomp+1),ea(iake+i,ipomp+1),
     &                        xa(iakx+i,ipomp+1),ya(iaky+i,ipomp+1),
     &                        za(iakz+i,ipomp+1),
     &                        ua(iaku+i,ipomp+1),va(iakv+i,ipomp+1),
     &                        wa(iakw+i,ipomp+1),
     &                        iblza(iakblz+i,ipomp+1),
     &                        nmeda(iaknmd+i,ipomp+1),
     &                        abs(t(ibkt+no,ipomp+1)),
     &                            wt(ibkwt+no,ipomp+1))

                     mark  = 1
                     markp = 0

               if( ipas .eq. 1 ) then

                     nfcsa(iaknfc+i,ipomp+1) = 1
                     xfcsa(iakxfc+i,ipomp+1) = 1.0d+10

                     r = unirn(dummy)
                     ap = abs( dfcl(ita) )

                  if( r .le. ap ) then

                     m = m + 1

                     if( m .gt. maxbn2 ) goto 9010

                     nfcsa(iaknfc+m,ipomp+1) = 2
                     xfcsa(iakxfc+m,ipomp+1) = disx

                     ntya(iaknty+m,ipomp+1) = ntya(iaknty+i,ipomp+1)
                     nkfa(iaknkf+m,ipomp+1) = nkfa(iaknkf+i,ipomp+1)

                     ea(iake+m,ipomp+1)   = ea(iake+i,ipomp+1)
                     ta(iakt+m,ipomp+1)   = ta(iakt+i,ipomp+1)
                     wta(iakwt+m,ipomp+1) = wta(iakwt+i,ipomp+1)
     &                            * ( 1.0d0 - prob ) / ap
                     if(probinv .gt. 0.d0)
     &                wta(iakwt+m,ipomp+1) =
     &                            wta(iakwt+i,ipomp+1) * probinv / ap

                     ua(iaku+m,ipomp+1) = ua(iaku+i,ipomp+1)
                     va(iakv+m,ipomp+1) = va(iakv+i,ipomp+1)
                     wa(iakw+m,ipomp+1) = wa(iakw+i,ipomp+1)

                     namea(iaknam+m,ipomp+1) = namea(iaknam+i,ipomp+1)

                     xa(iakx+m,ipomp+1) = xa(iakx+i,ipomp+1)
                     ya(iaky+m,ipomp+1) = ya(iaky+i,ipomp+1)
                     za(iakz+m,ipomp+1) = za(iakz+i,ipomp+1)

                     iblza(iakblz+m,ipomp+1) = iblza(iakblz+i,ipomp+1)
                     nmeda(iaknmd+m,ipomp+1) = nmeda(iaknmd+i,ipomp+1)
                     wtina(iakwin+m,ipomp+1) = wtina(iakwin+i,ipomp+1)
                     wtnza(iakwnz+m,ipomp+1) = wtnza(iakwnz+i,ipomp+1)

                     spxa(iakspx+m,ipomp+1) = spxa(iakspx+i,ipomp+1)
                     spya(iakspy+m,ipomp+1) = spya(iakspy+i,ipomp+1)
                     spza(iakspz+m,ipomp+1) = spza(iakspz+i,ipomp+1)

                     nzsta(iakzst+m,ipomp+1) = nzsta(iakzst+i,ipomp+1)
                     nsosa(iaksos+m,ipomp+1) = nsosa(iaksos+i,ipomp+1)

                     ncnta(iaknct+1,m,ipomp+1)=ncnta(iaknct+1,i,ipomp+1)
                     ncnta(iaknct+2,m,ipomp+1)=ncnta(iaknct+2,i,ipomp+1)
                     ncnta(iaknct+3,m,ipomp+1)=ncnta(iaknct+3,i,ipomp+1)
                     itetposa(ibtetposa+m,ipomp+1)=
     &                                     itetposa(ibtetposa+i,ipomp+1)

                  end if

                     wta(iakwt+i,ipomp+1) = wta(iakwt+i,ipomp+1) * prob

               end if

            end if

         end do

      end if

*-----------------------------------------------------------------------

         nabov = m

      return

*-----------------------------------------------------------------------
*     error message
*-----------------------------------------------------------------------

 9010 continue

         write(io,12) m,maxbn2,nocas
   12    format(/' *** error messge from s.fclsta ***'
     &   /' storage size (nabov) of particles is greater ',
     &    'than maxbn2.'
     &   /' nabov  =',i10
     &   /' maxbn2 =',i10
     &   /' nocas  =',i10)
         call parastop( 837 )

         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine fcldis(mark,markp,ipas,prob,probinv,disx,
     &                  ityp,kf,ein,
     &                  xi,yi,zi,ui,vi,wi,
     &                  iblzi,nmedi,timi,wti)
*                                                                      *
*       calculate pass through legth and probability                   *
*       for forcrd collisions                                          *
*       last modified by K.Niita on 24/07/2001                         *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       ipas    : 0=> no pass through, 1=> pass through                *
*       prob    : probability of pass through                          *
*       probinv : 1 - prob (for neutrinos)                             *
*       disx    : flight length sampling for forced collisions         *
*                                                                      *
*                                                                      *
************************************************************************
      use MEMBANKMOD !FURUTA
      use neutrino_mod, only : neutrino_Xsec
      use moddas_material

      use udm_Parameter
      use udm_Utility
      use udm_Manager

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( rpmass = 938.27, rnmass = 939.58 )
      parameter ( nfcdivmax = 1000 ) ! frtati 2023/12/22

*-----------------------------------------------------------------------

      common /eparm/  esmax, esmin, emin(20)
      common /ndemax/ dnmax(20)
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      common /kmat1g/ kmat(kvlmax)
      common /xgeosm/ ksig(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)
      common /regcm/  icmg(kvlmax)

      common /xinels/ ksige, ksign


*-----------------------------------------------------------------------
      common /muint/ imuint, imubrm, imuppd, imucap
      common /sigmuh/ sigmuv(200),sigmub(200),sigmup(200),
     &                nmuh,itzmuh(200),itamuh(200)
!$OMP THREADPRIVATE(/sigmuh/)

*-----------------------------------------------------------------------

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

      common /pnint/ ipnint
      common /gmuppd/ igmuppd ! S.Abe 2019/11/07
      common /xphotn/ totgam,totegs,totpni,totgmm ! y.sakaki 2019/09/26
!$OMP THREADPRIVATE(/xphotn/)

      data rcc / 2.997925d+10 /

c----------------------------------------------------------------------
      common /paraj/ mstz(300), parz(300)

      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)
      common /kmat1ka/ kmathd(kvlmax), kmathe(kvlmax)

      common /emode/  emodem, ge1, ge2, iemode

*-----------------------------------------------------------------------

      common /ccxsm/  icxsni, icxspi

*-----------------------------------------------------------------------
cfrtati 2023/12/22
      real*8, dimension(nfcdivmax) :: totcph, totcpn, pbbin

      integer,save:: iwarn=0 ! T.Sato 2024/09/23

*-----------------------------------------------------------------------
*     initial values and constants
*-----------------------------------------------------------------------

               ipas = 0
               prob = 0.0d0
               probinv = 0.0d0
               disx = 0.0d0

               mat  = nmedi
               ktyp = kf
               jtyp = ichgf(ityp,kf)
               mtyp = ibryf(ityp,kf)
               rtyp = rmtyp(ityp,kf)

*-----------------------------------------------------------------------
*        distance up to boundary
*-----------------------------------------------------------------------

               nmed = nmedi
               iblz = iblzi

               call gomdis(dist,xi,yi,zi,ui,vi,wi,
     &                     mark,markp,nmed,iblz)

                  if( mark .eq. -2 ) then

                     write(*,'(/''*** Geometry error in CG/GG *** ''/
     &               '' in sub. fcldis, cannot find boundary'')')
                     write(*,'( '' mark    ='',i6)') mark
                     write(*,'( '' region  ='',i6)') iblz
                     write(*,'( '' x, y, z :'',3e17.8)')  xi,yi,zi
                     write(*,'( '' u, v, w :'',3e17.8)')  ui,vi,wi

                     return

                  end if

*-----------------------------------------------------------------------
*        range and energy for charged particle
*-----------------------------------------------------------------------

               ecc = ein

            if( jtyp .ne. 0 .and. mat .ne. 0 ) then

               call rainge(ein,rng1,mat,ityp,ktyp,jtyp,rtyp)
cfrtati 2023/12/22
!               if( dist .ge. rng1 ) return
               if( dist .ge. rng1 ) then
!                 if( ityp.eq.1 .or. ityp.eq.15 .or. ityp.eq.18 ) then
                 if( ityp.eq.1 .or.(ityp.ge.15 .and. ityp.le.19)) then ! T.Sato 2024/09/23
                   dist = rng1
                 else
                 iwarn=iwarn+1
       if(iwarn.eq.1) then
        write(*,'("Ignore forced collision because the target,",es12.4,
     &  ", is thicker than the range of particle, ",es12.4)')
     &  dist,rng1
       endif
                   return
                 end if
               end if

               call ecol(ecc,dist,ein,rng1,
     &                   mat,ityp,ktyp,jtyp,rtyp)

               ecc = ( ein + ecc ) / 2.d0

            end if

*-----------------------------------------------------------------------
*     total mean free path from cross sections at average energy
*     Here we do not consider the Delta Ray production 2010/07/30
*-----------------------------------------------------------------------

               totttl = 0.d0
               tothyd = 0.d0
               totdcy = 0.d0

               totgam = 0.d0
               totegs = 0.d0
               totpni = 0.d0
               totgmm = 0.d0 ! y.sakaki 2019/09/26

               totcph = 0.d0 ! frtati 2023/12/22
               totcpn = 0.d0 ! frtati 2023/12/22
               nfcdiv = mstz(162)   ! T.Sato 2024/02/18

*-----------------------------------------------------------------------
*        totdcy : decay flight legth
*-----------------------------------------------------------------------

               totdcy = dcyfl(kf,ecc,rtyp)

               if( totdcy .gt. 1.0d+10 ) goto 1000

*-----------------------------------------------------------------------
*     for no void
*-----------------------------------------------------------------------

      if( mat .gt. 0 ) then

               lem   = nint( dnel_das(kmat0+mat) )
               hydro = denh_das(kmat0+mat)

*-----------------------------------------------------------------------
*        Negative muon capture in material
*-----------------------------------------------------------------------

         if( ityp .eq. 7 .and. imucap .gt. 0 ) then

            totdcy = 0.d0
            dentot = 0.d0

            do i = 1, lem

               itz = idnint( zz_das(kmat(mat)+i) )
               ita = idnint( a_das(kmat(mat)+i) )
               capr = dmax1(0.d0, caprmu(itz,ita)*1.d3)
     &                *den_das(kmat(mat)+i)
               dcyr = dmax1(0.d0, dcyrmu(itz)*1.d3)
     &                *den_das(kmat(mat)+i)
               dentot = dentot + den_das(kmat(mat)+i)
               totdcy = totdcy + capr + dcyr

            enddo

            if( hydro .gt. 0.d0 ) then

               itz = 1
               ita = 1

               capr = dmax1(0.d0, caprmu(itz,ita)*1.d3) * hydro
               dcyr = dmax1(0.d0, dcyrmu(itz)*1.d3) * hydro

               dentot = dentot + hydro
               totdcy = totdcy + capr + dcyr

            endif

            totdcy = totdcy / dentot
     &              / rcc / dsqrt(((ein/rtyp)+2.d0)*(ein/rtyp))

         endif


*-----------------------------------------------------------------------
*        user defined interaction
*-----------------------------------------------------------------------

         if(udm_int_num .gt. 0) then
          !-----------------------------
           totudm    =0d0
           totudm_mul=0d0
          !-----------------------------
           call fill_mat_info(mat)
           do i = 1, udm_int_num
             udm_Kin=ein
             udm_kf_incident=ktyp
             udm_sigt=0d0
             call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
             totudm_this = udm_sigt * denm(mat)
             totudm_(i) = totudm_this
             totudm     = totudm_this               + totudm
             totudm_mul = totudm_this * udm_bias(i) + totudm_mul
           enddo
           totttl = totttl + totudm
         endif

*-----------------------------------------------------------------------
*     non nucleus
*-----------------------------------------------------------------------

      if( ityp .lt. 15 ) then

               tme = timi / 10.0d0
               icl = idgr(iblzi)


*-----------------------------------------------------------------------
*        for proton
*-----------------------------------------------------------------------

            if( ityp .eq. 1 ) then

                  dmaxn = das_kmatd(kmatd(mat)+1)
                  dminn = das_kmatd(kmatd(mat)+2)

*-----------------------------------------------------------------------

                  if( hydro .gt. 0.0d0 ) then

                        dsmax = das_kmatd(kmatd(mat)+6)

                     if( ein .gt. dsmax ) then

                        kf1 = kf
                        kf2 = 2212

cfrtati 2023/12/22
                      deld = dist/dble(nfcdiv)
                      eint = ein
                      do idiv = 1, nfcdiv
                        call rainge(eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
                        call ecol(ecct,deld,eint,rng2,mat,
     &                            ityp,ktyp,jtyp,rtyp)
                        if(ecct .eq. 0.d0) exit ! 2025/3/11 ogawa. Infinite loop without this exit
!                        call sigjam(kf1,kf2,ein,sig,sigel,signo)
                        call sigjam(kf1,kf2,0.5d0*(eint+ecct),
     &                              sig,sigel,signo)

!                        tothyd = tothyd + sig * hydro
                        totcph(idiv) = sig * hydro

                        etint = ecct
                      end do
                     end if

                  end if

*-----------------------------------------------------------------------
*              mixed material
*-----------------------------------------------------------------------

               if( ein .gt. dminn .and. ein .le. dmaxn ) then

                     lema = 0

                     mk = mat
                     rh = denm(mat)

cfrtati 2023/12/22
                deld = dist/dble(nfcdiv)
                eint = ein
                do idiv = 1, nfcdiv
                  call rainge(eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
                  call ecol(ecct,deld,eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
                  if(ecct .eq. 0.d0) exit ! 2025/3/11 ogawa. Infinite loop without this exit
!                  call sig_tot(ityp,sigt,sigaa,ein,mk)
                  call sig_tot(ityp,sigt,sigaa,0.5d0*(eint+ecct),mk)

                  do i = 1, isigc(0)

                        dsmax = das_kmate(kmate(mk)+(i-1)*5+11)

                     if( ein .le. dsmax ) then


                        siggc(i) = siggc(i) * rh

                        lema = lema + 1
                        isigd(lema) = -i
                        siggd(lema) = siggc(i)

!                        totttl = totttl + siggc(i)
                        totcpn(idiv) = totcpn(idiv) + siggc(i)

                     end if

                  end do

*-----------------------------------------------------------------------

                     do i = 1, lem

                        dsmax = das_kmatd(kmatd(mat)+(i-1)*5+11)

                     if( ein .gt. dsmax ) then

                        itz = nint( zz_das(kmat(mat)+i) )
                        ita = nint( a_das(kmat(mat)+i) )

!                        call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)
                        call sigrc(ityp,0.5d0*(eint+ecct),ita,itz,
     &                             sigt,signe,sigel)

                        if( ielas .eq. 0 ) sigel = 0.0

                        signe = signe*den_das(kmat(mat)+i)
                        sigel = sigel*den_das(kmat(mat)+i)
                        sigt  = signe + sigel

                        siggn(ksign+i) = signe
                        sigge(ksige+i) = sigel

                        lema = lema + 1
                        isigd(lema) = i
                        siggd(lema) = sigt

!                        totttl = totttl + sigt
                        totcpn(idiv) = totcpn(idiv) + sigt

                     end if
                     end do

                        isigd(0) = lema

                  eint = ecct
                end do

*-----------------------------------------------------------------------
*              nuclear data
*-----------------------------------------------------------------------

               else if( ein .le. dminn ) then


                     mk = mat
                     rh = denm(mat)

cfrtati 2023/12/22
                   deld = dist/dble(nfcdiv)
                   eint = ein
                   do idiv = 1, nfcdiv
                     call rainge(eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
                     call ecol(ecct,deld,eint,rng2,mat,
     &                         ityp,ktyp,jtyp,rtyp)
!                     call sig_tot(ityp,sigt,sigaa,ein,mk)
                     if(ecct .eq. 0.d0) exit ! 2025/3/11 ogawa. Infinite loop without this exit
                     call sig_tot(ityp,sigt,sigaa,0.5d0*(eint+ecct),mk)

                     do i = 1, isigc(0)

                      if( iemode.eq.0 ) then ! S.H. 2024.12.16
                        siggc(i) = siggc(i) * rh
!                        totttl = totttl + siggc(i)
                        totcpn(idiv) = totcpn(idiv) + siggc(i)
                      else
                        siggcn(i) = siggcn(i) * rh
                        siggce(i) = siggce(i) * rh
                        totcpn(idiv) = totcpn(idiv)
     &                       + siggcn(i) + siggce(i)
                      end if

                     end do
                     eint = ecct
                   end do

*-----------------------------------------------------------------------
*              high energy by sigrc
*-----------------------------------------------------------------------

               else if( ein .gt. dmaxn ) then

cfrtati 2023/12/22
                   deld = dist/dble(nfcdiv)
                   eint = ein
                   do idiv = 1, nfcdiv
                     call rainge(eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
                     call ecol(ecct,deld,eint,rng2,mat,
     &                         ityp,ktyp,jtyp,rtyp)
                     if(ecct .eq. 0.d0) exit ! 2025/3/11 ogawa. Infinite loop without this exit
                     do i = 1, lem

                        itz = nint( zz_das(kmat(mat)+i) )
                        ita = nint( a_das(kmat(mat)+i) )

!                        call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)
                        if(eint+ecct .eq. 0.d0) then
                          signe = 0.d0
                          sigel = 0.d0
                        else
                          call sigrc(ityp,0.5d0*(eint+ecct),ita,itz,
     &                             sigt,signe,sigel)
                        endif

                        if( ielas .eq. 0 ) sigel = 0.0

                        signe = signe*den_das(kmat(mat)+i)
                        sigel = sigel*den_das(kmat(mat)+i)
                        sigt  = signe + sigel

                        siggn(ksign+i) = signe
                        sigge(ksige+i) = sigel

!                        totttl = totttl + sigt
                        totcpn(idiv) = totcpn(idiv) + sigt

                     end do
                     eint = ecct
                   end do

               end if

*-----------------------------------------------------------------------
*        for neutron
*-----------------------------------------------------------------------

            else if( ityp .eq.2 ) then

                  dmaxn = das_kmatd(kmatd(mat)+3)
                  dminn = das_kmatd(kmatd(mat)+4)

*-----------------------------------------------------------------------

                  if( hydro .gt. 0.0d0 ) then

                        dsmax = das_kmatd(kmatd(mat)+7)

                     if( ein .gt. dsmax ) then

                        kf1 = kf
                        kf2 = 2212

                        call sigjam(kf1,kf2,ein,sig,sigel,signo)

                        tothyd = tothyd + sig * hydro

                     end if

                  end if

*-----------------------------------------------------------------------
*              mixed material
*-----------------------------------------------------------------------

               if( ein .gt. dminn .and. ein .le. dmaxn ) then

                     lema = 0

                     mk = mat
                     rh = denm(mat)

                  call xstneu(0,sigt,sigaa,icl,ein,tme,mk)

                  do i = 1, isigc(0)

                        dsmax = das_kmate(kmate(mk)+(i-1)*5+12)

                     if( ein .le. dsmax ) then


                        siggc(i) = siggc(i) * rh

                        lema = lema + 1
                        isigd(lema) = -i
                        siggd(lema) = siggc(i)

                        totttl = totttl + siggc(i)

                     end if

                  end do

*-----------------------------------------------------------------------

                     do i = 1, lem

                        dsmax = das_kmatd(kmatd(mat)+(i-1)*5+12)

                     if( ein .gt. dsmax ) then

                        itz = nint( zz_das(kmat(mat)+i) )
                        ita = nint( a_das(kmat(mat)+i) )

                        call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)

                        if( ielas .eq. 0 ) sigel = 0.0

                        signe = signe*den_das(kmat(mat)+i)
                        sigel = sigel*den_das(kmat(mat)+i)
                        sigt  = signe + sigel

                        siggn(ksign+i) = signe
                        sigge(ksige+i) = sigel

                        lema = lema + 1
                        isigd(lema) = i
                        siggd(lema) = sigt

                        totttl = totttl + sigt

                     end if
                     end do

                        isigd(0) = lema

*-----------------------------------------------------------------------
*              nuclear data
*-----------------------------------------------------------------------

               else if( ein .le. dminn ) then


                     mk = mat
                     rh = denm(mat)

                  call xstneu(0,sigt,sigaa,icl,ein,tme,mk)

*-----------------------------------------------------------------------

                  do i = 1, isigc(0)

                     siggc(i) = siggc(i) * rh
                     totttl = totttl + siggc(i)

                  end do

*-----------------------------------------------------------------------
*              high energy by sigrc
*-----------------------------------------------------------------------

               else if( ein .gt. dmaxn ) then

                     do i = 1, lem

                        itz = nint( zz_das(kmat(mat)+i) )
                        ita = nint( a_das(kmat(mat)+i) )

                        call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)

                        if( ielas .eq. 0 ) sigel = 0.0

                        signe = signe*den_das(kmat(mat)+i)
                        sigel = sigel*den_das(kmat(mat)+i)
                        sigt  = signe + sigel

                        siggn(ksign+i) = signe
                        sigge(ksige+i) = sigel

                        totttl = totttl + sigt

                     end do

               end if

*-----------------------------------------------------------------------
*        photon
*-----------------------------------------------------------------------

            else if( ityp .eq. 14 .and. ein .le. dnmax(14) ) then

               mk = mat
               rh = denm(mat)

*-----------------------------------------------------------------------
*        photon for photonuclear reaction
*-----------------------------------------------------------------------

               if( ipnint .ge. 1 ) then

                  call xstpni_pn(0,sigt,ein,mk)

                  totpni = sigt * rh
                  totttl = totttl + sigt * rh

               endif

*-----------------------------------------------------------------------
*        photon for muon pair production
*-----------------------------------------------------------------------

               if( igmuppd .ge. 1 ) then

                  call xstgmm(0,sigt,ein,mk)

                  totgmm = sigt * rh
                  totttl = totttl + totgmm

               endif

*-----------------------------------------------------------------------
*        photon for nuclear data
*-----------------------------------------------------------------------

               if( iegsemi .eq. 0 ) then
                  call xstgam(0,sigt,sigaa,ein,mk)

                  totgam = sigt * rh
                  totttl = totttl + sigt * rh

*---------------------------------------------------------------------
*        photon by EGS5
*---------------------------------------------------------------------

               else if( iegsemi .ne. 0 ) then

                  call egs5pfpl(tstep,ein,wti,icl,
     &                          totpni,totegs)

                  totttl = totttl + totegs

               endif

*-----------------------------------------------------------------------
*        electron for nuclear data :  non
*-----------------------------------------------------------------------

            else if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &                 ein .le. dnmax(12) ) then

*-----------------------------------------------------------------------
*        for high energy part
*-----------------------------------------------------------------------

            else

*-----------------------------------------------------------------------
*              Hydrogen cross section by JAM
*-----------------------------------------------------------------------

                  sig   = 0.0
                  sigel = 0.0

               if( hydro .gt. 0.0d0 ) then

                  kf1 = kf
                  kf2 = 2212

                  call sigjam(kf1,kf2,ecc,sig,sigel,signo)

               end if

                  tothyd = sig * hydro

*-----------------------------------------------------------------------
*              total, elastic and nonelastic cross sections by sigrc
*-----------------------------------------------------------------------

               if( ityp .le. 2 ) then
*                    --- nuecleons ---

                  do i = 1, lem

                     itz = nint( zz_das(kmat(mat)+i) )
                     ita = nint( a_das(kmat(mat)+i) )

                     call sigrc(ityp,ecc,ita,itz,sigt,signe,sigel)

                     if( ( ityp .eq. 1 .and. ielas .le. 1 ) .or.
     &                   ( ityp .eq. 2 .and. ielas .eq. 0 ) )
     &                     sigel = 0.0

                     signe = signe*den_das(kmat(mat)+i)
                     sigel = sigel*den_das(kmat(mat)+i)
                     sigt  = signe + sigel

                     totttl = totttl + sigt

                  end do

*-----------------------------------------------------------------------
*                       --- muon ---
*-----------------------------------------------------------------------

               else if( ityp .eq. 6 .or. ityp .eq. 7 ) then

                  nmuh = 0

                  do i = 1, lem

                     itz = idnint( zz_das(kmat(mat)+i) )
                     ita = idnint( a_das(kmat(mat)+i) )

                     call sigmu(ein,ita,itz,sigvpi,sigbrm,sigppd)

                     sigvpi = sigvpi*den_das(kmat(mat)+i)
                     sigbrm = sigbrm*den_das(kmat(mat)+i)
                     sigppd = sigppd*den_das(kmat(mat)+i)

                     nmuh = nmuh + 1
                     itzmuh(nmuh) = itz
                     itamuh(nmuh) = ita
                     sigmuv(nmuh) = sigvpi
                     sigmub(nmuh) = sigbrm
                     sigmup(nmuh) = sigppd

                     totttl = totttl + sigvpi + sigbrm + sigppd

                  enddo

                  if( hydro .gt. 0.d0 ) then

                     itz = 1
                     ita = 1

                     call sigmu(ein,ita,itz,sigvpi,sigbrm,sigppd)

                     sigvpi = sigvpi * hydro
                     sigbrm = sigbrm * hydro
                     sigppd = sigppd * hydro

                     nmuh = nmuh + 1
                     itzmuh(nmuh) = itz
                     itamuh(nmuh) = ita
                     sigmuv(nmuh) = sigvpi
                     sigmub(nmuh) = sigbrm
                     sigmup(nmuh) = sigppd

                     totttl = totttl + sigvpi + sigbrm + sigppd

                  endif

*-----------------------------------------------------------------------

               else if( ityp .eq. 12 .or. ityp .eq. 13 ) then

*                    --- electron or positron ---

               else if( ityp .eq. 14 ) then

*                    --- photon ---

               else if( ityp .eq. 11 .and. abs(kf) .le. 16 ) then

                if( abs(kf) .eq. 15 ) then
*                    --- tau lepton ---
                else
*                    --- neutrino ---
                call neutrino_Xsec(kf, ein, hydro, mat, lem, sigt, sigh)
                totttl = totttl + sigt
                tothyd = tothyd + sigh
                endif

*-----------------------------------------------------------------------
*                       --- pion ---  (S.H. revised on 2016.7.11)
*-----------------------------------------------------------------------

               else if( ityp .ge. 3 .and. ityp .le. 5 ) then

                if ( icxspi .eq. 1 ) then ! PHITS original model

                 if ( ityp .eq. 3 ) then
                  pichg = 1d0
                 else if ( ityp .eq. 5 ) then
                  pichg = -1d0
                 else
                  pichg = 0d0
                 end if

*              --- for hydrogen ---
                 sigt = 0d0

                 if( hydro .gt. 0.0d0 ) then

                  zt = 1d0
                  at = 1d0

                  call pionXS(pichg,ecc,at,zt,sigt,signe,sigel,bmax)

                 end if

                 tothyd = sigt * hydro ! consider both elastic and inelastic

*              --- for nucleus ---
                 do i = 1, lem

                  zt = zz_das(kmat(mat)+i)
                  at = a_das(kmat(mat)+i)

                  call pionXS(pichg,ecc,at,zt,sigt,signe,sigel,bmax)

                  signe = signe*den_das(kmat(mat)+i)

                  totttl = totttl + signe ! consider only inelastic

                 end do

                else ! geometrical model

                 do i = 1, lem

                  totttl = totttl + sigg(ksig(mat)+i)

                 end do

                end if

*-----------------------------------------------------------------------

               else

                  do i = 1, lem

                     totttl = totttl + sigg(ksig(mat)+i)

                  end do

               end if

            end if

*-----------------------------------------------------------------------
*     transport particle is nucleus
*-----------------------------------------------------------------------

      else if( ityp .ge. 15 .and. ityp .le. 19 ) then

                  jta = mtyp
                  jtz = jtyp

cfrtati 2023/12/22
        iuselib = 0
        if( ityp.eq.15 .or. ityp.eq.18 ) then
         if(ityp.eq.15) then
          einmevpern = ein/2.0
          dmaxn = das_kmathd(kmathd(mat)+3)
          dminn = das_kmathd(kmathd(mat)+4)
         else if(ityp.eq.18) then
          einmevpern = ein/4.0
          dmaxn = das_kmathd(kmathd(mat)+5)
          dminn = das_kmathd(kmathd(mat)+6)
         endif
         deld = dist/dble(nfcdiv)
         eint = ein
         if( einmevpern .lt. dminn ) then
          mk = mat
          rh = denm(mat)
          do idiv = 1, nfcdiv
           call rainge(eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
           call ecol(ecct,deld,eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
           if(ecct .eq. 0.d0) exit ! 2025/3/11 ogawa. Infinite loop without this exit
           call sig_tot(ityp,sigt,sigaa,0.5d0*(eint+ecct),mk)
           do i = 1, isigc(0)
            siggc(i) = siggc(i) * rh
            totcpn(idiv) = totcpn(idiv) + siggc(i)
           end do
           eint = ecct
          end do
          iuselib = 1
         else if( einmevpern .gt. dminn .and.
     &            einmevpern .le. dmaxn ) then
          mk = mat
          rh = denm(mat)
          do idiv = 1, nfcdiv
           call rainge(eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
           call ecol(ecct,deld,eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
           if(ecct .eq. 0.d0) exit ! 2025/3/11 ogawa. Infinite loop without this exit
           call sig_tot(ityp,sigt,sigaa,0.5d0*(eint+ecct),mk)
           do i = 1, isigc(0)
            if( ityp.eq.15 ) then
             dsmax = das_kmate(kmate(mk)+(i-1)*5+14)
            else if ( ityp.eq.18 ) then
             dsmax = das_kmate(kmate(mk)+(i-1)*5+15)
            end if
            if( einpern.le.dsmax ) then
             siggc(i) = siggc(i) * rh
             totcpn(idiv) = totcpn(idiv) + siggc(i)
            end if
           end do
           do i = 1, lem
            if( ityp.eq.15 ) then
             dsmax = das_kmatd(kmatd(mk)+(i-1)*5+14)
            else if ( ityp.eq.18 ) then
             dsmax = das_kmatd(kmatd(mk)+(i-1)*5+15)
            end if
            if( einpern.gt.dsmax ) then
             zt = zz_das(kmat(mat)+i)
             at = a_das(kmat(mat)+i)
             call sighi(ap,zp,0.5d0*(eint+ecct),at,zt,signe,sigel,bmax)
             if( ielas .le. 2 ) sigel = 0.0
             signe = signe*den_das(kmat(mat)+i)
             sigel = sigel*den_das(kmat(mat)+i)
             sigt  = signe + sigel
             totcpn(idiv) = totcpn(idiv) + sigt
            end if
           end do
           eint = ecct
           iuselib = 1
          end do
         end if
        end if

*-----------------------------------------------------------------------
*              nucleus - proton
*-----------------------------------------------------------------------

                  signe = 0.0
                  sigel = 0.0

            if( hydro .gt. 0.0d0 ) then

cfrtati 2023/12/22
                  itypn = 1
!                  einn  = rpmass / rtyp * ecc

                deld = dist/dble(nfcdiv)
                eint = ein
                do idiv = 1, nfcdiv
                  call rainge(eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
                  call ecol(ecct,deld,eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
                  if(ecct .eq. 0.d0) exit ! 2025/3/11 ogawa. Infinite loop without this exit
                  einnt = 0.5d0*(eint+ecct) * rpmass / rtyp
                  call sigrc(itypn,einnt,jta,jtz,sigt,signe,sigel)
                  totcph(idiv) = signe * hydro
                  eint = ecct
                end do

            end if

!                  tothyd = signe * hydro

*-----------------------------------------------------------------------
*              nucleus - nucleus
*-----------------------------------------------------------------------

                  ap = dble( jta )
                  zp = dble( jtz )

cfrtati 2023/12/22
           if( iuselib.eq.0 ) then
             deld = dist/dble(nfcdiv)
             eint = ein
             do idiv = 1, nfcdiv
               call rainge(eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
               call ecol(ecct,deld,eint,rng2,mat,ityp,ktyp,jtyp,rtyp)
               if(ecct .eq. 0.d0) exit ! 2025/3/11 ogawa. Infinite loop without this exit
               do i = 1, lem

                  zt = zz_das(kmat(mat)+i)
                  at = a_das(kmat(mat)+i)

!                  call sighi(ap,zp,ecc,at,zt,signe,sigel,bmax)
                  call sighi(ap,zp,0.5d0*(eint+ecct),at,zt,
     &                       signe,sigel,bmax)

                  if( ielas .le. 2 ) sigel = 0.0

                  signe = signe*den_das(kmat(mat)+i)
                  sigel = sigel*den_das(kmat(mat)+i)
                  sigt  = signe + sigel

!                  totttl = totttl + sigt
                  totcpn(idiv) = totcpn(idiv) + sigt

               end do
               eint = ecct
             end do
           end if

*-----------------------------------------------------------------------

      end if
      end if

*-----------------------------------------------------------------------

 1000 continue

*-----------------------------------------------------------------------
cfrtati 2023/12/22 charged-particle forced collision for thick target
            if( ityp.eq.1 .or. (ityp.ge.15 .and. ityp.le.19) ) then
              prob = 1.d0
              do idiv = 1, nfcdiv
                prob = prob * exp( - deld*(totcpn(idiv)+totcph(idiv)) )
                pbbin(idiv) = 1.d0 - prob
              end do
              pbbin(:) = pbbin / ( 1.d0 - prob )
              r = unirn(dummy)
              do idiv = 1, nfcdiv
                if( r.le.pbbin(idiv) ) then
                  disx = deld * dble( idiv - 1 + unirn(dummy) )
                  exit
                end if
              end do
              ipas = 1
              return
            end if

*-----------------------------------------------------------------------
*     total mean free path
*-----------------------------------------------------------------------

            if( totttl + tothyd + totdcy .le. 0.0d0 ) then

               return

            else

               totmfp = 1.d0 / ( totttl + tothyd + totdcy )

            end if

*-----------------------------------------------------------------------
*     pass through probability
*-----------------------------------------------------------------------

            if( dist / totmfp .gt. 30.0 ) return

            prob = exp( - dist / totmfp )

            if(dist / totmfp .ne. 0.d0 .and. 1.d0 - prob .le. 1.d-8)then ! small reaction probability
              xx = dist/totmfp ! maclaurin expansion
              prob = 1.d0 - xx + xx**2/2.d0
              probinv = xx - xx**2/2.d0 + xx**3/6.d0
              ur = unirn(dummy)
              disx = totmfp * xx * ur * (1d0 - (1d0 - ur)/2d0 * xx) ! ~= -totmfp*log(1-ur*(1-exp(-xx)))
              ipas = 1
              if( mstz(198) .ne. 0  ) prob = 0.0d0
              return
            else

             if( prob .le. 1.d-8) return ! do not apply forced collision if reaction probability is high

            endif

*-----------------------------------------------------------------------
*     collision point
*-----------------------------------------------------------------------
cKN 2016/08/08 arg->argg

            argg = 1.d0 - unirn(dummy) * ( 1.d0 - prob )

         if( argg .le. 1.0d0 .and. argg .gt. 0.0d0 ) then

            disx = - totmfp  * log( argg )

         else

            return

         end if

            ipas = 1

*-----------------------------------------------------------------------
cKN 2015/12/19 for complete forced collision

      if( mstz(198) .ne. 0  ) prob = 0.0d0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function aimp(ityp,iblz,ilev,ilat,ii)
*                                                                      *
*     give the importance value of the cell                            *
*                                                                      *
*     2001/12/07 : last revised by K.Niita                             *
*                                                                      *
************************************************************************
      use moddas_region

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /impmsg/ iimpn, isimp, iswct, isstr, maxip,
     &                mnimp(7,0:20),
     &                kfimp(7), inimc(7), inimt(7)

      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)

      dimension     idas(1)
      equivalence ( das, idas )

      dimension ilat(5,10)
      dimension jlat(5,10)

      data jlat / 50 * 0 /
      save jlat !FURUTA
!$OMP THREADPRIVATE(jlat)
*-----------------------------------------------------------------------
*        give an importance factor
*-----------------------------------------------------------------------

            aimp = 1.0d0
            rimp = 1.0d0
            ii   = 0

            if( iwt .eq. 0 ) return

            do 100 j = 1, iimpn + 1

                        icp = 0

                        do l = 1, mnimp(j,20)

                           if( mnimp(j,l) .eq. ityp ) icp = 1

                        end do

                        if( icp .eq. 0 ) goto 100

                        ii = 0

                        idsm = inimt(j)
                        jdsm = 0

                        kdsm = kfimp(j)

                  do i = 1, mnimp(j,0) + 1

                        jj = 0

                        jdsm = jdsm + 1
                        ntrn = idas_inimt(idsm+jdsm)
                        jdsm = jdsm + 1
                        mtrn = idas_inimt(idsm+jdsm)

                        rimp = das_kfimp(kdsm-1+i)

                     do k = 1, ntrn

                           ii = ii + 1

                        call tregck(iblz,ilev,ilat,
     &                              mtrn,idas_inimt(idsm+jdsm+1),jj,icc)

                        if( icc .ne. 0 ) goto 300

                     end do

                        jdsm = jdsm + mtrn

                  end do

  100       continue

                  ErrCha = ''
                  ErrID = 'L:3521/R:aimp/F:celimp.f' !W06_001_001
                  call ErrWrite(ErrID,ErrCha)
                  write(*,'(''** Warning : '',
     &                      ''aimp could not fit any cell''/
     &                      '' iblz, ilev ='',2i4)') iblz, ilev
                  return

  300       continue

            aimp = aimp * rimp


            if( ilev .eq. 0 ) return

*-----------------------------------------------------------------------

         do m = 1, ilev

               kreg = ilat(1,m)
               jlev = 0

            do 110 j = 1, iimpn + 1

                        icp = 0

                        do l = 1, mnimp(j,20)

                           if( mnimp(j,l) .eq. ityp ) icp = 1

                        end do

                        if( icp .eq. 0 ) goto 110

                        idsm = inimt(j)
                        jdsm = 0

                        kdsm = kfimp(j)

                  do i = 1, mnimp(j,0) + 1

                        jj = 0

                        jdsm = jdsm + 1
                        ntrn = idas_inimt(idsm+jdsm)
                        jdsm = jdsm + 1
                        mtrn = idas_inimt(idsm+jdsm)

                        rimp = das_kfimp(kdsm-1+i)

                     do k = 1, ntrn

                        call tregck(kreg,jlev,jlat,
     &                              mtrn,idas_inimt(idsm+jdsm+1),jj,icc)

                        if( icc .ne. 0 ) goto 310

                     end do

                        jdsm = jdsm + mtrn

                  end do

  110       continue

                  ErrCha = ''
                  ErrID = 'L:3586/R:aimp/F:celimp.f' !W06_001_002
                  call ErrWrite(ErrID,ErrCha)
                  write(*,'(''** Warning : '',
     &                      ''aimp could not fit upper cell''/
     &                      '' kreg, ilev ='',2i4)') kreg, m


                  return

  310       continue

            aimp = aimp * rimp

         end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine awwd(ityp,iblz,ilev,ilat,ii,
     &                iwxp,iwyp,iwzp,
     &                enrg,timn,weit,cbup,cbsv,cbdn,iwit)
*                                                                      *
*     give the weight window value of the cell                         *
*                                                                      *
*        iwit =  1   : (default) zero weight -> weight cutoff          *
*        iwit =  0   : outside the window                              *
*        iwit = -1   : inside the window                               *
*        iwit = -10  : kill the particle                               *
*                                                                      *
*     2002/12/05 : last revised by K.Niita                             *
*                                                                      *
************************************************************************
      use moddas_region

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

cKN 2024/03/26
      common /wwindp/ wupn, wsurvn, mxspln, mwhere, mvoww
      common /wwindp0/ wupn0, wsurvn0, mxspln0, mwhere0, mvoww0
!$OMP THREADPRIVATE(/wwindp0/)

      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww

      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)

      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      dimension     idas(1)
      equivalence ( das, idas )

      dimension ilat(5,10)

*-----------------------------------------------------------------------
*        give an weight window factor
*-----------------------------------------------------------------------

            rimp = 0.0d0
            iwit = 1
            ii   = 0

            if( iwt .eq. 0 ) return

            do 100 j = 1, iwwdp

                        icp = 0

                        do l = 1, mnwwp(j,20)

                           if( mnwwp(j,l) .eq. ityp ) icp = 1

                        end do

                        if( icp .eq. 0 ) goto 100

                        ien = abs( ienww(j) )

                     if( ien .gt. 0 ) then

                        if( ienww(j) .gt. 0 ) then

                           epat = enrg

                        else

                           epat = timn

                        end if

                        do k = 1, ien

                           if( epat .lt. eenww(j,k) ) goto 200

                        end do

                           goto 100

  200                      ice = k

                     else

                           ice = 1

                     end if

*-----------------------------------------------------------------------
               if( iwmsh .eq. 1 ) then

                        ii = 0

                        idsm = inwwt(j)
                        jdsm = 0

                        kdsm = kfwwp(j)

                  do i = 1, mnwwp(j,0)

                        jj = 0

                        jdsm = jdsm + 1
                        ntrn = idas_inwwt(idsm+jdsm)
                        jdsm = jdsm + 1
                        mtrn = idas_inwwt(idsm+jdsm)

                        rimp = das_kfwwp(kdsm+(ice-1)*kvlmax+i-1)

                     do k = 1, ntrn

                           ii = ii + 1

                        call tregck(iblz,ilev,ilat,
     &                              mtrn,idas_inwwt(idsm+jdsm+1),jj,icc)

                        if( icc .ne. 0 ) goto 300

                     end do

                        jdsm = jdsm + mtrn

                  end do

*-----------------------------------------------------------------------

               else if( iwmsh .eq. 3 ) then

                     inx = iwxnm(j)
                     iny = iwynm(j)
                     inz = iwznm(j)

                     ixyz = inx * iny * inz

                  if( iwxp .lt. 1 .or. iwxp .gt. inx .or.
     &                iwyp .lt. 1 .or. iwyp .gt. iny .or.
     &                iwzp .lt. 1 .or. iwzp .gt. inz ) then

                     rimp = dvals
                     ii = ixyz + 1

                  else

                     icf = iwzp
     &                   + inz * ( iwyp - 1 )
     &                   + iny * inz * ( iwxp -1 )

                     ii = icf

                     rimp = das_kgwwp(kgwwp(j)+(ice-1)*ixyz+icf-1)

                  end if

                     goto 300

*-----------------------------------------------------------------------

               elseif( iwmsh .eq. 4 ) then

                        ic=idgr(iblz)
                        ihelem=iwxp
                        itet0=iwzp-10000
                        idsm = inwwt(j)
                        jdsm = idsm+1
                        nr=mnwwp(j,0)
                        ntrn=idas_inwwt(jdsm+2)
                        ntrn0=idas_inwwt(jdsm+3+ntrn)
                        mtrn=4+ntrn+ntrn0

                        call ttetck(ic,mtrn,idas_inwwt(inwwt(j)+1),
     &                       itet0,ihelem,icf,icc)

                        if(icc.eq.1)then
                         ii=icf
                         kdsm = kfwwp(j)
                         rimp = das_kfwwp(kdsm+(ice-1)*nr+icf-1)
                         goto 300
                        elseif(icc.lt.0)then
                         ii=nr+1
                         rimp = dvals
                         goto 300
                        endif

               end if

*-----------------------------------------------------------------------

  100       continue

                  return

  300       continue

*-----------------------------------------------------------------------

         if( rimp .lt. 0.0d0 ) then

               iwit = -10

         else if( rimp .gt. 0.0d0 ) then

               cbup = rimp * wupn0
               cbsv = rimp * wsurvn0
               cbdn = rimp

            if( weit .ge. cbdn .and. weit .le. cbup ) then

               iwit = -1

            else

               iwit = 0

            end if

         else

               iwit = 1

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine wwxdis(iwctl,dpx,dwwt,xwp,ywp,zwp,uwp,vwp,wwp,
     &                  nx,ny,nz,xm,ym,zm)
*                                                                      *
*       WW for xyz mesh,                                               *
*       iwctl=0 : give xyz id of the present mesh                      *
*       iwctl=1 : give the distance to next mesh and ID of next mesh   *
*       last modified by K.Niita on 2025/01/29                         *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
*-----------------------------------------------------------------------

      implicit double precision( a-h, o-z )

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww

      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)
      common /wwtrs/ iwwtr(6,4), rwwtr(6,13)

      common /wwxyzp/ iwxpp, iwypp, iwzpp
!$OMP THREADPRIVATE(/wwxyzp/)

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)

      dimension xm(nx+1)
      dimension ym(ny+1)
      dimension zm(nz+1)

      dimension ud(3)

      data dmax  /1.0d+19/
      data dmax0 /1.0d+18/

*-----------------------------------------------------------------------

         small = 1.d-14

*-----------------------------------------------------------------------
*        only one xyz mesh
*-----------------------------------------------------------------------

            m = 1

*-----------------------------------------------------------------------

            dpx = dmax

*-----------------------------------------------------------------------
*        final position
*-----------------------------------------------------------------------

            ud(1) = uwp
            ud(2) = vwp
            ud(3) = wwp

*-----------------------------------------------------------------------
*        transform positions
*-----------------------------------------------------------------------

            call trnsxx(xwp,ywp,zwp,xxa,yya,zza,iwwtr(m,4))

            xpp = xxa
            ypp = yya
            zpp = zza

*-----------------------------------------------------------------------
*     initial xyz position
*-----------------------------------------------------------------------

               ixm = 0
               ixc = 0

            if( ud(1) .gt. 0.0d0 ) ixk =  1
            if( ud(1) .eq. 0.0d0 ) ixk =  0
            if( ud(1) .lt. 0.0d0 ) ixk = -1

            if( ixk .ge. 0 ) then

               do i = 1, nx + 1
                  if( xm(i) .gt. xpp + small ) then
                         ixm = i - 1
                         ixc = i - 1
                         goto 32
                  end if
               end do

            else

               ixm = nx + 2
               ixc = nx + 1

               do i = 1, nx + 1
                  if( xm(i) .ge. xpp - small ) then
                         ixm = i
                         ixc = i - 1
                         goto 32
                  end if
               end do

            end if

   32       continue

*-----------------------------------------------------------------------

               iym = 0
               iyc = 0

            if( ud(2) .gt. 0.0d0 ) iyk =  1
            if( ud(2) .eq. 0.0d0 ) iyk =  0
            if( ud(2) .lt. 0.0d0 ) iyk = -1

            if( iyk .ge. 0 ) then

               do i = 1, ny + 1
                  if( ym(i) .gt. ypp + small ) then
                         iym = i - 1
                         iyc = i - 1
                         goto 35
                  end if
               end do

            else

               iym = ny + 2
               iyc = ny + 1

               do i = 1, ny + 1
                  if( ym(i) .ge. ypp - small ) then
                         iym = i
                         iyc = i - 1
                         goto 35
                  end if
               end do

            end if

   35       continue

*-----------------------------------------------------------------------

               izm = 0
               izc = 0

            if( ud(3) .gt. 0.0d0 ) izk =  1
            if( ud(3) .eq. 0.0d0 ) izk =  0
            if( ud(3) .lt. 0.0d0 ) izk = -1

            if( izk .ge. 0 ) then

               do i = 1, nz + 1
                  if( zm(i) .gt. zpp + small ) then
                         izm = i - 1
                         izc = i - 1
                         goto 38
                  end if
               end do

            else

               izm = nz + 2
               izc = nz + 1

               do i = 1, nz + 1
                  if( zm(i) .ge. zpp - small ) then
                         izm = i
                         izc = i - 1
                         goto 38
                  end if
               end do

            end if

   38       continue

*-----------------------------------------------------------------------
*        present postion
*-----------------------------------------------------------------------
cKN 2025/01/29

               iwxpp = ixc
               iwypp = iyc
               iwzpp = izc

            if( iwctl .eq. 0 ) return

*-----------------------------------------------------------------------
*        dx,dy,dz : distance to the nearest xyz mesh
*-----------------------------------------------------------------------

            dx = dmax
            dy = dmax
            dz = dmax

         if( ixk .ne. 0 ) then
            if( ixm + ixk .ge. 1 .and.
     &         ixm + ixk .le. nx + 1 ) then

               dx = abs( ( xm(ixm+ixk) - xpp ) / ud(1) )

            end if
         end if

         if( iyk .ne. 0 ) then
            if( iym + iyk .ge. 1 .and.
     &          iym + iyk .le. ny + 1 ) then

               dy = abs( ( ym(iym+iyk) - ypp ) / ud(2) )

            end if
         end if

         if( izk .ne. 0 ) then
            if( izm + izk .ge. 1 .and.
     &          izm + izk .le. nz + 1 ) then

               dz = abs( ( zm(izm+izk) - zpp ) / ud(3) )

            end if
         end if

*-----------------------------------------------------------------------
*        which boundary is the nearlist and distance to the boundary
*-----------------------------------------------------------------------

               dd = dx
               jj = 1

            if( dy .lt. dd ) then

               dd = dy
               jj = 2

            end if

            if( dz .lt. dd ) then

               dd = dz
               jj = 3

            end if

            if( dd .gt. dwwt ) return

               dpx = dd

*-----------------------------------------------------------------------
*        next position
*-----------------------------------------------------------------------

            if( jj .eq. 1 ) then
               ixc = ixc + ixk
            else if( jj .eq. 2 ) then
               iyc = iyc + iyk
            else if( jj .eq. 3 ) then
               izc = izc + izk
            end if

*-----------------------------------------------------------------------
cKN 2025/01/29

         if( ixc .lt. 1 .or. ixc .gt. nx .or.
     &       iyc .lt. 1 .or. iyc .gt. ny .or.
     &       izc .lt. 1 .or. izc .gt. nz ) then

               dpx = dmax
               return

         else

               iwxpp = ixc
               iwypp = iyc
               iwzpp = izc

         end if

*-----------------------------------------------------------------------

      return
      end

