************************************************************************
*                                                                      *
      subroutine inclin(ityp,eein,mmas,mchg)
*                                                                      *
*                                                                      *
*       control routine of incl calculation                            *
*       modified by S.Hashimoto on 2016/06/30                          *
*                                                                      *
*        call subroutine : init_incl, icascl, icascl45, parastop       *
*                          levdatset, levset1-5                        *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       inclg   : =1 : for p,n,pi+,pi-,d,t,3He,a                       *
*                 =2 : for p,n,pi+,pi-                                 *
*       inclv   : =0 : INCL version 4.6                                *
*       ityp    : particle type of projectile                          *
*       eein    : energy of projectile (MeV)                           *
*       mmas    : mass of target                                       *
*       mchg    : charge of target                                     *
*                                                                      *
*       conversion of the type of the projectile:                      *
*          ityp  --(kindp)-->  iprojl  --(kindj)-->  f(7)              *
*                                                                      *
*       table of particle type in phits, (lahet,) and incl             *
*          particle            ityp     iprojl     f(7)                *
*                             (phits)   (lahet)   (incl)               *
*          proton               1          1         1                 *
*          neutron              2          2         2                 *
*          pi+                  3          3         3                 *
*          pi0                  4          4         4                 *
*          pi-                  5          5         5                 *
*          mu+                  6          6         -                 *
*          mu-                  7          7         -                 *
*          K+                   8         13         -                 *
*          K0                   9         14(K0L)    -                 *
*          K-                  10         16         -                 *
*                                                                      *
*          other particles     11          -         -                 *
*                                                                      *
*          electron            12         19         -                 *
*          positron            13         20         -                 *
*          photon              14         12         -                 *
*                                                                      *
*          deuteron            15          8         6                 *
*          triton              16          9         7                 *
*          3He                 17         10         8                 *
*          Alpha               18         11         9                 *
*          nucleus             19         23         -                 *
*                                                                      *
*     output:                                                          *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclst   : total number of out going particles and nuclei      *
*                                                                      *
*        iclust(nclst)                                                 *
*                                                                      *
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclust(i,nclst)                                               *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ityp, see above                                *
*                  = 4,                                                *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        qclust(i,nclst)                                               *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, px (GeV/c)                                     *
*                  = 2, py (GeV/c)                                     *
*                  = 3, pz (GeV/c)                                     *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight change                                  *
*                  = 9, delay time                                     *
*                  = 10, x-displace                                    *
*                  = 11, y-displace                                    *
*                  = 12, z-displace                                    *
*                                                                      *
************************************************************************

      use levdat
      use NGSDATAMOD, only : bindeg

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param-physcnst.inc'

*-----------------------------------------------------------------------

      logical lflg
      common /cincl/  inclg, inclv
      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

*-----------------------------------------------------------------------

      integer zmat,amat,nbmat
      real*8 bmax_geo,abund
      common/incmat/zmat(500),amat(500),bmax_geo(6,500),abund(500),nbmat
!$OMP THREADPRIVATE(/incmat/)
      data amat(1)/0/
      data zmat(1)/0/
      real*8 r0,adif,rmaxws,drws,bmax,xfoisa,rlim
      integer nosurf, npaulstr
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,BMAX,XFOISA,RLIM(500),
     &          NPAULSTR,NOSURF
!$OMP THREADPRIVATE(/WS/)
      save init_graine
      data INIT_GRAINE/1/
!$OMP THREADPRIVATE(INIT_GRAINE)

*-----------------------------------------------------------------------

      integer :: kindp(20)
      DATA (kindp(lpt),lpt=1,20)/
     & 1,2,3,4,5,6,7,13,14,16,0,19,20,12,8,9,10,11,23,0/
      integer :: kindj(23)
      DATA (kindj(lpt),lpt=1,23)/
     & 1,2,3,4,5,0,0,6,7,8,9,0,0,0,0,0,0,0,0,0,0,0,0/

      integer ptype(300)
      real*8 :: epi(300),alphai(300),betai(300),gammai(300),tvt(300)
      integer zremi,aremi,jrem, jj
      real*8 esremi,erecremi,alremi,beremi,garemi,
     &     xjrem,yjrem,zjrem,t_fin
      real*8 bimpacti
      real*8 keproj,datarg,dztarg

      integer :: ncincl,nzincl
      common /clincl/ncincl,nzincl(165,112)
!$OMP THREADPRIVATE(/clincl/)
      real*8 f
      integer icoup
      common/calincl/f(30),icoup
!$OMP THREADPRIVATE(/calincl/)

      integer LBLNN,LNNEF
      real*8 cutnn,alph,tfl,pfl,v0l,s1nn,s2nn,v2p,v2n
      COMMON/OPTI/CUTNN,ALPH,
     s  TFL(4),PFL(4),V0L(4),S1NN,S2NN,V2P,V2N,LBLNN,LNNEF
!$OMP THREADPRIVATE(/OPTI/)

      common/inclphits/
     1     epi,alphai,betai,gammai,tvt,
     2     esremi,erecremi,alremi,beremi,garemi,
     3     bimpacti,xjrem,yjrem,zjrem,t_fin,
     4     ptype,zremi,aremi,jrem,ichcoul,nopart
!$OMP THREADPRIVATE(/inclphits/)
      data ncincl/0/
      logical generated
      common/crgauss/generated
!$OMP THREADPRIVATE(/crgauss/)

      common /crshi/  bplus, icrhi, ijudg, imadj, iqmax
      integer :: LCOUL_iadeu
      common /ciadeu/  LCOUL_iadeu
!$OMP THREADPRIVATE(/ciadeu/)

*-----------------------------------------------------------------------

      common /cidwba/ idwba

      integer iflag_incl45
      data iflag_incl45/0/
      save iflag_incl45
!$OMP THREADPRIVATE(iflag_incl45)
      integer icheckn1,icheckn2,icheckn3
      data icheckn1/0/,icheckn2/0/,icheckn3/0/
      save icheckn1,icheckn2,icheckn3
!$OMP THREADPRIVATE(icheckn1,icheckn2,icheckn3)
      integer icheckp1a,icheckp1b,icheckp2a,icheckp2b
      data icheckp1a/0/,icheckp1b/0/,icheckp2a/0/,icheckp2b/0/
      save icheckp1a,icheckp1b,icheckp2a,icheckp2b
!$OMP THREADPRIVATE(icheckp1a,icheckp1b,icheckp2a,icheckp2b)
      integer icheckp3a,icheckp3b
      data icheckp3a/0/,icheckp3b/0/
      save icheckp3a,icheckp3b
!$OMP THREADPRIVATE(icheckp3a,icheckp3b)
      integer icheckp1c
      data icheckp1c/0/
      save icheckp1c
!$OMP THREADPRIVATE(icheckp1c)

*-----------------------------------------------------------------------

      dimension rmspat(20)
      data rmspat/ 0.93827d0, 0.93958d0, 0.1396d0, 0.1350d0, 0.1396d0,
     &             0.0d0, 0.0d0, 0.4936d0, 0.4977d0, 0.4936d0, 4*0.0d0,
     &             1.8756d0, 2.8089d0, 2.8084d0, 3.7274d0, 0.0d0, 0.0d0/

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*        initial values
*-----------------------------------------------------------------------

      generated=.false.

      iprojl=kindp(ityp)
      keproj=eein
      datarg=real(mmas)
      dztarg=real(mchg)
      if ( iprojl .le. 0 ) go to 340

*-----------------------------------------------------------------------

      F(1)=datarg
      F(2)=dztarg
      F(3)=keproj
      F(4)=0d0
      F(5)=45.d0
      F(6)=1.d0
      F(7)=dble(kindj(iprojl))
      F(8)=0d0
      F(9)=1d0
      F(10)=0d0
      F(11)=3d0
      F(12)=6.83d0
      F(13)=3d0
      F(14)=0d0
      F(15)=1d0
      F(16)=1d0
      F(17)=0d0
      F(18)=3d0
      F(19)=1d0
      F(20)=2d0
      F(21)=1d0

*-----------------------------------------------------------------------
      iflag_incl45 = 0

      if ( idwba .eq.1 .and. ityp .eq. 15 ) then

* for 6,7Li(d,n)7,8Be and 6,7Li(d,p)7,8Li at 10 - 50 MeV
         if (mchg .eq. 3 .and. (mmas .eq. 6 .or. mmas .eq. 7)
     1        .and. eein .ge. 10d0 .and. eein .le. 50d0) then
            iflag_incl45 = 1

* for 9Be(d,n)10B and 9Be(d,p)10Be at 5 - 25 MeV
         else if (mchg .eq. 4 .and. mmas .eq. 9
     1           .and. eein .ge. 5d0 .and. eein .le. 25d0) then
            iflag_incl45 = 1

* for 12,13C(d,n)13,14N and 12,13C(d,p)13,14C at 10 - 50 MeV
         else if (mchg .eq. 6 .and. (mmas .eq. 12 .or. mmas .eq. 13)
     1           .and. eein .ge. 10d0 .and. eein .le. 50d0) then
            iflag_incl45 = 1

         end if

      end if
*-----------------------------------------------------------------------

      nosurf=-2
      xfoisa=8.d0
      NPAULSTR=0
      ichcoul=0
      cutnn=1880.d0
      s1nn=1910.d0
      s2nn=1911.d0
      NBMAT=1
      if ( amat(1) .ne. mmas ) then
         amat(1)=mmas
         INIT_GRAINE=1
      end if
      if ( zmat(1) .ne. mchg ) then
         zmat(1)=mchg
         INIT_GRAINE=1
      end if

*-----------------------------------------------------------------------

      eppin = eein / 1000.0d0
      rmsin = rmspat(ityp)
      if ( rmsin .le. 0.d0 ) go to 340

      if( ityp .eq. 1 ) then

         masim = mmas + 1
         mchim = mchg + 1
         einad = 0.0d0

      else if( ityp .eq. 2 ) then

         masim = mmas + 1
         mchim = mchg
         einad = 0.0d0

      else if( ityp .eq. 3 ) then

         masim = mmas
         mchim = mchg + 1
         einad = rstms(3)*1d3

      else if( ityp .eq. 4 ) then

         masim = mmas
         mchim = mchg
         einad = rstms(4)*1d3

      else if( ityp .eq. 5 ) then

         masim = mmas
         mchim = mchg - 1
         einad = rstms(5)*1d3

      else if( ityp .eq. 15 ) then

         masim = mmas + 2
         mchim = mchg + 1
         einad = -2.225d0

      else if( ityp .eq. 16 ) then

         masim = mmas + 3
         mchim = mchg + 1
         einad = -8.483d0

      else if( ityp .eq. 17 ) then

         masim = mmas + 3
         mchim = mchg + 2
         einad = -7.719d0

      else if( ityp .eq. 18 ) then

         masim = mmas + 4
         mchim = mchg + 2
         einad = -28.297d0

      else

         goto 340

      end if

*-----------------------------------------------------------------------
      if ( icrhi .eq. 3 .and. ityp .eq. 15 ) then

         Ecm = eein * datarg / ( 2d0 + datarg )

         rmsp  = radius(2d0,1d0)
         rmst  = radius(datarg,dztarg)
         factr = sqrt(5./3.)    ! = 1.29
         rp  = factr * rmsp
         rt  = factr * rmst

         Rcoul = rp + rt
         Bcoul = 1.44d0 * 1d0 * dztarg / Rcoul
         Eth = Bcoul * 0.8d0

         if ( Ecm .lt. Eth ) then
            LCOUL_iadeu = 1
         else
            LCOUL_iadeu = 0
         end if

      end if

*-----------------------------------------------------------------------
*        INCL
*-----------------------------------------------------------------------

      if ( INIT_GRAINE .eq. 1) then
         call INIT_INCL(INIT_GRAINE)
         INIT_GRAINE=0
      end if

      jj = 0
      ierrincl = 0

      if ( iflag_incl45 .eq. 1  ) then
         F(11)=1d0
         F(14)=1d0
         F(19)=0d0
         F(21)=1d0
         call icascl45(iprojl,keproj,datarg,dztarg,
     &                 ierrincl)

      else

         call icascl(iprojl,keproj,datarg,dztarg,jj,
     &                 ierrincl)

      end if

      if ( ierrincl .eq. 1 ) then
         nclst = -1
         return
      end if

*-----------------------------------------------------------------------
*        booking of the result of INCL
*-----------------------------------------------------------------------

               nbart = 0
               nchat = 0

               npipo = 0
               nping = 0
               npine = 0

               sume = 0.0d0
               tbene = 0.0d0

               poutx = 0.0d0
               pouty = 0.0d0
               poutz = 0.0d0

*-----------------------------------------------------------------------
*        pseudo collision
*-----------------------------------------------------------------------

               nclst = 0

         if( nopart .lt. 0 ) then

               nclst = -1

               return

         end if

         if( nopart .gt. 0 ) then
            do n = 1, nopart
               if( epi(n) .lt. 0.0d0 ) then
                  nclst = -1
                  return
               end if
            end do
         end if

*-----------------------------------------------------------------------


      nlevel = 0

      call levset(mmas , mchg, lflg)

      if (.not. lflg ) then
       elevel1st = 1d0 ! 1MeV is assumed in the unknown case

      else if ( ubound(elevel,1) .gt. 0 ) then
       elevel1st = elevel(1)*1d-3 ! keV -> MeV
       call levUNset
      else
       elevel1st = 1.0d0 ! No bound level TEMPORARY FIX 20201007 S.H.
      end if

      if ( nopart.eq.1 ) then
       if ( (ityp.le.5 .and. ityp.eq.ptype(1))
     &  .or. (ityp.ge.15 .and. ityp.le.18 .and. ityp.eq.ptype(1)+9) )
     &  then ! for p,n,pi+,pi0,pi-, or d,t,3he,alpha

        if( eein-epi(1) .lt. elevel1st ) then
         nclst = -1
         return
        end if

       end if
      end if

*-----------------------------------------------------------------------
*        real collision
*-----------------------------------------------------------------------

         if( nopart .gt. 0 ) then

            do n = 1, nopart

                  nclst = nclst + 1

                  lk = ptype(n)

               if( lk .eq. 1 ) then

                  kf    = 2212
                  ibary = 1
                  ipid  = 1
                  ippad = 1
                  ipprt = 1
                  ipneu = 0
                  ipchg = 1

                  bene = 0.0d0
                  rms = rstms(1)

               else if( lk .eq. 2 ) then

                  kf    = 2112
                  ibary = 1
                  ipid  = 2
                  ippad = 2
                  ipprt = 0
                  ipneu = 1
                  ipchg = 0

                  bene = 0.0d0
                  rms = rstms(2)

               else if( lk .eq. 3 ) then

                  npipo = npipo + 1

                  kf    = 211
                  ibary = 0
                  ipid  = 3
                  ippad = 3
                  ipprt = 0
                  ipneu = 0
                  ipchg = 1

                  bene = 0.0d0
                  rms = rstms(3)

               else if( lk .eq. 4 ) then

                  npine = npine + 1

                  kf    = 111
                  ibary = 0
                  ipid  = 3
                  ippad = 4
                  ipprt = 0
                  ipneu = 0
                  ipchg = 0

                  bene = 0.0d0
                  rms = rstms(4)

               else if( lk .eq. 5 ) then

                  nping = nping + 1

                  kf    = -211
                  ibary = 0
                  ipid  = 3
                  ippad = 5
                  ipprt = 0
                  ipneu = 0
                  ipchg = -1

                  bene = 0.0d0
                  rms = rstms(5)

               else if( lk .eq. 6 ) then

                  kf    = 1000002
                  ibary = 2
                  ipid  = 0
                  ippad = 15
                  ipprt = 1
                  ipneu = 1
                  ipchg = 1

                  bene = bindeg(ipprt,ipneu)
                  rms = 1.8756d0

               else if( lk .eq. 7 ) then

                  kf    = 1000003
                  ibary = 3
                  ipid  = 0
                  ippad = 16
                  ipprt = 1
                  ipneu = 2
                  ipchg = 1

                  bene = bindeg(ipprt,ipneu)
                  rms = 2.8089d0

               else if( lk .eq. 8 ) then

                  kf    = 2000003
                  ibary = 3
                  ipid  = 0
                  ippad = 17
                  ipprt = 2
                  ipneu = 1
                  ipchg = 2

                  bene = bindeg(ipprt,ipneu)
                  rms = 2.8084d0

               else if( lk .eq. 9 ) then

                  kf    = 2000004
                  ibary = 4
                  ipid  = 0
                  ippad = 18
                  ipprt = 2
                  ipneu = 2
                  ipchg = 2

                  rms = 3.7274d0
                  bene = bindeg(ipprt,ipneu)

               else if( lk .le. 128 ) then

                  ipid  = 0
                  ipprt = mod(lk,10)
                  ipneu = (lk-ipprt)/10 + 0.1 - ipprt
                  ippad = 19
                  ipchg = ipprt
                  ibary = ipprt + ipneu
                  kf    = 1000000 * ipprt + ibary

                  bene = bindeg(ipprt,ipneu)
                  rms = rstms(1) * ipprt + rstms(2) * ipneu
     &                - bene / 1000.0d0

               else

                  goto 340

               end if

*-----------------------------------------------------------------------

               if ( iflag_incl45 .eq. 1  ) then

               if ( ityp .eq. 15 .and. lk .eq. 2 ) then
                  derE = eein/2d0-epi(n)
                  if ( derE .ge. 1.1d0 .and. derE .le. 3.1d0 ) then

                     if ( gammai(n) .gt. 0.9999939d0 ) then
                        icheckn1=icheckn1+1
                        if ( icheckn1 .eq. 2500 ) then
                           icheckn1=0
                           goto 999
                        end if
                        nclst = -1
                        return
                     end if

                  else if ( derE .gt. 3.1d0 .and. derE .le. 4.2d0 ) then

                     if ( gammai(n) .gt. 0.9999939d0 ) then
                        icheckn2=icheckn2+1
                        if ( icheckn2 .eq. 10000 ) then
                           icheckn2=0
                           goto 999
                        end if
                        nclst = -1
                        return
                     end if

                  else if ( derE .gt. 4.2d0 .and. derE .lt. 6.0d0 ) then

                     if ( gammai(n) .gt. 0.9999939d0 ) then
                        icheckn3=icheckn3+1
                        if ( icheckn3 .eq. 500 ) then
                           icheckn3=0
                           goto 999
                        end if
                        nclst = -1
                        return
                     end if

                  end if
               end if
 999           continue
               if ( ityp .eq. 15 .and. lk .eq. 1 ) then
                  derE = eein/2d0-epi(n)

* for 6,7Li(d,p)7,8Li
                  if (mchg .eq. 3
     1                 .and. (mmas .eq. 6 .or. mmas .eq. 7) ) then

                     if ( derE .ge. 1.d0 .and. derE .le. 1.8d0 ) then
                        if ( gammai(n) .gt. 0.9998476d0 ) then
                           icheckp1a=icheckp1a+1
                           if ( icheckp1a .eq. 40 ) then
                              icheckp1a=0
                              goto 998
                           end if
                           nclst = -1
                           return
                        else if ( gammai(n) .gt. 0.9993908d0 ) then
                           icheckp1b=icheckp1b+1
                           if ( icheckp1b .eq. 8 ) then
                              icheckp1b=0
                              goto 998
                           end if
                           nclst = -1
                           return
                        end if

                     else if ( derE .ge. 3.6d0
     1                       .and. derE .le. 4.3d0 ) then
                        if ( gammai(n) .gt. 0.9998476d0 ) then
                           icheckp2a=icheckp2a+1
                           if ( icheckp2a .eq. 500 ) then
                              icheckp2a=0
                              goto 998
                           end if
                           nclst = -1
                           return
                        else if ( gammai(n) .gt. 0.9993908d0 ) then
                           icheckp2b=icheckp2b+1
                           if ( icheckp2b .eq. 50 ) then
                              icheckp2b=0
                              goto 998
                           end if
                           nclst = -1
                           return
                        end if

                     else if ( derE .gt. 1.8d0
     1                       .and. derE .lt. 3.6d0 ) then
                        if ( gammai(n) .gt. 0.9996573d0 ) then
                           icheckp3a=icheckp3a+1
                           if ( icheckp3a .eq. 8 ) then
                              icheckp3a=0
                              goto 998
                           end if
                           nclst = -1
                           return
                        else if ( gammai(n) .gt. 0.9993908d0 ) then
                           icheckp3b=icheckp3b+1
                           if ( icheckp3b .eq. 1 ) then
                              icheckp3b=0
                              goto 998
                           end if
                           nclst = -1
                           return
                        end if
                     end if


* for 9Be(d,p)10Be
                  else if (mchg .eq. 4 .and. mmas .eq. 9 ) then

                     if ( derE .ge. 2.1d0 .and. derE .le. 2.2d0 ) then
                        if ( gammai(n) .gt. 0.9993908d0 ) then
                           icheckp1a=icheckp1a+1
                           if ( icheckp1a .eq. 850 ) then
                              icheckp1a=0
                              goto 998
                           end if
                           nclst = -1
                           return
                        else if ( gammai(n) .gt. 0.9961947d0 ) then
                           icheckp1b=icheckp1b+1
                           if ( icheckp1b .eq. 400 ) then
                              icheckp1b=0
                              goto 998
                           end if
                           nclst = -1
                           return
                        else if ( gammai(n) .gt. 0.9659258d0 ) then
                           icheckp1c=icheckp1c+1
                           if ( icheckp1c .eq. 10 ) then
                              icheckp1c=0
                              goto 998
                           end if
                           nclst = -1
                           return
                        end if
                     end if


* for 12,13C(d,p)13,14C
                  else if (mchg .eq. 6
     1                 .and. (mmas .eq. 12 .or. mmas .eq. 13) ) then

                     if ( derE .ge. 1.1d0 .and. derE .le. 1.7d0 ) then
                        if ( gammai(n) .gt. 0.9993908d0 ) then
                           icheckp1a=icheckp1a+1
                           if ( icheckp1a .eq. 13 ) then
                              icheckp1a=0
                              goto 998
                           end if
                           nclst = -1
                           return
                        end if

                     else if ( derE .ge. 2.5d0
     1                       .and. derE .le. 2.8d0 ) then
                        if ( gammai(n) .gt. 0.9993908d0 ) then
                           icheckp2a=icheckp2a+1
                           if ( icheckp2a .eq. 450 ) then
                              icheckp2a=0
                              goto 998
                           end if
                           nclst = -1
                           return
                        end if

                     else if ( derE .gt. 1.7d0
     1                       .and. derE .lt. 2.3d0 ) then
                        if ( gammai(n) .gt. 0.9993908d0 ) then
                           icheckp3a=icheckp3a+1
                           if ( icheckp3a .eq. 3 ) then
                              icheckp3a=0
                              goto 998
                           end if
                           nclst = -1
                           return
                        end if
                     end if


                  end if
                  end if
 998              continue

                  end if

*-----------------------------------------------------------------------

                  nbart = nbart + ibary
                  nchat = nchat + ipchg
                  tbene = tbene + bene

                  sume = sume + epi(n)

                  epp = epi(n) / 1000.0d0

                  pouta = sqrt( epp**2 + 2.0d0 * epp * rms )

                  pxrv  = pouta * alphai(n)
                  pyrv  = pouta * betai(n)
                  pzrv  = pouta * gammai(n)

                  poutx = poutx + pxrv
                  pouty = pouty + pyrv
                  poutz = poutz + pzrv

*-----------------------------------------------------------------------
*        booking of outgoing particles
*-----------------------------------------------------------------------

                  iclust(nclst)    = ipid

                  jclust(0,nclst)  = 0
                  jclust(1,nclst)  = ipprt
                  jclust(2,nclst)  = ipneu
                  jclust(3,nclst)  = ippad
                  jclust(4,nclst)  = 0
                  jclust(5,nclst)  = ipchg
                  jclust(6,nclst)  = ibary
                  jclust(7,nclst)  = kf
                  jclust(8,nclst)  = 0

                  qclust(0,nclst)  = bimpacti
                  qclust(1,nclst)  = pxrv
                  qclust(2,nclst)  = pyrv
                  qclust(3,nclst)  = pzrv
                  qclust(4,nclst)  = epp + rms
                  qclust(5,nclst)  = rms
                  qclust(6,nclst)  = 0.0d0
                  qclust(7,nclst)  = epp * 1000.d0
                  qclust(8,nclst)  = 1.0d0
                  qclust(9,nclst)  = 0.0d0
                  qclust(10,nclst) = 0.0d0
                  qclust(11,nclst) = 0.0d0
                  qclust(12,nclst) = 0.0d0

            end do

         end if

*-----------------------------------------------------------------------
*           residual nucleus
*-----------------------------------------------------------------------

         masrs = masim - nbart
         mchrs = mchim - nchat

         if( masrs .gt. 0 .and. mchrs .ge. 0 ) then

C S.Hashimoto added IF statement for clusters composed purely of
C neutrons or protons. (2014.11.4)
          if ( ( masrs .gt. 1 .and. mchrs .eq. 0 ) .or.
     &           ( mchrs .gt. 1 .and. masrs .eq. mchrs ) ) then

             nresidue = max( masrs-mchrs, mchrs )

             piabs = dsqrt( eppin**2 + 2.0d0 * eppin * rmsin )
             presx = - poutx
             presy = - pouty
             presz = - poutz + piabs
             pabst = dsqrt( presx**2 + presy**2 + presz**2 )
             rsmas = rstms(1) * dble( mchrs )
     &            + rstms(2) * dble( masrs-mchrs )
             etota = dsqrt( pabst**2 + rsmas**2 )
             erres = ( etota - rsmas ) * 1000.0d0
             exres = eein + einad
     &            - sume
     &            - rstms(3)*1d3 * npipo - rstms(5)*1d3 * nping
     &            - rstms(4)*1d3 * npine
     &            - erres
     &            - bindeg(mchg,mmas-mchg)
     &            + tbene

           if ( exres + eein .lt. 0d0 ) then
              nclst = -1
              return
           end if

           absexres = dabs(exres)
           if ( absexres .gt. 0d0 ) then
              niteration = 10000
              epsiloni = 1d-5   !(10 keV)

            if ( exres .lt. 0d0 ) then
               facti = 1d0 - 1d-5
            else
               facti = 1d0 + 1d-5
            end if

            iteration = 0
            exresini = exres
            do while ( exresini*exres .gt. 0d0 .and.
     &           iteration .le. niteration )
               iteration = iteration + 1
               poutx = 0d0
               pouty = 0d0
               poutz = 0d0
               sume = 0d0

             do iclst=1,nclst
                qclust(1,iclst) = qclust(1,iclst) * facti
                qclust(2,iclst) = qclust(2,iclst) * facti
                qclust(3,iclst) = qclust(3,iclst) * facti
                qclust(4,iclst) =
     &               dsqrt( qclust(1,iclst)**2
     &               + qclust(2,iclst)**2 + qclust(3,iclst)**2
     &               + qclust(5,iclst)**2 )

                qclust(7,iclst) = 1000d0
     &               * ( qclust(4,iclst) - qclust(5,iclst) )
              if ( qclust(7,iclst) .lt. 0d0 ) then
                 nclst = -1
                 return
              end if

              poutx = poutx + qclust(1,iclst)
              pouty = pouty + qclust(2,iclst)
              poutz = poutz + qclust(3,iclst)
              sume = sume + qclust(7,iclst)

             end do

             piabs = dsqrt( eppin**2 + 2.0d0 * eppin * rmsin )
             presx = - poutx
             presy = - pouty
             presz = - poutz + piabs
             pabst = dsqrt( presx**2 + presy**2 + presz**2 )
             rsmas = rstms(1) * dble( mchrs )
     &            + rstms(2) * dble( masrs-mchrs )
             etota = dsqrt( pabst**2 + rsmas**2 )
             erres = ( etota - rsmas ) * 1000.0d0
             exres = eein + einad
     &            - sume
     &            - rstms(3)*1d3 * npipo - rstms(5)*1d3 * nping
     &            - rstms(4)*1d3 * npine
     &            - erres
     &            - bindeg(mchg,mmas-mchg)
     &            + tbene

            end do

            if( dabs(exres/1000d0) .le. epsiloni ) then
               exres = 0d0

            else
               nclst = -1
               return

            end if

           end if

           if( exres .lt. 0.0d0 ) then

              nclst = -1
              return

           end if

           exres = max( 0.0d0, exres )

*-----------------------------------------------------------------------
*           booking of the residual nucleons (only neutrons or protons)
*-----------------------------------------------------------------------

           if ( mchrs .eq. 0 ) then ! for only neutrons

            do iresidue = 1, nresidue

               nclst = nclst + 1

               iclust(nclst)    = 2

               jclust(0,nclst)  = 0
               jclust(1,nclst)  = 0
               jclust(2,nclst)  = 1
               jclust(3,nclst)  = 2
               jclust(4,nclst)  = 0
               jclust(5,nclst)  = 0
               jclust(6,nclst)  = 1
               jclust(7,nclst)  = 2112
               jclust(8,nclst)  = 0

               qclust(0,nclst)  = bimpacti
               qclust(1,nclst)  = presx / dble( nresidue )
               qclust(2,nclst)  = presy / dble( nresidue )
               qclust(3,nclst)  = presz / dble( nresidue )
               qclust(4,nclst)  = etota / dble( nresidue )
               qclust(5,nclst)  = rstms(2)
               qclust(6,nclst)  = 0.0d0
               fkinene = (qclust(4,nclst) - qclust(5,nclst)) *1000d0
               qclust(7,nclst)  = fkinene
               qclust(8,nclst)  = 1.0d0
               qclust(9,nclst)  = 0.0d0
               qclust(10,nclst) = 0.0d0
               qclust(11,nclst) = 0.0d0
               qclust(12,nclst) = 0.0d0

             if( fkinene .lt. 0.0d0 ) then
                nclst = -1
                return
             end if

            end do

*-----------------------------------------------------------------------

           else if  ( masrs .eq. mchrs ) then ! for only protons

            do iresidue = 1, nresidue

               nclst = nclst + 1

               iclust(nclst)    = 1

               jclust(0,nclst)  = 0
               jclust(1,nclst)  = 1
               jclust(2,nclst)  = 0
               jclust(3,nclst)  = 1
               jclust(4,nclst)  = 0
               jclust(5,nclst)  = 1
               jclust(6,nclst)  = 1
               jclust(7,nclst)  = 2212
               jclust(8,nclst)  = 0

               qclust(0,nclst)  = bimpacti
               qclust(1,nclst)  = presx / dble( nresidue )
               qclust(2,nclst)  = presy / dble( nresidue )
               qclust(3,nclst)  = presz / dble( nresidue )
               qclust(4,nclst)  = etota / dble( nresidue )
               qclust(5,nclst)  = rstms(1)
               qclust(6,nclst)  = 0.0d0
               fkinene = (qclust(4,nclst) - qclust(5,nclst)) *1000d0
               qclust(7,nclst)  = fkinene
               qclust(8,nclst)  = 1.0d0
               qclust(9,nclst)  = 0.0d0
               qclust(10,nclst) = 0.0d0
               qclust(11,nclst) = 0.0d0
               qclust(12,nclst) = 0.0d0

             if( fkinene .lt. 0.0d0 ) then
                nclst = -1
                return
             end if

            end do

*-----------------------------------------------------------------------

           else

              nclst = -1
              return

           end if

*-----------------------------------------------------------------------

          else ! except for only neutrons or only protons

             nclst = nclst + 1

             piabs = sqrt( eppin**2 + 2.0d0 * eppin * rmsin )

             presx = - poutx
             presy = - pouty
             presz = - poutz + piabs

             pabst = sqrt( presx**2 + presy**2 + presz**2 )
             rsmas = rstms(1) * dble( mchrs )
     &            + rstms(2) * dble( masrs-mchrs )
     &            - bindeg(mchrs,masrs-mchrs) / 1000.0d0

             etota = sqrt( pabst**2 + rsmas**2 )

             erres = ( etota - rsmas ) * 1000.0d0

             exres = eein + einad
     &            - sume
     &            - rstms(3)*1d3 * npipo - rstms(5)*1d3 * nping
     &            - rstms(4)*1d3 * npine
     &            - erres
     &            - bindeg(mchg,mmas-mchg)
     &            + bindeg(mchrs,masrs-mchrs)
     &            + tbene

           if( exres .lt. 0.0d0 ) then

              nclst = -1
              return

           end if

           exres = max( 0.0d0, exres )

*-----------------------------------------------------------------------
*           booking of the residual nucleus
*-----------------------------------------------------------------------

           iclust(nclst)    = 0

           jclust(0,nclst)  = jj
           jclust(1,nclst)  = mchrs
           jclust(2,nclst)  = masrs - mchrs
           jclust(3,nclst)  = 19
           jclust(4,nclst)  = 0
           jclust(5,nclst)  = mchrs
           jclust(6,nclst)  = masrs
           jclust(7,nclst)  = mchrs * 1000000 + masrs
           jclust(8,nclst)  = 0

           qclust(0,nclst)  = bimpacti
           qclust(1,nclst)  = presx
           qclust(2,nclst)  = presy
           qclust(3,nclst)  = presz
           qclust(4,nclst)  = etota
           qclust(5,nclst)  = rsmas
           qclust(6,nclst)  = exres
           qclust(7,nclst)  = ( etota - rsmas ) * 1000.d0
           qclust(8,nclst)  = 1.0d0
           qclust(9,nclst)  = 0.0d0
           qclust(10,nclst) = 0.0d0
           qclust(11,nclst) = 0.0d0
           qclust(12,nclst) = 0.0d0

          end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*        error in incl
*-----------------------------------------------------------------------

  340 continue

      write(6,1002) ityp,lk
 1002 format(/' *** error message from incl ***'
     &     /' invalid condition of itype or lk was found.'
     &     /' ityp =',i5,'  lk =',i5)
      call parastop( 844 )

*-----------------------------------------------------------------------

      return
      end
