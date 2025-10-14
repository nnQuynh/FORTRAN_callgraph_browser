c --------------------------------------------------------------------
c     2007/June version
c -------------------------------------------------------------------
c
c     Judgement of re-clustering after light complex-particle
c     coalescence
c
c     Taking into account the Coulomb barrier by a simple method
c     for all cases.
c
c     A method of picking up a nucleon is modified.
c
c     Input data has changed for a nuclear radius necessary for
c     judgement of leading nucleons
c
************************************************************************
*                                                                      *
*     L.P. Approx. (multi l.p. mode)                                   *
*     rclsmin < R_cm < rlpjdg : light cluster creation region          *
*                                                                      *
*     d + t + 3He + alpha production                                   *
*                                                                      *
*       p + n --> d                                                    *
*       d + n --> t,  d + p --> 3He                                    *
*       t + p --> alpha,  3He + n --> alpha                            *
*                                                                      *
************************************************************************
*                                                                      *
      subroutine lpinit
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

      parameter ( nnnn = 800 )

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson
!$OMP THREADPRIVATE(/vriab0/)

      common /grndc2/ r00, r01, saa, rada, radb
!$OMP THREADPRIVATE(/grndc2/)
      common /ldpart/ rlpjdg,ilpflg,idlp(nnnn)
!$OMP THREADPRIVATE(/ldpart/)
      common /lclsti/ eklcp(0:nnnn),ilclst(0:nnnn,0:nnnn)
!$OMP THREADPRIVATE(/lclsti/)
      common /clpara/ fdr,fdp,hhh,ddd,rcls00,rclsmin,iprnt
!$OMP THREADPRIVATE(/clpara/)

*-----------------------------------------------------------------------

      do k = 1, nnnn
         idlp(k) = 0
      enddo

      do i = 0, nnnn
         eklcp(i) = 0.d0
         do j = 0, nnnn
            ilclst(i,j) = 0
         enddo
      enddo

      ilpflg = 0


      if( rcls00 .lt. 0.d0 ) then
         radcls = -rcls00 * dble(massal)**(1.d0/3.d0)
      else
         radcls = rcls00
      endif

      rlpjdg = radcls + ddd

      if( rclsmin .lt. 0.d0 ) then

         rclsmin = radcls

      endif

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine lpjudge
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

      parameter ( nnnn = 800 )

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /const2/ dt, ntmax, iprun, iprun0
!$OMP THREADPRIVATE(/const2/)
      common /vriab0/ massal, massba, mmeson
!$OMP THREADPRIVATE(/vriab0/)
      common /vriab1/ b, llnow, ntnow
!$OMP THREADPRIVATE(/vriab1/)
      common /coodrp/ r(5,nnn),  p(6,nnn)
!$OMP THREADPRIVATE(/coodrp/)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)
      common /grndc2/ r00, r01, saa, rada, radb
!$OMP THREADPRIVATE(/grndc2/)
      common /const5/ pzpr, pxpr, rzpr, rxpr
!$OMP THREADPRIVATE(/const5/)
      common /const6/ pzta, pxta, rzta, rxta
!$OMP THREADPRIVATE(/const6/)
       common /coultr/ eccm, pzcc, rmax0, zeroz
!$OMP THREADPRIVATE(/coultr/)

      common /ldpart/ rlpjdg,ilpflg,idlp(nnnn)
!$OMP THREADPRIVATE(/ldpart/)
      common /lclsti/ eklcp(0:nnnn),ilclst(0:nnnn,0:nnnn)
!$OMP THREADPRIVATE(/lclsti/)
      common /qmdscmpar/ iesc(nnnn)
!$OMP THREADPRIVATE(/qmdscmpar/)

*-----------------------------------------------------------------------

      dimension idlp0(nnnn), itlp(nnnn)

*-----------------------------------------------------------------------

      do i = 1, nnnn
          idlp0(i) = 0
          itlp(i) = 0
      enddo

      nlp0 = 0

      do 100 j = 1, massal

         rpos = dsqrt( (r(1,j)-rxta)**2 + r(2,j)**2 + (r(3,j)-rzta)**2 )

         if( rpos .gt. rlpjdg .and.
     &       iesc(j) .eq. 1 .and. ihis(j) .gt. 0 ) then

            call epotlp(j,epot)

            ekin = p(4,j) - p(5,j)

            if( ekin+epot .gt. 0.d0 ) then
               nlp0 = nlp0 + 1
               idlp0(nlp0) = j
            endif

         endif

 100  continue

*-----------------------------------------------------------------------
*     case of elastc scatt.
*-----------------------------------------------------------------------

      if( ntnow .eq. ntmax ) then

         j = massal
         if( ilpflg .eq. 0 .and. ihis(j) .eq. 0 ) then

            rpos = dsqrt( (r(1,j)-rxta)**2 + r(2,j)**2
     &                   + (r(3,j)-rzta)**2 )

            call epotlp(j,epot)

            ekin = p(4,j) - p(5,j)

            if( ekin+epot .gt. 0.d0 ) then
               nlp0 = nlp0 + 1
               idlp0(nlp0) = j
            endif

         endif

      endif

*-----------------------------------------------------------------------
*     when a new l.p. appears
*-----------------------------------------------------------------------

      if( nlp0 .gt. ilpflg ) then

         nlp = 0

         do n1 = 1, nlp0

            isame = 0

            do n2 = 1, ilpflg
               if( idlp0(n1) .eq. idlp(n2) ) isame = 1
            enddo

            if( isame .eq. 0 ) then

               ilpflg = ilpflg + 1
               idlp(ilpflg) = idlp0(n1)

               nlp = nlp +1
               itlp(nlp) = idlp0(n1)

            endif

         enddo

         do n = 1, nlp
            call clsfrm(itlp(n))      ! Judgment of Cluster formation
         enddo

      endif

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine clsfrm(lp)
*                                                                      *
*                                                                      *
*        Last Revised:     2006 5 18  (by YW@Kyushu)                   *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to produce light ion clusters in terms of surface       *
*              coalescence model                                       *
*                                                                      *
*        Variables:  in common                                         *
*                                                                      *
************************************************************************
      use QMD_COOD2_MOD, only : rha, rhe, rhc, num

      implicit real*8(a-h,o-z)

      parameter ( nnnn = 800 )

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson
!$OMP THREADPRIVATE(/vriab0/)
      common /vriab1/ b, llnow, ntnow
!$OMP THREADPRIVATE(/vriab1/)

      common /coodrp/ r(5,nnn),  p(6,nnn)
!$OMP THREADPRIVATE(/coodrp/)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)

      common /clpara/ fdr,fdp,hhh,ddd,rcls00,rclsmin,iprnt
!$OMP THREADPRIVATE(/clpara/)
      common /ldpart/ rlpjdg,ilpflg,idlp(nnnn)
!$OMP THREADPRIVATE(/ldpart/)
      common /lclsti/ eklcp(0:nnnn),ilclst(0:nnnn,0:nnnn)
!$OMP THREADPRIVATE(/lclsti/)

      common /const5/ pzpr, pxpr, rzpr, rxpr
!$OMP THREADPRIVATE(/const5/)
      common /const6/ pzta, pxta, rzta, rxta
!$OMP THREADPRIVATE(/const6/)
      common /rpdist/ rrr2(nnnn),ir2srt(nnnn)
!$OMP THREADPRIVATE(/rpdist/)
      common /coultr/ eccm, pzcc, rmax0, zeroz
!$OMP THREADPRIVATE(/coultr/)

      common /qmdscmpar/ iesc(nnnn)
!$OMP THREADPRIVATE(/qmdscmpar/)

*-----------------------------------------------------------------------

      dimension mclst0(4), ekinlcp(4)

*-----------------------------------------------------------------------
*     identification of the cluster (deuteron)
*-----------------------------------------------------------------------

      do i = 1, nnnn
         num(i) = 0
      enddo

      do i = 1, 4
         ekinlcp(i) = 0.d0
      enddo

      i = lp         ! ID no. of the leading particle
      iclsflg = 0

*-----------------------------------------------------------------------
*     formation of possible deuterons
*-----------------------------------------------------------------------

      call r2sort(r(1,i),r(2,i),r(3,i))

      do 10 j = 1, massal

         jp = ir2srt(j)

         if( jp .eq. i .or. iesc(jp) .eq. 0 ) goto 10

         rdist2 = rrr2(jp)

c   Use of the definition of Letourneau (NPA 712, 2002)

         pdist2 = ( (p(1,i)-p(1,jp))**2 + (p(2,i)-p(2,jp))**2
     &             + (p(3,i)-p(3,jp))**2 ) / 4.d0

         rj = dsqrt( (r(1,jp)-rxta)**2 + r(2,jp)**2
     &              + (r(3,jp)-rzta)**2 )

         rdcmx = ( r(1,i) + r(1,jp) ) / 2.d0
         rdcmy = ( r(2,i) + r(2,jp) ) / 2.d0
         rdcmz = ( r(3,i) + r(3,jp) ) / 2.d0

         rccm = dsqrt((rdcmx-rxta)**2+rdcmy**2+(rdcmz-rzta)**2)

*-----------------------------------------------------------------------
*        deuteron coalescence
*-----------------------------------------------------------------------

         if( dsqrt(rdist2)*dsqrt(pdist2) .le. hhh .and.
     *       ichg(i)+ichg(jp) .eq. 1 .and.
     *       rccm .gt. rclsmin .and. rj .lt. rlpjdg ) then

            iclsflg = 1

            num(1) = i
            num(2) = jp
            mclst = 2
            mclst0(iclsflg) = mclst

            call ekincal(ekin,epot,iclsflg,num,mclst)
            ekinlcp(iclsflg) = ekin

            pdcmx = p(1,i) + p(1,jp)
            pdcmy = p(2,i) + p(2,jp)
            pdcmz = p(3,i) + p(3,jp)

            if(iprnt.eq.1) then
               write(51,'(2i7,4i3,4f10.4)')
     &          llnow, ntnow, i, jp, ichg(i),
     &          ichg(jp), rj, rccm, epot, ekin
            endif

            goto 20

         endif

 10   continue

*-----------------------------------------------------------------------

      goto 999   ! return

*-----------------------------------------------------------------------

 20   continue

      call r2sort(rdcmx,rdcmy,rdcmz)

      do 30 j = 1, massal

         jp = ir2srt(j)

         if( jp .eq. num(1) .or. jp .eq. num(2) .or.
     &       iesc(jp) .eq. 0 ) goto 30

         rdist2 = rrr2(jp)

         pdist2 = ( p(1,jp)*2.d0/3.d0 - pdcmx/3.d0 )**2
     &           + ( p(2,jp)*2.d0/3.d0 - pdcmy/3.d0 )**2
     &           + ( p(3,jp)*2.d0/3.d0 - pdcmz/3.d0 )**2

         rj = dsqrt( (r(1,jp)-rxta)**2 + r(2,jp)**2
     &                 +(r(3,jp)-rzta)**2 )

         i1 = num(1)
         i2 = num(2)
         i3 = jp

         rdcmx = ( r(1,i1) + r(1,i2) + r(1,i3) ) / 3.d0
         rdcmy = ( r(2,i1) + r(2,i2) + r(2,i3) ) / 3.d0
         rdcmz = ( r(3,i1) + r(3,i2) + r(3,i3) ) / 3.d0

         rccm = dsqrt( (rdcmx-rxta)**2 + rdcmy**2 + (rdcmz-rzta)**2 )

         if( dsqrt(rdist2) * dsqrt(pdist2) .le. hhh .and.
     &       rccm .gt. rclsmin .and. rj .lt. rlpjdg ) then

            if( ichg(jp) .eq. 0 ) iclsflg = 2      ! triton
            if( ichg(jp) .eq. 1 ) iclsflg = 3      ! He-3

            num(3) = jp
            mclst = 3
            mclst0(iclsflg) = mclst

            call ekincal(ekin,epot,iclsflg,num,mclst)
            ekinlcp(iclsflg) = ekin

            i1 = num(1)
            i2 = num(2)
            i3 = num(3)

            pdcmx = p(1,i1) + p(1,i2) + p(1,i3)
            pdcmy = p(2,i1) + p(2,i2) + p(2,i3)
            pdcmz = p(3,i1) + p(3,i2) + p(3,i3)

            if(iprnt.eq.1) then
               write(52,'(2i7,4i3,4f10.4)')
     &          llnow, ntnow, num(1), num(2), num(3),
     &          ichg(jp), rj, rccm, epot, ekin
            endif

            goto 40

         endif

 30   continue

*-----------------------------------------------------------------------

 40   continue

*-----------------------------------------------------------------------
*     alpha coalescence
*-----------------------------------------------------------------------

      if( iclsflg .eq. 2 .or. iclsflg .eq. 3 ) then

         call r2sort(rdcmx,rdcmy,rdcmz)

         do 50 j = 1, massal

            jp= ir2srt(j)
            if( jp .eq. num(1) .or. jp .eq. num(2) .or.
     &          jp .eq. num(3) .or. iesc(jp) .eq. 0 ) goto 50

            rdist2 = rrr2(jp)
            pdist2 = ( p(1,jp)*3.d0/4.d0 - pdcmx/4.d0 )**2
     &              + ( p(2,jp)*3.d0/4.d0 - pdcmy/4.d0 )**2
     &              + ( p(3,jp)*3.d0/4.d0 - pdcmz/4.d0 )**2

            rj = dsqrt( (r(1,jp)-rxta)**2 + r(2,jp)**2
     &                 + (r(3,jp)-rzta)**2 )

            i1 = num(1)
            i2 = num(2)
            i3 = num(3)
            i4 = jp

            rdcmx = ( r(1,i1) + r(1,i2) + r(1,i3) + r(1,i4) ) / 4.d0
            rdcmy = ( r(2,i1) + r(2,i2) + r(2,i3) + r(2,i4) ) / 4.d0
            rdcmz = ( r(3,i1) + r(3,i2) + r(3,i3) + r(3,i4) ) / 4.d0

            rccm = dsqrt( (rdcmx-rxta)**2 + rdcmy**2 + (rdcmz-rzta)**2 )

            if( dsqrt(rdist2)*dsqrt(pdist2) .le. hhh .and.
     &          rccm .gt. rclsmin .and. rj.lt.rlpjdg .and.
     &          ( iclsflg .eq. 2 .and. ichg(jp) .eq. 1 .or.
     &            iclsflg .eq. 3 .and .ichg(jp) .eq. 0 ) ) then

               iclsflg = 4    ! alpha
               num(4) = jp
               mclst = 4
               mclst0(iclsflg) = mclst

               call ekincal(ekin,epot,iclsflg,num,mclst)
               ekinlcp(iclsflg) = ekin

               if(iprnt.eq.1) then
                  write(53,'(2i7,5i3,4f10.4)')
     &             llnow, ntnow, num(1), num(2), num(3),num(4),
     &             ichg(jp), rj, rccm, epot, ekin
               endif

               goto 60

            endif

 50      continue

      endif

*-----------------------------------------------------------------------

 60   continue

      do icls = iclsflg, 1, -1

         if( ekinlcp(icls) .gt. 0.d0 ) then

            ilclst(0,0) = ilclst(0,0) + 1      ! total no. of light clusters
            i = ilclst(0,0)
            ilclst(i,0) = mclst0(icls)         ! number of nucleons in i-th cls
            do j = 1, mclst0(icls)
               ilclst(i,j) = num(j)
               iesc(num(j)) = 0                ! exclusion from composite system
            enddo

            eklcp(i) = ekinlcp(icls)  ! kinetic energy

            goto 100

         endif

      enddo

*-----------------------------------------------------------------------

 100  continue

*-----------------------------------------------------------------------

 999  return
      end


************************************************************************
*                                                                      *
      subroutine r2sort(rx0,ry0,rz0)
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

      parameter ( nnnn = 800 )

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson
!$OMP THREADPRIVATE(/vriab0/)
      common /coodrp/ r(5,nnn),  p(6,nnn)
!$OMP THREADPRIVATE(/coodrp/)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)

      common /rpdist/ rrr2(nnnn),ir2srt(nnnn)
!$OMP THREADPRIVATE(/rpdist/)

*-----------------------------------------------------------------------

      do i=1,massal
         rrr2(i) = (r(1,i)-rx0)**2 + (r(2,i)-ry0)**2 + (r(3,i)-rz0)**2
         ir2srt(i) = i
      enddo

*-----------------------------------------------------------------------
*     sorting
*-----------------------------------------------------------------------

 10   continue

      nexch = 0

      do i = 1, massal-1

         if( rrr2(ir2srt(i)) .gt. rrr2(ir2srt(i+1)) ) then

            nexch = nexch  + 1
            ibuf = ir2srt(i)
            ir2srt(i) = ir2srt(i+1)
            ir2srt(i+1) = ibuf

         endif

      enddo

      if( nexch .ne. 0 ) goto 10

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ekincal(ekinlcp,epot,iclsflg,inum,mclst)

*                                                                      *
*        Last Revised:     2006 5 18                                   *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate kinetic energy of emitted light cluster    *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              ekinlcp     : kinetic energy of lcp                     *
*              epot        : total potential energy                    *
*                                                                      *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg
      use QMD_COOD2_MOD, only : rha,rhe,rhc

      implicit real*8(a-h,o-z)

      parameter ( nnnn = 800 )

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson
!$OMP THREADPRIVATE(/vriab0/)
      common /vriab1/ b, llnow, ntnow
!$OMP THREADPRIVATE(/vriab1/)

      common /coodrp/ r(5,nnn),  p(6,nnn)
!$OMP THREADPRIVATE(/coodrp/)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)

      common /poten1/ gamm, c0, c3, cs, cl, wl
!$OMP THREADPRIVATE(/poten1/)

      common /qmdscmpar/ iesc(nnnn)
!$OMP THREADPRIVATE(/qmdscmpar/)

*-----------------------------------------------------------------------

      dimension inum(nnn), iprot(4), ineut(4)

*-----------------------------------------------------------------------

      data (iprot(i),i=1,4) /1,1,2,2/           ! d, t, 3He, alpha
      data (ineut(i),i=1,4) /1,2,1,2/

*-----------------------------------------------------------------------
*     Calculation of total potential energy for the cluster
*-----------------------------------------------------------------------

      epot = 0.d0

      do 100 n = 1, mclst

         i = inum(n)

         epot1 = 0.d0
         epot3 = 0.d0
         epots = 0.d0
         epotc = 0.d0

         do 110 j = 1, massal

            if( iesc(j) .eq. 0 ) goto 110

            epot1 = epot1 + rha(j,i)
            epotc = epotc + rhe(j,i)
            epots = epots + rha(j,i) * inuc(j) * inuc(i)
     &                     * ( 1.d0 - 2.d0 * dble(abs(ichg(j)-ichg(i))))

  110    continue

         epot3 = epot1 ** gamm

*-----------------------------------------------------------------------

         epot0 = c0 * epot1
     &          + c3 * epot3
     &          + cs * epots
     &          + cl * epotc

*-----------------------------------------------------------------------

         epot = epot + epot0

 100  continue

*-----------------------------------------------------------------------
*     Calculation of total kinetic energy
*-----------------------------------------------------------------------

      ekin = 0.d0

      do j = 1, mclst

         jp = inum(j)

         ekin = ekin + (p(4,jp) - p(5,jp))
      enddo

*-----------------------------------------------------------------------
*        ekinsum: total kinetic energy
*        epotcls: total potential + binding energy of the cluster
*        ekinlcp: emission energy of the cluster. It should be positive.
*-----------------------------------------------------------------------

      binlcp = bindeg(iprot(iclsflg),ineut(iclsflg)) / 1000.d0

      epotlcp = epot + binlcp
      ekinlcp = ekin + epotlcp

      masspar = 0
      nprot  = 0

      do 200 i = 1, massal

         masspar = masspar + iesc(i)
         nprot = nprot + ichg(i) * iesc(i)

 200  continue

      zlcp = dble( iprot(iclsflg) )
      zres = dble( nprot - iprot(iclsflg) )
      alcp = dble( iprot(iclsflg) + ineut(iclsflg) )
      ares = dble( masspar ) - alcp

      rcoul = 1.2d0 * (alcp**0.33333d0 + ares**0.33333d0)
      vcoul = 1.44d0 * zlcp * zres / rcoul /1000.d0         ! in GeV

      if( ekinlcp .ge. vcoul ) then

         fcoul = 1.d0 - vcoul / ekinlcp

         if( rn() .gt. fcoul ) ekinlcp = 0.d0

      else

         ekinlcp = 0.d0

      endif

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine cldist3
*                                                                      *
*                                                                      *
*        Last Revised:     2005 9 4  (by YW@Kyushu)                    *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine nuclear cluster.                           *
*                                                                      *
*        Variables:  in common                                         *
*                                                                      *
*              common /clusti/ itc(0:nnn,0:nnn)                        *
*              common /clustf/ nclst, iclust(nnn),                     *
*              common /clustg/ jclust(0:8,nnn),  qclust(0:12,nnn)      *
*                                                                      *
*              itc(0,0)    : total number of clusters plus particles   *
*              itc(i,0)    : number  of baryon in i-th cluster         *
*              itc(i,j)    : j-th id of baryon in i-th cluster         *
*                                                                      *
*              nclst       : total number of clusters plus particles   *
*                                                                      *
*              iclust(i)   : kind of i-th cluster or particle          *
*                                                                      *
*                     = 0  : Nucleus                                   *
*                     = 1  : proton                                    *
*                     = 2  : neutron                                   *
*                     = 3  : delta                                     *
*                     = 4  : N star                                    *
*                     = 5  : pions                                     *
*                     = 6  : gamma                                     *
*                     = 4, Gamma                                       *
*                     = 5, kaon                                        *
*                     = 6, muon                                        *
*                     = 7, others                                      *
*                                                                      *
*              jclust(k,i) : number of particles in i-th cluster       *
*                                                                      *
*                    (0,i) : jj, angular momentum                      *
*                    (1,i) : proton                                    *
*                    (2,i) : neutron                                   *
*                    (3,i) : delta                                     *
*                    (4,i) : N star                                    *
*                    (5,i) : charge                                    *
*                    (6,i) : collision history                         *
*                    (7,i) : sdm history                               *
*                   k = 0, angular momentum                            *
*                     = 1, proton number                               *
*                     = 2, neutron number                              *
*                     = 3, ip, see below                               *
*                     = 4, status of the particle 0: real, <0 : dead   *
*                     = 5, charge                                      *
*                     = 6, baryon number                               *
*                     = 7, kf code                                     *
*                                                                      *
*              qclust(k,i) : momentum of i-th cluster                  *
*                                                                      *
*                    (0,i) : b                                         *
*                    (1,i) : px                                        *
*                    (2,i) : py                                        *
*                    (3,i) : pz                                        *
*                    (4,i) : total energy, E = sqrt(m**2+p**2)         *
*                    (5,i) : rest mass                                 *
*                    (6,i) : excitation energy (MeV)                   *
*                                                                      *
*                   i = 0, impact parameter                            *
*                     = 1, px (GeV/c)                                  *
*                     = 2, py (GeV/c)                                  *
*                     = 3, pz (GeV/c)                                  *
*                     = 4, etot = sqrt( p**2 + rm**2 ) (GeV)           *
*                     = 5, rest mass (GeV)                             *
*                     = 6, excitation energy (MeV)                     *
*                     = 7, kinetic energy (MeV)                        *
*                     = 8, weight change                               *
*                     = 9, delay time                                  *
*                     = 10, x-displace                                 *
*                     = 11, y-displace                                 *
*                     = 12, z-displace                                 *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg
      use QMD_COOD2_MOD, only : rha, rhe, rhc, rr2, rbij, pp2,
     &                          mascl, num, isort, isorti, rhoa,
     &                          it, rs, ps

      implicit real*8(a-h,o-z)

      parameter ( nnnn = 800 )

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson
!$OMP THREADPRIVATE(/vriab0/)
      common /vriab1/ b, llnow, ntnow
!$OMP THREADPRIVATE(/vriab1/)

      common /coodrp/ r(5,nnn),  p(6,nnn)
!$OMP THREADPRIVATE(/coodrp/)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)
      common /poten1/ gamm, c0, c3, cs, cl, wl
!$OMP THREADPRIVATE(/poten1/)

      common /clusti/ itc(0:nnn,0:nnn)
!$OMP THREADPRIVATE(/clusti/)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn),  qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

      common /clpara/ fdr,fdp,hhh,ddd,rcls00,rclsmin,iprnt
!$OMP THREADPRIVATE(/clpara/)

      common /ldpart/ rlpjdg,ilpflg,idlp(nnnn)
!$OMP THREADPRIVATE(/ldpart/)
      common /lclsti/ eklcp(0:nnnn),ilclst(0:nnnn,0:nnnn)
!$OMP THREADPRIVATE(/lclsti/)
      common /qmdscmpar/ iesc(nnnn)
!$OMP THREADPRIVATE(/qmdscmpar/)

      common /const5/ pzpr, pxpr, rzpr, rxpr
!$OMP THREADPRIVATE(/const5/)
      common /const6/ pzta, pxta, rzta, rxta
!$OMP THREADPRIVATE(/const6/)

      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)

      common/jamevnt1/v(5,nnn),k(11,nnn)
!$OMP THREADPRIVATE(/jamevnt1/)

*-----------------------------------------------------------------------

      dimension       eklcp0(0:nnnn)

      dimension       num0(nnnn),num1(nnnn)
      dimension       mascl0(nnnn)

*-----------------------------------------------------------------------

      kcls = ilclst(0,0)      ! No. of produced light ion clusters
      kd = 0

      do i = 1, kcls          ! light clusters for kcls > 0

         mclst= ilclst(i,0)
         mascl(i)=mclst

         do j = 1, mclst
            kd = kd + 1
            num(kd)= ilclst(i,j)
         enddo

      enddo

*-----------------------------------------------------------------------
*     Re-idenification of clusters for the system except the light
*     clusters (d, t, 3He and alpha) --> whether leading nucleons
*     are re-trapped or not.
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*     Cluster distance
*-----------------------------------------------------------------------

      cpf2 = ( 1.5d0 * pi**2
     &             * ( 4.d0 * pi * wl )**(-1.5d0) )**(2.d0/3.d0)
     &      * hbc ** 2

      rcc2 = rclds ** 2

*-----------------------------------------------------------------------
*     Calculate overlap of the wave packets
*-----------------------------------------------------------------------

      do i = 1, massal

         rhoa(i) = 0.d0

         if( inds(i) .ne. 0 .and. iesc(i) .ne. 0 ) then

            do j = 1, massal

               rhoa(i) = rhoa(i) + rha(i,j)

            enddo

            rhoa(i) = ( rhoa(i) + 1.d0 ) ** (1.d0/3.d0)

         endif

      enddo

*-----------------------------------------------------------------------

      ij= 0

      do i = 1, massal

         if( inds(i) .ne. 0 .and. iesc(i) .ne. 0 ) then

            ij = ij + 1
            num1(ij) = i

         endif

      enddo

      nexlcp = ij      ! no. of nucleons which are not included in LCPs

      do i = 1, nexlcp

         num0(i) = i
         mascl0(i) = 1

      enddo

      nclst = 1
      ichek = 1

      do i = 1, nexlcp - 1

         j1 = ichek + 1
         id1 = inds( num1(num0(i)) )

         do j = j1, nexlcp

            id2 = inds( num1(num0(j)) )

            rdist2 = rr2( num1(num0(i)), num1(num0(j)) )
            pdist2 = pp2( num1(num0(i)), num1(num0(j)) )
            pcc2 = cpf2
     &            * ( rhoa( num1(num0(i)) ) + rhoa( num1(num0(j)) ) )**2

            if( rdist2 .lt. rcc2 .and. pdist2 .lt. pcc2 .and.
     &          id1 .eq. 1 .and. id2 .eq. 1) then

               ibuf = num0(ichek+1)
               num0(ichek+1) = num0(j)
               num0(j) = ibuf
               ichek = ichek + 1
               mascl0(nclst) = mascl0(nclst) + 1

            endif

         enddo

         if( ichek .eq. i ) then

            nclst = nclst + 1
            ichek = ichek + 1

         endif

      enddo

      do i = 1, nexlcp
         num(i+kd) = num1(num0(i))
      enddo

      do j = 1, nclst
         mascl(kcls+j) = mascl0(j)
      enddo

      nclst = nclst + kcls

*-----------------------------------------------------------------------
*     sort for summary
*-----------------------------------------------------------------------

      do i = 1, nclst

         isort(i) = i

      enddo

  510 continue

      nexch = 0

      do i = 1, nclst - 1

         if( mascl(isort(i)) .lt. mascl(isort(i+1)) ) then

            nexch = nexch + 1
            ibuf = isort(i)
            isort(i) = isort(i+1)
            isort(i+1) = ibuf

         endif

      enddo

      if( nexch .ne. 0 ) goto 510

      do i = 1, nclst

         isorti(isort(i)) = i

      enddo

      itc(0,0) = nclst

      inum = 0

      do i = 1, nclst

         itc(isorti(i),0) = mascl(i)
         eklcp0(isorti(i)) = eklcp(i)      ! kinetic energy for lcp

         do j = 1, mascl(i)

            inum = inum +1
            itc(isorti(i),j) = num(inum)

         enddo

      enddo


*-----------------------------------------------------------------------
*     add mesons
*-----------------------------------------------------------------------

      if( massal .gt. massba ) then

         do i = massba + 1, massal

            nclst = nclst + 1

            itc(nclst,0) = 1
            itc(nclst,1) = i

         enddo

         itc(0,0) = nclst

      endif

*-----------------------------------------------------------------------
*     Loop over clusters
*     Summarry of the cluster
*-----------------------------------------------------------------------

      mnucl = 0
      mprot = 0
      mneut = 0
      mpipo = 0
      mpine = 0
      mping = 0
      mmupo = 0
      mmune = 0
      mkapo = 0
      mkane = 0
      mkang = 0
      mothe = 0
      mgamm = 0

      ii = 0

*-----------------------------------------------------------------------

      do 1000 i = 1, nclst

*-----------------------------------------------------------------------

         mclst = itc(i,0)

         nchpa = 0

         nprot = 0
         nnuet = 0
         ndelt = 0
         nstar = 0

         jj = 0
         jc = 0
         js = 0

*-----------------------------------------------------------------------
*        This is cluster
*-----------------------------------------------------------------------

         if( mclst .gt. 1 ) then

            it(0) = mclst

            ipcst = 0

            mnucl = mnucl + 1

*-----------------------------------------------------------------------
*           Momentum of cluster.
*-----------------------------------------------------------------------

            pclx = 0.0
            pcly = 0.0
            pclz = 0.0

*-----------------------------------------------------------------------
*           Loop over mass of one cluster,
*           determine property of cluster.
*-----------------------------------------------------------------------

            do j = 1, mclst

               it(j) = itc(i,j)
               jp = it(j)

               pclx = pclx + p(1,jp)
               pcly = pcly + p(2,jp)
               pclz = pclz + p(3,jp)

               if( inds(jp) .eq. 1 ) then

                  if( ichg(jp) .eq. 1 ) then

                     nprot = nprot + 1

                  elseif( ichg(jp) .eq. 0 ) then

                     nnuet = nnuet + 1

                  endif

               elseif( inds(jp) .eq. 2 ) then

                     ndelt = ndelt + 1

               elseif( inds(jp) .eq. 3 ) then

                     nstar = nstar + 1

               endif

               nchpa = nchpa + ichg(jp)

            enddo

*-----------------------------------------------------------------------
*           Calculate binding energy of cluster.
*           Exciation energy and mass of the cluster.
*-----------------------------------------------------------------------

            if( nprot .gt. 0 .and. nnuet .gt. 0 ) then

               call etotal(1,ekin,epot,ebin,emas,epin,jj)

               texc = ebin * float( nprot + nnuet )
     &               + bindeg( nprot, nnuet ) / 1000.d0

               if( texc .lt. 0.d0 ) then
                  nclst = -1
                  return
               endif

               texc = dmax1( 0.d0, texc )

               tmas = ( emas + ebin ) * it(0)

               if( mclst.le.4 ) then          ! This is a light ion cluster
                  texc = 0.d0                 ! g.s. is assumed
                  tmas = 0.d0
                  do j = 1, mclst
                     jp = it(j)
                     tmas = tmas + p(5,jp)
                  enddo
               endif

               pabs = dsqrt( pclx**2 + pcly**2 + pclz**2 )
               etot = dsqrt( tmas**2
     &                      + pclx**2 + pcly**2 + pclz**2 )

               kf = nprot * 1000000 + nprot + nnuet
               ibary = nprot + nnuet

*-----------------------------------------------------------------------

               ii = ii + 1

               iclust(ii) = ipcst

               jclust(0,ii) = jj
               jclust(1,ii) = nprot
               jclust(2,ii) = nnuet
               jclust(3,ii) = 19
               jclust(4,ii) = 0
               jclust(5,ii) = nchpa
               jclust(6,ii) = ibary
               jclust(7,ii) = kf
               jclust(8,ii) = 0

               qclust(0,ii)  = b
               qclust(1,ii)  = pclx
               qclust(2,ii)  = pcly
               qclust(3,ii)  = pclz
               qclust(4,ii)  = etot
               qclust(5,ii)  = tmas
               qclust(6,ii)  = texc * 1000.d0
               qclust(7,ii)  = ( etot - tmas ) * 1000.d0
               qclust(8,ii)  = 1.d0
               qclust(9,ii)  = 0.d0
               qclust(10,ii) = 0.d0
               qclust(11,ii) = 0.d0
               qclust(12,ii) = 0.d0

               if( mclst .ge. 2 .and. mclst .le. 4 ) then      ! light ion cluster

                  ekin = eklcp0(i)      ! kinetic energy
                  pcl2p = ekin*(ekin+2.d0*tmas)

                  pcl2 = pclx**2 + pcly**2 + pclz**2
                  hlp = dsqrt(pcl2p/pcl2)

                  qclust(1,ii) = pclx * hlp
                  qclust(2,ii) = pcly * hlp
                  qclust(3,ii) = pclz * hlp
                  qclust(4,ii) = ekin + tmas

               endif

*-----------------------------------------------------------------------
*           strange nucleus
*-----------------------------------------------------------------------

            else

               do j = 1, mclst

                  ip = itc(i,j)

                  if( ichg(ip) .eq. 1 ) then

                     nprot = 1
                     nnuet = 0
                     ipcst = 1
                     ippad = 1
                     kf = 2212

                     mprot = mprot + 1

                  elseif( ichg(ip) .eq. 0 ) then

                     nprot = 0
                     nnuet = 1
                     ipcst = 2
                     ippad = 2
                     kf = 2112

                     mneut = mneut + 1

                  endif

                  pclx = p(1,ip)
                  pcly = p(2,ip)
                  pclz = p(3,ip)
                  pabs = dsqrt( pclx**2 + pcly**2 + pclz**2 )

                  etot = p(4,ip)
                  tmas = p(5,ip)

                  texc = 0.d0

                  ii = ii + 1

                  iclust(ii) = ipcst

                  jclust(0,ii) = jj
                  jclust(1,ii) = nprot
                  jclust(2,ii) = nnuet
                  jclust(3,ii) = ippad
                  jclust(4,ii) = 0
                  jclust(5,ii) = nprot
                  jclust(6,ii) = 1
                  jclust(7,ii) = kf
                  jclust(8,ii) = 0

                  qclust(0,ii)  = b
                  qclust(1,ii)  = pclx
                  qclust(2,ii)  = pcly
                  qclust(3,ii)  = pclz
                  qclust(4,ii)  = etot
                  qclust(5,ii)  = tmas
                  qclust(6,ii)  = texc * 1000.d0
                  qclust(7,ii)  = ( etot - tmas ) * 1000.d0
                  qclust(8,ii)  = 1.d0
                  qclust(9,ii)  = 0.d0
                  qclust(10,ii) = 0.d0
                  qclust(11,ii) = 0.d0
                  qclust(12,ii) = 0.d0

               enddo

            endif

*-----------------------------------------------------------------------
*        This is a baryon or meson
*-----------------------------------------------------------------------

         else

*-----------------------------------------------------------------------

            ip  = itc(i,1)

            if( inds(ip) .eq. 1 ) then

               if( ichg(ip) .eq. 1 ) then

                  nprot = 1
                  ipcst = 1

                  ippad = 1
                  ibary = 1
                  kf = 2212

                  mprot = mprot + 1

               elseif( ichg(ip) .eq. 0 ) then

                  nnuet = 1
                  ipcst = 2

                  nprot = 0
                  ippad = 2
                  ibary = 1
                  kf = 2112

                  mneut = mneut + 1

               endif

*-----------------------------------------------------------------------

            elseif( inds(ip) .eq. 2 ) then

               ndelt = 1

               ippad = 11
               ibary = 1
               nprot = 0
               nnuet = 0

               if( ichg(ip) .eq. 2 ) then

                  kf = 2224

               elseif( ichg(ip) .eq.  1 ) then

                  kf = 2214

               elseif( ichg(ip) .eq.  0 ) then

                  kf = 2114

               elseif( ichg(ip) .eq. -1 ) then

                  kf = 1114

               endif

               ipcst = 7
               mothe = mothe + 1

*-----------------------------------------------------------------------

            elseif( inds(ip) .eq. 3 ) then

               nstar = 1

               nstar = 1
               ippad = 11
               ibary = 1
               nprot = 0
               nnuet = 0

               if( ichg(ip) .eq. 1 ) then

                  kf = 12212

               elseif( ichg(ip) .eq. 0 ) then

                  kf = 12112

               endif

               ipcst = 7
               mothe = mothe + 1

*-----------------------------------------------------------------------

            elseif( inds(ip) .eq. 4 ) then

               if( ichg(ip) .eq. 1 ) then

                  kf = 211
                  ippad = 3
                  mpipo = mpipo + 1

               elseif( ichg(ip) .eq.  0 ) then

                  kf = 111
                  ippad = 4
                  mpine = mpine + 1

               elseif( ichg(ip) .eq. -1 ) then

                  kf = -211
                  ippad = 5
                  mping = mping + 1

               endif

               ipcst = 3
               ibary = 0
               nprot = 0
               nnuet = 0

*-----------------------------------------------------------------------

            elseif( inds(ip) .eq. 5 ) then

               kf = k(2,ip)
               ibary = k(9,ip) / 3
               nprot = 0
               nnuet = 0

               if( kf .eq. 22 ) then

                  ipcst = 4
                  ippad = 14

                  mgamm = mgamm + 1

               elseif( kf .eq. -321 .or.
     &                  kf .eq.  311 .or.
     &                  kf .eq.  321 ) then

                  ipcst = 5

                  if( kf .eq. 321 ) then

                     mkapo = mkapo + 1
                     ippad = 8

                  elseif( kf .eq. 311 ) then

                     mkane = mkane + 1
                     ippad = 9

                  elseif( kf .eq. -321 ) then

                     mkang = mkang + 1
                     ippad = 10

                  end if

               elseif( kf .eq. -13 .or.
     &                 kf .eq.  13 ) then

                  ipcst = 6

                  if( kf .eq. 13 ) then

                     mmupo = mmupo + 1
                     ippad = 6

                  elseif( kf .eq. -13 ) then

                     mmune = mmune + 1
                     ippad = 7

                  endif

               else

                  mothe = mothe + 1

                  ipcst = 7
                  ippad = 11

               end if

            end if

            nchpa = ichg(ip)

*-----------------------------------------------------------------------

            pclx = p(1,ip)
            pcly = p(2,ip)
            pclz = p(3,ip)

            pabs = dsqrt( pclx**2 + pcly**2 + pclz**2 )

            etot = p(4,ip)
            tmas = p(5,ip)

            texc = 0.d0

*-----------------------------------------------------------------------
*                 Multi-step history for one baryon case
*-----------------------------------------------------------------------

            jc = ihis(ip)

            if( jc .gt.  100 ) jc =  100
            if( jc .lt. -100 ) jc = -100

*-----------------------------------------------------------------------

            ii = ii + 1

            iclust(ii) = ipcst

            jclust(0,ii) = jj
            jclust(1,ii) = nprot
            jclust(2,ii) = nnuet
            jclust(3,ii) = ippad
            jclust(4,ii) = 0
            jclust(5,ii) = nchpa
            jclust(6,ii) = ibary
            jclust(7,ii) = kf
            jclust(8,ii) = 0

            qclust(0,ii)  = b
            qclust(1,ii)  = pclx
            qclust(2,ii)  = pcly
            qclust(3,ii)  = pclz
            qclust(4,ii)  = etot
            qclust(5,ii)  = tmas
            qclust(6,ii)  = texc * 1000.d0
            qclust(7,ii)  = ( etot - tmas ) * 1000.d0
            qclust(8,ii)  = 1.d0
            qclust(9,ii)  = 0.d0
            qclust(10,ii) = 0.d0
            qclust(11,ii) = 0.d0
            qclust(12,ii) = 0.d0

*-----------------------------------------------------------------------

         end if

 1000 continue

*-----------------------------------------------------------------------

      nclst = ii

*-----------------------------------------------------------------------

      numpat(0)  = mnucl
      numpat(1)  = mprot
      numpat(2)  = mneut
      numpat(3)  = mpipo
      numpat(4)  = mpine
      numpat(5)  = mping
      numpat(6)  = mmupo
      numpat(7)  = mmune
      numpat(8)  = mkapo
      numpat(9)  = mkane
      numpat(10) = mkang
      numpat(11) = mothe
      numpat(14) = mgamm

      rumpat(0)  = mnucl
      rumpat(1)  = mprot
      rumpat(2)  = mneut
      rumpat(3)  = mpipo
      rumpat(4)  = mpine
      rumpat(5)  = mping
      rumpat(6)  = mmupo
      rumpat(7)  = mmune
      rumpat(8)  = mkapo
      rumpat(9)  = mkane
      rumpat(10) = mkang
      rumpat(11) = mothe
      rumpat(14) = mgamm

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine epotlp(ilp,epot)
*                                                                      *
*        Last Revised:     2006 5 11                                   *
*                                                                      *
*        Purpose:                                                      *
*              to calculate potential energy of a leading particle     *
*                                                                      *
*        Variables:                                                    *
*              epot        : total potential energy                    *
*                                                                      *
************************************************************************
      use QMD_COOD2_MOD, only : rha, rhe, rhc

      implicit real*8(a-h,o-z)

      parameter ( nnnn = 800 )

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson
!$OMP THREADPRIVATE(/vriab0/)

      common /coodrp/ r(5,nnn),  p(6,nnn)
!$OMP THREADPRIVATE(/coodrp/)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)

      common /qmdscmpar/ iesc(nnnn)
!$OMP THREADPRIVATE(/qmdscmpar/)

      common /poten1/ gamm, c0, c3, cs, cl, wl

*-----------------------------------------------------------------------

      i = ilp

      epot1 = 0.d0
      epot3 = 0.d0
      epots = 0.d0
      epotc = 0.d0

      do 110 j = 1, massal

         if( j .eq. i .or. iesc(j) .eq. 0 ) goto 110

         epot1 = epot1 + rha(j,i)
         epotc = epotc + rhe(j,i)
         epots = epots + rha(j,i) * inuc(j) * inuc(i)
     &                  * ( 1.d0 - 2.d0 * dble(abs(ichg(j)-ichg(i))) )

  110 continue


      epot3 = epot1 ** gamm

*-----------------------------------------------------------------------

      epot = c0 * epot1 + c3 * epot3 + cs * epots + cl * epotc

*-----------------------------------------------------------------------

      return
      end
