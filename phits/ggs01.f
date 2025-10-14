************************************************************************
*                                                                      *
      subroutine initsk
*                                                                      *
*       initialize for the current task.                               *
*                                                                      *
*       Last modified by K.Niita on 2009/10/06                         *
*                                                                      *
************************************************************************

      use GGBANKMOD !FURUTA
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      kdb  = 0

      do 60 nlaj = 0,mlaj+mxa-1
   60    if(laj(nlaj+1).eq.0) goto 70 !FURUTA

   70    continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine chkcll(i1,m,j)
*                                                                      *
*       check the sense of xyz with respect to the surfaces of cell i1.*
*       m=0  simple cell:  return j=0 if xyz is in the cell.  otherwise*
*            return j=first bad surface.  ignore surface jsu.          *
*            complicated cell:  return j=0 if xyz is in the cell and,  *
*            if jsu is a surface of the cell, uvw is directed into the *
*            cell.  otherwise return j=1.                              *
*       m=1  for surface segmenting in tally.  return j=first good     *
*            surface.  if none, return j=0.                            *
*       m=2  like m=0, but cell is treated as complicated regardless.  *
*       m=3  preprocessor for track in places where lgc is not already *
*            loaded.  load logical elements into lgc.                  *
*       if jsu is a surface of the cell and either m>0 or the cell is  *
*       complicated, uvw must be appropriately defined.                *
*                                                                      *
*       Last modified by K.Niita on 2009/09/30                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

      real(8) xxx0,uuu0 !FURUTA

      if( i1.eq.0 ) return

*-----------------------------------------------------------------------
*     do all surfaces of cell i1 or until mm is matched.
*-----------------------------------------------------------------------

         j2 = abs(lca(i1))
         il = j2-lca(i1)+m
         mm = m
         if(il.ne.0.and.m.eq.0) mm = 2

      do 340 j = j2, abs(lca(i1+1))-1

         jk = abs(lja(j))
         if(il.eq.0) goto 30
         if(jk.gt.1000000) goto 20
         if(jk.ne.jsu) goto 40

   10    lgc(j-j2+1) = 1
         if(dble(lja(j))*angl(i1).lt.0.d0) lgc(j-j2+1) = 0 !FURUTA20230706
         goto 340

   20    lgc(j-j2+1) = jk
         goto 340

   30    if(jk.eq.jsu) goto 340

   40    k = kst(jk)
         i = lsc(jk)+1

         if ( k.ge.24 .and. k.le.26 .and. jtr(jk).ne.0 ) then
           if( jtr(jk).gt.0 ) then
             do jt = 1, mxtr
               if( abs(trf(1,jt)).eq.jtr(jk) ) then
                 m1 = jt
                 exit
               end if
             end do
           else
             m1 = -jtr(jk)/10000
             m2 = -jtr(jk) - m1*10000
           end if
           xxx_cp = xxx
           yyy_cp = yyy
           zzz_cp = zzz
           call trnsxx(xxx_cp,yyy_cp,zzz_cp,xxx,yyy,zzz,m1)
           if( jtr(jk).lt.0 ) then
             xxx_it = xxx
             yyy_it = yyy
             zzz_it = zzz
             call trnsxx(xxx_it,yyy_it,zzz_it,xxx,yyy,zzz,m2)
           end if
         end if
      goto (50,60,60,60,70,80,90,100,110,120,130,140,150,160,170,180,
     &      190,200,230,240,250,280,290,300,310,320) k

*-----------------------------------------------------------------------
* >>>>>  planes:  p, px, py, pz
*        if jsu<0, a plane may be coincident with jsu at a higher level.
*-----------------------------------------------------------------------

   50    t4 = scf(i)*xxx+scf(i+1)*yyy+scf(i+2)*zzz-scf(i+3)
         if(jsu.ge.0) goto 330
         if(ksc(jk).ne.ksc(-jsu)) goto 335

         t5 = scf(i)*uuu+scf(i+1)*vvv+scf(i+2)*www
         if(abs(t4).le.coincd*abs(t5)) t4 = t5
         goto 335

 60      select case(k)                             !FURUTA
         case(2)                                    !FURUTA
           xxx0 = xxx                               !FURUTA
           uuu0 = uuu                               !FURUTA
         case(3)                                    !FURUTA
           xxx0 = yyy                               !FURUTA
           uuu0 = vvv                               !FURUTA
         case(4)                                    !FURUTA
           xxx0 = zzz                               !FURUTA
           uuu0 = www                               !FURUTA
         case default                               !FURUTA
           write(*,*)'ERROR in chkcll: SELECT'      !FURUTA
           stop                                     !FURUTA
         end select                                 !FURUTA
         t4 = xxx0-scf(i)                           !FURUTA
         if(jsu.ge.0) goto 330                      !FURUTA
         if(ksc(jk).ne.ksc(-jsu)) goto 335 !FURUTA
         if(abs(t4).le.coincd*abs(uuu0)) t4 = uuu0  !FURUTA
         goto 335

*-----------------------------------------------------------------------
* >>>>>  spheres:  so, s, sx, sy, sz
*-----------------------------------------------------------------------

   70    t4 = xxx**2+yyy**2+zzz**2-scf(i)
         goto 330
   80    t4 = (xxx-scf(i))**2+(yyy-scf(i+1))**2
     &      + (zzz-scf(i+2))**2-scf(i+3)
         goto 330
   90    t4 = (xxx-scf(i))**2+yyy**2+zzz**2-scf(i+1)
         goto 330
  100    t4 = xxx**2+(yyy-scf(i))**2+zzz**2-scf(i+1)
         goto 330
  110    t4 = xxx**2+yyy**2+(zzz-scf(i))**2-scf(i+1)
         goto 330

*-----------------------------------------------------------------------
* >>>>>  cylinders:  c/x, c/y, c/z, cx, cy, cz
*-----------------------------------------------------------------------

  120    t4 = (yyy-scf(i))**2+(zzz-scf(i+1))**2-scf(i+2)
         goto 330

  130    t4 = (xxx-scf(i))**2+(zzz-scf(i+1))**2-scf(i+2)
         goto 330

  140    t4 = (xxx-scf(i))**2+(yyy-scf(i+1))**2-scf(i+2)
         goto 330

  150    t4 = yyy**2+zzz**2-scf(i)
         goto 330

  160    t4 = xxx**2+zzz**2-scf(i)
         goto 330

  170    t4 = xxx**2+yyy**2-scf(i)
         goto 330

*-----------------------------------------------------------------------
* >>>>>  cones:  k/x, k/y, k/z, kx, ky, kz
*-----------------------------------------------------------------------

  180    t1 = (yyy-scf(i+1))**2+(zzz-scf(i+2))**2
         t2 = xxx-scf(i)
         goto 210

  190    t1 = (xxx-scf(i))**2+(zzz-scf(i+2))**2
         t2 = yyy-scf(i+1)
         goto 210

  200    t1 = (xxx-scf(i))**2+(yyy-scf(i+1))**2
         t2 = zzz-scf(i+2)

  210    if(scf(i+4).eq.0.) goto 220
         t4 = t1-sign((scf(i+4)*t2)**2,scf(i+4)*t2)
         goto 330

  220    t4 = t1-scf(i+3)*t2**2
         goto 330

  230    t1 = yyy**2+zzz**2
         t2 = xxx-scf(i)
         goto 260

  240    t1 = xxx**2+zzz**2
         t2 = yyy-scf(i)
         goto 260

  250    t1 = xxx**2+yyy**2
         t2 = zzz-scf(i)

  260    if(scf(i+2).eq.0.) goto 270
         t4 = t1-sign((scf(i+2)*t2)**2,scf(i+2)*t2)
         goto 330

  270    t4 = t1-scf(i+1)*t2**2
         goto 330

*-----------------------------------------------------------------------
* >>>>>  special quadratic:  sq
*-----------------------------------------------------------------------

  280    t1 = xxx-scf(i+7)
         t2 = yyy-scf(i+8)
         t3 = zzz-scf(i+9)
         t4 = scf(i)*t1**2+scf(i+1)*t2**2+scf(i+2)*t3**2
     &      + 2.*(scf(i+3)*t1+scf(i+4)*t2+scf(i+5)*t3)+scf(i+6)
         goto 330

*-----------------------------------------------------------------------
* >>>>>  general quadratic:  gq
*-----------------------------------------------------------------------

  290    t4 = (scf(i)*xxx+scf(i+3)*yyy+scf(i+6))*xxx
     &      + (scf(i+1)*yyy+scf(i+4)*zzz+scf(i+7))*yyy
     &      + (scf(i+2)*zzz+scf(i+5)*xxx+scf(i+8))*zzz+scf(i+9)
         goto 330

*-----------------------------------------------------------------------
* >>>>>  tori:  tx, ty, tz
*-----------------------------------------------------------------------

  300    t4 = (sqrt((yyy-scf(i+1))**2+(zzz-scf(i+2))**2)-scf(i+3))**2
     &      + scf(i+6)*(xxx-scf(i))**2-scf(i+5)**2
         goto 325 ! frtati 2022/04/08

  310    t4 = (sqrt((xxx-scf(i))**2+(zzz-scf(i+2))**2)-scf(i+3))**2
     &      + scf(i+6)*(yyy-scf(i+1))**2-scf(i+5)**2
         goto 325 ! frtati 2022/04/08

  320    t4 = (sqrt((xxx-scf(i))**2+(yyy-scf(i+1))**2)-scf(i+3))**2
     &      + scf(i+6)*(zzz-scf(i+2))**2-scf(i+5)**2

  325    if( jtr(jk).ne.0 ) then
           xxx = xxx_cp
           yyy = yyy_cp
           zzz = zzz_cp
         end if

*-----------------------------------------------------------------------
*        set logical value into lgc.  return if it matches mm.
*-----------------------------------------------------------------------

 330     continue

         if(t4.eq.0.d0)then !FURUTA20251106
          jsu0=jsu
          jsu=abs(lja(j))
          lgc(j-j2+1) = 1
          if(dble(lja(j))*angl(i1).lt.0.d0) lgc(j-j2+1) = 0
          jsu=jsu0
          goto 340
         endif

  335    lgc(j-j2+1) = 1
         if(lja(j)*t4.lt.0.d0) lgc(j-j2+1) = 0
         if(lgc(j-j2+1).eq.mm) return

  340    continue

         j = 0
         if(mm.ne.2) return
         j = 1-lgeval(lgc,abs(lca(i1+1))-j2)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine fsurf(ih)
*                                                                      *
*       find the next bounding surface of cell ih or of a higher-level *
*       cell, and calculate its distance dls.                          *
*       if ih>mxa, for a set of tally segmenting surfaces, just return *
*       the intersections in dti.                                      *
*       this subroutine is used a lot.  some of the funny-looking      *
*       coding has been done that way for speed.                       *
*                                                                      *
*       Last modified by K.Niita on 2009/09/30                         *
*                                                                      *
************************************************************************
      use TETRAMOD, only: tetracald,tetrachk,itetragshow
      use GGBANKMOD !FURUTA
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      dimension dl(0:mxlv),jp(0:mxlv)

      real(8) xxx0,uuu0 !FURUTA
*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk
      integer ierr,ires,ic0,itet,iii0
      real(8) xx0(6),x0(3),u0(3)
*-----------------------------------------------------------------------
*        set up to do the cell in the current level.
*-----------------------------------------------------------------------


         ll = lev
         ic = ih

*-----------------------------------------------------------------------
*        find the distance to boundary, dl, in this cell.
*-----------------------------------------------------------------------

   10    dl(ll) = huge
         nlt = 0

*-----------------------------------------------------------------------
*        find tetrahedron surfaces (FURUTA20150714 TETRA)
*-----------------------------------------------------------------------
         if(junf.ne.0)then
          if(abs(lat(1,ic)).eq.3)then
           ierr=0
           call tetrabox(0,kkk,xx0,ierr)
           x0(1)=xxx
           x0(2)=yyy
           x0(3)=zzz
           u0(1)=uuu
           u0(2)=vvv
           u0(3)=www
           ic0=ic
           if(itetragshow.eq.0)then !FURUTA20221209
            itet=kkk-10000
            iii0=iii !FURUTA20171020
            call tetrachk(xx0,x0,u0,coincd,itet,ic,iii,ires,ierr)
            if(ierr.ne.0)mark=-2
            if(iii.lt.0)ic=itetcl(itet)
            if(ic0.eq.ic.and.ires.eq.1)ih=ic
            if(iii.ne.iii0)then
             icl=ic              !FURUTA20171020
             call tetragetmat(matnum) !FURUTA20200612
             call tetrasetmat(matnum) !FURUTA20200612
            endif
           endif
           call tetracald(xx0,x0,u0,kkk-10000,iii,d1,j,ierr)
           if(ierr.gt.0)mark=-2
           jp(ll) = j
           dl(ll) = d1
           goto 370
          endif
         endif

*-----------------------------------------------------------------------
*        calculate all intersections with all surfaces of this cell.
*-----------------------------------------------------------------------

         j2 = abs(lca(ic))
         il = j2-lca(ic)
         j3 = abs(lca(ic+1))-1
         do 20 j1 = j2, j3

*-----------------------------------------------------------------------
*        ignore logical operators and duplicate surfaces.
*-----------------------------------------------------------------------

         if(lcaj(j1).lt.0) goto 20 !FURUTA
         j = abs(lja(j1))
         k = kst(j)
         if(k.le.4.and.j.eq.jsu) goto 20
         i = lsc(j)+1
         d1 = -1.
         d2 = -1.
         tk = 0.
         t9 = 1.e10

         goto (30,40,40,40,50,60,70,80,90,100,110,120,130,140,150,160,
     &         170,180,190,200,210,220,230,240,240,240) k

*-----------------------------------------------------------------------
* >>>>>  planes:  p, px, py, pz
*        if jsu<0 check for coincident surface at higher level.
*-----------------------------------------------------------------------

   30    t1 = uuu*scf(i)+vvv*scf(i+1)+www*scf(i+2)
         if(t1.eq.0.) goto 20
         d1 = (scf(i+3)-xxx*scf(i)-yyy*scf(i+1)-zzz*scf(i+2))/t1
         if(d1.le.0.) goto 20
         if(jsu.ge.0) goto 300
         if(d1.gt.coincd) goto 300
         if(ksc(j).ne.ksc(-jsu)) goto 300
         if(t1*lja(j1).gt.0..or.il.ne.0) goto 20
         goto 300

 40      select case(k)                                 !FURUTA
         case(2)                                        !FURUTA
           xxx0 = xxx                                   !FURUTA
           uuu0 = uuu                                   !FURUTA
         case(3)                                        !FURUTA
           xxx0 = yyy                                   !FURUTA
           uuu0 = vvv                                   !FURUTA
         case(4)                                        !FURUTA
           xxx0 = zzz                                   !FURUTA
           uuu0 = www                                   !FURUTA
         case default                                   !FURUTA
           write(*,*)'ERROR in fsurf: SELECT'           !FURUTA
           stop                                         !FURUTA
         end select                                     !FURUTA
         if(uuu0.eq.0.) goto 20                         !FURUTA
         d1 = (scf(i)-xxx0)/uuu0                        !FURUTA
         if(d1.le.0.) goto 20                           !FURUTA
         if(jsu.ge.0) goto 300                          !FURUTA
         if(d1.gt.coincd) goto 300                      !FURUTA
         if(ksc(j).ne.ksc(-jsu)) goto 300      !FURUTA
         if(uuu0*lja(j1).gt.0..or.il.ne.0) goto 20 !FURUTA
         goto 300

*-----------------------------------------------------------------------
* >>>>>  spheres:  so, s, sx, sy, sz
*-----------------------------------------------------------------------

   50    a1 = xxx*uuu+yyy*vvv+zzz*www
         b1 = xxx**2+yyy**2+zzz**2-scf(i)
         goto 260

   60    t1 = xxx-scf(i)
         t2 = yyy-scf(i+1)
         t3 = zzz-scf(i+2)
         a1 = t1*uuu+t2*vvv+t3*www
         b1 = t1**2+t2**2+t3**2-scf(i+3)
         goto 260

   70    t1 = xxx-scf(i)
         a1 = t1*uuu+yyy*vvv+zzz*www
         b1 = t1**2+yyy**2+zzz**2-scf(i+1)
         goto 260

   80    t2 = yyy-scf(i)
         a1 = xxx*uuu+t2*vvv+zzz*www
         b1 = xxx**2+t2**2+zzz**2-scf(i+1)
         goto 260

   90    t3 = zzz-scf(i)
         a1 = xxx*uuu+yyy*vvv+t3*www
         b1 = xxx**2+yyy**2+t3**2-scf(i+1)
         goto 260

*-----------------------------------------------------------------------
* >>>>>  cylinders:  c/x, c/y, c/z, cx, cy, cz
*-----------------------------------------------------------------------

  100    t1 = vvv**2+www**2
         if(t1.eq.0.) goto 20
         t1 = 1./t1
         t2 = yyy-scf(i)
         t3 = zzz-scf(i+1)
         a1 = (t2*vvv+t3*www)*t1
         b1 = (t2**2+t3**2-scf(i+2))*t1
         goto 260

  110    t1 = uuu**2+www**2
         if(t1.eq.0.) goto 20
         t1 = 1./t1
         t2 = xxx-scf(i)
         t3 = zzz-scf(i+1)
         a1 = (t2*uuu+t3*www)*t1
         b1 = (t2**2+t3**2-scf(i+2))*t1
         goto 260

  120    t1 = uuu**2+vvv**2
         if(t1.eq.0.) goto 20
         t1 = 1./t1
         t2 = xxx-scf(i)
         t3 = yyy-scf(i+1)
         a1 = (t2*uuu+t3*vvv)*t1
         b1 = (t2**2+t3**2-scf(i+2))*t1
         goto 260

  130    t1 = vvv**2+www**2
         if(t1.eq.0.) goto 20
         t1 = 1./t1
         a1 = (yyy*vvv+zzz*www)*t1
         b1 = (yyy**2+zzz**2-scf(i))*t1
         goto 260

  140    t1=uuu**2+www**2
         if(t1.eq.0.)go to 20
         t1=1./t1
         a1=(xxx*uuu+zzz*www)*t1
         b1=(xxx**2+zzz**2-scf(i))*t1
         go to 260

  150    t1 = uuu**2+vvv**2
         if(t1.eq.0.) goto 20
         t1 = 1./t1
         a1 = (xxx*uuu+yyy*vvv)*t1
         b1 = (xxx**2+yyy**2-scf(i))*t1
         goto 260

*-----------------------------------------------------------------------
* >>>>>  cones:  k/x, k/y, k/z, kx, ky, kz
*-----------------------------------------------------------------------

  160    t4 = vvv**2+www**2-scf(i+3)*uuu**2
         if(abs(t4).gt.epss) t9 = 1./t4
         t1 = yyy-scf(i+1)
         t2 = zzz-scf(i+2)
         t3 = xxx-scf(i)
         uu = uuu
         tk = scf(i+4)
         a1 = (t1*vvv+t2*www-t3*uuu*scf(i+3))*t9
         b1 = (t1**2+t2**2-scf(i+3)*t3**2)*t9
         if(t9.ne.1.e10) goto 260
         goto 250

  170    t4 = uuu**2+www**2-scf(i+3)*vvv**2
         if(abs(t4).gt.epss) t9 = 1./t4
         t1 = xxx-scf(i)
         t2 = zzz-scf(i+2)
         t3 = yyy-scf(i+1)
         uu = vvv
         tk = scf(i+4)
         a1 = (t1*uuu+t2*www-t3*vvv*scf(i+3))*t9
         b1 = (t1**2+t2**2-scf(i+3)*t3**2)*t9
         if(t9.ne.1.e10) goto 260
         goto 250

  180    t4 = uuu**2+vvv**2-scf(i+3)*www**2
         if(abs(t4).gt.epss) t9 = 1./t4
         t1 = yyy-scf(i+1)
         t2 = xxx-scf(i)
         t3 = zzz-scf(i+2)
         uu = www
         tk = scf(i+4)
         a1 = (t1*vvv+t2*uuu-t3*www*scf(i+3))*t9
         b1 = (t1**2+t2**2-scf(i+3)*t3**2)*t9
         if(t9.ne.1.e10) goto 260
         goto 250

  190    t4 = vvv**2+www**2-scf(i+1)*uuu**2
         if(abs(t4).gt.epss) t9 = 1./t4
         t3 = xxx-scf(i)
         uu = uuu
         tk = scf(i+2)
         a1 = (zzz*www+yyy*vvv-t3*uuu*scf(i+1))*t9
         b1 = (zzz**2+yyy**2-scf(i+1)*t3**2)*t9
         if(t9.ne.1.e10) goto 260
         goto 250

  200    t4 = uuu**2+www**2-scf(i+1)*vvv**2
         if(abs(t4).gt.epss) t9 = 1./t4
         t3 = yyy-scf(i)
         uu = vvv
         tk = scf(i+2)
         a1 = (xxx*uuu+zzz*www-t3*vvv*scf(i+1))*t9
         b1 = (xxx**2+zzz**2-scf(i+1)*t3**2)*t9
         if(t9.ne.1.e10) goto 260
         goto 250

  210    t4 = uuu**2+vvv**2-scf(i+1)*www**2
         if(abs(t4).gt.epss) t9 = 1./t4
         t3 = zzz-scf(i)
         uu = www
         tk = scf(i+2)
         a1 = (xxx*uuu+yyy*vvv-t3*www*scf(i+1))*t9
         b1 = (xxx**2+yyy**2-scf(i+1)*t3**2)*t9
         if(t9.ne.1.e10) goto 260
         goto 250

*-----------------------------------------------------------------------
* >>>>>  special quadratic:  sq
*-----------------------------------------------------------------------

  220    t1 = scf(i)*uuu**2+scf(i+1)*vvv**2+scf(i+2)*www**2
         if(abs(t1).gt.epss) t9 = 1./t1
         t2 = xxx-scf(i+7)
         t3 = yyy-scf(i+8)
         t4 = zzz-scf(i+9)
         t5 = scf(i)*t2+scf(i+3)
         t6 = scf(i+1)*t3+scf(i+4)
         t7 = scf(i+2)*t4+scf(i+5)
         a1 = (t5*uuu+t6*vvv+t7*www)*t9
         b1 = ((t5+scf(i+3))*t2+(t6+scf(i+4))*t3+(t7+scf(i+5))*t4+
     &         scf(i+6))*t9
         if(t9.ne.1.e10) goto 260
         goto 250

*-----------------------------------------------------------------------
* >>>>>  general quadratic:  gq
*-----------------------------------------------------------------------

  230    t1 = scf(i)*uuu**2+scf(i+1)*vvv**2+scf(i+2)*www**2+
     &    (scf(i+3)*uuu+scf(i+4)*www)*vvv+scf(i+5)*uuu*www
         if(abs(t1).gt.epss) t9 = 1./t1
         a1 = ((2.*scf(i)*xxx+scf(i+6))*uuu
     &      + scf(i+3)*(xxx*vvv+yyy*uuu)
     &      + (2.*scf(i+1)*yyy+scf(i+7))*vvv
     &      + scf(i+4)*(zzz*vvv+yyy*www)
     &      + (2.*scf(i+2)*zzz+scf(i+8))*www
     &      + scf(i+5)*(xxx*www+zzz*uuu))*.5*t9
         b1 = ((scf(i)*xxx+scf(i+3)*yyy+scf(i+6))*xxx
     &      + (scf(i+1)*yyy+scf(i+4)*zzz+scf(i+7))*yyy
     &      + (scf(i+2)*zzz+scf(i+5)*xxx+scf(i+8))*zzz+scf(i+9))*t9
         if(t9.ne.1.e10) goto 260
         goto 250

*-----------------------------------------------------------------------
* >>>>>  tori:  tx, ty, tz
*-----------------------------------------------------------------------

  240     if( jtr(j).ne.0 ) then
           if( jtr(j).gt.0 ) then
             do jt = 1, mxtr
               if( abs(trf(1,jt)).eq.jtr(j) ) then
                 m1 = jt
                 exit
               end if
             end do
           else
             m1 = -jtr(j)/10000
             m2 = -jtr(j) - m1*10000
           end if
           xxx_cp = xxx
           yyy_cp = yyy
           zzz_cp = zzz
           uuu_cp = uuu
           vvv_cp = vvv
           www_cp = www
           call trnsxx(xxx_cp,yyy_cp,zzz_cp,xxx,yyy,zzz,m1)
           call trnsuu(uuu_cp,vvv_cp,www_cp,uuu,vvv,www,m1)
           if( jtr(j).lt.0 ) then
             xxx_it = xxx
             yyy_it = yyy
             zzz_it = zzz
             uuu_it = uuu
             vvv_it = vvv
             www_it = www
             call trnsxx(xxx_it,yyy_it,zzz_it,xxx,yyy,zzz,m2)
             call trnsuu(uuu_it,vvv_it,www_it,uuu,vvv,www,m2)
           end if
           call inttor(ic,i,j,d1)
           xxx = xxx_cp
           yyy = yyy_cp
           zzz = zzz_cp
           uuu = uuu_cp
           vvv = vvv_cp
           www = www_cp
         else
           call inttor(ic,i,j,d1)
         end if

         goto 290

*-----------------------------------------------------------------------
*        calculate a single distance.  (rare case)
*-----------------------------------------------------------------------

  250    if(abs(a1).le.epss.or.j.eq.jsu) goto 20
         d1 = -.5*b1/a1
         if(d1.le.0.) goto 20
         goto 280

*-----------------------------------------------------------------------
*        calculate a pair of distances.
*-----------------------------------------------------------------------

  260    t1=a1**2-b1
         if(t1.le.0.) goto 20
         if(j.ne.jsu) goto 270
         d1 = -2.*a1
         if(d1.le.0.) goto 20
         goto 280

  270    t2 = sqrt(t1)
         d1 = -a1+t2
         if(abs(d1).lt.t2*1.0d-12) d1=0.0d0 ! T.Sato 2023/07/16 a1 and -t2 are so close to, so should be 0. bug-fix due to optimization (sqrt is not precise enough?)
         if(d1.le.0.) goto 20
         d2 = -a1-t2
         if(abs(d2).lt.t2*1.0d-12) d2=0.0d0 ! T.Sato 2023/07/16 a1 and -t2 are so close to, so should be 0. bug-fix due to optimization (sqrt is not precise enough?)
*-----------------------------------------------------------------------
*        reject the distance to the wrong sheet of a cone.
*-----------------------------------------------------------------------

  280    if(tk.eq.0.) goto 300
         if(tk*(t3+d1*uu).lt.0.) d1 = -1.
         if(tk*(t3+d2*uu).lt.0.) d2 = -1.
         d1 = max(d1,d2)
  290    if(d1.le.0.) goto 20

*-----------------------------------------------------------------------
*        simple cell case:  select the smallest positive distance.
*-----------------------------------------------------------------------

  300    if(il.ne.0) goto 310
         if(d2.gt.0.) d1 = d2
         if(d1.ge.dl(ll)) goto 20
         jp(ll) = j
         dl(ll) = d1
         goto 20
*
*-----------------------------------------------------------------------
*        complicated cell case:  put positive roots in dti.
*-----------------------------------------------------------------------

  310    nlt = nlt+1
         dti(nlt) = d1
         iti(nlt) = j
         if(d2.le.0..or.d1.eq.d2) goto 20
         nlt = nlt+1
         dti(nlt) = d2
         iti(nlt) = j
         goto 20

   20 continue

*-----------------------------------------------------------------------
*        if complicated cell, find the distance to the cell boundary.
*-----------------------------------------------------------------------

  320    if(il.eq.0.and.junf.eq.0) goto 460
         if(il.eq.0) goto 380
         if(ih.gt.mxa) return
         if(nlt.eq.0) goto 370

  330    n1 = 1
      do 340 m = 1, nlt
  340    if(dti(m).lt.dti(n1)) n1 = m

         dl(ll) = dti(n1)
         if(dl(ll).eq.huge) goto 370

      do 350 j1 = j2, j3
  350    if(abs(lja(j1)).eq.iti(n1)) lgc(j1-j2+1) = 1-lgc(j1-j2+1)

         dti(n1) = huge
         if(lgeval(lgc,j3-j2+1).ne.0) goto 330

  360    jp(ll) = iti(n1)

*-----------------------------------------------------------------------
*        prepare to do the cell in the next higher level.
*        set jsu=-jsu for higher levels to check coincidence surfaces.
*-----------------------------------------------------------------------

  370    if(junf.eq.0) goto 460
  380    if(jun(ic).le.0) goto 420
         if(ll.ne.lev) goto 400

         udt(1,lev) = xxx !FURUTA
         udt(2,lev) = yyy !FURUTA
         udt(3,lev) = zzz !FURUTA
         udt(4,lev) = uuu !FURUTA
         udt(5,lev) = vvv !FURUTA
         udt(6,lev) = www !FURUTA

         js = jsu
         if(jsu.eq.0) goto 400
         if(jsu.le.mxj)then !FURUTA20221104
          if(ksc(jsu).ne.0) jsu = -js
         endif              !FURUTA20221104
  400    ll = ll - 1

         xxx = udt(1,ll) !FURUTA
         yyy = udt(2,ll) !FURUTA
         zzz = udt(3,ll) !FURUTA
         uuu = udt(4,ll) !FURUTA
         vvv = udt(5,ll) !FURUTA
         www = udt(6,ll) !FURUTA

         ic = udt(7,ll)
         if(lca(ic).lt.0) call chkcll(ic,3,j)
         goto 10

  420    if(ll.eq.lev) goto 440

         xxx = udt(1,lev) !FURUTA
         yyy = udt(2,lev) !FURUTA
         zzz = udt(3,lev) !FURUTA
         uuu = udt(4,lev) !FURUTA
         vvv = udt(5,lev) !FURUTA
         www = udt(6,lev) !FURUTA

         jsu = js

*-----------------------------------------------------------------------
*        find the level where the distance to boundary is least.
*-----------------------------------------------------------------------

  440    levp = ll

      do 450 i = ll+1,lev
         if(dl(i).ge.dl(levp)) goto 450
         if(dl(i)+coincd.lt.dl(levp)) goto 445
         if(jp(i).le.mxj)then !FURUTA20221104
          if(ksc(jp(i)).eq.0) goto 445
          if(ksc(jp(i)).eq.ksc(jp(levp))) goto 450
         endif                !FURUTA20221104
  445    levp = i
  450 continue

  460    dls = dl(levp)
         jap = jp(levp)

*-----------------------------------------------------------------------
*        return with error indicator if no intersection was found.
*-----------------------------------------------------------------------

         if(dls.eq.huge) kdb = 2

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine fcell(io,ierr,cs)
*                                                                      *
*       find the new cell, iap, into which the particle has crossed    *
*       from cell icl through surface jsu.  if cs .ne. 0, return the   *
*       angle between the track and the normal to surface jsu.         *
*                                                                      *
*       ierr  : error code                                             *
*         = 1 : the surface crossed is not a surface of this cell      *
*         = 2 : zero lattice element hit.                              *
*                                                                      *
*       Last modified by K.Niita on 2009/09/30                         *
*                                                                      *
************************************************************************
      use TETRAMOD, only: tetranext,tetrafnd
      use GGBANKMOD !FURUTA
      use LAFDATAMOD !FURUTA20201127
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      common /rcomon/ rcasc

*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk
      real(8) xx0(6),xxx0(3)
*-----------------------------------------------------------------------

         ierr = 0

*-----------------------------------------------------------------------
*      shift to the higher level if one was found in track.
*-----------------------------------------------------------------------

         if(levp.eq.lev) goto 10

         lev = levp
         xxx = udt(1,lev)
         yyy = udt(2,lev)
         zzz = udt(3,lev)
         uuu = udt(4,lev)
         vvv = udt(5,lev)
         www = udt(6,lev)
         icl = udt(7,lev)
         iii = udt(8,lev)
         jjj = udt(9,lev)
         kkk = udt(10,lev)

   10    if(cs.ne.0.) cs = angl(icl) !FURUTA20230706

*-----------------------------------------------------------------------
*     find surface crossed among the bounding surfaces of cell icl.
*-----------------------------------------------------------------------

      if(junf.ne.0)then
        if(abs(lat(1,icl)).eq.3)goto 35 !FURUTA20150714 TETRA
      end if

      do 20 i1 = abs(lca(icl)),abs(lca(icl+1))-1
   20    if(abs(lja(i1)).eq.jsu)go to 30


         ierr = 1
         return

   30    j2 = lcaj(i1)          !FURUTA
         j3 = abs(lcaj(i1+1))-1 !FURUTA

   35    continue !FURUTA20150714 TETRA

*-----------------------------------------------------------------------
*     move into the next lattice element if icl is a lattice.
*-----------------------------------------------------------------------

         if(junf.eq.0) goto 40
         if(lat(1,icl).eq.0) goto 40

cFURUTA20150714 TETRA---------------------------------------------------
         if(abs(lat(1,icl)).ne.3)then
          call intlat(i1)
         else
          j2=0 !FURUTA20160607 BUGFIX
          call tetrabox(io,kkk,xx0,ierr)
          xxx0(1)=xx0(1)
          xxx0(2)=xx0(3)
          xxx0(3)=xx0(5)
          call tetranext(xxx0,xxx,yyy,zzz,uuu,vvv,www,coincd,kkk-10000,
     &         jsu,iii,jjj,iap,ierr)
          if(jjj.lt.0)iap=itetcl(kkk-10000)
          if(jjj.eq.0)iap=-12345
         endif
*-----------------------------------------------------------------------

         if(iap.eq.-12345) goto 380
         j1 = j2-1
         goto 210

*-----------------------------------------------------------------------
*     look for the new cell in the laj list.
*     after the first 50 histories, omit the cell check if the cell
*     is simple and alone on the other side and if krflg indicates
*     that no geometry error has just occurred.
*-----------------------------------------------------------------------

   40    if(j3.lt.j2) goto 70
         if(laj(j2).gt.0) goto 50
         iap = abs(laj(j2))
         j1 = j2
         if( rcasc.gt.50.0 .and. krflg .lt. 2 ) goto 210

   50    iap = 0
      do 60 j1 = j2,j3
         if(iap.eq.1000001.and.laj(j1).ne.1000002) goto 60
         iap = abs(laj(j1))
         if(iap.gt.1000000) goto 60

         call chkcll(iap,0,j)

         if(j.eq.0) goto 210
   60 continue

*-----------------------------------------------------------------------
*     search the overflow region at the end of laj.
*-----------------------------------------------------------------------
   70 do 90 j1 = mlaj+1,nlaj !FURUTA
         iap = laj(j1)
         if(iap.eq.icl) goto 90
         if(junf.eq.0) goto 80
         if(abs(jun(iap)).ne.abs(jun(icl))) goto 90

   80    call chkcll(iap,0,j)

         if(j.eq.0) goto 210
   90 continue

*-----------------------------------------------------------------------
*     look for the new cell in the entire list of cells.
*-----------------------------------------------------------------------

      do 130 iap = 1,mxa
         if(iap.eq.icl) goto 130
         if(junf.eq.0) goto 100
         if(abs(jun(iap)).ne.abs(jun(icl))) goto 130
  100 do 110 i4 = abs(lca(iap)),abs(lca(iap+1))-1
  110    if(abs(lja(i4)).eq.jsu) goto 120
         goto 130

  120    call chkcll(iap,0,j)

         if(j.eq.0) goto 140
  130 continue
cFURUTA20250410 to avoid lost particle at contact point of surfaces-----
      outer: do iap=1,mxa
       if(iap.eq.icl)cycle
       if(junf.ne.0)then
          if(abs(jun(iap)).ne.abs(jun(icl)))cycle
       endif
       do i4 = abs(lca(iap)),abs(lca(iap+1))-1
        if(abs(lja(i4)).eq.jsu)cycle outer
       enddo
       call chkcll(iap,0,j)
       if(j.eq.0) goto  140
      enddo outer
c-----------------------------------------------------------------------      
         iap = 1
         kdb = 1
         return

*-----------------------------------------------------------------------
*     insert the new cell in laj and increment the lcaj pointers.
*-----------------------------------------------------------------------

  140    if(nlaj.ge.mlaj) goto 180
         j1 = j3+1
         j3 = j3+1

      do 150 iz = j1,nlaj !FURUTA
         i = j1+nlaj-iz   !FURUTA
  150    laj(i+1) = laj(i)

         nlaj = nlaj+1
         laj(j1) = iap

      do 160 j = i1,nlja
  160    lcaj(j+1) = sign(abs(lcaj(j+1))+1,lcaj(j+1)) !FURUTA

*-----------------------------------------------------------------------
*     mark the new other-side cell if it is simple and alone.
*-----------------------------------------------------------------------

         if(lca(iap).lt.0) goto 210

      do 170 j = 1,nlja
  170    if((j.lt.lca(iap).or.j.ge.abs(lca(iap+1))).and.
     &    lja(j).eq.lja(i4)) goto 210
         laj(j1) = -iap
         goto 210

*-----------------------------------------------------------------------
*     add the new cell to the overflow region at the end of laj.
*-----------------------------------------------------------------------

  180    j = 0
         nlaj = nlaj+1
         laj(nlaj) = iap !FURUTA

*-----------------------------------------------------------------------
*      move into the universe, if any, that fills the new cell.
*-----------------------------------------------------------------------

  210    if(junf.eq.0) return
         if(iap == 0) return !20220909
         if(mfl(1,iap).eq.0) return
         if(ksc(jsu).ne.0) jsu = -jsu

  215    j = -mfl(1,iap)
         if(j.gt.0) goto 220
         ju = -j

         call dwnlev(iap)

         if(mfl(2,iap).eq.0) goto 230
         iap = mfl(2,iap)
         goto 370

  220    n = 3+iii-laf(1,j+1)+laf(1,j+2)*(jjj-laf(2,j+1)+laf(2,j+2)
     &        * (kkk-laf(3,j+1)))

         if(laf(1,j+n).eq.abs(jun(iap))) goto 372
         if(laf(1,j+n).eq.0) goto 380
         ju = laf(1,j+n)

         call dwnlev(iap)

         if(laf(2,j+n).eq.0) goto 230
         iap = laf(2,j+n)

         goto 370
c-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        look for the new cell in the parentheses in the laj list.
*-----------------------------------------------------------------------

  230    if(j1.gt.mlaj) goto 250 !FURUTA
         if(laj(j1+1).ne.1000001) goto 270

      do 240 j4 = j1+2,j3-1
         if(laj(j4).eq.1000002) goto 270
         iap = abs(laj(j4))
         if(jun(iap).ne.ju) goto 240

         call chkcll(iap,0,j)

         if(j.eq.0) goto 370
  240 continue

*-----------------------------------------------------------------------
*        search the overflow region at the end of laj.
*-----------------------------------------------------------------------

  250 do 260 i = mlaj+1,nlaj
         iap = laj(i)
         if(jun(iap).ne.ju) goto 260

         call chkcll(iap,0,j)

         if(j.eq.0) goto 370
  260 continue

*-----------------------------------------------------------------------
*        add a pair of parentheses if necessary.
*-----------------------------------------------------------------------

  270    if(laj(j1+1).eq.1000001.and.j1.ne.j3) goto 300
         if(nlaj+1.ge.mlaj) goto 300

      do 280 iz = j1+1,nlaj !FURUTA
         i = j1+1+nlaj-iz   !FURUTA
  280    laj(i+2) = laj(i)

         nlaj = nlaj+2
         laj(j1+1) = 1000001
         laj(j1+2) = 1000002
         j3 = j3+2
         j4 = j1+2

      do 290 j = i1,nlja
  290    lcaj(j+1) = sign(abs(lcaj(j+1))+2,lcaj(+j+1)) !FURUTA

*-----------------------------------------------------------------------
*        look for the new cell among all cells of the universe.
*-----------------------------------------------------------------------

  300 do 310 iap = 1,mxa
         if(jun(iap).ne.ju) goto 310
         call chkcll(iap,0,j)
         if(j.eq.0) goto 320
  310 continue

         iap = 1
         kdb = 1
         return

*-----------------------------------------------------------------------
*        insert the new cell in laj and increment the lcaj pointers.
*-----------------------------------------------------------------------

  320    if(nlaj+1.ge.mlaj) goto 350

      do 330 iz = j4,nlaj !FURUTA
         i = j4+nlaj-iz   !FURUTA
  330    laj(i+1) = laj(i)

         nlaj = nlaj+1
         laj(j4) = iap

      do 340 j = i1,nlja
  340    lcaj(j+1) = sign(abs(lcaj(j+1))+1,lcaj(j+1)) !FURUTA

         j3 = j3+1
         goto 370

*-----------------------------------------------------------------------
*        add the new cell to the overflow region at the end of laj.
*-----------------------------------------------------------------------

  350    j = 0
         nlaj = max(nlaj,mlaj)+1
         laj(nlaj) = iap !FURUTA

*-----------------------------------------------------------------------
*        find the element if the cell is a lattice.
*-----------------------------------------------------------------------

  370    if(lat(1,iap).ne.0) then
          if(abs(lat(1,iap)).ne.3)then
           call fndlat(io,iap,ierr)
          else
           do itet=1,nlat3
            if(itetcl(itet).eq.iap)then
             kkk=10000+itet
             exit
            endif
           enddo
           call tetrabox(io,kkk,xx0,ierr)
           call tetrafnd(xx0,xxx,yyy,zzz,uuu,vvv,www,coincd,itet,
     &          jjj,iap,ierr)
          endif
          if( ierr .ne. 0 ) then
           ierr = 2
           return
          end if
         end if
         if(mfl(1,iap).ne.0) goto 215

*-----------------------------------------------------------------------
*        find coincident planar surface at lowest level, if any.
*-----------------------------------------------------------------------

  372    if(jsu.ge.0) return
         jsu = -jsu

      do 375 i = abs(lca(iap)),abs(lca(iap+1))-1
         j = abs(lja(i))
         if(j.gt.1000000) goto 375
         if(ksc(jsu).ne.ksc(j)) goto 375
         k = kst(j)
         l = lsc(j)

         select case(k)                                        !FURUTA
         case(1)                                               !FURUTA
           a = scf(l+1)*xxx+scf(l+2)*yyy+scf(l+3)*zzz-scf(l+4) !FURUTA
           b = scf(l+1)*uuu+scf(l+2)*vvv+scf(l+3)*www          !FURUTA
         case(2)                                               !FURUTA
           a = xxx-scf(l+1)                                    !FURUTA
           b = uuu                                             !FURUTA
         case(3)                                               !FURUTA
           a = yyy-scf(l+1)                                    !FURUTA
           b = vvv                                             !FURUTA
         case(4)                                               !FURUTA
           a = zzz-scf(l+1)                                    !FURUTA
           b = www                                             !FURUTA
         case default                                          !FURUTA
           write(*,*)'ERROR in fcell: SELECT'                  !FURUTA
           stop                                                !FURUTA
         end select                                            !FURUTA

         if(abs(a).gt.coincd*abs(b)) goto 375
         if(b*lja(i).gt.0.) jsu = j
  375 continue

         return

*-----------------------------------------------------------------------
*        track is exiting a finite lattice.  bad trouble if lev > 0.
*        real particles escape.  detector pseudoparticles get lost.
*-----------------------------------------------------------------------

  380    kdb = 1


         iap = -12345

         if(lev.ne.0) then
            ierr = 2
            return
         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine inttor(ih,i,j,d1)
*                                                                      *
*       calculate the intersections of the track with a torus.         *
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      parameter ( eps5 = 1.0d-08 )

*-----------------------------------------------------------------------

      dimension uu(5),xx(5)

*-----------------------------------------------------------------------

         ix = kst(j)-23

         uu(1) = uuu
         uu(4) = uuu
         uu(2) = vvv
         uu(5) = vvv
         uu(3) = www

         xx(1) = xxx-scf(i)
         xx(4) = xx(1)
         xx(2) = yyy-scf(i+1)
         xx(5) = xx(2)
         xx(3) = zzz-scf(i+2)

         tpp(1) = uu(ix)**2
         tpp(2) = 1.-tpp(1)
         tpp(3) = 1./(tpp(2)+tpp(1)*scf(i+6))
         tpp(4) = (2.*scf(i+3)*tpp(3))**2
         tpp(5) = 2.*(uu(ix+2)*xx(ix+2)+uu(ix+1)*xx(ix+1))
         tpp(6) = xx(ix+2)**2+xx(ix+1)**2

         t1 = (tpp(5)+2.*uu(ix)*xx(ix)*scf(i+6))*tpp(3)
         t2 = (tpp(6)+scf(i+6)*xx(ix)**2
     &        + scf(i+3)**2-scf(i+5)**2)*tpp(3)

*-----------------------------------------------------------------------
*        set up the quartic coefficients.  (uu(i),i=1,5)=1,b,c,d,e
*-----------------------------------------------------------------------

         uu(1) = 1.
         uu(2) = 2.*t1
         uu(3) = t1**2+2.*t2-tpp(4)*tpp(2)
         uu(4) = uu(2)*t2-tpp(4)*tpp(5)
         uu(5) = t2**2-tpp(4)*tpp(6)

*-----------------------------------------------------------------------
*        eliminate the smallest root if the particle is sitting on
*        the torus.  solve cubic instead of quartic.
*-----------------------------------------------------------------------

         n = 4
         if(jsu.eq.j) n = 3

*-----------------------------------------------------------------------
*        there are no positive roots if no coefficients are negative.
*-----------------------------------------------------------------------

      do 10 m = 2, n
   10    if(uu(m).lt.0.) goto 20

         if(uu(n+1).ge.0.) return

*-----------------------------------------------------------------------
*        solve the quartic or cubic equation for the intersections.
*-----------------------------------------------------------------------

   20    call quart(n,uu,xx,jj,0)
         if(jj.eq.0) return

*-----------------------------------------------------------------------
*        zero in on legitimate positive roots.
*-----------------------------------------------------------------------

      do 60 m = 1, jj

*-----------------------------------------------------------------------
*        eliminate negative roots.
*-----------------------------------------------------------------------

         if(xx(m).lt.0.) goto 60

*-----------------------------------------------------------------------
*        eliminate roots on the wrong branch of a degenerate torus.
*-----------------------------------------------------------------------

         if((t2+xx(m)*(t1+xx(m)))*scf(i+3).lt.0.) goto 60
         z0 = eps5*(1.+abs(xx(m)))

      do 30 nr = 1, 20
         t3 = uu(4)+xx(m)*(2.*uu(3)+xx(m)*(3.*uu(2)+4.*xx(m)))
         if(t3.eq.0.) goto 60
         t3 = (uu(5)+xx(m)*(uu(4)+xx(m)*(uu(3)
     &         + xx(m)*(uu(2)+xx(m)))))/t3
         xx(m) = xx(m) - t3
   30    if(abs(t3).le.z0) goto 40

         goto 60

*-----------------------------------------------------------------------
*        return with first (smallest) legitimate root if cell is simple.
*-----------------------------------------------------------------------

   40    if(lca(ih).lt.0) goto 50
         d1 = xx(m)
         return

*-----------------------------------------------------------------------
*        put the roots in dti if the cell is complicated.
*-----------------------------------------------------------------------

   50    nlt = nlt+1
         dti(nlt) = xx(m)
         iti(nlt) = j

   60 continue

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function angl(ic)
*                                                                      *
*       calculate the cosine of the angle between the direction of the *
*       particle and the positive normal, ang, to the surface jsu.    *
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************
      use TETRAMOD, only: tetraangl
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

      integer,intent(in) :: ic !FURUTA20230706
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

      if( jsu .eq. 0 ) then
         angl = 1.0
         return
      end if

*-----------------------------------------------------------------------
*     Angle to tetrahedron surfaces (FURUTA20150714 TETRA)
*-----------------------------------------------------------------------
      if(junf.ne.0)then
       if(abs(lat(1,ic)).eq.3)then !FURUTA20230706
        call tetraangl(jsu,iii,uuu,vvv,www,cs)
        angl=cs
        return
       end if
      endif

*-----------------------------------------------------------------------

         k = kst(jsu)
         i = lsc(jsu)+1
         ang(1) = 0.
         ang(2) = 0.
         ang(3) = 0.

         if ( k.ge.24 .and. k.le.26 .and. jtr(jsu).ne.0 ) then
           if( jtr(jsu).gt.0 ) then
             do jt = 1, mxtr
               if( abs(trf(1,jt)).eq.jtr(jsu) ) then
                 m1 = jt ! frtati 2022/10/19 bugfix
                 exit
               end if
             end do
           else
             m1 = -jtr(jsu)/10000
             m2 = -jtr(jsu) - m1*10000
           end if
           xxx_cp = xxx
           yyy_cp = yyy
           zzz_cp = zzz
           uuu_cp = uuu
           vvv_cp = vvv
           www_cp = www
           call trnsxx(xxx_cp,yyy_cp,zzz_cp,xxx,yyy,zzz,m1)
           call trnsuu(uuu_cp,vvv_cp,www_cp,uuu,vvv,www,m1)
           if( jtr(jsu).lt.0 ) then
             xxx_it = xxx
             yyy_it = yyy
             zzz_it = zzz
             uuu_it = uuu
             vvv_it = vvv
             www_it = www
             call trnsxx(xxx_it,yyy_it,zzz_it,xxx,yyy,zzz,m2)
             call trnsuu(uuu_it,vvv_it,www_it,uuu,vvv,www,m2)
           end if
         end if

         goto (10,20,20,20,30,40,50,50,50,60,70,80,90,100,110,120,
     &         130,140,150,160,170,180,190,200,210,220) k

*-----------------------------------------------------------------------
* >>>>>  planes:  p, px, py, pz
*-----------------------------------------------------------------------

   10    ang(1) = scf(i)
         ang(2) = scf(i+1)
         ang(3) = scf(i+2)
         goto 230

   20    ang(k-1) = 1.
         select case(k)                       !FURUTA
         case(2)                              !FURUTA
           angl = uuu                         !FURUTA
         case(3)                              !FURUTA
           angl = vvv                         !FURUTA
         case(4)                              !FURUTA
           angl = www                         !FURUTA
         case default                         !FURUTA
           write(0,*)'ERROR in ang(): SELECT' !FURUTA
           stop                               !FURUTA
         end select                           !FURUTA
         return

*-----------------------------------------------------------------------
* >>>>>  spheres:  so, s, sx, sy, sz
*-----------------------------------------------------------------------

   30    ang(1) = xxx
         ang(2) = yyy
         ang(3) = zzz
         goto 230

   40    ang(1) = xxx-scf(i)
         ang(2) = yyy-scf(i+1)
         ang(3) = zzz-scf(i+2)
         goto 230

   50    ang(1) = xxx
         ang(2) = yyy
         ang(3) = zzz
         ang(k-6) = ang(k-6)-scf(i)
         goto 230

*-----------------------------------------------------------------------
* >>>>>  cylinders:  c/x, c/y, c/z, cx, cy, cz
*-----------------------------------------------------------------------

   60    ang(2) = yyy-scf(i)
         ang(3) = zzz-scf(i+1)
         goto 230

   70    ang(1) = xxx-scf(i)
         ang(3) = zzz-scf(i+1)
         goto 230

   80    ang(1) = xxx-scf(i)
         ang(2) = yyy-scf(i+1)
         goto 230

   90    ang(2) = yyy
         ang(3) = zzz
         goto 230

  100    ang(1) = xxx
         ang(3) = zzz
         goto 230

  110    ang(1) = xxx
         ang(2) = yyy
         goto 230

*-----------------------------------------------------------------------
* >>>>>  cones:  k/x, k/y, k/z, kx, ky, kz
*-----------------------------------------------------------------------

  120    ang(1) = (scf(i)-xxx)*scf(i+3)
         ang(2) = yyy-scf(i+1)
         ang(3) = zzz-scf(i+2)
         goto 230

  130    ang(1) = xxx-scf(i)
         ang(2) = (scf(i+1)-yyy)*scf(i+3)
         ang(3) = zzz-scf(i+2)
         goto 230

  140    ang(1) = xxx-scf(i)
         ang(2) = yyy-scf(i+1)
         ang(3) = (scf(i+2)-zzz)*scf(i+3)
         goto 230

  150    ang(1) = (scf(i)-xxx)*scf(i+1)
         ang(2) = yyy
         ang(3) = zzz
         goto 230

  160    ang(1) = xxx
         ang(2) = (scf(i)-yyy)*scf(i+1)
         ang(3) = zzz
         goto 230

  170    ang(1) = xxx
         ang(2) = yyy
         ang(3) = (scf(i)-zzz)*scf(i+1)
         goto 230

*-----------------------------------------------------------------------
* >>>>>  special quadratic:  sq
*-----------------------------------------------------------------------

  180    ang(1) = scf(i)*(xxx-scf(i+7))+scf(i+3)
         ang(2) = scf(i+1)*(yyy-scf(i+8))+scf(i+4)
         ang(3) = scf(i+2)*(zzz-scf(i+9))+scf(i+5)
         goto 230

*-----------------------------------------------------------------------
* >>>>>  general quadratic:  gq
*-----------------------------------------------------------------------

  190    ang(1) = 2.*scf(i)*xxx+scf(i+3)*yyy+scf(i+5)*zzz+scf(i+6)
         ang(2) = 2.*scf(i+1)*yyy+scf(i+3)*xxx+scf(i+4)*zzz+scf(i+7)
         ang(3) = 2.*scf(i+2)*zzz+scf(i+4)*yyy+scf(i+5)*xxx+scf(i+8)
         goto 230

*-----------------------------------------------------------------------
* >>>>>  tori:  tx, ty, tz
*-----------------------------------------------------------------------

  200    t4 = 2.*(1.-scf(i+3)/sqrt((yyy-scf(i+1))**2+(zzz-scf(i+2))**2))
         ang(1) = 2.*scf(i+6)*(xxx-scf(i))
         ang(2) = (yyy-scf(i+1))*t4
         ang(3) = (zzz-scf(i+2))*t4
         goto 230

  210    t4 = 2.*(1.-scf(i+3)/sqrt((xxx-scf(i))**2+(zzz-scf(i+2))**2))
         ang(1) = (xxx-scf(i))*t4
         ang(2) = 2.*scf(i+6)*(yyy-scf(i+1))
         ang(3) = (zzz-scf(i+2))*t4
         goto 230

  220    t4 = 2.*(1.-scf(i+3)/sqrt((xxx-scf(i))**2+(yyy-scf(i+1))**2))
         ang(1) = (xxx-scf(i))*t4
         ang(2) = (yyy-scf(i+1))*t4
         ang(3) = 2.*scf(i+6)*(zzz-scf(i+2))

*-----------------------------------------------------------------------
*        normalize ang and calculate the cosine.
*-----------------------------------------------------------------------

  230    t4 = sqrt(ang(1)**2+ang(2)**2+ang(3)**2)
         if(t4.ne.0.) t4 = 1./t4
         ang(1) = ang(1)*t4
         ang(2) = ang(2)*t4
         ang(3) = ang(3)*t4
         angl = max(-one,min(one,ang(1)*uuu+ang(2)*vvv+ang(3)*www))
         if ( k.ge.24 .and. k.le.26 .and. jtr(jsu).ne.0 ) then
           xxx = xxx_cp
           yyy = yyy_cp
           zzz = zzz_cp
           uuu = uuu_cp
           vvv = vvv_cp
           www = www_cp
         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine deflat(io,ierr)
*                                                                      *
*       set up arrays for repeated structures / lattice universe map.  *
*       called from pass1 (to get nmaz and nmzu) and igeom.            *
*                                                                      *
*        each cell, ic, is filled with either a universe that is a     *
*        lattice or a universe that is a list of cells.                *
*                                                                      *
*        set up the following arrays:                                  *
*        nmzu = length of mazu array. (approximate in 1st pass)        *
*        nmzu = length of mazu array. (in 1st pass = max(mxa,n)        *
*               where n=max size of any finite lattice)                *
*        nmaz = length of maze array. ( = 0 in 1st pass)               *
*        mazp(1,ic) = i, index of cell ic in mazu(j) list.             *
*        mazp(2,ic) = universe address j of cell ic.                   *
*        mazp(3,ic) = address j of universe filling cell ic.           *
*        mazu(j-3) = i = universe name.                                *
*        mazu(j-2) = finite lattice cell filling universe i.           *
*        mazu(j-1) = total number of lowest level elements below u=i.  *
*        mazu(j)   = ne = number of cells/elements in universe i       *
*        mazu(j+k) = number of elements below kth cell/universe.       *
*                    later changed to offsets in maze array.           *
*        mazu(j+ne+k) = kth cell in universe i (repeated structures)   *
*        mazu(j+ne+k) = first cell of universe filling kth lat element *
*                                                                      *
*        set up all but pointers and subelements of mazp and mazu.     *
*                                                                      *
*        iu = current universe name.                                   *
*        lp = current universe address.                                *
*        ic = current cell.                                            *
*        kf = number of finished universes.                            *
*        ke = number of finished elements in this universe.            *
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
*       NOT USED ANYMORE confirmed by T.Furuta on 2020/11/27           *
************************************************************************

      use LAFDATAMOD !FURUTA20201127
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------
*     special initialization if called from pass1 (nmaz=0).
*-----------------------------------------------------------------------

         if(nmaz.ne.0) goto 20
         lmzu = llaf + mlaf
         nmzu = max(nmzu,mxa)
         lmaz = lmzu+3*max(nmzu,2000)

      do 10 i = lmzu+1, lmaz
   10    mazu(i) = 0

*-----------------------------------------------------------------------
*        loop over all cells and process all universes not already done.
*-----------------------------------------------------------------------

   20    kf = 0
         lp = 4

      do 140 ic = 1, mxa
         if(mazp(lmzp+2,ic).ne.0) goto 140
         iu = abs(jun(ljun+ic))

*-----------------------------------------------------------------------
*        get mazp(3,jc): cell jc is filled with this universe.
*-----------------------------------------------------------------------

         if(iu.eq.0) goto 40

      do 30 jc = 1, mxa
   30    if(mfl(lmfl+1,jc).eq.iu) mazp(lmzp+3,jc) = lp

*-----------------------------------------------------------------------
*        process universe iu.
*-----------------------------------------------------------------------

   40    mazu(lmzu+lp-3) = iu
         ke = 0
         if(mfl(lmfl+1,ic).ge.0) goto 70

*-----------------------------------------------------------------------
*        finite lattice.
*-----------------------------------------------------------------------

         mf = -mfl(lmfl+1,ic)
         mazu(lmzu+lp-2) = ic
         ne = max(1,laf(mf+1,2)*laf(mf+2,2)*laf(mf+3,2))
         mazp(lmzp+2,ic) = lp
         mazu(lmzu+lp) = ne

      do 60 ie = 1, ne
         ju = laf(mf+1,2+ie)
         if(ju.eq.0) goto 50
         if(ju.ne.iu) goto 60
         mazu(lmzu+lp+ne+ie) = ic
         mazu(lmzu+lp+ie) = 1
   50    ke = ke+1
   60 continue
         goto 100

*-----------------------------------------------------------------------
*        repeated structure, infinite lattice, or simple cell.
*-----------------------------------------------------------------------

   70    ne=0

      do 80 jc = ic, mxa
         if(abs(jun(ljun+jc)).ne.iu) goto 80
         ne = ne + 1
         mazu(lmzu+lp+ne) = jc
         mazp(lmzp+2,jc) = lp
         mazp(lmzp+1,jc) = ne
   80 continue
         mazu(lmzu+lp) = ne

*-----------------------------------------------------------------------
*        shift mazu(j+k) to mazu(j+ne+k) and fill mazu(j+k) in simple
*        (unfilled) case.
*-----------------------------------------------------------------------

      do 90 i = 1, ne
         jc = mazu(lmzu+lp+i)
         mazu(lmzu+lp+ne+i) = jc
         mazu(lmzu+lp+i) = 0
         if(mfl(lmfl+1,jc).ne.0) goto 90
         mazu(lmzu+lp+i) = 1
         ke = ke + 1
   90 continue
  100    if(ke.ne.ne) goto 120

*-----------------------------------------------------------------------
*        convert sizes to offsets for completely filled cells.
*-----------------------------------------------------------------------

         lz = 0
      do 110 i = 1, ne
         ls = mazu(lmzu+lp+i)
         mazu(lmzu+lp+i) = lz
  110    lz = lz+ls
         mazu(lmzu+lp-1) = lz
         kf = kf+1

  120    lp = lp+2*ne+4
         if(nmaz.gt.0) goto 140
         if(lmzu+lp+2*nmzu.lt.lmaz) goto 140
         n = 3*max(nmzu,(lmaz-lmzu)/10,2000)

      do 130 i = lmaz+1, lmaz+n
  130    mazu(i) = 0

         lmaz = lmaz+n

  140 continue

         if(nmaz.gt.0.and.lp-4.gt.nmzu) then

            write(io,'(''** Error: deflat nmzu dimension overflow'',
     &                 2i10)') nmzu, lp-4

            ierr = ierr + 1
            return

         end if

         if(nmaz.eq.0) nmzu = lp-4
         if(kf.eq.0.and.nmaz.gt.0) then

            write(io,'(''**Error: deflat, '',
     &                 ''no completely filled universe found.'')')

            ierr = ierr + 1
            return

         end if

*-----------------------------------------------------------------------
*        step through mazu array to determine how many elements or
*        cells fill each universe, mazu(j-1), and the number of
*        subelements, mazu(j+i) in each element i.
*
*        kl = number of elements finished this pass.
*-----------------------------------------------------------------------

  150    kl = 0
         lp = 4

*-----------------------------------------------------------------------
*        process next universe.
*-----------------------------------------------------------------------

  160    ne = mazu(lmzu+lp)
         if(mazu(lmzu+lp-1).ne.0) goto 300
         ke = 0
         ic = mazu(lmzu+lp-2)
         if(ic.eq.0) goto 250

*-----------------------------------------------------------------------
*        finite lattice universe.
*-----------------------------------------------------------------------

         mf = -mfl(lmfl+1,ic)
         jc = 1

      do 240 ie = 1, ne
         if(mazu(lmzu+lp+ie).ne.0) goto 230
         ju = laf(mf+1,2+ie)
         if(ju.eq.0) goto 230
         if(ju.ne.mazu(lmzu+lp-3)) goto 170
         mazu(lmzu+lp+ne+ie) = ic
         mazu(lmzu+lp+ie) = 1
         goto 220
  170    if(mazu(lmzu+lp+ne+ie).ne.0) goto 200
         if(abs(jun(ljun+jc)).eq.ju) goto 190

      do 180 jc = 1, mxa
  180    if(abs(jun(ljun+jc)).eq.ju) goto 190

         if(nmaz.gt.0) then

            write(io,'(''**Error: deflat, cannot find universe'',i7,
     &                 '' filling lattice cell'',i7)') ju, ncl(lncl+ic)

            ierr = ierr + 1
            return

         end if

         goto 230

  190    mazu(lmzu+lp+ne+ie) = jc
         goto 210

  200    jc = mazu(lmzu+lp+ne+ie)
  210    mp = mazp(lmzp+2,jc)
         if(mazu(lmzu+mp-1).eq.0) goto 240
         mazu(lmzu+lp+ie) = mazu(lmzu+mp-1)
  220    kl = kl+1
  230    ke = ke+1

  240    continue
         goto 280

*-----------------------------------------------------------------------
*        repeated structure or infinite lattice or simple universe.
*-----------------------------------------------------------------------

  250 do 270 ie = 1, ne
         if(mazu(lmzu+lp+ie).ne.0) goto 260
         jc = mazu(lmzu+lp+ne+ie)
         mp = mazp(lmzp+3,jc)
         if(mazu(lmzu+mp-1).eq.0) goto 270
         mazu(lmzu+lp+ie) = mazu(lmzu+mp-1)
         kl = kl+1
  260    ke = ke+1
  270 continue

*-----------------------------------------------------------------------
*        universe finished. set mazu(j-1) and convert subelement lengths
*        to cumulative.
*-----------------------------------------------------------------------

  280    if(ke.ne.ne) goto 300
         kf = kf+1
         lz = 0

      do 290 k = 1, ne
         ls = mazu(lmzu+lp+k)
         mazu(lmzu+lp+k) = lz
  290    lz = lz+ls
         mazu(lmzu+lp-1) = lz
         if(mazu(lmzu+lp-3).eq.0) goto 310

*-----------------------------------------------------------------------
*        go to next unresolved (unknown length) universe.
*-----------------------------------------------------------------------

  300    lp = lp+2*ne+4
         if(lp.lt.nmzu) goto 160
         if(kl.gt.0) goto 150

         if(nmaz.gt.0) then

            write(io,'(''**Error: improper geometry specification '',
     &                 ''detected in deflat.'')')

            ierr = ierr + 1
            return

         end if

*-----------------------------------------------------------------------
*        universe 0 finished.  thus all universes resolved.
*-----------------------------------------------------------------------

  310    if(lz.eq.0) goto 320
         if(nmaz.eq.0) goto 330


         if(lz.le.nmaz) return

            write(io,'(''**Warning: deflat nmaz dimension overflow.'',
     &                 2i10)') nmaz, lz

         lz = 0

  320    continue

  330    nmaz = lz

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function nsf(n)
*                                                                      *
*       return problem surface number of lnsf+j=n and facet number kfq *
*       if program surface j is a facet of a macrobody, then return    *
*       problem surface number of main macrobody and facet number.     *
*       if j is an ordinary surface, kfq is set to zero.               *
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

         nsf = nsfm(n)
         kfq = 0
         if(nsf.ne.0) return
         j = n
         ka = ksm(j)
         kfq = j-ka
         if( ka.eq.0 ) then
           nsf = 0
         else
           nsf = nsfm(ka)
         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine mbody1(io,ierr,ix,igkst,m1c,idsn)
*                                                                      *
*       set up  and check additional surfaces for macrobodies.         *
*       mbody1 call from oldcrd                                        *
*       mbody2 call from igeom                                         *
*       ix = surface index if called from oldcrd.                      *
*                                                                      *
*       a new array, ksm, is created with an entry for each surface.   *
*       if the surface is not part of a macrobody then ksm(lksm+mxj)=0 *
*       if the surface is a master surface of  the macrobody, then     *
*       ksm(lksm+mxj)=-macrobody surface type, that is incremented     *
*       by -1000 if it is infinite in one dimension.                   *
*       if the surface is one of the other surfaces of the macrobody,  *
*       then ksm(lksm+mxj) is set to the problem number of master      *
*       surface.                                                       *
*                                                                      *
*       Last modified by K.Niita on 2012/12/22                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

      parameter ( ibmt = 40 )
      parameter ( eps5 = 1.0d-08 )

*-----------------------------------------------------------------------

*     dimension sc(33)
*     dimension ia(3)
*     dimension mf(6,5)
      real(8),save:: sc(33)
      integer,save:: ia(3)
      integer,save:: mf(6,5)

*-----------------------------------------------------------------------

      character ksf(50)*3

      data ( ksf(i), i = 1, ibmt ) /
     &         'p  ','px ','py ','pz ','so ',
     &         's  ','sx ','sy ','sz ','c/x',
     &         'c/y','c/z','cx ','cy ','cz ',
     &         'k/x','k/y','k/z','kx ','ky ',
     &         'kz ','sq ','gq ','tx ','ty ',
     &         'tz ','x  ','y  ','z  ','box',
     &         'rpp','sph','rcc','rec','ell',
     &         'trc','wed','arb','rhp','hex'/

*-----------------------------------------------------------------------

   10    n = igkst

      do 20 k = 1, 27
         if(k.gt.n) sc(k) = 0.
   20    if(k.le.n) sc(k) = scf(ix+k)
         lsc(mxj) = ix
         iy = ix
         ksm(mxj) = -m1c
         kst(mxj) = m1c
         jn = 6
         if(m1c.eq.33.or.m1c.eq.34.or.m1c.eq.36) jn = 3
         if(m1c.eq.39) jn = 8
         if(m1c.eq.32.or.m1c.eq.35) jn = 0
         if(m1c.eq.37) jn = 5
         goto (40,150,190,240,830,910,960,1020,30,340) m1c-29

   30    continue

         write(io,'(''** Error in [surface] section,''/
     &   ''   macro-body = '',a3,'' is not yet implemented.''/
     &   ''   surface id ='',i7 )')
     &   ksf(m1c), idsn

         ierr = ierr + 1
         return

*-----------------------------------------------------------------------
*     set coefficients for each j surface of box
*-----------------------------------------------------------------------

   40    if(n.eq.9) jn = 4

      do 80 k = 1, jn/2
         x = sc(3+3*(k-1)+1)
         y = sc(3+3*(k-1)+2)
         z = sc(3+3*(k-1)+3)
         r = x**2+y**2+z**2

         if(r.eq.0.) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,'' : zero dimension.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

            ierr = ierr + 1
            return

         end if

         sc(21+k) = 0.
         if(r.eq.0.) goto 80
         r = sqrt(r)

         sc(13+3*(k-1)) = x/r
         sc(14+3*(k-1)) = y/r
         sc(15+3*(k-1)) = z/r
         sc(21+k) = r

         if(k.eq.1) goto 80
         if(k.eq.3) goto 60

*-----------------------------------------------------------------------
*        check for orthogonal vectors.
*-----------------------------------------------------------------------

         call mbodyo(1,13,16,0,0,ab,ac,sc,io,ierr0,m1c,idsm)

            if( ierr0 .gt. 0 ) then

              ierr = ierr + 1
              return

            end if

         goto 80

*-----------------------------------------------------------------------

   60    call mbodyo(2,13,16,25,19,ab,ac,sc,io,ierr0,m1c,idsm)

            if( ierr0 .gt. 0 ) then

              ierr = ierr + 1
              return

            end if

*-----------------------------------------------------------------------

   80    continue
         x = -1.

         if(jn.eq.4) then

            write(io,'(''** Warning in [surface] section,''/
     &      ''   macro-body = '',a3,'' : box is infinite in one'',
     &      '' dimension.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

         end if

*-----------------------------------------------------------------------
*        increment ksm if box is infinite in a dimension.
*-----------------------------------------------------------------------

         if(jn.eq.4) ksm(mxj) = -m1c-1000

      do 140 j = 1, jn
         jtr(mxj+j) = jtr(mxj)
         ksu(mxj+j) = ksu(mxj)
         ksm(mxj+j) = mxj
         nsfm(mxj+j) = 0
         lsc(mxj+j) = iy
         kst(mxj+j) = 1
         jj = (j+1)/2

         if(x.ge.0.) goto 90
         ij = 0
         if(sc(14+(jj-1)*3).eq.0..and.sc(15+(jj-1)*3).eq.0.) ij = 1
         if(sc(13+(jj-1)*3).eq.0..and.sc(15+(jj-1)*3).eq.0.) ij = 2
         if(sc(13+(jj-1)*3).eq.0..and.sc(14+(jj-1)*3).eq.0.) ij = 3

   90    call mbodyp(5,ii,ij,ik,iy,j,jj,jn,ia,x,h,sc)
  140    continue

         goto 780

*-----------------------------------------------------------------------
*        set coefficients for each j surface of rpp.
*-----------------------------------------------------------------------

  150    x = -1.
         ik = 0

*-----------------------------------------------------------------------
*        check for infinite in one dimension.
*-----------------------------------------------------------------------

         if(sc(1).eq.sc(2)) ik = 11
         if(sc(3).eq.sc(4)) ik = ik+12
         if(sc(5).eq.sc(6)) ik = ik+13

         if(sc(1).gt.sc(2).or.sc(3).gt.sc(4).or.
     &      sc(5).gt.sc(6).or.ik.gt.15) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,'' : invalid dimension for '',
     &      ''macrobody.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

            ierr = ierr + 1
            return

         end if

         if(ik.gt.10.and.ik.lt.14) jn = 4

         if(jn.eq.4) then

            write(io,'(''** Warning in [surface] section,''/
     &      ''   macro-body = '',a3,'' : rpp is infinite in one'',
     &      '' dimension.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

         end if

*-----------------------------------------------------------------------
*        increment ksm if rpp is infinite in a dimension.
*-----------------------------------------------------------------------

         if(jn.eq.4) ksm(mxj) = -m1c-1000

      do 180 j = 1,jn
         jtr(mxj+j) = jtr(mxj)
         ksu(mxj+j) = ksu(mxj)
         ksm(mxj+j) = mxj
         nsfm(mxj+j) = 0
         kst(mxj+j) = 1
         jj = (j+1)/2
         lsc(mxj+j) = iy
         ij=jj
         call mbodyp(6,ii,ij,ik,iy,j,jj,jn,ia,x,h,sc)
  180    continue
         goto 780

*-----------------------------------------------------------------------
*        set coefficients for surface of sph, set as normal surface.
*-----------------------------------------------------------------------

  190    if(sc(n).le.0.) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,'' : invalid dimension for '',
     &      ''macrobody.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

            ierr = ierr + 1
            return

         end if

         ksm(mxj) = 0
         if(n.eq.1) goto 200
         if(sc(1).ne.0..or.sc(2).ne.0..or.sc(3).ne.0.) goto 210

  200    kst(mxj) = 5
         r = sc(1)
         if(n.ne.1) r = sc(4)
         scf(iy+1) = r**2
         goto 230

  210    j = 0
         if(sc(2).eq.0..and.sc(3).eq.0.) j = 1
         if(sc(1).eq.0..and.sc(3).eq.0.) j = 2
         if(sc(1).eq.0..and.sc(2).eq.0.) j = 3
         kst(mxj) = 6+j
         if(j.eq.0) goto 220
         scf(iy+1) = sc(j)
         scf(iy+2) = sc(4)**2
         iy = iy+1
         goto 230

  220    scf(iy+1) = sc(1)
         scf(iy+2) = sc(2)
         scf(iy+3) = sc(3)
         scf(iy+4) = sc(4)**2
         iy = iy+3
  230    iy = iy+1
         goto 780

*-----------------------------------------------------------------------
*        set coefficients for each j surface of rcc.
*-----------------------------------------------------------------------

  240    x = -1.

         if((sc(4).eq.0..and.sc(5).eq.0..and.sc(6).eq.0.).or.
     &       sc(7).le.0.) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,'' : invalid dimension for '',
     &      ''macrobody.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

            ierr = ierr + 1
            return

         end if

         jj = 1
         ij = 0
         if(sc(5).eq.0..and.sc(6).eq.0.) ij = 1
         if(sc(4).eq.0..and.sc(6).eq.0.) ij = 2
         if(sc(4).eq.0..and.sc(5).eq.0.) ij = 3
         if(ij.ne.0) goto 250
         r = sqrt(sc(4)**2+sc(5)**2+sc(6)**2)
         sc(8) = sc(4)/r
         sc(9) = sc(5)/r
         sc(10) = sc(6)/r

  250    if(ij.eq.1) r = sc(2)**2+sc(3)**2
         if(ij.eq.2) r = sc(1)**2+sc(3)**2
         if(ij.eq.3) r = sc(1)**2+sc(2)**2
         if(ij.eq.0) r = sc(1)**2+sc(2)**2+sc(3)**2

      do 330 j = 1, jn
         jtr(mxj+j) = jtr(mxj)
         ksu(mxj+j) = ksu(mxj)
         ksm(mxj+j) = mxj
         nsfm(mxj+j) = 0
         kst(mxj+j) = 1
         lsc(mxj+j) = iy
         if(j.ne.1) goto 280

*-----------------------------------------------------------------------
*        set up cylinder as special surface or gq.
*-----------------------------------------------------------------------

         if(ij.eq.0) goto 270
         if(r.ne.0.) goto 260
         kst(mxj+j) = 12+ij
         scf(iy+1) = sc(7)**2
         iy = iy+1
         goto 330

  260    kst(mxj+j) = 9+ij
         scf(iy+1) = sc(1)
         scf(iy+2) = sc(3)
         if(ij.eq.1) scf(iy+1) = sc(2)
         if(ij.eq.3) scf(iy+2) = sc(2)
         scf(iy+3) = sc(7)**2
         iy = iy+3
         goto 330

  270    kst(mxj+j) = 23
         scf(iy+1) = 1.0-sc(8)**2
         scf(iy+2) = 1.0-sc(9)**2
         scf(iy+3) = 1.0-sc(10)**2
         scf(iy+4) = -2.0*sc(8)*sc(9)
         scf(iy+5) = -2.0*sc(9)*sc(10)
         scf(iy+6) = -2.0*sc(8)*sc(10)
         scf(iy+7) = -sc(2)*scf(iy+4)
     &               -sc(3)*scf(iy+6)-2.0*sc(1)*scf(iy+1)
         scf(iy+8) = -sc(1)*scf(iy+4)
     &               -sc(3)*scf(iy+5)-2.0*sc(2)*scf(iy+2)
         scf(iy+9) = -sc(1)*scf(iy+6)
     &               -sc(2)*scf(iy+5)-2.0*sc(3)*scf(iy+3)
         scf(iy+10) = sc(1)*sc(2)*scf(iy+4)
     &               +sc(2)*sc(3)*scf(iy+5)
     &               +sc(1)*sc(3)*scf(iy+6)
     &               +sc(1)**2*scf(iy+1)
     &               +sc(2)**2*scf(iy+2)+sc(3)**2*scf(iy+3)-sc(7)**2
         iy = iy+10
         goto 330

  280    call mbodyp(3,ii,ij,ik,iy,j,jj,jn,ia,x,h,sc)

  330 continue
         goto 780

*-----------------------------------------------------------------------
*        set coefficients for each j surface of rhp.
*
*        interpret defaults when data truncated.
*-----------------------------------------------------------------------

  340    x = -1.

         if((sc(4).eq.0..and.sc(5).eq.0..and.sc(6).eq.0.)) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,'' : invalid dimension for '',
     &      ''macrobody.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

            ierr = ierr + 1
            return

         end if

         if(n.gt.9.and.n.le.12) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,'' : no default allowed for '',
     &      ''third facet.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

            ierr = ierr + 1
            return

         end if

         ij = 0
         if(sc(5).eq.0..and.sc(6).eq.0.) ij = 1
         if(sc(4).eq.0..and.sc(6).eq.0.) ij = 2
         if(sc(4).eq.0..and.sc(5).eq.0.) ij = 3
         h = sqrt(sc(4)**2+sc(5)**2+sc(6)**2)
         sc(16) = sc(4)/h
         sc(17) = sc(5)/h
         sc(18) = sc(6)/h

*-----------------------------------------------------------------------
*        set facet for regular hexagon if default, and
*        check that facet is orthogonal to axis.
*-----------------------------------------------------------------------

         call mbodyr(ia,n,sc,io,ierr0,m1c,idsm)

         if(h.ge.1.e6) jn = 6

         if(jn.eq.6) then

            write(io,'(''** Warning in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : rhp is infinite in axial dimension.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

         end if

         if(jn.eq.6) ksm(mxj) = -m1c-1000

*-----------------------------------------------------------------------
*        create 6 or 8 surfaces.
*-----------------------------------------------------------------------

      do 770 j = 1, jn
         jtr(mxj+j) = jtr(mxj)
         ksu(mxj+j) = ksu(mxj)
         ksm(mxj+j) = mxj
         nsfm(mxj+j) = 0
         lsc(mxj+j) = iy
         kst(mxj+j) = 1
         jj = (j+1)/2
         if(j.gt.6) goto 720

         call mbodyp(7,ii,ij,ik,iy,j,jj,jn,ia,x,h,sc)

         goto 770

  720    call mbodyp(8,ii,ij,ik,iy,j,jj,jn,ia,x,h,sc)

  770 continue
         goto 780

*-----------------------------------------------------------------------
*     rec: set coefficients for each j surface.
*-----------------------------------------------------------------------

  830 x=-1.
      jj=1

      if((sc(4).eq.0..and.sc(5).eq.0..and.sc(6).eq.0.).or.
     & (sc(7).eq.0..and.sc(8).eq.0..and.sc(9).eq.0.).or.
     & (sc(10).eq.0..and.sc(11).eq.0..and.sc(12).eq.0.)) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid dimension for rec.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      if(n.eq.12.and.sc(10).eq.0..and.sc(11).eq.0..and.sc(12).eq.0.)
     & then
            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid dimension for rec.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

       if(n.eq.10.and.sc(10).le.0.) then
            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid dimension for rec.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

       end if

      if(n.eq.10) then
            write(io,'(''** Warning in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : alternate format used for rec macrobody.'')')
     &      ksf(m1c)
      end if


*-----------------------------------------------------------------------
*     determine square of major and minor radii, a2 and b2.
*-----------------------------------------------------------------------

      a2=sc(7)**2+sc(8)**2+sc(9)**2
      if(n.eq.12)b2=sc(10)**2+sc(11)**2+sc(12)**2
      if(n.eq.10)b2=sc(10)**2
      h2=sc(4)**2+sc(5)**2+sc(6)**2
      a=sqrt(a2)
      b=sqrt(b2)
      h=sqrt(h2)

*-----------------------------------------------------------------------
*     obtain r2 vector as cross product in alternate format.
*-----------------------------------------------------------------------

      if(n.eq.10)call crspro(sc(4),sc(7),sc(10))

      if(n.eq.10)bq=sqrt(sc(10)**2+sc(11)**2+sc(12)**2)

*-----------------------------------------------------------------------
*     first normalize.
*-----------------------------------------------------------------------

      do 840 kk=1,3
      sc(3+kk)=sc(3+kk)/h
      sc(6+kk)=sc(6+kk)/a
      if(n.eq.12)sc(9+kk)=sc(9+kk)/b
  840 if(n.eq.10)sc(9+kk)=sc(9+kk)/bq

*-----------------------------------------------------------------------
*     check for orthogonal vectors.
*-----------------------------------------------------------------------

      call mbodyo(1,4,7,0,0,ab,ac,sc,io,ierr0,m1c,idsm)

      if(abs(ab).gt..01)go to 860
      call mbodyo(2,4,7,25,10,ab,ac,sc,io,ierr0,m1c,idsm)

  860 af=sqrt(sc(7)**2+sc(8)**2+sc(9)**2)
      bf=sqrt(sc(10)**2+sc(11)**2+sc(12)**2)
      hf=sqrt(sc(4)**2+sc(5)**2+sc(6)**2)

*-----------------------------------------------------------------------
*     now restore lengths.
*-----------------------------------------------------------------------

      do 870 kk=1,3
      sc(3+kk)=sc(3+kk)*h/hf
      sc(6+kk)=sc(6+kk)*a/af
  870 sc(9+kk)=sc(9+kk)*b/bf
      ij=0
      if(sc(5).eq.0..and.sc(6).eq.0.)ij=1
      if(sc(4).eq.0..and.sc(6).eq.0.)ij=2
      if(sc(4).eq.0..and.sc(5).eq.0.)ij=3
      if(ij.eq.0)r=sc(1)**2+sc(2)**2+sc(3)**2
      if(ij.eq.1)r=sc(2)**2+sc(3)**2
      if(ij.eq.2)r=sc(1)**2+sc(3)**2
      if(ij.eq.3)r=sc(1)**2+sc(2)**2
      do 880 j=1,jn
      jtr(mxj+j)=jtr(mxj)
      ksu(mxj+j)=ksu(mxj)
      ksm(mxj+j)=mxj
      nsfm(mxj+j)=0
      kst(mxj+j)=1
      lsc(mxj+j)=iy
      if(j.ne.1)go to 890

*-----------------------------------------------------------------------
*     set up ellipse as special surface or gq.
*     sq must also have major axis along axis.
*-----------------------------------------------------------------------

      if(ij.eq.0)go to 885
      if(sc(7).ne.0..and.sc(8).ne.0.)go to 885
      if(sc(7).ne.0..and.sc(9).ne.0.)go to 885
      if(sc(8).ne.0..and.sc(9).ne.0.)go to 885
      kst(mxj+j)=22
      scf(iy+1)=a2*(1.-sc(4)**2/h2)-sc(7)**2*(1.-(b2/a2))
      scf(iy+2)=a2*(1.-sc(5)**2/h2)-sc(8)**2*(1.-(b2/a2))
      scf(iy+3)=a2*(1.-sc(6)**2/h2)-sc(9)**2*(1.-(b2/a2))
      if(ij.eq.1)scf(iy+1)=0.
      if(ij.eq.2)scf(iy+2)=0.
      if(ij.eq.3)scf(iy+3)=0.
      scf(iy+4)=0.
      scf(iy+5)=0.
      scf(iy+6)=0.
      if(ij.eq.1)scf(iy+7)=a2*r-(sc(2)**2*sc(8)**2+sc(3)**2*sc(9)**2)*
     & (1.-(b2/a2))-a2*b2-scf(iy+2)*sc(2)**2-scf(iy+3)*sc(3)**2
      if(ij.eq.2)scf(iy+7)=a2*r-(sc(1)**2*sc(7)**2+sc(3)**2*sc(9)**2)*
     & (1.-(b2/a2))-a2*b2-scf(iy+1)*sc(1)**2-scf(iy+3)*sc(3)**2
      if(ij.eq.3)scf(iy+7)=a2*r-(sc(1)**2*sc(7)**2+sc(2)**2*sc(8)**2)*
     & (1.-(b2/a2))-a2*b2-scf(iy+1)*sc(1)**2-scf(iy+2)*sc(2)**2
      scf(iy+8)=sc(1)
      scf(iy+9)=sc(2)
      scf(iy+10)=sc(3)
      iy=iy+10
      go to 880

*-----------------------------------------------------------------------
*     gq surface.
*-----------------------------------------------------------------------

  885 kst(mxj+j)=23
      scf(iy+1)=a2*(1.-sc(4)**2/h2)-sc(7)**2*(1.-(b2/a2))
      scf(iy+2)=a2*(1.-sc(5)**2/h2)-sc(8)**2*(1.-(b2/a2))
      scf(iy+3)=a2*(1.-sc(6)**2/h2)-sc(9)**2*(1.-(b2/a2))
      scf(iy+4)=-2.*a2*sc(4)*sc(5)/h2 -2.*sc(7)*sc(8)*(1.-(b2/a2))
      scf(iy+5)=-2.*a2*sc(5)*sc(6)/h2 -2.*sc(8)*sc(9)*(1.-(b2/a2))
      scf(iy+6)=-2.*a2*sc(4)*sc(6)/h2 -2.*sc(7)*sc(9)*(1.-(b2/a2))
      scf(iy+7)=a2*2.*((sc(4)/h2)*(sc(2)*sc(5)+sc(3)*sc(6))-
     & sc(1)*((1.-(sc(4)**2/h2))-((sc(7)**2/a2)*(1.-(b2/a2)))))
     & +2.*(1.-(b2/a2))*(sc(2)*sc(7)*sc(8)+sc(3)*sc(7)*sc(9))
      scf(iy+8)=a2*2.*((sc(5)/h2)*(sc(1)*sc(4)+sc(3)*sc(6)) -
     & sc(2)*((1.-(sc(5)**2/h2))-((sc(8)**2/a2)*(1.-(b2/a2)))))
     & +2.*(1.-(b2/a2))*(sc(1)*sc(7)*sc(8)+sc(3)*sc(8)*sc(9))
      scf(iy+9)=a2*2.*((sc(6)/h2)*(sc(1)*sc(4)+sc(2)*sc(5)) -
     & sc(3)*((1.-(sc(6)**2/h2))-((sc(9)**2/a2)*(1.-(b2/a2)))))
     & +2.*(1.-(b2/a2))*(sc(1)*sc(7)*sc(9)+sc(2)*sc(8)*sc(9))
      scf(iy+10)=a2*((-2./h2)*(sc(1)*sc(2)*sc(4)*sc(5)+sc(2)*sc(3)*
     & sc(5)*sc(6)+sc(1)*sc(3)*sc(4)*sc(6))+sc(1)**2*(1.-sc(4)**2/h2)+
     & sc(2)**2*(1.-sc(5)**2/h2)+sc(3)**2*(1.-sc(6)**2/h2))-
     & (sc(1)**2*sc(7)**2+sc(2)**2*sc(8)**2+sc(3)**2*sc(9)**2)*
     & (1.-b2/a2)-a2*b2
     & -2.*(sc(1)*sc(2)*sc(7)*sc(8)+sc(2)*sc(3)*sc(8)*sc(9)+
     & sc(1)*sc(3)*sc(7)*sc(9))*(1.-b2/a2)
      iy=iy+10
      go to 880

  890 call mbodyp(1,ii,ij,ik,iy,j,jj,jn,ia,x,h,sc)

  880 continue
      goto 780


*-----------------------------------------------------------------------
*     ell: set coefficients for surface as sq or gq surface.
*-----------------------------------------------------------------------

  910 if(sc(7).eq.0.) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid radius for ell.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      if(sc(7).ge.0.) goto 920

*-----------------------------------------------------------------------
*     using alternate format with center of ellipse, major axis
*     radius vector, and minor radius (as negative).
*-----------------------------------------------------------------------

      b2=sc(7)**2
      a2=sc(4)**2+sc(5)**2+sc(6)**2

      if(a2.eq.0.) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid radius for ell.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      a=sqrt(a2)
      u=sc(4)/a
      v=sc(5)/a
      w=sc(6)/a
      x=sc(1)
      y=sc(2)
      z=sc(3)

      goto 930

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

  920 a2=sc(7)**2
      c=(sc(4)-sc(1))**2+(sc(5)-sc(2))**2+(sc(6)-sc(3))**2

      if(c.le.0.) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid foci for ell.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      c=sqrt(c)
      u=(sc(4)-sc(1))/c
      v=(sc(5)-sc(2))/c
      w=(sc(6)-sc(3))/c
      c=c*0.5
      b2=a2-c**2

      if(b2.le.0.) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid foci for ell.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      x=.5*(sc(1)+sc(4))
      y=.5*(sc(2)+sc(5))
      z=.5*(sc(3)+sc(6))

*-----------------------------------------------------------------------
*     set up ellipse as sq or gq surface.
*-----------------------------------------------------------------------

  930 am=-a2+b2
      du=a2+u**2*am
      dv=a2+v**2*am
      dw=a2+w**2*am
      ksm(mxj)=0

      if((u.ne.0..and.v.ne.0.).or.(u.ne.0..and.w.ne.0.).or.
     & (v.ne.0..and.w.ne.0.)) goto 940

*-----------------------------------------------------------------------
*     sq.
*-----------------------------------------------------------------------

      kst(mxj)=22
      scf(iy+1)=du
      scf(iy+2)=dv
      scf(iy+3)=dw
      scf(iy+4)=0.
      scf(iy+5)=0.
      scf(iy+6)=0.
      scf(iy+7)=-a2*b2
      scf(iy+8)=x
      scf(iy+9)=y
      scf(iy+10)=z
      goto 950

  940 kst(mxj)=23

*-----------------------------------------------------------------------
*     gq.
*-----------------------------------------------------------------------

      scf(iy+1)=du
      scf(iy+2)=dv
      scf(iy+3)=dw
      scf(iy+4)=2.*u*v*am
      scf(iy+5)=2.*v*w*am
      scf(iy+6)=2.*u*w*am
      scf(iy+7)=-2.*x*du-2.*u*v*y*am-2.*u*w*z*am
      scf(iy+8)=-2.*y*dv-2.*u*v*x*am-2.*v*w*z*am
      scf(iy+9)=-2.*z*dw-2.*u*w*x*am-2.*v*w*y*am
      scf(iy+10)=x**2*du+y**2*dv+z**2*dw+x*y*2.*u*v*am+y*z*2.*v*w*am+
     & x*z*2.*u*w*am-a2*b2

  950 iy = iy + 10
      goto 780

*-----------------------------------------------------------------------
*     trc: set coefficients for each j surface.
*-----------------------------------------------------------------------

  960 x=-1.
      jj=1

      if((sc(4).eq.0..and.sc(5).eq.0..and.sc(6).eq.0.).or.
     & (sc(7).eq.0..and.sc(8).eq.0.).or.sc(7).lt.0..or.sc(8).lt.0.
     &  .or.sc(7).eq.sc(8))
     & then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid dimension for trc.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

       end if

      ij=0
      if(sc(5).eq.0..and.sc(6).eq.0.)ij=1
      if(sc(4).eq.0..and.sc(6).eq.0.)ij=2
      if(sc(4).eq.0..and.sc(5).eq.0.)ij=3
      rh=sqrt(sc(4)**2+sc(5)**2+sc(6)**2)
      sc(9)=sc(4)/rh
      sc(10)=sc(5)/rh
      sc(11)=sc(6)/rh
  970 if(ij.eq.1)r=sc(2)**2+sc(3)**2
      if(ij.eq.2)r=sc(1)**2+sc(3)**2
      if(ij.eq.3)r=sc(1)**2+sc(2)**2
      if(ij.eq.0)r=sc(1)**2+sc(2)**2+sc(3)**2

      do 1010 j=1,jn
      jtr(mxj+j)=jtr(mxj)
      ksu(mxj+j)=ksu(mxj)
      ksm(mxj+j)=mxj
      nsfm(mxj+j)=0
      kst(mxj+j)=1
      lsc(mxj+j)=iy
      if(j.ne.1) goto 1000

*-----------------------------------------------------------------------
*     set up cone as special surface or gq.
*-----------------------------------------------------------------------

      if(ij.eq.0) goto 990
      if(r.ne.0.) goto 980
      kst(mxj+j)=18+ij
      scf(iy+2)=(sc(7)-sc(8))**2/rh**2
      scf(iy+1)=sc(ij)+(sc(7)/sc(8))*sc(3+ij)/((sc(7)/sc(8))-1.)
      scf(iy+3)=0
      iy=iy+3
      goto 1010

  980 kst(mxj+j)=15+ij
      scf(iy+4)=(sc(7)-sc(8))**2/rh**2
      scf(iy+ij)=sc(ij)+(sc(7)/sc(8))*sc(3+ij)/((sc(7)/sc(8))-1.)
      scf(iy+5)=0
      if(ij.ne.1)scf(iy+1)=sc(1)
      if(ij.ne.2)scf(iy+2)=sc(2)
      if(ij.ne.3)scf(iy+3)=sc(3)
      iy=iy+5
      goto 1010

  990 kst(mxj+j)=23
      scf(iy+1)=1.-sc(9)**2
      scf(iy+2)=1.-sc(10)**2
      scf(iy+3)=1.-sc(11)**2
      scf(iy+4)=-2.*sc(9)*sc(10)
      scf(iy+5)=-2.*sc(10)*sc(11)
      scf(iy+6)=-2.*sc(9)*sc(11)
      scf(iy+7)=-sc(2)*scf(iy+4)-sc(3)*scf(iy+6)-2.*sc(1)*scf(iy+1)
      scf(iy+8)=-sc(1)*scf(iy+4)-sc(3)*scf(iy+5)-2.*sc(2)*scf(iy+2)
      scf(iy+9)=-sc(1)*scf(iy+6)-sc(2)*scf(iy+5)-2.*sc(3)*scf(iy+3)
      scf(iy+10)=sc(1)*sc(2)*scf(iy+4)+sc(2)*sc(3)*scf(iy+5)+
     & sc(1)*sc(3)*scf(iy+6)+sc(1)**2*scf(iy+1)+sc(2)**2*scf(iy+2)+
     & sc(3)**2*scf(iy+3)-sc(7)**2

*-----------------------------------------------------------------------
*     add terms to change cylinder to cone.
*-----------------------------------------------------------------------

      t2=(sc(7)-sc(8))**2/rh**2
      dp=-(sc(7)/sc(8))*rh/((sc(7)/sc(8))-1.)
      scf(iy+1)=scf(iy+1)-t2*sc(9)**2
      scf(iy+2)=scf(iy+2)-t2*sc(10)**2
      scf(iy+3)=scf(iy+3)-t2*sc(11)**2
      scf(iy+7)=scf(iy+7)-t2*
     &(sc(2)*scf(iy+4)+sc(3)*scf(iy+6)-2.*sc(1)*sc(9)**2+2.*dp*sc(9))
      scf(iy+8)=scf(iy+8)-t2*
     &(sc(1)*scf(iy+4)+sc(3)*scf(iy+5)-2.*sc(2)*sc(10)**2+2.*dp*sc(10))
      scf(iy+9)=scf(iy+9)-t2*
     &(sc(1)*scf(iy+6)+sc(2)*scf(iy+5)-2.*sc(3)*sc(11)**2+2.*dp*sc(11))
      scf(iy+10)=scf(iy+10)+sc(7)**2-t2*
     & (-sc(1)*sc(2)*scf(iy+4)-sc(2)*sc(3)*scf(iy+5)-
     & sc(1)*sc(3)*scf(iy+6)+sc(1)**2*sc(9)**2+sc(2)**2*sc(10)**2+
     & sc(3)**2*sc(11)**2+dp**2
     & -2.*dp*(sc(1)*sc(9)+sc(2)*sc(10)+sc(3)*sc(11)))
      scf(iy+4)=scf(iy+4)-2.*t2*sc(9)*sc(10)
      scf(iy+5)=scf(iy+5)-2.*t2*sc(10)*sc(11)
      scf(iy+6)=scf(iy+6)-2.*t2*sc(9)*sc(11)
      iy=iy+10

      goto 1010

 1000 call mbodyp(2,ii,ij,ik,iy,j,jj,jn,ia,x,h,sc)

 1010 continue

      goto 780

*-----------------------------------------------------------------------
*     wed: set coefficients for each j surface.
*-----------------------------------------------------------------------

 1020 if(n.ne.12) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid demension for wed.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

       end if

      do 1040 k=1,3
      x=sc(3+3*(k-1)+1)
      y=sc(3+3*(k-1)+2)
      z=sc(3+3*(k-1)+3)
      r=x**2+y**2+z**2

      if(r.eq.0.) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : zero demension for wed.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      sc(21+k)=0.
      if(r.eq.0.) goto 1040
      r=sqrt(r)
      sc(13+3*(k-1))=x/r
      sc(14+3*(k-1))=y/r
      sc(15+3*(k-1))=z/r
      sc(21+k)=r
      if(k.eq.1) goto 1040
      if(k.eq.3) goto 1030

*-----------------------------------------------------------------------
*     check for orthogonal vectors.
*-----------------------------------------------------------------------

      call mbodyo(1,13,16,0,0,ab,ac,sc,io,ierr0,m1c,idsm)
      go to 1040
 1030 call mbodyo(2,13,16,25,19,ab,ac,sc,io,ierr0,m1c,idsm)
 1040 continue

*-----------------------------------------------------------------------
*     find vector to slanted plane.
*-----------------------------------------------------------------------

      d=0
      do 1050 k=1,3
      sc(30+k)=sc(12+k)*sc(23)+sc(15+k)*sc(22)
 1050 d=d+sc(30+k)**2
      d=sqrt(d)
      do 1060 k=1,3
 1060 sc(30+k)=sc(30+k)/d
      d=sc(22)*sc(23)/sqrt(sc(22)**2+sc(23)**2)
      x=-1.
      do 1100 j=1,jn
      jtr(mxj+j)=jtr(mxj)
      ksu(mxj+j)=ksu(mxj)
      ksm(mxj+j)=mxj
      nsfm(mxj+j)=0
      lsc(mxj+j)=iy
      kst(mxj+j)=1
      jj=(j+1)/2
      if(j.gt.3)jj=(j+2)/2
      x=-x
      if(x.lt.0.) goto 1070
      ij=0
      if(sc(14+(jj-1)*3).eq.0..and.sc(15+(jj-1)*3).eq.0.)ij=1
      if(sc(13+(jj-1)*3).eq.0..and.sc(15+(jj-1)*3).eq.0.)ij=2
      if(sc(13+(jj-1)*3).eq.0..and.sc(14+(jj-1)*3).eq.0.)ij=3
 1070 if(j.eq.1) goto 1080
      call mbodyp(4,ii,ij,ik,iy,j,jj,jn,ia,x,h,sc)
      goto 1100

*-----------------------------------------------------------------------
*     do slanted plane.
*-----------------------------------------------------------------------

 1080 do 1090 k=1,4
       if(k.lt.4)then
        scf(iy+k)=sc(30+k)
       else
        scf(iy+k)=scf(iy+1)*sc(1)+scf(iy+2)*sc(2)+scf(iy+3)*sc(3)
        scf(iy+k)=scf(iy+k)+(scf(iy+1)*sc(31)+
     &       scf(iy+2)*sc(32)+scf(iy+3)*sc(33))*d
       endif
 1090 continue
      iy=iy+4
 1100 continue
      goto 780

*-----------------------------------------------------------------------
*     arb: set coefficients for each j surface.
*-----------------------------------------------------------------------

 1110 if(n.lt.28.or.n.gt.30) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid number entries for arb.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      nc=0
      im=0
      do 1120 ii=25,n
      i=sc(ii)

      if(i.lt.0.or.(i.eq.0.and.ii.eq.25)) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid four digit face entries for arb.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      if(i.le.0) goto 1120
      if(i.lt.1230.or.i.gt.8765) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid four digit face entries for arb.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      if(i.lt.1230.or.i.gt.8765) goto 1120
      nc=nc+1
      i1=i/1000

      if(i1.gt.8) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid four digit face entries for arb.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      mf(nc,1)=i1
      if(i1.gt.im)im=i1
      i2=i-i1*1000
      i2=i2/100

      if(i2.gt.8.or.i2.lt.1.or.i2.eq.i1) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid four digit face entries for arb.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      mf(nc,2)=i2
      if(i2.gt.im)im=i2
      i3=i-i1*1000-i2*100
      i3=i3/10

      if(i3.gt.8.or.i3.lt.1.or.i3.eq.i2.or.i3.eq.i1) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid four digit face entries for arb.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      mf(nc,3)=i3
      mf(nc,4)=0
      mf(nc,5)=3
      if(i3.gt.im)im=i3
      i4=i-i1*1000-i2*100-i3*10
      if(i4.eq.0)go to 1120

      if(i4.gt.8.or.i4.eq.i3.or.i4.eq.i2.or.i4.eq.i1) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : invalid four digit face entries for arb.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      mf(nc,4)=i4
      if(i4.gt.im)im=i4
      mf(nc,5)=4
 1120 continue

      if(nc.ge.4.and.nc.le.6) then

            write(io,'(''** Warning in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' has specified number faces.'',i3)')
     &      ksf(m1c), nc

      end if

      if(nc.lt.4.or.nc.gt.6) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' has specified number faces.'',i3)')
     &      ksf(m1c), nc

            ierr = ierr + 1
            return

      end if

      jn=nc

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

      nt=0
      do 1150 i=1,nc
      nt=nt+mf(i,5)
      do 1140 kk=1,4
      if(kk.gt.mf(i,5)) goto 1140
      j=mf(i,kk)
      do 1130 k=1,3
      if(i.eq.1.and.kk.eq.1)tpp(13+k)=0.
 1130 tpp(13+k)=tpp(13+k)+sc(k+(j-1)*3)
 1140 continue
 1150 continue
      do 1160 k=1,3
 1160 tpp(13+k)=tpp(13+k)/nt

*-----------------------------------------------------------------------
*    increment ksm by 1000*number faces.
*-----------------------------------------------------------------------

      ksm(mxj)=ksm(mxj)-1000*jn

*-----------------------------------------------------------------------
*     implement each plane surface of arb.
*-----------------------------------------------------------------------

      ie=0
      do 1260 j=1,jn
      jtr(mxj+j)=jtr(mxj)
      ksu(mxj+j)=ksu(mxj)
      ksm(mxj+j)=mxj
      nsfm(mxj+j)=0
      lsc(mxj+j)=iy
      kst(mxj+j)=1

*-----------------------------------------------------------------------
*     prepare to do cross product to obtain surface normal
*-----------------------------------------------------------------------

      j1=mf(j,1)
      j2=mf(j,2)
      j3=mf(j,3)
      j4=mf(j,4)

      do 1170 k=1,3
      tpp(k)=sc(k+(j2-1)*3)-sc(k+(j1-1)*3)
 1170 tpp(k+3)=sc(k+(j3-1)*3)-sc(k+(j2-1)*3)

      call crspro(tpp(1),tpp(4),tpp(10))
      s=0.
      do 1180 k=1,3
 1180 s=s+tpp(9+k)**2

      if(s.eq.0.) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' entries for arb not consistent.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

      end if

      s=sqrt(s)

*-----------------------------------------------------------------------
*     a,b,c coefficients are the normal to plane.
*     obtain d coefficient from plane eq. at first point.
*-----------------------------------------------------------------------

      scf(iy+4)=0.
      do 1190 k=1,3
      tpp(6+k)=tpp(9+k)/s
      scf(iy+k)=tpp(6+k)
 1190 scf(iy+4)=scf(iy+4)+scf(iy+k)*sc(k+(j1-1)*3)

*-----------------------------------------------------------------------
*     check sense and switch signs if necessary.
*-----------------------------------------------------------------------

      ss=-scf(iy+4)
      do 1200 k=1,3
 1200 ss=ss+scf(iy+k)*tpp(13+k)
      if(ss.lt.0.) goto 1220
      do 1210 k=1,4
 1210 scf(iy+k)=-scf(iy+k)

*-----------------------------------------------------------------------
*     if 4th point, determine how far from plane.
*-----------------------------------------------------------------------

 1220 if(j4.eq.0) goto 1260
      ds=scf(iy+4)
      do 1230 k=1,3
 1230 ds=ds-scf(iy+k)*sc(k+(j4-1)*3)
      if(abs(ds).lt.1.e-6) goto 1260
      ie=ie+1

      if(ie.eq.1) then

            write(io,'(''** Warning in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      ''fourth point more than 1.e-6 from plane for arb.'')')
     &      ksf(m1c)

      end if


      x=sc(1+(j4-1)*3)+scf(iy+1)*ds
      y=sc(2+(j4-1)*3)+scf(iy+2)*ds
      z=sc(3+(j4-1)*3)+scf(iy+3)*ds


 1260 iy=iy+4

*-----------------------------------------------------------------------
*     check for each surface of arb that mid-point of surface
*     based on 3 corners has negative sense with respect to
*     each of the other surfaces.
*-----------------------------------------------------------------------

      do 1290 j=1,jn
      x=0
      y=0
      z=0
      do 1270 k=1,3
      jk=mf(j,k)
      x=x+third*sc(1+(jk-1)*3)
      y=y+third*sc(2+(jk-1)*3)
 1270 z=z+third*sc(3+(jk-1)*3)
      do 1280 jj=1,jn
      if(jj.eq.j) goto 1280
      iq=iy-(jn-jj+1)*4
      ss=-scf(iq+4)
      ss=ss+scf(iq+1)*x
      ss=ss+scf(iq+2)*y
      ss=ss+scf(iq+3)*z
      if(ss.lt.0.) goto 1280

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      ''midpoint of arb surface,i3,20h has positive sense.'')')
     &      ksf(m1c)

            ierr = ierr + 1
            return

 1280 continue
 1290 continue
      goto 780

*-----------------------------------------------------------------------
*        finish macrocell.
*-----------------------------------------------------------------------

  780    if(jn.ne.0) jtr(mxj) = 0
         mxj = mxj + jn
         lsc(mxj+1) = iy

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine mbody2(io,ierr)
*                                                                      *
*       set up  and check additional surfaces for macrobodies.         *
*       mbody1 call from oldcrd                                        *
*       mbody2 call from igeom                                         *
*       ix = surface index if called from oldcrd.                      *
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

  800 do 890 ic = 1, mxa
         m1 = abs(lca(ic))
  810    m2 = abs(lca(ic+1))-1

      do 830 la = m1, m2
         if(lja(la).gt.1000000) goto 830
         if(la.gt.m1) then
           if(lja(la-1).eq.1000004) goto 830
         end if
         if(idna(la).lt.0) goto 830
         jb = namchg(2,abs(lja(la)))
         if(jb.gt.0.and.jb.le.mxj) goto 820
         if(abs(lja(la)).gt.1000) goto 830

            write(io,'(''** Error in checking cell (macrobody),''/
     &      ''   surface = '',i7,'' is not found for cell''/
     &      ''   cell id ='',i7 )')
     &      lja(la),ncl(ic)

            ierr = ierr + 1

         go to 830

  820    if(ksm(jb).lt.0) goto 840
  830 continue
         goto 890
  840    iz = 0

*-----------------------------------------------------------------------
*        deal with macrobody surface on cell ic -- expand to surfaces.
*-----------------------------------------------------------------------

         ix = lsc(jb)
         jz = 6
         ji = -ksm(jb)
         if(ji.eq.33.or.ji.eq.34.or.ji.eq.36) jz = 3
         if(ji.eq.39) jz = 8
         if(ji.eq.37) jz = 5
         if(ji.eq.1030.or.ji.eq.1031) jz = 4
         if(idna(la).ge.0) goto 850

         if(abs(idna(la)).gt.jz) then

            write(io,'(''** Error in checking cell (macrobody),''/
     &      ''   surface = '',i7,
     &      '' : the facet does not exist on cell''/
     &      ''   cell id ='',i7 )')
     &      lja(la),ncl(ic)

            ierr = ierr + 1

         end if

         m1 = la+1
         goto 810

*-----------------------------------------------------------------------
*        slot is already occupied,   lh=total additional.
*        include unions if outside macrobody -- positive sense.
*        move data up to make room for additional surfaces.
*-----------------------------------------------------------------------

  850    lh = jz-1
         if(lja(la).gt.0) lh = lh+(jz+1)

      do 860 ii = 1, nlja-la
         i = nlja-la+1-ii
         idna(la+lh+i) = idna(la+i)
  860    lja(la+lh+i) = lja(la+i)

         mq = lja(la)
         idna(la) = 1
         if(mq.gt.0) lja(la+1) = lja(la)
         if(mq.gt.0) lja(la) = 1000001
         if(mq.gt.0) lja(la+lh) = 1000002
         if(mq.gt.0) idna(la+1) = idna(la)
         if(mq.gt.0) idna(la) = 0
         if(mq.gt.0) idna(la+lh) = 0
         nlja = nlja+lh
         nljc = nljc+lh

      do 870 i = ic, mxa
  870    lca(i+1) = lca(i+1)+sign(lh,lca(i+1))
         if(mq.gt.0) lca(ic) = -abs(lca(ic))

*-----------------------------------------------------------------------
*        include surfaces for macrobody.
*-----------------------------------------------------------------------

         jj = la+1
         if(mq.gt.0) jj = jj+1

      do 880 i = 2, jz
         if(mq.gt.0) lja(jj) = 1000003
         if(mq.gt.0) idna(jj) = 0
         if(mq.gt.0) jj = jj+1
         lja(jj) = mq
         idna(jj) = i
         jj = jj+1
  880 continue

         if(mq.lt.0) jj = jj-1
         m1 = jj+1
         goto 810
  890    continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine mbodyo(mm,k2,k3,k4,k5,ab,ac,sc,io,ierr,m1c,idsn)
*                                                                      *
*        mm=1:  check orthogonality of vectors k2 and k3, adjust k3.   *
*        mm=2:  check cross product of k2 and k3 as k4, compare to     *
*               k5 and adjust k5.                                      *
*                                                                      *
*       Last modified by K.Niita on 2012/12/21                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension sc(33)

      parameter ( ibmt = 40 )
      parameter ( eps5 = 1.0d-08 )

*-----------------------------------------------------------------------

      character ksf(50)*3

      data ( ksf(i), i = 1, ibmt ) /
     &         'p  ','px ','py ','pz ','so ',
     &         's  ','sx ','sy ','sz ','c/x',
     &         'c/y','c/z','cx ','cy ','cz ',
     &         'k/x','k/y','k/z','kx ','ky ',
     &         'kz ','sq ','gq ','tx ','ty ',
     &         'tz ','x  ','y  ','z  ','box',
     &         'rpp','sph','rcc','rec','ell',
     &         'trc','wed','arb','rhp','hex'/

*-----------------------------------------------------------------------

         ierr = 0

*-----------------------------------------------------------------------

      if( mm .eq. 2 ) goto 20

*-----------------------------------------------------------------------
*     check for orthogonal vectors sc(k2) and sc(k3)
*-----------------------------------------------------------------------

      ab = sc(k2)*sc(k3) + sc(k2+1)*sc(k3+1) + sc(k2+2)*sc(k3+2)

      if( ab .eq. 0. ) return

         if(abs(ab).gt.0.01) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : second vector not orthogonal.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

            ierr = ierr + 1
            return

         end if

         if(abs(ab).gt.0.000001) then

            write(io,'(''** Warning in [surface] section,''/
     &      ''   macro-body = '',a3,'' : second vector adjusted '',
     &      ''to be orthogonal.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

         end if

*-----------------------------------------------------------------------
*     adjust slightly to be orthogonal
*-----------------------------------------------------------------------

      rb = sqrt(1.-ab**2)
      do 10 kk = 1, 3
   10   sc(k3+kk-1) = (sc(k3+kk-1)-sc(k2+kk-1)*ab) / rb

      return

*-----------------------------------------------------------------------
*     find vector k4 as cross product k2&k3, compare to k5
*-----------------------------------------------------------------------

   20 call crspro(sc(k2),sc(k3),sc(k4))

      ac = sc(k5)*sc(k4) + sc(k5+1)*sc(k4+1) + sc(k5+2)*sc(k4+2)

         if(abs(ac).lt.0.99995) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : third vector not orthogonal.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

            ierr = ierr + 1
            return

         end if

         if(abs(ac).lt.0.9999995) then

            write(io,'(''** Warning in [surface] section,''/
     &      ''   macro-body = '',a3,'' : third vector adjusted '',
     &      ''to be orthogonal.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

         end if

      do 30 kk=1,3
      if(ac.gt.0.)sc(k5+kk-1)=sc(k4+kk-1)
   30 if(ac.le.0.)sc(k5+kk-1)=-sc(k4+kk-1)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine mbodyp(mm,ii,ij,ik,iy,j,jj,jn,ia,x,h,sc)
*                                                                      *
*        provide planes for macrobodies.                               *
*        mm=1 to 8 treated differently according to the call.          *
*                                                                      *
*       Last modified by K.Niita on 2012/12/21                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      dimension sc(33)
      dimension ia(3)

*-----------------------------------------------------------------------

      x=-x
      if(mm.eq.4.and.j.ne.3)x=-x
      if(mm.ne.7.and.ij.eq.0)go to 30
      if(mm.eq.7)then
       if(ia(jj).eq.0)go to 30
      endif
      if(mm.lt.4.and.j.eq.2.and.sc(3+ij).lt.0.)go to 10
      if(mm.lt.4.and.j.eq.3.and.sc(3+ij).gt.0.)go to 10
      if((mm.eq.4.or.mm.eq.5).and.x*sc(12+(jj-1)*3+ij).lt.0.)go to 10
      if(mm.eq.6.and.x.lt.0.)go to 10
      if(mm.eq.7)then
       if(x*sc(15+jj*3+ia(jj)).lt.0.)go to 10
      endif
      if(mm.eq.8.and.x.lt.0..and.sc(15+ij).gt.0.)go to 10
      if(mm.eq.8.and.x.gt.0..and.sc(15+ij).lt.0.)go to 10
      if(mm.ne.7)kst(mxj+j)=1+ij
      if(mm.eq.7)kst(mxj+j)=1+ia(jj)
      if(mm.ne.6.and.mm.ne.7)scf(iy+1)=sc(ij)
      if(mm.lt.4.and.sc(3+ij).gt.0.)scf(iy+1)=sc(ij)+sc(3+ij)
      if((mm.eq.4.or.mm.eq.5).and.x.gt.0.)
     & scf(iy+1)=sc(ij)+sc(12+(jj-1)*3+ij)*sc(21+jj)
      if(mm.eq.6)scf(iy+1)=sc(j+1)
      if(mm.eq.6.and.jn.eq.4.and.ik-10.le.jj)scf(iy+1)=sc(j+3)
      if(mm.eq.7)scf(iy+1)=sc(ia(jj))+x*sc(15+jj*3+ia(jj))*sc(30+jj)
      if(mm.eq.8.and.x.gt.0.)scf(iy+1)=sc(ij)+sc(15+ij)*h
      iy=iy+1
      go to 50

*-----------------------------------------------------------------------
*     set up a general plane for - plane to keep negative sense.
*-----------------------------------------------------------------------

   10 do 20 k=1,4
      scf(iy+k)=0.
      if((mm.lt.6.or.mm.eq.8).and.k.le.3.and.ij.eq.k)scf(iy+k)=-1.
      if(mm.eq.6.and.k.le.3.and.jn.eq.6.and.jj.eq.k)scf(iy+k)=x
      if(mm.eq.6.and.k.le.3.and.jn.eq.4.and.ik-10.gt.jj.and.jj.eq.k)
     & scf(iy+k)=x
      if(mm.eq.6.and.k.le.3.and.jn.eq.4.and.ik-10.le.jj.and.jj+1.eq.k)
     & scf(iy+k)=x
      if(mm.eq.7)then
       if(k.le.3.and.ia(jj).eq.k)scf(iy+k)=-1.
      endif
      if(k.ne.4)go to 20
      if(mm.lt.4.and.sc(3+ij).gt.0.)scf(iy+k)=-sc(ij)
      if(mm.lt.4.and.sc(3+ij).lt.0.)scf(iy+k)=-(sc(ij)+sc(ij+3))
      if((mm.eq.4.or.mm.eq.5.or.mm.eq.8).and.x.lt.0.)scf(iy+k)=-sc(ij)
      if((mm.eq.4.or.mm.eq.5).and.x.gt.0.)
     & scf(iy+k)=-(sc(ij)+sc(12+(jj-1)*3+ij)*sc(21+jj))
      if(mm.eq.7)scf(iy+k)=-(sc(ia(jj))+x*sc(15+jj*3+ia(jj))*sc(30+jj))
      if(mm.eq.8.and.x.gt.0.)scf(iy+k)=-(sc(ij)+sc(15+ij)*h)
      if(mm.ne.6)go to 20
      scf(iy+k)=-sc(j-1)
      if(jn.eq.4.and.ik-10.le.jj)scf(iy+k)=-sc(j+1)
   20 continue
      iy=iy+4
      go to 50

*-----------------------------------------------------------------------
*     off-axis plane.
*-----------------------------------------------------------------------

   30 do 40 k=1,4
      if(mm.eq.1)scf(iy+k)=sc(3+k)*x/h
      if(mm.eq.2)scf(iy+k)=sc(8+k)*x
      if(mm.eq.3)scf(iy+k)=sc(7+k)*x
      if(mm.ge.4.and.mm.ne.7.and.mm.ne.8)scf(iy+k)=sc(12+(jj-1)*3+k)*x
      if(mm.eq.7)scf(iy+k)=sc(15+jj*3+k)*x
      if(mm.eq.8)scf(iy+k)=sc(15+k)*x
      if(k.ne.4)go to 40
      scf(iy+k)=scf(iy+1)*sc(1)+scf(iy+2)*sc(2)+scf(iy+3)*sc(3)
      if(mm.lt.4.and.j.eq.2)scf(iy+k)=scf(iy+k)+
     & scf(iy+1)*sc(4)+scf(iy+2)*sc(5)+scf(iy+3)*sc(6)
      if(mm.ge.4.and.mm.ne.7.and.mm.ne.8.and.x.gt.0.)
     & scf(iy+k)=scf(iy+k)+(scf(iy+1)*sc(13+(jj-1)*3)+
     & scf(iy+2)*sc(14+(jj-1)*3)+scf(iy+3)*sc(15+(jj-1)*3))*sc(21+jj)
      if(mm.eq.7)scf(iy+k)=scf(iy+k)+x*(scf(iy+1)*sc(16+jj*3)+
     & scf(iy+2)*sc(17+jj*3)+scf(iy+3)*sc(18+jj*3))*sc(30+jj)
      if(mm.eq.8.and.x.gt.0.)scf(iy+k)=scf(iy+k)+(scf(iy+1)*sc(16)+
     & scf(iy+2)*sc(17)+scf(iy+3)*sc(18))*h
   40 continue

      iy=iy+4

*-----------------------------------------------------------------------

   50 return
      end


************************************************************************
*                                                                      *
      subroutine mbodyr(ia,n,sc,io,ierr,m1c,idsm)
*                                                                      *
*       special vector computations for rhp.                           *
*                                                                      *
*       Last modified by K.Niita on 2012/12/21                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      dimension sc(33)
      dimension ia(3)

      parameter ( ibmt = 40 )

*-----------------------------------------------------------------------

      character ksf(50)*3

      data ( ksf(i), i = 1, ibmt ) /
     &         'p  ','px ','py ','pz ','so ',
     &         's  ','sx ','sy ','sz ','c/x',
     &         'c/y','c/z','cx ','cy ','cz ',
     &         'k/x','k/y','k/z','kx ','ky ',
     &         'kz ','sq ','gq ','tx ','ty ',
     &         'tz ','x  ','y  ','z  ','box',
     &         'rpp','sph','rcc','rec','ell',
     &         'trc','wed','arb','rhp','hex'/

*-----------------------------------------------------------------------
*     set facet for regular hexagon if default, and
*     check that facet is orthogonal to axis.
*-----------------------------------------------------------------------

      do 110 kk=1,3
      ia(kk)=0
      if(kk.ne.1.and.n.le.9)go to 110

*-----------------------------------------------------------------------
*     second facet vector is set orthogonal to axis and
*     60 degrees from first facet.
*     set ia(kk) to axis if facet kk along an axis and
*     obtain normalized vector.
*-----------------------------------------------------------------------

      if(sc(5+kk*3).eq.0..and.sc(6+kk*3).eq.0.)ia(kk)=1
      if(sc(4+kk*3).eq.0..and.sc(6+kk*3).eq.0.)ia(kk)=2
      if(sc(4+kk*3).eq.0..and.sc(5+kk*3).eq.0.)ia(kk)=3
      r=sqrt(sc(4+kk*3)**2+sc(5+kk*3)**2+sc(6+kk*3)**2)
      sc(16+kk*3)=sc(4+kk*3)/r
      sc(17+kk*3)=sc(5+kk*3)/r
      sc(18+kk*3)=sc(6+kk*3)/r
      sc(30+kk)=r
      if(kk.ne.1.and.n.le.9)go to 10

*-----------------------------------------------------------------------
*     orthogonality check of facet vector and axis
*-----------------------------------------------------------------------

      call mbodyo(1,16,16+3*kk,0,0,ab,ac,sc,io,ierr0,m1c,idsn)

            if( ierr0 .gt. 0 ) then

              ierr = ierr + 1
              return

            end if

   10 if(kk.ne.1.or.n.gt.9)go to 110
      sc(32)=r
      sc(33)=r

*-----------------------------------------------------------------------
*     second facet vector is set orthogonal to axis and
*     60 degrees from first facet in c-clockwise direction.
*-----------------------------------------------------------------------

      mq=0
      qc=0.
      qd=0.
      if(sc(16).eq.0.)go to 30
      m=1
      d=sc(20)-sc(19)*sc(17)/sc(16)
      if(abs(d).gt.1.e-7)go to 20
      mq=1
      d=sc(21)-sc(19)*sc(18)/sc(16)
      dz=.5/d
      qa=-sc(18)*dz/sc(16)
      qb=-sc(17)/sc(16)
      go to 70
   20 qa=.5/d
      qb=((sc(19)*sc(18)/sc(16))-sc(21))/d
      qc=-qa*sc(17)/sc(16)
      qd=(-qb*sc(17)/sc(16))-sc(18)/sc(16)
      go to 70
   30 if(sc(17).eq.0.)go to 50
      m=2
      d=sc(21)-sc(20)*sc(18)/sc(17)
      if(abs(d).gt.1.e-7)go to 40
      mq=1
      d=sc(19)-sc(20)*sc(16)/sc(17)
      dz=.5/d
      qa=-sc(16)*dz/sc(17)
      qb=-sc(18)/sc(17)
      go to 70
   40 qa=.5/d
      qb=((sc(20)*sc(16)/sc(17))-sc(19))/d
      qc=-qa*sc(18)/sc(17)
      qd=(-qb*sc(18)/sc(17))-sc(16)/sc(17)
      go to 70
   50 m=3
      d=sc(19)-sc(21)*sc(16)/sc(18)
      if(abs(d).gt.1.e-7)go to 60
      mq=1
      d=sc(20)-sc(21)*sc(17)/sc(18)
      dz=.5/d
      qa=-sc(17)*dz/sc(18)
      qb=-sc(16)/sc(18)
      go to 70
   60 qa=.5/d
      qb=((sc(21)*sc(17)/sc(18))-sc(20))/d
      qc=-qa*sc(16)/sc(18)
      qd=(-qb*sc(16)/sc(18))-sc(17)/sc(18)
   70 a=qd**2+qb**2+1
      b=2.*(qc*qd+qa*qb)
      c=qc**2+qa**2-1.0
      if(mq.eq.1)c=c+dz**2
      q=b**2-4.*a*c
      q=sqrt(q)

*-----------------------------------------------------------------------
*     first select root with plus sign.
*-----------------------------------------------------------------------

      call mbodys(1,m,mq,a,b,dz,q,qa,qb,sc)

*-----------------------------------------------------------------------
*     check if right root?
*-----------------------------------------------------------------------

      call crspro(sc(22),sc(19),sc(28))

      ba=(sc(16)*sc(28)+sc(17)*sc(29)+sc(18)*sc(30))/(0.5*sqrt(3.))

         if(abs(ba).lt.0.99995.or.abs(ba).gt.1.00005) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : failed to create second facet vector.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

            ierr = ierr + 1
            return

         end if

      if(ba.lt.0.)go to 90

*-----------------------------------------------------------------------
*     switch to other root.
*-----------------------------------------------------------------------

      q=-q
      call mbodys(1,m,mq,a,b,dz,q,qa,qb,sc)
      q=-q

*-----------------------------------------------------------------------
*     third facet vector is set orthogonal to axis and
*     120 degrees from first facet in c-clockwise direction.
*-----------------------------------------------------------------------

   90 qa=-qa
      qc=-qc
      b=-b
      if(mq.eq.1)dz=-dz

*-----------------------------------------------------------------------
*     first select root with plus sign.
*-----------------------------------------------------------------------

      call mbodys(2,m,mq,a,b,dz,q,qa,qb,sc)

  100 call crspro(sc(25),sc(19),sc(28))
      ba=(sc(16)*sc(28)+sc(17)*sc(29)+sc(18)*sc(30))/(0.5*sqrt(3.))

         if(abs(ba).lt.0.99995.or.abs(ba).gt.1.00005) then

            write(io,'(''** Error in [surface] section,''/
     &      ''   macro-body = '',a3,
     &      '' : failed to create third facet vector.''/
     &      ''   surface id ='',i7 )')
     &      ksf(m1c), idsn

            ierr = ierr + 1
            return

         end if

      if(ba.lt.0.)go to 110

*-----------------------------------------------------------------------
*     switch to other root.
*-----------------------------------------------------------------------

      q=-q
      call mbodys(2,m,mq,a,b,dz,q,qa,qb,sc)
      q=-q
  110 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine mbodys(mm,m,mq,a,b,dz,q,qa,qb,sc)
*                                                                      *
*       vector computations to set an sc vector for rhp.
*       mm controls vector set.
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension sc(33)

*-----------------------------------------------------------------------

      j=22
      if(mm.eq.2)j=25
      if(m.ne.1)go to 20
      if(mq.eq.0)go to 10
      sc(j+2)=dz
      sc(j+1)=(-b+q)/(2.*a)
      sc(j)=-(sc(17)*sc(j+1)+sc(18)*sc(j+2))/sc(16)
      go to 60
   10 sc(j+2)=(-b+q)/(2.*a)
      sc(j+1)=qa+qb*sc(j+2)
      sc(j)=-(sc(17)*sc(j+1)+sc(18)*sc(j+2))/sc(16)
      go to 60
   20 if(m.eq.3)go to 40
      if(mq.eq.0)go to 30
      sc(j)=dz
      sc(j+2)=(-b+q)/(2.*a)
      sc(j+1)=-(sc(18)*sc(j+2)+sc(16)*sc(j))/sc(17)
      go to 60
   30 sc(j)=(-b+q)/(2.*a)
      sc(j+2)=qa+qb*sc(j)
      sc(j+1)=-(sc(18)*sc(j+2)+sc(16)*sc(j))/sc(17)
      go to 60
   40 if(mq.eq.0)go to 50
      sc(j+1)=dz
      sc(j)=(-b+q)/(2.*a)
      sc(j+2)=-(sc(16)*sc(j)+sc(17)*sc(j+1))/sc(18)
      go to 60
   50 sc(j+1)=(-b+q)/(2.*a)
      sc(j)=qa+qb*sc(j+1)
      sc(j+2)=-(sc(16)*sc(j)+sc(17)*sc(j+1))/sc(18)

   60 return
      end


************************************************************************
*                                                                      *
      subroutine sufchk(io,ierr,ix,nx,m0c,m1c,idsn)
*                                                                      *
*       set up and check surface.                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      parameter ( ibmt = 40 )

*-----------------------------------------------------------------------

      character ksf(50)*3

      data ( ksf(i), i = 1, ibmt ) /
     &         'p  ','px ','py ','pz ','so ',
     &         's  ','sx ','sy ','sz ','c/x',
     &         'c/y','c/z','cx ','cy ','cz ',
     &         'k/x','k/y','k/z','kx ','ky ',
     &         'kz ','sq ','gq ','tx ','ty ',
     &         'tz ','x  ','y  ','z  ','box',
     &         'rpp','sph','rcc','rec','ell',
     &         'trc','wed','arb','rhp','hex'/

*-----------------------------------------------------------------------

   50          ii = 0

            if( m1c .ne. 1 .and.
     &          m1c .ne. 6 .and.
     &          abs(m1c-17) .gt. 1 ) goto 60

               if( scf(ix+1) .ne. 0. ) ii = 1
               if( scf(ix+2) .ne. 0. ) ii = 2 + 10 * ii
               if( scf(ix+3) .ne. 0. ) ii = 3 + 10 * ii

   60       continue

*-----------------------------------------------------------------------

            if( m1c .eq.  2 .or.
     &          m1c .eq.  3 .or.
     &          m1c .eq.  4 .or.
     &          m1c .eq. 22 ) goto 180

*-----------------------------------------------------------------------
*        p plane -- check for improper coefficients.
*-----------------------------------------------------------------------

            if( m1c .eq. 1 ) then

               if( ii .eq. 0 ) then

                  write(io,'(/''** ERROR : in [surface] section, '',
     &                        ''a, b, and c are all zero.''/
     &                        ''   surface id ='',i7)')  idsn

                  ierr = ierr + 1

                  goto 180

               end if

                  if( abs(ii-2) .gt. 1 )   goto 180
                  if( scf(ix+ii) .lt. 0. ) goto 180

                  scf(ix+1) = scf(nx) / scf(ix+ii)
                  if( ii .ne. 1 ) scf(ix+ii) = 0.
                  m1c = ii + 1

                  goto 120

            end if

*-----------------------------------------------------------------------
*        spheres and cylinders -- square the radius.
*-----------------------------------------------------------------------

            if( m1c .ge. 5 .and. m1c .le. 15 ) then

               scf(nx) = scf(nx)**2

               if( m1c .eq.  5 .or.
     &             m1c .gt. 12 .or.
     &             ii  .gt. 3 ) goto 180

               if( m1c .eq. 6 ) then

                  if( ii .eq. 0 ) scf(ix+1) = scf(nx)
                  if( ii .ne. 0 ) scf(ix+1) = scf(ix+ii)
                  if( ii .ne. 0 ) scf(ix+2) = scf(nx)

                  scf(ix+3) = 0.

                  if( ii .eq. 0 ) m1c = 5
                  if( ii .ne. 0 ) m1c = 6 + ii

                  goto 120

               end if

               if( scf(ix+1) .ne. 0. .or.
     &           ( m1c .gt. 9 .and. scf(ix+2) .ne. 0. ) ) goto 180

                  scf(ix+1) = scf(nx)
                  if( m1c .ge. 10 ) m1c = m1c + 3
                  if( m1c .le.  9 ) m1c = 5

                  goto 120

            end if

*-----------------------------------------------------------------------
*        cones -- set up signed slope for cone of one sheet.
*-----------------------------------------------------------------------

            if( m1c .ge. 16 .and. m1c .le. 21 ) then

               if( scf(nx) .ne. 0. )
     &             scf(nx) = sign( sqrt( abs(scf(nx-1)) ) , scf(nx) )

               if( scf(nx-1) .lt. 0. ) then

                  write(io,'(/''** ERROR : in [surface] section, '',
     &                        ''t**2 is negative.''/
     &                        ''   surface id ='',i7)')  idsn

                  ierr = ierr + 1

                  goto 180

               end if

               if( m1c .gt. 18 .or.
     &           ( ii .ne. 0 .and. ii .ne. m1c-15 ) ) goto 180

                  scf(ix+1) = scf(ix+ii)
                  scf(ix+2) = scf(nx-1)
                  scf(ix+3) = scf(nx)
                  scf(nx-1) = 0.
                  m1c = m1c + 3

                  goto 120

            end if

*-----------------------------------------------------------------------
*        gq -- convert symmetric gq to sq.
*-----------------------------------------------------------------------

            if( m1c .eq. 23 ) then

               if( scf(ix+4) .ne. 0. .or.
     &             scf(ix+5) .ne. 0. .or.
     &             scf(ix+6) .ne. 0. ) goto 180

                  scf(ix+4) = .5*scf(ix+7)
                  scf(ix+5) = .5*scf(ix+8)
                  scf(ix+6) = .5*scf(ix+9)
                  scf(ix+7) = scf(ix+10)
                  scf(ix+8) = 0.
                  scf(ix+9) = 0.

                  m1c = 22

                  goto 120

            end if

*-----------------------------------------------------------------------
*        tori -- set up square of ellipse axis ratio.
*-----------------------------------------------------------------------

            if( m1c .ge. 24 .and. m1c .le. 26 ) then

               if( scf(nx)   .le. 0. .or.
     &             scf(nx-1) .le. 0. ) then

                  write(io,'(/''** ERROR : in [surface] section, ''/
     &            ''   at least one ellipse semi-axis length is zero.''/
     &            ''   surface id ='',i7)')  idsn

                  ierr = ierr + 1

                  goto 180

               end if

               if( scf(nx-2) .ne. 0. ) then

                  ix = ix + 1
                  scf(nx+1) = ( scf(nx) / scf(nx-1) )**2

                  if( scf(nx) .gt. scf(nx-2) ) then

                     write(io,'(''** warning : in [surface] section, ''/
     &               ''   degenerate torus.'',
     &               ''   surface id ='',i7)')  idsn

                  end if

                  if( scf(nx) .eq. scf(nx-2) ) then

                     write(io,'(''** warning : in [surface] section, ''/
     &               ''   singular torus.  it can fail in plotting '',
     &                                   ''or tracking.''/
     &               ''   surface id ='',i7)')  idsn

                  end if

                  goto 180

               end if

*-----------------------------------------------------------------------
*           change a completely degenerate torus into an ellipsoid.
*-----------------------------------------------------------------------

                  scf(nx+1) = -scf(nx)**2
                  a = ( scf(nx) / scf(nx-1) )**2

                  do ii = 1, 3

                     scf(nx+1+ii) = scf(nx-6+ii)
                     scf(nx-6+ii) = 1.
                     scf(nx-3+ii) = 0.

                  end do

                  scf(nx+m1c-29) = a

                  m1c = 22

                  write(io,'(''** warning : in [surface] section, ''/
     &            ''   completely degenerate torus has been '',
     &                ''replaced by ellipsoid.''/
     &            ''   surface id ='',i7,
     &            '' :  surface type  '',a3,'' -> '',a3)')
     &                     idsn, ksf(m0c), ksf(m1c)

                  goto 180

            end if

*-----------------------------------------------------------------------
*        comment on surface simplifications.
*-----------------------------------------------------------------------

  120    continue

            scf(nx) = 0.

               write(io,'( ''** warning : in [surface] section, '',
     &                     ''surface has been replaced''/
     &                     ''   surface id ='',i7,
     &                     '' :  surface type  '',a3,'' -> '',a3)')
     &                     idsn, ksf(m0c), ksf(m1c)

            goto 180

*-----------------------------------------------------------------------

  180 return
      end


************************************************************************
*                                                                      *
      subroutine addsuf(io,ierr,ic)
*                                                                      *
*       create new surfaces for trcl cell ic.                          *
*       Last modified by K.Niita on 2010/01/13                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      dimension nc(26)
      data nc/4,4*1,4,3*2,3*3,3*1,3*5,3*3,2*10,3*7/

*-----------------------------------------------------------------------
*     loop over all surfaces of cell ic.
*-----------------------------------------------------------------------

      do 60 js = abs(lca(ic)), abs(lca(ic+1))-1

*-----------------------------------------------------------------------
*        generate a new surface for each unique surface of the cell.
*-----------------------------------------------------------------------

         if(lja(js).gt.1000000) goto 60
         if(js.gt.1) then
           if(lja(js-1).eq.1000004) goto 60
         end if

      do 10 k = abs(lca(ic)),js-1
         if(abs(lja(k)).eq.abs(lja(js)))then
           if(k.gt.1)then
             if(lja(k-1).ne.1000004) goto 60
           else
             goto 60
           endif
         endif
   10 continue
         mxj = mxj+1

*-----------------------------------------------------------------------
*        generate a unique name for the new surface.
*-----------------------------------------------------------------------

         if( ncl(ic) .le. 100 ) then
            ns = abs(lja(js))+9999*ncl(ic)
         else if( ncl(ic) .le. 1000 ) then
            ns = abs(lja(js))+999*ncl(ic)
         else if( ncl(ic) .le. 10000 ) then
            ns = abs(lja(js))+99*ncl(ic)
         else if( ncl(ic) .le. 100000 ) then
            ns = abs(lja(js))+9*ncl(ic)
         else
            ns = abs(lja(js))+ncl(ic)
         end if


         if( ns .gt. 1000000 ) ns = ns - ns/1000000*1000000

*-----------------------------------------------------------------------
         if(namchg(2,ns).ne.0) then

            insss = 0
  111       continue
            insss = insss + 1
            if( insss .gt. 10000 ) goto 112

            ns = ns + 1
            if(namchg(2,ns).ne.0) goto 111
            goto 113

  112       continue

            write(io,'(/''** ERROR : from s.addsuf, ''/
     &                  ''   generated surface name = '',i7/
     &                  ''   in cell'',i7,'' is not unique.'')')
     &                  ns, ncl(ic)

            ierr = ierr + 1
            return

         end if

  113    continue
*-----------------------------------------------------------------------

         nsfm(mxj) = ns
         jj = namchg(2,abs(lja(js)))

         if(jj.eq.0) then

            write(io,'(/''** ERROR : from s.addsuf, ''/
     &                  ''   surface = '',i7/
     &                  ''   of cell'',i7,'' is not defined.'')')
     &                  lja(js), ncl(ic)

            ierr = ierr + 1
            return

         end if

         nf = 1

*-----------------------------------------------------------------------
*        if macrobody, add the master macrobody surface.
*-----------------------------------------------------------------------

         if(idna(js).eq.0) goto 20
         jtr(mxj) = jtr(jj)
         kst(mxj) = kst(jj)
         ksu(mxj) = ksu(jj)
         ksm(mxj) = ksm(jj)
         lsc(mxj+1) = lsc(mxj)
         nf = 6
         ji = -ksm(jj)
         if(ji.eq.33.or.ji.eq.34.or.ji.eq.36) nf = 3
         if(ji.eq.39) nf = 8
         if(ji.eq.37) nf = 5
         if(ji.eq.1030.or.ji.eq.1031) nf = 4
         jj=jj+1
         mxj = mxj+1
         nsfm(mxj) = 0

*-----------------------------------------------------------------------
*        change the surface name.
*-----------------------------------------------------------------------

   20 do 30 kz = js, abs(lca(ic+1))-1
         k = js+abs(lca(ic+1))-1-kz
         if(abs(lja(k)).eq.abs(lja(js))) then
           if(k.gt.1) then
             if(lja(k-1).ne.1000004) lja(k) = sign(ns,lja(k))
           else
             lja(k) = sign(ns,lja(k))
           end if
         end if
   30 continue

*-----------------------------------------------------------------------
*        copy the attributes of the old surface.
*        if macrobody, do for each facet.
*-----------------------------------------------------------------------

      do 50 jf = 1, nf
         kst(mxj) = kst(jj)
         ksu(mxj) = ksu(jj)
         ksm(mxj) = ksm(jj)
         if(nf.ne.1) ksm(mxj) = mxj-jf
         lsc(mxj+1) = lsc(mxj)+nc(kst(mxj))

      do 40 i = 1, nc(kst(mxj))
   40    scf(lsc(mxj)+i) = scf(lsc(jj)+i)

*-----------------------------------------------------------------------
*        transform the surface coefficients.
*-----------------------------------------------------------------------

         jtr(mxj) = ktr(ic)

         if( kst(mxj).ge.24 .and. kst(mxj).le.26 ) then
           if( jtr(jj).ne.0 ) then
             do jt = 1, mxtr
               if( abs(trf(1,jt)).eq.ktr(ic) ) m1 = jt
               if( abs(trf(1,jt)).eq.jtr(jj) ) m2 = jt
             end do
             jtr(mxj) = -( 10000 * m1 + m2 )
           end if
         else
           call trfsuf(io,ierr,mxj)
         end if

            if( ierr .gt. 0 ) return

         if(jf.eq.nf) goto 50
         jj = jj + 1
         mxj = mxj + 1
         nsfm(mxj) = 0

   50 continue
   60 continue

*-----------------------------------------------------------------------
         if(kst(mxj) .eq. 27) then
          write(ErrCha,*)
     & 'trcl cannot be applied to cells containing p'
          ErrID = 'L:4469/R:addsuf/F:ggs01.f' !E07_001_001
          call ErrWrite(ErrID,ErrCha)
          stop
         endif
*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine callat(io,ierr)
*                                                                      *
*       calculate the lattice vectors and search constants.            *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************

      use LAFDATAMOD !FURUTA2020127
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      dimension sc(4,8),xp(3),tn(8)
      character ha*5

*-----------------------------------------------------------------------

      do 190 ic = 1, mxa
         if(lat(1,ic).eq.0.or.abs(lat(1,ic)).eq.3) goto 190
         !FURUTA20150714 TETRA
*-----------------------------------------------------------------------
*        put the lattice cell number into mfl or laf.
*-----------------------------------------------------------------------

         if(jun(ic).eq.0) goto 30

      do 20 i = 1, mxa
         if(i.eq.ic) goto 20
         if(mfl(1,i).eq.abs(jun(ic))) mfl(2,i) = ic
         if(mfl(1,i).ge.0) goto 20
         lp = -mfl(1,i)

      do 10 j = 3, 2+laf(1,lp+2)*laf(2,lp+2)*laf(3,lp+2)       !FURUTA20201127
   10    if(laf(1,lp+j).eq.abs(jun(ic))) laf(2,lp+j) = ic !FURUTA20201127

   20 continue

*-----------------------------------------------------------------------
*        put the coefficients of the surfaces into the same form.
*-----------------------------------------------------------------------

   30    n = abs(lca(ic+1))-lca(ic)
         ip = lat(2,ic)

      do 90 j = 1, n
         js = abs(lja(lca(ic)+j-1))
         ks = abs(kst(js))
         ls = lsc(js)
         if(ks.ne.1) goto 50
         tn(j) = sqrt(scf(ls+1)**2+scf(ls+2)**2+scf(ls+3)**2)

      do 40 i = 1, 4
   40    sc(i,j) = scf(ls+i)
         goto 70

   50 do 60 i = 1, 4
   60    sc(i,j) = 0.
         tn(j) = 1.
         sc(ks-1,j) = 1.
         sc(4,j) = scf(ls+1)
   70    if(mod(j,2).eq.1) goto 80
         a = sc(1,j)*sc(1,j-1)+sc(2,j)*sc(2,j-1)+sc(3,j)*sc(3,j-1)
         if(a.lt.0.) tn(j) = -tn(j)
         if(abs(a-tn(j)).gt..001) goto 200
   80 do 90 i = 1, 4
   90    sc(i,j) = sc(i,j)/tn(j)
         if(lat(1,ic).eq.2) goto 150

*-----------------------------------------------------------------------
*        hexahedron case
*-----------------------------------------------------------------------

         if(n.gt.2) goto 110

      do 100 i = 1, 3
         vcl(lvcl+i,5,ip) = sc(i,1)/(sc(4,1)-sc(4,2))
  100    vcl(lvcl+i,1,ip) = (sc(4,1)-sc(4,2))*sc(i,1)
         vcl(lvcl+1,4,ip) = .5*(sc(4,1)+sc(4,2))/(sc(4,1)-sc(4,2))
         goto 190

  110    if(n.eq.4) call crspro(sc(1,1),sc(1,3),sc(1,5))
                    call crspro(sc(1,3),sc(1,5),xp)

         t1 = xp(1)*sc(1,1)+xp(2)*sc(2,1)+xp(3)*sc(3,1)
         if(t1.eq.0.) goto 200
         t1 = (sc(4,1)-sc(4,2))/t1

      do 120 i = 1, 3
         vcl(lvcl+i,5,ip) = sc(i,1)/(sc(4,1)-sc(4,2))
  120    vcl(lvcl+i,1,ip) = t1*xp(i)
         vcl(lvcl+1,4,ip) = .5*(sc(4,1)+sc(4,2))/(sc(4,1)-sc(4,2))

         call crspro(sc(1,5),sc(1,1),xp)

         t1 = xp(1)*sc(1,3)+xp(2)*sc(2,3)+xp(3)*sc(3,3)
         if(t1.eq.0.) goto 200
         t1 = (sc(4,3)-sc(4,4))/t1

      do 130 i = 1, 3
         vcl(lvcl+i,6,ip) = sc(i,3)/(sc(4,3)-sc(4,4))
  130    vcl(lvcl+i,2,ip) = t1*xp(i)
         vcl(lvcl+2,4,ip) = .5*(sc(4,3)+sc(4,4))/(sc(4,3)-sc(4,4))
         if(n.eq.4) goto 190

         call crspro(sc(1,1),sc(1,3),xp)

         t1 = xp(1)*sc(1,5)+xp(2)*sc(2,5)+xp(3)*sc(3,5)
         if(t1.eq.0.) goto 200
         t1 = (sc(4,5)-sc(4,6))/t1

      do 140 i = 1, 3
         vcl(lvcl+i,7,ip) = sc(i,5)/(sc(4,5)-sc(4,6))
  140    vcl(lvcl+i,3,ip) = t1*xp(i)
         vcl(lvcl+3,4,ip) = .5*(sc(4,5)+sc(4,6))/(sc(4,5)-sc(4,6))
         goto 190

*-----------------------------------------------------------------------
*        hexagonal prism case.
*-----------------------------------------------------------------------

  150    if(n.eq.6) call crspro(sc(1,1),sc(1,3),sc(1,7))

      do 170 j = 1, 2
         call crspro(sc(1,6-2*j),sc(1,5),xp)

         t1 = sc(1,7)*xp(1)+sc(2,7)*xp(2)+sc(3,7)*xp(3)
         if(t1.eq.0.) goto 200

         call crspro(sc(1,6-2*j),sc(1,2*j),xp)

         s = (sc(1,7)*xp(1)+sc(2,7)*xp(2)+sc(3,7)*xp(3))/t1

         call crspro(sc(1,2*j),sc(1,5),xp)

         t = (sc(1,7)*xp(1)+sc(2,7)*xp(2)+sc(3,7)*xp(3))/t1
         bp = (sc(4,2*j-1)-s*sc(4,7-j))/t
         c = (sc(4,2*j)-t*sc(4,6-2*j))/s
         r = (bp-sc(4,6-2*j))/(c-sc(4,7-j))
         g = 1./((sc(4,4+j)-c)*r+sc(4,5-2*j)-sc(4,6-2*j))
         f = g*r

      do 160 i = 1, 3
  160    tpp(i) = f*sc(i,5)+g*sc(i,6-2*j)

         call crspro(tpp,sc(1,7),xp)

         t1 = sc(1,2*j)*xp(1)+sc(2,2*j)*xp(2)+sc(3,2*j)*xp(3)
         if(t1.eq.0.) goto 200
         t1 = (sc(4,2*j-1)-sc(4,2*j))/t1
         vcl(lvcl+3-j,4,ip) = .5*(f*(sc(4,4+j)+sc(4,7-j))
     &                      + g*(sc(4,5-2*j)+bp))
      do 170 i = 1, 3
         vcl(lvcl+i,7-j,ip) = tpp(i)
  170    vcl(lvcl+i,j,ip) = t1*xp(i)
         if(n.eq.6) goto 190

         call crspro(sc(1,1),sc(1,3),xp)

         t1 = (sc(4,7)-sc(4,8))
     &      / (xp(1)*sc(1,7)+xp(2)*sc(2,7)+xp(3)*sc(3,7))
      do 180 i = 1, 3
         vcl(lvcl+i,7,ip) = sc(i,7)/(sc(4,7)-sc(4,8))
  180    vcl(lvcl+i,3,ip) = t1*xp(i)
         vcl(lvcl+3,4,ip) = .5*(sc(4,7)+sc(4,8))/(sc(4,7)-sc(4,8))
  190    continue
         return

  200    continue

         write(io,'(/''** ERROR : in lattice setup,''/
     &   ''    wrong order of surfaces in lattice cell ='',i7)')
     &   ncl(ic)

         ierr = ierr + 1

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine celsuf(io,ierr)
*                                                                      *
*       finish up the cells and surfaces.                              *
*       Last modified by K.Niita on 2009/10/06                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      parameter ( ibmt = 40 )

*-----------------------------------------------------------------------

      character ksf(50)*3

      data ( ksf(i), i = 1, ibmt ) /
     &         'p  ','px ','py ','pz ','so ',
     &         's  ','sx ','sy ','sz ','c/x',
     &         'c/y','c/z','cx ','cy ','cz ',
     &         'k/x','k/y','k/z','kx ','ky ',
     &         'kz ','sq ','gq ','tx ','ty ',
     &         'tz ','x  ','y  ','z  ','box',
     &         'rpp','sph','rcc','rec','ell',
     &         'trc','wed','arb','rhp','hex'/

*-----------------------------------------------------------------------

      character hj(2,3)*10,hl*5,hp*130,hq*2,ht*6

*-----------------------------------------------------------------------
*        print the names and coefficients of the surfaces.
*-----------------------------------------------------------------------

         write(io,430)
  430    format(/
     &   '*** print the names and coefficients of the surfaces.'//
     &    9x,7hsurface,5x,34htrans  type   surface coefficients/)

      do 450 js = 1, mxj
         if(idnt(js).ne.0.and.idnt(js).ne.js) goto 450
         n = lsc(js+1)-lsc(js)

      if( n.ne.0 ) then

         do 440 i = 1, n
  440       tpp(i) = scf(lsc(js)+i)

         k = kst(js)
         if(k.ge.5.and.k.le.15) tpp(n) = sqrt(tpp(n))
         if(k.ge.16.and.k.le.21.and.tpp(n).ne.0.)
     &                          tpp(n) = sign(one,tpp(n))
         if(k.ge.16.and.k.le.21.and.tpp(n).eq.0..or.k.ge.24) n = n-1
         hl = ' '
         if(ksu(js).eq.-1) hl = 'refl.'
         if(ksu(js).eq.-2) hl = 'white'
         if(ksu(js).gt.0) hl = 'pbc  '
         ht = ' '
         if(jtr(js).ne.0) write(ht,'(i6)') jtr(js)
         j1 = nsf(js)
         hq = ' '
         if(kfq.ne.0) write(hq,'(1h.,i1)') kfq
         write(io,460) js,j1,hq,hl,ht,ksf(k),(tpp(i),i=1,n)

      else

         k = kst(js)
         hl = ' '
         if(ksu(js).eq.-1) hl = 'refl.'
         if(ksu(js).eq.-2) hl = 'white'
         if(ksu(js).gt.0) hl = 'pbc  '
         ht = ' '
         if(jtr(js).ne.0) write(ht,'(i6)') jtr(js)
         j1 = nsf(js)
         hq = ' '
         if(kfq.ne.0) write(hq,'(1h.,i1)') kfq
         write(io,460) js,j1,hq,hl,ht,ksf(k)

      end if

  450 continue
  460    format(i6,i7,a2,1x,a5,1x,a6,3x,a3,1p3e16.7/,4(32x,3e16.7,/))

*-----------------------------------------------------------------------
*      print out identical surface sets.
*-----------------------------------------------------------------------

         l = 1
         hp = ' '

      do 464 j = 1, idne(1)
         if(j.eq.1) j3 = 3
         if(j.eq.1) write(io,461)
  461    format(/'*** print out identical surface sets.'//
     &    3x,14hmaster surface,5x,18hidentical surfaces/)

      do 463 j2 = 1, idne(j3)
         write(hp(l+1:l+10),'(i8,2h  )') nsf(idne(j3+j2))
         if(kfq.ne.0) write(hp(l+9:l+10),'(1h.,i1)') kfq
         if(l+20.gt.80) call wonel(io,hp,l)
         if(l.eq.1) l = 9
  463    l = l + 10

         if(l.gt.19) call wonel(io,hp,l)
  464    j3 = j3+idne(j3)+1

*-----------------------------------------------------------------------
*        print out surface coefficients for any identical surfaces.
*-----------------------------------------------------------------------

      do 469 j = 1, idne(1)
         if(j.eq.1) j3 = 3
         if(j.eq.1) write(io,465)
  465    format(/'*** surface coefficients for identical surfaces not',
     &          ' used.'//9x,7hsurface,5x,
     &          34htrans  type   surface coefficients/)

      do 468 j4 = 2, idne(j3)
         js = idne(j3+j4)
         n = lsc(js+1)-lsc(js)

      do 466 i = 1, n
  466    tpp(i) = scf(lsc(js)+i)

         k = kst(js)
         if(k.ge.5.and.k.le.15) tpp(n) = sqrt(tpp(n))
         if(k.ge.16.and.k.le.21.and.tpp(n).ne.0.)
     &                          tpp(n) = sign(one,tpp(n))
         if(k.ge.16.and.k.le.21.and.tpp(n).eq.0..or.k.ge.24) n = n-1
         hl = ' '
         if(ksu(js).eq.-1) hl = 'refl.'
         if(ksu(js).eq.-2) hl = 'white'
         if(ksu(js).gt.0) hl = 'pbc  '
         ht = ' '
         if(jtr(js).ne.0) write(ht,'(i6)') jtr(js)
         j1 = nsf(js)
         hq = ' '
         if(kfq.ne.0) write(hq,'(1h.,i1)') kfq
  468    write(io,460) js,j1,hq,hl,ht,ksf(k),(tpp(i),i=1,n)
  469    j3 = j3+idne(j3)+1

*-----------------------------------------------------------------------
*        check for unused surfaces.
*-----------------------------------------------------------------------

  470    jp = 0

      do 500 j = 1, mxj
         if(idns(j).ne.0) goto 500
         if(ksm(j).ne.0) goto 500
         if(ksu(j).gt.0) jp = 1

      do 480 i = 1, mlja
  480    if(abs(lja(i)).eq.j) goto 500
      do 490 ic = 1, mxa
         k = nsf(j)+1000*ncl(ic)
      do 490 i = abs(lca(ic)),abs(lca(ic+1))-1
         if(lja(i).gt.1000000) goto 490
         if(nsf(abs(lja(i))).eq.k) goto 500
  490 continue

               write(io,'( ''** Warning : '',
     &                     '' surface = '',i7,
     &                     '' is not used for anything'')') nsf(j)

         ksu(j) = -3
  500 continue

*-----------------------------------------------------------------------
*     process periodic boundary conditions.
*-----------------------------------------------------------------------

         if(jp.ne.0) call intprb(io,ierr)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine wonel(io,ht,l)
*                                                                      *
*       flush the one-line print buffer ht and reset l to 1.           *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************

      character ht*80

*-----------------------------------------------------------------------

      write(io,'(a80)') ht
      ht = ' '
      l = 1

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine intprb(io,ierr)
*                                                                      *
*       initiation for periodic boundary conditions.                   *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      parameter ( eps5 = 1.0d-08 )
      parameter ( eps6 = 1.0d-12 )

*-----------------------------------------------------------------------

      common /regdm/  idmg(kvlmax)

*-----------------------------------------------------------------------
*        set up trf and ksu array pointers.
*-----------------------------------------------------------------------

         mx = mxtr+1

      do 20 j1 = 1, mxj

         if(ksu(j1).le.0) goto 20
         if(kst(j1).gt.4) then

            write(io,'(/''** ERROR : from s.intprb, ''/
     &                  ''   surface = '',i7,
     &                  ''  is nonplanar. cannot be periodic.'')')
     &                  nsf(j1)

            ierr = ierr + 1
            return

         end if

         j2 = namchg(2,ksu(j1))

         if(j2.eq.0) then

            write(io,'(/''** ERROR : from s.intprb, ''/
     &                  ''   surface = '',i7,
     &                  ''  cannot be periodic with surface ='',i7)')
     &                  nsf(j1), ksu(j1)

            ierr = ierr + 1
            return

         end if

         mxtr = mxtr+1
         trf(1,mxtr) = j2
         ksu(j1) = mxtr

*-----------------------------------------------------------------------
*        load surface coefficients into trf(ltrf+4+i,k), trf(ltrf+4,k).
*        trf(ltrf+4+i,mxtr) is also the surface normal vector.
*-----------------------------------------------------------------------

         ix = lsc(j1)
         if(kst(j1).eq.1)
     &         t1 = 1./sqrt(dotpro(scf(ix+1),scf(ix+1)))
      do 10 i = 1, 3
         trf(4+i,mxtr) = 0.
   10    if(kst(j1).eq.1) trf(4+i,mxtr) = t1*scf(ix+i)
         if(kst(j1).gt.1) trf(3+kst(j1),mxtr) = 1.
         if(kst(j1).gt.1) trf(4,mxtr) = scf(ix+1)
         if(kst(j1).eq.1) trf(4,mxtr) = t1*scf(ix+4)
   20    continue

*-----------------------------------------------------------------------
*        find tpp(i), periodic surface intersection and system vertical.
*-----------------------------------------------------------------------

         jd = 1
      do 60 m1 = mx, mxtr
      do 50 m2 = m1+1, mxtr

         call crspro(trf(5,m1),trf(5,m2),tpp(jd))

         if(abs(dotpro(tpp(1),tpp(jd))).lt.eps5) goto 50
         if(abs(abs(dotpro(tpp(jd),tpp(jd)))-1.).lt.eps5) goto 40
         t1 = 1./sqrt(abs(dotpro(tpp(jd),tpp(jd))))

      do 30 i = jd, jd+2
   30    tpp(i) = t1*tpp(i)
   40    if(abs(abs(dotpro(tpp(1),tpp(jd)))-1.).gt.eps5) goto 360
         jd = 4

   50 continue
   60 continue
         if(jd.eq.1) goto 360

*-----------------------------------------------------------------------
*        check reciprocity.
*-----------------------------------------------------------------------

         db = huge
         dt = -huge

      do 100 j1 = 1, mxj
         if(kst(j1).gt.4.or.ksu(j1).eq.-3) goto 100
         if(ksu(j1).le.0) goto 70
         j2 = nint(trf(1,ksu(j1)))
         j3 = nint(trf(1,ksu(j2)))

         if(j1.ne.j3) then

            write(io,'(/''** ERROR : from s.intprb, ''/
     &                  ''   surface = '',i7,
     &                  ''  is periodic with surface ='',i7)')
     &                  nsf(j2), nsf(j1)

            ierr = ierr + 1
            return

         end if

         goto 100

*-----------------------------------------------------------------------
*        find dt and db, geometry top and bottom along vertical, tpp(i).
*-----------------------------------------------------------------------

   70    ix = lsc(j1)
         if(kst(j1).gt.1) goto 80
         t1 = dotpro(scf(ix+1),scf(ix+1))
         t2 = dotpro(scf(ix+1),tpp(1))
         if(abs(t2**2-t1).gt.eps5*t1) goto 100
         t3 = sign(scf(ix+4)/sqrt(t1),t2)
         goto 90

   80    if(abs(abs(tpp(kst(j1)-1))-1.).gt.eps5) goto 100
         t3 = scf(ix+1)*tpp(kst(j1)-1)
   90    db = min(t3,db)
         dt = max(t3,dt)
  100    continue
         dh = .5*(dt+db)

*-----------------------------------------------------------------------
*        get trf(ltrf+7+i,k), perimeter vector in surface jsu
*        orthogonal to normal, trf(ltrf+4+i,k), and vertical, tpp(i).
*-----------------------------------------------------------------------

      do 270 jsu = 1, mxj

         if(ksu(jsu).le.0) goto 270
         k = ksu(jsu)

         call crspro(tpp(1),trf(5,k),trf(8,k))

         if(abs(dotpro(trf(8,k),trf(8,k))-1.).gt.eps5)
     &      goto 360

*-----------------------------------------------------------------------
*        set 3rd vector in jsu, trf(ltrf+10+i,k), to vertical, tpp(i).
*        find tpp(i+3), midpoint in jsu between top and bottom.
*-----------------------------------------------------------------------

      do 110 i = 1, 3
         trf(10+i,k) = tpp(i)
  110    tpp(i+3) = trf(4+i,k)*trf(4,k)+dh*tpp(i)

*-----------------------------------------------------------------------
*        find trf(ltrf+3,m), the distance of surface m from point
*        tpp(i+3) in perimeter direction trf(ltrf+7+i,k).
*-----------------------------------------------------------------------

         n = 1
      do 120 m = mx, mxtr
         trf(3,m) = huge
         t = dotpro(trf(5,m),trf(8,k))
         if(abs(t).ge.eps5) trf(3,m) = (trf(4,m)
     &                     - dotpro(trf(5,m),tpp(4)))/t
  120    if(trf(3,m).lt.trf(3,n)) n = m

         if(trf(3,n).eq.huge) then

            write(io,'(/''** ERROR : from s.intprb, ''/
     &                ''   not all periodic boundaries specified.'')')

            ierr = ierr + 1
            return

         end if

*-----------------------------------------------------------------------
*        move along trf(ltrf+7+i,k) perimeter direction in surface jsu.
*        gpblcm(i) is the center of each facet.
*        check inward/outward jsu surface normal, gpblcm(3+i).
*        set trf(ltrf+2,k) = +1/-1 for jsu normal in/out of cell.
*-----------------------------------------------------------------------

         d2 = huge
  130    d1 = trf(3,n)
         trf(3,n) = huge

      do 140 m = mx, mxtr
  140    if(trf(3,m).lt.trf(3,n)) n = m
         if(trf(3,n).eq.huge) goto 250
         dx = .5*(d1+trf(3,n))

         xxx = tpp(4)+dx*trf(8,k)  !FURUTA
         yyy = tpp(5)+dx*trf(9,k)  !FURUTA
         zzz = tpp(6)+dx*trf(10,k) !FURUTA
         uuu = trf(5,k)            !FURUTA
         vvv = trf(6,k)            !FURUTA
         www = trf(7,k)            !FURUTA

         jt = 0

      do 180 i1 = 1, mxa
         if(junf.ne.0) jt = abs(jun(i1))
         if(jt.ne.0) goto 180

      do 160 jk = abs(lca(i1)), abs(lca(i1+1))-1
  160    if(abs(lja(jk)).eq.jsu) goto 170
         goto 180

  170    call chkcll(i1,2,js)

         if(js.eq.0) goto 190
  180    continue
         goto 130

 190     uuu = -uuu !FURUTA
         vvv = -vvv !FURUTA
         www = -www !FURUTA
         jt = 0

      do 230 i2 = 1, mxa
         if(junf.ne.0) jt = abs(jun(i2))
         if(jt.ne.0) goto 230

      do 210 jk = abs(lca(i2)), abs(lca(i2+1))-1
  210    if(abs(lja(jk)).eq.jsu) goto 220
         goto 230

  220    call chkcll(i2,2,js)

         if(js.eq.0) goto 240
  230 continue
         goto 130
  240    iz = 0


         im1 = idmg(i1)
         im2 = idmg(i2)

         if(im1.ne.0.and.im2.eq.0) iz = 1.
         if(im1.eq.0.and.im2.ne.0) iz = -1.

         if(iz.eq.0) goto 130
         if(d2.eq.huge) d2 = d1
         d3 = trf(3,n)
         trf(2,k) = iz
         jz = max(-iz*i1,iz*i2)
         trf(17,k) = jz
         if(junf.eq.0) goto 130

         if(mfl(1,jz).ne.0) then

            write(io,'(/''** ERROR : from s.intprb, ''/
     &                  ''   zero importance cell= '',i7,
     &                  ''  cannot have fill.'')')
     &                  ncl(jz)

            ierr = ierr + 1
            return

         end if

         goto 130

  250    if(d2.eq.huge) then

            write(io,'(/''** ERROR : from s.intprb, ''/
     &                  ''   periodic surface = '',i7,
     &                  ''  must bound zero-importance cell.'')')
     &                  nsf(jsu)

            ierr = ierr + 1
            return

         end if

*-----------------------------------------------------------------------
*        trf(ltrf+13+i,k) is center of facet jsu.
*-----------------------------------------------------------------------

      do 260 i = 1, 3
  260    trf(13+i,k)=tpp(i+3)+.5*(d2+d3)*trf(7+i,k)
  270    continue

*-----------------------------------------------------------------------
*        convert trf array to translation/rotation matrix, r, where the
*        matrices r*trf(m1)=trf(m2).  thus r=trf(m2)*transpose(trf(m1)).
*        ss term forces m1 and m2 to be headed in and out of cell.
*-----------------------------------------------------------------------

      do 320 m1 = mx, mxtr
         m2 = ksu(nint(trf(1,m1)))
         if(m2.gt.m1) goto 320
         t1 = trf(17,m1)
         trf(17,m1) = trf(17,m2)
         trf(17,m2) = t1
         is = nint(-trf(2,m1)*trf(2,m2))

      do 290 i = 1, 3
      do 290 j = 1, 3
         tpp(3*i+j) = 0.

      do 280 k = 1, 3
         ss = max(2*k-5,is)
  280    tpp(3*i+j) = tpp(3*i+j)
     &              + ss*trf(1+3*k+j,m2)*trf(1+3*k+i,m1)
  290    if(abs(tpp(3*i+j)).lt.eps6) tpp(3*i+j) = 0.

      do 310 i = 1, 3
      do 300 j = 1, 3
         trf(3*i+j+1,m1) = tpp(3*i+j)
  300    trf(3*i+j+1,m2) = tpp(3*j+i)
         trf(1+i,m1) = trf(13+i,m2)
     &                    - dotpro(trf(3*i+2,m2),
     &                             trf(14,m1))
         trf(1+i,m2) = trf(13+i,m1)
     &                    - dotpro(trf(3*i+2,m1),
     &                             trf(14,m2))
  310    continue
  320    continue

*-----------------------------------------------------------------------
*        print pereiodic boundary transformation/rotation relations.
*-----------------------------------------------------------------------


         write(io,330)
  330    format('*** periodic boundary conditions'//
     &       13h surface pair,5x,31htransformation/rotation matrix./)

      do 340 j = 1, mxj
  340    if(ksu(j).gt.0) write(io,350) nsf(j),
     &              nsf(nint(trf(1,ksu(j)))),
     &              (trf(i,ksu(j)),i=2,13)
  350    format(i7,i6,1p3e14.6,3x,3e14.6/11x,3e14.6,3x,3e14.6)
         return

*-----------------------------------------------------------------------
*        error exit.
*-----------------------------------------------------------------------

  360    continue

            write(io,'(/''** ERROR : from s.intprb, ''/
     &                  ''   no common periodic boundary'',
     &                  '' intersection vector.'')')

            ierr = ierr + 1

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function dotpro(v1,v2)
*                                                                      *
*       calculate the dot product = v1 . v2                            *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension v1(3),v2(3)

*-----------------------------------------------------------------------

      dotpro = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)

      return
      end


************************************************************************
*                                                                      *
      subroutine refsuf(cs)
*                                                                      *
*       reflect the particle trajectory at surface jsu.                *
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      use LAFDATAMOD !FURUTA20201127

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

         l = lev
         if(levp.eq.l) goto 20
         udt(1,l) = xxx
         udt(2,l) = yyy
         udt(3,l) = zzz
         lev = levp

         xxx = udt(1,lev) !FURUTA
         yyy = udt(2,lev) !FURUTA
         zzz = udt(3,lev) !FURUTA
         uuu = udt(4,lev) !FURUTA
         vvv = udt(5,lev) !FURUTA
         www = udt(6,lev) !FURUTA

   20    cs = angl(int(udt(7,lev))) !FURUTA20230706

*-----------------------------------------------------------------------
*        white boundary (cosine distribution scatter off surface).
*-----------------------------------------------------------------------

         if(ksu(jsu).ne.-2) goto 22
         c = sign(sqrt(rang()),-cs)
         call dtcos(c,ang(1),uuu,0,irdm)
         goto 24

*-----------------------------------------------------------------------
*        reflecting boundary (specular reflection off surface).
*-----------------------------------------------------------------------

   22    uuu = uuu-2.*cs*ang(1)
         vvv = vvv-2.*cs*ang(2)
         www = www-2.*cs*ang(3)

*-----------------------------------------------------------------------
*        transfer reflection to other levels.
*-----------------------------------------------------------------------

   24    if(lev.ne.0) call higlev(uuu)
         if(lev.eq.l) goto 40

      do 30 lev = levp,l-1
         ic = udt(7,lev)
         udt(4,lev) = uuu
         udt(5,lev) = vvv
         udt(6,lev) = www
         j = -mfl(1,ic)
         if(j.lt.0) m = mfl(3,ic)
         if(j.gt.0) m = laf(3,j+3+int(udt(8,lev))              !FURUTA20201127
     &                - laf(1,j+1)+laf(1,j+2)                  !FURUTA20201127
     &                * (int(udt(9,lev))-laf(2,j+1)+laf(2,j+2) !FURUTA20201127
     &                * (int(udt(10,lev))-laf(3,j+1))))        !FURUTA20201127
         if(m.eq.0) goto 30
         uuu = udt(4,lev)*trf(5,m)+udt(5,lev)*trf(8,m)
     &       + udt(6,lev)*trf(11,m)
         vvv = udt(4,lev)*trf(6,m)+udt(5,lev)*trf(9,m)
     &       + udt(6,lev)*trf(12,m)
         www = udt(4,lev)*trf(7,m)+udt(5,lev)*trf(10,m)
     &       + udt(6,lev)*trf(13,m)
   30 continue

         lev = l
         xxx = udt(1,l)
         yyy = udt(2,l)
         zzz = udt(3,l)
   40    if(lca(icl).lt.0) call chkcll(icl,3,j)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine dtcos(c,a,d,l,irdm)
*                                                                      *
*       sample a direction d at an angle arccos(c) from axis a and     *
*       at an azimuthal angle sampled uniformly.                       *
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

      data ir /50/

*-----------------------------------------------------------------------

      character h*16
      dimension a(3),d(3)

*-----------------------------------------------------------------------

         if(abs(c).ge.1.) goto 30
   10    t1 = 2.*rang()-1.
         t2 = 2.*rang()-1.
         r = t1**2+t2**2
         if(r.gt.1.) goto 10
         r = sqrt((1.-c**2)/r)
         t1 = t1*r
         t2 = t2*r
         if(abs(a(3)).gt..9) goto 20
         s = sqrt(a(1)**2+a(2)**2)
         t = 1./s
         d(1) = a(1)*c+(t1*a(1)*a(3)-t2*a(2))*t
         d(2) = a(2)*c+(t1*a(2)*a(3)+t2*a(1))*t
         d(3) = a(3)*c-t1*s
         if(l.ne.0)call higlev(d)

*-----------------------------------------------------------------------
*        renormalize every 50 calls to prevent error buildup.
*-----------------------------------------------------------------------

         s = 1./sqrt(d(1)**2+d(2)**2+d(3)**2)
         d(1) = d(1)*s
         d(2) = d(2)*s
         d(3) = d(3)*s
         if(l.ne.0) call higlev(d)
         return

*-----------------------------------------------------------------------
*        special handling for the case of exceptionally large a(3).
*-----------------------------------------------------------------------

   20    s = sqrt(a(1)**2+a(3)**2)
         t = 1./s
         d(1) = a(1)*c+(t1*a(1)*a(2)+t2*a(3))*t
         d(2) = a(2)*c-t1*s
         d(3) = a(3)*c+(t1*a(3)*a(2)-t2*a(1))*t
         if(l.ne.0) call higlev(d)
         return

*-----------------------------------------------------------------------
*        special handling for abs(c)=1. or more.
*-----------------------------------------------------------------------

   30    if(abs(c).gt.1.d0) c = sign(1.d0,c)
         d(1) = c*a(1)
         d(2) = c*a(2)
         d(3) = c*a(3)
         if(l.ne.0) call higlev(d)

*-----------------------------------------------------------------------

      return
      end


