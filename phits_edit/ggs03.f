************************************************************************
*                                                                      *
      subroutine intlat(i1)
*                                                                      *
*       move into the lattice element beyond surface lja(llja+i1).     *
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************

      use LAFDATAMOD !FURUTA20201127
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      dimension ij(8),kj(8,2),ii(3)
      data ij/1,-1,1,-1,1,-1,1,-1/,kj/1,1,2,2,3,3,0,0,1,1,2,2,2,2,3,3/

*-----------------------------------------------------------------------

         iap = icl
         n = i1-lca(icl)+1
         l = lat(1,icl)
         i = ij(n)
         jsu = abs(lja(i1+i))
         k = kj(n,l)
         m = lat(2,icl)

         xxx = xxx-i*vcl(1,k,m)
         yyy = yyy-i*vcl(2,k,m)
         zzz = zzz-i*vcl(3,k,m)
         select case(k)                        !FURUTA
         case(1)                               !FURUTA
           iii=iii+i                           !FURUTA
           iik=iii                             !FURUTA
         case(2)                               !FURUTA
           jjj=jjj+i                           !FURUTA
           iik=jjj                             !FURUTA
         case(3)                               !FURUTA
           kkk=kkk+i                           !FURUTA
           iik=kkk                             !FURUTA
         case default                          !FURUTA
           write(*,*)'ERROR in INTLAT: SELECT' !FURUTA
           stop                                !FURUTA
         end select                            !FURUTA

         if(mfl(1,icl).ge.0) goto 10
         j = -mfl(1,icl)
         if(iik.lt.laf(k,j+1).or.iik.ge.laf(k,j+1)+laf(k,j+2)) !FURUTA20201127
     &        goto 20
   10    if(l.eq.1) return
         if(n.ne.5.and.n.ne.6) return

         xxx = xxx+i*vcl(1,1,m)
         yyy = yyy+i*vcl(2,1,m)
         zzz = zzz+i*vcl(3,1,m)
         iii=iii-i

         if(mfl(1,icl).ge.0) return
         j = -mfl(1,icl)
         if(iii.ge.laf(1,j+1).and.iii.lt.laf(1,j+1)+laf(1,j+2)) return !FURUTA20201127
   20    iap = -12345

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine dwnlev(ic)
*                                                                      *
*       transform coordinates to the next level down,                  *
*       into the universe filling cell ic.                             *
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************

      use LAFDATAMOD !FURUTA20201127
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

         udt(1,lev) = xxx
         udt(2,lev) = yyy
         udt(3,lev) = zzz
         udt(4,lev) = uuu
         udt(5,lev) = vvv
         udt(6,lev) = www
         udt(7,lev) = ic
         udt(8,lev) = iii
         udt(9,lev) = jjj
         udt(10,lev) = kkk

         j = -mfl(1,ic)
         if(j.lt.0) m = mfl(3,ic)
         if(j.gt.0) m = laf(3,j+3+iii-laf(1,j+1)      !FURUTA20201127
     &                + laf(1,j+2)*(jjj-laf(2,j+1)    !FURUTA20201127
     &                + laf(2,j+2)*(kkk-laf(3,j+1)))) !FURUTA20201127
         if(m.eq.0) goto 10

         xxx = udt(1,lev)*trf(5,m)
     &       + udt(2,lev)*trf(8,m)
     &       + udt(3,lev)*trf(11,m)+trf(2,m)
         yyy = udt(1,lev)*trf(6,m)
     &       + udt(2,lev)*trf(9,m)
     &       + udt(3,lev)*trf(12,m)+trf(3,m)
         zzz = udt(1,lev)*trf(7,m)
     &       + udt(2,lev)*trf(10,m)
     &       + udt(3,lev)*trf(13,m)+trf(4,m)
         uuu = udt(4,lev)*trf(5,m)
     &       + udt(5,lev)*trf(8,m)
     &       + udt(6,lev)*trf(11,m)
         vvv = udt(4,lev)*trf(6,m)
     &       + udt(5,lev)*trf(9,m)
     &       + udt(6,lev)*trf(12,m)
         www = udt(4,lev)*trf(7,m)
     &       + udt(5,lev)*trf(10,m)
     &       + udt(6,lev)*trf(13,m)
   10    lev = lev+1

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine higlev(u)
*                                                                      *
*       calculate direction cosines in all higher levels.              *
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************

      use LAFDATAMOD !FURUTA20201127
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      dimension u(3)

*-----------------------------------------------------------------------

         udt(4,lev) = u(1)
         udt(5,lev) = u(2)
         udt(6,lev) = u(3)

      do 20 lz = 0, lev-1
         l = lev-1-lz
         j = -mfl(1,int(udt(7,l)))
         if(j.lt.0) m = mfl(3,int(udt(7,l)))
         if(j.gt.0) m = laf(3,j+int(3+udt(8,l)-laf(1,j+1)    !FURUTA20201127
     &                + laf(1,j+2)*(udt(9,l)-laf(2,j+1)      !FURUTA20201127
     &                + laf(2,j+2)*(udt(10,l)-laf(3,j+1))))) !FURUTA20201127
         if(m.eq.0) goto 10
         udt(4,l) = udt(4,l+1)*trf(5,m)
     &            + udt(5,l+1)*trf(6,m)
     &            + udt(6,l+1)*trf(7,m)
         udt(5,l) = udt(4,l+1)*trf(8,m)
     &            + udt(5,l+1)*trf(9,m)
     &            + udt(6,l+1)*trf(10,m)
         udt(6,l) = udt(4,l+1)*trf(11,m)
     &            + udt(5,l+1)*trf(12,m)
     &            + udt(6,l+1)*trf(13,m)
         goto 20

   10    udt(4,l) = udt(4,l+1)
         udt(5,l) = udt(5,l+1)
         udt(6,l) = udt(6,l+1)

   20 continue

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine setsuf(ih,io,ierr)
*                                                                      *
*       set up icl, lev, and udt for ggmsor                            *
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                              *
*                                                                      *
************************************************************************
      use TETRAMOD, only: tetrafnd, ielem2icl
      use LAFDATAMOD !FURUTA20201127
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------
      common /regdc/  idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------

      common /igsherr/ igsher, icl01, icl02

*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk
      integer n,ks,ls,itet
      real(8) xx0(6)
*-----------------------------------------------------------------------

         ierr = 0

*-----------------------------------------------------------------------

         jt = 0
         mm = 0
         if(jsu.ne.0) mm = 2
         ju = 0
         if( ih .ne. 0 ) then
            js = 1
            icl0 = icl
            icl = ih
            if(junf.ne.0) jt = abs(jun(icl))
            if(jt.eq.ju) call chkcll(icl,mm,js)

            if(js.eq.0) then

               if( igsher .ne. 0 ) then
                     jss = 1
                  do icll = mxa, 1, -1
                     if( icll .eq. icl ) cycle
                     if(junf.ne.0) jt = abs(jun(icll))
                     if(jt.eq.ju) call chkcll(icll,mm,jss)
                     if(jss.eq.0) then
                        jss = 1
                        icl01 = idrg(icl)
                        icl02 = idrg(icll)
                        ierr = 2
                        return
                     end if
                  end do
               end if

               goto 50
            end if

            icl = icl0
         end if
*-----------------------------------------------------------------------
*        look for the cell in the list of previously found cells.
*        jt=cell icl universe;  ju=present location universe.
*-----------------------------------------------------------------------

   10    js = 1
      do 20 i = 1, nlse
         icl = lse(i)
         if(junf.ne.0) jt = abs(jun(icl))
         if(jt.eq.ju) call chkcll(icl,mm,js)

         if(js.eq.0) then

               if( igsher .ne. 0 ) then
                     jss = 1
                  do icll = mxa, 1, -1
                     if( icll .eq. icl ) cycle
                     if(junf.ne.0) jt = abs(jun(icll))
                     if(jt.eq.ju) call chkcll(icll,mm,jss)
                     if(jss.eq.0) then
                        jss = 1
                        icl01 = idrg(icl)
                        icl02 = idrg(icll)
                        ierr = 2
                        return
                     end if
                  end do
               end if

            goto 50
         end if

   20    continue

*-----------------------------------------------------------------------
*        look for the cell among all cells of the universe.
*-----------------------------------------------------------------------

      do 30 icl = 1, mxa
         if(junf.ne.0) jt = abs(jun(icl))
         if(jt.eq.ju) call chkcll(icl,mm,js)

         if(js.eq.0) then

               if( igsher .ne. 0 ) then
                     jss = 1
                  do icll = mxa, icl+1, -1
                     if(junf.ne.0) jt = abs(jun(icll))
                     if(jt.eq.ju) call chkcll(icll,mm,jss)
                     if(jss.eq.0) then
                        jss = 1
                        icl01 = idrg(icl)
                        icl02 = idrg(icll)
                        ierr = 2
                        return
                     end if
                  end do
               end if

            goto 40
         end if

   30    continue

*-----------------------------------------------------------------------
*     call expirx(0,'setsuf','source is not in any cell.')
*-----------------------------------------------------------------------

         ierr = 1
         return

*-----------------------------------------------------------------------
*        add the cell to the list of found cells.
*-----------------------------------------------------------------------

   40    nlse = nlse+1
         lse(nlse) = icl

*-----------------------------------------------------------------------
*        move into the universe, if any, that fills the cell.
*-----------------------------------------------------------------------

   50    if(junf.eq.0) return
   60    jsu = -abs(jsu)

         if(lat(1,icl).ne.0) then
          if(abs(lat(1,icl)).ne.3)then
           call fndlat(io,icl,ierr)
cFURUTA20150714 TETRA---------------------------------------------------
          else
           if(lat(1,icl).eq.3)then
            do itet=1,nlat3
             if(itetcl(itet).eq.icl)then
              kkk=10000+itet
              exit
             endif
            enddo
           endif
           if(iii.eq.0)then !FURUTA20160607
            call tetrabox(io,kkk,xx0,ierr)
            call tetrafnd(xx0,xxx,yyy,zzz,uuu,vvv,www,coincd,itet,
     &           iii,icl,ierr)
           else
            if(iii.gt.0)then
             icl=ielem2icl(iii)
            endif
           endif            !FURUTA20160607

*-----------------------------------------------------------------------
          endif
          if( ierr .ne. 0 ) return
         end if

         jsu = abs(jsu)
         if(mfl(1,icl).eq.0) return
         j = -mfl(1,icl)
         if(j.gt.0) goto 70
         ju = -j

         call dwnlev(icl)

         if(mfl(2,icl).eq.0) goto 10
         icl = mfl(2,icl)
         goto 60

   70    n = 3+iii-laf(1,j+1)+laf(1,j+2)*(jjj-laf(2,j+1) !FURUTA20201127
     &     + laf(2,j+2)*(kkk-laf(3,j+1)))                !FURUTA20201127
         if(laf(1,j+n).eq.abs(jun(icl))) return     !FURUTA20201127
         ju = laf(1,j+n)                                 !FURUTA20201127

         call dwnlev(icl)

         if(laf(2,j+n).eq.0) goto 10 !FURUTA20201127
         icl = laf(2,j+n)            !FURUTA20201127
         goto 60

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine fndlat(io,ic,ierr)
*                                                                      *
*       find the element in lattice ic that contains xxx,yyy,zzz and   *
*       adjust xxx,yyy,zzz to be relative to that element.             *
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************

      use LAFDATAMOD !FURUTA20201127
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      dimension ii(3),jx(14,2)
      data jx/1,1,2,2,3,3,1,1,2,2,2,2,3,3,
     &  1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1/

      dimension jj(3)

*-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
*        find the element.
*-----------------------------------------------------------------------

         m = lat(2,ic)
         ii(1)=iii !FURUTA
         ii(2)=0   !FURUTA
         ii(3)=0   !FURUTA
         i1 = (abs(lca(ic+1))-lca(ic))/2-lat(1,ic)+1

      do 10 i = 1, i1
         a = xxx*vcl(1,i+4,m)
     &     + yyy*vcl(2,i+4,m)
     &     + zzz*vcl(3,i+4,m)-vcl(i,4,m)
   10    ii(i) = nint(a)

*-----------------------------------------------------------------------
*        adjust xxx,yyy,zzz.
*-----------------------------------------------------------------------

 20      xxx = xxx-ii(1)*vcl(1,1,m)                         !FURUTA
     &             - ii(2)*vcl(1,2,m)-ii(3)*vcl(1,3,m) !FURUTA
         yyy = yyy-ii(1)*vcl(2,1,m)                         !FURUTA
     &             - ii(2)*vcl(2,2,m)-ii(3)*vcl(2,3,m) !FURUTA
         zzz = zzz-ii(1)*vcl(3,1,m)                         !FURUTA
     &             - ii(2)*vcl(3,2,m)-ii(3)*vcl(3,3,m) !FURUTA

*-----------------------------------------------------------------------
*        correct for coincident surfaces and hexagonal prism lattices.
*        check if new location completely inside lattice element.
*-----------------------------------------------------------------------

         jc = 0
   30    jc = jc+1
         if(jc.gt.10) goto 60

      do 50 jk = lca(ic), abs(lca(ic+1))-1
         j = abs(lja(jk))
         k = kst(j)
         l = lsc(j)
         select case(k)                                !FURUTA
         case(1)                                       !FURUTA
           t4 = scf(l+1)*xxx+scf(l+2)*yyy              !FURUTA
     &          + scf(l+3)*zzz-scf(l+4)                !FURUTA
           t5 = scf(l+1)*uuu+scf(l+2)*vvv+scf(l+3)*www !FURUTA
         case(2)                                       !FURUTA
           t4 = xxx-scf(l+1)                           !FURUTA
           t5 = uuu                                    !FURUTA
         case(3)                                       !FURUTA
           t4 = yyy-scf(l+1)                           !FURUTA
           t5 = vvv                                    !FURUTA
         case(4)                                       !FURUTA
           t4 = zzz-scf(l+1)                           !FURUTA
           t5 = www                                    !FURUTA
         case default                                  !FURUTA
           write(*,*)'ERROR in fndlat: SELECT'         !FURUTA
           stop                                        !FURUTA
         end select                                    !FURUTA
         if(jsu.ge.0) goto 40
         if(ksc(-jsu).ne.ksc(j)) goto 40
         if(abs(t4).le.coincd*abs(t5)) t4 = t5
   40    if(lja(jk)*t4.gt.0.) goto 50

*-----------------------------------------------------------------------
*        location coincident with or beyond surface j, index ix.
*-----------------------------------------------------------------------


         j1 = jk-lca(ic)+6*lat(1,ic)-5
         ix = jx(j1,1)
         iy = jx(j1,2)
         ii(ix) = ii(ix)+iy

         xxx = xxx-iy*vcl(1,ix,m)
         yyy = yyy-iy*vcl(2,ix,m)
         zzz = zzz-iy*vcl(3,ix,m)

*-----------------------------------------------------------------------
*        if on or beyond hex side 3, increment sides 1 and 2.
*-----------------------------------------------------------------------

         if(ix.ne.2.or.j1.lt.11) goto 30
         ii(1) = ii(1)-iy !FURUTA
         xxx = xxx+iy*vcl(1,1,m)
         yyy = yyy+iy*vcl(2,1,m)
         zzz = zzz+iy*vcl(3,1,m)
         goto 30
   50    continue

*-----------------------------------------------------------------------
*        final check that hexagonal element selected exists.
*-----------------------------------------------------------------------

         if(mfl(1,ic).ge.0)then !FURUTA
           iii=ii(1)                 !FURUTA
           jjj=ii(2)                 !FURUTA
           kkk=ii(3)                 !FURUTA
           return                    !FURUTA
         endif                       !FURUTA
         j = -mfl(1,ic)
         if(ii(1).lt.laf(1,j+1).or.ii(1).ge.laf(1,j+1)+laf(1,j+2))
     &       go to 60
         if(ii(2).lt.laf(2,j+1).or.ii(2).ge.laf(2,j+1)+laf(2,j+2))
     &       go to 60
         if(ii(3).lt.laf(3,j+1).or.ii(3).ge.laf(3,j+1)+laf(3,j+2))
     &       go to 60
c-----------------------------------------------------------------
         iii=ii(1) !FURUTA
         jjj=ii(2) !FURUTA
         kkk=ii(3) !FURUTA
         return

*-----------------------------------------------------------------------
*        cannot find lattice location.
*-----------------------------------------------------------------------

   60    kdb = 1

            ierr = 1


*-----------------------------------------------------------------------
            iii=ii(1)           !FURUTA
            jjj=ii(2)           !FURUTA
            kkk=ii(3)           !FURUTA
      return
      end

************************************************************************
*                                                                      *
      subroutine upposs(x,u,l,d,v,jf)
*                                                                      *
*       calculate the position and direction vectors, x and u, at      *
*       levels 0 thru l-1 for a particle travelling a distance d at    *
*       velocity v.  if jf.ne.0 an additional call to higlev is needed.*
*                                                                      *
*       Last modified by K.Niita on 2009/10/01                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      dimension x(3),u(3)

*-----------------------------------------------------------------------

      do 10 i = 1, 3
   10    x(i) = x(i)+d*u(i)

      if(l.eq.0) return

      call higlev(u)

      do 20 j = 0, l-1
      do 20 i = 1, 3
   20    udt(i,j) = udt(i,j)+d*udt(i+3,j)

*-----------------------------------------------------------------------
*     sometimes higlev must be called again to update direction cosine
*     as when uuu is used by knoelc, kxrelc, and brmelc because higlev
*     is previously called with uold.
*-----------------------------------------------------------------------

      if(jf.ne.0) call higlev(uuu)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine trnsxx(x1,y1,z1,x2,y2,z2,m)
*                                                                      *
*       transform the xyz coordinate                                   *
*                                                                      *
*       modified by K.Niita on 2003/08/13                              *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      if( m .eq. 0 ) then

         x2 = x1
         y2 = y1
         z2 = z1

      else

         x2 = x1 * trf( 5,m) +
     &        y1 * trf( 8,m) +
     &        z1 * trf(11,m) + trf(2,m)

         y2 = x1 * trf( 6,m) +
     &        y1 * trf( 9,m) +
     &        z1 * trf(12,m) + trf(3,m)

         z2 = x1 * trf( 7,m) +
     &        y1 * trf(10,m) +
     &        z1 * trf(13,m) + trf(4,m)

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine trnsxv(x1,y1,z1,x2,y2,z2,m)
*                                                                      *
*       inverse transform the xyz coordinate                           *
*                                                                      *
*       modified by K.Niita on 2003/08/15                              *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      if( m .eq. 0 ) then

         x2 = x1
         y2 = y1
         z2 = z1

      else

         x2 = ( x1 - trf(2,m) ) * trf( 5,m) +
     &        ( y1 - trf(3,m) ) * trf( 6,m) +
     &        ( z1 - trf(4,m) ) * trf( 7,m)

         y2 = ( x1 - trf(2,m) ) * trf( 8,m) +
     &        ( y1 - trf(3,m) ) * trf( 9,m) +
     &        ( z1 - trf(4,m) ) * trf(10,m)

         z2 = ( x1 - trf(2,m) ) * trf(11,m) +
     &        ( y1 - trf(3,m) ) * trf(12,m) +
     &        ( z1 - trf(4,m) ) * trf(13,m)

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine trnsuu(x1,y1,z1,x2,y2,z2,m)
*                                                                      *
*       transform the xyz direction                                    *
*                                                                      *
*       modified by K.Niita on 2003/09/03                              *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      if( m .eq. 0 ) then

         x2 = x1
         y2 = y1
         z2 = z1

      else

         x2 = x1 * trf( 5,m) +
     &        y1 * trf( 8,m) +
     &        z1 * trf(11,m)

         y2 = x1 * trf( 6,m) +
     &        y1 * trf( 9,m) +
     &        z1 * trf(12,m)

         z2 = x1 * trf( 7,m) +
     &        y1 * trf(10,m) +
     &        z1 * trf(13,m)

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine trnsuv(x1,y1,z1,x2,y2,z2,m)
*                                                                      *
*       inverse transform the xyz direction                            *
*                                                                      *
*       modified by K.Niita on 2003/09/03                              *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      if( m .eq. 0 ) then

         x2 = x1
         y2 = y1
         z2 = z1

      else

         x2 = x1 * trf( 5,m) +
     &        y1 * trf( 6,m) +
     &        z1 * trf( 7,m)

         y2 = x1 * trf( 8,m) +
     &        y1 * trf( 9,m) +
     &        z1 * trf(10,m)

         z2 = x1 * trf(11,m) +
     &        y1 * trf(12,m) +
     &        z1 * trf(13,m)

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gomsort(x2,y2,z2,u2,v2,w2,
     &                   nmed,iblz,mark,markp,ici,m)
*                                                                      *
*       check position after inverse transform the xyz coordinate      *
*                                                                      *
*       modified by K.Niita on 2003/08/13                              *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      if( m .eq. 0 ) then

            call gomsor(x2,y2,z2,u2,v2,w2,
     &                  nmed,iblz,mark,markp,ici)

      else

         x1 = ( x2 - trf(2,m) ) * trf( 5,m) +
     &        ( y2 - trf(3,m) ) * trf( 6,m) +
     &        ( z2 - trf(4,m) ) * trf( 7,m)

         y1 = ( x2 - trf(2,m) ) * trf( 8,m) +
     &        ( y2 - trf(3,m) ) * trf( 9,m) +
     &        ( z2 - trf(4,m) ) * trf(10,m)

         z1 = ( x2 - trf(2,m) ) * trf(11,m) +
     &        ( y2 - trf(3,m) ) * trf(12,m) +
     &        ( z2 - trf(4,m) ) * trf(13,m)

         u1 = u2 * trf( 5,m) +
     &        v2 * trf( 6,m) +
     &        w2 * trf( 7,m)

         v1 = u2 * trf( 8,m) +
     &        v2 * trf( 9,m) +
     &        w2 * trf(10,m)

         w1 = u2 * trf(11,m) +
     &        v2 * trf(12,m) +
     &        w2 * trf(13,m)

            call gomsor(x1,y1,z1,u1,v1,w1,
     &                  nmed,iblz,mark,markp,ici)

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gomprpt(mz,x2,y2,z2,xc2,yc2,zc2,u2,v2,w2,
     &                   nmed,iblz,mark,markp,m)
*                                                                      *
*       check boundary after inverse transform the xyz coordinate      *
*                                                                      *
*            modified by K.Niita on 2003/08/13                         *
*       Last modified by K.Niita on 2020/04/09                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      if( m .eq. 0 ) then

            call gomppg(mz,x2,y2,z2,xc2,yc2,zc2,u2,v2,w2,
     &                  nmed,iblz,mark,markp)

      else

         x1 = ( x2 - trf(2,m) ) * trf( 5,m) +
     &        ( y2 - trf(3,m) ) * trf( 6,m) +
     &        ( z2 - trf(4,m) ) * trf( 7,m)

         y1 = ( x2 - trf(2,m) ) * trf( 8,m) +
     &        ( y2 - trf(3,m) ) * trf( 9,m) +
     &        ( z2 - trf(4,m) ) * trf(10,m)

         z1 = ( x2 - trf(2,m) ) * trf(11,m) +
     &        ( y2 - trf(3,m) ) * trf(12,m) +
     &        ( z2 - trf(4,m) ) * trf(13,m)

         xc1 = ( xc2 - trf(2,m) ) * trf( 5,m) +
     &         ( yc2 - trf(3,m) ) * trf( 6,m) +
     &         ( zc2 - trf(4,m) ) * trf( 7,m)

         yc1 = ( xc2 - trf(2,m) ) * trf( 8,m) +
     &         ( yc2 - trf(3,m) ) * trf( 9,m) +
     &         ( zc2 - trf(4,m) ) * trf(10,m)

         zc1 = ( xc2 - trf(2,m) ) * trf(11,m) +
     &         ( yc2 - trf(3,m) ) * trf(12,m) +
     &         ( zc2 - trf(4,m) ) * trf(13,m)

         u1 = u2 * trf( 5,m) +
     &        v2 * trf( 6,m) +
     &        w2 * trf( 7,m)

         v1 = u2 * trf( 8,m) +
     &        v2 * trf( 9,m) +
     &        w2 * trf(10,m)

         w1 = u2 * trf(11,m) +
     &        v2 * trf(12,m) +
     &        w2 * trf(13,m)

            call gomppg(mz,x1,y1,z1,xc1,yc1,zc1,u1,v1,w1,
     &                  nmed,iblz,mark,markp)

         xc2 = xc1 * trf( 5,m) +
     &         yc1 * trf( 8,m) +
     &         zc1 * trf(11,m) + trf(2,m)

         yc2 = xc1 * trf( 6,m) +
     &         yc1 * trf( 9,m) +
     &         zc1 * trf(12,m) + trf(3,m)

         zc2 = xc1 * trf( 7,m) +
     &         yc1 * trf(10,m) +
     &         zc1 * trf(13,m) + trf(4,m)

      end if

*-----------------------------------------------------------------------

      return
      end

