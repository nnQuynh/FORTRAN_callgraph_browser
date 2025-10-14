************************************************************************
*                                                                      *
      subroutine gomsor(x,y,z,u,v,w,nmed,iblz,mark,markp,ici)
*                                                                      *
*       geometry control routine                                       *
*       modified by K.Niita on 2002/05/09                              *
*                                                                      *
*       ici = 0 ; first time from sorce                                *
*                 stop after not finding the cell                      *
*       ici =-1 ; same as 0 but return after error                     *
*       ici = 1 ; return with ierr=1 after not finding the cell        *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      common /ccggg/  icgg

*-----------------------------------------------------------------------
*     CG geometry
*-----------------------------------------------------------------------

      if( icgg .eq. 0 ) then

         call cggsor(x,y,z,u,v,w,nmed,iblz,mark,markp,ici)

*-----------------------------------------------------------------------
*     GG geometry
*-----------------------------------------------------------------------

      else

         call ggmsor(x,y,z,u,v,w,nmed,iblz,mark,markp,ici)

      end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine gomprp(icr,x,y,z,xc,yc,zc,u,v,w,
     &                  nmed,iblz,mark,markp)
*                                                                      *
*       geometry control routine                                       *
*       modified by K.Niita on 2001/11/29                              *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /ccggg/  icgg

*-----------------------------------------------------------------------
*     CG geometry
*-----------------------------------------------------------------------

      if( icgg .eq. 0 ) then

         call cggprp(icr,x,y,z,xc,yc,zc,u,v,w,
     &               nmed,iblz,mark,markp)

*-----------------------------------------------------------------------
*     GG geometry
*-----------------------------------------------------------------------

      else

         call ggmprp(icr,x,y,z,xc,yc,zc,u,v,w,
     &               nmed,iblz,mark,markp,icl0)

      end if


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gomppg(icr,x,y,z,xc,yc,zc,u,v,w,
     &                  nmed,iblz,mark,markp)
*                                                                      *
*       geometry control routine                                       *
*       modified by K.Niita on 2020/04/09                              *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /ccggg/  icgg

*-----------------------------------------------------------------------
*     CG geometry
*-----------------------------------------------------------------------

      if( icgg .eq. 0 ) then

         call cggprp(icr,x,y,z,xc,yc,zc,u,v,w,
     &               nmed,iblz,mark,markp)

*-----------------------------------------------------------------------
*     GG geometry
*-----------------------------------------------------------------------

      else

         call ggmpgp(icr,x,y,z,xc,yc,zc,u,v,w,
     &               nmed,iblz,mark,markp)

      end if


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gomdis(dist,x,y,z,u,v,w,
     &                  mark,markp,nmed,iblz)
*                                                                      *
*       get distance up to boundary                                    *
*       modified by K.Niita on 2002/05/02                              *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /ccggg/ icgg

*-----------------------------------------------------------------------
*     CG geometry
*-----------------------------------------------------------------------

      if( icgg .eq. 0 ) then

               call cggdis(dist,x,y,z,u,v,w,
     &                     nmed,iblz,mark,markp)

*-----------------------------------------------------------------------
*     GG geometry
*-----------------------------------------------------------------------

      else

               call ggmdis(dist,x,y,z,u,v,w,
     &                     nmed,iblz,mark,markp)

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gomnew(icr,dpr,x,y,z,xc,yc,zc,u,v,w,ec,
     &                  nmed,iblz,mark,markp)
*                                                                      *
*       geometry control routine                                       *
*       modified by K.Niita on 2002/05/02                              *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /ccggg/  icgg

*-----------------------------------------------------------------------

         iblz0 = iblz

*-----------------------------------------------------------------------
*     CG geometry
*-----------------------------------------------------------------------

      if( icgg .eq. 0 ) then

         call cggnew(icr,dpr,x,y,z,xc,yc,zc,u,v,w,ec,
     &               nmed,iblz,mark,markp)

*-----------------------------------------------------------------------
*     GG geometry
*-----------------------------------------------------------------------

      else

         call ggmnew(icr,dpr,x,y,z,xc,yc,zc,u,v,w,ec,
     &               nmed,iblz,mark,markp)

      end if


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gomupr(mark,markp,delt)
*                                                                      *
*       geometry control routine                                       *
*       modified by K.Niita on 2003/09/19                              *
*                                                                      *
*        markp = 1 for GG, 0 for CG                                    *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
*-----------------------------------------------------------------------

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /ccggg/  icgg

      common /paraj/ mstz(300), parz(300)
      real*8  parz
      integer mstz

*-----------------------------------------------------------------------

               xc(ibkxc+no,ipomp+1) = x(ibkx+no,ipomp+1) +
     &                                   u(ibku+no,ipomp+1) * delt
               yc(ibkyc+no,ipomp+1) = y(ibky+no,ipomp+1) +
     &                                   v(ibkv+no,ipomp+1) * delt
               zc(ibkzc+no,ipomp+1) = z(ibkz+no,ipomp+1) +
     &                                   w(ibkw+no,ipomp+1) * delt

*-----------------------------------------------------------------------

      if(mstz(85).eq.99)then ! debug output
         write(93,*)'-----------(take step here)-----------------'
         write(93,'(a10,f10.5)')'step=', delt
         write(93,'(6f10.5)') x(ibkx+no,ipomp+1),  y(ibky+no,ipomp+1),
     &                        z(ibkz+no,ipomp+1),
     $                        u(ibku+no,ipomp+1),  v(ibkv+no,ipomp+1),
     &                        w(ibkw+no,ipomp+1)
         write(93,'(6f10.5)')xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                       zc(ibkzc+no,ipomp+1),
     $                        u(ibku+no,ipomp+1),  v(ibkv+no,ipomp+1),
     &                        w(ibkw+no,ipomp+1)
      endif

*-----------------------------------------------------------------------
*        CG geometry
*-----------------------------------------------------------------------

         if( icgg .eq. 0 ) then

               mark  = 1
               markp = 0

*-----------------------------------------------------------------------
*        GG geometry
*-----------------------------------------------------------------------

         else

               mark  = 1
               markp = 1

               xxx = xxx + uuu * delt
               yyy = yyy + vvv * delt
               zzz = zzz + www * delt

               udt(1,lev) = xxx
               udt(2,lev) = yyy
               udt(3,lev) = zzz

            do l = 0, lev - 1

               udt(1,l) = udt(1,l) + udt(4,l) * delt
               udt(2,l) = udt(2,l) + udt(5,l) * delt
               udt(3,l) = udt(3,l) + udt(6,l) * delt

            end do

               jsu  = 0
               levp = 0

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gomupp(mark,markp,uus,vvs,wws)
*                                                                      *
*       geometry control routine                                       *
*       modified by K.Niita on 2003/09/19                              *
*                                                                      *
************************************************************************
      use LAFDATAMOD !FURUTA20201127
      use moddas_ggs !frtati20220905
*-----------------------------------------------------------------------

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      common /ccggg/  icgg

*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      integer iii0,jjj0,kkk0
*-----------------------------------------------------------------------
*     CG geometry
*-----------------------------------------------------------------------

      if( icgg .eq. 0 ) then

               mark  = 1
               markp = 0

*-----------------------------------------------------------------------
*     GG geometry
*-----------------------------------------------------------------------

      else

         if( mark  .eq. 0 .or.
     &       mark  .eq. 2 .or.
     &       markp .eq. 1 ) then

                  udt(4,0) = uus
                  udt(5,0) = vvs
                  udt(6,0) = wws
cFURUTA20150714 TETRA---------------------------------------------------
                  iii0=iii
                  jjj0=jjj
                  kkk0=kkk
*-----------------------------------------------------------------------
               do l = 1, lev

                  ic  = nint( udt(7,l-1) )
                  iii = udt(8,l-1)
                  jjj = udt(9,l-1)
                  kkk = udt(10,l-1)

                  j = -mfl(1,ic)
                  if(j.lt.0) m = mfl(3,ic)
                  if(j.gt.0) m = laf(3,j+3+iii-laf(1,j+1)+laf(1,j+2)*  !FURUTA20201127
     &                   (jjj-laf(2,j+1)+laf(2,j+2)*(kkk-laf(3,j+1)))) !FURUTA20201127

                  if( m .gt. 0 ) then

                     udt(4,l) = udt(4,l-1) * trf(5,m)
     &                        + udt(5,l-1) * trf(8,m)
     &                        + udt(6,l-1) * trf(11,m)

                     udt(5,l) = udt(4,l-1) * trf(6,m)
     &                        + udt(5,l-1) * trf(9,m)
     &                        + udt(6,l-1) * trf(12,m)

                     udt(6,l) = udt(4,l-1) * trf(7,m)
     &                        + udt(5,l-1) * trf(10,m)
     &                        + udt(6,l-1) * trf(13,m)
                  else

                     udt(4,l) = udt(4,l-1)
                     udt(5,l) = udt(5,l-1)
                     udt(6,l) = udt(6,l-1)

                  end if

               end do

                     uuu = udt(4,lev)
                     vvv = udt(5,lev)
                     www = udt(6,lev)

                     iii=iii0
                     jjj=jjj0
                     kkk=kkk0

         end if

      end if

*-----------------------------------------------------------------------

      return
      end


