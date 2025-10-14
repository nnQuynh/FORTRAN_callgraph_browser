      subroutine abend(ld)
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     92/06/22
c     time the module was last permanently updated:     08:39:50
c     programmer name:                                  j.t.west
c     module name:                                      maabend
c     current archiving level number:                   00001
c     current number of permanent updates:              00001
c     date of last access by librarian:                 93/03/15
c     dataset name:  x4s.scale4.master
c
c#######################################################################
      dimension ld(*)
c * * this routine called to get traceback for abnormal termination *
      call errtra
      return
      end
      subroutine abox(buf)
c * * this routine reads bpp and wpp body data and
c     converts it to box and wed body input formats. * * * *
      implicit real*8 (a-h,o-z)
      dimension buf(12), aij(3,3), bjk(3,3), cik(3,3),
     1          tr(3,3), ts(3,3), dx(3), cf(3)
      parameter ( pi = 3.14159 26535 89793 23846, radns = pi/180 )
      ue     = buf(7)*pi/180.0
      ve     = buf(8)*pi/180.0
      we     = buf(9)*pi/180.0
      cf(1)  = buf(1)
      cf(2)  = buf(3)
      cf(3)  = buf(5)
      dx(1)  = buf(2)-buf(1)
      dx(2)  = buf(4)-buf(3)
      dx(3)  = buf(6)-buf(5)
      do 20 i=1,3
        do 10 j=1,3
          aij(i,j) = 0.0
          bjk(i,j) = 0.0
          cik(i,j) = 0.0
          tr(i,j)  = 0.0
          ts(i,j)  = 0.0
  10      continue
        aij(i,i) = 1.0
        bjk(i,i) = 1.0
        cik(i,i) = 1.0
  20    continue
      if ( ue.ne.0.0 ) then
        aij(1,1) =  dcos(ue)
        aij(1,2) = -dsin(ue)
        aij(2,1) =  dsin(ue)
        aij(2,2) =  dcos(ue)
      end if
      if ( ve.ne.0.0 ) then
        bjk(2,2) =  dcos(ve)
        bjk(2,3) = -dsin(ve)
        bjk(3,2) =  dsin(ve)
        bjk(3,3) =  dcos(ve)
      end if
      if ( we.ne.0.0 ) then
        cik(1,1) =  dcos(we)
        cik(1,3) = -dsin(we)
        cik(3,1) =  dsin(we)
        cik(3,3) =  dcos(we)
      end if
      do 50 i=1,3
        do 40 j=1,3
          do 30 l=1,3
            ts(i,j) = ts(i,j)+bjk(i,l)*aij(l,j)
  30        continue
  40      continue
  50    continue
      do 80 i=1,3
        do 70 j=1,3
          do 60 l=1,3
            tr(i,j) = tr(i,j)+cik(i,l)*ts(l,j)
  60        continue
  70      continue
        buf(i) = cf(i)
  80    continue
      ii     = 3
      do 100 i=1,3
        do 90 j=1,3
          ii      = ii+1
          buf(ii) = tr(i,j)*dx(i)
  90      continue
 100    continue
      return
      end
      subroutine albert( f,ierr)
c * * this routine processes arb body data for storage in fpd array *
c * * it calculates unit normal vector and minimum distance to
c     origin for each plane containing a side of the arb.  * * * *
      implicit real*8 (a-h,o-z)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      dimension f(*),x(3,8),ix(4,6),v(3,4)
      do 10 i=1,8
      do 10 j=1,3
      k=3*(i-1)+j
      x(j,i)=f(k)
   10 f(k)=0.
c
c
      ipmax=0
      do  40 i=1,6
      k=0
      nside=i-1
      m=f(i+24)
      do  20  j=1,4
      ix(j,i)=m-(m/10)*10
      if( ix(j,i).ne.0) k=k+1
      if( ix(j,i).gt.ipmax) ipmax=ix(j,i)
   20 m=m/10
      if( k.eq.0 ) go to 50
      if( k.ge.3 ) go to 40
      write(ioe ,22) i,j,f(i+24)
   22 format(26h error in side description,2i10,f10.0)
      return
   40 continue
      nside=6
c
c  find  minimum distance between points
c
   50 dmin=1.0d+20
      ipmax=ipmax-1
      do  60 i=1,ipmax
      do  60 j=i,ipmax
      d=(x(1,i)-x(1,j+1))**2+(x(2,i)-x(2,j+1))**2+(x(3,i)-x(3,j+1))**2
      if((d.gt.0).and.(d.lt.dmin)) dmin=d
   60 continue
c
      dmin=dsqrt(dmin)
c
      ipmax=ipmax+1
      do   100   i=1,nside
      j1= ix( 3,i)
      do   62    j=2,4,2
      j2= ix( j,i)
      do   62    k=1,3
   62 v(k,j) = x(k,j1)- x(k,j2)
      a= v(2,2)*v(3,4) -  v(3,2)*v(2,4)
      b= v(3,2)*v(1,4) -  v(1,2)*v(3,4)
      c= v(1,2)*v(2,4) -  v(2,2)*v(1,4)
      d=-(a*x(1,j1)+ b*x(2,j1)+ c*x(3,j1) )
      eps=dsqrt(a*a + b*b + c*c)
      npl=0
      nmi=0
      do   80  j=1,ipmax
      ds=(a*x(1,j) + b*x(2,j)+ c*x(3,j) + d)/eps
      if(dabs(ds).lt.dmin*1.0d-6)  go to 80
      if(ds) 72,80,74
   72 nmi=nmi+1
      go to 80
   74 npl=npl+1
   80 continue
      if( (nmi.eq.0).and.(npl.gt.0)) go  to  90
      eps=-eps
      if(  (npl.eq.0).and.(nmi.gt.0))go  to  90
      ierr=ierr+1
   85 format(26h error in face description,3i10/1p,5e15.7)
      write(ioe ,85) i,nmi,npl,f(i+24),a,b,c,d
      return
   90 f( 4*i-3)=a/eps
      f( 4*i-2)=b/eps
      f( 4*i-1)=c/eps
      f( 4*i  )=d/eps
  100 continue
      f(25)=dmin
      f(26)=nside
      write( iot ,95) (f(i),i=1,26)
   95 format( 22x,1p,6e15.7 )
      do 30 i=1,8
      do 30 j=1,3
      k=3*(i-1)+j+26
  30  f(k)=x(j,i)
      do 45 i=1,6
      m=0
      do 42 j=1,4
  42  m=m+ix(j,i)*(10**(j-1))
      f(k+i)=dble(m)
  45  continue
      return
      end

*-----------------------------------------------------------------------

      subroutine argen(d,ld,ncmax,nbn,nbod,rby,ma,fpd,wlh,na)

*-----------------------------------------------------------------------

c   this routine determines the size of each array * * * * *
      implicit real*8 (a-h,o-z)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common/jomk/nbi,iri,jmk,iml,iop,imax

      dimension d(*),ld(*),ncmax(3),nbn(*),nbod(*),rby(*),ma(*),
     1          fpd(*),wlh(*),ixy(3)

      iz     = 3*(na-1)
      mr     = 1
      do 40 nd=1,3
        ndx     = nd
        rby(mr) = 0.0
        nmax    = ncmax(nd)
        do 10 i=1,3
  10      ixy(i) = 1
        if(nd.eq.1) ne     = 2
        if(nd.ge.2) ne     = 1
        do 30 l=1,nmax
          mr      = mr+1
          ixy(nd) = l
  20      continue
          ix     = ((ixy(3)-1)*ncmax(2) + (ixy(2)-1))*ncmax(1) + ixy(1)
          lm     = nbn(ix)
          if (lm.eq.0) then
            if (ixy(ne).ge.ncmax(ne)) then
              if (nd.eq.3) then
                ixy(2) = ixy(2)+1
                ixy(1) = 1
                if (ixy(2).le.ncmax(2)) go to 20
              else
                write(ioe,780) na
                call errtra
                return
              end if
            end if
            ixy(ne) = ixy(ne)+1
            go to 20
          end if
          ixy(ne) = 1
          call delta (ix,ndx,nbn,nbod,ma,fpd,dl,wlh)
          if (ndx.lt.0) then
            write(ioe,888) (ixy(j),j=1,3),ix,mr,na,ne,nmax
            call errtra
            return
          end if
          rby(mr) = dl+rby(mr-1)
  30      continue
        wlh(iz+nd) = rby(mr)
        mr     = mr+1
        if (iml.ne.0) then
          mb     = 7*nbi+7
          loc    = ma(mb)+2*nd+1
          fpd(loc) = rby(mr-1)
        end if
  40    continue
      return
 780  format(10x,' array no. ',i5,
     1       ' is improperly defined - fatal error in argen' )
 888  format(/5x,'argen  ixy=',3i6,'  ix=',i6,' mr=',i6,' na=',i6,
     1       ' ne=',i6,' nmax=',i6)
      end



      subroutine corner(nbo,xo,ma,fpd,d,ld,itype)
c * * this routine finds origin or vertex for reference body * *
      implicit real*8 (a-h,o-z)
      dimension xo(3),ma(*),fpd(*),d(*),ld(*)
      km=7*(nbo-1)+1
      itype=ma(km+2)
      kf=ma(7*nbo)+2
      if(itype.eq.9) go to 20
      xo(1)=fpd(kf)
      xo(2)=fpd(kf+1)
      xo(3)=fpd(kf+2)
      return
  20  continue
      xo(1)=fpd(kf)
      xo(2)=fpd(kf+2)
      xo(3)=fpd(kf+4)
      return
      end
      subroutine ctran(xb,ncmax,nx1,nx,ncn,rby,d,ld,ier)
c * * this routine determines particle coordinate relative to
c * lattice cell in array or coordinate relative to origin of array
      implicit real*8(a-h,o-z)
      real*4 d
      dimension xb(3),ncmax(*),rby(*),d(*),ld(*),nx1(3)
      dimension ncn(*)
      ier=0
      i=nx1(1)
      xb(1)=xb(1)+nx*rby(i)
      j=ncmax(1)+1+nx1(2)
      xb(2)=xb(2)+nx*rby(j)
      k=ncmax(1)+ncmax(2)+2+nx1(3)
      xb(3)=xb(3)+nx*rby(k)
      return
      end
c
      subroutine cubic(c,r,n)
c
c     solves a polynomial equation of the type
c    the coefficient of x**3 is assumed to be 1
c    r(3)contains the roots
c    n cintains the number of real roots
c    if there is one real root it will be in r(1)
c      with the complex roots r(2) +- r(3)*i
c
      implicit real*8 (a-h,o-z)
      dimension c(3),r(3)
      c1sq=c(1)*c(1)
      p=c(2)-c1sq/3.d0
      q=c(3)+c(1)*(2.*c1sq/27.d0-c(2)/3.d0)
      disc=4.*p*p*p+27.d0*q*q
      c3=c(1)/3.d0
      if(abs(disc).le.0.0001d0)disc=0.0d0
      if(disc.le.0.0d0)goto 10
c        1 real    2 complex
      n=1
      sqroot=sqrt(disc/108. )
      halfq=.5*q
      acu=-halfq+sqroot
      bcu=-halfq-sqroot
      a=sign(abs(acu)**.333333333333333d0,acu)
      b=sign(abs(bcu)**.333333333333333d0,bcu)
      ab=a+b
      r(1)=ab-c3
      r(2)=-.5*ab-c3
      r(3)=.866025404d0*(a-b)
      return
c        3 real roots
   10 n=3
      t=sqrt(abs(p)/3.d0)
      tt=t+t
      if(disc.eq.0.0)goto 20
      phi3=atan2(sqrt(-disc/27.d0),-q)/3.0d0
      r(1)=tt*cos(phi3)-c3
      r(2)=tt*cos(phi3+2.094395103d0)-c3
      r(3)=tt*cos(phi3-2.094395103d0)-c3
      return
c        2 equal roots or 3 equal roots
   20 r(1)=sign(tt,-q)-c3
      r(2)=sign(t,q)-c3
      r(3)=r(2)
      return
      end
      subroutine cwidth(ddx,ddy,ddz,rby,mx,my,mz,i,j,k)
c * * this routine called by finefi to calculate widths of cells
      implicit real*8 (a-h,o-z)
      dimension rby(*)
      ddx=rby(i+1)-rby(i)
      m=mx+1+j
      ddy=rby(m+1)-rby(m)
      n=mx+my+2+k
      ddz=rby(n+1)-rby(n)
      return
      end
      subroutine delta(i,nd,nbn,nbod,ma,fpd,dl,wlh)
c * * this routine computes size of dimension of lattice cell *
      implicit real*8 (a-h,o-z)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      dimension nbn(*),nbod(*),ma(*),fpd(*),wlh(*)
      ibod=nbn(i)
      if(ibod.gt.0) go to 10
      jbod=iabs(ibod)
      il=3*(jbod-1)+nd
      dl=wlh(il)
      return
  10  ibo=nbod(ibod)
  15  kma=7*(ibo-1)+1
      itype=ma(kma+2)
      kfpd=ma(kma+6)+2
      if(itype.eq.9) go to 20
      if(itype.ne.7) go to 40
      k7=kfpd+3*(nd-1)+3
      dl=dsqrt(fpd(k7)**2+fpd(k7+1)**2+fpd(k7+2)**2)
      return
  20  continue
      nn=nd
      k9=kfpd+2*(nn-1)
      dl=fpd(k9+1)-fpd(k9)
      return
  40  continue
      write(ioe,50) ibod,ibo,kma,itype,nd
      nd=-3
      return
  50  format(' in delta,universe reference body not rpp or box',5x,
     1 'ibod=',i5,' ibo=',i5,' ma index=',i5,' itype=',i6,' nd=',i5)
      end



      subroutine dipr(d,ip,ixd)
c * * this routine sets up dumps of mars data after error detection
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)

      dimension d(*),ip(*)

      if(nlev.le.0) return
      i3=ip(nq+3)
      i4=ip(nq+4)
      i6=ip(nq+6)
      i7=ip(nq+7)
      i9=ip(nq+9)
      i11=ip(nq+11)
      call dupr(d,d,ip,d(i3),d(i4),d(i6),d(i7),d(i9),
     1 d(i11),ixd)
      return
      end

      subroutine dxpr(fl,ntt,i,k,io)
c * * this routine prints array output when errors occur * * *
      implicit real*8 (a-h,o-z)
      dimension fl(*)
      do 10 j=1,ntt
      k1=k*(j-1)+1
      k2=k*j
      write(io,20) i,j,(fl(l),l=k1,k2)
  10  continue
  20  format(  5x,2i5,1p,5e15.5)
      return
      end

*-----------------------------------------------------------------------

      subroutine g1(s,ma,fpd,locreg,numbod,iror,ir1,ir2,nsto,distd)

c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     92/06/22
c     time the module was last permanently updated:     08:51:07
c     programmer name:                                  j.t.west
c     module name:                                      mag1
c     current archiving level number:                   00005
c     current number of permanent updates:              00005
c     date of last access by librarian:                 93/03/15
c     dataset name:  x4s.scale4.master
c
c#######################################################################

c * *  this is the control routine for combinatorial geometry.
c     in to g1, irprim=code zone
c     out of g1, ir=code zone

      implicit real*8 (a-h,o-z)

      dimension ma(*),fpd(*),locreg(*),numbod(*),iror(*),nsto(*)
      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      real*8 xb,wb,wp,xp,rin,rout,pinf,dist
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,ir,idbg ,
     1             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)
      common /dbg/ n,num,locat,isave,inext,irp,smin,inex
!$OMP THREADPRIVATE(/dbg/)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)

      real*8 dist0

      integer blzold
      common /orgi/ dist0,markg,nmedg,nblz,blzold,irpold
!$OMP THREADPRIVATE(/orgi/)

      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)

*-----------------------------------------------------------------------

      data eps/ 1.0d-5/
      save eps !FURUTA

      data ist/0/
      save ist !FURUTA
!$OMP THREADPRIVATE(ist)
*-----------------------------------------------------------------------

      lmax=ltma
      sp=0.0
      irp = irprim
      if(ll.ne.nsto(kbcz+irp-1)) go to 100
      if(nasc.gt.0) go to 110
      kloop=kloop+1

      dist = 0.0d0

  110 smin = pinf
      isave=0
  130 continue
      irr=irp
      n=locreg(irr) +1
      num=numbod(irr)*5+n-5
      if(idbg.ne.0) call pr(nsto,1)

*-----------------------------------------------------------------------
c the loop upto 300 finds the next body that the ray will intersect
*-----------------------------------------------------------------------

      do  300 i = n, num, 5

            nbo   = ma(i)
            locat = ma(i+1)

            call gg(locat,ma,fpd,nsto)

            if(ierror.ne.0) return

            if(idbg.ne.0) call pr(nsto,2)
            if( rout .le. 0.0 ) go to 300

            if(nbo) 147,300,197

  147       if( rout .le. dist )  goto 300
            if( rin  .gt. dist )  goto 1477
            go to 300

 1477       rms = rin - smin

            if(  rms  ) 149,151,148

  148       if(  rms - smin * eps ) 150,150,300
  149       if( -rms .ge. smin * eps ) goto 150

            k = ma(locat+6)
            fpd(k) = smin
            go to 151

  150       smin = rin

  151       nasc=-nbo
            lsurf=lri
            isave=i
            go to  300

  197       if( rout .gt. dist )  go to 1999
            goto 300

 1999       rms = rout - smin

            if(rms) 199,200,198

  198       if(  rms - smin * eps ) 200,200,300
  199       if( -rms .ge. smin * eps ) go to 200

            k=ma(locat+6)
            fpd(k+1)=smin
            go to 201

  200       smin = rout

  201       nasc=nbo
            lsurf=-lro
            isave=i

  300 continue

*-----------------------------------------------------------------------

      if( isave .eq. 0 ) then

         irprim = -3
         ierror = -1
         iect   = iect + 1

         if( iect .gt. 10 ) return

            write(ioe ,305) ir,xb,wb,dist
  305       format(' ****  warning from subroutine g1'/
     &      ' at region ',i5,' (x,y,z)=(',g10.3,',',g10.3,',',g10.3,')',
     &      ' direction =(',g10.3,',',g10.3,',',g10.3,')'/
     &      '  particle transport distanse=',g10.3)

         return

      end if

*-----------------------------------------------------------------------
*     now to find next region
*-----------------------------------------------------------------------

               s = smin - dist + sp

               dist = smin

               inext = isave + 2

*-----------------------------------------------------------------------
*     check for bugs by cKN
*-----------------------------------------------------------------------

               distd = abs( dist - dist0 )

*-----------------------------------------------------------------------
*     markg = 1 : no crossing
*-----------------------------------------------------------------------

         if( dist0 .lt. dist )  then

               s = s - dist + dist0

               dist = dist0

               irprim = ir
               ir = irp

               markg = 1

               return

         end if

*-----------------------------------------------------------------------
*     markg = 0 : crossing
*-----------------------------------------------------------------------

         markg = 0

*-----------------------------------------------------------------------

  310 continue

               irp = ma(inext)

                  if( idbg .ne. 0 ) call pr(nsto,3)
                  if( irp .eq. 0 ) goto  600
                  if( ll .ne. nsto(kbcz+irp-1) ) goto 100

               n   = locreg(irp) + 1
               num = numbod(irp) * 5 + n - 5

                  if( idbg .ne. 0 ) call pr(nsto,4)

*-----------------------------------------------------------------------
* the loop to 400 examines region irp to see if it is the next region
*-----------------------------------------------------------------------

         do 400 i = n, num, 5

            nbo   = ma(i)
            locat = ma(i+1)

            call gg(locat,ma,fpd,nsto)

               if( ierror .ne. 0 ) return

            dlta = eps * dist

               if( idbg .ne. 0 ) call pr(nsto,5)

            if( nbo .lt. 0 .and.
     &        ( ( rout - dist ) .gt. dlta .and.
     &          ( rin  - dist ) .le. dlta ) ) goto 500

            if( nbo .gt. 0 .and.
     &        ( ( rin  - dist ) .gt. dlta .or.
     &            dist .ge. rout ) ) goto 500

  400    continue

*-----------------------------------------------------------------------
*     not found a region
*-----------------------------------------------------------------------

            ma(inext+1) = ma(inext+1) + 1

            goto 750

*-----------------------------------------------------------------------
*     found a region
*        ( inext + 2 ) is due to new jtw component in ma array
*-----------------------------------------------------------------------

  500 continue

         inex = inext + 2

         inext = ma(inex)

            if( inext .gt. 0 ) go to 310

         inext = ldata

*-----------------------------------------------------------------------
*     sumitomo: treatment of nearest crossing point
*-----------------------------------------------------------------------

  600 continue

         if( idbg .ne. 0 ) call pr(nsto,6)

         kkloop = 0

  601 continue

*-----------------------------------------------------------------------
*        before looping thru  zone  irp see if it has been cked before
*        in table lookup loop preceding this loop
*-----------------------------------------------------------------------

      do  700 irr = 1, numr

            irp = irr

            if( ll .ne. nsto(kbcz+irp-1) ) goto 700

            if( kkloop .eq. 1 ) goto 615

            i2 = isave + 2

            if( ma(i2) .eq. 0 ) go to 615

  610    continue

            if( irp .eq. ma(i2) ) goto 700

            i2 = ma(i2+2)

            if( i2 .gt. 0 ) goto 610

  615    continue

            n   = locreg(irp) + 1
            num = numbod(irp) * 5 + n - 5

*-----------------------------------------------------------------------

         do 650 i = n, num, 5

            nbo   = ma(i)
            locat = ma(i+1)

            call gg(locat,ma,fpd,nsto)

               if( ierror .ne. 0 ) return

            dlta = eps * dist

               if( idbg .ne. 0 ) call pr(nsto,7)

            if( nbo .lt. 0 .and.
     &        ( ( rout - dist ) .gt. dlta .and.
     &          ( rin  - dist ) .le. dlta ) ) then

               if( kkloop .eq. 1 .and.
     &           ( ( rout - dist ) .le. dlta * 2.0 .or.
     &             ( rin  - dist ) .gt. dlta * 0.5 ) ) goto 650

               goto 700

            end if

            if( nbo .gt. 0 .and.
     &        ( ( rin - dist ) .gt. dlta .or.
     &            dist .ge. rout ) ) then

               if( kkloop .eq. 1 .and.
     &           ( ( rin - dist ) .le. dlta * 2.0 .and.
     &               dist .lt. rout ) ) goto 650

               goto 700

            end if

  650    continue

*-----------------------------------------------------------------------

      if( inext .ne. ldata )    goto 652
      if( inext .ge. lmax - 3 ) goto 651

         ma(inex) = inext
         ldata = ldata + 3

  652    ma(inext) = irp
         ma(inext+1) = 1

         goto 750

  651    if( ist .ne. 0 ) goto 750
         ist = 1

         write(ioe,655) lmax-3
  655 format(40x,'***********************************************'
     &      /40x,'  geometry search array full (limit=',i8,')'
     &      /40x,'***********************************************')

      goto 750

*-----------------------------------------------------------------------

  700 continue

*-----------------------------------------------------------------------

      irprim = 0

*-----------------------------------------------------------------------
*     for undefined region searching
*-----------------------------------------------------------------------

      ierror = 3
      iect = iect + 1
      if( iect .gt. 10 ) return

      kkloop = kkloop + 1

      if( kkloop .eq. 1 ) then

         write(ioe,1000) ir,iect,xb,wb,dist,
     &                   xb(1)+wb(1)*dist,
     &                   xb(2)+wb(2)*dist,
     &                   xb(3)+wb(3)*dist

         go to 601

      endif

      write(ioe,1001)

 1000 format(/'*** message in s.g1 ******************************'/
     &/' in searching region (',i4,') of point.   iect =',i4,
     &/'                       (',1p3e16.8,'),'
     &/' next point will be out of geometry or undefined region.'
     &/' particle direction    (',3e16.8,'),'
     &/' transport distance    =',e16.8
     &/' this trouble position (',3e16.8,')')

 1001 format(/' *** again error in s.g1 ***',
     &' :  the reason is above message.')

      do 710 irr = 1, numr

      if(ll.ne.nsto(kbcz+irr-1)) go to 710
      n=locreg(irr)+1
      num=numbod(irr)*5+n-5

      write(ioe,1002) irr,n,num
 1002 format(' reg=',i5,'   n,num=',2i7)

      do 711 i=n,num,5
      nbo=ma(i)
      locat=ma(i+1)
      call gg(locat,ma,fpd,nsto)

      dlta = eps * dist

      write(ioe,1003) i,nbo,rin,rout,dist,dlta
 1003 format('   i=',i5,'  nbo=',i4,'  rin,rout=',1p2e15.8,'  dist=',
     &       2e14.7)
      if(nbo) 712,711,713
  712 if((rout-dist).le.dlta.or.(rin-dist).gt.dlta) go to 711
      go to 710
  713 if((rin-dist).le.dlta.and.(dist.lt.rout)) go to 711
      go to 710
  711 continue
  710 continue

      call abend(nsto)

      return

*-----------------------------------------------------------------------

  750 continue

         irprim = iror(irp)

            if( idbg .ne. 0 ) call pr(nsto,8)

         if( irprim .ne. ir ) then
            ir = irp
            return
         end if

         sp = s
         blzold = irprim + 32768 * irp
         irprim = irp

      goto 110

*-----------------------------------------------------------------------

 100  write(ioe ,105)  irp,ll
 105  format(1h1,///,5x,'code zone',i5,' is not in universe ',i5,//)
      call pr(nsto,1)
      call abend(nsto)

      return
      end

*-----------------------------------------------------------------------


      subroutine geni(ma,fpd,    locreg,numbod,iror,mriz,mrcz,mmiz,
     1 mmcz,kr1,kr2,ngiz,ngcz,i1 ,iout,in,ibodt)

c * * this routine reads binary file written by jomin1 and
c     stores the data in combinatorial geometry array tables * * *

      implicit real*8 (a-h,o-z)
      dimension ngiz(*),ngcz(*)
      dimension ma(*),fpd(*),locreg(*),numbod(*),iror(*),mriz(*),
     1   mrcz(*),mmiz(*),mmcz(*),kr1(*),kr2(*)
      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
c
      real*8 xb,wb,wp,xp,rin,rout,pinf,dist
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,ir,idbg ,
     1             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)
c*mitu********* '91/4/25
C
CTAT/MOD(97/Oct) TO ADD C,P BODY TYPES ---------------------------(FROM)
Cc  1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
Cc arb sph rcc rec trc ell box wed rpp gel tor qua bpp wpp end
Cc 24  -2   1   6   2   1   6   6   0   6   3   4   3   3   0
C     1         4hrpp ,4hgel ,4htor ,4hqua ,4hbpp ,4hwpp ,4hend /
C
C-----------------------------------------------------------------
C (NOTE) IBIAS IS NOT USED IN CURRENT VERSION OF THIS ROUTINE
C
C BODYID       1   2   3   4   5   6   7   8   9  10  11  12  13  14
C MNEMONIC   ARB SPH RCC REC TRC ELL BOX WED RPP GEL TOR QUA BPP WPP
C IBIAS       24  -2   1   6   2   1   6   6   0   6   3   4   3   3
C
C BODYID      15  16  17  18  19  20  21  22  23  24
C MNEMONIC   P   PX  PY  PZ  PS  C   CX  CY  CZ  END
C IBIAS       -2  -5  -5  -5   1   1  -2  -2  -2
C
      dimension  ibias(23),  ity(24)
      character*4           city(24)
      equivalence (city(1),ity(1))
      data city/ 'arb ' ,'sph ' ,'rcc ' ,'rec ' ,'trc '
     1          ,'ell ' ,'box ' ,'wed ' ,'rpp ' ,'gel '
     1          ,'tor ' ,'qua ' ,'bpp ' ,'wpp ' ,'p   '
     1          ,'px  ' ,'py  ' ,'pz  ' ,'ps  ' ,'c   '
     1          ,'cx  ' ,'cy  ' ,'cz  ' ,'end '/
CTAT/MOD(97/Oct) TO ADD C,P BODY TYPES ---------------------------( TO )
C
c*mitu*********
      pinf=1.0d+20
      kloop=0
      do 80 i=1,ltma
   80 ma(i)=0
      write(iout,95)
   95 format(//50x, 9hbody data)
  105 format(2x,a3,1x,i4,1p,6e15.7/(10x,6e15.7))
  205 format(2x,a3,2i10,5x,9(2x,i5)/(30x,9(2x,i5)))
      ndx=0
      do 110 j=1,numb
      read(in)l,(fpd(ndx+i),i=1,l)
cst 98/10/30 changed for using free body number
      n=j
      n=7*n-6
      ma(n)=fpd(ndx+2)
      ma(n+2)=fpd(ndx+1)
      ma(n+6)=ndx+1
      kty=ity(ma(n+2))
      write(iout,105)kty    ,ma(n),(fpd(ndx+i),i=3,l)
      if(ma(n+2).eq.   1  )call albert(fpd(ndx+3),ierr)
      ndx=ndx+l
      if(ma(n+2).eq.1) ndx=ndx+26
c added by st informed from dr.asano for dos/v  95/08/29
      if(ma(n+2).eq.11) ndx=ndx+2
c added
  110 continue
      n=7*numb+1
      write(iout,170)numb,ndx
      l2=ndx
  170 format(1x,19hnumber of bodies   ,i5,
     .     / 1x,19hlength of fpd-array,i5)
c
c end of body data
c
      if(idbg.eq.0) go to 176
      write(iout,172)
  172 format( 50x,9hfpd array)
      do 174  i=1,l2,5
      k=i+4
  174 write( iout,175) i,(fpd(j),j=i,k),k
  175 format( i5,5e20.7,i5 )

  176 write(iout,180)
  180 format(//50x,15hinput zone data / 7x,10hinput zone,
     1 2x,10h code zone,8x,12hbody numbers )
      nx=n+4*ibodt-irtru-1
c     nczt is the total number of code zones
      nczt=0
      nczu=0
      do 300 j=1,irtru
      read(in  )l,(ma(nx+i),i=1,l)
      nf=3+ma(nx+3)
      nczu=nczu+1
      write(iout,205)ma(nx+1),j,nczu,(ma(nx+i),i=4,nf)
      nz=ma(nx+2)
      if(nz.le.1)go to 250
      do 240 kz=2,nz
      nczu=nczu+1
      ns=nf+2
      nf=nf+ma(nx+nf+1)+1
      write(iout,206)nczu,(ma(nx+i),i=ns,nf)
  240 continue
  250 continue
  206 format(15x,i10,5x,2hor,i5,8(2x,i5)/(30x,9(2x,i5)))
      izn=ma(nx+1)
      ncz=ma(nx+2)
      nx=nx+3
      do 280 nc=1,ncz
      nczt=nczt+1
      iror(nczt)=j
      mrcz(nczt)=mriz(j)
      mmcz(nczt)=mmiz(j)
      ngcz(nczt)=ngiz(j)
      nbo=ma(nx)
      nx=nx+1
      ma(n)=nczt
      locreg(nczt)=n
      n=n+1
      numbod(nczt)=nbo
      do 260 nb=1,nbo
      ma(n)=ma(nx)
      n=n+1
cst 98/10/30 changed for using free body number
cst below this line
      do 2100 k=1,numb
      if(iabs(ma(nx)).eq.ma(k*7-6)) goto 2110
 2100 continue
      call errtra
      return
 2110 ma(n)=7*k-6
cst above this line
      nx=nx+1
      ma(n+1)=0
      ma(n+2)=0
      ma(n+3)=0
      n=n+4
  260 continue
  280 continue
  300 continue
      ldata=n
      do 305 i=n,ltma
      ma(i)=0
  305 continue
      write(iout,310) irtru,numr,ltma
  310 format(  1x,23hnumber of input zones  ,i5
     .       / 1x,23hnumber of  code zones  ,i5
     .       / 1x,23hlength of integer array,i5)
      write(iout,315)
 315  format(//,' code zone    input zone    zone data loc.    no. of',
     1 ' bodies   region no.   media no.    box input zone   box code',
     2 ' zone ')
      write(iout,317) (i,iror(i),locreg(i),numbod(i),mrcz(i),
     1 mmcz(i),ngiz(iror(i)),ngcz(i),i=1,numr)
 317  format(i7,i13,i16,i18,i16,i13,i17,i18)
      kr1(1)=1
      kr2(1)=1
      l=1
      if(numr.le.1) go to 405
      do 400 i=2,numr
      if(iror(i-1)-iror(i) ) 390,380,390
  380 kr2(l)=i
      go to 400
  390 l=l+1
      kr1(l)=i
      kr2(l)=i
  400 continue
  405 write(iout,410) (i,kr1(i),kr2(i),i=1,l)
  410 format(//25h    i    kr1(i)    kr2(i)/(i5,i8,i10))
      if(idbg.eq.0) return
      write( iout,320 )
  320 format( 50x, 10h  ma-array)
      do  325  i=1,ldata,7
      k=i+6
  325 write( iout,326) i,(ma(j),j=i,k),k
  326 format(20x,i5,i10,6i5,i10)
c
c end of region data
c
      return
      end


*-----------------------------------------------------------------------

      subroutine gg(locat,ma,fpd,ld)

c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     92/06/22
c     time the module was last permanently updated:     08:53:41
c     programmer name:                                  j.t.west
c     module name:                                      magg
c     current archiving level number:                   00002
c     current number of permanent updates:              00002
c     date of last access by librarian:                 93/03/15
c     dataset name:  x4s.scale4.master
c
c#######################################################################

c      this routine is the workhorse of the combinatorial geometry pkg
c       it computes the distance to intersection for all body types  *

c  * * parameter ld added to gg so blank common unnecessary in pr * *

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      dimension ma(*),fpd(*),ld(*)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      dimension vhab(12),v(3),h(3),p(3),q(3)
      dimension asq(3),pv(3),g(3)
      real*8 xb,wb,wp,xp,rin,rout,pinf,dist
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,ir,idbg ,
     1             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)

      dimension am(3,3),axc(5),azz(6),azzz(4),axcq(4)

      dimension az(4)
      parameter  (pi = 3.14159265358979d0)

      dimension  ibox(3)
      equivalence (vhab(1),v(1)),(vhab(4),h(1)),(vhab(7),p(1)),
     1            (vhab(10),q(1))
      data ibox/2,1,3/

      data eps/1.0d-5/
      save eps !FURUTA
!$OMP THREADPRIVATE(eps)
      data pinf90/5.0d19/

*-----------------------------------------------------------------------

      l=locat
      loop=ma(l+1)
      itype=ma(l+2)
      k=ma(l+6)

      if( loop.eq.kloop ) then

        lri = ma(l+3)
        lro = ma(l+4)
        rin = fpd(k)
        rout= fpd(k+1)
        return

      endif

      if( loop.ne.kloop ) go to 1000
      lri=ma(l+3)
      lro=ma(l+4)
      rin=fpd(k)
      rout=fpd(k+1)
      return
 1000 rin=pinf
      rout=-pinf

      if(itype.eq.11) goto 1001

*-----------------------------------------------------------------------
C     CHECK OUT OF BOUNDS & BPP,WPP;
C     BECAUSE BOTH OF BPP AND WPP ARE ALREADY PROCESSED IN JOMIN1
*-----------------------------------------------------------------------

      if(  itype.le. 0                  ) goto 2011
      if(  itype.gt.23                  ) goto 2011
      if(  itype.ge.13.and.itype.le.14  ) goto 2011

      ma(l+1)=kloop
      ma(l+5) = ll
      eps = 1.0d-5

 1001 continue

*-----------------------------------------------------------------------
C             ARB  SPH  RCC  REC  TRC  ELL  BOX  WED  RPP
C             GEL  TOR  QUA (BPP  WPP)
C             P    PX   PY   PZ   PS   C    CX   CY   CZ
*-----------------------------------------------------------------------

      GO TO (1100,1200,1400,1300,1400,1600,1700,1800,1700,
     :       1900,2100,2200,2011,2011,
     :       3100,3200,3300,3400,3500,3600,3700,3800,3900 ),ITYPE


*-----------------------------------------------------------------------
C -----( P   )----- PLANE
*-----------------------------------------------------------------------

3100  CONTINUE
      LP=K+1
      WBWP=FPD(LP+1)*WB(1)
     1    +FPD(LP+2)*WB(2)
     1    +FPD(LP+3)*WB(3)
      XBWP=FPD(LP+4)
     1    -FPD(LP+1)*XB(1)
     1    -FPD(LP+2)*XB(2)
     1    -FPD(LP+3)*XB(3)

*-----------------------------------------------------------------------
C ----- COMMON PROCESS FOR P,PX,PY,PZ,PS -----
*-----------------------------------------------------------------------

3101  LRI =1
      LRO =1
      IF(WBWP) 3104,3102,3103
3102  ROUT=PINF90
      RIN =SIGN(PINF90,XBWP)
      IF(XBWP.EQ.0.) RIN=PINF90
      GOTO 2000
3103  ROUT=PINF90
      RIN =XBWP/WBWP
      GOTO 2000
3104  RIN =SIGN(PINF90,XBWP)
      ROUT=XBWP/WBWP
      IF (XBWP.GE.0.) ROUT = PINF90
      GOTO 2000

*-----------------------------------------------------------------------
C -----( PX  )----- X-PLANE
*-----------------------------------------------------------------------

3200  CONTINUE
      LP=K+1
      WBWP=WB(1)
      XBWP=FPD(LP+1)-XB(1)
      GOTO 3101

*-----------------------------------------------------------------------
C -----( PY  )----- Y-PLANE
*-----------------------------------------------------------------------

3300  CONTINUE
      LP=K+1
      WBWP=WB(2)
      XBWP=FPD(LP+1)-XB(2)
      GOTO 3101

*-----------------------------------------------------------------------
C -----( PZ  )----- Z-PLANE
*-----------------------------------------------------------------------

3400  CONTINUE
      LP=K+1
      WBWP=WB(3)
      XBWP=FPD(LP+1)-XB(3)
      GOTO 3101

*-----------------------------------------------------------------------
C -----( PS  )----- SHIFT PLANE
*-----------------------------------------------------------------------

3500  CONTINUE
      LP=K+1
      WBWP=FPD(LP+4)*WB(1)
     1    +FPD(LP+5)*WB(2)
     1    +FPD(LP+6)*WB(3)
      XBWP=FPD(LP+7)-(XB(1)-FPD(LP+1))*FPD(LP+4)
     1              -(XB(2)-FPD(LP+2))*FPD(LP+5)
     1              -(XB(3)-FPD(LP+3))*FPD(LP+6)
      GOTO 3101

*-----------------------------------------------------------------------
C -----( C   )----- CYLINDER
*-----------------------------------------------------------------------

3600  CONTINUE
      LP=K+1
      WP(1)=FPD(LP+4)
      WP(2)=FPD(LP+5)
      WP(3)=FPD(LP+6)
      R=FPD(LP+7)

*-----------------------------------------------------------------------
C ----- COMMON PROCESS FOR C,CX,CY,CZ -----
*-----------------------------------------------------------------------

3650  CONTINUE
      R2=R**2
      XPWP=0.
      WBWP=0.
      DO 3651 I=1,3
         XP(I)=XB(I)-FPD(LP+I)
         XPWP =XPWP +XP(I)*WP(I)
         WBWP =WBWP +WB(I)*WP(I)
3651  CONTINUE
      XRXR = 0.
      URUR = 0.
      XRUR = 0.
      DO 3652 I=1,3
         XRI  = XP(I)-WP(I)*XPWP
         URI  = WB(I)-WP(I)*WBWP
         XRXR = XRXR +XRI  *XRI
         URUR = URUR +URI  *URI
         XRUR = XRUR +XRI  *URI
3652  CONTINUE
      B   =XRXR-R2
      IF(URUR) 3653,3655,3653
3653  A=XRUR/URUR
      B=B/URUR
      DISC=A**2-B
      IF(DISC) 3655,3654,3654
3654  DISC=SQRT(DISC)
      RIN =-A - DISC
      ROUT=-A + DISC
      LRI = 1
      LRO = 1
      GOTO 2000
3655  ROUT = PINF90
      RIN  = SIGN(PINF90,B)
      LRI = 1
      LRO = 1
      GOTO 2000

*-----------------------------------------------------------------------
C -----( CX  )----- X-CYLINDER
*-----------------------------------------------------------------------

3700  CONTINUE
      WP(1)=1.
      WP(2)=0.
      WP(3)=0.
      LP=K+1
      R =FPD(LP+4)
      GOTO 3650

*-----------------------------------------------------------------------
C -----( CY  )----- Y-CYLINDER
*-----------------------------------------------------------------------

3800  CONTINUE
      WP(1)=0.
      WP(2)=1.
      WP(3)=0.
      LP=K+1
      R =FPD(LP+4)
      GOTO 3650

*-----------------------------------------------------------------------
C -----( CZ  )----- Z-CYLINDER
*-----------------------------------------------------------------------

3900  CONTINUE
      WP(1)=0.
      WP(2)=0.
      WP(3)=1.
      LP=K+1
      R =FPD(LP+4)
      GOTO 3650


*-----------------------------------------------------------------------
c arb  arbitrary polyhedron
*-----------------------------------------------------------------------

 1100 dmin=fpd(k+26)*1.0d-7
      nside=fpd(k+27)
      i1=k+2
      i2=i1+ 4*nside -1
      l=0
      do  1190  i=i1,i2,4
      dx=   fpd(i)*xb(1) + fpd(i+1)*xb(2)+fpd(i+2)*xb(3) + fpd(i+3)
      l=l+1
      dy=  -fpd(i)*wb(1)  -fpd(i+1)*wb(2)-fpd(i+2)*wb(3)
      if( dabs(dy ).le.1.0d-7 )  go to 1190
      dz=dx/dy
      pv(1)= xb(1) + dz*wb(1)
      pv(2)= xb(2) + dz*wb(2)
      pv(3)= xb(3) + dz*wb(3)
      do   1120  j=i1,i2,4
      if( j.eq.i) go to  1120
      dx=fpd(j)*pv(1)+fpd(j+1)*pv(2)+fpd(j+2)*pv(3)+fpd(j+3)
      if(  dx.ge.0 ) go to  1120
      if( -dx.lt.dmin ) go to  1120
      go  to 1190
 1120 continue
      if(dz .gt.rout) go to 1140
      rin = dz
      lri= l
      go to 2000
 1140 rin=rout
      lri= lro
      lro=l
      rout= dz
      if(rin .gt.-pinf) go to 2000

 1190 continue

      go to 2000

*-----------------------------------------------------------------------
c     sphere     sph
*-----------------------------------------------------------------------

 1200 continue
      dx= xb(1) -fpd(k+2)
      dy= xb(2) -fpd(k+3)
      dz= xb(3) -fpd(k+4)
      b= dx*wb(1) + dy*wb(2) + dz*wb(3)
      c= dx*dx + dy*dy + dz*dz -fpd(k+5)**2
      dx=b*b- c
      if( dx.lt.0.) go to 2000
      dy=dsqrt(dx)
      rin= -b-dy
      rout=-b+dy
      lri=1
      lro=1
      go to 2000

*-----------------------------------------------------------------------
c     rec   right elliptical cylinder
*-----------------------------------------------------------------------

 1300 do 1301 i=1,12
      j=k+1+i
 1301 vhab(i)=fpd(j)
      rin=-pinf
      rout=pinf
      lro=0
      lri=0

*-----------------------------------------------------------------------
c4    compute dot products of p.p and q.q
*-----------------------------------------------------------------------

      aa=p(1)*p(1)+p(2)*p(2)+p(3)*p(3)
      bb=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)

*-----------------------------------------------------------------------
c5    compute (v-xb) for x,y,z coordinates
*-----------------------------------------------------------------------

      v1xb1=v(1)-xb(1)
      v2xb2=v(2)-xb(2)
      v3xb3=v(3)-xb(3)

*-----------------------------------------------------------------------
c6    transform xb(x,y,z) to the coordinates of the rec
*-----------------------------------------------------------------------

      vpa=v1xb1*p(1)+v2xb2*p(2)+v3xb3*p(3)
      vpb=v1xb1*q(1)+v2xb2*q(2)+v3xb3*q(3)

*-----------------------------------------------------------------------
c7    transform wb(x,y,z) to the coordinates of the rec
*-----------------------------------------------------------------------

      wba=wb(1)*p(1)+wb(2)*p(2)+wb(3)*p(3)
      wbb=wb(1)*q(1)+wb(2)*q(2)+wb(3)*q(3)
      wbawba=wba*wba
      wbbwbb=wbb*wbb
      aaaa=aa*aa
      bbbb=bb*bb
      ambd=wba*vpa*bbbb+wbb*vpb*aaaa
      um=bbbb*vpa*vpa+aaaa*vpb*vpb-aaaa*bbbb
      den=wbawba*bbbb+wbbwbb*aaaa
      if(dabs(den).le.1.0d-6)go to 10
      ambda=ambd/den
      umu=um/den
      disc=ambda**2-umu
      if(disc.le.0.)goto 300
c
c8    compute the intersect points on the quadratic surface
c
      sd=dsqrt(disc)
      r1=ambda-sd
      r2=ambda+sd
      goto 20
   10 r1=-pinf
      r2=pinf
   20 hh=h(1)*h(1)+h(2)*h(2)+h(3)*h(3)
      wh=wb(1)*h(1)+wb(2)*h(2)+wb(3)*h(3)
      vph=v1xb1*h(1)+v2xb2*h(2)+v3xb3*h(3)
c
c9    determine if ray parallel to planar surfaces
c
      if(wh)40,70,50
   40 if(vph.ge.0.)goto 300
c
c10   compute the intersect points on the planar surfaces
c
      cp=vph/wh
      cm=(vph+hh)/wh
      lcp=1
      lcm=2
      goto100
   50 vphhh=vph+hh
      if(vphhh.le.0.)goto 300
      cp=vphhh/wh
      cm=vph/wh
      lcm=1
      lcp=2
      goto 100
   70 cp=pinf
      cm=-cp
  100 if(cm.gt.r1)goto 110
c
c11   rin for the quadratic surface
c
      rin=r1
      lri=3
      goto 120
c
c12   rin for a planar surface
c
  110 rin=cm
      lri=lcm
  120 if(cp.le.r2)goto 130
c
c13   rout for the quadratic surface
c
      rout=r2
      lro=3
      goto 200
c
c14   rout for a planar surface
c
  130 rout=cp
      lro=lcp
 200  if(rout.le.rin)goto300
      goto(210,210,220),lro
c
c15   determine if rout of planar surface occurs within elliptic
c     cross-section
c
  210 f1=den*rout**2-2.*ambd*rout+um
      if(f1)250,250,300
c
c16   determine if rout of quadratic occurs between planar surfaces
c
  220 f1=rout*wh-vph
      if(f1)300,250,230
  230 if(f1.gt.hh)goto 300
  250 goto(260,260,270),lri
c
c17   determine if rin of plane within elliptic cross section
c
  260 f1=den*rin**2-2.*ambd*rin+um
      if(f1.gt.0.0)go to 300
      go to 2000
c
c18   determine if rin of quadratic surface between planar surfaces
c
  270 f1=rin*wh-vph
      if(f1)300,2000,280
  280 if(f1.le.hh)goto 2000
c
c19   ray misses body
c
  300 rin=pinf
      rout=-pinf
      lri=0
      lro=0
      go to 2000

*-----------------------------------------------------------------------
c     rcc and trc  - right circular cyl and truncated cone
*-----------------------------------------------------------------------

 1400 rb=fpd(k+8)
      rt=fpd(k+9)

      if(itype.eq.3)rt=rb

      dx=fpd(k+2)-xb(1)
      dy=fpd(k+3)-xb(2)
      dz=fpd(k+4)-xb(3)

      h1=fpd(k+5)
      h2=fpd(k+6)
      h3=fpd(k+7)

      intsec=0
      intr1=0
      intr2=0

      pvpv=dx**2 + dy**2 + dz**2
      vpw=dx*wb(1) + dy*wb(2) + dz*wb(3)
      wh= h1*wb(1) + h2*wb(2) + h3*wb(3)

      vph=h1*dx + h2*dy + h3*dz
      hh =h1**2 + h2**2 + h3**2

      rtrb=rt-rb
      rrr= rb-rtrb/hh*vph
      vphhh=vph + hh
      um= hh*(pvpv-rrr**2) -vph**2
      ambd=hh*vpw -wh*(vph-rtrb*rrr)
      den= hh -wh**2*(1.0+rtrb**2/hh)


      if(dabs(1.-dabs(wh)/dsqrt(hh))) 1401,1401,1404

 1401 continue

      d = dmin1(rb,rt)
      if( pvpv - vpw**2 .lt. d**2 )  goto 1470

 1404 continue

      if(dabs(den).gt.1.0d-6 ) go to 1420

      if(rtrb.eq.0) goto 1470
      r2=um/(2.0*ambd)
      f1=r2*wh-vph
      if(f1.lt.0.0) goto 1470
      if((f1-hh).gt.0.0) goto 1470
      intsec=intsec+1
      if(wh.le.0.0)go to 1405
      if(rtrb)1410,1410,1415
 1405 if(rtrb.le.0.0)go to 1415
 1410 lro=3
      rout=r2
      go to 1480
 1415 lri=3
      rin=r2
      intsec=intsec+1
      go to 1472

 1420 ambda=ambd/den
      disc=ambda**2 -um/den
      if(disc)1498,1470,1422

 1422 continue

      sd=dsqrt(disc)
      r1=ambda-sd
      r2=ambda+sd
      f1=r2*wh-vph
      if(f1.lt.0.0)go to 1424
      if((f1-hh).gt.0.0)go to 1424
      intr2=intr2+1

 1424 f1=r1*wh-vph
      if(f1.lt.0.0)go to 1426
      if( (f1-hh).gt.0.) go to 1426
      intr1=intr1+1
      go to 1430

 1426 if(intr2.eq.0)go to 1470
      rout=r2
      rin=r2
      lro=3
      lri=3
      intsec=intsec+1
      go to 1470

 1430 if(intr2.gt.0)go to 1432
      rout=r1
      rin=r1
      lro=3
      lri=3
      intsec=intsec+1
      go to 1470

 1432 if(r1-r2)1434,1498,1436

 1434 rin=r1
      rout=r2
      lro=3
      lri=3
      go to 1496

 1436 rin=r2
      rout=r1
      lro=3
       lri=3
      go to 1496

 1470 continue

      if(wh)1472,1498,1480
 1472 if(vph.ge.0.0)go to 1498
      cp=vph/wh
      f1=cp**2 -2.0*cp*vpw + pvpv -rb**2
      if(f1.gt.0.0)go to 1474
      intsec=intsec+1
      rout=cp
      lro=1
      if(intsec.ge.2)go to 1496
 1474 cm=vphhh/wh
      f1=cm**2-2.0*((vpw+wh)*cm-vph)+hh+pvpv-rt**2
      if(f1.gt.0.0)go to 1498
      rin=cm
      lri=2
      go to 1496
 1480 if(vphhh.lt.0.0)go to 1498
      cp=vphhh/wh
      f1=cp**2-2.0*((vpw+wh)*cp-vph)+hh+pvpv-rt**2
      if(f1.gt.0.0)go to 1486
      intsec=intsec+1
      rout=cp
      lro=2
 1486 if(intsec.gt.1)go to 1496
      cm=vph/wh
      f1=cm**2-2.0*cm*vpw+pvpv-rb**2
      if(f1.gt.0.0)go to 1498
      rin=cm
      lri=1
 1496 go to 2000

 1498 continue

      rin  = 0.0
      rout = -pinf
      go to 2000

*-----------------------------------------------------------------------
c     ell    ellipsoid                      *           *
*-----------------------------------------------------------------------

 1600 continue
      a1=0.
      a2=0.
      b1=0.
      b2=0.
      ja=k+1
      do   1610  j=1,3
      ja=ja+1
      dx=xb(j)-fpd(ja)
      a1=a1+ dx*wb(j)
      b1=b1+ dx*dx
      dx=xb(j)-fpd(ja+3)
      a2=a2+ dx*wb(j)
 1610 b2=b2+ dx*dx
      a1=2.0*a1
      a2=2.0*a2
      c=fpd(k+8)
      c2=2.0*c
      a=(a2-a1)/c2
      b=( c**2+b2-b1)/c2
      alamd=a*a-1
      alam1=( a*b-a2*0.5)/alamd
      u=(b*b-b2)/alamd
      c=alam1*alam1-u
      if(c.lt.0 ) go to 2000
      c=dsqrt(c)
      rin= -alam1-c
      rout=-alam1+c
      lri=1
      lro=1
      go to 2000

*-----------------------------------------------------------------------
c     box and rpp     box     rectangular parallepiped
*-----------------------------------------------------------------------

 1700 continue
      rin =-pinf
      rout=+pinf
      do 1798   i=1,3
      if( itype- 9 )  1705,1701,1701
 1701 jv=k+2*i+1
      a =fpd(jv)-fpd(jv-1)
      vp= fpd(jv-1)-xb(i)
      w =wb(i)
      eps = 1.0d-6

      go to 1711
 1705 jv=k+1
      a=0.
      vp=0.
      w=0.
      ja=jv+3*i
      do  1710  j=1,3
      jv=jv+1
      ja=ja+1
      dx=fpd(ja)
      vp=vp+ (fpd(jv)  -xb(j))*dx
      w=w+ wb(j)*dx
 1710 a=a+dx*dx
 1711 if( w ) 1720,1712,1740
 1712 if(-vp.lt.0) go to 1799
      if(-vp-a )1798,1798,1799
 1720 dy=vp/w
      lo=2*ibox(i)-1
      if(dy.le.0 ) go to 1799
      dz=(vp+a)/w
      li=lo+1
      go to 1760
 1740 dy=(vp+a)/w
      lo=2*ibox(i)
      if( dy.le.0) go to 1799
      dz=vp/w
      li=lo-1
 1760 if(rout.le.dy) go to 1780
      rout=dy
      lro=lo
 1780 if( rin.ge.dz) go to 1798
      rin=dz
      lri=li
 1798 continue
      go to 2000
 1799 rin=pinf
      rout=-pinf
      go to 2000

*-----------------------------------------------------------------------
c     wed or raw    wedge            *            *
*-----------------------------------------------------------------------

 1800 continue
      rin=-pinf
      rout=pinf
      cm=-pinf
      cp=pinf
      l=0
      l1=0
      kk=0
      lri=0
      lro=0
      dx= xb(1)- fpd(k+2)
      dy= xb(2)- fpd(k+3)
      dz= xb(3)- fpd(k+4)
      do 1830  i=1,3
      jv=k+2+3*i
      asq(i)=fpd(jv)**2 + fpd(jv+1)**2 + fpd(jv+2)**2
      pv(i)=dx*fpd(jv) + dy*fpd(jv+1) + dz*fpd(jv+2)
      g(i)=wb(1)*fpd(jv)+wb(2)*fpd(jv+1) + wb(3)*fpd(jv+2)
      if( i.eq.3 ) go to 1801
      if( g(i) )  1810,1811,1860
 1810 if(-pv(i).ge.0) go to 1840
      temp=-pv(i)/g(i)
      if( temp.ge.cp) go to 1830
      cp=temp
      l=i
      if( i.gt.1 ) go to 1850
      lro=3
      go to 1830
 1850 lro=1
      go to 1830
 1860 if(-pv(i).le.0 ) go to 1830
      temp=-pv(i)/g(i)
      if(temp.le.cm) go to 1830
      cm=temp
      kk=i
      lri=3
      if( i.eq.1) go to 1830
      lri=1
      go to 1830
 1811 if( pv(i).le.0.) go to    1881
      if( pv(i).ge.asq(i)) go to 1881
 1830 l1=l1+i
 1801 if( g(3) ) 1815,1821,1823
 1815 temp=-pv(3)+asq(3)
      if(temp.ge.0.) go to 1818
      temp=temp/g(3)
      if( temp.le.cm) go to 1819
      cm=temp
      kk=3
      lri=6
 1818 if(-pv(3).ge.0.)  go to 1840
 1819 temp=-pv(3)/g(3)
      if(temp.ge.cp) go to 1829
      cp=temp
      l=3
      lro=5
      go to 1829
 1821 if(pv(3).le.0.) go to 1840
      if(pv(3).gt.asq(3) ) go to 1840
      go to 1829
 1823 if(-pv(3) .le. 0.) go to 1826
      temp = -pv(3)/g(3)
      if(temp.le.cm) go to 1826
      cm=temp
      kk=3
      lri=5
 1826 temp=-pv(3)+asq(3)
      if(temp.le.0) go to 1840
      temp=temp/g(3)
      if( temp.ge.cp ) go to 1829
      cp=temp
      l=3
      lro=6
 1829 ag=asq(2)*g(1) + asq(1)*g(2)
      pv4=pv(1)*asq(2) + pv(2)*asq(1)
      top=asq(1)*asq(2)-pv4
      if(ag) 1831,1835,1833
 1831 temp=top/ag
      if(temp.le.cm) go to 1838
      cm=temp
      kk=4
      lri=2
      go to 1838
 1833 if(top.lt.0.) go to 1840
      temp=top/ag
      if( temp-cp) 1837,1838,1838
 1835 if(pv4.le.0.) go to 1840
      if(-top) 1838,1840,1840
 1837 cp=temp
      l=4
      lro=2
 1838 if(l+kk.le.0) go to 1840
      rout=cp
      rin=cm
 1840 continue
      if( (rout.lt.pinf).and.(rout.gt.0.).and.(rout.gt.rin)) go to 2000
 1881 rout=-pinf
      rin=pinf
      lri=0
      lro=0
      go to 2000

*-----------------------------------------------------------------------
c*     gel    general ellipsoid
*-----------------------------------------------------------------------

 1900 continue
      dx    = xb(1) - fpd(k+2)
      dy    = xb(2) - fpd(k+3)
      dz    = xb(3) - fpd(k+4)

      r1 =  fpd(k+5)**2 + fpd(k+6)**2 + fpd(k+7)**2
      r2 =  fpd(k+8)**2 + fpd(k+9)**2 + fpd(k+10)**2
      r3 =  fpd(k+11)**2 + fpd(k+12)**2 + fpd(k+13)**2
      do 1250 i = 1,3
        am(1,i) = fpd(k+4+i)/r1
        am(2,i) = fpd(k+7+i)/r2
        am(3,i) = fpd(k+10+i)/r3
 1250 continue

      aa = 0.0
      bb = 0.0
      cc = -1.0
      do 1260 i = 1,3
        ax =  am(i,1)*dx + am(i,2)*dy + am(i,3)*dz
        ae =  am(i,1)*wb(1) + am(i,2)*wb(2) + am(i,3)*wb(3)
        aa = aa + ae**2
        bb = bb + ae*ax
        cc = cc + ax**2
 1260 continue

      dd = bb**2 - aa*cc
      if(dd.lt.0.0) then
        go to 2000
      else
        dd   =  sqrt(dd)
        rout = ( -bb + dd)/aa
        rin  = ( -bb - dd)/aa
        lri  = 1
        lro  = 1
        go to 2000
      endif


*-----------------------------------------------------------------------
c* torus
*-----------------------------------------------------------------------

 2100 continue
      fpd(k+11)= rin
      fpd(k+12)= rout
      dx = xb(1) - fpd(k+2)
      dy = xb(2) - fpd(k+3)
      dz = xb(3) - fpd(k+4)
      r1 = fpd(k+5)


*-----------------------------------------------------------------------
c* set the coefficients (ax(5)) of the quartic equation.
*-----------------------------------------------------------------------

      a2 = fpd(k+6)**2
      b2 = fpd(k+7)**2
      ixyz = fpd(k+8) + 0.1
      if( ixyz.eq.1) then
        aa = b2 + (a2-b2)*wb(1)**2
        bb = b2*wb(2)*dy + b2*wb(3)*dz + a2*wb(1)*dx
        cc = b2*dy**2 + b2*dz**2 + a2*dx**2 + b2*(r1**2-a2)
        dd = (2.0*r1*b2)**2*(1.0-wb(1)**2)
        ee = (2.0*r1*b2)**2*(wb(2)*dy+wb(3)*dz)
        ff = (2.0*r1*b2)**2*(dy**2+dz**2)
      endif
      if( ixyz.eq.2) then
        aa = b2 + (a2-b2)*wb(2)**2
        bb = b2*wb(1)*dx + b2*wb(3)*dz + a2*wb(2)*dy
        cc = b2*dx**2 + b2*dz**2 + a2*dy**2 + b2*(r1**2-a2)
        dd = (2.0*r1*b2)**2*(1.0-wb(2)**2)
        ee = (2.0*r1*b2)**2*(wb(1)*dx+wb(3)*dz)
        ff = (2.0*r1*b2)**2*(dx**2+dz**2)
      endif
      if(ixyz.eq.3) then
        aa = b2 + (a2-b2)*wb(3)**2
        bb = b2*wb(1)*dx + b2*wb(2)*dy + a2*wb(3)*dz
        cc = b2*dx**2 + b2*dy**2 + a2*dz**2 + b2*(r1**2-a2)
        dd = (2.0*r1*b2)**2*(1.0-wb(3)**2)
        ee = (2.0*r1*b2)**2*(wb(1)*dx+wb(2)*dy)
        ff = (2.0*r1*b2)**2*(dx**2+dy**2)
      endif

      axc(1) = aa*aa
      axc(2) = 4.0*aa*bb
      axc(3) = 4.0*bb*bb + 2.0*aa*cc - dd
      axc(4) = 4.0*bb*cc - 2.0*ee
      axc(5) = cc*cc - ff

*-----------------------------------------------------------------------
c* solve the quartic equation.
*-----------------------------------------------------------------------

      do 3954 i=1,4
         axcq(i)=axc(i+1)/axc(1)
 3954 continue

      call qrtic(axcq,az,n)

         if( n .eq. 0 ) goto 2000

c* select the real roots.

      ii = 0

      do 1270 i = 1,n
        ii = ii + 1
        azz(ii) = az(i)
 1270 continue

      if(ii.eq.0) go to 2000
      if(ii.eq.1.or.ii.eq.3) then
        write(6,*) '===== error in gg to find real roots ====='
        go to 2011
      endif

c* select only the real roots corresponding to the specified
c* angular range.
c* find the intersection points on the planar surface of the torus.

      theta1 = (pi/180.0d0)*fpd(k+ 9)
      theta2 = (pi/180.0d0)*fpd(k+10)

      if( theta1.gt.theta2) then
        write(6,*)  '***** error in gg (theta1>theta2) *****'
        go to 2011
      endif

      ipq = 0
      if( theta1.ne.theta2) then
        do 1280 i = 1, ii
          if(ixyz.eq.1) then
            x = dy + azz(i)*wb(2)
            y = dz + azz(i)*wb(3)
          endif
          if(ixyz.eq.2) then
            x = dz + azz(i)*wb(3)
            y = dx + azz(i)*wb(1)
          endif
          if(ixyz.eq.3) then
            x = dx + azz(i)*wb(1)
            y = dy + azz(i)*wb(2)
          endif
          r  = sqrt( x**2+y**2 )
          if( y.ge.0.d0) then
            theta = acos( x/r )
          else
            theta = 2.0*pi - acos( x/r )
          endif
          if( theta.le.theta2.and.theta.ge.theta1) then
            ipq = ipq + 1
            azzz(ipq) = azz(i)
          endif
 1280   continue
        ii    = ipq
        do 1281 i = 1,ii
          azz(i) = azzz(i)
 1281   continue
        if(ii.eq.4) go to 1282
        alpha = sin(theta1)
        beta  =-cos(theta1)
        alphap= sin(theta2)
        betap =-cos(theta2)
        if(ixyz.eq.1) then
          gama  = -fpd(k+3)*sin(theta1)+fpd(k+4)*cos(theta1)
          gamap = -fpd(k+3)*sin(theta2)+fpd(k+4)*cos(theta2)
          den   =  alpha*wb(2)+beta*wb(3)
          denp  =  alphap*wb(2)+betap*wb(3)
          if( abs(den ).gt.eps) then
            ii  = ii + 1
            azz(ii) = -(alpha *xb(2)+beta *xb(3)+gama )/den
            xx  = xb(1) + azz(ii)*wb(1)
            yy  = xb(2) + azz(ii)*wb(2)
            zz  = xb(3) + azz(ii)*wb(3)
            val = (yy-fpd(k+3))**2+(zz-fpd(k+4))**2
            val = (sqrt(val)-r1)**2/a2+(xx-fpd(k+2))**2/b2
            prod = -(yy-fpd(k+3))*beta+(zz-fpd(k+4))*alpha
            if(val.gt.1.0d0.or.prod.le.0.0) ii = ii - 1
          endif
          if( abs(denp).gt.eps) then
            ii  = ii + 1
            azz(ii) = -(alphap*xb(2)+betap*xb(3)+gamap)/denp
            xx  = xb(1) + azz(ii)*wb(1)
            yy  = xb(2) + azz(ii)*wb(2)
            zz  = xb(3) + azz(ii)*wb(3)
            val = (yy-fpd(k+3))**2+(zz-fpd(k+4))**2
            val = (sqrt(val)-r1)**2/a2+(xx-fpd(k+2))**2/b2
            prod = -(yy-fpd(k+3))*betap+(zz-fpd(k+4))*alphap
            if(val.gt.1.0d0.or.prod.le.0.0) ii = ii - 1
          endif
        endif
        if(ixyz.eq.2) then
          gama  = -fpd(k+4)*sin(theta1)+fpd(k+2)*cos(theta1)
          gamap = -fpd(k+4)*sin(theta2)+fpd(k+2)*cos(theta2)
          den   =  alpha*wb(3)+beta*wb(1)
          denp  =  alphap*wb(3)+betap*wb(1)
          if( abs(den ).gt.eps) then
            ii  = ii + 1
            azz(ii) = -(alpha *xb(3)+beta *xb(1)+gama )/den
            xx  = xb(1) + azz(ii)*wb(1)
            yy  = xb(2) + azz(ii)*wb(2)
            zz  = xb(3) + azz(ii)*wb(3)
            val = (zz-fpd(k+4))**2+(xx-fpd(k+2))**2
            val = (sqrt(val)-r1)**2/a2+(yy-fpd(k+3))**2/b2
            prod = -(zz-fpd(k+4))*beta+(xx-fpd(k+2))*alpha
            if(val.gt.1.0d0.or.prod.le.0.0) ii = ii - 1
          endif
          if( abs(denp).gt.eps) then
            ii  = ii + 1
            azz(ii) = -(alphap*xb(3)+betap*xb(1)+gamap)/denp
            xx  = xb(1) + azz(ii)*wb(1)
            yy  = xb(2) + azz(ii)*wb(2)
            zz  = xb(3) + azz(ii)*wb(3)
            val = (zz-fpd(k+4))**2+(xx-fpd(k+2))**2
            val = (sqrt(val)-r1)**2/a2+(yy-fpd(k+3))**2/b2
            prod = -(zz-fpd(k+4))*betap+(xx-fpd(k+2))*alphap
            if(val.gt.1.0d0.or.prod.le.0.0) ii = ii - 1
          endif
        endif
        if(ixyz.eq.3) then
          gama  = -fpd(k+2)*sin(theta1)+fpd(k+3)*cos(theta1)
          gamap = -fpd(k+2)*sin(theta2)+fpd(k+3)*cos(theta2)
          den   =  alpha*wb(1)+beta*wb(2)
          denp  =  alphap*wb(1)+betap*wb(2)
          if( abs(den ).gt.eps) then
            ii  = ii + 1
            azz(ii) = -(alpha *xb(1)+beta *xb(2)+gama )/den
            xx  = xb(1) + azz(ii)*wb(1)
            yy  = xb(2) + azz(ii)*wb(2)
            zz  = xb(3) + azz(ii)*wb(3)
            val = (xx-fpd(k+2))**2+(yy-fpd(k+3))**2
            val = (sqrt(val)-r1)**2/a2+(zz-fpd(k+4))**2/b2
            prod = -(xx-fpd(k+2))*beta+(yy-fpd(k+3))*alpha
            if(val.gt.1.0d0.or.prod.le.0.0) ii = ii - 1
          endif
          if( abs(denp).gt.eps) then
            ii  = ii + 1
            azz(ii) = -(alphap*xb(1)+betap*xb(2)+gamap)/denp
            xx  = xb(1) + azz(ii)*wb(1)
            yy  = xb(2) + azz(ii)*wb(2)
            zz  = xb(3) + azz(ii)*wb(3)
            val = (xx-fpd(k+2))**2+(yy-fpd(k+3))**2
            val = (sqrt(val)-r1)**2/a2+(zz-fpd(k+4))**2/b2
            prod = -(xx-fpd(k+2))*betap+(yy-fpd(k+3))*alphap
            if(val.gt.1.0d0.or.prod.le.0.0) ii = ii - 1
          endif
        endif
      endif

*-----------------------------------------------------------------------
c* find the entrance and exit points of the torus
*-----------------------------------------------------------------------

 1282 continue
      if(ii.eq.3.and.ipq.eq.2) ii = 2
      if(ii.eq.2) then
        if(azz(1).ge.azz(2)) then
          rout = azz(1)
          rin  = azz(2)
        else
          rin  = azz(1)
          rout = azz(2)
        endif
      elseif( ii.eq.4) then
        if(azz(3).ge.azz(4)) then
          save = azz(3)
          azz(3) = azz(4)
          azz(4) = save
        endif
        if(azz(1).ge.azz(2)) then
          save = azz(1)
          azz(1) = azz(2)
          azz(2) = save
        endif
        if(azz(2).ge.azz(4)) then
          save = azz(2)
          azz(2) = azz(4)
          azz(4) = save
        endif
        if(azz(1).ge.azz(3)) then
          save = azz(1)
          azz(1) = azz(3)
          azz(3) = save
        endif
        if(azz(2).ge.azz(3)) then
          save = azz(2)
          azz(2) = azz(3)
          azz(3) = save
        endif
c*
c ************ 1996.2.14 *********************
        if( azz(2)-dist.gt.eps  ) then
c ************ 1996.2.14 *********************
          rout = azz(2)
          rin  = azz(1)
          fpd(k+11) = azz(3)
          fpd(k+12) = azz(4)
        else
          rout = azz(4)
          rin  = azz(3)
          fpd(k+11) = azz(1)
          fpd(k+12) = azz(2)
        endif
      elseif( ii.le.1.or.ii.eq.3 ) then
        go to 2000
      endif
      lri = 1
      lro = 1
      go to 2000


*-----------------------------------------------------------------------
c* elliptical cone
*-----------------------------------------------------------------------

 2200 continue
      fpd2  = fpd(k+2)**2
      fpd3  = fpd(k+5)**2
      fpd4  = fpd(k+3)**2+fpd(k+6)**2-fpd(k+8)**2
      fpd6  = fpd(k+5)*fpd(k+6)
      fpd7  = fpd(k+2)*fpd(k+3)
      fpd8  = fpd(k+2)*fpd(k+4)
      fpd9  = fpd(k+5)*fpd(k+7)
      fpd10 = fpd(k+3)*fpd(k+4)+fpd(k+6)*fpd(k+7)-fpd(k+8)*fpd(k+9)
      fpd11 = fpd(k+4)**2+fpd(k+7)**2-fpd(k+9)**2
      zz1   = fpd(k+10)
      zz2   = fpd(k+11)
      waw   = fpd2*wb(1)**2+fpd3*wb(2)**2+fpd4*wb(3)**2
     >      + 2.0*fpd6*wb(2)*wb(3)
     >      + 2.0*fpd7*wb(3)*wb(1)
      xax   = fpd2*xb(1)**2+fpd3*xb(2)**2+fpd4*xb(3)**2
     >      + 2.0*fpd6*xb(2)*xb(3)
     >      + 2.0*fpd7*xb(3)*xb(1)
      pw    = fpd8*wb(1)+fpd9*wb(2)+fpd10*wb(3)
      px    = fpd8*xb(1)+fpd9*xb(2)+fpd10*xb(3)
      axw   = fpd2*xb(1)*wb(1)+fpd3*xb(2)*wb(2)
     >      + fpd4*xb(3)*wb(3)
     >      + fpd6*(xb(3)*wb(2)+xb(2)*wb(3))
     >      + fpd7*(xb(1)*wb(3)+xb(3)*wb(1))
      aa    = waw
      bb    = axw+pw
      cc    = xax+2.0*px+fpd11

      dd = bb**2 - aa*cc
      if(abs(aa).ge.eps.and.dd.lt.0.0d0) then
        go to 2000
      else

*-----------------------------------------------------------------------
c* find the intersection points on the quadratic surface
*-----------------------------------------------------------------------

        if(abs(aa).ge.eps) then
          dd   =  sqrt(dd)
          rr1  = ( -bb - dd)/aa
          rr2  = ( -bb + dd)/aa
          zt1  = xb(3)+rr1*wb(3)
          zt2  = xb(3)+rr2*wb(3)
        else
          if( abs(bb).lt.eps ) go to 2000
          rr1  = -cc/bb/2.0
          zt1  = xb(3)+rr1*wb(3)
          zt2 = zz2 + 0.1
        endif
        ii   = 0
        if( zt1.ge.zz1.and.zt1.le.zz2) then
          ii = ii + 1
          azz(ii) = rr1
        endif
        if( zt2.ge.zz1.and.zt2.le.zz2) then
          ii = ii + 1
          azz(ii) = rr2
        endif
*-----------------------------------------------------------------------
c* find the intersection points on the planar surface , if any
*-----------------------------------------------------------------------
        if( abs(wb(3)).ge.eps ) then
          fk1 = (zz1-xb(3))/wb(3)
          fk2 = (zz2-xb(3))/wb(3)
          xx1 = xb(1)+fk1*wb(1)
          xx2 = xb(1)+fk2*wb(1)
          yy1 = xb(2)+fk1*wb(2)
          yy2 = xb(2)+fk2*wb(2)
          test1 = (fpd(k+2)*xx1+fpd(k+3)*zz1+fpd(k+4))**2
     >          + (fpd(k+5)*yy1+fpd(k+6)*zz1+fpd(k+7))**2
     >          - (fpd(k+8)*zz1+fpd(k+9))**2
          test2 = (fpd(k+2)*xx2+fpd(k+3)*zz2+fpd(k+4))**2
     >          + (fpd(k+5)*yy2+fpd(k+6)*zz2+fpd(k+7))**2
     >          - (fpd(k+8)*zz2+fpd(k+9))**2
          if( test1.le.0.d0) then
            ii = ii + 1
            azz(ii) = fk1
          endif
          if( test2.le.0.d0) then
            ii = ii + 1
            azz(ii) = fk2
          endif
        endif
c*
        if(ii.le.1) go to 2000
        rout = dmax1(azz(1),azz(2))
        rin  = dmin1(azz(1),azz(2))
        lri  = 1
        lro  = 1
        go to 2000
      end if

*-----------------------------------------------------------------------
c     transfer here for all body types  *   *   * *
*-----------------------------------------------------------------------

 2000    ma(locat+3)=lri
         ma(locat+4)=lro

      if( dabs(rin-dist)  .lt. dist * eps ) rin=dist
      if( dabs(rout-dist) .lt. dist * eps ) rout=dist


      if( rin - rout ) 2003,2001,2002

 2001 rout = rout * 1.00001

      go to 2003

 2011 write(ioe ,2010) itype,ir,nbo
 2010 format(  13h in gg itype=,i5,5x,  3hir=,i5,5x,  4hnbo=,i5)

      call pr(ld,1)
      call abend(ld)

 2002 rout = -pinf

 2003 fpd(k)   = rin
      fpd(k+1) = rout


      return
      end


*-----------------------------------------------------------------------

      subroutine gtvlin(fpd,mriz,vnor,ivopt,nir,niz,i1,i0)

*-----------------------------------------------------------------------
c **  this routine reads or calculates volumes of each region in geom
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      dimension fpd(*),mriz(*),vnor(*)
c     assumes there are nir regions
c     i1 is input tape
c     i0 is output tape
c     ivopt is option for calculating volumes
c     ivopt = 0  set volume = 1
c     ivopt = 1  concentric spheres
c     ivopt = 2  slabs
c     ivopt = 3  read in volumes
      write(i0,1001) ivopt,nir
 1001 format(1h1,7h option,i2,37h was used in calculating volumes, for ,
     1 i3,8h regions,/,68h 0-set volumes = 1, 1-concentric spheres, 2-sl
     2abs, 3-inputvolumes.                      )
      iv=ivopt+1
      go to(500,200,400,100),iv
  100 read(i1     )(vnor(i),i=1,nir)
 1002 format(7e10.5)
      go to 700
  500 do 50 i=1,nir
   50 vnor(i)=1.0
      go to 700
c   con=4*pi/3
  200 con=4.18879
      j=1
      vols=0.0
      do 250 i=1,niz
      if(i.eq.niz) go to 225
      if(mriz(i).eq.mriz(i+1)) go to 250
c bug fixed ???? at scale-4.3 marslib
  225 js=(i-1)*6+6
      volb=con*fpd(js)**3
      vnor(j)=volb-vols
      vols=volb
      j=j+1
  250 continue
      jf=j-1
      if(jf.eq.nir) go to 700
      write(i0,1250) jf,nir
 1250 format(48h0**********error in volume calculation**********/
     1   10h      jf =,i3,10h     nir =,i3      )
      call errtra
      return
 400  continue
  600 continue
  700 write(i0,1070)
 1070 format(1h ,//,10x,71h volumes ( cc ) used in collisions density an
     1d track length estimators.    )
      iln=(nir-1)/10
      js=-9
 701  js=js+10
      jf=js+9
      if(jf.gt.nir) jf=nir
      iln=iln-1
      write(i0,1071) (j,j=js,jf)
 1071  format(1h ,6h   reg , 10(3x,i5,4x))
      write(i0,1072)(vnor(j),j=js,jf)
 1072 format(7h volume,10(1x,1pe10.3,1x)/)
      if(iln.ge.0) go to 701
      return
      end
cmit
      function ipack(b1,b2,b3)

c * * this routine called by jomin1 to pack 3 characters into 1 word

      logical*1  l(4),l1
      integer c,d,b1,b2,b3,blnk
      equivalence (c,l1), (d,l(1))
      character*4 cblnk
      equivalence (cblnk,blnk)
      data cblnk/ 4h    /
      d=blnk
      c=b1
      l(1)=l1
      c=b2
      l(2)=l1
      c=b3
      l(3)=l1
      ipack=d
      return
      end

************************************************************************
*                                                                      *
      subroutine jomin1(buf,iibuf,istr,nadd,lim)
*                                                                      *
*     this routine reads free-form combinatorial input for             *
*     body data, zone description, region, universe & media nos.       *
*                                                                      *
************************************************************************


      character*3 lor,iz,ity,itype,ibl
      character*4 title,aread
      character*1 prefix
      dimension ity(23),ibias(23),buf(*),iibuf(*),title(15)
      common/jomk/nbi,iri,jmk,iml,iop,imax
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)

*-----------------------------------------------------------------------

      double precision buf,dread

      character chinr*60
*-----------------------------------------------------------------------

      logical lread

      parameter ( ispcd = 2 )

      data ity/'arb','sph','rcc','rec','trc'
     1        ,'ell','box','wed','rpp','gel'
     1        ,'tor','qua','bpp','wpp','p  '
     1        ,'px ','py ','pz ','ps ','c  '
     1        ,'cx ','cy ','cz '/

      data ibias/ 30,    4,    7,   12,    8
     1        ,    7,   12,   12,    6,   12
     1        ,    9,   10,    9,    9,    4
     1        ,    1,    1,    1,    7,    7
     1        ,    4,    4,    4/

*-----------------------------------------------------------------------

      data i4/4/,m1/-1/
      data ibl,lor/'   ','or '/

      call ionums(intt,iot,i,j)

*-----------------------------------------------------------------------
*    ionums sets up i/o tape unit # for dread calls
*-----------------------------------------------------------------------

      j1     = 72

*-----------------------------------------------------------------------
*     rcrdln sets record length to j1
*-----------------------------------------------------------------------

      call rcrdln(j1,j2)
      iret   = 0
      if (jmk.eq.2) iout   = 15

*-----------------------------------------------------------------------
*     jmk=2  means keno 4-cg data
*-----------------------------------------------------------------------

      if (istr+33.gt.lim) go to 300

*-----------------------------------------------------------------------
*     read title card
*     write binary title on temp file
*-----------------------------------------------------------------------

      read (intt,10300) (title(i),i=1,15)
      write(iou2)      (title(i),i=1,15)

*-----------------------------------------------------------------------
*    read options card
*-----------------------------------------------------------------------

      call scanon

      ivopt  = dread(1,iret)
      if(ivopt.le.-1) ivopt=0

      idbg   = dread(0,iret)
      ibod   = dread(0,iret)
      naz    = dread(0,iret)

      lfpd   = 0
      ibodt  = 0
      numb   = 0
      irtru  = 0
      numr   = 0

*-----------------------------------------------------------------------
*     read body definition card
*-----------------------------------------------------------------------

  100 if ( iret.eq.0 ) then

        itype  = aread(i4,iret)

        do 110 j=1,23

        if (ity(j).eq.itype) go to 120

  110   continue

        write(ioe,10400) itype,ity

10400   format(/' *** error message from s.jomin1 ***'
     &  /' geometry type did not equal any of the following type.'
     &  /' itype (input) = ',a3
     &  /' available type = ',20(a3,2x))

        call errtra

        return

  120   continue
        numb   = numb+1
        ialp   = numb
        if (ibod.gt.0)      ialp   = dread(0,iret)
        buf(1) = j
        buf(2) = ialp
        n      = 2+ibias(j)
        nn     = 6
        if (ibias(j).gt.nn) nn     = ibias(j)
        lfpd   = lfpd+2+nn
        if (j.eq.1)         lfpd   = lfpd+26
        if (j.eq.11)        lfpd   = lfpd+2

        do 130 k=3,n
          buf(k) = dread(0,iret)
  130   continue

        if(j.ge.13.and.j.le.14) then

          j      = j-3
          buf(1) = j
          n      = 14
          lfpd   = lfpd+3
          call abox(buf(3))
        end if

*-----------------------------------------------------------------------
*       write binary body description block on temp file
*-----------------------------------------------------------------------

        write(iou2) n,(buf(i),i=1,n)
        go to 100
      end if
      call scanof
      itype  = aread(i4,iret)
      call rcrdln(j2,j1)
      iret   = 0
      call rstptr(j2)
      call scanon
  140 continue
      ltma   = 7*numb

*-----------------------------------------------------------------------
*     read zone description card
*-----------------------------------------------------------------------

      kn     = 3

  150 if ( iret.eq.0 ) then

        if ( .not.lread(i4,iret) ) then

          prefix=aread(m1,iret)
          iz=aread(i4,iret)
          if (prefix.ne.' ') then
            itype = prefix//iz
          else
            itype = iz
          end if
          iz = itype
          if (iz.eq.lor) then
            if (ibod.eq.0) go to 150
            if (n.eq.2)    go to 170
              iibuf(kn) = ibod
              ibodt  = ibodt+ibod
            go to 170
          end if

*-----------------------------------------------------------------------
*         a new input zone
*         write old input zone description block on temp file
*-----------------------------------------------------------------------

          if (irtru.eq.0)  go to 160
          iibuf(2)  = irts
          iibuf(kn) = ibod
          ibodt  = ibodt+ibod
          kn     = 3
          write(iou2) n,(iibuf(j),j=1,n)
  160     continue
          irts   = 0
          read (iz,'(a3)') iibuf(1)
          irtru  = irtru+1
          n      = 2

*-----------------------------------------------------------------------
*         check to see if first thing is an or
*-----------------------------------------------------------------------

  170     continue

*-----------------------------------------------------------------------
*         irts is number of code zones in this input zone
*-----------------------------------------------------------------------

          irts   = irts+1

*-----------------------------------------------------------------------
*         numr is total number of code zones
*-----------------------------------------------------------------------

          numr   = numr+1
          ibod   = 0
          n      = n+1
          kn     = n

        else

*-----------------------------------------------------------------------
*         should be a number
*         ibod is number of bodies in this code zone
*-----------------------------------------------------------------------

          ibod   = ibod+1
          n      = n+1
          if (istr+n+1.gt.lim) go to 300
          iibuf(n) = dread(0,iret)

        end if

        go to 150

      end if

*-----------------------------------------------------------------------
*     end of zone description blocks
*-----------------------------------------------------------------------

      call scanof
      itype  = aread(i4,iret)
      iret   = 0
      iibuf(2)  = irts
      iibuf(kn) = ibod
      ibodt  = ibodt+ibod
      write(iou2) n,(iibuf(j),j=1,n)
      ltma   = ltma+numr+5*ibodt
      rewind iou2
      rewind iout
      read(iou2 ) (title(i),i=1,15)
      write(iout) (title(i),i=1,15)
      nazt   = naz+5*irtru
      ltma   = ltma+3*nazt
      nadd   = 3+ltma+lfpd*ispcd+7*numr+7*irtru

*-----------------------------------------------------------------------
*    nadd is now the space used by cg data in jomin1
*-----------------------------------------------------------------------

      if (istr+nadd+1.gt.lim) go to 300
      iibuf(1)  = ivopt
      iibuf(2)  = idbg
      iibuf(3)  = irtru
      iibuf(4)  = numr
      iibuf(5)  = numb
      iibuf(6)  = nadd
      iibuf(7)  = nazt
      iibuf(8)  = ibodt
      iibuf(9)  = ltma
      iibuf(10) = lfpd

      write(iout) (iibuf(i),i=1,10)

*-----------------------------------------------------------------------
*     read region numbers
*-----------------------------------------------------------------------

      iibuf(1) = dread(1,iret)
      if (irtru.gt.1) then
        do 180 i=2,irtru
          iibuf(i)=dread(0,iret)
  180     continue
      end if

      write(iout) (iibuf(i),i=1,irtru)

*-----------------------------------------------------------------------
*     to calculate the number of input regions if ivopt=3
*-----------------------------------------------------------------------

      if (ivopt.eq.3) then
        nir    = 1
        do 200 i=2,irtru
          m      = i-1
          do 190 j=1,m
            if (iibuf(i).eq.iibuf(j)) go to 200
  190       continue
          nir    = nir+1
  200     continue
      end if

*-----------------------------------------------------------------------
*    read geometry box-space  type  ...
*-----------------------------------------------------------------------

      iibuf(1) = dread(1,iret)
      if (irtru.gt.1) then
        do 210 i=2,irtru
          iibuf(i) = dread(0,iret)
  210     continue
      end if

      write(iout) (iibuf(i),i=1,irtru)

*-----------------------------------------------------------------------
*     read media numbers
*-----------------------------------------------------------------------

      iibuf(1) = dread(1,iret)
      if (irtru.gt.1) then

        do 220 i=2,irtru
          iibuf(i)=dread(0,iret)
  220     continue
      end if

      write(iout) (iibuf(i),i=1,irtru)

*-----------------------------------------------------------------------
*     read temp blocks and put on regular file
*-----------------------------------------------------------------------

      do 230 j=1,numb
        read(iou2 ) l,(buf(i),i=1,l)
        write(iout) l,(buf(i),i=1,l)
  230   continue

      do 240 j=1,irtru
        read(iou2 ) l,(iibuf(i),i=1,l)
        write(iout) l,(iibuf(i),i=1,l)
  240   continue

      rewind iou2

*-----------------------------------------------------------------------
*       to read volumes of input regions if ivopt eq 3
*-----------------------------------------------------------------------

      if (ivopt.eq.3) then

        buf(1) = dread(1,iret)
        do 250 i=2,nir
          buf(i) = dread(0,iret)
  250     continue
        write(iout) (buf(i),i=1,nir)

      end if

      return

*-----------------------------------------------------------------------
*     if exceed size of buf array
*-----------------------------------------------------------------------

  300 continue


10300 format(15a4)
      write(ioe,10000) istr,n,nadd,lim,ltma,lfpd,numr,irtru,nazt,ibodt

10000 format(/' *** error message from s.jomin1 ***'
     &/' storage size of cg-geometry data was exceeded the limit.'
     &/' istr =',i8,'  n    =',i8,'  nadd =',i8,'  lim  =',i8
     &/' ltma =',i8,'  lfpd =',i8,'  numr =',i8,'  irtru=',i8,
     &'  nazt =',i8,'  ibodt=',i8)

      call errtra

      return
      end


*-----------------------------------------------------------------------

      subroutine level(nlv,ncmax,ni4,nbod,d,ld,kbz,mcz,numr)

*-----------------------------------------------------------------------

c *** this routine determines the highest level number an array is
c     referenced in and updates nlv table.  used on input only.

      implicit real*8 (a-h,o-z)

      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)

      dimension nlv(*),ncmax(*),ni4(*),nbod(*),d(*),ld(*),kbz(*),mcz(*)

      nlev=0
      ncl=1
  10  continue
      isp=0
      do 100 n=1,nar
      if(nlv(n).ne.ncl) go to 100
      m=3*(n-1)
      mx=ncmax(m+1)
      my=ncmax(m+2)
      mz=ncmax(m+3)
      nr=ni4(n)
      do 80 k=1,mz
      do 80 j=1,my
      do 80 i=1,mx
      ix=my*mx*(k-1)+mx*(j-1)+i
      lx=ix+nr-1
      lm=ld(lx)
      if(lm) 30,80,60
  30  jm=iabs(lm)
      nlv(jm)=ncl+1
      isp=1
      go to 80
  60  continue
      do 70 l=1,numr
      if(kbz(l).ne.lm) go to 70
      if(mcz(l).ge.0) go to 70
      if(mcz(l).le.-1000) go to 70
      jm=iabs(mcz(l))
      nlv(jm)=ncl+1
      isp=1
  70  continue
  80  continue
 100  continue
      if(isp.eq.0) go to 200
      ncl=ncl+1
      go to 10
 200  nlev=ncl
      return
      end

************************************************************************
*                                                                      *

      subroutine lookz(x,y,z,nsto ,ma,fpd,locreg,numbod,iror,nsor)


c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     92/06/22
c     time the module was last permanently updated:     08:55:32
c     programmer name:                                  j.t.west
c     module name:                                      malookz
c     current archiving level number:                   00004
c     current number of permanent updates:              00004
c     date of last access by librarian:                 93/03/15
c     dataset name:  x4s.scale4.master
c
c#######################################################################

c      this routine returns the combinatorial geometry zone of point  *
c          (x,y,z)  so tracking can be initialized.

      implicit real*8 (a-h,o-z)

      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)

      dimension ma(*),fpd(*),locreg(*),numbod(*),iror(*),nsor(*),nsto(*)
      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      real*8 dist0
      integer blzold
      common /orgi/ dist0,markg,nmedg,nblz,blzold,irpold
!$OMP THREADPRIVATE(/orgi/)
      real*8 xb,wb,wp,xp,rin,rout,pinf,dist
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,ir,idbg ,
     1             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)
      common /dbg/ n,num,locat,isave,inext,irp,smin,inex
!$OMP THREADPRIVATE(/dbg/)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)

      data eps /1.0d-5/
      save eps !FURUTA

*-----------------------------------------------------------------------

         if( ll .eq. 0 )     goto 10
         if( ll .ne. llold ) goto 20

 10   continue

         if( x .eq. xb(1) .and.
     &       y .eq. xb(2) .and.
     &       z .eq. xb(3) ) return

  20     llold = ll
         kloop = kloop + 1

*-----------------------------------------------------------------------

      xb(1) = x
      xb(2) = y
      xb(3) = z

      dist = 0.0d0

*-----------------------------------------------------------------------

      do 500 is = 1, numr

*-----------------------------------------------------------------------

               irp = nsor(is)

                  if( irp .eq. 0 ) goto  600

                  if( ll .ne. nsto(kbcz+irp-1) ) goto 500

               n   = locreg(irp) + 1
               num = numbod(irp) * 5 + n - 5

                  if(idbg.ne.0) call pr(nsto,4)

*-----------------------------------------------------------------------
*        the loop to 400 examines region irp to see
*        if it is the next region
*-----------------------------------------------------------------------

         do 400 i = n, num, 5

*-----------------------------------------------------------------------

               nbo   = ma(i)
               locat = ma(i+1)

            call gg(locat,ma,fpd,nsto)

               if(ierror.ne.0) return

            dlta = eps * dist

               if(idbg.ne.0) call pr(nsto,5)

            if(nbo) 320,400 ,330
  320       if((rout-dist).le.dlta.or.(rin-dist).gt.dlta) go to 400

            goto 500

  330       if( (rin-dist).le.dlta.and.(dist.lt.rout) )  go to 400

            goto 500

  400    continue

*-----------------------------------------------------------------------

            goto 750

  500 continue

*-----------------------------------------------------------------------

  600 continue

      inext = is
            if(idbg.ne.0) call pr(nsto,6)

      do  700 irp=1,numr

         if(ll.ne.nsto(kbcz+irp-1)) go to 700
         n = locreg(irp)+1
         num = numbod(irp)*5+n-5

         do  650  i = n, num, 5

            nbo=ma(i)
            locat=ma(i+1)

                  call gg(locat,ma,fpd,nsto)

                     if(ierror.ne.0) return

                  dlta = eps * dist

            if(idbg.ne.0) call pr(nsto,7)
            if(nbo) 620,650,630
  620       if((rout-dist).le.dlta.or.(rin-dist).gt.dlta) go to 650
            go to 700
  630       if((rin-dist).le.dlta.and.(dist.lt.rout)) go to 650
            go to 700

  650    continue

         nsor(inext) = irp

         goto 750

  700 continue

      iect = iect + 1

*-----------------------------------------------------------------------
*     for undefined region searching
*-----------------------------------------------------------------------

      ierror = 2

      if( iect .gt. 10 ) return

      write(ioe ,1000) ir,xb,wb,dist

 1000 format(' *** warning in lookz *** region(',i3,')  at point(',
     : 3(1pe10.3,','),')'/
     : ' next region will be out of geometry or undefined '/
     : '  particle direction =(',3(e10.3,','),')  '/
     : '  particle will transport through distance:',e12.4)


      irprim =0
      call abend(nsto)

      return

*-----------------------------------------------------------------------
*     irprim is the combinatorial input zone number - ir is code zone
*     this is consistent with definition in g1
*-----------------------------------------------------------------------

  750 continue

      irprim = iror(irp)
      ir = irp

         if(idbg.ne.0) call pr(nsto,8)

      nmedg = irprim
      nblz = 32768 * irp + irprim


      return
      end

*-----------------------------------------------------------------------

      function medno(ld,ir)

      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      dimension  ld(1)
      medno=ld(kmiz+ir-1)
      return
      end

*-----------------------------------------------------------------------

      integer function nregno(ld,ir)

c ----------------------------------------------------------------------
      common/gomloc/kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     &              kkr2,knsr,kvo,nadd,ldata,ltma,lfpd,numr,irtru,numb,
     &              nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      dimension ld(1)
c ----------------------------------------------------------------------
c
      nregno=ld(kriz+ir-1)
      return
      end

      subroutine norml(ma,fpd,locreg,numbod,ir1,ir2,un)
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     92/06/22
c     time the module was last permanently updated:     08:57:32
c     programmer name:                                  j.t.west
c     module name:                                      manorml
c     current archiving level number:                   00002
c     current number of permanent updates:              00002
c     date of last access by librarian:                 93/03/15
c     dataset name:  x4s.scale4.master
c
c#######################################################################

c      this routine calculates the unit vector to combgeom body nasc
c      *    *       *       *    at point xb+wb*dist    *    *     *

      implicit real*8 (a-h,o-z)
      dimension ma(*),fpd(*),locreg(*),numbod(*) ,ir1(*),ir2(*)
      dimension p(3),f1(3),f2(3)

      dimension un(3)

      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      real*8 xb,wb,wp,xp,rin,rout,pinf,dist
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,ir,idbg ,
     1             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)
      dimension x(3),h(3),ibox1(3),ibox2(3)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      data ibox1/10,7,4/
      data ibox2 /4,10,7/

*-----------------------------------------------------------------------
c  user beware       it is assumed that
c      xb+wb*dist  is on surface lsurf of body nasc.thus a normal call
c      must be proceeded by a call to g1 or lookz.
*-----------------------------------------------------------------------

      lsur=iabs(lsurf)
      do 25 i=1,3
   25 xp(i)=xb(i)+dist*wb(i)

      jr1= ir1(ir)
      jr2= ir2(ir)

      do 52 irr=jr1,jr2
      n=locreg(irr)+1
      num=numbod(irr)*5+n-5
      do 50 i=n,num,5
      if(nasc.eq.iabs(ma(i)))   go to 100
   50 continue
  52  continue

      write(ioe,1000) ir,nasc,(ma(i),i=n,num,5)
 1000 format(40h1invalid region or body in normal.   ir=,i7,5hnasc=,i7,
     .  /(10i7))
      iadr = 2-kma
      call abend(ma(iadr))
      return

  100 l=ma(i+1)
      itype=ma(l+2)
      k=ma(l+6)
c*mit******* '91/4/25
c            arb  sph  rcc  rec  trc  ell  box  wed  rpp
      go to (1100,1200,1400,1400,1400,1600,1700,1800,1900,
ctat 1998/2/26 by H.Row.Tateishi---from
ctat         gel  tor  qua
ctat :       2100,1400,1400),itype
c
c            gel  tor  qua  bpp  wpp  p    px   py   pz   ps
     :       2100,1400,1400,9999,9999,3000,3100,3200,3300,3400,
c            c    cx   cy   cz
     :       4000,4100,4200,4300),itype
ctat 1998/2/26 by H.Row.Tateishi---to
c*mit*******
c  arb
 1100 i=k+1+4*(lsur -1)
      do 1110 j=1,3
 1110 un(j)=fpd(i+j)
      go to 5000
c  sph
 1200 do 1210 i=1,3
 1210 un(i)=(xp(i)-fpd(k+1+i))/fpd(k+5)
      go to 5000
c
c  trc,rcc and rec
c
ctat 1998/2/26 by H.Row.Tateishi---from
 1400 continue
c
c process for tor & qua
      if(itype.eq.11.or.itype.eq.12) then
         do 1401 i=1,3
1401     un(i)=wb(i)
         goto 5000
      endif
c
      h2=0.0d0
ctat 1998/2/26 by H.Row.Tateishi---to
      do 1410 i=1,3
      h(i)=fpd(k+4+i)
 1410 h2=h2+h(i)**2
      go to (1415,1415,1435), lsur
 1415 h2=dsqrt(h2)
      do 1420 i=1,3
 1420 un(i)=h(i)/h2
      go to 5000
 1435 do 1440 i=1,3
 1440 x(i)=xp(i)-fpd(k+1+i)
      xmu=0.0d0
      a2=0.0d0
      do 1445 i=1,3
      j=mod(i,3)+1
      a2=a2+(h(i)*x(j)-h(j)*x(i))**2
 1445 xmu=xmu+fpd(k+4+i)*x(i)
      if(itype.eq.4) go to 1300
      r=fpd(k+8)-fpd(k+9)
      if(itype.eq.3) r=0.0d0
      b1=r/dsqrt(h2*(h2+r**2))
      b2=1.0d0/dsqrt(a2*(h2+r**2))
      do 1450 i=1,3
 1450 un(i)=b2*(x(i)*h2-h(i)*xmu)+b1*h(i)
      go to 5000
c
c   the following coding is for  surface 3 of the rec *****
 1300 ada=0.0d0
      bdb=0.0d0
      do 1305 i=1,3
      ada=ada+fpd(k+7+i)**2
 1305 bdb=bdb+fpd(k+10+i)**2
      c = dsqrt(ada-bdb)
c  c is the focus length  ****
      a=dsqrt(ada)
      do 1310 i=1,3
 1310 p(i) = c*fpd(k+7+i)/a
      twoa=a+a
      e=0.0d0
      do 1315 i=1,3
      temp=fpd(k+1+i)+h(i)*xmu/h2
      f1(i) =temp+p(i)
      f2(i) =temp-p(i)
 1315 e=e+(xp(i)-f1(i))**2
      e=dsqrt(e)
      e1= 0.0d0
      do 1320 i=1,3
      p(i) = f1(i)+twoa*(xp(i)-f1(i))/e
 1320 e1= (f2(i)-p(i))**2+e1
      e1= dsqrt(e1)
      do 1325 i=1,3
 1325 un(i) = (f2(i)-p(i))/e1
      go to 5000
c
c  ell
 1600 dx=0.0d0
      an=0.0d0
      c2=fpd(k+8)**2/4.0d0
      do 1620 i=1,3
      ia=i+k+1
      x(i)=xp(i)-.5*(fpd(ia)+fpd(ia+3))
      h(i)=.5*(fpd(ia+3)-fpd(ia))
 1620 dx=dx+x(i)*h(i)
      do 1630 i=1,3
      un(i)=c2*x(i)-dx*h(i)
 1630 an=an+un(i)**2
      if(an.lt.1.0d-12) go to 999
      an=dsqrt(an)
      do 1635 i=1,3
 1635 un(i)=un(i)/an
      go to 5000
  999 write(ioe,1002) nbo,lsurf,(xp(i),i=1,3)
 1002 format(31h round-off error in normal nbo=,i5,6hlsurf=,i5,3hxp=,
     1     3d12.5)
      iadr = 2-kma
      call abend(ma(iadr))
      return
c
c  box
 1700 l= 1+(lsur-1)/2
      i1=ibox1(l)
      i2=ibox2(l)
      do 1720 i=1,3
      x(i)=fpd(k+i1+i)
 1720 h(i)=fpd(k+i2+i)
 1721 an=0.0d0
      do 1730 i=1,3
      j=mod(i,3)+1
      jj=mod(i+1,3)+1
      un(i)=x(j)*h(jj)-x(jj)*h(j)
 1730 an=an+un(i)**2
      if(an.lt.1.0d-12) go to 999
      an=dsqrt(an)
      do 1735 i=1,3
 1735 un(i)=un(i)/an
      go to 5000
c
c wed or raw
 1800 if(lsur .ne.2) go to 1700
      do 1810 i=1,3
      x(i) = fpd(k+10+i)
 1810 h(i)= fpd(k+4+i)-fpd(k+7+i)
      go to 1721
c
c  rpp
 1900 do 1905 i=1,3
 1905 un(i)=0.0d0
      go to (1910,1910,1920,1920,1930,1930), lsur
 1910 un(2)=1.0d0
      go to 5000
 1920 un(1)=1.0d0
      go to 5000
 1930 un(3)=1.0d0
      go to 5000
c
c  gel
ctat 1998/2/26 by H.Row.Tateishi---from
2100  continue
      do 2110 i=1,3
2110  un(i)=wb(i)
ctat 1998/2/26 by H.Row.Tateishi---to
      go to 5000
c
ctat 1998/2/26 by H.Row.Tateishi---from
c
c    p
3000  continue
      do 3001 i=1,3
3001  un(i)=fpd(k+1+i)
      an=0.
      do 3002 i=1,3
3002  an=an+un(i)**2.
      if(an.lt.1.0d-12) go to 999
      an=dsqrt(an)
      do 3003 i=1,3
3003  un(i)=un(i)/an
      goto 5000
c
c   px
c   py
c   pz
3100  continue
3200  continue
3300  continue
      do 3301 i=1,3
3301  un(i)=0.
      if(itype.eq.16) un(1)=1.
      if(itype.eq.17) un(2)=1.
      if(itype.eq.18) un(3)=1.
      goto 5000
c
c  ps
3400  continue
      do 3401 i=1,3
3401  un(i)=fpd(k+4+i)
      an=0.
      do 3402 i=1,3
3402  an=an+un(i)**2.
      if(an.lt.1.0d-12) go to 999
      an=dsqrt(an)
      do 3403 i=1,3
3403  un(i)=un(i)/an
      goto 5000
c
c   c
4000  continue
      do 4001 i=1,3
      x(i)=fpd(k+1+i)
4001  h(i)=fpd(k+4+i)
      t=(h(1)*xp(1)+h(2)*xp(2)+h(3)*xp(3)
     1 -(h(1)*x (1)+h(2)*x (2)+h(3)*x (3)))
     1 /(h(1)**2.  +h(2)**2.  +h(3)**2.)
      do 4002 i=1,3
4002  un(i)=xp(i)-(x(i)+t*h(i))
      do 4004 i=1,3
4004  un(i)=un(i)/fpd(k+8)
      goto 5000
c
c  cx
c  cy
c  cz
4100  continue
4200  continue
4300  continue
      un(1)=xp(1)-fpd(k+2)
      un(2)=xp(2)-fpd(k+3)
      un(3)=xp(3)-fpd(k+4)
      if(itype.eq.21) un(1)=0.
      if(itype.eq.22) un(2)=0.
      if(itype.eq.23) un(3)=0.
      do 4302 i=1,3
4302  un(i)=un(i)/fpd(k+5)
      goto 5000
ctat 1998/2/26 by H.Row.Tateishi---to
c
c
c
 5000 xmu=0.0d0
      do 5005 i=1,3
 5005 xmu=xmu+wb(i)*un(i)
      an=dsign(1.0d0,-xmu)
      do 5010 i=1,3
 5010 un(i)=an*un(i)
ctat 1998/2/26 by H.Row.Tateishi---from
 9999 continue
ctat 1998/2/26 by H.Row.Tateishi---to
      return
      end
*-----------------------------------------------------------------------

      subroutine orthom(d,ld,ma,fpd,nasc1,tr)

*-----------------------------------------------------------------------
c **  this routine computes the rotational matrices of each box
c     which references an array or a universe. rota uses results.
      implicit real*8 (a-h,o-z)

      real*4 d

      dimension d(*),ld(*),ma(*),fpd(*),tr(3,3),dtr(3)

      mn=ma(nasc1*7)
      do 20 j1=1,3
      mo=mn+4+3*(j1-1)
      dtr(j1)=dsqrt(fpd(mo+1)**2+fpd(mo+2)**2+fpd(mo+3)**2)
      do 20 j2=1,3
      tr(j2,j1)=fpd(mo+j2)/dtr(j1)
  20  continue
      return
      end



      subroutine qrtic(c,r,n)
c
c    solves a polynomial equation of the type
c      x**4 + c(1)*x**3 + c(2)*x**2 + c(3)*x + c(4) = 0
c      using ferrari's solution of the quarti equation
c    the coefficient of x**4 is assumed to be 1
c    r(4) contains the roots
c    n contains the number of real roots
c    if there are 2 real roots they will be in r(1) and r(2)
c      with the complex roots r(3) +- r(4)*i
c    if there arare no real roots, the complex roots are
c      r(1) +- r(2)*i and r(3) +- r(4)*i
c    qrtic calls cubic to find the roots of
c      the resolvent cubic equation
c
      implicit real*8 (a-h,o-z)
      dimension c(4),cc(3),rr(3),r(4)
c   resolvent cubic
c
      c1sq=c(1)*c(1)
      cc(1)=-.5*c(2)
      cc(2)=.25*c(1)*c(3)-c(4)
      cc(3)=.125*( c(4)*(4.*c(2)-c1sq)-c(3)*c(3) )
      call cubic(cc,rr,nn)
c
      t=.25*c1sq-c(2)
      do 10 i=1,nn
      root=rr(i)
      asq=t+root+root
      if(abs(asq).le.0.0001d0)asq=0.
      if(asq.lt.0.0)goto 10
      bsq=root*root-c(4)
      if(abs(bsq).le.0.0001d0)bsq=0.
      if(bsq.ge.0.0)goto 20
   10 continue
      n=0
      return
c
   20 twoab=c(1)*root-c(3)
      a=sqrt(asq)
      b=sign(sqrt(bsq),twoab)
      n=0
      real=.25*(a+a-c(1))
      disc=real*real-root+b
      if(abs(disc).le.0.0001d0)disc=0.
      sqroot=sqrt(abs(disc))
      if(disc.lt.0.0)goto 30
c        disc.ge.0    2 real rootsc
c
      n=2
      r(1)=real+sqroot
      r(2)=real-sqroot
      goto 40
c        disc.lt.0     2 imaginary roots
   30 r(3)=real
      r(4)=sqroot
c
   40 real=real-a
      disc=real*real-root-b
      if(abs(disc).le.0.0001d0)disc=0.
      sqroot=sqrt(abs(disc))
      if(disc.lt.0.0)goto 50
c        disc.ge.0    2 real roots
c
      n=n+2
      r(n)=real-sqroot
      r(n-1)=real+sqroot
      return
c        disc.lt.0    2 imaginary roots
   50 r(n+1)=real
      r(n+2)=sqroot
      return
      end

*-----------------------------------------------------------------------

      subroutine reset(d,ld,nl,nlu)

*-----------------------------------------------------------------------
* * * this routine resets array tracking variables which may be
*     destroyed by point detector estimation. resets nesting table.
*-----------------------------------------------------------------------

      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      common /repeat/ jp(20)
      dimension d(*),ld(*)

      ll  = 0
      lm  = 0
      nl  = 0
      nlu = 0

      if( nlev .le. 0 ) return

      in = jp(12)
      nl = ld(in+6*nlev)
      nlu = ld(in+6*nlev+2)
      if(nl.eq.0) return
      nby = ld(6*(nl-1)+in)
      lm = ld(6*(nl-1)+in+1)
      ll = lm

      do 10 i=1,3
      il = in+6*(nl-1)+i+1
      nx1(i) = ld(il)
  10  continue

      return
      end

*-----------------------------------------------------------------------

      subroutine resetz(ix,i1,i2)

*-----------------------------------------------------------------------

      dimension ix(1)

      do 100 i=i1,i2
      ix(i)=0
 100  continue

      return
      end

*-----------------------------------------------------------------------

      subroutine restor(d,ld,nunt,l,naad,ndsn)

c * * this routine reads geometry data from binary file when
c     data was preprocessed in scale.  called by jomin
      dimension d(*),ld(*)
      common/repeat/jp(20)
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      common/parem/ljy(44)
!$OMP THREADPRIVATE(/parem/)
      dimension nqb(10),kaf(23)
      parameter ( ispcd = 2 )
      if (l.gt.0) go to 100
      read(nunt) jp,
     &     nby,nlev,nar,nq,iaw,iay,nf,nx1,                   !FURUTA
     &     kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1, !FURUTA
     1     kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,   !FURUTA
     2     numb,nir,kbiz,kbcz,                               !FURUTA
     &     ljy                                               !FURUTA
      ibod=jp(8)-jp(7)
      mus=jp(14)+2*ibod+2
      muz= mus+ibod
      naad=nadd
  100 continue
      if (l.le.0) return
      kma=l
      kma=2*(kma/2)+1
      if(kma .eq. l) go to 110
      nadd=nadd+1
      naad=nadd
  110 kfpd=kma+ltma
c......fpd array must start on double precision word boundary for ibm
      kfpd=2*(kfpd/2)+1
      klcr=kfpd+ispcd*lfpd
      knbd=klcr+numr
      kior=knbd+numr
      kriz=kior+numr
      krcz=kriz+irtru
      kbiz=krcz+numr
      kbcz=kbiz+irtru
      kmiz=kbcz+numr
      kmcz=kmiz+irtru
      kkr1=kmcz+numr
      kkr2=kkr1+irtru
      knsr=kkr2+irtru
      kvol=knsr+numr
      kvol=2*(kvol/2)+1
      n3=kvol-jp(1)+nir*ispcd
      n4=l+nadd-1
      m2=16
      if (ndsn.gt.0) m2=20
      do 120 j=1,m2
  120    jp(j)=jp(j)+n3
      read(nunt) (d(i),i=kma,n4)
      ll=0
      if (nlev.le.0) go to 160
      do 130 i=1,nar
         ix=jp(5)+i-1
         ld(ix)=ld(ix)+n3
         iy=jp(8)+3*(i-1)
         ld(iy)=ld(iy)+n3
         ld(iy+1)=ld(iy+1)+n3
         ld(iy+2)=ld(iy+2)+n3
  130    continue
      do 140 i=1,numb
         iz=jp(10)+i-1
         if (ld(iz).le.0) go to 140
         ld(iz)=ld(iz)+n3
  140    continue
      mus=mus+n3
      muz=muz+n3
      if (ndsn.le.0) go to 160
      ilm=jp(18)-jp(17)
      do 150 i=1,ilm
         ix=jp(17)+i-1
         ld(ix)=ld(ix)+n3
  150    continue
  160 continue
      return
      end
      subroutine rota(xa,xc,iflop,tr)
c * * this routine rotates particle and direction cosine when
c     entering or exiting coordinate system defined by a box.
      implicit real*8 (a-h,o-z)
      real*8 xb,wb,wp,xp,rin,rout,pinf,dist
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,ir,idbg ,
     1             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)
      dimension xa(3),xc(3),tr(3,3),wa(3)
      do 20 i=1,3
      xa(i)=0.0
      wa(i)=0.0
      do 20 j=1,3
      l=i
      m=j
      if(iflop.eq.1) go to 10
      l=j
      m=i
  10  continue
      xa(i)=tr(l,m)*xc(j)+xa(i)
      wa(i)=tr(l,m)*wb(j)+wa(i)
  20  continue
      do 25 i=1,3
  25  wb(i)=wa(i)
      return
      end
      subroutine sazar(ip,np,d,ld,ncmax,ndsn,nresp,nmost)
c *** this routine reads array edit input if ndsn gt 0
      dimension ip(20),d(*),ld(*),ncmax(*)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      ip(nq+16)=ip(nq+15)+ndsn*5
      ip(nq+17)=ip(nq+16)+ndsn+1
      icr=1
      iret=0
      iy=ip(nq+16)
      ld(iy)=ip(nq+17)
      nsmax=0
      do 130 i=1,ndsn
         ix=ip(nq+15)+5*(i-1)
         do 100 j=1,5
            iy=ix+j-1
            ld(iy)=0
            if (ld(iy-1).eq.0.and.j.ge.2) go to 100
            ld(iy)=iread(icr,iret)
            icr=0
  100       continue
         na=ld(ix)
         if (na.lt.0) go to 110
         nsz=1
         go to 120
  110    nlu=iabs(na)
         nx=3*(nlu-1)+1
         nsz=ncmax(nx)*ncmax(nx+1)*ncmax(nx+2)
  120    if (nsz.gt.nsmax) nsmax=nsz
         iy=ip(nq+16)+i
         ld(iy)=ld(iy-1)+nsz
  130    continue
      ip(nq+18)=ld(iy)+1
      ip(nq+19)=ip(nq+18)+ndsn*2
      do 140 i=1,ndsn
         ix=ip(nq+18)+i-1
         ld(ix)=iread(icr,iret)
         ld(ix+ndsn)=iread(icr,iret)
  140    continue
      nf=nsmax
      do 150 i=1,ndsn
         ix=ip(nq+15)+5*(i-1)
         ix4=ix+4
         iy=ip(nq+18)+i-1
         write(iot,10000) i,(ld(j),j=ix,ix4),ld(iy),ld(iy+ndsn)
  150    continue
      np=np+4
      if (2*nf.lt.20*nmost) return
      nmost=(nf/7)+1
      write(iot,10100) nmost
      return
10000 format(/,5x,'sazar',i5,5x,5i5,5x,i5,5x,i5)
10100 format(//,5x,' nmost increased to ',i6,' by subroutine sazar',//)
      end
      subroutine sorcer(ma,locreg,numbod)
c * * this routine sorts next zone of entry table in order of
c     most probable entry.  speeds up tracking.
      implicit real*8 (a-h,o-z)
      dimension ma(*),locreg(*),numbod(*)
      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
c
      do 800 i=1,numr
      n=locreg(i)+1
      num=numbod(i)*5+n-5
c
      do 800 j=n,num,5
      ict=j+2
  80  continue
      if(ma(ict+2).eq.0) go to 800
      i2=ict
      i1=ma(ict+2)
 100  continue
      if(ma(i1+1).gt.ma(i2+1)) i2=i1
      i1=ma(i1+2)
      if(i1.gt.0) go to 100
c
      if(i2.eq.ict) go to 170
      izt=ma(ict)
      iet=ma(ict+1)
      iz2=ma(i2)
      ie2=ma(i2+1)
      ma(ict)=iz2
      ma(ict+1)=ie2
      ma(i2)=izt
      ma(i2+1)=iet
c
 170  continue
      ict=ma(ict+2)
      if(ict.gt.0) go to 80
c
 800  continue
      return
      end
      subroutine sort(xb,ncmax,lu,rby,nxy,d,ld,wb,m,ier,ifl)
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     92/06/22
c     time the module was last permanently updated:     09:04:15
c     programmer name:                                  j.t.west
c     module name:                                      masort
c     current archiving level number:                   00002
c     current number of permanent updates:              00002
c     date of last access by librarian:                 93/03/15
c     dataset name:  x4s.scale4.master
c
c#######################################################################
c *** this routine determines cell position on each dimension **
c *** also, determines when particle is exiting an array lattice
c *******************
      implicit real*8 (a-h,o-z)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      dimension xb(3),wb(3),ncmax(3),rby(*),nxy(3),d(*),ld(*)
      data sml/1.0e-09/,osml/1.0000055/
      nmax=ncmax(m)+1
      ifl=0
      i=0
      lx1=1
      lx2=nmax
      small= dmax1(sml,sml*rby(lx2))
      if(xb(m).le.-small.or.rby(lx2)*osml.lt.xb(m)) go to 20
      if(igx.eq.0) go to 60
      lx1=nxy(m)-1
      if(lx1) 20,50,55
  50  lx1=1
  55  lx2=nxy(m)+1
      if(lx2-nmax) 60,60,20
  60  continue
      small= dmax1(sml,sml*rby(lx2))
      do 10 i=lx1,lx2
      if(rby(i)-small.gt.xb(m)) go to 15
      if(rby(i).gt.xb(m)-small) go to 12
  10  continue
      i=lx2
      if(igx.eq.1) go to 20
  12  continue
      lu=m
      nxy(m)=i-1
      if(wb(m).ge.0.0) nxy(m)=i
      if(dabs(wb(m)).gt.sml .and.igx.eq.1) go to 18
      if(nxy(m).eq.0) nxy(m)=1
      if(nxy(m).eq.nmax) nxy(m)=nmax-1
      if(nxy(m).ne.i) ifl=1
  18  continue
      xb(m)=rby(i)
      return
  15  nxy(m)=i-1
      return
  20  continue
c *** print error message - set error flag for pilot         *****
      ier=1
      write(ioe,668)
      write(ioe,666) m,nmax,igx,lx1,lx2,nxy,xb,rby(lx2 ),sml,small
      call abend(ld)
      return
 668  format(//,5x,'fatal error in sort  *** particle lost on',
     1 ' entering or exiting array.')
 666  format(//,' m,nmax,igx,lx1,lx2,nxy',3x,8i5,//,' xb,rby(lx2),',
     1 'sml,small',/,5x,1p,6e18.8/)
c  the following statments not reachable since 4/92 abend kills job
      end
      subroutine stora(d,ld,lm,lp,nl,nlu,nxy,ip,nbu)
c * * this routine updates nesting table during tracking  * * *
      implicit real*8 (a-h,o-z)
      real*4 d
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      real*8 xb,wb,wp,xp,rin,rout,pinf,dist
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,ir,idbg ,
     1             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)
      dimension d(*),ld(*),nxy(3),lp(*),ip(*)
      if(nl.eq.0) go to 60
      if(nl.lt.0) go to 300
      if(nl.gt.nlev) go to 300
      lp(6*(nl-1)+2)=lm
      nasc=lp(6*(nl-1)+1)
      if(nlo.eq.1) lp(6*(nl-1)+1)=nbu
      do 10 i=1,3
      il=6*(nl-1)+i+2
  10  lp(il)=nxy(i)
      lp(6*nl)=nlu
  60  continue
      if(nlo.ge.0) go to 40
      do 20 i=1,6
      il=6*nl+i
  20  lp(il)=0
  40  continue
      lp(6*nlev+1)=nl
      lp(6*nlev+2)=nlo
      lp(6*nlev+4)=lp(6*nlev+3)
      lp(6*nlev+3)=nlu
      nbu=0
      return
 300  continue
      write(ioe,310) nlo,nl,nlu,lm
 310  format(//,10x,' in stora ',
     1      //,10x,' either an illegal level transfer was attempted ',
     2 //,10x,'or an incorrect array request was made ','* nlo=',
     3 i5,' nl=',i5,' nlu=',i5,' lm=',i5 )
      ix=2
      call dipr(d,ip,ix)
      call errtra
c coment out by st 96/12/18
      return
      end
*-----------------------------------------------------------------------

      subroutine unis(np,ip,d,ld,ma,fpd,locreg,numr,mcz,kbz,ibod,numbod)


*-----------------------------------------------------------------------
c * * this routine reads, verifies, and edits universe type input * *

      implicit real*8 (a-h,o-z)

      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      common/jomk/nbi,iri,jmk,iml,iop,imax
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)

      dimension ip(*),d(*),ld(*),ma(*),fpd(*),locreg(*),mcz(*),kbz(*),
     1 numbod(*)

      iflow=0
      mus=ip(np+14)
      muz=mus+ibod
      ip(np+14)=muz+ibod
      icr=1
      ier=0
      write(iot,3)
   3  format(1h1,//,20x,'universe specifications',//,15x,'universe',
     1 10x,'1-simple/0-comjom',//)
      do 200 i=1,ibod
      ii=mus+i-1
      ld(ii)=1
      if(jmk.eq.1) go to 5
      ld(ii)=iread(icr,iret)
   5  continue
      icr=0
      ld(muz+i-1)=0
      if(ld(ii).ne.1) go to 100
      do 95 j=1,numr
      if(kbz(j).ne.i) go to 95
      if(mcz(j).ge.0) go to 20
      if(mcz(j).eq.-1000) go to 20
      ld(ii)=0
      iab=iabs(mcz(j))
      write(iot,25) i,j,iab
  25  format(//,5x,'warning -- universe',i5,' is not simple',/,
     1 5x,'code zone',i5,' contains array',i5,//)
  22  format(/,15x,i5,17x,i5)
      go to 100
  20  continue
      if(numbod(j).le.2) go to 15
      write(iot,10) j,i
  10  format(/,10x,'code zone',i4,' in universe',i4,' cannot be',
     1 ' a simple universe',//)
      ld(ii)=0
      go to 100
  15  continue
      if(mcz(j).ne.-1000) go to 95
      n=locreg(j)+1
      num=numbod(j)*5+n-5
      do 50 ni=n,num
      nbo=ma(ni)
      if(nbo.ge.0) go to 50
      nbp=-nbo
      do 45 m=1,numr
      if(kbz(m).ne.i) go to 45
      if(m.eq.j) go to 45
      im=locreg(m)+1
      mum=numbod(m)*5+im-5
      do 40 mi=im,mum
      nbq=ma(mi)
      if(nbq.ne.nbp) go to 40
      ld(muz+i-1)=m
      go to 100
  40  continue
  45  continue
      ier=1
      write(iot,55) nbp,i
  55  format(/,10x,'body',i4,' in universe',i4,' is not referenced',
     1 ' positive in the universe',//)
      go to 95
  50  continue
  95  continue
 100  continue
      write(iot,22) i,ld(ii)
 200  continue
      if(ier.eq.0) go to 110
      call errtra
      return
 110  continue
      return
      end
      subroutine zexits(d,ma,locreg,numbod)
c * * this routine prints the combinatorial summary of zone transfers
c     at the end of a run.
      implicit real*8 (a-h,o-z)
      dimension ma(*),locreg(*),numbod(*)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
c
      lc=57
      do 800 i=1,numr
      n=locreg(i)+1
      num=numbod(i)*5+n-5
      do 800 j=n,num,5
      ict=j+2
      ic=1
      go to 300
 100  continue
      i1=ma(ict+2)
      ic=2
      go to 300
 120  continue
      ict=i1
      if(i1.gt.0) go to 100
      go to 800
 300  continue
      iv=1
      if(lc.ge.54) iv=2
      go to (330,360),iv
 330  continue
      go to (370,380),ic
 360  continue
      write(iot,20)
      lc=8
 370  continue
      lc=lc+1
      write(iot,25) i,ma(j)
      if(ic.eq.1) go to 100
 380  continue
      lc=lc+1
      write(iot,30) ict,ma(ict),ma(ict+1),i1
      go to 120
 800  continue
  20  format(1h1,//,34x,'combinatorial geometry zone transfer summary',
     1' tables',//,32x,'exit',5x,'exit',5x,'loc. in',5x,'entering',5x,
     2 'no. of',5x,'loc. for',/,32x,'zone',5x,'body',5x,'ma array',6x,
     3 'zone',6x,'entries',5x,'next zone',/)
  25  format(31x,i4,6x,i4)
  30  format(48x,i6,8x,i3,7x,i6,6x,i6)
      return
      end
