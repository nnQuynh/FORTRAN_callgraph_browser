*-----------------------------------------------------------------------

      subroutine cali(nl,nlu,xa,dc,ip,ncmax,ni4,
     1 nbod,ncn,nba,fpd,ma,locreg,d,ld)

c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     92/06/22
c     time the module was last permanently updated:     08:43:28
c     programmer name:                                  j.t.west
c     module name:                                      macali
c     current archiving level number:                   00001
c     current number of permanent updates:              00001
c     date of last access by librarian:                 93/03/15
c     dataset name:  x4s.scale4.master
c
c#######################################################################

c * * this routine determines particle location . it finds universe
c       , array and combinatorial input zone number particle is in.

      implicit real*8 (a-h,o-z)

      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)

      real*4 d

      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      real*8 xb,wb,wp,xp,rin,rout,pinf,dist
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,ir,idbg ,
     1             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)
      common/ark/ngy,nll,nbb,nzy(3),xd(3)
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      dimension ma(*),fpd(*),xa(3),ip(*),ncmax(*),dc(3),
     1 ni4(*),nbod(*),ncn(*),nba(*),d(*),ld(*),   locreg(*)

*-----------------------------------------------------------------------

      ips = 0

   5  continue

      idb=idbg

      do 10 i=1,3

         wb(i) = dc(i)
         xb(i)=7777777.

  10  continue

  12  continue

*-----------------------------------------------------------------------

            call lookz(xa(1),xa(2),xa(3),d,ma,fpd,d(klcr),d(knbd),
     &                 d(kior),d(knsr))

               if( ierror .ne. 0 ) return
               if( irprim .le. 0 ) return

*-----------------------------------------------------------------------

            nmed = ld(kmiz+irprim-1)

            if(idb.ne.0) then
               write(iot,711) nl,nlu,ir,irprim,nmed,nlev
               call dipr(d,ip,-1)
 711           format(5x,'cali nl,nlu,ir,irprim,nmed,nlev',10i5)
            end if

            if( nmed .ge. 0 ) return

*-----------------------------------------------------------------------
c * * return to calling program if not in array * *
*-----------------------------------------------------------------------

      if( nlev .le. 0 ) return

*-----------------------------------------------------------------------

      if( nmed .eq. -1000 ) go to 40

c * * nmed=-1000 means particle exiting array * * *

      nlu=iabs(nmed)
      ia=locreg(ir)
      nby=ma(ia+1)
      ii=ip(nq+11)
      nbu=nby

      call trente(xa,ma,fpd,d,ld,ip,nq,nba,nby)

  15  nl=nl+1
      nlo=+1
      m=ni4(nlu)
      igx=0
*-----------------------------------------------------------------------

      call cell(xa,ma,fpd,nlu,ncmax,ni4,d(m),nbod,ncn,d,ld,nx1,
     1 lu,lm,ier,lz)

*-----------------------------------------------------------------------

ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.

      if(ierror.ne.0) return
      ii=ip(nq+11)

      call stora(d,ld,lm,ld(ii),nl,nlu,nx1,ip,nbu)

ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.

      if(ierror.ne.0) return
      nll=nlu
      nbb=nby
      do 22 mm=1,3
      nzy(mm)=nx1(mm)
      xd(mm)=xa(mm)
  22  continue
      if(lz+ier.eq.0) go to 21
      ips=0
      write(ioe,23) lz,lu,ier
  23  format(/,5x,' error in cali from cell * lz=',i5,' lu=',i5,
     1 ' ier=', i5)
      go to 60
  21  continue
      l1=-1
      if(idb.eq.1) call dipr(d,ip,l1)
      n=3*(nlu-1)+1
      nr=ncn(n)
      nx=-1
*-----------------------------------------------------------------------

      call ctran(xa,ncmax(n),nx1,nx,ncn(n),d(nr),d,ld,ier)


*-----------------------------------------------------------------------
      if(lm) 25,40,30
  25  continue
      nlu=iabs(lm)
      go to 15
  30  continue
      nlo=0
      ll=lm
      nby=nbod(lm)
      call rtexit(xa,ma,fpd,d,ld,ip,nq,nba,nby)
      go to 12
  40  continue
      do 50 i=1,3
      xa(i)=xa(i)-(1.0e-04)*wb(i)
  50  continue
      ips=ips+1
      if(ips.le.1) go to 5
  60  continue
      write(ioe,600) nmed,ips,lz,ier,lu,ir,irprim,lm,ll,nl
 600  format(5x,'error in cali * nmed,ips,lz,ier,lu,ir,irprim,lm',
     1 ',ll,nl',10i5)
      idbg=0
*-----------------------------------------------------------------------

      call pr(d,1)

*-----------------------------------------------------------------------

      call abend(ld)

      return
      end

*-----------------------------------------------------------------------
c     call finefi(ld(i3),ld(i4),ld(i6),ld(i7),d(i1),ma,fpd,numbod,
c    1 locreg,kbz,mcz,d,ld,ip,ibod,numr)

      subroutine finefi(ncmax,ni4,nbod,ncn,wlh,ma,fpd,numbod,locreg,
     1 kbz,mcz,d,ld,ip,ibod,numr)

*-----------------------------------------------------------------------

c * * this routine checks array input with comb. geom. input
c     to be sure arrays fit snugly * * * * *
      implicit real*8 (a-h,o-z)
      real*4 d
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      dimension ncmax(*),ni4(*),nbod(*),ncn(*),wlh(*),ma(*),fpd(*),
     1 numbod(*),locreg(*),kbz(*),mcz(*),d(*),ld(*),ip(*)

      data eps /1.0d-9/

      ife=0
      do 500 n=1,nar
      icl=0
      n2=3*(n-1)+1
      ly=ncn(n2)
      if(ly.gt.0) go to 10
      write(iot,12) n
  12  format(//5x,'warning - in finefi',
     1       ' no reference was made to array no.',i5//)
      go to 500
  10  continue
      m=3*(n-1)
      mx=ncmax(m+1)
      my=ncmax(m+2)
      mz=ncmax(m+3)
      nr=ni4(n)
      do 400 k=1,mz
      do 400 j=1,my
      do 400 i=1,mx
      call cwidth(ddx,ddy,ddz,d(ly),mx,my,mz,i,j,k)
      ix=my*mx*(k-1)+mx*(j-1)+i
      lx=ix+nr-1
      lm=ld(lx)
      if(lm) 20,400,200
  20  continue
      n1=iabs(lm)
      if(n1.le.nar) go to 50
      write(iot,40) n,i,j,k,n1
  40  format(//10x,'in finefi invalid array reference in array ',i5
     1,' position ',3i5,' to nonexistent array ',i5)
      ife=1
      go to 400
  50  continue
      icl=1
      nu=3*(n1-1)
      dx=wlh(nu+1)
      dy=wlh(nu+2)
      dz=wlh(nu+3)
      go to 300
 200  continue
      if(lm.le.ibod) go to 210
      ife=1
      write(ioe,205) n,i,j,k,lm
 205  format( //10x,'finefi fatal error - array ',i5,' position',3i5,
     1 ' refers to an invalid or undefined universe',i5  //)
      go to 400
 210  continue
      icl=1
      nby=nbod(lm)
      loc=ma(7*nby)+2
      itype=ma(7*nby-4)
      if(itype.eq.7) go to 250
      if(itype.eq.9) go to 240
      write(ioe,220) n,i,j,k,lm,nby
 220  format( //10x,'in finefi array no.',i5,' position ',3i5,
     1 ' refers to universe ',i5/10x,' defined by body ',i5,
     2 ' which is not an rpp or box - fatal error - '  //)
      ife=1
 240  continue
      dx=fpd(loc+1)-fpd(loc)
      dy=fpd(loc+3)-fpd(loc+2)
      dz=fpd(loc+5)-fpd(loc+4)
      go to 300
 250  continue
      dx=dsqrt(fpd(loc+3)**2+fpd(loc+4)**2+fpd(loc+5)**2)
      dy=dsqrt(fpd(loc+6)**2+fpd(loc+7)**2+fpd(loc+8)**2)
      dz=dsqrt(fpd(loc+9)**2+fpd(loc+10)**2+fpd(loc+11)**2)
 300  continue
      iff=0
      if(dabs(dx-ddx).gt.eps) iff=1
      if(dabs(dy-ddy).gt.eps) iff=1
      if(dabs(dz-ddz).gt.eps) iff=1
      if(iff.eq.0) go to 400
      ife=1
      write(ioe,444) n,i,j,k,lm,dx,ddx,dy,ddy,dz,ddz
 444  format( //10x,'finefi - fatal error in array ',i5,
     1       ' element position ',3i5/20x,' contains universe ',i5,
     2       ' misfit in position '/30x,1p,6e15.5  //)
 400  continue
      if(icl.eq.1) go to 500
      write(ioe,540) n
 540  format( //10x,'in finefi ',
     1  ' array ',i5,' has no valid elements -very fatal')
      ife=1
 500  continue
      ivl=0
      do 700 l=1,numr
      if(mcz(l).ge.0) go to 700
      if(mcz(l).eq.-1000) go to 700
      if(kbz(l).eq.0) ivl=1
      mn=locreg(l)+1
      nbd=numbod(l)
      if(nbd.le.1) go to 630
      write(iot,620) l,nbd,kbz(l),mcz(l)
 620  format(//10x,'finefi - code zone ',i5,'contains an array but ',
     1       ' is defined with more than one body - warning',3i5/
     2       ' first body will be reference body'//)
 630  continue
      nby=ma(mn)
      n=iabs(mcz(l))
      if(n.le.nar) go to 635
      ife=1
      write(ioe,633) l,n
 633  format(//10x,'finefi fatal error - code zone',i5,
     1       ' references an invalid array',i5//)
      go to 700
 635  continue
      loc=ma(7*nby)+2
      itype=ma(7*nby-4)
      if(itype.eq.7) go to 670
      if(itype.eq.9) go to 650
      ife=1
      write(ioe,640) l,n,nby
 640  format(//10x,
     1  'finefi  fatal error -code zone ',i5,' contains array',i5,
     2 ' inside body',i5/10x,' referenced body is not an rpp or box')
      go to 700
 650  continue
      dx=fpd(loc+1)-fpd(loc)
      dy=fpd(loc+3)-fpd(loc+2)
      dz=fpd(loc+5)-fpd(loc+4)
      go to 690
 670  continue
      dx=dsqrt(fpd(loc+3)**2+fpd(loc+4)**2+fpd(loc+5)**2)
      dy=dsqrt(fpd(loc+6)**2+fpd(loc+7)**2+fpd(loc+8)**2)
      dz=dsqrt(fpd(loc+9)**2+fpd(loc+10)**2+fpd(loc+11)**2)
 690  continue
      nu=3*(n-1)
      ddx=wlh(nu+1)
      ddy=wlh(nu+2)
      ddz=wlh(nu+3)
      iff=0
      if(dabs(dx-ddx).gt.eps) iff=1
      if(dabs(dy-ddy).gt.eps) iff=1
      if(dabs(dz-ddz).gt.eps) iff=1
      if(iff.eq.0) go to 700
      ife=1
      write(ioe,695) n,nby,l,dx,ddx,dy,ddy,dz,ddz
 695  format( //10x,'finefi fatal error - array',i5,' in body',i5,
     1' referenced in code zone',i5,' does not fit'/5x,1p,6e20.10//)
      write(ioe,777) (fpd(loc+i-1),i=1,6),
     1 ((wlh(3*(ii-1)+ij),ij=1,3),ii=1,nar)
 777  format(/5x,1p,6e20.10,99(/10x,1p,3e20.10))
 700  continue
      if(ivl.eq.1) go to 720
      ife=1
      write(ioe,740)
 740  format(//10x,'finefi fatal error - you have at least one ',
     1 'array but no array was referenced in a level zero zone '//)
 720  continue
      if(ife.eq.0) go to 800
      call dipr(d,ip,0)
      call pr(d,1)
      call errtra
      return
 800  continue
      return
      end

************************************************************************
*                                                                      *
      subroutine jomin(d,ld,istr,naad,i1,i0,ioer,n16,n17,
     &                 ndsn,nresp,nmost,itype,isig,lim)
*                                                                      *
*        this routine reads all geometry data                          *
*        combinatorial + array                                         *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)
      dimension kaf(23),nob(10)
      dimension d(*),ld(*)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      common/parem/ixb(44)
!$OMP THREADPRIVATE(/parem/)
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      common/repeat/jp(20)


      ioe = ioer
      intt = i1
      iot=i0
      iout=n16
      iou2=n17
      if(isig.eq.2) go to 200

c **  isig=2 means restore geometry data * * *

      lptj=istr
      if(isig.eq.3) go to 100

c *** isig=3 only in junbug - means keno type geometry ***
c *** when isig=3, jomick was called previous to jomin call ***

      call jomin1(d(istr),ld(istr),istr,naad,lim)

ctat/mod  FORTRAN STOP statement kills process that started from C main.
ctat/mod  So to avoid sudden deth, this flag ierror added.

      if(ierror.ne.0) return

ctat/mod
 100  call jomin2(d,istr,naad)
      left=lim-istr-naad
      iw=istr+naad
      np=1
*-----------------------------------------------------------------------
c     subroutine azip(ip,np,d,ld,ifwa,left,ma,fpd,numb,locreg,numr,mcz,
c    1 numbod,kbz)


      call azip(jp,np,d, d,iw,left,d(kma),d(kfpd),numb,d(klcr),numr,
     1 d(kmcz),d(knbd),d(kbcz))

*-----------------------------------------------------------------------
ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.

      if(ierror.ne.0) return

      if(ndsn.le.0) go to 120
      i4=jp(4)

      call sazar(jp,np,d,d,d(i4),ndsn,nresp,nmost)

 120  continue
      idm(1) = np
      naad= jp(np)-jp(1)+naad+1
      nadd=naad
      if(itype.eq.2) return

c *** itype=2 is for a junbug plotting run  ***

      rewind iout

c   unit iout will contain binary for all geom data - cg + array

      if(itype.eq.1) write(iout)jp,
     &     nby,nlev,nar,nq,iaw,iay,nf,nx1,                   !FURUTA
     &     kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1, !FURUTA
     1     kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,   !FURUTA
     2     numb,nir,kbiz,kbcz,                               !FURUTA
     &     ixb                                               !FURUTA

c     write all of geometry data on unit iout in binary

      lptk=jp(np)
      write(iout) (d(i),i=lptj,lptk)
      rewind iout
      return

 200  continue

c     read geom data into correct location in d array

      call restor(d,ld,n16,istr,naad,ndsn)

      return
      end


*-----------------------------------------------------------------------

      subroutine pilot(nl,nlu,xa,ip,ncmax,ni4,
     1 nbod,ncn,nba,fpd,ma,locreg,d,ld,distd)

c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     95/08/30
c     time the module was last permanently updated:     13:51:09
c     programmer name:                                  j.t.west
c     module name:                                      mapilot
c     current archiving level number:                   00003
c     current number of permanent updates:              00003
c     date of last access by librarian:                 95/08/30
c     dataset name:  x4s.scale4.master
c
c#######################################################################

c     this routine determines next input zone particle enters ***
c     on entering or leaving pilot, ir=code zone
c                                  irprim=input zone
c     inside  pilot when calling g1, irprim=code zone
c                                   ir=input zone
c ****** this makes pilot consistent with g1 ******

      implicit real*8 (a-h,o-z)

      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)

      real*4 d

      common/ark/ngy,nll,nbb,nzy(3),xd(3)
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      real*8 xb,wb,wp,xp,rin,rout,pinf,dist
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,ir,idbg ,
     1             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)
      real*8 dist0
      integer blzold
      common /orgi/ dist0,markg,nmedg,nblz,blzold,irpold
!$OMP THREADPRIVATE(/orgi/)
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common/ss/s1
!$OMP THREADPRIVATE(/ss/)

      dimension ma(*),fpd(*),xa(3),ip(*),ncmax(*),
     1 ni4(*),nbod(*),ncn(*),nba(*),d(*),ld(*),locreg(*),nx2(3)

      save lm
!$OMP THREADPRIVATE(lm)
*-----------------------------------------------------------------------

               idb = idbg
               ict = 0

               nmold = ld(kmiz+irprim-1)
               nrold = ld(kriz+irprim-1)

               ips   = 0
               iexit = 0
               s1    = 0.0
               in    = ip(nq+11)

               xb(1) = xa(1)
               xb(2) = xa(2)
               xb(3) = xa(3)

*-----------------------------------------------------------------------

  12  continue

*-----------------------------------------------------------------------

         if( irprim .eq. 0 .and. ir .eq. 0 ) then
            ierror = 1
            return
         end if

*-----------------------------------------------------------------------

               markg = 0

               itemp  = ir
               ir     = irprim
               irprim = itemp

*-----------------------------------------------------------------------

      call g1(s,ma,fpd,d(klcr),d(knbd),d(kior),ir1,ir2,d,distd)

*-----------------------------------------------------------------------

            if( ierror .ne. 0 .or. irprim .le. 0 ) return

*-----------------------------------------------------------------------

         s1 = s1 + s

 700  continue

         ict = ict + 1

*-----------------------------------------------------------------------

         if( ict .gt. 500 ) then
            write(ioe,640)
 640        format(' problem encountered in array processing,',
     &             ' particle will be terminated as escape.' )
            irprim = -3
            return
         end if

*-----------------------------------------------------------------------

         nmed = ld(kmcz+ir-1)

*-----------------------------------------------------------------------

         if( idb .ne. 0 ) then
            write(ioe,711) nl,nlu,ir,irprim,nmed,markg,nlev,dist0,nmold
 711        format(5x,'pilot crash',7i5,1pe15.5,i5)
            write(ioe,777) nasc,nbo,lsurf,lri,lro,kloop,loop,
     &                     itype,s,dist,rin,rout
 777        format(/,5x,8i5,/,5x,1p4e15.5,/)

            call dipr(d,ip,-1)

         end if

*-----------------------------------------------------------------------
*  *** if nmed<0, go to array processing routines
*-----------------------------------------------------------------------

      if( nmed .lt .0 ) goto 14

*-----------------------------------------------------------------------
*  **  if completed path; i.e.,  real collision - return
*-----------------------------------------------------------------------

      if( markg .eq. 1 ) then

               return

*-----------------------------------------------------------------------
*  *** if markg = 0 and nmed or nreg has changed,it's bdryx - return
*-----------------------------------------------------------------------

      else

               nreg = ld(krcz+ir-1)

*-----------------------------------------------------------------------
* ***  markg = 0, no change in nmed or nreg - continue in pilot
*-----------------------------------------------------------------------

            if( nmed .eq. nmold .and.
     &          nreg .eq. nrold .and.
     &          ipret .le. 0  ) then

               blzold = irprim + 32768 * ir
               irpold = irprim

               goto 12

            end if

               return

      end if

*-----------------------------------------------------------------------
* *** array geom calculations follow  **********
*-----------------------------------------------------------------------

  14  continue

      kloop=kloop+1
      do 17 i=1,3
      xa(i)=xa(i)+dist*wb(i)
  17  continue

      dist0 = dist0 - dist

      dist = 0.0

      if(nmed.eq.-1000) go to 500

c ******* particle entering array ************;

      nlu=iabs(nmed)
      ia=locreg(ir)
      nby=ma(ia+1)
      nbu=nby
      call trente(xa,ma,fpd,d,ld,ip,nq,nba,nby)
  20  nl=nl+1
      nlo=+1
      m=ni4(nlu)
      igx=0

      call cell(xa,ma,fpd,nlu,ncmax,ni4,d(m),nbod,ncn,d,ld,nx1,
     1 lu,lm,ier,lz)

      if(ierror.ne.0) return

      call stora(d,ld,lm,ld(in),nl,nlu,nx1,ip,nbu)

      if(ierror.ne.0) return

      nll=nlu
      nbb=nby
      do 22 mm=1,3
      nzy(mm)=nx1(mm)
      xd(mm)=xa(mm)
  22  continue
      l1=-1
      if(idb.eq.1) call dipr(d,ip,l1)
      if(ier.ne.0) go to 40

c **** lz gt 0  fatal on entering array ******

      if(lz.gt.0.or.lm.eq.0) go to 40
 300  continue
      n=3*(nlu-1)+1
      nr=ncn(n)
      nx=-1
      call ctran(xa,ncmax(n),nx1,nx,ncn(n),d(nr),d,ld,ier)
      if(lm) 25,40,30
  25  continue
      nlu=iabs(lm)
      go to 20
  30  continue
      nlo=0
      ll=lm
      nby=nbod(lm)
      call rtexit(xa,ma,fpd,d,ld,ip,nq,nba,nby)
      go to 580
 500  continue

c ***** particle exiting an array **********
c *****       nmed=-1000          **********

      do 510 i=1,3
 510  nx2(i)=nx1(i)
      nby=nbod(ll)
      n=3*(nlu-1)+1
      m=ni4(nlu)
      nr=ncn(n)
      call trente(xa,ma,fpd,d,ld,ip,nq,nba,nby)
      nx=1
      call ctran(xa,ncmax(n),nx2,nx,ncn(n),d(nr),d,ld,ier)
 515  continue
      igx=1
      call cell(xa,ma,fpd,nlu,ncmax,ni4,d(m),nbod,ncn,d,ld,nx1,
     1 lu,lm,ier,lz)
ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.
      if(ierror.ne.0) return
      if(ier.ne.0) go to 50
      if(lz.gt.0) go to 550
      if(lm.eq.0) go to 550
      iexit=0
 520  continue
      nx=-1
      call ctran(xa,ncmax(n),nx1,nx,ncn(n),d(nr),d,ld,ier)
      nlo=0
      if(lm.gt.0) go to 530
      call stora(d,ld,lm,ld(in),nl,nlu,nx1,ip,nbu)
ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.
      if(ierror.ne.0) return
      nlu=iabs(lm)
      go to 20
 530  continue
      nby=nbod(lm)
      ll=lm
      call rtexit(xa,ma,fpd,d,ld,ip,nq,nba,nby)
      go to 580
 550  nl1=nl
      if(iexit.eq.1) go to 50
      if(lm.eq.0) iexit=1
      nl=nl-1
      nlo=-1
      call clev(nl,nl,ld(in),ncmax,ni4,nbod,ncn,nba,xa,xa,
     1 d,ld,nx1,d(kma),d(kfpd),ip)
ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.
      if(ierror.ne.0) return
      nasc=ld(in+6*nl)
      nlu=nll
      lm=ngy
      n=3*(nlu-1)+1
      nr=ncn(n)
      m=ni4(nlu)
      call stora(d,ld,lm,ld(in),nl,nlu,nx1,ip,nbu)
ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.
      if(ierror.ne.0) return
      if(nl) 50,590,560
 560  continue
      if(lm) 515,50,520
 580  continue
      call stora(d,ld,lm,ld(in),nl,nlu,nx1,ip,nbu)
ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.
      if(ierror.ne.0) return
      nll=nlu
      nbb=nby
      do 24 mm=1,3
      nzy(mm)=nx1(mm)
      xd(mm)=xa(mm)
  24  continue
 590  continue
      l1=-1
      if(idb.eq.1) call dipr(d,ip,l1)

      markg=0

      if(ld(mus+ll-1).eq.1.and.ll.ge.1) go to 200
      do 595 i=1,3
 595  xb(i)=0.0
      if(lz.eq.0) go to 597
      xa(lz) =xa(lz)+ dsign(1.0d0,wb(lz))*1.0d-9

c ***** above stmt moves particle off bndry in correct direction ****

597   call lookz(xa(1),xa(2),xa(3),d,ma,fpd,d(klcr),d(knbd),

     1 d(kior),d(knsr))
      go to 210
 200  continue
      do 610 i=1,3
 610  xb(i) = xa(i)
      kloop=kloop+1
      ir=ld(muz+ll-1)
      irprim=ld(kior+ir-1)
      iflow=0
 210  continue
      if(irprim.le.0) return
      go to 700
  40  continue
      write(ioe,600) lz,ier,lm,iexit,nl,nmed,ir,irprim
 600  format(5x,'error in pilot on entering array',10i5)
      go to 55
  50  continue
      write(ioe,620) lz,ier,lm,iexit,nl,nmed,ir,irprim
 620  format(5x,'error in pilot on exiting array',10i5)
  55  idbg=0
      call pr(d,1)
      call abend(ld)
      return
      end


