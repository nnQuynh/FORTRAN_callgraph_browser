*-----------------------------------------------------------------------
c     call azip(jp,np,d, d,iw,left,d(kma),d(kfpd),numb,d(klcr),numr,
c    1 d(kmcz),d(knbd),d(kbcz))

      subroutine azip(ip,np,d,ld,ifwa,left,ma,fpd,numb,locreg,numr,mcz,
     1 numbod,kbz)

*-----------------------------------------------------------------------
c ** azip reads array geometry data **

      implicit real*8 (a-h,o-z)

ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.

      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)
*-----------------------------------------------------------------------
      real*4 d

*-----------------------------------------------------------------------
      parameter ( ispcd = 2 )
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      common/jomk/nbi,iri,jmk,iml,iop,imax
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)

      dimension d(*),ld(*),ip(*),ma(*),fpd(*),locreg(*),mcz(*),
     1 numbod(*),kbz(*)

      icr=1
      iret=0
      nar=0
      nlev=0
      ip(np)=ifwa
      nq=np
      id=1
   7  ix=ip(np)+3*(id-1)

      ld(ix)=iread(icr,iret)
         if(ld(ix).eq.0) go to 14
      ld(ix+1)=iread(0,iret)
      ld(ix+2)=iread(0,iret)

      icr=0
      id=id+1
      go to 7

  14  nar=id-1
      if(nar.le.0) go to 105
      if(jmk.le.1) go to 15
      if(iml.eq.1) go to 15
      nar=2
      id=id+1
      ld(ix)=1
      ld(ix+1)=1
      ld(ix+2)=1

  15  continue
      ip(np+1)=ip(np)+nar
      ip(np+1)=2*(ip(np+1)/2)+1
      ip(np+2)=ip(np+1)+ispcd*3*nar
      ip(np+3)=ip(np+2)
      ip(np+4)=ip(np+3)+3*nar
      n3=3*nar
      do 16 i=1,n3
      i1=ip(np)+i-1
      i2=ip(np+3)+i-1
      ld(i2)=ld(i1)
      ld(i1)=0
  16  continue
      nmax=0
      do 10 i=1,numr
      if(kbz(i).gt.nmax) nmax=kbz(i)
      if(mcz(i).ge.0) go to 10
      if(mcz(i).eq.-1000) go to 10
      ic=iabs(mcz(i))
      iz=ip(np)+ic-1
      ld(iz)=1
  10  continue
c.....array pointers....
c......ip(np+3) is the array size , nxmax by nymax by nzmax for each arr
c.....ip(np+4) is the lookup pointer for array cell definitions
c.... ip(np+5) is the cell name definition array for each array
c
      iaw=0
      iay=1
      do 30 n=1,nar
      ind=ip(np+3)+3*(n-1)
      iaz=ld(ind)*ld(ind+1)*ld(ind+2)
      iax=ld(ind)+ld(ind+1)+ld(ind+2)+3
      idn=ip(np+4)+n-1
      ld(idn)=iaz
      iay=iay+iaz
      iaw=iaw+iax
  30  continue
      ip(np+5)=ip(np+4)+nar+1
      ip(np+6)=ip(np+5)+iay-1
      ibod=1
      icr=1
      iret=0
      mm=0

      iop=iread(0,iret)

      do 40 n=1,nar
      idn=ip(np+4)+n-1
      ia=ld(idn)
      if(jmk.le.1) go to 125
      if(n.eq.1) go to 125
      ind=ip(np+5)+mm
      ld(ind)=imax
      if(imax.gt.ibod) ibod=imax
      go to 150

 125  continue

      if(iop.eq.1) go to 130
      if(iop.eq.2) go to 150

c iop = 0

      do 35 i=1,ia
      ind=ip(np+5)+mm+i-1
      ld(ind)=iread(icr,iret)
      icr=0
      if(ld(ind).gt.ibod) ibod=ld(ind)
  35  continue

      go to 150

c iop = 1
 130  continue
      ic=1
      iz=ip(np+3)+3*(n-1)
      nzmax=ld(iz+2)
      nymax=ld(iz+1)
      nxmax=ld(iz)

 135  continue
      lcx=iread(ic,ir)
      if(lcx.gt.ibod) ibod=lcx
      ic=0
      i1=iread(0,ir)
      i2=iread(0,ir)
      is=iread(0,ir)
      j1=iread(0,ir)
      j2=iread(0,ir)
      js=iread(0,ir)
      k1=iread(0,ir)
      k2=iread(0,ir)
      ks=iread(0,ir)
      if(is.le.0 .or. js.le.0 .or. ks.le.0) go to 136
      if(i1.lt.1.or.i1.gt.nxmax .or. i2.lt.1.or. i2.gt.nxmax)go to 136
      if(j1.lt.1.or.j1.gt.nymax .or. j2.lt.1.or. j2.gt.nymax)go to 136
      if(k1.lt.1.or.k1.gt.nzmax .or. k2.lt.1.or. k2.gt.nzmax)go to 136
      go to 138

 136  write(ioe, 120)i1,i2,is,j1,j2,js,k1,k2,ks
 120  format('0there is an error in the array data for loop indices'
     * ,' when iop=1', /
     1  '0**for i',3i5,5x,'for j',3i5,5x,'for k',3i5)
      call errtra
      return

 138  do 140 k=k1,k2,ks
      do 140 j=j1,j2,js
      do 140 i=i1,i2,is
      ind=nymax*nxmax*(k-1)+nxmax*(j-1)+i+ip(np+5)-1+mm
      ld(ind)=lcx
 140  continue
      iflag=iread(0,ir)
      if(iflag.eq.0) go to 135


c iop = 2
 150  continue
      ld(idn)=ip(np+5)+mm
      mm=mm+ia
  40  continue
      ld(idn+1)=ip(np+5)+mm
      if(iop.ne.2) go to 55
      ih=ip(np+4)
      ih1=-99
      ih2=1
      ibod=nmax
      call fidas(d,d,ld(ih),ih1,ih2,intt,iot)
  55  continue
      do 200 n=1,nar
      iz=ip(np+3)+3*(n-1)
      nzmax=ld(iz+2)
      nymax=ld(iz+1)
      nxmax=ld(iz)
      write(iot,185) n,nxmax,nymax,nzmax
      if(ibod.le.1) go to 200
 185  format(1h1,15x,'array no.',i5,5x,'array size is',i5,' by',
     1 i5,' by',i5   )
 187  format(6x,'x =',1x,30i3)
 188  format(4x,'y'  )
      idn=ip(np+4)+n-1
      mm=ld(idn)
      do 172 k=1,nzmax
      m1=1
      m2=30
 160  continue
      if(m2.gt.nxmax) m2=nxmax
      write(iot,175) k,n
      write(iot,187) (i,i=m1,m2)
      write(iot,188)
      do 170 jx=1,nymax
      j=nymax-jx+1
      ind=nymax*nxmax*(k-1)+nxmax*(j-1)-1+mm
      write(iot,180) j,(ld(ind+i),i=m1,m2)
 170  continue
      if(m2.ge.nxmax) go to 172
      m1=m1+30
      m2=m2+30
      go to 160
 172  continue
 200  continue
 175  format( //,15x,' level ',i5,' of array no. ',i5,/  )
 180  format(/,i5,5x,30i3,3(/,10x,30i3))
      ip(np+7)=ip(np+6)+ibod
      do 48 i=1,ibod
      do 46 l=1,numr
      if(kbz(l).ne.i) go to 46
      if(mcz(l).ne.-1000) go to 46
      n=locreg(l)+1
      num=numbod(l)*5+n-5
      do 45 j=n,num,5
      nby=ma(j)
      if(nby.gt.0) go to 45
      ix=ip(np+6)+i-1
      ld(ix)=iabs(nby)
      go to 48
  45  continue
  46  continue
      write(ioe,52) i
  52  format(//,10x,'fatal error - universe no.',i5,' is not defined')
      call errtra
      return
  48  continue
      ip(np+8)=ip(np+7)+3*nar
      ip(np+8)=2*(ip(np+8)/2)+1
      i0=ip(np)
      i3=ip(np+3)
      i4=ip(np+4)
      i6=ip(np+6)
*-----------------------------------------------------------------------
c     subroutine level(nlv,ncmax,ni4,nbod,d,ld,kbz,mcz,numr)

      call level(d(i0),d(i3),d(i4),d(i6),d,d,kbz,mcz,numr)

*-----------------------------------------------------------------------

      ind=ip(np+8)
      i1=ip(np+1)
      do 64 lv=1,nlev
      l=nlev-lv+1
      do 60 n=1,nar
      ix=ip(np)+n-1
      idx=ip(np+7)+3*(n-1)
      if(ld(ix).le.0) ld(idx)=0
      if(ld(ix).ne.l) go to 60
      i3=ip(np+3)+3*(n-1)
      i4=ip(np+4)+n-1
      i5=ld(i4)
      i6=ip(np+6)
      iax=ld(i3)+ld(i3+1)+ld(i3+2)+3
      ld(idx)=ind
      ld(idx+1)=ld(idx)+ispcd*(ld(i3)+1)
      ld(idx+2)=ld(idx+1)+ispcd*(ld(i3+1)+1)
*-----------------------------------------------------------------------

csubroutine argen(d,ld,ncmax, nbn,     nbod,    rby,ma,fpd,wlh,  na)
       call argen(d, d,ld(i3),ld(i5),ld(i6),ld(ind),ma,fpd,d(i1),n)

*-----------------------------------------------------------------------
ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.
      if(ierror.ne.0) return
      ind=ind+ispcd*iax
  60  continue
  64  continue
      ip(np+9)=ind
      ip(np+10)=ip(np+9)+numb
      ip(np+10)=ispcd*(ip(np+10)/ispcd)+1
      nf=ip(np+10)
      do 95 nby=1,numb
      l=ip(np+9)+nby-1
      ld(l)=0
      km=7*(nby-1)+1
      if(ma(km+2).ne.7) go to 95
      do 80 ib=1,numr
      if(mcz(ib).ge.0) go to 80
      if(mcz(ib).eq.-1000) go to 80
      ia=locreg(ib)
      if(nby.eq.ma(ia+1)) go to 87
  80  continue
      do 85 ib=1,ibod
      ii=ip(np+6)+ib-1
      if(nby.eq.ld(ii)) go to 87
  85  continue
      go to 95
  87  continue
      ld(l)=nf
*-----------------------------------------------------------------------
c     subroutine orthom(d,ld,ma,fpd,nasc1,tr)

      call orthom(d,d,ma,fpd,nby,d(nf))

*-----------------------------------------------------------------------

      nf=nf+ispcd*9
  95  continue
      ip(np+11)=nf
      ip(np+12)=ip(np+11)+6*nlev+4
      lt=ip(np+11)
      lu=lt+6*nlev+4
      do 100 i=lt,lu
      ld(i)=0
 100  continue
      i3=ip(np+3)
      i4=ip(np+4)
      i6=ip(np+6)
      i1=ip(np+1)
      i7=ip(np+7)
*-----------------------------------------------------------------------
c     subroutine finefi(ncmax,ni4,nbod,ncn,wlh,ma,fpd,numbod,locreg,
c    1 kbz,mcz,d,ld,ip,ibod,numr)

      call finefi(ld(i3),ld(i4),ld(i6),ld(i7),d(i1),ma,fpd,numbod,
     1 locreg,kbz,mcz,d,ld,ip,ibod,numr)


*-----------------------------------------------------------------------
ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.
      if(ierror.ne.0) return
      go to 108
 105  continue
      do 106 i=1,11
      ip(np+i)=ip(np+i-1)
 106  continue
      ip(np+12)=ip(np+11)+6*nlev+4
      ibod=0
 108  continue
      ip(np+13)=ip(np+12)+nar+1
      ip(np+14)=ip(np+13)+ispcd*ibod+2
      mus=ip(np+14)
      muz=mus
      iflow=0
      if(ibod.le.0) go to 109

*-----------------------------------------------------------------------
c     subroutine unis(np,ip,d,ld,ma,fpd,locreg,numr,mcz,kbz,ibod,numbod)

      call unis(np,ip,d,d,ma,fpd,locreg,numr,mcz,kbz,ibod,numbod)

*-----------------------------------------------------------------------

ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.
      if(ierror.ne.0) return
 109  continue
      ip(np+15)=ip(np+14)
      if(nar.ge.1) ip(np+15)=ip(np+14)+6*nlev+4
      np=np+15
      iit=ip(np)-ip(nq)+1
      left=left-iit
      write(iot,110) iit,left
 110  format(//,' array data requires',i10,' locations, leaving',
     1 i10,' locations',//)
      l1=0
      iaw=ibod
      if(nar.le.0) return
      return
      end
      subroutine bod(nxy,nbn,ll,ncmax,lz)
c * * this routine returns lattice cell contents * * *
c * * sets flag if cell is vacant * * *
      dimension nxy(3),nbn(*),ncmax(3)
      do 10 i=1,3
      if(nxy(i).le.0.or.nxy(i).gt.ncmax(i)) go to 20
  10  continue
      ind=ncmax(1)*ncmax(2)*(nxy(3)-1)+ncmax(1)*(nxy(2)-1)+nxy(1)
      ll=nbn(ind)
      return
  20  lz=i
      return
      end

*-----------------------------------------------------------------------

      subroutine cell(xm,ma,fpd,na,ncmax,ni4,nbn,nbod,ncn,d,ld,
     1 nxy,lu,lm,ier,lz)
c *  this routine calculates lattice cell position and content * *
c *  for current particle position when in or moving thru arrays *
      implicit real*8 (a-h,o-z)
ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.
      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)
      real*4 d
      real*8 xb,wb,wp,xp,rin,rout,pinf,dist
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,ir,idbg ,
     1             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      dimension xm(3),ma(*),fpd(*),ncmax(*),ni4(*),nbn(*),
     1 nbod(*),ncn(*),d(*),ld(*),nxy(3)
      dimension ix(3),ibox(6)
      dimension mj(3)
      data ibox/3,4,1,2,5,6/
      data ix/2,1,3/
      n=3*(na-1)+1
      lu=0
      ier=0
      mj(1)=1
      mj(2)=2
      mj(3)=3
      if(igx.ne.1) go to 5
      ib=ibox(lro)
      m=(ib+1)/2
      mj(1)=m
      mj(m)=1
   5  continue
      do 10 k=1,3
      j=mj(k)
      nr=ncn(n+j-1)

      call sort(xm,ncmax(n),lu,d(nr),nxy,d,ld,wb,j,ier,ifl)

ctat  fortran stop statement kills process that started from c main.
ctat  so to avoid sudden deth, this flag ierror added.
      if(ierror.ne.0) return
      if(igx*lu.gt.0) go to 14
  10  continue
  14  continue
      lz=0
      ips=0
  15  continue
      call bod(nxy,nbn,lm,ncmax(n),lz)
      if(lz.eq.0) go to 50
      lsurf=ix(lz)*2
      if(wb(lz).le.0.0) lsurf=ix(lz)*2-1
      lsurf=-lsurf
  50  continue
      if(lm.ne.0.and.lz.eq.0) return
      do 40 i=1,3
  40  if(dabs(wb(i)).le.1.0d-08) go to 25
      return
  25  ips=ips+1
      if(ips.gt.1) return
      lu=0
      nr=ncn(n+i-1)

      call sort(xm,ncmax(n),lu,d(nr),nxy,d,ld,wb,i,ier,ifl)

      if(lu.ne.i) go to 35
      if(ifl.eq.1) go to 35
  30  nxy(i)=nxy(i)-1
      go to 15
  35  continue
      return
      end
      subroutine clev(nl1,nl2,lp,ncmax,ni4,nbod,ncn,nba,xa,
     1 xb,d,ld,nxy,ma,fpd,ip)
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     92/06/22
c     time the module was last permanently updated:     08:47:48
c     programmer name:                                  j.t.west
c     module name:                                      maclev
c     current archiving level number:                   00001
c     current number of permanent updates:              00001
c     date of last access by librarian:                 93/03/15
c     dataset name:  x4s.scale4.master
c
c#######################################################################
c * this routine translates particle coordinate from higher level,*
c *  more local coord system to lower,more global geom level no.
      implicit real*8 (a-h,o-z)

      real*4 d

      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common/ark/ngy,nll,nbb,nzy(3),xd(3)
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      dimension lp(*),xa(3),xc(3),ip(*),ncmax(*),
     1 ni4(*),nbod(*),ncn(*),nba(*),nxy(3),d(*),ld(*),ma(*),fpd(*),xb(3)
c
c......subroutine clev translates coordinates in xa in level nl1
c....to coordinages xb in level nl2, where level nl2 is less than nl1, a
c.....level nl2 is a more global level (higher in the hieharchy) than le
c
c.....coordinate xa are in level nl1
c
      do 5 i=1,3
   5  xb(i)=xa(i)
      nlm=nl1-nl2+1
      if(nl1.lt.nl2) go to 777
      do 40 l=1,nlm
      nlt=nl1-l+1
      if(nlt.gt.0) go to 30
      nbb=lp(1)
      nll=lp(6)
      ngy=nbb
      ll=0
      call rtexit(xb,ma,fpd,d,ld,ip,nq,nba,nbb)
      return
  30  continue
      il=6*nlt
      ii=6*(nlt-1)
      if(nlt.eq.nlev) go to 15
      nby=lp(il+1)
      if(nby.eq.0) go to 15
      call rtexit(xb,ma,fpd,d,ld,ip,nq,nba,nby)
  15  continue
      ngy=lp(ii+2)
      do 10 i=1,3
  10  nxy(i)=lp(ii+i+2)
      nll=lp(il)
      if(ngy) 25,777,20
  20  continue
      ll=ngy
      nbb=nbod(ngy)
      call trente(xb,ma,fpd,d,ld,ip,nq,nba,nbb)
  25  continue
      n=3*(nll-1)+1
      nr=ncn(n)
      nx=+1

      call ctran(xb,ncmax(n),nxy,nx,ncn(n),d(nr),d,ld,ier)

  40  continue
      return
 777  continue
      write(ioe,600) nl1,nl2,ngy
 600  format(5x,'error in clev * nl1=',i5,' nl2=',i5,' ngy=',i5)
      call dipr(d,ip,2)
      call abend(ld)
      return
      end

      subroutine dupr(d,ld,ip,ncmax,ni4,nbod,ncn,nba,lp,ixd)
c * * this routine dumps mars array tables including nesting table
c * * it should be called only for debugging purposes * * * *
      implicit real*8 (a-h,o-z)
      real*4 d
      common/ark/ngy,nll,nbb,nzy(3),xd(3)
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      real*8 xb,wb,wp,xp,rin,rout,pinf,dist
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,ir,idbg ,
     1             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      dimension d(*),ld(*),ip(*),ncmax(*),ni4(*),
     1 nbod(*),ncn(*),nba(*),lp(*)
      parameter ( ispcd = 2 )
      if(ixd.lt.0) go to 80
      do 5 i=1,nar
      i1=ip(nq)+i-1
      i2=ip(nq+1)+ispcd*3*(i-1)
      write(iot,12) i,ld(i1)
      call dxpr(d(i2),1,i,3,iot)
   5  continue
      write(iot,10) nq,(ip(nq+i-1),i=1,13)
      write(iot,15) nby,nlev,nar,nq,iaw,iay,nf,nx1,ixd
      write(iot,17)
      do 20 i=1,nar
      n=3*(i-1)
      write(iot,25) i,(ncmax(n+j),j=1,3)
  20  continue
      write(iot,125)
      do 30 m=1,nar
      n=3*(m-1)
      write(iot,35) m,ni4(m),(ncn(n+i),i=1,3)
      ln=ni4(m)
      nx=ncmax(n+1)
      ny=ncmax(n+2)
      nz=ncmax(n+3)
      mm=0
      write(iot,105)
      do 47 i=1,3
      ntt=ncmax(n+i)+1
      lq=ncn(n+i)
      call dxpr(d(lq),ntt,i,1,iot)
  47  continue
  30  continue
      np=nq
      ibod=iaw
      do 200 n=1,nar
      iz=ip(np+3)+3*(n-1)
      nzmax=ld(iz+2)
      nymax=ld(iz+1)
      nxmax=ld(iz)
      write(iot,185) n,nxmax,nymax,nzmax
      if(ibod.le.1) go to 200
      idn=ip(np+4)+n-1
      mm=ld(idn)
      do 172 k=1,nzmax
      m1=1
      m2=30
 160  continue
      if(m2.gt.nxmax) m2=nxmax
      write(iot,175) k,n
      write(iot,187) (i,i=m1,m2)
      write(iot,188)
      do 170 jx=1,nymax
      j=nymax-jx+1
      ind=nymax*nxmax*(k-1)+nxmax*(j-1)-1+mm
      write(iot,180) j,(ld(ind+i),i=m1,m2)
 170  continue
      if(m2.ge.nxmax) go to 172
      m1=m1+30
      m2=m2+30
      go to 160
 172  continue
 200  continue
      write(iot,38)
      do 40 l=1,iaw
  40  write(iot,42) l,nbod(l)
      write(iot,45)
      do 60 i=1,numb
      write(iot,65) i,nba(i)
      if(nba(i).le.0) go to 60
      nrt=nba(i)
      call dxpr(d(nrt),3,i,3,iot)
  60  continue
      if(ixd.eq.0) return
  80  continue
      write(iot,85) nll,nbb,ll
      do 90 i=1,3
      write(iot,95) i,nzy(i),xd(i),xb(i),wb(i)
  90  continue
      write(iot,96)
      do 110 i=1,nlev
      ni=6*(i-1)
      write(iot,115) i,(lp(ni+j),j=1,6)
 110  continue
      write(iot,115) (lp(6*nlev+j),j=1,4)
      return
  12  format(//6x,'level',2x,'nlv' /,5x,2i5,10x,'array dimensions' )
  10  format(//,10x,'array pointers',//,1x,14i5,/)
  15  format(//,10x,' arar ',10i5, 5x,' ixd',i5)
  25  format(  i10,20x,3i10)
  17  format(//,7x,'level',5x,'array',5x,'array',8x,
     1 'array  dimensions',11x,'box',/,15x,'universe',2x,'reference',
     2 3x,'nxmax',5x,'nymax',5x,'nzmax   r-matrix',/,26x,'body',
     3 36x,'pointer',//)
 185  format(1h0    ,15x,'array no.',i5,5x,'array size is',i5,' by',
     1 i5,' by',i5   )
 187  format(6x,'x =',1x,30i3,3(//,10x,30i3))
 188  format(4x,'y',/)
 175  format( //,15x,' level ',i5,' of array no. ',i5,// )
 180  format(/,i5,5x,30i3,3(/,10x,30i3))
  42  format(  5x,2i5  )
  35  format(/,5x,5i10,/)
  65  format(  5x,2i10   )
  85  format(// ' debug',5x,'nll,nbb,ll',3i5 /4x,'i',3x,'nzy',8x,
     1 'xd(i)',10x,'xb(i)',10x,'wb(i)'  )
  95  format(2i5,3f15.5)
 115  format( 7i10)
  96  format(1h0,'  particle nesting array, lp '  )
 105  format(1h0,8x,'i',3x,'j',4x,'lattice cell bndry' )
 125  format(10x,'level',4x,'ary ptr',5x,'addresses - cell boundaries')
  38  format(5x,'univ',2x,'body' )
  45  format(8x,'body',7x,'nba' )
      end


************************************************************************
*                                                                      *
        subroutine jomin2(nstor,istr,naad)

c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     93/03/15
c     time the module was last permanently updated:     14:36:56
c     programmer name:                                  j.t.west
c     module name:                                      majomin2
c     current archiving level number:                   00003
c     current number of permanent updates:              00003
c     date of last access by librarian:                 93/03/15
c     dataset name:  x4s.scale4.master
c
c#######################################################################
c     jomin2 reads the binary geom file from unit iout and calls geni
      implicit real*8 (a-h,o-z)
      character*4 jty1(15)
      common/gomloc/ kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     1 kkr2,knsr,kvol,nadd,ldata,ltma,lfpd,numr,irtru,numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      real*8 xb,wb,wp,xp,rin,rout,pinf,dist
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,ir,idbg ,
     1             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)
      common/tape/intt,iot,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      dimension nstor(*)
      parameter ( ispcd = 2 )
      rewind iout
      read(iout)(jty1(i),i=1,15)
      read(iout)ivopt,idbg,irtru,numr,numb,ladd,nazt,ibodt,ltma,lfpd
      write(iot,25)(jty1(i),i=1,15),ivopt,idbg
 25   format (1h1//10x,15a4//20x,7hivopt =,i2,10x,6hidbg =,i2/)
      kma=istr
   10 kfpd=kma+ltma
c.....fpd array must start on double precision word boundary for ibm
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
      krf=krcz-1
      kmf=kmcz-1
      klf=kbcz-1
c     to read region data
      ndx=kriz-1
      read(iout)(nstor(ndx+j),j=1,irtru)
c ...read geometry box-space data      *************
      ndx=kbiz-1
      read(iout) (nstor(ndx+j),j=1,irtru)
c     to read media data
      ndx=kmiz-1
      read(iout)(nstor(ndx+j),j=1,irtru)
   40 format(1h0/10x,36hmorse region in input zone(i) array ,
     1   14h (mriz(i),i=1,,i3,1h)//(10x,22i5/))
   50 format(1h0/10x,34hmorse media in input zone(i) array,
     1   14h (mmiz(i),i=1,,i3,1h)//(10x,22i5/))
   55 format(1h0/10x,37hmorse universe in input zone(i) array,
     1   14h (nbiz(i),i=1,,i3,1h)//(10x,22i5/))

      call geni(nstor(kma),nstor(kfpd),            nstor(klcr),
     1 nstor(knbd),nstor(kior),nstor(kriz),nstor(krcz),nstor(kmiz),
     2 nstor(kmcz),nstor(kkr1),nstor(kkr2),nstor(kbiz),nstor(kbcz),
     3 intt,iot,iout,ibodt)

      write(iot,40) irtru, (nstor(i),i=kriz,krf)
      write(iot,50) irtru,(nstor(i),i=kmiz,kmf)
      write(iot,55) irtru,(nstor(i),i=kbiz,klf)
      nir=1
      do 70 i=2,irtru
      indx=kriz-1+i
      im1=i-1
      do 60 j=1,im1
      jndx=kriz-1+j
      if(nstor(jndx).eq.nstor(indx)) go to 70
   60 continue
      nir=nir+1
   70 continue
      naad=kvol+nir*ispcd-kma
c    nadd is space used by cg in jomin2  **** should be  same as jomin1
        nadd=naad
      do 80 i=1,numr
      j=knsr-1+i
   80 nstor(j)=0
      call gtvlin(nstor(kfpd),nstor(kriz),nstor(kvol),ivopt,
     1   nir,irtru,iout,iot)
      return
      end


      subroutine risk(nl1,nl2,lp,ncmax,ni4,nbod,ncn,nba,xa,
     1 xb,d,ld,nxy,ma,fpd,ip)
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     92/06/22
c     time the module was last permanently updated:     09:01:30
c     programmer name:                                  j.t.west
c     module name:                                      marisk
c     current archiving level number:                   00001
c     current number of permanent updates:              00001
c     date of last access by librarian:                 93/03/15
c     dataset name:  x4s.scale4.master
c
c#######################################################################
c * * this routine determines local coordinate of location in
c     level nl2 , given a global coord in level nl1. nl1 .le. nl2.
      implicit real*8 (a-h,o-z)
      real*4 d
      common/ark/ngy,nll,nbb,nzy(3),xd(3)
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      dimension lp(*),xa(3),xc(3),ip(*),ncmax(*),
     1 ni4(*),nbod(*),ncn(*),nba(*),nxy(3),d(*),ld(*),ma(*),fpd(*),xb(3)
c...this routine is the opposite of clev
      do 5 i=1,3
   5  xb(i)=xa(i)
      if(nl2.lt.nl1) go to 777
      do 40 nlt=nl1,nl2
      ii=6*(nlt-1)
      nby=lp(ii+1)
      ngy=lp(ii+2)
      do 10 i=1,3
  10  nxy(i)=lp(ii+i+2)
      nll=lp(ii+6)
      if(nby) 777,20,15
  15  continue
      call trente(xb,ma,fpd,d,ld,ip,nq,nba,nby)
  20  continue
      n=3*(nll-1)+1
      nr=ncn(n)
      nx=-1
      call ctran(xb,ncmax(n),nxy,nx,ncn(n),d(nr),d,ld,ier)
      if(ngy) 40,777,25
  25  nbb=nbod(ngy)
      call rtexit(xb,ma,fpd,d,ld,ip,nq,nba,nbb)
  40  continue
      return
 777  continue
      call dipr(d,ip,2)
      call abend(ld)
      return
      end

      subroutine rtexit(xb,ma,fpd,d,ld,ip,np,nba,nbo)
c * * this routine rotates, then translates coordinates of
c     reference body.
      implicit real*8 (a-h,o-z)
      real*4 d
      dimension xb(3),ma(*),fpd(*),d(*),ld(*),ip(*),xo(3),
     1 xa(3),nba(*)
      call corner(nbo,xo,ma,fpd,d,ld,itype)
      if(itype.eq.9) go to 15
      nr=nba(nbo)
      if(nr.le.0) go to 30
      call rota(xa,xb,1,d(nr))
      do 10 i=1,3
  10  xb(i)=xa(i)
  15  continue
      do 20 i=1,3
  20  xb(i)=xb(i)+xo(i)
      return
  30  continue
      call dipr(d,ip,2)
      call errtra
      return
      end

      subroutine trente(xb,ma,fpd,d,ld,ip,np,nba,nbo)
c * * this routine translates, then rotates particle coordinate
c     called on entering or exiting an array or universe.
      implicit real*8 (a-h,o-z)
      real*4 d
      dimension xa(3),xb(3),ma(*),fpd(*),d(*),ld(*),ip(*),
     1 xo(3),nba(*)
      call corner(nbo,xo,ma,fpd,d,ld,itype)
      do 10 i=1,3
  10  xb(i)=xb(i)-xo(i)
      if(itype.eq.9) return
      nr=nba(nbo)
      if(nr.le.0) go to 30
      call rota(xa,xb,0,d(nr))
      do 15 i=1,3
  15  xb(i)=xa(i)
      return
  30  continue
      call dipr(d,ip,2)
      call errtra
      return
      end
