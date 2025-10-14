! *******************************************************
function getCosSpec(icenv,ip,s,r,d,e,g)
! ip: particle ID
! s: W index
! r: cut-off rigidity in GV
! d: atmospheric depth in g/cm2
! e: energy in MeV/n
! g: local geometry effect
! *******************************************************
implicit real*8 (a-h, o-z)

getCosSpec=0.0
if(icenv.le.-1) then ! SEP or TP mode
 if(ip.ne.1) then ! only proton data is available
  write(*,*) 'proj should be proton for SEP or TP mode'
  stop
 endif
 if(icenv.eq.-1) then ! SEP mode
  ispe=nint(g) ! SPE ID
  getCosSpec=getSEPfluence(ispe,r,e) ! SEP fluence (/MeV/cm2)
 else  ! TP mode
  getCosSpec=getTPflux(s,d,e) ! TP fluence (/MeV/cm2/s)
 endif
elseif(icenv.eq.0) then ! GCR mode
 if(ip.eq.0.or.ip.ge.29) then
  write(*,*) 'proj should be proton or ions for GCR-in-space mode'
  stop
 else
  getCosSpec=getGCRSpace(ip,s,r,d,e)
 endif
else  ! PARMA mode
 getCosSpec=getSpec(ip,s,r,d,e,g)
endif

return

end

! *******************************************************
function getTPflux(s,alti,e) ! get TP fluence (/MeV/cm2/s), only in PHITS
! s: solar modulation
! alti: altitude in km between 340 and 420 (d is altitude in this program)
! e: energy in MeV
! *******************************************************
implicit real*8 (a-h, o-z)
parameter(nalti=5) ! number of altitude database
parameter (nsor=2) ! solar minimum & maximum
parameter(nebin=300)
real*8, save:: eTP(nebin)
real*8, save:: TPflu(nebin,nalti,nsor)
real*8, save:: altidata(nalti)
character chatmp1*1
integer*4, save:: ifirst
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn
data ifirst/0/
real*8, save:: spot(nsor)
character, save:: minmax(nsor)*3

data ifirst/0/
data spot/0.0,150.0/
data minmax/'min','max'/

if(ifirst.eq.0) then
 ifirst=1
 do is=1,nsor
  open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/leo/tpflux-'//minmax(is)//'.inp',status='old')
  read(28,'(a)') chatmp1
  read(28,*) (altidata(ia),ia=1,nalti)
  do ie=1,nebin
   read(28,*) eTP(ie),(TPflu(ie,ia,is),ia=1,nalti)
  enddo
  close(815)
 enddo
 if(alti.lt.altidata(1).or.alti.gt.altidata(nalti)) then
  write(*,'("Warning: alti",f10.1," should be between",f6.1," and",f6.1," km")') alti,altidata(1),altidata(nalti)
  if(alti.lt.altidata(1)) then
   write(*,'("orbit altitude is set to ",f6.1," km")') altidata(1)
  else
   write(*,'("orbit altitude is set to ",f6.1," km")') altidata(nalti)
  endif
 endif
endif

do ie=1,nebin
 if(e.lt.eTP(ie)) exit
enddo
if(ie.eq.1) then
 eratio=0.0
 ie=2
elseif(ie.eq.nebin+1) then
 eratio=1.0
 ie=nebin
else
 eratio=(e-eTP(ie-1))/(eTP(ie)-eTP(ie-1))
endif

do ia=1,nalti
 if(alti.lt.altidata(ia)) exit
enddo
if(ia.eq.1) then
 aratio=0.0
 ia=2
elseif(ia.eq.nalti+1) then
 aratio=1.0
 ia=nalti
else
 aratio=(alti-altidata(ia-1))/(altidata(ia)-altidata(ia-1))
endif

do is=1,nsor
 if(s.lt.spot(is)) exit
enddo
if(is.eq.1) then
 sratio=0.0
 is=2
elseif(is.eq.nsor+1) then
 sratio=1.0
 is=nsor
else
 sratio=(s-spot(is-1))/(spot(is)-spot(is-1))
endif

Ftmp1=TPflu(ie-1,ia-1,is-1)+eratio*(TPflu(ie,ia-1,is-1)-TPflu(ie-1,ia-1,is-1))
Ftmp2=TPflu(ie-1,ia  ,is-1)+eratio*(TPflu(ie,ia  ,is-1)-TPflu(ie-1,ia  ,is-1))
Ftmp3=Ftmp1+aratio*(Ftmp2-Ftmp1)
Ftmp4=TPflu(ie-1,ia-1,is  )+eratio*(TPflu(ie,ia-1,is  )-TPflu(ie-1,ia-1,is  ))
Ftmp5=TPflu(ie-1,ia  ,is  )+eratio*(TPflu(ie,ia  ,is  )-TPflu(ie-1,ia  ,is  ))
Ftmp6=Ftmp4+aratio*(Ftmp5-Ftmp4)

getTPflux=Ftmp3+sratio*(Ftmp6-Ftmp3)

return

end


! *******************************************************
function getSEPfluence(ispe,r,e) ! get SEP fluence (/cm2/(MeV/n)), only in PHITS
! ispe: SPE type (-1:Feb.1956, -2:Nov.1960, -3:Aug1972, -4:Oct.1989, -5:Jan 2005)
! e: energy in MeV
! *******************************************************
implicit real*8 (a-h, o-z)
parameter(nSPE=5) ! number of SPE type
parameter(nebin=82)
real*8, save:: eSEP(nebin)
real*8, save:: SEPflu(nebin,nSPE)
character chatmp1*1
integer*4, save:: ifirst
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn
data ifirst/0/
data Emp/938.27d0/ ! mass of proton, nucleus mass is simply assumed to be A*Emp

if(ifirst.eq.0) then
 ifirst=1
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/spe/spefluence.txt',status='old')
 read(28,'(a)') chatmp1
 read(28,'(a)') chatmp1
 do ie=1,nebin
  read(28,*) eSEP(ie),(SEPflu(ie,is),is=1,nSPE)
 enddo
 close(815)
endif

if(ispe.lt.1.or.ispe.gt.nSPE) then
 write(*,'(''environ (='',i4,'') should be between 1 to '',i2,'' for SPE mode'')')
 stop
endif

ia=1
iz=1

RGV=getRfromE(iz,ia,e,Emp*ia)*0.001 ! rigidity in GV
if(RGV.lt.r) then ! sharp cutoff
 getSEPfluence=0.0
 return
endif

do ie=1,nebin
 if(e.lt.eSEP(ie)) exit
enddo
if(ie.eq.1) then
 ratio=0.0
 ie=2
elseif(ie.eq.nebin+1) then
 ratio=1.0
 ie=nebin
else
 ratio=(e-eSEP(ie-1))/(eSEP(ie)-eSEP(ie-1))
endif

getSEPfluence=SEPflu(ie-1,ispe)+(SEPflu(ie,ispe)-SEPflu(ie-1,ispe))*ratio

return

end

! *******************************************************
function getGCRSpace(iz,s,r,d,e) ! get Free Space GCR flux (/cm2/s/(MeV/n)), only in PHITS
! iz: charge of particle
! s: W index
! r: cut-off rigidity in GV
! d: altitude of orbit (0 for free space)
! e: energy in MeV/n
! *******************************************************
implicit real*8 (a-h, o-z)
parameter(npart=28) ! upto Ni
integer*4, save:: iAnum(npart)
data ianum/ 1, 4, 7, 9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,59/
data Emp/938.27d0/ ! mass of proton, nucleus mass is simply assumed to be A*Emp

ia=ianum(iz)

RGV=getRfromE(iz,ia,e,Emp*ia)*0.001 ! rigidity in GV
if(d.eq.0) then ! free space mode
 if(RGV.gt.r) then ! sharp cutoff
  getGCRSpace=getTOAspec(iz,ia,e,s)*1.0e-4*4.0*acos(-1.0) ! TOA spec is in (/m2/(MeV/n)/sec/sr)
 else
  getGCRSpace=0
 endif
else ! LEO mode
 getGCRSpace=getTOAspec(iz,ia,e,s)*1.0e-4*4.0*acos(-1.0) ! TOA spec is in (/m2/(MeV/n)/sec/sr)
 getGCRSpace=getGCRSpace*getGTF(RGV,d) ! consider geomagnetic transmission function
endif

return

end

! *******************************************************
function getGTF(r,alti) ! get geomagnetic transmission function 
! r: rigidity in GV
! alti: altitude
! *******************************************************
implicit real*8 (a-h, o-z)
parameter(nalti=5) ! number of altitude database
parameter(nrbin=30) ! number of rigidity bin
real*8, save:: rGTF(nrbin)
real*8, save:: GTF(nrbin,nalti)
real*8, save:: altidata(nalti)
character chatmp1*1
integer*4, save:: ifirst
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn
data ifirst/0/

if(ifirst.eq.0) then
 ifirst=1
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/leo/GTF.inp',status='old')
 read(28,'(a)') chatmp1
 read(28,*) (altidata(ia),ia=1,nalti)
 do ir=1,nrbin
  read(28,*) rGTF(ir),(GTF(ir,ia),ia=1,nalti)
 enddo
 close(815)
 if(alti.lt.altidata(1).or.alti.gt.altidata(nalti)) then
  write(*,'("Warning: alti",f10.1," should be between",f6.1," and",f6.1," km")') alti,altidata(1),altidata(nalti)
  if(alti.lt.altidata(1)) then
   write(*,'("orbit altitude is set to ",f6.1," km")') altidata(1)
  else
   write(*,'("orbit altitude is set to ",f6.1," km")') altidata(nalti)
  endif
 endif
endif

do ir=1,nrbin
 if(r.lt.rGTF(ir)) exit
enddo
if(ir.eq.1) then
 rratio=0.0
 ir=2
elseif(ir.eq.nrbin+1) then
 rratio=1.0
 ir=nrbin
else
 rratio=(r-rGTF(ir-1))/(rGTF(ir)-rGTF(ir-1))
endif

do ia=1,nalti
 if(alti.lt.altidata(ia)) exit
enddo
if(ia.eq.1) then
 aratio=0.0
 ia=2
elseif(ia.eq.nalti+1) then
 aratio=1.0
 ia=nalti
else
 aratio=(alti-altidata(ia-1))/(altidata(ia)-altidata(ia-1))
endif

Ftmp1=GTF(ir-1,ia-1)+rratio*(GTF(ir,ia-1)-GTF(ir-1,ia-1))
Ftmp2=GTF(ir-1,ia  )+rratio*(GTF(ir,ia  )-GTF(ir-1,ia  ))

getGTF=Ftmp1+aratio*(Ftmp2-Ftmp1)

return

end

! ******************************************************
function getshadowcos(d) ! get shadow angle (cos) from atmospheric depth d g/cm^2
! ******************************************************
parameter(maxUS=75) ! number of altitude bin for US-Standard Air 
implicit real*8 (a-h, o-z)
real*8, save:: altUS(maxUS) ! altitude data for US-Standard Air 1976 
real*8, save:: depUS(maxUS) ! atmospheric depth data for US-Standard Air 1976
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn
integer*4, save:: ifirst
character chatmp1*1
data Rearth/6378.13/  ! radius of earh in km
data ifirst/0/

if(ifirst.eq.0) then ! come to this routine first, so read input data
 ifirst=1
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/AtomDepth.inp',status='old')
 read(28,'(a)') chatmp1
 do ia=1,maxUS ! read US standard air 1976 data
  read(28,*) altUS(ia),depUS(ia)
 enddo
 close(815)
endif

do ia=1,maxUS
 if(depUS(ia).lt.d) exit
enddo
if(ia.eq.1) then ! out of range
 ratio=0.0
 ia=2
elseif(ia.le.maxUS) then
 ratio=(depUS(ia-1)-d)/(depUS(ia-1)-depUS(ia))
else
 ratio=1.0
 ia=maxUS
endif

alti=altUS(ia-1)+(altUS(ia)-altUS(ia-1))*ratio ! altitude in km

getshadowcos=sqrt(1.0d0-(Rearth/(Rearth+alti))**2)

return

end

! *******************************************************
function getSpec(ip,s,r,d,e,g)
! ip: particle ID
! s: W index
! r: cut-off rigidity in GV
! d: atmospheric depth in g/cm2
! e: energy in MeV/n
! g: local geometry effect
! *******************************************************
implicit real*8 (a-h, o-z)

getSpec=0.0
if(ip.eq.0) then  ! neutron
 getSpec=getNeutSpec(s,r,d,g,e)
elseif(ip.ge.1.and.ip.le.28) then ! proton to Ni
 getSpec=getIonSpec(ip,s,r,d,e)
elseif(ip.eq.29.or.ip.eq.30) then ! Muon
 iptmp=ip-28
 getSpec=getMuonSpec(iptmp,s,r,d,e)
else
 iptmp=ip-28 ! 3:electron, 4:positron, 5:photon
 getSpec=getsecondary(iptmp,s,r,d,e)
endif

return

end

! *******************************************************
function getFl(ip,s,r,d)  ! get Fl value, s:solar modulation potential, r:Cut off rigidity, d:depth
! *******************************************************
implicit real*8 (a-h, o-z)
parameter (npart=11) ! number of particle type, 0:neutron, 1:proton, 2:alpha, 3:electron, 4:positron, 5:photon, 6-11: Li-O
parameter (nBdata=4) ! B(1) - B(4) : Fl= B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
parameter (nAdata=10) ! A(1) - A(10) : Bmin = A(1)+A(2)*r+A(3)/(1+exp((r-A(4))/A(5))), Bmin = A(6)+A(7)*r+A(8)/(1+exp((r-A(9))/A(10)))
parameter (nsor=2) ! solar minimum & maximum
character chatmp1*1,chatmp5*5
character, save:: pname(0:npart)*6
real*8, save:: A(0:npart,nBdata,nAdata),B(nBdata)
real*8, save:: spot(nsor),FL(nsor) 
integer*4, save:: ifirst
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

data ifirst/0/
data spot/0.0,150.0/
data pname/'neutro','proton','alphaa','elemag','elemag','elemag', & 
         & 'ions  ','ions  ','ions  ','ions  ','ions  ','ions  '/

if(ifirst.eq.0) then
 ifirst=1
 do i=0,5
  if(i.le.2) then ! neutron, proton, alpha
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/'//pname(i)//'/Rigid-Dep.inp',status='old')
  elseif(i.eq.3) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/'//pname(i)//'/Rigid-Dep-EL.inp',status='old')
  elseif(i.eq.4) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/'//pname(i)//'/Rigid-Dep-PO.inp',status='old')
  elseif(i.eq.5) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/'//pname(i)//'/Rigid-Dep-PH.inp',status='old')
  endif
  read(28,'(a)') chatmp1
  do ib=1,nBdata
   read(28,1010) chatmp5,(A(i,ib,ia),ia=1,nAdata)
  enddo
  close(28)
 enddo
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/ions/Rigid-Dep.inp',status='old')
 do i=npart-5,npart
  read(28,'(a)') chatmp1
  do ib=1,nBdata
   read(28,1010) chatmp5,(A(i,ib,ia),ia=1,nAdata)
  enddo
 enddo
 close(28)
endif
1010	format(a5,30es13.5)

do is=1,nsor ! solar minimum and maximum
 if(ip.eq.1.or.ip.eq.2.or.ip.ge.npart-5) then ! need not correction
  r1=r 
 else ! for neutron, electron, positron and photon, Rc for high-altitude should be corrected 
  if(ip.eq.0) then
   ipidx=ip ! for neutron
  else
   ipidx=ip-2 ! 1:electron, 2:positron, 3:photon
  endif
  r1=r*getBestR(is,r,d,ipidx)
 endif
 do ib=1,nBdata
  if(is.eq.1) then ! solar minimum 
   B(ib)=A(ip,ib,1)+A(ip,ib,2)*r1+A(ip,ib,3)/(1+exp((r1-A(ip,ib,4))/A(ip,ib,5)))
  else
   B(ib)=A(ip,ib,6)+A(ip,ib,7)*r1+A(ip,ib,8)/(1+exp((r1-A(ip,ib,9))/A(ip,ib,10)))
  endif
 enddo
 Fl(is)=B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
enddo

if(ip.le.5) then
 pow=getPow(ip,d,r) ! get Power index
else
 pow=getPow(ip+2,d,r) ! in getPow, ip should be +2 for ions
endif

A2=(Fl(1)-Fl(2))/(getFFPfromW(spot(1))**pow-getFFPfromW(spot(2))**pow)
A1=Fl(1)-A2*getFFPfromW(spot(1))**pow
getFl=a1+a2*getFFPfromW(s)**pow

return
end

! **********************************************************
function getFFPfromW(s) ! get FFP (MV) from W (sun spot number)
! **********************************************************
implicit real*8 (a-h, o-z)
if(s.ge.0) then
 getFFPfromW=370.0+3.0e-1*s**1.45 ! FFP in MV
else
 getFFPfromW=370.0-3.0e-1*abs(s)**1.45 ! FFP in MV
endif
return
end

! **********************************************************
function getRfromE(iz,ia,Ek,Em) ! get Rigidity in MV from Kinetic Energy (MeV/n)
! **********************************************************
implicit real*8 (a-h, o-z)
getRfromE=sqrt((ia*Ek)**2+2*Ek*ia*Em)/iz
return
end

! **********************************************************
function getEfromR(iz,Rm,COR) ! get Kinetic Energy (MeV) from Rigidity (MV)
! **********************************************************
implicit real*8 (a-h, o-z)
getEfromR=sqrt((iz*COR)**2+Rm**2)-Rm
return
end

! *******************************************************
function getPow(ip,d,r)  ! get Power of solar modulation dependence, r:Cut off rigidity, d:depth
! *******************************************************
implicit real*8 (a-h, o-z)
parameter (npart=13) ! number of particle type, neutron, proton, alpha, photon, electron positron, mu+, mu-, Li-O
parameter (nBdata=2) ! B(1) - B(4) : Pow = b1 + b2*d
parameter (nAdata=5) ! A(1) - A(5) : B = A(1)+A(2)*r+A(3)/(1+exp((r-A(4))/A(5)))
character chatmp1*1,chatmp5*5
character, save:: pname(0:npart)*6
real*8, save:: A(0:npart,nBdata,nAdata),B(nBdata)
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn
integer*4, save:: ifirst

data ifirst/0/
data pname/'neutro','proton','alphaa','elemag','elemag','elemag','muon--','muon--', &
         & 'ions  ','ions  ','ions  ','ions  ','ions  ','ions  '/

if(ifirst.eq.0) then
 ifirst=1
 do i=0,7
  if(i.le.2) then ! neutron, proton, alpha
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/'//pname(i)//'/solar-dep.inp',status='old')
  elseif(i.eq.3) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/'//pname(i)//'/solar-dep-EL.inp',status='old')
  elseif(i.eq.4) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/'//pname(i)//'/solar-dep-PO.inp',status='old')
  elseif(i.eq.5) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/'//pname(i)//'/solar-dep-PH.inp',status='old')
  elseif(i.eq.6) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/'//pname(i)//'/solar-dep.plus',status='old')
  elseif(i.eq.7) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/'//pname(i)//'/solar-dep.mins',status='old')
  endif
  read(28,'(a)') chatmp1
  do ib=1,nBdata
   read(28,*) (A(i,ib,ia),ia=1,nAdata)
  enddo
  close(28)
 enddo
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/ions/solar-dep.inp',status='old')
 read(28,'(a)') chatmp1 
 do i=npart-5,npart
  do ib=1,nBdata
   read(28,*) (A(i,ib,ia),ia=1,nAdata)
  enddo
 enddo
 close(28)
endif

do ib=1,nBdata
 b(ib)=a(ip,ib,1)+a(ip,ib,2)*r+a(ip,ib,3)/(1.0+exp((r-a(ip,ib,4))/a(ip,ib,5)))
enddo

getpow=b(1)+b(2)*d
return
end



! **********************************************************
function getBestR(is,r,d,ip) ! get best Rc data for high-altitude correction
! **********************************************************
parameter(nBdata=6)
parameter(ndep=26)
parameter(npart=3) ! high-altitude correction is necessary only for electron, positron, photon 
implicit real*8 (a-h, o-z)
real*8, save:: A(0:npart,nBdata,ndep),B(nBdata)
real*8, save:: dep(ndep)
character chatmp1*1,chatmp5*5
character pname(npart)*2
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

data ifirst/0/
data pname/'EL','PO','PH'/

if(ifirst.eq.0) then
 a(:,:,:)=0.0
 do i=0,npart
  if(i.eq.0) then ! neutron
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/neutro/bestR.inp',status='old')
  else ! electron, positron, photon
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/elemag/bestR-'//pname(i)//'.inp',status='old')
  endif
  read(28,'(A)') chatmp1
  do id=1,ndep
   read(28,*) dep(id),(A(i,ib,id),ib=1,nBdata)
  enddo
  close(28)
 enddo
 ifirst=1
endif

do id=1,ndep
 if(d.lt.dep(id)) exit
enddo
if(id.eq.1) then
 ratio=0.0
 id=2
elseif(id.eq.ndep+1) then
 ratio=1.0
 id=ndep
else
 ratio=(d-dep(id-1))/(dep(id)-dep(id-1))
endif

do ib=1,nBdata
 B(ib)=A(ip,ib,id-1)+(A(ip,ib,id)-A(ip,ib,id-1))*ratio
enddo

if(is.eq.1) then ! solar minimum
 getBestR=10**(b(1)+b(2)*r+b(3)/r)
else
 getBestR=10**(b(4)+b(5)*r+b(6)/r)
endif

return

end

! *******************************************************
function getNeutspec(s,r1,d1,g,e) ! get neutron flux
!     s:Wolf number
!     r:cut off rigidity (GV)
!     d:air depth (g/cm^2)
!     g:local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
!     e:neutron energy (MeV)
! *******************************************************
implicit real*8 (a-h, o-z)
parameter(nA=12) ! number of basic spectrum parameter
parameter(nG=6) ! number of geometry parameter
real*8, save:: A(nA) ! basic spectrum parameter
real*8, save:: geo(nG) ! geometry parameter
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

data ifirst/0/
data airbus/2.45/  ! weight of airbus340 (100t)
data Eth/2.5e-8/  ! themal energy

if(ifirst.eq.0) then ! first time, get universal parameter (i.e: independent of all parameters)
!      Read A parameter
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/neutro/fitting-lowspec.inp',status='old')
 read(28,'(a)') chatmp1
 read(28,1001) (A(ia),ia=1,nA) !A(4)&A(12) is s,r,d-dependence, so will be changed
 close(unit=28)
 ifirst=1
endif
1001	format(30e13.5)

r=max(1.0,r1) ! secondary particle fluxes are the same for Rc<1GV
d=max(0.15,d1) ! secondary particle fluxes are the same for d < 0.15 g/cm2

!     get condition dependent parameters
Fl=getFl(0,s,r,d)
A(12)=getA12(r,d)
A(4)=getA4(r,d)
call getGpara(g,geo) ! obtain G parameters
!     calculate flux (/cm^2/s/lethargy)
x=e  ! x is energy

evap=a(1)*(x/a(2))**a(3)*exp(-x/a(2))
gaus=a(4)*exp(-(log10(x)-log10(a(5)))**2/(2*log10(a(6))**2))
conti=a(7)*log10(x/a(8))*(1+tanh(a(9)*log10(x/a(10))))*(1-tanh(a(11)*log10(x/a(12))))

basic=conti+evap+gaus

basic=basic*CorrNeut(s,r,d,e) ! correction for high altitude

fG=geo(1)+geo(2)*log10(x/geo(3))*(1-tanh(geo(4)*log10(x/geo(5))))
if(g.lt.0.0) then
 if(g.lt.-10.0) then ! cabin
  gtmp=g+10.0
 else ! pilot
  gtmp=g
 endif
 fG=fG*(abs(gtmp)-10*int(abs(gtmp)/10.0))/airbus  ! consider size of aircraft
endif
geofactor=10.0**fG
ther=geo(6)*(x/Eth)**2*exp(-(x/Eth))

getNeutspec=Fl*(basic*geofactor+ther)/e

return
end

! *******************************************************
subroutine getGpara(g,geo) ! get surroudning environment parameters
! *******************************************************
implicit real*8 (a-h, o-z)
parameter(nG=6) ! number of geometry parameter
dimension geo(nG) ! geometry parameter
real*8, save:: P(24) ! Input parameters read from input file, P(1)-P(3) from Geo-Dep, P(4)-P(14) from Water-Dep, P(15)-P(24) from aircraft-dep
character chatmp1*1,chatmp4*4
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

1001	format(a4,1x,30e13.5)

if(P(1).eq.0.0) then
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/neutro/Geo-Dep.inp',status='old')
 read(28,'(a)') chatmp1
 read(28,*) p(1),p(2),p(3)
 close(28)
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/neutro/Water-Dep.inp',status='old')
 read(28,'(a)') chatmp1
 read(28,1001) chatmp4,p(4),p(5),p(6)
 read(28,1001) chatmp4,p(7),p(8),p(9)
 read(28,1001) chatmp4,p(10),p(11),p(12),p(13),p(14)
 close(28)
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/neutro/Aircraft-Dep.inp',status='old')
 read(28,'(a)') chatmp1
 read(28,*) (p(ip),ip=15,19)  ! for pilot location
 read(28,*) (p(ip),ip=20,24)  ! for passenger & small aircraft configuration
 close(28)
endif

if(g.ge.10.0) then ! in semi-infite atmosphere
 do ig=1,nG
  geo(ig)=0.0
 enddo
 geo(3)=1.0  ! if geo(3)=0, the value should be in NaN
 geo(5)=1.0  ! if geo(5)=0, the value should be in NaN
elseif(g.ge.0.0) then ! for normal ground case
 geo(1)=p(1)
 geo(2)=p(2)
 geo(3)=10.0**(p(4)+p(5)/(p(6)+g))
 geo(4)=p(3)
 geo(5)=p(7)+p(8)*g+p(9)*g**2
 geo(6)=(p(10)+p(11)*exp(-p(12)*g))/(1+p(13)*exp(-p(14)*g))
else   ! pilot or cabin location
 is=14
 if(g.le.-10.0) is=is+5 ! for passenger & small aircraft configuration, skip 5 more data
 do i=1,5
  geo(i)=p(is+i)
 enddo
 geo(6)=0.0  ! no thermal component
endif
return
end

! *******************************************************
function getA4(r,d)  ! get A4 value, s:solar modulation potential, r:Cut off rigidity, d:depth
! *******************************************************
implicit real*8 (a-h, o-z)
parameter (nBdata=4) ! B(1) - B(4) : Fl= B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
parameter (nAdata=6) ! A(1) - A(6) : B = A(1)+A(2)/(1+exp((r-A(3))/A(4))), A(5):A(1) for APmax, A(6):A(3) for APmax
character chatmp1*1,chatmp5*5
real*8, save:: A(nAdata),B(nBdata)
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

if(B(2).eq.0.0) then ! first time
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/neutro/Depth-Dep-mid.out',status='old')  ! read B(2)-B(4) (independent of s,r,d)
 read(28,'(a)') chatmp1
 read(28,'(a)') chatmp1
 read(28,*) tmp1,tmp2,B(2),B(3),B(4)
 close(28)
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/neutro/Rigid-Dep.inp',status='old')
 do i=1,5  ! skip 5 line, 1 header line + 4 Bdata
  read(28,'(a)') chatmp1
 enddo
 read(28,1010) chatmp5,(A(ia),ia=1,nAdata)
 close(28)
endif
1010	format(a5,30e13.5)	

B(1)=A(1)+A(2)*r+A(3)/(1+exp((r-A(4))/A(5)))
getA4=B(1)+B(2)*d/(1+B(3)*exp(B(4)*d))

return
end

! *******************************************************
function getA12(r,d)  ! get Fl value, s:solar modulation potential, r:Cut off rigidity, d:depth
! *******************************************************
implicit real*8 (a-h, o-z)
parameter (nBdata=4) ! B(1) - B(4) : Fl= B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
parameter (nAdata=6) ! A(1) - A(6) : B = A(1)+A(2)/(1+exp((r-A(3))/A(4))), A(5):A(1) for APmax, A(6):A(3) for APmax
character chatmp1*1,chatmp5*5
real*8, save:: A(nBdata,nAdata),B(nBdata)
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

if(B(4).eq.0.0) then ! first time
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/neutro/Depth-Dep-hig.out',status='old')  ! read B(2),B(4) (independent of s,r,d)
 read(28,'(a)') chatmp1
 read(28,'(a)') chatmp1
 read(28,*) tmp0,tmp1,tmp2,tmp3,B(4)
 close(28)
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/neutro/Rigid-Dep.inp',status='old')
 do i=1,6  ! skip 5 line, 1 header line + 5 Bdata
  read(28,'(a)') chatmp1
 enddo
 do ib=1,3 ! for B1 to B3
  read(28,1010) chatmp5,(A(ib,ia),ia=1,nAdata)
 enddo
 close(28)
endif
1010	format(a5,30es13.5)

do ib=1,3
 B(ib)=A(ib,1)+A(ib,2)*r+A(ib,3)/(1+exp((r-A(ib,4))/A(ib,5)))
enddo

getA12=B(1)*(exp(-B(2)*d)+B(3)*exp(-B(4)*d))

return
end

! *******************************************************
function CorrNeut(s,r,d,e)  ! get correction factor for high altitude data, only solar minimum and maximum
! *******************************************************
parameter(nhensu=9) ! number of hensu
parameter(mpara=5) ! number of parameters to COR dependence
parameter(ndep=26) ! number of depth
parameter(nsol=2) ! solar minimum and maximum
implicit real*8 (a-h,o-z)
dimension ainp(mpara,nhensu,ndep,nsol)
dimension dep(ndep)
dimension b(nhensu),c(nsol) ! temporary dimension
dimension spot(nsol)
character chatmp1*1
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

save ainp,spot,ifirst,dep

data ifirst/0/
data spot/0.0,150.0/

if(ifirst.eq.0) then ! first time call 
!     Read parameters for depth-independent parameters
! open(28,file='correction-depth-rigid-final.inp',status='old') ! only high altidue correction mode
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/neutro/correction-depth-rigid.inp',status='old') ! all altitude correction mode
 read(28,'(a)') chatmp1
 do ih=1,nhensu
  do id=1,ndep
   read(28,*) itmp,dep(id),((ainp(ip,ih,id,is),ip=1,mpara),is=1,nsol)
  enddo
 enddo
 close(28)
 ifirst=1
endif

! find depth ID
do id=1,ndep
 if(d.lt.dep(id)) exit
enddo
if(id.eq.1) then
 ratio=0.0
 id=2
elseif(id.eq.ndep+1) then
 ratio=1.0
 id=ndep
else
 ratio=(d-dep(id-1))/(dep(id)-dep(id-1))
endif

rc=max(1.0,r)
do is=1,nsol
 do ih=1,nhensu ! determine 9 hensu used in the correction equation
  d1=ainp(1,ih,id-1,is)+Ainp(2,ih,id-1,is)*rc+Ainp(3,ih,id-1,is)/(1+exp((rc-Ainp(4,ih,id-1,is))/Ainp(5,ih,id-1,is)))
  d2=ainp(1,ih,id-0,is)+Ainp(2,ih,id-0,is)*rc+Ainp(3,ih,id-0,is)/(1+exp((rc-Ainp(4,ih,id-0,is))/Ainp(5,ih,id-0,is)))
  b(ih)=d1+(d2-d1)*ratio
 enddo
 C(is)=10**(b(1)+(b(2)*log10(e)+b(3))*(1-tanh(b(4)*log10(e/b(5))))+(b(6)*log10(e)+b(7))*(1+tanh(b(8)*log10(e/b(9)))))
enddo

ip=0 ! always neutron 
pow=getPow(ip,d,r) ! get Power index
A2=(c(1)-c(2))/(getFFPfromW(spot(1))**pow-getFFPfromW(spot(2))**pow)
A1=c(1)-A2*getFFPfromW(spot(1))**pow
CorrNeut=a1+a2*getFFPfromW(s)**pow

return
end


! ******************************************************
function getMuonSpec(ip,s,r1,d1,e) ! get muon spectrum, ip=30:mu+, =31:mu-
! ******************************************************
implicit real*8 (a-h, o-z)	
parameter (nsol=2) ! solar minimum and maximum
parameter (nAdata=7) ! number of A parameter

dimension Acurr(nAdata)
dimension Fl(nsol),spot(nsol)

data restmass/105.6/
data spot/0.0,150.0/
data ethre/3.0e5/ ! threshold energy for high-energy muon correction

if(e.lt.1.0e-2) then
 getMuonSpec=0.0
 return
endif

r=max(1.0,r1) ! muon fluxes are the same for Rc < 1 GV
d=max(0.15,d1) ! secondary particle fluxes are the same for d < 0.15 g/cm2

beta=sqrt(1.0-(restmass/(restmass+e))**2)
tmp=max(2.0,log10(e)) ! below 100 MeV, this value should be constant

do is=1,nsol
 iptmp=ip ! 1 for mu+, 2 for mu+
 call getAmuon(Acurr,iptmp,is,d,r)
 if(e.gt.ethre) then
  acurr(5)=acurr(5)+0.4 ! /(1+exp((5.8088849-tmp)/0.24867240)) ! high energy correction, see fit/muon/PowerCorrection
  acurr(1)=acurr(1)*ethre**0.4 ! /(1+exp((5.8088849-tmp)/0.24867240))
 endif
 Fl(is)=acurr(1)*(e+(acurr(2)+acurr(4)*tmp)/beta**acurr(3))**(-acurr(5))*(1+exp(-acurr(6)*(log(e)+acurr(7))))
enddo

iptmp=ip+5 ! 6 for mu+, 7 for mu-
pow=getPow(iptmp,d,r) ! get Power index
A2=(Fl(1)-Fl(2))/(getFFPfromW(spot(1))**pow-getFFPfromW(spot(2))**pow)
A1=Fl(1)-A2*getFFPfromW(spot(1))**pow
getMuonSpec=a1+a2*getFFPfromW(s)**pow

return
end

! ******************************************************
subroutine getAmuon(Acurr,ip,is,d,r) ! get A parameter 
! ******************************************************
implicit real*8 (a-h, o-z)
parameter (npart=2) ! mu+ and mu-
parameter (nAdata=7) ! number of A parameter A(1) to A(7)
parameter (nBdata=10) ! A_min = B1+B2*r+B3/(1+exp((r-B4)/B5), A_max = B6+B7*r+B8/(1+exp((r-B9)/B10)
parameter (ndep=26) ! number of depth
dimension Acurr(nAdata)
dimension Bdata(nAdata,nBdata,npart,ndep),dep(ndep)
character chatmp1*1,chatmp4*4
character charge(npart)*4
save Bdata,dep
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

data charge/'plus','mins'/
data ifirst/0/

if(ifirst.eq.0) then  ! first time called this routine
 ifirst=1
 do ip2=1,npart
  open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/muon--/final135.'//charge(ip2),status='old') ! read A(1),A(3),A(5), they are solar independent
  read(28,*) chatmp1
  do id2=1,ndep
   read(28,*) dep(id2),Bdata(1,1,ip2,id2),Bdata(3,1,ip2,id2),Bdata(5,1,ip2,id2)
  enddo
  close(28)
  open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/muon--/final2467.'//charge(ip2),status='old') ! read A(2),A(4),A(6),A(7)
  read(28,*) chatmp1
  do id2=1,ndep
   read(28,*) dep(id2),(Bdata(2,ib,ip2,id2),ib=1,nBdata)
  enddo
  read(28,*) chatmp1
  do id2=1,ndep
   read(28,*) dep(id2),(Bdata(4,ib,ip2,id2),ib=1,nBdata)
  enddo
  read(28,*) chatmp1
  do id2=1,ndep
   read(28,*) dep(id2),(Bdata(6,ib,ip2,id2),ib=1,nBdata)
  enddo
  read(28,*) chatmp1
  do id2=1,ndep
   read(28,*) dep(id2),(Bdata(7,ib,ip2,id2),ib=1,nBdata)
  enddo
  close(28)
 enddo
endif

! find closest depth
do id=1,ndep
 if(d.lt.dep(id)) exit
enddo
if(id.eq.1) then
 ratio=0.0
 id=2
elseif(id.eq.ndep+1) then
 ratio=1.0
 id=ndep
else
 ratio=(d-dep(id-1))/(dep(id)-dep(id-1))
endif

rc=max(1.0,r) ! minimum Rc = 1.0GV

do ia=1,nAdata
 if(ia.eq.1) then
  Acurr(ia)=log(Bdata(ia,1,ip,id-1))+(log(Bdata(ia,1,ip,id))-log(Bdata(ia,1,ip,id-1)))*ratio ! log-interpolation
  Acurr(ia)=exp(Acurr(ia))
 elseif(ia.eq.3.or.ia.eq.5) then
  Acurr(ia)=Bdata(ia,1,ip,id-1)+(Bdata(ia,1,ip,id)-Bdata(ia,1,ip,id-1))*ratio
 else
  if(is.eq.1) then
   B1=Bdata(ia,1,ip,id-1)+Bdata(ia,2,ip,id-1)*rc+Bdata(ia,3,ip,id-1)/(1+exp((rc-Bdata(ia,4,ip,id-1))/Bdata(ia,5,ip,id-1)))
   B2=Bdata(ia,1,ip,id)+Bdata(ia,2,ip,id)*rc+Bdata(ia,3,ip,id)/(1+exp((rc-Bdata(ia,4,ip,id))/Bdata(ia,5,ip,id)))
  else
   B1=Bdata(ia,6,ip,id-1)+Bdata(ia,7,ip,id-1)*rc+Bdata(ia,8,ip,id-1)/(1+exp((rc-Bdata(ia,9,ip,id-1))/Bdata(ia,10,ip,id-1)))
   B2=Bdata(ia,6,ip,id)+Bdata(ia,7,ip,id)*rc+Bdata(ia,8,ip,id)/(1+exp((rc-Bdata(ia,9,ip,id))/Bdata(ia,10,ip,id)))
  endif 
  Acurr(ia)=B1+(B2-B1)*ratio
 endif
enddo

return
end


! *******************************************************
function getIonSpec(iz,s,r,d,e) ! get Ion flux
!     s:Wolf number
!     r:cut off rigidity (GV)
!     d:air depth (g/cm^2)
!     e:ion energy (MeV/n)
!     imode: include (=0) or exclude (=1) secondary particle
! *******************************************************
implicit real*8 (a-h, o-z)
parameter(nAdata=6)
parameter(npart=28)
parameter(ngroup=6) 

character chatmp1*1

real*8, save:: A(nAdata,ngroup) ! parameter used in combine.for
integer*4, save:: iAnum(npart)
integer*4, save:: ifirst

integer, save:: igidx(npart)    ! group index
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

data igidx/1,2,3,3,3,4,4,4,4,5 &
        & ,5,5,5,5,5,5,5,5,5,6 &  
        & ,6,6,6,6,6,6,6,6/

data ianum/ 1, 4, 7, 9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,59/

data ifirst/0/
data restmass/938.27d0/

if(ifirst.eq.0) then  ! first time call
 ifirst=1
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/ions/Combine.inp',status='old')
 read(28,*) chatmp1
 do i=1,ngroup
  read(28,*) (A(ia,i),ia=1,nAdata)
 enddo
 close(28)
endif

if(e.lt.1.0e-2) then  ! for lower energy, no output
 getIonSpec=0.0
 return
endif

x=e ! x is energy
ig=igidx(iz)

tmp=restmass*iAnum(iz)

Ecut=getEfromR(iZ,tmp,r*1000.0)/iAnum(iz)-a(6,ig)*d  ! MeV/n
EcPri=max(a(1,ig),Ecut*a(3,ig))
EcSec=max(a(2,ig),Ecut*a(3,ig))

if(iz.le.2) then
 ip=iz
else
 ip=iz+3 ! in getsecondary, ip=iz+3 (electron, positron, photon) for ions
endif

getIonSpec=getPrimary(iz,iAnum(iz),s,d,x)*0.5*(tanh(a(4,ig)*(x/EcPri-1))+1.0) &
        & +getsecondary(ip,s,r,d,x)*0.5*(tanh(a(5,ig)*(1-x/EcSec))+1.0)

return

end

! **********************************************************
function getsecondary(ip,s,r1,d1,e) ! get secondary particle flux (/cm^2/s/MeV)
! **********************************************************
parameter(nBdata=8)
parameter(ndep=26)
parameter(npart=11) ! proton, alpha, electron, positron, photon, Li, Be, B, C, N, O
implicit real*8 (a-h, o-z)
real*8, save:: A(nBdata,npart,ndep),B(nBdata)
real*8, save:: dep(ndep)
integer, save:: ifirst
character chatmp1*1,chatmp5*5
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

data ifirst/0/

r=max(1.0,r1) ! to get Fl, minimum Rc=1GV
d=max(0.15,d1) ! secondary particle fluxes are the same for d < 0.15 g/cm2

if(ifirst.eq.0) then
 do i=1,5 ! proton, alpha, electron, positron, photon
  if(i.eq.1) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/proton/fitting-lowspec.inp',status='old')
  elseif(i.eq.2) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/alphaa/fitting-lowspec.inp',status='old')
  elseif(i.eq.3) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/elemag/fitting-lowspec-EL.inp',status='old')
  elseif(i.eq.4) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/elemag/fitting-lowspec-PO.inp',status='old')
  elseif(i.eq.5) then
   open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/elemag/fitting-lowspec-PH.inp',status='old')
  endif  
  read(28,'(A)') chatmp1
  do id=1,ndep
   read(28,*) dep(id),(A(ib,i,id),ib=1,nBdata)
  enddo
  close(28)
 enddo
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/ions/fitting-lowspec.inp',status='old')
 read(28,'(a)') chatmp1
 do i=npart-5,npart
  read(28,*) (A(ib,i,1),ib=1,nBdata)
  do ib=1,nBdata
   do id=1,ndep
    A(ib,i,id)=A(ib,i,1) ! for Li to O, depth independent
   enddo
  enddo
 enddo
 close(28)
 ifirst=1
endif

if(e.lt.1.0e-2.or.ip.gt.npart) then  ! for lower energy or heavier ions, no output
 getsecondary=0.0
 return
endif

do id=1,ndep
 if(d.lt.dep(id)) exit
enddo
if(id.eq.1) then
 ratio=0.0
 id=2
elseif(id.eq.ndep+1) then
 ratio=1.0
 id=ndep
else
 ratio=(d-dep(id-1))/(dep(id)-dep(id-1))
endif

do ib=1,nBdata
 B(ib)=A(ib,ip,id-1)+(A(ib,ip,id)-A(ib,ip,id-1))*ratio
enddo

if(ip.eq.1.or.ip.eq.2) then ! proton or alpha
 getsecondary=getFl(ip,s,r,d)*(b(1)*e**b(2))/(1+b(3)*e**b(4))/(1+b(5)*e**b(6))*(1+exp(-b(7)*(log(e)+b(8))))
elseif(ip.eq.3.or.ip.eq.4) then ! electron or positron
 getsecondary=getFl(ip,s,r,d)*(b(1)*e**b(2))/(1+b(3)*e**b(4))/(1+b(5)*e**b(6))
elseif(ip.eq.5) then ! photon
 getsecondary=getFl(ip,s,r,d)*(b(1)*e**b(2))*(1+b(3)*e**b(4))/(1+b(5)*e**b(6))/(1+exp(-b(7)*(log(e)+b(8))))
else ! Li,Be,B,C,N,O
 getsecondary=getFl(ip,s,r,d)*(b(1)*e**b(2))/(1+b(3)*e**b(4))/(1+b(5)*e**b(6))/(1+exp(-b(7)*(log(e)+b(8))))
endif

return
	
end

! **********************************************************
function getTOAspec(iZ,iA,Ek,Spot) ! get TOA spectrum in (/(MeV/n)/s/m^2/sr)
!     Ek: Kinetic Energy in MeV/n
!     Spot: Wolf number estimated from count rate of neutron monitors
! **********************************************************
parameter(npart=28)
implicit real*8 (a-h,o-z)
real*8, save:: Dpara(npart),alpha(npart),gamma(npart),bpara(npart) ! Data for Nymmik Model
data Emp/938.27d0/ ! mass of proton, nucleus mass is simply assumed to be A*Emp

data Dpara/1.85e4,3.69e3,19.50,17.70,49.20,103.00,36.70,87.40, &
     &           3.19,16.40,4.43,19.30,4.17,13.40,1.15,3.06,1.30, &
     &           2.33,1.87,2.17,0.74,2.63,1.23,2.12,1.14,9.32,0.10,0.48/ ! ISO-Model taken from Matthia-ASR2013 

data alpha/2.85,3.12,3.41,4.30,3.93,3.18,3.77,3.11,4.05,3.11,3.14, &
     &           3.65,3.46,3.00,4.04,3.30,4.40,4.33,4.49,2.93,3.78,3.79, &
     &           3.50,3.28,3.29,3.01,4.25,3.52/

data gamma/2.74,2.77,2.82,3.05,2.96,2.76,2.89,2.70,2.82,2.76,2.84, &
     &           2.70,2.77,2.66,2.89,2.71,3.00,2.93,3.05,2.77,2.97,2.99, &
     &           2.94,2.89,2.74,2.63,2.63,2.63/
 
! ***** Determine LIS spectra based on DLR Model *****************
R=getRfromE(iz,ia,Ek,Emp*iA)*0.001 ! rigidity in GV
beta=sqrt(1-(Emp*iA/(Emp*iA+Ek*iA))**2)
dR2dE=0.001/iZ/beta*iA ! convert GV to MV, MeV to MeV/n
SpecLIS=Dpara(iz)*beta**alpha(iz)/R**gamma(iz)*dR2dE
! *******************************************************************

! consider solar modulation
R0=getFFPfromW(Spot)*0.001 ! FFP in GV from W, taken from Matthia-ASR2013
delta=0.02*Spot+4.7        ! taken from Matthia-ASR2013
getTOAspec=SpecLIS*(R/(R+R0))**delta

return
end


! *************************************************************
function getPrimary(iz,ia,s,d,e)  ! get Primary Flux
! *************************************************************
parameter(nAdata=3)
parameter(npart=28)
parameter(nepoint=6)
parameter(ngroup=6)
parameter(ndep=26)
implicit real*8 (a-h, o-z)	
real*8, save:: A(nAdata,ngroup,ndep)
real*8, save:: dEdxTable(npart,nepoint),epoint(nepoint) ! dE/dx for each particle
real*8, save:: dep(ndep) ! depth (g/cm2)
real*8, save:: down
integer, save:: igidx(npart)    ! group index
dimension B(nAdata) ! temporary used dimension

character gname(ngroup)*2
character chatmp1*1
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

data down/0.0/

data gname/'H-','He','Be','N-','Si','Fe'/

data igidx/1,2,3,3,3,4,4,4,4,5 &
        & ,5,5,5,5,5,5,5,5,5,6 &  
        & ,6,6,6,6,6,6,6,6/

if(down.eq.0.0) then ! first time call this routine
 down=1.720313d0
 do ig=1,ngroup
  open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/ions/primary-'//gname(ig)//'.inp',status='old')
  read(28,*) chatmp1
  do id=1,ndep
   read(28,*) dep(id),(a(i,ig,id),i=1,nAdata)
  enddo
  close(28)
 enddo
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/ions/dEdx-table.inp',status='old')
 read(28,'(a1)') chatmp1
 read(28,'(a1)') chatmp1
 read(28,*) (epoint(ie),ie=1,nepoint)
 do ip=1,npart
  read(28,*) itmp,itmp,(dEdxTable(ip,ie),ie=1,nepoint)
 enddo
 close(28)
endif

ig=igidx(iz)

do id=1,ndep
 if(d.lt.dep(id)) exit
enddo
if(id.eq.1) then
 ratio=0.0
 id=2
elseif(id.eq.ndep+1) then
 ratio=1.0
 id=ndep
else
 ratio=(d-dep(id-1))/(dep(id)-dep(id-1))
endif

do i=1,nAdata
 b(i)=a(i,ig,id-1)+(a(i,ig,id)-a(i,ig,id-1))*ratio
enddo

! find dEdx
do ie=2,nepoint-1
 if(e.le.epoint(ie)) exit
enddo
ratio=(log(e)-log(epoint(ie-1)))/(log(epoint(ie))-log(epoint(ie-1))) ! log-logg interpolation
ratio=max(0.0,min(1.0,ratio))
tmp=log(dEdxTable(iz,ie-1))+(log(dEdxTable(iz,ie))-log(dEdxTable(iz,ie-1)))*ratio ! log-log interpolation
dEdx=exp(tmp)

Eini=e+dEdx*d  ! Energy at the TOA
getPrimary=getTOAspec(iz,ia,Eini,s)*(b(1)*exp(-b(2)*d)+(1.0-b(1))*exp(-b(3)*d))
getPrimary=getPrimary*4.0*acos(-1.0)*1.0e-4/down  ! convert (/(MeV/n)/s/m^2/sr) to (/(MeV/n)/s/cm^2)

return

end

! *************************************************************
function get511flux(s,r,d)  ! get 511 keV photon flux in (/cm2/s)
! *************************************************************
parameter(ndep=26)
implicit real*8 (a-h, o-z)	
real*8, save:: F511(ndep),dep(ndep)
character chatmp1*1
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

data ifirst/0/

if(ifirst.eq.0) then
 ifirst=1
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/elemag/flux511keV.inp')
 read(28,'(a1)') chatmp1
 do id=1,ndep
  read(28,*) dep(id),F511(id)
 enddo
 close(28)
endif

do id=1,ndep
 if(d.lt.dep(id)) exit
enddo
if(id.eq.1) then
 ratio=0.0
 id=2
elseif(id.eq.ndep+1) then
 ratio=1.0
 id=ndep
else
 ratio=(d-dep(id-1))/(dep(id)-dep(id-1))
endif

Fratio=F511(id-1)+(F511(id)-F511(id-1))*ratio

iptmp=5 ! photon index
ene=0.511 ! 511 keV
get511flux=getsecondary(iptmp,s,r,d,ene)*Fratio

return

end


! ******************************************************
function getd(alti,cido) ! getd in g/cm^2, alti in km
! if -90 < cido < 90, use MSIS database
! else, use US standard air 1976
! ******************************************************
parameter(iMSIS=0)  ! 0:US standard atmosphere, 1:NRLMSISE database
parameter(maxUS=75) ! number of altitude bin for US-Standard Air 
parameter(maxMSIS=129) ! number of altitude bin for NRLMSISE-00 
parameter(maxlat=36) ! number of latitude bin for NRLMSISE-00
implicit real*8 (a-h, o-z)
real*8, save:: altUS(maxUS) ! altitude data for US-Standard Air 1976 
real*8, save:: altMSIS(maxMSIS) ! altitude data for NRLMSISE-00
real*8, save:: depUS(maxUS) ! atmospheric depth data for US-Standard Air 1976
real*8, save:: depMSIS(maxMSIS,maxlat) ! atmospheric depth data for each altitude & latitude for NRLMSISE-00
real*8, save:: glat(maxlat) ! latitude data
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn
integer*4, save:: ifirst
character chatmp1*1
data ifirst/0/

if(ifirst.eq.0) then ! come to this routine first, so read input data
 ifirst=1
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/AtomDepth.inp',status='old')
 read(28,'(a)') chatmp1
 do ia=1,maxUS ! read US standard air 1976 data
  read(28,*) altUS(ia),depUS(ia)
 enddo
 read(28,'(a)') chatmp1
 read(28,*) (glat(ido),ido=1,maxlat)
 do ia=1,maxMSIS ! read NRLMSISE-00 data
  read(28,*) altMSIS(ia),(depMSIS(ia,ido),ido=1,maxlat)
 enddo
 close(815)
endif

if((cido.lt.-90.01.or.cido.gt.90.01).or.iMSIS.eq.0) then ! US standard atmosphere 1976 mode
 do ia=1,maxUS
  if(altUS(ia).gt.alti) exit
 enddo
 if(ia.eq.1) then ! out of range
  write(6,*) 'Error in function getd'
  write(6,*) 'Altitude =',alti,' (km) should be higher than',altUS(1),' (km)'
  stop
 endif
 if(ia.eq.maxUS+1) then ! out of range
  write(6,*) 'Warning in function getd'
  write(6,*) 'Altitude =',alti,' (km) is too high. It is assumed to be',altUS(maxUS),' (km)'
  getd=depUS(maxUS)
  return
 endif
 ratio=(alti-altUS(ia-1))/(altUS(ia)-altUS(ia-1))
 getd=depUS(ia-1)+ratio*(depUS(ia)-depUS(ia-1))
else ! latitude is specified, so use NRLMSISE-00 data
 do ia=1,maxMSIS
  if(altMSIS(ia).gt.alti) exit
 enddo
 if(ia.eq.1) then ! out of range
  write(6,*) 'Error in function getd'
  write(6,*) 'Altitude =',alti,' (km) should be higher than',altMSIS(1),' (km)'
  stop
 endif
 if(ia.eq.maxMSIS+1) then ! out of range
  write(6,*) 'Error in function getd'
  write(6,*) 'Altitude =',alti,' (km) is too high. It is assumed to be',altMSIS(maxMSIS),' (km)'
  alti=altMSIS(maxMSIS)
  ia=maxMSIS
 endif
 ratio=(alti-altMSIS(ia-1))/(altMSIS(ia)-altMSIS(ia-1))
 do ido=2,maxlat-1 
  if(glat(ido).gt.cido) exit
 enddo
 ratio1=min(1.0,max(0.0,(cido-glat(ido-1))/(glat(ido)-glat(ido-1))))
 dep1=depMSIS(ia-1,ido-1)+ratio*(depMSIS(ia,ido-1)-depMSIS(ia-1,ido-1))
 dep2=depMSIS(ia-1,ido  )+ratio*(depMSIS(ia,ido  )-depMSIS(ia-1,ido  ))
 getd=dep1+ratio1*(dep2-dep1)
endif
return
end

! ******************************************************
function getr(cido,ckei) ! cido and ckei are center of ido&keido of each grid
! ******************************************************
implicit real*8 (a-h, o-z)
real*8, save:: cordata(181,361) ! maximum 1 deg step, +1 mean
real*8, save:: dpido(181),dpkei(361) ! data point ido & keido
character chatmp1*1
integer*4, save:: mkei,mido
real*8, save:: skei,sido
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn
integer*4, save:: ifirst
data ifirst/0/

if(ifirst.eq.0) then ! come to this routine first, so read input data
 ifirst=1
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/CORdata.inp',status='old')
 read(28,*) mkei,mido ! read step size
 read(28,*) chatmp1
 skei=360.0/(mkei-1)  ! keido step
 sido=180.0/(mido-1)  ! ido step
 do id=mido,1,-1     ! read from 90 to -90 deg
  do ik=mkei,1,-1      ! read from 180 to -180 deg
   read(28,*) dpkei(ik),dpido(id),cordata(id,ik)
  enddo
 enddo
 close(28)
endif

! **** Determine Cut-off Rigidity ***********
if(ckei.gt.180.and.ckei.lt.360) ckei=ckei-360 ! ckei should be given in Western longitude
id=min(mido-1,int((cido+90.0)/sido)+1)   ! lower ido bin
ik=min(mkei-1,int((ckei+180.0)/skei)+1)  ! lower keido bin
cor1=cordata(id,ik)*(dpido(id+1)-cido)/(dpido(id+1)-dpido(id))+cordata(id+1,ik)*(cido-dpido(id))/(dpido(id+1)-dpido(id))
cor2=cordata(id,ik+1)*(dpido(id+1)-cido)/(dpido(id+1)-dpido(id))+cordata(id+1,ik+1)*(cido-dpido(id))/(dpido(id+1)-dpido(id))
getr=cor1*(dpkei(ik+1)-ckei)/(dpkei(ik+1)-dpkei(ik))+cor2*(ckei-dpkei(ik))/(dpkei(ik+1)-dpkei(ik))

return

end



! subroutines for getting angular distribution
function getSpecAngFinal(ip,s,r,d,e,g,ang)
! ip: particle ID '1:neutro','2:proton','3:he---4','4:muon--','5:elepos','6:photon'
! s: W index
! r: cut-off rigidity in GV
! d: atmospheric depth in g/cm2
! e: energy in MeV/n
! g: local geometry effect
! ang: cos(theta)
implicit real*8 (a-h, o-z)
data coshigh/0.707106781186547d0/ ! 1/sqrt(2) (45 degree)
data coslow/0.50d0/               ! 1/2       (60 degree)

if(ip.eq.4.and.g.ge.0.0d0.and.g.le.1.0d0.and.ang.lt.coshigh) then ! ground level muon, condition for angle above 45 deg is introduced in 2024/11/10
 emin=1.1535d4 ! minimum energy for correction data
 if(e.ge.emin) then ! full correction
  getSpecAngFinal=getSpecAng(ip,s,r,d,e,g,1.0d0)*getGmuon(e,ang)
 else
  ratio=(getSpecAng(ip,s,r,d,emin,g,1.0d0)*getGmuon(emin,ang))/getSpecAng(ip,s,r,d,emin,g,ang)
  getSpecAngFinal=getSpecAng(ip,s,r,d,e,g,ang)*ratio
 endif
 if(ang.gt.coslow) then ! Partial correction for angle between 45 and 60 2024/11/10
  ratio=(1.0/ang-1.0/coshigh)/(1.0/coslow-1.0/coshigh)  ! interpolated with 1/cos
  getSpecAngFinal=getSpecAngFinal*ratio+getSpecAng(ip,s,r,d,e,g,ang)*(1.0-ratio)
 endif
else ! no correction
 getSpecAngFinal=getSpecAng(ip,s,r,d,e,g,ang)
endif

if(g.ge.100.0d0) getSpecAngFinal=getSpecAngFinal*BHfactor(ip,e,ang) ! Consider black hole

end

function getSpecAng(ip,s,r,d,e,g,ang)
! ip: particle ID
! s: W index
! r: cut-off rigidity in GV
! d: atmospheric depth in g/cm2
! e: energy in MeV/n
! g: local geometry effect
! ang: cos(theta)
parameter(npart=6) ! number of particle type (1:neutron, 2:proton, 3:heavy ion, 4:muon, 5:electron&positron, 6:photon
parameter(nsur=18)  ! number of surface, upto 52 km
parameter(ncor=7)   ! number of cut-off ridigity upto 20 GV
parameter(maxfit=8) ! number of parameters to express angular distribution
parameter(mpeach=10) ! number of parameters to express energy dependence of each parameter
parameter(ifour=4)  ! 4

implicit real*8 (a-h, o-z)

real*8, save:: ParaEdep(mpeach,maxfit,ncor,nsur,npart) ! parameter for expressing energy-differential energy dependence
real*8, save:: ParaEint(maxfit,ncor,nsur,npart) ! parameter for expressing energy-integrated energy dependence
real*8, save:: depth(nsur) ! depth (g/cm2) for parameters
real*8, save:: cor(ncor)   ! Rc (GV) for parameters
real*8, save:: ParaAdep(maxfit,ifour) ! parameter for expressing angular distribution, 1-4 is for each depth & Rc condition
real*8, save:: ratio1,ratio2 ! ratio must be saved
real*8, save:: emin(npart),emax(npart) ! maximum energy 
character chatmp1*1
character pname(npart)*6
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

data pname/'neutro','proton','he---4','muon--','elepos','photon'/
data emin/   1.0d-7,   1.0d0,   1.0d0,   1.0d1,  1.0d-1,  1.0d-2/      ! minimum energy
data emax/    1.0d4,   1.0d4,   1.0d4,   1.0d5,   1.0d4,   1.0d4/      ! maximum energy
data phimin/1.0e-3/ ! minimum value of phi

data ifirst/0/
data ipold/0/
data sold/0.0/
data rold/0.0/
data dold/0.0/
data eold/0.0/
data gold/0.0/

getSpecAng=1.0

! Read parameters
if(ifirst.eq.0) then ! first time call this routine
 ifirst=1
 do ip1=1,npart
  open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/angle/'//pname(ip1)//'.out',status='old')
  read(28,'(a1)') chatmp1
  do i=1,maxfit
   do is=1,nsur
    do ic=1,ncor
     read(28,*) itmp,depth(is),cor(ic),(ParaEdep(ii,i,ic,is,ip1),ii=1,mpeach)
    enddo
   enddo
  enddo
  close(28)
  open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/angle/'//pname(ip1)//'-Eint.out',status='old')
  read(28,'(a1)') chatmp1
  do is=1,nsur
   do ic=1,ncor
    read(28,*) tmp,tmp,(ParaEint(i,ic,is,ip1),i=1,maxfit)
   enddo
  enddo
  close(28)
 enddo
endif

! check previous condition
if(ip.eq.ipold.and.s.eq.sold.and.r.eq.rold.and.d.eq.dold.and.e.eq.eold.and.g.eq.gold) goto 10 ! same condition, need not determine parameter again
ipold=ip
sold=s
rold=r
dold=d
eold=e
gold=g
! Find depth
do is=1,nsur
 if(d.lt.depth(is)) exit
enddo
if(is.eq.1) then
 ratio1=0.0
 is=2
elseif(is.eq.nsur+1) then
 ratio1=1.0
 is=nsur
else
 ratio1=(d-depth(is-1))/(depth(is)-depth(is-1))
endif
! Find COR
do ic=1,ncor
 if(r.lt.cor(ic)) exit
enddo
if(ic.eq.1) then
 ratio2=0.0
 ic=2
elseif(ic.eq.ncor+1) then
 ratio2=1.0
 ic=ncor
else
 ratio2=(r-cor(ic-1))/(cor(ic)-cor(ic-1))
endif

idx=0
do is1=is-1,is
 do ic1=ic-1,ic
  idx=idx+1
  do i=1,maxfit
   if(e.eq.0.0) then ! energy integrated
    ParaAdep(i,idx)=ParaEint(i,ic1,is1,ip)
   else
    ene=max(min(e,emax(ip)),emin(ip))
    ParaAdep(i,idx)=getParaAdep(i,ene,ParaEdep(1,i,ic1,is1,ip))
    if(ip.eq.1.and.g.ge.0.0d0.and.g.le.1.0d0) ParaAdep(i,idx)=ParaAdep(i,idx)+getGneut(e,i) ! Ground level neutron correction, e instead of ene is used because data are down to 1e-8
    if(i.eq.2.and.ParaAdep(1,idx)+ParaAdep(2,idx).lt.phimin) ParaAdep(2,idx)=-ParaAdep(1,idx)+min(ParaAdep(1,idx),phimin) ! avoid negative value at cos(theta)=-1.0 
    if(i.eq.5.and.ParaAdep(4,idx)+ParaAdep(5,idx).lt.phimin) ParaAdep(5,idx)=-ParaAdep(4,idx)+min(ParaAdep(4,idx),phimin) ! avoid negative value at cos(theta)=1.0 
   endif
  enddo
 enddo
enddo

do idx=1,ifour
 call adjustParaAdep(ParaAdep(1,idx)) ! integration value is adjusted to 1
enddo  

10 continue ! if all conditions are same, jump to here

B1=funcAng(ang,ParaAdep(1,1))
B2=funcAng(ang,ParaAdep(1,2))
B3=funcAng(ang,ParaAdep(1,3))
B4=funcAng(ang,ParaAdep(1,4))

C1=B1+(B2-B1)*ratio2 ! Rc interpolation
C2=B3+(B4-B3)*ratio2 ! Rc interpolation

getSpecAng=C1+(C2-C1)*ratio1 ! Depth interpolation

return
end

function getParaAdep(i,x,A)
! i: parameter index
! x: energy
! A: ParaEdep
parameter(mpeach=10) ! number of parameters to express energy dependence of each parameter
implicit real*8 (a-h, o-z)
dimension A(mpeach)

a4=max(0.01,a(4))
a7=max(0.01,a(7))
a10=max(0.01,a(10))

getParaAdep=a(1)+a(2)/(1+exp((a(3)-log10(x))/a4))+a(5)/(1+exp((a(6)-log10(x))/a7))+a(8)/(1+exp((a(9)-log10(x))/a10))

if(i.eq.3.or.i.eq.6) getParaAdep=10.0**getParaAdep       ! a(3) & a(6) are determined by log10 fit
if(i.eq.7) getParaAdep=max(0.1d0,min(1.0d0,getParaAdep)) ! a(7) should be 0.1 < a7 1.0
if(i.eq.1.or.i.eq.4.or.i.eq.8) getParaAdep=max(0.0d0,getParaAdep)  ! a(1), a(4) & a(8) should be positive

return

end

subroutine AdjustParaAdep(A)
parameter (maxfit=8)
implicit real*8 (a-h, o-z)
dimension A(maxfit)
data one/1.0d0/
data atlow/-0.05/  ! lower threshold angle for 90 degree interpolation
data athig/ 0.05/  ! higher threshold angle for 90 degree interpolation

sum1=max(0.0,getint(abs(atlow),one,a(1),a(2),a(3)))                        ! -1.0 to atlow
sum2=max(0.0,(funcAng(atlow,a(1))+funcAng(athig,a(1)))*0.5*(athig-atlow)) ! atlow to athig
sum3=max(0.0,getint(athig,a(7),a(4),a(5),a(6)))                           ! athig to a(7)
if(a(7).lt.one) then
 sum4=max(0.0,(funcAng(a(7),a(1))+funcAng(one,a(1)))*0.5*(one-a(7)))           ! a(7) to 1.0
else
 sum4=0.0
endif
sum=max(1.0e-20,(sum1+sum2+sum3+sum4)*2.0*acos(-1.0))

A(1)=A(1)/sum
A(2)=A(2)/sum
A(4)=A(4)/sum
A(5)=A(5)/sum
A(8)=A(8)/sum

return
end

function getint(x0,x1,a1,a2,a3)
implicit real*8 (a-h,o-z)
getint=a1*(x1-x0)+a2*(x1**(a3+1)-x0**(a3+1))/(a3+1)
return
end

function funcAng(x,a)
parameter (maxfit=8)
implicit real*8 (a-h, o-z)
dimension A(maxfit)
data atlow/-0.05/  ! lower threshold angle for 90 degree interpolation
data athig/ 0.05/  ! higher threshold angle for 90 degree interpolation

if(x.le.atlow) then ! backward
 funcAng=a(1)+a(2)*abs(x)**a(3)
elseif(x.lt.athig) then ! intermediate
 tmp1=a(1)+a(2)*abs(atlow)**a(3)
 tmp2=a(4)+a(5)*athig**a(6)
 funcAng=tmp1+(tmp2-tmp1)*(x-atlow)/(athig-atlow) ! simply interpolate
elseif(x.le.a(7)) then ! forward
 funcAng=a(4)+a(5)*x**a(6)
else ! after peak
 tmp1=a(4)+a(5)*a(7)**a(6)
 tmp2=a(8)
 funcAng=tmp1+(tmp2-tmp1)*(x-a(7))/(1.0d0-a(7))
endif

return
end

function getGmuon(emid,ang)
implicit real*8 (a-h, o-z)

elog=max(4.062,log10(emid))
a2=9.9873006E-01+2.9141114E+00/(1.0+exp((5.8030900E+00-elog)/2.4585039E-01)) 
a3=max(0.01,2.4226042e1-1.5142933e1*elog+3.2012346*elog**2-2.2325286e-1*elog**3)
elog=min(6.0d0,elog) 
a4=1.4970401e1-5.3110524*elog+4.7458357e-1*elog**2 

if(ang.lt.0.0d0) then
 getGmuon=0.0d0
else
 x=1.0d0/max(0.001,ang)
 Scal=a2*sqrt(1-EXP(-((1.0/a2)**2)))
 getGmuon=(a2*sqrt(1-EXP(-((x/a2)**2)))-(1-EXP(-a4*(x-1)**a3)))/Scal
endif


return

end

function getGneut(emid,ID)
parameter(maxfit=8)   ! number of fitting parameter for angular dependence
parameter(maxEfit=3)  ! number of fitting parameter for energy dependence
implicit real*8 (a-h, o-z)
real*8, save:: Gneut(maxEfit,maxfit)
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

data ifirst/0/

if(ifirst.eq.0) then ! first time call this routine
! Read Ground Correction Factor
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/angle/NeutronGround.out',status='old')
 do i=1,maxfit
  read(28,*) (Gneut(ie,i),ie=1,maxEfit)
 enddo
 close(28)
 ifirst=1
endif 

elog=max(-8.0,log10(emid)) ! parameters are effective only above 0.01 eV
getGneut=Gneut(1,id)/(1+exp((elog-Gneut(2,id))/Gneut(3,id)))

return

end

function BHfactor(ip,e,ang) ! Black hole factor
implicit real*8 (a-h, o-z)
parameter(npart=6) ! only neutron, elepos, and photon
parameter(nBHpara=3) ! number of parameter
parameter(nBHeach=7) ! number of parameter to represent each parameter
real*8, save:: BHpara(nBHeach,nBHpara,npart)
dimension dimtmp(nBHpara) ! temporary used dimension
character chatmp1*1,chatmp40*40
character pname(npart)*6
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

data pname/'neutro','proton','he---4','muon--','elepos','photon'/
data ifirst/0/


if(ifirst.eq.0) then ! first time call this routine
! Read Black hole factor parameter (only for neutron, elepos, photon)
 do ip2=1,npart
  open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/angle/BkH-'//pname(ip2)//'.inp',status='old')
  read(28,'(a1)') chatmp1
  do i=1,nBHpara
   read(28,*) (BHpara(ii,i,ip2),ii=1,nBHeach)
  enddo
  close(28)
 enddo
 ifirst=1
endif 

if(ang.lt.0.0) then
 BHfactor=0.0 ! backward is always 0
else
 do i=1,nBHpara
  dimtmp(i)=doublesig(e,BHpara(1,i,ip))
 enddo
 BHfactor=dimtmp(1)+(dimtmp(2)-dimtmp(1))*ang**dimtmp(3)
endif

end

function doublesig(e,a) ! get double sigmoid
implicit real*8 (a-h,o-z)
parameter(nBHeach=7) ! number of parameter to represent each parameter
dimension a(nBHeach)
doublesig=a(1)+a(2)/(1+exp((a(3)-log10(e))/a(4)))+a(5)/(1+exp((a(6)-log10(e))/a(7)))
return
end

	function getHP(iy0,im0,id0,ic) ! get FFP from FFP tables
!   ic=1: obtain FFP from neutron monitor data
!	ic=2: obtain FFP from Wolf number
!	ic=3: suspected ground level event
!	ic=4: Too long time ago or future
!	ic=5: no such date
parameter(nmonth=12)
parameter(nday=31)
parameter(iymax=2100) ! maximum year
parameter(iymin=1614) ! earliest year
	
implicit real*8 (a-h, o-z)

real, save:: FFP(iymin:iymax,nmonth,nday)
real, save:: FFPuso(iymin:iymax) ! 0 & nmonth+1 data are used for interpolation
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

data ifirst/0/
integer*4, save:: iystart,iyend,iysUs,iyeUs

if(ifirst.eq.0) then ! first time call this subroutine
 do iy=iymin,iymax  ! Initialized
  FFPuso(iy)=0.0
  do im=1,nmonth
   do id=1,nday
    FFP(iy,im,id)=0.0
   enddo
  enddo
 enddo
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/FFPtable.day',status='old')
 read(28,*) iystart,iyend
 do im=1,nmonth
  do id=1,nday
   read(28,*) itmp1,itmp2,(FFP(iy,im,id),iy=iystart,iyend)
  enddo
 enddo
 close(28)
 open(28,file=chfn(1)(1:ilfn(1))//'/data/cosmicray/FFPtable.uso',status='old')
 read(28,*) iysUs,iyeUs
 do iy=iysUs,iyeUs
  read(28,*) itmp,FFPuso(iy)
 enddo
 close(28)
 ifirst=1
endif

! ****** Year, month, day Check **************
if(iy0.lt.iymin.or.iy0.gt.iymax) then  ! out of range
 ic=4
 getHP=0.0
 return
endif
if(im0.lt.1.or.im0.gt.nmonth.or.id0.lt.1.or.id0.gt.nday) then
 ic=5
 getHP=0.0
 return
endif

! ******Determine FFP from Neutron Monitor *************
if(FFP(iy0,im0,id0).gt.-99.0) then ! data exist
 if(iy0.ge.iystart.and.FFP(iy0,im0,id0).eq.-1000.0) then ! no such date
  ic=5
  getHP=0.0
 else  ! neutron monitor data exist
  getHP=FFP(iy0,im0,id0)
  if(getHP.gt.1000.0) then
   ic=3 ! GLE occurred	  
   write(*,*) 'GLE was occured on the day. The results might be wrong'
   getHP=getHP-10000.0
  else
   ic=1 ! normal data
  endif
 endif
 return
endif

!  ***** Determine FFP from Usoskin's data  *****
if(FFPuso(iy0).ne.0.0) then
 getHP=FFPuso(iy0)
 ic=2  ! determine FFP from Usoskin's data
 return
endif

!  **** No data **************************
ic=4
getHP=0.0
return

end
