************************************************************************
*                                                                      *
      subroutine setrnd(io,jo,ierr)
c ----------------------------------------------------------------------
c     controls for the pseudo-random number sequence.
c     When nrandgen=0, Linear Congruential Generator is used.
c     When nrandgen=1, xorshift is used.
c ----------------------------------------------------------------------
      implicit none
      integer :: io,jo,ierr
      integer :: nrandgen
      common /randn/ nrandgen
c ----------------------------------------------------------------------
c
      ierr = 0

      if ( nrandgen .eq. 0 ) then
         call setrnd_LCG(io,jo,ierr)
      else
         call setrnd_xorshift(io,jo,ierr)
      end if

      return
      end


************************************************************************
*                                                                      *
      subroutine setrnd_LCG(io,jo,ierr)
c ----------------------------------------------------------------------
c        randkk: rseed(>0) specified in [parameters]
c ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      include 'err.inc'
      parameter (fb=13008944d0,fs=170125d0,gb=1136868d0,gs=6328637d0,
     1 p=2d0**24,q=2d0**(-24),rm=5d0**19)
c
      integer*8 :: iranji64
      common /randm4/ rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,
     &                rnrtc,nstrid,inif, iranji64
      integer*8 :: iransb64
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)
      common /iradkk/ randkk,irskip
      dimension dbcn(14)
      data dbcn/14*0.d0/
c ----------------------------------------------------------------------
c
      nstrid=152917
      rnfb=fb
      rnfs=fs
      rngb=gb
      rngs=gs
      rnmult=rm
      inif=0

      if(randkk.gt.0.) dbcn(1) = randkk


c        set new random number multiplier, rnmult, if required.
c        rngb and rngs are the upper and lower 24 bits of rnmult.

      if(dbcn(13)+dbcn(14).le.0.)go to 30
      if(dbcn(14).gt.0.)rnmult=dbcn(14)
      if(aint((rnmult+.5)*.5).eq.aint((rnmult+1.5)*.5)) then
         write(io,*) 'ERROR: random number multiplier is even.'
         ErrCha = ''
         MsgID = 'L:69/R:setrnd_LCG/F:randmc.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,*) 'ERROR: random number multiplier is even.'
         ierr = 1
         return
      end if
      rngb=aint(rnmult*q)
      rngs=rnmult-rngb*p
      if(rngb+rngs.ge.p) then
         write(io,*) 'ERROR: random number multiplier rejected.'
         ErrCha = ''
         MsgID = 'L:80/R:setrnd_LCG/F:randmc.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,*) 'ERROR: random number multiplier rejected.'
         ierr = 1
         return
      end if
c
c        get rnfb (upper 24 bits) and rnfs (lower 24 bits) of
c        rnmult**nstrid which is used in advijk to advance the random
c        number by nstrid random numbers for each history.

      if(dbcn(13).gt.0.)nstrid=dbcn(13)
      rnfb=rngb
      rnfs=rngs
      ii=nstrid-1
      do 10 i=1,ii
      a=rngs*rnfs
      b=(rngb*rnfs-aint(rngb*rnfs*q)*p)
     1 +(rngs*rnfb-aint(rngs*rnfb*q)*p)+aint(a*q)
      rnfs=a-aint(a*q)*p
   10 rnfb=b-aint(b*q)*p

   20 format(/23h random number stride =,i19/
     1 27h random number multiplier =,f16.0,tl1,1h )
c
c        set the first random number, rijk, composed of
c        ranj (top 24 bits) and rani (bottom 24 bits).
   30 rijk=rnmult
      if(dbcn(1).gt.0.)rijk=dbcn(1)
      if(dbcn(1)+dbcn(8).ne.0.)inif=1

      rani=aint(rijk*q)
      ranj=rijk-rani*p

      if( inif .eq. 0 ) call advijk ! when rseed was not specified

      do 40 i=1,int(dbcn(8))
   40 call advijk

      if(dbcn(1)+dbcn(8).ne.0.) then
      end if

   50 format(' starting random number =',2x,f16.0,tl1,1h )
      dbcn(1)=0.
      dbcn(8)=0.

      return
      end


!***********************************************************************
      subroutine setrnd_xorshift(io,jo,ierr)
c ----------------------------------------------------------------------
c     when nrandgen=1; xorshift is used.
c
c        randkk: rseed(>0) specified in [parameters]
c ----------------------------------------------------------------------
      implicit none

      include 'err.inc'
      integer :: io,jo,ierr
      real*8 :: rijk,rans,ranb,randkk
      integer :: irskip
      integer*8 :: iransb64
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)
      common /iradkk/ randkk,irskip
      real*8 :: rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,rnrtc
      integer :: nstrid,inif
      integer*8 :: iranji64
      common /randm4/ rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,
     &                rnrtc,nstrid,inif, iranji64
      real*8 :: bitrseed
      common /randnbit/ bitrseed
c ----------------------------------------------------------------------

      rijk = bitrseed
      if ( rijk .eq. 0d0 ) then
 1000    format('** Error : random seed rijk is 0.')
         write(io,1000)
         ErrCha = ''
         MsgID = 'L:161/R:setrnd_xorshift/F:randmc.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,1000)
         ierr = ierr + 1
      end if

      iranji64 = transfer(rijk,iranji64)

c ----------------------------------------------------------------------
      return
      end subroutine
!***********************************************************************


************************************************************************
*                                                                      *
      subroutine advijk
c ----------------------------------------------------------------------
c     advance source random number, rijk, for the next history.
c ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /randn/ nrandgen
      integer*8 :: iransb64
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)

c ----------------------------------------------------------------------

      if ( nrandgen .eq. 0 ) then
         call advijk_LCG
      else
         call advijk_xorshift64
      end if

      return
      end


************************************************************************
*                                                                      *
      function rang()
c ----------------------------------------------------------------------
c        return the next pseudo-random number.
c ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /randn/ nrandgen
      integer*8 :: iransb64
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)
      common /ncall0/ ncall               !FURUTA
c ----------------------------------------------------------------------

      if ( nrandgen .eq. 0 ) then
         rang = rang_LCG()
      else
         rang = rang_xorshift64()
      end if

!$OMP ATOMIC
      ncall=ncall+1
      return
      end


************************************************************************
*                                                                      *
      subroutine advijk_LCG
c ----------------------------------------------------------------------
c     advijk_LCG and rang_LCG (below) are pseudo-random number generator
c     based on algorithm of Linear Congruential Generator.
c     Reference: R. Picard and T. Booth [LA-UR-08-06204]
c     period length is 2**46
c
c     ranj is the lower 24 bits of rijk, rani is the upper 24 bits.
c ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (p=2d0**24,q=2d0**(-24))
      integer*8 :: iranji64
      common /randm4/ rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,
     &                rnrtc,nstrid,inif, iranji64
      integer*8 :: iransb64
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)
c ----------------------------------------------------------------------
      a=rnfs*ranj
      b=(rnfb*ranj-aint(rnfb*ranj*q)*p)
     1 +(rnfs*rani-aint(rnfs*rani*q)*p)+aint(a*q)
      ranj=a-aint(a*q)*p
      rani=b-aint(b*q)*p
      rijk=rani*p+ranj
      return
      end


************************************************************************
*                                                                      *
      function rang_LCG()
c ----------------------------------------------------------------------
c        return the next pseudo-random number.
c ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (p=2d0**24,q=2d0**(-24),r=2d0**(-48))
c
      integer*8 :: iranji64
      common /randm4/ rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,
     &                rnrtc,nstrid,inif, iranji64
      integer*8 :: iransb64
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)
c ----------------------------------------------------------------------
c
c        rang()=mod(2**48*rang()*rnmult,2**48)
c        split rang() and rnmult into upper and lower 24-bit halves,
c        ranb,rans,rngb,rngs, respectively, to achieve 96-bit precision.
c        this expression for b is invalid unless rngb+rngs < 2**24
c        or unless the more elaborate expression of advijk is used.
c
      a=rngs*rans
      b=rngb*rans+rngs*ranb+aint(a*q)
      rans=a-aint(a*q)*p
      ranb=b-aint(b*q)*p
      rang_LCG=(ranb*p+rans)*r
      return
      end


!***********************************************************************
      subroutine advijk_xorshift64
c ----------------------------------------------------------------------
c     advijk_xorshift64 and rang_xorshift64 (below) are pseudo-random
c     number generator based on algorithm of xorshift.
c     These period lengths of xorshift64 and xorshift32 are 2**64-1.
c     Binary matrixes T=(I+L^a)(I+R^b)(I+L^c) with two sets
c     [a,b,c] =[13,7,17] and [29,27,37] are used for advijk and rang, respectively.
c
c     References:
c       G. Marsaglia, 2003. Xorshift RNGs. Journal of Statistical Software. 8, 14, 1-6,
c     and webpage for introduction and verification of this paper
c     [https://admiswalker.blogspot.com/2018/06/marsaglia-g-2003-xorshift-rngs.html
c     (in Japanese, 2020.1.29 accessed)]
c ----------------------------------------------------------------------
      implicit none
      real*8 :: rijk,rans,ranb
      integer*8 :: iransb64
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)
      integer*8 :: it
      real*8 :: rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,rnrtc
      integer :: nstrid,inif
      integer*8 :: iranji64
      common /randm4/ rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,
     &                rnrtc,nstrid,inif, iranji64

c ----------------------------------------------------------------------

      it = xor(iranji64,ishft(iranji64,13))
      it = xor(it,ishft(it,-7))
      iranji64 = xor(it,ishft(it,17))

      rijk = transfer(iranji64,rijk)

      return

c ----------------------------------------------------------------------
      end subroutine
!***********************************************************************


!***********************************************************************
      function rang_xorshift64()
      implicit none
      real*8 :: rang_xorshift64, r01
      real*8 :: rijk,rans,ranb
      integer*8 :: iransb64,it
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)
      integer*4 :: i32
      real*8 :: bit32, period
      parameter (bit32=2d0**32, period=2d0**32-1)

c ----------------------------------------------------------------------

      it = xor(iransb64,ishft(iransb64,29))
      it = xor(it,ishft(it,-27))
      iransb64 = xor(it,ishft(it,37))

      i32 = ibits(iransb64,32,32)
      if ( i32.eq.0 ) i32 = ibits(iransb64,0,32)
      r01 = dble(i32)
      if ( i32.lt.0 ) r01 = r01+bit32
      rang_xorshift64 = r01/period

      return

c ----------------------------------------------------------------------
      end function
!***********************************************************************
