!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas_source
      module moddas_source

      double precision, allocatable :: egmin(:)
      double precision, allocatable :: egmax(:)
      double precision, allocatable :: fegrp(:)
      double precision, allocatable :: rfe  (:)
      double precision, allocatable :: prw  (:)
      double precision, allocatable :: pwt  (:)

      double precision, allocatable :: agmin(:)
      double precision, allocatable :: agmax(:)
      double precision, allocatable :: fagrp(:)
      double precision, allocatable :: rfa  (:)
      double precision, allocatable :: paw  (:)
      double precision, allocatable :: pat  (:)

      double precision, allocatable :: tgmin(:)
      double precision, allocatable :: tgmax(:)
      double precision, allocatable :: ftgrp(:)
      double precision, allocatable :: rft  (:)
      double precision, allocatable :: ptw  (:)
      double precision, allocatable :: ptt  (:)

      double precision, allocatable :: slmin(:)
      double precision, allocatable :: slmax(:)
      double precision, allocatable :: flgrp(:)
      double precision, allocatable :: rgl  (:)
      double precision, allocatable :: rgw  (:)

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_source_initialize()

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      common /isorse/ ngrp(isrc), ngei(isrc), ngea(isrc), ngfe(isrc),
     &                ngft(isrc), ngll(isrc), ngpi(isrc), ngpw(isrc)
      common /isoraa/ narp(isrc), naei(isrc), naea(isrc), nafe(isrc),
     &                naft(isrc), nall(isrc), napi(isrc), napw(isrc)
      common /isortt/ ntrp(isrc), ntei(isrc), ntea(isrc), ntfe(isrc),
     &                ntft(isrc), ntll(isrc), ntpi(isrc), ntpw(isrc)
      common /isorsl/ nglp(isrc), ngli(isrc), ngla(isrc), ngfl(isrc),
     &                nglc(isrc), nglw(isrc)

      ngei(:) = 0
      ngea(:) = 0
      ngfe(:) = 0
      ngft(:) = 0
      ngpi(:) = 0
      ngpw(:) = 0

      naei(:) = 0
      naea(:) = 0
      nafe(:) = 0
      naft(:) = 0
      napi(:) = 0
      napw(:) = 0

      ntei(:) = 0
      ntea(:) = 0
      ntfe(:) = 0
      ntft(:) = 0
      ntpi(:) = 0
      ntpw(:) = 0

      ngli(:) = 0
      ngla(:) = 0
      ngfl(:) = 0
      nglc(:) = 0
      nglw(:) = 0
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module
