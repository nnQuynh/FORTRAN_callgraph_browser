!----------------------------cpcom.f------------------------------------
!  Version: 051219-1435
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12

! ----------------------------------------
! PEGS5 common file for CPCOM (Compton Profile Commons)
! ----------------------------------------

      integer nshell, icprof, mxshel, mxraw
      double precision cprof, avcprf, cprofi, qcap, qcap10, scprof,
     &                 scpsum, scproi, capin, capio, capils, elecni,
     &                 elecnj, elecno, cpimev
      COMMON/CPCOM/CPROF(31,102),AVCPRF(31),CPROFI(301),QCAP(31), 
     &             QCAP10(301),SCPROF(31,MXMXRAW),SCPSUM(31,MXMXRAW),
     &             SCPROI(301,MXMXRAW),CAPIN(MXMXRAW),CAPIO(MXMXRAW),
     &             CAPILS(MXMXRAW),
     &             ELECNI(MXMXRAW),ELECNJ(MXMXRAW),ELECNO(MXMXRAW),
     &             CPIMEV,
     &             NSHELL(MXMXRAW),ICPROF,MXSHEL,MXRAW

!--------------------------last line------------------------------------
