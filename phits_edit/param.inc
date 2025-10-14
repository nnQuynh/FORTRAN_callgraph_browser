************************************************************************
*                                                                      *
*     'param.inc'                                                      *
*                                                                      *
************************************************************************
!     when you increase mdas over 250000000, you have to define integer*8 instead of integer in the next line
      integer mdas,mcmx,mci,mmdas,mmmax,nbnds,mct,multmax,pnlmax
      integer kvlmax, kvmmax, itlmax, inevt, isrc, nbchmax, maxchapara
      integer mxprodxs
      real*8  das
      parameter ( mdas   = 80000000 )
      parameter ( kvlmax =     3000 )
      parameter ( kvmmax =  1000000 )
      parameter ( itlmax =      200 )
      parameter ( inevt  =       71 )
      parameter ( isrc   =     3000 )
      parameter ( nbchmax=    10000 )
      parameter ( multmax=      500 )
      parameter ( pnlmax =     3000 )
      parameter ( maxchapara = 1000 )
      parameter ( mxprodxs =  20000 )

*----------------------------------------------------------------------*
*                                                                      *
*       pnlmax : maximum number of material in photonuclear library    *
*                                                                      *
*----------------------------------------------------------------------*

      common /mdasa/ das( mdas )
      common /mdasb/ mmmax

*----------------------------------------------------------------------*
*                                                                      *
*       mdas   : total memory * 8 = byte                               *
*       mmmax  : maximum number of total array                         *
*                                                                      *
*       kvlmax : maximum number of regions, cell and material          *
*       kvmmax : maximum number of id for regions, cel and material    *
*                                                                      *
*       itlmax : number of maximum tally entry                         *
*       inevt  : number of collision type for summary                  *
*       isrc   : number of multi-source                                *
*       nbchmax: maximum number of batch assigned to parallel MPI node *
*       pnlmax : maximum number of material in photonuclear library    *
*       mxprodxs : maximum number of produced nuclides in prodxs       *
*                                                                      *
*----------------------------------------------------------------------*
