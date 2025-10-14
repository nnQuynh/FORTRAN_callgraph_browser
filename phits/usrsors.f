************************************************************************
*                                                                      *
      subroutine usrsors(x,y,z,u,v,w,e,wt,time,name,kf,nc1,nc2,nc3,
     &                   sx,sy,sz)
*                                                                      *
*        sample subroutine for user defined source.                    *
*                                                                      *
*        variables :                                                   *
*                                                                      *
*           x, y, z : position of the source.                          *
*           u, v, w : unit vector of the particle direction.           *
*           e       : kinetic energy of particle (MeV).                *
*           wt      : weight of particle.                              *
*           time    : initial time of particle. (ns)                   *
*           name    : usually = 1, for Coulmb spread.                  *
*           kf      : kf code of the particle.                         *
*           nc1     : initial value of counter 1                       *
*           nc2     : initial value of counter 2                       *
*           nc3     : initial value of counter 3                       *
*           sx,sy,sz : spin components                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        kf code table                                                 *
*                                                                      *
*           kf-code: ityp :  description                               *
*                                                                      *
*             2212 :   1  :  proton                                    *
*             2112 :   2  :  neutron                                   *
*              211 :   3  :  pion (+)                                  *
*              111 :   4  :  pion (0)                                  *
*             -211 :   5  :  pion (-)                                  *
*              -13 :   6  :  muon (+)                                  *
*               13 :   7  :  muon (-)                                  *
*              321 :   8  :  kaon (+)                                  *
*              311 :   9  :  kaon (0)                                  *
*             -321 :  10  :  kaon (-)                                  *
*                                                                      *
*           kf-code of the other transport particles                   *
*                                                                      *
*               12 :         nu_e                                      *
*               14 :         nu_mu                                     *
*              221 :         eta                                       *
*              331 :         eta'                                      *
*             -311 :         k0bar                                     *
*            -2112 :         nbar                                      *
*            -2212 :         pbar                                      *
*             3122 :         Lambda0                                   *
*             3222 :         Sigma+                                    *
*             3212 :         Sigma0                                    *
*             3112 :         Sigma-                                    *
*             3322 :         Xi0                                       *
*             3312 :         Xi-                                       *
*             3334 :         Omega-                                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        available function for random number                          *
*                                                                      *
*           unirn(dummy) : uniform random number from 0 to 1           *
*           gaurn(dummy) : gaussian random number                      *
*                          for exp( - x**2 / 2 / sig**2 ) : sig = 1.0  *
*                                                                      *
*----------------------------------------------------------------------*
*        use of parameter specified in your input                      *
*     You can use parameters c1 - c100 as cval(1) - cval(100)          *
*     If the parameter is defined more than once in your input file,   *
*     the value defined at last is used in this subroutine             *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      common /rval1/ cval(mxcval), aval(mxcval)

      data ifirst / 0 /
      save ifirst
!$OMP THREADPRIVATE(ifirst)
      character filenm*50
      logical   exex

*-----------------------------------------------------------------------
*        example of initialization
*-----------------------------------------------------------------------

         if( ifirst .eq. 0 ) then

c              filenm = 'input.dat'

c              inquire( file = filenm, exist = exex )
c              if( exex .eqv. .false. ) then
c                 write(*,*) 'file does not exist =>  ', filenm
c                 call parastop( 887 )
c                 end if

c              open(71, file = filenm, status = 'old' )



c              close(71)

               ifirst = 1

         end if

*-----------------------------------------------------------------------
*        example for 3 GeV proton with z-direction
*-----------------------------------------------------------------------

            x = 0.0
            y = 0.0
            z = 0.0

            u = 0.0
            v = 0.0
            w = 1.0

            e = 3000.0

           wt = 1.0
         time = 0.0
         name = 1

           kf = 2212

          nc1 = 0
          nc2 = 0
          nc3 = 0

           sx = 0.d0
           sy = 0.d0
           sz = 0.d0

*-----------------------------------------------------------------------

      return
      end

