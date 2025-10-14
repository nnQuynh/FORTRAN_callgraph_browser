************************************************************************
*                                                                      *
*     subroutine egs2phits(np,iq,e,u,v,w,wt)
      subroutine egs2phits(ncpls)
*                                                                      *
*       move all EGS5 particle information to PHITS
*       last modified by H.Iwase on 2012/3/6
*                                                                      *
*       ncpls = 1 : call from collision part, then nclsts
*       ncpls = 2 : call from decay part, then nclst
*       ncpls = 3 : call from annihilation part, then nclst
*       cKN 2014/09/09, T.Sato 2015/2/26
*                                                                      *
************************************************************************
      implicit none
!      save
      include 'include/egs5_h.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_useful.f'


      integer k
      integer nnn,nomp,nclsts,iclusts,jclusts
      include 'param00.inc'
      real*8 qclusts
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE( /clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE( /clustw/)

      integer nclst,iclust,jclust
      real*8 qclust
      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE( /clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE( /clustg/)

      integer kf,itype,id,charge
      real*8 rms,rmg,pr,etot,ekine,pxl,pyl,pzl

      integer ncpls

      integer numpal,i
      real*8 rumpal
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)

! T.Sato 2015/2/26 annihilation positron energy = emin(13)
      common /eparm/  esmax, esmin, emin(20)
      real*8 esmax,esmin,emin

      nclst  = 0
      nclsts = 0
      do i=0,20
       numpal(i)=0
       rumpal(i)=0
      enddo

       do k = np,1,-1

         if     ( iq(k) .eq. -1 ) then ! e-
            kf = 11             ! 11:e-, -11:e+, 22:photon
            itype = 7           ! 7:e-,e+, 4:photon
            id = 12             ! 12:e-, 13:e+, 14:photon
            charge = -1
            rms = RM
            rmg = rms / 1000d0
            e(k) = e(k) - RM
            pr  = sqrt( e(k) * ( e(k) + 2.0d0 * rms ) ) / 1000d0
            etot = sqrt( pr**2 + rmg**2 )
            ekine  = e(k)


         elseif ( iq(k) .eq.  1 ) then ! e+

            kf = -11            ! 11:e-, -11:e+, 22:photon
            itype = 7           ! 7:e-,e+, 4:photon
            id = 13             ! 12:e-, 13:e+, 14:photon
            charge = 1
            rms = RM
            rmg = rms / 1000d0
            e(k) = e(k) - RM
            pr  = sqrt( e(k) * ( e(k) + 2.0d0 * rms ) ) / 1000d0
            etot = sqrt( pr**2 + rmg**2 )
            ekine  = e(k)

         elseif ( iq(k) .eq.  0 ) then ! photon

            kf = 22             ! 11:e-, -11:e+, 22:photon
            itype = 4           ! 7:e-,e+, 4:photon
            id = 14             ! 12:e-, 13:e+, 14:photon
            charge = 0
            rms = RM
            rmg = 0d0
            e(k) = e(k)
            etot   = e(k)/1000d0
            ekine  = e(k)
            pr  = etot

         endif

         if( ncpls .eq. 1 ) then

            nclsts = nclsts + 1

            iclusts(nclsts)   = itype

            jclusts(0,nclsts) = 0
            jclusts(1,nclsts) = 0
            jclusts(2,nclsts) = 0
            jclusts(3,nclsts) = id
            jclusts(4,nclsts) = 0
cABE 2022/02/24, turn on the annihilation flag for positron with the energies below cutoff
            if( kf .eq. -11 .and. ekine .lt. emin(13) )
     &       jclusts(4,nclsts) = 2
            jclusts(5,nclsts) = charge
            jclusts(6,nclsts) = 0
            jclusts(7,nclsts) = kf
            jclusts(8,nclsts) = 0      ! isobar level

            pxl = pr * u(k)
            pyl = pr * v(k)
            pzl = pr * w(k)

            qclusts(0,nclsts)  = 0.0
            qclusts(1,nclsts)  = pxl
            qclusts(2,nclsts)  = pyl
            qclusts(3,nclsts)  = pzl
            qclusts(4,nclsts)  = etot              ! total E (GeV)
            qclusts(5,nclsts)  = rmg
            qclusts(6,nclsts)  = 0.0
            qclusts(7,nclsts)  = ekine             ! kinetic E
            qclusts(8,nclsts)  = 1.0d0

            qclusts(9,nclsts)  = 0.0
            qclusts(10,nclsts) = 0d0
            qclusts(11,nclsts) = 0d0
            qclusts(12,nclsts) = 0d0

        else

            nclst = nclst + 1

            iclust(nclst)   = itype

            jclust(0,nclst) = 0
            jclust(1,nclst) = 0
            jclust(2,nclst) = 0
            jclust(3,nclst) = id
            jclust(4,nclst) = 0
            jclust(5,nclst) = charge
            jclust(6,nclst) = 0
            jclust(7,nclst) = kf
            jclust(8,nclst) = 0      ! isobar level

            pxl = pr * u(k)
            pyl = pr * v(k)
            pzl = pr * w(k)

            qclust(0,nclst)  = 0.0
            qclust(1,nclst)  = pxl
            qclust(2,nclst)  = pyl
            qclust(3,nclst)  = pzl
            qclust(4,nclst)  = etot              ! total E (GeV)
            qclust(5,nclst)  = rmg
            qclust(6,nclst)  = 0.0
            qclust(7,nclst)  = ekine             ! kinetic E
            qclust(8,nclst)  = 1.0d0

            qclust(9,nclst)  = 0.0
            qclust(10,nclst) = 0d0
            qclust(11,nclst) = 0d0
            qclust(12,nclst) = 0d0

        end if

        enddo

        if( ncpls .eq. 3 ) then ! for annihilation, original positron information
            nclst = nclst + 1

            iclust(nclst)   = 7   ! 7:e-,e+, 4:photon

            jclust(0,nclst) = 0
            jclust(1,nclst) = 0
            jclust(2,nclst) = 0
            jclust(3,nclst) = 13  ! 12:e-, 13:e+, 14:photon
            jclust(4,nclst) = -1  ! dead particle
            jclust(5,nclst) = 1      ! charge
            jclust(6,nclst) = 0
            jclust(7,nclst) = -11    ! kf code
            jclust(8,nclst) = 0      ! isobar level

            qclust(0,nclst)  = 0.0
            qclust(1,nclst)  = 0.0
            qclust(2,nclst)  = 0.0
            qclust(3,nclst)  = 0.0
            qclust(4,nclst)  = 0.0              ! total E (GeV)
            qclust(5,nclst)  = RM*1.0d-3
            qclust(6,nclst)  = 0.0
            qclust(7,nclst)  = emin(13)         ! kinetic E = positron cut-off energy for EGS5
cABE 2022/02/24, e(np+1) is the positron energy saved in phits5annih
            if( e(np+1) .lt. emin(13) ) qclust(7,nclst)  = e(np+1)
            qclust(8,nclst)  = 1.0d0

            qclust(9,nclst)  = 0.0
            qclust(10,nclst) = 0d0
            qclust(11,nclst) = 0d0
            qclust(12,nclst) = 0d0
        endif

        if( ecapbind .ne. 0.0 ) then ! for EII, capbind energy is regarded as dead electron

            nclsts = nclsts + 1

            rms = RM
            rmg = rms / 1000d0
            pr = sqrt( ecapbind * ( ecapbind + 2.0d0 * rms ) ) / 1000d0
            etot = sqrt( pr**2 + rmg**2 )

            iclusts(nclsts)   = 7   ! 7:e-,e+, 4:photon

            jclusts(0,nclsts) = 0
            jclusts(1,nclsts) = 0
            jclusts(2,nclsts) = 0
            jclusts(3,nclsts) = 12  ! 12:e-, 13:e+, 14:photon
            jclusts(4,nclsts) = -11  ! dead particle
            jclusts(5,nclsts) = -1      ! charge
            jclusts(6,nclsts) = 0
            jclusts(7,nclsts) = 11    ! kf code
            jclusts(8,nclsts) = 0      ! isobar level

            qclusts(0,nclsts)  = 0.0
            qclusts(1,nclsts)  = 0.0
            qclusts(2,nclsts)  = 0.0
            qclusts(3,nclsts)  = 0.0
            qclusts(4,nclsts)  = etot              ! total E (GeV)
            qclusts(5,nclsts)  = rmg
            qclusts(6,nclsts)  = 0.0
            qclusts(7,nclsts)  = ecapbind             ! kinetic E
            qclusts(8,nclsts)  = 1.0d0

            qclusts(9,nclsts)  = 0.0
            qclusts(10,nclsts) = 0d0
            qclusts(11,nclsts) = 0d0
            qclusts(12,nclsts) = 0d0

            ecapbind = 0.0d0 ! Reset ecapbind

        endif


        return
        end

