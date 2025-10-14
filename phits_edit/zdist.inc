      parameter (mcsize=6)  ! number of input cell size
      parameter (mcsizeele=7)  ! number of input cell size for electron
      parameter (neinc=14)  ! number of input incident energy
      parameter (npart=4)
      parameter (mpara=7)   ! number of parameter
      parameter (mcst=2)    ! number of constant parameter
      parameter (npeak=9)   ! number of peak
      parameter (maxedist=70) ! number of z-bin used for normalization
      parameter (maxAuger=100) ! number of z-bin used for Auger calculation
      parameter (Wauger=100.0e-6)  ! Width of Auger Peak
      parameter (neele=10) ! number of electron incident energy

      common /fitcom/
     & Csize(mcsize),                  ! cell size
     & Csizeele(mcsizeele),            ! cell size for electron
     & Efit(neinc),                    ! incident kinetic energy (MeV/n)
     & Efitele(neele),                 ! incident electron energy (MeV/n)
     & Afit(mcsize,npart,neinc,mpara), ! Main Parameters for Nucleus
     & Pfit(mcsize,npart,neinc,npeak), ! Peak area for Nucleus
     & AfitEle(mcsizeele,neele,mpara), ! Main Parameters for Electron
     & PfitEle(mcsizeele,neele,npeak), ! Peak area for Electron
     & Acurr(mpara),Ccurr(mcst),Pcurr(npeak), ! current value of the parameters
     & Pmid(npeak),                    ! center of peak (MeV)
     & zlbin(maxedist+1),              ! lower boundary of z bin used for normalization
     & Auger(-maxAuger:maxAuger)       ! used for Auger peak width calculation
!$OMP THREADPRIVATE(/fitcom/)

      common /intfitcom/
     & iZfit(npart),  ! z number of particle
     & iponly         ! below this peak, no fitting data

