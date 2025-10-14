************************************************************************
*                                                                      *
      integer function nSmpSpNuDistData(isotope, nCf252option)
*                                                                      *
************************************************************************
c
c Copyright (c) 2006 The Regents of the University of California.
c All rights reserved.
c UCRL-CODE-224807
c
c This software was developed by Lawrence Livermore National Laboratory.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c 1. Redistributions of source code must retain the above copyright notice,
c    this list of conditions and the following disclaimer.
c 2. Redistributions in binary form must reproduce the above copyright notice,
c    this list of conditions and the following disclaimer in the documentation
c    and/or other materials provided with the distribution.
c 3. The name of the author may not be used to endorse or promote products
c    derived from this software without specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
c EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
c OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
c WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
c OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
c ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c
c  Description
c    Sample Number of Neutrons from spontaneous fission
c    (a) from the neutron multiplicity data for
c        U-238, Pu-238, Pu-240, Pu-242, Cm-242, Cm-244
c           using Holden and Zucker's tabulated data
c        Cf-252 using either Spencer's tabulated data or
c           Boldeman's data
c    (b) from Terrell's approximation using nubar for
c        Th-232,
c        U-232, U-233, U-234, U-235, U-236,
c        Np-237,
c        Pu-239, Pu-241,
c        Am-241,
c        Bk-249
c           using Ensslin's data.
c  Input
c    iso          - isotope
c    Cf252option  - 0 to use Spencer's tabulated data
c                   1 to use Boldeman's data
c  Output
c    SmpSpNuDistData - sampled multiplicity
c                      -1 is the isotope has
c                         no multiplicity data,
c                         nor any nubar data
c
c Convert from c++ to g77 for phits (liemph 2008)
c
      parameter (nSPfissIso = 8)
      parameter (nSPfissn = 11)
c
      integer isotope, nCf252option
      integer i, index, iflag
      real*8  sum, nubar
      real*8  r
      real*8  SmpSpNubarData
      real*8  fisslibrng, unirn
c
      real*8 sfnu(0:nSPfissIso,0:nSPfissn)
      data ((sfnu(i,j),j=0,nSPfissn-1),i=0,nSPfissIso-1) /
     & 0.0481677,0.2485215,0.4253044,0.2284094,0.0423438,0.0072533,
     & 0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
     & 0.0631852,0.2319644,0.3333230,0.2528207,0.0986461,0.0180199,
     & 0.0020407,0.0000000,0.0000000,0.0000000,0.0000000,
     & 0.0679423,0.2293159,0.3341228,0.2475507,0.0996922,0.0182398,
     & 0.0031364,0.0000000,0.0000000,0.0000000,0.0000000,
     & 0.0212550,0.1467407,0.3267531,0.3268277,0.1375090,0.0373815,
     & 0.0025912,0.0007551,0.0001867,0.0000000,0.0000000,
     & 0.0150050,0.1161725,0.2998427,0.3331614,0.1837748,0.0429780,
     & 0.0087914,0.0002744,0.0000000,0.0000000,0.0000000,
     & 0.0540647,0.2053880,0.3802279,0.2248483,0.1078646,0.0276366,
     & 0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,
     & 0.0021100,0.0246700,0.1229000,0.2714400,0.3076300,0.1877000,
     & 0.0677000,0.0140600,0.0016700,0.0001000,0.0000000,
     & 0.0020900,0.0262100,0.1262000,0.2752000,0.3018000,0.1846000,
     & 0.0668000,0.0150000,0.0021000,0.0000000,0.0000000 /
c
c  sample the spontaneous fission neutron number distribution
c
      index = -1
c
      if      (isotope .eq. 92238) then
        index = 0
      else if (isotope .eq. 94240) then
        index = 1
      else if (isotope .eq. 94242) then
        index = 2
      else if (isotope .eq. 96242) then
        index = 3
      else if (isotope .eq. 96244) then
        index = 4
      else if (isotope .eq. 94238) then
        index = 5
      else if (isotope .eq. 98252 .and. nCf252option .eq. 0) then
        index = 6
      else if (isotope .eq. 98252 .and. nCf252option .eq. 1) then
        index = 7
      endif
c
      if (index .ne. -1) then
c       CALL RANDOM_NUMBER(fisslibrng)
        fisslibrng = unirn(dummy)
        r=fisslibrng
        sum=0.0d0
c
c liemph
c to start do while loop correctly
        iflag = 1
        i = 0
c liemph
c
        do while ( iflag.eq.1 .and. i.lt.nSPfissn )
          sum=sum+sfnu(index,i)
          if (r .le. sum .or. sfnu(index,i+1) .eq. 0.0d0) then
              nSmpSpNuDistData = i
              iflag = 0
          endif
         i=i+1
        enddo
      else
c
c There is no full multiplicity distribution data available
c for that isotope, let's try to find a nubar for it in
c N. Ensslin, et.al., "Application Guide to Neutron
c Multiplicity Counting," LA-13422-M (November 1998)
c and use Terrell's approximation
c
        nubar = SmpSpNubarData(isotope)
        if (nubar .ne. -1.0d0) then
          nSmpSpNuDistData = nSmpTerrell(nubar)
        else
c There is no nubar information for that isotope, return -1,
c meaning no data available for that isotope
        nSmpSpNuDistData = -1
        endif
      endif
c
      end


************************************************************************
*                                                                      *
      real*8 function SmpWatt(ePart, iso)
*                                                                      *
************************************************************************
c
c Copyright (c) 2006 The Regents of the University of California.
c All rights reserved.
c UCRL-CODE-224807
c
c This software was developed by Lawrence Livermore National Laboratory.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c 1. Redistributions of source code must retain the above copyright notice,
c    this list of conditions and the following disclaimer.
c 2. Redistributions in binary form must reproduce the above copyright notice,
c    this list of conditions and the following disclaimer in the documentation
c    and/or other materials provided with the distribution.
c 3. The name of the author may not be used to endorse or promote products
c    derived from this software without specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
c EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
c OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
c WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
c OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
c ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c
c  Description
c    Sample Watt Spectrum as in TART (Kalos algorithm)
c
c  Input
c    ePart     - energy of incoming particle
c    iso       - isotope
c  Output
c              - energy of incoming particle
c
c     38 fissionable isotopes in ENDL + U-232 from ENDF7
c
c Convert from c++ to g77 for phits (liemph 2008)
      implicit real*8 (a-h,o-z)

      include 'err.inc'
c
      parameter (nZAfis=39)
      parameter (WATTEMIN=1.0d-6)
      parameter (WATTEMAX=20.0)
c
      integer nZA(0:nZAfis)
      data (nZA(i),i=0,nZAfis-1) /
     &   90231, 90232, 90233,
     &   91233,
     &   92232, 92233, 92234, 92235, 92236, 92237, 92238, 92239, 92240,
     &   93235, 93236, 93237, 93238,
     &   94237, 94238, 94239, 94240, 94241, 94242, 94243,
     &   95241, 95242, 95243,
     &   96242, 96243, 96244, 96245, 96246, 96247, 96248,
     &   97249,
     &   98249, 98250, 98251, 98252 /
c
      real*8 Watta (0:nZAfis,0:2)
      data ((Watta(i,j),j=0,2),i=0,nZAfis-1) /
     &   6.00949285d-05, -8.36695381d-03,  9.50939496d-01,
     &   6.54348443d-05, -8.86574327d-03,  9.55404490d-01,
     &   7.08173682d-05, -9.22676286d-03,  9.50088329d-01,
     &   6.35839062d-05, -8.63645973d-03,  9.24583535d-01,
     &   8.21929628d-05,  4.01922936d-03,  1.152121164d00,
     &   6.21335718d-05, -8.45651858d-03,  9.14717276d-01,
     &   6.81386135d-05, -8.99142394d-03,  9.21954824d-01,
     &   7.32627297d-05, -9.36908697d-03,  9.20107976d-01,
     &   8.06505279d-05, -9.95416671d-03,  9.27890410d-01,
     &   8.33208285d-05, -1.01073057d-02,  9.17691654d-01,
     &   8.96944680d-05, -1.06491070d-02,  9.25496030d-01,
     &   9.44608097d-05, -1.08940419d-02,  9.17795511d-01,
     &   1.01395704d-04, -1.15098159d-02,  9.29395462d-01,
     &   6.81110009d-05, -8.91619352d-03,  9.00047566d-01,
     &   7.21126359d-05, -9.20179363d-03,  8.95722889d-01,
     &   7.82371142d-05, -9.67050621d-03,  8.99574933d-01,
     &   8.27256297d-05, -9.99353009d-03,  8.97461897d-01,
     &   7.29458059d-05, -9.22415170d-03,  8.80996165d-01,
     &   8.02383914d-05, -9.78291439d-03,  8.88964070d-01,
     &   8.50641730d-05, -1.01099145d-02,  8.87304833d-01,
     &   9.10537157d-05, -1.05303084d-02,  8.89438514d-01,
     &   9.43014320d-05, -1.07133543d-02,  8.82632055d-01,
     &   1.02655616d-04, -1.13154691d-02,  8.91617174d-01,
     &   1.06118094d-04, -1.14971777d-02,  8.85181637d-01,
     &   9.08474473d-05, -1.04296303d-02,  8.71942958d-01,
     &   9.35633054d-05, -1.05612167d-02,  8.63930371d-01,
     &   1.01940441d-04, -1.11573929d-02,  8.73153437d-01,
     &   9.19501202d-05, -1.04229157d-02,  8.58681822d-01,
     &   9.42991674d-05, -1.05098872d-02,  8.49103546d-01,
     &   1.02747171d-04, -1.11371417d-02,  8.60434431d-01,
     &   1.05024967d-04, -1.12138980d-02,  8.51101942d-01,
     &   1.14130011d-04, -1.18692049d-02,  8.62838259d-01,
     &   1.15163673d-04, -1.18553822d-02,  8.51306646d-01,
     &   1.27169055d-04, -1.27033210d-02,  8.68623539d-01,
     &   1.24195213d-04, -1.24047085d-02,  8.48974077d-01,
     &   1.12616150d-04, -1.15135023d-02,  8.19708800d-01,
     &   1.23637465d-04, -1.22869889d-02,  8.35392018d-01,
     &   1.22724317d-04, -1.21677963d-02,  8.22569523d-01,
     &   1.33891595d-04, -1.29267762d-02,  8.37122909d-01 /
c
      real*8 ePart
      integer iso
c
      real*8 a
      real*8 b
      real*8 rand1,rand2, fisslibrng, dummy1, dummy2
      real*8 x,y,z
      real*8 eSmp
      integer i, iflag
c
      b=1.0d0
c
c   Find Watt parameters for isotope
c
      isoindex=-1

c    Iwamoto 2015/12/14
      do i=0,nZAfis-1
        if (iso .eq. nZA(i) ) isoindex = i
      enddo
      if (isoindex .eq. -1) then
       write(ErrCha,*) 'SmpWatt: No Watt spectrum available for iso ',iso
       ErrID = 'L:284/R:SmpWatt/F:fismul.f' !E03_001_001
       call ErrWrite(ErrID,ErrCha)
       stop
      endif
c
      a= Watta(isoindex,2) + ePart*( Watta(isoindex,1)
     &                     + ePart*Watta(isoindex,0) )
      x= 1.0d0 + (b/(8.d0*a))
      y= (x + dsqrt(x*x-1.0d0))/a
      z= a*y - 1.0d0
c
c
c liemph
c to start do while loop correctly
      iflag = 1
c liemph
      do while (iflag .eq.1 )
c       CALL RANDOM_NUMBER(fisslibrng)
        fisslibrng = unirn(dummy)
        rand1= -dlog(fisslibrng)
c       CALL RANDOM_NUMBER(fisslibrng)
        fisslibrng = unirn(dummy)
        rand2= -dlog(fisslibrng)
        eSmp= y*rand1
c
        dummy1 = (rand2-z*(rand1+1.d0))*(rand2-z*(rand1+1.d0))
        dummy2 = b*y*rand1
        if ( dummy1.gt.dummy2 .or.
     &       eSmp.lt.WATTEMIN .or.
     &       eSmp.gt.WATTEMAX) then
             iflag=1
         else
           iflag = 0
         endif
      enddo
c
      SmpWatt=eSmp
c
      return
      end


************************************************************************
*                                                                      *
      integer function nSmpNuDistDataPu239_241_MC(nubar)
*                                                                      *
************************************************************************
c
c Copyright (c) 2006 The Regents of the University of California.
c All rights reserved.
c UCRL-CODE-224807
c
c This software was developed by Lawrence Livermore National Laboratory.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c 1. Redistributions of source code must retain the above copyright notice,
c    this list of conditions and the following disclaimer.
c 2. Redistributions in binary form must reproduce the above copyright notice,
c    this list of conditions and the following disclaimer in the documentation
c    and/or other materials provided with the distribution.
c 3. The name of the author may not be used to endorse or promote products
c    derived from this software without specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
c EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
c OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
c WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
c OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
c ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c
c
c  Description
c    Sample Number of Neutrons from fission in Pu-239 and Pu-241 using
c    Zucker and Holden's tabulated data for Pu-239
c    The 11 P(nu) distributions are given as a function of nubar,
c    the average number of neutrons from induced fission for the
c    11 different energies (0 to 10 MeV), based on the Pu-239 data
c    from Zucker and Holden.
c
c  Input
c    nubar    - average number of neutrons per fission
c  Output
c    SmpNuDistDataPu239_241_MC  - sampled multiplicity
c
c Convert from c++ to g77 for phits (liemph 2008)
c
      real*8 nubar
      real*8  fisslibrng, unirn
c
      real*8 Pu239nu(0:11,0:9)
      data ((Pu239nu(i,j),j=0,8),i=0,10) /
     & .0108826, .0994916, .2748898, .3269196, .2046061,
     & .0726834, .0097282, .0006301, .0001685 ,
     & .0084842, .0790030, .2536175, .3289870, .2328111,
     & .0800161, .0155581, .0011760, .0003469 ,
     & .0062555, .0611921, .2265608, .3260637, .2588354,
     & .0956070, .0224705, .0025946, .0005205 ,
     & .0045860, .0477879, .1983002, .3184667, .2792811,
     & .1158950, .0301128, .0048471, .0007233 ,
     & .0032908, .0374390, .1704196, .3071862, .2948565,
     & .1392594, .0386738, .0078701, .0010046 ,
     & .0022750, .0291416, .1437645, .2928006, .3063902,
     & .1641647, .0484343, .0116151, .0014149 ,
     & .0014893, .0222369, .1190439, .2756297, .3144908,
     & .1892897, .0597353, .0160828, .0029917 ,
     & .0009061, .0163528, .0968110, .2558524, .3194566,
     & .2134888, .0729739, .0213339, .0020017 ,
     & .0004647, .0113283, .0775201, .2335926, .3213289,
     & .2356614, .0886183, .0274895, .0039531 ,
     & .0002800, .0071460, .0615577, .2089810, .3200121,
     & .2545846, .1072344, .0347255, .0054786 ,
     & .0002064, .0038856, .0492548, .1822078, .3154159,
     & .2687282, .1295143, .0432654, .0075217 /
c
      real*8 Pu239nubar(0:11)
      data (Pu239nubar(i),i=0,10) /
     & 2.8760000d0,
     & 3.0088800d0,
     & 3.1628300d0,
     & 3.3167800d0,
     & 3.4707300d0,
     & 3.6246800d0,
     & 3.7786300d0,
     & 3.9325800d0,
     & 4.0865300d0,
     & 4.2404900d0,
     & 4.3944400d0 /
c
      real*8 fraction, r, cum
      integer engind, nu
c
c  Check if nubar is within the range of experimental values
c
      if (nubar .ge. Pu239nubar(0) .and. nubar .le. Pu239nubar(10)) then
c
c     Use Zucker and Holden Data
c
        engind = 1
        do while (nubar .gt. Pu239nubar(engind) )
          engind = engind + 1
        enddo
        fraction = (nubar-Pu239nubar(engind-1))/
     &             (Pu239nubar(engind)-Pu239nubar(engind-1))
c       CALL RANDOM_NUMBER(fisslibrng)
        fisslibrng = unirn(dummy)
        if (fisslibrng .gt. fraction) engind = engind - 1
c       CALL RANDOM_NUMBER(fisslibrng)
        fisslibrng = unirn(dummy)
        r = fisslibrng
        nu = 0
        cum = Pu239nu(engind,0)
        do while ( r.gt.cum .and. nu.lt.8 )
          nu = nu + 1
          cum = cum + Pu239nu(engind,nu)
        enddo
        nSmpNuDistDataPu239_241_MC = nu
c
      else
c
c     Use Terrell's formula
c
        nSmpNuDistDataPu239_241_MC = nSmpTerrell(nubar)
      endif
c
      end




************************************************************************
*                                                                      *
      integer function nSmpNuDistDataU232_234_236_238_MC(nubar)
*                                                                      *
************************************************************************
c
c Copyright (c) 2006 The Regents of the University of California.
c All rights reserved.
c UCRL-CODE-224807
c
c This software was developed by Lawrence Livermore National Laboratory.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c 1. Redistributions of source code must retain the above copyright notice,
c    this list of conditions and the following disclaimer.
c 2. Redistributions in binary form must reproduce the above copyright notice,
c    this list of conditions and the following disclaimer in the documentation
c    and/or other materials provided with the distribution.
c 3. The name of the author may not be used to endorse or promote products
c    derived from this software without specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
c EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
c OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
c WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
c OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
c ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c
c  Description
c    Sample Number of Neutrons from fission in U-232, U-234, U-236,
c    and U-238 using Zucker and Holden's tabulated data for U-238
c    The 11 P(nu) distributions are given as a function of nubar,
c    the average number of neutrons from induced fission for the
c    11 different energies (0 to 10 MeV), based on the U-238 data
c    from Zucker and Holden.
c
c  Input
c    nubar    - average number of neutrons per fission
c  Output
c    SmpNuDistDataU232_234_236_238_MC  - sampled multiplicity
c
c Convert from c++ to g77 for phits (liemph 2008)
c
      real*8 nubar
      real*8  fisslibrng, unirn
c
      real*8 U238nu(0:11,0:9)
      data ((U238nu(i,j),j=0,8),i=0,10) /
     & .0396484, .2529541, .2939544, .2644470, .1111758,
     & .0312261, .0059347, .0005436, .0001158 ,
     & .0299076, .2043215, .2995886, .2914889, .1301480,
     & .0363119, .0073638, .0006947, .0001751 ,
     & .0226651, .1624020, .2957263, .3119098, .1528786,
     & .0434233, .0097473, .0009318, .0003159 ,
     & .0170253, .1272992, .2840540, .3260192, .1779579,
     & .0526575, .0130997, .0013467, .0005405 ,
     & .0124932, .0984797, .2661875, .3344938, .2040116,
     & .0640468, .0173837, .0020308, .0008730 ,
     & .0088167, .0751744, .2436570, .3379711, .2297901,
     & .0775971, .0225619, .0030689, .0013626 ,
     & .0058736, .0565985, .2179252, .3368863, .2541575,
     & .0933127, .0286200, .0045431, .0031316 ,
     & .0035997, .0420460, .1904095, .3314575, .2760413,
     & .1112075, .0355683, .0065387, .0031316 ,
     & .0019495, .0309087, .1625055, .3217392, .2943792,
     & .1313074, .0434347, .0091474, .0046284 ,
     & .0008767, .0226587, .1356058, .3076919, .3080816,
     & .1536446, .0522549, .0124682, .0067176 ,
     & .0003271, .0168184, .1111114, .2892434, .3160166,
     & .1782484, .0620617, .0166066, .0095665 /
c
      real*8 U238nubar(0:11)
      data (U238nubar(i),i=0,10) /
     & 2.2753781d0,
     & 2.4305631d0,
     & 2.5857481d0,
     & 2.7409331d0,
     & 2.8961181d0,
     & 3.0513031d0,
     & 3.2064881d0,
     & 3.3616731d0,
     & 3.5168581d0,
     & 3.6720432d0,
     & 3.8272281d0  /
c
      real*8 fraction, r, cum
      integer engind, nu
c
c  Check if nubar is within the range of experimental values
c
      if (nubar .ge. U238nubar(0) .and. nubar .le. U238nubar(10)) then
c
c     Use Zucker and Holden Data
c
        engind = 1
        do while (nubar .gt. U238nubar(engind) )
          engind = engind + 1
        enddo
        fraction = (nubar-U238nubar(engind-1))/
     &             (U238nubar(engind)-U238nubar(engind-1))
c       CALL RANDOM_NUMBER(fisslibrng)
        fisslibrng = unirn(dummy)
        if (fisslibrng .gt. fraction) engind = engind - 1
c       CALL RANDOM_NUMBER(fisslibrng)
        fisslibrng = unirn(dummy)
        r = fisslibrng
        nu = 0
        cum = U238nu(engind,0)
        do while ( r.gt.cum .and. nu.lt.8 )
          nu = nu + 1
          cum = cum + U238nu(engind,nu)
        enddo
        nSmpNuDistDataU232_234_236_238_MC = nu
c
      else
c
c     Use Terrell's formula
c
        nSmpNuDistDataU232_234_236_238_MC = nSmpTerrell(nubar)
      endif
c
      end




************************************************************************
*                                                                      *
      integer function nSmpNuDistDataU233_235_MC(nubar)
*                                                                      *
************************************************************************
c
c Copyright (c) 2006 The Regents of the University of California.
c All rights reserved.
c UCRL-CODE-224807
c
c This software was developed by Lawrence Livermore National Laboratory.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c 1. Redistributions of source code must retain the above copyright notice,
c    this list of conditions and the following disclaimer.
c 2. Redistributions in binary form must reproduce the above copyright notice,
c    this list of conditions and the following disclaimer in the documentation
c    and/or other materials provided with the distribution.
c 3. The name of the author may not be used to endorse or promote products
c    derived from this software without specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
c EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
c OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
c WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
c OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
c ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c
c  Description
c    Sample Number of Neutrons from fission in U-233 and U-235 using
c    (a) Gwin, Spencer and Ingle tabulated data at thermal
c        energies (0 MeV),
c    (b) Zucker and Holden's tabulated data for U-235 at 1 MeV and
c        higher.
c    The 11 P(nu) distributions are given as a function of nubar,
c    the average number of neutrons from induced fission for the
c    11 different energies (0 to 10 MeV), based on the U-235 data
c    above.
c
c
c  Input
c    nubar    - average number of neutrons per fission
c  Output
c    SmpNuDistDataU233_235_MC  - sampled multiplicity
c
c Convert from c++ to g77 for phits (liemph 2008)
c
      real*8 nubar
      real*8  fisslibrng, unirn
c
      real*8 U235nu(0:11,0:8)
      data ((U235nu(i,j),j=0,7),i=0,10) /
     & .0291000, .1660000, .3362000, .3074000, .1333000,
     & .0259000, .0021000, .0002000 ,
     & .0237898, .1555525, .3216515, .3150433, .1444732,
     & .0356013, .0034339, .0004546 ,
     & .0183989, .1384891, .3062123, .3217566, .1628673,
     & .0455972, .0055694, .0011093 ,
     & .0141460, .1194839, .2883075, .3266568, .1836014,
     & .0569113, .0089426, .0019504 ,
     & .0115208, .1032624, .2716849, .3283426, .2021206,
     & .0674456, .0128924, .0027307 ,
     & .0078498, .0802010, .2456595, .3308175, .2291646,
     & .0836912, .0187016, .0039148 ,
     & .0046272, .0563321, .2132296, .3290407, .2599806,
     & .1045974, .0265604, .0056322 ,
     & .0024659, .0360957, .1788634, .3210507, .2892537,
     & .1282576, .0360887, .0079244 ,
     & .0012702, .0216090, .1472227, .3083032, .3123950,
     & .1522540, .0462449, .0107009 ,
     & .0007288, .0134879, .1231200, .2949390, .3258251,
     & .1731879, .0551737, .0135376 ,
     & .0004373, .0080115, .1002329, .2779283, .3342611,
     & .1966100, .0650090, .0175099 /
c
      real*8 U235nubar(0:11)
      data (U235nubar(i),i=0,10) /
     & 2.4370000d0,
     & 2.5236700d0,
     & 2.6368200d0,
     & 2.7623400d0,
     & 2.8738400d0,
     & 3.0386999d0,
     & 3.2316099d0,
     & 3.4272800d0,
     & 3.6041900d0,
     & 3.7395900d0,
     & 3.8749800d0 /
      real*8 fraction, r, cum
      integer engind, nu
c
c  Check if nubar is within the range of experimental values
c
      if (nubar .ge. U235nubar(0) .and. nubar .le. U235nubar(10)) then
c
c     Use Zucker and Holden Data
c
        engind = 1
        do while (nubar .gt. U235nubar(engind) )
          engind = engind + 1
        enddo
        fraction = (nubar-U235nubar(engind-1))/
     &             (U235nubar(engind)-U235nubar(engind-1))
c       CALL RANDOM_NUMBER(fisslibrng)
        fisslibrng = unirn(dummy)
        if (fisslibrng .gt. fraction) engind = engind - 1
c       CALL RANDOM_NUMBER(fisslibrng)
        fisslibrng = unirn(dummy)
        r = fisslibrng
        nu = 0
        cum = U235nu(engind,0)
        do while ( r.gt.cum .and. nu.lt.8 )
          nu = nu + 1
          cum = cum + U235nu(engind,nu)
        enddo
        nSmpNuDistDataU233_235_MC = nu
c
      else
c
c     Use Terrell's formula
c
        nSmpNuDistDataU233_235_MC = nSmpTerrell(nubar)
      endif
c
      end




************************************************************************
*                                                                      *
      integer function nSmpTerrell(nubar)
*                                                                      *
************************************************************************
c
c Copyright (c) 2006 The Regents of the University of California.
c All rights reserved.
c UCRL-CODE-224807
c
c This software was developed by Lawrence Livermore National Laboratory.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c 1. Redistributions of source code must retain the above copyright notice,
c    this list of conditions and the following disclaimer.
c 2. Redistributions in binary form must reproduce the above copyright notice,
c    this list of conditions and the following disclaimer in the documentation
c    and/or other materials provided with the distribution.
c 3. The name of the author may not be used to endorse or promote products
c    derived from this software without specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
c EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
c OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
c WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
c OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
c ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c
c  Description
c    Sample Fission Number from Terrell's modified Gaussian distribution
c
c    method uses Red Cullen's algoritm UCRL-TR-222526
c
c  Input
c    nubar       - average number of neutrons per fission
c  Output
c    SmpTerrell  - sampled multiplicity
c
c Convert from c++ to g77 for phits (liemph 2008)
c
      include 'err.inc'

      parameter (TWOPI=6.283185307)
      parameter (SQRT2=1.414213562)
      parameter (BSHIFT=-0.43287)
      parameter (WIDTH_C=1.079)
c
      real*8 nubar
      real*8 fisslibrng
c
      real*8 width
      real*8 temp1, temp2, expo, cshift
      real*8 rw, theta, sampleg, unirn
c
      if (nubar .lt. WIDTH_C) then
       write(ErrCha,*) 'SmpTerrell: fission nubar out of range'
       ErrID = 'L:785/R:nSmpTerrell/F:fismul.f' !E03_002_001
       call ErrWrite(ErrID,ErrCha)
      endif
c
      width = SQRT2 * WIDTH_C
      temp1 = nubar + 0.5
      temp2 = temp1/width
      temp2 = temp2*temp2
      expo  = dexp(-temp2)
      cshift = temp1 + BSHIFT * WIDTH_C * expo/(1. - expo)
c
c liemph
c to start do while loop correctly
      sampleg = -1.0d0
c liemph
c
      do while (sampleg .lt. 0d0)
c       CALL RANDOM_NUMBER(fisslibrng)
        fisslibrng = unirn(dummy)
        rw = dsqrt( -dlog(fisslibrng) )
c       CALL RANDOM_NUMBER(fisslibrng)
        fisslibrng = unirn(dummy)
        theta = TWOPI * fisslibrng
        sampleg = width * rw * dcos(theta) + cshift
      enddo
c
      nSmpTerrell = int(sampleg)
c
      return
      end

************************************************************************
*                                                                      *
      real*8 function SmpSpNubarData(isotope)
*                                                                      *
************************************************************************

c Copyright (c) 2006 The Regents of the University of California.
c All rights reserved.
c Copyright (c) 2006 The Regents of the University of California.
c All rights reserved.
c UCRL-CODE-224807
c
c This software was developed by Lawrence Livermore National Laboratory.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c 1. Redistributions of source code must retain the above copyright notice,
c    this list of conditions and the following disclaimer.
c 2. Redistributions in binary form must reproduce the above copyright notice,
c    this list of conditions and the following disclaimer in the documentation
c    and/or other materials provided with the distribution.
c 3. The name of the author may not be used to endorse or promote products
c    derived from this software without specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
c EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
c OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
c WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
c OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
c ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c UCRL-CODE-224807
c
c This software was developed by Lawrence Livermore National Laboratory.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c 1. Redistributions of source code must retain the above copyright notice,
c    this list of conditions and the following disclaimer.
c 2. Redistributions in binary form must reproduce the above copyright notice,
c    this list of conditions and the following disclaimer in the documentation
c    and/or other materials provided with the distribution.
c 3. The name of the author may not be used to endorse or promote products
c    derived from this software without specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
c EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
c OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
c WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
c OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
c ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c
c   Description
c     Determine average number of neutrons from spontaneous fission for
c         Th-232,
c         U-232, U-233, U-234, U-235, U-236, U-238
c         Np-237,
c         Pu-239, Pu-240, Pu-241,  Pu-242
c         Am-241,
c         Cm-242, Cm-244,
c         Bk-249,
c         Cf-252
c     Based on Ensslin's data.
c     N. Ensslin, et.al., "Application Guide to Neutron Multiplicity Counting,"
c     LA-13422-M (November 1998)
c
c   Input
c     iso          - isotope
c   Output
c     SmpSpNubarData - average number of neutrons
c                      -1. is the isotope has
c                          no nubar data
c
c Convert from c++ to g77 for phits (liemph 2008)
c

      parameter (nSPfissn = 18)

      integer isotope, i
      integer spzaid (nSPfissn)
      real*8  spnubar(nSPfissn)

      data (spzaid(i),i=1,nSPfissn) /
     & 90232, 92232, 92233, 92234, 92235,
     & 92236, 92238, 93237, 94238, 94239,
     & 94240, 94241, 94242, 95241, 96242,
     & 96244, 97249, 98252/

      data (spnubar(i),i=1,nSPfissn) /
     & 2.14,  1.71, 1.76,  1.81, 1.86,
     & 1.91,  2.01, 2.05,  2.21, 2.16,
     & 2.156, 2.25, 2.145, 3.22, 2.54,
     & 2.72,  3.40, 3.757/

            SmpSpNubarData = -1.0d0

      do i = 1, nSPfissn

         if( isotope .eq. spzaid(i) ) then

            SmpSpNubarData = spnubar(i)
            return

         end if

      end do

      return
      end

