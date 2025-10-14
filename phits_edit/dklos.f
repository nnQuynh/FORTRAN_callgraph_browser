************************************************************************
*                                                                      *
      subroutine dklos(iprj,kprj,eein)
*                                                                      *
*                                                                      *
*       particle decay                                                 *
*       modified by K.Niita on 2013/08/22                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       iprj    : particle type of projectile                          *
*       kprj    : kf code of projectile                                *
*       eein    : energy of projectile (MeV)                           *
*                                                                      *
*     output :  in common from dklos2, dklos3 and c9decay              *
*                                                                      *
*                                                                      *
************************************************************************

      use udm_Parameter
      use udm_Utility
      use udm_Manager
      use decay_utility

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)

      common /othdcp/ nodcp(2000,6), nodcn
      common /dcayp/  adcayp(20), bdcayp(20)

*-----------------------------------------------------------------------
! T.Sato, add neutron decay 2016/03/02
         if( iprj .le. 1 ) return

*-----------------------------------------------------------------------

         kdcay = 0

*-----------------------------------------------------------------------
*     user defined particle decay
*
*     "user defined particle" function can modify any particle decay
*     pattern. So, this function is written at the beginning of dklos.
*-----------------------------------------------------------------------

         if(udm_part_num .gt. 0) then
           udm_kf_for_21=kprj
           udm_Kin=eein
           set_decay_success=.false.
           call user_defined_particle(21)
           if(set_decay_success) goto 1000
           if( 900000 .le. abs(kprj) .and. abs(kprj) .le. 999999) then
             goto 1000
           endif
         endif

*-----------------------------------------------------------------------
*     Neutron decay to proton, electron and antielectron neutrino
*-----------------------------------------------------------------------

         if( iprj .eq. 2 ) then

                  kdcay = 3

                  ip1 = 1
                  kp1 = 2212

                  ip2 = 12
                  kp2 = 11

                  ip3 = 11
                  kp3 = -12

                  call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)

*-----------------------------------------------------------------------
*        9C decay to two alpha and proton through 8Be or 5Li
*-----------------------------------------------------------------------

         elseif( kprj .eq. 6000009 ) then

                     kdcay = 3

                     ip1 = 1
                     kp1 = 2212

                     ip2 = 18
                     kp2 = 2000004

                     ip3 = 18
                     kp3 = 2000004

                  call c9decay(eein)

*-----------------------------------------------------------------------
*        positron to two photons
*-----------------------------------------------------------------------

         else if( iprj .eq. 13 ) then

                     kdcay = 2

                     ip1 = 14
                     kp1 = 22

                     ip2 = 14
                     kp2 = 22

                  call egs5annih(eein)

*-----------------------------------------------------------------------
*        pi0 decay to two photons
*-----------------------------------------------------------------------

         else if( iprj .eq. 4 ) then

                     kdcay = 2

                     ip1 = 14
                     kp1 = 22

                     ip2 = 14
                     kp2 = 22

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*        pi+ and pi- decay to mu+, mu- and neutrino
*-----------------------------------------------------------------------

         else if( iprj .eq. 3 .or. iprj .eq. 5 ) then

                     kdcay = 2

                  if( iprj .eq. 3 ) then

                     ip1 = 6
                     kp1 = -13

                     ip2 = 11
                     kp2 = 14

                  else if( iprj .eq. 5 ) then

                     ip1 = 7
                     kp1 = 13

                     ip2 = 11
                     kp2 = -14

                  end if


                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*        muon : e + nu_e_bar + nu_mu
*-----------------------------------------------------------------------

         else if( iprj .eq. 6 .or. iprj .eq. 7 ) then

                     kdcay = 3

                  if( iprj .eq. 6 ) then

                     ip1 = 11
                     kp1 = 12

                     ip2 = 11
                     kp2 = -14

                     ip3 = 13
                     kp3= -11

                  else if( iprj .eq. 7 ) then

                     ip1 = 11
                     kp1 = -12

                     ip2 = 11
                     kp2 = 14

                     ip3 = 12
                     kp3 = 11

                  end if


                  call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)

*-----------------------------------------------------------------------
*        K_L0 decay (kf=130)
*        decay channel : (1) 20.275 % : pi+ e-  nu_e_bar
*                        (2) 20.275 % : pi- e+  nu_e
*                        (3) 13.52  % : pi+ mu- nu_mu_bar
*                        (4) 13.52  % : pi- mu+ nu_mu
*                        (5) 19.52  % : pi0 pi0 pi0
*                        (6) 12.55  % : pi+ pi- pi0
*                        sum = 99.66
*        [Ref] https://ccwww.kek.jp/pdg/2020/tables/rpp2020-tab-mesons-strange.pdf
*-----------------------------------------------------------------------
         else if( kprj .eq. 130 ) then

             Br1=20.275
             Br2=20.275
             Br3=13.52
             Br4=13.52
             Br5=19.52
             Br6=12.55
             BrTOT=Br1+Br2+Br3+Br4+Br5+Br6
             Br1=Br1/BrTOT
             Br2=Br2/BrTOT
             Br3=Br3/BrTOT
             Br4=Br4/BrTOT
             Br5=Br5/BrTOT
             Br6=Br6/BrTOT

             prob = rn(0)

             ! ----------------------------------
             ! (1) 20.275 % : pi+ e-  nu_e_bar
             if( prob .lt. Br1 ) then
                 kdcay = 3
                 ip1 =  3;  kp1 =  211 ! pi+
                 ip2 = 12;  kp2 =   11 ! e-
                 ip3 = 11;  kp3 =  -12 ! nu_e_bar
                 call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)
             ! ----------------------------------
             ! (2) 20.275 % : pi- e+  nu_e
             else if( prob .lt. Br1+Br2 ) then
                 kdcay = 3
                 ip1 =  5;  kp1 = -211 ! pi-
                 ip2 = 13;  kp2 =  -11 ! e+
                 ip3 = 11;  kp3 =   12 ! nu_e
                 call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)
             ! ----------------------------------
             ! (3) 13.52  % : pi+ mu- nu_mu_bar
             else if( prob .lt. Br1+Br2+Br3 ) then
                 kdcay = 3
                 ip1 =  3;  kp1 =  211 ! pi+
                 ip2 =  7;  kp2 =   13 ! mu-
                 ip3 = 11;  kp3 =  -14 ! nu_mu_bar
                 call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)
             ! ----------------------------------
             ! (4) 13.52  % : pi- mu+ nu_mu
             else if( prob .lt. Br1+Br2+Br3+Br4 ) then
                 kdcay = 3
                 ip1 =  5;  kp1 = -211 ! pi-
                 ip2 =  6;  kp2 =  -13 ! mu+
                 ip3 = 11;  kp3 =   14 ! nu_mu
                 call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)
             ! ----------------------------------
             ! (5) 19.52  % : pi0 pi0 pi0
             else if( prob .lt. Br1+Br2+Br3+Br4+Br5 ) then
                 kdcay = 3
                 ip1 =  4;  kp1 =  111 ! pi0
                 ip2 =  4;  kp2 =  111 ! pi0
                 ip3 =  4;  kp3 =  111 ! pi0
                 call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)
             ! ----------------------------------
             ! (6) 12.55  % : pi+ pi- pi0
             else
                 kdcay = 3
                 ip1 =  3;  kp1 =  211 ! pi+
                 ip2 =  5;  kp2 = -211 ! pi-
                 ip3 =  4;  kp3 =  111 ! pi0
                 call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)
             end if

*-----------------------------------------------------------------------
*        K_S0 decay (kf=310)
*        decay channel : (1) 30.69 % : pi0 pi0
*                        (2) 69.20 % : pi+ pi-
*                        sum = 99.89
*        [Ref] https://ccwww.kek.jp/pdg/2020/tables/rpp2020-tab-mesons-strange.pdf
*-----------------------------------------------------------------------
         else if( kprj .eq. 310 ) then

             Br1=30.69
             Br2=69.20
             BrTOT=Br1+Br2
             Br1=Br1/BrTOT
             Br2=Br2/BrTOT

             prob = rn(0)

             ! ----------------------------------
             ! (1) 30.69 % : pi0 pi0
             if( prob .lt. Br1 ) then
                 kdcay = 2
                 ip1 =  4;  kp1 =  111 ! pi0
                 ip2 =  4;  kp2 =  111 ! pi0
                 call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)
             ! ----------------------------------
             ! (2) 69.20 % : pi+ pi-
             else
                 kdcay = 2
                 ip1 =  3;  kp1 =  211 ! pi+
                 ip2 =  5;  kp2 = -211 ! pi-
                 call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)
             end if

*-----------------------------------------------------------------------
*        k0, k0bar decay
*        decay channel : (1) 68.61 % : pi+,  pi-
*                        (2) 31.39 % : pi0, pi0
*                        (3) other   : two photons
*-----------------------------------------------------------------------

         else if(( iprj .eq.  9 .and. kprj .eq.  311 ) .or. ! y.sakaki 2021/08
     &           ( iprj .eq. 11 .and. kprj .eq. -311 ) ) then

                     kdcay = 2

                     prob = rn(0)

                  if( prob .lt. 0.6861 ) then

                     idec = 1

                  else if( prob .lt. 0.6861 + 0.3139 ) then

                     idec = 2

                  else

                     idec = 3

                  end if

*-----------------------------------------------------------------------
*           (1) pi+ , pi-
*-----------------------------------------------------------------------

            if( idec .eq. 1 ) then

                     ip1 = 3
                     kp1 = 211

                     ip2 = 5
                     kp2 = -211

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*           (2) pi0 , pi0
*-----------------------------------------------------------------------

            else if( idec .eq. 2 ) then

                     ip1 = 4
                     kp1 = 111

                     ip2 = 4
                     kp2 = 111

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*           (3) two photons
*-----------------------------------------------------------------------

            else if( idec .eq. 3 ) then

                     ip1 = 14
                     kp1 = 22

                     ip2 = 14
                     kp2 = 22

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

            end if

*-----------------------------------------------------------------------
*        kaon+ (iprj=8), kaon- (iprj=10) decay
*        decay channel (kaon+) : (1) 63.56  % : mu+     nu_mu
*                                (2)  5.07  % : pi0 e+  nu_e
*                                (3)  3.352 % : pi0 mu+ nu_mu
*                                (4) 20.67  % : pi+ pi0
*                                (5)  1.760 % : pi+ pi0 pi0
*                                (6)  5.583 % : pi+ pi+ pi-
*                                sum = 99.995
*        [Ref] https://pdg.lbl.gov/2022/tables/rpp2022-tab-mesons-strange.pdf
C ip= 3 (pi+)
C ip= 4 (pi0)
C ip= 5 (pi-)
C ip= 6 (mu+)
C ip= 7 (mu-)
C ip=12 (e-)
C ip=13 (e+)
*-----------------------------------------------------------------------
         else if( iprj .eq. 8 .or. iprj .eq. 10 ) then

             Br1=63.56d0
             Br2= 5.07d0
             Br3= 3.352d0
             Br4=20.67d0
             Br5= 1.760d0
             Br6= 5.583d0
             BrTOT=Br1+Br2+Br3+Br4+Br5+Br6
             Br1=Br1/BrTOT
             Br2=Br2/BrTOT
             Br3=Br3/BrTOT
             Br4=Br4/BrTOT
             Br5=Br5/BrTOT
             Br6=Br6/BrTOT

             prob = rn(0)
             ! ----------------------------------
             ! (1) 63.56  % : mu+     nu_mu
             if( prob .lt. Br1 ) then
                 kdcay = 2
                 if( iprj .eq. 8 ) then
                     ip1 =  6;  kp1 =  -13 ! mu+
                     ip2 = 11;  kp2 =   14 ! nu_mu
                 else
                     ip1 =  7;  kp1 =   13 ! mu-
                     ip2 = 11;  kp2 =  -14 ! nu_mu_bar
                 endif
                 call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)
             ! ----------------------------------
             ! (2)  5.07  % : pi0 e+  nu_e
             else if( prob .lt. Br1+Br2 ) then
                 kdcay = 3
                 if( iprj .eq. 8 ) then
                     ip1 =  4;  kp1 =  111 ! pi0
                     ip2 = 13;  kp2 =  -11 ! e+
                     ip3 = 11;  kp3 =   12 ! nu_e
                 else
                     ip1 =  4;  kp1 =  111 ! pi0
                     ip2 = 12;  kp2 =   11 ! e-
                     ip3 = 11;  kp3 =  -12 ! nu_e_bar
                 endif
                 call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)
             ! ----------------------------------
             ! (3)  3.352 % : pi0 mu+ nu_mu
             else if( prob .lt. Br1+Br2+Br3 ) then
                 kdcay = 3
                 if( iprj .eq. 8 ) then
                     ip1 =  4;  kp1 =  111 ! pi0
                     ip2 =  6;  kp2 =  -13 ! mu+
                     ip3 = 11;  kp3 =   14 ! nu_mu
                 else
                     ip1 =  4;  kp1 =  111 ! pi0
                     ip2 =  7;  kp2 =   13 ! mu-
                     ip3 = 11;  kp3 =  -14 ! nu_mu_bar
                 endif
                 call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)
             ! ----------------------------------
             ! (4) 20.67  % : pi+ pi0
             else if( prob .lt. Br1+Br2+Br3+Br4 ) then
                 kdcay = 2
                 if( iprj .eq. 8 ) then
                     ip1 =  3;  kp1 =  211 ! pi+
                     ip2 =  4;  kp2 =  111 ! pi0
                 else
                     ip1 =  5;  kp1 = -211 ! pi-
                     ip2 =  4;  kp2 =  111 ! pi0
                 endif
                 call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)
             ! ----------------------------------
             ! (5)  1.760 % : pi+ pi0 pi0
             else if( prob .lt. Br1+Br2+Br3+Br4+Br5 ) then
                 kdcay = 3
                 if( iprj .eq. 8 ) then
                     ip1 =  3;  kp1 =  211 ! pi+
                     ip2 =  4;  kp2 =  111 ! pi0
                     ip3 =  4;  kp3 =  111 ! pi0
                 else
                     ip1 =  5;  kp1 = -211 ! pi-
                     ip2 =  4;  kp2 =  111 ! pi0
                     ip3 =  4;  kp3 =  111 ! pi0
                 endif
                 call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)
             ! ----------------------------------
             ! (6)  5.583 % : pi+ pi+ pi-
             else
                 kdcay = 3
                 if( iprj .eq. 8 ) then
                     ip1 =  3;  kp1 =  211 ! pi+
                     ip2 =  3;  kp2 =  211 ! pi+
                     ip3 =  5;  kp3 = -211 ! pi-
                 else
                     ip1 =  5;  kp1 = -211 ! pi-
                     ip2 =  5;  kp2 = -211 ! pi-
                     ip3 =  3;  kp3 =  211 ! pi+
                 endif
                 call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)
             end if

*-----------------------------------------------------------------------
*        eta(221) decay
*        decay channel : (1) 38.9 % : two photons
*                        (2) 31.9 % : three pi0
*                        (3) 23.7 % : pi+ pi- pi0
*                        (4) other  : pi+ pi- photon
*-----------------------------------------------------------------------

         else if( iprj .eq. 11 .and. kprj .eq. 221 ) then

                     prob = rn(0)

                  if( prob .lt. 0.389 ) then

                     idec = 1
                     kdcay = 2

                  else if( prob .lt. 0.389 + 0.319 ) then

                     idec = 2
                     kdcay = 3

                  else if( prob .lt. 0.389 + 0.319 + 0.237 ) then

                     idec = 3
                     kdcay = 3

                  else

                     idec = 4
                     kdcay = 3

                  end if

*-----------------------------------------------------------------------
*           (1) two photons
*-----------------------------------------------------------------------

               if( idec .eq. 1 ) then

                     ip1 = 14
                     kp1 = 22

                     ip2 = 14
                     kp2 = 22

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*           (2) three pi0
*-----------------------------------------------------------------------

               else if( idec .eq. 2 ) then

                     ip1 = 4
                     kp1 = 111

                     ip2 = 4
                     kp2 = 111

                     ip3 = 4
                     kp3 = 111

                  call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)

*-----------------------------------------------------------------------
*           (3) pi+ pi- pi0
*-----------------------------------------------------------------------

               else if( idec .eq. 3 ) then

                     ip1 = 3
                     kp1 = 211

                     ip2 = 4
                     kp2 = 111

                     ip3 = 5
                     kp3 = -211

                  call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)

*-----------------------------------------------------------------------
*           (3) pi+ pi- photon
*-----------------------------------------------------------------------

               else

                     ip1 = 3
                     kp1 = 211

                     ip2 = 5
                     kp2 = -211

                     ip3 = 14
                     kp3 = 22

                  call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)

               end if

*-----------------------------------------------------------------------
*        eta'(331) decay
*        decay channel : (1) 44.1 % : pi+, pi-, eta
*                        (2) 20.5 % : pi0, pi0, eta
*                        (3) 30.1 % : pi+, pi-, photon
*                        (4) other  : two photons
*-----------------------------------------------------------------------

         else if( iprj .eq. 11 .and. kprj .eq. 331 ) then

                     prob = rn(0)

                  if( prob .lt. 0.441 ) then

                     idec = 1
                     kdcay = 3

                  else if( prob .lt. 0.441 + 0.205 ) then

                     idec = 2
                     kdcay = 3

                  else if( prob .lt. 0.441 + 0.205 + 0.301 ) then

                     idec = 3
                     kdcay = 3

                  else

                     idec = 4
                     kdcay = 2

                  end if

*-----------------------------------------------------------------------
*           (1) pi+, pi-, eta
*-----------------------------------------------------------------------

               if( idec .eq. 1 ) then

                     ip1 = 3
                     kp1 = 211

                     ip2 = 5
                     kp2 = -211

                     ip3 = 11
                     kp3 = 221

                  call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)

*-----------------------------------------------------------------------
*           (2) pi0, pi0, eta
*-----------------------------------------------------------------------

               else if( idec .eq. 2 ) then

                     ip1 = 4
                     kp1 = 111

                     ip2 = 4
                     kp2 = 111

                     ip3 = 11
                     kp3 = 221

                  call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)

*-----------------------------------------------------------------------
*           (3) pi+, pi-, photon
*-----------------------------------------------------------------------

               else if( idec .eq. 3 ) then

                     ip1 = 3
                     kp1 = 211

                     ip2 = 5
                     kp2 = -211

                     ip3 = 14
                     kp3 = 22

                  call dklos3(iprj,kprj,eein,ip1,kp1,ip2,kp2,ip3,kp3)

*-----------------------------------------------------------------------
*           (4) two photons
*-----------------------------------------------------------------------

               else

                     ip1 = 14
                     kp1 = 22

                     ip2 = 14
                     kp2 = 22

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

               end if

*-----------------------------------------------------------------------
*        Lambda0(3122) decay
*        decay channel : (1) 64.1 % : p, pi-
*                        (2) other  : n, pi0
*-----------------------------------------------------------------------

         else if( iprj .eq. 11 .and. kprj .eq. 3122 ) then

                     prob = rn(0)
                     kdcay = 2

                  if( prob .lt. 0.641 ) then

                     idec = 1

                  else

                     idec = 2

                  end if

*-----------------------------------------------------------------------
*           (1) p, pi-
*-----------------------------------------------------------------------

               if( idec .eq. 1 ) then

                     ip1 = 1
                     kp1 = 2212

                     ip2 = 5
                     kp2 = -211

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*           (2) n, pi0
*-----------------------------------------------------------------------

               else if( idec .eq. 2 ) then

                     ip1 = 2
                     kp1 = 2112

                     ip2 = 4
                     kp2 = 111

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

               end if

*-----------------------------------------------------------------------
*        Sigma+(3222) decay
*        decay channel : (1) 51.57 % : p, pi0
*                        (2) other   : n, pi+
*-----------------------------------------------------------------------

         else if( iprj .eq. 11 .and. kprj .eq. 3222 ) then

                     prob = rn(0)
                     kdcay = 2

                  if( prob .lt. 0.5157 ) then

                     idec = 1

                  else

                     idec = 2

                  end if

*-----------------------------------------------------------------------
*           (1) p, pi0
*-----------------------------------------------------------------------

               if( idec .eq. 1 ) then

                     ip1 = 1
                     kp1 = 2212

                     ip2 = 4
                     kp2 = 111

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*           (2) n, pi+
*-----------------------------------------------------------------------

               else if( idec .eq. 2 ) then

                     ip1 = 2
                     kp1 = 2112

                     ip2 = 3
                     kp2 = 211

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

               end if

*-----------------------------------------------------------------------
*        Sigma0(3212) decay
*        decay channel : (1) 100 % : Lambda0 + photon
*-----------------------------------------------------------------------

         else if( iprj .eq. 11 .and. kprj .eq. 3212 ) then

                     prob = rn(0)
                     kdcay = 2

                     ip1 = 11
                     kp1 = 3122

                     ip2 = 14
                     kp2 = 22

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*        Sigma-(3112) decay
*        decay channel : (1) 100 % : n + pi-
*-----------------------------------------------------------------------

         else if( iprj .eq. 11 .and. kprj .eq. 3112 ) then

                     prob = rn(0)
                     kdcay = 2

                     ip1 = 2
                     kp1 = 2112

                     ip2 = 5
                     kp2 = -211

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*        Xi0(3322) decay
*        decay channel : (1) 100 % : Lambda0 + pi0
*-----------------------------------------------------------------------

         else if( iprj .eq. 11 .and. kprj .eq. 3322 ) then

                     prob = rn(0)
                     kdcay = 2

                     ip1 = 11
                     kp1 = 3122

                     ip2 = 4
                     kp2 = 111

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*        Xi-(3312) decay
*        decay channel : (1) 100 % : Lambda0 + pi-
*-----------------------------------------------------------------------

         else if( iprj .eq. 11 .and. kprj .eq. 3312 ) then

                     prob = rn(0)
                     kdcay = 2

                     ip1 = 11
                     kp1 = 3122

                     ip2 = 5
                     kp2 = -211

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*        Omega-(3334) decay
*        decay channel : (1) 67.8 % : Lambda0 + K-
*                        (2) 23.6 % : Xi0 + pi-
*                        (2) other  : Xi- + pi0
*-----------------------------------------------------------------------

         else if( iprj .eq. 11 .and. kprj .eq. 3334 ) then

                     prob = rn(0)
                     kdcay = 2

                  if( prob .lt. 0.678 ) then

                     idec = 1

                  else if( prob .lt. 0.678 + 0.236 ) then

                     idec = 2

                  else

                     idec = 3

                  end if

*-----------------------------------------------------------------------
*           (1) Lambda0 + K-
*-----------------------------------------------------------------------

               if( idec .eq. 1 ) then

                     ip1 = 11
                     kp1 = 3122

                     ip2 = 10
                     kp2 = -321

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*           (2) Xi0 + pi-
*-----------------------------------------------------------------------

               else if( idec .eq. 2 ) then

                     ip1 = 11
                     kp1 = 3322

                     ip2 = 5
                     kp2 = -211

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

*-----------------------------------------------------------------------
*           (3) Xi- + pi0
*-----------------------------------------------------------------------

               else if( idec .eq. 3 ) then

                     ip1 = 11
                     kp1 = 3312

                     ip2 = 4
                     kp2 = 111

                  call dklos2(iprj,kprj,eein,ip1,kp1,ip2,kp2)

               end if

*-----------------------------------------------------------------------
*        other particles
*-----------------------------------------------------------------------

         else if( kprj ==     15) then; call decay_with_data(kprj, eein) !      tau- (1777.0)
         else if( kprj ==    -15) then; call decay_with_data(kprj, eein) !      tau+ (1777.0)
         else if( kprj ==    113) then; call decay_with_data(kprj, eein) !      rho0 (768.50)
         else if( kprj ==    213) then; call decay_with_data(kprj, eein) !      rho+ (766.90)
         else if( kprj ==   -213) then; call decay_with_data(kprj, eein) !      rho- (766.90)
         else if( kprj ==    223) then; call decay_with_data(kprj, eein) !      omega(781.94)
         else if( kprj ==  20213) then; call decay_with_data(kprj, eein) !      a_1+ (1230.0)
         else if( kprj == -20213) then; call decay_with_data(kprj, eein) !      a_1- (1230.0)
         else if( kprj ==    313) then; call decay_with_data(kprj, eein) !      K*0  (896.10)
         else if( kprj ==   -313) then; call decay_with_data(kprj, eein) ! anti-K*0  (896.10)
         else if( kprj ==    323) then; call decay_with_data(kprj, eein) !      K*+  (891.60)
         else if( kprj ==   -323) then; call decay_with_data(kprj, eein) !      K*-  (891.60)
         else if( kprj ==    333) then; call decay_with_data(kprj, eein) !      phi  (1019.4)
         else if( kprj ==    411) then; call decay_with_data(kprj, eein) !      D+   (1869.3)
         else if( kprj ==   -411) then; call decay_with_data(kprj, eein) !      D-   (1869.3)
         else if( kprj ==    421) then; call decay_with_data(kprj, eein) !      D0   (1864.5)
         else if( kprj ==   -421) then; call decay_with_data(kprj, eein) ! anti-D0   (1864.5)
         else if( kprj ==    431) then; call decay_with_data(kprj, eein) !      D_s+ (1968.5)
         else if( kprj ==   -431) then; call decay_with_data(kprj, eein) !      D_s- (1968.5)
         else if( kprj ==    511) then; call decay_with_data(kprj, eein) !      B0   (5279.2)
         else if( kprj ==   -511) then; call decay_with_data(kprj, eein) ! anti-B0   (5279.2)
         else if( kprj ==    521) then; call decay_with_data(kprj, eein) !      B+   (5278.9)
         else if( kprj ==   -521) then; call decay_with_data(kprj, eein) !      B-   (5278.9)
         else if( kprj ==    531) then; call decay_with_data(kprj, eein) !      B_s0 (5369.3)
         else if( kprj ==   -531) then; call decay_with_data(kprj, eein) ! anti-B_s0 (5369.3)
         else if( kprj ==    541) then; call decay_with_data(kprj, eein) !      B_c+ (6594.0)
         else if( kprj ==   -541) then; call decay_with_data(kprj, eein) !      B_c- (6594.0)

*-----------------------------------------------------------------------
*        no data
*-----------------------------------------------------------------------

         else

            if(udm_counter_dklos .eq. 0)
     &      write(6,'(''*** no data in dklos.f kf = '',i7)') kprj
            if(iudmodel .gt. 0) udm_counter_dklos=1

            return

         end if

*-----------------------------------------------------------------------
 1000 continue

*-----------------------------------------------------------------------
*     sumary decay mode
*-----------------------------------------------------------------------

            if( kdcay .eq. 0 ) return
!$OMP CRITICAL (dcayp_crit)
                  bdcayp(iprj) = bdcayp(iprj) + 1.0
                  adcayp(iprj) = adcayp(iprj) + oldwt
!$OMP END CRITICAL (dcayp_crit)

               if( iprj .eq. 11 ) then
!$OMP CRITICAL (dcp_crit)
                  if( nodcn .eq. 0 ) then

                        nodcn = nodcn + 1
                        nodcp(nodcn,1) = 1
                        nodcp(nodcn,2) = kprj

                     if( kdcay .eq. 2 ) then

                        nodcp(nodcn,3) = 2
                        nodcp(nodcn,4) = kp1
                        nodcp(nodcn,5) = kp2

                     else if( kdcay .eq. 3 ) then

                        nodcp(nodcn,3) = 3
                        nodcp(nodcn,4) = kp1
                        nodcp(nodcn,5) = kp2
                        nodcp(nodcn,6) = kp3

                     end if

                  else

                     do m = 1, nodcn

                        if( kdcay .eq. 2 ) then

                           if( kprj  .eq. nodcp(m,2) .and.
     &                         kdcay .eq. nodcp(m,3) .and.
     &                         kp1   .eq. nodcp(m,4) .and.
     &                         kp2   .eq. nodcp(m,5) ) then

                              nodcp(m,1) = nodcp(m,1) + 1
                              goto 310

                           end if

                        else if( kdcay .eq. 3 ) then

                           if( kprj  .eq. nodcp(m,2) .and.
     &                         kdcay .eq. nodcp(m,3) .and.
     &                         kp1   .eq. nodcp(m,4) .and.
     &                         kp2   .eq. nodcp(m,5) .and.
     &                         kp3   .eq. nodcp(m,6) ) then

                              nodcp(m,1) = nodcp(m,1) + 1
                              goto 310

                           end if

                        end if

                     end do

                        nodcn = nodcn + 1
                        nodcp(nodcn,1) = 1
                        nodcp(nodcn,2) = kprj

                     if( kdcay .eq. 2 ) then

                        nodcp(nodcn,3) = 2
                        nodcp(nodcn,4) = kp1
                        nodcp(nodcn,5) = kp2

                     else if( kdcay .eq. 3 ) then

                        nodcp(nodcn,3) = 3
                        nodcp(nodcn,4) = kp1
                        nodcp(nodcn,5) = kp2
                        nodcp(nodcn,6) = kp3

                     end if

  310                continue

                  end if
!$OMP END CRITICAL (dcp_crit)
               end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine c9decay(eein)
*                                                                      *
*                                                                      *
*       9C decay through beta decay to 9B                              *
*          data from E. Gete et.al. Phys. Rev. C61, 064310 (2000)      *                                                                      *
*                                                                      *
*                          total   8Be+p   5Li+alpha                   *
*          ground  of 9B:  49.1%   49.10%    0.00%                     *
*          2.34MeV of 9B:  36.8%    0.19%   36.60%                     *
*          2.80MeV of 9B:   7.0%    6.30%    0.73%                     *
*         12.16MeV of 9B:   7.1%    1.88%    5.00%                     *
*                   total 100.0    57.47%   42.53%                     *
*                                                                      *
*          57.5%  9B -> 8Be + p -> alpha + alpha + p                   *
*          42.5%  9B -> 5Li + alpha -> alpha + alpha + p               *
*                                                                      *
*       modified by K.Niita on 2006/10/17                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       eein    : kinetic energy of 9C (MeV)                           *
*                                                                      *
*     output :  in common                                              *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( pi = 3.1415926535898 )

*-----------------------------------------------------------------------

      include 'param00.inc'

*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)

*-----------------------------------------------------------------------

      dimension brach0(8), brach(8), bexc(8)

      data  brach0 /49.10,  0.19,  6.30,  1.88,
     &               0.00, 36.60,  0.70,  5.00/

      data  brach  /0.491,  0.4929,  0.5559,  0.5747,
     &              0.5747, 0.9407,  0.9480,  1.00/

      data  bexc   / 0.0,  2.34,  2.8,  12.16,
     &               0.0,  2.34,  2.8,  12.16/

*-----------------------------------------------------------------------
*        which channel ?
*-----------------------------------------------------------------------

            rms = rn(0)

         do i = 1, 8

            if( brach(i) .ge. rms ) goto 100

         end do

  100    ibc = i

            ex9b = bexc(ibc)

            ip0 = 19
            kp0 = 5000009

*-----------------------------------------------------------------------

         if( ibc .le. 4 ) then

            ip1 = 19
            id1 = 0
            kp1 = 4000008

            ip2 = 1
            id2 = 1
            kp2 = 2212

            ip3 = 18
            id3 = 18
            kp3 = 2000004

            ip4 = 18
            id4 = 18
            kp4 = 2000004

         else

            ip1 = 19
            id1 = 0
            kp1 = 3000005

            ip2 = 18
            id2 = 18
            kp2 = 2000004

            ip3 = 1
            id3 = 1
            kp3 = 2212

            ip4 = 18
            id4 = 18
            kp4 = 2000004

         end if

*-----------------------------------------------------------------------
*        initial values (MeV)
*-----------------------------------------------------------------------

               ein   = eein
               rmas0 = rmtyp(ip0,kp0) + ex9b
               rmas1 = rmtyp(ip1,kp1)
               rmas2 = rmtyp(ip2,kp2)
               rmas3 = rmtyp(ip3,kp3)
               rmas4 = rmtyp(ip4,kp4)

*-----------------------------------------------------------------------
*        outgoing momentum in cm system from kp0 -> kp1 + kp2
*-----------------------------------------------------------------------

               pabs  = sqrt( ( rmas0**2 - ( rmas1 + rmas2 )**2 )
     &                     * ( rmas0**2 - ( rmas1 - rmas2 )**2 ) )
     &               / 2.0 / rmas0

               etc1 = sqrt( rmas1**2 + pabs**2 )
               etc2 = sqrt( rmas2**2 + pabs**2 )

               cos1 = 1.0 - 2.0 * rn(0)
               sin1 = sqrt( 1.0 - cos1**2 )
               phi1 = 2.0 * pi * rn(0)

               pxc1 = pabs * sin1 * cos(phi1)
               pyc1 = pabs * sin1 * sin(phi1)
               pzc1 = pabs * cos1

*-----------------------------------------------------------------------
*        Lorentz transform to lab system
*-----------------------------------------------------------------------

               plab = sqrt( ein**2 + 2.0 * rmas0 * ein )
               elab = sqrt( plab**2 + rmas0**2 )

               betl = plab / elab
               gaml = elab / rmas0

               pxl1 =  pxc1
               pyl1 =  pyc1
               pzl1 =  pzc1 * gaml + etc1 * gaml * betl
               pal1 = sqrt( pxl1**2 + pyl1**2 + pzl1**2 )
               exl1 = sqrt( rmas1**2 + pal1**2 )

               pxl2 = -pxc1
               pyl2 = -pyc1
               pzl2 = -pzc1 * gaml + etc2 * gaml * betl
               pal2 = sqrt( pxl2**2 + pyl2**2 + pzl2**2 )
               exl2 = sqrt( rmas2**2 + pal2**2 )

*-----------------------------------------------------------------------
*        outgoing momentum in cm system from kp1 -> kp3 + kp4
*-----------------------------------------------------------------------

               pabs  = sqrt( ( rmas1**2 - ( rmas3 + rmas4 )**2 )
     &                     * ( rmas1**2 - ( rmas3 - rmas4 )**2 ) )
     &               / 2.0 / rmas1

               etc3 = sqrt( rmas3**2 + pabs**2 )
               etc4 = sqrt( rmas4**2 + pabs**2 )

               cos3 = 1.0 - 2.0 * rn(0)
               sin3 = sqrt( 1.0 - cos1**2 )
               phi3 = 2.0 * pi * rn(0)

               pxc3 = pabs * sin3 * cos(phi3)
               pyc3 = pabs * sin3 * sin(phi3)
               pzc3 = pabs * cos3

*-----------------------------------------------------------------------
*        Lorentz transform to lab system
*-----------------------------------------------------------------------

               plab = pal1
               elab = exl1

               betx = pxl1 / elab
               bety = pyl1 / elab
               betz = pzl1 / elab
               gaml = elab / rmas1

               p1beta = pxc3 * betx + pyc3 * bety + pzc3 * betz
               e3cm  = sqrt (rmas3**2 + pxc3**2 + pyc3**2 + pzc3**2)
               transf = gaml * ( gaml * p1beta / (gaml+1) - e3cm )

               pxl3 =  pxc3 + betx * transf
               pyl3 =  pyc3 + bety * transf
               pzl3 =  pzc3 + betz * transf

               pal3 = sqrt( pxl3**2 + pyl3**2 + pzl3**2 )
               exl3 = sqrt( rmas3**2 + pal3**2 )

               e4cm  = sqrt (rmas4**2 + pxc3**2 + pyc3**2 + pzc3**2)
               transf = gaml * ( -gaml * p1beta / (gaml+1) + e4cm )

               pxl4 = -pxc3 + betx * transf
               pyl4 = -pyc3 + bety * transf
               pzl4 = -pzc3 + betz * transf

               pal4 = sqrt( pxl4**2 + pyl4**2 + pzl4**2 )
               exl4 = sqrt( rmas4**2 + pal4**2 )

*-----------------------------------------------------------------------
*        booking of outgoing two particles
*-----------------------------------------------------------------------

               nclst = 3

               numpat(id2)  = numpat(id2) + 1
               rumpat(id2)  = rumpat(id2) + 1

               iclust(1)    = ipatf(ip2,kp2)

               jclust(0,1)  = 0
               jclust(1,1)  = ichgf(ip2,kp2)
               jclust(2,1)  = ibryf(ip2,kp2) - ichgf(ip2,kp2)
               jclust(3,1)  = id2
               jclust(4,1)  = 0
               jclust(5,1)  = ichgf(ip2,kp2)
               jclust(6,1)  = ibryf(ip2,kp2)
               jclust(7,1)  = kp2
               jclust(8,1)  = 0

               qclust(0,1)  = 0.0
               qclust(1,1)  = pxl2 / 1000.d0
               qclust(2,1)  = pyl2 / 1000.d0
               qclust(3,1)  = pzl2 / 1000.d0
               qclust(4,1)  = exl2 / 1000.d0
               qclust(5,1)  = rmas2 / 1000.d0
               qclust(6,1)  = 0.0
               qclust(7,1)  = ( exl2 - rmas2 )
               qclust(8,1)  = 1.0
               qclust(9,1)  = 0.0
               qclust(10,1) = 0.0d0
               qclust(11,1) = 0.0d0
               qclust(12,1) = 0.0d0

               numpat(id3)  = numpat(id3) + 1
               rumpat(id3)  = rumpat(id3) + 1

               iclust(2)    = ipatf(ip3,kp3)

               jclust(0,2)  = 0
               jclust(1,2)  = ichgf(ip3,kp3)
               jclust(2,2)  = ibryf(ip3,kp3) - ichgf(ip3,kp3)
               jclust(3,2)  = id3
               jclust(4,2)  = 0
               jclust(5,2)  = ichgf(ip3,kp3)
               jclust(6,2)  = ibryf(ip3,kp3)
               jclust(7,2)  = kp3
               jclust(8,2)  = 0

               qclust(0,2)  = 0.0
               qclust(1,2)  = pxl3 / 1000.d0
               qclust(2,2)  = pyl3 / 1000.d0
               qclust(3,2)  = pzl3 / 1000.d0
               qclust(4,2)  = exl3 / 1000.d0
               qclust(5,2)  = rmas3 / 1000.d0
               qclust(6,2)  = 0.0
               qclust(7,2)  = ( exl3 - rmas3 )
               qclust(8,2)  = 1.0
               qclust(9,2)  = 0.0
               qclust(10,2) = 0.0d0
               qclust(11,2) = 0.0d0
               qclust(12,2) = 0.0d0

               numpat(id4)  = numpat(id4) + 1
               rumpat(id4)  = rumpat(id4) + 1

               iclust(3)    = ipatf(ip4,kp4)

               jclust(0,3)  = 0
               jclust(1,3)  = ichgf(ip4,kp4)
               jclust(2,3)  = ibryf(ip4,kp4) - ichgf(ip4,kp4)
               jclust(3,3)  = id4
               jclust(4,3)  = 0
               jclust(5,3)  = ichgf(ip4,kp4)
               jclust(6,3)  = ibryf(ip4,kp4)
               jclust(7,3)  = kp4
               jclust(8,3)  = 0

               qclust(0,3)  = 0.0
               qclust(1,3)  = pxl4 / 1000.d0
               qclust(2,3)  = pyl4 / 1000.d0
               qclust(3,3)  = pzl4 / 1000.d0
               qclust(4,3)  = exl4 / 1000.d0
               qclust(5,3)  = rmas4 / 1000.d0
               qclust(6,3)  = 0.0
               qclust(7,3)  = ( exl4 - rmas4 )
               qclust(8,3)  = 1.0
               qclust(9,3)  = 0.0
               qclust(10,3) = 0.0d0
               qclust(11,3) = 0.0d0
               qclust(12,3) = 0.0d0

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine dklos2(ip0,kp0,eein,ip1,kp1,ip2,kp2)
*                                                                      *
*                                                                      *
*       two-body particle decay                                        *
*       modified by K.Niita on 01/05/2000                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       ip0     : particle type of projectile                          *
*       kp0     : kf code of projectile                                *
*       eein    : energy of projectile (MeV)                           *
*       ip1,2   : particle type of outgoing particles                  *
*       kp1,2   : kf code of outgoing prticles                         *
*                                                                      *
*     output :                                                         *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclst  = 2  : total number of out going particles and nuclei  *
*                                                                      *
*        iclust(nclst)                                                 *
*                                                                      *
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclust(i,nclst)                                               *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4,                                                *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        qclust(i,nclst)                                               *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, px (GeV/c)                                     *
*                  = 2, py (GeV/c)                                     *
*                  = 3, pz (GeV/c)                                     *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight change                                  *
*                  = 9, delay time                                     *
*                  = 10, x-displace                                    *
*                  = 11, y-displace                                    *
*                  = 12, z-displace                                    *
*                                                                      *
*        numpat(i) : total number of out going particles or nuclei     *
*                                                                      *
*                i =  0, nuclei                                        *
*                  =  1, proton                                        *
*                  =  2, neutron                                       *
*                  =  3, pi+                                           *
*                  =  4, pi0                                           *
*                  =  5, pi-                                           *
*                  =  6, mu+                                           *
*                  =  7, mu-                                           *
*                  =  8, K+                                            *
*                  =  9, K0                                            *
*                  = 10, K-                                            *
*                  = 11, other particles                               *
*                  = 14, gammma                                        *
*                  = 15, deuteron                                      *
*                  = 16, triton                                        *
*                  = 17, 3He                                           *
*                  = 18, Alpha                                         *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( pi = 3.1415926535898 )

*-----------------------------------------------------------------------

      include 'param00.inc'

*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)

*-----------------------------------------------------------------------
*        initial values
*-----------------------------------------------------------------------

               ein   = eein / 1000.0
               rmas0 = rmtyp(ip0,kp0) / 1000.0
               rmas1 = rmtyp(ip1,kp1) / 1000.0
               rmas2 = rmtyp(ip2,kp2) / 1000.0

*-----------------------------------------------------------------------
*        outgoing momentum in cm system
*-----------------------------------------------------------------------

               pabs  = sqrt( ( rmas0**2 - ( rmas1 + rmas2 )**2 )
     &                     * ( rmas0**2 - ( rmas1 - rmas2 )**2 ) )
     &               / 2.0 / rmas0

               etc1 = sqrt( rmas1**2 + pabs**2 )
               etc2 = sqrt( rmas2**2 + pabs**2 )

               cos1 = 1.0 - 2.0 * rn(0)
               sin1 = sqrt( 1.0 - cos1**2 )
               phi1 = 2.0 * pi * rn(0)

               pxc1 = pabs * sin1 * cos(phi1)
               pyc1 = pabs * sin1 * sin(phi1)
               pzc1 = pabs * cos1

*-----------------------------------------------------------------------
*        Lorentz transform to lab system
*-----------------------------------------------------------------------

               plab = sqrt( ein**2 + 2.0 * rmas0 * ein )
               elab = sqrt( plab**2 + rmas0**2 )

               betl = plab / elab
               gaml = elab / rmas0

               pxl1 =  pxc1
               pyl1 =  pyc1
               pzl1 =  pzc1 * gaml + etc1 * gaml * betl
               pal1 = sqrt( pxl1**2 + pyl1**2 + pzl1**2 )
               exl1 = sqrt( rmas1**2 + pal1**2 )

               pxl2 = -pxc1
               pyl2 = -pyc1
               pzl2 = -pzc1 * gaml + etc2 * gaml * betl
               pal2 = sqrt( pxl2**2 + pyl2**2 + pzl2**2 )
               exl2 = sqrt( rmas2**2 + pal2**2 )

*-----------------------------------------------------------------------
*        booking of outgoing two particles
*-----------------------------------------------------------------------

               nclst = 2

               numpat(ip1)  = numpat(ip1) + 1
               rumpat(ip1)  = rumpat(ip1) + 1

               iclust(1)    = ipatf(ip1,kp1)

               jclust(0,1)  = 0
               jclust(1,1)  = 0
               jclust(2,1)  = 0
               jclust(3,1)  = ip1
               jclust(4,1)  = 0
               jclust(5,1)  = ichgf(ip1,kp1)
               jclust(6,1)  = ibryf(ip1,kp1)
               jclust(7,1)  = kp1
               jclust(8,1)  = 0

               qclust(0,1)  = 0.0
               qclust(1,1)  = pxl1
               qclust(2,1)  = pyl1
               qclust(3,1)  = pzl1
               qclust(4,1)  = exl1
               qclust(5,1)  = rmas1
               qclust(6,1)  = 0.0
               qclust(7,1)  = ( exl1 - rmas1 ) * 1000.
               qclust(8,1)  = 1.0
               qclust(9,1)  = 0.0
               qclust(10,1) = 0.0d0
               qclust(11,1) = 0.0d0
               qclust(12,1) = 0.0d0

               numpat(ip2)  = numpat(ip2) + 1
               rumpat(ip2)  = rumpat(ip2) + 1

               iclust(2)    = ipatf(ip2,kp2)

               jclust(0,2)  = 0
               jclust(1,2)  = 0
               jclust(2,2)  = 0
               jclust(3,2)  = ip2
               jclust(4,2)  = 0
               jclust(5,2)  = ichgf(ip2,kp2)
               jclust(6,2)  = ibryf(ip2,kp2)
               jclust(7,2)  = kp2
               jclust(8,2)  = 0

               qclust(0,2)  = 0.0
               qclust(1,2)  = pxl2
               qclust(2,2)  = pyl2
               qclust(3,2)  = pzl2
               qclust(4,2)  = exl2
               qclust(5,2)  = rmas2
               qclust(6,2)  = 0.0
               qclust(7,2)  = ( exl2 - rmas2 ) * 1000.
               qclust(8,2)  = 1.0
               qclust(9,2)  = 0.0
               qclust(10,2) = 0.0d0
               qclust(11,2) = 0.0d0
               qclust(12,2) = 0.0d0

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine dklos3(ip0,kp0,eein,ip1,kp1,ip2,kp2,ip3,kp3)
*                                                                      *
*                                                                      *
*       three-body particle decay                                      *
*       modified by K.Niita on 01/05/2000                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       ip0     : particle type of projectile                          *
*       kp0     : kf code of projectile                                *
*       eein    : energy of projectile (MeV)                           *
*       ip1,2,3 : particle type of outgoing particles                  *
*       kp1,2,3 : kf code of outgoing prticles                         *
*                                                                      *
*     output :                                                         *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclst  = 3  : total number of out going particles and nuclei  *
*                                                                      *
*        iclust(nclst)                                                 *
*                                                                      *
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclust(i,nclst)                                               *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4,                                                *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        qclust(i,nclst)                                               *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, px (GeV/c)                                     *
*                  = 2, py (GeV/c)                                     *
*                  = 3, pz (GeV/c)                                     *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight change                                  *
*                  = 9, delay time                                     *
*                  = 10, x-displace                                    *
*                  = 11, y-displace                                    *
*                  = 12, z-displace                                    *
*                                                                      *
*        numpat(i) : total number of out going particles or nuclei     *
*                                                                      *
*                i =  0, nuclei                                        *
*                  =  1, proton                                        *
*                  =  2, neutron                                       *
*                  =  3, pi+                                           *
*                  =  4, pi0                                           *
*                  =  5, pi-                                           *
*                  =  6, mu+                                           *
*                  =  7, mu-                                           *
*                  =  8, K+                                            *
*                  =  9, K0                                            *
*                  = 10, K-                                            *
*                  = 11, other particles                               *
*                  = 14, gammma                                        *
*                  = 15, deuteron                                      *
*                  = 16, triton                                        *
*                  = 17, 3He                                           *
*                  = 18, Alpha                                         *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( pi = 3.1415926535898 )

*-----------------------------------------------------------------------

      include 'param00.inc'

*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)

*-----------------------------------------------------------------------
*        initial values
*-----------------------------------------------------------------------

               ein   = eein / 1000.0
               rmas0 = rmtyp(ip0,kp0) / 1000.0
               rmas1 = rmtyp(ip1,kp1) / 1000.0
               rmas2 = rmtyp(ip2,kp2) / 1000.0
               rmas3 = rmtyp(ip3,kp3) / 1000.0

*-----------------------------------------------------------------------
*        determin rm12 and rm13
*-----------------------------------------------------------------------

               srm12 = rmas1 + rmas2
               drm12 = rmas0 - rmas3

               srm13 = rmas1 + rmas3
               drm13 = rmas0 - rmas2

  100       continue

               x1 = rn(0)
               x2 = rn(0)

               rm12 = srm12**2 + ( drm12**2 - srm12**2 ) * x1
               rm13 = srm13**2 + ( drm13**2 - srm13**2 ) * x2

               e1s = ( rm12 + rmas1**2 - rmas2**2 ) / 2.0 / sqrt(rm12)
               e3s = ( rmas0**2 - rm12 - rmas3**2 ) / 2.0 / sqrt(rm12)

               sr13 = ( e1s + e3s )**2
     &              - ( sqrt( e1s**2 - rmas1**2 )
     &                + sqrt( e3s**2 - rmas3**2 ) )**2

               dr13 = ( e1s + e3s )**2
     &              - ( sqrt( e1s**2 - rmas1**2 )
     &                - sqrt( e3s**2 - rmas3**2 ) )**2

               if( rm13 .ge. dr13 .or. rm13 .le. sr13 ) goto 100

                  rm12 = sqrt(rm12)
                  rm13 = sqrt(rm13)

*-----------------------------------------------------------------------
*        outgoing momentum of particle 3 in cm system
*-----------------------------------------------------------------------

               pab3  = sqrt( ( rmas0**2 - ( rm12 + rmas3 )**2 )
     &                     * ( rmas0**2 - ( rm12 - rmas3 )**2 ) )
     &               / 2.0 / rmas0

               etc3 = sqrt( rmas3**2 + pab3**2 )

               cos3 = 1.0 - 2.0 * rn(0)
               sin3 = sqrt( 1.0 - cos3**2 )
               phi3 = 2.0 * pi * rn(0)

               pxc3 = pab3 * sin3 * cos(phi3)
               pyc3 = pab3 * sin3 * sin(phi3)
               pzc3 = pab3 * cos3

*-----------------------------------------------------------------------
*        outgoing momentum of particle 1, 2 in 1-2 cm system
*-----------------------------------------------------------------------

               pac1  = sqrt( ( rm12**2 - ( rmas1 + rmas2 )**2 )
     &                     * ( rm12**2 - ( rmas1 - rmas2 )**2 ) )
     &               / 2.0 / rm12

               ecc1 = sqrt( rmas1**2 + pac1**2 )
               ecc2 = sqrt( rmas2**2 + pac1**2 )

               cos1 = 1.0 - 2.0 * rn(0)
               sin1 = sqrt( 1.0 - cos1**2 )
               phi1 = 2.0 * pi * rn(0)

               pxs1 = pac1 * sin1 * cos(phi1)
               pys1 = pac1 * sin1 * sin(phi1)
               pzs1 = pac1 * cos1

*-----------------------------------------------------------------------
*        Lorentz transform to cm system from 1-2 cm frame
*-----------------------------------------------------------------------

               e12 = sqrt( rm12**2 + pxc3**2 + pyc3**2 + pzc3**2 )

               betx = pxc3 / e12
               bety = pyc3 / e12
               betz = pzc3 / e12
               bet2 = betx**2 + bety**2 + betz**2
               gamm = e12 / rm12

               p1bet = pxs1 * betx + pys1 * bety + pzs1 * betz
               tran1 = gamm * (  gamm * p1bet / ( gamm + 1.0 ) - ecc1 )
               tran2 = gamm * ( -gamm * p1bet / ( gamm + 1.0 ) - ecc2 )

               pxc1 = betx * tran1 + pxs1
               pyc1 = bety * tran1 + pys1
               pzc1 = betz * tran1 + pzs1

               etc1 = sqrt( pxc1**2 + pyc1**2 + pzc1**2 + rmas1**2 )

               pxc2 = betx * tran2 - pxs1
               pyc2 = bety * tran2 - pys1
               pzc2 = betz * tran2 - pzs1

               etc2 = sqrt( pxc2**2 + pyc2**2 + pzc2**2 + rmas2**2 )

*-----------------------------------------------------------------------
*        Lorentz transform to lab system
*-----------------------------------------------------------------------

               plab = sqrt( ein**2 + 2.0 * rmas0 * ein )
               elab = sqrt( plab**2 + rmas0**2 )

               betl = plab / elab
               gaml = elab / rmas0

               pxl1 = pxc1
               pyl1 = pyc1
               pzl1 = pzc1 * gaml + etc1 * gaml * betl
               pal1 = sqrt( pxl1**2 + pyl1**2 + pzl1**2 )
               exl1 = sqrt( rmas1**2 + pal1**2 )

               pxl2 = pxc2
               pyl2 = pyc2
               pzl2 = pzc2 * gaml + etc2 * gaml * betl
               pal2 = sqrt( pxl2**2 + pyl2**2 + pzl2**2 )
               exl2 = sqrt( rmas2**2 + pal2**2 )

               pxl3 = pxc3
               pyl3 = pyc3
               pzl3 = pzc3 * gaml + etc3 * gaml * betl
               pal3 = sqrt( pxl3**2 + pyl3**2 + pzl3**2 )
               exl3 = sqrt( rmas3**2 + pal3**2 )

*-----------------------------------------------------------------------
*        booking of outgoing two particles
*-----------------------------------------------------------------------

               nclst = 3

               numpat(ip1)  = numpat(ip1) + 1
               rumpat(ip1)  = rumpat(ip1) + 1

               iclust(1)    = ipatf(ip1,kp1)
               jclust(0,1)  = 0
               jclust(1,1)  = 0
               jclust(2,1)  = 0
               jclust(3,1)  = ip1
               jclust(4,1)  = 0
               jclust(5,1)  = ichgf(ip1,kp1)
               jclust(6,1)  = ibryf(ip1,kp1)
               jclust(7,1)  = kp1
               jclust(8,1)  = 0

               qclust(0,1)  = 0.0
               qclust(1,1)  = pxl1
               qclust(2,1)  = pyl1
               qclust(3,1)  = pzl1
               qclust(4,1)  = exl1
               qclust(5,1)  = rmas1
               qclust(6,1)  = 0.0
               qclust(7,1)  = ( exl1 - rmas1 ) * 1000.
               qclust(8,1)  = 1.0
               qclust(9,1)  = 0.0
               qclust(10,1) = 0.0d0
               qclust(11,1) = 0.0d0
               qclust(12,1) = 0.0d0

               numpat(ip2)  = numpat(ip2) + 1
               rumpat(ip2)  = rumpat(ip2) + 1

               iclust(2)    = ipatf(ip2,kp2)

               jclust(0,2)  = 0
               jclust(1,2)  = 0
               jclust(2,2)  = 0
               jclust(3,2)  = ip2
               jclust(4,2)  = 0
               jclust(5,2)  = ichgf(ip2,kp2)
               jclust(6,2)  = ibryf(ip2,kp2)
               jclust(7,2)  = kp2
               jclust(8,2)  = 0

               qclust(0,2)  = 0.0
               qclust(1,2)  = pxl2
               qclust(2,2)  = pyl2
               qclust(3,2)  = pzl2
               qclust(4,2)  = exl2
               qclust(5,2)  = rmas2
               qclust(6,2)  = 0.0
               qclust(7,2)  = ( exl2 - rmas2 ) * 1000.
               qclust(8,2)  = 1.0
               qclust(9,2)  = 0.0
               qclust(10,2) = 0.0d0
               qclust(11,2) = 0.0d0
               qclust(12,2) = 0.0d0

               numpat(ip3)  = numpat(ip3) + 1
               rumpat(ip3)  = rumpat(ip3) + 1

               iclust(3)    = ipatf(ip3,kp3)

               jclust(0,3)  = 0
               jclust(1,3)  = 0
               jclust(2,3)  = 0
               jclust(3,3)  = ip3
               jclust(4,3)  = 0
               jclust(5,3)  = ichgf(ip3,kp3)
               jclust(6,3)  = ibryf(ip3,kp3)
               jclust(7,3)  = kp3
               jclust(8,3)  = 0

               qclust(0,3)  = 0.0
               qclust(1,3)  = pxl3
               qclust(2,3)  = pyl3
               qclust(3,3)  = pzl3
               qclust(4,3)  = exl3
               qclust(5,3)  = rmas3
               qclust(6,3)  = 0.0
               qclust(7,3)  = ( exl3 - rmas3 ) * 1000.
               qclust(8,3)  = 1.0
               qclust(9,3)  = 0.0
               qclust(10,3) = 0.0d0
               qclust(11,3) = 0.0d0
               qclust(12,3) = 0.0d0

*-----------------------------------------------------------------------

      return
      end

