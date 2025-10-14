************************************************************************
*                                                                      *
      module neutrino_mod
*                                                                      *
*                                                                      *
*        Last Revised:     2017 11 10  by T.Ogawa                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              Calculate neutrino reactions                            *
*                                                                      *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)
      integer ntrnore

      real*8, dimension(136), private:: sigden, sigdpp, sigdnn,
     & sigdnupn, sigdnbpn
      real*8, dimension(45), private:: sigpen, sigpn, enelm, costelm

      real*8, private::
     & sig2H_nu_e_b_cc(2,0:122,79,37),
     & sig2H_nu_e_cc  (2,0:122,79,37),
     & sig2H_nu_nc    (2,0:218,72,37),
     & sig2H_nu_bar_nc(2,0:218,72,37),
     & eng2H_nu(79)

      data eng2H_nu/
     &   1.5d0, 1.6d0, 1.8d0, 2.0d0, 2.2d0, 2.4d0, 2.6d0, 2.8d0,
     &   3.0d0, 3.2d0, 3.4d0, 3.6d0, 3.8d0, 4.0d0, 4.2d0, 4.4d0, 4.6d0,
     &   4.8d0, 5.0d0, 5.2d0, 5.4d0, 5.6d0, 5.8d0, 6.0d0, 6.2d0, 6.4d0,
     &   6.6d0, 6.8d0, 7.0d0, 7.2d0, 7.4d0, 7.6d0, 7.8d0, 8.0d0, 8.2d0,
     &   8.4d0, 8.6d0, 8.8d0, 9.0d0, 9.2d0, 9.4d0, 9.6d0, 9.8d0, 10.0d0,
     &  10.5d0, 11.0d0, 11.5d0, 12.0d0, 12.5d0, 13.0d0, 13.5d0, 14.0d0,
     &  14.5d0, 15.0d0, 15.5d0, 16.0d0, 16.5d0, 17.0d0, 17.5d0, 18.0d0,
     &  18.5d0, 19.0d0, 19.5d0,  20.d0,  30.d0,  40.d0,  50.d0,  60.d0,
     &   70.d0,  80.d0,  90.d0, 100.d0, 110.d0, 120.d0, 130.d0, 140.d0,
     &  150.d0, 160.d0, 170.d0/

      integer, private:: iover, ideux_read
      data iover /0/

      data ideux_read /0/

c S. NAKAMURA, T. SATO, V. GUDKOV, AND K. KUBODERA, PHYSICAL REVIEW C, 63, 0346 (2001)
c d + nu reaction energy bin (Unit : MeV )
      data sigden /
     &  1.4d0,  1.6d0,  1.8d0,  2.0d0,  2.2d0,  2.4d0,  2.6d0,  2.8d0,
     &  3.0d0,  3.2d0,  3.4d0,  3.6d0,  3.8d0,  4.0d0,  4.2d0,  4.4d0,
     &  4.6d0,  4.8d0,  5.0d0,  5.2d0,  5.4d0,  5.6d0,  5.8d0,  6.0d0,
     &  6.2d0,  6.4d0,  6.6d0,  6.8d0,  7.0d0,  7.2d0,  7.4d0,  7.6d0,
     &  7.8d0,  8.0d0,  8.2d0,  8.4d0,  8.6d0,  8.8d0,  9.0d0,  9.2d0,
     &  9.4d0,  9.6d0,  9.8d0, 10.0d0, 10.2d0, 10.4d0, 10.6d0, 10.8d0,
     & 11.0d0, 11.2d0, 11.4d0, 11.6d0, 11.8d0, 12.0d0, 12.2d0, 12.4d0,
     & 12.6d0, 12.8d0, 13.0d0, 13.5d0, 14.0d0, 14.5d0, 15.0d0, 15.5d0,
     & 16.0d0, 16.5d0, 17.0d0, 17.5d0, 18.0d0, 18.5d0, 19.0d0, 19.5d0,
     & 20.0d0, 20.5d0, 21.0d0, 21.5d0, 22.0d0, 22.5d0, 23.0d0, 23.5d0,
     & 24.0d0, 24.5d0,  25.d0,  26.d0,  27.d0,  28.d0,  29.d0,  30.d0,
     &  31.d0,  32.d0,  33.d0,  34.d0,  35.d0,  36.d0,  37.d0,  38.d0,
     &  39.d0,  40.d0,  41.d0,  42.d0,  43.d0,  44.d0,  45.d0,  46.d0,
     &  47.d0,  48.d0,  49.d0,  50.d0,  51.d0,  52.d0,  53.d0,  54.d0,
     &  55.d0,  60.d0,  65.d0,  70.d0,  75.d0,  80.d0,  85.d0,  90.d0,
     &  95.d0, 100.d0, 105.d0, 110.d0, 115.d0, 120.d0, 125.d0, 130.d0,
     & 135.d0, 140.d0, 145.d0, 150.d0, 155.d0, 160.d0, 165.d0, 170.d0/

c d + nu_e -> e- + p + p Cross section (Unit : 1e-40 cm^2 = 1e-16 b )
      data sigdpp /
     & 4.680d-8, 1.147d-6, 1.147d-5, 3.603d-5, 7.833d-5, 1.404d-4,
     & 2.242d-4, 3.315d-4, 4.639d-4, 6.228d-4, 8.095d-4, 1.025d-3,
     & 1.271d-3, 1.547d-3, 1.855d-3, 2.196d-3, 2.570d-3, 2.978d-3,
     & 3.420d-3, 3.897d-3, 4.410d-3, 4.959d-3, 5.544d-3, 6.166d-3,
     & 6.825d-3, 7.522d-3, 8.258d-3, 9.031d-3, 9.843d-3, 1.069d-2,
     & 1.159d-2, 1.252d-2, 1.349d-2, 1.450d-2, 1.555d-2, 1.664d-2,
     & 1.777d-2, 1.894d-2, 2.016d-2, 2.141d-2, 2.271d-2, 2.405d-2,
     & 2.544d-2, 2.686d-2, 2.833d-2, 2.984d-2, 3.139d-2, 3.299d-2,
     & 3.463d-2, 3.631d-2, 3.804d-2, 3.981d-2, 4.163d-2, 4.349d-2,
     & 4.539d-2, 4.734d-2, 4.933d-2, 5.137d-2, 5.346d-2, 5.887d-2,
     & 6.456d-2, 7.054d-2, 7.681d-2, 8.338d-2, 9.024d-2, 9.740d-2,
     & 1.049d-1, 1.126d-1, 1.207d-1, 1.291d-1, 1.378d-1, 1.468d-1,
     & 1.561d-1, 1.657d-1, 1.757d-1, 1.859d-1, 1.965d-1, 2.074d-1,
     & 2.187d-1, 2.303d-1, 2.422d-1, 2.545d-1, 2.671d-1, 2.933d-1,
     & 3.209d-1, 3.499d-1, 3.803d-1, 4.121d-1, 4.454d-1, 4.802d-1,
     & 5.164d-1, 5.541d-1, 5.934d-1, 6.342d-1, 6.765d-1, 7.204d-1,
     & 7.659d-1, 8.130d-1, 8.617d-1, 9.120d-1, 9.639d-1,  1.018d0,
     &  1.073d0,  1.130d0,  1.188d0,  1.248d0,  1.310d0,  1.374d0,
     &  1.440d0,  1.507d0,  1.575d0,  1.646d0,  1.718d0,  2.107d0,
     &  2.540d0,  3.018d0,  3.540d0,  4.108d0,  4.721d0,  5.378d0,
     &  6.079d0,  6.824d0,  7.612d0,  8.440d0,  9.307d0, 1.021d+1,
     & 1.116d+1, 1.214d+1, 1.315d+1, 1.420d+1, 1.528d+1, 1.639d+1,
     & 1.753d+1, 1.870d+1, 1.989d+1, 2.111d+1/

c d + nu_e_bar -> e+ + n + n Cross section (Unit : 1e-40 cm^2 = 1e-16 b )
      data sigdnn /
     &  0.000d0,  0.000d0,  0.000d0,  0.000d0,  0.000d0,  0.000d0,
     &  0.000d0,  0.000d0,  0.000d0,  0.000d0,  0.000d0,  0.000d0,
     &  0.000d0,  0.000d0, 1.115d-5, 4.554d-5, 1.010d-4, 1.787d-4,
     & 2.799d-4, 4.059d-4, 5.578d-4, 7.364d-4, 9.427d-4, 1.177d-3,
     & 1.441d-3, 1.733d-3, 2.056d-3, 2.409d-3, 2.792d-3, 3.206d-3,
     & 3.652d-3, 4.127d-3, 4.635d-3, 5.175d-3, 5.746d-3, 6.349d-3,
     & 6.984d-3, 7.652d-3, 8.351d-3, 9.082d-3, 9.846d-3, 1.064d-2,
     & 1.147d-2, 1.233d-2, 1.322d-2, 1.415d-2, 1.510d-2, 1.609d-2,
     & 1.712d-2, 1.817d-2, 1.925d-2, 2.037d-2, 2.152d-2, 2.270d-2,
     & 2.392d-2, 2.516d-2, 2.644d-2, 2.775d-2, 2.909d-2, 3.258d-2,
     & 3.626d-2, 4.015d-2, 4.422d-2, 4.849d-2, 5.295d-2, 5.760d-2,
     & 6.244d-2, 6.747d-2, 7.268d-2, 7.809d-2, 8.367d-2, 8.944d-2,
     & 9.539d-2, 1.015d-1, 1.078d-1, 1.143d-1, 1.210d-1, 1.278d-1,
     & 1.348d-1, 1.420d-1, 1.494d-1, 1.569d-1, 1.646d-1, 1.805d-1,
     & 1.971d-1, 2.143d-1, 2.322d-1, 2.507d-1, 2.698d-1, 2.896d-1,
     & 3.099d-1, 3.309d-1, 3.525d-1, 3.746d-1, 3.973d-1, 4.206d-1,
     & 4.445d-1, 4.689d-1, 4.938d-1, 5.193d-1, 5.453d-1, 5.718d-1,
     & 5.988d-1, 6.264d-1, 6.544d-1, 6.829d-1, 7.119d-1, 7.413d-1,
     & 7.712d-1, 8.016d-1, 8.324d-1, 8.636d-1, 8.953d-1,  1.060d0,
     &  1.233d0,  1.415d0,  1.606d0,  1.802d0,  2.004d0,  2.212d0,
     &  2.424d0,  2.640d0,  2.859d0,  3.081d0,  3.306d0,  3.532d0,
     &  3.760d0,  3.990d0,  4.220d0,  4.452d0,  4.684d0,  4.918d0,
     &  5.151d0,  5.385d0,  5.621d0,  5.856d0/

c d + nu_* -> p + n + nu_* Cross section (Unit : 1e-40 cm^2 = 1e-16 b )
      data sigdnupn /
     &   0.000d0,  0.000d0,  0.000d0,  0.000d0,  0.000d0, 4.279d-7,
     &  4.258d-6, 1.457d-5, 3.355d-5, 6.286d-5, 1.038d-4, 1.574d-4,
     &  2.246d-4, 3.060d-4, 4.024d-4, 5.142d-4, 6.420d-4, 7.860d-4,
     &  9.468d-4, 1.125d-3, 1.320d-3, 1.533d-3, 1.763d-3, 2.012d-3,
     &  2.279d-3, 2.564d-3, 2.868d-3, 3.191d-3, 3.532d-3, 3.893d-3,
     &  4.273d-3, 4.672d-3, 5.091d-3, 5.529d-3, 5.987d-3, 6.464d-3,
     &  6.961d-3, 7.479d-3, 8.016d-3, 8.573d-3, 9.150d-3, 9.747d-3,
     &  1.036d-2, 1.100d-2, 1.166d-2, 1.234d-2, 1.304d-2, 1.376d-2,
     &  1.450d-2, 1.526d-2, 1.604d-2, 1.684d-2, 1.767d-2, 1.851d-2,
     &  1.938d-2, 2.026d-2, 2.117d-2, 2.210d-2, 2.305d-2, 2.551d-2,
     &  2.811d-2, 3.084d-2, 3.371d-2, 3.671d-2, 3.984d-2, 4.311d-2,
     &  4.651d-2, 5.006d-2, 5.374d-2, 5.755d-2, 6.151d-2, 6.560d-2,
     &  6.984d-2, 7.421d-2, 7.872d-2, 8.338d-2, 8.817d-2, 9.311d-2,
     &  9.819d-2, 1.034d-1, 1.088d-1, 1.143d-1, 1.199d-1, 1.317d-1,
     &  1.440d-1, 1.569d-1, 1.704d-1, 1.845d-1, 1.992d-1, 2.145d-1,
     &  2.304d-1, 2.469d-1, 2.640d-1, 2.817d-1, 3.001d-1, 3.190d-1,
     &  3.386d-1, 3.588d-1, 3.796d-1, 4.011d-1, 4.232d-1, 4.459d-1,
     &  4.692d-1, 4.932d-1, 5.178d-1, 5.430d-1, 5.689d-1, 5.954d-1,
     &  6.226d-1, 6.504d-1, 6.788d-1, 7.079d-1, 7.376d-1, 8.957d-1,
     &   1.070d0,  1.260d0,  1.465d0,  1.686d0,  1.922d0,  2.172d0,
     &   2.437d0,  2.715d0,  3.007d0,  3.313d0,  3.630d0,  3.958d0,
     &   4.298d0,  4.648d0,  5.009d0,  5.378d0,  5.756d0,  6.143d0,
     &   6.539d0,  6.941d0,  7.350d0,  7.765d0/


c d + nu_*_bar -> p + n + nu_*_bar Cross section (Unit : 1e-40 cm^2 = 1e-16 b )
      data sigdnbpn /
     &  0.000d0,  0.000d0,  0.000d0,  0.000d0,  0.000d0, 4.248d-7,
     & 4.222d-6, 1.443d-5, 3.320d-5, 6.213d-5, 1.025d-4, 1.553d-4,
     & 2.213d-4, 3.012d-4, 3.956d-4, 5.049d-4, 6.297d-4, 7.702d-4,
     & 9.267d-4, 1.100d-3, 1.289d-3, 1.495d-3, 1.718d-3, 1.958d-3,
     & 2.215d-3, 2.490d-3, 2.782d-3, 3.092d-3, 3.419d-3, 3.764d-3,
     & 4.126d-3, 4.506d-3, 4.904d-3, 5.320d-3, 5.754d-3, 6.206d-3,
     & 6.676d-3, 7.163d-3, 7.669d-3, 8.193d-3, 8.735d-3, 9.294d-3,
     & 9.872d-3, 1.047d-2, 1.108d-2, 1.171d-2, 1.236d-2, 1.303d-2,
     & 1.372d-2, 1.442d-2, 1.514d-2, 1.588d-2, 1.664d-2, 1.741d-2,
     & 1.821d-2, 1.902d-2, 1.985d-2, 2.069d-2, 2.156d-2, 2.379d-2,
     & 2.614d-2, 2.860d-2, 3.117d-2, 3.385d-2, 3.663d-2, 3.953d-2,
     & 4.253d-2, 4.564d-2, 4.886d-2, 5.218d-2, 5.561d-2, 5.915d-2,
     & 6.279d-2, 6.653d-2, 7.038d-2, 7.434d-2, 7.839d-2, 8.255d-2,
     & 8.681d-2, 9.117d-2, 9.564d-2, 1.002d-1, 1.049d-1, 1.145d-1,
     & 1.245d-1, 1.350d-1, 1.458d-1, 1.570d-1, 1.685d-1, 1.805d-1,
     & 1.928d-1, 2.055d-1, 2.186d-1, 2.320d-1, 2.458d-1, 2.600d-1,
     & 2.745d-1, 2.893d-1, 3.045d-1, 3.200d-1, 3.359d-1, 3.521d-1,
     & 3.686d-1, 3.854d-1, 4.026d-1, 4.201d-1, 4.379d-1, 4.559d-1,
     & 4.743d-1, 4.930d-1, 5.120d-1, 5.313d-1, 5.509d-1, 6.528d-1,
     & 7.612d-1, 8.757d-1, 9.959d-1,  1.121d0,  1.250d0,  1.383d0,
     &  1.520d0,  1.660d0,  1.803d0,  1.949d0,  2.097d0,  2.247d0,
     &  2.397d0,  2.549d0,  2.702d0,  2.855d0,  3.009d0,  3.163d0,
     &  3.318d0,  3.472d0,  3.627d0,  3.781d0/

************************************************************************
c d + nu_e -> e- + p + p   e+ double-differential distribution

************************************************************************

c p + nu_e-bar -> n + e+  Energy bin (Unit : MeV )
      data sigpen /
     & 1.806d0, 2.01d0, 2.25d0, 2.51d0, 2.80d0, 3.12d0, 3.48d0, 3.89d0,
     &  4.33d0, 4.84d0, 5.40d0, 6.02d0, 6.72d0, 7.49d0, 8.36d0, 8.83d0,
     &  9.85d0, 11.0d0, 12.3d0, 13.7d0, 15.3d0, 17.0d0, 19.0d0, 21.2d0,
     &  23.6d0, 26.4d0, 29.4d0, 32.8d0, 36.6d0, 40.9d0, 43.2d0, 48.2d0,
     &  53.7d0, 59.9d0, 66.9d0, 74.6d0, 83.2d0, 92.9d0, 104.d0, 116.d0,
     &  129.d0, 144.d0, 160.d0, 179.d0, 200.d0/

c p + nu_e-bar -> n + e+  Cross section (Unit : 1e-40 cm^2 = 1e-16 b )
      data sigpn /
     &  0.00d0, 3.51d-4, 7.35d-4, 1.27d-3, 2.02d-3, 3.04d-3, 4.40d-3,
     & 6.19d-3, 8.54d-3, 1.16d-2, 1.55d-2, 2.05d-2, 2.69d-2, 3.49d-2,
     & 4.51d-2, 5.11d-2, 6.54d-2, 8.32d-2, 1.05d-1, 1.33d-1, 1.67d-1,
     & 2.09d-1, 2.61d-1, 3.24d-1, 4.01d-1, 4.95d-1, 6.08d-1, 7.44d-1,
     & 9.08d-1,  1.10d0,  1.21d0,  1.47d0,  1.76d0,  2.10d0,  2.50d0,
     &  2.96d0,  3.48d0,  4.07d0,  4.73d0,  5.46d0,  6.27d0,  7.15d0,
     &  8.10d0,  9.13d0,  10.2d0/

c p + nu_e-bar -> n + e+  mean electron energy
      data enelm /0.00d0,
     & 0.719d0, 0.952d0, 1.210d0, 1.500d0, 1.820d0, 2.180d0, 2.580d0,
     & 3.030d0, 3.520d0, 4.080d0, 4.690d0, 5.380d0, 6.150d0, 7.000d0,
     & 7.460d0, 8.470d0, 9.580d0, 10.80d0, 12.20d0, 13.70d0, 15.50d0,
     & 17.40d0, 19.50d0, 21.80d0, 24.40d0, 27.30d0, 30.50d0, 34.10d0,
     & 38.00d0, 40.20d0, 44.80d0, 49.90d0, 55.60d0, 61.80d0, 68.80d0,
     & 76.50d0, 85.00d0, 94.50d0, 105.0d0, 117.0d0, 130.0d0, 144.0d0,
     & 161.0d0, 179.0d0/

c p + nu_e-bar -> n + e+  mean cosine of electron
      data costelm /0.00d0,
     & -0.021d0, -0.025d0, -0.027d0, -0.027d0, -0.027d0, -0.027d0,
     & -0.026d0, -0.025d0, -0.024d0, -0.023d0, -0.022d0, -0.020d0,
     & -0.018d0, -0.016d0, -0.015d0, -0.013d0, -0.010d0, -0.007d0,
     & -0.003d0, 0.0006d0,  0.005d0,  0.010d0,  0.015d0,  0.021d0,
     &  0.028d0,  0.036d0,  0.044d0,  0.054d0,  0.065d0,  0.070d0,
     &  0.083d0,  0.097d0,  0.113d0,  0.131d0,  0.151d0,  0.173d0,
     &  0.198d0,  0.225d0,  0.255d0,  0.288d0,  0.323d0,  0.361d0,
     &  0.400d0,  0.442d0/

      contains

************************************************************************
*                                                                      *
      subroutine neutrino_Xsec(kf,ein,hydro,mat,lem, sigt, sigh)
*                                                                      *
*     calculation of reaction cross-sections                           *
*     modified by T.Ogawa     on 2017/11/10                            *
*                                                                      *
*     input:                                                           *
*        kf     : incident particle kf-code (+12 or -12)               *
*        ein    : incident nu energy (MeV)                             *
*      hydro    : hydrogen density (1.d24 atoms/cm^3)                  *
*        mat    : pointer offset                                       *
*        lem    : # of elements except hydrogen                        *
*                                                                      *
*     output:                                                          *
*       sigh   : hydrogen cross-section (b)                            *
*       sigt   : other total cross-section (b)                         *
*                                                                      *
************************************************************************
      use MEMBANKMOD, only : sigge, siggn
      use moddas_material

      implicit doubleprecision(a-h,o-z)

      include 'param.inc'

      common /kmat1g/ kmat(kvlmax)
      common /xinels/ ksige, ksign


      sigt = 0.d0
      sigh = 0.d0
      signe = 0.d0
      sigel = 0.d0
      signe1 = 0.d0
      sigel1 = 0.d0

      if(ntrnore .eq. 0) return ! Disregard neutrino reaction by default

      if(hydro .gt. 0.d0 ) then

        if(kf .eq. -12) then ! Charge current reaction

          call signuCC(kf,ein,1,1,signe1,sigel1)

          signe = signe1 * hydro
          sigel = sigel1 * hydro
          sigh  = signe  + sigel
        endif

          call signue_e(kf,ein,1,1,signe1,sigel1)

          signe = signe + signe1 * hydro
          sigel = sigel + sigel1 * hydro
          sigh  = signe + sigel

      endif

c      return ! 2018/1/4 Currently, nu+d reactions are not implemented

      signe = 0.d0 ! reset temporal variable
      sigel = 0.d0 !

      do i = 1, lem

         itz = nint( zz_das(kmat(mat)+i) )
         ita = nint( a_das(kmat(mat)+i) )

         if(itz .eq. 1 .and. ita .eq. 2) then ! deuteron breakup

!$OMP CRITICAL (ideux_read_crit)
         if(ideux_read .eq. 0) then ! read deuteron-nu cross section
          ideux_read = 1
          call nu_deu_X_read
         endif
!$OMP END CRITICAL (ideux_read_crit)

          call signuNC(kf,ein,ita,itz,signe1,sigel1)

          signe = signe + signe1*den_das(kmat(mat)+i)
          sigel = sigel + sigel1*den_das(kmat(mat)+i)
         endif

          call signue_e(kf,ein,ita,itz,signe1,sigel1)

          signe = signe + signe1*den_das(kmat(mat)+i)
          sigel = sigel + sigel1*den_das(kmat(mat)+i)
          sigt  = signe + sigel

          siggn(ksign+i) = signe
          sigge(ksige+i) = sigel

          if(abs(kf) .eq. 12) then ! deuteron inverse beta decay
           call signuCC(kf,ein,ita,itz,signe1,sigel1)

           signe = signe + signe1*den_das(kmat(mat)+i)
           sigel = sigel + sigel1*den_das(kmat(mat)+i)
           sigt  = signe + sigel

           siggn(ksign+i) = signe
           sigge(ksige+i) = sigel
          endif

      end do

      end subroutine



************************************************************************
*                                                                      *
      subroutine signue_e(ktyp,ein,iat,izt,signe,sigel)
*                                                                      *
*     calculation of neutrino-electron scattering X-section            *
*     modified by T.Ogawa     on 2018/09/03                            *
*                                                                      *
*     input:                                                           *
*        ktyp   : incident particle kf-code (+-12,+-14,+-16)           *
*        ein    : incident nu energy (MeV)                             *
*        ita    : target mass number                                   *
*        itz    : target charge number                                 *
*                                                                      *
*     output:                                                          *
*       signe   : cross-section (b)                                    *
*       sigel   : dummy                                                *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)
      include 'param-physcnst.inc'

* X-section parameters**********************
      gf    = 1.16637d-11            ! Fermi coupling constant (/MeV^2)
      hbc   = physc(3)               ! MeV fm
      stwsq = 0.2324d0               ! sine of Weinberg angle
********************************************

      sigc = 2.d0 * gf**2 * rstms(12) * 1.d3 / physc(1) * hbc**2 * 1.d-2 ! 0.01 b/fm^2
      signe = 0.d0
      sigel = 0.d0

      if(abs(ktyp) .eq. 12) then ! (anti)electron neutrino + e-
       gl = stwsq + 0.5d0
       gr = stwsq
      else  ! (anti)tau or (anti)mu neutrino + e-
       gl = stwsq - 0.5d0
       gr = stwsq
      endif

      if(ktyp .le. 0) then ! anti particle. invert gl and gr
       temp = gl
       gl = gr
       gr = temp
      endif

      signe = sigc * ((gl**2 + gr**2/3.d0) * ein - gl * gr * rstms(12) *
     & 1.d3 /2.d0)

      if(ktyp .eq. 14) signe = signe + sigc * ein * max(0.d0 ,(1.d0 -
     & rstms(7)**2/(2.d0 * rstms(12) * ein * 1.d-3)))

      signe = signe * izt ! multiply by number of electron/atom
      sigel = sigel * izt ! because X-sections are given per electron

      end subroutine

************************************************************************
*                                                                      *
      subroutine signuCC(kf,ein,ita,itz,signe,sigel)
*                                                                      *
*     calculation of chage current reaction cross-sections             *
*     modified by T.Ogawa     on 2017/11/10                            *
*                                                                      *
*        calculates total, nonelastic and elastic cross-sections       *
*                                                                      *
*     input:                                                           *
*        kf     : incident particle kf-code (+12 or -12)               *
*        ein    : incident nu energy (MeV)                             *
*        ita    : target mass number                                   *
*        itz    : target charge number                                 *
*                                                                      *
*     output:                                                          *
*       signe   : nonelastic cross-section (b)                         *
*       sigel   : sigt-sigr=elastic scattering cross-section (b)       *
*                                                                      *
************************************************************************
      implicit doubleprecision(a-h,o-z)
      include 'err.inc'

       sigel = 0.d0
       signe = 0.d0

******  proton ******
      if(itz .eq. 1 .and. ita .eq. 1 .and. kf .eq. -12) then ! proton + nu_e-bar = neutron + positron
	
       if(ein .lt. 1.806d0) then ! Below threshold
        signe = 0.d0
        return
       endif

       inden = 1
       do while (ein .ge. sigpen(inden) )
        inden = inden + 1
       enddo

c attention! Infinitely interpolated.
        signe = sigpn(inden) + ( sigpn(inden + 1) - sigpn(inden) )
     & * (ein - sigpen(inden))/( sigpen(inden+1)  - sigpen(inden))
        signe = signe * 1.d-16
******  proton end ******

******  deuteron ******
      elseif(itz .eq. 1 .and. ita .eq. 2 .and. abs(kf) .eq. 12) then ! d + nu_e_bar = n + n + e+ or d + nu_e = p + p + e-

       if(ein .lt. 4.026783d0 .and. kf .eq. -12 .or.
     &    ein .lt. 1.440117d0 .and. kf .eq. 12 ) then ! Below threshold
        signe = 0.d0
        return
       elseif(ein .lt. 13.d0) then ! above threshold
        inden = int( ein / 0.2d0) - 6
       elseif(ein .lt. 25.d0) then
        inden = int( ein / 0.5d0) + 33
       elseif(ein .lt. 55.d0) then
        inden = int( ein ) + 58
       elseif(ein .lt. 170.d0) then
        inden = int( ein / 5.d0 ) + 102
       else
        signe = 0.d0
        if(iover .eq. 0) then
         write(ErrCha,*) 'Neutrino energy is
     & too high. reaction disregarded'
        ErrID = 'L:450/R:signuCC/F:neutrinomod.f' !W03_001_001
        call ErrWrite(ErrID,ErrCha)
        endif

        iover = 1
        return
       endif

c attention! Infinitely interpolated.
       if(kf .eq. -12) then ! d + nu_e_bar = n + n + e+
        signe = sigdnn(inden) + ( sigdnn(inden + 1) - sigdnn(inden) )
     &   * (ein - sigden(inden))/(sigden(inden+1) - sigden(inden))
        signe = signe * 1.d-16
       else ! d + nu_e = p + p + e-
        signe = sigdpp(inden) + ( sigdpp(inden + 1) - sigdpp(inden) )
     &   * (ein - sigden(inden))/(sigden(inden+1) - sigden(inden))
        signe = signe * 1.d-16
       endif


      endif
******  deuteron end ******
*-----------------------------------------------------------------------

      return

      end subroutine


************************************************************************
*                                                                      *
      subroutine signuNC(kf,ein,iat,izt,signe,sigel)
*                                                                      *
*     calculation of reaction cross-sections for deuteron beak-up      *
*     reaction cross section                                           *
*     created by T.Ogawa     on 2017/11/20                             *
*                                                                      *
*     input:                                                           *
*        kf     : incident particle kf-code (all types of neutrinos)   *
*        ein    : incident nu energy (MeV)                             *
*                                                                      *
*     output:                                                          *
*       sigt    : total cross-section (b)                              *
*       signe   : nonelastic cross-section (b)                         *
*       sigel   : sigt-sigr=elastic scattering cross-section (b)       *
*                                                                      *
************************************************************************
      implicit doubleprecision(a-h,o-z)
      include 'err.inc'

      sigel = 0.d0
      signe = 0.d0

      if(izt .eq. 1 .and. iat .eq. 2) then ! d + nu_* = n + p + nu_*

       if(ein .lt. 2.2d0) then ! Below threshold
        signe = 0.d0
        return
       elseif(ein .lt. 13.d0) then ! above threshold
        inden = int( ein / 0.2d0) - 6
       elseif(ein .lt. 25.d0) then
        inden = int( ein / 0.5d0) + 33
       elseif(ein .lt. 55.d0) then
        inden = int( ein ) + 58
       elseif(ein .lt. 170.d0) then
        inden = int( ein / 5.d0 ) + 102
       else
        signe = 0.d0
        if(iover .eq. 0) then
         write(ErrCha,*) 'Neutrino energy is
     & too high. reaction disregarded'
        ErrID = 'L:521/R:signuNC/F:neutrinomod.f' !W03_001_002
        call ErrWrite(ErrID,ErrCha)
        endif

        iover = 1
        return
       endif

c attention! Infinitely interpolated.
       if(kf .lt. 0) then ! d + nu_*_bar = n + p + nu_*_bar
        signe = sigdnbpn(inden) + (sigdnbpn(inden+1) - sigdnbpn(inden))
     &   * (ein - sigden(inden))/(sigden(inden+1) - sigden(inden))
        signe = signe * 1.d-16
       else ! d + nu_* = n + p + nu_*
        signe = sigdnupn(inden) + (sigdnupn(inden+1) - sigdnupn(inden))
     &   * (ein - sigden(inden))/(sigden(inden+1) - sigden(inden))
        signe = signe * 1.d-16
       endif


      endif

      return

*-----------------------------------------------------------------------

      end subroutine

*-----------------------------------------------------------------------
*-------------Upto here   Cross section calculation part----------------
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*-------------Below       Kinematics calculation part ------------------
*-----------------------------------------------------------------------

************************************************************************
*                                                                      *
      subroutine neutrino_kineml(ktyp,ein,ilept)
*                                                                      *
*     calculation of electron-neutrino scattering                      *
*     eletron + nu_arbitrary = lepton + nu_arbitrary                 *
*     created by T.Ogawa     on 2018/09/03                             *
*                                                                      *
*     input:                                                           *
*        ktyp   : incident particle kf-code (definitely -12)           *
*        ein    : incident nu energy (MeV)                             *
*        ilept  : product lepton id number                             *
*                                                                      *
*     output:                                                          *
*       iclust, jclust, qclust                                         *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)
      include 'param00.inc'
      include 'param02.inc'
      include 'param-physcnst.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

*------ calculate probability--------------------------------------------
        if(ktyp .eq. 14) then
         stwsq = 0.2324d0               ! sine of Weinberg angle
         gl = stwsq - 0.5d0
         gr = stwsq
         sigmu = ein * max(0.d0 ,(1.d0 - rstms(7)**2/(2.d0 * rstms(12) *
     &    ein * 1.d-3)))
         sigel = ((gl**2 + gr**2/3.d0) * ein - gl * gr * rstms(12) *
     &    1.d3 /2.d0)
         if( sigmu/(sigmu+sigel) .gt. unirn(dummy) ) ilept = 7 ! muon production
        endif

*------ calculate kinematics--------------------------------------------
        flmass = rstms(ilept) * 1.d3 ! final lepton mass in MeV
        stwsq = 0.2324d0               ! sine of Weinberg angle

        if(ktyp .eq. 12) then
        gl = stwsq + 0.5d0
        gr = stwsq
        elseif(ktyp .eq. -12) then
        gl = stwsq
        gr = stwsq + 0.5d0
        elseif(ktyp .eq. 14 .or. ktyp .eq. 16) then
        gl = stwsq - 0.5d0
        gr = stwsq
        elseif(ktyp .eq. -14 .or. ktyp .eq. -16) then
        gl = stwsq
        gr = stwsq - 0.5d0
        endif

c Determine electron recoil energy in lab frame ( = t )
        sigmax = max(gl**2 + gr**2, gl**2 - gl*gr*flmass/ein) ! max of differential X-section
        icont = 1
        do icont = 1, 1000
        t = unirn(dummy) * ein /(1.d0 + rstms(12) * 1.d3 / 2.d0 / ein) ! max is less than ein because of kinematic limit
        if(gl**2 + gr**2 * (1.d0 - t/ein)**2 - gl * gr * flmass
     &  * t/ein**2 .ge. sigmax * unirn(dummy)) exit
        enddo

        pel = sqrt(t**2 + 2.d0 * t * flmass)

        costh = (ein**2 + pel**2 - (ein-t)**2 ) ! electron scattering angle
     &  / (2.d0 * ein * pel )

        sinth = Sqrt(1.d0 - costh**2 )

        phi = 2.d0 * unirn(dummy) * physc(1)
*------ finalization ---------------------------------------------------

c e- + nu  -> lepton + nu
        nclst  =  2
        iclust(1)    = 7
        iclust(2)    = 7

c particle #2   Electron
        jclust(0,2)  = 0
        jclust(1,2)  = 0
        jclust(2,2)  = 0
        jclust(3,2)  = ilept
        jclust(4,2)  = 0
        jclust(5,2)  = ichgf(ilept,kfft(ilept))
        jclust(6,2)  = 0
        jclust(7,2)  = kfft(ilept)
        jclust(8,2)  = 0

        qclust(0,2)  = 0.d0
        qclust(1,2)  = pel * sinth * cos(phi) * 1.d-3
        qclust(2,2)  = pel * sinth * sin(phi) * 1.d-3
        qclust(3,2)  = pel * costh * 1.d-3
        qclust(4,2)  = sqrt(pel**2 + flmass **2) * 1.d-3
        qclust(5,2)  = flmass * 1.d-3
        qclust(6,2)  = 0.d0
        qclust(7,2)  = t ! qclust(4,2) * 1.d3 - rstms(12)*1.d3
        qclust(8,2)  = 1.0d0
        qclust(9,2)  = 0.0d0
        qclust(10,2) = 0.0d0
        qclust(11,2) = 0.0d0
        qclust(12,2) = 0.0d0

c particle #1   Neutrino
        jclust(0,1)  = 0
        jclust(1,1)  = 0
        jclust(2,1)  = 0
        jclust(3,1)  = 11
        jclust(4,1)  = 0
        jclust(5,1)  = 0
        jclust(6,1)  = 0
        jclust(7,1)  = ktyp
        jclust(8,1)  = 0

        qclust(0,1)  = 0.d0
        qclust(1,1)  = -qclust(1,2)
        qclust(2,1)  = -qclust(2,2)
        qclust(3,1)  = ein - qclust(3,2)
        qclust(4,1)  = (ein - t) * 1.d-3 ! sqrt(pne**2 + rneumass **2) * 1.d-3
        qclust(5,1)  = 0.d0 ! rneumass * 1.d-3
        qclust(6,1)  = 0.d0
        qclust(7,1)  = qclust(4,1) * 1.d3 ! qclust(4,1) * 1.d3 - rneumass
        qclust(8,1)  = 1.0d0
        qclust(9,1)  = 0.0d0
        qclust(10,1) = 0.0d0
        qclust(11,1) = 0.0d0
        qclust(12,1) = 0.0d0

      return

*-----------------------------------------------------------------------

      end subroutine


************************************************************************
*                                                                      *
      subroutine neutrino_kinem1(ktyp,iat,izt,ein)
*                                                                      *
*     calculation of proton charge current reaction                    *
*     i.e.  inverse reaction of proton + beta decay                    *
*     proton + nu_e-bar = neutron + positron                           *
*     created by T.Ogawa     on 2017/11/21                             *
*                                                                      *
*     input:                                                           *
*        ktyp   : incident particle kf-code (definitely -12)           *
*        iat    : target mass number                                   *
*        izt    : target charge number                                 *
*        ein    : incident nu energy (MeV)                             *
*                                                                      *
*     output:                                                          *
*       iclust, jclust, qclust                                         *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)
      include 'param00.inc'
      include 'param02.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      parameter ( elmass = 0.511d0)

       sCc = 0.d0
       s_e = 0.d0
       call signuCC(ktyp,ein,iat,izt,sCc,sigel)
       call signue_e(ktyp,ein,iat,izt,s_e,sigel)

       p = unirn(dummy)
       if((p .gt. sCc/(sCc+s_e) .and. ktyp .eq. -12) .or.
     &  ktyp .ne. -12) then ! If ktyp .eq. -12, both inverse beta decay and electron scattering are possible
        call neutrino_kineml(ktyp,ein,12) ! Neutral current reaction with electron
        return
       endif

*------ Inverse beta decay ---------------------------------------------
*------ calculate kinematics--------------------------------------------

      inden = 1
      do while (ein .ge. sigpen(inden) )
       inden = inden + 1
      enddo

c Mean electron energy and angle
      enelme = enelm(inden) + ( enelm(inden + 1) - enelm(inden) )
     & * (ein - sigpen(inden))/( sigpen(inden+1)  - sigpen(inden))
     & - elmass

      costelme = costelm(inden) + (costelm(inden + 1) - costelm(inden))
     & * (ein - sigpen(inden))/( sigpen(inden+1)  - sigpen(inden))

c Gaussian distribution around enelme
      ewid = (ein + rpmass - rnmass - elmass - enelme) * 0.5d0 ! width

      do i = 1, 100
       do ii = 1, 100
       enel = enelme - 6.d0 * ewid + 8.d0 * ewid * unirn(dummy)
       if( (1.d0 - exp(- (enelme - enel)/ewid )) * exp(- (enelme - enel)
     &  /ewid ) .gt. unirn(dummy) ) exit
       enddo


       eneu = ein + rpmass - rnmass - elmass - enel

       if(eneu .le. 0.d0) cycle

       pel = sqrt(enel**2 + 2.d0 * elmass * enel)
       pne = sqrt(eneu**2 + 2.d0 * rnmass * eneu)

       if(eneu .gt. 0.d0 .and. ein - pel - pne .le. 0.d0 ) exit
       if(i .eq. 100) then
        enel = enelme
        exit
       endif
      enddo

c Solve momentum conservation
      pelr = sqrt(2.d0 * (ein**2*pel**2 + ein**2*pne**2 + pel**2*pne**2)
     & - pne**4 - ein**4 - pel**4 )/(2.d0 * ein)
      the  = 2.d0 * pi * unirn(dummy)
      pelx = pelr * cos(the)
      pely = pelr * sin(the)
      pelz = Sqrt(pel**2 - pelr**2) * (-1)**nint(unirn(dummy))
      pnez = ein - pelz

*-----------------------------------------------------------------------

*------ finalization ---------------------------------------------------

c p + nu_e-bar -> n + e+
        nclst  =  2
        iclust(1)    = 2
        iclust(2)    = 7

c particle #1   Neutron
        jclust(0,1)  = 0
        jclust(1,1)  = 0
        jclust(2,1)  = 1
        jclust(3,1)  = 2
        jclust(4,1)  = 0
        jclust(5,1)  = 0
        jclust(6,1)  = 1
        jclust(7,1)  = 2112
        jclust(8,1)  = 0

        qclust(0,1)  = 0.d0
        qclust(1,1)  = -pelx * 1.d-3
        qclust(2,1)  = -pely * 1.d-3
        qclust(3,1)  = pnez * 1.d-3
        qclust(4,1)  = sqrt(pne**2 + rnmass**2) * 1.d-3
        qclust(5,1)  = rnmass * 1.d-3
        qclust(6,1)  = 0.d0
        qclust(7,1)  = qclust(4,1) * 1.d3 - rnmass
        qclust(8,1)  = 1.0d0
        qclust(9,1)  = 0.0d0
        qclust(10,1) = 0.0d0
        qclust(11,1) = 0.0d0
        qclust(12,1) = 0.0d0

c particle #2   Positron
        jclust(0,2)  = 0
        jclust(1,2)  = 0
        jclust(2,2)  = 0
        jclust(3,2)  = 13
        jclust(4,2)  = 0
        jclust(5,2)  = 1
        jclust(6,2)  = 0
        jclust(7,2)  = -11
        jclust(8,2)  = 0

        qclust(0,2)  = 0.d0
        qclust(1,2)  = pelx * 1.d-3
        qclust(2,2)  = pely * 1.d-3
        qclust(3,2)  = pelz * 1.d-3
        qclust(4,2)  = sqrt(pel**2 + elmass**2) * 1.d-3
        qclust(5,2)  = elmass * 1.d-3
        qclust(6,2)  = 0.d0
        qclust(7,2)  = qclust(4,2) * 1.d3 - elmass
        qclust(8,2)  = 1.0d0
        qclust(9,2)  = 0.0d0
        qclust(10,2) = 0.0d0
        qclust(11,2) = 0.0d0
        qclust(12,2) = 0.0d0

      return

*-----------------------------------------------------------------------

      end subroutine

************************************************************************
*                                                                      *
      subroutine neutrino_kinem2(ktyp,iat,izt,ein)
*                                                                      *
*     calculation of deuteron + nutrino reaction                       *
*     created by T.Ogawa     on 2017/11/21                             *
*                                                                      *
*     input:                                                           *
*        ktyp   : incident particle kf-code                            *
*        iat    : target mass number                                   *
*        izt    : target charge number                                 *
*        ein    : incident nu energy (MeV)                             *
*                                                                      *
*     output:                                                          *
*       iclust, jclust, qclust                                         *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)

       sCc = 0.d0
       sNc = 0.d0
       s_e = 0.d0
       call signuCC(ktyp,ein,iat,izt,sCc,sigel)
       call signuNC(ktyp,ein,iat,izt,sNc,sigel)
       call signue_e(ktyp,ein,iat,izt,s_e,sigel)

      if(abs(ktyp) .eq. 12) then ! Both inverse beta decay and breakup are possible
       p = unirn(dummy)
       if(p .lt. sCc/(sCc+sNc+s_e) ) then
        call neutrino_kinem21(ktyp,iat,izt,ein) ! Charge current reaction
       elseif(p .lt. (sCc+sNc)/(sCc+sNc+s_e)) then
        call neutrino_kinem20(ktyp,iat,izt,ein) ! Neutral current reaction
       else
        call neutrino_kineml(ktyp,ein,12) ! Neutral current reaction with electron
       endif
      else ! No charge current reactions below 100 MeV for nu_mu, nu_tau.
       if(p .lt. sNc/(sNc+s_e) ) then
        call neutrino_kinem20(ktyp,iat,izt,ein) ! Neutral current reaction
       else
        call neutrino_kineml(ktyp,ein,12) ! Neutral current reaction with electron
       endif
      endif

      return

*-----------------------------------------------------------------------

      end subroutine

************************************************************************
*                                                                      *
      subroutine neutrino_kinem20(ktyp,iat,izt,ein)
*                                                                      *
*     calculation of neutral current reactions with neutrons           *
*     deuteron + nu_* = neutron + proton + nu_*                        *
*     created by T.Ogawa     on 2017/11/21                             *
*                                                                      *
*     input:                                                           *
*        ktyp   : incident particle kf-code                            *
*        iat    : target mass number                                   *
*        izt    : target charge number                                 *
*        ein    : incident nu energy (MeV)                             *
*                                                                      *
*     output:                                                          *
*       iclust, jclust, qclust                                         *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)
      include 'param00.inc'
      include 'param02.inc'
      include 'param-physcnst.inc'

      dimension eng(2), ang(2)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

*------ sample neutron energy-------------------------------------------

       if(ein .lt. 13.d0) then ! above threshold
        inden = int( ein / 0.2d0) - 6
       elseif(ein .lt. 25.d0) then
        inden = int( ein / 0.5d0) + 33
       elseif(ein .lt. 55.d0) then
        inden = int( ein ) + 58
       elseif(ein .lt. 170.d0) then
        inden = int( ein / 5.d0 ) + 102
       endif

100    sigsum = 0.d0
       do i = 1, 2 ! sample neutron momentum and proton momentum
       if(ktyp .gt. 0) then
        do j = 1, 36 ! determine angle from
         sigsum = sigsum + sig2H_nu_nc(2,0,inden,j) *
     &    (cos(dble(j-1) * 5.d0 /180.d0 ) - cos(dble(j) * 5.d0 /180.d0))
        enddo
        sigsum = sigsum * unirn(dummy)
        j = 0
        do while(sigsum .gt. 0.d0)
         j = j + 1
         sigsum = sigsum - sig2H_nu_nc(2,0,inden,j) *
     &    (cos(dble(j-1) * 5.d0 /180.d0 ) - cos(dble(j) * 5.d0 /180.d0))
        enddo
        djang = sigsum / sig2H_nu_nc(2,0,inden,j) /
     &    (cos(dble(j-1) * 5.d0 /180.d0 ) - cos(dble(j) * 5.d0 /180.d0))
        jang = j ! determine angle here

        sigsum1 = sig2H_nu_nc(2,0,inden,jang) * unirn(dummy)
        l = 0
        do while(sigsum1 .gt. 0.d0) ! determine outgoing energy
         l = l + 1
         sigsum1 = sigsum1 - sig2H_nu_nc(2,l,inden,jang)
        enddo
        dl = sigsum1/sig2H_nu_nc(2,l,inden,jang)

        eng(i) = sig2H_nu_nc(1,l,inden,jang) + dl *
     &   (sig2H_nu_nc(1,l,inden,jang) - sig2H_nu_nc(1,l-1,inden,jang))
        ang(i) = 5.d0 * (1.d0 - djang + dble(jang))/180.d0 * pi

       else
        do j = 1, 36 ! determine angle from
         sigsum = sigsum + sig2H_nu_bar_nc(2,0,inden,j) *
     &    (cos(dble(j-1) * 5.d0 /180.d0 ) - cos(dble(j) * 5.d0 /180.d0))
        enddo
        sigsum = sigsum * unirn(dummy)
        j = 0
        do while(sigsum .gt. 0.d0)
         j = j + 1
         sigsum = sigsum - sig2H_nu_bar_nc(2,0,inden,j) *
     &    (cos(dble(j-1) * 5.d0 /180.d0 ) - cos(dble(j) * 5.d0 /180.d0))
        enddo
        djang = sigsum / sig2H_nu_bar_nc(2,0,inden,j) /
     &    (cos(dble(j-1) * 5.d0 /180.d0 ) - cos(dble(j) * 5.d0 /180.d0))
        jang = j ! determine angle here

        sigsum1 = sig2H_nu_bar_nc(2,0,inden,jang) * unirn(dummy)
        l = 0
        do while(sigsum1 .gt. 0.d0) ! determine outgoing energy
         l = l + 1
         sigsum1 = sigsum1 - sig2H_nu_bar_nc(2,l,inden,jang)
        enddo
        dl = sigsum1/sig2H_nu_bar_nc(2,l,inden,jang)

        eng(i) = sig2H_nu_bar_nc(1,l,inden,jang) + dl *
     &   (sig2H_nu_bar_nc(1,l,inden,jang)
     &  - sig2H_nu_bar_nc(1,l-1,inden,jang))
        ang(i) = 5.d0 * (1.d0 - djang + dble(jang))/180.d0 * pi
       endif
       enddo

       eng = eng * 1.d-3 ! conversion from MeV to GeV

*-----------------------------------------------------------------------

*------ calculate kinematics--------------------------------------------


*-----------------------------------------------------------------------

*------ finalization ---------------------------------------------------

c d + nu_* -> p + n + nu_*  ! deuteron break-up
        nclst  =  3
        iclust(1)    = 7
        iclust(2)    = 1
        iclust(3)    = 2

c particle #1   Neutrino
        jclust(0,1)  = 0
        jclust(1,1)  = 0
        jclust(2,1)  = 0
        jclust(3,1)  = 11
        jclust(4,1)  = 0
        jclust(5,1)  = 0
        jclust(6,1)  = 0
        jclust(7,1)  = ktyp
        jclust(8,1)  = 0

c particle #2   Proton
        jclust(0,2)  = 0
        jclust(1,2)  = 1
        jclust(2,2)  = 0
        jclust(3,2)  = 1
        jclust(4,2)  = 0
        jclust(5,2)  = 1
        jclust(6,2)  = 1
        jclust(7,2)  = 2212
        jclust(8,2)  = 0

c particle #3   Neutron
        jclust(0,3)  = 0
        jclust(1,3)  = 0
        jclust(2,3)  = 1
        jclust(3,3)  = 2
        jclust(4,3)  = 0
        jclust(5,3)  = 0
        jclust(6,3)  = 1
        jclust(7,3)  = 2112
        jclust(8,3)  = 0

c proton
        phi = 2.d0 * pi * unirn(dummy)
        qclust(1,2) = sqrt(eng(2)**2 + 2.d0 * rstms(1) * eng(2)) *
     &   sin(ang(2)) * Cos(phi)
        qclust(2,2) = sqrt(eng(2)**2 + 2.d0 * rstms(1) * eng(2)) *
     &   sin(ang(2)) * Sin(phi)
        qclust(3,2) = sqrt(eng(2)**2 + 2.d0 * rstms(1) * eng(2)) *
     &   cos(ang(2))
        qclust(4,2) = rstms(1) + eng(2)
        qclust(5,2) = rstms(1)
        qclust(7,2) = eng(2) * 1.d3

c neutron
        phi = 2.d0 * pi * unirn(dummy)
        qclust(1,3) = sqrt(eng(1)**2 + 2.d0 * rstms(2) * eng(1)) *
     &   sin(ang(1)) * Cos(phi)
        qclust(2,3) = sqrt(eng(1)**2 + 2.d0 * rstms(2) * eng(1)) *
     &   sin(ang(1)) * Sin(phi)
        qclust(3,3) = sqrt(eng(1)**2 + 2.d0 * rstms(2) * eng(1)) *
     &   cos(ang(1))
        qclust(4,3) = rstms(2) + eng(1)
        qclust(5,3) = rstms(2)
        qclust(7,3) = eng(1) * 1.d3

c netrino*        p and Ekin are determined by conservation law
        qclust(1,1)  = - qclust(1,2) - qclust(1,3)
        qclust(2,1)  = - qclust(2,2) - qclust(2,3)
        qclust(3,1)  = - qclust(3,2) - qclust(3,3) + ein * 1.d-3
        qclust(4,1)  = sqrt(qclust(1,1)**2 + qclust(2,1)**2
     &   + qclust(3,1)**2 )
        qclust(5,1)  = 0.d0
        qclust(7,1)  = qclust(4,1) * 1.d3

        iloop = 0
        do while(abs(log((2.2d0 + sum(qclust(7,1:3))) / ein)) .gt.1.d-4)
         iloop = iloop + 1

         if(abs(log((2.2d0 + sum(qclust(7,1:3)))/ein)) .gt. 2.5d0)
     &    goto 100

         if(iloop .gt. 100) goto 100 !
c When energy is too low, scale up everything to adjust
         if( (2.2d0 + sum(qclust(7,1:3))) .lt. ein ) then
          do i = 1, 3
           do j = 1, 2
            qclust(j,i) = qclust(j,i)/((2.2d0 + sum(qclust(7,1:3)))/ein)
           enddo
           qclust(4,i) = sqrt(qclust(1,i)**2 + qclust(2,i)**2 +
     &      qclust(3,i)**2 + qclust(5,i)**2 )
           qclust(7,i) = (qclust(4,i) - qclust(5,i)) * 1.d3
          enddo
         else
c When energy is too high, it is nu to blame. Scale up momentum of p and n to reduce nu momentum
          if( qclust(3,2) + qclust(3,3) .lt. 0.d0 ) goto 100 ! No adjustable scenario
          do i = 2, 3
           do j = 1, 3
            qclust(j,i) = qclust(j,i)*((2.2d0 + sum(qclust(7,1:3)))/ein)
           enddo
           qclust(4,i) = sqrt(qclust(1,i)**2 + qclust(2,i)**2 +
     &      qclust(3,i)**2 + qclust(5,i)**2 )
           qclust(7,i) = (qclust(4,i) - qclust(5,i)) * 1.d3
          enddo
c netrino*        p and Ekin are determined by conservation law
          qclust(1,1)  = - qclust(1,2) - qclust(1,3)
          qclust(2,1)  = - qclust(2,2) - qclust(2,3)
          qclust(3,1)  = - qclust(3,2) - qclust(3,3) + ein * 1.d-3
          qclust(4,1)  = sqrt(qclust(1,1)**2 + qclust(2,1)**2
     &     + qclust(3,1)**2 )
          qclust(5,1)  = 0.d0
          qclust(7,1)  = qclust(4,1) * 1.d3
         endif
        enddo



        qclust(6,1:3) = 0.d0
        qclust(8,1:3) = 1.d0
        qclust(9,1:3) = 0.d0
        qclust(10,1:3) = 0.d0
        qclust(11,1:3) = 0.d0
        qclust(12,1:3) = 0.d0

      return

*-----------------------------------------------------------------------

      end subroutine

************************************************************************
*                                                                      *
      subroutine neutrino_kinem21(ktyp,iat,izt,ein)
*                                                                      *
*     calculation of charge current reaction                           *
*     deuteron + nu_e = proton + proton + electron                     *
*     created by T.Ogawa     on 2017/11/21                             *
*                                                                      *
*     input:                                                           *
*        ktyp   : incident particle kf-code (definitely 12 or -12)     *
*        iat    : target mass number                                   *
*        izt    : target charge number                                 *
*        ein    : incident nu energy (MeV)                             *
*                                                                      *
*     output:                                                          *
*       iclust, jclust, qclust                                         *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)
      include 'param00.inc'
      include 'param02.inc'
      include 'param-physcnst.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

*-----------------------------------------------------------------------


*------ sample lepton energy-------------------------------------------

       if(ein .lt. 13.d0) then ! above threshold
        inden = int( ein / 0.2d0) - 6
       elseif(ein .lt. 25.d0) then
        inden = int( ein / 0.5d0) + 33
       elseif(ein .lt. 55.d0) then
        inden = int( ein ) + 58
       elseif(ein .lt. 170.d0) then
        inden = int( ein / 5.d0 ) + 102
       endif

100    sigsum = 0.d0
       if(ktyp .gt. 0) then
        do j = 1, 36 ! determine angle from
         sigsum = sigsum + sig2H_nu_e_cc(2,0,inden,j) *
     &    (cos(dble(j-1) * 5.d0 /180.d0 ) - cos(dble(j) * 5.d0 /180.d0))
        enddo
        sigsum = sigsum * unirn(dummy)
        j = 0
        do while(sigsum .gt. 0.d0)
         j = j + 1
         sigsum = sigsum - sig2H_nu_e_cc(2,0,inden,j) *
     &    (cos(dble(j-1) * 5.d0 /180.d0 ) - cos(dble(j) * 5.d0 /180.d0))
        enddo
        djang = sigsum / sig2H_nu_e_cc(2,0,inden,j) /
     &    (cos(dble(j-1) * 5.d0 /180.d0 ) - cos(dble(j) * 5.d0 /180.d0))
        jang = j ! determine angle here

        sigsum1 = sig2H_nu_e_cc(2,0,inden,jang) * unirn(dummy)
        l = 0
        do while(sigsum1 .gt. 0.d0) ! determine outgoing energy
         l = l + 1
         sigsum1 = sigsum1 - sig2H_nu_e_cc(2,l,inden,jang)
        enddo
        dl = sigsum1/sig2H_nu_e_cc(2,l,inden,jang)

        eng = sig2H_nu_e_cc(1,l,inden,jang) + dl *
     & (sig2H_nu_e_cc(1,l,inden,jang) - sig2H_nu_e_cc(1,l-1,inden,jang))
        ang = 5.d0 * (1.d0 - djang + dble(jang))/180.d0 * pi

       else
        do j = 1, 36 ! determine angle from
         sigsum = sigsum + sig2H_nu_e_b_cc(2,0,inden,j) *
     &    (cos(dble(j-1) * 5.d0 /180.d0 ) - cos(dble(j) * 5.d0 /180.d0))
        enddo
        sigsum = sigsum * unirn(dummy)
        j = 0
        do while(sigsum .gt. 0.d0)
         j = j + 1
         sigsum = sigsum - sig2H_nu_e_b_cc(2,0,inden,j) *
     &    (cos(dble(j-1) * 5.d0 /180.d0 ) - cos(dble(j) * 5.d0 /180.d0))
        enddo
        djang =  sigsum / sig2H_nu_e_b_cc(2,0,inden,j) /
     &    (cos(dble(j-1) * 5.d0 /180.d0 ) - cos(dble(j) * 5.d0 /180.d0))
        jang = j ! determine angle here

        sigsum1 = sig2H_nu_e_b_cc(2,0,inden,jang) * unirn(dummy)
        l = 0
        do while(sigsum1 .gt. 0.d0) ! determine outgoing energy
         l = l + 1
         sigsum1 = sigsum1 - sig2H_nu_e_b_cc(2,l,inden,jang)
        enddo
        dl = sigsum1/sig2H_nu_e_b_cc(2,l,inden,jang)

        eng = sig2H_nu_e_b_cc(1,l,inden,jang) + dl *
     &   (sig2H_nu_e_b_cc(1,l,inden,jang)
     &  - sig2H_nu_e_b_cc(1,l-1,inden,jang))
        ang = 5.d0 * (1.d0 - djang + dble(jang))/180.d0 * pi
       endif

       eng = eng * 1.d-3 ! conversion from MeV to GeV

*	inden	incident energy bin
*	l ~ l+1	outgoing energy bin
*     dl      outgoing energy in-bin correction
*	jang
*      ~ jang+1	outgoing angular bin
*     djang   outgoing angle in-bin correction
*-----------------------------------------------------------------------


*------ finalization ---------------------------------------------------

        nclst  =  3
        if(ktyp .gt. 0) then ! c d + nu_e -> p + p + e-    ! deuteron inverse beta decay
         iclust(1)    = 1
         iclust(2)    = 1
         iclust(3)    = 7

c particle #1   Proton
         jclust(0,1)  = 0
         jclust(1,1)  = 1
         jclust(2,1)  = 0
         jclust(3,1)  = 1
         jclust(4,1)  = 0
         jclust(5,1)  = 1
         jclust(6,1)  = 1
         jclust(7,1)  = 2212
         jclust(8,1)  = 0

c particle #2   Proton
         jclust(0,2)  = 0
         jclust(1,2)  = 1
         jclust(2,2)  = 0
         jclust(3,2)  = 1
         jclust(4,2)  = 0
         jclust(5,2)  = 1
         jclust(6,2)  = 1
         jclust(7,2)  = 2212
         jclust(8,2)  = 0

c particle #3   Electron
         jclust(0,3)  = 0
         jclust(1,3)  = 0
         jclust(2,3)  = 0
         jclust(3,3)  = 12
         jclust(4,3)  = 0
         jclust(5,3)  = -1
         jclust(6,3)  = 0
         jclust(7,3)  = 11
         jclust(8,3)  = 0

        else  ! c d + nu_e_bar -> n + n + e+    ! deuteron inverse beta decay
         iclust(1)    = 2
         iclust(2)    = 2
         iclust(3)    = 7

c particle #1   Neutron
         jclust(0,1)  = 0
         jclust(1,1)  = 0
         jclust(2,1)  = 1
         jclust(3,1)  = 2
         jclust(4,1)  = 0
         jclust(5,1)  = 0
         jclust(6,1)  = 1
         jclust(7,1)  = 2112
         jclust(8,1)  = 0

c particle #2   Neutron
         jclust(0,2)  = 0
         jclust(1,2)  = 0
         jclust(2,2)  = 1
         jclust(3,2)  = 2
         jclust(4,2)  = 0
         jclust(5,2)  = 0
         jclust(6,2)  = 1
         jclust(7,2)  = 2112
         jclust(8,2)  = 0

c particle #3   Positron
         jclust(0,3)  = 0
         jclust(1,3)  = 0
         jclust(2,3)  = 0
         jclust(3,3)  = 13
         jclust(4,3)  = 0
         jclust(5,3)  = 1
         jclust(6,3)  = 0
         jclust(7,3)  = -11
         jclust(8,3)  = 0

        endif

c sample initial neutron energy
        if(unirn(dummy) .lt. 2.564d-2 ** 2 ) then ! D wave
         ur = 0.625177d0 * unirn(dummy)
         enn1 = 0.04440d0 + 13.45d0 * ur - 296.5d0 * ur**2
     & + 3619.d0 * ur**3 - 2.415d4 * ur**4 + 9.232d4 * ur**5
     & - 1.985d5 * ur**6 + 2.032d5 * ur**7 + 1.538d3 * ur**8
     & - 1.853d5 * ur**9 + 1.133d5 * ur**10 ! Sample energy in GeV from integral function
        else ! S wave
         ur = 356.799d0 * unirn(dummy)
         if(ur .lt. 56.54d0) then
         enn1 = 1.d-2/56.54d0 * ur
         else
         enn1 = 3.243d-172 * Exp(ur) - 0.0202d0 * ur + 0.00118d0 * ur**2
     &  - 2.853d-5  * ur**3 + 3.808d-7  * ur**4 - 3.101d-9  * ur**5
     &  + 1.606d-11 * ur**6 - 5.313d-14 * ur**7 + 1.088d-16 * ur**8
     &  - 1.256d-19 * ur**9 + 6.254d-23 * ur**10 ! Sample energy in GeV from integral function
         endif
        endif

        enn1 = enn1 * 1.d-3

        p1 = sqrt(enn1**2 + 2.d0 * rstms(2) * enn1)
        costhe = unirn(dummy)
        sinthe = sqrt(1.d0 - costhe**2 )
        phi = 2.d0 * pi * unirn(dummy)

c positron or electron
        qclust(1,3) = sqrt(eng ** 2 + 2.d0 * eng* rstms(12)) * sin(ang)
     &   * Cos(phi)
        qclust(2,3) = sqrt(eng ** 2 + 2.d0 * eng* rstms(12)) * sin(ang)
     &   * Sin(phi)
        qclust(3,3) = sqrt(eng ** 2 + 2.d0 * eng* rstms(12)) * cos(ang)
        qclust(4,3) = rstms(12) + eng
        qclust(5,3) = rstms(12)
        qclust(7,3) = eng * 1.d3

        if(ktyp .gt. 0) then ! switch mass
         pmas = rstms(1)
        else
         pmas = rstms(2)
        endif
c neutron or proton 1
        qclust(1,1) = p1 * sinthe * cos(phi)
        qclust(2,1) = p1 * sinthe * sin(phi)
        qclust(3,1) = p1 * costhe
        qclust(4,1) = sqrt(p1**2 + pmas**2)
        qclust(5,1) = pmas
        qclust(7,1) = enn1 * 1.d3

c neutron or proton 2
        qclust(1,2) = - qclust(1,1) - qclust(1,3)
        qclust(2,2) = - qclust(2,1) - qclust(2,3)
        qclust(3,2) = - qclust(3,1) - qclust(3,3) + ein * 1.d-3
        qclust(4,2) = sqrt(qclust(1,2)**2 + qclust(2,2)**2 + qclust(3,2)
     &   **2 + pmas**2 )
        qclust(5,2) = pmas
        qclust(7,2) = (qclust(4,2) - qclust(5,2)) * 1.d3

        r = (2.2d-3 + 2.d0 * pmas - rstms(1) - rstms(2) + rstms(12)
     &   + sum(qclust(7,1:3)) * 1.d-3)/(ein * 1.d-3)
        if( r  .gt. 2.5d0)  goto 100 ! reject. Energy is hardly rebalanced.
        iloop = 0
        do while(abs(log( r )) .gt.1.d-4)
        iloop = iloop + 1
        if(iloop .gt. 100) goto 100 !
         do i = 1, 3
          do j = 1, 2
           qclust(j,i) = qclust(j,i)/ r
          enddo
          qclust(4,i) = sqrt(qclust(1,i)**2 + qclust(2,i)**2 +
     &     qclust(3,i)**2 + qclust(5,i)**2 )
          qclust(7,i) = (qclust(4,i) - qclust(5,i)) * 1.d3
         enddo
        r = (2.2d-3 + 2.d0 * pmas - rstms(1) - rstms(2) + rstms(12)
     &   + sum(qclust(7,1:3)) * 1.d-3)/(ein * 1.d-3)
        enddo

        qclust(6,1:3) = 0.d0
        qclust(8,1:3) = 1.d0
        qclust(9,1:3) = 0.d0
        qclust(10,1:3) = 0.d0
        qclust(11,1:3) = 0.d0
        qclust(12,1:3) = 0.d0

      return

*-----------------------------------------------------------------------

      end subroutine

************************************************************************
*                                                                      *
      subroutine neutrino_kinem3(ktyp,iat,izt,ein)
*                                                                      *
*     calculation of arbitrary nucleus + nutrino reaction              *
*     created by T.Ogawa     on 2018/09/06                             *
*                                                                      *
*     input:                                                           *
*        ktyp   : incident particle kf-code                            *
*        iat    : target mass number                                   *
*        izt    : target charge number                                 *
*        ein    : incident nu energy (MeV)                             *
*                                                                      *
*     output:                                                          *
*       iclust, jclust, qclust                                         *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)
      include 'err.inc'

       sCc = 0.d0
       sNc = 0.d0
       s_e = 0.d0
       call signuCC(ktyp,ein,iat,izt,sCc,sigel) ! currently 0
       call signuNC(ktyp,ein,iat,izt,sNc,sigel) ! currently 0
       call signue_e(ktyp,ein,iat,izt,s_e,sigel)

       p = unirn(dummy)
       if(p .lt. sCc/(sCc+sNc+s_e) ) then
        Write(ErrCha,*) "Charge current neutrino reaction with A
     & > 2 nucleus. Error."
       ErrID = 'L:1463/R:neutrino_kinem3/F:neutrinomod.f' !E00_008_001
       call ErrWrite(ErrID,ErrCha)

       elseif(p .lt. (sCc+sNc)/(sCc+sNc+s_e)) then
        Write(ErrCha,*) "Charge current neutrino reaction with A
     & > 2 nucleus. Error."
       ErrID = 'L:1469/R:neutrino_kinem3/F:neutrinomod.f' !E00_008_002
       call ErrWrite(ErrID,ErrCha)
       else
        call neutrino_kineml(ktyp,ein,12) ! Neutral current reaction with electron
       endif

      return

*-----------------------------------------------------------------------

      end subroutine

************************************************************************
*                                                                      *
      subroutine nu_deu_X_read
*                                                                      *
*     Read cross sections for deutron-induced reactions                *
*     Neutral current (nu + d -> n + p, nu_bar + d -> n+ p)            *
*     Charge current (nu_e + d -> n + n + e+, nu_e_bar + d -> p + p + e-) *
*     created by T.Ogawa     on 2018/01/18                             *
*                                                                      *
*     output:                                                          *
*                                                *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)
      dimension ne(4)
      character      hnp(34)*11
      character*230   xsdir, hdpath, hdpth, klin
      character hl*281, htn*10, yen*1
      common /dircha/ idirch
      common /gm007/ xsdir, hdpath, hdpth, klin, hnp, htn

*-----------------------------------------------------------------------

      yen  = char(92)
      inp = 26
      sig2H_nu_e_b_cc = 0.d0
      sig2H_nu_e_cc   = 0.d0
      sig2H_nu_nc     = 0.d0
      sig2H_nu_bar_nc = 0.d0


          if( idirch .eq. 0 ) then
             if( leng(hdpth) .gt. 0 )
     &         hl = hdpth(1:leng(hdpth))//'/'//'ddx'//'/'//'nu'//'/'//
     &'CC'//'/'//'nu_e_bar+2H.dat'! //'/'//cheng(k)//'_'//chang(j)
          else
             if( leng(hdpth) .gt. 0 )
     &         hl = hdpth(1:leng(hdpth))//yen//'ddx'//yen//'nu'//yen//
     &'CC'//yen//'nu_e_bar+2H.dat'! //yen//cheng(k)//'_'//chang(j)
          end if
          open(inp,file=hl,status='old')
       do j = 1, 37! emission angle
        do k = 15, 79 ! incident energy
         l = 0
         sigac = 0.d0
         do
          read(inp,*,end=100) sig2H_nu_e_b_cc(1,l+1,k,j),
     &     sig2H_nu_e_b_cc(2,l+1,k,j) ! read energy and X section
          if(l.ge.1) sigac = sigac + sig2H_nu_e_b_cc(2,l+1,k,j) *
     &     (sig2H_nu_e_b_cc(1,l+1,k,j) - sig2H_nu_e_b_cc(1,l,k,j))
          if(sig2H_nu_e_b_cc(1,l+1,k,j) .ne. 0.d0 .and.
     &       sig2H_nu_e_b_cc(2,l+1,k,j) .eq. 0.d0) exit ! end of one data section
          l = l + 1
         enddo
  100    sig2H_nu_e_b_cc(1,0,k,j) = l ! valid table length
         sig2H_nu_e_b_cc(2,0,k,j) = sigac
        enddo
       enddo


          if( idirch .eq. 0 ) then
             if( leng(hdpth) .gt. 0 )
     &         hl = hdpth(1:leng(hdpth))//'/'//'ddx'//'/'//'nu'//'/'//
     &'CC'//'/'//'nu_e+2H.dat'! //'/'//cheng(k)//'_'//chang(j)
          else
             if( leng(hdpth) .gt. 0 )
     &         hl = hdpth(1:leng(hdpth))//yen//'ddx'//yen//'nu'//yen//
     &'CC'//yen//'nu_e+2H.dat'! //yen//cheng(k)//'_'//chang(j)
          end if
          open(inp,file=hl,status='old')
       do j = 1, 37! emission angle
        do k = 1, 79 ! incident energy
         l = 0
         sigac = 0.d0
         do
          read(inp,*,end=101) sig2H_nu_e_cc(1,l+1,k,j),
     &     sig2H_nu_e_cc(2,l+1,k,j) ! read energy and X section
          if(l.ge.1) sigac = sigac + sig2H_nu_e_cc(2,l+1,k,j) *
     &     (sig2H_nu_e_cc(1,l+1,k,j) - sig2H_nu_e_cc(1,l,k,j))
          if(sig2H_nu_e_cc(1,l+1,k,j) .ne. 0.d0 .and.
     &       sig2H_nu_e_cc(2,l+1,k,j) .eq. 0.d0) exit ! end of one data section
          l = l + 1
         enddo
  101    sig2H_nu_e_cc(1,0,k,j) = l ! valid table length
         sig2H_nu_e_cc(2,0,k,j) = sigac
        enddo
       enddo


          if( idirch .eq. 0 ) then
             if( leng(hdpth) .gt. 0 )
     &         hl = hdpth(1:leng(hdpth))//'/'//'ddx'//'/'//'nu'//'/'//
     &'NC'//'/'//'nu+2H.dat'! //'/'//cheng(k)//'_'//chang(j)
          else
             if( leng(hdpth) .gt. 0 )
     &         hl = hdpth(1:leng(hdpth))//yen//'ddx'//yen//'nu'//yen//
     &'NC'//yen//'nu+2H.dat'! //yen//cheng(k)//'_'//chang(j)
          end if
          open(inp,file=hl,status='old')
       do j = 1, 37! emission angle
        do k = 5, 72 ! incident energy
         l = 0
         sigac = 0.d0
         do
          read(inp,*,end=102) sig2H_nu_nc(1,l+1,k,j),
     &     sig2H_nu_nc(2,l+1,k,j) ! read energy and X section
          if(l.ge.1) sigac = sigac + sig2H_nu_nc(2,l+1,k,j) *
     &     (sig2H_nu_nc(1,l+1,k,j) - sig2H_nu_nc(1,l,k,j))
          if(sig2H_nu_nc(1,l+1,k,j) .ne. 0.d0 .and.
     &       sig2H_nu_nc(2,l+1,k,j) .eq. 0.d0) exit ! end of one data section
          l = l + 1
         enddo
  102    sig2H_nu_nc(1,0,k,j) = l ! valid table length
         sig2H_nu_nc(2,0,k,j) = sigac
        enddo
       enddo


          if( idirch .eq. 0 ) then
             if( leng(hdpth) .gt. 0 )
     &         hl = hdpth(1:leng(hdpth))//'/'//'ddx'//'/'//'nu'//'/'//
     &'NC'//'/'//'nu_bar+2H.dat'! //'/'//cheng(k)//'_'//chang(j)
          else
             if( leng(hdpth) .gt. 0 )
     &         hl = hdpth(1:leng(hdpth))//yen//'ddx'//yen//'nu'//yen//
     &'NC'//yen//'nu_bar+2H.dat'! //yen//cheng(k)//'_'//chang(j)
          end if
          open(inp,file=hl,status='old')
       do j = 1, 37! emission angle
        do k = 5, 72 ! incident energy
         l = 0
         sigac = 0.d0
         do
          read(inp,*,end=103) sig2H_nu_bar_nc(1,l+1,k,j),
     &     sig2H_nu_bar_nc(2,l+1,k,j) ! read energy and X section
          if(l.ge.1) sigac = sigac + sig2H_nu_bar_nc(2,l+1,k,j) *
     &     (sig2H_nu_bar_nc(1,l+1,k,j) - sig2H_nu_bar_nc(1,l,k,j))
          if(sig2H_nu_bar_nc(1,l+1,k,j) .ne. 0.d0 .and.
     &       sig2H_nu_bar_nc(2,l+1,k,j) .eq. 0.d0) exit ! end of one data section
          l = l + 1
         enddo
  103    sig2H_nu_bar_nc(1,0,k,j) = l ! valid table length
         sig2H_nu_bar_nc(2,0,k,j) = sigac
        enddo
       enddo

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

      end subroutine

      end module

