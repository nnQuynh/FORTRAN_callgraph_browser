************************************************************************
*                                                                      *
      block data nntable
*                                                                      *
*              percent abundance of natural nuclei                     *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include 'param01.inc'

*-----------------------------------------------------------------------

      common /natnuc/ natnn(maxpt), natnm(maxpt,10), patnn(maxpt,10)

*-----------------------------------------------------------------------

      data   natnn(1), ( natnm(1,i), patnn(1,i), i = 1, 2 ) / ! H
     &       2,
     &       1,  99.985,
     &       2,   0.015 /

      data   natnn(2), ( natnm(2,i), patnn(2,i), i = 1, 2 ) / ! He
     &       2,
     &       3,   0.000137,
     &       4,  99.999863 /

      data   natnn(3), ( natnm(3,i), patnn(3,i), i = 1, 2 ) / ! Li
     &       2,
     &       6,   7.59,
     &       7,  92.41 /

      data   natnn(4), ( natnm(4,i), patnn(4,i), i = 1, 1 ) / ! Be
     &       1,
     &       9, 100.0 /

      data   natnn(5), ( natnm(5,i), patnn(5,i), i = 1, 2 ) / ! B
     &       2,
     &      10,  19.9,
     &      11,  80.1 /

      data   natnn(6), ( natnm(6,i), patnn(6,i), i = 1, 2 ) / ! C
     &       2,
     &      12,  98.89,
     &      13,   1.11 /

      data   natnn(7), ( natnm(7,i), patnn(7,i), i = 1, 2 ) / ! N
     &       2,
     &      14,  99.634,
     &      15,   0.366 /

      data   natnn(8), ( natnm(8,i), patnn(8,i), i = 1, 1 ) / ! O
     &       1,
     &      16,  100.0 /

      data   natnn(9), ( natnm(9,i), patnn(9,i), i = 1, 1 ) / ! F
     &       1,
     &      19, 100.0 /

      data   natnn(10), ( natnm(10,i), patnn(10,i), i = 1, 1 ) / ! Ne
     &       1,
     &      20,  100.0 /

      data   natnn(11), ( natnm(11,i), patnn(11,i), i = 1, 1 ) / ! Na
     &       1,
     &      23, 100.0 /

      data   natnn(12), ( natnm(12,i), patnn(12,i), i = 1, 3 ) / ! Mg
     &       3,
     &      24,  78.99,
     &      25,  10.0,
     &      26,  11.01 /

      data   natnn(13), ( natnm(13,i), patnn(13,i), i = 1, 1 ) / ! Al
     &       1,
     &      27, 100.0 /

      data   natnn(14), ( natnm(14,i), patnn(14,i), i = 1, 3 ) / ! Si
     &       3,
     &      28,  92.23,
     &      29,   4.67,
     &      30,   3.10 /

      data   natnn(15), ( natnm(15,i), patnn(15,i), i = 1, 1 ) / ! P
     &       1,
     &      31, 100.0 /

      data   natnn(16), ( natnm(16,i), patnn(16,i), i = 1, 4 ) / ! S
     &       4,
     &      32,  95.02,
     &      33,   0.75,
     &      34,   4.21,
     &      36,   0.02 /

      data   natnn(17), ( natnm(17,i), patnn(17,i), i = 1, 2 ) / ! Cl
     &       2,
     &      35,  75.77,
     &      37,  24.23 /

      data   natnn(18), ( natnm(18,i), patnn(18,i), i = 1, 1 ) / ! Ar
     &       1,
     &      40,  100.0 /

      data   natnn(19), ( natnm(19,i), patnn(19,i), i = 1, 3 ) / ! K
     &       3,
     &      39,  93.2581,
     &      40,   0.0117,
     &      41,   6.7302 /

      data   natnn(20), ( natnm(20,i), patnn(20,i), i = 1, 6 ) / ! Ca
     &       6,
     &      40,  96.941,
     &      42,   0.647,
     &      43,   0.135,
     &      44,   2.086,
     &      46,   0.004,
     &      48,   0.187 /

      data   natnn(21), ( natnm(21,i), patnn(21,i), i = 1, 1 ) / ! Sc
     &       1,
     &      45, 100.0 /

      data   natnn(22), ( natnm(22,i), patnn(22,i), i = 1, 5 ) / ! Ti
     &       5,
     &      46,   8.25,
     &      47,   7.44,
     &      48,  73.72,
     &      49,   5.41,
     &      50,   5.18 /

      data   natnn(23), ( natnm(23,i), patnn(23,i), i = 1, 2 ) / ! V
     &       2,
     &      50,   0.25,
     &      51,  99.75 /

      data   natnn(24), ( natnm(24,i), patnn(24,i), i = 1, 4 ) / ! Cr
     &       4,
     &      50,   4.345,
     &      52,  83.789,
     &      53,   9.501,
     &      54,   2.365 /

      data   natnn(25), ( natnm(25,i), patnn(25,i), i = 1, 1 ) / ! Mn
     &       1,
     &      55, 100.0 /

      data   natnn(26), ( natnm(26,i), patnn(26,i), i = 1, 4 ) / ! Fe
     &       4,
     &      54,   5.845,
     &      56,  91.754,
     &      57,   2.119,
     &      58,   0.282 /

      data   natnn(27), ( natnm(27,i), patnn(27,i), i = 1, 1 ) / ! Co
     &       1,
     &      59, 100.0 /

      data   natnn(28), ( natnm(28,i), patnn(28,i), i = 1, 5 ) / ! Ni
     &       5,
     &      58,  68.077,
     &      60,  26.223,
     &      61,   1.140,
     &      62,   3.634,
     &      64,   0.926 /

      data   natnn(29), ( natnm(29,i), patnn(29,i), i = 1, 2 ) / ! Cu
     &       2,
     &      63,  69.17,
     &      65,  30.83 /

      data   natnn(30), ( natnm(30,i), patnn(30,i), i = 1, 5 ) / ! Zn
     &       5,
     &      64,  48.6,
     &      66,  27.9,
     &      67,   4.1,
     &      68,  18.8,
     &      70,   0.6 /

      data   natnn(31), ( natnm(31,i), patnn(31,i), i = 1, 2 ) / ! Ga
     &       2,
     &      69,  60.108,
     &      71,  39.892 /

      data   natnn(32), ( natnm(32,i), patnn(32,i), i = 1, 5 ) / ! Ge
     &       5,
     &      70,  21.23,
     &      72,  27.66,
     &      73,   7.73,
     &      74,  35.94,
     &      76,   7.44 /

      data   natnn(33), ( natnm(33,i), patnn(33,i), i = 1, 1 ) / ! As
     &       1,
     &      75, 100.0 /

      data   natnn(34), ( natnm(34,i), patnn(34,i), i = 1, 6 ) / ! Se
     &       6,
     &      74,   0.89,
     &      76,   9.36,
     &      77,   7.63,
     &      78,  23.78,
     &      80,  49.61,
     &      82,   8.73 /

      data   natnn(35), ( natnm(35,i), patnn(35,i), i = 1, 2 ) / ! Br
     &       2,
     &      79,  50.69,
     &      81,  49.31 /

      data   natnn(36), ( natnm(36,i), patnn(36,i), i = 1, 6 ) / ! Kr
     &       6,
     &      78,   0.35,
     &      80,   2.25,
     &      82,  11.6,
     &      83,  11.5,
     &      84,  57.0,
     &      86,  17.3 /

      data   natnn(37), ( natnm(37,i), patnn(37,i), i = 1, 2 ) / ! Rb
     &       2,
     &      85,  72.165,
     &      87,  27.835 /

      data   natnn(38), ( natnm(38,i), patnn(38,i), i = 1, 4 ) / ! Sr
     &       4,
     &      84,   0.56,
     &      86,   9.86,
     &      87,   7.00,
     &      88,  82.58 /

      data   natnn(39), ( natnm(39,i), patnn(39,i), i = 1, 1 ) / ! Y
     &       1,
     &      89, 100.0 /

      data   natnn(40), ( natnm(40,i), patnn(40,i), i = 1, 5 ) / ! Zr
     &       5,
     &      90,  51.45,
     &      91,  11.22,
     &      92,  17.15,
     &      94,  17.38,
     &      96,   2.80 /

      data   natnn(41), ( natnm(41,i), patnn(41,i), i = 1, 1 ) / ! Nb
     &       1,
     &      93, 100.0 /

      data   natnn(42), ( natnm(42,i), patnn(42,i), i = 1, 7 ) / ! Mo
     &       7,
     &      92,  14.84,
     &      94,   9.25,
     &      95,  15.92,
     &      96,  16.68,
     &      97,   9.55,
     &      98,  24.13,
     &     100,   9.63 /

      data   natnn(43) / 0 / ! Tc

      data   natnn(44), ( natnm(44,i), patnn(44,i), i = 1, 7 ) / ! Ru
     &       7,
     &      96,   5.52,
     &      98,   1.88,
     &      99,  12.7,
     &     100,  12.6,
     &     101,  17.0,
     &     102,  31.6,
     &     104,  18.7 /

      data   natnn(45), ( natnm(45,i), patnn(45,i), i = 1, 1 ) / ! Rh
     &       1,
     &     103, 100.0 /

      data   natnn(46), ( natnm(46,i), patnn(46,i), i = 1, 6 ) / ! Pd
     &       6,
     &     102,   1.02,
     &     104,  11.14,
     &     105,  22.33,
     &     106,  27.33,
     &     108,  26.46,
     &     110,  11.72 /

      data   natnn(47), ( natnm(47,i), patnn(47,i), i = 1, 2 ) / ! Ag
     &       2,
     &     107,  51.839,
     &     109,  48.161 /

      data   natnn(48), ( natnm(48,i), patnn(48,i), i = 1, 8 ) / ! Cd
     &       8,
     &     106,   1.25,
     &     108,   0.89,
     &     110,  12.49,
     &     111,  12.80,
     &     112,  24.13,
     &     113,  12.22,
     &     114,  28.73,
     &     116,   7.49 /

      data   natnn(49), ( natnm(49,i), patnn(49,i), i = 1, 2 ) / ! In
     &       2,
     &     113,   4.29,
     &     115,  95.71 /

      data   natnn(50), ( natnm(50,i), patnn(50,i), i = 1, 10 ) / ! Sn
     &      10,
     &     112,   0.97,
     &     114,   0.65,
     &     115,   0.34,
     &     116,  14.54,
     &     117,   7.68,
     &     118,  24.23,
     &     119,   8.58,
     &     120,  32.59,
     &     122,   4.63,
     &     124,   5.79 /

      data   natnn(51), ( natnm(51,i), patnn(51,i), i = 1, 2 ) / ! Sb
     &       2,
     &     121,  57.21,
     &     123,  42.79 /

      data   natnn(52), ( natnm(52,i), patnn(52,i), i = 1, 8 ) / ! Te
     &       8,
     &     120,   0.096,
     &     122,   2.603,
     &     123,   0.908,
     &     124,   4.816,
     &     125,   7.139,
     &     126,  18.952,
     &     128,  31.687,
     &     130,  33.799 /

      data   natnn(53), ( natnm(53,i), patnn(53,i), i = 1, 1 ) / ! I
     &       1,
     &     127, 100.0 /

      data   natnn(54), ( natnm(54,i), patnn(54,i), i = 1, 9 ) / ! Xe
     &       9,
     &     124,   0.10,
     &     126,   0.09,
     &     128,   1.91,
     &     129,  26.4,
     &     130,   4.1,
     &     131,  21.2,
     &     132,  26.9,
     &     134,  10.4,
     &     136,   8.9 /

      data   natnn(55), ( natnm(55,i), patnn(55,i), i = 1, 1 ) / ! Cs
     &       1,
     &     133, 100.0 /

      data   natnn(56), ( natnm(56,i), patnn(56,i), i = 1, 7 ) / ! Ba
     &       7,
     &     130,   0.106,
     &     132,   0.101,
     &     134,   2.417,
     &     135,   6.592,
     &     136,   7.854,
     &     137,  11.23,
     &     138,  71.7 /

      data   natnn(57), ( natnm(57,i), patnn(57,i), i = 1, 2 ) / ! La
     &       2,
     &     138,   0.0902,
     &     139,  99.9098 /

      data   natnn(58), ( natnm(58,i), patnn(58,i), i = 1, 2 ) / ! Ce, revised by T.Sato 2023/10/10 because no nuclear data for Ce138
     &       2,
     &     140,  88.84,
     &     142,  11.16 /

      data   natnn(59), ( natnm(59,i), patnn(59,i), i = 1, 1 ) / ! Pr
     &       1,
     &     141, 100.0 /

      data   natnn(60), ( natnm(60,i), patnn(60,i), i = 1, 7 ) / ! Nd
     &       7,
     &     142,  27.13,
     &     143,  12.18,
     &     144,  23.80,
     &     145,   8.30,
     &     146,  17.19,
     &     148,   5.76,
     &     150,   5.64 /

      data   natnn(61) / 0 / ! Pm

      data   natnn(62), ( natnm(62,i), patnn(62,i), i = 1, 7 ) / ! Sm
     &       7,
     &     144,   3.1,
     &     147,  15.0,
     &     148,  11.3,
     &     149,  13.8,
     &     150,   7.4,
     &     152,  26.7,
     &     154,  22.7 /

      data   natnn(63), ( natnm(63,i), patnn(63,i), i = 1, 2 ) / ! Eu
     &       2,
     &     151,  47.8,
     &     153,  52.2 /

      data   natnn(64), ( natnm(64,i), patnn(64,i), i = 1, 7 ) / ! Gd
     &       7,
     &     152,   0.20,
     &     154,   2.18,
     &     155,  14.80,
     &     156,  20.47,
     &     157,  15.65,
     &     158,  24.84,
     &     160,  21.86 /

      data   natnn(65), ( natnm(65,i), patnn(65,i), i = 1, 1 ) / ! Tb
     &       1,
     &     159, 100.0 /

      data   natnn(66), ( natnm(66,i), patnn(66,i), i = 1, 7 ) / ! Dy
     &       7,
     &     156,   0.06,
     &     158,   0.10,
     &     160,   2.34,
     &     161,  18.9,
     &     162,  25.5,
     &     163,  24.9,
     &     164,  28.2 /

      data   natnn(67), ( natnm(67,i), patnn(67,i), i = 1, 1 ) / ! Ho
     &       1,
     &     165, 100.0 /

      data   natnn(68), ( natnm(68,i), patnn(68,i), i = 1, 6 ) / ! Er
     &       6,
     &     162,   0.14,
     &     164,   1.61,
     &     166,  33.6,
     &     167,  22.95,
     &     168,  26.8,
     &     170,  14.9 /

      data   natnn(69), ( natnm(69,i), patnn(69,i), i = 1, 1 ) / ! Tm
     &       1,
     &     169, 100.0 /

      data   natnn(70), ( natnm(70,i), patnn(70,i), i = 1, 7 ) / ! Yb
     &       7,
     &     168,   0.13,
     &     170,   3.05,
     &     171,  14.3,
     &     172,  21.9,
     &     173,  16.12,
     &     174,  31.8,
     &     176,  12.7 /

      data   natnn(71), ( natnm(71,i), patnn(71,i), i = 1, 1 ) / ! Lu
     &       1,
     &     175, 100.0 /

      data   natnn(72), ( natnm(72,i), patnn(72,i), i = 1, 6 ) / ! Hf
     &       6,
     &     174,   0.162,
     &     176,   5.206,
     &     177,  18.606,
     &     178,  27.297,
     &     179,  13.629,
     &     180,  35.100 /

      data   natnn(73), ( natnm(73,i), patnn(73,i), i = 1, 1 ) / ! Ta
     &       1,
     &     181, 100.0 /

      data   natnn(74), ( natnm(74,i), patnn(74,i), i = 1, 5 ) / ! W
     &       5,
     &     180,   0.120,
     &     182,  26.498,
     &     183,  14.314,
     &     184,  30.642,
     &     186,  28.426 /

      data   natnn(75), ( natnm(75,i), patnn(75,i), i = 1, 2 ) / ! Re
     &       2,
     &     185,  37.40,
     &     187,  62.60 /

      data   natnn(76), ( natnm(76,i), patnn(76,i), i = 1, 7 ) / ! Os
     &       7,
     &     184,   0.020,
     &     186,   1.58,
     &     187,   1.6,
     &     188,  13.3,
     &     189,  16.1,
     &     190,  26.4,
     &     192,  41.0 /

      data   natnn(77), ( natnm(77,i), patnn(77,i), i = 1, 2 ) / ! Ir
     &       2,
     &     191,  37.3,
     &     193,  62.7 /

      data   natnn(78), ( natnm(78,i), patnn(78,i), i = 1, 6 ) / ! Pt
     &       6,
     &     190,   0.01,
     &     192,   0.79,
     &     194,  32.9,
     &     195,  33.8,
     &     196,  25.3,
     &     198,   7.2 /

      data   natnn(79), ( natnm(79,i), patnn(79,i), i = 1, 1 ) / ! Au
     &       1,
     &     197, 100.0 /

      data   natnn(80), ( natnm(80,i), patnn(80,i), i = 1, 7 ) / ! Hg
     &       7,
     &     196,   0.15,
     &     198,   9.97,
     &     199,  16.87,
     &     200,  23.10,
     &     201,  13.18,
     &     202,  29.86,
     &     204,   6.87 /

      data   natnn(81), ( natnm(81,i), patnn(81,i), i = 1, 2 ) / ! Tl
     &       2,
     &     203,  29.524,
     &     205,  70.476 /

      data   natnn(82), ( natnm(82,i), patnn(82,i), i = 1, 4 ) / ! Pb
     &       4,
     &     204,   1.4,
     &     206,  24.1,
     &     207,  22.1,
     &     208,  52.4 /

      data   natnn(83), ( natnm(83,i), patnn(83,i), i = 1, 1 ) / ! Bi
     &       1,
     &     209, 100.0 /

      data   ( natnn(i), i = 84, 89 ) / 6*0 /

      data   natnn(90), ( natnm(90,i), patnn(90,i), i = 1, 1 ) / ! Th
     &       1,
     &     232, 100.0 /

      data   natnn(91) / 0 / ! Pa

      data   natnn(92), ( natnm(92,i), patnn(92,i), i = 1, 3 ) / ! U
     &       3,
     &     234,   0.0055,
     &     235,   0.7200,
     &     238,  99.2745 /

      data   ( natnn(i), i = 93, maxpt ) / 44*0 /

*-----------------------------------------------------------------------

      end

************************************************************************
      subroutine readnatural
*     T.Sato 2023/12/26 read natural abundance data if data exist      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param01.inc'
      include 'err.inc'
*-----------------------------------------------------------------------

      common /natnuc/ natnn(maxpt), natnm(maxpt,10), patnn(maxpt,10)

      common /nlibcom/ nlibdef ! T.Sato 2023/12/26, read from natural_abundance.dat
      character*2 nlibdef

      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

      character chatmp1*1,chatmp9*9,chatmp200*200

      nlibdef='50' ! default nlibdef, for JENDL4.0

      io=29

!  for making default table
!      open(io,file='natural_abundance.dat')
!      iexist=1 ! always exist
!      do ip=1,maxpt
!       do is=1,natnn(ip)
!        write(io,'(i1,2i4,f15.10)') iexist,ip,natnm(ip,is),patnn(ip,is)
!       enddo
!      enddo
!      close(io)
      open(io,file=chfn(7)(1:ilfn(7)),iostat=ios,status='old')
      if(ios.eq.0) then ! read only when data exist. If file does not exist error occurs anyway later
       do i=1,100000
        read(io,'(a9)',end=5) chatmp9
        if(chatmp9.eq.'directory') exit
       enddo
       do i=1,100000
        read(io,'(a200)',end=5) chatmp200
        do ii=1,197
         if(chatmp200(ii:ii).eq.'.') then
          if(chatmp200(ii+1:ii+3).eq.'20c'.
     &    or.chatmp200(ii+1:ii+3).eq.'50c') then
           read(chatmp200(ii+1:ii+2),'(a2)') nlibdef ! first neutron data written in xsdir
           goto 5
          endif
         endif
        enddo
       enddo
  5    close(io)
      endif

      open(io,file=chfn(1)(1:ilfn(1))//'/data/natural_abundance_'//
     &nlibdef//'c.dat',iostat=ios,status='old')

      if(ios.ne.0) then ! data do not exist
       write(ErrCha,'("natural_abundance_",a2,"c.dat does not exist ",
     & "in data folder. Original data adjusted for JENDL-4 are used",
     & " for natural isotope expansion")')
     & nlibdef
       ErrID = 'L:625/R:readnatural/F:ndata01.f'
       call ErrWrite(ErrID,ErrCha)
      else  ! data exist so update the natural abundance
       natnn(:)=0
       natnm(:,:)=0.0
       patnn(:,:)=0.0
       read(io,'(a)') chatmp1 ! skip header
       read(io,'(a)') chatmp1 ! skip header
       do i=1,100000
        read(io,*,end=10) iexist,ip,ia,frac
        if(iexist.eq.1) then ! nuclear data exist
         natnn(ip)=natnn(ip)+1
         natnm(ip,natnn(ip))=ia
         patnn(ip,natnn(ip))=frac
        endif
       enddo
  10   close(io)
       do ip=1,maxpt  ! normalized to 100%
        if(natnn(ip).ne.0) then ! natural isotope exist
         sum=0.0
         do is=1,natnn(ip)
          sum=sum+patnn(ip,is)
         enddo
         do is=1,natnn(ip)
          patnn(ip,is)=patnn(ip,is)*100.0d0/sum
         enddo
        endif
       enddo
      endif

      return

      end
