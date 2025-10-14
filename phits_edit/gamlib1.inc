      parameter (maxph=  200)
      parameter (len1 =41600)
      parameter (len2 = 1600)
      parameter (len3 =71000)
      parameter (lms  =   20)
C
      common/ big1   /izao(len2),ipt(len2),iref(105)
      common/ big2   /elo(len1),isp(len1)              ! Caution. Common blocks "big1" and "big2" appear in gem.f, too.
      common/ big3   /idk(len1),tlev(len1)
      common/ gross1 /pmode(len3)
      common/ gross2 /lfin(len3)
