************************************************************************
*                                                                      *
*     angel01.inc                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      parameter ( mh = 42 )  ! up to 20 column
*                 max column terms in one h: section

      integer:: mfc
      parameter ( mfc = 100 )
*                 mfc : max terms of equation
      parameter ( mxcval = 999 )
*                 mxcval : max number of user-defined variables

*-----------------------------------------------------------------------

      parameter ( icolm = 1024 )
*                 icolm : max column number, should be less than ichrl

      parameter ( ichrl = 2048 )
*                 ichrl : length of character for dimension
*                         and max column number for ps file
*     When you change this parameter, you have to change
*     "character dum*2048" appeared twice in a-wtext.f
*     Search "ichrl =" in a-wtext.f and revise them
*-----------------------------------------------------------------------

      parameter ( inig = 500 )
*                 inig : length of character for dimension

      parameter ( mc = 1000 )
*                 mc : max number of line comments (total)
*                      max number of arrows
*                      max number of multi records

      parameter ( numtic = 1000 )
*                 numtic : max number of tic number

      parameter ( ipsm = 100 )
*                   max number of include ps file in one graph

      parameter ( nr = 50 )
      parameter ( nc = 50 )

*     nr: maxmum number of row in one table
*     nc: maxmum number of column in one table

*-----------------------------------------------------------------------

      parameter ( r3max = 1.7976931348623158d+308 ) ! S.H. xorshift (2020.2.6)
*                   max real number

c S.H. extended these values (2022.12.14
      parameter ( r0max = 1.0d+90,  r1max = 1.0d+91 )
      parameter ( r2max = -1.0d+91 )
      parameter ( r9max = 1.0d+89 )
      parameter ( r9min = 1.0d-89 )

      parameter ( pi = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      parameter (jol=18)
*                   output temporary file of lines
      parameter (jil=19)
*                   output temporary file of region line
      parameter (jhc=20)
*                   output temporary file of 3D color plot
      parameter (jhl=21)
*                   output temporary file of 3D cluster plot
      parameter (jer=22)
*                   for error message
      parameter (jps=23)
*                   input ps file
      parameter (jab=24)
*                   angelbat for psview only for unix
      parameter (jht=25)
*                   final eps output for angel
      parameter (jhs=26)
*                   temporary file of gauss smearing fot contour plot
*                   temporary new input file
*                   temporary file for wtab
      parameter (jwt=71)
*                   temporary file for wtab:
      parameter (jhm=72)
*                   temporary file for w: wt:
      parameter (jha=73)
*                   temporary file for aw:
      parameter (jhy=74)
*                   temporary file for y-comment
      parameter (jhb=75)
*                   temporary file for bitmap
      parameter (jhd=76)
*                   temporary file for bitmap
      parameter (jhp=77)
*                   temporary file for bitmap
      parameter (jhq=78)
*                   temporary file for bitmap
      parameter (jhr=79)
*                   temporary file for bitmap

*-----------------------------------------------------------------------

