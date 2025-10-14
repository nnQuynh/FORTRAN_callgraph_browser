************************************************************************
*                                                                      *
      block data angelin
*                                                                      *
*       block data for angel
*       modified by K.Niita on 2002/08/02                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

      common /rval1/ cval(mxcval), aval(mxcval)

*-----------------------------------------------------------------------

      data cval /mxcval*r2max/
      data aval /mxcval*r2max/

      end


************************************************************************
*                                                                      *
      subroutine openf(jsi,jsn,dsin)
*                                                                      *
*              purpose : open include file                             *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      character dsin(0:9)*200

*-----------------------------------------------------------------------

         jsn = jsn + 1
         jsi = jsn + 30

            open(jsi, file = dsin(jsn), status = 'OLD' )

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine closef(jsi,jsn)
*                                                                      *
*              purpose : close include file                            *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

            if( jsi .ge. 31 ) then

                  close(jsi)

                  jsi = jsi - 1
                  jsn = jsn - 1

               if( jsi .eq. 30 ) jsi = 5

            else

                  jsn = 0

            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine inclf(jsn,jsi,dsin,dum,k,ilst,ill,ilf,idsi,ierr)
*                                                                      *
*              purpose : define include file                           *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character dum1*1
      character dsin(0:9)*200
      dimension ill(0:9), ilf(0:9)
      dimension idsi(0:9)

      character m_err*200
      common /error/ m_err, l_err, k_err

      logical exex

*-----------------------------------------------------------------------

         ierr = 0

            if( jsn .ge. 9 ) then

               m_err = 'Number of Include File is too deep. '//
     &                 ' Max Number is 8.'
               ErrCha = ''
               ErrID = 'L:121/R:inclf/F:a-utl00.f'
               l_err = ill(jsn)
               k_err = jsn

               ierr = 1

               return

            end if

*-----------------------------------------------------------------------

   10    k = k + 1

               if( k .gt. ilst ) goto 999

            if( dum(k) .eq. '{' ) then

               ini = k + 1

            else if( dum(k) .eq. '}' ) then

               inf = k - 1

               goto 11

            end if

               goto 10

*-----------------------------------------------------------------------

   11          continue

               if( ini .gt. inf ) goto 999


            do 20 i = 1, 200

                  dsin(jsn+1)(i:i) = ' '

   20       continue

            do 30 i = ini, inf

                  j = i - ini + 1

                  dsin(jsn+1)(j:j) = dum(i)

   30       continue

                  idsi(jsn+1) = inf - ini + 1

*-----------------------------------------------------------------------
*     open include file :  input file dsin exist ?
*-----------------------------------------------------------------------

            inquire( file = dsin(jsn+1), exist = exex )

            if( exex .eqv. .false. ) then

               m_err = 'Input File Name Error in INFL: '//
     &                 ' File not exist ->> '//
     &                           dsin(jsn+1)(1:idsi(jsn+1))
               ErrCha = ''
               ErrID = 'L:186/R:inclf/F:a-utl00.f'
               l_err = ill(jsn)
               k_err = jsn

               ierr  = 1

               return

            end if

               call openf(jsi,jsn,dsin)

*-----------------------------------------------------------------------
*     read line number description
*-----------------------------------------------------------------------

            ill(jsn) = 0
            ilf(jsn) = 100000000 !FURUTA20130416

            imi = -1
            imf = -1
            imm = -1

   40       k = k + 1

            if( k .gt. ilst ) then

               if( imi .eq. -1 ) return

               goto 998

            end if

            if( dum(k) .eq. '[' ) then

               imi = k

            else if( dum(k) .eq. ']' ) then

               imf = k

               if( imi .eq. -1 ) goto 998
               if( imm .eq. -1 ) goto 998

               goto 41

            else if( dum(k) .eq. '-' ) then

               if( imi .eq. -1 ) goto 998

               imm = k

            end if

               goto 40

*-----------------------------------------------------------------------

   41          continue

            if( imi .ne. -1 .and.
     &          imm .ne. -1 .and.
     &          imf .ne. -1 ) then

               if( imf .lt. imm .or. imm .lt. imi )  goto 998

*-----------------------------------------------------------------------

                  icc = imi

                  call qnum(dum,icc,ilst,rnni,jerr,'[','-')

                  if( jerr .ne. 0 ) goto 998

                  ilin = nint(rnni)

                  if( ilin .gt. 0 ) ill(jsn) = ilin - 1

               do 50 i = 1, ill(jsn)

                  read(jsi,'(a1)') dum1

   50          continue

*-----------------------------------------------------------------------

                  icc = imm

                  call qnum(dum,icc,ilst,rnnf,jerr,'-',']')

                  if( jerr .ne. 0 ) goto 998

                  ilfi = nint(rnnf)

                  if( ilfi .gt. 0 ) ilf(jsn) = ilfi

*-----------------------------------------------------------------------

            else

                  goto 998

            end if

                  if( ilf(jsn) .lt. ill(jsn) ) goto 998

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

  999 continue

               m_err = 'File Name is wrong in INFL: '//
     &                 ' File Name is specified by {filename}'
               ErrCha = ''
               ErrID = 'L:303/R:inclf/F:a-utl00.f'
               l_err = ill(jsn)
               k_err = jsn

               ierr = 1
               return

  998 continue

               m_err = 'INFL: Line Number specification is wrong. '//
     &                 ' Format is [3-6].'
               ErrCha = ''
               ErrID = 'L:315/R:inclf/F:a-utl00.f'
               l_err = ill(jsn-1)
               k_err = jsn-1

               ierr = 1

      return
      end


************************************************************************
*                                                                      *
      subroutine pnum(dum,ic,icf,prn,ierr)
*                                                                      *
*        purpose : character to real number                            *
*                  one number from ic to icf                           *
*                  start from (, [, { and end ), ], }                  *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1

      character c2*1

      logical deqn3

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

            prn  = 0.0
            ierr = 1

*-----------------------------------------------------------------------

      if( dum(ic) .eq. '(' ) then

            c2 = ')'

         ici = ic + 1

  100    ic  = ic + 1

            if( ic .gt. icf ) return

            if( dum(ic) .eq. ' ' .or. dum(ic) .eq. tub .or.
     &          deqn3(dum(ic)) ) goto 100

            if( dum(ic) .ne. c2 ) return

               call rnum(prn,dum,ici,ic-1,ierr)

                  if(ierr.ne.0) return

                  ierr = 0

*-----------------------------------------------------------------------

      else if( dum(ic) .eq. '[' .or. dum(ic) .eq. '{' ) then

         if( dum(ic) .eq. '[' ) then

            c2 = ']'

         else if( dum(ic) .eq. '{' ) then

            c2 = '}'

         end if

            call func12(dum,ic,ierr,icf,c2,prn,xval)

                  if( dum(ic) .ne. c2 ) return
                  if(ierr.ne.0) return

                  ierr = 0

      end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine qnum(dum,ic,iclm,prn,ierr,c1,c2)
*                                                                      *
*        purpose : integer only                                        *
*                  character to real number from 'c1' to 'c2'          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character c1*1, c2*1

      character dum(ichrl)*1

      logical deqn1

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

            prn  = 0.0
            ierr = 1

            if( dum(ic) .ne. c1 ) return

         ici = ic + 1

  100    ic  = ic + 1

            if( ic .gt. iclm ) return

            if( dum(ic) .eq. ' ' .or. dum(ic) .eq. tub .or.
     &          deqn1(dum(ic)) ) goto 100

            if( dum(ic) .ne. c2 ) return

               call rnum(prn,dum,ici,ic-1,ierr)

                  if(ierr.ne.0) return

                  ierr = 0

      return
      end

************************************************************************
*                                                                      *
      subroutine rnum(rnm,dum,ini,ifi,ierr)
*                                                                      *
*       purpose:       characters to real number rnm                   *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character rdum*2048

*-----------------------------------------------------------------------

      ierr = 0

            iall = ifi - ini + 1

         do i = ini, ifi
            j = i - ini + 1
            rdum(j:j) = dum(i)
         end do
         do j = iall + 1, 2048
            rdum(j:j) = ' '
         end do

            if( iall .lt. 20 ) then

               read(rdum,'(G20.0)',err=999) rnmd

            else if( iall .lt. 30 ) then

               read(rdum,'(G30.0)',err=999) rnmd

            else if( iall .lt. 40 ) then

               read(rdum,'(G40.0)',err=999) rnmd

            else

               goto 999

            end if

               rnm = rnmd

      return

  999 ierr = 1

      return
      end


************************************************************************
*                                                                      *
      subroutine idd2(dum,ic,iclm,idd,ierr)
*                                                                      *
*       PURPOSE:       CHARACTERS TO INTEGER NUMBER IDD                *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1,d2*2,d3*3,d4*4

      logical deqn1

      data maxdigit/4/ ! T.Sato 2023/04/08 applicable to more than 2 digit

*-----------------------------------------------------------------------

            idd = -1

      do i=0,maxdigit

       if( deqn1(dum(ic+i)).and.ic+i.le.iclm ) then
        continue
       else
        exit
       endif

      enddo

      if(i.eq.0) then ! not number
       return
      elseif(i.eq.1) then ! 1 digit
       read(dum(ic),*) idd
      elseif(i.eq.2) then ! 2 digit
       d2=dum(ic)//dum(ic+1)
       read(d2,*) idd
      elseif(i.eq.3) then ! 3 digit
       d3=dum(ic)//dum(ic+1)//dum(ic+2)
       read(d3,*) idd
      elseif(i.eq.4) then ! 4 digit
       d4=dum(ic)//dum(ic+1)//dum(ic+2)//dum(ic+3)
       read(d4,*) idd
      else
       ierr = 1
      endif
      ic=ic+i

      return
      end

************************************************************************
*                                                                      *
      subroutine chlow(dum,lum)
*                                                                      *
*          Convert a text into lower letters.                          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      character*1 c

      do i = 1, ichrl
         c = dum(i)
         if( c .ge. 'A' .and. c .le. 'Z' )
     &       c = char(ichar(c)+ichar('a')-ichar('A'))
         lum(i) = c
      end do

      return
      end

************************************************************************
*                                                                      *
      function dnen1(dum)
*                                                                      *
************************************************************************
      character dum*1
      logical dnen1
*-----------------------------------------------------------------------

      dnen1  =  dum.ne.'0'.and.
     &          dum.ne.'1'.and.dum.ne.'2'.and.
     &          dum.ne.'3'.and.dum.ne.'4'.and.
     &          dum.ne.'5'.and.dum.ne.'6'.and.
     &          dum.ne.'7'.and.dum.ne.'8'.and.
     &          dum.ne.'9'

      return
      end

************************************************************************
*                                                                      *
      function dnen2(dum)
*                                                                      *
************************************************************************
      character dum*1
      logical dnen2
*-----------------------------------------------------------------------

      dnen2  =  dum.ne.'0'.and.
     &          dum.ne.'1'.and.dum.ne.'2'.and.
     &          dum.ne.'3'.and.dum.ne.'4'.and.
     &          dum.ne.'5'.and.dum.ne.'6'.and.
     &          dum.ne.'7'.and.dum.ne.'8'.and.
     &          dum.ne.'9'.and.dum.ne.'.'.and.
     &          dum.ne.'+'.and.dum.ne.'-'

      return
      end

************************************************************************
*                                                                      *
      function dnen3(dum)
*                                                                      *
************************************************************************
      character dum*1
      logical dnen3
*-----------------------------------------------------------------------

      dnen3  =  dum.ne.'0'.and.
     &          dum.ne.'1'.and.dum.ne.'2'.and.
     &          dum.ne.'3'.and.dum.ne.'4'.and.
     &          dum.ne.'5'.and.dum.ne.'6'.and.
     &          dum.ne.'7'.and.dum.ne.'8'.and.
     &          dum.ne.'9'.and.dum.ne.'.'.and.
     &          dum.ne.'+'.and.dum.ne.'-'.and.
     &          dum.ne.'D'.and.dum.ne.'d'.and.
     &          dum.ne.'E'.and.dum.ne.'e'

      return
      end

************************************************************************
*                                                                      *
      function dnen4(dum)
*                                                                      *
************************************************************************
      character dum*1
      logical dnen4
*-----------------------------------------------------------------------

      dnen4  =  dum.ne.'0'.and.
     &          dum.ne.'1'.and.dum.ne.'2'.and.
     &          dum.ne.'3'.and.dum.ne.'4'.and.
     &          dum.ne.'5'.and.dum.ne.'6'.and.
     &          dum.ne.'7'.and.dum.ne.'8'.and.
     &          dum.ne.'9'.and.
     &          dum.ne.'+'.and.dum.ne.'-'

      return
      end

************************************************************************
*                                                                      *
      function deqn0(dum)
*                                                                      *
************************************************************************
      character dum*1
      logical deqn0
*-----------------------------------------------------------------------

      deqn0  =  dum.eq.'0'.or.
     &          dum.eq.'1'.or.dum.eq.'2'.or.
     &          dum.eq.'3'.or.dum.eq.'4'.or.
     &          dum.eq.'5'.or.dum.eq.'6'.or.
     &          dum.eq.'7'.or.dum.eq.'8'.or.
     &          dum.eq.'9'.or.dum.eq.'.'
      return
      end

************************************************************************
*                                                                      *
      function deqn1(dum)
*                                                                      *
************************************************************************
      character dum*1
      logical deqn1
*-----------------------------------------------------------------------

      deqn1  =  dum.eq.'0'.or.
     &          dum.eq.'1'.or.dum.eq.'2'.or.
     &          dum.eq.'3'.or.dum.eq.'4'.or.
     &          dum.eq.'5'.or.dum.eq.'6'.or.
     &          dum.eq.'7'.or.dum.eq.'8'.or.
     &          dum.eq.'9'
      return
      end

************************************************************************
*                                                                      *
      function deqn2(dum)
*                                                                      *
************************************************************************
      character dum*1
      logical deqn2
*-----------------------------------------------------------------------

      deqn2  =  dum.eq.'0'.or.
     &          dum.eq.'1'.or.dum.eq.'2'.or.
     &          dum.eq.'3'.or.dum.eq.'4'.or.
     &          dum.eq.'5'.or.dum.eq.'6'.or.
     &          dum.eq.'7'.or.dum.eq.'8'.or.
     &          dum.eq.'9'.or.dum.eq.'.'.or.
     &          dum.eq.'D'.or.dum.eq.'d'.or.
     &          dum.eq.'E'.or.dum.eq.'e'
      return
      end

************************************************************************
*                                                                      *
      function deqn3(dum)
*                                                                      *
************************************************************************
      character dum*1
      logical deqn3
*-----------------------------------------------------------------------

      deqn3  =  dum.eq.'0'.or.
     &          dum.eq.'1'.or.dum.eq.'2'.or.
     &          dum.eq.'3'.or.dum.eq.'4'.or.
     &          dum.eq.'5'.or.dum.eq.'6'.or.
     &          dum.eq.'7'.or.dum.eq.'8'.or.
     &          dum.eq.'9'.or.dum.eq.'.'.or.
     &          dum.eq.'+'.or.dum.eq.'-'.or.
     &          dum.eq.'D'.or.dum.eq.'d'.or.
     &          dum.eq.'E'.or.dum.eq.'e'
      return
      end

************************************************************************
*                                                                      *
      function deqn4(dum)
*                                                                      *
************************************************************************
      character dum*1
      logical deqn4
*-----------------------------------------------------------------------

      deqn4  =  dum.eq.'0'.or.
     &          dum.eq.'1'.or.dum.eq.'2'.or.
     &          dum.eq.'3'.or.dum.eq.'4'.or.
     &          dum.eq.'5'.or.dum.eq.'6'.or.
     &          dum.eq.'7'.or.dum.eq.'8'.or.
     &          dum.eq.'9'.or.
     &          dum.eq.'+'.or.dum.eq.'-'
      return
      end

************************************************************************
*                                                                      *
      function deqn5(dum)
*                                                                      *
************************************************************************
      character dum*1
      logical deqn5
*-----------------------------------------------------------------------

      deqn5  =  dum.eq.'0'.or.
     &          dum.eq.'1'.or.dum.eq.'2'.or.
     &          dum.eq.'3'.or.dum.eq.'4'.or.
     &          dum.eq.'5'.or.dum.eq.'6'.or.
     &          dum.eq.'7'.or.dum.eq.'8'.or.
     &          dum.eq.'9'.or.dum.eq.'.'.or.
     &          dum.eq.'+'.or.dum.eq.'-'
      return
      end

************************************************************************
*                                                                      *
      function dcom2(dum)
*                                                                      *
*     separate characters : ' ', '	', (, [, {                     *
*                                                                      *
************************************************************************
      character dum*1
      logical dcom2
      character tub*1
      tub = char(9)
*-----------------------------------------------------------------------

      dcom2  =  dum.eq.' '.or.dum.eq.tub.or.
     &          dum.eq.'('.or.dum.eq.'{'.or.
     &          dum.eq.'['
      return
      end

************************************************************************
*                                                                      *
      function dcom3(dum) ! T.Sato 2021/02/12 for reading multiplier subsection
*                                                                      *
*     find logical characters used in multiplier subsection
*                                                                      *
************************************************************************
      character dum*1
      logical dcom3
      character tub*1
      tub = char(9)
*-----------------------------------------------------------------------

      dcom3  =  dum.eq.' '.or.dum.eq.tub.or.
     &          dum.eq.'('.or.dum.eq.')'.or.
     &          dum.eq.':'

      return
      end

************************************************************************
*                                                                      *
      block data colordf
*                                                                      *
*              definition of colors                                    *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter ( icnm = 36 )
      parameter ( iccx = 15 )
*                 icnm : max number color definition
*                 iccx : max character length of color name

      common /coldef/ coldf(3,icnm), lcoln(icnm), colnm(icnm), icolr(30)
      character colnm*15

*-----------------------------------------------------------------------

      data icolr / 11,  5, 23, 26, 21,
     &              1, 15,  8, 20,  6,
     &             24, 10, 18,  2, 16,
     &             25,  9,  4, 12, 22,
     &             17, 29, 13,  3, 27,
     &              7, 14, 30, 19, 28/

      data ( colnm(j),
     &     ( coldf(i,j), i = 1, 3 ), lcoln(j), j = 1, icnm ) /
c          123456789012345
     &    'darkred        ', 1.000,1.000,0.600, 7, !  1
     &    'pink           ', 1.000,0.500,1.000, 4, !  2
     &    'red            ', 1.000,1.000,1.000, 3, !  3
     &    'pastelpink     ', 0.900,0.500,1.000,10, !  4
     &    'orangeyellow   ', 0.867,1.000,1.000,12, !  5
     &    'orange         ', 0.933,1.000,1.000, 6, !  6
     &    'darkbrown      ', 0.900,1.000,0.300, 9, !  7
     &    'pastelbrown    ', 0.900,0.600,0.500,11, !  8
     &    'brown          ', 0.900,1.000,0.500, 5, !  9
     &    'camel          ', 0.800,0.700,0.700, 5, ! 10
     &    'yellowgreen    ', 0.700,1.000,1.000,11, ! 11
     &    'yellow         ', 0.800,1.000,1.000, 6, ! 12
     &    'darkgreen      ', 0.600,1.000,0.500, 9, ! 13
     &    'mossgreen      ', 0.500,1.000,0.300, 9, ! 14
     &    'green          ', 0.600,1.000,1.000, 5, ! 15
     &    'pastelyellow   ', 0.800,0.700,1.000,12, ! 16
     &    'pastelgreen    ', 0.700,0.600,1.000,11, ! 17
     &    'pastelcyan     ', 0.400,0.400,1.000,10, ! 18
     &    'cyanblue       ', 0.400,1.000,0.500, 8, ! 19
     &    'cyan           ', 0.400,1.000,1.000, 4, ! 20
     &    'bluegreen      ', 0.500,1.000,1.000, 9, ! 21
     &    'blue           ', 0.200,1.000,1.000, 4, ! 22
     &    'pastelblue     ', 0.250,0.400,1.000,10, ! 23
     &    'pastelviolet   ', 0.133,0.400,1.000,12, ! 24
     &    'pastelpurple   ', 0.100,0.400,0.500,12, ! 25
     &    'pastelmagenta  ', 0.067,0.600,1.000,13, ! 26
     &    'violet         ', 0.133,1.000,1.000, 6, ! 27
     &    'purple         ', 0.100,1.000,0.500, 6, ! 28
     &    'magenta        ', 0.067,1.000,1.000, 7, ! 29
     &    'winered        ', 0.002,0.800,0.700, 7, ! 30
c
     &    'matblack       ',-0.200,1.000,1.000, 8,
     &    'black          ',-0.000,1.000,1.000, 5,
     &    'darkgray       ',-0.400,1.000,1.000, 8,
     &    'gray           ',-0.600,1.000,1.000, 4,
     &    'lightgray      ',-0.800,1.000,1.000, 9,
     &    'white          ',-1.000,1.000,1.000, 5/
c          123456789012345

      end

************************************************************************
*                                                                      *
      subroutine dcols(icc,lum,ic,tcol,ierr,iclm,c1,c2)
*                                                                      *
*      PURPOSE  :   READ COLOR AND GRAY SACLE                          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character lum(ichrl)*1
      character c1*1, c2*1
      character clnam*15
      logical deqn1,deqn3
      dimension tcol(3)

*-----------------------------------------------------------------------

      parameter ( icnm = 36 )
      parameter ( iccx = 15 )

      common /coldef/ coldf(3,icnm), lcoln(icnm), colnm(icnm), icolr(30)
      character colnm*15

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

            ierr = 0
            icol = 0
            icob = 0
            icoc = 0

            tcol(1) = -r1max
            tcol(2) = 1.0
            tcol(3) = 1.0

*-----------------------------------------------------------------------

            ic = ic + 1

            if( lum(ic) .ne. c1 ) goto 999

  711       ic = ic + 1

            if( ic .gt. iclm ) goto 999

            if( lum(ic) .eq. ' ' .or. lum(ic) .eq . tub ) goto 711
            if( lum(ic) .eq. c2 ) goto 730

*-----------------------------------------------------------------------
*        color by numeric
*-----------------------------------------------------------------------

            if( deqn1(lum(ic)) .or. lum(ic) .eq. '.' .or.
     &          lum(ic) .eq. '-' .or. lum(ic) .eq. '+' ) then

               if( icol .ne. 0 .or. icoc .ne. 0 ) goto 999

*-----------------------------------------------------------------------

                  ici = ic
  710             ic  = ic + 1

               if( ic .gt. iclm ) goto 999
               if( deqn3(lum(ic)) ) goto 710

                     icob = icob + 1
                     if( icob .gt. 3 ) goto 999

                     ic  = ic - 1
                     icf = ic

                     call rnum(tcol0,lum,ici,icf,ierr)
                        if(ierr.ne.0) goto 999

                        if( tcol0 .le. 0.0 ) then

                           if( icob .eq. 1 ) then
                              tcol1 = max(-1.0d0,tcol0)
                              tcol2 = min(-0.0d0,tcol1)
                              tcol(1) = -tcol2 - 2
                           else
                              tcol(icob) = 0.0
                           end if

                        else if( tcol0 .gt. 0.0 ) then
                              tcol1 = max(0.0d0,tcol0)
C                             Approximately convert to standard HSB hue.
                              if( icob .eq. 1 .and. 1d0 .lt. tcol1 .and.
     &                                            tcol1 .le. 360d0) then
                                tcol1=(1d0-tcol1/360d0)**1.4
                              endif
                              tcol2 = min(1.0d0,tcol1)
                           if( icob .eq. 1 ) then
                              tcol(1) = 2.5 * ( 1.0 - tcol2 ) + 1.0
                           else
                              tcol(icob) = tcol2
                           end if

                        end if

*-----------------------------------------------------------------------
*        color by symbol
*-----------------------------------------------------------------------

            else

                  if( icob .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*              color or gray scal by color name
*-----------------------------------------------------------------------

                  do k = 1, iccx

                     clnam(k:k) = lum(ic+k-1)

                  end do

                  do i = 1, icnm

                     if( clnam(1:lcoln(i)) .eq. colnm(i)(1:lcoln(i)) )
     &               goto 100

                  end do

                     goto 200

  100             continue

                     ipm = i

                     icoc = icoc + 1

                  if( coldf(1,ipm) .gt. 0.0 ) then

                     tcol(1) = 2.5 * ( 1.0 - coldf(1,ipm) ) + 1.0

                  else

                     tcol(1) = -coldf(1,ipm) - 2.0

                  end if

                     tcol(2) =  coldf(2,ipm)
                     tcol(3) =  coldf(3,ipm)

                     ic = ic + lcoln(ipm) - 1

                     goto 711

*-----------------------------------------------------------------------
*              color or gray scal by symbol
*-----------------------------------------------------------------------

  200             continue

                  if( icoc .ne. 0 ) goto 999
                  icol = icol + 1

*-----------------------------------------------------------------------
*                 symbl color
*-----------------------------------------------------------------------

                  if(lum(ic).eq.'r') then
                     if( tcol(1) .gt. 0.0 ) then
                         tcol(1) = tcol(1) + 0.16667
                     else
                         tcol(1) = 1.0
                     end if
                  else if(lum(ic).eq.'y') then
                     if( tcol(1) .gt. 0.0 ) then
                         tcol(1) = tcol(1) + 0.16667
                     else
                         tcol(1) = 1.5
                     end if
                  else if(lum(ic).eq.'g') then
                     if( tcol(1) .gt. 0.0 ) then
                         tcol(1) = tcol(1) + 0.16667
                     else
                         tcol(1) = 2.0
                     end if
                  else if(lum(ic).eq.'c') then
                     if( tcol(1) .gt. 0.0 ) then
                         tcol(1) = tcol(1) + 0.16667
                     else
                         tcol(1) = 2.5
                     end if
                  else if(lum(ic).eq.'b') then
                     if( tcol(1) .gt. 0.0 ) then
                         tcol(1) = tcol(1) + 0.16667
                     else
                         tcol(1) = 3.0
                     end if

*-----------------------------------------------------------------------
*                 gray scal
*-----------------------------------------------------------------------

                  else if(lum(ic).eq.'e') then
                     if( tcol(1) .gt. -r0max .and.
     &                   tcol(1) .le. 0.0 ) then
                        tcol(1) = tcol(1) - 0.05
                     else
                        tcol(1) = -2.0
                     end if
                  else if(lum(ic).eq.'f') then
                     if( tcol(1) .gt. -r0max .and.
     &                   tcol(1) .le. 0.0 ) then
                        tcol(1) = tcol(1) - 0.05
                     else
                        tcol(1) = -1.8
                     end if
                  else if(lum(ic).eq.'j') then
                     if( tcol(1) .gt. -r0max .and.
     &                   tcol(1) .le. 0.0 ) then
                        tcol(1) = tcol(1) - 0.05
                     else
                        tcol(1) = -1.6
                     end if
                  else if(lum(ic).eq.'k') then
                     if( tcol(1) .gt. -r0max .and.
     &                   tcol(1) .le. 0.0 ) then
                        tcol(1) = tcol(1) - 0.05
                     else
                        tcol(1) = -1.4
                     end if
                  else if(lum(ic).eq.'o') then
                     if( tcol(1) .gt. -r0max .and.
     &                  tcol(1) .le. 0.0 ) then
                        tcol(1) = tcol(1) - 0.05
                     else
                        tcol(1) = -1.2
                     end if
                  else if(lum(ic).eq.'w') then
                     if( tcol(1) .gt. -r0max .and.
     &                   tcol(1) .le. 0.0 ) then
                        tcol(1) = tcol(1) - 0.05
                     else
                        tcol(1) = -1.0
                     end if
                  else
                     goto 999
                  end if

*-----------------------------------------------------------------------

            end if

                  goto 711

*-----------------------------------------------------------------------

  730       ic = ic + 1

            if( icol .eq. 0 .and. icob .eq. 0 .and. icoc .eq. 0 )
     &      goto 999
            if( icc .eq. 1 ) then
c T.Sato for debug mode
                  tmp=tcol(1)
                  tcol(1) = chue(tmp)
            end if

      return

*-----------------------------------------------------------------------

  999 ierr = 1

      return
      end

************************************************************************
*                                                                      *
      function chue(rin0)
*                                                                      *
*     TAKEN FROM SCIPLOT FOR COLOR DEFINITION   1-2                    *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

         if( rin0 .le. 0.0 ) then

            rin = max(-2.0d0,rin0)
            rin = min(-1.0d0,rin)

            chue = sqrt( rin + 2.0 ) - 2.0

         else

            rin = max(1.0d0,rin0)
            rin = min(3.49d0,rin)

            if( rin .lt. 3.5000 ) rut  =   (rin - 2.25)**2 + 39.0/16.0
            if( rin .lt. 2.5000 ) rut  = - (rin - 2.75)**2 + 41.0/16.0
            if( rin .lt. 2.0000 ) rut  =   (rin - 1.25)**2 + 23.0/16.0
            if( rin .lt. 1.5000 ) rut  = - (rin - 1.75)**2 + 25.0/16.0

            chue = 2.0 * ( rut - 1.0 ) / 2.0 / 3.0  +  1.0

         end if

      return
      end

************************************************************************
*                                                                      *
      subroutine ctomo(rcol,clmo)
*                                                                      *
*       color to mono                                                  *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension rcol(3)

*-----------------------------------------------------------------------

         redu = 0.2

         rcol(1) = ( abs( rcol(1) - 1.0 ) )**0.4 * ( 1.0 - redu )
     &           + redu

         bri = ( 1.0 - rcol(1) ) * ( 1.0 - rcol(2) ) * 1.0
         drk = ( rcol(1) - rcol(1) * rcol(3) ) * 0.8

         rcol(1) = min( 1.0d0, max( 0.0d0, rcol(1) + bri - drk ) )

         rcol(1) = rcol(1)**clmo

         rcol(1) = rcol(1) - 2.0

*-----------------------------------------------------------------------

      return
      end

