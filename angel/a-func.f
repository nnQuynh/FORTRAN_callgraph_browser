************************************************************************
*                                                                      *
      subroutine func01(dum,ici,iclm,ierr,duml,
     &                  rk,kd,kp,ki,jd,mg,md,me,mp,id,iaac)
*                                                                      *
*      purpose  : read function description                            *
*                                                                      *
*      variables :                                                     *
*             id : element index                                       *
*             ip : class index                                         *
*                                                                      *
*         kp(id) : class of id element                                 *
*         ki(id) : int. func. index or column number                   *
*                                                                      *
*             kind of elements  kd(i) jd(i)                            *
*             'X'         :        1     0                             *
*             'Y'         :        2     0                             *
*             'C'         :        3     0                             *
*             number      :        4     0                             *
*             '*'         :        5     2                             *
*             '/'         :        6     2                             *
*             '+'         :        7     3                             *
*             '-'         :        8     3                             *
*             '**'        :        9     1                             *
*             int. func.  :       10    -2                             *
*             '('         :       11    -1                             *
*             ')'         :       12    -1                             *
*             ','         :       13    -3                             *
*             'A'         :       14     0                             *
*             'V'         :       15     0                             *
*             'DY'        :       16     0                             *
*             'DX'        :       17     0                             *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character duml*1

      dimension rk(0:mfc)
      dimension kd(0:mfc), kp(0:mfc), ki(0:mfc), jd(0:mfc)
      dimension md(0:mfc), mg(0:mfc), me(0:mfc), mp(0:mfc)

      dimension mq(0:mfc)
      dimension lp(0:mfc)
      dimension jp(0:mfc)

      logical deqn0
      logical deqn4
      logical func02

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------
*     INITIALIZATION
*-----------------------------------------------------------------------

            ierr = 0

            iaac = 0

            ip = 0
            id = 0

         do 10 i = 0, mfc

            md(i) = 0
            mg(i) = 0
            me(i) = 0
            mp(i) = 0

            kd(i) = 0
            kp(i) = 0
            ki(i) = 0
            jd(i) = -5

            mq(i) = 0
            lp(i) = 0
            jp(i) = 0

   10    continue

*-----------------------------------------------------------------------
*     READ TEXT
*-----------------------------------------------------------------------

            ic = ici

  100       ic = ic + 1

               if( ic .gt. iclm ) goto 999

               if( dum(ic) .eq. ' ' .or. dum(ic) .eq. tub ) goto 100

*-----------------------------------------------------------------------


               if( dum(ic) .eq. duml ) then

                  ici = ic

                  goto 980

               end if

*-----------------------------------------------------------------------
*        START
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*     ---- INTRINSIC FUNCTION ----
*-----------------------------------------------------------------------

         if( func02(dum,ic,ifn,imc) ) then

               if( imc .eq. 1 ) then

                        mp(id+1) = mp(id) + 1
                        mg(id+1) = 1
                        md(id+1) = id + 1
                        mq(mp(id+1)) = id + 1

                  if( mp(id) .ge. 1 ) then

                     if( kp(id) .eq. kp(md(id)) ) then

                        ip = ip + 1

                     end if

                  end if

               else

                  if( mp(id) .ge. 1 ) then

                        mp(id+1) = mp(id)
                        md(id+1) = md(id)
                        mg(id+1) = mg(id)
                        mq(mp(id+1)) = id + 1

                     if( kp(id) .eq. kp(md(id)) ) then

                        ip = ip + 1

                     end if

                  end if

               end if


               ip = ip + 1
               id = id + 1

               kd(id) = 10
               kp(id) = ip

               lp(ip) = id
               jp(id) = lp(ip)

               ki(id) = ifn


               jd(id) = -2

*-----------------------------------------------------------------------

         else if( dum(ic) .eq. '(' ) then

            if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               if( kp(id) .eq. kp(md(id)) ) then

                  ip = ip + 1
                  lp(ip) = id + 1

               end if

            end if

               ip = ip + 1
               id = id + 1

               kd(id) = 11
               kp(id) = ip

               lp(ip) = id
               jp(id) = lp(ip)

               jd(id) = -1

*-----------------------------------------------------------------------

         else if( dum(ic) .eq. ')' ) then


            if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)

                  mq(mp(id+1)) = id + 1

               if( md(id) .eq. lp(ip-1) ) then

                  ip = ip - 1

               end if

            end if

               id = id + 1

               kd(id) = 12
               kp(id) = ip
               jp(id) = lp(ip)

               ip = ip - 1

               if( ip .lt. 0 ) goto 999

            if( mg(id) .ge. 1 .and.
     &          md(id) .eq. jp(id) ) then

               if( mg(id) .lt. 2 ) goto 999

               if( ( ki(md(id)) .eq.  9 .or.
     &               ki(md(id)) .eq. 11 .or.
     &               ki(md(id)) .eq. 19 ) .and.
     &               mg(id) .gt. 2 ) goto 999


               me(md(id)) = mg(id)


               mp(id) = max(0,mp(id)-1)

               mg(id) = mg(mq(mp(id)))
               md(id) = md(mq(mp(id)))

            end if

               jd(id) = -1

*-----------------------------------------------------------------------

         else if( dum(ic) .eq. ',' ) then

            if( mg(id) .lt. 1 ) goto 999

               mp(id+1) = mp(id)
               md(id+1) = md(id)
               mg(id+1) = mg(id) + 1
               mq(mp(id+1)) = id + 1

               ip = ip - 1

               id = id + 1

               kd(id) = 13
               kp(id) = ip

               jp(id) = lp(ip)

               jd(id) = -3

*-----------------------------------------------------------------------

         else if( dum(ic) .eq. '*' .and. dum(ic+1) .ne. '*' ) then

               if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               end if

               id = id + 1

               kd(id) = 5
               kp(id) = ip

               jp(id) = lp(ip)

               jd(id) = 2

*-----------------------------------------------------------------------

         else if( dum(ic) .eq. '/' ) then

               if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               end if

               id = id + 1

               kd(id) = 6
               kp(id) = ip

               jp(id) = lp(ip)

               jd(id) = 2

*-----------------------------------------------------------------------

         else if( dum(ic) .eq. '+' ) then

            if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               if( kp(id) .eq. kp(md(id)) ) then

                  ip = ip + 1
                  lp(ip) = id + 1

               end if

            end if

               id = id + 1

               kd(id) = 7
               kp(id) = ip

               jp(id) = lp(ip)

               jd(id) = 3

*-----------------------------------------------------------------------

         else if( dum(ic) .eq. '-' ) then

            if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               if( kp(id) .eq. kp(md(id)) ) then

                  ip = ip + 1
                  lp(ip) = id + 1

               end if

            end if

               id = id + 1

               kd(id) = 8
               kp(id) = ip

               jp(id) = lp(ip)

               jd(id) = 3

*-----------------------------------------------------------------------

         else if( dum(ic) .eq. '*' .and. dum(ic+1) .eq. '*' ) then

               if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               end if

               id = id + 1

               kd(id) = 9
               kp(id) = ip

               jp(id) = lp(ip)

               ic = ic + 1

               jd(id) = 1

*-----------------------------------------------------------------------
*     ---- V-VALUES ----
*-----------------------------------------------------------------------

         else if( dum(ic) .eq. 'V' .or. dum(ic) .eq. 'v' ) then

            if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               if( kp(id) .eq. kp(md(id)) ) then

                  ip = ip + 1
                  lp(ip) = id + 1

               end if

            end if

               id = id + 1

               kd(id) = 15
               kp(id) = ip

               jp(id) = lp(ip)

               jd(id) = 0

*-----------------------------------------------------------------------
*     ---- X ----
*-----------------------------------------------------------------------

         else if( dum(ic) .eq. 'X' .or. dum(ic) .eq. 'x' ) then

            if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               if( kp(id) .eq. kp(md(id)) ) then

                  ip = ip + 1
                  lp(ip) = id + 1

               end if

            end if

               id = id + 1

               kd(id) = 1
               kp(id) = ip

               jp(id) = lp(ip)

               jd(id) = 0

*-----------------------------------------------------------------------
*     ---- DX-VALUES ----
*-----------------------------------------------------------------------

         else if( ( dum(ic) .eq. 'D' .or. dum(ic) .eq. 'd' ) .and.
     &            ( dum(ic+1) .eq. 'X' .or. dum(ic+1) .eq. 'x' )
     &            ) then

            if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               if( kp(id) .eq. kp(md(id)) ) then

                  ip = ip + 1
                  lp(ip) = id + 1

               end if

            end if

               id = id + 1

               kd(id) = 17
               kp(id) = ip

               jp(id) = lp(ip)

               ki(id) = mh

               icc = ic + 2

            if( dum(icc) .eq. 'P' .or. dum(icc) .eq. 'p' ) then

                  ic = icc

            else if( dum(icc) .eq. 'M' .or. dum(icc) .eq. 'm' ) then

                  ki(id) = - ki(id)

                  ic = icc

            else

                  ic = icc - 1

            end if

               jd(id) = 0

*-----------------------------------------------------------------------
*     ---- Y VALUABLE SECTION ----
*-----------------------------------------------------------------------

         else if( dum(ic) .eq. 'Y' .or. dum(ic) .eq. 'y' ) then

            if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               if( kp(id) .eq. kp(md(id)) ) then

                  ip = ip + 1
                  lp(ip) = id + 1

               end if

            end if

               id = id + 1

               kd(id) = 2
               kp(id) = ip

               jp(id) = lp(ip)

               icc = ic + 1

            call idd2(dum,icc,iclm,idd,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( idd .lt. 0 ) goto 999

                  ki(id) = idd


               ic = icc - 1

               jd(id) = 0

*-----------------------------------------------------------------------
*     ---- DY-VALUES ----
*-----------------------------------------------------------------------

         else if( ( dum(ic) .eq. 'D' .or. dum(ic) .eq. 'd' ) .and.
     &            ( dum(ic+1) .eq. 'Y' .or. dum(ic+1) .eq. 'y' )
     &            ) then

            if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               if( kp(id) .eq. kp(md(id)) ) then

                  ip = ip + 1
                  lp(ip) = id + 1

               end if

            end if

               id = id + 1

               kd(id) = 16
               kp(id) = ip

               jp(id) = lp(ip)

               icc = ic + 2

            call idd2(dum,icc,iclm,idd,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( idd .lt. 0 ) goto 999

               ki(id) = idd


            if( dum(icc) .eq. 'P' .or. dum(icc) .eq. 'p' ) then

                  ic = icc

            else if( dum(icc) .eq. 'M' .or. dum(icc) .eq. 'm' ) then

                  ki(id) = - ki(id)

                  ic = icc

            else

                  ic = icc - 1

            end if

               jd(id) = 0

*-----------------------------------------------------------------------
*     ---- CONSTANTS ----
*-----------------------------------------------------------------------

         else if( dum(ic) .eq. 'C' .or. dum(ic) .eq. 'c' ) then

            if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               if( kp(id) .eq. kp(md(id)) ) then

                  ip = ip + 1
                  lp(ip) = id + 1

               end if

            end if

               id = id + 1

               kd(id) = 3
               kp(id) = ip

               jp(id) = lp(ip)

               icc = ic + 1

               call idd2(dum,icc,iclm,idd,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( idd .lt. 0 ) goto 999

               ki(id) = idd

               ic = icc - 1

               jd(id) = 0

*-----------------------------------------------------------------------
*     ---- PARAMETERS ----
*-----------------------------------------------------------------------

         else if( dum(ic) .eq. 'A' .or. dum(ic) .eq. 'a' ) then

            if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               if( kp(id) .eq. kp(md(id)) ) then

                  ip = ip + 1
                  lp(ip) = id + 1

               end if

            end if

               id = id + 1

               kd(id) = 14
               kp(id) = ip

               jp(id) = lp(ip)

               icc = ic + 1

               call idd2(dum,icc,iclm,idd,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( idd .lt. 0 ) goto 999

               ki(id) = idd

               ic = icc - 1

               jd(id) = 0

               iaac = 1

*-----------------------------------------------------------------------
*     ---- REAL NUMBER, PI AND RANDOM NUMBER RAN ----
*-----------------------------------------------------------------------

         else if( deqn0(dum(ic)) .or.
     &          ( ( dum(ic) .eq. 'P' .or. dum(ic) .eq. 'p' ) .and.
     &            ( dum(ic+1) .eq. 'I' .or.
     &              dum(ic+1) .eq. 'i' ) ) .or.
     &          ( ( dum(ic) .eq. 'R' .or. dum(ic) .eq. 'r' ) .and.
     &            ( dum(ic+1) .eq. 'A' .or.
     &              dum(ic+1) .eq. 'a' ) .and.
     &            ( dum(ic+2) .eq. 'N' .or.
     &              dum(ic+2) .eq. 'n' ) ) ) then

            if( mp(id) .ge. 1 ) then

                  mp(id+1) = mp(id)
                  md(id+1) = md(id)
                  mg(id+1) = mg(id)
                  mq(mp(id+1)) = id + 1

               if( kp(id) .eq. kp(md(id)) ) then

                  ip = ip + 1
                  lp(ip) = id + 1

               end if

            end if

               id = id + 1

               kd(id) = 4
               kp(id) = ip

               jp(id) = lp(ip)

               jd(id) = 0

*-----------------------------------------------------------------------
*           ---- REAL NUMBER ----
*-----------------------------------------------------------------------

            if( deqn0(dum(ic)) ) then

                  ici = ic

                  ice = 0

  200          ic = ic + 1

                  if( ic .gt. iclm ) goto 999

                  if( dum(ic) .eq. duml ) goto 300

                  if( ice .eq. 0 .and. deqn0(dum(ic)) ) goto 200

                  if( ice .eq. 1 .and. deqn4(dum(ic)) ) then

                     ice = 0

                     goto 200

                  end if

                  if( dum(ic) .eq. 'E' .or. dum(ic) .eq. 'e' .or.
     &                dum(ic) .eq. 'D' .or. dum(ic) .eq. 'd' ) then

                     if( ice .ne. 0 ) goto 999

                     ice = 1

                     goto 200

                  end if


  300          icf = ic - 1

                  call rnum(rk(id),dum,ici,icf,ierr)

                  if( ierr .ne. 0 ) goto 999

               ic = icf

*-----------------------------------------------------------------------
*           ---- PI ----
*-----------------------------------------------------------------------

            else if( dum(ic) .eq. 'P' .or. dum(ic) .eq. 'p' ) then

               ic = ic + 1

               rk(id) = pi

*-----------------------------------------------------------------------
*           ---- RANDOM NUMBER ----
*-----------------------------------------------------------------------

            else if( dum(ic) .eq. 'R' .or. dum(ic) .eq. 'r' ) then

               ic = ic + 2

               rk(id) = r1max

            end if


*-----------------------------------------------------------------------
*     ---- ERROR ----
*-----------------------------------------------------------------------

         else

            goto 999

         end if

*-----------------------------------------------------------------------

            goto 100

*-----------------------------------------------------------------------
*     LAST CHACK
*-----------------------------------------------------------------------

  980 continue


      return

*-----------------------------------------------------------------------
*     ERROR END
*-----------------------------------------------------------------------

  999 ierr = 1

      return
      end


************************************************************************
*                                                                      *
      function func02(dum,ic,ifn,imc)
*                                                                      *
*       PURPOSE:       IDENTIFY THE INTRINSIC FUNCTION                 *
*                                                                      *
*       INTTRINSIC FUNCTION    ID                                      *
*                                                                      *
*           FLOAT          ;   1                                       *
*           DBLE           ;   1                                       *
*           INT            ;   2                                       *
*           ABS            ;   3                                       *
*           EXP            ;   4                                       *
*           LOG            ;   5                                       *
*           LOG10          ;   6                                       *
*           MAX( , , , )   ;   7                                       *
*           MIN( , , , )   ;   8                                       *
*           MOD( , )       ;   9                                       *
*           NINT           ;   10                                      *
*           SIGN( , )      ;   11                                      *
*           SQRT           ;   12                                      *
*           ACOS           ;   13                                      *
*           ASIN           ;   14                                      *
*           ATAN           ;   15                                      *
*        X  ACOSH          ;   16                                      *
*        X  ASINH          ;   17                                      *
*        X  ATANH          ;   18                                      *
*           ATAN2( , )     ;   19                                      *
*           COS            ;   20                                      *
*           COSH           ;   21                                      *
*           SIN            ;   22                                      *
*           SINH           ;   23                                      *
*           TAN            ;   24                                      *
*           TANH           ;   25                                      *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      logical func02

*-----------------------------------------------------------------------

            func02 = .true.

            ifn = 0
            imc = 0

*-----------------------------------------------------------------------

         if(
     &       ( dum(ic  ) .eq. 'F' .or. dum(ic  ) .eq. 'f' ) .and.
     &       ( dum(ic+1) .eq. 'L' .or. dum(ic+1) .eq. 'l' ) .and.
     &       ( dum(ic+2) .eq. 'O' .or. dum(ic+2) .eq. 'o' ) .and.
     &       ( dum(ic+3) .eq. 'A' .or. dum(ic+3) .eq. 'a' ) .and.
     &       ( dum(ic+4) .eq. 'T' .or. dum(ic+4) .eq. 't' ) .and.
     &         dum(ic+5) .eq. '(' ) then

               ifn = 1
               ic  = ic + 5

         else if(
     &       ( dum(ic  ) .eq. 'D' .or. dum(ic  ) .eq. 'd' ) .and.
     &       ( dum(ic+1) .eq. 'B' .or. dum(ic+1) .eq. 'b' ) .and.
     &       ( dum(ic+2) .eq. 'L' .or. dum(ic+2) .eq. 'l' ) .and.
     &       ( dum(ic+3) .eq. 'E' .or. dum(ic+3) .eq. 'e' ) .and.
     &         dum(ic+4) .eq. '(' ) then

               ifn = 1
               ic  = ic + 4

         else if(
     &       ( dum(ic  ) .eq. 'I' .or. dum(ic  ) .eq. 'i' ) .and.
     &       ( dum(ic+1) .eq. 'N' .or. dum(ic+1) .eq. 'n' ) .and.
     &       ( dum(ic+2) .eq. 'T' .or. dum(ic+2) .eq. 't' ) .and.
     &         dum(ic+3) .eq. '(' ) then

               ifn = 2
               ic  = ic + 3


         else if(
     &       ( dum(ic  ) .eq. 'A' .or. dum(ic  ) .eq. 'a' ) .and.
     &       ( dum(ic+1) .eq. 'B' .or. dum(ic+1) .eq. 'b' ) .and.
     &       ( dum(ic+2) .eq. 'S' .or. dum(ic+2) .eq. 's' ) .and.
     &         dum(ic+3) .eq. '(' ) then

               ifn = 3
               ic  = ic + 3


         else if(
     &       ( dum(ic  ) .eq. 'E' .or. dum(ic  ) .eq. 'e' ) .and.
     &       ( dum(ic+1) .eq. 'X' .or. dum(ic+1) .eq. 'x' ) .and.
     &       ( dum(ic+2) .eq. 'P' .or. dum(ic+2) .eq. 'p' ) .and.
     &         dum(ic+3) .eq. '(' ) then

               ifn = 4
               ic  = ic + 3


         else if(
     &       ( dum(ic  ) .eq. 'L' .or. dum(ic  ) .eq. 'l' ) .and.
     &       ( dum(ic+1) .eq. 'O' .or. dum(ic+1) .eq. 'o' ) .and.
     &       ( dum(ic+2) .eq. 'G' .or. dum(ic+2) .eq. 'g' ) .and.
     &         dum(ic+3) .eq. '(' ) then

               ifn = 5
               ic  = ic + 3


         else if(
     &       ( dum(ic  ) .eq. 'L' .or. dum(ic  ) .eq. 'l' ) .and.
     &       ( dum(ic+1) .eq. 'O' .or. dum(ic+1) .eq. 'o' ) .and.
     &       ( dum(ic+2) .eq. 'G' .or. dum(ic+2) .eq. 'g' ) .and.
     &       ( dum(ic+3) .eq. '1' ) .and.
     &       ( dum(ic+4) .eq. '0' ) .and.
     &         dum(ic+5) .eq. '(' ) then

               ifn = 6
               ic  = ic + 5


         else if(
     &       ( dum(ic  ) .eq. 'M' .or. dum(ic  ) .eq. 'm' ) .and.
     &       ( dum(ic+1) .eq. 'A' .or. dum(ic+1) .eq. 'a' ) .and.
     &       ( dum(ic+2) .eq. 'X' .or. dum(ic+2) .eq. 'x' ) .and.
     &         dum(ic+3) .eq. '(' ) then

               ifn = 7
               imc = 1
               ic  = ic + 3


         else if(
     &       ( dum(ic  ) .eq. 'M' .or. dum(ic  ) .eq. 'm' ) .and.
     &       ( dum(ic+1) .eq. 'I' .or. dum(ic+1) .eq. 'i' ) .and.
     &       ( dum(ic+2) .eq. 'N' .or. dum(ic+2) .eq. 'n' ) .and.
     &         dum(ic+3) .eq. '(' ) then

               ifn = 8
               imc = 1
               ic  = ic + 3


         else if(
     &       ( dum(ic  ) .eq. 'M' .or. dum(ic  ) .eq. 'm' ) .and.
     &       ( dum(ic+1) .eq. 'O' .or. dum(ic+1) .eq. 'o' ) .and.
     &       ( dum(ic+2) .eq. 'D' .or. dum(ic+2) .eq. 'd' ) .and.
     &         dum(ic+3) .eq. '(' ) then

               ifn = 9
               imc = 1
               ic  = ic + 3


         else if(
     &       ( dum(ic  ) .eq. 'N' .or. dum(ic  ) .eq. 'n' ) .and.
     &       ( dum(ic+1) .eq. 'I' .or. dum(ic+1) .eq. 'i' ) .and.
     &       ( dum(ic+2) .eq. 'N' .or. dum(ic+2) .eq. 'n' ) .and.
     &       ( dum(ic+3) .eq. 'T' .or. dum(ic+3) .eq. 't' ) .and.
     &         dum(ic+4) .eq. '(' ) then

               ifn = 10
               ic  = ic + 4


         else if(
     &       ( dum(ic  ) .eq. 'S' .or. dum(ic  ) .eq. 's' ) .and.
     &       ( dum(ic+1) .eq. 'I' .or. dum(ic+1) .eq. 'i' ) .and.
     &       ( dum(ic+2) .eq. 'G' .or. dum(ic+2) .eq. 'g' ) .and.
     &       ( dum(ic+3) .eq. 'N' .or. dum(ic+3) .eq. 'n' ) .and.
     &         dum(ic+4) .eq. '(' ) then

               ifn = 11
               imc = 1
               ic  = ic + 4


         else if(
     &       ( dum(ic  ) .eq. 'S' .or. dum(ic  ) .eq. 's' ) .and.
     &       ( dum(ic+1) .eq. 'Q' .or. dum(ic+1) .eq. 'q' ) .and.
     &       ( dum(ic+2) .eq. 'R' .or. dum(ic+2) .eq. 'r' ) .and.
     &       ( dum(ic+3) .eq. 'T' .or. dum(ic+3) .eq. 't' ) .and.
     &         dum(ic+4) .eq. '(' ) then

               ifn = 12
               ic  = ic + 4


         else if(
     &       ( dum(ic  ) .eq. 'A' .or. dum(ic  ) .eq. 'a' ) .and.
     &       ( dum(ic+1) .eq. 'C' .or. dum(ic+1) .eq. 'c' ) .and.
     &       ( dum(ic+2) .eq. 'O' .or. dum(ic+2) .eq. 'o' ) .and.
     &       ( dum(ic+3) .eq. 'S' .or. dum(ic+3) .eq. 's' ) .and.
     &         dum(ic+4) .eq. '(' ) then

               ifn = 13
               ic  = ic + 4


         else if(
     &       ( dum(ic  ) .eq. 'A' .or. dum(ic  ) .eq. 'a' ) .and.
     &       ( dum(ic+1) .eq. 'S' .or. dum(ic+1) .eq. 's' ) .and.
     &       ( dum(ic+2) .eq. 'I' .or. dum(ic+2) .eq. 'i' ) .and.
     &       ( dum(ic+3) .eq. 'N' .or. dum(ic+3) .eq. 'n' ) .and.
     &         dum(ic+4) .eq. '(' ) then

               ifn = 14
               ic  = ic + 4


         else if(
     &       ( dum(ic  ) .eq. 'A' .or. dum(ic  ) .eq. 'a' ) .and.
     &       ( dum(ic+1) .eq. 'T' .or. dum(ic+1) .eq. 't' ) .and.
     &       ( dum(ic+2) .eq. 'A' .or. dum(ic+2) .eq. 'a' ) .and.
     &       ( dum(ic+3) .eq. 'N' .or. dum(ic+3) .eq. 'n' ) .and.
     &         dum(ic+4) .eq. '(' ) then

               ifn = 15
               ic  = ic + 4










         else if(
     &       ( dum(ic  ) .eq. 'A' .or. dum(ic  ) .eq. 'a' ) .and.
     &       ( dum(ic+1) .eq. 'T' .or. dum(ic+1) .eq. 't' ) .and.
     &       ( dum(ic+2) .eq. 'A' .or. dum(ic+2) .eq. 'a' ) .and.
     &       ( dum(ic+3) .eq. 'N' .or. dum(ic+3) .eq. 'n' ) .and.
     &       ( dum(ic+4) .eq. '2' ) .and.
     &         dum(ic+5) .eq. '(' ) then

               ifn = 19
               imc = 1
               ic  = ic + 5


         else if(
     &       ( dum(ic  ) .eq. 'C' .or. dum(ic  ) .eq. 'c' ) .and.
     &       ( dum(ic+1) .eq. 'O' .or. dum(ic+1) .eq. 'o' ) .and.
     &       ( dum(ic+2) .eq. 'S' .or. dum(ic+2) .eq. 's' ) .and.
     &         dum(ic+3) .eq. '(' ) then

               ifn = 20
               ic  = ic + 3


         else if(
     &       ( dum(ic  ) .eq. 'C' .or. dum(ic  ) .eq. 'c' ) .and.
     &       ( dum(ic+1) .eq. 'O' .or. dum(ic+1) .eq. 'o' ) .and.
     &       ( dum(ic+2) .eq. 'S' .or. dum(ic+2) .eq. 's' ) .and.
     &       ( dum(ic+3) .eq. 'H' .or. dum(ic+3) .eq. 'h' ) .and.
     &         dum(ic+4) .eq. '(' ) then

               ifn = 21
               ic  = ic + 4


         else if(
     &       ( dum(ic  ) .eq. 'S' .or. dum(ic  ) .eq. 's' ) .and.
     &       ( dum(ic+1) .eq. 'I' .or. dum(ic+1) .eq. 'i' ) .and.
     &       ( dum(ic+2) .eq. 'N' .or. dum(ic+2) .eq. 'n' ) .and.
     &         dum(ic+3) .eq. '(' ) then

               ifn = 22
               ic  = ic + 3


         else if(
     &       ( dum(ic  ) .eq. 'S' .or. dum(ic  ) .eq. 's' ) .and.
     &       ( dum(ic+1) .eq. 'I' .or. dum(ic+1) .eq. 'i' ) .and.
     &       ( dum(ic+2) .eq. 'N' .or. dum(ic+2) .eq. 'n' ) .and.
     &       ( dum(ic+3) .eq. 'H' .or. dum(ic+3) .eq. 'h' ) .and.
     &         dum(ic+4) .eq. '(' ) then

               ifn = 23
               ic  = ic + 4


         else if(
     &       ( dum(ic  ) .eq. 'T' .or. dum(ic  ) .eq. 't' ) .and.
     &       ( dum(ic+1) .eq. 'A' .or. dum(ic+1) .eq. 'a' ) .and.
     &       ( dum(ic+2) .eq. 'N' .or. dum(ic+2) .eq. 'n' ) .and.
     &         dum(ic+3) .eq. '(' ) then

               ifn = 24
               ic  = ic + 3


         else if(
     &       ( dum(ic  ) .eq. 'T' .or. dum(ic  ) .eq. 't' ) .and.
     &       ( dum(ic+1) .eq. 'A' .or. dum(ic+1) .eq. 'a' ) .and.
     &       ( dum(ic+2) .eq. 'N' .or. dum(ic+2) .eq. 'n' ) .and.
     &       ( dum(ic+3) .eq. 'H' .or. dum(ic+3) .eq. 'h' ) .and.
     &         dum(ic+4) .eq. '(' ) then

               ifn = 25
               ic  = ic + 4

         else

            func02 = .false.

         end if

*-----------------------------------------------------------------------

      return
      end



************************************************************************
*                                                                      *
      subroutine func03(ierr,
     &                  kd,kp,ki,jd,mg,md,me,mp,id,
     &                  ivi,ivn,ifp,ipr)
*                                                                      *
*      PURPOSE  : FROM FUNCTION DESCRIPTION                            *
*                 DETERMINE THE PROCEDURE                              *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      dimension kd(0:mfc), kp(0:mfc), ki(0:mfc), jd(0:mfc)
      dimension md(0:mfc), mg(0:mfc), me(0:mfc), mp(0:mfc)
      dimension ivi(0:mfc,0:mfc), ivn(0:mfc), ifp(0:mfc)

      dimension jf(0:mfc)
      dimension idv(0:mfc)

*-----------------------------------------------------------------------
*        initialization
*-----------------------------------------------------------------------

               ierr = 0

               do 10 i = 0, mfc

                  jf(i) = 0
                  idv(i) = 0

   10          continue

               ipr = 0
               idd = id

               do 20 i = 1, idd

                  idv(i) = i

   20          continue

               icons = 1

*-----------------------------------------------------------------------

   90 continue

*-----------------------------------------------------------------------
*           last idd
*-----------------------------------------------------------------------

            if( idd .eq. 1 ) then

               if( icons .eq. 1 ) then

                     i = idv(1)

                     if( jd(i) .ne. 0 ) goto 999

                     ipr = ipr + 1

                     ivn(ipr) = 1
                     ifp(ipr) = -1

                  if( jf(i) .eq. 0 ) then

                     ivi(1,ipr) = -i

                  else

                     ivi(1,ipr) = jf(i)

                  end if

               end if

                  goto 1000

            end if

*-----------------------------------------------------------------------
*        find one highest region
*-----------------------------------------------------------------------

               ikp = 0
               iid = 1

            do 100 j = 1, idd

                  i = idv(j)

               if( kp(i) .gt. ikp ) then

                  ikp = kp(i)

                  iid = j

                  goto 100

               end if

               if( kp(i) .lt. ikp ) then

                  ifd = j - 1

                  goto 200

               end if

  100       continue

                  ifd = idd

  200       continue

*-----------------------------------------------------------------------
*        FIND ONE PROCEDURE IN THE ABOVE REGION
*-----------------------------------------------------------------------

               ii1 = iid

  290    continue

*-----------------------------------------------------------------------

            if( iid .eq. ifd ) then

                  kp(idv(iid)) = max(0,kp(idv(iid)) - 1)

                  goto 90

            end if

*-----------------------------------------------------------------------
*           SEARCH FIRST VARIABLE
*-----------------------------------------------------------------------

                  icons = 0

            do 300 j = iid, ifd

                  i = idv(j)

*-----------------------------------------------------------------------
*              ----( **, *, / )----
*-----------------------------------------------------------------------

               if( jd(i) .eq. 1 .or. jd(i) .eq. 2 )  goto 999

*-----------------------------------------------------------------------
*              ----( + )----
*-----------------------------------------------------------------------

               if( kd(i) .eq. 7 .and. ( j+1 .gt. ifd .or.
     &           ( j+1 .le. ifd .and. jd(idv(j+1)) .ne. 0 ) ) ) goto 999

               if( kd(i) .eq. 7 ) then

                  do 351 k = j, idd-1

                     idv(k) = idv(k+1)

  351             continue

                     ifd = ifd - 1
                     idd = idd - 1

                     icons = 1

                  goto 290

               end if

*-----------------------------------------------------------------------
*              ----( - )----
*-----------------------------------------------------------------------

               if( kd(i) .eq. 8 .and. ( j+1 .gt. ifd .or.
     &           ( j+1 .le. ifd .and. jd(idv(j+1)) .ne. 0 ) ) ) goto 999

               if( kd(i) .eq. 8 .and. ( j+1 .eq. ifd .or.
     &           ( j+2 .le. ifd .and. jd(idv(j+2)) .ne. 1 ) ) ) then

                  ipr = ipr + 1

                  ivi(1,ipr) = 0
                  ifp(ipr)   = 0
                  ivn(ipr)   = 2

                  if( jf(idv(j+1)) .eq. 0 ) then

                     ivi(2,ipr) = - idv(j+1)

                  else

                     ivi(2,ipr) = jf(idv(j+1))

                  end if

                     jf(idv(j+1)) = ipr

                  do 352 k = j, idd-1

                     idv(k) = idv(k+1)

  352             continue

                     ifd = ifd - 1
                     idd = idd - 1

                  goto 290

               end if

*-----------------------------------------------------------------------
*              ----( ( ONE VARIABLE ) )----
*-----------------------------------------------------------------------

               if( jd(i) .eq. 0 .and.
     &             ifd .eq. iid + 2 .and.
     &             jd(idv(iid)) .eq. -1 .and.
     &             jd(idv(ifd)) .eq. -1 ) then

                     kp(i) = max(0,kp(i)-1)

                     idv(iid) = idv(iid+1)

                  do 355 k = iid+1, idd-2

                     idv(k) = idv(k+2)

  355             continue

                     idd = idd - 2

                     icons = 1

                     goto 90

               end if

*-----------------------------------------------------------------------
*              ----( FUNCTION( ONE VARIABLE ) )----
*-----------------------------------------------------------------------

               if( jd(i) .eq. 0 .and.
     &             ifd .eq. iid + 2 .and.
     &             jd(idv(iid)) .eq. -2 .and.
     &             jd(idv(ifd)) .eq. -1 ) then

                     ipr = ipr + 1

                  if( jf(idv(j)) .eq. 0 ) then

                     ivi(1,ipr) = - idv(j)

                  else

                     ivi(1,ipr) = jf(idv(j))

                  end if

                     ifp(ipr) = idv(iid)
                     ivn(ipr) = 1

                     jf(i) = ipr


                     kp(i) = max(0,kp(i)-1)

                     idv(iid) = idv(iid+1)

                  do 356 k = iid+1, idd-2

                     idv(k) = idv(k+2)

  356             continue

                     idd = idd - 2

                     goto 90

               end if

*-----------------------------------------------------------------------
*              ----( FUNCTION( A,B,C SOME VARIABLES ) )----
*-----------------------------------------------------------------------

*                    MAX(A,B,C)
*                    MIN(A,B,C)
*                    ATAN2(A,B)

*-----------------------------------------------------------------------

               if( jd(i) .eq. 0 .and.
     &             jd(idv(iid)) .eq. -2 .and.
     &           ( ki(idv(iid)) .eq.  7 .or.
     &             ki(idv(iid)) .eq.  8 .or.
     &             ki(idv(iid)) .eq.  9 .or.
     &             ki(idv(iid)) .eq. 11 .or.
     &             ki(idv(iid)) .eq. 19 ) ) then


                     ipr = ipr + 1

                     ifp(ipr) = idv(iid)
                     ivn(ipr) = me(idv(iid))

                     if( ivn(ipr) .lt. 2 ) goto 999

                  do 358 k = 1, ivn(ipr)

                        kk = iid + ( k - 1) * 2 + 1

                        if( k .gt. 1 .and.
     &                      jd(idv(kk-1)) .ne. -3 ) goto 999

                        if( k .eq. ivn(ipr) .and.
     &                      jd(idv(kk+1)) .ne. -1 ) goto 999

                     if( jf(idv(kk)) .eq. 0 ) then

                        ivi(k,ipr) = - idv(kk)

                     else

                        ivi(k,ipr) = jf(idv(kk))

                     end if

  358             continue

                        jf(i) = ipr

                        kp(i) = max(0,kp(i)-1)

                        idv(iid) = idv(iid+1)

                     do 357 k = iid+1, idd-2*ivn(ipr)

                        idv(k) = idv(k+2*ivn(ipr))

  357                continue

                        idd = idd - 2*ivn(ipr)

                        goto 90

               end if

*-----------------------------------------------------------------------
*              ----( FIRST VARIABLE )----
*-----------------------------------------------------------------------

               if( jd(i) .eq. 0 ) then

                  ii1 = j

                  goto 301

               end if

*-----------------------------------------------------------------------

  300       continue

*-----------------------------------------------------------------------
*              ----( () )----
*-----------------------------------------------------------------------

               if( ifd .eq. iid + 1 .and.
     &             jd(idv(iid)) .eq. -1 .and.
     &             jd(idv(ifd)) .eq. -1 ) then

                  do 353 k = iid, idd-2

                     idv(k) = idv(k+2)

  353             continue

                     idd = idd - 2

                     goto 90

               end if

*-----------------------------------------------------------------------
*              ---- ERROR ----
*-----------------------------------------------------------------------

                  goto 999

  301    continue

*-----------------------------------------------------------------------
*           FOR FUNC(A,B,C)
*-----------------------------------------------------------------------

            if( ii1 .eq. ifd ) then

               if( mp(idv(iid)) .ge. 1 ) then

                  kp(idv(iid)) = max(0,kp(idv(iid)) - 1)

                  goto 90

               else

                  goto 999

               end if

            end if

*-----------------------------------------------------------------------
*           SEARCH SECOND VARIABLE
*-----------------------------------------------------------------------

               ii2 = ii1
               if1 = 0
               if2 = 0
               in2 = 0
               ifj = 0

            do 302 j = ii1 + 1, ifd
! T.Sato 2015/3/5 to avoid bug due to intel fortran 15.0
               if(j.eq.ifd+1) exit

               i = idv(j)

               if( if1 .eq. 0 .and. jd(i) .eq. 0 ) goto 999

               if( jd(i) .eq. 0 ) then

                  ii2 = j
                  in2 = 1

               end if

               if( in2 .eq. 0 .and.
     &           ( jd(i) .eq. 1 .or.
     &             jd(i) .eq. 2 .or.
     &             jd(i) .eq. 3 ) ) then

                  if1 = jd(i)
                  ifj = j

               end if

               if( in2 .eq. 1 .and.
     &           ( jd(i) .eq. 1 .or.
     &             jd(i) .eq. 2 .or.
     &             jd(i) .eq. 3 ) ) then

                  if2 = jd(i)

                  if( if2 .lt. if1 ) then

                     ii1 = ii2

                     goto 301

                  else

                     goto 303

                  end if

               end if

  302       continue

               if( in2 .eq. 0 ) goto 999

*-----------------------------------------------------------------------
*           ---- DEFINE ONE PROCEDURE ----
*-----------------------------------------------------------------------

  303       continue

                  ipr = ipr + 1

               if( jf(idv(ii1)) .eq. 0 ) then

                  ivi(1,ipr) = - idv(ii1)

               else

                  ivi(1,ipr) = jf(idv(ii1))

               end if

               if( jf(idv(ii2)) .eq. 0 ) then

                  ivi(2,ipr) = - idv(ii2)

               else

                  ivi(2,ipr) = jf(idv(ii2))

               end if

                  ifp(ipr) = idv(ifj)
                  ivn(ipr) = 2


                  jf(idv(ii1)) = ipr
                  jd(idv(ii1)) = 0

                  do 354 k = ii1+1, idd-2

                     idv(k) = idv(k+2)

  354             continue

                     ifd = ifd - 2
                     idd = idd - 2

                  goto 290


*-----------------------------------------------------------------------
*        FINAL CHECK
*-----------------------------------------------------------------------

 1000    continue

            if( jd(idv(idd)) .ne. 0 ) goto 999

      return

*-----------------------------------------------------------------------

  999 continue
      ierr = 1

      return
      end


************************************************************************
*                                                                      *
      subroutine func04(dum,ic0,iclm,ierr,
     &                  xvin,xvdd,nnxv)
*                                                                      *
*      PURPOSE  : READ V-VALUE DESCRIPTION                             *
*                                                                      *
*      FORMAT   : V=[ XVIN, XVFN, NNXV ]                               *
*                                                                      *
*                     XVIN ; INITIAL VALUE OF V                        *
*                     XVFN ; FINAL VALUE OF V                          *
*                     NNXV ; NUMBER OF VPOINT - 1                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1

      common /rval1/ cval(mxcval), aval(mxcval)

      dimension xvv(3)

*-----------------------------------------------------------------------
*     INITIALIZATION
*-----------------------------------------------------------------------

            ierr = 0

            ic  = ic0
            ici = ic0

            icn = 0

*-----------------------------------------------------------------------
*        READ TEXT AND EVALUATE FUNCTION
*-----------------------------------------------------------------------

  100       ic  = ic + 1

               if( ic .gt. iclm ) goto 999

               if( dum(ic) .ne. ',' .and. dum(ic) .ne. ']' ) goto 100

*-----------------------------------------------------------------------

               if( dum(ic) .eq. ',' ) then

                     icn = icn + 1

                     if( icn .gt. 2 ) goto 999

                  call func12(dum,ici,ierr,iclm,
     &                        ',',cvvv,xval)

                     if( ierr .ne. 0 ) goto 999

                     xvv(icn) = cvvv

                     ici = ici

                     goto 100

               else if( dum(ic) .eq. ']' ) then

                     icn = icn + 1

                     if( icn .ne. 3 ) goto 999


                  call func12(dum,ici,ierr,iclm,
     &                        ']',cvvv,xval)

                     if( ierr .ne. 0 ) goto 999

                     xvv(icn) = cvvv

                     goto 300

               end if

*-----------------------------------------------------------------------

            goto 999

*-----------------------------------------------------------------------
*        FINAL CHECK
*-----------------------------------------------------------------------

  300    ic0 = ici

            xvin = xvv(1)
            xvfn = xvv(2)
            nnxv = nint( xvv(3) ) + 1

            if( nnxv .lt. 2 ) goto 999

            xvdd = ( xvfn - xvin ) / dble( nnxv - 1 )

         return

*-----------------------------------------------------------------------
*        ERROR
*-----------------------------------------------------------------------

  999 continue

      ierr = 1

      return
      end


************************************************************************
*                                                                      *
      subroutine func05(nfunc,nnfun,neval,eval,yvalp,ierr)
*                                                                      *
*      PURPOSE  : EVALUATE ONE PROCEDURE                               *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      dimension eval(0:mfc)

*-----------------------------------------------------------------------

         ierr = 0

         if( nfunc .le. 4 .or. nfunc .ge. 11 ) then

               yvalp = eval(1)

         else if( nfunc .eq. 5 ) then

               yvalp = eval(1) * eval(2)

         else if( nfunc .eq. 6 ) then

               yvalp = eval(1) / eval(2)

         else if( nfunc .eq. 7 ) then

               yvalp = eval(1) + eval(2)

         else if( nfunc .eq. 8 ) then

               yvalp = eval(1) - eval(2)

         else if( nfunc .eq. 9 ) then

            if( eval(1) .ge. 0.0 ) then

               yvalp = eval(1) ** eval(2)

            else

               if( abs(eval(2) - real(nint(eval(2)))) .lt. 1.e-4 ) then

                  yvalp = eval(1) ** nint(eval(2))

               else

                  goto 999

               end if

            end if


         else if( nfunc .eq. 10 ) then

            if( nnfun .eq. 1 ) then

               yvalp = dble( nint( eval(1) ) )

            else if( nnfun .eq.  2 ) then

               yvalp = int( eval(1) )

            else if( nnfun .eq.  3 ) then

               yvalp = abs( eval(1) )

            else if( nnfun .eq.  4 ) then

               yvalp = exp( eval(1) )

            else if( nnfun .eq.  5 ) then

               if( eval(1) .le. 0.0 ) goto 999

               yvalp = log( eval(1) )

            else if( nnfun .eq.  6 ) then

               if( eval(1) .le. 0.0 ) goto 999

               yvalp = log10( eval(1) )

            else if( nnfun .eq.  7 ) then

               yvalp = max( eval(1), eval(2) )

               if( neval .ge. 3 ) then

                  do 10 i = 3, neval

                     yvalp = max( yvalp, eval(i) )

   10             continue

               end if

            else if( nnfun .eq.  8 ) then

               yvalp = min( eval(1), eval(2) )

               if( neval .ge. 3 ) then

                  do 20 i = 3, neval

                     yvalp = min( yvalp, eval(i) )

   20             continue

               end if

            else if( nnfun .eq.  9 ) then

               yvalp = mod( eval(1), eval(2) )

            else if( nnfun .eq. 10 ) then

               yvalp = nint( eval(1) )

            else if( nnfun .eq. 11 ) then

               yvalp = sign( eval(1), eval(2) )

            else if( nnfun .eq. 12 ) then

               if( eval(1) .le. 0 ) goto 999

               yvalp = sqrt( eval(1) )

            else if( nnfun .eq. 13 ) then

               yvalp = acos( eval(1) )

            else if( nnfun .eq. 14 ) then

               yvalp = asin( eval(1) )

            else if( nnfun .eq. 15 ) then

               yvalp = atan( eval(1) )

*-----------------------------------------------------------------------







*-----------------------------------------------------------------------

            else if( nnfun .eq. 19 ) then

               yvalp = atan2( eval(1), eval(2) )

            else if( nnfun .eq. 20 ) then

               yvalp = cos( eval(1) )

            else if( nnfun .eq. 21 ) then

               yvalp = cosh( eval(1) )

            else if( nnfun .eq. 22 ) then

               yvalp = sin( eval(1) )

            else if( nnfun .eq. 23 ) then

               yvalp = sinh( eval(1) )

            else if( nnfun .eq. 24 ) then

               yvalp = tan( eval(1) )

            else if( nnfun .eq. 25 ) then

               yvalp = tanh( eval(1) )

            end if

         end if

      return

*-----------------------------------------------------------------------

  999 ierr = 1

      return
      end

************************************************************************
*                                                                      *
      subroutine func07(idn,kd,ki,rk,ivi,ivn,ifp,ipr,id,
     &                  kdjj,kijj,rkjj,ivijj,ivnjj,ifpjj,
     &                  iprjj,idrjj)
*                                                                      *
*      PURPOSE  : STORE THE ONE VALUABLE PROCEDURE                     *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      dimension rk(0:mfc)
      dimension kd(0:mfc), kp(0:mfc), ki(0:mfc), jd(0:mfc)
      dimension ivi(0:mfc,0:mfc), ivn(0:mfc), ifp(0:mfc)

      dimension rkjj(mh,0:mfc)
      dimension iprjj(mh), idrjj(mh)
      dimension kdjj(mh,0:mfc), kijj(mh,0:mfc)
      dimension ivnjj(mh,0:mfc), ifpjj(mh,0:mfc)
      dimension ivijj(mh,0:mfc,0:mfc)

*-----------------------------------------------------------------------

            iprjj(idn) = ipr

         do 100 i = 1, ipr

            ifpjj(idn,i) = ifp(i)
            ivnjj(idn,i) = ivn(i)

         do 200 j = 1, ivn(i)

            ivijj(idn,j,i) = ivi(j,i)

  200    continue
  100    continue


            idrjj(idn) = id

         do 300 i = 1, id

            kdjj(idn,i) = kd(i)
            kijj(idn,i) = ki(i)
            rkjj(idn,i) = rk(i)

  300    continue

*-----------------------------------------------------------------------

         return
         end


************************************************************************
*                                                                      *
      subroutine func08(jj,ila,daxy,
     &                  ifdc,ivdc,idyd,ided,idef,idvd,idxd,
     &                  ipr,ifp,kd,ki,rk,ivn,ivi,
     &                  ierr,jsn,illi,
     &                  nevk,nevn)
*                                                                      *
*      PURPOSE  : EVALUATE ONE COLUMN VALUE                            *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      dimension daxy(mh+2,*)
      dimension dval(mh+1)

      dimension idyd(-1:mh),ided(-1:mh),idef(-1:mh)

      common /rval1/ cval(mxcval), aval(mxcval)

      dimension rk(0:mfc)
      dimension kd(0:mfc), kp(0:mfc), ki(0:mfc), jd(0:mfc)
      dimension ivi(0:mfc,0:mfc), ivn(0:mfc), ifp(0:mfc)

      dimension ifdc(mh), ivdc(mh)

      dimension nfun(0:mfc), nnfn(0:mfc), nevl(0:mfc)
      dimension nevk(0:mfc,0:mfc), nevn(0:mfc,0:mfc)

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------
*        CHECK FOR PROCEDURE
*-----------------------------------------------------------------------

            do 311 i = 1, ipr

*-----------------------------------------------------------------------
*                 DEFINITION OF FUNCTION
*-----------------------------------------------------------------------

                  if( ifp(i) .eq. 0 ) then

                        nfunc = 5
                        nnfun = 0

                  else if( ifp(i) .eq. -1 ) then

                        nfunc = 1
                        nnfun = 0

                  else if( ifp(i) .gt. 0 ) then

                        nfunc = kd(ifp(i))

                     if( nfunc .ne. 10 ) then

                        nnfun = 0

                     else

                        nnfun = ki(ifp(i))

                     end if

                  else

                        goto 998

                  end if

                        nfun(i) = nfunc
                        nnfn(i) = nnfun

*-----------------------------------------------------------------------
*                 DEFINITION OF VARIABLES
*-----------------------------------------------------------------------

                        neval = ivn(i)
                        nevl(i) = neval

                  do 321 j = 1, neval

                           iiv = ivi(j,i)

                     if( iiv .lt. 0 ) then

                           iiv = - iiv

*-----------------------------------------------------------------------
*                          --- X-VALUE ---
*-----------------------------------------------------------------------

                        if( kd(iiv) .eq. 1 ) then

                           if( ivdc(idxd) .eq. 0 ) goto 998

                           nevk(i,j) = 1
                           nevn(i,j) = idxd

*-----------------------------------------------------------------------
*                          --- Y-VALUE ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 2 ) then

                           kdn = idyd(ki(iiv))

                           if( kdn .eq. -1 ) goto 998

                           if( ivdc(kdn) .eq. 0 ) goto 998

                           nevk(i,j) = 1
                           nevn(i,j) = kdn

*-----------------------------------------------------------------------
*                          --- CONSTANTS ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 3 ) then

                           if( cval(ki(iiv)) .lt. -r0max )
     &                     goto 998

                           nevk(i,j) = 2
                           nevn(i,j) = ki(iiv)

*-----------------------------------------------------------------------
*                          --- REAL NUMBERS ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 4 ) then

                           nevk(i,j) = 3
                           nevn(i,j) = iiv

*-----------------------------------------------------------------------
*                          --- V-VALUES ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 15 ) then

                           if( idvd .ne. 1 ) goto 998

                           if( ivdc(1) .eq. 0 ) stop 106 !goto 998

                           nevk(i,j) = 1
                           nevn(i,j) = 1

*-----------------------------------------------------------------------
*                          --- DY-VALUES ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 16 ) then

                           kid  = ki(iiv)

                           if( kid .ge. 0 ) then

                              kdn = ided(kid)

                              if( kdn .eq. -1 ) goto 998

                           else

                              kid  = -kid

                              kdn = idef(kid)

                              if( kdn .eq. -1 ) goto 998

                           end if

                           if( ivdc(kdn) .eq. 0 ) goto 998

                           nevk(i,j) = 1
                           nevn(i,j) = kdn

*-----------------------------------------------------------------------
*                          --- DX-VALUES ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 17 ) then

                           kid  = ki(iiv)

                           if( kid .ge. 0 ) then

                              kdn = ided(kid)

                              if( kdn .eq. -1 ) goto 998

                           else

                              kid  = -kid

                              kdn = idef(kid)

                              if( kdn .eq. -1 ) goto 998

                           end if

                           if( ivdc(kdn) .eq. 0 ) goto 998

                           nevk(i,j) = 1
                           nevn(i,j) = kdn

*-----------------------------------------------------------------------
*                          --- PARAMETERS ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 14 ) then

                           if( aval(ki(iiv)) .lt. -r0max )
     &                     goto 998

                           nevk(i,j) = 6
                           nevn(i,j) = ki(iiv)

*-----------------------------------------------------------------------

                        else

                           goto 998

                        end if


                     else if( iiv .gt. 0 ) then

                           if( iiv .ge. i ) goto 998

                           nevk(i,j) = 4
                           nevn(i,j) = iiv

                     else if( iiv .eq. 0 ) then

                           nevk(i,j) = 5
                           nevn(i,j) = 0

                     end if

  321             continue

  311       continue


*-----------------------------------------------------------------------
*        DO LOOP FOR X
*-----------------------------------------------------------------------

         do 300 ix = 1, ila

                     do 310 j = 1, mh + 1

                        dval(j) = daxy(j,ix)

  310                continue


                     call func10(rk,dval,
     &                           ipr,ierr,dyval,
     &                           nfun,nnfn,nevl,nevk,nevn)

                        if( ierr .ne. 0 ) goto 998

                     daxy(jj,ix) = dyval

  300    continue


            ierr = 0

            return

*-----------------------------------------------------------------------

  998    continue

            m_err = 'H: Detail of Function Description is Wrong.'
            ErrCha = ''
            ErrID = 'L:2408/R:func08/F:a-func.f'
            l_err = illi
            k_err = jsn

            ierr = 1

            return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine func09(jj,idn,ila,daxy,
     &                  ifdc,ivdc,idyd,ided,idef,idvd,idxd,
     &                  ipr,ifp,kd,ki,rk,ivn,ivi,
     &                  ierr,jsn,illi,iffy,iffd,iyyd,iydd,
     &                  nevk,nevn)
*                                                                      *
*      PURPOSE  : FITTING THE FUNCTION                                 *
*                                                                      *
*        IFFY(JJ) = FITTING Y ID                                       *
*        IFFD(JJ) = 1 -> WITHOUT ERROR                                 *
*        IFFD(JJ) = 2 -> WITH ERROR                                    *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      dimension daxy(mh+2,1)
      dimension dval(mh+1)

      common /rval1/ cval(mxcval), aval(mxcval)

      dimension ifdc(mh), ivdc(mh)
      dimension iffy(-1:mh), iffd(-1:mh)
      dimension idyd(-1:mh),ided(-1:mh),idef(-1:mh)
      dimension iyyd(mh,2),iydd(mh)

      dimension rk(0:mfc)
      dimension kd(0:mfc), kp(0:mfc), ki(0:mfc), jd(0:mfc)
      dimension ivi(0:mfc,0:mfc), ivn(0:mfc), ifp(0:mfc)

      dimension nfun(0:mfc), nnfn(0:mfc), nevl(0:mfc)
      dimension nevk(0:mfc,0:mfc), nevn(0:mfc,0:mfc)

      dimension ianm(0:mfc), ikta(0:mfc), iatk(0:mfc), iata(0:mfc)
      dimension nhit(0:mfc), rdev(0:mfc)

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------
*        INITIAL VALUE
*-----------------------------------------------------------------------

               ian = 0

*-----------------------------------------------------------------------
*        CHECK OF FITTING Y-VALUE
*-----------------------------------------------------------------------

               kdy = idyd(iffy(jj))

               if( kdy .eq. -1 ) goto 998

               if( ivdc(kdy) .eq. 0 ) goto 998

*-----------------------------------------------------------------------
*        CHECK OF FITTING Y-VALUE WITH ERROR BAR
*-----------------------------------------------------------------------

            if( iffd(jj) .eq. 2 ) then

               if( iydd(kdy) .eq. 0 ) goto 998

                  do 90 ix = 1, ila

                     daxy(mh+2,ix) = 0.0

   90             continue

               if( iyyd(kdy,1) .ne. idn+1 ) then

                  if( ivdc(iyyd(kdy,1)) .eq. 0 ) goto 998

                  do 100 ix = 1, ila

                     if( daxy(iyyd(kdy,1),ix) .le. 0.0 ) goto 997

                     daxy(mh+2,ix) = daxy(mh+2,ix)
     &                             + daxy(iyyd(kdy,1),ix)

  100             continue

               end if

               if( iyyd(kdy,2) .ne. idn+1 ) then

                  if( ivdc(iyyd(kdy,2)) .eq. 0 ) goto 998

                  do 101 ix = 1, ila

                     if( daxy(iyyd(kdy,2),ix) .le. 0.0 ) goto 997

                     daxy(mh+2,ix) = daxy(mh+2,ix)
     &                             + daxy(iyyd(kdy,1),ix)

  101             continue

               end if

            end if

*-----------------------------------------------------------------------
*        CHECK FOR PROCEDURE
*-----------------------------------------------------------------------

            do 311 i = 1, ipr

*-----------------------------------------------------------------------
*                 DEFINITION OF FUNCTION
*-----------------------------------------------------------------------

                  if( ifp(i) .eq. 0 ) then

                        nfunc = 5
                        nnfun = 0

                  else if( ifp(i) .eq. -1 ) then

                        nfunc = 1
                        nnfun = 0

                  else if( ifp(i) .gt. 0 ) then

                        nfunc = kd(ifp(i))

                     if( nfunc .ne. 10 ) then

                        nnfun = 0

                     else

                        nnfun = ki(ifp(i))

                     end if

                  else

                        goto 998

                  end if

                        nfun(i) = nfunc
                        nnfn(i) = nnfun

*-----------------------------------------------------------------------
*                 DEFINITION OF VARIABLES
*-----------------------------------------------------------------------

                        neval = ivn(i)
                        nevl(i) = neval

                  do 321 j = 1, neval

                           iiv = ivi(j,i)

                     if( iiv .lt. 0 ) then

                           iiv = - iiv

*-----------------------------------------------------------------------
*                          --- X-VALUE ---
*-----------------------------------------------------------------------

                        if( kd(iiv) .eq. 1 ) then

                           if( ivdc(idxd) .eq. 0 ) goto 998

                           nevk(i,j) = 1
                           nevn(i,j) = idxd

*-----------------------------------------------------------------------
*                          --- Y-VALUE ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 2 ) then

                           kdn = idyd(ki(iiv))

                           if( kdn .eq. -1 ) goto 998

                           if( ivdc(kdn) .eq. 0 ) goto 998

                           nevk(i,j) = 1
                           nevn(i,j) = kdn

*-----------------------------------------------------------------------
*                          --- CONSTANTS ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 3 ) then

                           if( cval(ki(iiv)) .lt. -r0max )
     &                     goto 998

                           nevk(i,j) = 2
                           nevn(i,j) = ki(iiv)

*-----------------------------------------------------------------------
*                          --- REAL NUMBERS ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 4 ) then

                           nevk(i,j) = 3
                           nevn(i,j) = iiv

*-----------------------------------------------------------------------
*                          --- V-VALUES ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 15 ) then

                           if( idvd .ne. 1 ) goto 998

                           if( ivdc(1) .eq. 0 ) goto 998

                           nevk(i,j) = 1
                           nevn(i,j) = 1

*-----------------------------------------------------------------------
*                          --- DY-VALUES ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 16 ) then

                           kid  = ki(iiv)

                           if( kid .ge. 0 ) then

                              kdn = ided(kid)

                              if( kdn .eq. -1 ) goto 998

                           else

                              kid  = -kid

                              kdn = idef(kid)

                              if( kdn .eq. -1 ) goto 998

                           end if

                           if( ivdc(kdn) .eq. 0 ) goto 998

                           nevk(i,j) = 1
                           nevn(i,j) = kdn

*-----------------------------------------------------------------------
*                          --- DX-VALUES ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 17 ) then

                           kid  = ki(iiv)

                           if( kid .ge. 0 ) then

                              kdn = ided(kid)

                              if( kdn .eq. -1 ) goto 998

                           else

                              kid  = -kid

                              kdn = idef(kid)

                              if( kdn .eq. -1 ) goto 998

                           end if

                           if( ivdc(kdn) .eq. 0 ) goto 998

                           nevk(i,j) = 1
                           nevn(i,j) = kdn

*-----------------------------------------------------------------------
*                          --- PARAMETERS ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 14 ) then

                           nevk(i,j) = 6
                           nevn(i,j) = ki(iiv)

                           ian = ian + 1
                           ianm(ian) = ki(iiv)

*-----------------------------------------------------------------------

                        else

                           goto 998

                        end if


                     else if( iiv .gt. 0 ) then

                           if( iiv .ge. i ) goto 998

                           nevk(i,j) = 4
                           nevn(i,j) = iiv

                     else if( iiv .eq. 0 ) then

                           nevk(i,j) = 5
                           nevn(i,j) = 0

                     end if

  321             continue

  311       continue

*=======================================================================
*        SUMMARY OF PARAMETERS
*-----------------------------------------------------------------------

            if( ian .eq. 0 ) goto 998

            do 200 i = 1, mfc

               iatk(i) = 0
               ikta(i) = 0

  200       continue

               iaa = 0

            do 201 i = 1, ian

               if( ikta(ianm(i)) .eq. 0 ) then

                  iaa = iaa + 1

                  iatk(iaa) = ianm(i)

                  ikta(ianm(i)) = iaa

               end if

  201       continue

*=======================================================================
*        INITIAL VALUE OF PARAMETERS
*-----------------------------------------------------------------------

                  eps1 = 1.0e-7
                  eps2 = 1.0e-5

                  metro = 3000

                  metri = metro * iaa

                  nhic = 10

                  rdevi = 1.0

                  rdevr = 0.8

                  avali = 1.0 + rn(0)

                  mhit = 0
                  lhit = 0

*-----------------------------------------------------------------------

               do 400 i = 1, iaa

                  aval(iatk(i)) = avali

                  rdev(i) = rdevi

                  nhit(i) = 0

  400          continue

*-----------------------------------------------------------------------
*           INITIAL VALUE OF FITTING RKAI**2
*-----------------------------------------------------------------------

                  rkai2 = 0.0

                  do 300 ix = 1, ila

                     do 310 j = 1, mh + 1

                        dval(j) = daxy(j,ix)

  310                continue


                     call func10(rk,dval,
     &                           ipr,ierr,dyval,
     &                           nfun,nnfn,nevl,nevk,nevn)

                                 if( ierr .ne. 0 ) goto 998

                     daxy(jj,ix) = dyval


                     rkai1 = daxy(kdy,ix) - dyval

                     if( iffd(jj) .eq. 2 ) rkai1 = rkai1 / daxy(mh+2,ix)

                     rkai2 = rkai2 + rkai1**2

  300             continue

                     rkai2i = rkai2


*-----------------------------------------------------------------------
*           METROPOLIS FOR RKAI**2
*-----------------------------------------------------------------------

            do 500 imet = 1, metro

            do 501 iran = 1, iaa

                  if( rdev(iran) .lt. eps2 ) goto 501

                     lhit = lhit + 1

                  if( nhit(iran) .gt. nhic ) then

                     rdev(iran) = rdev(iran) * rdevr

                     nhit(iran) = 0

                  end if


                  saval = aval(iatk(iran))

                  skai2 = rkai2

                  aval(iatk(iran)) = saval
     &                             + max( eps1, abs( saval ) )
     &                             * ( 2.0 * rn(0) - 1.0 )
     &                             * rdev(iran)

*-----------------------------------------------------------------------

                  rkai2 = 0.0

                  do 301 ix = 1, ila

                     do 331 j = 1, mh + 1

                        dval(j) = daxy(j,ix)

  331                continue


                     call func10(rk,dval,
     &                           ipr,ierr,dyval,
     &                           nfun,nnfn,nevl,nevk,nevn)

                                 if( ierr .ne. 0 ) goto 998

                     daxy(jj,ix) = dyval


                     rkai1 = daxy(kdy,ix) - dyval

                     if( iffd(jj) .eq. 2 ) rkai1 = rkai1 / daxy(mh+2,ix)

                     rkai2 = rkai2 + rkai1**2

  301             continue

*-----------------------------------------------------------------------

                     rdif = rkai2 - skai2

                  if( rdif .le. 0.0 ) then

                     mhit = mhit +1

                     nhit(iran) = 0

                  else if( rdif .gt. 0.0 ) then

                     aval(iatk(iran)) = saval

                     rkai2 = skai2

                     nhit(iran) = nhit(iran) + 1

                  end if

  501       continue
  500       continue


*-----------------------------------------------------------------------
*        END OF METROPOLIS AND WRITE PARAMETERS ON *
*-----------------------------------------------------------------------

                  write(*,*)
                  write(*,'(''  METRO('',i5,'')  HIT = '',
     &                      i5,'' / '',i5)') metri, mhit, lhit

                  write(*,*)

                  write(*,'(''  KAI**2 = '',g16.8)') rkai2

                  write(*,*)


                  iass = 0

               do 542 j = 1, iaa

                  iakk = mfc + 1

                  do 541 i = 1, iaa

                 if( iatk(i) .lt. iakk .and.
     &               iatk(i) .gt. iass ) then

                    iakk = iatk(i)
                    iacr = iakk

                 end if

  541             continue

                    iass = iakk
                    iata(j) = iacr

  542          continue


               do 540 i = 1, iaa

                  write(*,'(''     A'',i2.2,'' = '',g16.8)')
     &                  iata(i), aval(iata(i))

  540          continue

*-----------------------------------------------------------------------
*        FINAL EVALUATE THE FUNCTION
*-----------------------------------------------------------------------

                  do 302 ix = 1, ila

                     do 312 j = 1, mh + 1

                        dval(j) = daxy(j,ix)

  312                continue

                     call func10(rk,dval,
     &                           ipr,ierr,dyval,
     &                           nfun,nnfn,nevl,nevk,nevn)

                                 if( ierr .ne. 0 ) goto 998

                     daxy(jj,ix) = dyval

  302             continue

*-----------------------------------------------------------------------

            ierr = 0

            return

*-----------------------------------------------------------------------

  997    continue

            m_err = 'H: Fitting Function Has Negative or Zero Error.'
            ErrCha = ''
            ErrID = 'L:3008/R:func09/F:a-func.f'
            l_err = illi
            k_err = jsn

            ierr = 1

            return

*-----------------------------------------------------------------------

  998    continue

            m_err = 'H: Fitting Function Description is Wrong.'
            ErrCha = ''
            ErrID = 'L:3022/R:func09/F:a-func.f'
            l_err = illi
            k_err = jsn

            ierr = 1

            return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine func10(rk,dval,
     &                  ipr,ierr,dyval,
     &                  nfun,nnfn,nevl,nevk,nevn)
*                                                                      *
*      PURPOSE  : EVALUATE ONE COLUMN                                  *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      dimension dval(mh+1)

      common /rval1/ cval(mxcval), aval(mxcval)
      dimension rk(0:mfc)

      dimension eval(0:mfc)
      dimension yval(0:mfc)

      dimension nfun(0:mfc), nnfn(0:mfc), nevl(0:mfc)
      dimension nevk(0:mfc,0:mfc), nevn(0:mfc,0:mfc)

*-----------------------------------------------------------------------
*        DO LOOP FOR ONE PROCEDURE
*-----------------------------------------------------------------------

            do 310 i = 1, ipr

*-----------------------------------------------------------------------
*                 DEFINITION OF FUNCTION
*-----------------------------------------------------------------------

                        nfunc = nfun(i)
                        nnfun = nnfn(i)

*-----------------------------------------------------------------------
*                 DEFINITION OF VARIABLES
*-----------------------------------------------------------------------

                           neval = nevl(i)

                  do 320 j = 1, neval

                     if( nevk(i,j) .eq. 1 ) then

                           eval(j) = dval(nevn(i,j))

                     else if( nevk(i,j) .eq. 2 ) then

                           eval(j) = cval(nevn(i,j))

                     else if( nevk(i,j) .eq. 3 ) then

                        if( rk(nevn(i,j)) .le. r3max ) then ! S.H. xorshift (2020.2.6)

                           eval(j) = rk(nevn(i,j))

                        else

                           eval(j) = rn(0)

                        end if

                     else if( nevk(i,j) .eq. 4 ) then

                           eval(j) = yval(nevn(i,j))

                     else if( nevk(i,j) .eq. 5 ) then

                           eval(j) = -1.0

                     else if( nevk(i,j) .eq. 6 ) then

                           eval(j) = aval(nevn(i,j))

                     end if

  320             continue

*-----------------------------------------------------------------------
*                 EVALUATE ONE PROCEDURE YVAL(I)
*-----------------------------------------------------------------------

               call func05(nfunc,nnfun,neval,eval,yvalp,ierr)

                     if( ierr .ne. 0 ) goto 998

                     yval(i) = yvalp

*-----------------------------------------------------------------------

  310       continue

                     dyval = yval(ipr)


            ierr = 0

            return

*-----------------------------------------------------------------------

  998    continue

            ierr = 1

            return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine func11(idn,kd,ki,rk,ivi,ivn,ifp,ipr,
     &                  kdjj,kijj,rkjj,ivijj,ivnjj,ifpjj,
     &                  iprjj,idrjj)
*                                                                      *
*      PURPOSE  : RE-STORE THE ONE VALUABLE PROCEDURE                  *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      dimension kd(0:mfc), ki(0:mfc)
      dimension rk(0:mfc)
      dimension ivi(0:mfc,0:mfc), ivn(0:mfc), ifp(0:mfc)

      dimension iprjj(mh), idrjj(mh)

      dimension kdjj(mh,0:mfc), kijj(mh,0:mfc)
      dimension rkjj(mh,0:mfc)
      dimension ivnjj(mh,0:mfc), ifpjj(mh,0:mfc)
      dimension ivijj(mh,0:mfc,0:mfc)

*-----------------------------------------------------------------------

            ipr = iprjj(idn)

         do 100 i = 1, ipr

            ifp(i) = ifpjj(idn,i)
            ivn(i) = ivnjj(idn,i)

         do 200 j = 1, ivn(i)

            ivi(j,i) = ivijj(idn,j,i)

  200    continue
  100    continue


            id = idrjj(idn)

         do 300 i = 1, id

            kd(i) = kdjj(idn,i)
            ki(i) = kijj(idn,i)
            rk(i) = rkjj(idn,i)

  300    continue

*-----------------------------------------------------------------------

         return
         end


************************************************************************
*                                                                      *
      subroutine func12(dum,ic,ierr,iclm,duml,cvvv,xval)
*                                                                      *
*                                                                      *
*      PURPOSE  : EVALUATE ONE FUNCTION WITHOUT VARIABLES              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character duml*1

      common /rval1/ cval(mxcval), aval(mxcval)

      dimension dval(mh+1)

*-----------------------------------------------------------------------

      dimension ivi(0:mfc,0:mfc)
      dimension nevk(0:mfc,0:mfc), nevn(0:mfc,0:mfc)

*-----------------------------------------------------------------------

      dimension rk(0:mfc)
      dimension kd(0:mfc), kp(0:mfc), ki(0:mfc), jd(0:mfc)
      dimension md(0:mfc), mg(0:mfc), me(0:mfc), mp(0:mfc)
      dimension ivn(0:mfc), ifp(0:mfc)
      dimension nfun(0:mfc), nnfn(0:mfc), nevl(0:mfc)

*-----------------------------------------------------------------------
*     INITIALIZATION
*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------
*           READ FUNCTION
*-----------------------------------------------------------------------

                  call func01(dum,ic,iclm,ierr,duml,
     &                        rk,kd,kp,ki,jd,mg,md,me,mp,id,iaac)

                        if( ierr .ne. 0 ) goto 911

                  call func03(ierr,
     &                        kd,kp,ki,jd,mg,md,me,mp,id,
     &                        ivi,ivn,ifp,ipr)

                     if( ierr .ne. 0 ) goto 911

*-----------------------------------------------------------------------
*           CHECK FOR PROCEDURE
*-----------------------------------------------------------------------

               do 311 i = 1, ipr

*-----------------------------------------------------------------------
*                 DEFINITION OF FUNCTION
*-----------------------------------------------------------------------

                  if( ifp(i) .eq. 0 ) then

                        nfunc = 5
                        nnfun = 0

                  else if( ifp(i) .eq. -1 ) then

                        nfunc = 1
                        nnfun = 0

                  else if( ifp(i) .gt. 0 ) then

                        nfunc = kd(ifp(i))

                     if( nfunc .ne. 10 ) then

                        nnfun = 0

                     else

                        nnfun = ki(ifp(i))

                     end if

                  else

                        goto 911

                  end if

                        nfun(i) = nfunc
                        nnfn(i) = nnfun

*-----------------------------------------------------------------------
*                 DEFINITION OF VARIABLES
*-----------------------------------------------------------------------

                        neval = ivn(i)
                        nevl(i) = neval

                  do 321 j = 1, neval

                           iiv = ivi(j,i)

*-----------------------------------------------------------------------

                     if( iiv .lt. 0 ) then

                           iiv = - iiv

*-----------------------------------------------------------------------
*                          --- X-VALUE ---
*-----------------------------------------------------------------------

                        if( kd(iiv) .eq. 1 ) then

                           nevk(i,j) = 1
                           nevn(i,j) = 1

*-----------------------------------------------------------------------
*                          --- CONSTANTS ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 3 ) then

                           if ( ki(iiv) .gt. mxcval ) goto 911
                           if( cval(ki(iiv)) .lt. -r0max )
     &                     goto 911

                           nevk(i,j) = 2
                           nevn(i,j) = ki(iiv)

*-----------------------------------------------------------------------
*                          --- REAL NUMBERS ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 4 ) then

                           nevk(i,j) = 3
                           nevn(i,j) = iiv

*-----------------------------------------------------------------------
*                          --- PARAMETERS ---
*-----------------------------------------------------------------------

                        else if( kd(iiv) .eq. 14 ) then

                           if( aval(ki(iiv)) .lt. -r0max )
     &                     goto 911

                           nevk(i,j) = 6
                           nevn(i,j) = ki(iiv)

*-----------------------------------------------------------------------

                        else

                           goto 911

                        end if

*-----------------------------------------------------------------------

                     else if( iiv .gt. 0 ) then

                           if( iiv .ge. i ) goto 911

                           nevk(i,j) = 4
                           nevn(i,j) = iiv

*-----------------------------------------------------------------------

                     else if( iiv .eq. 0 ) then

                           nevk(i,j) = 5
                           nevn(i,j) = 0

                     end if

*-----------------------------------------------------------------------

  321             continue

  311          continue

*-----------------------------------------------------------------------
*           EVALUATE FUNCTION
*-----------------------------------------------------------------------

                     dval(1) = xval

                  call func10(rk,dval,
     &                        ipr,ierr,cvvv,
     &                        nfun,nnfn,nevl,nevk,nevn)

                     if( ierr .ne. 0 ) goto 911

*-----------------------------------------------------------------------
*        ERROR
*-----------------------------------------------------------------------

         return

  911 continue

         ierr = 1

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine setcv(lum,ic0,iclm,dsin,idsi,ill,jsn,ierr)
*                                                                      *
*      PURPOSE  : READ SET CONSTANTS DESCRIPTION                       *
*                                                                      *
*      FORMAT   : SET: C1[23] C2[PI]                                   *
*                                                                      *
*                                                                      *
************************************************************************

      use cvaloutmod ! S.H. 2023.10.27
      use charvarmod, only : ErrLine_Adjust,irwt ! S.H. 2023.11.27

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character lum(ichrl)*1

      common /rval1/ cval(mxcval), aval(mxcval)
      common /paraj/ mstz(300), parz(300)

      dimension ill(0:9)

      character dsin(0:9)*200
      dimension idsi(0:9)

      character m_err*200
      common /error/ m_err, l_err, k_err

      integer l_errtmp
      character ctmp*200

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------
*     INITIALIZATION
*-----------------------------------------------------------------------

            ierr = 0

            ic = ic0

*-----------------------------------------------------------------------
*     READ TEXT
*-----------------------------------------------------------------------

  100       ic  = ic + 1

               if( ic .gt. iclm ) goto 900

               if( lum(ic) .eq. ' ' .or. lum(ic) .eq. tub ) goto 100


*-----------------------------------------------------------------------

            if( lum(ic) .eq. 'c' ) then

*-----------------------------------------------------------------------

               ic = ic + 1

                  call idd2(lum,ic,iclm,idd,ierr)

                  if( ierr .ne. 0 ) then

                      m_err = 'SET: Description of C Specific'//
     &                        ' Number is Wrong.'
                      ErrCha = ''
                      ErrID = 'L:3515/R:setcv/F:a-func.f'
                      l_err = ill(jsn)
                      k_err = jsn

                      goto 999

                  end if

                  if( idd .le. 0 .or. idd .gt. mxcval ) then

                      m_err = 'SET: C Specific Number is Wrong.'
                      ErrCha = ''
                      ErrID = 'L:3527/R:setcv/F:a-func.f'
                      l_err = ill(jsn)
                      k_err = jsn

                      goto 999

                  end if

*-----------------------------------------------------------------------

                  icnum = idd

*-----------------------------------------------------------------------

                  if( lum(ic) .ne. '[' ) then

                      m_err = 'SET: Description of C is Wrong.'//
     &                        ' It should be C2[pi*4.0] for example.'
                      ErrCha = ''
                      ErrID = 'L:3546/R:setcv/F:a-func.f'
                      l_err = ill(jsn)
                      k_err = jsn

                      goto 999

                  end if

*-----------------------------------------------------------------------
*              EVALUATE ONE FUNCTION
*-----------------------------------------------------------------------

                  call func12(lum,ic,ierr,iclm,
     &                        ']',cvvv,xval)

                        if( ierr .ne. 0 ) goto 911


                  cval(icnum) = cvvv

                  if ( mstz(161) .ne. 0 .and.
     &                 imemorycval .ne. 0 .and.
     &                 icvalout .le. iabs(mstz(161)) ) then
                   icvalout = icvalout + 1
                   cvalout(1,icvalout) = icnum
                   cvalout(2,icvalout) = cval(icnum)
                   if( irwt .eq. 0) then ! $RWT=0
                    cvalout(3,icvalout) = ill(jsn)
                    cvalout(4,icvalout) = len_trim(dsin(jsn))
                    cfncvalout(icvalout) = dsin(jsn)(1:idsi(jsn))
                   else ! $RWT>0
                    l_errtmp = ill(jsn)
                    call ErrLine_Adjust(ctmp,ill(jsn))
                    cvalout(3,icvalout) = ill(jsn)
                    cvalout(4,icvalout) = len_trim(ctmp)
                    cfncvalout(icvalout) = ctmp
                    ill(jsn) = l_errtmp
                   end if
                  end if

*-----------------------------------------------------------------------

                  goto 100

*-----------------------------------------------------------------------

            else

                  goto 998

            end if

*-----------------------------------------------------------------------

  900    continue

         return

*-----------------------------------------------------------------------
*        ERROR
*-----------------------------------------------------------------------

  911 continue

                      ierr = 1

                      m_err = 'SET: Description of C[ ] is Wrong.'
                      ErrCha = ''
                      ErrID = 'L:3614/R:setcv/F:a-func.f'
                      l_err = ill(jsn)
                      k_err = jsn

      return

  998 continue

                      ierr = 1

                      m_err = 'SET: Description of C is Wrong.'
                      ErrCha = ''
                      ErrID = 'L:3626/R:setcv/F:a-func.f'
                      l_err = ill(jsn)
                      k_err = jsn

      return

  999 continue
                     ierr = 1
      return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine spln00(noned,noner,syss,ila,ierr,
     &                  idc,ixdcs,iydcs,jtyp,iwi,ipols,
     &                  imars,rcol,rcob,icoma,facsiz,
     &                  jof,jif,ityp,
     &                  dspl)

*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter ( eps = 1.0e-5 )

      dimension dspl(15,*)
      dimension rcol(3), rcob(3)

*-----------------------------------------------------------------------

            ierr = 0

            isys = nint( syss )

*-----------------------------------------------------------------------
*        CHECK CLOSED PATH
*-----------------------------------------------------------------------

            icw = 0

            if( abs( dspl(1,1) - dspl(1,ila) ) .lt. eps .and.
     &          abs( dspl(2,1) - dspl(2,ila) ) .lt. eps ) icw = 1

*-----------------------------------------------------------------------
*           CHECK THE MONOTONICTY OF THE FUNCTION
*-----------------------------------------------------------------------

               if( dspl(1,2) .gt. dspl(1,1) ) then

                  imono = 1

               else if( dspl(1,2) .lt. dspl(1,1) ) then

                  imono = -1

               else

                  imono = 0

               end if

               do 100 i = 3, ila

                  if( imono .eq. 1 ) then

                     if( dspl(1,i) .lt. dspl(1,i-1) ) then

                        imono = 0
                        goto 110

                     end if

                  else if( imono .eq. -1 ) then

                     if( dspl(1,i) .gt. dspl(1,i-1) ) then

                        imono = 0
                        goto 110

                     end if

                  else if( imono .eq. 0 ) then

                     if( dspl(1,i) .gt. dspl(1,i-1) ) then

                        imono = 1

                     else if( dspl(1,i) .lt. dspl(1,i-1) ) then

                        imono = -1

                     end if

                  end if

  100          continue

  110          continue

*-----------------------------------------------------------------------
*        NORMAL OR DISTANCE SPLINE
*-----------------------------------------------------------------------
*           DISTANCE: X, Y ARE SPLINED BY THE DISTANCE
*-----------------------------------------------------------------------

            if( imono .eq. 0 ) then

                     ila0 = ila

                        dspl(3,1) = dspl(1,1)
                        dspl(4,1) = dspl(2,1)

                        dspl(5,1) = 0.0

                        xylng = 0.0

                        ij = 1

                  do 730 ik = 2, ila

                        dist = sqrt(
     &                        + ( dspl(1,ik) - dspl(1,ik-1) )**2
     &                        + ( dspl(2,ik) - dspl(2,ik-1) )**2 )

                     if( dist .gt. 0.0 ) then

                        ij = ij + 1

                        xylng = xylng + dist

                        dspl(5,ij) = xylng

                        dspl(3,ij) = dspl(1,ik)
                        dspl(4,ij) = dspl(2,ik)

                     else

                        ila0 = ila0 - 1

                     end if

  730             continue

                        ilas = ila0 + ( ila0 - 1 ) * abs(isys)

*-----------------------------------------------------------------------
*                 X-VALUE
*-----------------------------------------------------------------------

                  do 300 j = 1, ila0

                     dspl(10,j) = dspl(5,j)
                     dspl(11,j) = dspl(3,j)

  300             continue

                     call spline(dspl,ila0,isys,ierr,icw)

                           if( ierr .ne. 0 ) goto 999

                  do 310 j = 1, ilas

                     dspl(6,j) = dspl(9,j)

  310             continue

*-----------------------------------------------------------------------
*                 Y-VALUE
*-----------------------------------------------------------------------

                  do 320 j = 1, ila0

                     dspl(10,j) = dspl(5,j)
                     dspl(11,j) = dspl(4,j)

  320             continue

                     call spline(dspl,ila0,isys,ierr,icw)

                           if( ierr .ne. 0 ) goto 999

                  do 330 j = 1, ilas

                     dspl(7,j) = dspl(9,j)

  330             continue

*-----------------------------------------------------------------------
*           NORMAL SPLINE
*-----------------------------------------------------------------------

            else

                     ila0 = ila

                        dspl(3,1) = dspl(1,1)
                        dspl(4,1) = dspl(2,1)

                        ij = 1

                  do 740 ik = 2, ila

                     if( dspl(1,ik) .ne. dspl(1,ik-1) ) then

                        ij = ij + 1

                        dspl(3,ij) = dspl(1,ik)
                        dspl(4,ij) = dspl(2,ik)

                     else

                        ila0 = ila0 - 1

                     end if

  740             continue

                        ilas = ila0 + ( ila0 - 1 ) * abs(isys)

                  do 340 j = 1, ila0

                     dspl(10,j) = dspl(3,j)
                     dspl(11,j) = dspl(4,j)

  340             continue

                     call spline(dspl,ila0,isys,ierr,icw)

                           if( ierr .ne. 0 ) goto 999

                  do 350 j = 1, ilas

                     dspl(6,j) = dspl(8,j)
                     dspl(7,j) = dspl(9,j)

  350             continue

            end if

*-----------------------------------------------------------------------
*           WRITE ON THE FILE
*-----------------------------------------------------------------------

                        if( ityp .ne. 11 ) then

                           noned = noned + 1

                           write(jof) ilas, idc, ixdcs, iydcs

                           write(jof)
     &                     jtyp,iwi,ipols,imars,rcol,rcob,icoma,facsiz

                        else

                           noner = noner + 1

                           write(jif) ilas, idc, ixdcs, iydcs

                           write(jif)
     &                     jtyp,iwi,ipols,imars,rcol,rcob,icoma,facsiz

                        end if

                        ddxm = 0.0
                        ddxp = 0.0
                        ddym = 0.0
                        ddyp = 0.0

                     do 726 ij = 1, ilas

                        if( ityp .ne. 11 ) then

                           write(jof)
     &                     dspl(6,ij),dspl(7,ij),ddxm,ddxp,ddym,ddyp

                        else

                           write(jif)
     &                     dspl(6,ij),dspl(7,ij),ddxm,ddxp,ddym,ddyp

                        end if

  726                continue

*-----------------------------------------------------------------------

      return

  999 ierr = 1

      return
      end


************************************************************************
*                                                                      *
      subroutine spline(dspl,ila,isys,ierr,icw)
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      dimension dspl(15,*)

*-----------------------------------------------------------------------

         ierr = 0

         ilas = ila + ( ila - 1 ) * abs(isys)

         sysn = dble(abs(isys)) + 1.0

*-----------------------------------------------------------------------
*     SPLINE CUBE
*-----------------------------------------------------------------------

         if( isys .gt. 0 ) then

*-----------------------------------------------------------------------

                     call spcoef(dspl,ila,icw)

               do 400 i = 1, ila - 1

                     j0 = ( i - 1 ) * ( isys + 1 ) + 1

                     dspl(8,j0) = dspl(10,i)
                     dspl(9,j0) = dspl(11,i)

                     dist = ( dspl(10,i+1) - dspl(10,i) ) / sysn

                  do 500 j = 1, isys

                     ii = j0 + j

                     xin = dspl(10,i) + dble(j) * dist

                     dspl(8,ii) = xin
                     dspl(9,ii) = splcub(dspl,ila,xin)

  500             continue

  400          continue

                     dspl(8,ilas) = dspl(10,ila)
                     dspl(9,ilas) = dspl(11,ila)

*-----------------------------------------------------------------------
*     LAGRANGE 4 TH
*-----------------------------------------------------------------------

         else if( isys .lt. 0 ) then

               isys = - isys

               if( ila .lt. 5 ) goto 999

               do 100 i = 1, ila - 1

                     j0 = ( i - 1 ) * ( isys + 1 ) + 1

                     dspl(8,j0) = dspl(10,i)
                     dspl(9,j0) = dspl(11,i)

                     dist = ( dspl(10,i+1) - dspl(10,i) ) / sysn

                  do 200 j = 1, isys

                     ii = j0 + j

                     xin = dspl(10,i) + dble(j) * dist

                     call splag4(dspl,i,ila,xin,yout,ierr)

                        if( ierr .ne. 0 ) goto 999

                     dspl(8,ii) = xin
                     dspl(9,ii) = yout

  200             continue

  100          continue

                     dspl(8,ilas) = dspl(10,ila)
                     dspl(9,ilas) = dspl(11,ila)

            end if

*-----------------------------------------------------------------------

      return

  999 ierr = 1

      return
      end


************************************************************************
*                                                                      *
      subroutine spcoef(dspl,n,icw)
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      dimension dspl(15,*)

*-----------------------------------------------------------------------
*        RE-ORDERING THE X-VALUE
*-----------------------------------------------------------------------

         do 100 i = 1, n

            dspl(15,i) = dble( i )

  100    continue

         do 300 i = 1, n - 1

               ip1 = i + 1

            do 200 j = ip1, n

               ii = nint( dspl(15,i) )
               ij = nint( dspl(15,j) )

               if( dspl(10,ii) .le. dspl(10,ij) ) goto 200

               ritemp     = dspl(15,i)
               dspl(15,i) = dspl(15,j)
               dspl(15,j) = ritemp

  200       continue

  300    continue

*-----------------------------------------------------------------------
*        FIRST DERIVATIVE IS GIVEN OR NOT
*-----------------------------------------------------------------------

            if( icw .eq. 0 ) then

               yp1 = 0.0
               ypn = 0.0

               dspl(12,1) = 0.0
               dspl(13,1) = 0.0

            else

               yp1 = ( dspl(11,2) - dspl(11,n-1) )
     &             / ( dspl(10,2) - dspl(10,1)
     &               + dspl(10,n) - dspl(10,n-1) )

               ypn = yp1

               dspl(12,1) = -0.5
               dspl(13,1) = ( 3.0 / ( dspl(10,2) - dspl(10,1) ) )
     &                    * ( ( dspl(11,2) - dspl(11,1) )
     &                    / ( dspl(10,2) - dspl(10,1) ) - yp1 )

            end if

*-----------------------------------------------------------------------

            do 400 i = 2, n - 1

               ii = nint( dspl(15,i) )
               im = nint( dspl(15,i-1) )
               ip = nint( dspl(15,i+1) )

               hip = ( dspl(10,ii) - dspl(10,im) )
     &             / ( dspl(10,ip) - dspl(10,im) )

               temp = hip * dspl(12,i-1) + 2.0

               dspl(12,i) = ( hip - 1.0 ) / temp

               dspl(13,i) = ( 6.0 * ( ( dspl(11,ip) - dspl(11,ii) )
     &                              / ( dspl(10,ip) - dspl(10,ii) )
     &                              - ( dspl(11,ii) - dspl(11,im) )
     &                              / ( dspl(10,ii) - dspl(10,im) ) )
     &                              / ( dspl(10,ip) - dspl(10,im) )
     &                              - hip * dspl(13,i-1) ) / temp

  400       continue

*-----------------------------------------------------------------------

               nn = nint( dspl(15,n) )
               nm = nint( dspl(15,n-1) )

            if( icw .eq. 0 ) then

               qn = 0.0
               un = 0.0

            else

               qn = 0.5
               un = ( 3.0 / ( dspl(10,nn) - dspl(10,nm) ) )
     &            * ( ypn - ( dspl(11,nn) - dspl(11,nm) )
     &            / ( dspl(10,nn) - dspl(10,nm) ) )

            end if

               dspl(12,n) = ( un - qn * dspl(13,n-1) )
     &                    / ( qn * dspl(12,n-1) + 1.0 )

*-----------------------------------------------------------------------

            do 500 i = n - 1, 1, -1

               dspl(12,i) = dspl(12,i) * dspl(12,i+1) + dspl(13,i)

  500       continue

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function splcub(dspl,n,x)
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      dimension dspl(15,*)

*-----------------------------------------------------------------------

               i1 = nint( dspl(15, 1 ) )

         if( x .lt. dspl(10,i1) ) then

               i2 = nint( dspl(15, 2 ) )
               h1 = dspl(10,i2) - dspl(10,i1)

               splcub = dspl(11,i1)
     &                + ( x - dspl(10,i1) )
     &                * ( ( dspl(11,i2) - dspl(11,i1) )
     &                / h1 - h1 * dspl(12,2) / 6.0 )

               return

         end if

*-----------------------------------------------------------------------

               in = nint( dspl(15, n ) )

         if( x .gt. dspl(10,in) ) then

               inm1 = nint( dspl(15, n - 1 ) )
               hnm1 = dspl(10,in) - dspl(10,inm1)

               splcub = dspl(11,in)
     &                + ( x -dspl(10,in) )
     &                * ( ( dspl(11,in) - dspl(11,inm1) ) / hnm1
     &                   + hnm1 * dspl(12,n-1) / 6.0 )

               return

         end if

*-----------------------------------------------------------------------

         do 800 i = 2, n

               ii = nint( dspl(15, i ) )

               if( x .le. dspl(10,ii) ) goto 900

  800    continue

  900          l = i - 1

               il   = nint( dspl(15, l ) )
               ilp1 = nint( dspl(15, l + 1 ) )

               a = dspl(10,ilp1) - x
               b = x - dspl(10,il)

               hl = dspl(10,ilp1) - dspl(10,il)

               splcub = a * dspl(12,l)
     &                * ( a**2 / hl - hl ) / 6.0
     &                + b * dspl(12,l+1) * ( b**2 / hl - hl ) / 6.0
     &                + ( a * dspl(11,il) + b * dspl(11,ilp1) ) / hl

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine splag4(dspl,in,ila,xin,yout,ierr)
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      dimension dspl(15,*)

      dimension isel(5)
      dimension v02(5), v03(5), v04(5), v05(5)

*-----------------------------------------------------------------------
*     LAGRANGE 4 TH
*-----------------------------------------------------------------------

         ierr = 0

         do 50 i = 1, ila

            dspl(14,i) = dspl(10,i)

   50    continue

*-----------------------------------------------------------------------
*     SEARCH MINIMUM 5 TH POINTS
*-----------------------------------------------------------------------

               if( dspl(10,in) .eq. dspl(10,in+1) ) goto 999

               if( abs( dspl(10,in)   - xin ) .lt.
     &             abs( dspl(10,in+1) - xin ) ) then

                  isel(1) = in
                  isel(2) = in + 1

                  isec = 1

               else

                  isel(1) = in + 1
                  isel(2) = in

                  isec = 2

               end if

*-----------------------------------------------------------------------

               if( isec .eq. 1 ) then

                  if( in - 1 .lt. 1 ) then

                              isel(3) = in + 2
                              isel(4) = in + 3
                              isel(5) = in + 4

                  else

                              isel(3) = in - 1

                     if( in + 2 .gt. ila ) then

                              isel(4) = in - 2
                              isel(5) = in - 3

                     else

                              isel(4) = in + 2

                        if( in - 2 .lt. 1 ) then

                              isel(5) = in + 3

                        else

                              isel(5) = in - 2

                        end if

                     end if

                  end if

               else if( isec .eq. 2 ) then

                  if( in + 2 .gt. ila ) then

                              isel(3) = in - 1
                              isel(4) = in - 2
                              isel(5) = in - 3

                  else

                              isel(3) = in + 2

                     if( in - 1 .lt. 1 ) then

                              isel(4) = in + 3
                              isel(5) = in + 4

                     else

                              isel(4) = in - 1

                        if( in + 3 .gt. ila ) then

                              isel(5) = in - 2

                        else

                              isel(5) = in + 3

                        end if

                     end if

                  end if

               end if

*-----------------------------------------------------------------------
*     LAGRANGE 4 TH
*-----------------------------------------------------------------------

      do 300 i = 2, 5

            v02(i) = ( dspl(11,isel(1)) * ( dspl(14,isel(i)) - xin )
     &               - dspl(11,isel(i)) * ( dspl(14,isel(1)) - xin ) )
     &             / ( dspl(14,isel(i)) - dspl(14,isel(1)) )

  300 continue

*-----------------------------------------------------------------------

      do 301 i = 3, 5

            v03(i) = ( v02(2) * ( dspl(14,isel(i)) - xin )
     &               - v02(i) * ( dspl(14,isel(2)) - xin ) )
     &             / ( dspl(14,isel(i)) - dspl(14,isel(2)) )

  301 continue

*-----------------------------------------------------------------------

      do 302 i = 4, 5

            v04(i) = ( v03(3) * ( dspl(14,isel(i)) - xin )
     &               - v03(i) * ( dspl(14,isel(3)) - xin ) )
     &             / ( dspl(14,isel(i)) - dspl(14,isel(3)) )

  302 continue

*-----------------------------------------------------------------------

            v05(5) = ( v04(4) * ( dspl(14,isel(5)) - xin )
     &               - v04(5) * ( dspl(14,isel(4)) - xin ) )
     &             / ( dspl(14,isel(5)) - dspl(14,isel(4)) )

            yout = v05(5)

*-----------------------------------------------------------------------

            if( abs(yout) .gt.
     &          100.0 * max(1.0d0,abs(dspl(11,in)),abs(dspl(11,in+1))) )
     &      goto 999

      return

  999 ierr = 1

      return
      end


