************************************************************************
*                                                                      *
      subroutine jpncode(chaf,icn,ifon)
*                                                                      *
*       PURPOSE : DETERMINE THE JAPANESE KANJI CODE                    *
*                                                                      *
*                 IFON = 1 : SJIS : 128 < ONE BITE < 161               *
*                 IFON = 2 : JIS  : 27 36 66(64),,66(74)               *
*                 IFON = 3 : EUC  : ONE BITE > 161                     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1

*-----------------------------------------------------------------------

            if( ifon .eq. 1 .or. ifon .eq. 2 ) return

                  i = 0

  101       i = i + 1

                  irr = ichar( chaf(i) )

                  if( irr .gt. 128 .and. irr .lt. 161 ) then

                     ifon = 1

                     return

                  else if( irr .eq. 27 .and.
     &                     ichar(chaf(i+1)) .eq. 36 .and.
     &                   ( ichar(chaf(i+2)) .eq. 66  .or.
     &                     ichar(chaf(i+2)) .eq. 64 ) ) then

                     ifon = 2

                     return

                  else if( irr .gt. 128 ) then

                     ifon = 3
                     i = i + 1

                  end if

            if( i .lt. icn ) goto 101

      return
      end

************************************************************************
*                                                                      *
      subroutine jpnprep(chaf,icn,ifon)
*                                                                      *
*              PRE-PROCEDURE FOR JAPANESE KANJI                        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(ichrl)*1

*-----------------------------------------------------------------------

      character yen*1

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------

            yen = char(92)

*-----------------------------------------------------------------------

               if( ifon .le. 0 ) return

*-----------------------------------------------------------------------
*     INSERT \377\001 AND \377\000 BEFORE AND AFTER THE KANJI
*-----------------------------------------------------------------------

                  i = 0
                  j = 0

                  k = 0

  101       i = i + 1

                  if( i .gt. icn ) goto 1000

                  irr = ichar( chaf(i) )

         if( ifon .eq. 1 .or. ifon .eq. 3 ) then

            if( k .eq. 0 .and. irr .gt. 128 ) then

                  if( j + 9 .gt. ichrl ) goto 1000

                  chag(j+1) = yen
                  chag(j+2) = '3'
                  chag(j+3) = '7'
                  chag(j+4) = '7'
                  chag(j+5) = yen
                  chag(j+6) = '0'
                  chag(j+7) = '0'
                  chag(j+8) = '1'

                  chag(j+9) = chaf(i)

                  j = j + 9
                  i = i + 1

                  k = 1

            else if( k .eq. 1 .and. irr .gt. 128 ) then

                  if( j + 1 .gt. ichrl ) goto 1000

                  chag(j+1) = chaf(i)

                  j = j + 1
                  i = i + 1

            else if( k .eq. 1 .and. irr .le. 128 ) then

                  if( j + 8 .gt. ichrl ) goto 1000

                  chag(j+1) = yen
                  chag(j+2) = '3'
                  chag(j+3) = '7'
                  chag(j+4) = '7'
                  chag(j+5) = yen
                  chag(j+6) = '0'
                  chag(j+7) = '0'
                  chag(j+8) = '0'

                  j = j + 8

                  k = 0

            end if

         else if( ifon .eq. 2 ) then

            if( irr .eq. 27 .and.
     &          ichar(chaf(i+1)) .eq. 36 .and.
     &        ( ichar(chaf(i+2)) .eq. 66 .or.
     &          ichar(chaf(i+2)) .eq. 64 ) ) then

                  if( j + 8 .gt. ichrl ) goto 1000

                  chag(j+1) = yen
                  chag(j+2) = '3'
                  chag(j+3) = '7'
                  chag(j+4) = '7'
                  chag(j+5) = yen
                  chag(j+6) = '0'
                  chag(j+7) = '0'
                  chag(j+8) = '1'

                  j = j + 8
                  i = i + 3

                  k = 1

            else
     &      if( irr .eq. 27 .and.
     &          ichar(chaf(i+1)) .eq. 40 .and.
     &        ( ichar(chaf(i+2)) .eq. 66 .or.
     &          ichar(chaf(i+2)) .eq. 74 ) ) then

                  if( j + 8 .gt. ichrl ) goto 1000

                  chag(j+1) = yen
                  chag(j+2) = '3'
                  chag(j+3) = '7'
                  chag(j+4) = '7'
                  chag(j+5) = yen
                  chag(j+6) = '0'
                  chag(j+7) = '0'
                  chag(j+8) = '0'

                  j = j + 8
                  i = i + 3

                  k = 0

            end if

         end if

*-----------------------------------------------------------------------

            if( j + 1 .gt. ichrl ) goto 1000

            if( i .le. icn ) then

                  j = j + 1
                  chag(j) = chaf(i)

                  goto 101

            end if

*-----------------------------------------------------------------------
*     REWRITE THE TEXT
*-----------------------------------------------------------------------

 1000       continue


            do 200 i = 1, j

                  chaf( i ) = chag( i )

  200       continue

                  icn = j

      return
      end


************************************************************************
*                                                                      *
      function spjpn(i,icn,icf,chaf,ifon,njpn)
*                                                                      *
*              SKIP JAPANESE TWO BITE CHARACTERS                       *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1

      logical spjpn

*-----------------------------------------------------------------------

      character yen*1

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------

            yen = char(92)

*-----------------------------------------------------------------------

               icf  = i
               j    = i

               njpn = 0

*-----------------------------------------------------------------------

            if( ifon .gt. 0 .and.
     &          icn .ge. 7 .and. i .le. icn - 6 ) then

               if( chaf(j)   .eq. yen .and.
     &             chaf(j+1) .eq. '3' .and.
     &             chaf(j+2) .eq. '7' .and.
     &             chaf(j+3) .eq. '7' .and.
     &             chaf(j+4) .eq. yen .and.
     &             chaf(j+5) .eq. '0' .and.
     &             chaf(j+6) .eq. '0' .and.
     &             chaf(j+7) .eq. '1' ) then

                   j = j + 7

  100              j = j + 1

                   njpn = njpn + 1

                  if( chaf(j)   .eq. yen .and.
     &                chaf(j+1) .eq. '3' .and.
     &                chaf(j+2) .eq. '7' .and.
     &                chaf(j+3) .eq. '7' .and.
     &                chaf(j+4) .eq. yen .and.
     &                chaf(j+5) .eq. '0' .and.
     &                chaf(j+6) .eq. '0' .and.
     &                chaf(j+7) .eq. '0' ) then

                      icf = j + 7

                      njpn = ( njpn - 1 ) / 2

                     spjpn = .true.
                     return

                  end if

                  if( j .lt. icn ) goto 100

                      icf = icn

                      njpn = njpn / 2

                     spjpn = .true.
                     return

               end if

            end if

                  spjpn = .false.

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine prtex0(chaf,icn,ifon)
*                                                                      *
*              PRE-PRE-PROCEDURE FOR TEX MODE                          *
*              TEX MODE IS SPECIFIED BETWEEN CHARACTER $               *
*              SELECT STRING BETWEEN $ AND $                           *
*                                                                      *
*              INSERT '\' BEFORE '(' AND ')'                           *
*                                                                      *
*              TAB TO 8 SPACE                                          *
*                                                                      *
*              INSERT '{' AND '}' FOR SPECIAL CHARACTER AFTER '\'      *
*              LISTED ABOVE WITHOUT '{' JUST BEFORE THE '\'            *
*                                                                      *
*              TRANSLATE NEW MACRO FOR TEX                             *
*                 \TeX, \LaTeX, \ANGEL, \PHITS, \ll, \gg               *
*                 \stackrel{ }{ } -> \mathop{ }^{ }                    *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chat(ichrl)*1
      character chag(ichrl)*1
      dimension ie(0:inig)

      logical spjpn
      logical prepc
      logical premc
      logical preqc

*-----------------------------------------------------------------------

      character yen*1
      character tub*1
      yen = char(92)
      tub = char(9)

*-----------------------------------------------------------------------

               do 20 i = 0, inig

                  ie( i ) = -1

   20          continue

*-----------------------------------------------------------------------
*        INSERT '\' BEFORE '(' OR ')' SKIPPING THE JAPANESE CODE
*        CHANGE TAB TO SPACE
*        INSERT { AND } FOR \large and so on WITHOUT { JUST BEFORE \
*-----------------------------------------------------------------------

               i  = 0
               j  = 0
               ib = 0

  100       continue

               ipsp = 0
               if( i .gt. 0 ) then
                  if( chaf(i) .eq. '{' ) ipsp = 1
               end if

               i = i + 1

*-----------------------------------------------------------------------

            if( spjpn(i,icn,icf,chaf,ifon,njpn) ) then

*              SKIP JAPANESE TWO BITE CHARACTERS

*-----------------------------------------------------------------------

                  do k = i, icf

                        j = j + 1
                        chag( j ) = chaf( k )

                  end do

                        i = icf

*-----------------------------------------------------------------------
*           SPECIAL CASE OF '(' AND ')'
*-----------------------------------------------------------------------

            else if( chaf(i)   .eq. yen .and.
     &               chaf(i+1) .eq. 'r' .and.
     &               chaf(i+2) .eq. 'i' .and.
     &               chaf(i+3) .eq. 'g' .and.
     &               chaf(i+4) .eq. 'h' .and.
     &               chaf(i+5) .eq. 't' .and.
     &               chaf(i+6) .eq. ')' ) then

                     do m = 0, 6

                        j = j + 1
                        chag(j) = chaf( i + m )

                     end do

                        i = i + 6

            else if( chaf(i)   .eq. yen .and.
     &               chaf(i+1) .eq. 'l' .and.
     &               chaf(i+2) .eq. 'e' .and.
     &               chaf(i+3) .eq. 'f' .and.
     &               chaf(i+4) .eq. 't' .and.
     &               chaf(i+5) .eq. '(' ) then

                     do m = 0, 5

                        j = j + 1
                        chag(j) = chaf( i + m )

                     end do

                        i = i + 5

            else if( chaf(i)   .eq. yen .and.
     &             ( chaf(i+1) .eq. 'B' .or.
     &               chaf(i+1) .eq. 'b' ) .and.
     &               chaf(i+2) .eq. 'i' .and.
     &               chaf(i+3) .eq. 'g' .and.
     &             ( chaf(i+4) .eq. '(' .or.
     &               chaf(i+4) .eq. ')' ) ) then

                     do m = 0, 4

                        j = j + 1
                        chag(j) = chaf( i + m )

                     end do

                        i = i + 4

            else if( chaf(i)   .eq. yen .and.
     &             ( chaf(i+1) .eq. '(' .or. chaf(i+1) .eq. ')' .or.
     &               chaf(i+1) .eq. '{' .or. chaf(i+1) .eq. '}' )
     &              ) then

                     do m = 0, 1

                        j = j + 1
                        chag(j) = chaf( i + m )

                     end do

                        i = i + 1

*-----------------------------------------------------------------------
*           INSERT '\' BEFORE '(' AND ')'
*-----------------------------------------------------------------------

            else if( chaf(i) .eq. '(' .or. chaf(i) .eq. ')' ) then

                        j = j + 1
                        chag( j ) = yen
                        j = j + 1
                        chag( j ) = chaf(i)

*-----------------------------------------------------------------------
*           TAB TO 8 SPACE
*-----------------------------------------------------------------------

            else if( chaf(i) .eq. tub ) then

                     do m = 1, 8

                        j = j + 1
                        chag( j ) = ' '

                     end do

*-----------------------------------------------------------------------
*           CHACK PARENTHESIS
*-----------------------------------------------------------------------

            else if( chaf(i) .eq. '{' ) then

                        ib = ib + 1

                        j = j + 1
                        chag( j ) = chaf(i)

            else if( chaf(i) .eq. '}' ) then

  300                   continue

                     if( ie(ib) .eq. ib ) then

                        j = j + 1
                        chag( j ) = '}'

                        ie(ib) = -1

                        ib = ib - 1

                        goto 300

                     end if

                        ib = ib - 1

                        if( ib .lt. 0 ) ib = 0

                        j = j + 1
                        chag( j ) = chaf(i)

*-----------------------------------------------------------------------

*           INSERT '{' AND '}' FOR SPECIAL CHARACTER AFTER '\'
*           LISTED ABOVE WITHOUT '{' JUST BEFORE THE '\'

*-----------------------------------------------------------------------

            else if( prepc(i,j,chaf,chag,icn,ie,ib,ipsp) ) then

*-----------------------------------------------------------------------

*           INSERT '{' AND '}' FOR SPECIAL CHARACTER AFTER '\'
*           LISTED ABOVE WITHOUT '{' JUST AFTER THE '\'

*-----------------------------------------------------------------------

            else if( preqc(i,j,chaf,chag,icn) ) then

*-----------------------------------------------------------------------

*           TRANSLATE NEW MACRO FOR TEX
*                 \TeX, \LaTeX, \ANGEL, \PHITS, \ll, \gg
*                 \stackrel{ }{ } -> \mathop{ }^{ }

*-----------------------------------------------------------------------

            else if( premc(i,j,chaf,chag,icn,ifon) ) then

*-----------------------------------------------------------------------

            else

                        j = j + 1
                        chag( j ) = chaf(i)

            end if

               if( i .lt. icn ) goto 100

*-----------------------------------------------------------------------

               icn = j

               do 110 i = 1, icn

                        chaf(i) = chag(i)

  110          continue

*-----------------------------------------------------------------------
*        CHANGE TEX MODE SKIPPING THE JAPANESE CHARACTERS
*-----------------------------------------------------------------------

               idt = 0
               idi = 1
               i   = 0
               j   = 0

  200          i = i + 1

*-----------------------------------------------------------------------

            if( spjpn(i,icn,icf,chaf,ifon,njpn) ) then

*              SKIP JAPANESE TWO BITE CHARACTERS

*-----------------------------------------------------------------------

                  do 510 k = i, icf

                        j = j + 1
                        chag( j ) = chaf( k )

  510             continue

                        i = icf

*-----------------------------------------------------------------------

            else if( chaf(i)   .eq. yen .and.
     &               chaf(i+1) .eq. '$' ) then

                     do m = 0, 1

                        j = j + 1
                        chag(j) = chaf( i + m )

                     end do

                        i = i + 1

            else if( idt .eq. 0 .and. chaf(i) .eq. '$' ) then

                  idt = 1
                  idi = j + 1


            else if( idt .eq. 1 .and.
     &             ( chaf(i) .eq. '$' .or. i .eq. icn ) ) then

                     if( i .eq. icn .and. chaf(i) .ne. '$' ) then

                        j = j + 1
                        chag( j ) = chaf( i )

                     end if

                  idt = 0
                  idf = j

                  jcn = idf - idi + 1

                  do 220 k = idi, idf

                     chat( k - idi + 1 ) = chag(k)

  220             continue


                     call prtex(chat,jcn,ifon)


                  do 230 k = 1, jcn

                     chag( idi - 1 + k ) = chat(k)

  230             continue

                     j = idi - 1 + jcn

            else

                     j = j + 1
                     chag( j ) = chaf( i )

            end if

               if( i .lt. icn ) goto 200

*-----------------------------------------------------------------------

               icn = j

               do 210 i = 1, icn

                        chaf(i) = chag(i)

  210          continue


*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function premc(k,l,chaf,chag,icn,ifon)
*                                                                      *
*           TRANSLATE NEW MACRO FOR TEX                                *
*                 \TeX, \LaTeX, \ANGEL, \PHITS, \ll, \gg               *
*                 \stackrel{ }{ } -> \mathop{ }^{ }                    *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter ( icmax = 9, iscn = 5 )

      parameter ( ic02 = 4 )
      parameter ( ic03 = 1 )
      parameter ( ic05 = 3 )
      parameter ( ic09 = 1 )

      parameter ( ic01  = 0 )
      parameter ( ic04  = 0 )
      parameter ( ic06  = 0 )
      parameter ( ic07  = 0 )
      parameter ( ic08  = 0 )

*-----------------------------------------------------------------------

      logical premc
      logical spjpn

      character chaf(ichrl)*1
      character chag(ichrl)*1

*-----------------------------------------------------------------------

      character  isex(icmax,iscn)*20
      character  isch(icmax,iscn)*200
      dimension  iscl(icmax,iscn)
      dimension  nic(icmax)

      character dum20*20

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic02),i=2,2)
     &         /'ll',
     &          'gg',
     &          'AA',
     &          'aa'/
      data ((iscl(i,j),j=1,ic02),i=2,2)
     &         /8,
     &          8,
     &          7,
     &          7/
      data ((isch(i,j),j=1,ic02),i=2,2)
     &         /'@lt@!@lt',
     &          '@gt@!@gt',
     &          '@Ang{A}',
     &          '@ang{a}'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic03),i=3,3)
     &         /'TeX'/
      data ((iscl(i,j),j=1,ic03),i=3,3)
     &         /23/
      data ((isch(i,j),j=1,ic03),i=3,3)
     &         /'{@rmT_{@!{@LARGEE}@!}X}'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic05),i=5,5)
     &         /'LaTeX',
     &          'ANGEL',
     &          'PHITS'/
      data ((iscl(i,j),j=1,ic05),i=5,5)
     &         /97,
     &          75,
     &          55/
      data ((isch(i,j),j=1,ic05),i=5,5)
     &         /'{@rmL_{{@040@>@!@!@!@!@!@!}^{@tiny{@small{@040}^{{@Huge
     &{@Huge{@LARGEA@!}}}}}}}T_{@!{@LARGEE}@!}X}',
     &          '{@Large{@rm{@itA{@TINY{@tiny@!}{@TINYN{@tiny@!}}G}_{@!{
     &@largeE}}{@TINYL}}}}',
     &       '{@Large{@rm{@itP{@TINY{@tiny@!}HI}T{@TINY{@tiny@!}S}}}}'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic09),i=9,9)
     &         /'stackrel{'/
      data ((iscl(i,j),j=1,ic09),i=9,9)
     &         /8/
      data ((isch(i,j),j=1,ic09),i=9,9)
     &         /'@mathop{'/

*-----------------------------------------------------------------------

      data nic/ ic01, ic02, ic03, ic04, ic05, ic06, ic07, ic08, ic09 /

*-----------------------------------------------------------------------

      data iyenn / 0 /
      save iyenn

      character yen*1

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------

            yen = char(92)

*-----------------------------------------------------------------------

      if( iyenn .eq. 0 ) then

         do i = 1, icmax
         do j = 1, iscn
         do m = 1, 200

            if( isch(i,j)(m:m) .eq. '@' ) isch(i,j)(m:m) = yen

         end do
         end do
         end do

         iyenn = 1

      end if

*-----------------------------------------------------------------------
*           INSERT '{' AND '}' FOR SPECIAL CHARACTER AFTER '\'
*           LISTED ABOVE WITHOUT '\' JUST BEFORE THE '\'
*-----------------------------------------------------------------------

      if( chaf(k) .eq. yen ) then

                  ichl = 0
                  ichn = 0

                  jfin = min( icn, k+icmax )

                  dum20 = chaf(k+1)

*-----------------------------------------------------------------------

               do 102 j = k + 2, jfin

                     jj = j - k

                     dum20 = dum20(1:jj-1)//chaf(j)

                  if( nic(jj) .gt. 0 ) then

                     do 201 i = 1, nic(jj)

                        if( dum20(1:jj) .eq. isex(jj,i) ) then

                           ichl = jj
                           ichn = i

                        end if

  201                continue

                  end if


  102          continue

*-----------------------------------------------------------------------
*        WRITE NEW MACRO
*-----------------------------------------------------------------------

         if( ichl .gt. 0 ) then

*-----------------------------------------------------------------------
*           FOR \stackrel{ }{ }  to  \mathop{ }^{ }
*-----------------------------------------------------------------------

            if( ichl .eq. 9 .and. ichn .eq. 1 ) then

*-----------------------------------------------------------------------

                     kk = k + ichl

                     icc = 0

                     is1 = kk + 1

                     is2 = -1
                     is3 = -1
                     is4 = -1

                     i = kk

  100                i = i + 1

                     if( i .gt. icn ) then

                           icc = -1
                           goto 101

                     end if

*-----------------------------------------------------------------------

                     if( spjpn(i,icn,icf,chaf,ifon,njpn) ) then

*                       SKIP JAPANESE TWO BITE CHARACTERS

                           i = icf

                     else if( chaf(i) .eq. yen .and.
     &                      ( chaf(i+1) .eq. '{' .or.
     &                        chaf(i+1) .eq. '}' ) ) then

                           i = i + 1

                     else if( chaf(i) .eq. '}' ) then

                        if( is3 .lt. 0 ) then

                           if( is2 .ge. 0 ) then

                              icc = -1
                              goto 101

                           end if

                           is2 = i - 1

                        else if( is3 .ge. 0 ) then

                           is4 = i - 1

                           goto 101

                        end if

                     else if( chaf(i) .eq. '{' ) then

                        if( is3 .ge. 0 ) then

                           icc = -1
                           goto 101

                        else if( is2 .lt. 0 ) then

                           icc = -1
                           goto 101

                        else if( is2 .ge. 0 ) then

                           is3 = i + 1

                        end if

                     end if

                           goto 100

*-----------------------------------------------------------------------

  101             continue

                  if( icc .lt. 0 ) then

                     goto 109

                  else

                     if( is2 .lt. is1 .or. is4 .lt. is3 ) goto 109

                     do 301 ll = 1, iscl(ichl,ichn)

                        l = l + 1
                        chag( l ) = isch(ichl,ichn)(ll:ll)

  301                continue

                     do 302 ll = is3, is4

                        l = l + 1
                        chag( l ) = chaf(ll)

  302                continue

                        l = l + 1
                        chag( l ) = '}'
                        l = l + 1
                        chag( l ) = '^'
                        l = l + 1
                        chag( l ) = '{'

                     do 303 ll = is1, is2

                        l = l + 1
                        chag( l ) = chaf(ll)

  303                continue

                        l = l + 1
                        chag( l ) = '}'

                        k = is4 + 1

                  end if

                     goto 110

*-----------------------------------------------------------------------
*              ERROR
*-----------------------------------------------------------------------

  109             continue

                        l = l + 1
                        chag( l ) = yen
                        l = l + 1
                        chag( l ) = yen

  110             continue

*-----------------------------------------------------------------------
*           END OF \stackrel
*-----------------------------------------------------------------------

            else

                     do 300 ll = 1, iscl(ichl,ichn)

                        l = l + 1
                        chag( l ) = isch(ichl,ichn)(ll:ll)

  300                continue

                        k = k + ichl

                        if( chaf(k+1) .eq. ' ' ) k = k + 1

            end if

*-----------------------------------------------------------------------

                  premc = .true.
                  return

         end if

*-----------------------------------------------------------------------

      end if

                  premc = .false.

      return
      end

************************************************************************
*                                                                      *
      function preqc(k,l,chaf,chag,icn)
*                                                                      *
*           INSERT '{' AND '}' FOR SPECIAL CHARACTER AFTER '\'         *
*           LISTED ABOVE WITHOUT '{' JUST AFTER THE '\'                *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      logical preqc

      parameter ( icmax = 14, iscn = 6 )
      parameter ( ic03 = 4 )
      parameter ( ic04 = 3 )
      parameter ( ic05 = 6 )
      parameter ( ic06 = 1 )
      parameter ( ic08 = 2 )
      parameter ( ic09 = 2 )
      parameter ( ic013 = 1 )
      parameter ( ic014 = 1 )

      parameter ( ic01  = 0 )
      parameter ( ic02  = 0 )
      parameter ( ic07  = 0 )
      parameter ( ic010 = 0 )
      parameter ( ic011 = 0 )
      parameter ( ic012 = 0 )

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(ichrl)*1

*-----------------------------------------------------------------------

      character isex(icmax,iscn)*20
      dimension nic(icmax)

      character dum20*20

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic03),i=3,3)
     &         /'hat','bar','dot','vec'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic04),i=4,4)
     &         /'ddot','sqrt','frac'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic05),i=5,5)
     &         /'check','breve','acute','grave','tilde','fracs'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic06),i=6,6)
     &         /'widehat'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic08),i=8,8)
     &         /'overline','stackrel'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic09),i=9,9)
     &         /'widetilde','underline'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic013),i=13,13)
     &         /'overleftarrow'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic014),i=14,14)
     &         /'overrightarrow'/

*-----------------------------------------------------------------------

      data nic/ ic01, ic02, ic03, ic04, ic05, ic06, ic07, ic08,
     &          ic09, ic010,ic011,ic012,ic013,ic014/

*-----------------------------------------------------------------------

      character yen*1

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------

            yen = char(92)

*-----------------------------------------------------------------------
*           INSERT '{' AND '}' FOR SPECIAL CHARACTER AFTER '\'
*           LISTED ABOVE WITHOUT '{' JUST AFTER THE '\'
*-----------------------------------------------------------------------

         if( chaf(k) .eq. yen ) then

                  ichl = 0
                  ichn = 0

                  jfin = min( icn, k+icmax )

                  dum20 = chaf(k+1)

*-----------------------------------------------------------------------

               do 102 j = k + 2, jfin

                     jj = j - k

                     dum20 = dum20(1:jj-1)//chaf(j)

                  if( nic(jj) .gt. 0 ) then

                     do 201 i = 1, nic(jj)

                        if( dum20(1:jj) .eq. isex(jj,i) ) then

                           ichl = jj
                           ichn = i

                        end if

  201                continue

                  end if


  102          continue

*-----------------------------------------------------------------------

                  if( ichl .eq. 0             .or.
     &                chaf(k+ichl+1) .eq. '{' ) goto 999


                        l = l + 1
                        chag( l ) = yen

                     do 300 ll = 1, ichl

                        l = l + 1
                        chag( l ) = isex(ichl,ichn)(ll:ll)

  300                continue


                        k = k + ichl + 1

                        l = l + 1
                        chag( l ) = '{'

                     if( chaf(k) .eq. ' ' ) k = k + 1

                     if( chaf(k) .eq. '}' ) then

                        l = l + 1
                        chag( l ) = ' '
                        l = l + 1
                        chag( l ) = '}'
                        l = l + 1
                        chag( l ) = chaf(k)

                     else if( chaf(k) .eq. yen ) then

                        l = l + 1
                        chag( l ) = chaf(k)

  310                   k = k + 1

                        if( k .le. icn .and.
     &                      chaf(k) .ne. ' ' .and.
     &                      chaf(k) .ne. yen .and.
     &                      chaf(k) .ne. '}' ) then

                              l = l + 1
                              chag( l ) = chaf(k)

                              goto 310

                        end if

                           l = l + 1
                           chag( l ) = '}'

                           k = k - 1

                     else if( chaf(k) .ne. '{' ) then

                        l = l + 1
                        chag( l ) = chaf(k)
                        l = l + 1
                        chag( l ) = '}'

                     end if


                  if( ( ichl .eq. 4 .and. ichn .eq. 3 ) .or.
     &                ( ichl .eq. 5 .and. ichn .eq. 6 ) .or.
     &                ( ichl .eq. 8 .and. ichn .eq. 2 ) ) then

                        k = k + 1

                     if( chaf(k) .eq. ' ' ) k = k + 1

                        l = l + 1
                        chag( l ) = '{'

                     if( chaf(k) .eq. '}' ) then

                        l = l + 1
                        chag( l ) = ' '
                        l = l + 1
                        chag( l ) = '}'
                        l = l + 1
                        chag( l ) = chaf(k)

                     else if( chaf(k) .eq. yen ) then

                        l = l + 1
                        chag( l ) = chaf(k)

  320                   k = k + 1

                        if( k .le. icn .and.
     &                      chaf(k) .ne. ' ' .and.
     &                      chaf(k) .ne. yen .and.
     &                      chaf(k) .ne. '}' ) then

                              l = l + 1
                              chag( l ) = chaf(k)

                              goto 320

                        end if

                           l = l + 1
                           chag( l ) = '}'

                           k = k - 1

                     else if( chaf(k) .ne. '{' ) then

                        l = l + 1
                        chag( l ) = chaf(k)
                        l = l + 1
                        chag( l ) = '}'

                     end if


                  end if

*-----------------------------------------------------------------------

                     preqc = .true.
                     return

*-----------------------------------------------------------------------

         end if

  999                preqc = .false.

      return
      end

************************************************************************
*                                                                      *
      function prepc(k,l,chaf,chag,icn,ie,ib,ipsp)
*                                                                      *
*           INSERT '{' AND '}' FOR SPECIAL CHARACTER AFTER '\'         *
*           LISTED ABOVE WITHOUT '{' JUST BEFORE THE '\'               *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      logical prepc

      parameter ( icmax = 13, iscn = 13 )
      parameter ( ic02 = 12 )
      parameter ( ic04 = 5 )
      parameter ( ic05 = 5 )
      parameter ( ic07 = 4 )
      parameter ( ic09 = 12 )
      parameter ( ic010 = 1 )
      parameter ( ic012 = 1 )
      parameter ( ic013 = 4 )

      parameter ( ic01  = 0 )
      parameter ( ic03  = 0 )
      parameter ( ic06  = 0 )
      parameter ( ic08  = 0 )
      parameter ( ic011 = 0 )

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(ichrl)*1
      dimension ie(0:inig)

*-----------------------------------------------------------------------

      character isex(icmax,iscn)*20
      dimension nic(icmax)

      character dum20*20

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic02),i=2,2)
     &         /'sf', 'hv', 'rm', 'tm', 'tt', 'cr', 'it', 'bf', 'ib',
     &          'sb', 'mn', 'gt'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic04),i=4,4)
     &         /'tiny','Tiny','TINY','huge','Huge'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic05),i=5,5)
     &         /'small','large','Large','LARGE','color'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic07),i=7,7)
     &         /'ovalbox','Ovalbox','ovalBox','OvalBox'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic09),i=9,9)
     &         /'singlebox','Singlebox','singleBox','SingleBox',
     &          'doublebox','Doublebox','doubleBox','DoubleBox',
     &          'shadowbox','Shadowbox','shadowBox','ShadowBox'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic010),i=10,10)
     &         /'scriptsize'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic012),i=12,12)
     &         /'footnotesize'/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic013),i=13,13)
     &         /'ovalshadowbox','Ovalshadowbox','ovalshadowBox',
     &          'OvalshadowBox'/

*-----------------------------------------------------------------------

      data nic/ ic01, ic02, ic03, ic04, ic05, ic06, ic07, ic08,
     &          ic09, ic010,ic011,ic012,ic013/

*-----------------------------------------------------------------------

      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------
*           INSERT '{' AND '}' FOR SPECIAL CHARACTER AFTER '\'
*           LISTED ABOVE WITHOUT '{' JUST BEFORE THE '\'
*-----------------------------------------------------------------------

         if( chaf(k) .eq. yen ) then

                  ichl = 0
                  ichn = 0

                  jfin = min( icn, k+icmax )

                  dum20 = chaf(k+1)

*-----------------------------------------------------------------------

               do 102 j = k + 2, jfin

                     jj = j - k

                     dum20 = dum20(1:jj-1)//chaf(j)

                  if( nic(jj) .gt. 0 ) then

                     do 201 i = 1, nic(jj)

                        if( dum20(1:jj) .eq. isex(jj,i) ) then

                           ichl = jj
                           ichn = i

                        end if

  201                continue

                  end if


  102          continue

*-----------------------------------------------------------------------

               if( ichl .gt. 0 ) then

                  if( ib .gt. 0 .and. ipsp .ne. 1 ) then

                        ib = ib + 1

                        ie(ib) = ib

                        l = l + 1
                        chag( l ) = '{'

                        l = l + 1
                        chag( l ) = chaf(k)

                        prepc = .true.
                        return

                  end if

               end if

*-----------------------------------------------------------------------

         end if

                        prepc = .false.

      return
      end

************************************************************************
*                                                                      *
      subroutine prtex(chaf,icn,ifon)
*                                                                      *
*              PRE-PROCEDURE FOR TEX MODE                              *
*                                                                      *
*              ADD {\it        }                                       *
*                                                                      *
*              +      -> \+                                            *
*              -      -> \-                                            *
*              |      -> \|                                            *
*              =      -> \eq                                           *
*              number -> \number                                       *
*              [, ]   -> \left[, \right]                               *
*              \(, \) -> \left(, \right)                               *
*              <, >   -> \left<, \right>                               *
*              \{, \} -> \left\{, \right\}                             *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(ichrl)*1

      character tex01*5

      logical prtex1, numec
      logical spjpn

*-----------------------------------------------------------------------

      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

            tex01 = '{'//yen//'it '

*-----------------------------------------------------------------------
*     ADD '{' AND '}' FOR ONE CHARACTER AFTER '^' OR '_'
*-----------------------------------------------------------------------

               i   = 0
               j   = 0

  101       i = i + 1

*-----------------------------------------------------------------------

            if( spjpn(i,icn,icf,chaf,ifon,njpn) ) then

*              SKIP JAPANESE TWO BITE CHARACTERS

*-----------------------------------------------------------------------

                  do 500 k = i, icf

                        j = j + 1
                        chag( j ) = chaf( k )

  500             continue

                        i = icf

*-----------------------------------------------------------------------

            else if( chaf(i) .eq. yen .and.
     &             ( chaf(i+1) .eq. '^' .or. chaf(i) .eq. '_' ) ) then

                     do m = 0, 1

                        j = j + 1
                        chag( j ) = chaf( i + m )

                     end do

                        i = i + 1

            else if( ( chaf(i) .eq. '^' .or. chaf(i) .eq. '_' ) .and.
     &                 chaf(i+1) .ne. '{' .and. i .lt. icn ) then

                              j = j + 1
                        chag( j ) = chaf(i)
                              j = j + 1
                        chag( j ) = '{'
                              j = j + 1
                        chag( j ) = chaf( i + 1 )
                              j = j + 1
                        chag( j ) = '}'

                        i = i + 1

            else

                        j = j + 1
                        chag( j ) = chaf( i )

            end if


            if( i .lt. icn ) goto 101


            do 201 k = 1, j

                  chaf( k ) = chag( k )

  201       continue

                  icn = j

*-----------------------------------------------------------------------
*     PREPROCEDURE FOR TEX MODE
*-----------------------------------------------------------------------

            do 10 k = 1, 5

               chag(k) = tex01(k:k)

   10       continue

               i   = 0
               j   = 5

  100       i = i + 1

*-----------------------------------------------------------------------

            if( spjpn(i,icn,icf,chaf,ifon,njpn) ) then

*              SKIP JAPANESE TWO BITE CHARACTERS

*-----------------------------------------------------------------------

                  do 510 k = i, icf

                        j = j + 1
                        chag( j ) = chaf( k )

  510             continue

                        i = icf

*-----------------------------------------------------------------------

            else if( chaf(i)   .eq. yen .and.
     &             ( chaf(i+1) .eq. 'h' .or. chaf(i+1) .eq. 'v' ) .and.
     &               chaf(i+2) .eq. 's' .and.
     &               chaf(i+3) .eq. 'p' .and.
     &               chaf(i+4) .eq. 'a' .and.
     &               chaf(i+5) .eq. 'c' .and.
     &               chaf(i+6) .eq. 'e' .and.
     &               chaf(i+7) .eq. '{' ) then

                  do 511 k = i + 8, icn

                     if( chaf(k) .eq. '}' ) goto 512

  511             continue

  512                icf = k

                  do 513 k = i, icf

                        j = j + 1
                        chag( j ) = chaf( k )

  513             continue

                        i = icf

            else if( prtex1(i,j,chaf,chag,icn) ) then


            else if( numec(chaf(i)) ) then

                  chag(j+1) = yen
                  chag(j+2) = chaf(i)

                  j = j + 2

            else if( chaf(i) .eq. '+' .or.
     &               chaf(i) .eq. '-' .or.
     &               chaf(i) .eq. '|' ) then

                  chag(j+1) = yen
                  chag(j+2) = chaf(i)

                  j = j + 2

            else if( chaf(i) .eq. '=' ) then

                  chag(j+1) = yen
                  chag(j+2) = 'e'
                  chag(j+3) = 'q'

                  j = j + 3

               if( chaf(i+1) .eq. ' ' ) then

                  j = j + 1
                  chag(j) = ' '

               end if

            else if( chaf(i) .eq. '[' .or.
     &               chaf(i) .eq. '<' ) then

                  chag(j+1) = yen
                  chag(j+2) = 'l'
                  chag(j+3) = 'e'
                  chag(j+4) = 'f'
                  chag(j+5) = 't'
                  chag(j+6) = chaf(i)

                  j = j + 6

            else if( chaf(i) .eq. ']' .or.
     &               chaf(i) .eq. '>' ) then

                  chag(j+1) = yen
                  chag(j+2) = 'r'
                  chag(j+3) = 'i'
                  chag(j+4) = 'g'
                  chag(j+5) = 'h'
                  chag(j+6) = 't'
                  chag(j+7) = chaf(i)

                  j = j + 7

            else if( chaf(i)   .eq. yen .and.
     &               chaf(i+1) .eq. '(' ) then

                  chag(j+1) = yen
                  chag(j+2) = 'l'
                  chag(j+3) = 'e'
                  chag(j+4) = 'f'
                  chag(j+5) = 't'
                  chag(j+6) = chaf(i+1)

                  j = j + 6
                  i = i + 1

            else if( chaf(i)   .eq. yen .and.
     &               chaf(i+1) .eq. ')' ) then

                  chag(j+1) = yen
                  chag(j+2) = 'r'
                  chag(j+3) = 'i'
                  chag(j+4) = 'g'
                  chag(j+5) = 'h'
                  chag(j+6) = 't'
                  chag(j+7) = chaf(i+1)

                  j = j + 7
                  i = i + 1

            else if( chaf(i)   .eq. yen .and.
     &               chaf(i+1) .eq. '{' ) then

                  chag(j+1) = yen
                  chag(j+2) = 'l'
                  chag(j+3) = 'e'
                  chag(j+4) = 'f'
                  chag(j+5) = 't'
                  chag(j+6) = yen
                  chag(j+7) = chaf(i+1)

                  j = j + 7
                  i = i + 1

            else if( chaf(i)   .eq. yen .and.
     &               chaf(i+1) .eq. '}' ) then

                  chag(j+1) = yen
                  chag(j+2) = 'r'
                  chag(j+3) = 'i'
                  chag(j+4) = 'g'
                  chag(j+5) = 'h'
                  chag(j+6) = 't'
                  chag(j+7) = yen
                  chag(j+8) = chaf(i+1)

                  j = j + 8
                  i = i + 1

            else

                        j = j + 1
                        chag( j ) = chaf( i )

            end if


            if( i .lt. icn ) goto 100

*-----------------------------------------------------------------------

                  chag(j+1) = '}'

                  j = j + 1

            do 200 k = 1, j

                  chaf( k ) = chag( k )

  200       continue

                  icn = j


      return
      end

************************************************************************
*                                                                      *
      function prtex1(i,j,chaf,chag,icn)
*                                                                      *
*              \+CHARACHTER FOR SYMBOL                                 *
*                                                                      *
*              \+                ->  \+                                *
*              \-                ->  \-                                *
*              \|                ->  \|                                *
*              \>                ->  \>                                *
*              \eq               ->  \eq                               *
*              \equal            ->  \equal                            *
*              \123              ->  \123                              *
*                                                                      *
*              \number           ->  number                            *
*                                                                      *
*              \left[, \right]   ->  \left[, \right]                   *
*              \left(, \right)   ->  \left(, \right)                   *
*              \left<, \right>   ->  \left<, \right>                   *
*              \left\{, \right\} ->  \left\{, \right\}                 *
*              \Big,big[(\{,])\} ->  \Big,big[(\{,])\}                 *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(ichrl)*1

      logical prtex1, numec

*-----------------------------------------------------------------------

      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

         if( chaf(i) .eq. yen ) then

               prtex1 = .true.

            if( chaf(i+1) .eq. '+' .or.
     &          chaf(i+1) .eq. '-' .or.
     &          chaf(i+1) .eq. '|' .or.
     &          chaf(i+1) .eq. '>' ) then

               chag(j+1) = chaf(i)
               chag(j+2) = chaf(i+1)

               j = j + 2
               i = i + 1

            else if( chaf(i+1) .eq. 'e' .and. chaf(i+2) .eq. 'q' ) then

               chag(j+1) = chaf(i)
               chag(j+2) = chaf(i+1)
               chag(j+3) = chaf(i+2)

               j = j + 3
               i = i + 2

            else if( numec(chaf(i+1)) .and.
     &               numec(chaf(i+2)) .and.
     &               numec(chaf(i+3)) ) then

               chag(j+1) = chaf(i)
               chag(j+2) = chaf(i+1)
               chag(j+3) = chaf(i+2)
               chag(j+4) = chaf(i+3)

               j = j + 4
               i = i + 3

            else if( numec(chaf(i+1)) .and.
     &               numec(chaf(i+2)) ) then

               chag(j+1) = chaf(i)
               chag(j+2) = chaf(i+1)
               chag(j+3) = chaf(i+2)

               j = j + 3
               i = i + 2

            else if( numec(chaf(i+1)) ) then

               chag(j+1) = chaf(i+1)

               j = j + 1
               i = i + 1

            else if( ( chaf(i+1) .eq. 'B' .or.
     &                 chaf(i+1) .eq. 'b' ) .and.
     &               ( chaf(i+4) .eq. '[' .or.
     &                 chaf(i+4) .eq. ']' .or.
     &                 chaf(i+4) .eq. '(' .or.
     &                 chaf(i+4) .eq. ')' ) ) then

               chag(j+1) = chaf(i)
               chag(j+2) = chaf(i+1)
               chag(j+3) = chaf(i+2)
               chag(j+4) = chaf(i+3)
               chag(j+5) = chaf(i+4)

               j = j + 5
               i = i + 4

            else if( ( chaf(i+1) .eq. 'B' .or.
     &                 chaf(i+1) .eq. 'b' ) .and.
     &               ( chaf(i+4) .eq. yen ) .and.
     &               ( chaf(i+5) .eq. '{' .or.
     &                 chaf(i+5) .eq. '}' ) ) then

               chag(j+1) = chaf(i)
               chag(j+2) = chaf(i+1)
               chag(j+3) = chaf(i+2)
               chag(j+4) = chaf(i+3)
               chag(j+5) = chaf(i+4)
               chag(j+6) = chaf(i+5)

               j = j + 6
               i = i + 5

            else if( chaf(i+1) .eq. 'e' .and.
     &               chaf(i+2) .eq. 'q' .and.
     &               chaf(i+3) .eq. 'u' .and.
     &               chaf(i+4) .eq. 'a' .and.
     &               chaf(i+5) .eq. 'l' ) then

               chag(j+1) = chaf(i)
               chag(j+2) = chaf(i+1)
               chag(j+3) = chaf(i+2)
               chag(j+4) = chaf(i+3)
               chag(j+5) = chaf(i+4)
               chag(j+6) = chaf(i+5)

               j = j + 6
               i = i + 5

            else if(   chaf(i+1) .eq. 'l' .and.
     &                 chaf(i+2) .eq. 'e' .and.
     &               ( chaf(i+5) .eq. '[' .or.
     &                 chaf(i+5) .eq. '(' .or.
     &                 chaf(i+5) .eq. '<' ) ) then

               chag(j+1) = chaf(i)
               chag(j+2) = chaf(i+1)
               chag(j+3) = chaf(i+2)
               chag(j+4) = chaf(i+3)
               chag(j+5) = chaf(i+4)
               chag(j+6) = chaf(i+5)

               j = j + 6
               i = i + 5

            else if(   chaf(i+1) .eq. 'r' .and.
     &                 chaf(i+2) .eq. 'i' .and.
     &               ( chaf(i+6) .eq. ']' .or.
     &                 chaf(i+6) .eq. ')' .or.
     &                 chaf(i+6) .eq. '>' ) ) then

               chag(j+1) = chaf(i)
               chag(j+2) = chaf(i+1)
               chag(j+3) = chaf(i+2)
               chag(j+4) = chaf(i+3)
               chag(j+5) = chaf(i+4)
               chag(j+6) = chaf(i+5)
               chag(j+7) = chaf(i+6)

               j = j + 7
               i = i + 6

            else if(   chaf(i+1) .eq. 'l' .and.
     &                 chaf(i+2) .eq. 'e' .and.
     &                 chaf(i+5) .eq. yen .and.
     &                 chaf(i+6) .eq. '{' ) then

               chag(j+1) = chaf(i)
               chag(j+2) = chaf(i+1)
               chag(j+3) = chaf(i+2)
               chag(j+4) = chaf(i+3)
               chag(j+5) = chaf(i+4)
               chag(j+6) = chaf(i+5)
               chag(j+7) = chaf(i+6)

               j = j + 7
               i = i + 6

            else if(   chaf(i+1) .eq. 'r' .and.
     &                 chaf(i+2) .eq. 'i' .and.
     &                 chaf(i+6) .eq. yen .and.
     &                 chaf(i+7) .eq. '}' ) then

               chag(j+1) = chaf(i)
               chag(j+2) = chaf(i+1)
               chag(j+3) = chaf(i+2)
               chag(j+4) = chaf(i+3)
               chag(j+5) = chaf(i+4)
               chag(j+6) = chaf(i+5)
               chag(j+7) = chaf(i+6)
               chag(j+8) = chaf(i+7)

               j = j + 8
               i = i + 7

            else

               prtex1 = .false.

            end if

         else

               prtex1 = .false.

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine psufx(chaf,icn,isuf,isup,isdn,ifon)
*                                                                      *
*              PRE-PROCEDURE FOR SUFFIXES                              *
*                                                                      *
*             {^{  }_{  }}                                             *
*                                                                      *
*             |__ ISUF(K,3) = END IF FIRST }                           *
*             |__ ISUF(K,2) = END OF }                                 *
*             |__ ISUF(K,1) = 1, 2, 3, 4                               *
*                                                                      *
*             ISUF(K,1) = 1 A^B or A_B                                 *
*             ISUF(K,1) = 2 A^B_C or A_B^C                             *
*             ISUF(K,1) = 3 ^AB or _AB                                 *
*             ISUF(K,1) = 4 ^A_BC or _A^BC                             *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(ichrl)*1

      dimension isu(0:500), isd(0:500), ib(0:500)

      logical psuf1,psuf2,psuf3,psuf4,psuf5,psuf6
      logical spjpn

      dimension it(0:500,2)
      dimension isuf(0:ichrl,3)

*-----------------------------------------------------------------------

      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

               do 10 i = 1, 3
               do 10 j = 0, ichrl
                  isuf(j,i) = 0
   10          continue

               do 20 i = 1, ichrl

                  chag(i) = ' '

   20          continue

*-----------------------------------------------------------------------

               isup = 0
               isdn = 0

               isst = 1

*-----------------------------------------------------------------------

               i   = 0
               j   = 0
               is  = 0
               isn = 0

               ib(is) = 0

  100       continue

            ipsp = 0
            if( i .gt. 0 ) then
               if( chaf(i) .eq. '^' .or. chaf(i) .eq. '_' ) ipsp = 1
            end if

            i = i + 1

*-----------------------------------------------------------------------

            if( chaf(i) .eq. yen .and. i .lt. icn .and.
     &        ( chaf(i+1) .eq. '^' .or. chaf(i+1) .eq. '_' .or.
     &          chaf(i+1) .eq. '{' .or. chaf(i+1) .eq. '}' )  ) then

                        j = j + 1
                        chag( j ) = chaf( i )
                        j = j + 1
                        chag( j ) = chaf( i + 1 )
                        i = i + 1

*-----------------------------------------------------------------------

            else if( psuf5(i,j,chaf,chag,icn,is,isn,isu,isd,ib,
     &                     it,isuf,isup,isdn,isst) ) then

*              END OF ONE SUFFIX BLOCK : IS

*-----------------------------------------------------------------------

            else if( psuf4(i,j,chaf,chag,icn,is,isn,
     &                     ib,ifon,ipsp) ) then

*              ADD '{' AND '}' FOR ONE CHARACTER AFTER '^' OR '_'

*-----------------------------------------------------------------------

            else if( spjpn(i,icn,icf,chaf,ifon,njpn) ) then

*              SKIP JAPANESE TWO BITE CHARACTERS

*-----------------------------------------------------------------------

                  do 500 k = i, icf

                        j = j + 1

                        chag( j ) = chaf( k )

  500             continue

                        i = icf

*-----------------------------------------------------------------------

            else if( psuf1(i,j,chaf,chag,icn,is,isn,isu,isd,ib,
     &                     it,isuf,isup,isdn,isst) ) then

*              FIRST SUFIX SIGNAL : ISN = 0, '^' OR '_'

*-----------------------------------------------------------------------

            else if( psuf2(i,j,chaf,chag,icn,is,isn,isu,isd,ib,
     &                     it,isuf,isup,isdn,isst) ) then

*              SECOND SUFIX SIGNAL : ISN = 1, '^' OR '_'

*-----------------------------------------------------------------------

            else if( psuf3(i,j,chaf,chag,icn,is,isn,isu,isd,ib,
     &                     it,isuf,isup,isdn,isst) ) then

*              DOUBLE SUFFIXES : ENTER ONE BLANK

*-----------------------------------------------------------------------

            else if( psuf6(i,j,chaf,chag,icn,is,isn,ib) ) then

*              '{' AND '}'

*-----------------------------------------------------------------------

            else

                        j = j + 1
                        chag( j ) = chaf( i )

            end if


            if( i .lt. icn ) goto 100

*-----------------------------------------------------------------------

                  is = is + 1

            do 220 k = 1, 500

                  is = is - 1

               if( is .ge. 0 ) then

                     ib(is) = ib(is) + 1

                     do 230 l = 1, 500

                        ib(is) = ib(is) - 1

                        if( ib(is) .gt. 0 ) then

                                    j = j + 1
                              chag( j ) = '}'

                        else

                           goto 240

                        end if

  230                continue
  240                continue

                  if( is .gt. 0 ) then

                           j = j + 1
                     chag( j ) = '}'

                     isuf( it(is,1), 2 ) = j
                     isst = it(is,2)

                  end if

               else

                     goto 210

               end if


  220       continue
  210       continue

*-----------------------------------------------------------------------

            do 200 k = 1, j

                  chaf( k ) = chag( k )

  200       continue

                  icn = j


*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function psuf1(i,j,chaf,chag,icn,is,isn,isu,isd,ib,
     &               it,isuf,isup,isdn,isst)
*                                                                      *
*              FIRST SUFIX SIGNAL : ISN = 0, '^' OR '_'                *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(ichrl)*1

      dimension isu(0:500), isd(0:500), ib(0:500)

      dimension it(0:500,2)
      dimension isuf(0:ichrl,3)

      logical psuf1

*-----------------------------------------------------------------------

      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

                  ipsq = 0

            if( chaf(i) .eq. ' ' .and.
     &        ( chaf(i+1) .eq. '^' .or. chaf(i+1) .eq. '_' ) ) then

                  j = j + 1
                  chag( j ) = ' '

                  ipsq = 1
                  i = i + 1

            end if

*-----------------------------------------------------------------------

            if( isn .eq. 0 .and. i .lt. icn .and.
     &        ( chaf(i) .eq. '^' .or. chaf(i) .eq. '_' )  ) then

                              j = j + 1
                        chag( j ) = '{'

                        is = is + 1

                        it( is, 1 ) = j
                        it( is, 2 ) = isst

                        if( isst .eq. 1 ) then

                           if( chaf(i) .eq. '^' ) then

                              isup = 1

                           else if( chaf(i) .eq. '_' ) then

                              isdn = 1

                           end if

                        end if

                        isst = isst + 1

                     if( i .eq. 1 .or. ipsq .eq. 1 ) then

                        isuf( j, 1 ) = 3

                     else

                        isuf( j, 1 ) = 1

                     end if


                     if( chaf(i) .eq. '^' ) then

                        isu(is) = 1
                        isd(is) = 0

                     else

                        isu(is) = 0
                        isd(is) = 1

                     end if

                        ib(is) = 0

                              j = j + 1
                        chag( j ) = chaf( i )

*-----------------------------------------------------------------------

               if(  chaf(i+1) .eq. '^' .or. chaf(i+1) .eq. '_' ) then

                              j = j + 1
                        chag( j ) = '{'
                              j = j + 1
                        chag( j ) = yen
                              j = j + 1
                        chag( j ) = chaf( i + 1 )
                              j = j + 1
                        chag( j ) = '}'
                              i = i + 1

                        isn = 1

               end if

*-----------------------------------------------------------------------

               psuf1 = .true.

*-----------------------------------------------------------------------

            else

               psuf1 = .false.

            end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function psuf2(i,j,chaf,chag,icn,is,isn,isu,isd,ib,
     &               it,isuf,isup,isdn,isst)
*                                                                      *
*              SECOND SUFIX SIGNAL : ISN = 1, '^' OR '_'               *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(ichrl)*1

      dimension isu(0:500), isd(0:500), ib(0:500)

      dimension it(0:500,2)
      dimension isuf(0:ichrl,3)

      logical psuf2

*-----------------------------------------------------------------------

      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

         if( ( isn .eq. 1 ) .and.
     &      (   ( chaf(i) .eq. '^' .or. chaf(i) .eq. '_' ) .and.
     &          ( i .lt. icn ) .and.
     &        ( ( chaf(i) .eq. '^' .and. isu(is) .eq. 0 ) .or.
     &          ( chaf(i) .eq. '_' .and. isd(is) .eq. 0 ) ) ) ) then

                     if( chaf(i) .eq. '^' ) then

                        isu(is) = isu(is) + 1

                     else

                        isd(is) = isd(is) + 1

                     end if

                        ib(is) = 0

                        isn = 0

                              j = j + 1

                        chag(j) = chaf(i)


                        isuf( it(is,1), 1 ) = isuf( it(is,1), 1 ) + 1
                        isuf( it(is,1), 3 ) = j - 1


                        if( isst .eq. 2 ) then

                           if( chaf(i) .eq. '^' ) then

                              isup = 1

                           else if( chaf(i) .eq. '_' ) then

                              isdn = 1

                           end if

                        end if

               psuf2 = .true.

*-----------------------------------------------------------------------

            else if( ( isn .eq. 1 ) .and.
     &      (   ( chaf(i) .eq. '^' .or. chaf(i) .eq. '_' ) .and.
     &          ( i .eq. icn ) .and.
     &        ( ( chaf(i) .eq. '^' .and. isu(is) .eq. 0 ) .or.
     &          ( chaf(i) .eq. '_' .and. isd(is) .eq. 0 ) ) ) ) then

                              j = j + 1
                        chag( j ) = yen
                              j = j + 1
                        chag( j ) = chaf( i )

*-----------------------------------------------------------------------

               psuf2 = .true.

         else

               psuf2 = .false.

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function psuf3(i,j,chaf,chag,icn,is,isn,isu,isd,ib,
     &               it,isuf,isup,isdn,isst)
*                                                                      *
*              DOUBLE SUFFIXES : ENTER ONE BLANK                       *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(ichrl)*1

      dimension isu(0:500), isd(0:500), ib(0:500)

      dimension it(0:500,2)
      dimension isuf(0:ichrl,3)

      logical psuf3

*-----------------------------------------------------------------------

      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

         if( ( isn .eq. 1 ) .and.
     &      (   ( chaf(i) .eq. '^' .or. chaf(i) .eq. '_' ) .and.
     &          ( i .lt. icn ) .and. ( chaf(i+1) .ne. '}' ) .and.
     &        ( ( chaf(i) .eq. '^' .and. isu(is) .ne. 0 ) .or.
     &          ( chaf(i) .eq. '_' .and. isd(is) .ne. 0 ) ) ) ) then

                        j = j + 1
                  chag( j ) = '}'

                  isuf( it(is,1), 2 ) = j
                  isst = it(is,2)

                        j = j + 1
                  chag( j ) = ' '
                        j = j + 1
                  chag( j ) = '{'


                  it( is, 1 )  = j
                  isuf( j, 1 ) = 1


                     if( chaf(i) .eq. '^' ) then

                        isu(is) = 1
                        isd(is) = 0

                     else

                        isu(is) = 0
                        isd(is) = 1

                     end if

                        ib(is) = 0

                        j = j + 1
                  chag( j ) = chaf( i )

                        isn = 0

                     psuf3 = .true.

*-----------------------------------------------------------------------

         else if( ( isn .eq. 1 ) .and.
     &       ( chaf(i) .eq. '^' .or. chaf(i) .eq. '_' ) .and.
     &       ( i .eq. icn .or. chaf(i+1) .eq. '}' ) ) then

                        j = j + 1
                  chag( j ) = '}'

                  isuf( it(is,1), 2 ) = j
                  isst = it(is,2)

  100             is = is - 1

               if( ib(is) .eq. 0 .and. is .gt. 0 ) then

                       j = j + 1
                  chag(j) = '}'

                  isuf( it(is,1), 2 ) = j
                  isst = it(is,2)

                  goto 100

               end if

                        j = j + 1
                  chag( j ) = yen
                        j = j + 1
                  chag( j ) = chaf( i )

                     psuf3 = .true.

*-----------------------------------------------------------------------

         else

               psuf3 = .false.

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function psuf4(i,j,chaf,chag,icn,is,isn,ib,ifon,ipsp)
*                                                                      *
*              ADD '{' AND '}' FOR ONE CHARACTER AFTER '^' OR '_'      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(ichrl)*1

      dimension ib(0:500)

      logical psuf4

*-----------------------------------------------------------------------

      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

            if( is .ne. 0 .and. ib(is) .eq. 0 .and. isn .eq. 0 .and.
     &          ipsp .eq. 1 .and. chaf(i) .ne. '{' ) then

*-----------------------------------------------------------------------

                  if( ifon .gt. 0 .and.
     &                icn .ge. 7 .and. i .le. icn - 6 .and.
     &                chaf(i)   .eq. yen .and.
     &                chaf(i+1) .eq. '3' .and.
     &                chaf(i+2) .eq. '7' .and.
     &                chaf(i+3) .eq. '7' .and.
     &                chaf(i+4) .eq. yen .and.
     &                chaf(i+5) .eq. '0' .and.
     &                chaf(i+6) .eq. '0' .and.
     &                chaf(i+7) .eq. '1' ) then

                              j = j + 1
                        chag( j ) = '{'

                              j = j + 1
                        chag( j ) = chaf( i )

                              i = i + 1
                              j = j + 1
                        chag( j ) = chaf( i )
                              i = i + 1
                              j = j + 1
                        chag( j ) = chaf( i )
                              i = i + 1
                              j = j + 1
                        chag( j ) = chaf( i )
                              i = i + 1
                              j = j + 1
                        chag( j ) = chaf( i )
                              i = i + 1
                              j = j + 1
                        chag( j ) = chaf( i )
                              i = i + 1
                              j = j + 1
                        chag( j ) = chaf( i )
                              i = i + 1
                              j = j + 1
                        chag( j ) = chaf( i )

                              i = i + 1
                              j = j + 1
                        chag( j ) = chaf( i )
                              i = i + 1
                              j = j + 1
                        chag( j ) = chaf( i )

                              j = j + 1
                        chag( j ) = yen
                              j = j + 1
                        chag( j ) = '3'
                              j = j + 1
                        chag( j ) = '7'
                              j = j + 1
                        chag( j ) = '7'
                              j = j + 1
                        chag( j ) = yen
                              j = j + 1
                        chag( j ) = '0'
                              j = j + 1
                        chag( j ) = '0'
                              j = j + 1
                        chag( j ) = '0'
                              j = j + 1
                        chag( j ) = '}'


                     if( chaf(i+1) .eq. yen .and.
     &                   chaf(i+2) .eq. '3' .and.
     &                   chaf(i+3) .eq. '7' .and.
     &                   chaf(i+4) .eq. '7' .and.
     &                   chaf(i+5) .eq. yen .and.
     &                   chaf(i+6) .eq. '0' .and.
     &                   chaf(i+7) .eq. '0' .and.
     &                   chaf(i+8) .eq. '0' ) then

                        i = i + 8

                     else if( i + 1 .lt. icn ) then

                        i = i - 8

                        chaf(i+1) = yen
                        chaf(i+2) = '3'
                        chaf(i+3) = '7'
                        chaf(i+4) = '7'
                        chaf(i+5) = yen
                        chaf(i+6) = '0'
                        chaf(i+7) = '0'
                        chaf(i+8) = '1'

                     end if

*-----------------------------------------------------------------------

                  else

                              j = j + 1
                        chag( j ) = '{'
                              j = j + 1
                        chag( j ) = chaf( i )
                              j = j + 1
                        chag( j ) = '}'

                  end if

                        isn = 1

*-----------------------------------------------------------------------

               psuf4 = .true.

            else

               psuf4 = .false.

            end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function psuf5(i,j,chaf,chag,icn,is,isn,isu,isd,ib,
     &               it,isuf,isup,isdn,isst)
*                                                                      *
*              END OF ONE SUFFIX BLOCK : IS                            *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(ichrl)*1

      dimension isu(0:500), isd(0:500), ib(0:500)

      dimension it(0:500,2)
      dimension isuf(0:ichrl,3)

      logical psuf5

*-----------------------------------------------------------------------

      character yen*1

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------

            yen = char(92)

*-----------------------------------------------------------------------

         if( isn .eq. 1 .and.
     &       chaf(i) .ne. '^' .and.
     &       chaf(i) .ne. '_' ) then

                        j = j + 1
                  chag( j ) = '}'

                  isuf( it(is,1), 2 ) = j
                  isst = it(is,2)


  100             is = is - 1

               if( ib(is) .eq. 0 .and. is .gt. 0 ) then

                        j = j + 1
                  chag( j ) = '}'

                  isuf( it(is,1), 2 ) = j
                  isst = it(is,2)

                  goto 100

               end if

                        j = j + 1
                  chag( j ) = chaf( i )

                  isn = 0

*-----------------------------------------------------------------------
*CAUTION

               if( chaf(i) .eq. '}' ) then

                     ib(is) = ib(is) - 1

                     if( ib(is) .lt. 0 ) ib(is) = 0

                     if( ib(is) .eq. 0 .and. is .gt. 0 )  isn = 1

               end if

*-----------------------------------------------------------------------

                  psuf5 = .true.

         else

               psuf5 = .false.

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function psuf6(i,j,chaf,chag,icn,is,isn,ib)
*                                                                      *
*              '{' AND '}'                                             *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(ichrl)*1

      dimension ib(0:500)

      logical psuf6

*-----------------------------------------------------------------------

      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

            if( i .gt. 1 .and. chaf(i) .eq. '{' ) then

                     ib(is) = ib(is) + 1

                              j = j + 1

                        chag(j) = chaf(i)

                        psuf6 = .true.

            else if( chaf(i) .eq. '}' ) then

                             j = j + 1
                        chag(j) = chaf(i)

                     ib(is) = ib(is) - 1

                     if( ib(is) .lt. 0 ) ib(is) = 0

                     if( ib(is) .eq. 0 .and. is .gt. 0 )  isn = 1

                        psuf6 = .true.

*-----------------------------------------------------------------------

            else

               psuf6 = .false.

            end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function alphb(dum)
*                                                                      *
*              DUM IS ALPHABET OR NOT                                  *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      character dum*1

      logical alphb

*-----------------------------------------------------------------------

      alphb  =  dum.eq.'A'.or.dum.eq.'B'.or.dum.eq.'C'.or.
     &          dum.eq.'D'.or.dum.eq.'E'.or.dum.eq.'F'.or.
     &          dum.eq.'G'.or.dum.eq.'H'.or.dum.eq.'I'.or.
     &          dum.eq.'J'.or.dum.eq.'K'.or.dum.eq.'L'.or.
     &          dum.eq.'M'.or.dum.eq.'N'.or.dum.eq.'O'.or.
     &          dum.eq.'P'.or.dum.eq.'Q'.or.dum.eq.'R'.or.
     &          dum.eq.'S'.or.dum.eq.'T'.or.dum.eq.'U'.or.
     &          dum.eq.'V'.or.dum.eq.'W'.or.dum.eq.'X'.or.
     &          dum.eq.'Y'.or.dum.eq.'Z'.or.
     &          dum.eq.'a'.or.dum.eq.'b'.or.dum.eq.'c'.or.
     &          dum.eq.'d'.or.dum.eq.'e'.or.dum.eq.'f'.or.
     &          dum.eq.'g'.or.dum.eq.'h'.or.dum.eq.'i'.or.
     &          dum.eq.'j'.or.dum.eq.'k'.or.dum.eq.'l'.or.
     &          dum.eq.'m'.or.dum.eq.'n'.or.dum.eq.'o'.or.
     &          dum.eq.'p'.or.dum.eq.'q'.or.dum.eq.'r'.or.
     &          dum.eq.'s'.or.dum.eq.'t'.or.dum.eq.'u'.or.
     &          dum.eq.'v'.or.dum.eq.'w'.or.dum.eq.'x'.or.
     &          dum.eq.'y'.or.dum.eq.'z'


      return
      end


************************************************************************
*                                                                      *
      function numec(dum)
*                                                                      *
*              DUM IS NUMBER OR NOT                                   *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      character dum*1

      logical numec

*-----------------------------------------------------------------------

      numec  =  dum.eq.'1'.or.dum.eq.'2'.or.dum.eq.'3'.or.
     &          dum.eq.'4'.or.dum.eq.'5'.or.dum.eq.'6'.or.
     &          dum.eq.'7'.or.dum.eq.'8'.or.dum.eq.'9'.or.
     &          dum.eq.'0'


      return
      end


************************************************************************
*                                                                      *
      subroutine wtext(jhf,chag,idbg,clal,clmo,ipin,xpo,ypo,
     &                 charo,icno,ixc,iyc,iang,ang,if00,ifon,
     &                 imull,indt,fts0,clc0)
*                                                                      *
*      PURPOSE  : WRITE CHARACTER                                      *
*                                                                      *
*          jhf  : OUTPUT FILE NUMBER                                   *
*         IDBG  : =1 WRITE COMENT IN OUTPUT FILE                       *
*         CLAL  : COLOR ALL                                            *
*         CLMO  : MONO COLOR                                           *
*                                                                      *
*         IPIN  : =0 INITIAL POSITION IS ALREADY DECLARED,             *
*                                                                      *
*                 =1 DECLARES HERE,                                    *
*                 =2 BASE LINE SKIP,                                   *
*                 =3 END OF BASE LINE SKIP,                            *
*                                                                      *
*                 >4 ONE LINE WITH BOX FROM WT:  <= IBOX(IW)           *
*                                                                      *
*                 = -1 FOR TABULAR                                     *
*                                                                      *
*                 = -2 FOR LEGEND BOX                                  *
*                                                                      *
*        XPO    : INITIAL X POSITION                                   *
*        YPO    : INITIAL Y POSITION                                   *
*        ICOLM  : MAX LENGTH OF ONE LINE                               *
*        CHARO  : CHARACTER ARRAY                                      *
*         ICNO  : NUMBER OF CHARACTERS                                 *
*          IXC  : START POSITION SELECTION FOR X 1 - 3                 *
*          IYC  : START POSITION SELECTION FOR Y 0 - 3                 *
*         IANG  : =0 NO ROTATION,  =1 ROTATION                         *
*          ANG  : ANGLE OF ROTAITION                                   *
*         IF00  : DEFAULT FONT                                         *
*         IFON  : JAPANESE CODE                                        *
*        IMULL  : > 0 ; NUMBER OF MULTI LINES                          *
*                 < 0 ; IMULL = - ( ILN *1000 + ICM )                  *
*                       ILN = ROW NUMBER                               *
*                       ICM = COLUMN NUMBER                            *
*                                                                      *
*         INDT  : INDENT CONTROL VALIABLE FOR MULTILINE COMMENTS       *
*                                                                      *
*         FTS0  : INITIAL FONT SIZE                                    *
*         CLC0  : INITIAL FONT COLOR                                   *
*                                                                      *
*       IFD(I)  : ID FOR JAPANESE KANJI INITIALIZATION                 *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chag(inig,0:ichrl)*1

*-----------------------------------------------------------------------

      common /fdfn/ ifd(0:13),ikan1,ikan2,ikan3,ikan4,
     &              idsm(9), idln(6), idxe, idye, ibd1, ibd2, ibd3,
     &              ibd4, ibd5, ibd6, ibd7, ibd8, ibd9, ibd10, ibd11

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      dimension ifn(0:inig), ifs(0:inig), ilt(0:inig), icl(0:inig)
      dimension clp(0:inig,3,3)
      dimension ixp(0:inig), iyp(0:inig), isi(0:inig), iap(0:inig,3)
      dimension ibf(0:inig), ibs(0:inig), igs(0:inig),
     &          ibb(0:inig), ibc(0:inig)
      dimension iwt(0:inig), ibt(0:inig,2), iah(0:inig,6)
      dimension isuf(0:ichrl,3)
      character iacd(0:inig)*3
*
      character charo(ichrl)*1

*-----------------------------------------------------------------------

      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb

*-----------------------------------------------------------------------

      common /wtval1/ strl0, strh0, strb0
      common /wtval2/ clc(0:inig,3), fts(0:inig),
     &                stw(0:inig), sth(0:inig), stb(0:inig),
     &                xps(0:inig), yps(0:inig),
     &                xpp(0:inig), ypp(0:inig),
     &                ixs(0:inig), iys(0:inig), jfn(0:inig)
      common /wtval3/ strla(0:inig), strha(0:inig), strba(0:inig)
      common /wtval4/ stwh(0:inig),  stws(0:inig),
     &                stwu(0:inig),  stwd(0:inig),
     &                strlu(0:inig), strld(0:inig)
      common /wtval5/ rsid(0:inig),  rlng(0:inig),  rsca(0:inig)
      common /wtval6/ stwfu(0:inig), sthfu(0:inig), stbfu(0:inig),
     &                stwfd(0:inig), sthfd(0:inig), stbfd(0:inig),
     &                strhb(0:inig), strhh(0:inig), strbb(0:inig)
      common /wtval7/ strhtp, strbtp, strhbt, strbbt, strhall
      common /wtval8/ xpsi(0:inig), ypsi(0:inig), strl(0:inig),
     &                xpsl(0:inig), xpsr(0:inig), rind(0:inig)

*-----------------------------------------------------------------------

      common /lbox/  bcb(3), bcl(3), bcs(3), bds, bdd, bdc

*-----------------------------------------------------------------------

      common /box/  bxws, bxwl, bxss, bxsl, bxds, bxdl, bxls, bxll
      common /suf/  ftht, ftss, ftks, ftsv,
     &              yplu, ypld, ypku, ypkd, yplh, ypls
      common /con/  cm, dd
      common /pap/  a4w, a4h, wmg, hmg, wct, hct

*-----------------------------------------------------------------------

      dimension clal(3), clc0(3)

      logical chars1, chars2, chars3, chars4, chars5, chars6, chars7
      logical spjpn

*-----------------------------------------------------------------------

      character yen*1

            yen = char(92)

*-----------------------------------------------------------------------
*     SAVE ORIGINAL STRING
*-----------------------------------------------------------------------

            do 112 i = 1, ichrl

                chaf(i) = ' '

  112       continue

            do 111 i = 1, icno

                chaf(i) = charo(i)

  111       continue

               icn = icno

*-----------------------------------------------------------------------
*     PRE PROCEDURE FOR TEX MODE BETWEEN $ AND $
*     AND INSERT '\' BEFORE '(' OR ')'
*-----------------------------------------------------------------------

            call prtex0(chaf,icn,ifon)

*-----------------------------------------------------------------------
*     PRE PROCEDURE FOR SUFIXES
*-----------------------------------------------------------------------

            call psufx(chaf,icn,isuf,isup,isdn,ifon)

*-----------------------------------------------------------------------
*     INITIAL VALUES
*-----------------------------------------------------------------------

               strl0  = 0.0
               strh0  = fts0 * ftht
               strb0  = 0.0

               fts(0) = fts0
               stw(0) = 0.0

               clc(0,1) = clc0(1)
               clc(0,2) = clc0(2)
               clc(0,3) = clc0(3)

            if( ipin .eq.  0 .or.
     &          ipin .gt.  4 .or.
     &          ipin .eq. -1 .or.
     &          ipin .eq. -2 ) then

               xpsin = xpo
               ypsin = ypo

            end if

               xps(0) = 0.0
               yps(0) = 0.0

               xpp(0) = 0.0
               ypp(0) = 0.0

*-----------------------------------------------------------------------

               iftdn    = if00

               ig      = 1
               ilt(ig) = 0
               icl(ig) = 0
               ifn(ig) = iftdn
               ifs(ig) = 0
               ixp(ig) = 1
               iyp(ig) = 0
               isi(ig) = 0

               ids        = 0
               iap(ids,1) = 0
               iap(ids,2) = 0
               iap(ids,3) = 0

               ib      = 0
               ibf(ib) = iftdn

            do 99 i = 0, inig

               ibs(i)   = 0      ! S.Abe 2015/06/18
               igs(i)   = 0
               iwt(i)   = 0
               ibb(i)   = -1
               ibc(i)   = -1
               ibt(i,1) = 0
               ibt(i,2) = 0

               iah(i,1) = 0
               iah(i,2) = 0
               iah(i,3) = 0
               iah(i,4) = 0
               iah(i,5) = 0

               clp(i,1,1) = -r1max
               clp(i,1,2) = 1.0
               clp(i,1,3) = 1.0
               clp(i,2,1) = -r1max
               clp(i,2,2) = 1.0
               clp(i,2,3) = 1.0
               clp(i,3,1) = -r1max
               clp(i,3,2) = 1.0
               clp(i,3,3) = 1.0

   99       continue

               idh        = 0
               idv        = 0

*-----------------------------------------------------------------------
*     INITIALIZATION OF REARRANGE STARTING POSITION
*-----------------------------------------------------------------------

                  idv = 1

                  call write00(8,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                         ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                         idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------

         k = 0

  100 continue

            ipsp = 0
         if( k .gt. 0 ) then
            if( chaf(k) .eq. yen ) ipsp = 1
         end if

         k = k + 1

*-----------------------------------------------------------------------

         if( spjpn(k,icn,icf,chaf,ifon,njpn) ) then

*              SKIP JAPANESE TWO BITE CHARACTERS

*-----------------------------------------------------------------------

               if( ifn(ig) .ne. iftdn ) then

                     if( ilt(ig) .gt. 0 ) then

                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                     else

                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                     end if

                        ilt(ig) = 0
                        ifn(ig) = iftdn

                        ifs(ig) = ipfs
                        icl(ig) = ipcl
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi

               end if

*-----------------------------------------------------------------------

                  do 500 kk = k, icf

                        ilt(ig) = ilt(ig) + 1

                        chag(ig,ilt(ig)) = chaf(kk)

  500             continue

                        k = icf

*-----------------------------------------------------------------------

         else if( chars1(jhf,idbg,k,chaf,icn,chag,ig,
     &                   ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                   ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                   idh,iacd,idv,iah,iwt,ibt,ifon,imull) ) then

*              SUFFIXES AFTER '^', '_'

*-----------------------------------------------------------------------

         else if( chars2(jhf,idbg,k,chaf,icn,chag,ig,
     &                   ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                   ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                   idh,iacd,idv,iah,iwt,ibt,ifon,imull,
     &                   ipsp) ) then

*              CHANGE FONT '{\sf','{\hv','{\rm','{\tm','{\tt','{\cr',
*                          '{\it','{\bf','{\ib','{\sb','{\mn','{\gt'

*-----------------------------------------------------------------------

         else if( chars3(jhf,idbg,k,chaf,icn,chag,ig,
     &                   ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                   ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                   idh,iacd,idv,iah,iwt,ibt,ifon,imull) ) then

*              CHANGE FONT '{\tiny','{\scriptsize','{\footnotesize',
*                                   '{\Tiny',      '{\TINY',
*                          '{\small',
*                          '{\large','{\Large','{\LARGE',
*                          '{\huge', '{\Huge',

*-----------------------------------------------------------------------

         else if( chars4(jhf,idbg,clal,clmo,k,chaf,icn,chag,ig,
     &                   ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                   ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                   idh,iacd,idv,iah,iwt,ibt,ifon,imull) ) then

*             CHANGE COLOR '{\color{r} text}'

*-----------------------------------------------------------------------

         else if( chars5(jhf,idbg,clal,clmo,k,chaf,icn,chag,ig,
     &                   ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                   ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                   idh,iacd,idv,iah,iwt,ibt,ifon,imull) ) then

*             BOX AND COLOR BOX
*                          '{\singlebox{c}{b}{l}     text}'
*                          '{\ovalbox{c}{b}{l}       text}'
*                          '{\doublebox{c}{b}{l}     text}'
*                          '{\shadowbox{c}{b}{s}     text}'
*                          '{\ovalshadowbox{c}{b}{s} text}'

*-----------------------------------------------------------------------

         else if( chars6(jhf,idbg,k,chaf,icn,chag,ig,
     &                   ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                   ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                   idh,iacd,idv,iah,iwt,ibt,ifon,imull) ) then

*              SPECIAL CHARACTER AFTER '\'

*-----------------------------------------------------------------------

         else if( chars7(jhf,idbg,k,chaf,icn,chag,ig,
     &                   ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                   ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                   idh,iacd,idv,iah,iwt,ibt,ifon,imull,
     &                   ipsp) ) then

*              NEGLECT NON-SPECIAL '{' AND '}'
*              OR END OF CHANGE FONT, COLOR, FONTSIZE ETC.

*-----------------------------------------------------------------------

         else if( chaf(k) .eq. '&' .and. ipsp .ne. 1 ) then

            if( imull .gt. 0 .and.
     &          ids .eq. 0 .and. idh .eq. 0 .and.
     &        ( ( ( ipin .eq. 1 .or. ipin .eq. 2 .or.
     &              ipin .eq. 3 ) .and. idv .eq. 1 ) ) ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        rind(imull) = strla(0)

                        indt = indt + 1

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = iftdn

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0


            end if

*-----------------------------------------------------------------------

         else

               if( ifn(ig) .ne. iftdn ) then

                     if( ilt(ig) .gt. 0 ) then

                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                     else

                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                     end if

                        ilt(ig) = 0
                        ifn(ig) = iftdn

                        ifs(ig) = ipfs
                        icl(ig) = ipcl
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi

               else if( ilt(ig) .gt. ichrl - 100 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = iftdn

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

               end if

                  ilt(ig) = ilt(ig) + 1

                  chag(ig,ilt(ig)) = chaf(k)

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

         if( k .lt. icn ) goto 100

*-----------------------------------------------------------------------

               if( ilt(ig) .gt. 0 ) then

                  call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                        ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                        idh,idv,iah,iwt,ibt,ifon,imull,phs)

               else

                  ig = ig - 1

               end if

*-----------------------------------------------------------------------
*     END OF WRITE EACH GROUP
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*     REARRANGE STARTING POSITION
*-----------------------------------------------------------------------

                  strl0 = strla(0)
                  strh0 = strha(0) - strba(0)
                  strb0 = strba(0)

               if( ipin .ge. 1 .and. ipin .le. 3 ) then

                  strl(imull) = strl0

               end if

               if( ipin .eq. 1 ) then

                  strhtp = strh0
                  strbtp = strb0

               else if( ipin .eq. 3 ) then

                  strhbt = strh0
                  strbbt = strb0

               end if

*-----------------------------------------------------------------------
*     ROTATE THE COORDINATE FOR ONE LINE COMMENT
*-----------------------------------------------------------------------

               if( ( ipin .eq. 0 .or. ipin .gt. 4 ) .and.
     &               iang .ne. 0 ) then

                  if( idbg .eq. 1 )  write(jhf,'()')

                  write(jhf,'(3g14.5,'' rotin'')')
     &            xpsin, ypsin, ang

                  xpsain = xpsin
                  ypsain = ypsin
                  angain = ang

                  call bbox(1,1,xpsain,ypsain,angain)

                  xpsin = 0.0
                  ypsin = 0.0

               end if

*-----------------------------------------------------------------------
*     DETERMINE STARTING POSITION
*-----------------------------------------------------------------------

            if( ipin .eq. 0 .or. ipin .eq. -2 ) then

*-----------------------------------------------------------------------
*              NORMAL ONE LINE COMMENT
*-----------------------------------------------------------------------

                  if( ixc .eq. 1 ) then

                     xpss = xpsin

                  else if( ixc .eq. 2 ) then

                     xpss = xpsin - strl0 / 2.0

                  else if( ixc .eq. 3 ) then

                     xpss = xpsin - strl0

                  end if


                  if( iyc .eq. 0 ) then

                     ypss = ypsin

                  else if( iyc .eq. 1 ) then

                     ypss = ypsin - strb0

                  else if( iyc .eq. 2 ) then

                     ypss = ypsin - strb0 - strh0 / 2.0

                  else if( iyc .eq. 3 ) then

                     ypss = ypsin - strb0 - strh0

                  end if

                     call bbox(1,0,xpss,ypss+strha(0),0.d0)
                     call bbox(1,0,xpss,ypss+strba(0),0.d0)
                     call bbox(1,0,xpss+strla(0),ypss+strha(0),0.d0)
                     call bbox(1,0,xpss+strla(0),ypss+strba(0),0.d0)

               if( ipin .eq. -2 ) then

                  write(jhf,'(''} B /mts'',i3.3,'' {'')') imull

               end if

*-----------------------------------------------------------------------

            else if( ipin .gt. 4 ) then

*-----------------------------------------------------------------------
*              ONE LINE COMMENT WITH BOX FROM WT:
*-----------------------------------------------------------------------

                  if( ixc .eq. 1 ) then

                     xpss = xpsin + bds

                  else if( ixc .eq. 2 ) then

                     xpss = xpsin - strl0 / 2.0

                  else if( ixc .eq. 3 ) then

                     xpss = xpsin - bds - strl0

                  end if


                  if( iyc .eq. 0 ) then

                     ypss = ypsin + bds

                  else if( iyc .eq. 1 ) then

                     ypss = ypsin + bds - strb0

                  else if( iyc .eq. 2 ) then

                     ypss = ypsin - strb0 - strh0 / 2.0

                  else if( iyc .eq. 3 ) then

                     ypss = ypsin - bds - strb0 - strh0

                  end if


                     b1x = xpss - bds
                     b2x = xpss + strl0 + bds

                     b1y = ypss + strb0 - bds
                     b2y = ypss + strh0 + strb0 + bds

*-----------------------------------------------------------------------

               iccb = 0

               call wrbox(jhf,idbg,0,ipin,iccb,
     &                    b1x, b2x, b1y, b2y,
     &                    bdd, bcb, bcl, bdc, bcs)

*-----------------------------------------------------------------------

                  if( ipin .ge. 40 ) then

                     b2x = b2x + bdc
                     b1y = b1y - bdc

                  end if

                     call bbox(1,0,b1x,b1y,0.d0)
                     call bbox(1,0,b1x,b2y,0.d0)
                     call bbox(1,0,b2x,b1y,0.d0)
                     call bbox(1,0,b2x,b2y,0.d0)

*-----------------------------------------------------------------------

            else if( ipin .gt. 0 .and. ipin .lt. 4 ) then

*-----------------------------------------------------------------------
*              MULTI LINE COMMENT
*-----------------------------------------------------------------------

               if( ixc .eq. 1 ) then

                  xpp(0) = xpsi(imull)

               else if( ixc .eq. 2 ) then

                  xpp(0) = xpsi(imull) - strl0 / 2.0

               else if( ixc .eq. 3 ) then

                  xpp(0) = xpsi(imull) - strl0

               end if

                  ypp(0) = ypsi(imull)

                  xpsl(imull) = xpp(0)
                  xpsr(imull) = xpp(0) + strl0


                  write(jhf,'(''} B /mts'',i3.3,'' {'')') imull


*-----------------------------------------------------------------------

            else if( ipin .eq. -1 ) then

*-----------------------------------------------------------------------
*              FOR TABLE
*-----------------------------------------------------------------------

                  write(jhf,'(''} B /mts'',i6.6,'' {'')') -imull

            end if

*-----------------------------------------------------------------------
*        STARTING POINT
*-----------------------------------------------------------------------

            if( ipin .eq. 0 .or. ipin .gt. 4 .or. ipin .eq. -2 ) then

                  if( idbg .eq. 1 .and. ipin .ne. -2 )  write(jhf,'()')

               write(jhf,'(''/xps '',g14.5,'' N /yps '',g14.5,'' N'')')
     &                       xpss, ypss

            end if

*-----------------------------------------------------------------------
*     WRITE THE STRINGS
*-----------------------------------------------------------------------


            do 511 i = 1, ig

                  xpp(i) = xpp(ixs(i)) + xps(i)
                  ypp(i) = ypp(iys(i)) + yps(i)

  511       continue


            do 510 i = 1, ig

                  if( ifd(ifn(i)) .eq. 0 .and. ilt(i) .gt. 0 )
     &            call ftdfn(jhf,ifon,ifn(i),idbg)

                  call write00(10,idbg,jhf,i,chag,igs,ib,icl,clp,
     &                         ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                         idh,idv,iah,iwt,ibt,ifon,imull,phs)

  510       continue


*-----------------------------------------------------------------------

            if( ( ipin .gt. 0 .and. ipin .lt. 4 ) .or.
     &            ipin .eq. -1 .or. ipin .eq. -2 ) then

                  write(jhf,'(''} B'')')

            end if

*-----------------------------------------------------------------------
*     UN ROTATE THE COORDINATE FOR ONE LINE COMMENT
*-----------------------------------------------------------------------

               if( ( ipin .eq. 0 .or. ipin .gt. 4 ) .and.
     &               iang .ne. 0 ) then

                  if( idbg .eq. 1 )  write(jhf,'()')

                  write(jhf,'(''rotout'')')

                  call bbox(1,1,0.d0,0.d0,-angain)
                  call bbox(1,1,-xpsain,-ypsain,0.d0)

               end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function chars1(jhf,idbg,k,chaf,icn,chag,ig,
     &                ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                idh,iacd,idv,iah,iwt,ibt,ifon,imull)
*                                                                      *
*              SUFFIXES AFTER '^', '_'                                 *
*                                                                      *
*              {^{  }_{  }}                                            *
*                                                                      *
*                    |__IAP(IDS,2)                                     *
*               |__IAP(IDS,1)                                          *
*                                                                      *
*              |__IGS(IB); IG                                          *
*              |__IBS(IB); K                                           *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(inig,0:ichrl)*1
      dimension ifn(0:inig), ifs(0:inig), ilt(0:inig), icl(0:inig)
      dimension clp(0:inig,3,3)
      dimension ixp(0:inig), iyp(0:inig), isi(0:inig), iap(0:inig,3)
      dimension ibf(0:inig), ibs(0:inig), igs(0:inig),
     &          ibb(0:inig), ibc(0:inig)
      dimension iwt(0:inig), ibt(0:inig,2), iah(0:inig,6)
      dimension isuf(0:ichrl,3)
      character iacd(0:inig)*3

      logical chars1
      save ipfn

*-----------------------------------------------------------------------

         if( isuf(k,1) .ne. 0 ) then

*-----------------------------------------------------------------------

               if(  ( isuf(k,1) .eq. 3 .or. isuf(k,1) .eq. 4 ) .and.
     &                k .ne. 1 .and. ilt(ig) .gt. 0 ) then

                     ilt(ig) = ilt(ig) - 1


               end if

*-----------------------------------------------------------------------

                        ipfn    = ibf(ib)
                        ib      = ib + 1

*-----------------------------------------------------------------------
*           CAUTION !!!  08.10.97
*-----------------------------------------------------------------------

                  if( ilt(ig) .gt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1


                  else if( ( iah(idv,4) .ne. 60 .and.
     &                       iah(idv,4) .ne. 61 ) .and.
     &                   ( ifs(ig) .ne. ig - 1 .or.
     &                     icl(ig) .ne. ig - 1 .or.
     &                     ixp(ig) .ne. ig     .or.
     &                     iyp(ig) .ne. ig - 1 .or.
     &                     isi(ig) .ne. 0 .or. idh .eq. 2 ) ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                  end if

                        ibf(ib) = ipfn
                        ibs(ib) = k
                        igs(ib) = ig

                        ilt(ig) = 0
                        ifn(ig) = ibf(ib)

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

*-----------------------------------------------------------------------

                  if( chaf(k+1) .eq. '^' ) then

                     if( ibt(ig,1) .eq. 4 ) then

                        isi(ig) = 110

                     else

                        isi(ig) = 1

                     end if

                  else

                     if( ibt(ig,1) .eq. 4 ) then

                        isi(ig) = -110

                     else

                        isi(ig) = -1

                     end if

                  end if


               if( isuf(k,1) .eq. 3 .or. isuf(k,1) .eq. 4 ) then

                        ibt(ig,1)  = 0
                        ibt(ig,2)  = 0

               end if


               if( isuf(k,1) .eq. 2 .or. isuf(k,1) .eq. 4 ) then

                        ids = ids + 1

                        iap(ids,2) = 0

                  if( isuf(k,1) .eq. 2 ) then

                        iap(ids,3) = 2

                  else if( isuf(k,1) .eq. 4 ) then

                        iap(ids,3) = 4

                  end if

                  if( chaf(k+1) .eq. '^' ) then

                        iap(ids,1) = ig

                  else

                        iap(ids,1) = - ig

                  end if

                        call write00(24,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------

               else if( ibt(ig,1) .ge. 1 ) then

                        ids = ids + 1

                        iap(ids,2) = 0
                        iap(ids,3) = 5

                     if( chaf(k+1) .eq. '^' ) then

                        iap(ids,1) = ig

                     else

                        iap(ids,1) = - ig

                     end if

                        call write00(24,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

               end if

*-----------------------------------------------------------------------

               if( ibt(ig,1) .eq. 3 .or. ibt(ig,1) .eq. 4 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        idv = idv + 1

                        iah(idv,1) = ig
                        iah(idv,2) = ids
                        iah(idv,3) = ipfn

                        iah(idv,4) = 6

                        call write00(8,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

               end if

*-----------------------------------------------------------------------


                        k = k + 1


                        chars1 = .true.


*-----------------------------------------------------------------------

         else

                  chars1 = .false.

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function chars2(jhf,idbg,k,chaf,icn,chag,ig,
     &                ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                idh,iacd,idv,iah,iwt,ibt,ifon,imull,
     &                ipsp)
*                                                                      *
*              CHANGE FONT '{\sf','{\hv','{\rm','{\tm','{\tt','{\cr',  *
*                          '{\it','{\bf','{\ib','{\sb','{\mc','{\gt'   *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(inig,0:ichrl)*1
      dimension ifn(0:inig), ifs(0:inig), ilt(0:inig), icl(0:inig)
      dimension clp(0:inig,3,3)
      dimension ixp(0:inig), iyp(0:inig), isi(0:inig), iap(0:inig,3)
      dimension ibf(0:inig), ibs(0:inig), igs(0:inig),
     &          ibb(0:inig), ibc(0:inig)
      dimension iwt(0:inig), ibt(0:inig,2), iah(0:inig,6)
      dimension isuf(0:ichrl,3)
      character iacd(0:inig)*3

      logical chars2
      save ipfn

*-----------------------------------------------------------------------

      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

         if( chaf(k) .eq. '{' .and. chaf(k+1) .eq. yen .and.
     &       ipsp .ne. 1 .and.
     &     ( chaf(k+2)//chaf(k+3) .eq. 'rm' .or.
     &       chaf(k+2)//chaf(k+3) .eq. 'tm' .or.
     &       chaf(k+2)//chaf(k+3) .eq. 'sf' .or.
     &       chaf(k+2)//chaf(k+3) .eq. 'hv' .or.
     &       chaf(k+2)//chaf(k+3) .eq. 'tt' .or.
     &       chaf(k+2)//chaf(k+3) .eq. 'cr' .or.
     &       chaf(k+2)//chaf(k+3) .eq. 'it' .or.
     &       chaf(k+2)//chaf(k+3) .eq. 'bf' .or.
     &       chaf(k+2)//chaf(k+3) .eq. 'ib' .or.
     &       chaf(k+2)//chaf(k+3) .eq. 'mc' .or.
     &       chaf(k+2)//chaf(k+3) .eq. 'gt' .or.
     &       chaf(k+2)//chaf(k+3) .eq. 'sb' )  ) then

                     ipfn    = ibf(ib)

                     ib      = ib + 1

                     ibs(ib) = 0


               if( chaf(k+2)//chaf(k+3) .eq. 'sf' .or.
     &             chaf(k+2)//chaf(k+3) .eq. 'hv' .or.
     &             chaf(k+2)//chaf(k+3) .eq. 'gt' ) then

                     ibf(ib) = 0

               else if( chaf(k+2)//chaf(k+3) .eq. 'rm' .or.
     &                  chaf(k+2)//chaf(k+3) .eq. 'tm' .or.
     &                  chaf(k+2)//chaf(k+3) .eq. 'mc' ) then

                     ibf(ib) = 5

               else if( chaf(k+2)//chaf(k+3) .eq. 'tt' .or.
     &                  chaf(k+2)//chaf(k+3) .eq. 'cr' ) then

                     ibf(ib) = 9

               else if( chaf(k+2)//chaf(k+3) .eq. 'it' ) then

                  if( ipfn .ge. 0 .and. ipfn .le. 3 ) then

                     ibf(ib) = 1

                  else if( ipfn .ge. 5 .and. ipfn .le. 8 ) then

                     ibf(ib) = 6

                  else if( ipfn .ge. 9 .and. ipfn .le. 12 ) then

                     ibf(ib) = 10

                  end if

               else if( chaf(k+2)//chaf(k+3) .eq. 'bf' ) then

                  if( ipfn .ge. 0 .and. ipfn .le. 3 ) then

                     ibf(ib) = 2

                  else if( ipfn .ge. 5 .and. ipfn .le. 8 ) then

                     ibf(ib) = 7

                  else if( ipfn .ge. 9 .and. ipfn .le. 12 ) then

                     ibf(ib) = 11

                  end if

               else if( chaf(k+2)//chaf(k+3) .eq. 'ib' ) then

                  if( ipfn .ge. 0 .and. ipfn .le. 3 ) then

                     ibf(ib) = 3

                  else if( ipfn .ge. 5 .and. ipfn .le. 8 ) then

                     ibf(ib) = 8

                  else if( ipfn .ge. 9 .and. ipfn .le. 12 ) then

                     ibf(ib) = 12

                  end if

               else if( chaf(k+2)//chaf(k+3) .eq. 'sb' ) then

                     ibf(ib) = 4

               end if

*-----------------------------------------------------------------------

               iftdn = ibf(ib)

               if( ibf(ib) .ne. ipfn ) then

                     if( ilt(ig) .gt. 0 ) then

                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                     else

                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                     end if

                        ilt(ig) = 0
                        ifn(ig) = ibf(ib)

                        ifs(ig) = ipfs
                        icl(ig) = ipcl
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi

               end if

*-----------------------------------------------------------------------

               if( chaf(k+4) .eq. ' ' ) then

                     k = k + 4

               else

                     k = k + 3

               end if

*-----------------------------------------------------------------------


                  chars2 = .true.

         else

                  chars2 = .false.

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function chars3(jhf,idbg,k,chaf,icn,chag,ig,
     &                ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                idh,iacd,idv,iah,iwt,ibt,ifon,imull)
*                                                                      *
*              CHANGE FONT '{\tiny','{\scriptsize','{\footnotesize',   *
*                                   '{\Tiny',      '{\TINY',           *
*                          '{\small',                                  *
*                          '{\large','{\Large','{\LARGE',              *
*                          '{\huge', '{\Huge'                          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(inig,0:ichrl)*1
      dimension ifn(0:inig), ifs(0:inig), ilt(0:inig), icl(0:inig)
      dimension clp(0:inig,3,3)
      dimension ixp(0:inig), iyp(0:inig), isi(0:inig), iap(0:inig,3)
      dimension ibf(0:inig), ibs(0:inig), igs(0:inig),
     &          ibb(0:inig), ibc(0:inig)
      dimension iwt(0:inig), ibt(0:inig,2), iah(0:inig,6)
      dimension isuf(0:ichrl,3)
      character iacd(0:inig)*3

      character dum*2048
*         ( ichrl = 2048 ) should be the same as dum

      logical chars3
      save ipfn

*-----------------------------------------------------------------------

      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

            do 200 i = 1, ichrl

               dum(i:i) = chaf(i)

  200       continue

*-----------------------------------------------------------------------

         if( chaf(k) .eq. '{' .and. chaf(k+1) .eq. yen .and.
     &     ( dum(k+2:k+ 5) .eq. 'tiny' .or.
     &       dum(k+2:k+ 5) .eq. 'Tiny' .or.
     &       dum(k+2:k+ 5) .eq. 'TINY' .or.
     &       dum(k+2:k+11) .eq. 'scriptsize' .or.
     &       dum(k+2:k+13) .eq. 'footnotesize' .or.
     &       dum(k+2:k+ 6) .eq. 'small' .or.
     &       dum(k+2:k+ 6) .eq. 'large' .or.
     &       dum(k+2:k+ 6) .eq. 'Large' .or.
     &       dum(k+2:k+ 6) .eq. 'LARGE' .or.
     &       dum(k+2:k+ 5) .eq. 'huge' .or.
     &       dum(k+2:k+ 5) .eq. 'Huge' ) ) then


                     ipfn    = ibf(ib)

                     ib      = ib + 1

                     ibs(ib) = 0
                     ibf(ib) = ipfn


               if( dum(k+2:k+5) .eq. 'tiny' ) then

                     ibbs    = 101
                     ibbl    = 4

               else if( dum(k+2:k+5) .eq. 'Tiny' ) then

                     ibbs    = 102
                     ibbl    = 4

               else if( dum(k+2:k+11) .eq. 'scriptsize' ) then

                     ibbs    = 102
                     ibbl    = 10

               else if( dum(k+2:k+5) .eq. 'TINY' ) then

                     ibbs    = 103
                     ibbl    = 4

               else if( dum(k+2:k+13) .eq. 'footnotesize' ) then

                     ibbs    = 103
                     ibbl    = 12

               else if( dum(k+2:k+ 6) .eq. 'small' ) then

                     ibbs    = 104
                     ibbl    = 5

               else if( dum(k+2:k+ 6) .eq. 'large' ) then

                     ibbs    = 105
                     ibbl    = 5

               else if( dum(k+2:k+ 6) .eq. 'Large' ) then

                     ibbs    = 106
                     ibbl    = 5

               else if( dum(k+2:k+ 6) .eq. 'LARGE' ) then

                     ibbs    = 107
                     ibbl    = 5

               else if( dum(k+2:k+ 5) .eq. 'huge' ) then

                     ibbs    = 108
                     ibbl    = 4

               else if( dum(k+2:k+ 5) .eq. 'Huge' ) then

                     ibbs    = 109
                     ibbl    = 4

               end if

*-----------------------------------------------------------------------

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                        ibb(ib) = ig - 1

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = ibbs

*-----------------------------------------------------------------------

               if( chaf(k+ibbl+2) .eq. ' ' ) then

                     k = k + ibbl + 2

               else

                     k = k + ibbl + 1

               end if

*-----------------------------------------------------------------------


                  chars3 = .true.

         else

                  chars3 = .false.

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function chars4(jhf,idbg,clal,clmo,k,chaf,icn,chag,ig,
     &                ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                idh,iacd,idv,iah,iwt,ibt,ifon,imull)
*                                                                      *
*             CHANGE COLOR '{\color{r} text}'                          *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(inig,0:ichrl)*1
      dimension ifn(0:inig), ifs(0:inig), ilt(0:inig), icl(0:inig)
      dimension clp(0:inig,3,3)
      dimension ixp(0:inig), iyp(0:inig), isi(0:inig), iap(0:inig,3)
      dimension ibf(0:inig), ibs(0:inig), igs(0:inig),
     &          ibb(0:inig), ibc(0:inig)
      dimension iwt(0:inig), ibt(0:inig,2), iah(0:inig,6)
      dimension isuf(0:ichrl,3)
      character iacd(0:inig)*3

      logical chars4
      save ipfn

      dimension clal(3), colc(3)

*-----------------------------------------------------------------------

      character yen*1

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------

            yen = char(92)

*-----------------------------------------------------------------------

               colc(1) = -r1max
               colc(2) = 1.0
               colc(3) = 1.0

               ierr = 0

*-----------------------------------------------------------------------

         if( chaf(k) .eq. '{' .and. chaf(k+1) .eq. yen .and.
     &     ( chaf(k+2)//chaf(k+3)//chaf(k+4)//chaf(k+5)//chaf(k+6)
     &      .eq. 'color' ) ) then

                     ki = k + 6
                     kf = k + 6

                        call dcols(1,chaf,ki,colc,ierr,icn,'{','}')

                  if( ierr .eq. 0 ) then

                     kf = ki

                  else

                     kf = kf

                  end if

*-----------------------------------------------------------------------

                     ipfn    = ibf(ib)

                     ib      = ib + 1

                     ibs(ib) = 0
                     ibf(ib) = ipfn


                  if( ierr .ne. 0 ) then

                     k = kf

                     chars4 = .true.

                     return

                  end if

*-----------------------------------------------------------------------

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                        ibc(ib) = ig - 1

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

*-----------------------------------------------------------------------
*        COLOR DEFINITION
*-----------------------------------------------------------------------

                  if( colc(1) .gt. -r0max .and. clmo .lt. -r0max ) then

                     if( clal(1) .gt. -r0max ) then
                        colc(1) = clal(1)
                        colc(2) = clal(2)
                        colc(3) = clal(3)
                     end if

                        clp(ig,1,1) = colc(1)
                        clp(ig,1,2) = colc(2)
                        clp(ig,1,3) = colc(3)

                        icl(ig) = ig

                  else

                        icl(ig) = ig - 1

                  end if

*-----------------------------------------------------------------------

               if( chaf(kf) .eq. ' ' ) then

                     k = kf

               else

                     k = kf - 1

               end if

*-----------------------------------------------------------------------


                  chars4 = .true.

         else

                  chars4 = .false.

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function chars5(jhf,idbg,clal,clmo,k,chaf,icn,chag,ig,
     &                ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                idh,iacd,idv,iah,iwt,ibt,ifon,imull)
*                                                                      *
*             BOX AND COLOR BOX                                        *
*                          '{\singlebox{c}{b}{l}     text}'            *
*                          '{\ovalbox{c}{b}{l}       text}'            *
*                          '{\doublebox{c}{b}{l}     text}'            *
*                          '{\shadowbox{c}{b}{s}     text}'            *
*                          '{\ovalshadowbox{c}{b}{s} text}'            *
*                                                                      *
*           IAH(IDV,4) = 10 ; singlebox                                *
*                        11 ; Singlebox                                *
*                        12 ; singleBox                                *
*                        13 ; SingleBox                                *
*                        20 ; ovalbox                                  *
*                        21 ; Ovalbox                                  *
*                        22 ; ovalBox                                  *
*                        23 ; OvalBox                                  *
*                        30 ; doublebox                                *
*                        31 ; Doublebox                                *
*                        32 ; doubleBox                                *
*                        33 ; DoubleBox                                *
*                        40 ; shadowbox                                *
*                        41 ; Shadowbox                                *
*                        42 ; shadowBox                                *
*                        43 ; ShadowBox                                *
*                        50 ; ovalshadowbox                            *
*                        51 ; Ovalshadowbox                            *
*                        52 ; ovalshadowBox                            *
*                        53 ; OvalshadowBox                            *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(inig,0:ichrl)*1
      dimension ifn(0:inig), ifs(0:inig), ilt(0:inig), icl(0:inig)
      dimension clp(0:inig,3,3)
      dimension ixp(0:inig), iyp(0:inig), isi(0:inig), iap(0:inig,3)
      dimension ibf(0:inig), ibs(0:inig), igs(0:inig),
     &          ibb(0:inig), ibc(0:inig)
      dimension iwt(0:inig), ibt(0:inig,2), iah(0:inig,6)
      dimension isuf(0:ichrl,3)
      character iacd(0:inig)*3

      character dum*2048
*         ( ichrl = 2048 ) should be the same as dum

      logical chars5
      save ipfn

      dimension clal(3), colc(3)

*-----------------------------------------------------------------------

      character yen*1

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------

            yen = char(92)

*-----------------------------------------------------------------------

            do 200 i = 1, ichrl

               dum(i:i) = chaf(i)

  200       continue

*-----------------------------------------------------------------------

         if( chaf(k) .eq. '{' .and. chaf(k+1) .eq. yen .and.
     &     ( ( dum(k+2:k+10) .eq. 'singlebox' ) .or.
     &       ( dum(k+2:k+10) .eq. 'Singlebox' ) .or.
     &       ( dum(k+2:k+10) .eq. 'singleBox' ) .or.
     &       ( dum(k+2:k+10) .eq. 'SingleBox' ) .or.
     &       ( dum(k+2:k+ 8) .eq. 'ovalbox'   ) .or.
     &       ( dum(k+2:k+ 8) .eq. 'Ovalbox'   ) .or.
     &       ( dum(k+2:k+ 8) .eq. 'ovalBox'   ) .or.
     &       ( dum(k+2:k+ 8) .eq. 'OvalBox'   ) .or.
     &       ( dum(k+2:k+10) .eq. 'doublebox' ) .or.
     &       ( dum(k+2:k+10) .eq. 'Doublebox' ) .or.
     &       ( dum(k+2:k+10) .eq. 'doubleBox' ) .or.
     &       ( dum(k+2:k+10) .eq. 'DoubleBox' ) .or.
     &       ( dum(k+2:k+10) .eq. 'shadowbox' ) .or.
     &       ( dum(k+2:k+10) .eq. 'Shadowbox' ) .or.
     &       ( dum(k+2:k+10) .eq. 'shadowBox' ) .or.
     &       ( dum(k+2:k+10) .eq. 'ShadowBox' ) .or.
     &       ( dum(k+2:k+14) .eq. 'ovalshadowbox' ) .or.
     &       ( dum(k+2:k+14) .eq. 'Ovalshadowbox' ) .or.
     &       ( dum(k+2:k+14) .eq. 'ovalshadowBox' ) .or.
     &       ( dum(k+2:k+14) .eq. 'OvalshadowBox' ) )
     &     ) then

               if( dum(k+2:k+10) .eq. 'singlebox' ) then

                     iboxd = 10
                     iboxl = 9

               else if( dum(k+2:k+10) .eq. 'Singlebox' ) then

                     iboxd = 11
                     iboxl = 9

               else if( dum(k+2:k+10) .eq. 'singleBox' ) then

                     iboxd = 12
                     iboxl = 9

               else if( dum(k+2:k+10) .eq. 'SingleBox' ) then

                     iboxd = 13
                     iboxl = 9

               else if( dum(k+2:k+ 8) .eq. 'ovalbox' ) then

                     iboxd = 20
                     iboxl = 7

               else if( dum(k+2:k+ 8) .eq. 'Ovalbox' ) then

                     iboxd = 21
                     iboxl = 7

               else if( dum(k+2:k+ 8) .eq. 'ovalBox' ) then

                     iboxd = 22
                     iboxl = 7

               else if( dum(k+2:k+ 8) .eq. 'OvalBox' ) then

                     iboxd = 23
                     iboxl = 7

               else if( dum(k+2:k+10) .eq. 'doublebox' ) then

                     iboxd = 30
                     iboxl = 9

               else if( dum(k+2:k+10) .eq. 'Doublebox' ) then

                     iboxd = 31
                     iboxl = 9

               else if( dum(k+2:k+10) .eq. 'doubleBox' ) then

                     iboxd = 32
                     iboxl = 9

               else if( dum(k+2:k+10) .eq. 'DoubleBox' ) then

                     iboxd = 33
                     iboxl = 9

               else if( dum(k+2:k+10) .eq. 'shadowbox' ) then

                     iboxd = 40
                     iboxl = 9

               else if( dum(k+2:k+10) .eq. 'Shadowbox' ) then

                     iboxd = 41
                     iboxl = 9

               else if( dum(k+2:k+10) .eq. 'shadowBox' ) then

                     iboxd = 42
                     iboxl = 9

               else if( dum(k+2:k+10) .eq. 'ShadowBox' ) then

                     iboxd = 43
                     iboxl = 9

               else if( dum(k+2:k+14) .eq. 'ovalshadowbox' ) then

                     iboxd = 50
                     iboxl = 13

               else if( dum(k+2:k+14) .eq. 'Ovalshadowbox' ) then

                     iboxd = 51
                     iboxl = 13

               else if( dum(k+2:k+14) .eq. 'ovalshadowBox' ) then

                     iboxd = 52
                     iboxl = 13

               else if( dum(k+2:k+14) .eq. 'OvalshadowBox' ) then

                     iboxd = 53
                     iboxl = 13

               end if

*-----------------------------------------------------------------------

                        ipfn = ibf(ib)

                     if( ilt(ig). gt. 0 ) then

                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                     else

                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                     end if

*-----------------------------------------------------------------------

                        iwt(ig) = iboxd

*-----------------------------------------------------------------------

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ipfs
                        icl(ig) = ipcl
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi


                     if( iboxd .eq. 10 .or. iboxd .eq. 11 .or.
     &                   iboxd .eq. 20 .or. iboxd .eq. 21 .or.
     &                   iboxd .eq. 30 .or. iboxd .eq. 31 .or.
     &                   iboxd .eq. 40 .or. iboxd .eq. 41 .or.
     &                   iboxd .eq. 50 .or. iboxd .eq. 51 ) then

                        idh = -10

                     else

                        idh = -11

                     end if

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        idh = 0

*-----------------------------------------------------------------------

                        ig      = ig + 1

                        ib      = ib + 1

                        ibs(ib) = 0
                        ibf(ib) = ipfn

                        igs(ib) = ig - 1

                        idv     = idv + 1

                        iah(idv,1) = ig - 1
                        iah(idv,2) = ids
                        iah(idv,3) = ipfn

                        iah(idv,4) = iboxd

                        call write00(8,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------
*        COLOR DEFINITION
*-----------------------------------------------------------------------

                        ki = k + 1 + iboxl
                        kf = ki + 1

                  do 100 jj = 1, 3

                        call dcols(1,chaf,ki,colc,ierr,icn,'{','}')

                        if( ierr .ne. 0 ) then

                           goto 101

                        end if

                           clp(ig,jj,1) = colc(1)
                           clp(ig,jj,2) = colc(2)
                           clp(ig,jj,3) = colc(3)

                           kf = ki
                           ki = ki - 1

  100             continue

  101             continue

*-----------------------------------------------------------------------

                  if( clp(ig,1,1) .gt. -r0max .and.
     &                clal(1) .gt. -r0max ) then
                        clp(ig,1,1) = clal(1)
                        clp(ig,1,2) = clal(2)
                        clp(ig,1,3) = clal(3)
                  end if

                  if( clp(ig,3,1) .gt. -r0max .and.
     &                clal(1) .gt. -r0max ) then
                        clp(ig,3,1) = clal(1)
                        clp(ig,3,2) = clal(2)
                        clp(ig,3,3) = clal(3)
                  end if

                  if( clp(ig,2,1) .lt. -r0max .or.
     &                clmo .gt. -r0max .or.
     &                clal(1) .gt. -r0max ) clp(ig,2,1) = -1.0

                  if( clp(ig,3,1) .lt. -r0max .or.
     &                clmo .gt. -r0max ) clp(ig,3,1) = -3.0

*-----------------------------------------------------------------------

                        ibc(ib) = ig - 1

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

                     if( clp(ig,1,1) .lt. -r0max .or.
     &                   clmo .gt. -r0max ) then

                        icl(ig) = ig - 1

                     else

                        icl(ig) = ig

                     end if

*-----------------------------------------------------------------------

               if( chaf(kf) .eq. ' ' ) then

                     k = kf

               else

                     k = kf - 1

               end if

*-----------------------------------------------------------------------

                  chars5 = .true.

         else

                  chars5 = .false.

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function chars6(jhf,idbg,k,chaf,icn,chag,ig,
     &                ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                idh,iacd,idv,iah,iwt,ibt,ifon,imull)
*                                                                      *
*           FOR SPECIAL CHARACTER AFTER '\'                            *
*                                                                      *
*           IAH(IDV,1) = IG                                            *
*           IAH(IDV,2) = IDS                                           *
*           IAH(IDV,3) = IFN(IG)                                       *
*                                                                      *
*           IAH(IDV,4) = 0  ; NORMAL ACCENT GROUP                      *
*                      = 1  ; VECTOR                                   *
*                      = 2  ; OVERLINE                                 *
*                      = 3  ; UNDERLINE                                *
*                      = 4  ; OVER RIGHT ARROW                         *
*                      = 5  ; OVER LEFT ARROW                          *
*                      = 6  ; SUM, PROD, LIM                           *
*                      = 7  ; FRAC                                     *
*                      = 8  ; FRACS                                    *
*                      = 9  ; SQRT                                     *
*                      = 60 ; MATHOP                                   *
*                      = 61 ; KS ( KANJI SUFIX, FURIGANA )             *
*                                                                      *
*           IAH(IDV,5) = IAH(IDV,1) FOR \frac{ }{ }                    *
*                                                                      *
*           ISI(I)     = |1| ; SUFFIX                                  *
*                      =  2  ; HAT                                     *
*                      =  3  ; VECTOR                                  *
*                      =  4  ; HBAR                                    *
*                      =  5  ; UNDERLINE                               *
*                      =  6 7 8   ; Big                                *
*                      =  9 10 11 ; big                                *
*                      =  21 ; fracs                                   *
*                                                                      *
*                      >  100 ; CHANGE FONT SIZE                       *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      parameter ( icmax = 15, iscn = 44 )
      parameter ( ic01 = 38 )
      parameter ( ic02 = 26 )
      parameter ( ic03 = 44 )
      parameter ( ic04 = 39 )
      parameter ( ic05 = 41 )
      parameter ( ic06 = 34 )
      parameter ( ic07 = 15 )
      parameter ( ic08 = 9 )
      parameter ( ic09 = 10 )
      parameter ( ic010 = 5 )
      parameter ( ic011 = 1 )
      parameter ( ic013 = 1 )
      parameter ( ic014 = 3 )
      parameter ( ic015 = 1 )
      parameter ( ic012 = 0 )

*-----------------------------------------------------------------------

      common /wtval2/ clc(0:inig,3), fts(0:inig),
     &                stw(0:inig), sth(0:inig), stb(0:inig),
     &                xps(0:inig), yps(0:inig),
     &                xpp(0:inig), ypp(0:inig),
     &                ixs(0:inig), iys(0:inig), jfn(0:inig)
      common /wtval7/ strhtp, strbtp, strhbt, strbbt, strhall
      common /wtval8/ xpsi(0:inig), ypsi(0:inig), strl(0:inig),
     &                xpsl(0:inig), xpsr(0:inig), rind(0:inig)

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

      character chaf(ichrl)*1
      character chag(inig,0:ichrl)*1
      dimension ifn(0:inig), ifs(0:inig), ilt(0:inig), icl(0:inig)
      dimension clp(0:inig,3,3)
      dimension ixp(0:inig), iyp(0:inig), isi(0:inig), iap(0:inig,3)
      dimension ibf(0:inig), ibs(0:inig), igs(0:inig),
     &          ibb(0:inig), ibc(0:inig)
      dimension iwt(0:inig), ibt(0:inig,2), iah(0:inig,6)
      dimension isuf(0:ichrl,3)
      character iacd(0:inig)*3

      character cdum*10

      logical chars6
      logical numec
      logical alphb

*-----------------------------------------------------------------------

      character     today1*11, timestr1*5, today*100, timestr*100
      character     infnm(ichrl)*1, infn(ichrl)*1, avers*5
      common /dfil1/today, timestr, today1, timestr1, infnm, infn, avers
      common /dfil2/jtoday, jtime, jnm, mpage, iname, javers

*-----------------------------------------------------------------------

      character isex(icmax,iscn)*20
      character isch(icmax,iscn)*3
      dimension isfn(icmax,iscn)
      dimension jsfn(icmax,iscn)
      dimension nic(icmax)

      character dum20*20

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic01),i=1,1)
     &         /'+',  '-',  '#',  '%',  '&',  'S',  'P',  'o',  'O',
     &          'l',  'L',  'i',  '1',  '2',  '3',  '4',  '5',  '6',
     &          '7',  '8',  '9',  '0',  ',',  '>',  ';',  '!',  '$',
     &          '"',  '`',  "'",  '~',  '=',  '.',  '^',  '_',  '{',
     &          '}',  '|'/
      data ((isch(i,j),j=1,ic01),i=1,1)
     &         /'053','055','043','045','046','247','266','371','351',
     &          '370','350','365','061','062','063','064','065','066',
     &          '067','070','071','060','000','000','000','000','044',
     &          '310','301','302','304','305','307','136','137','173',
     &          '175','174'/
      data ((isfn(i,j),j=1,ic01),i=1,1)
     &         /4,     4,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
     &         -1,    -1,   -1,   -4,   -4,   -4,   -4,   -4,   -4,
     &         -4,    -4,   -4,   -4,   -1,   -1,   -1,   -1,   -1,
     &         -1,    -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
     &         -1,     4/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic02),i=2,2)
     &         /'mu', 'nu', 'xi', 'pi', 'Xi', 'Pi', 'pm', 'mp', 'le',
     &          'in', 'ge', 'ne', 'to', 'wp', 'Re', 'Im', 'oe', 'OE',
     &          'ae', 'AE', 'ss', 'eq', 'ln', 'Pr', 'lt', 'gt'/
      data ((isch(i,j),j=1,ic02),i=2,2)
     &         /'155','156','170','160','130','120','261','261','243',
     &          '316','263','271','256','303','302','301','372','352',
     &          '361','341','373','075','000','000','074','076'/
      data ((isfn(i,j),j=1,ic02),i=2,2)
     &         /13,   13,   13,   13,    4,    4,    4,    4,    4,
     &          4,     4,    4,    4,    4,    4,    4,   -1,   -1,
     &         -1,    -1,   -1,    4,   -4,   -4,    4,    4/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic03),i=3,3)
     &         /'eta','rho','tau','phi','chi','psi','Phi','Psi',
     &          'div','ast','cap','cup','vee','lor','sim','neq',
     &          'mid','bot','neg','sum','int','dag','yen','lim',
     &          'arg','cos','cot','csc','deg','det','dim','exp',
     &          'gcd','hom','inf','ker','log','max','min','sec',
     &          'sin','sup','tan','ks{'/
      data ((isch(i,j),j=1,ic03),i=3,3)
     &         /'150','162','164','146','143','171','106','131',
     &          '270','052','307','310','332','332','176','271',
     &          '275','136','330','345','362','262','245','lim',
     &          '000','000','000','000','000','000','000','000',
     &          '000','000','000','000','000','000','000','000',
     &          '000','000','000','000'/
      data ((isfn(i,j),j=1,ic03),i=3,3)
     &         /13,   13,   13,   13,   13,   13,    4,    4,
     &          4,     4,    4,    4,    4,    4,    4,    4,
     &          4,     4,    4,    4,    4,   -1,   -1,   -3,
     &         -4,    -4,   -4,   -4,   -4,   -4,   -4,   -4,
     &         -4,    -4,   -4,   -4,   -4,   -4,   -4,   -4,
     &         -4,    -4,   -4,   -1/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic04),i=4,4)
     &         /'beta','zeta','iota','cdot','land','plus','cong','perp',
     &          'gets','surd','lnot','prod','ddag','hat{','bar{','dot{',
     &          'vec{','Big[','big[','Big]','big]','Big(','big(','Big)',
     &          'big)','sums','lims','quad','cosh','coth','sinh','tanh',
     &          'circ','time','file','page','hbar','Ang{','ang{'/
      data ((isch(i,j),j=1,ic04),i=4,4)
     &         /'142', '172', '151', '264', '331', '053', '100', '136',
     &          '254', '326', '330', '325', '263', '303', '305', '307',
     &          '256', '352', '352', '372', '372', '347', '347', '367',
     &          '367', '345', 'lim', '000', '000', '000', '000', '000',
     &          '260', '000', '000', '000', '305', '260', '260'/
      data ((isfn(i,j),j=1,ic04),i=4,4)
     &         /13,    13,    13,    -1,     4,     4,     4,     4,
     &          4,      4,     4,     4,    -1,    -2,    -2,    -2,
     &          4,      4,     4,     4,     4,     4,     4,     4,
     &          4,      4,    -3,    -1,    -4,    -4,    -4,    -4,
     &          4,     -1,    -1,    -1,    -1,     4,     4/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic05),i=5,5)
     &         /'alpha','gamma','delta','theta','kappa','omega',
     &          'varpi','Gamma','Delta','Theta','Sigma','Omega',
     &          'times','wedge','oplus','minus','notin','equiv',
     &          'aleph','imath','infty','prime','nabla','angle',
     &          'lceil','rceil','ddot{','Big@{','big@{','Big@}',
     &          'big@}','prods','frac{','sqrt{','equal','left[',
     &          'left(','qquad','left<','sigma','today'/
      data ((isch(i,j),j=1,ic05),i=5,5)
     &         /'141',  '147',  '144',  '161',  '153',  '167',
     &          '166',  '107',  '104',  '121',  '123',  '127',
     &          '264',  '331',  '305',  '055',  '317',  '272',
     &          '300',  '365',  '245',  '242',  '321',  '320',
     &          '351',  '371',  '310',  '355',  '355',  '375',
     &          '375',  '325',  '055',  '326',  '075',  '133',
     &          '050',  '000',  '341',  '163',  '000'/
      data ((isfn(i,j),j=1,ic05),i=5,5)
     &         /13,     13,     13,     13,     13,     13,
     &          13,      4,      4,      4,      4,      4,
     &          4,       4,      4,      4,      4,      4,
     &          4,      -1,      4,      4,      4,      4,
     &          4,       4,     -2,      4,      4,      4,
     &          4,       4,      4,      4,      4,      4,
     &          4,      -1,      4,     13,     -1/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic06),i=6,6)
     &         /'lambda','varrho','varphi','Lambda','bullet','otimes',
     &          'dagger','subset','supset','approx','propto','forall',
     &          'exists','pounds','lfloor','rfloor','langle','rangle',
     &          'check{','breve{','acute{','grave{','tilde{','Bigint',
     &          'bigint','fracs{','right]','right)','left@{','right>',
     &          'arccos','arcsin','arctan','oslash'/
      data ((isch(i,j),j=1,ic06),i=6,6)
     &         /'154',   '162',   '152',   '114',  '267',   '304',
     &          '262',   '314',   '311',   '273',  '265',   '042',
     &          '044',   '243',   '353',   '373',  '341',   '361',
     &          '317',   '306',   '302',   '301',  '304',   '364',
     &          '364',   '055',   '135',   '051',  '173',   '361',
     &          '000',   '000',   '000',   '306'/
      data ((isfn(i,j),j=1,ic06),i=6,6)
     &         /13,      13,      13,       4,      4,       4,
     &         -1,        4,       4,       4,      4,       4,
     &          4,       -1,       4,       4,      4,       4,
     &         -2,       -2,      -2,      -2,     -2,       4,
     &          4,        4,       4,       4,      4,       4,
     &         -4,       -4,      -4,       4/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic07),i=7,7)
     &         /'epsilon','upsilon','Upsilon','diamond','ddagger',
     &          'partial','uparrow','Uparrow','right@}','lim inf',
     &          'lim sup','mathop{','version','hspace{','vspace{'/
      data ((isch(i,j),j=1,ic07),i=7,7)
     &         /'145',    '165',    '125',    '340',    '263',
     &          '266',    '255',    '335',    '175',    '000',
     &          '000',    '000',    '000',    '000',    '000'/
      data ((isfn(i,j),j=1,ic07),i=7,7)
     &         /13,       13,        4,        4,       -1,
     &          4,         4,        4,        4,       -4,
     &         -4,        -1,       -1,       -1,       -1/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic08),i=8,8)
     &         /'vartheta','varsigma','setminus','subseteq','supseteq',
     &          'emptyset','clubsuit','register','widehat{'/
      data ((isch(i,j),j=1,ic08),i=8,8)
     &         /'112',     '126',     '134',     '315',     '312',
     &          '306',     '247',     '322',     '303'/
      data ((isfn(i,j),j=1,ic08),i=8,8)
     &         /13,        13,        -1,         4,         4,
     &          4,          4,         4,        -2/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic09),i=9,9)
     &         /'leftarrow','Leftarrow','heartsuit','spadesuit',
     &          'copyright','trademark','backslash','downarrow',
     &          'Downarrow','overline{'/
      data ((isch(i,j),j=1,ic09),i=9,9)
     &         /'254',      '334',      '251',      '252',
     &          '323',      '324',      '134',      '257',
     &          '337',      '055'/
      data ((isfn(i,j),j=1,ic09),i=9,9)
     &         /4,           4,          4,          4,
     &          4,           4,         -1,          4,
     &          4,           4/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic010),i=10,10)
     &         /'varepsilon','rightarrow','Rightarrow','widetilde{',
     &          'underline{'/
      data ((isch(i,j),j=1,ic010),i=10,10)
     &         /'145',       '256',       '336',       '304',
     &          '055'/
      data ((isfn(i,j),j=1,ic010),i=10,10)
     &         /13,           4,           4,          -2,
     &          4/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic011),i=11,11)
     &         /'diamondsuit'/
      data ((isch(i,j),j=1,ic011),i=11,11)
     &         /'250'/
      data ((isfn(i,j),j=1,ic011),i=11,11)
     &         /4/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic013),i=13,13)
     &         /'hookleftarrow'/
      data ((isch(i,j),j=1,ic013),i=13,13)
     &         /'277'/
      data ((isfn(i,j),j=1,ic013),i=13,13)
     &         /4/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic014),i=14,14)
     &         /'leftrightarrow','Leftrightarrow','overleftarrow{'/
      data ((isch(i,j),j=1,ic014),i=14,14)
     &         /'253',           '333',           '254'/
      data ((isfn(i,j),j=1,ic014),i=14,14)
     &         /4,                4,               4/

*-----------------------------------------------------------------------

      data ((isex(i,j),j=1,ic015),i=15,15)
     &         /'overrightarrow{'/
      data ((isch(i,j),j=1,ic015),i=15,15)
     &         /'256'/
      data ((isfn(i,j),j=1,ic015),i=15,15)
     &         /4/

*-----------------------------------------------------------------------

      data nic/ ic01, ic02, ic03, ic04, ic05, ic06, ic07, ic08,
     &          ic09, ic010,ic011,ic012,ic013,ic014,ic015/

*-----------------------------------------------------------------------

      character bigc(3)*3
      save ipfn

*-----------------------------------------------------------------------

      character yen*1

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------

            yen = char(92)

            isex(5,28)(4:4) = yen
            isex(5,29)(4:4) = yen
            isex(5,30)(4:4) = yen
            isex(5,31)(4:4) = yen
            isex(6,29)(5:5) = yen
            isex(7, 9)(6:6) = yen

*-----------------------------------------------------------------------

      if( chaf(k) .ne. yen ) then

            chars6 = .false.

      else

*-----------------------------------------------------------------------

            chars6 = .true.

            ichl = 0
            ichn = 0

*-----------------------------------------------------------------------
*     FOR SPECIAL CHARACTER AFTER '\'
*-----------------------------------------------------------------------

               jfin = min( icn, k+icmax )

               dum20 = chaf(k+1)

*-----------------------------------------------------------------------

               jj = 1

            if( nic(jj) .gt. 0 ) then

               do 202 i = 1, nic(jj)

                  if( dum20(1:jj) .eq. isex(jj,i) ) then

                     ichl = jj
                     ichn = i

                  end if

  202          continue

            end if

*-----------------------------------------------------------------------

         do 100 j = k + 2, jfin

               jj = j - k

               dum20 = dum20(1:jj-1)//chaf(j)

            if( nic(jj) .gt. 0 ) then

               do 201 i = 1, nic(jj)

                  if( dum20(1:jj) .eq. isex(jj,i) ) then

                     ichl = jj
                     ichn = i

                  end if

  201          continue

            end if


  100    continue

*-----------------------------------------------------------------------
*        \123 or \12 IS EXCLUDED IN TEX EXPRESION
*        THESE ARE THE POSTSCRIPT CODE
*-----------------------------------------------------------------------

            if( ichl .eq. 1 .and.
     &          ichn .ge. 13 .and. ichn .le. 22 .and.
     &          numec(chaf(k+2)) ) then

                  ichl = 0

            end if

*-----------------------------------------------------------------------
*        \+STRING IS NOT INCLUDED IN TeX
*-----------------------------------------------------------------------

            if( ichl .eq. 0 ) then

               if( ifn(ig) .ne. iftdn ) then

                     if( ilt(ig) .gt. 0 ) then

                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                     else

                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                     end if

                        ilt(ig) = 0
                        ifn(ig) = iftdn

                        ifs(ig) = ipfs
                        icl(ig) = ipcl
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi

               end if

*-----------------------------------------------------------------------
*        SPECIAL CHARACTER FOR PS ( last \,  \\ \( \) \12 \123 )
*-----------------------------------------------------------------------

               if( k .eq. icn ) then

                  chag( ig, ilt(ig) + 1 ) = yen
                  chag( ig, ilt(ig) + 2 ) = yen

                  ilt(ig) = ilt(ig) + 2

               else if( chaf(k+1) .eq. yen .or.
     &                  chaf(k+1) .eq. '(' .or.
     &                  chaf(k+1) .eq. ')' ) then

                  chag( ig, ilt(ig) + 1 ) = chaf(k)
                  chag( ig, ilt(ig) + 2 ) = chaf(k+1)

                  k       = k + 1
                  ilt(ig) = ilt(ig) + 2

               else if( numec(chaf(k+1)) .and.
     &                  numec(chaf(k+2)) .and.
     &                  numec(chaf(k+3)) ) then

                  chag( ig, ilt(ig) + 1 ) = chaf(k)
                  chag( ig, ilt(ig) + 2 ) = chaf(k+1)
                  chag( ig, ilt(ig) + 3 ) = chaf(k+2)
                  chag( ig, ilt(ig) + 4 ) = chaf(k+3)

                  k       = k + 3
                  ilt(ig) = ilt(ig) + 4

               else if( numec(chaf(k+1)) .and.
     &                  numec(chaf(k+2)) ) then

                  chag( ig, ilt(ig) + 1 ) = chaf(k)
                  chag( ig, ilt(ig) + 2 ) = chaf(k+1)
                  chag( ig, ilt(ig) + 3 ) = chaf(k+2)

                  k       = k + 2
                  ilt(ig) = ilt(ig) + 3

*-----------------------------------------------------------------------
*           \+STRING IS NOT SUPPORTED
*-----------------------------------------------------------------------

               else

                  chag( ig, ilt(ig) + 1 ) = yen
                  chag( ig, ilt(ig) + 2 ) = yen

                  ilt(ig) = ilt(ig) + 2

               end if

*-----------------------------------------------------------------------
*           FOR SPECIAL CHARACTERS
*-----------------------------------------------------------------------

            else if( ichl .gt. 0 ) then

*-----------------------------------------------------------------------
*           DETERMINE THE FONT
*-----------------------------------------------------------------------

               jsfn(ichl,ichn) = isfn(ichl,ichn)

               if( isfn(ichl,ichn) .lt. 0 ) then

                  if( ifn(ig) .eq. 4 .or. ifn(ig) .eq. 13 ) then

                     if( ( ifn(ig) .eq.  4 .and. iftdn .eq.  4 ) .or.
     &                   ( ifn(ig) .eq. 13 .and. iftdn .eq. 13 ) ) then

                        jsfn(ichl,ichn) = 0

                     else

                        jsfn(ichl,ichn) = iftdn

                     end if

                  else

                        jsfn(ichl,ichn) = ifn(ig)

                  end if

               end if


               if( isfn(ichl,ichn) .le. -3 ) then

                  if( jsfn(ichl,ichn) .eq. 1 .or.
     &                jsfn(ichl,ichn) .eq. 2 .or.
     &                jsfn(ichl,ichn) .eq. 3 ) then

                        jsfn(ichl,ichn) = 0

                  else if( jsfn(ichl,ichn) .eq. 6 .or.
     &                     jsfn(ichl,ichn) .eq. 7 .or.
     &                     jsfn(ichl,ichn) .eq. 8 ) then

                        jsfn(ichl,ichn) = 5

                  else if( jsfn(ichl,ichn) .eq. 10 .or.
     &                     jsfn(ichl,ichn) .eq. 11 .or.
     &                     jsfn(ichl,ichn) .eq. 12 ) then

                        jsfn(ichl,ichn) = 9

                  end if

               end if

*-----------------------------------------------------------------------
*     ACCENT GROUP OF \"A TYPE AND \hbar
*-----------------------------------------------------------------------

               if( ( ichl .eq. 1 .and. ichn .eq. 28 ) .or.
     &             ( ichl .eq. 1 .and. ichn .eq. 29 ) .or.
     &             ( ichl .eq. 1 .and. ichn .eq. 30 ) .or.
     &             ( ichl .eq. 1 .and. ichn .eq. 31 ) .or.
     &             ( ichl .eq. 1 .and. ichn .eq. 32 ) .or.
     &             ( ichl .eq. 1 .and. ichn .eq. 33 ) .or.
     &             ( ichl .eq. 4 .and. ichn .eq. 37 ) ) then

                        isiu = 2
                        hbar = 0

                     if( ichl .eq. 4 .and. ichn .eq. 37 ) then

                        isiu = 4
                        hbar = 1

                     end if

                        iah(  idv + 1, 3 ) = jsfn(ichl,ichn)

                        iacd( idv + 1 )    = isch(ichl,ichn)

*-----------------------------------------------------------------------

                     if( ilt(ig) .gt. 0 ) then

                        ipfn = ibf(ib)

                     else

                        ipfn = ifn(ig)

                     end if

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

*-----------------------------------------------------------------------

                        idv = idv + 1

                        iah(idv,1) = ig - 1
                        iah(idv,2) = ids
                        iah(idv,4) = 0

                        call write00(8,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------

                        ilt(ig) = 1
                        ifn(ig) = iah(idv,3)

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

                        chag(ig,1) = chaf(k+ichl+1)

                        if( hbar .eq. 1 ) chag(ig,1) = 'h'

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------

                        ig = ig + 1

                        ilt(ig) = 4
                        ifn(ig) = iah(idv,3)

                        ifs(ig) = ig - 2
                        icl(ig) = ig - 2
                        ixp(ig) = ig
                        iyp(ig) = ig - 2
                        isi(ig) = isiu

                        chag(ig,1) = yen
                        chag(ig,2) = isch(ichl,ichn)(1:1)
                        chag(ig,3) = isch(ichl,ichn)(2:2)
                        chag(ig,4) = isch(ichl,ichn)(3:3)

                        idh     = 1

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 3
                        icl(ig) = ig - 3
                        ixp(ig) = ig - 2
                        iyp(ig) = ig - 3

                        idh     = 0

                        isi(ig) = isiu

                        call write00(9,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        idh     = 2

                        isi(ig) = 0

                        iah(idv,1) = 0
                        iah(idv,2) = 0
                        iah(idv,3) = 0
                        iah(idv,4) = 0
                        iah(idv,5) = 0

                        idv = idv - 1

                        k   = k + ichl

                        if( hbar .eq. 0 .or.
     &                    ( hbar .ne. 0 .and. chaf(k+1) .eq. ' ' ) )
     &                      k = k + 1


*-----------------------------------------------------------------------
*     ACCENT GROUP OF \***{***} TYPE
*-----------------------------------------------------------------------

               else if( ( ichl .eq. 4 .and. ichn .eq. 14 ) .or.
     &                  ( ichl .eq. 6 .and. ichn .eq. 19 ) .or.
     &                  ( ichl .eq. 6 .and. ichn .eq. 20 ) .or.
     &                  ( ichl .eq. 6 .and. ichn .eq. 21 ) .or.
     &                  ( ichl .eq. 6 .and. ichn .eq. 22 ) .or.
     &                  ( ichl .eq. 6 .and. ichn .eq. 23 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 15 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 16 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 17 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 38 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 39 ) .or.
     &                  ( ichl .eq. 9 .and. ichn .eq. 10 ) .or.
     &                  ( ichl .eq. 8 .and. ichn .eq.  9 ) .or.
     &                  ( ichl .eq.10 .and. ichn .eq.  4 ) .or.
     &                  ( ichl .eq.10 .and. ichn .eq.  5 ) .or.
     &                  ( ichl .eq.15 .and. ichn .eq.  1 ) .or.
     &                  ( ichl .eq.14 .and. ichn .eq.  3 ) .or.
     &                  ( ichl .eq. 5 .and. ichn .eq. 27 ) ) then


                        iah(  idv + 1, 3 ) = jsfn(ichl,ichn)

                        iacd( idv + 1 )    = isch(ichl,ichn)

*-----------------------------------------------------------------------

                     if( ilt(ig) .gt. 0 ) then

                        ipfn = ibf(ib)

                     else

                        ipfn = ifn(ig)

                     end if

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

*-----------------------------------------------------------------------

                        k  = k + ichl - 1

                        igs(ib+1)  = ig - 1

                        idv        = idv + 1

                        iah(idv,1) = ig - 1
                        iah(idv,2) = ids

                     if( ichl .eq. 4 .and. ichn .eq. 17 ) then

                        iah(idv,4) = 1

                     else if( ichl .eq. 4 .and. ichn .eq. 38 ) then

                        iah(idv,4) = -1

                     else if( ichl .eq. 4 .and. ichn .eq. 39 ) then

                        iah(idv,4) = -2

                     else if( ( ichl .eq. 9 .and. ichn .eq. 10 ) .or.
     &                        ( ichl .eq. 8 .and. ichn .eq.  9 ) .or.
     &                        ( ichl .eq.10 .and. ichn .eq.  4 ) ) then

                        iah(idv,4) = 2

                     else if( ( ichl .eq.10 .and. ichn .eq.  5 ) ) then

                        iah(idv,4) = 3

                     else if( ( ichl .eq.15 .and. ichn .eq.  1 ) ) then

                        iah(idv,4) = 4

                     else if( ( ichl .eq.14 .and. ichn .eq.  3 ) ) then

                        iah(idv,4) = 5

                     else

                        iah(idv,4) = 0

                     end if


*-----------------------------------------------------------------------

                        call write00(8,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0



*-----------------------------------------------------------------------
*     BIG CHARACTER CASE
*-----------------------------------------------------------------------

               else if( ( ichl .eq. 6 .and. ichn .eq. 24 ) .or.
     &                  ( ichl .eq. 6 .and. ichn .eq. 25 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 18 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 19 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 20 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 21 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 22 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 23 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 24 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 25 ) .or.
     &                  ( ichl .eq. 5 .and. ichn .eq. 28 ) .or.
     &                  ( ichl .eq. 5 .and. ichn .eq. 29 ) .or.
     &                  ( ichl .eq. 5 .and. ichn .eq. 30 ) .or.
     &                  ( ichl .eq. 5 .and. ichn .eq. 31 ) ) then

                        bigc(1)(1:1) = isch(ichl,ichn)(1:1)
                        bigc(1)(2:2) = isch(ichl,ichn)(2:2)
                        bigc(1)(3:3) = isch(ichl,ichn)(3:3)

                  if( ( ichl .eq. 6 .and. ichn .eq. 24 ) .or.
     &                ( ichl .eq. 6 .and. ichn .eq. 25 ) ) then

                        bigc(2)(1:1) = '3'
                        bigc(2)(2:2) = '6'
                        bigc(2)(3:3) = '3'

                        bigc(3)(1:1) = '3'
                        bigc(3)(2:2) = '6'
                        bigc(3)(3:3) = '5'

                  else if( ( ichl .eq. 4 .and. ichn .eq. 18 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 19 ) ) then

                        bigc(2)(1:1) = '3'
                        bigc(2)(2:2) = '5'
                        bigc(2)(3:3) = '1'

                        bigc(3)(1:1) = '3'
                        bigc(3)(2:2) = '5'
                        bigc(3)(3:3) = '3'

                  else if( ( ichl .eq. 4 .and. ichn .eq. 20 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 21 ) ) then

                        bigc(2)(1:1) = '3'
                        bigc(2)(2:2) = '7'
                        bigc(2)(3:3) = '1'

                        bigc(3)(1:1) = '3'
                        bigc(3)(2:2) = '7'
                        bigc(3)(3:3) = '3'

                  else if( ( ichl .eq. 4 .and. ichn .eq. 22 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 23 ) ) then

                        bigc(2)(1:1) = '3'
                        bigc(2)(2:2) = '4'
                        bigc(2)(3:3) = '6'

                        bigc(3)(1:1) = '3'
                        bigc(3)(2:2) = '5'
                        bigc(3)(3:3) = '0'

                  else if( ( ichl .eq. 4 .and. ichn .eq. 24 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 25 ) ) then

                        bigc(2)(1:1) = '3'
                        bigc(2)(2:2) = '6'
                        bigc(2)(3:3) = '6'

                        bigc(3)(1:1) = '3'
                        bigc(3)(2:2) = '7'
                        bigc(3)(3:3) = '0'

                  else if( ( ichl .eq. 5 .and. ichn .eq. 28 ) .or.
     &                     ( ichl .eq. 5 .and. ichn .eq. 29 ) ) then

                        bigc(2)(1:1) = '3'
                        bigc(2)(2:2) = '5'
                        bigc(2)(3:3) = '4'

                        bigc(3)(1:1) = '3'
                        bigc(3)(2:2) = '5'
                        bigc(3)(3:3) = '6'

                  else if( ( ichl .eq. 5 .and. ichn .eq. 30 ) .or.
     &                     ( ichl .eq. 5 .and. ichn .eq. 31 ) ) then

                        bigc(2)(1:1) = '3'
                        bigc(2)(2:2) = '7'
                        bigc(2)(3:3) = '4'

                        bigc(3)(1:1) = '3'
                        bigc(3)(2:2) = '7'
                        bigc(3)(3:3) = '6'

                  end if

                        ibig = 0

                  if(   ichl .eq. 6 .and. ichn .eq. 25   ) then

                        ibig = 3

                  end if

                  if( ( ichl .eq. 4 .and. ichn .eq. 19 ) .or.
     &                ( ichl .eq. 4 .and. ichn .eq. 21 ) .or.
     &                ( ichl .eq. 4 .and. ichn .eq. 23 ) .or.
     &                ( ichl .eq. 4 .and. ichn .eq. 25 ) .or.
     &                ( ichl .eq. 5 .and. ichn .eq. 29 ) .or.
     &                ( ichl .eq. 5 .and. ichn .eq. 31 ) ) then

                        ibig = 6

                  end if


*-----------------------------------------------------------------------

                     if( ilt(ig) .gt. 0 ) then

                        ipfn = ibf(ib)

                     else

                        ipfn = ifn(ig)

                     end if

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

*-----------------------------------------------------------------------

                     k = k + ichl

                     if( alphb(isex(ichl,ichn)(ichl:ichl)) .and.
     &                   chaf(k+1) .eq. ' ' ) k = k + 1

*-----------------------------------------------------------------------

                        idv = idv + 1

                        iah(idv,1) = ig - 1
                        iah(idv,2) = ids
                        iah(idv,3) = ig - 1
                        iah(idv,4) = 0

                        call write00(8,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------

                        idh = 3

                     do 500 kk = 1, 3

                        ilt(ig) = 4
                        ifn(ig) = jsfn(ichl,ichn)

                        ifs(ig) = ig - kk
                        icl(ig) = ig - kk
                        ixp(ig) = ig - kk + 1
                        iyp(ig) = ig - kk
                        isi(ig) = kk + 5 + ibig

                        chag(ig,1) = yen
                        chag(ig,2) = bigc(kk)(1:1)
                        chag(ig,3) = bigc(kk)(2:2)
                        chag(ig,4) = bigc(kk)(3:3)

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

  500                continue

*-----------------------------------------------------------------------

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 4
                        icl(ig) = ig - 4
                        ixp(ig) = ig - 3
                        iyp(ig) = ig - 4
                        isi(ig) = 0

                        idh = 0

                  if(      ( ichl .eq. 6 .and. ichn .eq. 24 ) .or.
     &                     ( ichl .eq. 6 .and. ichn .eq. 25 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 18 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 19 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 22 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 23 ) .or.
     &                     ( ichl .eq. 5 .and. ichn .eq. 28 ) .or.
     &                     ( ichl .eq. 5 .and. ichn .eq. 29 ) .or.
     &                     ( ichl .eq. 5 .and. ichn .eq. 30 ) .or.
     &                     ( ichl .eq. 5 .and. ichn .eq. 31 ) ) then

                        call write00(11,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                  else if( ( ichl .eq. 4 .and. ichn .eq. 20 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 21 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 24 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 25 ) ) then

                        call write00(12,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                  end if

*-----------------------------------------------------------------------

                        idh = 2


                        iah(idv,1) = 0
                        iah(idv,2) = 0
                        iah(idv,3) = 0
                        iah(idv,4) = 0
                        iah(idv,5) = 0

                        idv = idv - 1

                  if( ( ichl .eq. 6 .and. ichn .eq. 24 ) .or.
     &                ( ichl .eq. 6 .and. ichn .eq. 25 ) ) then

                        ibt(ig+1,1) = 1
                        ibt(ig+1,2) = ig

                  else if( ( ichl .eq. 4 .and. ichn .eq. 18 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 19 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 20 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 21 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 22 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 23 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 24 ) .or.
     &                     ( ichl .eq. 4 .and. ichn .eq. 25 ) .or.
     &                     ( ichl .eq. 5 .and. ichn .eq. 28 ) .or.
     &                     ( ichl .eq. 5 .and. ichn .eq. 29 ) .or.
     &                     ( ichl .eq. 5 .and. ichn .eq. 30 ) .or.
     &                     ( ichl .eq. 5 .and. ichn .eq. 31 ) ) then

                        ibt(ig+1,1) = 2
                        ibt(ig+1,2) = ig

                  end if


*-----------------------------------------------------------------------
*     \MATHOP{ } AND \KS{ }
*-----------------------------------------------------------------------

               else if( ( ichl .eq. 7 .and. ichn .eq. 12 ) .or.
     &                  ( ichl .eq. 3 .and. ichn .eq. 44 ) ) then


                        iah(  idv + 1, 3 ) = jsfn(ichl,ichn)

*-----------------------------------------------------------------------

                     if( ilt(ig) .gt. 0 ) then

                        ipfn = ibf(ib)

                     else

                        ipfn = ifn(ig)

                     end if

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

*-----------------------------------------------------------------------

                        k  = k + ichl - 1

                        igs(ib+1)  = ig - 1

                        idv        = idv + 1

                        iah(idv,1) = ig - 1
                        iah(idv,2) = ids

                     if( ichl .eq. 7 .and. ichn .eq. 12 ) then

                        iah(idv,4) = 60

                     else if( ichl .eq. 3 .and. ichn .eq. 44 ) then

                        iah(idv,4) = 61

                     end if

*-----------------------------------------------------------------------

                        call write00(8,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0



*-----------------------------------------------------------------------
*     \sum, \prod and \lim cases
*-----------------------------------------------------------------------

               else if( ( ( ichl .eq. 3 .and. ichn .eq. 20 ) .or.
     &                    ( ichl .eq. 3 .and. ichn .eq. 24 ) .or.
     &                    ( ichl .eq. 4 .and. ichn .eq. 12 ) ) .and.
     &                      isuf(k+ichl+1,1) .ne. 0 ) then

                     if( ilt(ig) .gt. 0 ) then

                        ipfn = ibf(ib)

                     else

                        ipfn = ifn(ig)

                     end if

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

*-----------------------------------------------------------------------

                        idv = idv + 1

                        iah(idv,1) = ig - 1
                        iah(idv,2) = ids
                        iah(idv,3) = ig - 1

                        iah(idv,4) = 6

                        call write00(8,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------

                        ilt(ig) = 0
                        ifn(ig) = jsfn(ichl,ichn)

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig -1
                        isi(ig) = 0

                     chag( ig, ilt(ig) + 1 ) = yen
                     chag( ig, ilt(ig) + 2 ) = isch(ichl,ichn)(1:1)
                     chag( ig, ilt(ig) + 3 ) = isch(ichl,ichn)(2:2)
                     chag( ig, ilt(ig) + 4 ) = isch(ichl,ichn)(3:3)

                     ilt(ig) = ilt(ig) + 4

                     k = k + ichl

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0


                        ibt(ig,1) = 3
                        ibt(ig,2) = ig - 2


*-----------------------------------------------------------------------
*     \frac and \fracs
*-----------------------------------------------------------------------

               else if( ( ichl .eq. 5 .and. ichn .eq. 33 ) .or.
     &                  ( ichl .eq. 6 .and. ichn .eq. 26 ) ) then


                        iah(  idv + 1, 3 ) = jsfn(ichl,ichn)

                        iacd( idv + 1 )    = isch(ichl,ichn)

*-----------------------------------------------------------------------

                     if( ilt(ig) .gt. 0 ) then

                        ipfn = ibf(ib)

                     else

                        ipfn = ifn(ig)

                     end if

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

*-----------------------------------------------------------------------

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1

                     if( ichl .eq. 5 .and. ichn .eq. 33 ) then

                        isi(ig) = 0

                     else if( ( ichl .eq. 6 .and. ichn .eq. 26 ) ) then

                        isi(ig) = 21

                     end if

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------

                        k  = k + ichl - 1

                        igs(ib+1)  = ig

                        idv        = idv + 1

                        iah(idv,1) = ig
                        iah(idv,2) = ids

                     if( ichl .eq. 5 .and. ichn .eq. 33 ) then

                        iah(idv,4) = 7

                     else if( ( ichl .eq. 6 .and. ichn .eq. 26 ) ) then

                        iah(idv,4) = 8

                     end if


*-----------------------------------------------------------------------

                        call write00(8,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0


*-----------------------------------------------------------------------
*     SQUARE ROOT
*-----------------------------------------------------------------------

               else if( ichl .eq. 5 .and. ichn .eq. 34 ) then


                        iah(  idv + 1, 3 ) = jsfn(ichl,ichn)

*-----------------------------------------------------------------------

                     if( ilt(ig) .gt. 0 ) then

                        ipfn = ibf(ib)

                     else

                        ipfn = ifn(ig)

                     end if

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

*-----------------------------------------------------------------------

                        k  = k + ichl - 1

                        ilt(ig) = 4
                        ifn(ig) = jsfn(ichl,ichn)

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

                        chag( ig, 1 ) = yen
                        chag( ig, 2 ) = isch(ichl,ichn)(1:1)
                        chag( ig, 3 ) = isch(ichl,ichn)(2:2)
                        chag( ig, 4 ) = isch(ichl,ichn)(3:3)

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        ig = ig + 1

                        ilt(ig) = 4
                        ifn(ig) = jsfn(ichl,ichn)

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig - 1
                        iyp(ig) = ig - 1
                        isi(ig) = 0

                        iwt(ig) = 1

                        chag( ig, 1 ) = yen
                        chag( ig, 2 ) = '1'
                        chag( ig, 3 ) = '4'
                        chag( ig, 4 ) = '0'

                        idh = 3

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        idh = 0

                        ig = ig + 1

                        igs(ib+1)  = ig - 3

                        idv        = idv + 1

                        iah(idv,1) = ig - 3
                        iah(idv,2) = ids

                        iah(idv,4) = 9

*-----------------------------------------------------------------------

                        call write00(18,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        call write00(8,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 3
                        icl(ig) = ig - 3
                        ixp(ig) = ig - 1
                        iyp(ig) = ig - 3
                        isi(ig) = 0



*-----------------------------------------------------------------------
*     \vspace{ }
*-----------------------------------------------------------------------

               else if( ichl .eq. 7 .and. ichn .eq. 15 ) then

                        icin = k + ichl

                        call pnum(chaf,icin,icolm,phv,ierr)

                        if( ierr .ne. 0 ) then

                           m_err = 'Number of '//yen//
     &                             'vspace{ } is wrong in '//
     &                             ' Mult Comment WT: or Tablur WTAB:'
                           ErrCha = ''
                           ErrID = 'L:6559/R:chars6/F:a-wtext.f'

                           goto 999

                        end if

                        k = icin

                     if( imull .gt. 0 ) then

                        ypsi(imull) = ypsi(imull) - fts(0) * phv

                        strhall = strhall + fts(0) * phv

                     end if

*-----------------------------------------------------------------------
*     \quad, \qquad, \,, \>, \;, \!, and \hspace{ }
*-----------------------------------------------------------------------

               else if( ( ichl .eq. 4 .and. ichn .eq. 28 ) .or.
     &                  ( ichl .eq. 5 .and. ichn .eq. 38 ) .or.
     &                  ( ichl .eq. 1 .and. ichn .eq. 23 ) .or.
     &                  ( ichl .eq. 1 .and. ichn .eq. 24 ) .or.
     &                  ( ichl .eq. 1 .and. ichn .eq. 25 ) .or.
     &                  ( ichl .eq. 1 .and. ichn .eq. 26 ) .or.
     &                  ( ichl .eq. 7 .and. ichn .eq. 14 ) ) then

                        ipfn = ifn(ig)

                     if( ilt(ig). gt. 0 ) then

                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                     else

                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                     end if

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ipfs
                        icl(ig) = ipcl
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi

                  if( ichl .eq. 4 .and. ichn .eq. 28 ) then

                        idh = -1

                  else if( ichl .eq. 5 .and. ichn .eq. 38 ) then

                        idh = -2

                  else if( ichl .eq. 1 .and. ichn .eq. 23 ) then

                        idh = -3

                  else if( ichl .eq. 1 .and. ichn .eq. 24 ) then

                        idh = -4

                  else if( ichl .eq. 1 .and. ichn .eq. 25 ) then

                        idh = -5

                  else if( ichl .eq. 1 .and. ichn .eq. 26 ) then

                        idh = -6

                  else if( ichl .eq. 7 .and. ichn .eq. 14 ) then

                        idh = -7

                        icin = k + ichl

                        call pnum(chaf,icin,icolm,phs,ierr)

                        if( ierr .ne. 0 ) then

                           m_err = 'Number of '//yen//
     &                             'hspace{ } is wrong in '//
     &                             ' Mult Comment WT:'
                           ErrCha = ''
                           ErrID = 'L:6660/R:chars6/F:a-wtext.f'

                           goto 999

                        end if

                        ichl = icin - k

                  end if


                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        idh = 0

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0


                        k = k + ichl

                     if( alphb(isex(ichl,ichn)(ichl:ichl)) .and.
     &                   chaf(k+1) .eq. ' ' ) k = k + 1

*-----------------------------------------------------------------------
*     NORMAL CASE FOR sin, cos AND so on
*-----------------------------------------------------------------------

               else if( isfn(ichl,ichn) .eq. -4 ) then

                        ipfn = ifn(ig)

                  if( jsfn(ichl,ichn) .ne. ifn(ig) ) then

                     if( ilt(ig) .gt. 0 ) then

                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                     else

                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                     end if

                        ilt(ig) = 0
                        ifn(ig) = jsfn(ichl,ichn)

                        ifs(ig) = ipfs
                        icl(ig) = ipcl
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi

                  end if

                  do 300 j = 1, ichl

                     chag( ig, ilt(ig) + j ) = isex(ichl,ichn)(j:j)

  300             continue

                     ilt(ig) = ilt(ig) + ichl

                        k = k + ichl

                     if( alphb(isex(ichl,ichn)(ichl:ichl)) .and.
     &                   chaf(k+1) .eq. ' ' ) k = k + 1

*-----------------------------------------------------------------------
*     today, time, page, version and file cases
*-----------------------------------------------------------------------

               else if( ( ichl .eq. 5 .and. ichn .eq. 41 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 34 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 35 ) .or.
     &                  ( ichl .eq. 4 .and. ichn .eq. 36 ) .or.
     &                  ( ichl .eq. 7 .and. ichn .eq. 13 ) ) then

                        ipfn = ifn(ig)

                  if( jsfn(ichl,ichn) .ne. ifn(ig) ) then

                     if( ilt(ig) .gt. 0 ) then

                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                     else

                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                     end if

                        ilt(ig) = 0
                        ifn(ig) = jsfn(ichl,ichn)

                        ifs(ig) = ipfs
                        icl(ig) = ipcl
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi

                  end if


                  if( ichl .eq. 5 .and. ichn .eq. 41 ) then

                        istl = jtoday

                     do 200 ist = 1, istl

                        chag( ig, ilt(ig) + ist ) = today(ist:ist)

  200                continue

                  else if( ichl .eq. 4 .and. ichn .eq. 34 ) then

                        istl = jtime

                     do 210 ist = 1, istl

                        chag( ig, ilt(ig) + ist ) = timestr(ist:ist)

  210                continue

                  else if( ichl .eq. 4 .and. ichn .eq. 35 ) then

                        istl = jnm

                     do 220 ist = 1, istl

                        chag( ig, ilt(ig) + ist ) = infnm(ist)

  220                continue

                  else if( ichl .eq. 7 .and. ichn .eq. 13 ) then

                        istl = javers

                     do 240 ist = 1, istl

                        chag( ig, ilt(ig) + ist ) = avers(ist:ist)

  240                continue

                  else if( ichl .eq. 4 .and. ichn .eq. 36 ) then

                        if( mpage .lt.  10 ) then
                           istl = 1
                        else if( mpage .lt. 100 ) then
                           istl = 2
                        else if( mpage .lt. 1000 ) then
                           istl = 3
                        else
                           istl = 4
                        end if

                           write(cdum,'(I4)') mpage

                     do 230 ist = 1, istl

                        chag( ig, ilt(ig) + ist ) =
     &                          cdum(4-istl+ist:4-istl+ist)

  230                continue

                  end if


                     ilt(ig) = ilt(ig) + istl

                     k = k + ichl

                     if( alphb(isex(ichl,ichn)(ichl:ichl)) .and.
     &                   chaf(k+1) .eq. ' ' ) k = k + 1

*-----------------------------------------------------------------------
*     NORMAL CASE
*-----------------------------------------------------------------------

               else

                        ipfn = ifn(ig)

                  if( jsfn(ichl,ichn) .ne. ifn(ig) ) then

                     if( ilt(ig) .gt. 0 ) then

                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                     else

                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                     end if

                        ilt(ig) = 0
                        ifn(ig) = jsfn(ichl,ichn)

                        ifs(ig) = ipfs
                        icl(ig) = ipcl
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi

                  end if

                     chag( ig, ilt(ig) + 1 ) = yen
                     chag( ig, ilt(ig) + 2 ) = isch(ichl,ichn)(1:1)
                     chag( ig, ilt(ig) + 3 ) = isch(ichl,ichn)(2:2)
                     chag( ig, ilt(ig) + 4 ) = isch(ichl,ichn)(3:3)

                     ilt(ig) = ilt(ig) + 4

                     k = k + ichl

                     if( alphb(isex(ichl,ichn)(ichl:ichl)) .and.
     &                   chaf(k+1) .eq. ' ' ) k = k + 1

*-----------------------------------------------------------------------

               end if

*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

  999       continue

            open(30,file='error.ang',status='UNKNOWN')

            write(30,'('' ***** Error Message in '',
     &                 ''ANGEL *****''/)')
            write( 6,'('' ***** Error Message in '',
     &                 ''ANGEL *****''/)')
            call ErrWrite(ErrID, ErrCha)
            write(30,'('' ERROR = '',a200)') m_err
            write( 6,'('' ERROR = '',a200)') m_err

            close(30)


      stop
      end

************************************************************************
*                                                                      *
      function chars7(jhf,idbg,k,chaf,icn,chag,ig,
     &                ib,iftdn,ibf,isuf,ibs,igs,ibb,ibc,
     &                ilt,ifn,ifs,ixp,iyp,isi,iap,ids,icl,clp,
     &                idh,iacd,idv,iah,iwt,ibt,ifon,imull,
     &                ipsp)
*                                                                      *
*              NEGLECT NON-SPECIAL '{' AND '}'                         *
*              OR END OF CHANGE FONT                                   *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chaf(ichrl)*1
      character chag(inig,0:ichrl)*1
      dimension ifn(0:inig), ifs(0:inig), ilt(0:inig), icl(0:inig)
      dimension clp(0:inig,3,3)
      dimension ixp(0:inig), iyp(0:inig), isi(0:inig), iap(0:inig,3)
      dimension ibf(0:inig), ibs(0:inig), igs(0:inig),
     &          ibb(0:inig), ibc(0:inig)
      dimension iwt(0:inig), ibt(0:inig,2), iah(0:inig,6)
      dimension isuf(0:ichrl,3)
      character iacd(0:inig)*3

      logical chars7
      save ipfn

*-----------------------------------------------------------------------

      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

         if( chaf(k) .eq. '{' .and. ipsp .ne. 1 ) then

                  ipfn    = ibf(ib)
                  ib      = ib + 1
                  ibf(ib) = ipfn
                  ibs(ib) = 0

                  chars7  = .true.

*-----------------------------------------------------------------------
*     END OF '}'
*-----------------------------------------------------------------------

         else if( chaf(k) .eq. '}' .and. ipsp .ne. 1 ) then

            if( ib .gt. 0 ) then

*-----------------------------------------------------------------------
*           END OF FIRST PART OF DOUBLE SUFIX
*-----------------------------------------------------------------------

               if( ( chaf(k+1) .eq. '^' .or. chaf(k+1) .eq. '_' ) .and.
     &             ( isuf(ibs(ib-1),1) .eq. 2 .or.
     &               isuf(ibs(ib-1),1) .eq. 4 ) .and.
     &             ( isuf(ibs(ib-1),3) .eq. k ) ) then


                     if( ibf(ib) .ne. ibf(ib-1) )
     &               iftdn = ibf(ib-1)

                     if( ilt(ig) .gt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                     end if


                  if( ibt(igs(ib-1),1) .eq. 3 .or.
     &                ibt(igs(ib-1),1) .eq. 4 ) then

                       call write00(5,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                     if( ibt(igs(ib-1),1) .eq. 3 ) then

                       call write00(7,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                     else if( ibt(igs(ib-1),1) .eq. 4 ) then

                       call write00(107,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                     end if

                  end if


                        ilt(ig) = 0
                        ifn(ig) = ifn(igs(ib-1)-1)

                        ifs(ig) = igs(ib-1) - 1
                        icl(ig) = igs(ib-1) - 1
                        ixp(ig) = igs(ib-1)
                        iyp(ig) = igs(ib-1) - 1


                  if( chaf(k+1) .eq. '^' ) then

                     if( ibt(igs(ib-1),1) .eq. 4 ) then

                        isi(ig) = 110

                     else

                        isi(ig) = 1

                     end if

                  else

                     if( ibt(igs(ib-1),1) .eq. 4 ) then

                        isi(ig) = -110

                     else

                        isi(ig) = -1

                     end if

                  end if


                        iap(ids,2) = ig

                     if( chaf(k+1) .eq. '^' ) then

                        iap(ids,1) =   igs(ib-1)

                     else

                        iap(ids,1) = - igs(ib-1)

                     end if


                     if( isuf(ibs(ib-1),1) .eq. 2 ) then

                        ibt(ig,1) = ibt(abs(iap(ids,1)),1)
                        ibt(ig,2) = ibt(abs(iap(ids,1)),2)

                     end if

                        call write00(24,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        k = k + 1

                     if( ibt(ig,1) .eq. 3 .or. ibt(ig,1) .eq. 4 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        iah(idv,1) = ig
                        iah(idv,2) = ids
cKN 2003/07/07  ipfs ->ipfn
                        iah(idv,3) = ipfn
                        iah(idv,4) = 6

                        call write00(8,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0


                     end if

*-----------------------------------------------------------------------
*         END OF DOUBLE SUFIX
*-----------------------------------------------------------------------

               else if( ( isuf(ibs(ib),1) .eq. 2 .or.
     &                    isuf(ibs(ib),1) .eq. 4 ) .and.
     &                    isuf(ibs(ib),2) .eq. k ) then

                  if( ilt(ig) .gt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                  end if

                     if( isuf(ibs(ib),1) .eq. 2 ) then

                        call write00(2,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                     else if( isuf(ibs(ib),1) .eq. 4 ) then

                        call write00(4,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                     end if

                  if( ibt(igs(ib), 1 ) .eq. 3 ) then

                       call write00(67,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                              ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                              idh,idv,iah,iwt,ibt,ifon,imull,phs)

                  else if( ibt(igs(ib), 1 ) .eq. 4 ) then
*
                       call write00(167,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                              ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                              idh,idv,iah,iwt,ibt,ifon,imull,phs)

                  end if


                     if( ilt(ig) .gt. 0 ) ig = ig + 1

                        iap(ids,1) = 0
                        iap(ids,2) = 0
                        iap(ids,3) = 0

                        ids = ids - 1

                        ilt(ig) = 0
                        ifn(ig) = ifn(igs(ib)-1)

                        ifs(ig) = igs(ib) - 1
                        icl(ig) = igs(ib) - 1
                        ixp(ig) = - igs(ib)
                        iyp(ig) = igs(ib) - 1
                        isi(ig) = 0


                  if( ibt(igs(ib), 1 ) .eq. 3 .or.
     &                ibt(igs(ib), 1 ) .eq. 4 ) then

                        idh = 0

                        call write00(69,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        iah(idv,1) = 0
                        iah(idv,2) = 0
                        iah(idv,3) = 0
                        iah(idv,4) = 0
                        iah(idv,5) = 0

                        idv = idv - 1

                        idh = 2

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        idh = 0

                        call write00(70,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ifn(ig-1)

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

                        iah(idv,1) = 0
                        iah(idv,2) = 0
                        iah(idv,3) = 0
                        iah(idv,4) = 0
                        iah(idv,5) = 0

                        idv = idv - 1

                  end if

*-----------------------------------------------------------------------
*        END OF SINGLE SUFIXE
*-----------------------------------------------------------------------

               else if( ( isuf(ibs(ib),1) .eq. 1 .or.
     &                    isuf(ibs(ib),1) .eq. 3 ) .and.
     &                    isuf(ibs(ib),2) .eq. k ) then

                     if( ilt(ig) .gt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                     end if

                  if( ibt(igs(ib), 1 ) .ge. 1 ) then

                        call write00(5,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                     if( ibt(igs(ib), 1 ) .eq. 3 ) then

                       call write00(6,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                              ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                              idh,idv,iah,iwt,ibt,ifon,imull,phs)

                     else if( ibt(igs(ib), 1 ) .eq. 4 ) then

                       call write00(106,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                              ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                              idh,idv,iah,iwt,ibt,ifon,imull,phs)

                     end if

                        iap(ids,1) = 0
                        iap(ids,2) = 0
                        iap(ids,3) = 0

                        ids = ids - 1

                  end if


                     if( ilt(ig) .gt. 0 ) ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ifn(igs(ib)-1)

                        ifs(ig) = igs(ib) - 1
                        icl(ig) = igs(ib) - 1
                        iyp(ig) = igs(ib) - 1
                        isi(ig) = 0

                        ixp(ig) = ig

                     if( ibt(igs(ib), 1 ) .ge. 1 ) then

                        ixp(ig) = igs(ib)

                     end if


                  if( ibt(igs(ib), 1 ) .eq. 3 .or.
     &                ibt(igs(ib), 1 ) .eq. 4 ) then

                        idh = 0

                        call write00(69,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        iah(idv,1) = 0
                        iah(idv,2) = 0
                        iah(idv,3) = 0
                        iah(idv,4) = 0
                        iah(idv,5) = 0

                        idv = idv - 1

                        idh = 2

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        idh = 0


                        call write00(70,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ifn(ig-1)

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

                        iah(idv,1) = 0
                        iah(idv,2) = 0
                        iah(idv,3) = 0
                        iah(idv,4) = 0
                        iah(idv,5) = 0

                        idv = idv - 1

                  end if

*-----------------------------------------------------------------------
*        END OF FIRST PART OF FRAC AND FRACS
*-----------------------------------------------------------------------

               else if( igs(ib) .ne. 0 .and.
     &                  igs(ib) .eq. iah(idv,1) .and.
     &                ( iah(idv,4) .eq. 7 .or.
     &                  iah(idv,4) .eq. 8 ) .and.
     &                  chaf(k+1) .eq. '{' ) then

                  if( ilt(ig) .gt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                           ig = ig + 1

                  else

                     if( ixp(ig) .lt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                           ig = ig + 1

                     end if

                  end if


                        call write00(90,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        ilt(ig) = 0
                        ifn(ig) = ifn(igs(ib))

                        if( ibf(ib) .ne. ibf(ib-1) )
     &                  iftdn = ibf(ib-1)

                        ifs(ig) = ifs(igs(ib))
                        icl(ig) = icl(igs(ib))
                        ixp(ig) = ixp(igs(ib))
                        iyp(ig) = iyp(igs(ib))

                     if( iah(idv,4) .eq. 7 ) then

                        isi(ig) = 0

                     else if( iah(idv,4) .eq. 8 ) then

                        isi(ig) = -21

                     end if

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        igs(ib)    = ig

                        iah(idv,5) = iah(idv,1)

                        iah(idv,1) = ig
                        iah(idv,2) = ids


                        call write00(8,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0


*-----------------------------------------------------------------------
*        END OF FRAC AND FRACS
*-----------------------------------------------------------------------

               else if( igs(ib) .ne. 0 .and.
     &                  igs(ib) .eq. iah(idv,1) .and.
     &                ( iah(idv,4) .eq. 7 .or.
     &                  iah(idv,4) .eq. 8 ) ) then

                  if( ilt(ig) .gt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                           ig = ig + 1

                  else

                     if( ixp(ig) .lt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                           ig = ig + 1

                     end if

                  end if


                        call write00(91,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        ilt(ig) = 4
                        ifn(ig) = 4

                        ifs(ig) = ifs(igs(ib))
                        icl(ig) = icl(igs(ib))
                        ixp(ig) = ixp(igs(ib))
                        iyp(ig) = iyp(igs(ib))
                        isi(ig) = 0

                        chag(ig,1) = yen
                        chag(ig,2) = '0'
                        chag(ig,3) = '5'
                        chag(ig,4) = '5'

                        iwt(ig)    = 1

                        idh = 1

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        idh = 0

                        call write00(92,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ifn(igs(ib))

                        ifs(ig) = iah(idv,5) - 1
                        icl(ig) = iah(idv,5) - 1
                        ixp(ig) = iah(idv,5)
                        iyp(ig) = iah(idv,5) - 1
                        isi(ig) = 0

                        idh = 2

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        idh = 0

                        iah(idv,1) = iah(idv,5) - 1

                        call write00(70,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ifn(igs(ib))

                        if( ibf(ib) .ne. ibf(ib-1) )
     &                  iftdn = ibf(ib-1)

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

                        iah(idv,1) = 0
                        iah(idv,2) = 0
                        iah(idv,3) = 0
                        iah(idv,4) = 0
                        iah(idv,5) = 0

                        idv     = idv - 1

                        igs(ib) = 0


*-----------------------------------------------------------------------
*        END OF SQRT
*-----------------------------------------------------------------------

               else if( igs(ib) .ne. 0 .and.
     &                  igs(ib) .eq. iah(idv,1) .and.
     &                  iah(idv,4) .eq. 9 ) then

                  if( ilt(ig) .gt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                           ig = ig + 1

                  else

                     if( ixp(ig) .lt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                           ig = ig + 1

                     end if

                  end if


                        call write00(80,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        call write00(70,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)


                        ilt(ig) = 0
                        ifn(ig) = ifn(iah(idv,1))

                        ifs(ig) = iah(idv,1)
                        icl(ig) = iah(idv,1)
                        ixp(ig) = iah(idv,1) + 1
                        iyp(ig) = iah(idv,1)
                        isi(ig) = 0

                        idh = 2

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        idh = 0

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ifn(iah(idv,1))

                        if( ibf(ib) .ne. ibf(ib-1) )
     &                  iftdn = ibf(ib-1)

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

                        iah(idv,1) = 0
                        iah(idv,2) = 0
                        iah(idv,3) = 0
                        iah(idv,4) = 0
                        iah(idv,5) = 0

                        idv     = idv - 1

                        igs(ib) = 0


*-----------------------------------------------------------------------
*        END OF \MATHOP{  }, \KS{  }
*-----------------------------------------------------------------------

               else if( igs(ib) .ne. 0 .and.
     &                  igs(ib) .eq. iah(idv,1) .and.
     &                ( iah(idv,4) .eq. 60 .or.
     &                  iah(idv,4) .eq. 61 ) ) then

                     if( iah(idv,4) .eq. 60 ) then

                        ipbt1 = 3

                     else if( iah(idv,4) .eq. 61 ) then

                        ipbt1 = 4

                     end if

                  if( ilt(ig) .gt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                           ig = ig + 1

                  else

                     if( ixp(ig) .lt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                           ig = ig + 1

                     end if

                  end if


                        ilt(ig) = 0
                        ifn(ig) = ifn(iah(idv,1))

                        if( ibf(ib) .ne. ibf(ib-1) )
     &                  iftdn = ibf(ib-1)

                        ifs(ig) = iah(idv,1)
                        icl(ig) = iah(idv,1)
                        ixp(ig) = iah(idv,1) + 1
                        iyp(ig) = iah(idv,1)
                        isi(ig) = 0


                  if( isuf(k+1,1) .eq. 0 ) then

                        call write00(70,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ixp(ig) = ig

                        iah(idv,1) = 0
                        iah(idv,2) = 0
                        iah(idv,3) = 0
                        iah(idv,4) = 0
                        iah(idv,5) = 0

                        idv = idv - 1

                  else if( isuf(k+1,1) .ne. 0 ) then

                        ibt(ig,1) = ipbt1
                        ibt(ig,2) = ig - 2

                  end if

*-----------------------------------------------------------------------
*        END OF ACCENT GROUP OF \***{***} TYPE
*-----------------------------------------------------------------------

               else if( igs(ib) .ne. 0 .and.
     &                  igs(ib) .eq. iah(idv,1) .and.
     &                  iah(idv,4) .le. 5 ) then

                  if( ilt(ig) .gt. 0 ) then

                        jpfn = ifn(ig)

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                           ig = ig + 1

                  else

                           jpfn = ifn(ig)

                     if( ixp(ig) .lt. 0 ) then

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                           ig = ig + 1

                     end if

                  end if

                              ipfn = jpfn

                        if( ipfn .eq. 4 ) then

                           if( iah(idv,3) .ne. 4 ) then

                              ipfn = 0

                           else if( iah(idv,3) .eq. 4 ) then

                              ipfn = 4

                           end if

                        else if( ipfn .ne. 4 ) then

                           if( iah(idv,3) .ne. 4 ) then

                              ipfn = 0

                           else if( iah(idv,3) .eq. 4 ) then

                              ipfn = 4

                           end if

                        end if

*-----------------------------------------------------------------------

                        ilt(ig) = 4
                        ifn(ig) = ipfn

                        ifs(ig) = igs(ib)
                        icl(ig) = igs(ib)
                        ixp(ig) = igs(ib) + 1
                        iyp(ig) = igs(ib)
                        isi(ig) = 2

                        chag(ig,1) = yen
                        chag(ig,2) = iacd(idv)(1:1)
                        chag(ig,3) = iacd(idv)(2:2)
                        chag(ig,4) = iacd(idv)(3:3)


                     if( iah(idv,4) .eq. 1 .or.
     &                   iah(idv,4) .eq. 4 .or.
     &                   iah(idv,4) .eq. 5 ) then

                        isi(ig) = 3

                     else if( iah(idv,4) .eq. -1 ) then

                        isi(ig) = 31

                     else if( iah(idv,4) .eq. -2 ) then

                        isi(ig) = 32

                     else if( iah(idv,4) .eq. 3 ) then

                        isi(ig) = 5

                     end if

                     if( iah(idv,4) .eq. 2 .or.
     &                   iah(idv,4) .eq. 3 ) then

                        iwt(ig) = 1

                     end if


                        idh     = 1

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ifn(igs(ib))

                        ifs(ig) = igs(ib)
                        icl(ig) = igs(ib)
                        ixp(ig) = igs(ib) + 1
                        iyp(ig) = igs(ib)
                        isi(ig) = 0

                        idh     = 0

                  if( iah(idv,4) .le. 1 ) then

                        call write00(9,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                  else if( iah(idv,4) .eq. 2 ) then

                        call write00(19,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                  else if( iah(idv,4) .eq. 3 ) then

                        call write00(29,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                  else if( iah(idv,4) .eq. 4 .or.
     &                     iah(idv,4) .eq. 5 ) then

                     if( iah(idv,4) .eq. 4 ) then

                        call write00(39,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                     else if( iah(idv,4) .eq. 5 ) then

                        call write00(49,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                     end if

                        ilt(ig) = 4
                        ifn(ig) = ifn(ig-1)

                        ifs(ig) = ig - 1
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

                        chag(ig,1) = yen
                        chag(ig,2) = '2'
                        chag(ig,3) = '7'
                        chag(ig,4) = '6'

                        iwt(ig) = 1

                     if( iah(idv,4) .eq. 4 ) then

                        call write00(40,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                     else if( iah(idv,4) .eq. 5 ) then

                        call write00(50,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                     end if


                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ifn(igs(ib))

                        ifs(ig) = igs(ib)
                        icl(ig) = igs(ib)
                        ixp(ig) = igs(ib) + 1
                        iyp(ig) = igs(ib)
                        isi(ig) = 0

                  end if

*-----------------------------------------------------------------------

                        idh     = 2

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------

                        ig = ig + 1

                        ilt(ig) = 0
                        ifn(ig) = ifn(igs(ib))

                        if( ibf(ib) .ne. ibf(ib-1) )
     &                  iftdn = ibf(ib-1)

                        ifs(ig) = igs(ib)
                        icl(ig) = ig - 1
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

                        idh     = 0

                        iah(idv,1) = 0
                        iah(idv,2) = 0
                        iah(idv,3) = 0
                        iah(idv,4) = 0
                        iah(idv,5) = 0

                        idv     = idv - 1

                        igs(ib) = 0

*-----------------------------------------------------------------------
*        END OF BOX
*-----------------------------------------------------------------------

               else if( igs(ib) .ne. 0 .and.
     &                  igs(ib) .eq. iah(idv,1) .and.
     &                  iah(idv,4) .ge. 10 ) then

                        iboxd = iah(idv,4)

                     if( ilt(ig). gt. 0 ) then

                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                     else

                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                     end if

*-----------------------------------------------------------------------

                        call write00(70,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------

*                       IAH(IDV,1) = 0 <- NOT NOW, WHICH IS USED BELOW

*-----------------------------------------------------------------------

                        iah(idv,2) = 0
                        iah(idv,3) = 0
                        iah(idv,4) = 0
                        iah(idv,5) = 0

                        idv     = idv - 1

                        igs(ib) = 0

*-----------------------------------------------------------------------

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        if( ibf(ib) .ne. ibf(ib-1) )
     &                  iftdn = ibf(ib-1)

                        ifs(ig) = ipfs
                        icl(ig) = ipcl
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi

*-----------------------------------------------------------------------

                     if( iboxd .eq. 40 .or. iboxd .eq. 41 .or.
     &                   iboxd .eq. 50 .or. iboxd .eq. 51 ) then

                              idh = -14

                     else if( iboxd .eq. 42 .or. iboxd .eq. 43 .or.
     &                        iboxd .eq. 52 .or. iboxd .eq. 53 ) then

                              idh = -15

                     else if( iboxd .eq. 10 .or. iboxd .eq. 11 .or.
     &                        iboxd .eq. 20 .or. iboxd .eq. 21 .or.
     &                        iboxd .eq. 30 .or. iboxd .eq. 31 ) then

                              idh = -12

                     else

                              idh = -13

                     end if

*-----------------------------------------------------------------------

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

*-----------------------------------------------------------------------
*                       HERE : IAH(IDV+1,1) = 0
*-----------------------------------------------------------------------

                        iah(idv+1,1) = 0

*-----------------------------------------------------------------------

                        idh = 0

                        ig = ig + 1

*-----------------------------------------------------------------------

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        ifs(ig) = ig - 1
                        icl(ig) = ibc(ib)
                        ixp(ig) = ig
                        iyp(ig) = ig - 1
                        isi(ig) = 0

*-----------------------------------------------------------------------
*     END OF FONT SIZE CHANGE
*-----------------------------------------------------------------------

               else if( ibb(ib) .ne. -1 ) then

                  if( ilt(ig) .gt. 0 ) then

                        ipfn = ifn(ig)
                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                  else

                        ipfn = ifn(ig)
                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                  end if

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        if( ibf(ib) .ne. ibf(ib-1) )
     &                  iftdn = ibf(ib-1)

                        ifs(ig) = ibb(ib)
                        icl(ig) = ipcl
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi

                        ibb(ib) = -1

*-----------------------------------------------------------------------
*     END OF CHARACTER COLOR CHANGE
*-----------------------------------------------------------------------

               else if( ibc(ib) .ne. -1 ) then

                  if( ilt(ig) .gt. 0 ) then

                        ipfn = ifn(ig)
                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                  else

                        ipfn = ifn(ig)
                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                  end if

                        ilt(ig) = 0
                        ifn(ig) = ipfn

                        if( ibf(ib) .ne. ibf(ib-1) )
     &                  iftdn = ibf(ib-1)

                        ifs(ig) = ipfs
                        icl(ig) = ibc(ib)
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi

                        ibc(ib) = -1

*-----------------------------------------------------------------------
*     END OF FONT NAME CHANGE
*-----------------------------------------------------------------------

               else if( ibf(ib) .ne. ibf(ib-1) ) then

                        iftdn = ibf(ib-1)

                     if( ilt(ig) .gt. 0 ) then

                        ipfs = ig
                        ipcl = ig
                        ipxp = ig + 1
                        ipyp = ig
                        ipsi = 0

                        call write00(0,idbg,jhf,ig,chag,igs,ib,icl,clp,
     &                               ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                               idh,idv,iah,iwt,ibt,ifon,imull,phs)

                        ig = ig + 1

                     else

                        ipfs = ifs(ig)
                        ipcl = icl(ig)
                        ipxp = ixp(ig)
                        ipyp = iyp(ig)
                        ipsi = isi(ig)

                     end if

                        ilt(ig) = 0
                        ifn(ig) = ibf(ib-1)

                        ifs(ig) = ipfs
                        icl(ig) = ipcl
                        ixp(ig) = ipxp
                        iyp(ig) = ipyp
                        isi(ig) = ipsi


*-----------------------------------------------------------------------

               end if


            end if

                        ib = ib - 1

                        if( ib .lt. 0 ) ib = 0

*-----------------------------------------------------------------------

                  chars7 = .true.

         else

                  chars7 = .false.

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine write00(icwt,idbg,jhf,i,chag,igs,ib,icl,clp,
     &                   ifn,ifs,ilt,ixp,iyp,isi,iap,ids,
     &                   idh,idv,iah,iwt,ibt,ifon,imull,phs)
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chag(inig,0:ichrl)*1
      dimension ifn(0:inig), ifs(0:inig), ilt(0:inig), icl(0:inig)
      dimension clp(0:inig,3,3)
      dimension ixp(0:inig), iyp(0:inig), isi(0:inig), iap(0:inig,3)
      dimension iwt(0:inig), ibt(0:inig,2), iah(0:inig,6), igs(0:inig)

      character chaf(ichrl)*1

      dimension bcb(3), bcl(3), bcs(3)

*-----------------------------------------------------------------------

      common /wtval1/ strl0, strh0, strb0
      common /wtval2/ clc(0:inig,3), fts(0:inig),
     &                stw(0:inig), sth(0:inig), stb(0:inig),
     &                xps(0:inig), yps(0:inig),
     &                xpp(0:inig), ypp(0:inig),
     &                ixs(0:inig), iys(0:inig), jfn(0:inig)
      common /wtval3/ strla(0:inig), strha(0:inig), strba(0:inig)
      common /wtval4/ stwh(0:inig),  stws(0:inig),
     &                stwu(0:inig),  stwd(0:inig),
     &                strlu(0:inig), strld(0:inig)
      common /wtval5/ rsid(0:inig),  rlng(0:inig),  rsca(0:inig)
      common /wtval6/ stwfu(0:inig), sthfu(0:inig), stbfu(0:inig),
     &                stwfd(0:inig), sthfd(0:inig), stbfd(0:inig),
     &                strhb(0:inig), strhh(0:inig), strbb(0:inig)
      common /wtval7/ strhtp, strbtp, strhbt, strbbt, strhall
      common /wtval8/ xpsi(0:inig), ypsi(0:inig), strl(0:inig),
     &                xpsl(0:inig), xpsr(0:inig), rind(0:inig)

*-----------------------------------------------------------------------

      common /box/  bxws, bxwl, bxss, bxsl, bxds, bxdl, bxls, bxll
      common /suf/  ftht, ftss, ftks, ftsv,
     &              yplu, ypld, ypku, ypkd, yplh, ypls
      common /con/  cm, dd
      common /pap/  a4w, a4h, wmg, hmg, wct, hct

*-----------------------------------------------------------------------

      common /stwrl/  stwrl4(3,4)

*-----------------------------------------------------------------------


         if( icwt .eq. 4 ) then

                     stws(abs(iap(ids,1))-1)
     &                                 = max( strlu(abs(iap(ids,1))-1),
     &                                        strld(abs(iap(ids,1))-1) )

                     stwu(abs(iap(ids,1))-1) = stws(abs(iap(ids,1))-1)
     &                                       - strlu(abs(iap(ids,1))-1)

                     stwd(abs(iap(ids,1))-1) = stws(abs(iap(ids,1))-1)
     &                                       - strld(abs(iap(ids,1))-1)

            if( iap(ids,1) .gt. 0 ) then

                     xps(abs(iap(ids,1))) = stw(abs(iap(ids,1))-1)
     &                                    + stwd(abs(iap(ids,1))-1)
                     ixs(abs(iap(ids,1))) = abs(iap(ids,1)) - 1
                     xpp(abs(iap(ids,1))) = xpp(ixs(abs(iap(ids,1))))
     &                                    + xps(abs(iap(ids,1)))


                     xps(iap(ids,2)) = stw(abs(iap(ids,1))-1)
     &                                    + stwu(abs(iap(ids,1))-1)
                     ixs(iap(ids,2)) = abs(iap(ids,1)) - 1
                     xpp(iap(ids,2)) = xpp(ixs(iap(ids,2)))
     &                               + xps(iap(ids,2))

            else if( iap(ids,1) .lt. 0 ) then

                     xps(abs(iap(ids,1))) = stw(abs(iap(ids,1))-1)
     &                                    + stwu(abs(iap(ids,1))-1)
                     ixs(abs(iap(ids,1))) = abs(iap(ids,1)) - 1
                     xpp(abs(iap(ids,1))) = xpp(ixs(abs(iap(ids,1))))
     &                                    + xps(abs(iap(ids,1)))

                     xps(iap(ids,2)) = stw(abs(iap(ids,1))-1)
     &                                    + stwd(abs(iap(ids,1))-1)
                     ixs(iap(ids,2)) = abs(iap(ids,1)) - 1
                     xpp(iap(ids,2)) = xpp(ixs(iap(ids,2)))
     &                               + xps(iap(ids,2))

            end if

               if( ids - 1 .le. iah(idv,2) ) then

                     strla(iah(idv,1)) = strla(iah(idv,1))
     &                                 + stws(abs(iap(ids,1))-1)

               else

                  if( iap(ids-1,1) .gt. 0 ) then

                     strlu(iap(ids-1,1)-1) = strlu(iap(ids-1,1)-1)
     &                                     + stws(abs(iap(ids,1))-1)

                  else if( iap(ids-1,1) .lt. 0 ) then

                     strld(-iap(ids-1,1)-1) = strld(-iap(ids-1,1)-1)
     &                                      + stws(abs(iap(ids,1))-1)

                  end if

               end if

         else if( icwt .eq.  6 .or. icwt .eq. 106 .or.
     &            icwt .eq.  7 .or. icwt .eq. 107 .or.
     &            icwt .eq. 67 .or. icwt .eq. 167 ) then

                  if( iap(ids,1) .gt. 0 ) then

                     strlu(abs(iap(ids,1))-1) = strla(iah(idv,1))

                     strla(iah(idv-1,1)) = max( strla(iah(idv-1,1)),
     &                                     strlu(abs(iap(ids,1))-1) )

                  else if( iap(ids,1) .lt. 0 ) then

                     strld(abs(iap(ids,1))-1) = strla(iah(idv,1))

                     strla(iah(idv-1,1)) = max( strla(iah(idv-1,1)),
     &                                     strld(abs(iap(ids,1))-1) )

                  end if

                  if( icwt .eq.  6 .or. icwt .eq. 106 .or.
     &                icwt .eq. 67 .or. icwt .eq. 167 ) then

                     xps(abs(iap(ids,1))-1) = stw(abs(iap(ids,1))-2)
     &                                    + ( strla(iah(idv-1,1))
     &                                      - stw(abs(iap(ids,1))-1) )
     &                                    / 2.0
                     ixs(abs(iap(ids,1))-1) = abs(iap(ids,1)) - 2
                     xpp(abs(iap(ids,1))-1) =
     &                                      xpp(ixs(abs(iap(ids,1))-1))
     &                                    + xps(abs(iap(ids,1))-1)

                  end if


               if( iap(ids,1) .gt. 0 ) then

                  if( icwt .eq. 6 .or. icwt .eq. 106 ) then

                     xps(abs(iap(ids,1))) = stw(abs(iap(ids,1))-2)
     &                                  + ( strla(iah(idv-1,1))
     &                                    - strlu(abs(iap(ids,1))-1) )
     &                                  / 2.0
                     ixs(abs(iap(ids,1))) = abs(iap(ids,1)) - 2
                     xpp(abs(iap(ids,1))) = xpp(ixs(abs(iap(ids,1))))
     &                                    + xps(abs(iap(ids,1)))

                  end if

                  if( icwt .eq. 67 .or. icwt .eq. 167 ) then

                     xps(abs(iap(ids,1))) = stw(abs(iap(ids,1))-2)
     &                                  + ( strla(iah(idv-1,1))
     &                                    - strld(abs(iap(ids,1))-1) )
     &                                  / 2.0
                     ixs(abs(iap(ids,1))) = abs(iap(ids,1)) - 2
                     xpp(abs(iap(ids,1))) = xpp(ixs(abs(iap(ids,1))))
     &                                    + xps(abs(iap(ids,1)))

                     xps(abs(iap(ids,2))) = stw(abs(iap(ids,1))-2)
     &                                  + ( strla(iah(idv-1,1))
     &                                    - strlu(abs(iap(ids,1))-1) )
     &                                  / 2.0
                     ixs(abs(iap(ids,2))) = abs(iap(ids,1)) - 2
                     xpp(abs(iap(ids,2))) = xpp(ixs(abs(iap(ids,2))))
     &                                    + xps(abs(iap(ids,2)))

                  end if

                  if( icwt .eq. 67 ) then

                     yps(abs(iap(ids,2))) = sth(abs(iap(ids,1))-1)
     &                                    + fts(abs(iap(ids,1))-1)
     &                                    * ypls

                     iys(abs(iap(ids,2))) = abs(iap(ids,1)) - 2

                     ypp(abs(iap(ids,2))) = ypp(iys(abs(iap(ids,2))))
     &                                    + yps(abs(iap(ids,2)))

                  end if

                  if( icwt .eq. 6 .or. icwt .eq. 7 ) then

                     yps(abs(iap(ids,1))) = sth(abs(iap(ids,1))-1)
     &                                    + fts(abs(iap(ids,1))-1)
     &                                    * ypls

                     iys(abs(iap(ids,1))) = abs(iap(ids,1)) - 2

                     ypp(abs(iap(ids,1))) = ypp(iys(abs(iap(ids,1))))
     &                                    + yps(abs(iap(ids,1)))

                  end if

               else if( iap(ids,1) .lt. 0 ) then

                  if( icwt .eq. 6 .or. icwt .eq. 106 ) then

                     xps(abs(iap(ids,1))) = stw(abs(iap(ids,1))-2)
     &                                    + ( strla(iah(idv-1,1))
     &                                      - strld(abs(iap(ids,1))-1) )
     &                                    / 2.0

                     ixs(abs(iap(ids,1))) = abs(iap(ids,1)) - 2

                  end if

                  if( icwt .eq. 67 .or. icwt .eq. 167 ) then

                     xps(abs(iap(ids,1))) = stw(abs(iap(ids,1))-2)
     &                                    + ( strla(iah(idv-1,1))
     &                                      - strlu(abs(iap(ids,1))-1) )
     &                                    / 2.0

                     ixs(abs(iap(ids,1))) = abs(iap(ids,1)) - 2

                     xps(abs(iap(ids,2))) = stw(abs(iap(ids,1))-2)
     &                                    + ( strla(iah(idv-1,1))
     &                                      - strld(abs(iap(ids,1))-1) )
     &                                    / 2.0

                     ixs(abs(iap(ids,2))) = abs(iap(ids,1)) - 2

                  end if

                  if( icwt .eq. 67 ) then

                     yps(abs(iap(ids,2))) = stb(abs(iap(ids,1))-1)
     &                                    - fts(abs(iap(ids,1))-1)
     &                                    * ypls
     &                                    - strha(abs(iap(ids,2)))

                     iys(abs(iap(ids,2))) = abs(iap(ids,1)) - 2

                     ypp(abs(iap(ids,2))) = ypp(iys(abs(iap(ids,2))))
     &                                    + yps(abs(iap(ids,2)))

                  end if

                  if( icwt .eq. 6 .or. icwt .eq. 7 ) then

                     yps(abs(iap(ids,1))) = stb(abs(iap(ids,1))-1)
     &                                    - fts(abs(iap(ids,1))-1)
     &                                    * ypls
     &                                    - strha(abs(iap(ids,1)))

                     iys(abs(iap(ids,1))) = abs(iap(ids,1)) - 2

                     ypp(abs(iap(ids,1))) = ypp(iys(abs(iap(ids,1))))
     &                                    + yps(abs(iap(ids,1)))

                  end if

               end if

                  if( icwt .eq. 6 .or. icwt .eq. 106 ) then

                     stwh(abs(iap(ids,1))-1) = ( strla(iah(idv-1,1))
     &                                       - stw(abs(iap(ids,1))-1) )
     &                                       / 2.0

                  end if


                  if( icwt .eq. 67 .or. icwt .eq. 167 ) then

                     stwh(abs(iap(ids,1))-1) = ( strla(iah(idv-1,1))
     &                                       - stw(abs(iap(ids,1))-1) )
     &                                       / 2.0

                  end if

                  if( ( icwt .eq. 7 .or. icwt .eq. 107 ) .and.
     &                  idv .gt. 1 ) then

                     strha(iah(idv-1,1)) = max( strha(iah(idv-1,1)),
     &                                     strha(iah(idv,1))
     &                                   + ( ypp(iah(idv,1))
     &                                     - ypp(iah(idv-1,1)) ) )

                     strba(iah(idv-1,1)) = min( strba(iah(idv-1,1)),
     &                                     strba(iah(idv,1))
     &                                   + ( ypp(iah(idv,1))
     &                                     - ypp(iah(idv-1,1)) ) )

                  end if


         else if( icwt .eq. 2 .or.
     &            icwt .eq. 5 ) then

                  if( ibt(abs(iap(ids,1)),1) .ge. 1 .or.
     &                icwt .eq. 5 ) then

                     stws(abs(iap(ids,1))-1) = 0.0

                  else

                     stws(abs(iap(ids,1))-1) = max(
     &                                         strlu(abs(iap(ids,1))-1),
     &                                         strld(abs(iap(ids,1))-1))

                  end if


               if( ids - 1 .le. iah(idv,2) ) then

                     strla(iah(idv,1)) = strla(iah(idv,1))
     &                                 + stws(abs(iap(ids,1))-1)

               else

                  if( iap(ids-1,1) .gt. 0 ) then

                     strlu(iap(ids-1,1)-1) = strlu(iap(ids-1,1)-1)
     &                                     + stws(abs(iap(ids,1))-1)

                  else if( iap(ids-1,1) .lt. 0 ) then

                     strld(-iap(ids-1,1)-1) = strld(-iap(ids-1,1)-1)
     &                                      + stws(abs(iap(ids,1))-1)

                  end if

               end if

         else if( icwt .eq. 24 ) then

                  if( iap(ids,1) .gt. 0 ) then

                     strlu(iap(ids,1)-1) = 0.0

                  else if( iap(ids,1) .lt. 0 ) then

                     strld(-iap(ids,1)-1) = 0.0

                  end if

         else if( icwt .eq. 18 ) then

                     rsid(i-2) = stw(i-2)
     &                         - stwrl4(3,3) * fts(i-2) / 1000.0

                     rsid(i-1) = stwrl4(2,4) * fts(i-1) / 1000.0

                     rlng(i-1) = stwrl4(3,4) * fts(i-1) / 1000.0
     &                         - stwrl4(2,4) * fts(i-1) / 1000.0

         else if( icwt .eq. 80 ) then

                     rsca(iah(idv,1)+2) = ( strla(iah(idv,1)) +
     &                                      rsid(iah(idv,1)+1) )
     &                                  /   rlng(iah(idv,1)+2)

                     xps(iah(idv,1)+2) = stw(iah(idv,1))
     &                                 - rsid(iah(idv,1)+2)
     &                                 * ( rsca(iah(idv,1)+2) - 1.0 )

                     ixs(iah(idv,1)+2) = iah(idv,1)

                     xpp(iah(idv,1)+2) = xpp(ixs(iah(idv,1)+2))
     &                                 + xps(iah(idv,1)+2)


                     sth(iah(idv,1)+2) = sth(iah(idv,1)+2)
     &                                 + ypp(iah(idv,1)+2)
     &                                 - ypp(iah(idv,1))

                     stb(iah(idv,1)+2) = stb(iah(idv,1)+2)
     &                                 + ypp(iah(idv,1)+2)
     &                                 - ypp(iah(idv,1))

                     strha(iah(idv,1)) = max(strha(iah(idv,1)),
     &                                       sth(iah(idv,1)+2))

                     strba(iah(idv,1)) = max(strba(iah(idv,1)),
     &                                       stb(iah(idv,1)+2))

                     stwh(iah(idv,1)) = stw(iah(idv,1)+1)
     &                                + strla(iah(idv,1))

         else if( icwt .eq. 8 ) then

                     strla(iah(idv,1)) = 0.0
                     strha(iah(idv,1)) = 0.0
                     strba(iah(idv,1)) = 0.0

         else if( icwt .eq. 90 ) then

                     stwfu(i) = strla(iah(idv,1))
                     sthfu(i) = strha(iah(idv,1))
                     stbfu(i) = strba(iah(idv,1))

         else if( icwt .eq. 91 ) then

                     stwfd(iah(idv,1)) = strla(iah(idv,1))
                     sthfd(iah(idv,1)) = strha(iah(idv,1))
                     stbfd(iah(idv,1)) = strba(iah(idv,1))

                     stwh(iah(idv,5)-1) = max( stwfu(iah(idv,1)),
     &                                         stwfd(iah(idv,1)) )

         else if( icwt .eq. 92 ) then

                     strhb(i) = ( sth(i) + stb(i) ) / 2.0
                     rsid(i)  = stwrl4(2,1) * fts(i) / 1000.0

                     rsca(i) = stwh(iah(idv,5)-1)
     &                       / ( stw(i) - rsid(i) * 2.0 )


                     xps(i) = stw(iah(idv,5)-1)
     &                      - rsid(i) * rsca(i)
                     ixs(i) = iah(idv,5) - 1
                     xpp(i) = xpp(ixs(i)) + xps(i)


                     xps(iah(idv,5)) = stw(iah(idv,5)-1)
     &                               + ( stwh(iah(idv,5)-1)
     &                                 - stwfu(iah(idv,1)) ) / 2.0
                     ixs(iah(idv,5)) = iah(idv,5) - 1
                     xpp(iah(idv,5)) = xpp(ixs(iah(idv,5)))
     &                               + xps(iah(idv,5))


                     xps(iah(idv,1)) = stw(iah(idv,5)-1)
     &                               + ( stwh(iah(idv,5)-1)
     &                                 - stwfd(iah(idv,1)) ) / 2.0
                     ixs(iah(idv,1)) = iah(idv,5) - 1
                     xpp(iah(idv,1)) = xpp(ixs(iah(idv,1)))
     &                               + xps(iah(idv,1))


                  if( iah(idv,4) .eq. 7 ) then

                     yps(iah(idv,5)) = strhb(i) - stbfu(iah(idv,1))
     &                               + fts(iah(idv,5)) * ypls

                     iys(iah(idv,5)) = iah(idv,5) - 1

                     ypp(iah(idv,5)) = ypp(iys(iah(idv,5)))
     &                               + yps(iah(idv,5))


                     yps(iah(idv,1)) = strhb(i) - sthfd(iah(idv,1))
     &                               - fts(iah(idv,1)) * ypls

                     iys(iah(idv,1)) = iah(idv,5) - 1

                     ypp(iah(idv,1)) = ypp(iys(iah(idv,1)))
     &                               + yps(iah(idv,1))


                  else if( iah(idv,4) .eq. 8 ) then

                     yps(iah(idv,5)) = strhb(i) - stbfu(iah(idv,1))
     &                               + fts(iah(idv,5)) * ftss * ypls

                     iys(iah(idv,5)) = iah(idv,5) - 1

                     ypp(iah(idv,5)) = ypp(iys(iah(idv,5)))
     &                               + yps(iah(idv,5))


                     yps(iah(idv,1)) = strhb(i) - sthfd(iah(idv,1))
     &                               - fts(iah(idv,1)) * ftss * ypls

                     iys(iah(idv,1)) = iah(idv,5) - 1

                     ypp(iah(idv,1)) = ypp(iys(iah(idv,1)))
     &                               + yps(iah(idv,1))


                  end if

                     strla(iah(idv,5)-1) = stwh(iah(idv,5)-1)

                     strha(iah(idv,5)-1) = sthfu(iah(idv,1))
     &                                   + ( ypp(iah(idv,5))
     &                                     - ypp(iah(idv,5)-1) )

                     strba(iah(idv,5)-1) = stbfd(iah(idv,1))
     &                                   + ( ypp(iah(idv,1))
     &                                     - ypp(iah(idv,5)-1) )


         else if( icwt .eq. 9  .or.
     &            icwt .eq. 11 .or.
     &            icwt .eq. 12 .or.
     &            icwt .eq. 19 .or.
     &            icwt .eq. 29 .or.
     &            icwt .eq. 39 .or.
     &            icwt .eq. 49 .or.
     &            icwt .eq. 69 .or.
     &            icwt .eq. 70 ) then


               if( icwt .eq.  9 .or.
     &             icwt .eq. 39 .or.
     &             icwt .eq. 49 ) then

                     stwh(abs(ixp(i))-1) = max( strla(iah(idv,1)),
     &                                          stw(i-1) )

               else if( icwt .eq. 19 .or.
     &                  icwt .eq. 29 ) then

                     stwh(abs(ixp(i))-1) = strla(iah(idv,1))


               else if( icwt .eq. 11 ) then

                     strha(iah(idv,1)) = sth(i-2)
                     strba(iah(idv,1)) = stb(i-1)

                     stwh(abs(ixp(i))-1) = stw(i-3)


               else if( icwt .eq. 12 ) then

                     strha(iah(idv,1)) = sth(i-2)
                     strba(iah(idv,1)) = stb(i-1)

                     stwh(abs(ixp(i))-1) = stw(i-3) + stw(i-3) * 0.5


               end if

            if( icwt .eq. 69 .and. idv .gt. 1 ) then

                     strha(iah(idv-1,1)) = max( strha(iah(idv-1,1)),
     &                                          strha(iah(idv,1))
     &                                      + ( ypp(iah(idv,1))
     &                                        - ypp(iah(idv-1,1)) ) )


                     strba(iah(idv-1,1)) = min( strba(iah(idv-1,1)),
     &                                          strba(iah(idv,1))
     &                                      + ( ypp(iah(idv,1))
     &                                        - ypp(iah(idv-1,1)) ) )


            end if


            if( icwt .eq. 70 ) then

               if( idv .eq. 1 ) then

                  if( iap(ids,1) .gt. 0 ) then

                     strlu(iap(ids,1)-1) = strlu(iap(ids,1)-1)
     &                                   + strla(iah(idv,1))

                  else if( iap(ids,1) .lt. 0 ) then

                     strld(-iap(ids,1)-1) = strld(-iap(ids,1)-1)
     &                                    + strla(iah(idv,1))

                  end if

               else if( ids .le. iah(idv-1,2) ) then

                     strla(iah(idv-1,1)) = strla(iah(idv-1,1))
     &                                   + strla(iah(idv,1))

                     strha(iah(idv-1,1)) = max( strha(iah(idv-1,1)),
     &                                          strha(iah(idv,1))
     &                                      + ( ypp(iah(idv,1))
     &                                        - ypp(iah(idv-1,1)) ) )

                     strba(iah(idv-1,1)) = min( strba(iah(idv-1,1)),
     &                                          strba(iah(idv,1))
     &                                      + ( ypp(iah(idv,1))
     &                                        - ypp(iah(idv-1,1)) ) )


               else if( ids .gt. iah(idv-1,2) ) then

                  if( iap(ids,1) .gt. 0 ) then

                     strlu(iap(ids,1)-1) = strlu(iap(ids,1)-1)
     &                                   + strla(iah(idv,1))

                  else if( iap(ids,1) .lt. 0 ) then

                     strld(-iap(ids,1)-1) = strld(-iap(ids,1)-1)
     &                                    + strla(iah(idv,1))

                  end if

                     strha(iah(idv-1,1)) = max( strha(iah(idv-1,1)),
     &                                          strha(iah(idv,1))
     &                                      + ( ypp(iah(idv,1))
     &                                        - ypp(iah(idv-1,1)) ) )

                     strba(iah(idv-1,1)) = min( strba(iah(idv-1,1)),
     &                                          strba(iah(idv,1))
     &                                      + ( ypp(iah(idv,1))
     &                                        - ypp(iah(idv-1,1)) ) )

               end if


            end if


            if( icwt .ne. 69 .and. icwt .ne. 70 ) then

               if( idv .eq. 1 ) then

                  if( iap(ids,1) .gt. 0 ) then

                     strlu(iap(ids,1)-1) = strlu(iap(ids,1)-1)
     &                                   + stwh(abs(ixp(i))-1)

                  else if( iap(ids,1) .lt. 0 ) then

                     strld(-iap(ids,1)-1) = strld(-iap(ids,1)-1)
     &                                    + stwh(abs(ixp(i))-1)

                  end if

               else if( ids .le. iah(idv-1,2) ) then

                     strla(iah(idv-1,1)) = strla(iah(idv-1,1))
     &                                   + stwh(abs(ixp(i))-1)

               else if( ids .gt. iah(idv-1,2) ) then

                  if( iap(ids,1) .gt. 0 ) then

                     strlu(iap(ids,1)-1) = strlu(iap(ids,1)-1)
     &                                   + stwh(abs(ixp(i))-1)

                  else if( iap(ids,1) .lt. 0 ) then

                     strld(-iap(ids,1)-1) = strld(-iap(ids,1)-1)
     &                                    + stwh(abs(ixp(i))-1)

                  end if

               end if

            end if


            if( icwt .eq. 9 ) then

                  if( isi(i) .eq. 4 ) then

                     xps(i-1) = stw(abs(ixp(i))-1)
                     ixs(i-1) = abs(ixp(i)) - 1
                     xpp(i-1) = xpp(ixs(i-1)) + xps(i-1)

                  else if( isi(i) .ne. 4 ) then

                     xps(abs(ixp(i))) = stw(abs(ixp(i))-1)
     &                                + ( stwh(abs(ixp(i))-1)
     &                                  - strla(iah(idv,1)) ) / 2.0

                     ixs(abs(ixp(i))) = abs(ixp(i)) - 1

                     xpp(abs(ixp(i))) = xpp(ixs(abs(ixp(i))))
     &                                + xps(abs(ixp(i)))


                     xps(i-1) = stw(abs(ixp(i))-1)
     &                        + ( stwh(abs(ixp(i)) - 1)
     &                          - stw(i-1) ) / 2.0

                     ixs(i-1) = abs(ixp(i)) - 1

                     xpp(i-1) = xpp(ixs(i-1)) + xps(i-1)

                  end if


            else if( icwt .eq. 39 .or.
     &               icwt .eq. 49 ) then

                     xps(abs(ixp(i))) = stw(abs(ixp(i))-1)
     &                                + ( stwh(abs(ixp(i))-1)
     &                                  - strla(iah(idv,1)) ) / 2.0

                     ixs(abs(ixp(i))) = abs(ixp(i)) - 1

                     xpp(abs(ixp(i))) = xpp(ixs(abs(ixp(i))))
     &                                + xps(abs(ixp(i)))


                  if( icwt .eq. 39 ) then

                     xps(i-1) = stw(abs(ixp(i))-1)
     &                        + stwh(abs(ixp(i)) - 1)
     &                        - stw(i-1)

                     ixs(i-1) = abs(ixp(i)) - 1

                     xpp(i-1) = xpp(ixs(i-1)) + xps(i-1)


                  else

                     xps(i-1) = stw(abs(ixp(i))-1)
                     ixs(i-1) = abs(ixp(i)) - 1
                     xpp(i-1) = xpp(ixs(i-1)) + xps(i-1)


                  end if

            else if( icwt .eq. 19 .or. icwt .eq. 29 ) then

                     rsid(i-1) = stwrl4(2,1) * fts(i-1) / 1000.0

                     rsca(i-1) = stwh(abs(ixp(i))-1)
     &                         / ( stw(i-1) - rsid(i-1) * 2.0 )

                     xps(i-1) = stw(abs(ixp(i))-1)
     &                        - rsid(i-1) * rsca(i-1)

                     ixs(i-1) = abs(ixp(i)) - 1

                     xpp(i-1) = xpp(ixs(i-1)) + xps(i-1)


            else if( ( icwt .eq. 11 .or. icwt .eq. 12 ) .and.
     &                 idv .gt. 1 ) then

                     strha(iah(idv-1,1)) = max( strha(iah(idv-1,1)),
     &                                          strha(iah(idv,1))
     &                                      + ( ypp(iah(idv,1))
     &                                        - ypp(iah(idv-1,1)) ) )

                     strba(iah(idv-1,1)) = min( strba(iah(idv-1,1)),
     &                                          strba(iah(idv,1))
     &                                      + ( ypp(iah(idv,1))
     &                                        - ypp(iah(idv-1,1)) ) )

            end if


*-----------------------------------------------------------------------

         else if( icwt .eq. 40 .or.
     &            icwt .eq. 50 ) then

*-----------------------------------------------------------------------

                  if( i .ne. icl(i) ) then

                     clc(i,1) = clc(icl(i),1)
                     clc(i,2) = clc(icl(i),2)
                     clc(i,3) = clc(icl(i),3)

                  else if( i .eq. icl(i) ) then

                     clc(i,1) = clp(i,1,1)
                     clc(i,2) = clp(i,1,2)
                     clc(i,3) = clp(i,1,3)

                  end if

*-----------------------------------------------------------------------

                     jfn(i) = ifn(i)

                     fts(i) = fts(ifs(i))

*-----------------------------------------------------------------------

                     icn = ilt(i)

                     do 801 k = 1, icn

                        chaf(k) = chag(i,k)

  801                continue

                     call strwhb(icn,chaf,fts(i),jfn(i),
     &                           strwf,strhf,strbf,ifon)

                     stw(i) = strwf

*-----------------------------------------------------------------------

                     yps(i) = 0.0
                     iys(i) = abs(iyp(i))
                     ypp(i) = ypp(iys(i)) + yps(i)

                     rsid(i) = stwrl4(2,2) * fts(i) / 1000.0

                     rsca(i) = ( stwh(iah(idv,1)) + rsid(i) * 2.0 )
     &                       / ( stw(i) - rsid(i) * 2.0 )


                  if( icwt .eq. 40 ) then

                     xps(i) = stw(iah(idv,1)) - rsid(i) * rsca(i)
                     ixs(i) = iah(idv,1)
                     xpp(i) = xpp(ixs(i)) + xps(i)

                  else if( icwt .eq. 50 ) then

                     xps(i) = stw(iah(idv,1)) - rsid(i) * rsca(i)
     &                      - rsid(i) * 2.0
                     ixs(i) = iah(idv,1)
                     xpp(i) = xpp(ixs(i)) + xps(i)


                  end if

*-----------------------------------------------------------------------


                     call wstrg(jhf,i,chag,ilt,ifon)


*-----------------------------------------------------------------------

         else if( icwt .eq. 0 ) then

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*              COLOR OF CHARACTERS
*-----------------------------------------------------------------------

                  if( i .ne. icl(i) ) then

                     clc(i,1) = clc(icl(i),1)
                     clc(i,2) = clc(icl(i),2)
                     clc(i,3) = clc(icl(i),3)

                  else if( i .eq. icl(i) ) then

                     clc(i,1) = clp(i,1,1)
                     clc(i,2) = clp(i,1,2)
                     clc(i,3) = clp(i,1,3)

                  end if

*-----------------------------------------------------------------------
*              FONT OF CHARACTER
*-----------------------------------------------------------------------

                     jfn(i) = ifn(i)

*-----------------------------------------------------------------------
*              FONT SIZE
*-----------------------------------------------------------------------

                  if( abs(isi(i)) .ne. 1 .and. isi(i) .ne. 3 .and.
     &                abs(isi(i)) .ne. 21 .and.
     &                abs(isi(i)) .ne. 31 .and.
     &                abs(isi(i)) .ne. 32 .and.
     &                abs(isi(i)) .lt. 100 ) then

                        fts(i) = fts(ifs(i))

                  else if( abs(isi(i)) .eq. 1 .or.
     &                     abs(isi(i)) .eq. 21 ) then

                        fts(i) = fts(ifs(i)) * ftss

                  else if( isi(i) .eq. 3 ) then

                        fts(i) = fts(ifs(i)) * ftsv

                  else if( isi(i) .eq. 101 ) then

                        fts(i) = fts(ifs(i)) * 0.5

                  else if( isi(i) .eq. 31 .or.
     &                     isi(i) .eq. 32 ) then

                        fts(i) = fts(ifs(i)) * 0.6

                  else if( isi(i) .eq. 102 ) then

                        fts(i) = fts(ifs(i)) * ftss

                  else if( isi(i) .eq. 103 ) then

                        fts(i) = fts(ifs(i)) * 0.8

                  else if( isi(i) .eq. 104 ) then

                        fts(i) = fts(ifs(i)) * 0.9

                  else if( isi(i) .eq. 105 ) then

                        fts(i) = fts(ifs(i)) * 1.1

                  else if( isi(i) .eq. 106 ) then

                        fts(i) = fts(ifs(i)) * 1.2

                  else if( isi(i) .eq. 107 ) then

                        fts(i) = fts(ifs(i)) * 1.4

                  else if( isi(i) .eq. 108 ) then

                        fts(i) = fts(ifs(i)) * 1.7

                  else if( isi(i) .eq. 109 ) then

                        fts(i) = fts(ifs(i)) * 2.0

                  else if( abs(isi(i)) .eq. 110 ) then

                        fts(i) = fts(ifs(i)) * ftks

                  end if

*-----------------------------------------------------------------------
*           STRING WIDTH
*-----------------------------------------------------------------------

               if( idh .ge. 0 ) then

*-----------------------------------------------------------------------

                     icn = ilt(i)

                  if( icn .eq. 0 ) then

                     strwf = 0.0
                     strhf = 0.0
                     strbf = 0.0

                  else

                     do 800 k = 1, icn

                        chaf(k) = chag(i,k)

  800                continue

                     call strwhb(icn,chaf,fts(i),jfn(i),
     &                           strwf,strhf,strbf,ifon)

                  end if

                     stw(i) = strwf
                     sth(i) = strhf
                     stb(i) = strbf

*-----------------------------------------------------------------------

               else if( idh .eq. -1 ) then

                     stw(i) = fts(i)
                     sth(i) = 0.0
                     stb(i) = 0.0

               else if( idh .eq. -2 ) then

                     stw(i) = fts(i) * 2.0
                     sth(i) = 0.0
                     stb(i) = 0.0

               else if( idh .eq. -3 ) then

                     stw(i) = fts(i) * 3.0 / 18.0
                     sth(i) = 0.0
                     stb(i) = 0.0

               else if( idh .eq. -4 ) then

                     stw(i) = fts(i) * 4.0 / 18.0
                     sth(i) = 0.0
                     stb(i) = 0.0

               else if( idh .eq. -5 ) then

                     stw(i) = fts(i) * 5.0 / 18.0
                     sth(i) = 0.0
                     stb(i) = 0.0

               else if( idh .eq. -6 ) then

                     stw(i) = - fts(i) * 3.0 / 18.0
                     sth(i) = 0.0
                     stb(i) = 0.0

               else if( idh .eq. -7 ) then

                     stw(i) = fts(i) * phs
                     sth(i) = 0.0
                     stb(i) = 0.0

               else if( idh .eq. -10 ) then

                     stw(i) = fts(i) * bxws
                     sth(i) = 0.0
                     stb(i) = 0.0

               else if( idh .eq. -11 ) then

                     stw(i) = fts(i) * bxwl
                     sth(i) = 0.0
                     stb(i) = 0.0

               else if( idh .eq. -12 ) then

                     stw(i) = fts(i) * bxws
                     sth(i) = strha(iah(idv+1,1)) + fts(i) * bxws
     &                      + ypp(iah(idv+1,1)) - ypp(iah(idv,1))
                     stb(i) = strba(iah(idv+1,1)) - fts(i) * bxws
     &                      + ypp(iah(idv+1,1)) - ypp(iah(idv,1))

               else if( idh .eq. -13 ) then

                     stw(i) = fts(i) * bxwl
                     sth(i) = strha(iah(idv+1,1)) + fts(i) * bxwl
     &                      + ypp(iah(idv+1,1)) - ypp(iah(idv,1))
                     stb(i) = strba(iah(idv+1,1)) - fts(i) * bxwl
     &                      + ypp(iah(idv+1,1)) - ypp(iah(idv,1))

               else if( idh .eq. -14 ) then

                     stw(i) = fts(i) * ( bxws + bxss )
                     sth(i) = strha(iah(idv+1,1)) + fts(i) * bxws
     &                      + ypp(iah(idv+1,1)) - ypp(iah(idv,1))
                     stb(i) = strba(iah(idv+1,1))
     &                      - fts(i) * ( bxws + bxss )
     &                      + ypp(iah(idv+1,1)) - ypp(iah(idv,1))

               else if( idh .eq. -15 ) then

                     stw(i) = fts(i) * ( bxwl + bxsl )
                     sth(i) = strha(iah(idv+1,1)) + fts(i) * bxwl
     &                      + ypp(iah(idv+1,1)) - ypp(iah(idv,1))
                     stb(i) = strba(iah(idv+1,1))
     &                      - fts(i) * ( bxwl + bxsl )
     &                      + ypp(iah(idv+1,1)) - ypp(iah(idv,1))

               end if

*-----------------------------------------------------------------------
*           WRITE STRING
*-----------------------------------------------------------------------

               if( ilt(i) .gt. 0 ) then

                     call wstrg(jhf,i,chag,ilt,ifon)

               end if

*-----------------------------------------------------------------------
*           X POSITION
*-----------------------------------------------------------------------

            if( ( idh .le. 0 .or. idh .ge. 2 ) .and.
     &        ( ( ids .eq. 0 ) .or.
     &          ( iap(ids,3) .ne. 4 .or.
     &          ( iap(ids,1) .ne. i .and. iap(ids,2) .ne. i ) ) ) ) then


               if( idh .eq. 2 ) then

                     xps(i) = stw(abs(ixp(i))-1) + stwh(abs(ixp(i))-1)
                     ixs(i) = abs(ixp(i)) - 1
                     xpp(i) = xpp(ixs(i)) + xps(i)

               else if( ixp(i) .gt. 0 ) then

                  if( ibt(i,1) .eq. 1 .and. isi(i) .eq. 1 ) then

                     xps(i) = stw(ixp(i)-1) + stw(ibt(i,2)-2) / 5.0
                     ixs(i) = ixp(i)-1
                     xpp(i) = xpp(ixs(i)) + xps(i)

                  else if( ibt(i,1) .eq. 1 .and. isi(i) .eq. -1 ) then

                     xps(i) = stw(ixp(i)-1) - stw(ibt(i,2)-1) / 4.0
                     ixs(i) = ixp(i)-1
                     xpp(i) = xpp(ixs(i)) + xps(i)

                  else

                     xps(i) = stw(ixp(i)-1)
                     ixs(i) = ixp(i)-1
                     xpp(i) = xpp(ixs(i)) + xps(i)

                  end if

               else if( ixp(i) .lt. 0 ) then

                     xps(i) = stw(abs(ixp(i))-1) + stws(abs(ixp(i))-1)
                     ixs(i) = abs(ixp(i))-1
                     xpp(i) = xpp(ixs(i)) + xps(i)

               end if

            end if

*-----------------------------------------------------------------------
*           Y POSITION
*-----------------------------------------------------------------------

               if( abs(isi(i)) .eq. 0 .or.
     &             abs(isi(i)) .eq. 21 .or.
     &           ( isi(i) .ge. 100 .and. isi(i) .le. 109 ) ) then

                     yps(i) = 0.0
                     iys(i) = abs(iyp(i))
                     ypp(i) = ypp(iys(i)) + yps(i)

               else if( isi(i) .eq. 1 .or. isi(i) .eq. 110 ) then

                  if( ibt(i,1) .eq. 1 .or.
     &                ibt(i,1) .eq. 2 ) then

                     yps(i) = fts(ibt(i,2)-2) * 0.65
                     iys(i) = ibt(i,2) - 2
                     ypp(i) = ypp(iys(i)) + yps(i)

                  else if( ibt(i,1) .eq. 4 ) then

                     yps(i) = fts(ifs(i)) * ypku
                     iys(i) = abs(iyp(i))
                     ypp(i) = ypp(iys(i)) + yps(i)

                  else

                     yps(i) = fts(ifs(i)) * yplu
                     iys(i) = abs(iyp(i))
                     ypp(i) = ypp(iys(i)) + yps(i)

                  end if

               else if( isi(i) .eq. -1 .or. isi(i) .eq. -110 ) then

                  if( ibt(i,1) .eq. 1 .or.
     &                ibt(i,1) .eq. 2 ) then

                     yps(i) = stb(ibt(i,2)-1)
     &                      - fts(ibt(i,2)-1) * 0.2
                     iys(i) = 0
                     ypp(i) = ypp(iys(i)) + yps(i)

                  else if( ibt(i,1) .eq. 4 ) then

                     yps(i) = fts(ifs(i)) * ypkd
                     iys(i) = abs(iyp(i))
                     ypp(i) = ypp(iys(i)) + yps(i)

                  else

                     yps(i) = fts(ifs(i)) * ypld
                     iys(i) = abs(iyp(i))
                     ypp(i) = ypp(iys(i)) + yps(i)

                  end if

               else if( isi(i) .eq. 2 .or.
     &                  isi(i) .eq. 3 .or.
     &                  isi(i) .eq. 31.or.
     &                  isi(i) .eq. 32 ) then

                     strhh(i) = ( sth(i) + stb(i) ) / 2.0

                     yplhh = yplh
                     if( isi(i) .eq. 31 ) yplhh = 0.1

                     yps(i) = strha(iah(idv,1))
     &                      + fts(iah(idv,1)) * yplhh
     &                      - strhh(i)
                     iys(i) = iah(idv,1)
                     ypp(i) = ypp(iys(i)) + yps(i)

                  if( idh .eq. 1 .and. idv .gt. 1 ) then

                     sth(i) = ( sth(i) + stb(i) ) / 2.0
     &                      + ypp(i) - ypp(iah(idv-1,1))

                     strha(iah(idv-1,1)) = max( strha(iah(idv-1,1)),
     &                                          sth(i) )

                     strba(iah(idv-1,1)) = min( strba(iah(idv-1,1)),
     &                                          strba(iah(idv,1))
     &                                      + ( ypp(iah(idv,1))
     &                                        - ypp(iah(idv-1,1)) ) )

                  end if

               else if( isi(i) .eq. 5 ) then

                     strbb(i) = ( sth(i) + stb(i) ) / 2.0

                     yps(i) = strba(iah(idv,1))
     &                      - fts(iah(idv,1)) * yplh
     &                      - strbb(i)
                     iys(i) = iah(idv,1)
                     ypp(i) = ypp(iys(i)) + yps(i)

                  if( idh .eq. 1 .and. idv .gt. 1 ) then

                     stb(i) = ( sth(i) + stb(i) ) / 2.0
     &                      + ypp(i) - ypp(iah(idv-1,1))

                     strha(iah(idv-1,1)) = max( strha(iah(idv-1,1)),
     &                                          strha(iah(idv,1))
     &                                      + ( ypp(iah(idv,1))
     &                                        - ypp(iah(idv-1,1)) ) )

                     strba(iah(idv-1,1)) = min( strba(iah(idv-1,1)),
     &                                          stb(i) )

                  end if

               else if( isi(i) .eq. 4 ) then

                     yps(i) = 0.0
                     iys(i) = abs(iyp(i))
                     ypp(i) = ypp(iys(i)) + yps(i)

                  if( idh .eq. 1 .and. idv .gt. 1 ) then

                     strha(iah(idv-1,1)) = max( strha(iah(idv-1,1)),
     &                                          strha(iah(idv,1))
     &                                      + ( ypp(iah(idv,1))
     &                                        - ypp(iah(idv-1,1)) ) )

                  end if

               else if( isi(i) .eq.  6 .or.
     &                  isi(i) .eq.  9 .or.
     &                  isi(i) .eq. 12 ) then

                     strhh(i) = sth(i)
                     strbb(i) = stb(i)

                     yps(i) = ( fts(iah(idv,1)) * ftht
     &                        - ( strhh(i) + strbb(i) ) ) / 2.0
                     iys(i) = iah(idv,1)
                     ypp(i) = ypp(iys(i)) + yps(i)

                     sth(i) = sth(i) + ypp(i) - ypp(iah(idv,1))
                     stb(i) = stb(i) + ypp(i) - ypp(iah(idv,1))

               else if( isi(i) .eq.  7 .or.
     &                  isi(i) .eq. 10 .or.
     &                  isi(i) .eq. 13 ) then

                     strbb(i) = stb(i)

                  if( isi(i) .eq. 7 ) then

                     yps(i) = ( sth(i-1) - strbb(i) ) * 0.99
                     iys(i) = iah(idv,1)
                     ypp(i) = ypp(iys(i)) + yps(i)

                  else if( isi(i) .eq. 10 ) then

                     yps(i) = ( sth(i-1) - strbb(i) ) * 0.55
                     iys(i) = iah(idv,1)
                     ypp(i) = ypp(iys(i)) + yps(i)

                  else if( isi(i) .eq. 13 ) then

                     yps(i) = ( sth(i-1) - strbb(i) ) * 0.72
                     iys(i) = iah(idv,1)
                     ypp(i) = ypp(iys(i)) + yps(i)

                  end if

               else if( isi(i) .eq.  8 .or.
     &                  isi(i) .eq. 11 .or.
     &                  isi(i) .eq. 14 ) then

                     strhh(i) = sth(i)

                  if( isi(i) .eq. 8 ) then

                     yps(i) = ( stb(i-2) - strhh(i) ) * 0.99
                     iys(i) = iah(idv,1)
                     ypp(i) = ypp(iys(i)) + yps(i)

                  else if( isi(i) .eq. 11 ) then

                     yps(i) = ( stb(i-2) - strhh(i) ) * 0.55
                     iys(i) = iah(idv,1)
                     ypp(i) = ypp(iys(i)) + yps(i)

                  else if( isi(i) .eq. 14 ) then

                     yps(i) = ( stb(i-2) - strhh(i) ) * 0.72
                     iys(i) = iah(idv,1)
                     ypp(i) = ypp(iys(i)) + yps(i)

                  end if

               end if

*-----------------------------------------------------------------------


            if( ( ilt(i) .gt. 0 .and.
     &          ( idh .eq. 0 .or. idh .eq. 2 ) ) .or.
     &            idh .lt. 0 ) then

               if( ids .le. iah(idv,2) ) then


                  strla(iah(idv,1)) = strla(iah(idv,1)) + stw(i)


               else if( ids .gt. iah(idv,2) ) then

                  if( iap(ids,1) .gt. 0 ) then

                     strlu(iap(ids,1)-1) = strlu(iap(ids,1)-1)
     &                                   + stw(i)

                  else if( iap(ids,1) .lt. 0 ) then

                     strld(-iap(ids,1)-1) = strld(-iap(ids,1)-1)
     &                                    + stw(i)

                  end if

               end if

                  if( idh .gt. -12 ) then

                     sth(i) = sth(i) + ypp(i) - ypp(iah(idv,1))
                     stb(i) = stb(i) + ypp(i) - ypp(iah(idv,1))

                  end if

                     strha(iah(idv,1)) = max(strha(iah(idv,1)),sth(i))
                     strba(iah(idv,1)) = min(strba(iah(idv,1)),stb(i))

            end if

*-----------------------------------------------------------------------
*        WRITE STRING AND POSITION
*-----------------------------------------------------------------------

         else if( icwt .eq. 10 ) then

*-----------------------------------------------------------------------

            if( iwt(i) .ne. 1 ) then

                  if( ilt(i) .gt. 0 ) then

                     write(jhf,'(3f7.3,'' sc str'',I3.3,'' xps '',
     &                           g14.5,'' add'',/''yps '',g14.5,
     &                           '' add ft'',i3.3,'' '',
     &                           g14.5,'' stsw'')')
     &               clc(i,1), clc(i,2), clc(i,3),
     &               i, xpp(i), ypp(i), jfn(i), fts(i)

                  end if

*-----------------------------------------------------------------------
*              STRING BOX OUTPUT
*-----------------------------------------------------------------------

               if( iwt(i) .ge. 10 .and. iwt(i) .le. 53 ) then

*-----------------------------------------------------------------------

                        b1x = xpp(i)

                  if( iwt(i) .eq. 10 .or. iwt(i) .eq. 11 .or.
     &                iwt(i) .eq. 20 .or. iwt(i) .eq. 21 .or.
     &                iwt(i) .eq. 30 .or. iwt(i) .eq. 31 .or.
     &                iwt(i) .eq. 40 .or. iwt(i) .eq. 41 .or.
     &                iwt(i) .eq. 50 .or. iwt(i) .eq. 51 ) then

                        b2x = xpp(i) + strla(i) + fts(i) * bxws * 2.0
                        b1y = ypp(i) + strba(i) - fts(i) * bxws
                        b2y = ypp(i) + strha(i) + fts(i) * bxws

                  else

                        b2x = xpp(i) + strla(i) + fts(i) * bxwl * 2.0
                        b1y = ypp(i) + strba(i) - fts(i) * bxwl
                        b2y = ypp(i) + strha(i) + fts(i) * bxwl

                  end if

                  if( iwt(i)/2*2 .eq. iwt(i) ) then

                        bdd = fts(i) * bxls

                  else

                        bdd = fts(i) * bxll

                  end if

                        bcb(1) = clp(i+1,2,1)
                        bcb(2) = clp(i+1,2,2)
                        bcb(3) = clp(i+1,2,3)

                  if( iwt(i) .ge. 10 .and. iwt(i) .le. 33 ) then

                     if( clp(i+1,3,1) .gt. -2.5 ) then

                        bcl(1) = clp(i+1,3,1)
                        bcl(2) = clp(i+1,3,2)
                        bcl(3) = clp(i+1,3,3)

                     else

                        bcl(1) = clc(i,1)
                        bcl(2) = clc(i,2)
                        bcl(3) = clc(i,3)

                     end if

                  else if( iwt(i) .ge. 40 .and. iwt(i) .le. 53 ) then

                        bcl(1) = clc(i,1)
                        bcl(2) = clc(i,2)
                        bcl(3) = clc(i,3)

                  end if

                  if( iwt(i) .eq. 20 .or. iwt(i) .eq. 21 .or.
     &                iwt(i) .eq. 40 .or. iwt(i) .eq. 41 .or.
     &                iwt(i) .eq. 50 .or. iwt(i) .eq. 51 ) then

                        bdc = fts(i) * bxss

                  else if( iwt(i) .eq. 22 .or. iwt(i) .eq. 23 .or.
     &                     iwt(i) .eq. 42 .or. iwt(i) .eq. 43 .or.
     &                     iwt(i) .eq. 52 .or. iwt(i) .eq. 53 ) then

                        bdc = fts(i) * bxsl

                  else if( iwt(i) .eq. 30 .or. iwt(i) .eq. 31 ) then

                        bdc = fts(i) * bxds

                  else if( iwt(i) .eq. 32 .or. iwt(i) .eq. 33 ) then

                        bdc = fts(i) * bxdl

                  end if

                  if( iwt(i) .ge. 40 .and. iwt(i) .le. 53 ) then

                     if( clp(i+1,3,1) .gt. -2.5 ) then

                        bcs(1) = clp(i+1,3,1)
                        bcs(2) = clp(i+1,3,2)
                        bcs(3) = clp(i+1,3,3)

                     else

                        bcs(1) = clc(i,1)
                        bcs(2) = clc(i,2)
                        bcs(3) = clc(i,3)

                     end if

                  end if

*-----------------------------------------------------------------------
*                 WRITE BOX
*-----------------------------------------------------------------------

                  iccb = 0

                  call wrbox(jhf,idbg,1,iwt(i),iccb,
     &                       b1x, b2x, b1y, b2y,
     &                       bdd, bcb, bcl, bdc, bcs)

*-----------------------------------------------------------------------

               end if

*-----------------------------------------------------------------------
*        WITH SCALE STRING
*-----------------------------------------------------------------------

            else if( iwt(i) .eq. 1 ) then

                  if( ilt(i) .gt. 0 ) then

                     write(jhf,'(3f7.3,'' sc str'',i3.3,'' xps '',
     &                           g14.5,'' add'',/''yps '',g14.5,
     &                           '' add '',
     &                           g14.5,'' 1 ft'',i3.3,'' '',
     &                           g14.5,'' stss'')')
     &               clc(i,1), clc(i,2), clc(i,3),
     &               i, xpp(i), ypp(i), rsca(i), jfn(i), fts(i)

                  end if

            end if

         end if

*-----------------------------------------------------------------------
         idh = 0
*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine wstrg(jhf,i,chag,ilt,ifon)
*                                                                      *
*        PURPOSE    :  DELETE THE LAST CONTROL CHARACTER FOR KANJI     *
*                                                                      *
*                      INSERT '\' BEFORE PS SPECIAL CHARACTER          *
*                      '(', ')', AND '\' IN THE JAPANESE MODE          *
*                                                                      *
*                      WRITE STRING ON FILE FOR TEXT                   *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character chag(inig,0:ichrl)*1
      character chaf(ichrl)*1
      dimension ilt(0:inig)

      character cendef(3)*1
      data      cendef/')',' ','N'/

      character yen*1

      logical numec

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------

            yen = char(92)

*-----------------------------------------------------------------------
*     DELETE THE LAST CONTROL CHARACTER FOR JAPANESE KANJI
*     AND INSERT '\' BEFORE PS SPECIAL CHARACTER '(', ')', AND '\'
*     IN THE JAPANESE MODE
*-----------------------------------------------------------------------

         if( ifon .gt. 0 ) then

*-----------------------------------------------------------------------
*     DELETE THE LAST CONTROL CHARACTER FOR JAPANESE KANJI
*-----------------------------------------------------------------------
cKN 2024/01/24

               if( ilt(i) .gt. 8 ) then
                  if( chag(i,ilt(i)-7) .eq. yen .and.
     &                chag(i,ilt(i)-6) .eq. '3' .and.
     &                chag(i,ilt(i)-5) .eq. '7' .and.
     &                chag(i,ilt(i)-4) .eq. '7' .and.
     &                chag(i,ilt(i)-3) .eq. yen .and.
     &                chag(i,ilt(i)-2) .eq. '0' .and.
     &                chag(i,ilt(i)-1) .eq. '0' ) ilt(i) = ilt(i) - 8
               end if

*-----------------------------------------------------------------------
*     INSERT '\' BEFORE PS SPECIAL CHARACTER '(', ')', AND '\'
*     IN THE JAPANESE MODE
*-----------------------------------------------------------------------

            if( ilt(i) .gt. 7 ) then

                  j = 0
                  k = 0
                  l = 0

  101          j = j + 1

                  if( chag(i,j)   .eq. yen .and.
     &                chag(i,j+1) .eq. '3' .and.
     &                chag(i,j+2) .eq. '7' .and.
     &                chag(i,j+3) .eq. '7' .and.
     &                chag(i,j+4) .eq. yen .and.
     &                chag(i,j+5) .eq. '0' .and.
     &                chag(i,j+6) .eq. '0' .and.
     &                chag(i,j+7) .eq. '1' ) then

                     k = 1

                        l = l + 1

                        chaf(l)   = chag(i,j)
                        chaf(l+1) = chag(i,j+1)
                        chaf(l+2) = chag(i,j+2)
                        chaf(l+3) = chag(i,j+3)
                        chaf(l+4) = chag(i,j+4)
                        chaf(l+5) = chag(i,j+5)
                        chaf(l+6) = chag(i,j+6)
                        chaf(l+7) = chag(i,j+7)

                        l = l + 7
                        j = j + 7

                  else
     &            if( chag(i,j)   .eq. yen .and.
     &                chag(i,j+1) .eq. '3' .and.
     &                chag(i,j+2) .eq. '7' .and.
     &                chag(i,j+3) .eq. '7' .and.
     &                chag(i,j+4) .eq. yen .and.
     &                chag(i,j+5) .eq. '0' .and.
     &                chag(i,j+6) .eq. '0' .and.
     &                chag(i,j+7) .eq. '0' ) then

                     k = 0

                        l = l + 1

                        chaf(l)   = chag(i,j)
                        chaf(l+1) = chag(i,j+1)
                        chaf(l+2) = chag(i,j+2)
                        chaf(l+3) = chag(i,j+3)
                        chaf(l+4) = chag(i,j+4)
                        chaf(l+5) = chag(i,j+5)
                        chaf(l+6) = chag(i,j+6)
                        chaf(l+7) = chag(i,j+7)

                        l = l + 7
                        j = j + 7

                  else
     &            if( k .eq. 1 .and.
     &              ( chag(i,j) .eq. '(' .or.
     &                chag(i,j) .eq. ')' .or.
     &                chag(i,j) .eq. yen ) ) then

                     l = l + 1
                     chaf(l) = yen
                     l = l + 1
                     chaf(l) = chag(i,j)

                  else

                     l = l + 1
                     chaf(l) = chag(i,j)

                  end if


                  if( j .lt. ilt(i) ) goto 101


                  ilt(i) = l

               do 102 l = 1, ilt(i)

                  chag(i,l) = chaf(l)

  102          continue


            end if

*-----------------------------------------------------------------------

         end if



*-----------------------------------------------------------------------
*     CUT LONG STRING
*-----------------------------------------------------------------------

                  icut  =  0
                  icont =  0
                  icut1 = 60

               if( ilt(i) .gt. icut1 + 3 ) then

  300             icont = icont + 1
                  icut2 = icut + icut1
                  idec1 =  0

                  if( chag(i,icut2)   .eq. yen .and.
     &              ( chag(i,icut2+1) .eq. yen .or.
     &                chag(i,icut2+1) .eq. '(' .or.
     &                chag(i,icut2+1) .eq. ')' .or.
     &              ( numec(chag(i,icut2+1))   .and.
     &                numec(chag(i,icut2+2)) .and.
     &                numec(chag(i,icut2+3)) ) ) ) then

                        idec1 = 1

                  else if( chag(i,icut2-1) .eq. yen .and.
     &                numec(chag(i,icut2))   .and.
     &                numec(chag(i,icut2+1)) .and.
     &                numec(chag(i,icut2+2)) )  then

                        idec1 = 2

                  else if(  chag(i,icut2-2) .eq. yen .and.
     &                numec(chag(i,icut2-1)) .and.
     &                numec(chag(i,icut2))   .and.
     &                numec(chag(i,icut2+1)) ) then

                        idec1 = 3

                  else if(  chag(i,icut2-3) .eq. yen .and.
     &                numec(chag(i,icut2-2)) .and.
     &                numec(chag(i,icut2-1)) .and.
     &                numec(chag(i,icut2)) )  then

                        idec1 = 4

                  end if

                     icut = icut2 - idec1

                  if( icont .eq. 1 ) then

                     write(jhf,'(''/str'',i3.3,'' ('',100a1)')
     &               i, ( chag(i,n), n = 1, icut ), yen

                     icut1 = icut1 + 8

                     goto 300

                  else if( ilt(i) - icut2 + icut1 .gt. icut1 + 3 ) then

                     write(jhf,'(100A1)')
     &               ( chag(i,n), n = icut2 - icut1 + 1, icut ),
     &                 yen

                     goto 300

                  else

                     write(jhf,'(100A1)')
     &               ( chag(i,n), n = icut2 - icut1 + 1, ilt(i) ),
     &                 cendef

                  end if

               else

                     write(jhf,'(''/str'',i3.3,'' ('',100a1)')
     &               i, ( chag(i,n), n = 1, ilt(i) ), cendef

               end if

*-----------------------------------------------------------------------

      return
      end

