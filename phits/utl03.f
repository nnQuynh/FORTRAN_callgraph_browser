************************************************************************
*                                                                      *
      subroutine readl(jsn,jsi,dsin,idsi,ill,ilf,cmmt,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
*                                                                      *
*       read lines of input files                                      *
*       modified by K.Niita on 2001/02/26                              *
*                                                                      *
************************************************************************
      use CHARVARMOD, only: ErrLine_Adjust,irwt
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err
      integer l_errtmp

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chins*200, chlws*200, chcms*200
      character ctmp*200
      parameter ( icolms = 200 )

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      character yen*1

      character*(*) cmmt

      integer kfind, ifind, ninit, nlast

*-----------------------------------------------------------------------

            yen  = char(92)
            isql = 0
            i3s  = 1
            i4s  = 1

*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140       ill(jsn) = ill(jsn) + 1

*-----------------------------------------------------------------------

            read(jsi,'(a200)', iostat=ios ) chin
            call ZspcReplace(chin)
C ----------------------------------
            if( ios .eq. -1 ) goto 251

*-----------------------------------------------------------------------

  240       continue

               if( ill(jsn) .le. ilf(jsn) ) goto 250

*-----------------------------------------------------------------------
*           close files
*-----------------------------------------------------------------------

  251       continue

                  call closef(jsi,jsn)

               if( jsn .le. 0 ) then

                  jpn = 3
                  return

               else if( jsn .gt. 0 ) then

                  goto 140

               end if

*-----------------------------------------------------------------------
*        one line edit
*-----------------------------------------------------------------------

  250    continue

                  call chlngt(chin,200,i1,i2)

                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,cmmt)

                  chcm = chlw
                  call chcomp(chcm,i1,i3,i4)

                  k = i1

*-----------------------------------------------------------------------
*        comment line
*        normal line  iskip = 0
*        blank line   iskip = 1
*        comment line iskip = 2
*             ( comment character or 'c ' within 5 column )
*-----------------------------------------------------------------------

            if( i1 .eq. 0 .and. i2 .eq. 0 ) then

               iskip = 1
               return

            else if( index(cmmt,chin(i1:i1)) .ne. 0 ) then

               iskip = 2
               return

            else if( chlw(i1:i1+1) .eq. 'c ' .and. i1 .le. 5  ) then

               iskip = 2
               return

            else

               iskip = 0

            end if

*-----------------------------------------------------------------------
*        sequential line
*-----------------------------------------------------------------------

            if( chlw(i3:i3) .eq. yen ) then

               isql = isql + 1

  150          continue

               if( i3s + i3 - i1 + 1 .gt. 200 ) goto 998

               chin(i3:i3) = ' '
               chlw(i3:i3) = ' '
               chcm(i4:i4) = ' '

               chins(i3s:i3s+i3-i1) = chin(i1:i3)
               chlws(i3s:i3s+i3-i1) = chlw(i1:i3)
               chcms(i4s:i4s+i4-i1) = chlw(i1:i4)

               i3s = i3s + i3 - i1 + 1
               i4s = i4s + i4 - i1 + 1

                  ill(jsn) = ill(jsn) + 1

                  read(jsi,'(a200)', iostat = ios ) chin
                  if( ios .eq. -1 ) goto 999

                  if( ill(jsn) .gt. ilf(jsn) ) goto 999

                  call chlngt(chin,200,i1,i2)

                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,cmmt)

                  chcm = chlw
                  call chcomp(chcm,i1,i3,i4)

                  if( chlw(i3:i3) .eq. yen ) goto 150

               if( i3s + i3 - i1 + 1 .gt. 200 ) goto 998

               chins(i3s:i3s+i3-i1) = chin(i1:i3)
               chlws(i3s:i3s+i3-i1) = chlw(i1:i3)
               chcms(i4s:i4s+i4-i1) = chlw(i1:i4)

               i3s = i3s + i3 - i1 + 1
               i4s = i4s + i4 - i1 + 1

            end if

            if( isql .gt. 0  ) then

               i1 = 1
               i2 = i3s - 1
               i3 = i3s - 1
               i4 = i4s - 1

               chin(i1:i3) = chins(i1:i3)
               chlw(i1:i3) = chlws(i1:i3)
               chcm(i1:i4) = chcms(i1:i4)

               k = i1

            end if

*-----------------------------------------------------------------------
*        infl: {filename} [1-15]  ; include file
*-----------------------------------------------------------------------

            if( jpn .ne. 2 .and. jpn .ne. 4
     &           .and. chlw(i1:i1+4) .eq. 'infl:' ) then

                  k = k + 4

c  *** T.Sato 2014/2/20, avoid error when -, [, or ] are used in the comment
               do iii=i3+1,i2
                chin(iii:iii)=' '
               enddo
c  *****************************

               call inclf(jsn,jsi,dsin,chin,k,200,ill,ilf,idsi,ierr)

                  if( ierr .ne. 0 ) return

                  goto 140

            else if( jpn .eq. 2 .and. chlw(i1:i1+4) .eq. 'infl:' ) then

               kfind = i1
               ifind = 0
               ninit = i1
               nlast = i3
               do while ( kfind .le. i3 .and. ifind .ne. 1 )
                  kfind = kfind + 1
                  if ( chlw(kfind:kfind) .eq. '{' ) then
                     ninit = kfind + 1
                  else if ( chlw(kfind:kfind) .eq. '}' ) then
                     nlast = kfind - 1
                     ifind = 1
                  end if
               end do

               if ( ifind .eq. 1 ) then
                  if( irwt /= 0) then
                   l_errtmp = ill(jsn)
                   call ErrLine_Adjust(ctmp,ill(jsn))
                  write(*,*)
     &                 '*** Warning: PHITS skipped infl:{ ',
     &                 trim(adjustl(chlw(ninit:nlast))),' } at line',
     &                 ill(jsn),' of ',trim(ctmp)
                   ill(jsn) = l_errtmp
                  else
                  write(*,*)
     &                 '*** Warning: PHITS skipped infl:{ ',
     &                 trim(adjustl(chlw(ninit:nlast))),' } at line',
     &                 ill(jsn),' of ',dsin(jsn)(1:idsi(jsn))
                  endif

               else

                  if( irwt /= 0) then
                   l_errtmp = ill(jsn)
                   call ErrLine_Adjust(ctmp,ill(jsn))
                  write(*,*)
     &                 '*** Warning: PHITS skipped infl: at line',
     &                 ill(jsn),' of ',trim(ctmp)
                   ill(jsn) = l_errtmp
                  else
                  write(*,*)
     &                 '*** Warning: PHITS skipped infl: at line',
     &                 ill(jsn),' of ',dsin(jsn)(1:idsi(jsn))
                  endif

               end if

            else if( jpn .eq. 4 .and. chlw(i1:i1+4) .eq. 'infl:' ) then

cFURUTA20200515 infl: is skipped after fill *:* *:* *:* by ivoxel=1

               kfind = i1
               ifind = 0
               ninit = i1
               nlast = i3
               do while ( kfind .le. i3 .and. ifind .ne. 1 )
                  kfind = kfind + 1
                  if ( chlw(kfind:kfind) .eq. '{' ) then
                     ninit = kfind + 1
                  else if ( chlw(kfind:kfind) .eq. '}' ) then
                     nlast = kfind - 1
                     ifind = 1
                  end if
               end do

               if ( ifind .eq. 1 ) then

                  if( irwt /= 0) then
                   l_errtmp = ill(jsn)
                   call ErrLine_Adjust(ctmp,ill(jsn))
                  write(*,*)
     &                 '*** infl:{ ',
     &                 trim(adjustl(chlw(ninit:nlast))),' } at line',
     &               ill(jsn),' of ',trim(ctmp),
     &               ' is skipped by ivoxel=1 or icells=1'
                   ill(jsn) = l_errtmp
                  else
                  write(*,*)
     &                 '*** infl:{ ',
     &                 trim(adjustl(chlw(ninit:nlast))),' } at line',
     &               ill(jsn),' of ',dsin(jsn)(1:idsi(jsn)),
     &               ' is skipped by ivoxel=1 or icells=1'
                  endif
               else

                  if( irwt /= 0) then
                   l_errtmp = ill(jsn)
                   call ErrLine_Adjust(ctmp,ill(jsn))
                  write(*,*)
     &                 '*** infl: at line',
     &                 ill(jsn),' of ',trim(ctmp),
     &               ' is skipped by ivoxel=1 or icells=1'
                   ill(jsn) = l_errtmp
                  else
                  write(*,*)
     &                 '*** infl: at line',
     &                 ill(jsn),' of ',dsin(jsn)(1:idsi(jsn)),
     &               ' is skipped by ivoxel=1 or icells=1'
                  endif
               end if

               goto 140

            end if

*-----------------------------------------------------------------------
*        set: c1[5.0] c23[pi*4.0]  ; set constants
*-----------------------------------------------------------------------

            if( jpn .ne. 2 .and. chlw(i1:i1+3) .eq. 'set:' ) then

                  k = k + 3

               call setcvp(chlw,k,icolms,dsin,idsi,ill,jsn,ierr)

                  if( ierr .ne. 0 ) return

                  goto 140

            else if( jpn .eq. 2 .and. chlw(i1:i1+3) .eq. 'set:' ) then

             if( irwt /= 0) then
             l_errtmp = ill(jsn)
              call ErrLine_Adjust(ctmp,ill(jsn))
               write(*,*)
     &              '*** Warning: PHITS skipped set: at line',
     &              ill(jsn),' of ',trim(ctmp)
             ill(jsn) = l_errtmp
             else
               write(*,*)
     &              '*** Warning: PHITS skipped set: at line',
     &              ill(jsn),' of ',dsin(jsn)(1:idsi(jsn))
             endif

            end if

*-----------------------------------------------------------------------
*       stop reading or stop reading the sections
*-----------------------------------------------------------------------

            if( jpn .le. 1 ) then

               if( chlw(i1:i1+1) .eq. 'q:' ) then

                  jsn0 = jsn

                  do j = jsn0, 1, -1

                     call closef(jsi,jsn)

                  end do

                     jpn = 3
                     return

               else if( chlw(i1:i1+2) .eq. 'qp:' ) then

                     jpn = 2
                     goto 140

               end if

            end if

*-----------------------------------------------------------------------
*       skip lines under qp: up to [sections]
*-----------------------------------------------------------------------

         if( jpn .eq. 2 ) then

            if( i1 .gt. 5 .or. chlw(i1:i1) .ne. '[' ) goto 140

            jpn = 1

         end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Sequential line include too more characters: max(200)'
         ErrCha = ''
         ErrID = 'L:411/R:readl/F:utl03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Sequential line meets the end of file.'
         ErrCha = ''
         ErrID = 'L:423/R:readl/F:utl03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine readm(jsn,jsi,dsin,idsi,ill,ilf,cmmt,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
*                                                                      *
*       read lines of input files for only [material]                  *
*       in order to avoid "c " within 5 column                         *
*       modified by K.Niita on 2016/11/25                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chins*200, chlws*200, chcms*200
      parameter ( icolms = 200 )

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      character yen*1

      character*(*) cmmt

*-----------------------------------------------------------------------

            yen  = char(92)
            isql = 0
            i3s  = 1
            i4s  = 1

*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140       ill(jsn) = ill(jsn) + 1

*-----------------------------------------------------------------------

            read(jsi,'(a200)', iostat=ios ) chin
            call ZspcReplace(chin)
C ----------------------------------
            if( ios .eq. -1 ) goto 251

*-----------------------------------------------------------------------

  240       continue

               if( ill(jsn) .le. ilf(jsn) ) goto 250

*-----------------------------------------------------------------------
*           close files
*-----------------------------------------------------------------------

  251       continue

                  call closef(jsi,jsn)

               if( jsn .le. 0 ) then

                  jpn = 3
                  return

               else if( jsn .gt. 0 ) then

                  goto 140

               end if

*-----------------------------------------------------------------------
*        one line edit
*-----------------------------------------------------------------------

  250    continue

                  call chlngt(chin,200,i1,i2)

                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,cmmt)

                  chcm = chlw
                  call chcomp(chcm,i1,i3,i4)

                  k = i1

*-----------------------------------------------------------------------
*        comment line
*        normal line  iskip = 0
*        blank line   iskip = 1
*        comment line iskip = 2
*             ( comment character or 'c ' within 5 column )
*-----------------------------------------------------------------------

            if( i1 .eq. 0 .and. i2 .eq. 0 ) then

               iskip = 1
               return

            else if( index(cmmt,chin(i1:i1)) .ne. 0 ) then

               iskip = 2
               return



            else

               iskip = 0

            end if

*-----------------------------------------------------------------------
*        sequential line
*-----------------------------------------------------------------------

            if( chlw(i3:i3) .eq. yen ) then

               isql = isql + 1

  150          continue

               if( i3s + i3 - i1 + 1 .gt. 200 ) goto 998

               chin(i3:i3) = ' '
               chlw(i3:i3) = ' '
               chcm(i4:i4) = ' '

               chins(i3s:i3s+i3-i1) = chin(i1:i3)
               chlws(i3s:i3s+i3-i1) = chlw(i1:i3)
               chcms(i4s:i4s+i4-i1) = chlw(i1:i4)

               i3s = i3s + i3 - i1 + 1
               i4s = i4s + i4 - i1 + 1

                  ill(jsn) = ill(jsn) + 1

                  read(jsi,'(a200)', iostat = ios ) chin
                  if( ios .eq. -1 ) goto 999

                  if( ill(jsn) .gt. ilf(jsn) ) goto 999

                  call chlngt(chin,200,i1,i2)

                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,cmmt)

                  chcm = chlw
                  call chcomp(chcm,i1,i3,i4)

                  if( chlw(i3:i3) .eq. yen ) goto 150

               if( i3s + i3 - i1 + 1 .gt. 200 ) goto 998

               chins(i3s:i3s+i3-i1) = chin(i1:i3)
               chlws(i3s:i3s+i3-i1) = chlw(i1:i3)
               chcms(i4s:i4s+i4-i1) = chlw(i1:i4)

               i3s = i3s + i3 - i1 + 1
               i4s = i4s + i4 - i1 + 1

            end if

            if( isql .gt. 0  ) then

               i1 = 1
               i2 = i3s - 1
               i3 = i3s - 1
               i4 = i4s - 1

               chin(i1:i3) = chins(i1:i3)
               chlw(i1:i3) = chlws(i1:i3)
               chcm(i1:i4) = chcms(i1:i4)

               k = i1

            end if

*-----------------------------------------------------------------------
*        infl: {filename} [1-15]  ; include file
*-----------------------------------------------------------------------

            if( jpn .ne. 2. and. chlw(i1:i1+4) .eq. 'infl:' ) then

                  k = k + 4

c  *** T.Sato 2014/2/20, avoid error when -, [, or ] are used in the comment
               do iii=i3+1,i2
                chin(iii:iii)=' '
               enddo
c  *****************************

               call inclf(jsn,jsi,dsin,chin,k,200,ill,ilf,idsi,ierr)

                  if( ierr .ne. 0 ) return

                  goto 140

            end if

*-----------------------------------------------------------------------
*        set: c1[5.0] c23[pi*4.0]  ; set constants
*-----------------------------------------------------------------------

            if( jpn .ne. 2 .and. chlw(i1:i1+3) .eq. 'set:' ) then

                  k = k + 3

               call setcvp(chlw,k,icolms,dsin,idsi,ill,jsn,ierr)

                  if( ierr .ne. 0 ) return

                  goto 140

            end if

*-----------------------------------------------------------------------
*       stop reading or stop reading the sections
*-----------------------------------------------------------------------

            if( jpn .le. 1 ) then

               if( chlw(i1:i1+1) .eq. 'q:' ) then

                  jsn0 = jsn

                  do j = jsn0, 1, -1

                     call closef(jsi,jsn)

                  end do

                     jpn = 3
                     return

               else if( chlw(i1:i1+2) .eq. 'qp:' ) then

                     jpn = 2
                     goto 140

               end if

            end if

*-----------------------------------------------------------------------
*       skip lines under qp: up to [sections]
*-----------------------------------------------------------------------

         if( jpn .eq. 2 ) then

            if( i1 .gt. 5 .or. chlw(i1:i1) .ne. '[' ) goto 140

            jpn = 1

         end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Sequential line include too more characters: max(200)'
         ErrCha = ''
         ErrID = 'L:716/R:readm/F:utl03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Sequential line meets the end of file.'
         ErrCha = ''
         ErrID = 'L:728/R:readm/F:utl03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine setcvp(chlw,ic0,iclm,dsin,idsi,ill,jsn,ierr)
*                                                                      *
*      purpose  : preprocedure of setcv from phits to angel            *
*                 change chlw to dum                                   *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character chlw*(*)
      dimension ill(0:9)

      character dsin(0:9)*200
      dimension idsi(0:9)

*-----------------------------------------------------------------------

            do i = ic0, iclm

                  dum(i) = chlw(i:i)

            end do

               call setcv(dum,ic0,iclm,dsin,idsi,ill,jsn,ierr)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine snum(chlw,ic0,ic1,ic2,prn,ierr)
*                                                                      *
*        purpose : character to real number                            *
*                  one number from ic0 to ic1 by separator ' '         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character chlw*(*)

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

            ic   = ic0
            icl  = ic1

*-----------------------------------------------------------------------

            do i = ic, icl

               if( chlw(i:i) .ne. ' ' .and.
     &             chlw(i:i) .ne. tub ) goto 100

            end do

               ierr = 1
               return

  100          ic = i

*-----------------------------------------------------------------------

               ii = 0

            do i = ic, icl

               ii = ii + 1

                  dum(ii+1) = chlw(i:i)

            end do

              jc  = 2
              jcl = ii + 1

*-----------------------------------------------------------------------

         if( chlw(ic:ic) .ne. '[' .and.
     &       chlw(ic:ic) .ne. '{' ) then

            do i = ic + 1, icl

               if( chlw(i:i) .eq. ' ' .or.
     &             chlw(i:i) .eq. tub ) goto 200

            end do

                  i = jcl + 1 + ic - 2

  200             jcl = i - ic + 2
                  jc  = 1

                  dum(jc)  = '['
                  dum(jcl) = ']'

         end if

*-----------------------------------------------------------------------

               call pnum(dum,jc,jcl,prn,ierr)

               ic2 = jc + ic - 2 + 1

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine unum(chlw,ic0,ic1,ic2,prn,ierr)
*                                                                      *
*        purpose : character to real number                            *
*                  one number from ic0 to ic1 by separator             *
*                  ),},] or non numerical character dnen3              *
*                  this unum neglect the operaters                     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      logical dnen3

      character dum(ichrl)*1
      character chlw*(*)

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

            ic   = ic0
            icl  = ic1

*-----------------------------------------------------------------------

            do i = ic, icl

               if( chlw(i:i) .ne. ' ' .and.
     &             chlw(i:i) .ne. tub ) goto 100

            end do

               ierr = 1
               return

  100          ic = i

*-----------------------------------------------------------------------

               ii = 0

            do i = ic, icl

               ii = ii + 1

                  dum(ii+1) = chlw(i:i)

            end do

              jc  = 2
              jcl = ii + 1

*-----------------------------------------------------------------------

               ipars = 0

         if( chlw(ic:ic) .ne. '[' .and.
     &       chlw(ic:ic) .ne. '{' ) then

               ipars = 1

            do i = ic + 1, icl

               if( chlw(i:i) .eq. ' ' .or.
     &             chlw(i:i) .eq. tub ) goto 200

               if( chlw(i:i) .eq. ']' .or.
     &             chlw(i:i) .eq. '}' ) goto 200

               if( dnen3( chlw(i:i) ) ) goto 200

            end do

                  i = jcl + 1 + ic - 2

  200             jcl = i - ic + 2
                  jc  = 1

                  dum(jc)  = '['
                  dum(jcl) = ']'

         end if

*-----------------------------------------------------------------------

               call pnum(dum,jc,jcl,prn,ierr)

               ic2 = jc + ic - 2 + 1

               if( ipars .eq. 1 ) ic2 = ic2 - 1

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine onum(chlw,ic0,ic1,prn,ierr)
*                                                                      *
*        purpose : character to real number                            *
*                  one number from ic0 to ic1                          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character chlw*(*)

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

            ic   = ic0
            icl  = ic1

*-----------------------------------------------------------------------

            do i = ic, icl

               if( chlw(i:i) .ne. ' ' .and.
     &             chlw(i:i) .ne. tub ) goto 100

            end do

               ierr = 1
               return

  100          ic = i

*-----------------------------------------------------------------------

               ii = 0

            do i = ic, icl

               ii = ii + 1

                  dum(ii+1) = chlw(i:i)

            end do

              jc  = 2
              jcl = ii + 1

*-----------------------------------------------------------------------

         if( chlw(ic:ic) .ne. '[' .and.
     &       chlw(ic:ic) .ne. '{' ) then

                  jc  = 1
                  jcl = jcl + 1

                  dum(jc)  = '['
                  dum(jcl) = ']'

         end if

*-----------------------------------------------------------------------

               call pnum(dum,jc,jcl,prn,ierr)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gnum(chlw,ic0,ic1,prn,xval,ierr)
*                                                                      *
*        purpose : function to real number                             *
*                  function from ic0 to ic1                            *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character chlw*(*)
      character c2*1

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

            prn  = 0.0
            ierr = 1

*-----------------------------------------------------------------------

            ic   = ic0
            icl  = ic1

*-----------------------------------------------------------------------

            do i = ic, icl

               if( chlw(i:i) .ne. ' ' .and.
     &             chlw(i:i) .ne. tub ) goto 100

            end do

               ierr = 1
               return

  100          ic = i

*-----------------------------------------------------------------------

               ii = 0

            do i = ic, icl

               ii = ii + 1

                  dum(ii+1) = chlw(i:i)

            end do

              jc  = 2
              jcl = ii + 1

*-----------------------------------------------------------------------

         if( chlw(ic:ic) .ne. '[' .and.
     &       chlw(ic:ic) .ne. '{' ) then

                  jc  = 1
                  jcl = jcl + 1

                  dum(jc)  = '['
                  dum(jcl) = ']'

         end if

*-----------------------------------------------------------------------

         if( dum(jc) .eq. '[' .or. dum(jc) .eq. '{' ) then

            if( dum(jc) .eq. '[' ) then

               c2 = ']'

            else if( dum(jc) .eq. '{' ) then

               c2 = '}'

            end if

            call func12(dum,jc,ierr,jcl,c2,prn,xval)

                  if( dum(jc) .ne. c2 ) return
                  if(ierr.ne.0) return

                  ierr = 0

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine chcomp(t,i1,i3,i4)
*                                                                      *
*          Convert a text into the compress form removing the blank    *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)
      character*(*) t
      character dum*2048
      character tub*1
      tub = char(9)

         i4 = 0
         if( i1 .eq. 0 .and. i3 .eq. 0 ) return

         j = i1 - 1
      do i = i1, i3
         if( t(i:i) .ne. ' ' .and. t(i:i) .ne. tub ) then
            j = j + 1
            dum(j:j) = t(i:i)
         end if
      end do

         i4 = j

      do i = i1, i4
         t(i:i) = dum(i:i)
      end do

      if( i4 .lt. i3 ) then
         do i = i4+1, i3
            t(i:i) = ' '
         end do
      end if


      return
      end

************************************************************************
*                                                                      *
      subroutine chcaps(t,i1,i2,i3,cmmt)
*                                                                      *
*          Convert a text into lower letters.                          *
*          and delete the comment and blank in the last part           *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)
      character*(*) t
      character*1 c
      character*(*) cmmt
      character tub*1
      tub = char(9)

         i3 = 0
         if( i1 .eq. 0 .and. i2 .eq. 0 ) return

      do i = i1, i2

         if( index(cmmt,t(i:i)) .ne. 0 ) goto 100
         c = t(i:i)
         if( c .ge. 'A' .and. c .le. 'Z' )
     &       c = char(ichar(c)+ichar('a')-ichar('A'))
         t(i:i) = c
      end do

  100 continue
         i4 = i - 1
         do i = i4, i1, -1
            if( t(i:i) .ne. ' ' .and. t(i:i) .ne. tub ) goto 200
         end do

         i = i1 - 1

  200 continue

         i3 = i

      do i = i3+1, i2
         t(i:i) = ' '
      end do

      return
      end

************************************************************************
*                                                                      *
      subroutine chcptl(t,i1,i2)
*                                                                      *
*          Convert a text into upper letters.                          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)
      character*(*) t
      character*1 c

      do i = i1, i2
         c = t(i:i)
         if( c .ge. 'a' .and. c .le. 'z' )
     &       c = char(ichar(c)+ichar('A')-ichar('a'))
         t(i:i) = c
      end do

      return
      end

************************************************************************
*                                                                      *
      function knump(chin,i1,i2)
*                                                                      *
*          get column number of final ' ) '                            *
*                                                                      *
************************************************************************

      character chin*200

         knump = i2 + 1

         if( chin(i1:i1) .ne. '(' ) return

                  ibra = 0

         do i = i1, i2

            if( chin(i:i) .eq. ')' ) then

                  ibra = ibra - 1

                  if( ibra .eq. 0 ) goto 100

            else if( chin(i:i) .eq. '(' ) then

                  ibra = ibra + 1

            end if

         end do

            return

  100       knump = i

      return
      end

************************************************************************
*                                                                      *
      function inumc(chin,i1,i2,cc)
*                                                                      *
*          get column number of character cc                           *
*                                                                      *
************************************************************************

      character chin*200
      character cc*1

         do i = i1, i2

            if( chin(i:i) .eq. cc ) then

               inumc = i
               return

            end if

         end do

               inumc = i2 + 1

      return
      end

************************************************************************
*                                                                      *
      function jnumc(chin,i1,i2)
*                                                                      *
*          get column number of non space character                    *
*                                                                      *
************************************************************************

      character chin*200
      character tub*1
      tub = char(9)


         if( i1 .gt. i2 ) goto 130

         do i = i1, i2

            if( chin(i:i) .ne. ' ' .and.
     &          chin(i:i) .ne. tub ) goto 120

         end do

  130    continue

               jnumc = i2 + 1
               return

  120    continue

               jnumc = i

      return
      end

************************************************************************
*                                                                      *
      subroutine chlngt(hl,nn,icl1,icl2)
*                                                                      *
*     supply the starting and ending column number of non-blank        *
*     character. We replace the dub by 8 blank.                        *
*     hl : searching character data                                    *
*     nn : total number of character data                              *
*     icl1 : starting column                                           *
*     icl2 : ending column                                             *
*     (created by k.kosako and modified by k.niita)                    *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      character hl*(*)
      character dum*10000
      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------
*     replace tub by 8 blanks
*-----------------------------------------------------------------------


         k = 0

      do i = 1, nn

         if( hl(i:i) .eq. tub ) then

            do j = 1, 8

               k = k + 1
               dum(k:k) = ' '

            end do

         else

            k = k + 1
            dum(k:k) = hl(i:i)

         end if

      end do

      do i = 1, nn

         hl(i:i) = dum(i:i)

      end do

*-----------------------------------------------------------------------

      icl1 = 0
      icl2 = 0
      do 110 i = 1, nn
  110 if( hl(i:i) .ne. ' ' .and. hl(i:i) .ne. tub ) goto 120
      goto 900
  120 icl1 = i
      do 130 i = nn, icl1, -1
  130 if( hl(i:i) .ne. ' ' .and. hl(i:i) .ne. tub ) goto 140
  140 icl2 = i
  900 continue

      return
      end


************************************************************************
*                                                                      *
      subroutine mapnum(hl,nn,ic,icl1,icl2,ierr)
*                                                                      *
*     supply the starting and ending column number of non-blank        *
*     character. We replace the dub by 8 blank.                        *
*     hl : searching character data                                    *
*     nn : total number of character data                              *
*     icl1 : starting column                                           *
*     icl2 : ending column                                             *
*     (created by ASTOM)                                               *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      character hl*(*)
*      character dum*100000
      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

      ierr = 0
      icl1 = 0
      icl2 = 0
      do 110 i = 1, nn
  110 if( hl(i:i) .ne. ' ' .and. hl(i:i) .ne. tub
     &    .and. i .ge. ic ) goto 120
      ierr = 1
      goto 900
  120 icl1 = i
      do 130 i = 1, nn
  130 if( (hl(i:i) .eq. ' ' .or. hl(i:i) .eq. tub )
     &     .and. i .gt. icl1 ) goto 140
  140 icl2 = i
  900 continue

      return
      end

************************************************************************
      subroutine dichotomy(ix,ixmax,y,ydata) ! T.Sato 2020/08/10, Nibunho
!  ix: return value (ydata(ix) >= y)
!  ixmax: dimension size
!  y: reference y value
!  ydata: y data
      implicit double precision (a-h,o-z)

      dimension ydata(ixmax)

      if(ixmax.eq.1) then
       ix=1
       return
      endif

      ix0=1
      ix1=ixmax/2 ! initial condition
      ix2=ixmax

      do i=1,ixmax
       if(ix2-ix0.le.3) exit ! enough close
       if(ydata(ix1).ge.y) then ! already too high y value
        ix2=ix1
        ix1=(ix0+ix1)/2
       else
        ix0=ix1
        ix1=(ix1+ix2)/2
       endif
      enddo

      do ix1=ix0,ix2
       if(ydata(ix1).ge.y) exit
      enddo

      ix=min(ix1,ixmax)

      return

      end
