
      subroutine calc_stdev(m,X,sig,A,B,conv)

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /taliin/ rsouin, nzztin, nrgnin
        common /cparm/  maxbch, maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
        common /restart/ resc2(itlmax), resc3(itlmax)
        common /stat / istdev, irestart, ireschk

        common /tcntl/ icntl, inucr

        common /sumtal01/ isistdev, ismaxcas, rsijklst
        integer   :: isistdev(itlmax), ismaxcas(itlmax)
        dimension rsijklst(itlmax)

        common /talout/ itall
        common /mpi00/ npe, me

        common /tall21/ rtfac(itlmax)  ! kitamura23/04/30

*-----------------------------------------------------------------------

        if ( icntl.eq.13 .or. icntl.eq.17 ) then

           istdev = isistdev(m)
           maxcas = ismaxcas(m)

        end if

*-----------------------------------------------------------------------

        if ( istdev .eq. 2 ) then
          c2 = rsouin + resc2(m)
          c3 = dble(nobch) * dble(maxcas) + resc3(m)
        else
          c2 = rsouin / maxcas + resc2(m)
          c3 = dble(nobch) + resc3(m)
        end if

        if ( itall.eq.4 .and. istdev.eq.2 ) then
          if ( npe .gt. 1 ) then
            c2 = rsouin/dble(nobch/(npe-1))+resc2(m)
            c3 = dble(npe-1)*dble(maxcas)+resc3(m)
          else
            c2 = rsouin/dble(nobch)+resc2(m)
            c3 = dble(maxcas)+resc3(m)
          end if
        end if

*-----------------------------------------------------------------------

        X = A / c2

        if (c3 .gt. 1) then
          rer = ( B / ((c2/c3)**2) - c3*(X**2)) / (c3**2 - c3)
          if ( rer .lt. 0 ) rer = 0
        else
          rer = X ** 2
        endif

        sig = sqrt(rer) / X
        X = X * conv

        if(rtfac(m) .lt. 0.d0) X = X * c2  ! kitamura23/04/30

*-----------------------------------------------------------------------

      end subroutine


      subroutine invert_stdev(m,A,B,X,sig,conv)

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /restart/ resc2(itlmax), resc3(itlmax)

        common /tall21/ rtfac(itlmax)  ! kitamura23/04/30

*-----------------------------------------------------------------------

        c2 = resc2(m) !! c2 means W.
        c3 = resc3(m) !! c3 means N.

*-----------------------------------------------------------------------

        Xc = X * conv


        if( rtfac(m) .ge. 0.d0 ) then

           rer = ( Xc * sig ) ** 2
           B = ( rer * (c3**2 - c3) + c3 * Xc**2 ) * ((c2/c3)**2)
           A = Xc * c2

        else

           rer = ( Xc/c2 * sig ) ** 2
           B = ( rer * (c3**2 - c3) + c3 * (Xc/c2)**2 ) * ((c2/c3)**2)
           A = Xc

        end if

*-----------------------------------------------------------------------

      end subroutine

      subroutine calc_heat_stdev(m,X,sig,A,B,conv)

        use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
cABE 2015/12/03, add tall09
        common /tall21/ rtfac(itlmax)

        common /taliin/ rsouin, nzztin, nrgnin
        common /cparm/  maxbch, maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
        common /restart/ resc2(itlmax), resc3(itlmax)
        common /stat / istdev, irestart, ireschk

        common /tcntl/ icntl, inucr

        common /sumtal01/ isistdev, ismaxcas, rsijklst
        integer   :: isistdev(itlmax), ismaxcas(itlmax)
        dimension rsijklst(itlmax)

        common /talout/ itall
        common /mpi00/ npe, me

*-----------------------------------------------------------------------

        if ( icntl .eq. 13 ) then

           istdev = isistdev(m)
           maxcas = ismaxcas(m)

        end if

*-----------------------------------------------------------------------

        if ( istdev .eq. 2 ) then
          c2 = rsouin + resc2(m)
          c3 = dble(nobch) * dble(maxcas) + resc3(m)
        else if ( istdev .ne. 2 .and. itout(m) .ge. 4) then
          c2 = rsouin + resc2(m) * maxcas
          c3 = ( dble(nobch) + resc3(m) ) * maxcas
        else
          c2 = rsouin / maxcas + resc2(m)
          c3 = dble(nobch) + resc3(m)
        end if

        if ( itall.eq.4 .and. istdev.eq.2 ) then
          if ( npe .gt. 1 ) then
            c2 = rsouin/dble(nobch/(npe-1))+resc2(m)
            c3 = dble(npe-1)*dble(maxcas)+resc3(m)
          else
            c2 = rsouin/dble(nobch)+resc2(m)
            c3 = dble(maxcas)+resc3(m)
          end if
        end if

*-----------------------------------------------------------------------


        X = A / c2


          if (c3 .gt. 1) then
            rer = ( B / ((c2/c3)**2) - c3*(X**2)) / (c3**2 - c3)
            if ( rer .lt. 0 ) rer = 0
          else
            rer = X ** 2
          endif
          sig = sqrt(rer) / X

          X = X * conv  ! kitamura23/03/31

        if(rtfac(m) .lt. 0.d0) X = X * c2  ! kitamura23/04/30


*-----------------------------------------------------------------------

      end subroutine

      subroutine invert_heat_stdev(m,A,B,X,sig,conv)

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall21/ rtfac(itlmax)

        common /cparm/  maxbch, maxcas

        common /restart/ resc2(itlmax), resc3(itlmax)
        common /stat / istdev, irestart, ireschk

        common /tcntl/ icntl, inucr

        common /sumtal01/ isistdev, ismaxcas, rsijklst
        integer   :: isistdev(itlmax), ismaxcas(itlmax)
        dimension rsijklst(itlmax)

*-----------------------------------------------------------------------

        if ( icntl .eq. 13 ) then

           istdev = isistdev(m)
           maxcas = ismaxcas(m)

        end if

*-----------------------------------------------------------------------

        c2 = resc2(m) !! c2 means W.
        c3 = resc3(m) !! c3 means N.

        if ( istdev .ne. 2 .and. itout(m) .ge. 4) then
          c2 = resc2(m) * maxcas
          c3 = resc3(m) * maxcas
        end if

*-----------------------------------------------------------------------


        Xc = X * conv


        if( rtfac(m) .ge. 0.d0 ) then

           rer = ( Xc * sig ) ** 2
           B = ( rer * (c3**2 - c3) + c3 * Xc**2 ) * ((c2/c3)**2)
           A = Xc * c2

        else

           rer = ( Xc/c2 * sig ) ** 2
           B = ( rer * (c3**2 - c3) + c3 * (Xc/c2)**2 ) * ((c2/c3)**2)
           A = Xc

        end if

*-----------------------------------------------------------------------

      end subroutine


      subroutine calc_deposit_stdev(m,X,sig,A,B,conv,ip) ! S.Abe 2015/12/03, add ip

        use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
cABE 2015/12/03, add tall09

        common /taliin/ rsouin, nzztin, nrgnin
        common /cparm/  maxbch, maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
        common /restart/ resc2(itlmax), resc3(itlmax)
        common /stat / istdev, irestart, ireschk

        common /tcntl/ icntl, inucr

        common /sumtal01/ isistdev, ismaxcas, rsijklst
        integer   :: isistdev(itlmax), ismaxcas(itlmax)
        dimension rsijklst(itlmax)

        common /talout/ itall
        common /mpi00/ npe, me

        common /tall21/ rtfac(itlmax)  ! kitamura23/04/30

*-----------------------------------------------------------------------

        if ( icntl .eq. 13 ) then

           istdev = isistdev(m)
           maxcas = ismaxcas(m)

        end if

*-----------------------------------------------------------------------

        if ( istdev .eq. 2 ) then
          c2 = rsouin + resc2(m)
          c3 = dble(nobch) * dble(maxcas) + resc3(m)
        else
          if ( itout(m) .ge. 2 .or. itals(m) .eq. 34 ) then
            c2 = rsouin + resc2(m) * maxcas
            c3 = ( dble(nobch) + resc3(m) ) * maxcas
          else
            c2 = rsouin / maxcas + resc2(m)
            c3 = dble(nobch) + resc3(m)
          endif
        end if

        if ( itall.eq.4 .and. istdev.eq.2 ) then
          if ( npe .gt. 1 ) then
            c2 = rsouin/dble(nobch/(npe-1))+resc2(m)
            c3 = dble(npe-1)*dble(maxcas)+resc3(m)
          else
            c2 = rsouin/dble(nobch)+resc2(m)
            c3 = dble(maxcas)+resc3(m)
          end if
        end if

*-----------------------------------------------------------------------

        X = A / c2
        if ( (itout(m).ge.2 .or. itals(m) .eq. 34) .and.
     &        itpat(m,ip,1) .eq. 20 ) then

          sig = 1.0 / sqrt(A)  ! probability mode

        else

          if (c3 .gt. 1) then
            rer = ( B / ((c2/c3)**2) - c3*(X**2)) / (c3**2 - c3)
            if ( rer .lt. 0 ) rer = 0
          else
            rer = X ** 2
          endif

          sig = sqrt(rer) / X

        endif

        X = X * conv

        if(rtfac(m) .lt. 0.d0) X = X * c2  ! kitamura23/04/30

*-----------------------------------------------------------------------

      end subroutine


      subroutine invert_deposit_stdev(m,A,B,X,sig,conv,ip)

        use partmod, only: itpat ! T.Sato 2022/09/09

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)

        common /tall21/ rtfac(itlmax)  ! kitamura23/04/30

*-----------------------------------------------------------------------

        common /restart/ resc2(itlmax), resc3(itlmax)

        common /cparm/ maxbch, maxcas
        common /stat / istdev, irestart, ireschk

        common /tcntl/ icntl, inucr

        common /sumtal01/ isistdev, ismaxcas, rsijklst
        integer   :: isistdev(itlmax), ismaxcas(itlmax)
        dimension rsijklst(itlmax)

*-----------------------------------------------------------------------

        if ( icntl .eq. 13 ) then

           istdev = isistdev(m)
           maxcas = ismaxcas(m)

        end if

*-----------------------------------------------------------------------

        if( istdev .ne. 2 .and.
     &     ( itout(m) .ge. 2 .or. itals(m) .eq. 34 )
     &    .and. itpat(m,ip,1) .eq. 20 ) then  ! T.Sato 2022/09/09

          c2 = resc2(m) * maxcas !! c2 means W.
          c3 = resc3(m) * maxcas !! c3 means N.
        else
          c2 = resc2(m) !! c2 means W.
          c3 = resc3(m) !! c3 means N.
        endif

*-----------------------------------------------------------------------

        Xc = X * conv

        if ( (itout(m).ge.2 .or. itals(m) .eq. 34) .and.
     &        itpat(m,ip,1) .eq. 20 ) then  ! T.Sato 2022/09/09
         if(sig.gt.0.0) then
          B = (1.0d0/sig)**2
         else
          B = 0.0
         endif
        else

           if( rtfac(m) .ge. 0.d0 ) then
             rer = ( Xc * sig ) ** 2
             B = ( rer * (c3**2 - c3) + c3 * Xc**2 ) * ((c2/c3)**2)
           else
             rer = ( Xc/c2 * sig ) ** 2
             B = ( rer * (c3**2 - c3) + c3 * (Xc/c2)**2 ) * ((c2/c3)**2)
           end if
        endif

        if( rtfac(m) .ge. 0.d0 ) then
           A = Xc * c2
        else
           A = Xc
        endif

*-----------------------------------------------------------------------

      end subroutine

      subroutine calc_prob_stdev(m,X,sig,A,B,conv)

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /taliin/ rsouin, nzztin, nrgnin
        common /cparm/  maxbch, maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
        common /restart/ resc2(itlmax), resc3(itlmax)
        common /stat / istdev, irestart, ireschk

        common /tall71/ itmorp(itlmax),itname(itlmax)

        common /tcntl/ icntl, inucr

        common /sumtal01/ isistdev, ismaxcas, rsijklst
        integer   :: isistdev(itlmax), ismaxcas(itlmax)
        dimension rsijklst(itlmax)

        common /talout/ itall
        common /mpi00/ npe, me

        common /tall21/ rtfac(itlmax)  ! kitamura23/04/30

*-----------------------------------------------------------------------
        if( icntl .eq. 13 ) then

           istdev = isistdev(m)
           maxcas = ismaxcas(m)

        endif

*-----------------------------------------------------------------------

        if( istdev .eq. 2 ) then

           c2 = rsouin + resc2(m)
           c3 = dble(nobch) * dble(maxcas) + resc3(m)

        else

           if( itmorp(m) .eq. 0 ) then

              c2 = rsouin / dble(maxcas) + resc2(m)
              c3 = dble(nobch) + resc3(m)

           elseif( itmorp(m) .ne. 0 ) then

              c2 = rsouin + resc2(m) * dble(maxcas)
              c3 = ( dble(nobch) + resc3(m) ) * dble(maxcas)

           endif

        endif

cfrtati 2021/03/06, copy & paste to here by Ogawa & Sato
        if ( itall.eq.4 .and. istdev.eq.2 ) then
          if ( npe .gt. 1 ) then
            c2 = rsouin/dble(nobch/(npe-1))+resc2(m)
            c3 = dble(npe-1)*dble(maxcas)+resc3(m)
          else
            c2 = rsouin/dble(nobch)+resc2(m)
            c3 = dble(maxcas)+resc3(m)
          end if
        end if

*-----------------------------------------------------------------------

        X = A / c2
        if( itmorp(m) .eq. 0 ) then

           if( c3 .gt. 1 ) then

              rer = ( B / ((c2/c3)**2) - c3*(X**2)) / (c3**2 - c3)
              if( rer .lt. 0.d0 ) rer = 0.d0

           else

              rer = X ** 2

           endif

           sig = dsqrt(rer) / X

        else

           sig = 1.d0 / dsqrt(B)

        endif

        X = X * conv

        if(rtfac(m) .lt. 0.d0) X = X * c2  ! kitamura23/04/30

*-----------------------------------------------------------------------

      end subroutine


      subroutine invert_prob_stdev(m,A,B,X,sig,conv)

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /cparm/ maxbch, maxcas
        common /restart/ resc2(itlmax), resc3(itlmax)
        common /stat / istdev, irestart, ireschk

        common /tall71/ itmorp(itlmax),itname(itlmax)

        common /tcntl/ icntl, inucr

        common /sumtal01/ isistdev, ismaxcas, rsijklst
        integer   :: isistdev(itlmax), ismaxcas(itlmax)
        dimension rsijklst(itlmax)

        common /tall21/ rtfac(itlmax)  ! kitamura23/04/30

*-----------------------------------------------------------------------

        if( icntl .eq. 13 ) then

           istdev = isistdev(m)
           maxcas = ismaxcas(m)

        endif

*-----------------------------------------------------------------------

        if( istdev .eq. 2 ) then
           c2 = resc2(m) !! c2 means W.
           c3 = resc3(m) !! c3 means N.
        else
           if( itmorp(m) .eq. 0 ) then
              c2 = resc2(m) !! c2 means W.
              c3 = resc3(m) !! c3 means N.
           elseif( itmorp(m) .ne. 0 ) then
              c2 = resc2(m) * maxcas !! c2 means W.
              c3 = resc3(m) * maxcas !! c3 means N.
           endif
        endif

*-----------------------------------------------------------------------

        Xc = X * conv

        if( itmorp(m) .eq. 0 ) then

           if( rtfac(m) .ge. 0.d0 ) then
             rer = ( Xc * sig ) ** 2
             B = ( rer * (c3**2 - c3) + c3 * Xc**2 ) * ((c2/c3)**2)
           else
             rer = ( Xc/c2 * sig ) ** 2
             B = ( rer * (c3**2 - c3) + c3 * (Xc/c2)**2 ) * ((c2/c3)**2)
           end if
        else
         if(sig.gt.0.0) then
          B = (1.0d0/sig)**2
         else
          B = 0.0
         endif
        endif

        if( rtfac(m) .ge. 0.d0 ) then
           A = Xc * c2
        else
           A = Xc
        endif

*-----------------------------------------------------------------------

      end subroutine


      subroutine calc_stdev_ext(m,X1,X2,X3,X4,A1,A2,A3,A4,conv)

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /taliin/ rsouin, nzztin, nrgnin
        common /cparm/  maxbch, maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
        common /restart/ resc2(itlmax), resc3(itlmax)
        common /stat / istdev, irestart, ireschk

        common /tcntl/ icntl, inucr

        common /sumtal01/ isistdev, ismaxcas, rsijklst
        integer   :: isistdev(itlmax), ismaxcas(itlmax)
        dimension rsijklst(itlmax)

        common /talout/ itall
        common /mpi00/ npe, me

*-----------------------------------------------------------------------

        if ( icntl.eq.13 .or. icntl.eq.17 ) then

           istdev = isistdev(m)
           maxcas = ismaxcas(m)

        end if

*-----------------------------------------------------------------------

        if ( istdev .eq. 2 ) then
          c2 = rsouin + resc2(m)
          c3 = dble(nobch) * dble(maxcas) + resc3(m)
        else
          c2 = rsouin / maxcas + resc2(m)
          c3 = dble(nobch) + resc3(m)
        end if

        if ( itall.eq.4 .and. istdev.eq.2 ) then
          if ( npe .gt. 1 ) then
            c2 = rsouin/dble(nobch/(npe-1))+resc2(m)
            c3 = dble(npe-1)*dble(maxcas)+resc3(m)
          else
            c2 = rsouin/dble(nobch)+resc2(m)
            c3 = dble(maxcas)+resc3(m)
          end if
        end if

*-----------------------------------------------------------------------

        X1 = A1 / c2

        if (c3 .gt. 1) then
          rer = ( A2 / ((c2/c3)**2) - c3*(X1**2)) / (c3**2 - c3)
          if ( rer .lt. 0 ) rer = 0
        else
          rer = X1 ** 2
        endif

        X2 = sqrt(rer) / X1
        X1 = X1 * conv

*-----------------------------------------------------------------------

        B1 = A1 / (c2/c3)
        B2 = A2 / (c2/c3)**2
        B3 = A3 / (c2/c3)**3
        B4 = A4 / (c2/c3)**4
        if (c3 .gt. 1) then
           X3= B3
           X4= B4 -4*B1*B3/c3 +8*B2*B1**2/c3**2 -4*B1**4/c3**3 -B2**2/c3
           Xsecond = (B2-B1**2/c3)**2
           if ( Xsecond .gt. 0d0 ) then
              X4= X4/Xsecond
           else
              X4= 0d0
           end if
        else
           X3 = A3
           X4 = 0d0
        end if

*-----------------------------------------------------------------------

      end subroutine
