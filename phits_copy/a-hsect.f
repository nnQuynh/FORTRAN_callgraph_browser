************************************************************************
*                                                                      *
      subroutine hone(dum,lum,in,jsi,jsn,ill,ilf,dsin,idsi,
     &                ymax,ymin,ymin2,
     &                xmax,xmin,xmin2,noned,noner,ncom,ierr,jof,jif,
     &                iycm,rycm,ifon)
*                                                                      *
*         MH:      MAX NUMBER OF COLUMN IN ONE SECTION                 *
*         ML:      MAX NUMBER OF ROW IN ONE SECTION                    *
*                                                                      *
*         MC:      MAX NUMBER OF LINE COMMENTS (TOTAL)                 *
*                                                                      *
*         JOF:     OUTPUT FILE OF LINES                                *
*                                                                      *
*         JIL:     OUTPUT FILE OF REGION                               *
*                                                                      *
*         ILL(JSN):  PRESENT LINE                                      *
*         ILF(JSN):  FINAL LINE                                        *
*                                                                      *
*         NONED:   NUMBER OF LINES WRITTEN IN JOF                      *
*         NONER:   NUMBER OF REGIONS WRITTEN IN JIL                    *
*         NCOM:    NUMBER OF LINE COMMENTS                             *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      integer,allocatable:: ivi(:,:)
      integer,allocatable:: nevk(:,:), nevn(:,:)

      real(8),allocatable:: rkjj(:,:)
      dimension iprjj(mh), idrjj(mh)
      integer,allocatable:: kdjj(:,:), kijj(:,:)
      integer,allocatable:: ivnjj(:,:), ifpjj(:,:)
      integer,allocatable:: ivijj(:,:,:)

      character jdum(ichrl)*1

      character dum(ichrl)*1
      character lum(ichrl)*1
      dimension dfac(5),fac(mh,5)
      dimension iycm(mc,8),rycm(mc,7)
      dimension idc(mh),idyd(-1:mh),iddy(mh),jycm(mh),
     &          iylc(mh),iyss(mh),syss(mh),iylt(mh),ycol(mh,3),
     &          iysb(mh),iyhh(mh),iyyd(mh,2),iybb(mh),ycob(mh,3),
     &          ided(-1:mh),idde(mh),ides(mh),
     &          idef(-1:mh),iffy(-1:mh), iffd(-1:mh), iydd(mh)

      dimension rk(0:mfc)
      dimension ifdc(mh), ivdc(mh)
      dimension kd(0:mfc), kp(0:mfc), ki(0:mfc), jd(0:mfc)
      dimension md(0:mfc), mg(0:mfc), me(0:mfc), mp(0:mfc)
      dimension ivn(0:mfc), ifp(0:mfc)

      common /rval1/ cval(mxcval), aval(mxcval)

      character m_err*200
      common /error/ m_err, l_err, k_err

      dimension ill(0:9), ilf(0:9)
      character dsin(0:9)*200
      dimension idsi(0:9)

      logical deqn1

      dimension coln(3)

      real(8),allocatable:: daxy(:,:)
      logical isv

*-----------------------------------------------------------------------
      logical, save :: second_call

*-----------------------------------------------------------------------

      character yen*1
      character tub*1
      yen = char(92)
      tub = char(9)

*-----------------------------------------------------------------------

      allocate(ivi(0:mfc,0:mfc),
     &         nevk(0:mfc,0:mfc), nevn(0:mfc,0:mfc),
     &         rkjj(mh,0:mfc),
     &         kdjj(mh,0:mfc), kijj(mh,0:mfc),
     &         ivnjj(mh,0:mfc), ifpjj(mh,0:mfc),
     &         ivijj(mh,0:mfc,0:mfc))

*-----------------------------------------------------------------------
cHM 2017/11/29 add new variable to judge weather this function already have been called.
      second_call = .false.

*-----------------------------------------------------------------------
*     VARIABLES AND INITIALIZATION
*-----------------------------------------------------------------------

*        INV: =0, ;X IS GIVEN BY ARRAY
*             =1, ;X IS GIVEN BY NUMBERS ( default )

*        IX: CHECK OF X-AXIS NUMBER
*        IV: CHECK OF V-COLUMN

*        IDN : ID NUMBER OF COLUMN
*        IDD : ID NUMBER OF Y-COLUMN

*        IDDY : IDN -> IDD
*        IDYD : IDD -> IDN

*        IYYD(I,1[2]) = IDN + 1 -> NO ERROR +[-]
*        IYYD(I,1[2]) = I       -> WITH ERORR +[-]

*        IYDD(I) = 1 -> WITH ERROR -
*                = 2 -> WITH ERROR +
*                = 3 -> WITH ERROR +-

*        NO NUMBER  -> IDD = -1
*        NX   :     -> IDD =  0

*        IFDC(IDN) = 0 -> NUMBER
*                  > 1 -> FUNCTION

*        IVDC(IDN) = 0 -> IDN COLUMN HAS NOT BEEN EVALUATED
*                    1 -> IDN COLUMN HAS BEEN ALREADY EVALUATED

*        IFFY(IDN) = FITTING Y ID
*        IFFD(IDN) = 1 -> WITHOUT ERROR
*        IFFD(IDN) = 2 -> WITH ERROR

*        IDC(IDN) : IDENTIFICATION OF COLUMN

*        NNVV : NUMBER OF POINTS FOR FUNCTION DEFINE BY V=[1,2,3]

*-----------------------------------------------------------------------

*        IDD : ID NUMBER OF Y-COLUMN

*        IDDE : IDN -> IDD
*        IDED : IDD -> IDN ( for D+ )
*        IDEF : IDD -> IDN ( for D- )

*        NO NUMBER  -> IDD = -1
*        NO NUMBER  -> IDD = MH FOR X

*-----------------------------------------------------------------------

         do 50 i = -1, mh

            idyd(i) = -1
            ided(i) = -1
            idef(i) = -1
            iffy(i) = -1
            iffd(i) = -1

   50    continue

         do 51 i = 1, mh

            ifdc(i) = 0
            ivdc(i) = 0

   51    continue

*-----------------------------------------------------------------------

            nnvv = 0

            inv  = 1
            icol = 0

            ix = 0
            iv = 0

            ycob(1,1) = -r1max !FURUTA ???
            ycob(1,2) = 1.0
            ycob(1,3) = 1.0

*-----------------------------------------------------------------------
*        IC: POINTER OF COLUMN
*-----------------------------------------------------------------------

            ic = in - 1

*-----------------------------------------------------------------------

*     IDN: ID NUMBER

*        IFDC(IDN) = 0;  REAL NUMBER COLUMN
*                  = 1;  V-VALUE COLUMN
*                  = 2;  FUNCTION COLUMN
*                  = 3;  FUNCTION COLUMN WITH A12    (IAAC)
*                  = 4;  FUNCTION COLUMN FOR FITTING (IFFC)

*        IVDC(IDN) = 0;  NOT EVALUATED YET
*                  = 1;  ALREADY EVALUATED

*          IDY: ID OF Y FOR DEFAULT SYMBOL ( = 3 + n INITIAL )
*          IDE: ID OF ERROR
*          IDX: ID OF X-ERROR

*          IDXD: COLUMN OF X
*          IDVD: COLUMN OF V

*-----------------------------------------------------------------------

            idn = 0
            idy = 3
            ide = 0

            idxd = -1
            idvd = -1

*-----------------------------------------------------------------------
*     FAC:  FACTORS OF 1; /, 2; *, 3; -, 4; +, 5; **
*     JYCM(I): NUMBER OF CHRACTERS OF COMMENTS
*-----------------------------------------------------------------------

         do 10 i = 1, mh

            fac(i,1) =  1.0
            fac(i,2) =  1.0
            fac(i,3) =  0.0
            fac(i,4) =  0.0
            fac(i,5) = -r1max

            jycm(i)  =  0

   10    continue

*-----------------------------------------------------------------------
*     IDC(I) : IDENTIFICATION OF COLUMN
*-----------------------------------------------------------------------

*     IDC(I) = 1  :  X           - COLUMN
*     IDC(I) = 2  :  Y           - COLUMN
*     IDC(I) = 3  :
*     IDC(I) = 4  :  N,  NY      - COLUMN
*     IDC(I) = 5  :
*     IDC(I) = 6  :  D,  DY      - COLUMN
*     IDC(I) = 7  :  ND, NDY     - COLUMN
*     IDC(I) = 8  :  DX          - COLUMN
*     IDC(I) = 9  :  NDX         - COLUMN
*     IDC(I) = 10 :  V=[1,2,3]   - COLUMN

*-----------------------------------------------------------------------

  100    ic = ic + 1

*-----------------------------------------------------------------------

         if( ic .gt .ichrl ) goto 500

         if( lum(ic) .eq. ' ' .or .lum(ic) .eq. tub ) goto 100

         if( lum(ic) .ne. 'x' .and. lum(ic) .ne. 'y' .and.
     &       lum(ic) .ne. 'n' .and. lum(ic) .ne. 'v' .and.
     &       lum(ic) .ne. 'd') then

            m_err = 'H: Parameter Should be Started X, Y, V, N, or D.'
            ErrCha = ''
            ErrID = 'L:268/R:hone/F:a-hsect.f'
            l_err = ill(jsn)
            k_err = jsn

            goto 999

         end if

*-----------------------------------------------------------------------

*     IN FIRST LINE OF ONE SECTION, THE FIRST CHARACTER OF WORDS
*     SHOULD BE STARTED BY X, Y, N, OR D.

*           X: X-COLUMN
*           Y: Y-COLUMN
*           V: V-COLUMN
*           N: NY, ND, DUMY COLUMN, NOT PLOTTED
*           D: ERROR COLUMN

*=======================================================================

      isv = ( lum(ic) .eq. 'v' )
      if( isv ) then

*=======================================================================

                  if( lum(ic+1) .ne. '=' .or. lum(ic+2) .ne. '[' ) then

                     m_err = 'H: V-Column Should be V=[1,2,3].'
                     ErrCha = ''
                     ErrID = 'L:298/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

*-----------------------------------------------------------------------

               idn = idn + 1

                  if( idn .ne. 1 ) then

                     m_err = 'H: V-Column Should be Used'//
     &                       ' at First Column in H: line.'
                     ErrCha = ''
                     ErrID = 'L:315/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

*-----------------------------------------------------------------------

               iv = iv + 1

                  if( iv .gt. 1 ) then

                     m_err = 'H: Number of V-Column Should be One.'
                     ErrCha = ''
                     ErrID = 'L:331/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

*-----------------------------------------------------------------------
*           READ V-VALUE DESCRIPTION
*-----------------------------------------------------------------------

               ic = ic + 2

                  call func04(lum,ic,ichrl,ierr,
     &                        vvin,vvdd,nnvv)

                  if( ierr .ne. 0 ) then

                     m_err = 'H: Form of V-Column V=[1,2,3] is Wrong.'
                     ErrCha = ''
                     ErrID = 'L:352/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

      end if

*-----------------------------------------------------------------------

cHM 2017/11/29 judge to do initialise or not.
      if(.not. second_call) then
cHM 2017/11/29 if this is the first call, deallocate daxy and allocate daxy.
         if (allocated(daxy)) then
            deallocate(daxy)
         endif

         allocate(daxy(mh+2,nnvv+1))
         daxy(:,:) = 0.d0
      endif

*-----------------------------------------------------------------------
*           EVALUATE THE V-VALUES
*-----------------------------------------------------------------------

      if( isv ) then

cHM 2017/11/29 in the case of V-column, make second_call true not to initialise daxy variable.
         second_call = .true.

            do 60 i = 1, nnvv

               daxy(1,i) = vvin + dble(i-1) * vvdd

   60       continue

               ifdc(idn) = 1
               ivdc(idn) = 1

               idc(idn)  = 10
               idvd      = idn

               inv       = 0

*-----------------------------------------------------------------------

                  if( lum(ic+1) .ne. ' ' .and.
     &                lum(ic+1) .ne. tub ) then

                     m_err = 'H: Form of V-Column V=[1,2,3] is Wrong.'
                     ErrCha = ''
                     ErrID = 'L:405/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

*-----------------------------------------------------------------------

            goto 100

      end if

*-----------------------------------------------------------------------
*     END OF V-COLUMN DESCRIPTION
*-----------------------------------------------------------------------


*=======================================================================

      if( lum(ic) .eq. 'x' ) then

*=======================================================================

               ix = ix + 1

                  if( ix .gt. 1 ) then

                     m_err = 'H: Number of X-Column Should be One.'
                     ErrCha = ''
                     ErrID = 'L:436/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

*-----------------------------------------------------------------------

               idn = idn + 1

                  if( idn .gt. mh - 1 ) then

                     m_err = 'H: Number of Column Exceeds '//
     &                       'the Limit.'
                     ErrCha = ''
                     ErrID = 'L:453/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

*-----------------------------------------------------------------------

               idc(idn) = 1
               idxd     = idn
               ycol(idn,1) = -r1max
               ycol(idn,2) = 1.0
               ycol(idn,3) = 1.0
               ycob(idn,1) = -r1max
               ycob(idn,2) = 1.0
               ycob(idn,3) = 1.0

*-----------------------------------------------------------------------
*        FOR FUNCTION
*-----------------------------------------------------------------------

            if( lum(ic+1) .eq. '=' .and. lum(ic+2) .eq. '[' ) then

                  ic = ic + 2

               call func01(lum,ic,ichrl,ierr,']',
     &                     rk,kd,kp,ki,jd,mg,md,me,mp,id,iaac)

                     if( ierr .ne. 0 ) goto 910

               call func03(ierr,
     &                     kd,kp,ki,jd,mg,md,me,mp,id,
     &                     ivi,ivn,ifp,ipr)

                     if( ierr .ne. 0 ) goto 910


               call func07(idn,kd,ki,rk,ivi,ivn,ifp,ipr,id,
     &                     kdjj,kijj,rkjj,ivijj,ivnjj,ifpjj,
     &                     iprjj,idrjj)


                  ifdc(idn) = 2 + iaac

                  goto 920


  910             continue

                     m_err = 'H: Description of X=[  ] '//
     &                       'is Wrong.'
                     ErrCha = ''
                     ErrID = 'L:507/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

  920          continue

            else

                  if( inv .eq. 0 ) then

                     m_err = 'H: X Should be a Function of V=[ ] '//
     &                       'i.e. X=[  ] .'
                     ErrCha = ''
                     ErrID = 'L:522/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

               ivdc(idn) = 1

            end if

*-----------------------------------------------------------------------
*           FACTOR OF X-VALUE
*-----------------------------------------------------------------------

               ic = ic + 1

                  if( lum(ic) .eq. ' ' .or. lum(ic) .eq. tub .or.
     &                ic .gt. ichrl ) goto 100

               ici = ic

*-----------------------------------------------------------------------

  110          ic = ic + 1

                  if( ic .le .ichrl .and.
     &                lum(ic) .ne. ' ' .and.
     &                lum(ic) .ne. tub ) goto 110

                  if( ic - 1 .eq. ici ) then

                     m_err = 'H: X-Factor Description is Wrong,'
                     ErrCha = ''
                     ErrID = 'L:557/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

               call a_fact(lum,ichrl,ici,ic-1,dfac,ierr)

                  if( ierr .ne. 0 ) then

                     m_err = 'H: X-Factor Description is Wrong.'
                     ErrCha = ''
                     ErrID = 'L:571/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

                  do 112 i = 1, 5

                     fac(idn,i) = dfac(i)

  112             continue

*-----------------------------------------------------------------------

            goto 100

      end if

*-----------------------------------------------------------------------
*     END OF X-COLUMN
*-----------------------------------------------------------------------

*=======================================================================
      if( lum(ic) .eq. 'y' .or.
     &  ( lum(ic) .eq. 'n' .and. lum(ic+1) .ne. 'd' )) then

*=======================================================================
*        Y, N, NY
*-----------------------------------------------------------------------

               idn = idn + 1

                if( idn .gt. mh - 1 ) then

                   write(ErrCha,*)
     &            'Warning: ignore more than 20 column data in ANGEL'
                   ErrID = 'L:609/R:hone/F:a-hsect.f'
                   call ErrWrite(ErrID,ErrCha)
                   goto 500  ! change error to warning T.Sato 2023/04/06

                end if

*-----------------------------------------------------------------------

               if( lum(ic) .eq. 'y' ) then

                  idy = idy + 1

                  idyb = idy

               else

                  idyb = 1

               end if

               iylc(idn) = 0

               iyss(idn) = 0

               syss(idn) = 4.0

               iylt(idn) = 4
               iyhh(idn) = 0
               iybb(idn) = 0

               ycol(idn,1) = -r1max
               ycol(idn,2) = 1.0
               ycol(idn,3) = 1.0

               ycob(idn,1) = -r1max
               ycob(idn,2) = 1.0
               ycob(idn,3) = 1.0

               iysb(idn) = idyb

               if( iysb(idn) .gt. 16 ) then

                  idy = 1

                  iysb(idn) = 1

               end if

               iffc = 0
               icol = 0

*-----------------------------------------------------------------------

               if( lum(ic) .eq. 'y' ) then

                     idc(idn) = 2

               else if( lum(ic) .eq. 'n' ) then

                  if( lum(ic+1) .eq. 'y' ) then

                     idc(idn) = 4
                     ic = ic + 1

                  else if( lum(ic+1) .eq. 'x' ) then

                      m_err = 'H: NX is Error. Please Use NY'
                      ErrCha = ''
                      ErrID = 'L:677/R:hone/F:a-hsect.f'
                      l_err = ill(jsn)
                      k_err = jsn

                      goto 999

                  else

                   if( lum(ic+1) .ne. ' ' .and.
     &                 lum(ic+1) .ne. '/' .and. lum(ic+1) .ne. '*' .and.
     &                 lum(ic+1) .ne. '+' .and. lum(ic+1) .ne. '-' .and.
     &                 lum(ic+1) .ne. '(' .and. lum(ic+1) .ne. ',' .and.
     &                 lum(ic+1) .ne. '0' .and. lum(ic+1) .ne. '1' .and.
     &                 lum(ic+1) .ne. '2' .and. lum(ic+1) .ne. '3' .and.
     &                 lum(ic+1) .ne. '4' .and. lum(ic+1) .ne. '5' .and.
     &                 lum(ic+1) .ne. '6' .and. lum(ic+1) .ne. '7' .and.
     &                 lum(ic+1) .ne. '8' .and. lum(ic+1) .ne. '9') then

                        m_err = 'H: X, Y, N Description is Wrong.'
                        ErrCha = ''
                        ErrID = 'L:697/R:hone/F:a-hsect.f'
                        l_err = ill(jsn)
                        k_err = jsn

                        goto 999

                   end if

                     idc(idn) = 4

                  end if

               end if

*-----------------------------------------------------------------------

*     IDD : ID NUMBER OF Y-COLUMN

*       IDDY : IDN -> IDD
*       IDYD : IDD -> IDN

*       NO NUMBER  -> IDD = -1

*-----------------------------------------------------------------------

            ic = ic + 1

                  call idd2(lum,ic,ichrl,idd,ierr)

                  if( ierr .ne. 0 ) then

                      m_err = 'H: Description of Y Specific'//
     &                        ' Number is Wrong.'
                      ErrCha = ''
                      ErrID = 'L:731/R:hone/F:a-hsect.f'
                      l_err = ill(jsn)
                      k_err = jsn

                      goto 999

                  end if


                  iddy(idn) = idd
                  idyd(idd) = idn

*-----------------------------------------------------------------------
*        FOR FUNCTION
*-----------------------------------------------------------------------

            if( ( lum(ic) .eq. '=' .and. lum(ic+1) .eq. '[' ) .or.
     &          ( lum(ic) .eq. '=' .and. lum(ic+1) .eq. 'f' .and.
     &            lum(ic+2) .eq. '{' ) .or.
     &          ( lum(ic) .eq. '=' .and. lum(ic+1) .eq. 'f' .and.
     &            lum(ic+2) .eq. 'd' .and. lum(ic+3) .eq. '{' ) ) then

*-----------------------------------------------------------------------

               if( lum(ic+1) .eq. 'f' ) then

*-----------------------------------------------------------------------

                  if( lum(ic+2) .eq. '{' ) then

                     iffd(idn) = 1

                     ic = ic + 2

                  else if( lum(ic+3) .eq. '{' ) then

                     iffd(idn) = 2

                     ic = ic + 3

                  end if


  400             ic = ic + 1

                     if( lum(ic) .eq. ' ' .or .lum(ic) .eq. tub )
     &               goto 400

                     if( lum(ic) .ne. 'y' ) goto 922

                  ic = ic + 1

                  call idd2(lum,ic,ichrl,idd,ierr)

                     if( ierr .ne. 0 ) goto 922

                     iffy(idn) = idd

                     ic = ic - 1

  401             ic = ic + 1

                     if( lum(ic) .eq. ' ' .or .lum(ic) .eq. tub )
     &               goto 401

                     if( lum(ic) .ne. '}' .and.
     &                   lum(ic+1) .ne. '[' ) goto 922

                     iffc = 1


                  goto 923

  922             continue

                        m_err = 'H: Description of Y3=F{Y2}[ ] is '//
     &                          ' Wrong.'
                        ErrCha = ''
                        ErrID = 'L:809/R:hone/F:a-hsect.f'
                         l_err = ill(jsn)
                         k_err = jsn

                         goto 999

  923             continue

               end if

*-----------------------------------------------------------------------

                  ic = ic + 1

               call func01(lum,ic,ichrl,ierr,']',
     &                     rk,kd,kp,ki,jd,mg,md,me,mp,id,iaac)

                     if( ierr .ne. 0 ) goto 911

                     if( iffc .eq. 1 .and. iaac .eq. 0 ) goto 911

               call func03(ierr,
     &                     kd,kp,ki,jd,mg,md,me,mp,id,
     &                     ivi,ivn,ifp,ipr)

                     if( ierr .ne. 0 ) goto 911


               call func07(idn,kd,ki,rk,ivi,ivn,ifp,ipr,id,
     &                     kdjj,kijj,rkjj,ivijj,ivnjj,ifpjj,
     &                     iprjj,idrjj)


                  ifdc(idn) = 2 + iaac + iffc

                  goto 921


  911             continue

                     m_err = 'H: Description of Y=[  ] '//
     &                       'is Wrong.'
                     ErrCha = ''
                     ErrID = 'L:852/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999


  921          ic = ic + 1

*-----------------------------------------------------------------------

            else

                  if( inv .eq. 0 ) then

                     m_err = 'H: Y Should be a Function of V=[ ] '//
     &                       'i.e. Y=[  ] .'
                     ErrCha = ''
                     ErrID = 'L:870/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

               ivdc(idn) = 1

            end if

*-----------------------------------------------------------------------

         if( ic .gt. ichrl .or.
     &       lum(ic) .eq. ' '. or. lum(ic) .eq. tub ) goto 100


  200  continue

*-----------------------------------------------------------------------
*        FACTOR OF Y-VALUE
*-----------------------------------------------------------------------

            if( lum(ic) .eq. '/' .or.
     &          lum(ic) .eq. '*' .or.
     &          lum(ic) .eq. '-' .or.
     &          lum(ic) .eq. '+' ) then

               ici = ic

  120          ic = ic + 1

               if( ic .le. ichrl .and.
     &             lum(ic) .ne. ' ' .and .lum(ic) .ne. tub .and.
     &             lum(ic) .ne. ',' .and. lum(ic) .ne. '(' ) goto 120

               if( ic .le. ichrl - 2 .and. lum(ic) .eq. '(' .and.
     &           ( lum(ic-1) .eq. '*' .or.
     &             lum(ic-1) .eq. '/' .or.
     &             lum(ic-1) .eq. '+' .or.
     &             lum(ic-1) .eq. '-' ) ) goto 120

               if( ic-1 .eq. ici ) then

                  m_err = 'H: Y-Factor Description is Wrong,'
                  ErrCha = ''
                  ErrID = 'L:917/R:hone/F:a-hsect.f'
                  l_err = ill(jsn)
                  k_err = jsn

                  goto 999

               end if

                  call a_fact(lum,ichrl,ici,ic-1,dfac,ierr)

                     if( ierr .ne. 0 ) then

                        m_err = 'H: Y-Factor Description is Wrong.'
                        ErrCha = ''
                        ErrID = 'L:931/R:hone/F:a-hsect.f'
                        l_err = ill(jsn)
                        k_err = jsn

                        goto 999

                     end if

                     do 122 i =1, 5

                          fac(idn,i) = dfac(i)

  122                continue

            end if

*-----------------------------------------------------------------------
*        JYCM, YJCM : COMMENT OF Y
*-----------------------------------------------------------------------

            if(lum(ic).eq.'(') then

               ici = ic + 1

  130          ic = ic + 1

                  if( ic .le. ichrl ) then

                     if( lum(ic) .eq. ')' .and.
     &                   lum(ic-1) .ne. yen ) goto 131

                     goto 130

                  end if

  131             if( ic .gt. ichrl .or. ici .eq. ic ) then

                     m_err = 'H: Y-Comment Description is Wrong.'
                     ErrCha = ''
                     ErrID = 'L:970/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

                  ijycm = 0

                  do 132 i = ici, ic-1

                     ijycm = ijycm + 1

                     jdum(ijycm) = dum(i)

  132             continue


                  call jpncode(jdum,ijycm,ifon)
                  call jpnprep(jdum,ijycm,ifon)


                  jycm(idn) = ijycm

               if( idc(idn) .eq. 2 ) then

                  write(jhy) ( jdum(i), i = 1, ijycm )

               end if

               ic = ic + 1

            end if

*-----------------------------------------------------------------------

         if( lum(ic) .eq. '/' .or. lum(ic) .eq. '*' .or.
     &       lum(ic) .eq. '-' .or. lum(ic) .eq. '+' .or.
     &       lum(ic) .eq. '(' ) goto 200

         if( lum(ic) .eq. ' ' .or. lum(ic) .eq. tub .or.
     &       ic .gt. ichrl ) goto 100

         if( lum(ic) .ne. ',' ) then

            m_err = 'H: Y Description is Wrong.'
            ErrCha = ''
            ErrID = 'L:1018/R:hone/F:a-hsect.f'
            l_err = ill(jsn)
            k_err = jsn

            goto 999

         end if

*-----------------------------------------------------------------------

  140    ic = ic + 1

         if( ic .gt. ichrl .or.
     &       lum(ic) .eq. ' ' .or. lum(ic) .eq. tub ) goto 100

*-----------------------------------------------------------------------
*     IYLC : LINE STYLE

*         'N' IYLC = -1 :  NO LINE
*         'L' IYLC =  0 :  SOLID
*         'M' IYLC =  1 :  DOT-DASHED
*         'D' IYLC =  2 :  SHORT-DASHED
*         'U' IYLC =  3 :  LONG-DASHED
*         'P' IYLC =  4 :  DOTTED
*         'Q' IYLC =  5 :  DOT-DOT-DASHED
*         'V' IYLC =  6 :  DOT-DOT-DOT-DASHED
*         'I' IYLC = 11 :  INTERIA

*-----------------------------------------------------------------------

            if( lum(ic) .eq. 'n' ) then

               iylc(idn) = -1
               goto 140

            end if

            if( lum(ic) .eq. 'l' ) then

               iylc(idn) = 0
               goto 140

            end if

            if( lum(ic) .eq. 'm' ) then

               iylc(idn) = 1
               goto 140

            end if

            if( lum(ic) .eq. 'd' ) then

               iylc(idn) = 2
               goto 140

            end if

            if( lum(ic) .eq. 'u' ) then

               iylc(idn) = 3
               goto 140

            end if

            if( lum(ic) .eq. 'p' ) then

               iylc(idn) = 4
               goto 140

            end if

            if( lum(ic) .eq. 'q' ) then

               iylc(idn) = 5
               goto 140

            end if

            if( lum(ic) .eq. 'v' ) then

               iylc(idn) = 6
               goto 140

            end if

            if( lum(ic) .eq. 'i' .and. lum(ic+1) .ne. 'i' ) then

               iylc(idn) = 11
               goto 140

            end if

            if( lum(ic) .eq. 'i' .and. lum(ic+1) .eq. 'i' ) then

               iylc(idn) = 12
               ic = ic + 1

               goto 140

            end if

*-----------------------------------------------------------------------
*     YCOL : COLOR OF REGIONS, LINES AND SYMBOLS
*     YCOB : COLOR OF INNER OF SYMBOL FOR 3, 5, 7, 9, 11
*-----------------------------------------------------------------------

*            'R' YCOL = 1.0 :  RED
*            'Y' YCOL = 1.5 :  YELLOW
*            'G' YCOL = 2.0 :  GREEN
*            'C' YCOL = 2.5 :  CYAN
*            'B' YCOL = 3.0 :  BLUE

*-----------------------------------------------------------------------
*           COLOR BY NUMERIC
*-----------------------------------------------------------------------

            if( lum(ic) .eq. 'c' .and. lum(ic+1) .eq. '[' ) then

                  call dcols(0,lum,ic,coln,ierr,ichrl,'[',']')

                     if( ierr .ne. 0 ) then

                        m_err = 'H: Number of C[ ] is Wrong'
                        ErrCha = ''
                        ErrID = 'L:1143/R:hone/F:a-hsect.f'
                        l_err = ill(jsn)
                        k_err = jsn

                        goto 999

                     end if

                  ycol(idn,1) = coln(1)
                  ycol(idn,2) = coln(2)
                  ycol(idn,3) = coln(3)

                  ic = ic - 1
                  icol = icol + 1

                  goto 140

            end if

*-----------------------------------------------------------------------
*           COLOR BY SYMBOL
*-----------------------------------------------------------------------

            if( lum(ic) .eq. 'r' ) then
               if( icol .eq. 0 ) then
                  if( ycol(idn,1) .gt. 0.0 ) then
                    ycol(idn,1) = ycol(idn,1) + 0.16667
                  else
                    ycol(idn,1) = 1.0
                  end if
               end if
                  goto 140
            end if

            if( lum(ic) .eq. 'y' ) then
               if( icol .eq. 0 ) then
                  if( ycol(idn,1) .gt. 0.0 ) then
                    ycol(idn,1) = ycol(idn,1) + 0.16667
                  else
                    ycol(idn,1) = 1.5
                  end if
               end if
                  goto 140
            end if

            if( lum(ic) .eq. 'g' ) then
               if( icol .eq. 0 ) then
                  if( ycol(idn,1) .gt. 0.0 ) then
                    ycol(idn,1) = ycol(idn,1) + 0.16667
                  else
                    ycol(idn,1) = 2.0
                  end if
               end if
                  goto 140
            end if

            if( lum(ic) .eq. 'c' ) then
               if( icol .eq. 0 ) then
                  if( ycol(idn,1) .gt. 0.0 ) then
                    ycol(idn,1) = ycol(idn,1) + 0.16667
                  else
                    ycol(idn,1) = 2.5
                  end if
               end if
                  goto 140
            end if

            if( lum(ic) .eq. 'b' ) then
               if( icol .eq. 0 ) then
                  if( ycol(idn,1) .gt. 0.0 ) then
                    ycol(idn,1) = ycol(idn,1) + 0.16667
                  else
                    ycol(idn,1) = 3.0
                  end if
               end if
                  goto 140
            end if

*-----------------------------------------------------------------------
*     YCOL : GRAY SCALE OF REGIONS, LINES AND SYMBOLS
*-----------------------------------------------------------------------

*            'W' YCOL = -1.0 :  WHITE
*            'O' YCOL = -1.1 :  GRAY TO WHITE
*            'K' YCOL = -1.3 :  GRAY TO WHITE
*            'J' YCOL = -1.5 :  GRAY TO WHITE
*            'F' YCOL = -1.7 :  GRAY TO WHITE
*            'E' YCOL = -2.0 :  BLACK (DEFAULT)

*-----------------------------------------------------------------------

            if( lum(ic) .eq. 'w' ) then
               if( icol .eq. 0 ) then
                  if( ycol(idn,1) .gt. -r0max .and.
     &                ycol(idn,1) .le. 0.0 ) then
                     ycol(idn,1) = ycol(idn,1) - 0.05
                  else
                     ycol(idn,1) = -1.0
                  end if
                     ycol(idn,1) = max( -2.0d0, ycol(idn,1) )
               end if
                  goto 140
            end if

            if( lum(ic) .eq. 'e' ) then
               if( icol .eq. 0 ) then
                  if( ycol(idn,1) .gt. -r0max .and.
     &                ycol(idn,1) .le. 0.0 ) then
                     ycol(idn,1) = ycol(idn,1) - 0.05
                  else
                     ycol(idn,1) = -2.0
                  end if
                     ycol(idn,1) = max( -2.0d0, ycol(idn,1) )
               end if
                  goto 140
            end if

            if( lum(ic) .eq. 'f' ) then
               if( icol .eq. 0 ) then
                  if( ycol(idn,1) .gt. -r0max .and.
     &                ycol(idn,1) .le. 0.0 ) then
                     ycol(idn,1) = ycol(idn,1) - 0.05
                  else
                     ycol(idn,1) = -1.8
                  end if
                     ycol(idn,1) = max( -2.0d0, ycol(idn,1) )
               end if
                  goto 140
            end if

            if( lum(ic) .eq. 'j' ) then
               if( icol .eq. 0 ) then
                  if( ycol(idn,1) .gt. -r0max .and.
     &                ycol(idn,1) .le. 0.0 ) then
                     ycol(idn,1) = ycol(idn,1) - 0.05
                  else
                     ycol(idn,1) = -1.6
                  end if
                     ycol(idn,1) = max( -2.0d0, ycol(idn,1) )
               end if
                  goto 140
            end if

            if( lum(ic) .eq. 'k' ) then
               if( icol .eq. 0 ) then
                  if( ycol(idn,1) .gt. -r0max .and.
     &                ycol(idn,1) .le. 0.0 ) then
                     ycol(idn,1) = ycol(idn,1) - 0.05
                  else
                     ycol(idn,1) = -1.4
                  end if
                     ycol(idn,1) = max( -2.0d0, ycol(idn,1) )
               end if
                  goto 140
            end if

            if( lum(ic) .eq. 'o' ) then
               if( icol .eq. 0 ) then
                  if( ycol(idn,1) .gt. -r0max .and.
     &                ycol(idn,1) .le. 0.0 ) then
                     ycol(idn,1) = ycol(idn,1) - 0.05
                  else
                     ycol(idn,1) = -1.2
                  end if
                     ycol(idn,1) = max( -2.0d0, ycol(idn,1) )
               end if
                  goto 140
            end if

*-----------------------------------------------------------------------
*     IYSS : 'S' SPLINE [ ]
*     IYLT : 'T' THICK LINE, (4), T(7), TT(10), TTT(13)
*     IYLT : 'Z' THIN LINE,  (4), Z(3), ZZ(2), ZZZ(1)
*     IYBB : 'A' BIGER SIMBOL A(1.5), AA(2.0), AAA(2.5)
*     IYBB : 'X' SMALL SIMBOL X(0.75), XX(0.5), XXX(0.25)
*     IYHH : 'H' HISTGRAM
*-----------------------------------------------------------------------

            if( lum(ic) .eq. 's' ) then

                  iyss(idn) = 1

               if( lum(ic+1) .eq. '[' ) then

                  ic = ic + 1

                  call pnum(lum,ic,ichrl,syrn,ierr)

                     if( ierr .ne. 0 ) then

                        m_err = 'H: Number of S[ ] is Wrong'
                        ErrCha = ''
                        ErrID = 'L:1335/R:hone/F:a-hsect.f'
                        l_err = ill(jsn)
                        k_err = jsn

                        goto 999

                     end if

                  syss(idn) = syrn

                     isys = nint( syss(idn) )

                     if( isys .eq. 0 ) then

                        m_err = 'H: Number of S[ ] Should not be Zero'
                        ErrCha = ''
                        ErrID = 'L:1351/R:hone/F:a-hsect.f'
                        l_err = ill(jsn)
                        k_err = jsn

                        goto 999

                     end if

               end if

               goto 140

            end if

            if( lum(ic) .eq. 't' ) then

               iylt(idn) = iylt(idn) + 3

               goto 140

            end if

            if( lum(ic) .eq. 'z' ) then

               iylt(idn) = iylt(idn) - 1

               if( iylt(idn) .eq. 0 ) iylt(idn) = 1

               goto 140

            end if

            if( lum(ic) .eq. 'a' ) then

               iybb(idn) = iybb(idn)+1

               goto 140

            end if

            if( lum(ic) .eq. 'x' ) then

               iybb(idn) = iybb(idn)-1

               goto 140

            end if

            if( lum(ic) .eq. 'h' ) then

               iyhh(idn) = iyhh(idn) + 1

               goto 140

            end if

*-----------------------------------------------------------------------
*     MARKER TYPE

*        '0'   : NO MARKER
*        '1'   : DOT
*        '2'   : PLUS
*        '3'   : OPEN  CIRCLE
*        '4'   : SOLID CIRCLE
*        '5'   : OPEN  SQUARE
*        '6'   : SOLID SQUARE
*        '7'   : OPEN  TRIANGLE (UP)
*        '8'   : SOLID TRIANGLE (UP)
*        '9'   : OPEN  TRIANGLE (DOWN)
*        '10'  : SOLID TRIANGLE (DOWN)
*        '11'  : OPEN  DIAMOND
*        '12'  : SOLID DIAMOND
*-----------------------------------------------------------------------

            if( deqn1(lum(ic)) ) then


               call idd2(lum,ic,ichrl,idd,ierr)


                  if( ierr .ne. 0 ) then

                     m_err = 'H: Y, Marker-Type Description is Wrong.'
                     ErrCha = ''
                     ErrID = 'L:1435/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if


                  if( idd .lt. 0 .or. idd .gt. 34 .or.
     &              ( idd .ge. 17 .and. idd .le. 22 ) ) then

                     m_err = 'H: Y, This Marker-Type is Not Supported.'
                     ErrCha = ''
                     ErrID = 'L:1449/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

                  if( idd .ge. 23 ) then

                     if( idd/2*2 .ne. idd ) then

                        idd = ( idd - 23 ) / 2 + 17

                     else

                     m_err = 'H: Y, This Marker-Type is Not Supported.'
                     ErrCha = ''
                     ErrID = 'L:1467/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                     end if

                  end if


                     iysb(idn) = idd


               if( lum(ic) .eq. '[' ) then

                     ic = ic - 1

                  call dcols(1,lum,ic,coln,ierr,icolm,'[',']')

                  if( ierr .ne. 0 ) then

                     m_err = 'H: Y, Color of Inside the Symbol is'//
     &                       ' Wrong. 23[r] for example.'
                     ErrCha = ''
                     ErrID = 'L:1492/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

                     ycob(idn,1) = coln(1)
                     ycob(idn,2) = coln(2)
                     ycob(idn,3) = coln(3)

               end if

                  ic = ic - 1

               goto 140


            end if

*-----------------------------------------------------------------------

               m_err = 'H: Y, Marker-Type, Line-Style, '//
     &                 'Description is Wrong.'
               ErrCha = ''
               ErrID = 'L:1518/R:hone/F:a-hsect.f'
               l_err = ill(jsn)
               k_err = jsn

               goto 999

*-----------------------------------------------------------------------
*     END OF Y SECTION
*-----------------------------------------------------------------------

      end if

*=======================================================================

      if( lum(ic) .eq. 'd' .or.
     &  ( lum(ic) .eq. 'n' .and. lum(ic+1) .eq. 'd' ) ) then

*=======================================================================

            idn = idn + 1

                if( idn .gt. mh - 1 ) then

                   write(ErrCha,*)
     &            'Warning: ignore more than 20 column data in ANGEL'
                   ErrID = 'L:1543/R:hone/F:a-hsect.f'
                   call ErrWrite(ErrID,ErrCha)
                   goto 500  ! change error to warning T.Sato 2023/04/06

                end if

*-----------------------------------------------------------------------

               ide = ide + 1
               ides(idn) = 0

            if( lum(ic) .eq. 'd' ) then

               if( lum(ic+1) .eq. 'x' ) then

                  idc(idn) = 8
                  ic = ic + 1

               else
     &         if( lum(ic+1) .eq. 'y' ) then

                  idc(idn) = 6
                  ic = ic + 1

               else

                  idc(idn) = 6

               end if

            else

               if( lum(ic+2) .eq. 'x' ) then

                  idc(idn) = 9
                  ic = ic + 2

               else
     &         if( lum(ic+2) .eq. 'y' ) then

                  idc(idn) = 7
                  ic = ic + 2

               else

                  idc(idn) = 7
                  ic = ic + 1

               end if

            end if

*-----------------------------------------------------------------------
*           IDD : ID NUMBER OF Y-COLUMN

*             IDDE : IDN -> IDD
*             IDED : IDD -> IDN ( for D+ )
*             IDEF : IDD -> IDN ( for D- )

*             NO NUMBER  -> IDD = -1
*             NO NUMBER  -> IDD = MH FOR X
*-----------------------------------------------------------------------

         ic = ic + 1

            call idd2(lum,ic,ichrl,idd,ierr)

                  if( ierr .ne. 0 ) then

                    m_err = 'H: Description of D, ND Specific '//
     &                      'Number is Wrong.'
                    ErrCha = ''
                    ErrID = 'L:1615/R:hone/F:a-hsect.f'
                    l_err = ill(jsn)
                    k_err = jsn

                    goto 999

                  end if

                  if( idd .ne. -1 .and. lum(ic-1) .eq. 'x' ) then

                      m_err = 'H: DX Specific'//
     &                        ' Number is not necessary.'
                      ErrCha = ''
                      ErrID = 'L:1628/R:hone/F:a-hsect.f'
                      l_err = ill(jsn)
                      k_err = jsn

                      goto 999

                  end if

                  if( idd .eq. -1 .and. lum(ic-1) .eq. 'x' ) idd = mh


               idde(idn) = idd

*-----------------------------------------------------------------------
*        FOR FUNCTION
*-----------------------------------------------------------------------

            if( lum(ic) .eq. '=' .and. lum(ic+1) .eq. '[' ) then

                  ic = ic + 1

               call func01(lum,ic,ichrl,ierr,']',
     &                     rk,kd,kp,ki,jd,mg,md,me,mp,id,iaac)

                     if( ierr .ne. 0 ) goto 931

               call func03(ierr,
     &                     kd,kp,ki,jd,mg,md,me,mp,id,
     &                     ivi,ivn,ifp,ipr)

                     if( ierr .ne. 0 ) goto 931


               call func07(idn,kd,ki,rk,ivi,ivn,ifp,ipr,id,
     &                     kdjj,kijj,rkjj,ivijj,ivnjj,ifpjj,
     &                     iprjj,idrjj)


                  ifdc(idn) = 2 + iaac

                  goto 941


  931             continue

                     m_err = 'H: Description of D=[  ] '//
     &                       'is Wrong.'
                     ErrCha = ''
                     ErrID = 'L:1676/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999


  941          ic = ic + 1

            else

                  if( inv .eq. 0 ) then

                     m_err = 'H: D Should be a Function of V=[ ] '//
     &                       'i.e. D=[  ] .'
                     ErrCha = ''
                     ErrID = 'L:1692/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

               ivdc(idn) = 1

            end if

*-----------------------------------------------------------------------

         if( ic .gt. ichrl .or.
     &       lum(ic) .eq. ' ' .or. lum(ic) .eq. tub ) goto 100


  300    continue

*-----------------------------------------------------------------------

               if( lum(ic) .eq. '+' ) then

                 ides(idn) = ides(idn) + 1
                 ided(idd) = idn

                 ic = ic + 1

               end if

               if( lum(ic) .eq. '-' ) then

                 ides(idn) = ides(idn) - 1
                 idef(idd) = idn

                 ic = ic + 1

               end if

*-----------------------------------------------------------------------
*     FACTOR OF ERROR COLUMN
*-----------------------------------------------------------------------

          if( lum(ic) .eq. '/' .or. lum(ic) .eq. '*' ) then

             ici = ic

  320        ic = ic + 1

             if( ic .le. ichrl .and.
     &         ( ( lum(ic) .ne. '+' .or. lum(ic-1) .eq. 'e' .or.
     &             lum(ic-1) .eq. '*' ) .and.
     &           ( lum(ic) .ne. '-' .or. lum(ic-1) .eq. 'e' .or.
     &             lum(ic-1) .eq. '*' ) .and.
     &           ( lum(ic) .ne. ' ' .and. lum(ic) .ne. tub ) ) )
     &       goto 320

                 if( ic-1 .eq. ici ) goto 999

               call a_fact(lum,ichrl,ici,ic-1,dfac,ierr)

                  if(ierr.ne.0) then

                     m_err = 'H: Factor of Error Column is Wrong.'
                     ErrCha = ''
                     ErrID = 'L:1758/R:hone/F:a-hsect.f'
                     l_err = ill(jsn)
                     k_err = jsn

                     goto 999

                  end if

                  do 322 i = 1, 5

                     fac(idn,i) = dfac(i)

  322             continue

         end if

*-----------------------------------------------------------------------

         if( lum(ic) .eq. '+' .or. lum(ic) .eq. '-' .or.
     &       lum(ic) .eq. '/' .or. lum(ic) .eq. '*' ) goto 300

*-----------------------------------------------------------------------

         if( ided(idd) .eq. -1 .and. idef(idd) .eq. -1 ) then

                  ided(idd) = idn
                  idef(idd) = idn

         end if

*-----------------------------------------------------------------------

         if( ic .gt. ichrl. or.
     &       lum(ic) .eq.' ' .or. lum(ic) .eq. tub ) goto 100


            m_err = 'H: Parameter Description of Error (D, ND)is Wrong.'
            ErrCha = ''
            ErrID = 'L:1796/R:hone/F:a-hsect.f'
            l_err = ill(jsn)
            k_err = jsn

            goto 999

*-----------------------------------------------------------------------
*     END OF D SECTION
*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

         goto 100

*-----------------------------------------------------------------------
*        END OF H: LINE
*-----------------------------------------------------------------------

  500    continue

*-----------------------------------------------------------------------
*     CHECK X-COLUMN
*-----------------------------------------------------------------------

            if( ix .eq. 0 ) then

               m_err = 'H: There is No X-Column in This Section.'
               ErrCha = ''
               ErrID = 'L:1826/R:hone/F:a-hsect.f'
               l_err = ill(jsn)
               k_err = jsn

               goto 999

            end if

*-----------------------------------------------------------------------
*     SUMMARY OF ERROR BAR
*-----------------------------------------------------------------------

                  jdy = 0
                  jdx = 0

               do 505 i = 1, idn

                  iyyd(i,1) = idn + 1
                  iyyd(i,2) = idn + 1

                  if( idc(i) .eq. 6 .or. idc(i) .eq. 7 ) jdy = jdy + 1
                  if( idc(i) .eq. 8 .or. idc(i) .eq. 9 ) jdx = jdx + 1

  505          continue

*-----------------------------------------------------------------------

            if( jdy + jdx .gt. 0 ) then

                  ly = 0
                  lx = 0

               do 720 i = 1, idn

                  if( idc(i) .eq. 1 ) lx = i

                  if( idc(i) .eq. 2 .or. idc(i) .eq. 4 ) ly = i

                  if( idc(i) .eq. 6 .or. idc(i) .eq. 8 ) then

                     if( idc(i) .eq. 8 ) ly = lx

                     if( idde(i) .eq. -1 .or. idde(i) .eq. mh ) then

                        if( ly .eq. 0 ) then

                           m_err = 'H: Error Column Should be '//
     &                             'Just After the Data Column, '//
     &                             'or, Use ID Number'
                           ErrCha = ''
                           ErrID = 'L:1876/R:hone/F:a-hsect.f'
                           l_err = ill(jsn)
                           k_err = jsn

                           goto 999

                        end if

                           iddd = ly

                     else

                        if( idyd(idde(i)) .eq. -1 ) then

                           m_err = 'H: ID Mismatch of Error Column, '//
     &                             'There is no Y-Data Column.'
                           ErrCha = ''
                           ErrID = 'L:1893/R:hone/F:a-hsect.f'
                           l_err = ill(jsn)
                           k_err = jsn

                           goto 999

                        end if

                        if( idc(idyd(idde(i))) .eq. 2 .or.
     &                      idc(idyd(idde(i))) .eq. 4 ) then

                            iddd = idyd(idde(i))

                        else

                           m_err = 'H: ID Mismatch of Error Column '//
     &                             'and the Data Column.'
                           ErrCha = ''
                           ErrID = 'L:1911/R:hone/F:a-hsect.f'
                           l_err = ill(jsn)
                           k_err = jsn

                           goto 999

                        end if

                     end if


                        if( ides(i) .eq. 0 )  then

                           iyyd(iddd,1) = i
                           iyyd(iddd,2) = i

                        else if( ides(i) .eq.  1 ) then

                           iyyd(iddd,2) = i

                        else if( ides(i) .eq. -1 ) then

                           iyyd(iddd,1) = i

                        else

                           m_err = 'H: Error Description is Wrong.'
                           ErrCha = ''
                           ErrID = 'L:1939/R:hone/F:a-hsect.f'
                           l_err = ill(jsn)
                           k_err = jsn

                           goto 999

                        end if

                  end if

  720          continue

            end if

*-----------------------------------------------------------------------

               do 721 i = 1, idn

                     iydd(i) = 0

                     if( iyyd(i,1) .ne. idn+1 ) iydd(i) = 1
                     if( iyyd(i,2) .ne. idn+1 ) iydd(i) = 2

                     if( iyyd(i,1) .ne. idn+1 .and.
     &                   iyyd(i,2) .ne. idn+1 ) iydd(i) = 3

  721          continue


*-----------------------------------------------------------------------
*     SUMMARY AND READ DATA
*-----------------------------------------------------------------------

            call hone2(dum,lum,jsi,jsn,ill,ilf,dsin,idsi,daxy,
     &                 ymax,ymin,ymin2,
     &                 xmax,xmin,xmin2,noned,noner,ncom,ierr,jof,jif,
     &                 iycm,rycm,iddy,ided,idef,iffy,iffd,iydd,
     &                 idn,idc,fac,idde,idyd,iyyd,ides,iyhh,iylc,
     &                 iylt,iyss,syss,iysb,iybb,jycm,ycol,ycob,
     &                 inv,ifdc,ivdc,nnvv,idvd,idxd,
     &                 iprjj,idrjj,kdjj,kijj,rkjj,ivijj,ivnjj,ifpjj,
     &                 nevk,nevn,ivi,rk,kd,kp,ki,jd,md,mg,me,mp,
     &                 ivn,ifp)

            deallocate(daxy)
            deallocate(ivi, nevk, nevn, rkjj,
     &                 kdjj, kijj, ivnjj, ifpjj, ivijj)
      return

*-----------------------------------------------------------------------

  999 ierr = 1

      if (allocated(daxy)) deallocate(daxy)
      deallocate(ivi, nevk, nevn, rkjj,
     &           kdjj, kijj, ivnjj, ifpjj, ivijj)
      return
      end


************************************************************************
*                                                                      *
      subroutine hone2(dum,lum,jsi,jsn,ill,ilf,dsin,idsi,daxyin,
     &                 ymax,ymin,ymin2,
     &                 xmax,xmin,xmin2,noned,noner,ncom,ierr,jof,jif,
     &                 iycm,rycm,iddy,ided,idef,iffy,iffd,iydd,
     &                 idn,idc,fac,idde,idyd,iyyd,ides,iyhh,iylc,
     &                 iylt,iyss,syss,iysb,iybb,jycm,ycol,ycob,
     &                 inv,ifdc,ivdc,nnvv,idvd,idxd,
     &                 iprjj,idrjj,kdjj,kijj,rkjj,ivijj,ivnjj,ifpjj,
     &                 nevk,nevn,ivi,rk,kd,kp,ki,jd,md,mg,me,mp,
     &                 ivn,ifp)
*                                                                      *
*         MH:      MAX NUMBER OF COLUMN IN ONE SECTION                 *
*         ML:      MAX NUMBER OF ROW IN ONE SECTION                    *
*                                                                      *
*         MC:      MAX NUMBER OF LINE COMMENTS (TOTAL)                 *
*                                                                      *
*         JOF:     OUTPUT FILE OF LINES                                *
*                                                                      *
*         JIL:     OUTPUT FILE OF REGION                               *
*                                                                      *
*         ILL:     PRESENT LINE                                        *
*         ILF:     FINAL LINE OF THE PRESENT INPUT FILE                *
*                                                                      *
*         NONED:   NUMBER OF LINES WRITTEN IN JOF                      *
*         NONER:   NUMBER OF REGIONS WRITTEN IN JIL                    *
*         NCOM:    NUMBER OF LINE COMMENTS                             *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      real(8),intent(in)::  daxyin(mh+2,nnvv+1)
      real(8),allocatable:: daxy(:,:),daxyt(:,:)
      integer ndaxy2, nndaxy2  ! size of 2nd rank in daxy
      integer,parameter:: ipnd2 = 10

      real(8),allocatable:: dspl(:,:)

      dimension ivi(0:mfc,0:mfc)

      dimension rkjj(mh,0:mfc)
      dimension iprjj(mh), idrjj(mh)
      dimension kdjj(mh,0:mfc), kijj(mh,0:mfc)
      dimension ivnjj(mh,0:mfc), ifpjj(mh,0:mfc)
      dimension ivijj(mh,0:mfc,0:mfc)

      dimension nevk(0:mfc,0:mfc), nevn(0:mfc,0:mfc)

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      dimension fac(mh,5)
      dimension iycm(mc,8),rycm(mc,7)
      dimension idc(mh),idyd(-1:mh),iddy(mh),jycm(mh),
     &          iylc(mh),iyss(mh),syss(mh),iylt(mh),ycol(mh,3),
     &          iysb(mh),iyhh(mh),iyyd(mh,2),iybb(mh),ycob(mh,3),
     &          ided(-1:mh),idde(mh),ides(mh),
     &          idef(-1:mh),iffy(-1:mh), iffd(-1:mh), iydd(mh)

      dimension rk(0:mfc)
      dimension ifdc(mh), ivdc(mh)
      dimension kd(0:mfc), kp(0:mfc), ki(0:mfc), jd(0:mfc)
      dimension md(0:mfc), mg(0:mfc), me(0:mfc), mp(0:mfc)
      dimension ivn(0:mfc), ifp(0:mfc)

      common /rval1/ cval(mxcval), aval(mxcval)

      character m_err*200
      common /error/ m_err, l_err, k_err

      dimension ill(0:9), ilf(0:9)
      character dsin(0:9)*200
      dimension idsi(0:9)

      logical dnen2,deqn3

      dimension rcol(3), rcob(3)

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------
*     allocate and initialize daxy
*-----------------------------------------------------------------------

      ndaxy2  = nnvv+1
      nndaxy2 = nnvv+ipnd2
      allocate(daxy(mh+2,nndaxy2))
      daxy(:,1:ndaxy2) = daxyin(:,1:ndaxy2)

*-----------------------------------------------------------------------
*     INITIAL LINE AND X-COLUMN
*-----------------------------------------------------------------------

            illi = ill(jsn)

            ix = idxd

*-----------------------------------------------------------------------
*     READ DATA
*-----------------------------------------------------------------------

      if( inv .eq. 1 ) then

*-----------------------------------------------------------------------
*           ILA; CURRENT NUMBER OF LINES, ILA0; OF REAL ROW
*-----------------------------------------------------------------------

            ila0 = 0
            ila  = 0

  510       ila0 = ila0 + 1
            ila  = ila  + 1

*-----------------------------------------------------------------------


            read(jsi,'(10000a1)', iostat = ios ) (dum(ic),ic=1,icolm)
            if( ios .eq. -1 ) goto 410

            call chlow(dum,lum)

            if( ill(jsn) + ila .gt. ilf(jsn) ) goto 410


*-----------------------------------------------------------------------
*        SKIP THE BLANK
*-----------------------------------------------------------------------

            k = 1

               do 145 l = 1, icolm

                  if(lum(l) .ne. ' ' .and. lum(l) .ne. tub) goto 146

  145          continue

               goto 410

  146       k = l

*-----------------------------------------------------------------------
*        INCLFL: INCLUDE FILE ONLY FOR ILA0 = 1
*-----------------------------------------------------------------------

            if( ila0 .eq. 1 .and.
     &          lum(k)//lum(k+1)//lum(k+2)//lum(k+3)//lum(k+4) .eq.
     &          'infl:' ) then

                  ill(jsn) = ill(jsn) + ila

                  k = k + 4

               call inclf(jsn,jsi,dsin,dum,k,icolm,ill,ilf,idsi,ierr)

                  if( ierr .ne. 0 ) goto 999

               ila0 = 0
               ila  = 0

               goto 510

            end if


*-----------------------------------------------------------------------
*        COMENT LINE STARTED '#'
*-----------------------------------------------------------------------

            if( lum(k) .eq. '#' ) then

               ila0 = ila0 - 1

               goto 510

            end if

*-----------------------------------------------------------------------
*        END OF DATA BY _:, __:, OR NON NUMERIC CHARACTERS
*-----------------------------------------------------------------------

            if( lum(k+1) .eq. ':' .or.
     &          lum(k+2) .eq. ':' .or.
     &          dnen2(lum(k)) ) goto 410

*-----------------------------------------------------------------------
*        CHECK OF ONE LINE AND READ DATA
*-----------------------------------------------------------------------

            ic  = k - 1
            icn = 0

  530       ic = ic + 1

*-----------------------------------------------------------------------

               if( ic .gt. icolm ) then

  570             continue

                  if( icn + 1 .le. idn .and.
     &                ifdc(icn+1) .ne. 0 ) then

                     icn = icn + 1

                     goto 570

                  end if

                     goto 540

               end if

*-----------------------------------------------------------------------

               if( lum(ic) .eq. ' ' .or. lum(ic) .eq. tub ) goto 530

               if( dnen2(lum(ic)) ) then

                  m_err = 'Some Part of the Numeric Date is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:2230/R:hone2/F:a-hsect.f'
                  l_err = ill(jsn) + ila
                  k_err = jsn
                  goto 999

               end if

*-----------------------------------------------------------------------
*        FIND ONE NUMBER
*-----------------------------------------------------------------------

  560          icn = icn + 1

            if( icn .gt. idn ) then

               icn = icn - 1

               goto 540

            end if

            if( ifdc(icn) .ne. 0 ) goto 560

               ici = ic

*-----------------------------------------------------------------------

  550       ic = ic + 1

               if( ic .gt. icolm ) then

                  call rnum(rrnm,lum,ici,ic-1,ierr)

                  if( ierr .ne. 0 ) then

                     m_err = 'Some Part of the Numeric Date is Wrong.'
                     ErrCha = ''
                     ErrID = 'L:2267/R:hone2/F:a-hsect.f'
                     l_err = ill(jsn) + ila
                     k_err = jsn

                     goto 999

                  end if

*-----------------------------------------------------------------------
*                 size check and resize daxy
*-----------------------------------------------------------------------

                  if (ndaxy2 .lt. ila0) then
                    nndaxy2 = ila0+ipnd2
                    allocate(daxyt(mh+2,ndaxy2))
                    daxyt(:,1:ndaxy2) = daxy(:,1:ndaxy2)
                    deallocate(daxy)
                    allocate(daxy(mh+2,nndaxy2))
                    daxy(:,1:ndaxy2) = daxyt(:,1:ndaxy2)
                    deallocate(daxyt)
                    ndaxy2 = nndaxy2
                  end if

*-----------------------------------------------------------------------

                  daxy(icn,ila0) = rrnm

                  goto 540

               end if


               if( lum(ic) .eq. ' ' .or. lum(ic) .eq. tub ) then

                  call rnum(rrnm,lum,ici,ic-1,ierr)

                  if( ierr .ne. 0 ) then

                     m_err = 'Some Part of the Numeric Date is Wrong.'
                     ErrCha = ''
                     ErrID = 'L:2307/R:hone2/F:a-hsect.f'
                     l_err = ill(jsn) + ila
                     k_err = jsn

                     goto 999

                  end if

*-----------------------------------------------------------------------
*                 check size and resize daxy
*-----------------------------------------------------------------------

                  if (ndaxy2 .lt. ila0) then
                    nndaxy2 = ila0+ipnd2
                    allocate(daxyt(mh+2,ndaxy2))
                    daxyt(:,1:ndaxy2) = daxy(:,1:ndaxy2)
                    deallocate(daxy)
                    allocate(daxy(mh+2,nndaxy2))
                    daxy(:,1:ndaxy2) = daxyt(:,1:ndaxy2)
                    deallocate(daxyt)
                    ndaxy2 = nndaxy2
                  end if

*-----------------------------------------------------------------------

                  daxy(icn,ila0) = rrnm

                  goto 530

               end if


               if( deqn3(lum(ic)) ) goto 550

                  m_err = 'Some Part of the Numeric Date is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:2343/R:hone2/F:a-hsect.f'
                  l_err = ill(jsn) + ila
                  k_err = jsn

                  goto 999


*-----------------------------------------------------------------------
*        BLANK LINE IS THE END OF THIS SECTION
*-----------------------------------------------------------------------

  540       if( icn .eq. 0 ) goto 410

*-----------------------------------------------------------------------

               if( icn .ne. idn ) then

                  m_err = 'H: Mismatch of Column Description and Data.'
                  ErrCha = ''
                  ErrID = 'L:2362/R:hone2/F:a-hsect.f'
                  l_err = ill(jsn) + ila
                  k_err = jsn

                  goto 999

               end if

*-----------------------------------------------------------------------

            goto 510

*-----------------------------------------------------------------------
*        CHECK OF THE END OF DATA
*-----------------------------------------------------------------------

  410       ila0 = ila0 - 1
            ila  = ila  - 1

               if( ila0 .eq. 0 ) then

                  m_err = 'H: There is No Data in This H: Section.'
                  ErrCha = ''
                  ErrID = 'L:2385/R:hone2/F:a-hsect.f'
                  l_err = ill(jsn)
                  k_err = jsn

                  goto 999

               end if

*-----------------------------------------------------------------------

            ill(jsn) = ill(jsn) + ila + 1

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*        END OF READ DATA AND FOR FUNCTION V=[1,2,3]
*-----------------------------------------------------------------------

      else

            read(jsi,'(10000a1)', iostat = ios ) (dum(ic),ic=1,icolm)
            if( ios .eq. -1 ) goto 800

            call chlow(dum,lum)

            if( ill(jsn) + 1 .gt. ilf(jsn) ) goto 800

               goto 810

  800       continue

               dum(1) = 'C'
               dum(2) = ':'
               lum(1) = 'c'
               lum(2) = ':'

  810       ill(jsn) = ill(jsn) + 1

            ila0 = nnvv

      end if


*=======================================================================
*     SUMMARY :
*=======================================================================

            ila = ila0

*-----------------------------------------------------------------------
*           size check and resize daxy
*-----------------------------------------------------------------------

            if (ndaxy2 .lt. ila0) then
              nndaxy2 = ila0
              allocate(daxyt(mh+2,ndaxy2))
              daxyt(:,1:ndaxy2) = daxy(:,1:ndaxy2)
              deallocate(daxy)
              allocate(daxy(mh+2,nndaxy2))
              daxy(:,1:ndaxy2) = daxyt(:,1:ndaxy2)
              deallocate(daxyt)
              ndaxy2 = nndaxy2
            end if

*=======================================================================

            do 690 i = 1, ila

                     daxy(idn+1,i) = 0.0

  690       continue

*-----------------------------------------------------------------------
*     SUMMARY : FACTOR FOR NORMAL COLUMN
*-----------------------------------------------------------------------

            do 700 i = 1, ila

               do 710 j = 1, idn

                  if( ivdc(j) .eq. 1 ) then

                     if( fac(j,5) .gt. -r0max )
     &               daxy(j,i) = daxy(j,i)**fac(j,5)

                     daxy(j,i) = daxy(j,i) / fac(j,1) * fac(j,2)
     &                                     - fac(j,3) + fac(j,4)

                  end if

  710          continue

  700       continue

*=======================================================================
*     SUMMARY : FUNCTION WITHOUT A12
*=======================================================================

            do 750 j = 1, idn

               if( ifdc(j) .eq. 2 ) then

                  jdn = j

                  call func11(jdn,kd,ki,rk,ivi,ivn,ifp,ipr,
     &                        kdjj,kijj,rkjj,ivijj,ivnjj,ifpjj,
     &                        iprjj,idrjj)

                  call func08(jdn,ila,daxy,
     &                        ifdc,ivdc,idyd,ided,idef,idvd,idxd,
     &                        ipr,ifp,kd,ki,rk,ivn,ivi,
     &                        ierr,jsn,illi,
     &                        nevk,nevn)

                        if( ierr .ne. 0 ) goto 999

                  ivdc(j) = 1

*-----------------------------------------------------------------------

                  do 720 i = 1, ila

                        if( fac(j,5) .gt. -r0max )
     &                  daxy(j,i) = daxy(j,i)**fac(j,5)

                        daxy(j,i) = daxy(j,i) / fac(j,1) * fac(j,2)
     &                                        - fac(j,3) + fac(j,4)

  720             continue

               end if

  750       continue


*=======================================================================
*     SUMMARY : FUNCTION FITTING
*=======================================================================

            do 751 j = 1, idn

               if( ifdc(j) .eq. 4 ) then

                  jdn = j

                  call func11(jdn,kd,ki,rk,ivi,ivn,ifp,ipr,
     &                        kdjj,kijj,rkjj,ivijj,ivnjj,ifpjj,
     &                        iprjj,idrjj)

                  call func09(jdn,idn,ila,daxy,
     &                        ifdc,ivdc,idyd,ided,idef,idvd,idxd,
     &                        ipr,ifp,kd,ki,rk,ivn,ivi,
     &                        ierr,jsn,illi,iffy,iffd,iyyd,iydd,
     &                        nevk,nevn)

                        if( ierr .ne. 0 ) goto 999

                  ivdc(j) = 1

*-----------------------------------------------------------------------

                  do 740 i = 1, ila

                        if( fac(j,5) .gt. -r0max )
     &                  daxy(j,i) = daxy(j,i)**fac(j,5)

                        daxy(j,i) = daxy(j,i) / fac(j,1) * fac(j,2)
     &                                        - fac(j,3) + fac(j,4)

  740             continue

               end if

  751       continue


*=======================================================================
*     SUMMARY : FUNCTION WITH A12
*=======================================================================

            do 752 j = 1, idn

               if( ifdc(j) .eq. 3 ) then

                  jdn = j

                  call func11(jdn,kd,ki,rk,ivi,ivn,ifp,ipr,
     &                        kdjj,kijj,rkjj,ivijj,ivnjj,ifpjj,
     &                        iprjj,idrjj)

                  call func08(jdn,ila,daxy,
     &                        ifdc,ivdc,idyd,ided,idef,idvd,idxd,
     &                        ipr,ifp,kd,ki,rk,ivn,ivi,
     &                        ierr,jsn,illi,
     &                        nevk,nevn)

                        if( ierr .ne. 0 ) goto 999

                  ivdc(j) = 1

*-----------------------------------------------------------------------

                  do 760 i = 1, ila

                        if( fac(j,5) .gt. -r0max )
     &                  daxy(j,i) = daxy(j,i)**fac(j,5)

                        daxy(j,i) = daxy(j,i) / fac(j,1) * fac(j,2)
     &                                        - fac(j,3) + fac(j,4)

  760             continue

               end if

  752       continue


*=======================================================================
*     SUMMARY : MAXIMUM AND MINMUM VALUE OF X AND Y
*=======================================================================

            do 716 i = 1, ila

               do 717 j = 1, idn

                  if( idc(j) .eq. 1 ) then

                    if( daxy(j,i) .gt. xmax ) xmax = daxy(j,i)
                    if( daxy(j,i) .lt. xmin ) xmin = daxy(j,i)
                    if( daxy(j,i) .lt. xmin2 .and.
     &                  daxy(j,i) .gt. 0.0 ) xmin2 = daxy(j,i)

                  end if

                  if( idc(j) .eq. 2 ) then

                    if( daxy(j,i) .gt. ymax ) ymax = daxy(j,i)
                    if( daxy(j,i) .lt. ymin ) ymin = daxy(j,i)
                    if( daxy(j,i) .lt. ymin2 .and.
     &                  daxy(j,i) .gt. 0.0 ) ymin2 = daxy(j,i)

                  end if


  717          continue

  716       continue


*-----------------------------------------------------------------------
*     WRITE ONEDIM TO jof OR jif AND Y-COMMENTS
*-----------------------------------------------------------------------

            zerr = 0.0

      do 721 i = 1, idn

               if( ycol(i,1) .ge. -r0max ) ycol(i,1) = chue(ycol(i,1))
               if( ycol(i,1) .lt. -r0max ) ycol(i,1) = -3.0
               if( ycob(i,1) .lt. -r0max ) ycob(i,1) = -1.0

*-----------------------------------------------------------------------

         if( idc(i) .eq. 2 .or. idc(i) .eq. 3 ) then

*-----------------------------------------------------------------------
*     NORMAL LINES
*-----------------------------------------------------------------------

            if( iyhh(i) .eq. 0 ) then

*-----------------------------------------------------------------------

                  ityp  = iylc(i)
                    if(iylc(i).eq.-1) ityp = 0

                  iwi   = iylt(i)

                  ipol  = iyss(i)
                  ipols = 0

                  imar  = iysb(i)
                  imars = 0

                  rcol(1) = ycol(i,1)
                  rcol(2) = ycol(i,2)
                  rcol(3) = ycol(i,3)
                  rcob(1) = ycob(i,1)
                  rcob(2) = ycob(i,2)
                  rcob(3) = ycob(i,3)

                  icoma  = 0
                  icomas = 1

                     if(iylc(i).eq.-1) icoma = 1

                  facsiz = 1.0

                    if( iybb(i) .ge. 1 )
     &                 facsiz = 1.5 + dble( iybb(i) - 1 ) * 0.5

                    if( iybb(i) .eq. -1 ) facsiz = 0.8
                    if( iybb(i) .eq. -2 ) facsiz = 0.6
                    if( iybb(i) .eq. -3 ) facsiz = 0.5
                    if( iybb(i) .eq. -4 ) facsiz = 0.4
                    if( iybb(i) .lt. -4 ) facsiz = 0.0001

                     ixdc  = iydd(ix)
                     iydc  = iydd(i)

                     ixdcs = 0
                     iydcs = 0


                     if( ityp .eq. 12 ) then

                        jtyp = 11

                     else

                        jtyp = ityp

                     end if

                     if( jtyp .eq. 11 ) then

                        imar  = 0
                        icoma = 0
                        ixdc  = 0
                        iydc  = 0

                        ddxm  = 0.0
                        ddxp  = 0.0
                        ddym  = 0.0
                        ddyp  = 0.0

                     end if


*-----------------------------------------------------------------------
*           MINIMUM POINTS FOR SPLINE
*-----------------------------------------------------------------------

               if( ipol .ne. 0 .and. ila .lt. 5 .and.
     &             syss(i) .lt. 0.0 ) ipol = 0

               if( ipol .ne. 0 .and. ila .lt. 3 .and.
     &             syss(i) .gt. 0.0 ) ipol = 0

*-----------------------------------------------------------------------
               if( ipol .eq. 1 ) then
*-----------------------------------------------------------------------
*                 **** SPLINE
*-----------------------------------------------------------------------

                  if( icoma .ne. 1 ) then

                     allocate(dspl(15,ila))

                     do 600 k = 1, ila

                        dspl(1,k) = daxy(ix,k)
                        dspl(2,k) = daxy(i,k)

  600                continue

                     call spln00(noned,noner,syss(i),ila,ierr,
     &                        idc(i),ixdcs,iydcs,jtyp,iwi,ipols,
     &                        imars,rcol,rcob,icoma,facsiz,
     &                        jof,jif,ityp,
     &                        dspl)

                     deallocate(dspl)

                     if( ierr .ne. 0 ) goto 910

                  end if

*-----------------------------------------------------------------------
*                 **** WRITE SYMBOL WITH ERROR BAR
*-----------------------------------------------------------------------

                  if( imar .ne. 0 .or.
     &              ( imar .eq. 0 .and.
     &              ( ixdc .ne. 0 .or. iydc .ne. 0 ) ) ) then

                     noned = noned + 1

                     write(jof) ila, idc(i), ixdc, iydc

                     write(jof)
     &               jtyp,iwi,ipols,imar,rcol,rcob,icomas,facsiz

                     do 728 ij = 1, ila

                          ddxm = daxy(iyyd(ix,1),ij)
                          ddxp = daxy(iyyd(ix,2),ij)

                          ddym = daxy(iyyd(i,1),ij)
                          ddyp = daxy(iyyd(i,2),ij)

                          write(jof)
     &                    daxy(ix,ij),daxy(i,ij),ddxm,ddxp,ddym,ddyp

  728                continue

                  end if


                  goto 900

*-----------------------------------------------------------------------
               end if
*-----------------------------------------------------------------------

  910             continue

*-----------------------------------------------------------------------
*                  **** NORMAL
*-----------------------------------------------------------------------

                     ipol = 0

                        if( ityp .ne. 11 ) then

                           noned = noned + 1

                           write(jof) ila, idc(i), ixdc, iydc

                           write(jof)
     &                     jtyp,iwi,ipol,imar,rcol,rcob,icoma,facsiz

                        else

                           noner = noner + 1

                           write(jif) ila, idc(i), ixdc, iydc

                           write(jif)
     &                     jtyp,iwi,ipol,imar,rcol,rcob,icoma,facsiz

                        end if


                     do 727 ij = 1, ila

                        if( jtyp .ne. 11 ) then

                           ddxm = daxy(iyyd(ix,1),ij)
                           ddxp = daxy(iyyd(ix,2),ij)

                           ddym = daxy(iyyd(i,1),ij)
                           ddyp = daxy(iyyd(i,2),ij)

                        end if

                        if( ityp .ne. 11 ) then

                           write(jof)
     &                     daxy(ix,ij),daxy(i,ij),ddxm,ddxp,ddym,ddyp

                        else

                           write(jif)
     &                     daxy(ix,ij),daxy(i,ij),ddxm,ddxp,ddym,ddyp

                        end if

  727                continue

*-----------------------------------------------------------------------

  900             continue


*-----------------------------------------------------------------------
*     HISTOGRAM
*-----------------------------------------------------------------------

            else if( iyhh(i) .gt. 0 ) then

*-----------------------------------------------------------------------

                  ityp  = iylc(i)
                    if(iylc(i).eq.-1) ityp = 0

                  iwi   = iylt(i)
                  ipol  = 0
                  imar  = 0

                  rcol(1) = ycol(i,1)
                  rcol(2) = ycol(i,2)
                  rcol(3) = ycol(i,3)

                  icoma = 0
                  facsiz = 0.0

                  ixdc = 0
                  iydc = 0

                noned = noned + 1


              if(iyhh(i) .le. 2 ) then

                write(jof) 2*ila-1,idc(i), ixdc, iydc

              else if(iyhh(i) .ge. 3 ) then

                write(jof) 3*ila-2,idc(i), ixdc, iydc

              end if


                write(jof) ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz

                write(jof) daxy(ix,1),daxy(i,1),zerr,zerr,zerr,zerr

             do 723 ij=2,ila

               if( iyhh(i) .eq. 1 ) then

                write(jof) daxy(ix,ij),daxy(i,ij-1),
     &                      zerr,zerr,zerr,zerr
                write(jof) daxy(ix,ij),daxy(i,ij),
     &                      zerr,zerr,zerr,zerr

               else if( iyhh(i) .eq. 2 ) then

                write(jof) daxy(ix,ij-1),daxy(i,ij),
     &                     zerr,zerr,zerr,zerr
                write(jof) daxy(ix,ij),daxy(i,ij),
     &                     zerr,zerr,zerr,zerr

               else if( iyhh(i) .ge. 3 ) then

                     xmidl = ( daxy(ix,ij-1) + daxy(ix,ij) ) / 2.0

                write(jof) xmidl,daxy(i,ij-1),
     &                     zerr,zerr,zerr,zerr
                write(jof) xmidl,daxy(i,ij),
     &                     zerr,zerr,zerr,zerr
                write(jof) daxy(ix,ij),daxy(i,ij),
     &                     zerr,zerr,zerr,zerr

               end if

  723        continue

*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------
*     Y COMMENTS
*-----------------------------------------------------------------------

           if( jycm(i) .ne. 0 ) then

                  ncom = ncom + 1

                  iycm(ncom,1) = jycm(i)
                  iycm(ncom,6) = 0

                  if(iyhh(i).ne.0) iycm(ncom,6) = 1

                  iycm(ncom,7) = iydd(i)
                  iycm(ncom,8) = iydd(ix)

                  iycm(ncom,2) = ityp
                  iycm(ncom,3) = iwi
                  iycm(ncom,4) = imar
                  iycm(ncom,5) = icoma

                  rycm(ncom,1) = facsiz
                  rycm(ncom,2) = rcol(1)
                  rycm(ncom,3) = rcol(2)
                  rycm(ncom,4) = rcol(3)
                  rycm(ncom,5) = rcob(1)
                  rycm(ncom,6) = rcob(2)
                  rycm(ncom,7) = rcob(3)

           end if

*-----------------------------------------------------------------------
*     HISTOGRAM + ERROR BARS
*-----------------------------------------------------------------------

            if( iyhh(i) .gt. 0 ) then

              if( iyyd(i,1) .ne. idn+1 .or.
     &            iyyd(i,2) .ne. idn+1 ) then

                     ityp  = iylc(i)
                       if(iylc(i).eq.-1) ityp = 0

                     iwi   = iylt(i)
                     ipol  = 0

*                 -----------------
                     imar  = -1
*                 -----------------

                     rcol(1) = ycol(i,1)
                     rcol(2) = ycol(i,2)
                     rcol(3) = ycol(i,3)
                     rcob(1) = ycob(i,1)
                     rcob(2) = ycob(i,2)
                     rcob(3) = ycob(i,3)

                     icoma = 1

                     facsiz = 1.0

                     ixdc = 0
                     iydc = iydd(i)

                     noned = noned + 1


                     write(jof) ila-1,idc(i), ixdc, iydc

                     write(jof)
     &               ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz


                  do 724 ij = 2, ila


                     if( iyhh(i) .eq. 1 ) then

                        ijj = ij - 1
                        xmidl = (daxy(ix,ij)+daxy(ix,ij-1))/2.0

                     else if( iyhh(i) .eq. 2 ) then

                        ijj = ij
                        xmidl = (daxy(ix,ij)+daxy(ix,ij-1))/2.0

                     else if( iyhh(i) .ge. 3 ) then

                        ijj = ij
                        xmidl = daxy(ix,ij)

                     end if

                        ddym = daxy(iyyd(i,1),ijj)
                        ddyp = daxy(iyyd(i,2),ijj)

                     if( daxy(ix,ij) .ne. daxy(ix,ij-1) ) then

                        write(jof)
     &                  xmidl, daxy(i,ijj), zerr,zerr,ddym,ddyp

                     else

                        write(jof)
     &                  xmidl, daxy(i,ijj), zerr,zerr,zerr,zerr

                     end if

  724              continue


              end if

            end if

*-----------------------------------------------------------------------

         end if

  721 continue

*-----------------------------------------------------------------------

      deallocate(daxy)
      return

  999 ierr = 1

      if(allocated(daxy)) deallocate(daxy)
      return
      end


