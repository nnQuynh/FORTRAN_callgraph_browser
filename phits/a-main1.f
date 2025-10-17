************************************************************************
*                                                                      *
      function zonep(xmin,xmax,xps,ilog)
*                                                                      *
*              TRANSLATE THE COODINATE TO ZERO-ONE VALUE OF AXIS       *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

            if( ilog .eq. 0 ) then

               zonep = ( xps - xmin ) / ( xmax - xmin )

            else if( ilog .eq. 1 ) then

               if( xps .gt. 0.0 ) then

               zonep = 1.0d0
     &               * ( log10(xps)  - log10(xmin) )
     &               / ( log10(xmax) - log10(xmin) )

               else
               zonep = -1.0d0
               end if

            end if


      return
      end

************************************************************************
*                                                                      *
      subroutine para(dum,lum,in,form,xpmin,xpmax,ypmin,ypmax,
     &                afac,iylog,ixlog,atxs,
     &                ispac,ispax,ispay,icut,rcut,ncut,
     &                ixtic,iytic,nocm,notc,notn,notf,
     &                itic,iwd2,ipd2,ipds,ipdc,
     &                scal,noms,nolg,notl,ibsw,
     &                jxtic,jytic,iccr,
     &                regx,regy,regs,regd,irot,xorg,yorg,xfac,nofr,
     &                noxt,noxn,noyt,noyn,ierr,izdo,icct,ibfon,itfon,
     &                clal,clax,cltx,clnm,cltl,cllg,clfr,clms,clhd,
     &                clbg,clin,clmo,izlog,isec,noaf,angz,
     &                ixstd,iystd,ixltd,iyltd,xsd,xld,ysd,yld,
     &                mdecx,mdecy,xstv,xltv,ystv,yltv,rlptl,sybw,idbg,
     &                ixnum,iynum,ixtxt,iytxt,xtxp,ytxp,tlxp,ifon,idat,
     &                clgb,clgl,clgs,ilbox,
     &                nopa,parx,pary,pars,clpa,ipbox,clpb,clpl,clps,
     &                nocn,conx,cony,cons,clcn,isbox,clnb,clnl,clns,
     &                xmul,ymul,zmul,icmcg,itncg,icmyy,itnyy,isccg,
     &                erwd,ixexp,iyexp,
cKN 2024/01/24
     &                ih2fs)
*                                                                      *
*                                                                      *
*        FORMAT: P: NOSP CUTS(1.5, 2.5, 3.0) FORM(1.2)                 *
*                                                                      *
*          IZDO     : NUMBER OF MULTI RECORD                           *
*                                                                      *
*          IBSW     : SHOW BOUNDING BOX FOR EACH PAGE                  *
*                                                                      *
*          XMIN( ), XMAX( ), YMIN( ), YMAX( )                          *
*                                                                      *
*          IROT     : =1; LANDSCAPE(DEFAULT) =0; PORTRATE              *
*                     -1; SMALL LANDSCAPE                              *
*          XORG( )  : POSITION OF THE X-ORIGIN, UNIT IS X-AXIS         *
*          YORG( )  : POSITION OF THE Y-ORIGIN, UNIT IS Y-AXIS         *
*          ANGL( )  : ANGLE OF THE GRAPH (DEGREE)                      *
*          XFAC( )  : LENGTH OF THE X-AXIS REL.TO DEFAULT              *
*                                                                      *
*         IBFON( )  : BASIC FONT                                       *
*         ITFON( )  : TIC FONT                                         *
*                                                                      *
*          LEGX( )  : X-POSITION OF LEGEND                             *
*          LEGY( )  : Y-POSITION OF LEGEND                             *
*          LEGS( )  : SIZE OF LEGEND ( DEFAULT =1.0 )                  *
*                                                                      *
*          LBOX( )  : BOX FOR LEGEND                                   *
*          LBCB( )  : COLOR OF BACKGROUND OF LEGEND BOX                *
*          LBCL( )  : COLOR OF LINE OF LEGEND BOX                      *
*          LBCS( )  : COLOR OF SHADOW OF LEGEND BOX                    *
*                                                                      *
*          CLLG( )  : COLOR OF LEGEND TEXT                             *
*                                                                      *
*          NOLG     : NO LEGEND IF DEFAULT POSITION                    *
*          NOTL     : NO TITLE IF EXIST                                *
*                                                                      *
*          NOXT     : NO X-AXIS TEXT                                   *
*          NOXN     : NO X-AXIS NUMBER                                 *
*          NOYT     : NO Y-AXIS TEXT                                   *
*          NOYN     : NO Y-AXIS NUMBER                                 *
*                                                                      *
*          MONO     : CHANGE COLOR TO MONO                             *
*          CLAL( )  : COLOR OF ALL                                     *
*          CLAX( )  : COLOR OF AXIS                                    *
*          CLTX( )  : COLOR OF AXIS TEXT                               *
*          CLNM( )  : COLOR OF AXIS NUMBER                             *
*          CLTL( )  : COLOR OF GRAPH TITLE                             *
*          CLFR( )  : COLOR OF FRAME                                   *
*          CLMS( )  : COLOR OF MESSAGE OF FILE NAME AND DATE           *
*          CLHD( )  : COLOR OF CLUSTER PLOT                            *
*          CLBG( )  : COLOR OF BACK GROUND                             *
*          CLIN( )  : COLOR OF INSIDE THE AXIS                         *
*                                                                      *
*          SCAL( )  : SCALING OF THE FIGURE FOR SPECIAL.DAT            *
*          FORM( )  : RATIO Y-AXIS TO X-AXIS                           *
*          AFAC( )  : HIGHT OF THE AXIS TEXT REL TO NORMAL (MAX 2.0)   *
*                                                                      *
*          ATXS( )  : SIZE OF AXIS TEXT (DEFAULT= 1)                   *
*                                                                      *
*          YLOG, YLIN, XLOG, XLIN,                                     *
*                                                                      *
*          SPAC, NOSP  : MARGIN OR NO MARGIN                           *
*          NOMS     :  NO MSSAGES OF FILE NAME AND DATE                *
*          NOCM     :  NO COMMENTS IF EXIT                             *
*          NOFR     :  NO FRAME                                        *
*          FRAM     :  FRAME (DEFAULT)                                 *
*                                                                      *
*          CUTS(A, B, C) : CUTS OF CONTOUR PLOT                        *
*          COLS(A, B, C) : COLORS  OF CONTOUR PLOT                     *
*          COLN          : NO COLOR REGION OF CONTOUR PLOT             *
*          COLC          : COLOR REGION  OF CONTOUR PLOT               *
*          ICUT          : DEFAULT NUMBER OF CUTS                      *
*                                                                      *
*          XDTC( )  : 0=NORMAL, 1=LARGE, -1=SMALL  FOR AUTO            *
*          YDTC( ),                                                    *
*                                                                      *
*          IXSTD    : =0 AUTO, =1 MANUAL                               *
*          IYSTD, IXLTD, IYLTD                                         *
*                                                                      *
*          XSTD( )  : >0 SHORT TIC DISTANCE; <0 NO SHORT TIC           *
*          YSTD( )    XSD, YSD                                         *
*                                                                      *
*          XLTD( )  : >0 LONG TIC DISTANCE; <0 NO SHORT AND LONG TIC   *
*          YLTD( )    XLD, YLD                 AND NO NUMBER           *
*                                                                      *
*          MDECX( ) : = -10 AUTO, OTHER;  MANUAL                       *
*          MDECY( )   XTDC, YTDC                                       *
*                                                                      *
*          XSTV( )  : = EXPLICIT ONE VALUE OF THE TICS                 *
*          XLTV( )    YSTV( ) YLTV( )                                  *
*                                                                      *
*                                                                      *
*          ITIC( )  : 1=INSIDE (DEFAULT), -1=OUTSIDE                   *
*                                                                      *
*          BFON( )  : BASIC FONT                                       *
*                     1: Helveticsa, 2: Times-Roman, 3: Courier        *
*          TFON( )  : TIC FONT                                         *
*                     1: Helveticsa, 2: Times-Roman, 3: Courier        *
*                                                                      *
*          CLUS( )  : STORENGTH OF CLUSTER PLOT                        *
*                                                                      *
*          IWD2( )  : WIDTH OF CONTOUR LINE 1/300 INCH, DEFAULT=4      *
*                                                                      *
*          IPD2     : SPLINE FOR CONTOUR PLOT, DEFAULT IS NO           *
*          IPDS     : smoothing FOR CONTOUR PLOT, DEFAULT IS NO        *
*                                                                      *
*          IPDC     : DUPLICATE OF COLOR AND MONOCLE CLUSTERS          *
*                                                                      *
*         CMAX( )   : MAX VALUE OF CLUSTER AND COLOR PLOT              *
*         CMIN( )   : MIN VALUE OF CLUSTER AND COLOR PLOT              *
*                                                                      *
*         DMAX( )   : MAX CUT OFF VALUE OF CLUSTER AND COLOR PLOT      *
*         DMIN( )   : MIN CUT OFF VALUE OF CLUSTER AND COLOR PLOT      *
*                                                                      *
*         LPTL( )   : LINE PATTERN LENGTH                              *
*                     < 0 -> DEFAULT * SQRT(IWI/4)                     *
*                     > 0 -> DEFAULT * LPTL( )                         *
*                                                                      *
*         SYBW( )   : LINE WIDTH OF SYMBOLS                            *
*                     < 0 -> ACCORDING TO THE SIZE OF SYMBOL           *
*                            X ; -1,  A ; +1(5), AA ; +2(7)            *
*                     > 0 -> LINE WIDTH IS SYBW( ) DD                  *
*                                                                      *
*         IDBG      : 0 ; DEFAULT, NE 0 ; DEBUG MODE                   *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /h2/ xdd,ydd,
     &            xymax,xymin,clus,cmax,cmin,cxymax,cxymin,xyabs,
     &            dmax,dmin,ixnm,iynm

      common /h3/ smax(3), smin(3)
      common /pap/  a4w, a4h, wmg, hmg, wct, hct

      character dum(ichrl)*1
      character lum(ichrl)*1
      character dum4*4
      dimension rcut(mc,5)
      dimension iwd2(mc)
      logical dnen2,deqn1,deqn3

      character c1*1,c2*1,c3*1

      character rum*20

      common /error/ m_err, l_err, k_err
      character m_err*200

      dimension  clal(3), clax(3), cltx(3), clnm(3), cltl(3), cllg(3),
     &           clpa(3), clcn(3), clfr(3), clms(3), clgb(3), clgl(3),
     &           clgs(3), clpb(3), clpl(3), clps(3), clnb(3), clnl(3),
     &           clns(3), clhd(3), clbg(3), clin(3)

      dimension rcol(3)

      common /cmap/ ndis, icrev, cmap
      character cmap*99

      character tub*1
      tub = char(9)

*----------------------------------------------------------------------*

         c1 = '('
         c2 = ')'

         cm = 28.346457

*----------------------------------------------------------------------*

      ic=in-1

  100 ic=ic+1
        if(ic.gt.icolm) goto 500
        if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 100
        if(ic+3.gt.icolm) goto 999

         do k = ic, ic + 3

            dum4(k-ic+1:k-ic+1) = lum(k)

         end do

        ic = ic + 3

*----------------------------------------------------------------------*

        if(dum4.eq.'date') then
              ic=ic+1
           call pkan(dum,ic,icolm,idat,ierr)
              if(ierr.ne.0) goto 999

               if( ifon .eq. -1 .and. idat .gt. 0 ) ifon = idat

              ic=ic+1

        else if(dum4.eq.'tlxp') then
              ic=ic+1
           call pnum(lum,ic,icolm,tlxp,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'xtxp') then
              ic=ic+1
           call pnum(lum,ic,icolm,xtxp,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'ytxp') then
              ic=ic+1
           call pnum(lum,ic,icolm,ytxp,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'xstv') then
              ic=ic+1
           call pnum(lum,ic,icolm,xstv,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'xltv') then
              ic=ic+1
           call pnum(lum,ic,icolm,xltv,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'ystv') then
              ic=ic+1
           call pnum(lum,ic,icolm,ystv,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'yltv') then
              ic=ic+1
           call pnum(lum,ic,icolm,yltv,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'xstd') then
              ic=ic+1
              ixstd=1
           call pnum(lum,ic,icolm,xsd,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'ystd') then
              ic=ic+1
              iystd=1
           call pnum(lum,ic,icolm,ysd,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'xltd') then
              ic=ic+1
              ixltd=1
           call pnum(lum,ic,icolm,xld,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'yltd') then
              ic=ic+1
              iyltd=1
           call pnum(lum,ic,icolm,yld,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

*----------------------------------------------------------------------*
        else if(dum4.eq.'h2fs') then
              ih2fs = 1
              ic=ic+1

        else if(dum4.eq.'h2nr') then
              ih2fs = 0
              ic=ic+1

        else if(dum4.eq.'cmyy') then
              icmyy = 1
              ic=ic+1

        else if(dum4.eq.'nsyy') then
              itnyy = 1
              ic=ic+1

        else if(dum4.eq.'scsc') then
              isccg = -1
              ic=ic+1

        else if(dum4.eq.'scmn') then
              isccg = 1
              ic=ic+1

        else if(dum4.eq.'schr') then
              isccg = 2
              ic=ic+1

        else if(dum4.eq.'scdy') then
              isccg = 3
              ic=ic+1

        else if(dum4.eq.'scyr') then
              isccg = 4
              ic=ic+1

        else if(dum4.eq.'cmnm') then
              icmcg = 7
              ic=ic+1

        else if(dum4.eq.'cmum') then
              icmcg = 4
              ic=ic+1

        else if(dum4.eq.'cmmm') then
              icmcg = 1
              ic=ic+1

        else if(dum4.eq.'cmmt') then
              icmcg = -2
              ic=ic+1

        else if(dum4.eq.'cmkm') then
              icmcg = -5
              ic=ic+1

        else if(dum4.eq.'nsps') then
              itncg = 3
              ic=ic+1

        else if(dum4.eq.'nsus') then
              itncg = -3
              ic=ic+1

        else if(dum4.eq.'nsms') then
              itncg = -6
              ic=ic+1

        else if(dum4.eq.'nssc') then
              itncg = -9
              ic=ic+1

*----------------------------------------------------------------------*
cKN 2024/01/24
        else if(dum4.eq.'erwd') then
              ic=ic+1
           call pnum(lum,ic,icolm,erwd,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'xdec') then
              ixexp = 0
              ic=ic+1
        else if(dum4.eq.'xexp') then
              ixexp = 1
              ic=ic+1
        else if(dum4.eq.'ydec') then
              iyexp = 0
              ic=ic+1
        else if(dum4.eq.'yexp') then
              iyexp = 1
              ic=ic+1

*----------------------------------------------------------------------*
        else if(dum4.eq.'xmul') then
              ic=ic+1
           call pnum(lum,ic,icolm,xmul,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'ymul') then
              ic=ic+1
           call pnum(lum,ic,icolm,ymul,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'zmul') then
              ic=ic+1
           call pnum(lum,ic,icolm,zmul,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

*----------------------------------------------------------------------*

        else if(dum4.eq.'icut') then
              ic=ic+1
           call pnum(lum,ic,icolm,rncut,ierr)
              if(ierr.ne.0) goto 999
              ncut=nint(rncut)
              if(ncut.gt.mc) goto 999
              ic=ic+1

        else if(dum4.eq.'sybw') then
              ic=ic+1
           call pnum(lum,ic,icolm,sybw,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'lptl') then
              ic=ic+1
           call pnum(lum,ic,icolm,rlptl,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'legx') then
              ic=ic+1
           call pnum(lum,ic,icolm,regx,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'legy') then
              ic=ic+1
           call pnum(lum,ic,icolm,regy,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'legs') then
              ic=ic+1
           call pnum(lum,ic,icolm,regs,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'legd') then
              ic=ic+1
           call pnum(lum,ic,icolm,regd,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'parx') then
              ic=ic+1
           call pnum(lum,ic,icolm,parx,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'pary') then
              ic=ic+1
           call pnum(lum,ic,icolm,pary,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'pars') then
              ic=ic+1
           call pnum(lum,ic,icolm,pars,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'conx') then
              ic=ic+1
           call pnum(lum,ic,icolm,conx,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'cony') then
              ic=ic+1
           call pnum(lum,ic,icolm,cony,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'cons') then
              ic=ic+1
           call pnum(lum,ic,icolm,cons,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'xfac') then
              ic=ic+1
           call pnum(lum,ic,icolm,xfac,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'xorg') then
              ic=ic+1
           call pnum(lum,ic,icolm,xorg,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'yorg') then
              ic=ic+1
           call pnum(lum,ic,icolm,yorg,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'angl') then
              ic=ic+1
           call pnum(lum,ic,icolm,angz,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'scal') then
              ic=ic+1
           call pnum(lum,ic,icolm,scal,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'form') then
              ic=ic+1
           call pnum(lum,ic,icolm,form,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'xmin') then
              ic=ic+1
           call pnum(lum,ic,icolm,xpmin,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'xmax') then
              ic=ic+1
           call pnum(lum,ic,icolm,xpmax,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'ymax') then
              ic=ic+1
           call pnum(lum,ic,icolm,ypmax,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'ymin') then
              ic=ic+1
           call pnum(lum,ic,icolm,ypmin,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'afac') then
              ic=ic+1
           call pnum(lum,ic,icolm,afac,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'atxs') then
              ic=ic+1
           call pnum(lum,ic,icolm,atxs,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'clus') then
              ic=ic+1
           call pnum(lum,ic,icolm,clus,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'cmax') then
              ic=ic+1
           call pnum(lum,ic,icolm,cmax,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'cmin') then
              ic=ic+1
           call pnum(lum,ic,icolm,cmin,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'dmax') then
              ic=ic+1
           call pnum(lum,ic,icolm,dmax,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'dmin') then
              ic=ic+1
           call pnum(lum,ic,icolm,dmin,ierr)
              if(ierr.ne.0) goto 999
              ic=ic+1

        else if(dum4.eq.'bfon') then
              ic=ic+1
           call pnum(lum,ic,icolm,bfon,ierr)
              if(ierr.ne.0) goto 999
              ibfon=nint(bfon)

              if( ibfon .le. -1 .or. ibfon .ge. 13 ) then
                  ibfon = 0
              end if

              ic=ic+1

        else if(dum4.eq.'tfon') then
              ic=ic+1
           call pnum(lum,ic,icolm,tfon,ierr)
              if(ierr.ne.0) goto 999
              itfon=nint(tfon)

              if( itfon .le. -1 .or. itfon .ge. 13 ) then
                  itfon = 0
              end if

              ic=ic+1

        else if(dum4.eq.'itic') then
              ic=ic+1
           call pnum(lum,ic,icolm,ztic,ierr)
              if(ierr.ne.0) goto 999
              itic=nint(ztic)

              if(itic.eq.0) then
              else if(itic.ne.1.and.itic.ne.-1) then
                  itic = 1
              end if

              ic=ic+1

        else if(dum4.eq.'xnum') then
              ic=ic+1
           call pnum(lum,ic,icolm,zxnum,ierr)
              if(ierr.ne.0) goto 999
              ixnum=nint(zxnum)
              if(ixnum.lt.-1.or.ixnum.gt.2)
     &        ixnum=1
              if(ixnum.eq.0) then
                  noxn=0
              else
                  noxn=1
              end if
              ic=ic+1

        else if(dum4.eq.'ynum') then
              ic=ic+1
           call pnum(lum,ic,icolm,zynum,ierr)
              if(ierr.ne.0) goto 999
              iynum=nint(zynum)
              if(iynum.lt.-1.or.iynum.gt.2)
     &        iynum=1
              if(iynum.eq.0) then
                  noyn=0
              else
                  noyn=1
              end if
              ic=ic+1

        else if(dum4.eq.'xtxt') then
              ic=ic+1
           call pnum(lum,ic,icolm,zxtxt,ierr)
              if(ierr.ne.0) goto 999
              ixtxt=nint(zxtxt)
              if(ixtxt.lt.-1.or.ixtxt.gt.2)
     &        ixtxt=1
              if(ixtxt.eq.0) then
                  noxt=0
              else
                  noxt=1
              end if
              ic=ic+1

        else if(dum4.eq.'ytxt') then
              ic=ic+1
           call pnum(lum,ic,icolm,zytxt,ierr)
              if(ierr.ne.0) goto 999
              iytxt=nint(zytxt)
              if(iytxt.lt.-1.or.iytxt.gt.2)
     &        iytxt=1
              if(iytxt.eq.0) then
                  noyt=0
              else
                  noyt=1
              end if
              ic=ic+1

        else if(dum4.eq.'xtdc') then
              ic=ic+1
           call pnum(lum,ic,icolm,xtdc,ierr)
              if(ierr.ne.0) goto 999
              mdecx=nint(xtdc)
              ic=ic+1

        else if(dum4.eq.'ytdc') then
              ic=ic+1
           call pnum(lum,ic,icolm,ytdc,ierr)
              if(ierr.ne.0) goto 999
              mdecy=nint(ytdc)
              ic=ic+1

        else if(dum4.eq.'xdtc') then
              ic=ic+1
           call pnum(lum,ic,icolm,xtic,ierr)
              if(ierr.ne.0) goto 999
              ixtic=nint(xtic)
              if(ixtic.ne.0.and.ixtic.ne.1.and.ixtic.ne.-1)
     &        ixtic=0
              ic=ic+1

        else if(dum4.eq.'ydtc') then
              ic=ic+1
           call pnum(lum,ic,icolm,ytic,ierr)
              if(ierr.ne.0) goto 999
              iytic=nint(ytic)
              if(iytic.ne.0.and.iytic.ne.1.and.iytic.ne.-1)
     &        iytic=0
              ic=ic+1

        else if(dum4.eq.'xtic') then
              ic=ic+1
           call pnum(lum,ic,icolm,xtic,ierr)
              if(ierr.ne.0) goto 999
              jxtic=nint(xtic)
              if(jxtic.gt.2.or.jxtic.lt.-1)
     &        jxtic=1
              ic=ic+1

        else if(dum4.eq.'ytic') then
              ic=ic+1
           call pnum(lum,ic,icolm,ytic,ierr)
              if(ierr.ne.0) goto 999
              jytic=nint(ytic)
              if(jytic.gt.2.or.jytic.lt.-1)
     &        jytic=0
              ic=ic+1

*----------------------------------------------------------------------*

        else if(dum4.eq.'a4us') then
              a4w = 21.59  * cm
              a4h = 27.94  * cm
              ic=ic+1

        else if(dum4.eq.'a3pp') then
              a4w = 29.7  * cm
              a4h = 42.0  * cm
              ic=ic+1

        else if(dum4.eq.'a4pp') then
              a4w = 21.0  * cm
              a4h = 29.7  * cm
              ic=ic+1

        else if(dum4.eq.'a5pp') then
              a4w = 14.8  * cm
              a4h = 21.0  * cm
              ic=ic+1

        else if(dum4.eq.'b3pp') then
              a4w = 35.3  * cm
              a4h = 50.0  * cm
              ic=ic+1

        else if(dum4.eq.'b4pp') then
              a4w = 25.0  * cm
              a4h = 35.3  * cm
              ic=ic+1

        else if(dum4.eq.'b5pp') then
              a4w = 17.6  * cm
              a4h = 25.0  * cm
              ic=ic+1

*----------------------------------------------------------------------*

        else if(dum4.eq.'bbsw') then
              if(izdo.eq.1) ibsw=1
              ic=ic+1

        else if(dum4.eq.'nobb') then
              if(izdo.eq.1) ibsw=0
              ic=ic+1

        else if(dum4.eq.'debg') then
              if(izdo.eq.1) idbg=1
              ic=ic+1

        else if(dum4.eq.'nodb') then
              if(izdo.eq.1) idbg=0
              ic=ic+1

        else if(dum4.eq.'notf') then
              notf=0
              ic=ic+1

        else if(dum4.eq.'notn') then
              notn=0
              ic=ic+1

        else if(dum4.eq.'notc') then
              notc=0
              ic=ic+1

        else if(dum4.eq.'noyn') then
              noyn=0
              ic=ic+1

        else if(dum4.eq.'noyt') then
              noyt=0
              ic=ic+1

        else if(dum4.eq.'noxn') then
              noxn=0
              ic=ic+1

        else if(dum4.eq.'noxt') then
              noxt=0
              ic=ic+1

        else if(dum4.eq.'nocm') then
              nocm=0
              ic=ic+1

        else if(dum4.eq.'notl') then
              notl=0
              ic=ic+1

        else if(dum4.eq.'titl') then
              notl=1
              ic=ic+1

        else if(dum4.eq.'nolg') then
              nolg=0
              ic=ic+1

        else if(dum4.eq.'legn') then
              nolg=1
              ic=ic+1

        else if(dum4.eq.'nopa') then
              nopa=0
              ic=ic+1

        else if(dum4.eq.'para') then
              nopa=1
              ic=ic+1

        else if(dum4.eq.'nocn') then
              nocn=0
              ic=ic+1

        else if(dum4.eq.'cnst') then
              nocn=1
              ic=ic+1

        else if(dum4.eq.'noms') then
              noms=1
              ic=ic+1

        else if(dum4.eq.'mssg') then
              noms=0
              ic=ic+1

        else if(dum4.eq.'fram') then
              nofr=0
              ic=ic+1

        else if(dum4.eq.'nofr') then
              nofr=1
              ic=ic+1

        else if(dum4.eq.'irot') then
              if( izdo .eq. 1 ) irot = 0
              ic=ic+1

        else if(dum4.eq.'land') then
              if( izdo .eq. 1 ) irot = 1
              ic=ic+1

        else if(dum4.eq.'port') then
              if( izdo .eq. 1 ) irot = 0
              ic=ic+1

        else if(dum4.eq.'slnd') then
              if( izdo .eq. 1 ) irot = -1
              ic=ic+1


        else if(dum4.eq.'noaf') then
              noaf = 0
              ic=ic+1

        else if(dum4.eq.'afrm') then
              noaf = 1
              ic=ic+1

        else if(dum4.eq.'secp') then

            if( izdo .eq. 1 ) then

               ic = ic + 1

                  scal  = 1.0
                  xpmin = 0.0
                  xpmax = 16.0
                  ypmin = 0.0
                  ypmax = 24.0
                  irot  = 0
                  nofr  = 1
                  noms  = 1
                  noxt  = 0
                  noyt  = 0
                  ispac = 1
                  xfac  = 1.142857
                  form  = 1.5
                  xorg  = -0.158
                  yorg  = -0.0878
                  itic  = -1
                  afac  = 0.4
                  ixltd = 1
                  iyltd = 1
                  ixstd = 1
                  iystd = 1
                  xld   = 1.0
                  yld   = 1.0
                  xsd   = 0.1
                  ysd   = 0.1
                  ixnum = 2
                  iynum = 2

                  isec  = 1

            end if


        else if(dum4.eq.'secl') then

            if( izdo .eq. 1 ) then

               ic = ic + 1

                  scal  = 1.0
                  xpmin = 0.0
                  xpmax = 24.0
                  ypmin = 0.0
                  ypmax = 16.0
                  irot  = 1
                  nofr  = 1
                  noms  = 1
                  noxt  = 0
                  noyt  = 0
                  ispac = 1
                  xfac  = 1.71429
                  form  = 0.6666667
                  xorg  = -0.129
                  yorg  = -0.095
                  itic  = -1
                  afac  = 0.4
                  ixltd = 1
                  iyltd = 1
                  ixstd = 1
                  iystd = 1
                  xld   = 1.0
                  yld   = 1.0
                  xsd   = 0.1
                  ysd   = 0.1
                  ixnum = 2
                  iynum = 2

                  isec  = 2

            end if


        else if(dum4.eq.'ipd2') then

              ipd2 = 4

              ic = ic + 1

               if( lum(ic) .eq. '[' ) then

                  call pnum(lum,ic,ichrl,syrn,ierr)

                     if( ierr .ne. 0 ) goto 999

                  ipd2 = nint(syrn)

                  if( ipd2 .eq. 0 ) goto 999

                  ic = ic + 1

               end if

        else if(dum4.eq.'ipds') then

              ipds = 1

              ic = ic + 1

               if( lum(ic) .eq. '[' ) then

                  call pnum(lum,ic,ichrl,syrn,ierr)

                     if( ierr .ne. 0 ) goto 999

                  ipds = nint(syrn)

                  if( ipds .eq. 0 ) goto 999

                  ic = ic + 1

               end if


        else if(dum4.eq.'npdc') then
              ipdc=0
              ic=ic+1

        else if(dum4.eq.'ipdc') then
              ipdc=1
              ic=ic+1

        else if(dum4.eq.'ylog') then
              iylog=1
              ic=ic+1

        else if(dum4.eq.'ylin') then
              iylog=0
              ic=ic+1

        else if(dum4.eq.'xlog') then
              ixlog=1
              ic=ic+1

        else if(dum4.eq.'xlin') then
              ixlog=0
              ic=ic+1

        else if(dum4.eq.'zlog') then
              izlog=1
              ic=ic+1

        else if(dum4.eq.'zlin') then
              izlog=0
              ic=ic+1

        else if(dum4.eq.'spac') then
              ispac=0
              ic=ic+1

        else if(dum4.eq.'nosp') then
              ispac=1
              ic=ic+1

        else if(dum4.eq.'nosx') then
              ispax=1
              ic=ic+1

        else if(dum4.eq.'nosy') then
              ispay=1
              ic=ic+1

        else if(dum4.eq.'coln') then

              ic = ic + 1

               if( lum(ic) .eq. '(' ) then

                  call pnum(lum,ic,ichrl,vvv,ierr)

                     if( ierr .ne. 0 ) goto 999

                  iccr = nint( vvv )

                  iccr = min(3,max(1,iccr))

                  ic = ic + 1

               end if


        else if(dum4.eq.'mono') then

              clmo = 1.0

              ic = ic + 1

               if( lum(ic) .eq. '(' ) then

                  call pnum(lum,ic,ichrl,vvv,ierr)

                     if( ierr .ne. 0 ) goto 999

                  clmo = min(1.0d0,max(0.0d0,vvv))

                  ic = ic + 1

               end if

        else if(dum4.eq.'colo') then

              clmo = -r1max

              ic=ic+1

*----------------------------------------------------------------------*

        else if(dum4.eq.'cuts') then

           do 200 i=1,mc
             rcut(i,1)=0.0
  200      continue

           icn=0
           ic=ic+1

           if(lum(ic).ne.'(') goto 999
  211       ic=ic+1
            if(ic.gt.icolm) goto 999
            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 211
            if(dnen2(lum(ic))) goto 999

  220         ici=ic
  210         ic=ic+1
              if(ic.gt.icolm) goto 999
              if(lum(ic).eq.' '.or.lum(ic).eq.tub.or.deqn3(lum(ic)))
     &                        goto 210

               if(lum(ic).eq.','.or.lum(ic).eq.')') then
                  icf=ic-1
                  call rnum(rc,lum,ici,icf,ierr)
                  if(ierr.ne.0) goto 999
                  icn=icn+1
                  if(icn.gt.mc) goto 999
                  rcut(icn,1)=rc
                  if(lum(ic).eq.')') goto 300
                  if(lum(ic).eq.',') then
                     ic=ic+1
                     goto 220
                  end if
               else
                  goto 999
               end if
  300     icut=icn
          ic = ic+1

*----------------------------------------------------------------------*

        else if(dum4.eq.'cllg') then
             call dcols(1,lum,ic,cllg,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'clpa') then
             call dcols(1,lum,ic,clpa,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'clcn') then
             call dcols(1,lum,ic,clcn,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'lbcb') then
             call dcols(1,lum,ic,clgb,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'lbcl') then
             call dcols(1,lum,ic,clgl,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'lbcs') then
             call dcols(1,lum,ic,clgs,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'pbcb') then
             call dcols(1,lum,ic,clpb,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'pbcl') then
             call dcols(1,lum,ic,clpl,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'pbcs') then
             call dcols(1,lum,ic,clps,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'cbcb') then
             call dcols(1,lum,ic,clnb,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'cbcl') then
             call dcols(1,lum,ic,clnl,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'cbcs') then
             call dcols(1,lum,ic,clns,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

*----------------------------------------------------------------------*

        else if(dum4.eq.'lbox') then

            ic = ic + 1

               if( lum(ic) .ne. '(' ) goto 999

  600       ic = ic + 1

            if( ic .gt. icolm ) goto 999

            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 600

                  do 601 k = 1, 13

                     rum(k:k) = dum(ic+k-1)

  601             continue

               if( rum(1:9) .eq. 'singlebox' ) then

                     ilbox = 10
                     iboxl = 9

               else if( rum(1:9) .eq. 'Singlebox' ) then

                     ilbox = 11
                     iboxl = 9

               else if( rum(1:9) .eq. 'singleBox' ) then

                     ilbox = 12
                     iboxl = 9

               else if( rum(1:9) .eq. 'SingleBox' ) then

                     ilbox = 13
                     iboxl = 9

               else if( rum(1:7) .eq. 'ovalbox' ) then

                     ilbox = 20
                     iboxl = 7

               else if( rum(1:7) .eq. 'Ovalbox' ) then

                     ilbox = 21
                     iboxl = 7

               else if( rum(1:7) .eq. 'ovalBox' ) then

                     ilbox = 22
                     iboxl = 7

               else if( rum(1:7) .eq. 'OvalBox' ) then

                     ilbox = 23
                     iboxl = 7

               else if( rum(1:9) .eq. 'doublebox' ) then

                     ilbox = 30
                     iboxl = 9

               else if( rum(1:9) .eq. 'Doublebox' ) then

                     ilbox = 31
                     iboxl = 9

               else if( rum(1:9) .eq. 'doubleBox' ) then

                     ilbox = 32
                     iboxl = 9

               else if( rum(1:9) .eq. 'DoubleBox' ) then

                     ilbox = 33
                     iboxl = 9

               else if( rum(1:9) .eq. 'shadowbox' ) then

                     ilbox = 40
                     iboxl = 9

               else if( rum(1:9) .eq. 'Shadowbox' ) then

                     ilbox = 41
                     iboxl = 9

               else if( rum(1:9) .eq. 'shadowBox' ) then

                     ilbox = 42
                     iboxl = 9

               else if( rum(1:9) .eq. 'ShadowBox' ) then

                     ilbox = 43
                     iboxl = 9

               else if( rum(1:13) .eq. 'ovalshadowbox' ) then

                     ilbox = 50
                     iboxl = 13

               else if( rum(1:13) .eq. 'Ovalshadowbox' ) then

                     ilbox = 51
                     iboxl = 13

               else if( rum(1:13) .eq. 'ovalshadowBox' ) then

                     ilbox = 52
                     iboxl = 13

               else if( rum(1:13) .eq. 'OvalshadowBox' ) then

                     ilbox = 53
                     iboxl = 13

               else

                     goto 999

               end if

            ic = ic + iboxl - 1

  602       ic = ic + 1

            if( ic .gt. icolm ) goto 999

            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 602

            if( lum(ic) .ne. ')' ) goto 999

            ic = ic + 1

*----------------------------------------------------------------------*

        else if(dum4.eq.'pbox') then

            ic = ic + 1

               if( lum(ic) .ne. '(' ) goto 999

  610       ic = ic + 1

            if( ic .gt. icolm ) goto 999

            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 610

                  do 611 k = 1, 13

                     rum(k:k) = dum(ic+k-1)

  611             continue

               if( rum(1:9) .eq. 'singlebox' ) then

                     ipbox = 10
                     iboxl = 9

               else if( rum(1:9) .eq. 'Singlebox' ) then

                     ibox = 11
                     iboxl = 9

               else if( rum(1:9) .eq. 'singleBox' ) then

                     ipbox = 12
                     iboxl = 9

               else if( rum(1:9) .eq. 'SingleBox' ) then

                     ipbox = 13
                     iboxl = 9

               else if( rum(1:7) .eq. 'ovalbox' ) then

                     ipbox = 20
                     iboxl = 7

               else if( rum(1:7) .eq. 'Ovalbox' ) then

                     ipbox = 21
                     iboxl = 7

               else if( rum(1:7) .eq. 'ovalBox' ) then

                     ipbox = 22
                     iboxl = 7

               else if( rum(1:7) .eq. 'OvalBox' ) then

                     ipbox = 23
                     iboxl = 7

               else if( rum(1:9) .eq. 'doublebox' ) then

                     ipbox = 30
                     iboxl = 9

               else if( rum(1:9) .eq. 'Doublebox' ) then

                     ipbox = 31
                     iboxl = 9

               else if( rum(1:9) .eq. 'doubleBox' ) then

                     ipbox = 32
                     iboxl = 9

               else if( rum(1:9) .eq. 'DoubleBox' ) then

                     ipbox = 33
                     iboxl = 9

               else if( rum(1:9) .eq. 'shadowbox' ) then

                     ipbox = 40
                     iboxl = 9

               else if( rum(1:9) .eq. 'Shadowbox' ) then

                     ipbox = 41
                     iboxl = 9

               else if( rum(1:9) .eq. 'shadowBox' ) then

                     ipbox = 42
                     iboxl = 9

               else if( rum(1:9) .eq. 'ShadowBox' ) then

                     ipbox = 43
                     iboxl = 9

               else if( rum(1:13) .eq. 'ovalshadowbox' ) then

                     ipbox = 50
                     iboxl = 13

               else if( rum(1:13) .eq. 'Ovalshadowbox' ) then

                     ipbox = 51
                     iboxl = 13

               else if( rum(1:13) .eq. 'ovalshadowBox' ) then

                     ipbox = 52
                     iboxl = 13

               else if( rum(1:13) .eq. 'OvalshadowBox' ) then

                     ipbox = 53
                     iboxl = 13

               else

                     goto 999

               end if

            ic = ic + iboxl - 1

  612       ic = ic + 1

            if( ic .gt. icolm ) goto 999

            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 612

            if( lum(ic) .ne. ')' ) goto 999

            ic = ic + 1

*----------------------------------------------------------------------*

        else if(dum4.eq.'cbox') then

            ic = ic + 1

               if( lum(ic) .ne. '(' ) goto 999

  620       ic = ic + 1

            if( ic .gt. icolm ) goto 999

            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 620

                  do 621 k = 1, 13

                     rum(k:k) = dum(ic+k-1)

  621             continue

               if( rum(1:9) .eq. 'singlebox' ) then

                     isbox = 10
                     iboxl = 9

               else if( rum(1:9) .eq. 'Singlebox' ) then

                     isbox = 11
                     iboxl = 9

               else if( rum(1:9) .eq. 'singleBox' ) then

                     isbox = 12
                     iboxl = 9

               else if( rum(1:9) .eq. 'SingleBox' ) then

                     isbox = 13
                     iboxl = 9

               else if( rum(1:7) .eq. 'ovalbox' ) then

                     isbox = 20
                     iboxl = 7

               else if( rum(1:7) .eq. 'Ovalbox' ) then

                     isbox = 21
                     iboxl = 7

               else if( rum(1:7) .eq. 'ovalBox' ) then

                     isbox = 22
                     iboxl = 7

               else if( rum(1:7) .eq. 'OvalBox' ) then

                     isbox = 23
                     iboxl = 7

               else if( rum(1:9) .eq. 'doublebox' ) then

                     isbox = 30
                     iboxl = 9

               else if( rum(1:9) .eq. 'Doublebox' ) then

                     isbox = 31
                     iboxl = 9

               else if( rum(1:9) .eq. 'doubleBox' ) then

                     isbox = 32
                     iboxl = 9

               else if( rum(1:9) .eq. 'DoubleBox' ) then

                     isbox = 33
                     iboxl = 9

               else if( rum(1:9) .eq. 'shadowbox' ) then

                     isbox = 40
                     iboxl = 9

               else if( rum(1:9) .eq. 'Shadowbox' ) then

                     isbox = 41
                     iboxl = 9

               else if( rum(1:9) .eq. 'shadowBox' ) then

                     isbox = 42
                     iboxl = 9

               else if( rum(1:9) .eq. 'ShadowBox' ) then

                     isbox = 43
                     iboxl = 9

               else if( rum(1:13) .eq. 'ovalshadowbox' ) then

                     isbox = 50
                     iboxl = 13

               else if( rum(1:13) .eq. 'Ovalshadowbox' ) then

                     isbox = 51
                     iboxl = 13

               else if( rum(1:13) .eq. 'ovalshadowBox' ) then

                     isbox = 52
                     iboxl = 13

               else if( rum(1:13) .eq. 'OvalshadowBox' ) then

                     isbox = 53
                     iboxl = 13

               else

                     goto 999

               end if

            ic = ic + iboxl - 1

  622       ic = ic + 1

            if( ic .gt. icolm ) goto 999

            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 622

            if( lum(ic) .ne. ')' ) goto 999

            ic = ic + 1

*----------------------------------------------------------------------*

        else if(dum4.eq.'clal') then
             call dcols(1,lum,ic,clal,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'clax') then
             call dcols(1,lum,ic,clax,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'cltx') then
             call dcols(1,lum,ic,cltx,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'clnm') then
             call dcols(1,lum,ic,clnm,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'cltl') then
             call dcols(1,lum,ic,cltl,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'clfr') then
             call dcols(1,lum,ic,clfr,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'clms') then
             call dcols(1,lum,ic,clms,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'clhd') then
             call dcols(1,lum,ic,clhd,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'clbg') then
             call dcols(1,lum,ic,clbg,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

        else if(dum4.eq.'clin') then
             call dcols(1,lum,ic,clin,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

*----------------------------------------------------------------------*

        else if(dum4.eq.'smax') then
             call dcols(0,lum,ic,rcol,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

               smax(1) = rcol(1)
               smax(2) = rcol(2)
               smax(3) = rcol(3)

        else if(dum4.eq.'smin') then
             call dcols(0,lum,ic,rcol,ierr,icolm,c1,c2)
             if(ierr.eq.1) goto 999

               smin(1) = rcol(1)
               smin(2) = rcol(2)
               smin(3) = rcol(3)

*----------------------------------------------------------------------*

        else if(dum4.eq.'cols') then

           do 400 i = 1, mc

             rcut(i,2) = -2.0
             rcut(i,3) =  1.0
             rcut(i,4) =  1.0

  400      continue

            icol = 0
            icn  = 0
            ic = ic + 1

          if(lum(ic).eq.' '.or.lum(ic).eq.tub) then

             icct = -1

          else

            if(lum(ic).ne.'(') goto 999

  411       ic = ic + 1

            if( ic .gt. icolm ) goto 999
            if( icol .eq. 0 .and.
     &        ( lum(ic) .eq. ' ' .or. lum(ic) .eq. tub ) ) goto 411

               if( icol .eq. 0 ) ici = ic - 2

            if( lum(ic) .eq. ',' .or. lum(ic) .eq. ')' ) then

                  if( icol .eq. 0 ) goto 999

               call dcols(1,lum,ici,rcol,ierr,icolm,lum(ici+1),lum(ic))

                  if( ierr .ne. 0 ) goto 999

                      icn  = icn + 1
                      icol = 0

                      if( icn. gt. mc ) goto 999

                      rcut(icn,2) = rcol(1)
                      rcut(icn,3) = rcol(2)
                      rcut(icn,4) = rcol(3)

                  if( lum(ic) .eq. ')' ) goto 430

            else

                  icol = icol + 1

            end if

               goto 411


  430       icct = icn
            ic   = ic + 1

          end if

*----------------------------------------------------------------------*

        else if(dum4.eq.'conl') then

           do i = 1, mc
             rcut(i,5) =  0.0
           end do

            icol = 0
            icn  = 0
            ic = ic + 1

            if(lum(ic).ne.'(') goto 999

  441       ic = ic + 1

            if( ic .gt. icolm ) goto 999
            if( lum(ic) .eq. ' ' .or. lum(ic) .eq. tub ) goto 441

            if( lum(ic) .eq. ',' .or. lum(ic) .eq. ')' ) then

                  if( icol .eq. 0 ) goto 999
                      icol = 0

                  if( lum(ic) .eq. ')' ) goto 450

            else

                  icol = icol + 1

                  icn  = icn  + 1
                  if( icn. gt. mc ) goto 999

               if( lum(ic) .eq. 'l' ) then
                  rlin =  0.0
               else if( lum(ic) .eq. 'm' ) then
                  rlin =  1.0
               else if( lum(ic) .eq. 'd' ) then
                  rlin =  2.0
               else if( lum(ic) .eq. 'u' ) then
                  rlin =  3.0
               else if( lum(ic) .eq. 'p' ) then
                  rlin =  4.0
               else if( lum(ic) .eq. 'q' ) then
                  rlin =  5.0
               else if( lum(ic) .eq. 'v' ) then
                  rlin =  6.0
               else

                  goto 999

               end if

                  do icl = icn, mc
                     rcut(icl,5) =  rlin
                  end do

            end if

               goto 441

  450       ic = ic + 1

*----------------------------------------------------------------------*

        else if(dum4.eq.'iwd2') then

           do i = 1, mc
             iwd2(i) =  4
           end do

           icn=0
           ic=ic+1

           if(lum(ic).ne.'(') goto 999
  711       ic=ic+1
            if(ic.gt.icolm) goto 999
            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 711
            if(dnen2(lum(ic))) goto 999

  720         ici=ic
  710         ic=ic+1
              if(ic.gt.icolm) goto 999
              if(lum(ic).eq.' '.or.lum(ic).eq.tub.or.deqn3(lum(ic)))
     &                        goto 710

               if(lum(ic).eq.','.or.lum(ic).eq.')') then
                  icf=ic-1
                  call rnum(rc,lum,ici,icf,ierr)
                  if(ierr.ne.0) goto 999
                  icn=icn+1
                  if(icn.gt.mc) goto 999
                  iwd2(icn)=nint(rc)
                  if(lum(ic).eq.')') goto 700
                  if(lum(ic).eq.',') then
                     ic=ic+1
                     goto 720
                  end if
               else
                  goto 999
               end if
  700     ic = ic+1

          do icl = icn, mc
             iwd2(icl) = iwd2(icn)
          end do

*-----------------------------------------------------------------------

        else if(dum4.eq.'cmap') then
            itmp=0
  751       ic=ic+1
            if( lum(ic) .eq. ' ' ) goto 751
            if( lum(ic) .eq. tub ) goto 751
            if( lum(ic) .eq. '(' ) goto 751
            if( lum(ic) .eq. '"' ) goto 751
            if( lum(ic) .eq. "'" ) goto 751
            if( lum(ic) .eq. ")" ) goto 752 ! exit
            itmp=itmp+1
            cmap(itmp:itmp)=lum(ic)
            goto 751
  752       if(itmp>1) then
               if(cmap(itmp-1:itmp)=="_r") then
                  icrev=1 ! reverse colormap
                  cmap=cmap(1:itmp-2)
               endif
            endif
            ic=ic+1

        else if(dum4.eq.'ndis') then !
              ic=ic+1
           call pnum(lum,ic,icolm,xndis,ierr)
              if(ierr.ne.0) goto 999
              ndis=nint(xndis)
              ic=ic+1

*----------------------------------------------------------------------*

        else

            goto 999

        end if

*-----------------------------------------------------------------------

      if(lum(ic).ne.' '.and.lum(ic).ne.tub) goto 999

      goto 100

  500 continue

      return

*-----------------------------------------------------------------------

  999 ierr = 1

      m_err = 'P: Parameter Description is Wrong in [ '//dum4//' ]'
      ErrCha = ''
      ErrID = 'L:1933/R:para/F:a-main1.f'

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine prehtwo(lum,in,jsn,ill,
     &                   xfin,yfin,xint,yint,ixy,ic,
     &                   ierr)
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /h2/ xdd,ydd,
     &            xymax,xymin,clus,cmax,cmin,cxymax,cxymin,xyabs,
     &            dmax,dmin,ixnm,iynm

      character lum(ichrl)*1

      common /error/ m_err, l_err, k_err
      character m_err*200

      dimension ill(0:9)

      logical dnen2,deqn3

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

      ixy=0

      ic=in-1

      do
  100 ic=ic+1
        if(ic.gt.icolm) goto 500
        if(lum(ic).eq.' '.or.lum(ic).eq.tub.or.lum(ic).eq.';'.or.
     &     lum(ic).eq.':'.or.lum(ic).eq.',') goto 100
        if(lum(ic).ne.'y'.and.lum(ic).ne.'x') goto 998

      if(lum(ic).eq.'x') then
         if(ixy.eq.0) ixy=-1

  200    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(lum(ic).eq.' '.or.lum(ic).eq.tub.or.lum(ic).eq.'=')
     &                      goto 200
            if(dnen2(lum(ic))) goto 998
         ici=ic

  210    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(deqn3(lum(ic))) goto 210
            if(lum(ic).ne.' '.and.lum(ic).ne.tub) goto 998
            call rnum(xint,lum,ici,ic-1,ierr)
            if(ierr.ne.0) goto 999

  220    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 220
            if(lum(ic)//lum(ic+1).ne.'to') goto 998
            ic=ic+1

  230    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 230
            if(dnen2(lum(ic))) goto 998
         ici=ic

  240    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(deqn3(lum(ic))) goto 240
            if(lum(ic).ne.' '.and.lum(ic).ne.tub) goto 998
            call rnum(xfin,lum,ici,ic-1,ierr)
            if(ierr.ne.0) goto 999

  250    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 250
            if(lum(ic)//lum(ic+1).ne.'by') goto 998
            ic=ic+1

  260    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 260
            if(dnen2(lum(ic))) goto 998
         ici=ic

  270    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(deqn3(lum(ic))) goto 270
            if(lum(ic).ne.' '.and.lum(ic).ne.tub
     &               .and.lum(ic).ne.':'.and.
     &         lum(ic).ne.';'.and.lum(ic).ne.',') goto 998
            call rnum(xdd,lum,ici,ic-1,ierr)
            if(ierr.ne.0) goto 999
            if(xdd.le.0.0) xdd=-xdd
      end if


      if(lum(ic).eq.'y') then
         if(ixy.eq.0) ixy=1

  300    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(lum(ic).eq.' '.or.lum(ic).eq.tub.or.lum(ic).eq.'=')
     &                      goto 300
            if(dnen2(lum(ic))) goto 998
         ici=ic

  310    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(deqn3(lum(ic))) goto 310
            if(lum(ic).ne.' '.and.lum(ic).ne.tub) goto 998
            call rnum(yint,lum,ici,ic-1,ierr)
            if(ierr.ne.0) goto 999

  320    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 320
            if(lum(ic)//lum(ic+1).ne.'to') goto 998
            ic=ic+1

  330    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 330
            if(dnen2(lum(ic))) goto 998
         ici=ic

  340    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(deqn3(lum(ic))) goto 340
            if(lum(ic).ne.' '.and.lum(ic).ne.tub) goto 998
            call rnum(yfin,lum,ici,ic-1,ierr)
            if(ierr.ne.0) goto 999

  350    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 350
            if(lum(ic)//lum(ic+1).ne.'by') goto 998
            ic=ic+1

  360    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 360
            if(dnen2(lum(ic))) goto 998
         ici=ic

  370    ic=ic+1
            if(ic.gt.icolm) goto 999
            if(deqn3(lum(ic))) goto 370
            if(lum(ic).ne.' '.and.lum(ic).ne.tub
     &               .and.lum(ic).ne.':'.and.
     &         lum(ic).ne.';'.and.lum(ic).ne.',') goto 998
            call rnum(ydd,lum,ici,ic-1,ierr)
            if(ierr.ne.0) goto 999
            if(ydd.le.0.0) ydd=-ydd
      end if

      enddo
  500 continue

*-----------------------------------------------------------------------
*       IXNM : NUMBER OF X-POINTS
*       IYNM : NUMBER OF Y-POINTS
*-----------------------------------------------------------------------

      ixnm = nint(abs(xfin-xint)/xdd)+1
      iynm = nint(abs(yfin-yint)/ydd)+1

      return

  998 m_err = 'H2: HC: HD: X, Y Descriptions are Wrong.'
      ErrCha = ''
      ErrID = 'L:2122/R:prehtwo/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn

  999 ierr  = 1

      return
      end

************************************************************************
*                                                                      *
      subroutine htwo(dum,lum,in,jsi,jsn,ill,ilf,dsin,idsi,
     &                ymax,ymin,xmax,xmin,
     &                nhl,nh2,ierr,mhc2,ipds,ipdc,clhd,
     &                izlog,daxy,dax,day,
     &                xfin,yfin,xint,yint,ixy,ic,
     &                ih2fs)
*                                                                      *
*        PURPOSE:  FOR 3D CLUSTER PLOT                                 *
*                                                                      *
*         NHC:    NUMBER OF COLORS  POINTS NOT EQUAL TO ZERO           *
*                                                                      *
*         NHL:    NUMBER OF CLUSTER POINTS NOT EQUAL TO ZERO           *
*                                                                      *
*         NH2:    FOR CONTOUR PLOT                                     *
*                                                                      *
*         MHC2:   MHC2=1 : CLUSTER PLOT, MHC2=2 : CONTOUR PLOT         *
*                 MHC2=3 : COLOR PLOT                                  *
*                                                                      *
*         IPDC:   DUPLICATE OF COLOR AND MONOCLE CLUSTERS              *
*         IPDs:   smoothing for contour plot                           *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /h2/ xdd,ydd,
     &            xymax,xymin,clus,cmax,cmin,cxymax,cxymin,xyabs,
     &            dmax,dmin,ixnm,iynm

      common /h3/ smax(3), smin(3)
      common /h5/ nhc, mhc(100)
      common /h6/ dhgh(100), dwih(100)

      dimension daxy(ixnm*iynm)
      dimension dax(ixnm), day(iynm)

      character dum(ichrl)*1
      character lum(ichrl)*1

      common /error/ m_err, l_err, k_err
      character m_err*200

      dimension ill(0:9), ilf(0:9)
      character dsin(0:9)*200
      dimension idsi(0:9)

      logical dnen2,deqn3

      dimension dfct(-50:50,-50:50)
      dimension clhd(3)
      dimension rcolp(3), rcoln(3), rcol(3)

      common /cmap/ ndis, icrev, cmap
      character cmap*99

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------
*     set memory
*-----------------------------------------------------------------------

         if(xfin.gt.xint) then
           xin=xint
           if(xfin.gt.xmax) xmax=xfin
           if(xint.lt.xmin) xmin=xint
         else
           xin=xfin
           if(xint.gt.xmax) xmax=xint
           if(xfin.lt.xmin) xmin=xfin
         end if

         do 400 i=1,ixnm
           dax(i)=xin+dble(i-1)*xdd
  400    continue

         if(yfin.gt.yint) then
           yin=yint
           if(yfin.gt.ymax) ymax=yfin
           if(yint.lt.ymin) ymin=yint
         else
           yin=yfin
           if(yint.gt.ymax) ymax=yint
           if(yfin.lt.ymin) ymin=yfin
         end if
         do 401 i=1,iynm
           day(i)=yin+dble(i-1)*ydd
  401    continue

*-----------------------------------------------------------------------
*     SEQUENTIAL CONTROL OF THE DATA
*-----------------------------------------------------------------------

               if( xfin .gt. xint ) then

                  ixi = 1
                  ixf = ixnm
                  ixd = 1

               else

                  ixf = 1
                  ixi = ixnm
                  ixd = -1

               end if

               if( yfin .gt. yint ) then

                  iyi = 1
                  iyf = iynm
                  iyd = 1

               else

                  iyf = 1
                  iyi = iynm
                  iyd = -1

               end if

                  iiin = ixi
                  iifn = ixf
                  iidd = ixd

                  jjin = iyi
                  jjfn = iyf
                  jjdd = iyd

               if( ixy .eq. 1 ) then

                  jjal = abs( iifn - iiin ) + 1

               else

                  jjal = abs( jjfn - jjin ) + 1

               end if

                  itot = ixnm * iynm

*-----------------------------------------------------------------------

      if( ih2fs .eq. 0 ) then

*-----------------------------------------------------------------------
*     CHECK OF TOTAL NUMBER
*-----------------------------------------------------------------------

         ila0 = 0
         ila  = 0
         idn  = 0

*-----------------------------------------------------------------------

  510    ila0 = ila0 + 1
         ila  = ila  + 1

         iln = 0

*-----------------------------------------------------------------------

         read(jsi,'(10000a1)', iostat = ios ) (dum(ic),ic=1,icolm)
         if( ios .eq. -1 ) goto 560

         call chlow(dum,lum)

            if( ill(jsn) + ila .gt. ilf(jsn) ) goto 560

*-----------------------------------------------------------------------
*        SKIP BLANCK AND BLANK LINE IS THE END OF THE DATA
*-----------------------------------------------------------------------

            k = 1
            do 145 l = 1, icolm
               if(lum(l).ne.' '.and.lum(l).ne.tub) goto 146
  145       continue

            goto 560

  146       k = l

*-----------------------------------------------------------------------
*     INCLFL: INCLUDE FILE ONLY FOR ILA0 = 1
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
*        END OF THE H2: SECTION
*-----------------------------------------------------------------------

         if( lum(k+1) .eq. ':' .or.
     &       lum(k+2) .eq. ':' .or.
     &       dnen2(lum(k)) ) goto 560

*-----------------------------------------------------------------------
*     CHECK OF NUMERICS
*-----------------------------------------------------------------------

         ic = k - 1

  530    ic = ic + 1

            if( ic .gt. icolm .and. iln .eq. 0 ) goto 560
            if( ic .gt. icolm ) goto 510
            if( lum(ic) .eq. ' ' .or. lum(ic) .eq. tub ) goto 530

            if( dnen2(lum(ic)) ) goto 997

*-----------------------------------------------------------------------
*        FIND ONE NUMBER AND ASSIGN THE SUFIXES
*-----------------------------------------------------------------------

            idn = idn + 1
            iln = iln + 1

            ici = ic

         if( ixy .eq. 1 ) then

            jorg = ( idn - 1 ) / jjal + 1
            iorg = idn - ( idn - 1 ) / jjal * jjal

         else

            iorg = ( idn - 1 ) / jjal + 1
            jorg = idn - ( idn - 1 ) / jjal * jjal

         end if

            iix  = iiin + ( iorg - 1 ) * iidd
            iiy  = jjin + ( jorg - 1 ) * jjdd

*-----------------------------------------------------------------------

  550    ic = ic + 1

            if( ic .gt. icolm .and. iln .eq. 0 ) goto 560

            if( ic .gt. icolm ) then

               call rnum(rrnm,lum,ici,ic-1,ierr)
               if(ierr.ne.0) goto 997

               daxy(iix+(iiy-1)*ixnm) = rrnm

               goto 510

            end if


            if(  lum(ic) .eq. ' ' .or. lum(ic) .eq. tub ) then

               call rnum(rrnm,lum,ici,ic-1,ierr)
               if(ierr.ne.0) goto 997

               daxy(iix+(iiy-1)*ixnm) = rrnm

               goto 530

            end if

            if( deqn3(lum(ic)) ) goto 550

            goto 997

*-----------------------------------------------------------------------
*        END OF CHECK OF DATA AND CHECK TOTAL NUMBERS OF DATA
*-----------------------------------------------------------------------

  560       if( idn .ne. itot ) then

               m_err = 'H2:, HD: or HC: Number of Data Mismatches '//
     &                 'with the Description.'
               ErrCha = ''
               ErrID = 'L:2446/R:htwo/F:a-main1.f'
               l_err = ill(jsn)
               k_err = jsn
               goto 999

            end if

*-----------------------------------------------------------------------

            ill(jsn) = ill(jsn) + ila

*-----------------------------------------------------------------------

      else if( ih2fs .eq. 1 ) then

*-----------------------------------------------------------------------

         if( ixy .eq. 1 ) then

            read(jsi,*, iostat = ios )
     &       ((daxy(iix+(iiy-1)*ixnm),iix=ixi,ixf,ixd),iiy=iyi,iyf,iyd)

         else

            read(jsi,*, iostat = ios )
     &       ((daxy(iiy+(iix-1)*iynm),iiy=iyi,iyf,iyd),iix=ixi,ixf,ixd)

         end if

            if( ios .eq. -1 ) then

               m_err = 'H2:, HD: or HC: Data is something wrong '
               ErrCha = ''
               ErrID = 'L:2479/R:htwo/F:a-main1.f'
               l_err = ill(jsn)
               k_err = jsn
               goto 999

            end if

*-----------------------------------------------------------------------
cKN for PHITS only ???

          if( itot .eq. 100 ) then

            ila = ( itot - 1 ) / 20 + 2

          else

            ila = ( itot - 1 ) / 10 + 2

          end if

*-----------------------------------------------------------------------

  670    continue


         read(jsi,'(10000a1)', iostat = ios ) (dum(ic),ic=1,icolm)
         if( ios .eq. -1 ) goto 660

         call chlow(dum,lum)

  660    continue

         ill(jsn) = ill(jsn) + ila

      end if

*-----------------------------------------------------------------------
*     SUMMARY OF THE DATA
*-----------------------------------------------------------------------
cKN 2020/10/24 43,41,42  -> 83,81,82

      if( izlog .eq. 1 ) then

           xymax = -83.0
           xymin =  81.0

         do 705 i = 1, ixnm
         do 715 j = 1, iynm

         if(daxy(i+(j-1)*ixnm).lt.1.d-82) daxy(i+(j-1)*ixnm)=1.d-82 ! T.Sato 2020/06/19, to avoid error when negative value


              daxy(i+(j-1)*ixnm) = log10(daxy(i+(j-1)*ixnm))

              if(daxy(i+(j-1)*ixnm).gt.xymax)
     &           xymax = daxy(i+(j-1)*ixnm)

              if(daxy(i+(j-1)*ixnm).lt.xymin.and.
     &           daxy(i+(j-1)*ixnm).gt.-41.0)
     &           xymin=daxy(i+(j-1)*ixnm)


  715    continue
  705    continue

            cxymax = xymax
            cxymin = xymin

            if(cmax.gt.-r0max) cxymax = log10(cmax)
            if(cmin.lt.+r0max) cxymin = log10(cmin)

            dxymax = xymax
            dxymin = xymin

            if(dmax.gt.-r0max) dxymax = log10(dmax)
            if(dmin.lt.+r0max) dxymin = log10(dmin)

      else

           xymax = - r1max
           xymin =   r1max

         do 700 i=1,ixnm
         do 710 j=1,iynm

            if(daxy(i+(j-1)*ixnm).gt.xymax) xymax=daxy(i+(j-1)*ixnm)
            if(daxy(i+(j-1)*ixnm).lt.xymin) xymin=daxy(i+(j-1)*ixnm)

  710    continue
  700    continue

            cxymax = xymax
            cxymin = xymin

            if(cmax.gt.-r0max) cxymax = cmax
            if(cmin.lt.+r0max) cxymin = cmin

            dxymax = xymax
            dxymin = xymin

            if(dmax.gt.-r0max) dxymax = dmax
            if(dmin.lt.+r0max) dxymin = dmin

      end if

*-----------------------------------------------------------------------
*     WRITE DATA TO JHL FOR CLUSTER PLOT
*-----------------------------------------------------------------------

      if( mhc2 .eq. 1 .or. mhc2 .eq. 4 )  then

                     xyabs = dxymax - dxymin

               if( abs( xyabs ) .lt. 1.0d-41 ) then

                  if( abs( dxymax ) .lt. 1.0d-41 ) then

                     xyabs =  1.0
                     ddmin = -1.0

                  else

                     xyabs = abs( dxymax )
                     ddmin = 0.0

                  end if

               else

                  if( dxymax * dxymin .lt. 0.0 ) then

                     xyabs = max( abs(dxymax), abs(dxymin) )
                     ddmin = 0.0

                  end if

               end if

                  rcolp(1)  = clhd(1)
                  rcolp(2)  = clhd(2)
                  rcolp(3)  = clhd(3)
                  rcoln(1)  = clhd(1)
                  rcoln(2)  = clhd(2)
                  rcoln(3)  = clhd(3)

         if( clhd(1) .le. 0.0 ) then

                  if( clhd(1) .lt. -r0max ) clhd(1) = -2.0

                  rcolp(1) = clhd(1)

               if( clhd(1) .ge. -1.6 ) then

                  rcoln(1) = -2.0

               else

                  rcoln(1) = -1.4

               end if

         else

            if( rcolp(1) .le. 1.5 ) then

               rcoln(1) = 1.5 + rcolp(1) - 1.0

            else if( rcolp(1) .gt. 1.5 ) then

               rcoln(1) = 0.5 + rcolp(1) - 1.0

            end if

         end if


      do 720 i=1,ixnm
      do 730 j=1,iynm

        if(daxy(i+(j-1)*ixnm).ne.0.0) then
        if(daxy(i+(j-1)*ixnm).ge.dxymin.and.
     &     daxy(i+(j-1)*ixnm).le.dxymax) then

          dazz = max(dble(daxy(i+(j-1)*ixnm)),cxymin)
          dazz = min(dazz,cxymax)

          if( dazz .ge. 0.0 ) then
                 rcol(1) = rcolp(1)
                 rcol(2) = rcolp(2)
                 rcol(3) = rcolp(3)
          else
                 rcol(1) = rcoln(1)
                 rcol(2) = rcoln(2)
                 rcol(3) = rcoln(3)
          end if

          if( izlog .eq. 1 ) then
                 rate = sqrt( abs(dazz - ddmin) / xyabs )
                 rcol(1) = rcolp(1)
                 rcol(2) = rcolp(2)
                 rcol(3) = rcolp(3)
          else
                 rate = sqrt( abs(dazz - ddmin) / xyabs )
          end if

             write(jhl)
     &         dax(i),day(j),xdd*rate,ydd*rate,rcol

             nhl=nhl+1

        end if
        end if

  730 continue
  720 continue

      end if

*-----------------------------------------------------------------------
*     WRITE DATA TO JHC FOR COLOR PLOT
*-----------------------------------------------------------------------

      if( mhc2 .eq. 3 .or. mhc2 .eq. 5 )  then

         nhc = nhc + 1
         if( nhc .gt. 100 ) goto 995

                     xyabs = cxymax - cxymin

               if( abs( xyabs ) .lt. 1.0d-41 ) then

                  if( abs( cxymax ) .lt. 1.0d-41 ) then

                     xyabs  =  1.0
                     cxymin = -1.0

                  else

                     xyabs = abs( cxymax ) / 2.0
                     cxymin = cxymax - xyabs

                  end if

               end if

                  if(ipdc.eq.0) then

                     xdp = xdd * 1.03
                     ydp = ydd * 1.03

                  else

                     xdp = xdd / 2.0 * 1.03
                     ydp = ydd / 2.0 * 1.03

                  end if

                  if( smax(1) * smin(1) .le. 0.0d0 ) then

                     smin(1) = 3.0
                     smax(1) = 1.0

                  end if

         write(jhc) xdp, ydp

      do 721 i = 1, ixnm

             dax0 = dax(i)
             dax1 = dax0

      do 731 j = 1, iynm

       if(daxy(i+(j-1)*ixnm).ge.dxymin.and.
     &    daxy(i+(j-1)*ixnm).le.dxymax ) then

             dhi0 = max(dble(daxy(i+(j-1)*ixnm)),cxymin)
             dhi0 = min(dhi0,cxymax)
             day0 = day(j)
             rcl1 = smin(1) + ( dhi0 - cxymin )
     &            / xyabs * ( smax(1) - smin(1) )
             rcl2 = smin(2) + ( dhi0 - cxymin )
     &            / xyabs * ( smax(2) - smin(2) )
             rcl3 = smin(3) + ( dhi0 - cxymin )
     &            / xyabs * ( smax(3) - smin(3) )

             rcol(1) = chue(rcl1)
             rcol(2) = rcl2
             rcol(3) = rcl3

             call applyColorMap(cmap,ndis,icrev,rcol,rcl1)

             write(jhc) dax0,day0,rcol

             mhc(nhc) = mhc(nhc) + 1

        if(ipdc.eq.1) then

          if(i.ne.ixnm.and.j.ne.iynm) then

             dhi1 = ( dhi0 + daxy(i+(j+1-1)*ixnm) ) / 2.0
             dhi2 = ( dhi0 + daxy(i+1+(j-1)*ixnm) ) / 2.0
             dhi3 = ( dhi0 + daxy(i+(j+1-1)*ixnm)
     &                     + daxy(i+1+(j-1)*ixnm)
     &                     + daxy(i+1+(j+1-1)*ixnm) ) / 4.0

             dhi1 = max(dhi1,cxymin)
             dhi1 = min(dhi1,cxymax)
             dhi2 = max(dhi2,cxymin)
             dhi2 = min(dhi2,cxymax)
             dhi3 = max(dhi3,cxymin)
             dhi3 = min(dhi3,cxymax)

             day1 = ( day0 + day(j+1) ) / 2.0
             day2 =   day0
             day3 =   day1

             dax2 = ( dax0 + dax(i+1) ) / 2.0
             dax3 =   dax2

             rcl1 = smin(1) + ( dhi1 - cxymin )
     &            / xyabs * ( smax(1) - smin(1) )
             rcl2 = smin(2) + ( dhi1 - cxymin )
     &            / xyabs * ( smax(2) - smin(2) )
             rcl3 = smin(3) + ( dhi1 - cxymin )
     &            / xyabs * ( smax(3) - smin(3) )

             rcol(1) = chue(rcl1)
             rcol(2) = rcl2
             rcol(3) = rcl3

             write(jhc) dax1,day1,rcol

             mhc(nhc) = mhc(nhc) + 1

             rcl1 = smin(1) + ( dhi2 - cxymin )
     &            / xyabs * ( smax(1) - smin(1) )
             rcl2 = smin(2) + ( dhi2 - cxymin )
     &            / xyabs * ( smax(2) - smin(2) )
             rcl3 = smin(3) + ( dhi2 - cxymin )
     &            / xyabs * ( smax(3) - smin(3) )

             rcol(1) = chue(rcl1)
             rcol(2) = rcl2
             rcol(3) = rcl3

             write(jhc) dax2,day2,rcol

             mhc(nhc) = mhc(nhc) + 1

             rcl1 = smin(1) + ( dhi3 - cxymin )
     &            / xyabs * ( smax(1) - smin(1) )
             rcl2 = smin(2) + ( dhi3 - cxymin )
     &            / xyabs * ( smax(2) - smin(2) )
             rcl3 = smin(3) + ( dhi3 - cxymin )
     &            / xyabs * ( smax(3) - smin(3) )

             rcol(1) = chue(rcl1)
             rcol(2) = rcl2
             rcol(3) = rcl3

             write(jhc) dax3,day3,rcol

             mhc(nhc) = mhc(nhc) + 1

          else if(i.ne.ixnm.and.j.eq.iynm) then

             dhi2 = ( dhi0 + daxy(i+1+(j-1)*ixnm) ) / 2.0
             dhi2 = max(dhi2,cxymin)
             dhi2 = min(dhi2,cxymax)

             day2 =   day0
             dax2 = ( dax0 + dax(i+1) ) / 2.0

             rcl1 = smin(1) + ( dhi2 - cxymin )
     &            / xyabs * ( smax(1) - smin(1) )
             rcl2 = smin(2) + ( dhi2 - cxymin )
     &            / xyabs * ( smax(2) - smin(2) )
             rcl3 = smin(3) + ( dhi2 - cxymin )
     &            / xyabs * ( smax(3) - smin(3) )

             rcol(1) = chue(rcl1)
             rcol(2) = rcl2
             rcol(3) = rcl3

             write(jhc) dax2,day2,rcol

             mhc(nhc) = mhc(nhc) + 1

          else if(i.eq.ixnm.and.j.ne.iynm) then

             dhi1 = ( dhi0 + daxy(i+(j+1-1)*ixnm) ) / 2.0
             dhi1 = max(dhi1,cxymin)
             dhi1 = min(dhi1,cxymax)

             day1 = ( day0 + day(j+1) ) / 2.0

             rcl1 = smin(1) + ( dhi1 - cxymin )
     &            / xyabs * ( smax(1) - smin(1) )
             rcl2 = smin(2) + ( dhi1 - cxymin )
     &            / xyabs * ( smax(2) - smin(2) )
             rcl3 = smin(3) + ( dhi1 - cxymin )
     &            / xyabs * ( smax(3) - smin(3) )

             rcol(1) = chue(rcl1)
             rcol(2) = rcl2
             rcol(3) = rcl3

             write(jhc) dax1,day1,rcol

             mhc(nhc) = mhc(nhc) + 1

          end if

        end if

       end if

  731 continue
  721 continue

      end if

*-----------------------------------------------------------------------
*     FOR CONTOUR PLOT (ONLY COUNT THE CASE)
*-----------------------------------------------------------------------

      if( mhc2 .eq. 2 .or. mhc2 .eq. 4 .or. mhc2 .eq. 5 ) then

            if( mhc2 .eq. 2 .or. mhc2 .eq. 4 ) then

                     xyabs = cxymax - cxymin

               if( abs( xyabs ) .lt. 1.0d-41 ) then

                  if( abs( cxymax ) .lt. 1.0d-41 ) then

                     xyabs  =  1.0
                     cxymin = -1.0

                  else

                     xyabs = abs( cxymax ) / 2.0
                     cxymin = cxymax - xyabs

                  end if

               end if

            end if

*-----------------------------------------------------------------------

            nh2 = nh2 + 1

*-----------------------------------------------------------------------

         if( ipds .gt. 0 ) then

               open(jhs,form='unformatted',status='scratch')

               if( ipds .gt. 50 ) ipds = 50

               wid = sqrt( dble(ipds) * 0.5 )

            do i=-ipds, ipds
            do j=-ipds, ipds

               dfct(i,j) = exp(-(dble(i)**2+dble(j)**2)/(2.0*wid**2))

            end do
            end do

            do i=1,ixnm
            do j=1,iynm

                  sekf = 0.0
                  sekv = 0.0

               do ii = -ipds, ipds
               do jj = -ipds, ipds

                  k = i + ii
                  m = j + jj

                  if( k .ge. 1 .and. k .le. ixnm .and.
     &                m .ge. 1 .and. m .le. iynm ) then

                     sekf = sekf + dfct(ii,jj)
                     sekv = sekv + daxy(k+(m-1)*ixnm) * dfct(ii,jj)

                  end if

               end do
               end do

                  dbxy = sekv / sekf
                  write(jhs) sekv / sekf

            end do
            end do

                  rewind(jhs)

            do i=1,ixnm
            do j=1,iynm

                  read(jhs) dbxy
                  daxy(i+(j-1)*ixnm) = dbxy

            end do
            end do

            close(jhs)

         end if

*-----------------------------------------------------------------------

         return

      end if

*-----------------------------------------------------------------------

      return

  995 m_err = 'HC: section is exceed 100.'
      ErrCha = ''
      ErrID = 'L:3008/R:htwo/F:a-main1.f'
      l_err = ill(jsn)+ila
      k_err = jsn

      goto 999

  996 m_err = 'H2: HC: HD: memory is larger than mdas in ANGEL'
      ErrCha = ''
      ErrID = 'L:3016/R:htwo/F:a-main1.f'
      l_err = ill(jsn)+ila
      k_err = jsn

      goto 999

  997 m_err = 'H2: HC: HD: Numerical Data is Something Wrong.'
      ErrCha = ''
      ErrID = 'L:3024/R:htwo/F:a-main1.f'
      l_err = ill(jsn)+ila
      k_err = jsn

      goto 999

  998 m_err = 'H2: HC: HD: X, Y Descriptions are Wrong.'
      ErrCha = ''
      ErrID = 'L:3032/R:htwo/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn

  999 ierr  = 1

      return
      end


************************************************************************
*                                                                      *
      subroutine bmap(dum,lum,in,jsi,jsn,ill,ilf,dsin,idsi,
     &                ibmap,ymax,ymin,xmax,xmin,ierr)
*                                                                      *
*        purpose:  bitmap for hb: clip: bmap:                          *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      character c1*1, c2*1

      common /error/ m_err, l_err, k_err
      character m_err*200

      dimension ill(0:9), ilf(0:9)
      character dsin(0:9)*200
      dimension idsi(0:9)

      logical dnen2,deqn3

      dimension val(5), rcol(3), bcol(3), pcol(3)

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

      ierr = 0

      iline = 0
      ilizt = 4

      iclip = 0
      ihsb  = 1

      mclip = 0

      nfram = 0
      nclip = 0
      nbmap = 0
      npath = 0

      nfrmm = 0
      nclmm = 0
      nptmm = 0
      nbtmm = 0

      rcol(1) = -2.0
      rcol(2) =  1.0
      rcol(3) =  1.0

      bcol(1) = -r1max
      bcol(2) =  1.0
      bcol(3) =  1.0

      pcol(1) = -r1max
      pcol(2) =  1.0
      pcol(3) =  1.0

*-----------------------------------------------------------------------
*     read initial line for parameter
*-----------------------------------------------------------------------

      ic = in - 1

  100 ic = ic + 1

            if( ic .gt. icolm ) goto 200
            if( lum(ic) .eq. ' ' .or. lum(ic).eq.tub ) goto 100

         if( lum(ic)//lum(ic+1)//lum(ic+2)//lum(ic+3) .eq. 'clip' ) then

            iclip = 1
            ic = ic + 3

         else if( lum(ic)//lum(ic+1)//lum(ic+2)//lum(ic+3)//
     &            lum(ic+4)//lum(ic+5)  .eq. 'noclip' ) then

            iclip = 0
            ic = ic + 5

         else if( lum(ic)//lum(ic+1)//lum(ic+2)//lum(ic+3) .eq.
     &          'line' ) then

            iline = 1
            ic = ic + 3

         else if( lum(ic)//lum(ic+1)//lum(ic+2)//lum(ic+3)//
     &            lum(ic+4)//lum(ic+5)  .eq. 'noline' ) then

            iline = 0
            ic = ic + 5

         else if(lum(ic)//lum(ic+1).eq.'c(') then

            c1 = '('
            c2 = ')'

           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

            if(ierr.eq.1) goto 986

         else if(lum(ic)//lum(ic+1)//lum(ic+2).eq.'cb(') then

            ic = ic + 1

            c1 = '('
            c2 = ')'

           call dcols(1,lum,ic,bcol,ierr,icolm,c1,c2)

            if(ierr.eq.1) goto 986

         else if( lum(ic)//lum(ic+1)//lum(ic+2) .eq. 'hsb' ) then

            ihsb = 1
            ic = ic + 2

         else if( lum(ic)//lum(ic+1)//lum(ic+2) .eq. 'rgb' ) then

            ihsb = 0
            ic = ic + 2

         else if(lum(ic) .eq. 't' ) then

            ilizt = ilizt + 3

         else if( lum(ic) .eq. 'z' ) then

            ilizt = ilizt - 1

            if( ilizt .le. 0.0 ) ilizt = 1

         end if

         goto 100

*-----------------------------------------------------------------------

  200 continue

         ibmap = ibmap + 1

*-----------------------------------------------------------------------

         ila  = 0
         ila0 = 0
         iln  = 0

  510    ila  = ila  + 1
         ila0 = ila0 + 1

         if( mclip .eq. 1 .and. iln .eq. 2 .and. iclip .eq. 1 ) then

               nclmm = nclmm + 1
               write(jhb) val(1), val(2)

               if( val(1) .lt. xmin ) xmin = val(1)
               if( val(1) .gt. xmax ) xmax = val(1)
               if( val(2) .lt. ymin ) ymin = val(2)
               if( val(2) .gt. ymax ) ymax = val(2)

         else if( mclip .eq. 2 .and. iln .eq. 2 ) then

               nptmm = nptmm + 1
               write(jhq) val(1), val(2)

            if( iclip .eq. 0 ) then

               if( val(1) .lt. xmin ) xmin = val(1)
               if( val(1) .gt. xmax ) xmax = val(1)
               if( val(2) .lt. ymin ) ymin = val(2)
               if( val(2) .gt. ymax ) ymax = val(2)

            end if

         else if( mclip .eq. 3 .and. iln .eq. 5 ) then

               nbtmm = nbtmm + 1
               write(jhq) val(1), val(2), val(3), val(4), val(5)

            if( iclip .eq. 0 ) then

               if( val(1) - wid / 2.0 .lt. xmin )
     &               xmin = val(1) - wid / 2.0
               if( val(1) + wid / 2.0 .gt. xmax )
     &               xmax = val(1) + wid / 2.0
               if( val(2) - hgt / 2.0 .lt. ymin )
     &               ymin = val(2) - hgt / 2.0
               if( val(2) + hgt / 2.0 .gt. ymax )
     &               ymax = val(2) + hgt / 2.0

            end if

         else if( mclip .eq. 4 .and. iln .eq. 2 ) then

               nfrmm = nfrmm + 1
               write(jhb) val(1), val(2)

         end if

         iln  = 0

*-----------------------------------------------------------------------

         read(jsi,'(10000a1)', iostat = ios ) (dum(ic),ic=1,icolm)
         if( ios .eq. -1 ) goto 560

         call chlow(dum,lum)

            if( ill(jsn) + ila .gt. ilf(jsn) ) goto 560

*-----------------------------------------------------------------------
*        skip blanck and blank line is the end of the data
*-----------------------------------------------------------------------

            k = 1
            do 145 l = 1, icolm
               if(lum(l).ne.' '.and.lum(l).ne.tub) goto 146
  145       continue

            goto 560

  146       k = l

*-----------------------------------------------------------------------
*     INCLFL: INCLUDE FILE ONLY FOR ILA0 = 1
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

*-----------------------------------------------------------------------
*        frame: subsection
*-----------------------------------------------------------------------

         else if( lum(k)//lum(k+1)//lum(k+2)//lum(k+3)//
     &            lum(k+4)//lum(k+5) .eq. 'frame:' ) then

            if( nbmap .gt. 0 ) goto 981
            if( npath .gt. 0 ) goto 981
            if( nclip .gt. 0 ) goto 981

            nfram = nfram + 1
            mclip = 4
            if( nfram .gt. 1 ) goto 980

            ila0 = ila0 - 1
            goto 510

*-----------------------------------------------------------------------
*        clip: subsection
*-----------------------------------------------------------------------

         else if( lum(k)//lum(k+1)//lum(k+2)//lum(k+3)//lum(k+4) .eq.
     &            'clip:' ) then

            if( nbmap .gt. 0 ) goto 991
            if( npath .gt. 0 ) goto 991

            nclip = nclip + 1
            mclip = 1
            if( nclip .gt. 1 ) goto 998

            ila0 = ila0 - 1
            goto 510

*-----------------------------------------------------------------------
*        path: subsection
*-----------------------------------------------------------------------

         else if( lum(k)//lum(k+1)//lum(k+2)//lum(k+3)//lum(k+4) .eq.
     &            'path:' ) then

            if( iclip .eq. 1 .and. nclip .eq. 0 ) goto 990

*-----------------------------------------------------------------------

            if( mclip .eq. 2 ) then

                  if( nptmm .eq. 0 ) goto 978

                  rewind(jhq)

                  write(jhp) nptmm, pcol

               do j = 1, nptmm

                  read(jhq) xv, yv
                  write(jhp) xv, yv

               end do

            else if( mclip .eq. 3 ) then

                  if( nbtmm .eq. 0 ) goto 977

                  rewind(jhq)

                  write(jhr) nbtmm, wid, hgt

               do j = 1, nbtmm

                  read(jhq) xv, yv, cl1, cl2, cl3
                  write(jhr) xv, yv, cl1, cl2, cl3

               end do

            end if

*-----------------------------------------------------------------------

               npath = npath + 1
               mclip = 2

               rewind(jhq)

               nptmm = 0

            ic =  k + 4
  110       ic = ic + 1

            if( ic .gt. icolm ) goto 983
            if( lum(ic) .eq. ' ' .or. lum(ic).eq.tub ) goto 110

            if(lum(ic)//lum(ic+1).eq.'c(') then

               c1 = '('
               c2 = ')'

              call dcols(1,lum,ic,pcol,ierr,icolm,c1,c2)

               if(ierr.eq.1) goto 984

            end if

            ila0 = ila0 - 1
            goto 510

*-----------------------------------------------------------------------
*        bmap: subsection
*-----------------------------------------------------------------------

         else if( lum(k)//lum(k+1)//lum(k+2)//lum(k+3)//lum(k+4) .eq.
     &            'bmap:' ) then

            if( iclip .eq. 1 .and. nclip .eq. 0 ) goto 990

*-----------------------------------------------------------------------

            if( mclip .eq. 2 ) then

                  if( nptmm .eq. 0 ) goto 978

                  rewind(jhq)

                  write(jhp) nptmm, pcol

               do j = 1, nptmm

                  read(jhq) xv, yv
                  write(jhp) xv, yv

               end do

            else if( mclip .eq. 3 ) then

                  if( nbtmm .eq. 0 ) goto 977

                  rewind(jhq)

                  write(jhr) nbtmm, wid, hgt

               do j = 1, nbtmm

                  read(jhq) xv, yv, cl1, cl2, cl3
                  write(jhr) xv, yv, cl1, cl2, cl3

               end do

            end if

*-----------------------------------------------------------------------

            nbmap = nbmap + 1
            mclip = 3

               rewind(jhq)

               nbtmm = 0
               iwid  = 0
               ihgt  = 0

            ic =  k + 4
  120       ic = ic + 1

            if( ic .gt. icolm ) goto 130
            if( lum(ic) .eq. ' ' .or. lum(ic).eq.tub ) goto 120

         if( lum(ic)//lum(ic+1)//lum(ic+2)//lum(ic+3)//
     &       lum(ic+4) .eq. 'width' ) then

            ic = ic + 3
  101       ic = ic + 1
            if( ic .gt. icolm ) goto 979
            if( lum(ic) .ne. '=' ) goto 101
  102       ic = ic + 1
            if( ic .gt. icolm ) goto 979
            if( lum(ic) .eq. ' ' .or. lum(ic).eq.tub ) goto 102
               ici = ic
  103       ic = ic + 1
            if( ic .gt. icolm ) goto 979
            if(deqn3(lum(ic))) goto 103
            if(lum(ic).ne.' '.and.lum(ic).ne.tub) goto 995
               icf = ic
               call rnum(wid,lum,ici,icf,ierr)
               if(ierr.ne.0) goto 995
               iwid = iwid + 1

         else if( lum(ic)//lum(ic+1)//lum(ic+2)//lum(ic+3)//
     &            lum(ic+4) .eq. 'hight' ) then

            ic = ic + 3
  111       ic = ic + 1
            if( ic .gt. icolm ) goto 979
            if( lum(ic) .ne. '=' ) goto 111
  112       ic = ic + 1
            if( ic .gt. icolm ) goto 979
            if( lum(ic) .eq. ' ' .or. lum(ic).eq.tub ) goto 112
               ici = ic
  113       ic = ic + 1
            if( ic .gt. icolm ) goto 979
            if(deqn3(lum(ic))) goto 113
            if(lum(ic).ne.' '.and.lum(ic).ne.tub) goto 994
               icf = ic
               call rnum(hgt,lum,ici,icf,ierr)
               if(ierr.ne.0) goto 994
               ihgt = ihgt + 1

         end if

            goto 120

  130    continue

            if( iwid .eq. 0 ) goto 993
            if( ihgt .eq. 0 ) goto 992

            ila0 = ila0 - 1
            goto 510

*-----------------------------------------------------------------------
*        end of the HB: section
*-----------------------------------------------------------------------

         else if( lum(k+1) .eq. ':' .or.
     &            lum(k+2) .eq. ':' .or.
     &            dnen2(lum(k)) ) then

            goto 560

         end if

*-----------------------------------------------------------------------
*     CHECK OF NUMERICS
*-----------------------------------------------------------------------

         ic  = k - 1
  530    ic = ic + 1

            if( ic .gt. icolm .and. iln .eq. 0 ) goto 560
            if( ic .gt. icolm ) goto 510
            if( lum(ic) .eq. ' ' .or. lum(ic) .eq. tub ) goto 530

            if( dnen2(lum(ic)) ) goto 997

            ici = ic

*-----------------------------------------------------------------------

  550    ic = ic + 1

            if( ic .gt. icolm .and. iln .eq. 0 ) goto 560

            if( ic .gt. icolm ) then

               call rnum(rrnm,lum,ici,ic-1,ierr)
               if(ierr.ne.0) goto 997

               iln = iln + 1
               if( mclip .eq. 1 .and. iln .gt. 2 ) goto 989
               if( mclip .eq. 2 .and. iln .gt. 2 ) goto 985
               if( mclip .eq. 3 .and. iln .gt. 5 ) goto 988
               if( mclip .eq. 4 .and. iln .gt. 2 ) goto 985

               val(iln) = rrnm

               goto 510

            end if

            if(  lum(ic) .eq. ' ' .or. lum(ic) .eq. tub ) then

               call rnum(rrnm,lum,ici,ic-1,ierr)
               if(ierr.ne.0) goto 997

               iln = iln + 1
               if( mclip .eq. 1 .and. iln .gt. 2 ) goto 989
               if( mclip .eq. 2 .and. iln .gt. 2 ) goto 985
               if( mclip .eq. 3 .and. iln .gt. 5 ) goto 988
               if( mclip .eq. 4 .and. iln .gt. 2 ) goto 985

               val(iln) = rrnm

               goto 530

            end if

            if( deqn3(lum(ic)) ) goto 550

            goto 997

*-----------------------------------------------------------------------
*        end of check of data
*-----------------------------------------------------------------------

  560    continue

*-----------------------------------------------------------------------

            if( mclip .eq. 2 ) then

                  if( nptmm .eq. 0 ) goto 978

                  rewind(jhq)

                  write(jhp) nptmm, pcol

               do j = 1, nptmm

                  read(jhq) xv, yv
                  write(jhp) xv, yv

               end do

            else if( mclip .eq. 3 ) then

                  if( nbtmm .eq. 0 ) goto 977

                  rewind(jhq)

                  write(jhr) nbtmm, wid, hgt

               do j = 1, nbtmm

                  read(jhq) xv, yv, cl1, cl2, cl3
                  write(jhr) xv, yv, cl1, cl2, cl3

               end do

            end if

*-----------------------------------------------------------------------

         if( iclip .eq. 1 .and. nclmm .eq. 0 ) goto 990
         if( iclip .eq. 0 .and. nbmap .eq. 0 .and.
     &                          npath .eq. 0 ) goto 987

         write(jhd) iclip, nfrmm, nclmm, nbmap, npath,
     &              ihsb, iline, ilizt, rcol, bcol

*-----------------------------------------------------------------------

            ill(jsn) = ill(jsn) + ila

*-----------------------------------------------------------------------

      return

  977 m_err = 'bmap: there is no data'
      ErrCha = ''
      ErrID = 'L:3660/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  978 m_err = 'path: there is no data'
      ErrCha = ''
      ErrID = 'L:3667/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  979 m_err = 'bmap: width= or hight= is missing.'
      ErrCha = ''
      ErrID = 'L:3674/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  980 m_err = 'hb: fram: subsection is appeared twice.'
      ErrCha = ''
      ErrID = 'L:3681/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  981 m_err = 'hb: frame: should be before clip: bmap: or path:'
      ErrCha = ''
      ErrID = 'L:3688/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  983 m_err = 'path: C( ) is missing.'
      ErrCha = ''
      ErrID = 'L:3695/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  984 m_err = 'path: Number in C( ) is Wrong.'
      ErrCha = ''
      ErrID = 'L:3702/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  985 m_err = 'hb: frame: or path: should be 2 numbers: x y'
      ErrCha = ''
      ErrID = 'L:3709/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  986 m_err = 'hb: Number in C( ) is Wrong.'
      ErrCha = ''
      ErrID = 'L:3716/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  987 m_err = 'hb: bmap: and path: is missing'
      ErrCha = ''
      ErrID = 'L:3723/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  988 m_err = 'hb: bmap: should be 5 numbers: x y c1 c2 c3'
      ErrCha = ''
      ErrID = 'L:3730/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  989 m_err = 'hb: clip: should be two numbers: x y'
      ErrCha = ''
      ErrID = 'L:3737/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  990 m_err = 'hb: clip: is missing'
      ErrCha = ''
      ErrID = 'L:3744/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  991 m_err = 'hb: clip: should be before bmap: or path:'
      ErrCha = ''
      ErrID = 'L:3751/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  992 m_err = 'bmap: hight = is missing'
      ErrCha = ''
      ErrID = 'L:3758/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  993 m_err = 'bmap: width = is missing'
      ErrCha = ''
      ErrID = 'L:3765/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  994 m_err = 'bmap: number of hight = is wrong'
      ErrCha = ''
      ErrID = 'L:3772/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  995 m_err = 'bmap: number of width = is wrong'
      ErrCha = ''
      ErrID = 'L:3779/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  997 m_err = 'hb: Numerical Data is Something Wrong.'
      ErrCha = ''
      ErrID = 'L:3786/R:bmap/F:a-main1.f'
      l_err = ill(jsn)+ila
      k_err = jsn

      goto 999

  998 m_err = 'hb: clip: subsection is appeared twice.'
      ErrCha = ''
      ErrID = 'L:3794/R:bmap/F:a-main1.f'
      l_err = ill(jsn)
      k_err = jsn
      goto 999

  999 ierr  = 1

      return
      end

************************************************************************
*                                                                      *
      subroutine wcmt(dum,lum,in0,iw,wxys,iwx,iwy,iwn,ierr,ifon,
     &                ill,ilf,jsi,jsn,iend,iwf,iwb,ibox,cbox)
*                                                                      *
*      PURPOSE  :   READ COMMENTS FOR LONG TEXT                        *
*                                                                      *
*      FORMAT   :  WT: X( ) Y( ) S( ) IX( ) IY( ) C( ) F( ) B( )       *
*                      BOX(singlebox) CB( ) CL( ) CS( ) A( )           *
*                  text                                                *
*                  text                                                *
*                  E:   END MARKER                                     *
*                                                                      *
*           IW    :   NUMBER OF COMMENTS                               *
*                                                                      *
*           WXYS(IW,9) :   1-> X-COORDINATE   X( )                     *
*                          2-> Y-CORDINATE    Y( )                     *
*                          3-> SIZE OF CHARACTER REL TO NORMAL  S( )   *
*                          4,5,6 -> COLOR OF CHARACTERS  C( )          *
*                          7-> BASE LINE SKIP REL TO FONT SIZE B( )    *
*                              DEFAULT IS 1.5 FOR JAPANESE             *
*                              ==>only    1.2 FOR ENGLISH              *
*                          8-> ANGLE  A( )  DEGREE                     *
*                          9-> BASE LINE SKIP FOR BLANK LINE BB( )     *
*                                                                      *
*           IWX(IW)     :  X JUSTIFICATION OF TEXT IX( )               *
*                         = 1 ; LEFT     * DEFAUT                      *
*                           2 ; CENTER                                 *
*                           3 ; RIGHT                                  *
*                                                                      *
*           IWY(IW)     :  Y JUSTIFICATION OF THE FIRST LINE IY( )     *
*                         = 1 ; BOTTOM                                 *
*                           2 ; CENTER                                 *
*                           3 ; TOP      * DEFAULT                     *
*                                                                      *
*           WCM(IW,ICHRL) :  TEXT                                      *
*                                                                      *
*           IWN(IW)     :  LENGTH OF COMMENTS                          *
*                                                                      *
*           IWF(IW)     :  FONT OF COMMENTS                            *
*                                                                      *
*           IWB(IW)     :  BASELINESKIP OF COMMENTS            IPIN    *
*                          -1   -> ONE LINE                       0    *
*                          -2   -> INITIAL LINE OF MULTI LINES    1    *
*                          IW-1 -> MULTI LINES                    2    *
*                          -3   -> END LINE OF MULTI LINES        3    *
*                                                                      *
*           CBOX(IW,3)  : COLOR OF BOX                                 *
*                         CBOX(IW,1) : COLOR OF BACK GROUND            *
*                         CBOX(IW,2) : COLOR OF LINE                   *
*                         CBOX(IW,3) : COLOR OF SHADOW                 *
*                                                                      *
*           IBOX(IW)    : 0-> NO BOX                                   *
*                                                                      *
*                        10 ; singlebox                                *
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
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      dimension wxys(mc,9),iwx(mc),iwy(mc),iwn(mc),iwf(mc),iwb(mc)
      dimension ibox(mc),cbox(mc,3,3)
      character wcm(ichrl)*1
      character m_err*200
      common /error/ m_err, l_err, k_err


      character rum*2048
*         ( ichrl = 2048 )

      character c1*1, c2*1

      dimension ill(0:9), ilf(0:9)
      dimension rcol(3)

      character yen*1
      character tub*1
      yen = char(92)
      tub = char(9)

*-----------------------------------------------------------------------

               c1 = '('
               c2 = ')'

               iwx(iw) = 1
               iwy(iw) = 3

               iwn(iw) = 0
               iwf(iw) = -1
               iwb(iw) = -2

               wxys(iw,1) = -r1max
               wxys(iw,2) = -r1max
               wxys(iw,3) = 1.0
               wxys(iw,4) = -r1max
               wxys(iw,5) = 1.0
               wxys(iw,6) = 1.0
               wxys(iw,7) = -r1max
               wxys(iw,8) = -r1max
               wxys(iw,9) = -r1max

               ibox(iw)   = 0
               cbox(iw,1,1) = -r1max
               cbox(iw,2,1) = -r1max
               cbox(iw,3,1) = -r1max
               cbox(iw,1,2) = 1.0
               cbox(iw,1,3) = 1.0
               cbox(iw,2,2) = 1.0
               cbox(iw,2,3) = 1.0
               cbox(iw,3,2) = 1.0
               cbox(iw,3,3) = 1.0

               iend = 0

               ic   = in0 - 1

*-----------------------------------------------------------------------

  200 ic = ic + 1

         if( ic .gt. icolm ) goto 700

         if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 200

      if(lum(ic).ne.'x'.and.lum(ic).ne.'y'.and.lum(ic).ne.'s'.and.
     &   lum(ic).ne.'c'.and.lum(ic).ne.'f'.and.lum(ic).ne.'b'.and.
     &   lum(ic).ne.'a'.and.lum(ic).ne.'i') then

         m_err = 'WT: Text/ X() Y() S() IX() IY() C() F() B() A() '//
     &           'BOX() CB() CL() CS() BB(); Unexpected Parameter.'
         ErrCha = ''
         ErrID = 'L:3966/R:wcmt/F:a-main1.f'
         goto 999

      end if


      if(lum(ic).eq.'x') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WT: Number in X( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:3979/R:wcmt/F:a-main1.f'
               goto 999
            end if
            wxys(iw,1)=rnm
            ic=ic+1

      else if(lum(ic).eq.'y') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WT: Number in Y( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:3992/R:wcmt/F:a-main1.f'
               goto 999
            end if
            wxys(iw,2)=rnm
            ic=ic+1

      else if(lum(ic).eq.'a') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WT: Number in A( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4005/R:wcmt/F:a-main1.f'
               goto 999
            end if
            wxys(iw,8)=rnm
            ic=ic+1

      else if(lum(ic).eq.'s') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WT: Number in S( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4018/R:wcmt/F:a-main1.f'
               goto 999
            end if
            wxys(iw,3)=rnm
            ic=ic+1

      else if( lum(ic)//lum(ic+1)//lum(ic+2) .eq. 'box' ) then

               ic = ic + 2

               if(lum(ic+1) .ne. '(') then
                  m_err = 'WT: Description of BOX( ) is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:4031/R:wcmt/F:a-main1.f'
                  goto 999
               end if

               ic = ic + 1

  550          ic = ic + 1

               if( ic .gt. icolm ) then
                  m_err = 'WT: Description of BOX( ) is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:4042/R:wcmt/F:a-main1.f'
                  goto 999
               end if

               if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 550


               do 500 k = ic, icolm

                  kk = k - ic + 1

                  rum(kk:kk) = dum(k)

  500          continue

         if( ( rum(1: 9) .eq. 'singlebox' ) .or.
     &       ( rum(1: 9) .eq. 'Singlebox' ) .or.
     &       ( rum(1: 9) .eq. 'singleBox' ) .or.
     &       ( rum(1: 9) .eq. 'SingleBox' ) .or.
     &       ( rum(1: 7) .eq. 'ovalbox'   ) .or.
     &       ( rum(1: 7) .eq. 'Ovalbox'   ) .or.
     &       ( rum(1: 7) .eq. 'ovalBox'   ) .or.
     &       ( rum(1: 7) .eq. 'OvalBox'   ) .or.
     &       ( rum(1: 9) .eq. 'doublebox' ) .or.
     &       ( rum(1: 9) .eq. 'Doublebox' ) .or.
     &       ( rum(1: 9) .eq. 'doubleBox' ) .or.
     &       ( rum(1: 9) .eq. 'DoubleBox' ) .or.
     &       ( rum(1: 9) .eq. 'shadowbox' ) .or.
     &       ( rum(1: 9) .eq. 'Shadowbox' ) .or.
     &       ( rum(1: 9) .eq. 'shadowBox' ) .or.
     &       ( rum(1: 9) .eq. 'ShadowBox' ) .or.
     &       ( rum(1:13) .eq. 'ovalshadowbox' ) .or.
     &       ( rum(1:13) .eq. 'Ovalshadowbox' ) .or.
     &       ( rum(1:13) .eq. 'ovalshadowBox' ) .or.
     &       ( rum(1:13) .eq. 'OvalshadowBox' )
     &     ) then

               if( rum(1:9) .eq. 'singlebox' ) then

                     iboxd = 10
                     iboxl = 9

               else if( rum(1:9) .eq. 'Singlebox' ) then

                     iboxd = 11
                     iboxl = 9

               else if( rum(1:9) .eq. 'singleBox' ) then

                     iboxd = 12
                     iboxl = 9

               else if( rum(1:9) .eq. 'SingleBox' ) then

                     iboxd = 13
                     iboxl = 9

               else if( rum(1:7) .eq. 'ovalbox' ) then

                     iboxd = 20
                     iboxl = 7

               else if( rum(1:7) .eq. 'Ovalbox' ) then

                     iboxd = 21
                     iboxl = 7

               else if( rum(1:7) .eq. 'ovalBox' ) then

                     iboxd = 22
                     iboxl = 7

               else if( rum(1:7) .eq. 'OvalBox' ) then

                     iboxd = 23
                     iboxl = 7

               else if( rum(1:9) .eq. 'doublebox' ) then

                     iboxd = 30
                     iboxl = 9

               else if( rum(1:9) .eq. 'Doublebox' ) then

                     iboxd = 31
                     iboxl = 9

               else if( rum(1:9) .eq. 'doubleBox' ) then

                     iboxd = 32
                     iboxl = 9

               else if( rum(1:9) .eq. 'DoubleBox' ) then

                     iboxd = 33
                     iboxl = 9

               else if( rum(1:9) .eq. 'shadowbox' ) then

                     iboxd = 40
                     iboxl = 9

               else if( rum(1:9) .eq. 'Shadowbox' ) then

                     iboxd = 41
                     iboxl = 9

               else if( rum(1:9) .eq. 'shadowBox' ) then

                     iboxd = 42
                     iboxl = 9

               else if( rum(1:9) .eq. 'ShadowBox' ) then

                     iboxd = 43
                     iboxl = 9

               else if( rum(1:13) .eq. 'ovalshadowbox' ) then

                     iboxd = 50
                     iboxl = 13

               else if( rum(1:13) .eq. 'Ovalshadowbox' ) then

                     iboxd = 51
                     iboxl = 13

               else if( rum(1:13) .eq. 'ovalshadowBox' ) then

                     iboxd = 52
                     iboxl = 13

               else if( rum(1:13) .eq. 'OvalshadowBox' ) then

                     iboxd = 53
                     iboxl = 13

               end if

         else

                  m_err = 'WT: Description of BOX( ) is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:4185/R:wcmt/F:a-main1.f'
                  goto 999

         end if

               ic = ic + iboxl - 1

  560          ic = ic + 1

               if( ic .gt. icolm ) then
                  m_err = 'WT: Description of BOX( ) is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:4197/R:wcmt/F:a-main1.f'
                  goto 999
               end if

               if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 560

               if( lum(ic) .ne. ')' ) then
                  m_err = 'WT: Description of BOX( ) is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:4206/R:wcmt/F:a-main1.f'
                  goto 999
               end if


               ibox(iw) = iboxd

*-----------------------------------------------------------------------

      else if( lum(ic)//lum(ic+1).eq.'bb') then

            ic=ic+2
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WT: Number in BB( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4222/R:wcmt/F:a-main1.f'
               goto 999
            end if
            wxys(iw,9)=rnm
            ic=ic+1

      else if(lum(ic).eq.'b') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WT: Number in B( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4235/R:wcmt/F:a-main1.f'
               goto 999
            end if
            wxys(iw,7)=rnm
            ic=ic+1

      else if(lum(ic).eq.'f') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WT: Number in F( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4248/R:wcmt/F:a-main1.f'
               goto 999
            end if
            iwf(iw)=nint(rnm)
            if( iwf(iw) .lt. 0 .or. iwf(iw) .gt. 12 ) iwf(iw) = -1
            ic=ic+1

      else if( lum(ic)//lum(ic+1) .eq. 'cb' ) then

                  ic = ic + 1

            call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

                  cbox(iw,1,1) = rcol(1)
                  cbox(iw,1,2) = rcol(2)
                  cbox(iw,1,3) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'WT: Number in CB( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4268/R:wcmt/F:a-main1.f'
               goto 999
            end if

      else if( lum(ic)//lum(ic+1) .eq. 'cl' ) then

                  ic = ic + 1

            call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

                  cbox(iw,2,1) = rcol(1)
                  cbox(iw,2,2) = rcol(2)
                  cbox(iw,2,3) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'WT: Number in CL( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4285/R:wcmt/F:a-main1.f'
               goto 999
            end if

      else if( lum(ic)//lum(ic+1) .eq. 'cs' ) then

                  ic = ic + 1

            call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

                  cbox(iw,3,1) = rcol(1)
                  cbox(iw,3,2) = rcol(2)
                  cbox(iw,3,3) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'WT: Number in CS( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4302/R:wcmt/F:a-main1.f'
               goto 999
            end if

      else if(lum(ic).eq.'c') then

            call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               wxys(iw,4) = rcol(1)
               wxys(iw,5) = rcol(2)
               wxys(iw,6) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'WT: Number in C( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4317/R:wcmt/F:a-main1.f'
               goto 999
            end if

      else if( lum(ic)//lum(ic+1) .eq. 'ix' ) then

            ic=ic+2
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WT: Number in IX( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4328/R:wcmt/F:a-main1.f'
               goto 999
            end if
            if(nint(rnm).ne.1.and.nint(rnm).ne.2.and.
     &         nint(rnm).ne.3) then
               m_err = 'WT: Number in IX() is 1(Left), '//
     &                 '2(Center), or 3(Right).'
               ErrCha = ''
               ErrID = 'L:4336/R:wcmt/F:a-main1.f'
               goto 999
            end if
            iwx(iw)=nint(rnm)
            ic=ic+1

      else if( lum(ic)//lum(ic+1) .eq. 'iy' ) then

            ic=ic+2
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WT: Number in IY( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4349/R:wcmt/F:a-main1.f'
               goto 999
            end if
            if(nint(rnm).ne.1.and.nint(rnm).ne.2.and.
     &         nint(rnm).ne.3) then
               m_err = 'WT: Number in IY() is 1(Bottom), '//
     &                 '2(Center), or 3(Top).'
               ErrCha = ''
               ErrID = 'L:4357/R:wcmt/F:a-main1.f'
               goto 999
            end if
            iwy(iw)=nint(rnm)
            ic=ic+1

      end if

            goto 200

*-----------------------------------------------------------------------

  700    continue

*-----------------------------------------------------------------------
*        READ NEW LINE
*-----------------------------------------------------------------------

  750       ill(jsn)  = ill(jsn) + 1


            read(jsi,'(10000a1)', iostat = ios ) (dum(ic),ic=1,icolm)
            if( ios .eq. -1 ) goto 350

            call chlow(dum,lum)

*----------------------------------------------------------------------*
*        END OF WT: SECTION BY E:
*----------------------------------------------------------------------*

            do 103 i = 1, icolm - 3

               if(lum(i).ne.' '.and.lum(i).ne.tub) goto 113

  103       continue

  113       continue

            if( lum(i)//lum(i+1) .eq. 'e:' ) then

                  iw = iw - 1

                  iwb(iw) = -3

                  return

            end if

*----------------------------------------------------------------------*
*        CHARACTER LENGTH
*----------------------------------------------------------------------*

            do 101 i = icolm, 1, -1

               if(lum(i).ne.' '.and.lum(i).ne.tub) goto 111

  101       continue

*----------------------------------------------------------------------*
*        BLANK LINE
*----------------------------------------------------------------------*

               iwn(iw) = -1

               goto 755

*----------------------------------------------------------------------*
*        SEQUENTIAL LINE
*----------------------------------------------------------------------*

  111       icf = i
            ict = 0

                  iseql = 0
               if( icf .eq. 1 .and. dum(icf) .eq. yen )
     &            iseql = 1
               if( icf .gt. 1 ) then
                  if( dum(icf) .eq. yen .and. dum(icf-1) .ne. yen .and.
     &                                  ichar(dum(icf-1)) .le. 128 )
     &            iseql = 2
               end if

         if( iseql .gt. 0 ) then

  410          icf = icf - 1

               do 401 i = 1, icf

                  wcm( ict + i ) = dum( i )

  401          continue

                  iwn(iw) = ict + icf
                  ict = ict + icf

                  ill(jsn)  = ill(jsn) + 1


                  read(jsi,'(10000a1)',iostat=ios) (dum(ic),ic=1,icolm)
                  if( ios .eq. -1 ) goto 350
                  call chlow(dum,lum)

               do 102 i = icolm, 1, -1

                  if(dum(i).ne.' '.and.dum(i).ne.tub) goto 112

  102          continue
  112          icf = i

                  jseql = 0
               if( icf .eq. 1 .and. dum(icf) .eq. yen )
     &            jseql = 1
               if( icf .gt. 1 ) then
                  if( dum(icf) .eq. yen .and. dum(icf-1) .ne. yen .and.
     &                                  ichar(dum(icf-1)) .le. 128 )
     &            jseql = 2
               end if

               if( jseql .gt. 0 ) goto 410

         end if

*----------------------------------------------------------------------*
*        NORMAL LINE
*----------------------------------------------------------------------*

               do 400 i = 1, icf

                  wcm( ict + i ) = dum( i )

  400          continue

                  iwn(iw) = ict + icf

*----------------------------------------------------------------------*
*        CHECK OF LANGAGE CODE
*----------------------------------------------------------------------*

               call jpncode(wcm,iwn(iw),ifon)
               call jpnprep(wcm,iwn(iw),ifon)

               write(jhm) (wcm(i),i=1,iwn(iw))

*----------------------------------------------------------------------*
*        NEXT LINE
*----------------------------------------------------------------------*

  755          iw = iw + 1

               if( iw .gt. mc ) then
                  m_err = 'Number of Comment is larger than mc'
                  ErrCha = ''
                  ErrID = 'L:4509/R:wcmt/F:a-main1.f'
                  goto 999
               end if

               iwx(iw) = iwx(iw-1)
               iwy(iw) = iwy(iw-1)

               iwn(iw) = 0
               iwf(iw) = iwf(iw-1)
               iwb(iw) = iw - 1

               wxys(iw,1) = wxys(iw-1,1)
               wxys(iw,2) = wxys(iw-1,2)
               wxys(iw,3) = wxys(iw-1,3)
               wxys(iw,4) = wxys(iw-1,4)
               wxys(iw,5) = wxys(iw-1,5)
               wxys(iw,6) = wxys(iw-1,6)
               wxys(iw,7) = wxys(iw-1,7)
               wxys(iw,8) = wxys(iw-1,8)
               wxys(iw,9) = wxys(iw-1,9)

               ibox(iw)   = ibox(iw-1)

               cbox(iw,1,1) = cbox(iw-1,1,1)
               cbox(iw,1,2) = cbox(iw-1,1,2)
               cbox(iw,1,3) = cbox(iw-1,1,3)
               cbox(iw,2,1) = cbox(iw-1,2,1)
               cbox(iw,2,2) = cbox(iw-1,2,2)
               cbox(iw,2,3) = cbox(iw-1,2,3)
               cbox(iw,3,1) = cbox(iw-1,3,1)
               cbox(iw,3,2) = cbox(iw-1,3,2)
               cbox(iw,3,3) = cbox(iw-1,3,3)

               goto 750

*----------------------------------------------------------------------*

  350          iend = 1
               return

  999          ierr = 1
               l_err = ill(jsn)
               k_err = jsn
               return
               end

************************************************************************
*                                                                      *
      subroutine wcmm(dum,lum,in0,iw,wxys,iwx,iwy,iwn,ierr,ifon,
     &                ill,ilf,jsi,jsn,iend,iwf,iwb)
*                                                                      *
*      PURPOSE  :   READ COMMENTS FOR ONE LINE COMMENT                 *
*                                                                      *
*      FORMAT   :  W: TEXT / X( ) Y( ) S( ) IX( ) IY( ) C( ) F( )      *
*                            A( )                                      *
*                                                                      *
*           IW    :   NUMBER OF COMMENTS                               *
*                                                                      *
*           WXYS(IW,9) :   1-> X-COORDINATE   X( )                     *
*                          2-> Y-CORDINATE    Y( )                     *
*                          3-> SIZE OF CHARACTER REL TO NORMAL  S( )   *
*                          4,5,6 -> COLOR OF CHARACTERS  C( )          *
*                          7-> BASE LINE SKIP REL TO FONT SIZE B( )    *
*                              DEFAULT IS 1.5 FOR JAPANESE             *
*                              ==>only    1.2 FOR ENGLISH              *
*                          8-> ANGLE  A( )  DEGREE                     *
*                          9-> BASE LINE SKIP FOR BLANK LINE BB( )     *
*                                                                      *
*           IWX(IW)     :  X JUSTIFICATION OF TEXT     IX( )           *
*                         = 1 ; LEFT     * DEFAUT                      *
*                           2 ; CENTER                                 *
*                           3 ; RIGHT                                  *
*                                                                      *
*           IWY(IW)     :  Y JUSTIFICATION OF TEXT     IY( )           *
*                         = 1 ; BOTTOM                                 *
*                           2 ; CENTER                                 *
*                           3 ; TOP                                    *
*                           0 ; BASELINE * DEFAULT                     *
*                                                                      *
*           WCM(IW,ICHRL) :  TEXT               W: TEXT /              *
*                                                                      *
*           IWN(IW)     :  LENGTH OF COMMENTS                          *
*                                                                      *
*           IWF(IW)     :  FONT OF COMMENTS                            *
*                                                                      *
*           IWB(IW)     :  BASELINESKIP OF COMMENTS (NOT USED HERE)    *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      dimension wxys(mc,9),iwx(mc),iwy(mc),iwn(mc),iwf(mc),iwb(mc)
      character wcm(ichrl)*1
      character m_err*200
      common /error/ m_err, l_err, k_err

      character c1*1, c2

      dimension ill(0:9), ilf(0:9)
      dimension rcol(3)

      character yen*1
      character tub*1
      yen = char(92)
      tub = char(9)

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------


*----------------------------------------------------------------------*

               c1 = '('
               c2 = ')'

               iwx(iw) = 1
               iwy(iw) = 0

               iwn(iw) = 0
               iwf(iw) = -1
               iwb(iw) = -1

               wxys(iw,1) = -r1max
               wxys(iw,2) = -r1max
               wxys(iw,3) = 1.0
               wxys(iw,4) = -r1max
               wxys(iw,5) = 1.0
               wxys(iw,6) = 1.0
               wxys(iw,7) = -r1max
               wxys(iw,8) = -r1max
               wxys(iw,9) = -r1max

               in   = in0
               ict  = 0
               iend = 0

*----------------------------------------------------------------------*

            do 101 i = icolm, in, -1

               if(lum(i).ne.' '.and.lum(i).ne.tub) goto 111

  101       continue

               return

  111       continue

*----------------------------------------------------------------------*

                  iseql = 0
               if( i .eq. 1 .and. dum(i) .eq. yen )
     &            iseql = 1
               if( i .gt. 1 ) then
                  if( dum(i) .eq. yen .and. dum(i-1) .ne. yen .and.
     &                                ichar(dum(i-1)) .le. 128 )
     &            iseql = 2
               end if

         if( iseql .gt. 0 ) then

  410          icf = i - 1

               do 401 i = 1, icf - in + 1

                  wcm( ict + i ) = dum( in + i - 1 )

  401          continue

                  ict = ict + icf - in + 1

                  ill(jsn)  = ill(jsn) + 1
                  in  = 1


                  read(jsi,'(10000a1)',iostat=ios) (dum(ic),ic=1,icolm)
                  if( ios .eq. -1 ) goto 350
                  call chlow(dum,lum)

               do 102 i = icolm, 1, -1

                  if(dum(i).ne.' '.and.dum(i).ne.tub) goto 112

  102          continue
  112          continue

                  jseql = 0
               if( icf .eq. 1 .and. dum(i) .eq. yen )
     &            jseql = 1
               if( icf .gt. 1 ) then
                  if( dum(i) .eq. yen .and. dum(i-1) .ne. yen .and.
     &                                ichar(dum(i-1)) .le. 128 )
     &            jseql = 2
               end if

               if( jseql .gt. 0 ) goto 410

         end if

*----------------------------------------------------------------------*

            goto 351
  350       continue

               iend = 1

  351       continue

               icms = 0

            do 100 i = icolm, in, -1

               if( lum(i) .eq. ']' ) icms = 1
               if( lum(i) .eq. '[' ) icms = 0

               if( icms .eq. 0 .and. lum(i) .eq. '/' ) goto 110

  100       continue

               m_err = 'W: / is Missing Between Text and Parameters.'
               ErrCha = ''
               ErrID = 'L:4740/R:wcmm/F:a-main1.f'
               goto 999

  110    icf = i - 1

            if( icf .eq. in - 1 ) then

               if( ict .eq. 0 ) then

                  m_err = 'W: There is No Text Before /.'
                  ErrCha = ''
                  ErrID = 'L:4751/R:wcmm/F:a-main1.f'
                  goto 999

               else

                  iwn(iw) = ict

               end if

            else

               do 400 i = 1, icf - in + 1

                  wcm( ict + i ) = dum( in + i - 1 )

  400          continue

                  iwn(iw) = ict + icf - in + 1

            end if

*----------------------------------------------------------------------*

               call jpncode(wcm,iwn(iw),ifon)
               call jpnprep(wcm,iwn(iw),ifon)

               write(jhm) (wcm(i),i=1,iwn(iw))

*----------------------------------------------------------------------*

      ic = icf + 1

  200 ic = ic + 1

         if( ic .gt. icolm ) return

         if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 200

      if(lum(ic).ne.'x'.and.lum(ic).ne.'y'.and.lum(ic).ne.'s'.and.
     &   lum(ic).ne.'c'.and.lum(ic).ne.'f'.and.lum(ic).ne.'a'.and.
     &   lum(ic).ne.'i') then

         m_err = 'W: Text/ X() Y() S() IX() IY() C() F() A(); '//
     &           'Unexpected Parameter.'
         ErrCha = ''
         ErrID = 'L:4796/R:wcmm/F:a-main1.f'
         goto 999

      end if


      if(lum(ic).eq.'x') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'W: Number in X( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4809/R:wcmm/F:a-main1.f'
               goto 999
            end if
            wxys(iw,1)=rnm
            ic=ic+1

      else if(lum(ic).eq.'y') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'W: Number in Y( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4822/R:wcmm/F:a-main1.f'
               goto 999
            end if
            wxys(iw,2)=rnm
            ic=ic+1

      else if(lum(ic).eq.'a') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'W: Number in A( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4835/R:wcmm/F:a-main1.f'
               goto 999
            end if
            wxys(iw,8)=rnm
            ic=ic+1

      else if(lum(ic).eq.'s') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'W: Number in S( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4848/R:wcmm/F:a-main1.f'
               goto 999
            end if
            wxys(iw,3)=rnm
            ic=ic+1

      else if(lum(ic).eq.'f') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'W: Number in F( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4861/R:wcmm/F:a-main1.f'
               goto 999
            end if
            iwf(iw)=nint(rnm)
            if( iwf(iw) .lt. 0 .or. iwf(iw) .gt. 12 ) iwf(iw) = -1
            ic=ic+1

      else if(lum(ic).eq.'c') then

            call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               wxys(iw,4) = rcol(1)
               wxys(iw,5) = rcol(2)
               wxys(iw,6) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'W: Number in C( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4879/R:wcmm/F:a-main1.f'
               goto 999
            end if

      else if(lum(ic)//lum(ic+1) .eq. 'ix' ) then

            ic=ic+2
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'W: Number in IX( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4890/R:wcmm/F:a-main1.f'
               goto 999
            end if
            if(nint(rnm).ne.1.and.nint(rnm).ne.2.and.
     &         nint(rnm).ne.3) then
               m_err = 'W: Number in IX() is 1(Left), '//
     &                 '2(Center), or 3(Right).'
               ErrCha = ''
               ErrID = 'L:4898/R:wcmm/F:a-main1.f'
               goto 999
            end if
            iwx(iw)=nint(rnm)
            ic=ic+1

      else if(lum(ic)//lum(ic+1) .eq. 'iy' ) then

            ic=ic+2
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'W: Number in IY( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:4911/R:wcmm/F:a-main1.f'
               goto 999
            end if
            if(nint(rnm).ne.1.and.nint(rnm).ne.2.and.
     &         nint(rnm).ne.3.and.nint(rnm).ne.0) then
               m_err = 'W: Number in IY() is 1(Bottom), '//
     &                 '2(Center), 3(Top), or 0(BaseLine).'
               ErrCha = ''
               ErrID = 'L:4919/R:wcmm/F:a-main1.f'
               goto 999
            end if
            iwy(iw)=nint(rnm)
            ic=ic+1

      end if


      goto 200

*-----------------------------------------------------------------------

  999 ierr = 1
      l_err = ill(jsn)
      k_err = jsn
      return
      end


************************************************************************
*                                                                      *
      subroutine wctb(dum,lum,in0,iv,ierr,ifon,ill,ilf,jsi,jsn,iend)
*                                                                      *
*      PURPOSE  :   READ COMMENTS FOR TABULAR                          *
*                                                                      *
*      FORMAT   :  WTAB: TAB or tab{|l|c|c|r|}                         *
*                        X( ) Y( ) S( ) IX( ) IY( ) C( ) F( )          *
*                        CB( ) CL( ) A( )                              *
*                  text & text & text & text                           *
*                  \hline \cline{2-3} \vline{|l|c|c|r|}                *
*                  \vspace{ } \tabtopsp{ }                             *
*                  text & text & text & text \selcolor{ }              *
*                  text & text & text & text                           *
*                  E:   END MARKER                                     *
*                                                                      *
*           IV    :   NUMBER OF TABLE                                  *
*                                                                      *
*           VXYS(13)    :  1-> X-COORDINATE   X( )                     *
*                          2-> Y-CORDINATE    Y( )                     *
*                          3-> ANGLE  A( )  DEGREE                     *
*                          4-> SIZE OF CHARACTER REL TO NORMAL  S( )   *
*                          5,6,7    -> COLOR OF CHARACTERS             *
*                          8,9,10   -> COLOR OF BACKGROUND             *
*                          11,12,13 -> COLOR OF LINE                   *
*                                                                      *
*                                                                      *
*           IVX(IV)     :  X JUSTIFICATION OF TABLE IX( )              *
*                         = 1 ; LEFT     * DEFAUT                      *
*                           2 ; CENTER                                 *
*                           3 ; RIGHT                                  *
*                                                                      *
*           IVY(IV)     :  Y JUSTIFICATION OF THE TABLE IY( )          *
*                         = 1 ; BOTTOM                                 *
*                           2 ; CENTER                                 *
*                           3 ; TOP      * DEFAULT                     *
*                                                                      *
*           IVR(IV)     :  NUMBER OF ROW ELEMENTS                      *
*           IVC(IV)     :  NUMBER OF COLUMN ELEMENTS                   *
*                                                                      *
*           IVS(IV)     :  SPACEING OF THE TABLE                       *
*                         = 0 SMALL                                    *
*                         = 1 LARGE                                    *
*                                                                      *
*           IVF(IV)     :  FONT OF COMMENTS                            *
*                                                                      *
*           VCM(IV,IR,IC,ICHRL) :  TEXT                                *
*           LVM(IV,IR,IC)       :  LENGTH OF TEXT                      *
*                                                                      *
*           VSPT(IV,1,IR)       :  VSPACE OF THE LINE                  *
*           VSPT(IV,2,IR)       :  VSPACE OF TEXT (TABTOPSPACE)        *
*                                                                      *
*           ITTX(IV,IR,IC)      :  X POSITION JUSTIFICATION            *
*                                  = 1 ; LEFT     * DEFAUT             *
*                                    2 ; CENTER                        *
*                                    3 ; RIGHT                         *
*                                                                      *
*           ITTH(IV,IR,IC)      :  HORIZONTAL LINE                     *
*           ITTV(IV,IR,IC)      :  VERTICAL LINE                       *
*                                                                      *
*                                 = 0 NO LINE                          *
*                                 = 1 ONE LINE                         *
*                                 = 2 DOUBLE LINE                      *
*                                                                      *
*                                 = -1 ONE LINE ( THICK LINE )         *
*                                 = -2 DOUBLE LINE ( THICK LINE )      *
*                                                                      *
*           ITTL(IV)            : 0-> NO LINE                          *
*                                 1-> WITH LINE                        *
*                                                                      *
*           CSL(IV,IR,IC)       :  COLOR OF SEL                        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      character jdum(ichrl)*1
      character jlum(ichrl)*1

*-----------------------------------------------------------------------

      dimension vxys(13)
      dimension vspt(2,0:nr)
      dimension ittx(0:nr,0:nc), ittv(0:nr,0:nc), itth(0:nr,0:nc)
      dimension csl(0:nr,0:nc,3)
      dimension lvm(0:nr,0:nc)
      dimension rcol(3)

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

      character c1*1, c2*1

      dimension ill(0:9), ilf(0:9)

      logical spjpn

      character yen*1
      character tub*1
      yen = char(92)
      tub = char(9)

*-----------------------------------------------------------------------

            open(jhs,form='unformatted',status='scratch')

*-----------------------------------------------------------------------
*        INITIALIZATION
*-----------------------------------------------------------------------

               c1 = '('
               c2 = ')'

               ivx = 1
               ivy = 3

               ivf = -1

               ivs = 0

               ivr = 0
               ivc = 0

            do 10 i = 0, nr

               vspt(1,i) = 0.0
               vspt(2,i) = 0.0

            do 10 j = 1, nc

               lvm(i,j) = 0
               csl(i,j,1) = -r1max
               csl(i,j,2) = 1.0
               csl(i,j,3) = 1.0

   10       continue

            do 20 i = 0, nc
            do 20 j = 0, nr

               ittx(j,i) = 1
               ittv(j,i) = 0
               itth(j,i) = 0

   20       continue

               vxys( 1) = -r1max
               vxys( 2) = -r1max
               vxys( 3) = -r1max
               vxys( 4) = 1.0
               vxys( 5) = -r1max
               vxys( 6) = 1.0
               vxys( 7) = 1.0
               vxys( 8) = -r1max
               vxys( 9) = 1.0
               vxys(10) = 1.0
               vxys(11) = -r1max
               vxys(12) = 1.0
               vxys(13) = 1.0

               l_err = ill(jsn)
               k_err = jsn

*-----------------------------------------------------------------------

               iend = 0

               ic   = in0 - 1

*-----------------------------------------------------------------------

  200 ic = ic + 1

         if( ic .gt. icolm ) goto 700

         if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 200

      if(lum(ic).ne.'x'.and.lum(ic).ne.'y'.and.lum(ic).ne.'s'.and.
     &   lum(ic).ne.'c'.and.lum(ic).ne.'f'.and.lum(ic).ne.'t'.and.
     &   lum(ic).ne.'a'.and.lum(ic).ne.'i') then

         m_err = 'WTAB: X() Y() S() IX() IY() C() F() A() '//
     &           'tab{} CB() CL() ; Unexpected Parameter.'
         ErrCha = ''
         ErrID = 'L:5134/R:wctb/F:a-main1.f'
         goto 999

      end if


      if(lum(ic).eq.'x') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WTAB: Number in X( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:5147/R:wctb/F:a-main1.f'
               goto 999
            end if
            vxys(1)=rnm
            ic=ic+1

      else if(lum(ic).eq.'y') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WTAB: Number in Y( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:5160/R:wctb/F:a-main1.f'
               goto 999
            end if
            vxys(2)=rnm
            ic=ic+1

      else if(lum(ic).eq.'a') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WTAB: Number in A( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:5173/R:wctb/F:a-main1.f'
               goto 999
            end if
            vxys(3)=rnm
            ic=ic+1

      else if(lum(ic).eq.'s') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WTAB: Number in S( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:5186/R:wctb/F:a-main1.f'
               goto 999
            end if
            vxys(4)=rnm
            ic=ic+1

      else if(lum(ic).eq.'f') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WTAB: Number in F( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:5199/R:wctb/F:a-main1.f'
               goto 999
            end if
            ivf=nint(rnm)
            if( ivf .lt. 0 .or. ivf .gt. 12 ) ivf = -1
            ic=ic+1

      else if(lum(ic)//lum(ic+1).eq.'cb') then

                  ic = ic + 1

            call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               vxys( 8) = rcol(1)
               vxys( 9) = rcol(2)
               vxys(10) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'WTAB: Number in CB( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:5219/R:wctb/F:a-main1.f'
               goto 999
            end if

      else if(lum(ic)//lum(ic+1).eq.'cl') then

                  ic = ic + 1

            call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               vxys(11) = rcol(1)
               vxys(12) = rcol(2)
               vxys(13) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'WTAB: Number in CL( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:5236/R:wctb/F:a-main1.f'
               goto 999
            end if

      else if(lum(ic).eq.'c') then

            call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               vxys(5) = rcol(1)
               vxys(6) = rcol(2)
               vxys(7) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'WTAB: Number in C( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:5251/R:wctb/F:a-main1.f'
               goto 999
            end if

      else if(lum(ic)//lum(ic+1).eq.'ix') then

            ic=ic+2
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WTAB: Number in IX( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:5262/R:wctb/F:a-main1.f'
               goto 999
            end if
            if(nint(rnm).ne.1.and.nint(rnm).ne.2.and.
     &         nint(rnm).ne.3) then
               m_err = 'WTAB: Number in IX() is 1(Left), '//
     &                 '2(Center), or 3(Right).'
               ErrCha = ''
               ErrID = 'L:5270/R:wctb/F:a-main1.f'
               goto 999
            end if
            ivx=nint(rnm)
            ic=ic+1

      else if(lum(ic)//lum(ic+1).eq.'iy') then

            ic=ic+2
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'WTAB: Number in IY( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:5283/R:wctb/F:a-main1.f'
               goto 999
            end if
            if(nint(rnm).ne.1.and.nint(rnm).ne.2.and.
     &         nint(rnm).ne.3) then
               m_err = 'WTAB: Number in IY() is 1(Bottom), '//
     &                 '2(Center), or 3(Top).'
               ErrCha = ''
               ErrID = 'L:5291/R:wctb/F:a-main1.f'
               goto 999
            end if
            ivy=nint(rnm)
            ic=ic+1

*-----------------------------------------------------------------------

      else if(lum(ic)//lum(ic+1)//lum(ic+2)//lum(ic+3) .eq.
     &        'tab{' ) then

                  if( dum(ic) .eq. 'T' ) then

                     ivs = 1

                  end if

            icc = 0

            ic = ic + 3

  560       ic = ic + 1

               if( ic .gt. icolm ) then
                  m_err = 'WTAB: Description of tab{} is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:5317/R:wctb/F:a-main1.f'
                  goto 999
               end if

               if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 560


               if( lum(ic) .ne. '|' .and. lum(ic) .ne. 'l' .and.
     &             lum(ic) .ne. 'r' .and. lum(ic) .ne. 'c' .and.
     &             lum(ic) .ne. '}' .and. lum(ic) .ne. '!' ) then

                  m_err = 'WTAB: Description of tab{} is Wrong. '//
     &                    'for example {!l|c|r!}'
                  ErrCha = ''
                  ErrID = 'L:5331/R:wctb/F:a-main1.f'
                  goto 999

               end if


               if( lum(ic) .eq. '}' ) then

                  if( icc .eq. 0 ) then
                     m_err = 'WTAB: Description of tab{} is Wrong. '//
     &                       'for example {!l|c|r!}'
                     ErrCha = ''
                     ErrID = 'L:5343/R:wctb/F:a-main1.f'
                     goto 999
                  end if

                  goto 570

               else if( lum(ic) .eq. '|' ) then

                  do 134 k = 0, nr

                     ittv(k,icc) = abs(ittv(k,icc)) + 1

                     ittl = 1

  134             continue

               else if( lum(ic) .eq. '!' ) then

                  do 135 k = 0, nr

                     ittv(k,icc) = - abs(ittv(k,icc)) - 1

                     ittl = 1

  135             continue

               else if( lum(ic) .eq. 'l' ) then

                  icc = icc + 1

                  do 131 k = 0, nr

                     ittx(k,icc) = 1

  131             continue

               else if( lum(ic) .eq. 'c' ) then

                  icc = icc + 1

                  do 132 k = 0, nr

                     ittx(k,icc) = 2

  132             continue

               else if( lum(ic) .eq. 'r' ) then

                  icc = icc + 1

                  do 133 k = 0, nr

                     ittx(k,icc) = 3

  133             continue

               end if

                  goto 560

  570          continue

               ivc = icc

      end if

            goto 200

*-----------------------------------------------------------------------

  700    continue

*-----------------------------------------------------------------------
*        READ NEW LINE
*-----------------------------------------------------------------------

            icr = 0

  750       ill(jsn)  = ill(jsn) + 1

            icr = icr + 1

            read(jsi,'(10000a1)', iostat = ios ) (dum(ic),ic=1,icolm)
            if( ios .eq. -1 ) goto 350
            call chlow(dum,lum)

*-----------------------------------------------------------------------
*        END OF WTAB: SECTION BY E:
*-----------------------------------------------------------------------

            do 103 i = 1, icolm - 3

               if(lum(i).ne.' '.and.lum(i).ne.tub) goto 113

  103       continue

  113       continue


            if( lum(i)//lum(i+1) .eq. 'e:' ) then

                  ivr = icr - 1

                  goto 351

            end if

*-----------------------------------------------------------------------
*        FINAL CHARACTER
*-----------------------------------------------------------------------

            do 101 i = icolm, 1, -1

               if(lum(i).ne.' '.and.lum(i).ne.tub) then

                  icf = i

                  goto 111

               end if

  101       continue

*-----------------------------------------------------------------------
*        BLANK LINE
*-----------------------------------------------------------------------

               icr = icr - 1

               goto 750

*-----------------------------------------------------------------------
*        INITIAL CHARACTER
*-----------------------------------------------------------------------

  111    continue

            do 104 i = 1, icf

               if(lum(i).ne.' '.and.lum(i).ne.tub) then

                  ici = i

                  goto 105

               end if

  104       continue

  105       continue


*-----------------------------------------------------------------------
*        CONTROL LINE STARTED BY \
*           \hline, \cline{2-3}, \vline{!r|c|l!}, \vspace{ }
*-----------------------------------------------------------------------

            if( lum(ici) .eq. yen ) then

               icr = icr - 1

               i = ici - 1

  106          i = i + 1

                  if(lum(i).eq.' '.or.lum(i).eq.tub) goto 106

                  if(   lum(i)   .eq. yen .and.
     &              (   lum(i+1) .eq. 'h' .or.
     &                  lum(i+1) .eq. 'c' .or.
     &                  lum(i+1) .eq. 'v' ) .and.
     &                ( lum(i+2) .eq. 'l' .and.
     &                  lum(i+3) .eq. 'i' .and.
     &                  lum(i+4) .eq. 'n' .and.
     &                  lum(i+5) .eq. 'e'       )  .or.
     &              (   lum(i+1) .eq. 'v' .and.
     &                  lum(i+2) .eq. 's' .and.
     &                  lum(i+3) .eq. 'p' .and.
     &                  lum(i+4) .eq. 'a' .and.
     &                  lum(i+5) .eq. 'c' .and.
     &                  lum(i+6) .eq. 'e'       )  .or.
     &              (   lum(i+1) .eq. 't' .and.
     &                  lum(i+2) .eq. 'a' .and.
     &                  lum(i+3) .eq. 'b' .and.
     &                  lum(i+4) .eq. 't' .and.
     &                  lum(i+5) .eq. 'o' .and.
     &                  lum(i+6) .eq. 'p' .and.
     &                  lum(i+7) .eq. 's' .and.
     &                  lum(i+8) .eq. 'p'       )  ) then

*-----------------------------------------------------------------------

                     if( dum(i+1) .eq. 'h' .or.
     &                   dum(i+1) .eq. 'H' ) then

*-----------------------------------------------------------------------

                        do 251 k = 0, ivc

                           if( dum(i+1) .eq. 'h' ) then

                              itth(icr,k) = abs(itth(icr,k)) + 1

                              ittl = 1

                           else if( dum(i+1) .eq. 'H' ) then

                              itth(icr,k) = - abs(itth(icr,k)) - 1

                              ittl = 1

                           end if

  251                   continue

                           i = i + 5

*-----------------------------------------------------------------------

                     else if( dum(i+1) .eq. 'C' .or.
     &                        dum(i+1) .eq. 'c' ) then

*-----------------------------------------------------------------------

                        if( dum(i+1) .eq. 'c' ) then

                           icln = 1

                        else if( dum(i+1) .eq. 'C' ) then

                           icln = -1

                        end if


                        if( lum(i+6) .ne. '{' ) then
                           m_err = 'WTAB: Description of '//yen//
     &                             'cline{} is Wrong. for example '//
     &                             yen//'cline{2-3}'
                           ErrCha = ''
                           ErrID = 'L:5583/R:wctb/F:a-main1.f'
                           l_err = ill(jsn)
                           k_err = jsn
                           goto 999
                        end if

                           i = i + 6

                           call qnum(lum,i,icolm,qrn,ierr,'{','-')

                              if( ierr .ne. 0 ) then
                                 m_err = 'WTAB: Description of '//yen//
     &                                   'cline{} is Wrong. '//
     &                                   'for example '//yen//
     &                                   'cline{2-3}'
                                 ErrCha = ''
                                 ErrID = 'L:5599/R:wctb/F:a-main1.f'
                                 l_err = ill(jsn)
                                 k_err = jsn
                                 goto 999
                              end if

                           iinc = nint(qrn)

                              if( iinc .lt. 1 .or.
     &                            iinc .gt. ivc ) then
                                 m_err = 'WTAB: Description of '//yen//
     &                                   'cline{} is Wrong. '//
     &                                   'for example '//yen//
     &                                   'cline{2-3}'
                                 ErrCha = ''
                                 ErrID = 'L:5614/R:wctb/F:a-main1.f'
                                 l_err = ill(jsn)
                                 k_err = jsn
                                 goto 999
                              end if

                           call qnum(lum,i,icolm,qrn,ierr,'-','}')

                              if( ierr .ne. 0 ) then
                                 m_err = 'WTAB: Description of '//yen//
     &                                   'cline{} is Wrong. '//
     &                                   'for example '//yen//
     &                                   'cline{2-3}'
                                 ErrCha = ''
                                 ErrID = 'L:5628/R:wctb/F:a-main1.f'
                                 l_err = ill(jsn)
                                 k_err = jsn
                                 goto 999
                              end if

                           ifnc = nint(qrn)

                              if( iinc .lt. 0 .or.
     &                            ifnc .lt. iinc .or.
     &                            ifnc .gt. ivc ) then
                                 m_err = 'WTAB: Description of '//yen//
     &                                   'cline{} is Wrong. '//
     &                                   'for example '//yen//
     &                                   'cline{2-3}'
                                 ErrCha = ''
                                 ErrID = 'L:5644/R:wctb/F:a-main1.f'
                                 l_err = ill(jsn)
                                 k_err = jsn
                                 goto 999
                              end if

                        do 252 k = iinc, ifnc

                           if( abs( itth(icr,k) ) .le. 1 ) then

                              itth(icr,k) = icln

                              ittl = 1

                           end if

  252                   continue

*-----------------------------------------------------------------------

                     else if( lum(i+1)//lum(i+2) .eq. 'vs' ) then

*-----------------------------------------------------------------------

                        if( lum(i+7) .ne. '{' ) then
                           m_err = 'WTAB: Description of '//yen//
     &                             'vspace{} is Wrong.'
                           ErrCha = ''
                           ErrID = 'L:5672/R:wctb/F:a-main1.f'
                           l_err = ill(jsn)
                           k_err = jsn
                           goto 999
                        end if

                           i = i + 7

                           call pnum(lum,i,icolm,prn,ierr)

                              if( ierr .ne. 0 ) then
                                 m_err = 'WTAB: Description of '//yen//
     &                                   'vspace{} is Wrong. '
                                 ErrCha = ''
                                 ErrID = 'L:5686/R:wctb/F:a-main1.f'
                                 l_err = ill(jsn)
                                 k_err = jsn
                                 goto 999
                              end if

                           vspt(1,icr) = prn

*-----------------------------------------------------------------------

                     else if( lum(i+1)//lum(i+2) .eq. 'ta' ) then

*-----------------------------------------------------------------------

                        if( lum(i+9) .ne. '{' ) then
                           m_err = 'WTAB: Description of '//yen//
     &                             'tabtopsp{} is Wrong.'
                           ErrCha = ''
                           ErrID = 'L:5704/R:wctb/F:a-main1.f'
                           l_err = ill(jsn)
                           k_err = jsn
                           goto 999
                        end if

                           i = i + 9

                           call pnum(lum,i,icolm,prn,ierr)

                              if( ierr .ne. 0 ) then
                                 m_err = 'WTAB: Description of '//yen//
     &                                   'tabtopsp{} is Wrong. '
                                 ErrCha = ''
                                 ErrID = 'L:5718/R:wctb/F:a-main1.f'
                                 l_err = ill(jsn)
                                 k_err = jsn
                                 goto 999
                              end if

                           vspt(2,icr+1) = prn

*-----------------------------------------------------------------------

                     else if( lum(i+1)//lum(i+2) .eq. 'vl' ) then

*-----------------------------------------------------------------------

                        if( lum(i+6) .ne. '{' ) then
                           m_err = 'WTAB: Description of '//yen//
     &                             'vline{} is Wrong. for example '//
     &                             yen//'vline{!r|c|l!}'
                           ErrCha = ''
                           ErrID = 'L:5737/R:wctb/F:a-main1.f'
                           l_err = ill(jsn)
                           k_err = jsn
                           goto 999
                        end if

                           icc = 0

                           i = i + 6

  561                      i = i + 1

                           if( i .gt. icf ) then
                              m_err = 'WTAB: Description of '//yen//
     &                                'vline is Wrong.'
                              ErrCha = ''
                              ErrID = 'L:5753/R:wctb/F:a-main1.f'
                              l_err = ill(jsn)
                              k_err = jsn
                              goto 999
                           end if

                        if( lum(i) .eq. ' ' .or. lum(i) .eq. tub )
     &                  goto 561


                        if( lum(i) .ne. '|' .and.
     &                      lum(i) .ne. '!' .and.
     &                      lum(i) .ne. 'l' .and.
     &                      lum(i) .ne. 'r' .and.
     &                      lum(i) .ne. 'c' .and.
     &                      lum(i) .ne. '}' ) then

                           m_err = 'WTAB: Description of '//yen//
     &                             'vline{} '//
     &                             'is Wrong. for example {!l|c|r!}'
                           ErrCha = ''
                           ErrID = 'L:5774/R:wctb/F:a-main1.f'
                           l_err = ill(jsn)
                           k_err = jsn
                           goto 999

                        end if


                        if( lum(i) .eq. '}' ) then

                           if( icc .ne. ivc ) then

                              m_err = 'WTAB: Description of '//yen//
     &                                'vline{} '//
     &                                'is Wrong. for example '//yen//
     &                                'vline{!l|c|r!}. Column Number'//
     &                                ' should be the same as in tab{ }'
                              ErrCha = ''
                              ErrID = 'L:5792/R:wctb/F:a-main1.f'
                              l_err = ill(jsn)
                              k_err = jsn
                              goto 999

                           end if

                           goto 575

                        else if( lum(i) .eq. '|' ) then

                           if( abs( ittv(icr,icc) ) .le. 1 ) then

                              ittv(icr,icc) = 1

                              ittl = 1

                           end if

                        else if( lum(i) .eq. '!' ) then

                           if( abs( ittv(icr,icc) ) .le. 1 ) then

                              ittv(icr,icc) = -1

                              ittl = 1

                           end if

                        else if( lum(i) .eq. 'l' ) then

                           icc = icc + 1

                              ittx(icr,icc) = 1

                        else if( lum(i) .eq. 'c' ) then

                           icc = icc + 1

                              ittx(icr,icc) = 2

                        else if( lum(i) .eq. 'r' ) then

                           icc = icc + 1

                              ittx(icr,icc) = 3

                        end if

                           goto 561

*-----------------------------------------------------------------------

                     else

                        m_err = 'WTAB: Description of Control line '//
     &                          'is Wrong. for example '//yen//
     &                          'hline, '//yen//
     &                          'cline{2-3}, '//yen//
     &                          'vline{!r|c|l!} '//yen//'vspace{}'
                        ErrCha = ''
                        ErrID = 'L:5853/R:wctb/F:a-main1.f'
                        l_err = ill(jsn)
                        k_err = jsn
                        goto 999

                     end if

*-----------------------------------------------------------------------

                  else

                        icr = icr + 1

                        goto 580

                  end if

*-----------------------------------------------------------------------

  575          continue

               if( i .lt. icf ) goto 106

               goto 750

            end if

*-----------------------------------------------------------------------
*        SEQUENTIAL TEXT LINES
*-----------------------------------------------------------------------

  580    continue

               ict = 0

                  iseql = 0
               if( icf .eq. 1 .and. dum(icf) .eq. yen )
     &            iseql = 1
               if( icf .gt. 1 ) then
                  if( dum(icf) .eq. yen .and. dum(icf-1) .ne. yen .and.
     &                                  ichar(dum(icf-1)) .le. 128 )
     &            iseql = 2
               end if

         if( iseql .gt. 0 ) then

  410          icf = icf - 1

               do 401 i = 1, icf

                  jdum( ict + i ) = dum( i )

  401          continue

                  ict = ict + icf

                  ill(jsn)  = ill(jsn) + 1

                  read(jsi,'(10000a1)',iostat=ios) (dum(ic),ic=1,icolm)
                  if( ios .eq. -1 ) goto 350
                  call chlow(dum,lum)

               do 102 i = icolm, 1, -1

                  if(lum(i).ne.' '.and.lum(i).ne.tub) goto 112

  102          continue
  112          icf = i

                  jseql = 0
               if( icf .eq. 1 .and. dum(icf) .eq. yen )
     &            jseql = 1
               if( icf .gt. 1 ) then
                  if( dum(icf) .eq. yen .and. dum(icf-1) .ne. yen .and.
     &                                  ichar(dum(icf-1)) .le. 128 )
     &            jseql = 2
               end if

               if( jseql .gt. 0 ) goto 410

         end if

*-----------------------------------------------------------------------
*        NORMAL LINE
*-----------------------------------------------------------------------

               do 400 i = 1, icf

                  jdum( ict + i ) = dum( i )

  400          continue

                  ict = ict + icf

*-----------------------------------------------------------------------
*           NEGLECT LAST \\ FOR TEX EXPRESSION
*-----------------------------------------------------------------------

            if( ict .gt. 2 ) then

               if( jdum(ict-1) .eq. yen .and. jdum(ict) .eq. yen .and.
     &                                  ichar(jdum(ict-2)) .le. 128 )
     &            ict = ict - 2

            end if

*-----------------------------------------------------------------------
*        CHECK OF LANGAGE CODE AND
*        PRE-PROCEDURE FOR JAPANESE CODE
*-----------------------------------------------------------------------

                  call jpncode(jdum,ict,ifon)
                  call jpnprep(jdum,ict,ifon)

*-----------------------------------------------------------------------
*        CHECK EACH ELEMENTS
*-----------------------------------------------------------------------

                  icc = 0
                  ibi = 0

                  k = 0

  100             k = k + 1

                  if( k .gt. ict ) goto 110

*-----------------------------------------------------------------------
*           SEL COLOR
*-----------------------------------------------------------------------

                  if( ibi .eq. 0 .and.
     &                jdum(k)   .eq. yen .and.
     &                jdum(k+1) .eq. 's' .and.
     &                jdum(k+2) .eq. 'e' .and.
     &                jdum(k+3) .eq. 'l' .and.
     &                jdum(k+4) .eq. 'c' .and.
     &                jdum(k+5) .eq. 'o' .and.
     &                jdum(k+6) .eq. 'l' .and.
     &                jdum(k+7) .eq. 'o' .and.
     &                jdum(k+8) .eq. 'r' .and.
     &                jdum(k+9) .eq. '{' ) then

                     k = k + 8

                        call chlow(jdum,jlum)

                     call dcols(1,jlum,k,rcol,ierr,1000,'{','}')

                     if(ierr.eq.1) then

                        m_err = yen//'selcolor{ } is wrong number'
                        ErrCha = ''
                        ErrID = 'L:6006/R:wctb/F:a-main1.f'
                        l_err = ill(jsn)
                        k_err = jsn
                        goto 999
                     end if

                     csl(icr,icc+1,1) = rcol(1)
                     csl(icr,icc+1,2) = rcol(2)
                     csl(icr,icc+1,3) = rcol(3)

                  end if

*-----------------------------------------------------------------------

                  if( ibi .eq. 0 .and.
     &              ( jdum(k) .eq. ' '.or.
     &                jdum(k). eq. tub ) ) goto 100

                     if( ibi .eq. 0 ) ibi = k

*-----------------------------------------------------------------------
*              SKIP JAPANESE TWO BITE CHARACTERS AND INITIAL BLANCK
*-----------------------------------------------------------------------

                  if( spjpn(k,ict,icf,jdum,ifon,njpn) ) then

                           k = icf

*-----------------------------------------------------------------------

                  else if( jdum(k) .eq. yen .and. k .lt. ict .and.
     &                     jdum(k+1) .eq. '&' ) then

                           k = k + 1

*-----------------------------------------------------------------------

                  else if( jdum(k) .eq. '&' ) then

                           icc = icc + 1

                     if( k .gt. ibi ) then

                        do 120 l = k-1, ibi, -1

                           if( jdum(l) .ne. ' ' .and.
     &                         jdum(l) .ne. tub ) goto 121

  120                   continue

  121                      ibf = l

                           write(jhs) (jdum(l),l=ibi,ibf)

                           lvm(icr,icc) = ibf - ibi + 1

                     else

                           lvm(icr,icc) = 0

                     end if

                           ibi = 0

                  end if

                           goto 100

*-----------------------------------------------------------------------

  110             continue

                           icc = icc + 1

                           if( icc .ne. ivc ) then

                              m_err = 'WTAB: Number of Column '//
     &                                'is Wrong. Column Number'//
     &                                ' should be the same as in tab{ }'
                              ErrCha = ''
                              ErrID = 'L:6086/R:wctb/F:a-main1.f'
                              l_err = ill(jsn)
                              k_err = jsn
                              goto 999

                           end if

                     if( ibi .gt. 0 .and. k .gt. ibi ) then

                        do 123 l = k-1, ibi, -1

                           if( jdum(l) .ne. ' ' .and.
     &                         jdum(l) .ne. tub ) goto 124

  123                   continue

  124                      ibf = l

                           write(jhs) (jdum(l),l=ibi,ibf)

                           lvm(icr,icc) = ibf - ibi + 1

                     else

                           lvm(icr,icc) = 0

                     end if

*-----------------------------------------------------------------------

               goto 750

*-----------------------------------------------------------------------

  350          iend = 1
  351          continue

*-----------------------------------------------------------------------
*        write information on file
*-----------------------------------------------------------------------

            rewind(jhs)

            write(jwt) ivx,ivy,ivr,ivc,ivf,ivs,ittl

         do i = 1, 13
            write(jwt) vxys(i)
         end do

         do i = 0, ivr
            write(jwt) vspt(1,i), vspt(2,i)
         do j = 0, ivc
            write(jwt) ittx(i,j), ittv(i,j), itth(i,j)
         end do
         end do

         do i = 1, ivr
         do j = 1, ivc
            write(jwt) csl(i,j,1), csl(i,j,2), csl(i,j,3)
            write(jwt) lvm(i,j)
         end do
         end do

         do i = 1, ivr
         do j = 1, ivc

            if( lvm(i,j) .gt. 0 ) then

                read(jhs) (jdum(k),k=1,lvm(i,j))
               write(jwt) (jdum(k),k=1,lvm(i,j))

            end if

         end do
         end do

*-----------------------------------------------------------------------

         close(jhs)

      return

*-----------------------------------------------------------------------

  999          ierr = 1
               return
               end


************************************************************************
*                                                                      *
      subroutine arrcm(dum,lum,in0,ia,axys,iax,ian,ierr,ifon,
     &                 ill,ilf,jsi,jsn,iend,iaf)
*                                                                      *
*      PURPOSE  :   READ ARROW AND COMMENTS                            *
*                                                                      *
*      FORMAT   :  AW: TEXT / X( ) Y( ) AX( ) AY( ) S( ) C( ) A( )     *
*                             T or Z, IR or IL                         *
*                                                                      *
*           IA    :   NUMBER OF COMMENTS                               *
*                                                                      *
*           AXYS(IW,11) :  1-> X-COORDINATE OF ORIGIN    X( )          *
*                          2-> Y-COORDINATE OF ORIGIN    Y( )          *
*                          3-> X-COORDINATE OF ARROW    AX( )          *
*                          4-> Y-COORDINATE OF ARROW    AY( )          *
*                          5-> X-COORDINATE OF COMMENT                 *
*                          6-> SIZE OF COMMENT           S( )          *
*                          7,8,9 -> COLORS OF CHARACTERS C( )          *
*                         10-> WIDTH OF LINE  ( 4 dd default ) T or Z  *
*                         11-> ANGLE OF ARROW            A( )          *
*                                                                      *
*           IAX(IA)     :  X JUSTIFICATION OF TEXT    IR or IL         *
*                         = 1 ; LEFT                                   *
*                           3 ; RIGHT                                  *
*                                                                      *
*           ACM(IA,ICHRL) :  TEXT               A: TEXT /              *
*                                                                      *
*           IAN(IA)     :  LENGTH OF COMMENTS                          *
*                                                                      *
*           IAF(IA)     :  FONT OF COMMENTS                            *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      dimension axys(mc,11),iax(mc),ian(mc),iaf(mc)
      character acm(ichrl)*1

      common /error/ m_err, l_err, k_err
      character m_err*200

      dimension ill(0:9), ilf(0:9)
      dimension rcol(3)

      character c1*1, c2*2

      character yen*1
      character tub*1
      tub = char(9)
      yen = char(92)

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

                  c1 = '('
                  c2 = ')'

                  iax(ia)    = 0

                  ian(ia)    = 0
                  iaf(ia)    = -1

                  axys(ia,1) = r1max
                  axys(ia,2) = r1max
                  axys(ia,3) = r1max
                  axys(ia,4) = r1max
                  axys(ia,6) = 1.0
                  axys(ia,7) = -r1max
                  axys(ia,8) = 1.0
                  axys(ia,9) = 1.0
                  axys(ia,10) = 4.0
                  axys(ia,11) = 1.0

                  in   = in0
                  ict  = 0
                  iend = 0

*-----------------------------------------------------------------------

            do 101 i = icolm, in, -1

               if(lum(i).ne.' '.and.lum(i).ne.tub) goto 111

  101       continue

               return

  111       continue

*----------------------------------------------------------------------*

                  iseql = 0
               if( i .eq. 1 .and. dum(i) .eq. yen )
     &            iseql = 1
               if( i .gt. 1 ) then
                  if( dum(i) .eq. yen .and. dum(i-1) .ne. yen .and.
     &                                ichar(dum(i-1)) .le. 128 )
     &            iseql = 2
               end if

         if( iseql .gt. 0 ) then

  410          icf = i - 1

               do 401 i = 1, icf - in + 1

                  acm( ict + i ) = dum( in + i - 1 )

  401          continue

                  ict = ict + icf - in + 1

                  ill(jsn)  = ill(jsn) + 1
                  in  = 1

                  read(jsi,'(10000a1)',iostat=ios) (dum(ic),ic=1,icolm)
                  if( ios .eq. -1 ) goto 350
                  call chlow(dum,lum)

               do 102 i = icolm, 1, -1

                  if(dum(i).ne.' '.and.dum(i).ne.tub) goto 112

  102          continue
  112          continue

                  jseql = 0
               if( icf .eq. 1 .and. dum(i) .eq. yen )
     &            jseql = 1
               if( icf .gt. 1 ) then
                  if( dum(i) .eq. yen .and. dum(i-1) .ne. yen .and.
     &                                ichar(dum(i-1)) .le. 128 )
     &            jseql = 2
               end if

               if( jseql .gt. 0 ) goto 410

         end if

*----------------------------------------------------------------------*

               goto 351
  350          continue

                  iend = 1

  351          continue

               icms = 0

            do 100 i = icolm, in, -1

               if( lum(i) .eq. ']' ) icms = 1
               if( lum(i) .eq. '[' ) icms = 0

               if( icms .eq. 0 .and. lum(i) .eq. '/' ) goto 110

  100          continue

               m_err = 'AW: / is Missing Between Text and Parameters.'
               ErrCha = ''
               ErrID = 'L:6350/R:arrcm/F:a-main1.f'
               goto 999

  110       icf = i - 1

            if( icf .eq. in - 1 ) then

               if( ict .eq. 0 ) then

                  m_err = 'W: There is No Text Before /.'
                  ErrCha = ''
                  ErrID = 'L:6361/R:arrcm/F:a-main1.f'
                  goto 999

               else

                  ian(ia) = ict

               end if

            else

               do 400 i = 1, icf - in + 1

                  acm( ict + i ) = dum( in + i - 1 )

  400          continue

                  ian(ia) = ict + icf - in + 1

            end if

*-----------------------------------------------------------------------

                call jpncode(acm,ian(ia),ifon)
                call jpnprep(acm,ian(ia),ifon)

               write(jha) (acm(i),i=1,ian(ia))

*-----------------------------------------------------------------------

      ic = icf + 1

  200 ic = ic + 1

         if(ic.gt.icolm) goto 500

         if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 200

      if(lum(ic).ne.'x'.and.lum(ic).ne.'y'.and.lum(ic).ne.'i'.and.
     &   lum(ic).ne.'s'.and.lum(ic).ne.'f'.and.lum(ic).ne.'c'.and.
     &   lum(ic).ne.'t'.and.lum(ic).ne.'z'.and.lum(ic).ne.'a') then

         m_err = 'AW: Text/ X() Y() AX() AY() S() F() A()'//
     &           ' T Z IR IL; Unexpected Parameter.'
         ErrCha = ''
         ErrID = 'L:6406/R:arrcm/F:a-main1.f'

         goto 999

      end if


      if(lum(ic).eq.'x') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'AW: Number in X( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:6420/R:arrcm/F:a-main1.f'
               goto 999
            end if
            axys(ia,1)=rnm

      else if(lum(ic).eq.'y') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'AW: Number in Y( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:6432/R:arrcm/F:a-main1.f'
               goto 999
            end if
            axys(ia,2)=rnm

      else if(lum(ic)//lum(ic+1) .eq. 'ax' ) then

            ic=ic+2
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'AW: Number in AX( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:6444/R:arrcm/F:a-main1.f'
               goto 999
            end if
            axys(ia,3)=rnm

      else if(lum(ic)//lum(ic+1) .eq. 'ay' ) then

            ic=ic+2
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'AW: Number in AY( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:6456/R:arrcm/F:a-main1.f'
               goto 999
            end if
            axys(ia,4)=rnm

      else if(lum(ic)//lum(ic+1) .eq. 'ir' ) then

            iax(ia) = 3
            ic=ic+1

      else if(lum(ic)//lum(ic+1) .eq. 'il' ) then

            iax(ia) = 1
            ic=ic+1

      else if(lum(ic).eq.'s') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'AW: Number in S( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:6478/R:arrcm/F:a-main1.f'
               goto 999
            end if
            axys(ia,6)=rnm

      else if(lum(ic).eq.'f') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'AW: Number in F( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:6490/R:arrcm/F:a-main1.f'
               goto 999
            end if

            iaf(ia)=nint(rnm)

            if( iaf(ia) .lt. 0 .or. iaf(ia) .gt. 12 ) iaf(ia) = -1


      else if(lum(ic).eq.'c') then

            call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               axys(ia,7) = rcol(1)
               axys(ia,8) = rcol(2)
               axys(ia,9) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'AW: Number in C( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:6510/R:arrcm/F:a-main1.f'
               goto 999
            end if

      else if(lum(ic).eq.'t') then

            axys(ia,10) = axys(ia,10) + 3.0

      else if(lum(ic).eq.'z') then

            axys(ia,10) = axys(ia,10) - 1.0

            if( axys(ia,10) .le. 0.0 ) axys(ia,10) = 1.0

      else if(lum(ic).eq.'a') then

            ic=ic+1
            call pnum(lum,ic,icolm,rnm,ierr)
            if(ierr.ne.0) then
               m_err = 'AW: Number in A( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:6531/R:arrcm/F:a-main1.f'
               goto 999
            end if
            axys(ia,11)=rnm

      end if


      goto 200


  500 continue

      do 300 i = 1, 4
        if(axys(ia,i).gt.1.0d+41) then
           m_err = 'AW: You need X() Y() AX() AY() Four '//
     &             'Position at Least.'
           ErrCha = ''
           ErrID = 'L:6549/R:arrcm/F:a-main1.f'
           goto 999
        end if
  300 continue

      return
  999 ierr=1
      l_err = ill(jsn)
      k_err = jsn
      return
      end

************************************************************************
*                                                                      *
      subroutine arrow(dum,lum,in,ill,ilf,jsn,ib,bxys,ierr)
*                                                                      *
*      PURPOSE  :   READ ARROW ONLY                                    *
*                                                                      *
*      FORMAT   :  A: X( ) Y( ) AX( ) AY( ) C( ) A( ) N, T or Z        *
*                                                                      *
*                  N: NOLINE                                           *
*                  Z: THIN LINE WIDTH                                  *
*                  T: THICK LINE WIDTH                                 *
*               A( ): ANGLE OF THE ARROW RELATIVE TO DEFAULT           *
*                                                                      *
*           BXYS(IB,10  :  1-> X-COORDINATE OF ORIGIN    X( )          *
*                          2-> Y-COORDINATE OF ORIGIN    Y( )          *
*                          3-> X-COORDINATE OF ARROW    AX( )          *
*                          4-> Y-COORDINATE OF ARROW    AY( )          *
*                          5,6,7 -> COLOR OF ARROW       C( )          *
*                          8-> WIDTH OF LINE  ( 4 dd default ) T or Z  *
*                          9-> -1 ; NO LINE              N             *
*                         10-> ANGLE OF ARROW            A( )          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      dimension bxys(mc,10)

      common /error/ m_err, l_err, k_err
      character m_err*200

      dimension ill(0:9), ilf(0:9)
      dimension rcol(3)

      character c1*1, c2*1

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

               c1 ='('
               c2 =')'

               bxys(ib, 1) = r1max
               bxys(ib, 2) = r1max
               bxys(ib, 3) = r1max
               bxys(ib, 4) = r1max
               bxys(ib, 5) = -r1max
               bxys(ib, 6) = 1.0
               bxys(ib, 7) = 1.0
               bxys(ib, 8) = 4.0
               bxys(ib, 9) = 1.0
               bxys(ib,10) = 1.0

*-----------------------------------------------------------------------

      ic = in - 1

  200 ic = ic + 1

         if(ic.gt.icolm) goto 500
         if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 200

      if(lum(ic).ne.'x'.and.lum(ic).ne.'y'.and.lum(ic).ne.'c'.and.
     &   lum(ic).ne.'n'.and.lum(ic).ne.'t'.and.lum(ic).ne.'z'.and.
     &   lum(ic).ne.'a') then
         m_err = 'A: X() Y() AX() AY() C() A() T Z N; '//
     &           'Unexpected Parameter.'
         ErrCha = ''
         ErrID = 'L:6640/R:arrow/F:a-main1.f'
         goto 999
      end if

      if(lum(ic).eq.'x') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'A: Number in X( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:6651/R:arrow/F:a-main1.f'
            goto 999
         end if
         bxys(ib,1)=rnm


      else if(lum(ic).eq.'y') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'A: Number in Y( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:6664/R:arrow/F:a-main1.f'
            goto 999
         end if
         bxys(ib,2)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'ax' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'A: Number in AX( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:6677/R:arrow/F:a-main1.f'
            goto 999
         end if
         bxys(ib,3)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'ay' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'A: Number in AY( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:6690/R:arrow/F:a-main1.f'
            goto 999
         end if
         bxys(ib,4)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'cb' ) then

         ic=ic+1
           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

            if(ierr.eq.1) then
               m_err = 'A: Number in CB( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:6704/R:arrow/F:a-main1.f'
               goto 999
            end if


      else if(lum(ic).eq.'c') then

           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               bxys(ib,5) = rcol(1)
               bxys(ib,6) = rcol(2)
               bxys(ib,7) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'A: Number in C( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:6720/R:arrow/F:a-main1.f'
               goto 999
            end if


      else if(lum(ic).eq.'t') then

               bxys(ib,8) = bxys(ib,8) + 3.0


      else if(lum(ic).eq.'z') then

               bxys(ib,8) = bxys(ib,8) - 1.0

               if( bxys(ib,8) .le. 0.0 ) bxys(ib,8) = 1.0


      else if(lum(ic).eq.'n') then

               bxys(ib,9) = - 1.0

      else if(lum(ic).eq.'a') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'A: Number in A( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:6748/R:arrow/F:a-main1.f'
            goto 999
         end if
         bxys(ib,10)=rnm


      end if

*----------------------------------------------------------------------*

      goto 200

  500 continue

      do 300 i=1,4
        if(bxys(ib,i).gt.1.0d+41) then
           m_err = 'A: You need X() Y() AX() AY() Four '//
     &             'Position at Least.'
           ErrCha = ''
           ErrID = 'L:6767/R:arrow/F:a-main1.f'
           goto 999
        end if
  300 continue

      return
  999 ierr=1
      l_err = ill(jsn)
      k_err = jsn
      return
      end

************************************************************************
*                                                                      *
      subroutine arrowc(dum,lum,in,ill,ilf,jsn,ib,cxys,ierr)
*                                                                      *
*      PURPOSE  :   READ ARROW ONLY                                    *
*                                                                      *
*      FORMAT   :  AB: X( ) Y( ) AX( ) AY( ) C( ) CB( ) A( )           *
*                      N, T or Z                                       *
*                                                                      *
*                  N: NOLINE                                           *
*                  Z: THIN LINE WIDTH                                  *
*                  T: THICK LINE WIDTH                                 *
*               A( ): ANGLE OF THE ARROW RELATIVE TO DEFAULT           *
*                                                                      *
*           CXYS(IB,11) :  1-> X-COORDINATE OF ORIGIN    X( )          *
*                          2-> Y-COORDINATE OF ORIGIN    Y( )          *
*                          3-> X-COORDINATE OF ARROW    AX( )          *
*                          4-> Y-COORDINATE OF ARROW    AY( )          *
*                          5,6,7 -> COLOR OF ARROW LINE  C( )          *
*                          8-> WIDTH OF LINE  ( 4 dd default ) T or Z  *
*                          9-> -1 ; NO LINE              N             *
*                         10-> ANGLE OF ARROW            A( )          *
*                         11,12,13 -> COLOR OF ARROW INTERIA  CB( )    *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      dimension cxys(mc,13)

      common /error/ m_err, l_err, k_err
      character m_err*200

      dimension ill(0:9), ilf(0:9)
      dimension rcol(3)

      character c1*1, c2*1

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

               c1 ='('
               c2 =')'

               cxys(ib, 1) = r1max
               cxys(ib, 2) = r1max
               cxys(ib, 3) = r1max
               cxys(ib, 4) = r1max
               cxys(ib, 5) = -r1max
               cxys(ib, 6) = 1.0
               cxys(ib, 7) = 1.0
               cxys(ib, 8) = 8.0
               cxys(ib, 9) = 1.0
               cxys(ib,10) = 1.0
               cxys(ib,11) = -r1max
               cxys(ib,12) = 1.0
               cxys(ib,13) = 1.0

*-----------------------------------------------------------------------

      ic = in - 1

  200 ic = ic + 1

         if(ic.gt.icolm) goto 500
         if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 200

      if(lum(ic).ne.'x'.and.lum(ic).ne.'y'.and.lum(ic).ne.'c'.and.
     &   lum(ic).ne.'n'.and.lum(ic).ne.'t'.and.lum(ic).ne.'z'.and.
     &   lum(ic).ne.'a') then
         m_err = 'AB: X() Y() AX() AY() C() CB() A() T Z N; '//
     &           'Unexpected Parameter.'
         ErrCha = ''
         ErrID = 'L:6863/R:arrowc/F:a-main1.f'
         goto 999
      end if

      if(lum(ic).eq.'x') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'AB: Number in X( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:6874/R:arrowc/F:a-main1.f'
            goto 999
         end if
         cxys(ib,1)=rnm


      else if(lum(ic).eq.'y') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'AB: Number in Y( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:6887/R:arrowc/F:a-main1.f'
            goto 999
         end if
         cxys(ib,2)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'ax' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'AB: Number in AX( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:6900/R:arrowc/F:a-main1.f'
            goto 999
         end if
         cxys(ib,3)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'ay' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'AB: Number in AY( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:6913/R:arrowc/F:a-main1.f'
            goto 999
         end if
         cxys(ib,4)=rnm


      else if(lum(ic).eq.'t') then

               cxys(ib,8) = cxys(ib,8) + 4.0


      else if(lum(ic).eq.'z') then

               cxys(ib,8) = cxys(ib,8) - 2.0

               if( cxys(ib,8) .le. 2.0 ) cxys(ib,8) = 2.0


      else if(lum(ic).eq.'n') then

               cxys(ib,9) = - 1.0

      else if(lum(ic).eq.'a') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'AB: Number in A( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:6942/R:arrowc/F:a-main1.f'
            goto 999
         end if
         cxys(ib,10)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'cb' ) then

         ic=ic+1
           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               cxys(ib,11) = rcol(1)
               cxys(ib,12) = rcol(2)
               cxys(ib,13) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'AB: Number in CB( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:6960/R:arrowc/F:a-main1.f'
               goto 999
            end if

      else if(lum(ic).eq.'c') then

           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               cxys(ib,5) = rcol(1)
               cxys(ib,6) = rcol(2)
               cxys(ib,7) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'AB: Number in C( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:6975/R:arrowc/F:a-main1.f'
               goto 999
            end if


      end if

*----------------------------------------------------------------------*

      goto 200

  500 continue

      do 300 i=1,4
        if(cxys(ib,i).gt.1.0d+41) then
           m_err = 'AB: You need X() Y() AX() AY() Four '//
     &             'Position at Least.'
           ErrCha = ''
           ErrID = 'L:6993/R:arrowc/F:a-main1.f'
           goto 999
        end if
  300 continue

      return
  999 ierr=1
      l_err = ill(jsn)
      k_err = jsn
      return
      end

************************************************************************
*                                                                      *
      subroutine polg(dum,lum,in,ill,ilf,jsn,ib,pxys,ierr)
*                                                                      *
*      PURPOSE  :   READ POLYGON                                       *
*                                                                      *
*      FORMAT   :  POLG: X( ) Y( ) S( ) SX( ) SY( ) CL( ) CB( ) A( )   *
*                        PL( )                                         *
*                                                                      *
*           PXYS(IB,14) :  1-> X-COORDINATE OF ORIGIN    X( )          *
*                          2-> Y-COORDINATE OF ORIGIN    Y( )          *
*                          3-> NUMBER OF CORNAER        PL( )          *
*                          4-> SCAL OF POLYGON           S( )          *
*                          5-> X SCAL OF POLYGON        SX( )          *
*                          6-> Y SCAL OF POLYGON        SY( )          *
*                          7-> ANGLE OF POLYGON          A( )          *
*                          8,9.10 -> COLOR OF POLYGON LINE  CL( )      *
*                         11,12,13 -> COLOR OF POLYGON INTERIA CB( )   *
*                         14-> WIDTH OF LINE  ( 4 dd default ) T or Z  *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      dimension pxys(mc,14)

      common /error/ m_err, l_err, k_err
      character m_err*200

      dimension ill(0:9), ilf(0:9)
      dimension rcol(3)

      character c1*1, c2*1

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

               c1 ='('
               c2 =')'

               pxys(ib, 1) = -r1max
               pxys(ib, 2) = -r1max
               pxys(ib, 3) = 5.0
               pxys(ib, 4) = 1.0
               pxys(ib, 5) = 1.0
               pxys(ib, 6) = 1.0
               pxys(ib, 7) = -r1max
               pxys(ib, 8) = -r1max
               pxys(ib, 9) = 1.0
               pxys(ib,10) = 1.0
               pxys(ib,11) = -r1max
               pxys(ib,12) = 1.0
               pxys(ib,13) = 1.0
               pxys(ib,14) = 4.0

*-----------------------------------------------------------------------

      ic = in - 1

  200 ic = ic + 1

         if(ic.gt.icolm) goto 500
         if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 200

      if(lum(ic).ne.'x'.and.lum(ic).ne.'y'.and.lum(ic).ne.'c'.and.
     &   lum(ic).ne.'s'.and.lum(ic).ne.'p'.and.lum(ic).ne.'t'.and.
     &   lum(ic).ne.'z'.and.lum(ic).ne.'a') then
         m_err = 'POLG: X() Y() PL() S() SX() SY() A() '//
     &           'CL() CB() T Z ; Unexpected Parameter.'
         ErrCha = ''
         ErrID = 'L:7086/R:polg/F:a-main1.f'
         goto 999
      end if

      if(lum(ic).eq.'x') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'POLG: Number in X( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7097/R:polg/F:a-main1.f'
            goto 999
         end if
         pxys(ib,1)=rnm


      else if(lum(ic).eq.'y') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'POLG: Number in Y( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7110/R:polg/F:a-main1.f'
            goto 999
         end if
         pxys(ib,2)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'pl' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'POLG: Number in PL( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7123/R:polg/F:a-main1.f'
            goto 999
         end if
         pxys(ib,3)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 's(' ) then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'POLG: Number in S( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7136/R:polg/F:a-main1.f'
            goto 999
         end if
         pxys(ib,4)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'sx' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'POLG: Number in SX( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7149/R:polg/F:a-main1.f'
            goto 999
         end if
         pxys(ib,5)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'sy' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'POLG: Number in SY( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7162/R:polg/F:a-main1.f'
            goto 999
         end if
         pxys(ib,6)=rnm


      else if(lum(ic).eq.'a') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'POLG: Number in A( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7175/R:polg/F:a-main1.f'
            goto 999
         end if
         pxys(ib,7)=rnm


      else if(lum(ic).eq.'t') then

               pxys(ib,14) = pxys(ib,14) + 3.0


      else if(lum(ic).eq.'z') then

               pxys(ib,14) = pxys(ib,14) - 1.0

               if( pxys(ib,14) .le. 1.0 ) pxys(ib,14) = 1.0


      else if(lum(ic)//lum(ic+1) .eq. 'cl' ) then

         ic=ic+1
           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               pxys(ib, 8) = rcol(1)
               pxys(ib, 9) = rcol(2)
               pxys(ib,10) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'POLG: Number in CL( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:7205/R:polg/F:a-main1.f'
               goto 999
            end if


      else if(lum(ic)//lum(ic+1) .eq. 'cb' ) then

         ic=ic+1
           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               pxys(ib,11) = rcol(1)
               pxys(ib,12) = rcol(2)
               pxys(ib,13) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'POLG: Number in CB( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:7222/R:polg/F:a-main1.f'
               goto 999
            end if

      end if

*----------------------------------------------------------------------*

      goto 200

  500 continue

      do 300 i=1,3
        if(pxys(ib,i).lt.-r0max) then
           m_err = 'POLG: You need X() Y() '//
     &             'Position of Polygon at Least.'
           ErrCha = ''
           ErrID = 'L:7239/R:polg/F:a-main1.f'
           goto 999
        end if
  300 continue

      return
  999 ierr=1
      l_err = ill(jsn)
      k_err = jsn
      return
      end

************************************************************************
*                                                                      *
      subroutine rbox(dum,lum,in,ill,ilf,jsn,ib,pxys,ierr)
*                                                                      *
*      PURPOSE  :  READ BOX                                            *
*                                                                      *
*      FORMAT   :  BOX: X( ) Y( ) S( ) SX( ) SY( ) CL( ) CB( ) A( )    *
*                       BOX(singlebox) CS( )                           *
*                                                                      *
*           PXYS(IB,16) :  1-> X-COORDINATE OF ORIGIN    X( )          *
*                          2-> Y-COORDINATE OF ORIGIN    Y( )          *
*                          3-> SCAL OF BOX               S( )          *
*                          4-> X SCAL OF BOX            SX( )          *
*                          5-> Y SCAL OF BOX            SY( )          *
*                          6-> ANGLE OF BOX              A( )          *
*                          7-> TYPE OF BOX singlebox, ......           *
*                          8.9.10-> COLOR OF BOX SHADOW CS( )          *
*                        11,12,13-> COLOR OF BOX LINE   CL( )          *
*                        14,15,16-> COLOR OF BOX INTERIA CB( )         *
*                                                                      *
*           PXYS(IB,7)  =  10 ; singlebox                              *
*                          11 ; Singlebox                              *
*                          12 ; singleBox                              *
*                          13 ; SingleBox                              *
*                          20 ; ovalbox                                *
*                          21 ; Ovalbox                                *
*                          22 ; ovalBox                                *
*                          23 ; OvalBox                                *
*                          30 ; doublebox                              *
*                          31 ; Doublebox                              *
*                          32 ; doubleBox                              *
*                          33 ; DoubleBox                              *
*                          40 ; shadowbox                              *
*                          41 ; Shadowbox                              *
*                          42 ; shadowBox                              *
*                          43 ; ShadowBox                              *
*                          50 ; ovalshadowbox                          *
*                          51 ; Ovalshadowbox                          *
*                          52 ; ovalshadowBox                          *
*                          53 ; OvalshadowBox                          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      dimension pxys(mc,16)

      character rum*2048
*         ( ichrl = 2048 )

      common /error/ m_err, l_err, k_err
      character m_err*200

      dimension ill(0:9), ilf(0:9)
      dimension rcol(3)

      character c1*1, c2*1

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

               c1 ='('
               c2 =')'

               pxys(ib, 1)  = -r1max
               pxys(ib, 2)  = -r1max
               pxys(ib, 3)  = 1.0
               pxys(ib, 4)  = 1.0
               pxys(ib, 5)  = 1.0
               pxys(ib, 6)  = -r1max
               pxys(ib, 7)  = 10.0
               pxys(ib, 8)  = -r1max
               pxys(ib, 9)  = 1.0
               pxys(ib,10)  = 1.0
               pxys(ib,11)  = -r1max
               pxys(ib,12)  = 1.0
               pxys(ib,13)  = 1.0
               pxys(ib,14)  = -r1max
               pxys(ib,15)  = 1.0
               pxys(ib,16)  = 1.0

*-----------------------------------------------------------------------

      ic = in - 1

  200 ic = ic + 1

         if(ic.gt.icolm) goto 600
         if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 200


      if(lum(ic).ne.'x'.and.lum(ic).ne.'y'.and.lum(ic).ne.'c'.and.
     &   lum(ic).ne.'s'.and.lum(ic).ne.'b'.and.lum(ic).ne.'a') then
         m_err = 'BOX: BOX( ) X() Y() S() SX() SY() A() '//
     &           'CL() CB() CS() ; Unexpected Parameter.'
         ErrCha = ''
         ErrID = 'L:7358/R:rbox/F:a-main1.f'
         goto 999
      end if

      if(lum(ic).eq.'x') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'BOX: Number in X( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7369/R:rbox/F:a-main1.f'
            goto 999
         end if
         pxys(ib,1)=rnm


      else if(lum(ic).eq.'y') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'BOX: Number in Y( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7382/R:rbox/F:a-main1.f'
            goto 999
         end if
         pxys(ib,2)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 's(' ) then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'BOX: Number in S( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7395/R:rbox/F:a-main1.f'
            goto 999
         end if
         pxys(ib,3)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'sx' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'BOX: Number in SX( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7408/R:rbox/F:a-main1.f'
            goto 999
         end if
         pxys(ib,4)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'sy' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'BOX: Number in SY( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7421/R:rbox/F:a-main1.f'
            goto 999
         end if
         pxys(ib,5)=rnm


      else if(lum(ic).eq.'a') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'BOX: Number in A( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7434/R:rbox/F:a-main1.f'
            goto 999
         end if
         pxys(ib,6)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'cl' ) then

         ic=ic+1
           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               pxys(ib,11) = rcol(1)
               pxys(ib,12) = rcol(2)
               pxys(ib,13) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'BOX: Number in CL( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:7452/R:rbox/F:a-main1.f'
               goto 999
            end if


      else if(lum(ic)//lum(ic+1) .eq. 'cb' ) then

         ic=ic+1
           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               pxys(ib,14) = rcol(1)
               pxys(ib,15) = rcol(2)
               pxys(ib,16) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'BOX: Number in CB( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:7469/R:rbox/F:a-main1.f'
               goto 999
            end if


      else if(lum(ic)//lum(ic+1) .eq. 'cs' ) then

         ic=ic+1
           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               pxys(ib, 8) = rcol(1)
               pxys(ib, 9) = rcol(2)
               pxys(ib,10) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'BOX: Number in CS( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:7486/R:rbox/F:a-main1.f'
               goto 999
            end if

*-----------------------------------------------------------------------

      else if(lum(ic)//lum(ic+1)//lum(ic+2) .eq. 'box' ) then

               ic = ic + 2

               if(lum(ic+1) .ne. '(') then
                  m_err = 'BOX: Description of BOX( ) is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:7499/R:rbox/F:a-main1.f'
                  goto 999
               end if

               ic = ic + 1

  550          ic = ic + 1

               if( ic .gt. icolm ) then
                  m_err = 'BOX: Description of BOX( ) is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:7510/R:rbox/F:a-main1.f'
                  goto 999
               end if

               if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 550


               do 500 k = ic, icolm

                  kk = k - ic + 1

                  rum(kk:kk) = dum(k)

  500          continue

         if( ( rum(1: 9) .eq. 'singlebox' ) .or.
     &       ( rum(1: 9) .eq. 'Singlebox' ) .or.
     &       ( rum(1: 9) .eq. 'singleBox' ) .or.
     &       ( rum(1: 9) .eq. 'SingleBox' ) .or.
     &       ( rum(1: 7) .eq. 'ovalbox'   ) .or.
     &       ( rum(1: 7) .eq. 'Ovalbox'   ) .or.
     &       ( rum(1: 7) .eq. 'ovalBox'   ) .or.
     &       ( rum(1: 7) .eq. 'OvalBox'   ) .or.
     &       ( rum(1: 9) .eq. 'doublebox' ) .or.
     &       ( rum(1: 9) .eq. 'Doublebox' ) .or.
     &       ( rum(1: 9) .eq. 'doubleBox' ) .or.
     &       ( rum(1: 9) .eq. 'DoubleBox' ) .or.
     &       ( rum(1: 9) .eq. 'shadowbox' ) .or.
     &       ( rum(1: 9) .eq. 'Shadowbox' ) .or.
     &       ( rum(1: 9) .eq. 'shadowBox' ) .or.
     &       ( rum(1: 9) .eq. 'ShadowBox' ) .or.
     &       ( rum(1:13) .eq. 'ovalshadowbox' ) .or.
     &       ( rum(1:13) .eq. 'Ovalshadowbox' ) .or.
     &       ( rum(1:13) .eq. 'ovalshadowBox' ) .or.
     &       ( rum(1:13) .eq. 'OvalshadowBox' )
     &     ) then

               if( rum(1:9) .eq. 'singlebox' ) then

                     iboxd = 10
                     iboxl = 9

               else if( rum(1:9) .eq. 'Singlebox' ) then

                     iboxd = 11
                     iboxl = 9

               else if( rum(1:9) .eq. 'singleBox' ) then

                     iboxd = 12
                     iboxl = 9

               else if( rum(1:9) .eq. 'SingleBox' ) then

                     iboxd = 13
                     iboxl = 9

               else if( rum(1:7) .eq. 'ovalbox' ) then

                     iboxd = 20
                     iboxl = 7

               else if( rum(1:7) .eq. 'Ovalbox' ) then

                     iboxd = 21
                     iboxl = 7

               else if( rum(1:7) .eq. 'ovalBox' ) then

                     iboxd = 22
                     iboxl = 7

               else if( rum(1:7) .eq. 'OvalBox' ) then

                     iboxd = 23
                     iboxl = 7

               else if( rum(1:9) .eq. 'doublebox' ) then

                     iboxd = 30
                     iboxl = 9

               else if( rum(1:9) .eq. 'Doublebox' ) then

                     iboxd = 31
                     iboxl = 9

               else if( rum(1:9) .eq. 'doubleBox' ) then

                     iboxd = 32
                     iboxl = 9

               else if( rum(1:9) .eq. 'DoubleBox' ) then

                     iboxd = 33
                     iboxl = 9

               else if( rum(1:9) .eq. 'shadowbox' ) then

                     iboxd = 40
                     iboxl = 9

               else if( rum(1:9) .eq. 'Shadowbox' ) then

                     iboxd = 41
                     iboxl = 9

               else if( rum(1:9) .eq. 'shadowBox' ) then

                     iboxd = 42
                     iboxl = 9

               else if( rum(1:9) .eq. 'ShadowBox' ) then

                     iboxd = 43
                     iboxl = 9

               else if( rum(1:13) .eq. 'ovalshadowbox' ) then

                     iboxd = 50
                     iboxl = 13

               else if( rum(1:13) .eq. 'Ovalshadowbox' ) then

                     iboxd = 51
                     iboxl = 13

               else if( rum(1:13) .eq. 'ovalshadowBox' ) then

                     iboxd = 52
                     iboxl = 13

               else if( rum(1:13) .eq. 'OvalshadowBox' ) then

                     iboxd = 53
                     iboxl = 13

               end if

         else

                  m_err = 'BOX: Description of BOX( ) is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:7653/R:rbox/F:a-main1.f'
                  goto 999

         end if

               ic = ic + iboxl - 1

  560          ic = ic + 1

               if( ic .gt. icolm ) then
                  m_err = 'BOX: Description of BOX( ) is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:7665/R:rbox/F:a-main1.f'
                  goto 999
               end if

               if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 560

               if( lum(ic) .ne. ')' ) then
                  m_err = 'BOX: Description of BOX( ) is Wrong.'
                  ErrCha = ''
                  ErrID = 'L:7674/R:rbox/F:a-main1.f'
                  goto 999
               end if


               pxys(ib,7) = dble(iboxd)

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

      goto 200

  600 continue

      do 300 i=1,2
        if(pxys(ib,i).lt.-r0max) then
           m_err = 'BOX: You need X() Y() '//
     &             'Position of BOX at Least.'
           ErrCha = ''
           ErrID = 'L:7696/R:rbox/F:a-main1.f'
           goto 999
        end if
  300 continue

      return
  999 ierr=1
      l_err = ill(jsn)
      k_err = jsn
      return
      end

************************************************************************
*                                                                      *
      subroutine ribn(dum,lum,in,ill,ilf,jsn,ib,pxys,ierr)
*                                                                      *
*      PURPOSE  :   READ RIBBON                                        *
*                                                                      *
*      FORMAT   :  RIBN: X( ) Y( ) S( ) SX( ) SY( ) CL( ) CB( ) A( )   *
*                        CS( ) T or Z                                  *
*                                                                      *
*           PXYS(IB,10) :  1-> X-COORDINATE OF ORIGIN    X( )          *
*                          2-> Y-COORDINATE OF ORIGIN    Y( )          *
*                          3-> SCAL OF RIBBON           S( )           *
*                          4-> X SCAL OF RIBBON        SX( )           *
*                          5-> Y SCAL OF RIBBON        SY( )           *
*                          6-> ANGLE OF RIBBON          A( )           *
*                          7-> WIDTH OF LINE  ( 4 dd default ) T or Z  *
*                          8,9,10-> COLOR OF RIBBON SHADOW  CS( )      *
*                        11,12,13-> COLOR OF RIBBON LINE    CL( )      *
*                        14,15,16-> COLOR OF RIBBON INTERIA CB( )      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      dimension pxys(mc,16)

      common /error/ m_err, l_err, k_err
      character m_err*200

      dimension ill(0:9), ilf(0:9)
      dimension rcol(3)

      character c1*1, c2*1

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

               c1 ='('
               c2 =')'

               pxys(ib, 1)  = -r1max
               pxys(ib, 2)  = -r1max
               pxys(ib, 3)  = 1.0
               pxys(ib, 4)  = 1.0
               pxys(ib, 5)  = 1.0
               pxys(ib, 6)  = -r1max
               pxys(ib, 7)  = 4.0
               pxys(ib, 8)  = -r1max
               pxys(ib, 9)  = 1.0
               pxys(ib,10)  = 1.0
               pxys(ib,11)  = -r1max
               pxys(ib,12)  = 1.0
               pxys(ib,13)  = 1.0
               pxys(ib,14)  = -r1max
               pxys(ib,15)  = 1.0
               pxys(ib,16)  = 1.0

*-----------------------------------------------------------------------

      ic = in - 1

  200 ic = ic + 1

         if(ic.gt.icolm) goto 500
         if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 200

      if(lum(ic).ne.'x'.and.lum(ic).ne.'y'.and.lum(ic).ne.'c'.and.
     &   lum(ic).ne.'s'.and.lum(ic).ne.'t'.and.lum(ic).ne.'z'.and.
     &   lum(ic).ne.'a') then
         m_err = 'RIBN: X() Y() S() SX() SY() A() '//
     &           'CL() CB() CS() T Z ; Unexpected Parameter.'
         ErrCha = ''
         ErrID = 'L:7791/R:ribn/F:a-main1.f'
         goto 999
      end if

      if(lum(ic).eq.'x') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'RIBN: Number in X( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7802/R:ribn/F:a-main1.f'
            goto 999
         end if
         pxys(ib,1)=rnm


      else if(lum(ic).eq.'y') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'RIBN: Number in Y( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7815/R:ribn/F:a-main1.f'
            goto 999
         end if
         pxys(ib,2)=rnm


cKN 2024/01/24
c     else if(lum(ic)//lum(ic+1) .eq. 's(' ) then
      else if(lum(ic) .eq. 's' ) then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'RIBN: Number in S( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7830/R:ribn/F:a-main1.f'
            goto 999
         end if
         pxys(ib,3)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'sx' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'RIBN: Number in SX( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7843/R:ribn/F:a-main1.f'
            goto 999
         end if
         pxys(ib,4)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'sy' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'RIBN: Number in SY( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7856/R:ribn/F:a-main1.f'
            goto 999
         end if
         pxys(ib,5)=rnm


      else if(lum(ic).eq.'a') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'RIBN: Number in A( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:7869/R:ribn/F:a-main1.f'
            goto 999
         end if
         pxys(ib,6)=rnm


      else if(lum(ic).eq.'t') then

               pxys(ib,7) = pxys(ib,7) + 3.0


      else if(lum(ic).eq.'z') then

               pxys(ib,7) = pxys(ib,7) - 1.0

               if( pxys(ib,7) .le. 1.0 ) pxys(ib,7) = 1.0


      else if(lum(ic)//lum(ic+1) .eq. 'cl' ) then

         ic=ic+1
           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               pxys(ib,11) = rcol(1)
               pxys(ib,12) = rcol(2)
               pxys(ib,13) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'RIBN: Number in CL( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:7899/R:ribn/F:a-main1.f'
               goto 999
            end if


      else if(lum(ic)//lum(ic+1) .eq. 'cb' ) then

         ic=ic+1
           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               pxys(ib,14) = rcol(1)
               pxys(ib,15) = rcol(2)
               pxys(ib,16) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'RIBN: Number in CB( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:7916/R:ribn/F:a-main1.f'
               goto 999
            end if


      else if(lum(ic)//lum(ic+1) .eq. 'cs' ) then

         ic=ic+1
           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               pxys(ib, 8) = rcol(1)
               pxys(ib, 9) = rcol(2)
               pxys(ib,10) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'RIBN: Number in CS( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:7933/R:ribn/F:a-main1.f'
               goto 999
            end if

      end if

*----------------------------------------------------------------------*

      goto 200

  500 continue

      do 300 i=1,2
        if(pxys(ib,i).lt.-r0max) then
           m_err = 'RIBN: You need X() Y() '//
     &             'Position and Number of conner of Star at Least.'
           ErrCha = ''
           ErrID = 'L:7950/R:ribn/F:a-main1.f'
           goto 999
        end if
  300 continue

      return
  999 ierr=1
      l_err = ill(jsn)
      k_err = jsn
      return
      end

************************************************************************
*                                                                      *
      subroutine star(dum,lum,in,ill,ilf,jsn,ib,pxys,ierr)
*                                                                      *
*      PURPOSE  :   READ STAR                                          *
*                                                                      *
*      FORMAT   :  STAR: X( ) Y( ) S( ) SX( ) SY( ) CL( ) CB( ) A( )   *
*                        PL( ) V( )  T or Z                            *
*                                                                      *
*           PXYS(IB,11) :  1-> X-COORDINATE OF ORIGIN    X( )          *
*                          2-> Y-COORDINATE OF ORIGIN    Y( )          *
*                          3-> NUMBER OF CORNER         PL( )          *
*                          4-> SCAL OF POLYGON           S( )          *
*                          5-> X SCAL OF POLYGON        SX( )          *
*                          6-> Y SCAL OF POLYGON        SY( )          *
*                          7-> ANGLE OF POLYGON          A( )          *
*                          8-> WIDTH OF LINE  ( 4 dd default ) T or Z  *
*                          9-> VALLEY OF STAR  0.5 IS DEFAULT          *
*                         10,11,12-> COLOR OF POLYGON LINE    CL( )    *
*                         13,14,15-> COLOR OF POLYGON INTERIA CB( )    *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      dimension pxys(mc,15)

      common /error/ m_err, l_err, k_err
      character m_err*200

      dimension ill(0:9), ilf(0:9)
      dimension rcol(3)

      character c1*1, c2*1

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

               c1 ='('
               c2 =')'

               pxys(ib,1)  = -r1max
               pxys(ib,2)  = -r1max
               pxys(ib,3)  = 5.0
               pxys(ib,4)  = 1.0
               pxys(ib,5)  = 1.0
               pxys(ib,6)  = 1.0
               pxys(ib,7)  = -r1max
               pxys(ib,8) = 4.0
               pxys(ib,9) = -r1max
               pxys(ib,10)  = -r1max
               pxys(ib,11)  = 1.0
               pxys(ib,12)  = 1.0
               pxys(ib,13)  = -r1max
               pxys(ib,14)  = 1.0
               pxys(ib,15)  = 1.0

*-----------------------------------------------------------------------

      ic = in - 1

  200 ic = ic + 1

         if(ic.gt.icolm) goto 500
         if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 200

      if(lum(ic).ne.'x'.and.lum(ic).ne.'y'.and.lum(ic).ne.'c'.and.
     &   lum(ic).ne.'s'.and.lum(ic).ne.'p'.and.lum(ic).ne.'t'.and.
     &   lum(ic).ne.'v'.and.lum(ic).ne.'z'.and.lum(ic).ne.'a') then
         m_err = 'STAR: X() Y() PL() V() S() SX() SY() A() '//
     &           'CL() CB() T Z ; Unexpected Parameter.'
         ErrCha = ''
         ErrID = 'L:8045/R:star/F:a-main1.f'
         goto 999
      end if

      if(lum(ic).eq.'x') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'STAR: Number in X( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:8056/R:star/F:a-main1.f'
            goto 999
         end if
         pxys(ib,1)=rnm


      else if(lum(ic).eq.'y') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'STAR: Number in Y( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:8069/R:star/F:a-main1.f'
            goto 999
         end if
         pxys(ib,2)=rnm


      else if(lum(ic).eq.'v') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'STAR: Number in V( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:8082/R:star/F:a-main1.f'
            goto 999
         end if
         pxys(ib,9)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'pl' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'STAR: Number in PL( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:8095/R:star/F:a-main1.f'
            goto 999
         end if
         pxys(ib,3)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 's(' ) then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'STAR: Number in S( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:8108/R:star/F:a-main1.f'
            goto 999
         end if
         pxys(ib,4)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'sx' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'STAR: Number in SX( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:8121/R:star/F:a-main1.f'
            goto 999
         end if
         pxys(ib,5)=rnm


      else if(lum(ic)//lum(ic+1) .eq. 'sy' ) then

         ic=ic+2
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'STAR: Number in SY( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:8134/R:star/F:a-main1.f'
            goto 999
         end if
         pxys(ib,6)=rnm


      else if(lum(ic).eq.'a') then

         ic=ic+1
         call pnum(lum,ic,icolm,rnm,ierr)
         if(ierr.ne.0) then
            m_err = 'STAR: Number in A( ) is Wrong.'
            ErrCha = ''
            ErrID = 'L:8147/R:star/F:a-main1.f'
            goto 999
         end if
         pxys(ib,7)=rnm


      else if(lum(ic).eq.'t') then

               pxys(ib,8) = pxys(ib,8) + 3.0


      else if(lum(ic).eq.'z') then

               pxys(ib,8) = pxys(ib,8) - 1.0

               if( pxys(ib,8) .le. 1.0 ) pxys(ib,8) = 1.0


      else if(lum(ic)//lum(ic+1) .eq. 'cl' ) then

         ic=ic+1
           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               pxys(ib,10) = rcol(1)
               pxys(ib,11) = rcol(2)
               pxys(ib,12) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'STAR: Number in CL( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:8177/R:star/F:a-main1.f'
               goto 999
            end if


      else if(lum(ic)//lum(ic+1) .eq. 'cb' ) then

         ic=ic+1
           call dcols(1,lum,ic,rcol,ierr,icolm,c1,c2)

               pxys(ib,13) = rcol(1)
               pxys(ib,14) = rcol(2)
               pxys(ib,15) = rcol(3)

            if(ierr.eq.1) then
               m_err = 'STAR: Number in CB( ) is Wrong.'
               ErrCha = ''
               ErrID = 'L:8194/R:star/F:a-main1.f'
               goto 999
            end if

      end if

*----------------------------------------------------------------------*

      goto 200

  500 continue

      do 300 i=1,2
        if(pxys(ib,i).lt.-r0max) then
           m_err = 'STAR: You need X() Y() '//
     &             'Position of Star at Least.'
           ErrCha = ''
           ErrID = 'L:8211/R:star/F:a-main1.f'
           goto 999
        end if
  300 continue

      return
  999 ierr=1
      l_err = ill(jsn)
      k_err = jsn
      return
      end

************************************************************************
*                                                                      *
      subroutine lgline(ixlog,iylog,xmax,xmin,ymax,ymin,
     &                  regx,regy,regs,regd,idbg,xfac,
     &                  ncom,iycm,rycm,regxl,regxc,regyd,jhf,txr,
     &                  clal,clmo,rlptl,sybw,erwd)
cKN 2024/01/24
*                                                                      *
*      PURPOSE  :   WRITE LEGEND LINES                                 *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter ( ss0 = 0.254 )
      parameter ( regx0  = 0.6, regxl0 = 1.4 )
      parameter ( regxc0 = 0.4, regyd0 = 0.9 )

      common /frm/  xal, yal
      common /con/  cm, dd

      dimension iycm(mc,8),rycm(mc,7)

      dimension xyp(6,5)
      dimension rcol(3), rcob(3)
      dimension clal(3)

*-----------------------------------------------------------------------

      do 450 i = 1, ncom

                      ityp   = iycm(i,2)
                      iwi    = iycm(i,3)
                      ipol   = 0
                      imar   = iycm(i,4)
                      rcol(1) = rycm(i,2)
                      rcol(2) = rycm(i,3)
                      rcol(3) = rycm(i,4)
                      rcob(1) = rycm(i,5)
                      rcob(2) = rycm(i,6)
                      rcob(3) = rycm(i,7)

         if( ityp .ne. 11 ) then

               if( clal(1) .gt. -r0max ) then
                  rcol(1) = clal(1)
                  rcol(2) = clal(2)
                  rcol(3) = clal(3)
               end if

               if( clmo .gt. -r0max .and. rcol(1) .gt. 0.0 )
     &             rcol(1) = -2.0
               if( rcol(1) .lt. -2.5 ) rcol(1) = -2.0

            if( imar .eq.  3 .or. imar .eq.  5 .or.
     &          imar .eq.  7 .or. imar .eq.  9 .or.
     &          imar .eq. 11 .or. imar .eq. 13 ) then

               if( clal(1) .gt. -r0max ) rcob(1) = -1.0

               if( clmo .gt. -r0max .and. rcob(1) .gt. 0.0 )
     &             rcob(1) = -1.0

               if( clmo .gt. -r0max .and. rcob(1) .lt. -r0max )
     &             rcob(1) = -1.0

               if( rcob(1) .lt. -2.5 ) rcob(1) = -1.0

            end if

         else if( ityp .eq. 11 ) then

               if( clmo .gt. -r0max .and. rcol(1) .gt. 0.0 )
     &         call ctomo(rcol,clmo)

               if( rcol(1) .lt. -2.5 ) rcol(1) = -1.0

         end if

                      icoma  = iycm(i,5)
                      facsiz = rycm(i,1)

*-----------------------------------------------------------------------

         if( imar .le. 0 ) then

            if( iycm(i,6) .eq. 0 ) then

               if( ityp .ne. 11 ) then

                  if( icoma .ne. 1 .and. iycm(i,8) .eq. 0 ) then

*                   ---- NORMAL LINE WITHOUT SYMBOL ----

                     ild   = 2
                     ixdc  = 0
                     iydc  = 0
                     rleng = regxl * xal

                     xyp(1,1) = regx
                     xyp(1,2) = regx + regxl
                     xyp(2,1) = regy - dble( i - 1 ) * regyd
                     xyp(2,2) = xyp(2,1)

                        call wline01(2,jhf,ild,ixdc,iydc,rleng,idbg,
     &                       ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                       xyp,xfac,rlptl)

                  end if

                  if( iycm(i,7) .ne. 0 .or. iycm(i,8) .ne. 0 ) then

*                    ---- DRAW ERROR BARS ----

                           ild   = 1

                           xyp(1,1) = regx + regxl / 2.0
                           xyp(2,1) = regy - dble( i - 1 ) * regyd

                           ixdc  = iycm(i,8)
                           iydc  = iycm(i,7)

                           xddl  = max( ss0 * facsiz
     &                                / xal, regxl * 0.2 )
                           yddl  = max( ss0 * facsiz
     &                                / yal, regyd * 0.3 )

                           xyp(3,1) = xyp(1,1) - xddl
                           xyp(4,1) = xyp(1,1) + xddl
                           xyp(5,1) = xyp(2,1) - yddl
                           xyp(6,1) = xyp(2,1) + yddl

                        call weror01(2,jhf,ild,ixdc,iydc,rleng,idbg,
     &                       ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                       xyp,xfac,erwd)
cKN 2024/01/24

                  end if

               else

*                    ---- REGION ----

                     ild   = 4
                     ixdc  = 0
                     iydc  = 0
                     rleng = regxl * xal * 0.8

                     xyp(1,1) = regx + regxl * 0.1
                     xyp(1,2) = regx + regxl * 0.9
                     xyp(2,1) = regy - dble( i - 1 ) * regyd
     &                        - regyd / 4.0
                     xyp(2,2) = xyp(2,1)
                     xyp(1,3) = xyp(1,2)
                     xyp(1,4) = xyp(1,1)
                     xyp(2,3) = regy - dble( i - 1 ) * regyd
     &                        + regyd / 4.0
                     xyp(2,4) = xyp(2,3)

                        call wline01(2,jhf,ild,ixdc,iydc,rleng,idbg,
     &                       ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                       xyp,xfac,rlptl)

               end if

            else

*                ---- HISTOGRAM ----

                  ild   = 5
                  ixdc  = 0
                  iydc  = 0

                  xyp(1,1) = regx
                  xyp(1,2) = xyp(1,1)
                  xyp(1,3) = regx + regxl / 2.0
                  xyp(1,4) = xyp(1,3)
                  xyp(1,5) = regx + regxl
                  yone1    = regy - dble( i - 1 ) * regyd
                  xyp(2,1) = yone1 - regyd * 0.3
                  xyp(2,2) = yone1 + regyd * 0.3
                  xyp(2,3) = xyp(2,2)
                  xyp(2,4) = yone1
                  xyp(2,5) = yone1

                     rleng = 0.0

                  do 100 k = 2, ild

                     rleng = rleng + sqrt(
     &                       ( xyp(1,k) - xyp(1,k-1) )**2 * xal**2
     &                     + ( xyp(2,k) - xyp(2,k-1) )**2 * yal**2 )

  100             continue

                     call wline01(2,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,rlptl)

            end if

*-----------------------------------------------------------------------

         else if( imar .gt. 0 ) then


            if( imar .eq. 17 .or. imar .eq. 18 .or.
     &          imar .eq. 19 .or. imar .eq. 20 .or.
     &          imar .eq. 21 .or. imar .eq. 22 ) then

               if( icoma .ne. 1 .and. iycm(i,8) .eq. 0 ) then

*                 ---- CLIP MARKER AND DRAW LINE ----

                        ild   = 1
                        ixdc  = 0
                        iydc  = 0

                        xyp(1,1) = regx + regxl / 2.0
                        xyp(2,1) = regy - dble( i - 1 ) * regyd

                     call wsymb01(3,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,sybw)


                        ild   = 2
                        ixdc  = 0
                        iydc  = 0
                        rleng = regxl * xal

                        xyp(1,1) = regx
                        xyp(1,2) = regx + regxl
                        xyp(2,1) = regy - dble( i - 1 ) * regyd
                        xyp(2,2) = xyp(2,1)

                     call wline01(3,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,rlptl)

               end if

               if( iycm(i,7) .ne. 0 .or. iycm(i,8) .ne. 0 ) then

*                 ---- CLIP MARKER AND DRAW ERROR BARS ----

                        ild   = 1

                        xyp(1,1) = regx + regxl / 2.0
                        xyp(2,1) = regy - dble( i - 1 ) * regyd

                        ixdc  = iycm(i,8)
                        iydc  = iycm(i,7)

                        if( icoma .ne. 1 ) ixdc = 0

                        xddl  = max( ss0 * facsiz / xal, regxl * 0.2 )
                        yddl  = max( ss0 * facsiz / yal, regyd * 0.3 )

                        xyp(3,1) = xyp(1,1) - xddl
                        xyp(4,1) = xyp(1,1) + xddl
                        xyp(5,1) = xyp(2,1) - yddl
                        xyp(6,1) = xyp(2,1) + yddl


                     call weror01(3,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,erwd)
cKN 2024/01/24

               end if

            else

               if( icoma .ne. 1 .and. iycm(i,8) .eq. 0 ) then

*                 ---- DRAW LINE ----

                        ild   = 2
                        ixdc  = 0
                        iydc  = 0
                        rleng = regxl * xal

                        xyp(1,1) = regx
                        xyp(1,2) = regx + regxl
                        xyp(2,1) = regy - dble( i - 1 ) * regyd
                        xyp(2,2) = xyp(2,1)

                     call wline01(2,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,rlptl)

               end if

               if( iycm(i,7) .ne. 0 .or. iycm(i,8) .ne. 0 ) then

*                 ---- DRAW ERROR BARS ----

                        ild   = 1

                        xyp(1,1) = regx + regxl / 2.0
                        xyp(2,1) = regy - dble( i - 1 ) * regyd

                        ixdc  = iycm(i,8)
                        iydc  = iycm(i,7)

                        xddl  = max( ss0 * facsiz / xal, regxl * 0.2 )
                        yddl  = max( ss0 * facsiz / yal, regyd * 0.3 )

                        xyp(3,1) = xyp(1,1) - xddl
                        xyp(4,1) = xyp(1,1) + xddl
                        xyp(5,1) = xyp(2,1) - yddl
                        xyp(6,1) = xyp(2,1) + yddl

                     call weror01(2,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,erwd)
cKN 2024/01/24

               end if

            end if

            if( imar .gt. 0 ) then

*                 ---- DRAW MAKERS ----

                        ild   = 1
                        ixdc  = 0
                        iydc  = 0

                        xyp(1,1) = regx + regxl / 2.0
                        xyp(2,1) = regy - dble( i - 1 ) * regyd

                     call wsymb01(2,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,sybw)

            end if


         end if

  450 continue

      return
      end

************************************************************************
*                                                                      *
      subroutine zmult(lum,in,xorg,yorg,angz,
     &                 noxt,noxn,noyt,noyn,ierr)
*                                                                      *
*                                                                      *
*        PURPOSE    :  DEFINE NEW AXISIS                               *
*        FORMAT:   Z: XORG( ) YORG( ) ANGL( ) NOXT,NOYT,NOXN,NOYN      *
*                                                                      *
*          XORG     : POSITION OF THE X-ORIGIN, UNIT IS PREVIUS X-AXIS *
*          YORG     : POSITION OF THE Y-ORIGIN, UNIT IS PREVIUS Y-AXIS *
*          ANGZ     : ANGLE OF THE NEW GRAPH REL TO PREVIUS GRAPH      *
*                                                                      *
*          NOXT     : NO X-AXIS TEXT                                   *
*          NOXN     : NO X-AXIS NUMBER                                 *
*          NOYT     : NO Y-AXIS TEXT                                   *
*          NOYN     : NO Y-AXIS NUMBER                                 *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character lum(ichrl)*1,dum4*4

      common /error/ m_err, l_err, k_err
      character m_err*200

      character c1*1,c2*1,c3*1

      character tub*1
      tub = char(9)

         c1 = '('
         c2 = ')'

*-----------------------------------------------------------------------

      xorg = 0.0
      yorg = 0.0
      angz = -r1max

      noxt = 1
      noxn = 1
      noyt = 1
      noyn = 1

      ic=in-1

  100 ic=ic+1
        if(ic.gt.icolm) goto 500
        if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 100

        if(lum(ic).ne.'x'.and. lum(ic).ne.'y'.and.
     &     lum(ic).ne.'a'.and. lum(ic).ne.'n') then
           m_err = 'Z: XORG() YORG() ANGL() NOXT NOXN NOYT NOYN;'//
     &             ' Unexpected Parameter.'
           ErrCha = ''
           ErrID = 'L:8642/R:zmult/F:a-main1.f'
           goto 999
        end if
        if(ic+3.gt.icolm) goto 999

         do k = ic, ic + 3
            dum4(k-ic+1:k-ic+1) = lum(k)
         end do

        ic = ic + 3

        if(dum4.eq.'xorg') then
              ic=ic+1
           call pnum(lum,ic,icolm,xorg,ierr)
              if(ierr.ne.0) then
                 m_err = 'Z: Number in XORG() is Wrong.'
                 ErrCha = ''
                 ErrID = 'L:8659/R:zmult/F:a-main1.f'
                 goto 999
              end if
              ic=ic+1
              goto 100

        else if(dum4.eq.'yorg') then
              ic=ic+1
           call pnum(lum,ic,icolm,yorg,ierr)
              if(ierr.ne.0) then
                 m_err = 'Z: Number in YORG() is Wrong.'
                 ErrCha = ''
                 ErrID = 'L:8671/R:zmult/F:a-main1.f'
                 goto 999
              end if
              ic=ic+1
              goto 100

        else if(dum4.eq.'angl') then
              ic=ic+1
           call pnum(lum,ic,icolm,angz,ierr)
              if(ierr.ne.0) then
                 m_err = 'Z: Number in ANGL() is Wrong.'
                 ErrCha = ''
                 ErrID = 'L:8683/R:zmult/F:a-main1.f'
                 goto 999
              end if
              ic=ic+1
              goto 100

        else if(dum4.eq.'noyn') then
              noyn=0
              ic=ic+1
              goto 100

        else if(dum4.eq.'noyt') then
              noyt=0
              ic=ic+1
              goto 100

        else if(dum4.eq.'noxn') then
              noxn=0
              ic=ic+1
              goto 100

        else if(dum4.eq.'noxt') then
              noxt=0
              ic=ic+1
              goto 100

        else

           m_err = 'Z: XORG() YORG() ANGL();'//
     &             ' Unexpected Parameter.'
           ErrCha = ''
           ErrID = 'L:8714/R:zmult/F:a-main1.f'
           goto 999

        end if

  500 continue

      return

*-----------------------------------------------------------------------

  999 ierr=1

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine warrw01(icw,jhf,idbg,clal,clmo,
     &                   x1,y1,x2,y2,xl,yl,rcol,dd,rnl,ara)
*                                                                      *
*                                                                      *
*        PURPOSE    :  WRITE ARROWS ON jhf                             *
*                                                                      *
*              ICW  :  0 -> ARROW ONLY, 1-> ARROW + LINE FOR COMMENT   *
*                                                                      *
*       X1,Y1,X2,Y2 : COORDINATE OF THE INITIAL AND FINAL POINT        *
*              RCOL : COLOR OF ARROWS                                  *
*                DD : WIDTH OF ARROW LINE                              *
*               RNL : -1 -> NO ARROW LINE, 1 -> WITH ARROW LINE        *
*               ARA : ANGLE OF THE ARROW RELATIVE TO DEFAULT           *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter (ddcm = 8.4666667e-3)

      parameter ( rang0 =  5.0 )
      parameter ( rlng  = 12.0 )

      common /frm/  xal, yal
      dimension clal(3), rcol(3)

*-----------------------------------------------------------------------

            rang = rang0 / ara

            ahed = ( rang + 1.0 ) / 2.0 * dd * ddcm
            aang = rang * dd * ddcm
            alng = rlng * dd * ddcm
            aalg = alng / rang

            x1 = x1 * xal
            x2 = x2 * xal
            y1 = y1 * yal
            y2 = y2 * yal

            rleng = sqrt( ( x2 - x1 )**2 + ( y2 - y1 )**2 )

            if( rleng .lt. alng ) rnl = -1.0

*-----------------------------------------------------------------------

         if( rnl .gt. 0.0 ) then

               x3 = x2 + ahed / rleng * ( x1 - x2 )
               y3 = y2 + ahed / rleng * ( y1 - y2 )

               x6 = x2 + alng / rleng * ( x1 - x2 )
               y6 = y2 + alng / rleng * ( y1 - y2 )

         else

               x3 = x1 + 2.0 / 3.0 * alng / rleng * ( x2 - x1 )
               y3 = y1 + 2.0 / 3.0 * alng / rleng * ( y2 - y1 )

               x6 = x1 - 1.0 / 3.0 * alng / rleng * ( x2 - x1 )
               y6 = y1 - 1.0 / 3.0 * alng / rleng * ( y2 - y1 )

               x2 = x3
               y2 = y3

         end if


         if( x2 .eq. x6 ) then

               x4 = x6 - aalg
               x5 = x6 + aalg
               y4 = y6
               y5 = y6

         else if( y2 .eq. y6 ) then

               y4 = y6 - aalg
               y5 = y6 + aalg
               x4 = x6
               x5 = x6

         else

               yxf = ( y2 - y6 ) / ( x2 - x6 )
               yyf = sqrt( yxf**2 + 1.0 )
               xxf = sqrt( 1.0 / yxf**2 + 1.0 )

            if( yxf .gt. 0.0 ) then

               x4 = x6 - aalg / xxf
               y4 = y6 + aalg / yyf
               x5 = x6 + aalg / xxf
               y5 = y6 - aalg / yyf

            else

               x4 = x6 + aalg / xxf
               y4 = y6 + aalg / yyf
               x5 = x6 - aalg / xxf
               y5 = y6 - aalg / yyf

            end if

         end if

               x1 = x1 / xal
               x2 = x2 / xal
               x3 = x3 / xal
               x4 = x4 / xal
               x5 = x5 / xal

               y1 = y1 / yal
               y2 = y2 / yal
               y3 = y3 / yal
               y4 = y4 / yal
               y5 = y5 / yal

*-----------------------------------------------------------------------

            iwi = nint(dd)


               if( clal(1) .gt. -r0max ) then
                  rcol(1) = clal(1)
                  rcol(2) = clal(2)
                  rcol(3) = clal(3)
               end if

               if( clmo .gt. -r0max .and. rcol(1) .gt. 0.0 )
     &            rcol(1) = -2.0
               if( rcol(1) .lt. -r0max ) rcol(1) = -2.0


                  write(jhf,'(3f7.3,'' sc sd0 '',I2,'' dd lw'')')
     &                       rcol, iwi

*-----------------------------------------------------------------------

               write(jhf,'(2g14.5,'' m '',2g14.5,'' l ''/
     &                     2g14.5,'' l cp fl'')') x2, y2, x4, y4, x5, y5

                  call bbox(3,0,x2,y2,0.d0)
                  call bbox(3,0,x4,y4,0.d0)
                  call bbox(3,0,x5,y5,0.d0)

*-----------------------------------------------------------------------

         if( rnl .gt. 0.0 ) then

            if( icw .eq. 0 ) then

               write(jhf,'(2g14.5,'' m '',
     &                     2g14.5,'' l st'')') x3, y3, x1, y1

                  call bbox(3,0,x1,y1,0.d0)

            else

               write(jhf,'(2g14.5,'' m '',2g14.5,'' l ''/
     &                     2g14.5,'' l st'')') x3, y3, x1, y1, xl, yl

                  call bbox(3,0,x1,y1,0.d0)
                  call bbox(3,0,xl,yl,0.d0)

            end if

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine wsbox(jhf,idbg,clal,clmo,pxys,ib,
     &                 xmin,xmax,ymin,ymax,ixlog,iylog)
*                                                                      *
*                                                                      *
*        PURPOSE    :  WRITE BOX ON jhf                                *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter (ddcm = 8.4666667e-3)

      parameter ( tfs0  = 17.0 )
      parameter ( unit0 =  2.0 )

      common /frm/  xal, yal
      common /con/  cm, dd

      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb

      common /box/  bxws, bxwl, bxss, bxsl, bxds, bxdl, bxls, bxll

      dimension pxys(mc,16)
      dimension clal(3), clgb(3), clgl(3), clgs(3)
      dimension bcb(3), bcl(3), bcs(3)

*-----------------------------------------------------------------------
*           TRANSLATE AND ROTATE
*-----------------------------------------------------------------------

                  x1 = zonep(xmin,xmax,pxys(ib,1),ixlog) * xal * cm
                  y1 = zonep(ymin,ymax,pxys(ib,2),iylog) * yal * cm

                     xpsain = x1
                     ypsain = y1

                     write(jhf,'(2g14.5,'' TR'')') x1, y1

                  if( pxys(ib,6) .gt. -r0max ) then

                     iang   = 1
                     angain = pxys(ib,6)

                     write(jhf,'(g14.5,'' rotate'')') angain

                  else

                     iang   = 0
                     angain = 0.0

                  end if

                     call bbox(1,1,xpsain,ypsain,angain)

                     xpsin = 0.0
                     ypsin = 0.0


*-----------------------------------------------------------------------
*           TYPE AND COLOR, SIZE, DIMENSION  DEFINITION
*-----------------------------------------------------------------------

                     ilbox = nint( pxys(ib,7) )

                     clgb(1) = pxys(ib,14)
                     clgb(2) = pxys(ib,15)
                     clgb(3) = pxys(ib,16)

                     clgl(1) = pxys(ib,11)
                     clgl(2) = pxys(ib,12)
                     clgl(3) = pxys(ib,13)

                     clgs(1) = pxys(ib,8)
                     clgs(2) = pxys(ib,9)
                     clgs(3) = pxys(ib,10)

                     tfs  = tfs0 * pxys(ib,3)

                     b1x  = - tfs / 2.0 * pxys(ib,4)
                     b2x  =   tfs / 2.0 * pxys(ib,4)
                     b1y  = - tfs / 2.0 * pxys(ib,5)
                     b2y  =   tfs / 2.0 * pxys(ib,5)

*-----------------------------------------------------------------------
*           WRITE BOX ON jhf
*-----------------------------------------------------------------------

                        if( clgb(1) .lt. -r0max .or.
     &                      clmo .gt. -r0max .or.
     &                      clal(1) .gt. -r0max ) clgb(1) = -1.0

                        if( clal(1) .gt. -r0max ) then

                           clgl(1) = clal(1)
                           clgl(2) = clal(2)
                           clgl(3) = clal(3)
                           clgs(1) = clal(1)
                           clgs(2) = clal(2)
                           clgs(3) = clal(3)

                        end if

                        if( clgl(1) .lt. -r0max .or.
     &                      clmo .gt. -r0max ) clgl(1) = -2.0

                        if( clgs(1) .lt. -r0max .or.
     &                      clmo .gt. -r0max ) clgs(1) = -2.0

                        bcb(1) = clgb(1)
                        bcb(2) = clgb(2)
                        bcb(3) = clgb(3)
                        bcl(1) = clgl(1)
                        bcl(2) = clgl(2)
                        bcl(3) = clgl(3)
                        bcs(1) = clgs(1)
                        bcs(2) = clgs(2)
                        bcs(3) = clgs(3)


                  if( ilbox .eq. 10 .or. ilbox .eq. 11. or.
     &                ilbox .eq. 20 .or. ilbox .eq. 21. or.
     &                ilbox .eq. 30 .or. ilbox .eq. 31. or.
     &                ilbox .eq. 40 .or. ilbox .eq. 41. or.
     &                ilbox .eq. 50 .or. ilbox .eq. 51. ) then

                        bds = tfs * bxwl

                  else

                        bds = tfs * bxwl * 2.0

                  end if

                  if( ilbox .eq. 11 .or. ilbox .eq. 13. or.
     &                ilbox .eq. 21 .or. ilbox .eq. 23. or.
     &                ilbox .eq. 31 .or. ilbox .eq. 33. or.
     &                ilbox .eq. 41 .or. ilbox .eq. 43. or.
     &                ilbox .eq. 51 .or. ilbox .eq. 53. ) then

                        bdd = tfs * bxll

                  else

                        bdd = tfs * bxls

                  end if

                  if( ilbox .eq. 20 .or. ilbox .eq. 21. or.
     &                ilbox .eq. 40 .or. ilbox .eq. 41. or.
     &                ilbox .eq. 50 .or. ilbox .eq. 51. ) then

                        bdc = tfs * bxsl

                  else
     &            if( ilbox .eq. 22 .or. ilbox .eq. 23. or.
     &                ilbox .eq. 42 .or. ilbox .eq. 43. or.
     &                ilbox .eq. 52 .or. ilbox .eq. 53. ) then

                        bdc = tfs * bxsl * 2.0

                  end if

                  if( ilbox .eq. 30 .or. ilbox .eq. 31. ) then

                        bdc = tfs * bxdl

                  else
     &            if( ilbox .eq. 32 .or. ilbox .eq. 33. ) then

                        bdc = tfs * bxdl * 2.0

                  end if

                     b1x = b1x - bds
                     b2x = b2x + bds

                     b1y = b1y - bds
                     b2y = b2y + bds

*-----------------------------------------------------------------------

                     iccb = 0

                     call wrbox(jhf,idbg,0,ilbox,iccb,
     &                          b1x, b2x, b1y, b2y,
     &                          bdd, bcb, bcl, bdc, bcs)

*-----------------------------------------------------------------------

                  if( ilbox .ge. 40 ) then

                     b2x = b2x + bdc
                     b1y = b1y - bdc

                  end if

                     call bbox(1,0,b1x,b1y,0.d0)
                     call bbox(1,0,b1x,b2y,0.d0)
                     call bbox(1,0,b2x,b1y,0.d0)
                     call bbox(1,0,b2x,b2y,0.d0)


*-----------------------------------------------------------------------
*           TRANSLATE AND ROTATE
*-----------------------------------------------------------------------

                  if( iang .eq. 0 ) then

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' TR'')') -xpsain, -ypsain

                     call bbox(1,1,-xpsain,-ypsain,0.d0)

                  else

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(g14.5,'' rotate'')') -angain

                     write(jhf,'(2g14.5,'' TR'')') -xpsain, -ypsain

                     call bbox(1,1,0.d0,0.d0,-angain)
                     call bbox(1,1,-xpsain,-ypsain,0.d0)

                  end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine wpolg(jhf,idbg,clal,clmo,pxys,ib,
     &                 xmin,xmax,ymin,ymax,ixlog,iylog)
*                                                                      *
*                                                                      *
*        PURPOSE    :  WRITE POLYGON ON jhf                            *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter (ddcm = 8.4666667e-3)

      parameter ( unit0 =  24.0 )

      common /frm/  xal, yal
      common /con/  cm, dd

      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb

      dimension pxys(mc,14)
      dimension clal(3), rcoll(3), rcolb(3)

*-----------------------------------------------------------------------
*           TRANSLATE AND ROTATE
*-----------------------------------------------------------------------

                  x1 = zonep(xmin,xmax,pxys(ib,1),ixlog) * xal * cm
                  y1 = zonep(ymin,ymax,pxys(ib,2),iylog) * yal * cm

                     xpsain = x1
                     ypsain = y1

                     write(jhf,'(2g14.5,'' TR'')') x1, y1

                  if( pxys(ib,7) .gt. -r0max ) then

                     iang   = 1
                     angain = pxys(ib,7)

                     write(jhf,'(g14.5,'' rotate'')') angain

                  else

                     iang   = 0
                     angain = 0.0

                  end if

                     call bbox(1,1,xpsain,ypsain,angain)

                     xpsin = 0.0
                     ypsin = 0.0


*-----------------------------------------------------------------------
*           SCALE IN
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     pxys(ib,5) = pxys(ib,4) * pxys(ib,5)
                     pxys(ib,6) = pxys(ib,4) * pxys(ib,6)

                     write(jhf,'(2g14.5,'' scale'')')
     &                           pxys(ib,5), pxys(ib,6)

                     sxg0 = sxg0 * pxys(ib,5)
                     syg0 = syg0 * pxys(ib,6)

*-----------------------------------------------------------------------
*           CORNNER OF POLYGONS PATH
*-----------------------------------------------------------------------

                  if( pxys(ib,3) .lt. 0.0 ) then

                     npic = -1

                  else

                     npic =  1

                  end if

                     npol = nint( abs( pxys(ib,3) ) )

                  if( npol .lt. 3 ) then

                     npol = 0

                  else if( npol .gt. 20 ) then

                     npol = 20

                  end if

                  if( idbg .eq. 1 )  write(jhf,'()')

*-----------------------------------------------------------------------
*           CIRCLE
*-----------------------------------------------------------------------

               if( npol .eq. 0 ) then

                  write(jhf,'(g14.5,'' 0 M 0 0 '',g14.5,
     &                  '' 0 360 arc'')') unit0, unit0


                     dang = 2.0 * pi / 10.0

                  do 200 j = 1, 20

                     angn = dble( j - 1 ) * dang

                     xpos = unit0 * cos( angn )
                     ypos = unit0 * sin( angn )

                     call bbox(1,0,xpos,ypos,0.d0)

  200             continue

               end if

*-----------------------------------------------------------------------
*           POLYGONS
*-----------------------------------------------------------------------

               if( npol .gt. 0 ) then

                     dang = 2.0 * pi / dble( npol )

                  if( npic .gt. 0 ) then

                     anin = 90.0 * pi / 180.0

                  else

                     anin = 90.0 * pi / 180.0 - dang / 2.0

                  end if

                  do 100 j = 1, npol

                     angn = anin + dble( j - 1 ) * dang

                     xpos = unit0 * cos( angn )
                     ypos = unit0 * sin( angn )

                     call bbox(1,0,xpos,ypos,0.d0)


                     if( j .eq. 1 ) then

                        write(jhf,'(2g14.5,'' M'')') xpos, ypos

                     else if( j .eq. npol ) then

                        write(jhf,'(2g14.5,'' L cp'')') xpos, ypos

                     else

                        write(jhf,'(2g14.5,'' L'')') xpos, ypos

                     end if

  100             continue

               end if

*-----------------------------------------------------------------------
*           SCALE OUT
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' scale'')')
     &                  1./pxys(ib,5), 1./pxys(ib,6)

                     sxg0 = sxg0 / pxys(ib,5)
                     syg0 = syg0 / pxys(ib,6)

*-----------------------------------------------------------------------
*           FILL AND STRING POLYGONS PATH
*-----------------------------------------------------------------------

               iwi = nint( pxys(ib,14) )

               rcoll(1) = pxys(ib,8)
               rcoll(2) = pxys(ib,9)
               rcoll(3) = pxys(ib,10)

               rcolb(1) = pxys(ib,11)
               rcolb(2) = pxys(ib,12)
               rcolb(3) = pxys(ib,13)

               if( clal(1) .gt. -r0max ) then

                     rcoll(1) = clal(1)
                     rcoll(2) = clal(2)
                     rcoll(3) = clal(3)

                  if( rcolb(1) .lt. -r0max ) then

                     rcolb(1) = -3.0

                  else

                     rcolb(1) = clal(1)
                     rcolb(2) = clal(2)
                     rcolb(3) = clal(3)

                  end if

               end if

               if( clmo .gt. -r0max .and. rcoll(1) .gt. 0.0 )
     &             rcoll(1) = -2.0
               if( clmo .gt. -r0max .and. rcolb(1) .gt. 0.0 )
     &             rcolb(1) = -1.0

               if( rcoll(1) .lt. -r0max ) rcoll(1) = -2.0
               if( rcolb(1) .lt. -r0max ) rcolb(1) = -3.0

               if( rcolb(1) .gt. -2.5 ) then

                  write(jhf,'(''gs '',3f7.3,'' sc fl gr sd0 '',i2,
     &                        '' dd lw '',3f7.3,'' sc st'')')
     &                       rcolb, iwi, rcoll

               else

                  write(jhf,'(3f7.3,'' sc sd0 '',i2,'' dd lw st'')')
     &                       rcoll, iwi

               end if


*-----------------------------------------------------------------------
*           TRANSLATE AND ROTATE
*-----------------------------------------------------------------------

                  if( iang .eq. 0 ) then

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' TR'')') -xpsain, -ypsain

                     call bbox(1,1,-xpsain,-ypsain,0.d0)

                  else

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(g14.5,'' rotate'')') -angain

                     write(jhf,'(2g14.5,'' TR'')') -xpsain, -ypsain

                     call bbox(1,1,0.d0,0.d0,-angain)
                     call bbox(1,1,-xpsain,-ypsain,0.d0)

                  end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine wstar(jhf,idbg,clal,clmo,pxys,ib,
     &                 xmin,xmax,ymin,ymax,ixlog,iylog)
*                                                                      *
*                                                                      *
*        PURPOSE    :  WRITE STAR ON jhf                               *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter (ddcm = 8.4666667e-3)

      parameter ( unit0 =  24.0 )

      common /frm/  xal, yal
      common /con/  cm, dd

      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb

      dimension pxys(mc,15)
      dimension clal(3), rcoll(3), rcolb(3)

*-----------------------------------------------------------------------
*           TRANSLATE AND ROTATE
*-----------------------------------------------------------------------

                  x1 = zonep(xmin,xmax,pxys(ib,1),ixlog) * xal * cm
                  y1 = zonep(ymin,ymax,pxys(ib,2),iylog) * yal * cm

                     xpsain = x1
                     ypsain = y1

                     write(jhf,'(2g14.5,'' TR'')') x1, y1

                  if( pxys(ib,7) .gt. -r0max ) then

                     iang   = 1
                     angain = pxys(ib,7)

                     write(jhf,'(g14.5,'' rotate'')') angain

                  else

                     iang   = 0
                     angain = 0.0

                  end if

                     call bbox(1,1,xpsain,ypsain,angain)

                     xpsin = 0.0
                     ypsin = 0.0


*-----------------------------------------------------------------------
*           SCALE IN
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     pxys(ib,5) = pxys(ib,4) * pxys(ib,5)
                     pxys(ib,6) = pxys(ib,4) * pxys(ib,6)

                     write(jhf,'(2g14.5,'' scale'')')
     &                           pxys(ib,5), pxys(ib,6)

                     sxg0 = sxg0 * pxys(ib,5)
                     syg0 = syg0 * pxys(ib,6)

*-----------------------------------------------------------------------
*           CORNNER OF POLYGONS PATH
*-----------------------------------------------------------------------

                  if( pxys(ib,3) .lt. 0.0 ) then

                     npic = -1

                  else

                     npic =  1

                  end if

                     npol = nint( abs( pxys(ib,3) ) )

                  if( npol .lt. 3 ) then

                     npol = 3

                  else if( npol .gt. 40 ) then

                     npol = 40

                  end if

                  if( idbg .eq. 1 )  write(jhf,'()')

*-----------------------------------------------------------------------
*           POLYGONS
*-----------------------------------------------------------------------

                     dang = 2.0 * pi / dble( npol ) / 2.0

                     sunit = unit0 * cos( dang ) / 2.0

                  if( pxys(ib,9) .gt. -r0max ) then

                     sunit = sunit * pxys(ib,9)

                  end if

                  if( npic .gt. 0 ) then

                     anin = 90.0 * pi / 180.0

                  else

                     anin = 90.0 * pi / 180.0 - dang

                  end if

                  do 100 j = 1, 2 * npol

                     angn = anin + dble( j - 1 ) * dang

                     if( j / 2 * 2 .ne. j ) then

                        xpos = unit0 * cos( angn )
                        ypos = unit0 * sin( angn )

                     else

                        xpos = sunit * cos( angn )
                        ypos = sunit * sin( angn )

                     end if

                        call bbox(1,0,xpos,ypos,0.d0)


                     if( j .eq. 1 ) then

                        write(jhf,'(2g14.5,'' M'')') xpos, ypos

                     else if( j .eq. 2 * npol ) then

                        write(jhf,'(2g14.5,'' L cp'')') xpos, ypos

                     else

                        write(jhf,'(2g14.5,'' L'')') xpos, ypos

                     end if

  100             continue

*-----------------------------------------------------------------------
*           SCALE OUT
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' scale'')')
     &                  1./pxys(ib,5), 1./pxys(ib,6)

                     sxg0 = sxg0 / pxys(ib,5)
                     syg0 = syg0 / pxys(ib,6)

*-----------------------------------------------------------------------
*           FILL AND STRING POLYGONS PATH
*-----------------------------------------------------------------------

               iwi = nint( pxys(ib,8) )

               rcoll(1) = pxys(ib,10)
               rcoll(2) = pxys(ib,11)
               rcoll(3) = pxys(ib,12)
               rcolb(1) = pxys(ib,13)
               rcolb(2) = pxys(ib,14)
               rcolb(3) = pxys(ib,15)

               if( clal(1) .gt. -r0max ) then

                     rcoll(1) = clal(1)
                     rcoll(2) = clal(2)
                     rcoll(3) = clal(3)

                  if( rcolb(1) .lt. -r0max ) then

                     rcolb(1) = -3.0

                  else

                     rcolb(1) = clal(1)
                     rcolb(2) = clal(2)
                     rcolb(3) = clal(3)

                  end if

               end if

               if( clmo .gt. -r0max .and. rcoll(1) .gt. 0.0 )
     &             rcoll(1) = -2.0
               if( clmo .gt. -r0max .and. rcolb(1) .gt. 0.0 )
     &             rcolb(1) = -1.0

               if( rcoll(1) .lt. -r0max ) rcoll(1) = -2.0
               if( rcolb(1) .lt. -r0max ) rcolb(1) = -3.0

               if( rcolb(1) .gt. -2.5 ) then

                  write(jhf,'(''gs '',3f7.3,'' sc fl gr sd0 '',i2,
     &                        '' dd lw '',3f7.3,'' sc st'')')
     &                       rcolb, iwi, rcoll

               else

                  write(jhf,'(3f7.3,'' sc sd0 '',i2,'' dd lw st'')')
     &                       rcoll, iwi

               end if


*-----------------------------------------------------------------------
*           TRANSLATE AND ROTATE
*-----------------------------------------------------------------------

                  if( iang .eq. 0 ) then

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' TR'')') -xpsain, -ypsain

                     call bbox(1,1,-xpsain,-ypsain,0.d0)

                  else

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(g14.5,'' rotate'')') -angain

                     write(jhf,'(2g14.5,'' TR'')') -xpsain, -ypsain

                     call bbox(1,1,0.d0,0.d0,-angain)
                     call bbox(1,1,-xpsain,-ypsain,0.d0)

                  end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine wribn(jhf,idbg,clal,clmo,pxys,ib,
     &                 xmin,xmax,ymin,ymax,ixlog,iylog)
*                                                                      *
*                                                                      *
*        PURPOSE    :  WRITE RIBN ON jhf                               *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter (ddcm = 8.4666667e-3)

      parameter ( unit0 =  24.0 )

      common /frm/  xal, yal
      common /con/  cm, dd

      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb

      dimension pxys(mc,16)
      dimension clal(3), rcoll(3), rcolb(3), rcols(3)

*-----------------------------------------------------------------------
*           TRANSLATE AND ROTATE
*-----------------------------------------------------------------------

                  x1 = zonep(xmin,xmax,pxys(ib,1),ixlog) * xal * cm
                  y1 = zonep(ymin,ymax,pxys(ib,2),iylog) * yal * cm

                     xpsain = x1
                     ypsain = y1

                     write(jhf,'(2g14.5,'' TR'')') x1, y1

                  if( pxys(ib,6) .gt. -r0max ) then

                     iang   = 1
                     angain = pxys(ib,6)

                     write(jhf,'(g14.5,'' rotate'')') angain

                  else

                     iang   = 0
                     angain = 0.0

                  end if

                     call bbox(1,1,xpsain,ypsain,angain)

                     xpsin = 0.0
                     ypsin = 0.0


*-----------------------------------------------------------------------
*           SCALE AND UNIT
*-----------------------------------------------------------------------

                  ru = unit0 / 2.0

                  pxys(ib,4) = pxys(ib,3) * pxys(ib,4)
                  pxys(ib,5) = pxys(ib,3) * pxys(ib,5)

*-----------------------------------------------------------------------
*           COLOR AND LINE WIDTH DEFINITION
*-----------------------------------------------------------------------

               iwi = nint( pxys(ib,7) )

               rcoll(1) = pxys(ib,11)
               rcoll(2) = pxys(ib,12)
               rcoll(3) = pxys(ib,13)
               rcolb(1) = pxys(ib,14)
               rcolb(2) = pxys(ib,15)
               rcolb(3) = pxys(ib,16)
               rcols(1) = pxys(ib, 8)
               rcols(2) = pxys(ib, 9)
               rcols(3) = pxys(ib,10)

               if( clal(1) .gt. -r0max ) then

                     rcoll(1) = clal(1)
                     rcoll(2) = clal(2)
                     rcoll(3) = clal(3)
                     rcols(1) = -1.5

                  if( rcolb(1) .lt. -r0max ) then

                     rcolb(1) = -1.0

                  else

                     rcolb(1) = clal(1)
                     rcolb(2) = clal(2)
                     rcolb(3) = clal(3)

                  end if

               end if

               if( clmo .gt. -r0max .and. rcoll(1) .gt. 0.0 )
     &             rcoll(1) = -2.0
               if( clmo .gt. -r0max .and. rcolb(1) .gt. 0.0 )
     &             rcolb(1) = -1.0
               if( clmo .gt. -r0max .and. rcols(1) .gt. 0.0 )
     &             rcols(1) = -1.5

               if( rcoll(1) .lt. -r0max ) rcoll(1) = -2.0
               if( rcolb(1) .lt. -r0max ) rcolb(1) = -1.0
               if( rcols(1) .lt. -r0max ) rcols(1) = -1.5

*-----------------------------------------------------------------------
* (1)       SCALE IN
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' scale'')')
     &                           pxys(ib,4), pxys(ib,5)

                     sxg0 = sxg0 * pxys(ib,4)
                     syg0 = syg0 * pxys(ib,5)

*-----------------------------------------------------------------------
*           RIBBON PART(1)
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                        write(jhf,'(2g14.5,'' M'')')
     &                  -4.5*ru,  1.75*ru

                        write(jhf,'(4g14.5,/2g14.5,'' curveto'')')
     &                  -4.5*ru,  2.35*ru, -6.5*ru, 2.00*ru,
     &                  -7.5*ru,  3.15*ru

                        call bbox(1,0,-7.5*ru,3.15*ru,0.d0)

                        write(jhf,'(4g14.5,/2g14.5,'' curveto'')')
     &                  -7.2*ru, 2.15*ru, -6.8*ru,  1.8*ru,
     &                  -6.5*ru, 1.40*ru

                        write(jhf,'(2g14.5,'' L'')')
     &                  -7.5*ru,  1.25*ru

                        call bbox(1,0,-7.5*ru,1.25*ru,0.d0)

                        write(jhf,'(4g14.5,/2g14.5,'' curveto cp'')')
     &                  -6.5*ru, -0.10*ru, -4.5*ru,  0.15*ru,
     &                  -4.5*ru, -0.25*ru

                        write(jhf,'(2g14.5,'' M'')')
     &                  4.5*ru,  1.75*ru

                        write(jhf,'(4g14.5,/2g14.5,'' curveto'')')
     &                  4.5*ru,  2.35*ru, 6.5*ru, 2.00*ru,
     &                  7.5*ru,  3.15*ru

                        call bbox(1,0,7.5*ru,3.15*ru,0.d0)

                        write(jhf,'(4g14.5,/2g14.5,'' curveto'')')
     &                  7.2*ru, 2.15*ru, 6.8*ru,  1.8*ru,
     &                  6.5*ru, 1.40*ru

                        write(jhf,'(2g14.5,'' L'')')
     &                  7.5*ru,  1.25*ru

                        call bbox(1,0,7.5*ru,1.25*ru,0.d0)

                        write(jhf,'(4g14.5,/2g14.5,'' curveto cp'')')
     &                  6.5*ru, -0.10*ru, 4.5*ru,  0.15*ru,
     &                  4.5*ru, -0.25*ru

*-----------------------------------------------------------------------
*           SCALE OUT
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' scale'')')
     &                  1./pxys(ib,4), 1./pxys(ib,5)

                     sxg0 = sxg0 / pxys(ib,4)
                     syg0 = syg0 / pxys(ib,5)

*-----------------------------------------------------------------------
*           FILL AND STRING RIBBON PATH
*-----------------------------------------------------------------------

                  write(jhf,'(''gs '',3f7.3,'' sc fl gr sd0 '',i2,
     &                        '' dd lw '',3f7.3,'' sc st'')')
     &                       rcolb, iwi, rcoll

*-----------------------------------------------------------------------
* (2)       SCALE IN
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' scale'')')
     &                           pxys(ib,4), pxys(ib,5)

                     sxg0 = sxg0 * pxys(ib,4)
                     syg0 = syg0 * pxys(ib,5)

*-----------------------------------------------------------------------
*           RIBBON PART(2)
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                        write(jhf,'(2g14.5,'' M'')')
     &                  -6.0*ru,  1.35*ru

                        write(jhf,'(4g14.5,/2g14.5,'' curveto'')')
     &                  -6.0*ru,  1.95*ru, -4.5*ru, 1.15*ru,
     &                  -4.5*ru,  1.75*ru

                        write(jhf,'(2g14.5,'' L'')')
     &                  -4.5*ru,  -0.25*ru

                        write(jhf,'(4g14.5,/2g14.5,'' curveto cp'')')
     &                  -4.5*ru, -0.85*ru, -6.0*ru,  -0.05*ru,
     &                  -6.0*ru, -0.65*ru

                        write(jhf,'(2g14.5,'' M'')')
     &                  6.0*ru,  1.35*ru

                        write(jhf,'(4g14.5,/2g14.5,'' curveto'')')
     &                  6.0*ru,  1.95*ru,  4.5*ru, 1.15*ru,
     &                  4.5*ru,  1.75*ru

                        write(jhf,'(2g14.5,'' L'')')
     &                  4.5*ru,  -0.25*ru

                        write(jhf,'(4g14.5,/2g14.5,'' curveto cp'')')
     &                  4.5*ru, -0.85*ru, 6.0*ru,  -0.05*ru,
     &                  6.0*ru, -0.65*ru

*-----------------------------------------------------------------------
*           SCALE OUT
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' scale'')')
     &                  1./pxys(ib,4), 1./pxys(ib,5)

                     sxg0 = sxg0 / pxys(ib,4)
                     syg0 = syg0 / pxys(ib,5)

*-----------------------------------------------------------------------
*           FILL AND STRING RIBBON PATH
*-----------------------------------------------------------------------

                  write(jhf,'(''gs '',3f7.3,'' sc fl gr sd0 '',i2,
     &                        '' dd lw '',3f7.3,'' sc st'')')
     &                       rcols, iwi, rcoll

*-----------------------------------------------------------------------
* (3)       SCALE IN
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' scale'')')
     &                           pxys(ib,4), pxys(ib,5)

                     sxg0 = sxg0 * pxys(ib,4)
                     syg0 = syg0 * pxys(ib,5)

*-----------------------------------------------------------------------
*           RIBBON PART(3)
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                        ru = unit0 / 2.0

                        write(jhf,'(2g14.5,'' M'')')
     &                  0.0*ru,  -1.0*ru

                        write(jhf,'(4g14.5,/2g14.5,'' curveto'')')
     &                  -2.0*ru, -ru, -4.5*ru, -2.0*ru,
     &                  -6.0*ru, -0.65*ru

                        call bbox(1,0,-4.5*ru,-1.70*ru,0.d0)
                        call bbox(1,0,-6.0*ru,-0.65*ru,0.d0)

                        write(jhf,'(2g14.5,'' L'')')
     &                  -6.0*ru,  1.35*ru

                        write(jhf,'(4g14.5,/2g14.5,'' curveto'')')
     &                  -4.5*ru, 0.0*ru, -2.0*ru,  ru,
     &                  0.0*ru,  1.0*ru

                        write(jhf,'(4g14.5,/2g14.5,'' curveto'')')
     &                  2.0*ru, ru, 4.5*ru, 0.0*ru,
     &                  6.0*ru, 1.35*ru

                        write(jhf,'(2g14.5,'' L'')')
     &                  6.0*ru,  -0.65*ru

                        write(jhf,'(4g14.5,/2g14.5,'' curveto cp'')')
     &                  4.5*ru, -2.0*ru,  2.0*ru,  -1.0*ru,
     &                  0.0*ru, -1.0*ru

                        call bbox(1,0,4.5*ru,-1.70*ru,0.d0)
                        call bbox(1,0,6.0*ru,-0.65*ru,0.d0)


*-----------------------------------------------------------------------
*           SCALE OUT
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' scale'')')
     &                  1./pxys(ib,4), 1./pxys(ib,5)

                     sxg0 = sxg0 / pxys(ib,4)
                     syg0 = syg0 / pxys(ib,5)

*-----------------------------------------------------------------------
*           FILL AND STRING RIBBON PATH
*-----------------------------------------------------------------------

                  write(jhf,'(''gs '',3f7.3,'' sc fl gr sd0 '',i2,
     &                        '' dd lw '',3f7.3,'' sc st'')')
     &                       rcolb, iwi, rcoll

*-----------------------------------------------------------------------
*           TRANSLATE AND ROTATE
*-----------------------------------------------------------------------

                  if( iang .eq. 0 ) then

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' TR'')') -xpsain, -ypsain

                     call bbox(1,1,-xpsain,-ypsain,0.d0)

                  else

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(g14.5,'' rotate'')') -angain

                     write(jhf,'(2g14.5,'' TR'')') -xpsain, -ypsain

                     call bbox(1,1,0.d0,0.d0,-angain)
                     call bbox(1,1,-xpsain,-ypsain,0.d0)

                  end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine warrw02(jhf,idbg,clal,clmo,rcoll,rcolb,
     &                   x1,y1,x2,y2,xl,yl,dd,rnl,ara)
*                                                                      *
*                                                                      *
*        PURPOSE    :  WRITE ARROWS ON jhf                             *
*                                                                      *
*       X1,Y1,X2,Y2 : COORDINATE OF THE INITIAL AND FINAL POINT        *
*             RCOLL : COLOR OF ARROWS LINE                             *
*             RCOLB : COLOR OF ARROWS INTERIA                          *
*                DD : WIDTH OF ARROW LINE                              *
*               RNL : -1 -> NO ARROW LINE, 1 -> WITH ARROW LINE        *
*               ARA : ANGLE OF THE ARROW RELATIVE TO DEFAULT           *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter (ddcm = 8.4666667e-3)

      parameter ( rang0 =  2.0 )
      parameter ( rlng  = 12.0 )

      common /frm/  xal, yal

      dimension clal(3)
      dimension rcoll(3), rcolb(3)

*-----------------------------------------------------------------------

            rang = rang0 / ara

            ahed = ( rang + 1.0 ) / 2.0 * dd * ddcm
            aang = rang * dd * ddcm
            alng = rlng * dd * ddcm
            aalg = alng / rang

            x1 = x1 * xal
            x2 = x2 * xal
            y1 = y1 * yal
            y2 = y2 * yal

            rleng = sqrt( ( x2 - x1 )**2 + ( y2 - y1 )**2 )

            if( rleng .lt. alng ) rnl = -1.0

*-----------------------------------------------------------------------

         if( rnl .gt. 0.0 ) then

               x3 = x2 + ahed / rleng * ( x1 - x2 )
               y3 = y2 + ahed / rleng * ( y1 - y2 )

               x6 = x2 + alng / rleng * ( x1 - x2 )
               y6 = y2 + alng / rleng * ( y1 - y2 )

         else

               x3 = x1 + 2.0 / 3.0 * alng / rleng * ( x2 - x1 )
               y3 = y1 + 2.0 / 3.0 * alng / rleng * ( y2 - y1 )

               x6 = x1 - 1.0 / 3.0 * alng / rleng * ( x2 - x1 )
               y6 = y1 - 1.0 / 3.0 * alng / rleng * ( y2 - y1 )

               x2 = x3
               y2 = y3

         end if


         if( x2 .eq. x6 ) then

               x4 = x6 - aalg
               x5 = x6 + aalg
               y4 = y6
               y5 = y6

               x7  = x4 + ( x5 - x4 ) * 1.0 / 3.0
               x10 = x4 + ( x5 - x4 ) * 2.0 / 3.0
               y7  = y6
               y10 = y6

               x8  = x7
               x9  = x10
               y8  = y1
               y9  = y1

         else if( y2 .eq. y6 ) then

               y4 = y6 - aalg
               y5 = y6 + aalg
               x4 = x6
               x5 = x6

               y7  = y4 + ( y5 - y4 ) * 1.0 / 3.0
               y10 = y4 + ( y5 - y4 ) * 2.0 / 3.0
               x7  = x6
               x10 = x6

               x8  = x1
               x9  = x1
               y8  = y7
               y9  = y10

         else

               yxf = ( y2 - y6 ) / ( x2 - x6 )
               yyf = sqrt( yxf**2 + 1.0 )
               xxf = sqrt( 1.0 / yxf**2 + 1.0 )

            if( yxf .gt. 0.0 ) then

               x4 = x6 - aalg / xxf
               y4 = y6 + aalg / yyf
               x5 = x6 + aalg / xxf
               y5 = y6 - aalg / yyf

            else

               x4 = x6 + aalg / xxf
               y4 = y6 + aalg / yyf
               x5 = x6 - aalg / xxf
               y5 = y6 - aalg / yyf

            end if

            if( rnl .gt. 0.0 ) then

               x7  = x4 + ( x5 - x4 ) * 1.0 / 3.0
               y7  = y4 + ( y5 - y4 ) * 1.0 / 3.0

               x10 = x4 + ( x5 - x4 ) * 2.0 / 3.0
               y10 = y4 + ( y5 - y4 ) * 2.0 / 3.0

               a26 = ( y2 - y6 ) / ( x2 - x6 )

               a78 = a26
               b78 = y7 - a78 * x7

               a90 = a26
               b90 = y10 - a90 * x10

               a89 = ( y4 - y5 ) / ( x4 - x5 )
               b89 = y1 - a89 * x1

               x8 = ( b89 - b78 ) / ( a78 - a89 )
               y8 = a89 * x8 + b89

               x9 = ( b89 - b90 ) / ( a90 - a89 )
               y9 = a89 * x9 + b89

            end if

         end if

               x1 = x1 / xal
               x2 = x2 / xal
               x3 = x3 / xal
               x4 = x4 / xal
               x5 = x5 / xal

               y1 = y1 / yal
               y2 = y2 / yal
               y3 = y3 / yal
               y4 = y4 / yal
               y5 = y5 / yal

            if( rnl .gt. 0.0 ) then

               x7  = x7  / xal
               x8  = x8  / xal
               x9  = x9  / xal
               x10 = x10 / xal

               y7  = y7  / yal
               y8  = y8  / yal
               y9  = y9  / yal
               y10 = y10 / yal

            end if

*-----------------------------------------------------------------------

            iwi = nint( dd / 2.0 )

               if( clal(1) .gt. -r0max ) then

                     rcoll(1) = clal(1)
                     rcoll(2) = clal(2)
                     rcoll(3) = clal(3)

                  if( rcolb(1) .lt. -r0max ) then

                     rcolb(1) = -3.0

                  else

                     rcolb(1) = clal(1)
                     rcolb(2) = clal(2)
                     rcolb(3) = clal(3)

                  end if

               end if

               if( clmo .gt. -r0max .and. rcoll(1) .gt. 0.0 )
     &             rcoll(1) = -2.0
               if( clmo .gt. -r0max .and. rcolb(1) .gt. 0.0 )
     &             rcolb(1) = -1.0

               if( rcoll(1) .lt. -r0max ) rcoll(1) = -2.0
               if( rcolb(1) .lt. -r0max ) rcolb(1) = -3.0

               if( rcolb(1) .gt. -2.5 ) then

                  write(jhf,'(3f7.3,'' sc sd0 '',i2,'' dd lw'')')
     &                       rcolb, iwi

               else

                  write(jhf,'(''sd0 '',i2,'' dd lw'')')
     &                       iwi

               end if

*-----------------------------------------------------------------------

         if( rnl .gt. 0.0 ) then

               write(jhf,'(2g14.5,'' m '',2g14.5,'' l '',/
     &                     2g14.5,'' l'')') x2, y2, x4, y4, x7, y7
               write(jhf,'(2g14.5,'' l '',2g14.5,'' l '',/
     &                     2g14.5,'' l'')') x8, y8, x9, y9, x10, y10

            if( rcolb(1) .gt. -2.5 ) then

               write(jhf,'(2g14.5,'' l cp gs fl gr '',3f7.3,
     &                     '' sc st'')') x5, y5, rcoll

            else

               write(jhf,'(2g14.5,'' l cp '',3f7.3,
     &                     '' sc st'')') x5, y5, rcoll

            end if

                  call bbox(3,0,x2,y2,0.d0)
                  call bbox(3,0,x4,y4,0.d0)
                  call bbox(3,0,x5,y5,0.d0)
                  call bbox(3,0,x8,y8,0.d0)
                  call bbox(3,0,x9,y9,0.d0)

         else

               write(jhf,'(2g14.5,'' m '',
     &                     2g14.5,'' l'')') x2, y2, x4, y4

            if( rcolb(1) .gt. -2.5 ) then

               write(jhf,'(2g14.5,'' l cp gs fl gr '',3f7.3,
     &                     '' sc st'')') x5, y5, rcoll

            else

               write(jhf,'(2g14.5,'' l cp '',3f7.3,
     &                     '' sc st'')') x5, y5, rcoll

            end if

                  call bbox(3,0,x2,y2,0.d0)
                  call bbox(3,0,x4,y4,0.d0)
                  call bbox(3,0,x5,y5,0.d0)

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine pkan(dum,ic,icf,idat,ierr)
*                                                                      *
*        PURPOSE : JUDGE THE LANGUAGE FOR DATE FROM (   )              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character chaf(ichrl)*1

*-----------------------------------------------------------------------

            prn  = 0.0
            ierr = 1

            if( dum(ic) .ne. '(' ) return

         ici = ic + 1
  100    ic  = ic + 1

            if( ic .gt. icf ) return

            if( dum(ic) .eq. ')' ) then

               icn = ic - ici

               do 200 i = 1, icn

                     j = ici + i - 1

                     chaf(i) = dum(j)

  200          continue

                     call jpncode(chaf,icn,idat)

                     ierr = 0

                     return

            end if

            goto 100


      end


************************************************************************
*                                                                      *
      subroutine a_fact(lum,iclm,ici,icf,dfac,ierr)
*                                                                      *
*       PURPOSE:       READ FACTOR (CHARACTERS TO REAL NUMBER)         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character lum(ichrl)*1
      character dmm*1, dm1
      dimension dfac(5)

      logical dnen2, dnen1, dnen3, deqn2

*-----------------------------------------------------------------------
*     INITIAL VALUES
*-----------------------------------------------------------------------

               xnd = 1.0
               xnt = 1.0
               xnm = 0.0
               xnp = 0.0
               xnu = -r1max

               ic = ici
               ip = 0

*-----------------------------------------------------------------------

         if( lum(ic) .eq. '*' .and. lum(ic+1) .eq. '*' ) then

               ic = ic + 2

               if( lum(ic) .eq. '(' ) then

                  ic = ic + 1
                  ip = 1

               end if

                  dm1 = lum(ic)

                  if( dnen2(dm1) ) goto 999

                  if( ic .eq. icf .and.
     &              ( dm1 .eq. '.' .or.
     &                dm1 .eq. '+' .or.
     &                dm1 .eq. '-' ) ) goto 999

            icn = ic

  246       ic = ic + 1

               if( ic .gt. icf ) goto 250

               dm1 = lum(ic)

               if( dm1 .ne. '*' .and. dm1 .ne. '/' .and.
     &             dnen3(dm1) .and.
     &           ( ip .eq. 0 .and. dm1 .eq. ')' ) ) goto 999


               if( deqn2(dm1) ) goto 246

               if( ( dm1 .eq. '+' .or. dm1 .eq. '-' ) .and.
     &               lum(ic-1) .eq. 'e' ) goto 246

  250          call rnum(xnu,lum,icn,ic-1,ierr)

                  if( ierr .ne. 0 ) goto 999

               if( ip .eq. 1 ) then

                  if( lum(ic) .ne. ')' ) goto 999
                  ic = ic + 1

               end if

               if( ic .gt. icf ) goto 800

         end if

*-----------------------------------------------------------------------

               ip = 0

               dmm = lum(ic)

               if( dmm .ne. '+' .and. dmm .ne. '-' .and.
     &             dmm .ne. '*' .and. dmm .ne. '/' ) goto 999

               ic = ic + 1

                  if( lum(ic) .eq. '(' ) then
                      ic = ic + 1
                      ip = 1
                  end if

               dm1 = lum(ic)

                   if( ( ip .eq. 0 .and. dnen1(dm1) .and.
     &                   dm1 .ne. '.' ) .or.
     &                 ( ip .eq. 1 .and. dnen2(dm1) ) ) goto 999

                   if( ic .eq. icf .and. dm1 .eq. '.' ) goto 999

               icn = ic

  146          ic = ic + 1

                  if( ic .gt. icf ) goto 150

               dm1 = lum(ic)

                  if( dnen3(dm1) .and.
     &              ( ip .eq. 0 .and. dm1 .eq. ')' ) ) goto 999

                  if( deqn2(dm1) ) goto 146

                  if( ( dm1 .eq. '+' .or. dm1 .eq. '-' ) .and.
     &                  lum(ic-1) .eq. 'e' ) goto 146

                  if( ( ( dm1 .eq. '+' .or. dm1 .eq. '-' ) .and.
     &                    lum(ic-1) .ne. 'e' ) .and.
     &                  ( dmm .eq. '+' .or. dmm .eq. '-' ) )  goto 999

*-----------------------------------------------------------------------

  150    continue

*-----------------------------------------------------------------------

            if( dmm .eq. '+' .or. dmm .eq. '-' ) then

               if( dmm .eq. '+' ) then

                  call rnum(xnp,lum,icn,ic-1,ierr)

                     if( ierr .ne. 0 ) goto 999

               end if

               if( dmm .eq. '-' ) then

                  call rnum(xnm,lum,icn,ic-1,ierr)

                     if( ierr .ne. 0 ) goto 999

               end if

               if( ip .eq. 1 ) then

                  if( lum(ic) .ne. ')' ) goto 999
                  ic =ic + 1

               end if

            end if

*-----------------------------------------------------------------------

         if( dmm .eq. '/' .or. dmm. eq. '*' ) then

                  if( dmm .eq. '/' ) then

                     call rnum(xnd,lum,icn,ic-1,ierr)

                        if( ierr .ne. 0 ) goto 999

                  end if

                  if( dmm .eq. '*' ) then

                     call rnum(xnt,lum,icn,ic-1,ierr)

                        if( ierr .ne. 0 ) goto 999

                  end if

                  if( ip .eq. 1 ) then

                     if( lum(ic) .ne. ')' ) goto 999
                     ic = ic + 1

                  end if


            if( ic .le. icf ) then

                  if( ic .eq. icf ) goto 999

                  dmm = lum(ic)

                  if( dmm .ne. '+' .and. dmm .ne. '-' ) goto 999

                  ic = ic + 1
                  ip = 0

                  if( lum(ic) .eq. '(' ) then

                     ic = ic + 1
                     ip = 1

                  end if

                  if( ic .gt. iclm - 2 ) goto 999

                  dm1 = lum(ic)

                  if( ( ip .eq. 0 .and.
     &                  dnen1(dm1) .and. dm1 .ne. '.' ) .or.
     &                ( ip .eq. 1 .and. dnen2(dm1) ) )  goto 999

                  if( ic .eq. icf .and. dm1 .eq. '.' ) goto 999

                  icn = ic

  147             ic = ic + 1

                     if( ic .gt. icf ) goto 160

                  dm1 = lum(ic)

                  if( dnen3(dm1) .and.
     &              ( ip .eq. 0 .and. dm1 .eq. ')' ) )  goto 999

                  if( deqn2(dm1) ) goto 147

                  if( dm1 .eq. '+' .or. dm1 .eq. '-' ) then

                     if( lum(ic-1) .eq. 'e' ) goto 147

                     goto 999

                  end if

  160          continue

                  if( dmm .eq. '+' ) then

                     call rnum(xnp,lum,icn,ic-1,ierr)

                        if( ierr .ne. 0 ) goto 999

                  end if

                  if( dmm .eq. '-' ) then

                     call rnum(xnm,lum,icn,ic-1,ierr)

                        if( ierr .ne. 0 ) goto 999

                  end if

                  if( ip .eq. 1 ) then

                     if( lum(ic) .ne. ')' ) goto 999
                     ic = ic + 1

                  end if

            end if

         end if

*-----------------------------------------------------------------------

  800                dfac(1) = xnd
                     dfac(2) = xnt
                     dfac(3) = xnm
                     dfac(4) = xnp
                     dfac(5) = xnu

       return

*-----------------------------------------------------------------------

  999  ierr=1

       return
       end


************************************************************************
*                                                                      *
      subroutine atitle(ill,jsi,jsn,dum,lum,in1,rtt,itt,ifon,iend)
*                                                                      *
*        PURPOSE;                                                      *
*                      READ TITLE COMMENTS                             *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      character dumn(ichrl)*1
      character rtt(ichrl)*1

      dimension ill(0:9)

      character yen*1
      character tub*1
      tub = char(9)
      yen = char(92)

*-----------------------------------------------------------------------

         iend = 0
         in   = in1


         ici = 0
         icf = 0

         do 100 i = in, icolm

            if(dum(i).ne.' '.and.dum(i).ne.tub) goto 110

  100    continue

            goto 300

  110    ici = i

         do 200 i = icolm, in, -1

            if(dum(i).ne.' '.and.dum(i).ne.tub) goto 210

  200    continue

            goto 300

  210    icf = i

         itt = icf - ici + 1

         do 400 i = 1, itt

            rtt(i) = dum(ici+i-1)

  400    continue

*-----------------------------------------------------------------------

                  iseql = 0
               if( icf .eq. 1 .and. dum(icf) .eq. yen )
     &            iseql = 1
               if( icf .gt. 1 ) then
                  if( dum(icf) .eq. yen .and. dum(icf-1) .ne. yen .and.
     &                                  ichar(dum(icf-1)) .le. 128 )
     &            iseql = 2
               end if

      if( iseql .gt. 0 ) then

  360    ill(jsn) = ill(jsn) + 1

         read(jsi,'(10000a1)', iostat = ios ) (dum(ic),ic=1,icolm)
         if( ios .eq. -1 ) goto 350

         icf = 0

         do 201 i = icolm, 1, -1

            if(dum(i).ne.' '.and.dum(i).ne.tub) goto 211

  201    continue

            goto 351

  211    icf = i

         do 401 i = 1, icf

            rtt(itt-1+i) = dum(i)

  401    continue

            itt = itt + icf - 1

                  jseql = 0
               if( icf .eq. 1 .and. dum(icf) .eq. yen )
     &            jseql = 1
               if( icf .gt. 1 ) then
                  if( dum(icf) .eq. yen .and. dum(icf-1) .ne. yen .and.
     &                                  ichar(dum(icf-1)) .le. 128 )
     &            jseql = 2
               end if

               if( jseql .gt. 0 ) goto 360

      end if

*-----------------------------------------------------------------------

               goto 351
  350          iend = 1
  351          continue

            call jpncode(rtt,itt,ifon)
            call jpnprep(rtt,itt,ifon)

            return

*-----------------------------------------------------------------------

  300    itt = 0

      return
      end

************************************************************************
cKN 2024/01/24
*                                                                      *
      subroutine xyrmax(xmin,xmax,ixlog,ispac,ierr,
     &                  xpmin,xpmax)
*                                                                      *
*       PURPOSE : DETERMINE MAX AND MIN VALUE OF AXIS                  *
*                                                                      *
*       AXPAR ; % OF (XMAX-XMIN) TO (XMIN-0.0)  ALL>0                  *
*       AXSPC ; % OF EXTRA SPACE OF THE AXIS                           *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

        parameter( axpar = 0.26,  axspc = 0.07 )

*-----------------------------------------------------------------------
*          X-AXIS LINEAR CASE
*-----------------------------------------------------------------------

           if( ixlog .eq. 0 ) then

               xd = xmax - xmin

                 if( xmin .ge. 0.0 ) then

                    if( ispac .eq. 0 ) then

                       if( xmin / xd .lt. axpar ) then

                          xmin = 0.0
                          xmax = xmax * ( 1.0 + axspc )

                       else

                          xmin = xmin - xd * axspc
                          xmax = xmax + xd * axspc

                          if( xmin .lt. 0.0 ) xmin = 0.0

                       end if

                    end if

                 else

                    if( xmax .gt. 0.0 ) then

                       if( ispac .eq . 0 ) then

                             xmin = xmin - xd * axspc
                             xmax = xmax + xd * axspc

                       end if

                    else

                       if( ispac .eq. 0 ) then

                          if( xmax / ( xmin - xmax ) .lt. axpar ) then

                             xmax = 0.0
                             xmin = xmin * ( 1.0 + axspc )

                           else

                             xmin = xmin - xd * axspc
                             xmax = xmax + xd * axspc

                           end if

                       end if

                    end if

                 end if

*-----------------------------------------------------------------------
*          X-AXIS LOG CASE
*-----------------------------------------------------------------------

           else

              if( xmin .lt. 0.0 ) goto 999

              if( xmin .eq. 0.0 ) xmin = r9min

              if( ispac .eq. 0 ) then

                 xd = log10(xmax) - log10(xmin)

                 xmin = 10.0**( log10(xmin) - xd * axspc )

                 xmax = 10.0**( log10(xmax) + xd * axspc )

              end if

*-----------------------------------------------------------------------
cKN 2024/01/24

               if( xpmin .ge.  r0max .and.
     &             xpmax .le. -r0max ) then

                  if( xd .lt. 1.0 ) then

                     xmin = xmin / 5.0
                     xmax = xmax * 5.0

                  end if

               end if

*-----------------------------------------------------------------------

           end if

      return

*-----------------------------------------------------------------------

  999 ierr = 1


      return
      end


************************************************************************
*                                                                      *
      subroutine xytic(ixlog,xxmin,xxmax,ixtic,xld,xsd,
     &              nxlta,nxsta,xlta,xsta,xcta,nxtch,
     &              idecx,ixstd,ixltd,mdecx,xstv,xltv,ierr,
     &              xxmul)
*                                                                      *
*           purpose : determine tic point and tic number               *
*                                                                      *
*            ixlog       : 1 log scale, 0 linear scale                 *
*            xmin, xmax  : min and max values of axis                  *
*            ixtic       : 1 increase long tics,                       *
*                          0 default, -1 decrease                      *
*            xld, xsd    : mesh of short and long tics                 *
*            nxlta, nxsta: number of long or short tics                *
*            xlta, xsta  : position of long or short tics              *
*            xcta        : character of tics number                    *
*            nxtch       : number of character of each tics            *
*            idecx       : number of decimals of axis                  *
*            ixstd       : =0 auto, =1 manual                          *
*            ixltd       : =0 auto, =1 manual                          *
*            mdecx       : manual number of decimals of axis           *
*                          if mdecx=-10 auto                           *
*            xstv, xltv  : explicit one value of the tics              *
*                                                                      *
*                                                                      *
*    ####    ic1, ic2, ic3 : critical values of range    ####          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter ( ic1 = 15, ic2 = 30, ic3 = 60 )
      parameter ( snb = 1.0e-5 )

*-----------------------------------------------------------------------

      character dum(30)*1,d1*1,d2*2,d3*3,d4*4
      character rdum*100

      dimension xlta(numtic), xsta(numtic), nxtch(numtic)
      character xcta(numtic,30)*1

*----------------------------------------------------------------------*
                             xmin = xxmin
      if( xxmul .ne. 0.0d0 ) xmin = xxmin * xxmul
                             xmax = xxmax
      if( xxmul .ne. 0.0d0 ) xmax = xxmax * xxmul

*----------------------------------------------------------------------*

      nxlta = 0
      nxsta = 0

      idecx = -1

      xd = xmax - xmin


*-----------------------------------------------------------------------
*     LOG SCALE
*-----------------------------------------------------------------------

      if( ixlog .eq. 1 ) then

         write(rdum,'(E11.4e3)') xmin
cKN 2024/01/03
c        write(rdum,'(E10.4)') xmin

         d1=rdum(3:3)
         read(d1,'(I1)') ixd

c        write(*,'(e11.4e3)') xmin
c        write(*,'(e10.4)') xmin

cKN 2024/01/03
         d4=rdum(8:8)//rdum(9:9)//rdum(10:10)//rdum(11:11)
         read(d4,'(I4)') iee

         gmax = log10(xmax)
         gmin = log10(xmin)
         gd   = gmax - gmin

*----------------------------------------------------------------------*

         if( ixstd .eq. 0 .or. ( ixstd .eq. 1 .and. xsd .ge. 0.0 )) then

                  xsmin = dble(ixd) * 10.**(iee-1)

               if( xsmin .lt. xmin ) then

                  ixd = ixd + 1

               end if

                  idd = iee - 1

                  if( ixd .eq. 10 ) then

                     ixd = 1
                     idd = idd + 1

                  end if

            do 320 n = 1, numtic

               xstva =  dble( ixd ) * 10.d0**dble(idd)
               xstan = ( log10( xstva ) - gmin ) / gd

               if( xstan .gt. 1.00001 ) goto 330

               xsta(n) = xstan

                  ixd = ixd + 1

                  if( ixd .eq. 10 ) then

                     ixd = 1
                     idd = idd + 1

                  end if

  320       continue

  330          nxsta = n - 1

         end if

*----------------------------------------------------------------------*

         if( ixltd .eq. 0 .or. ( ixltd .eq. 1 .and. xld .ge. 0.0 )) then

               xlmin = 10.**dble(iee-1)
               iemin = iee - 1

               if( xlmin .lt. xmin ) then

                  xlmin = 10.**dble(iee)
                  iemin = iee

               end if

               xlmin = ( log10( xlmin ) - gmin ) / gd


            do 300 n = 1, numtic

               xltan = dble( n - 1 ) / gd + xlmin
               if( xltan .gt. 1.00001 ) goto 310

               xlta(n) = xltan

*-----------------------------------------------------------------------
*           Get Tic Number
*-----------------------------------------------------------------------

               ienum = ( n - 1 ) + iemin

               write(rdum,'(I4)') ienum

                  j = 0

               do 500 i = 1, 4

                  if( rdum(i:i) .eq. ' ' ) goto 500

                  j = j + 1

                  xcta(n,j) = rdum(i:i)

  500          continue

                  nxtch(n) = j

cKN 2024/01/24 w
c            write(*,*) (xcta(n,j),j=1,nxtch(n)), idecx, nxtch(n)

*-----------------------------------------------------------------------

  300       continue

  310          nxlta = n - 1

         end if


*----------------------------------------------------------------------*
*     LINEAR SCAL
*----------------------------------------------------------------------*

      else

         write(rdum,'(E11.4e3)') xd
cKN 2024/01/03
c        write(rdum,'(E10.4)') xd

         d2=rdum(3:3)//rdum(4:4)
         read(d2,'(I2)') ixd
         if(ixd.lt.10) goto 999

c       write(*,'(e11.4e3)') xd
c       write(*,'(e10.4)') xd

cKN 2024/01/03
c        d3=rdum(8:8)//rdum(9:9)//rdum(10:10)
c        read(d3,'(I3)') iee
         d4=rdum(8:8)//rdum(9:9)//rdum(10:10)//rdum(11:11)
         read(d4,'(I4)') iee

*----------------------------------------------------------------------*
*                                                                      *
*  tic formation      1993/05/21  by koji niita                        *
*                                                                      *
*                                                                      *
*     10 < ixd < 15        ic1 = 15                                    *
*                                                                      *
*                  default   xld = 2,  xsd = 1                         *
*                    small   xld = 1,  xsd = 0.5                       *
*                    large   xld = 5,  xsd = 1                         *
*                                                                      *
*     15 < ixd < 30        ic2 = 30                                    *
*                                                                      *
*                  default   xld = 5,  xsd = 1                         *
*                    small   xld = 2,  xsd = 1                         *
*                    large   xld = 10, xsd = 5                         *
*                                                                      *
*     30 < ixd < 60        ic3 = 60                                    *
*                                                                      *
*                  default   xld = 10, xsd = 5                         *
*                    small   xld = 5,  xsd = 1                         *
*                    large   xld = 20, xsd = 10                        *
*                                                                      *
*     60 < ixd < 99        ic3 = 60                                    *
*                                                                      *
*                  default   xld = 20, xsd = 10                        *
*                    small   xld = 10, xsd = 5                         *
*                    large   xld = 50, xsd = 10                        *
*                                                                      *
*----------------------------------------------------------------------*

cKN 2024/12/3

      if(ixd.lt.ic1) then
         if(ixtic.eq.0) then

            xlda = 2.0d0 * 10.d0**(iee-2)
            xsda = 1.0d0 * 10.d0**(iee-2)

            if(iee.lt.2) idecx=iabs(iee-2)

         else if(ixtic.eq.-1) then

            xlda = 1.0d0 * 10.d0**(iee-2)
            xsda = 0.5d0 * 10.d0**(iee-2)

            if(iee.lt.2) idecx=iabs(iee-2)

         else if(ixtic.eq.1 ) then

            xlda = 5.0d0 * 10.d0**(iee-2)
            xsda = 1.0d0 * 10.d0**(iee-2)

            if(iee.lt.2) idecx=iabs(iee-2)

         end if
      end if

      if(ixd.ge.ic1.and.ixd.lt.ic2) then
         if(ixtic.eq.0) then

            xlda = 5.0d0 * 10.d0**(iee-2)
            xsda = 1.0d0 * 10.d0**(iee-2)

            if(iee.lt.2) idecx=iabs(iee-2)

         else if(ixtic.eq.-1) then

            xlda = 2.0d0 * 10.d0**(iee-2)
            xsda = 1.0d0 * 10.d0**(iee-2)

            if(iee.lt.2) idecx=iabs(iee-2)

         else if(ixtic.eq.1 ) then

            xlda = 10.0d0 * 10.d0**(iee-2)
            xsda = 5.0d0 * 10.d0**(iee-2)

            if(iee.lt.1) idecx=iabs(iee-1)

         end if
      end if

      if(ixd.ge.ic2.and.ixd.lt.ic3) then
         if(ixtic.eq.0) then

            xlda = 10.0d0 * 10.d0**(iee-2)
            xsda = 5.0d0 * 10.d0**(iee-2)

            if(iee.lt.1) idecx=iabs(iee-1)

         else if(ixtic.eq.-1) then

            xlda = 5.0d0 * 10.d0**(iee-2)
            xsda = 1.0d0 * 10.d0**(iee-2)

            if(iee.lt.2) idecx=iabs(iee-2)

         else if(ixtic.eq.1 ) then

            xlda = 20.0d0 * 10.d0**(iee-2)
            xsda = 10.0d0 * 10.d0**(iee-2)

            if(iee.lt.1) idecx=iabs(iee-1)

         end if
      end if

      if(ixd.ge.ic3) then
         if(ixtic.eq.0) then

            xlda = 20.0d0 * 10.d0**(iee-2)
            xsda = 10.0d0 * 10.d0**(iee-2)

            if(iee.lt.1) idecx=iabs(iee-1)

         else if(ixtic.eq.-1) then

            xlda = 10.0d0 * 10.d0**(iee-2)
            xsda = 5.0d0 * 10.d0**(iee-2)

            if(iee.lt.1) idecx=iabs(iee-1)

         else if(ixtic.eq.1 ) then

            xlda = 50.0d0 * 10.d0**(iee-2)
            xsda = 10.0d0 * 10.d0**(iee-2)

            if(iee.lt.1) idecx=iabs(iee-1)

         end if
      end if

cKN 2024/12/3

*-----------------------------------------------------------------------


      if(ixltd.eq.0)    xld = xlda
      if(ixstd.eq.0)    xsd = xsda

      if(mdecx.ne.-10)  idecx = mdecx

*-----------------------------------------------------------------------

      xlmi=real(nint((xmin-xltv)/xld))*xld+xltv
           if(xlmi.lt.xmin-xld*snb) xlmi=xlmi+xld

      xsmi=real(int((xmin-xstv)/xsd))*xsd+xstv
           if(xsmi.lt.xmin-xsd*snb) xsmi=xsmi+xsd

*-----------------------------------------------------------------------

      if( xld .gt. 0.0 ) then

         do 100 n = 1, numtic

            xltan = ( ( n - 1 ) * xld + xlmi - xmin ) / xd
            if( xltan .gt. 1.00001 ) goto 101

            xlta(n) = xltan

*-----------------------------------------------------------------------
*           Get Tic Number
*-----------------------------------------------------------------------

            xctan = ( n - 1 ) * xld + xlmi


            if( idecx .gt. 0 ) then

               xctan = ( dble( nint( xctan * 10.**idecx ) )
     &                 + sign( 0.1d0, xctan ) ) / 10.**idecx

            else

               xctan = dble( nint( xctan , 8) ) + sign( 0.1d0, xctan )

            end if


cKN 2024/01/24 --->       ***
            write(rdum,'(f100.50)') xctan
c           write(rdum,'(f26.10)') xctan
c           write(rdum,'(f22.10)') xctan
c           read(rdum,'(22a1)') (dum(i),i=1,22)

            j  = 0
            k  = 0
            kk = -1

cKN 2024/01/24 <---       ***
            do 400 i = 1, 100

               if( rdum(i:i) .eq. ' ' ) goto 400

               j = j + 1

               if( rdum(i:i) .eq. '.' ) then

                  k = 1

                  if( i .gt. 1 ) then
                     if( rdum(i-1:i-1) .eq. ' ' .or.
     &                   rdum(i-1:i-1) .eq. '-' .or.
     &                   rdum(i-1:i-1) .eq. '+' ) then

                        xcta(n,j) = '0'
                        j = j + 1

                     end if
                  else if( i .eq. 1 ) then

                        xcta(n,j) = '0'
                        j = j + 1

                  end if

               end if

               if( k .eq. 1 ) kk = kk + 1

               if( kk .gt. idecx ) goto 410

                  xcta(n,j) = rdum(i:i)

  400       continue

  410       nxtch(n) = j - 1

cKN 2024/01/24 w
c             write(*,*) (xcta(n,j),j=1,nxtch(n)), idecx, nxtch(n)

*-----------------------------------------------------------------------


  100    continue

  101       nxlta = n - 1

      end if

*-----------------------------------------------------------------------

      if( xsd .gt. 0.0 ) then

         do 200 n = 1, numtic

            xstan = ( ( n - 1 ) * xsd + xsmi - xmin ) / xd
            if( xstan .gt. 1.00001 ) goto 201

            xsta(n) = xstan

  200    continue

  201       nxsta = n - 1

      end if



*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

      return
  999 ierr=1
      return
      end


************************************************************************
*                                                                      *
      subroutine contjol(jof,icut,ncut,rcut,icct,iccr,
     &                   iwd2,ipd2,noned,izlog,ibmap,
     &                   daxy,ixn,iyn,dax,day,ierr)
*                                                                      *
*      PURPOSE  :  WRITE CONTOUR ON JOL                                *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      dimension rcut(mc,5)
      dimension iwd2(mc)

      common /h2/ xdd,ydd,
     &            xymax,xymin,clus,cmax,cmin,cxymax,cxymin,xyabs,
     &            dmax,dmin,ixnm,iynm

      common /h3/ smax(3), smin(3)

      dimension dax(ixn), day(ixn)

      dimension daxy(ixn,iyn)
      dimension itrace(ixn,iyn)

*-----------------------------------------------------------------------
*           WRITE CONTOUR ON JOF
*-----------------------------------------------------------------------

*           ICUT :  NUMBER OF CUTS
*           IWD2 ;  WIDTH OF LINES (4, 7, 10, 13)

*-----------------------------------------------------------------------

            ierr = 0

                  if( icut .eq. 0 ) then

                     icut = ncut
                     dcut = ( xymax - xymin ) / dble(icut+1)

                     do 510 i = 1, icut

                        rcut(i,1) = xymin + dble(i) * dcut

  510                continue

                  else if( izlog .ne. 0 ) then

                     do 511 i = 1, icut

                        rcut(i,1) = log10( rcut(i,1) )

  511                continue

                  end if

                  if( icct .eq. 0 .and. iccr .eq. 1 ) then

                     icct = icut

                     do 520 i = 1, icut + 1

                        rcut(i,2) = -2.0

  520                continue

                     goto 560

                  end if

                  if( icct .eq. 0 .and. iccr .gt. 1 ) icct = -1

                  if( icct .gt. 0 .and. icct .lt. icut ) then

                     do 530 i = icct + 1, icut + 1

                        rcut(i,2) = rcut(icct,2)
                        rcut(i,3) = rcut(icct,3)
                        rcut(i,4) = rcut(icct,4)

  530                continue

                     icct = icut

                     goto 560

                  end if


                  if( icct .lt. 0 ) then

                     if( smax(1) * smin(1) .le. 0.0d0 ) then
                        smin(1) = 3.0
                        smax(1) = 1.0
                     end if

                     if( icut .gt. 0 ) then

                         do 550 ik = 1, icut

                            dhh = max(rcut(ik,1),cxymin)
                            dhh = min(dhh,cxymax)

                            rcut(ik,2) = smin(1) + ( dhh - cxymin )
     &                                 / xyabs * ( smax(1) - smin(1) )
                            rcut(ik,2) = chue(rcut(ik,2))
                            rcut(ik,3) = smin(2) + ( dhh - cxymin )
     &                                 / xyabs * ( smax(2) - smin(2) )
                            rcut(ik,4) = smin(3) + ( dhh - cxymin )
     &                                 / xyabs * ( smax(3) - smin(3) )


  550                    continue

                            ik = icut + 1

                            rcut(ik,2) = smax(1)
                            rcut(ik,2) = chue(rcut(ik,2))
                            rcut(ik,3) = smax(2)
                            rcut(ik,4) = smax(3)

                     end if

                  end if

  560             continue

*-----------------------------------------------------------------------
c
*-----------------------------------------------------------------------
*        ono's program
*-----------------------------------------------------------------------

                  i = 0

  500             if( i + 1 .gt. icut ) goto 555

                  i = i + 1

                  call cntdrw(ixnm,iynm,xdd,ydd,
     &                        dax(1),day(1),daxy(1,1),
     &                        rcut(i,1),jof,
     &                        iwd2(i),ipd2,noned,
     &                        rcut(i,2),rcut(i,3),rcut(i,4),
     &                        rcut(i,5),
     &                        itrace(1,1),ierr)

                     if( ierr .ne. 0 ) return

                  goto 500

  555             continue

*-----------------------------------------------------------------------
*              INITIALIZE AGAIN FOR CONTOUR PLOT
*-----------------------------------------------------------------------

                  icut  =  0
                  ncut  =  8
                  icct  =  0

*-----------------------------------------------------------------------


      return
      end

************************************************************************
*                                                                      *
      subroutine cntdrwf(iccr,ibmap,
     &                  icut,rcut,iwd2,izlog,
     &                  nx,ny,vx,vy,imat,daxy,
     &                  dmax,dmin,xdel,ydel)
*                                                                      *
*       sub program of output counter plot                             *
*       last modified by K.Niita on 2004/12/16                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

      parameter ( eps = 1.0d-6 )

*-----------------------------------------------------------------------

      dimension rcut(mc,5)
      dimension iwd2(mc)
      dimension gcut(mc,5)
      dimension nord(mc)

*-----------------------------------------------------------------------

      dimension   imat(nx,ny)
      dimension   daxy(nx,ny)
      dimension   vx(nx), vy(ny)
      dimension   rcol(3), bcol(3)

*-----------------------------------------------------------------------
*     initialization
*-----------------------------------------------------------------------

            sddel = min( xdel/dble(nx), ydel/dble(ny) )**2

            iclip = 1
            nfrmm = 0
            nbmap = 0
            npath = 0
            ihsb  = 1
            ilizt = iwd2(1)

*-----------------------------------------------------------------------
*     reordering of rcut to gcut
*-----------------------------------------------------------------------

            if( icut .gt. 1 ) then

               do i = 1, icut
                  nord(i) = i
               end do

               do i = 1, icut - 1
                     bnmin = rcut(nord(i),1)
                  do j = i + 1, icut
                     if( rcut(nord(j),1) .lt. bnmin ) then
                        bnmin = rcut(nord(j),1)
                        ii = nord(i)
                        nord(i) = nord(j)
                        nord(j) = ii
                     end if
                  end do
               end do

               do i = 1, icut
                  gcut(i+2,1) = rcut(nord(i),1)
                  gcut(i+2,2) = rcut(nord(i),2)
                  gcut(i+2,3) = rcut(nord(i),3)
                  gcut(i+2,4) = rcut(nord(i),4)
               end do

            else
                  gcut(3,1) = rcut(1,1)
                  gcut(3,2) = rcut(1,2)
                  gcut(3,3) = rcut(1,3)
                  gcut(3,4) = rcut(1,4)
            end if

*-----------------------------------------------------------------------

               if( izlog .eq. 1 ) then

                  fxymax =  41.0
                  dxymax =  41.0
                  dxymin = -41.0
                  fxymin = -41.0

                  if(dmax.gt.-r0max) dxymax = log10(dmax)
                  if(dmin.lt.+r0max) dxymin = log10(dmin)

               else

                  fxymax =  r1max
                  dxymax =  r1max
                  dxymin = -r1max
                  fxymin = -r1max

                  if(dmax.gt.-r0max) dxymax = dmax
                  if(dmin.lt.+r0max) dxymin = dmin

               end if

                  gcut(1,1) = fxymin
                  gcut(1,2) = -1.0
                  gcut(1,3) =  1.0
                  gcut(1,4) =  1.0

                  gcut(2,1) = dxymin
                  gcut(2,2) = -1.0
                  gcut(2,3) =  1.0
                  gcut(2,4) =  1.0

                  gcut(icut+3,1) = dxymax
                  gcut(icut+3,2) = rcut(icut+1,2)
                  gcut(icut+3,3) = rcut(icut+1,3)
                  gcut(icut+3,4) = rcut(icut+1,4)

                  gcut(icut+4,1) = fxymax
                  gcut(icut+4,2) = -1.0
                  gcut(icut+4,3) =  1.0
                  gcut(icut+4,4) =  1.0

                  ncut = icut + 4

*-----------------------------------------------------------------------

               do j = 1, ny
               do i = 1, nx

                  do k = 2, ncut

                     if( daxy(i,j) .le. gcut(k,1) ) goto 210

                  end do

                  k = ncut

  210             continue

                  imat(i,j) = k

               end do
               end do

*-----------------------------------------------------------------------
*     start check each cell
*-----------------------------------------------------------------------

  100 continue

               do j = 1, ny
               do i = 1, nx

                  if( imat(i,j) .lt. 0 ) imat(i,j) = 0

               end do
               end do

               do j = 1, ny
               do i = 1, nx

                  if( imat(i,j) .gt. 0 ) goto 200

               end do
               end do

               goto 900

*-----------------------------------------------------------------------
*        start new region
*-----------------------------------------------------------------------

  200    continue

               kk = imat(i,j)

               iline = 0
            if( iccr .eq. 1 .or. iccr .eq. 3 ) iline = 1

            if( iccr .eq. 1 ) then
               rcol(1) = gcut(kk,2)
               rcol(2) = gcut(kk,3)
               rcol(3) = gcut(kk,4)
               bcol(1) = -r1max
               bcol(2) =  1.0
               bcol(3) =  1.0
            else
               rcol(1) = -2.0
               rcol(2) = 1.0
               rcol(3) = 1.0
               bcol(1) = gcut(kk,2)
               bcol(2) = gcut(kk,3)
               bcol(3) = gcut(kk,4)
            end if

               vmax = gcut(kk,1)
               vmin = gcut(kk-1,1)

               rewind(jhq)
               ip = 0

               i0 = i
               j0 = j


*-----------------------------------------------------------------------
*           ic = 5 : first left check
*-----------------------------------------------------------------------

            ic = 5

                  ii = i - 1
                  jj = j

            if( ii .gt. 0 ) then

                     vval = vmax
                     if( daxy(ii,jj) .le. vmin ) vval = vmin

                  if( abs( ( daxy(ii,jj) - daxy(i,j) ) / vval )
     &                .gt. eps ) then

                     fac = ( daxy(ii,jj) - vval )
     &                   / ( daxy(ii,jj) - daxy(i,j) )

                  else

                     fac = 0.0

                  end if

                     x0 = vx(ii) + xdel * fac
                     y0 = vy(jj)

            else

                     x0 = vx(i)
                     y0 = vy(j)

            end if

                     x1 = x0
                     y1 = y0


*-----------------------------------------------------------------------
*           ic = 6 : first left-down check
*-----------------------------------------------------------------------

            ic = 6

                  ii = i - 1
                  jj = j - 1

                  icm = 0

            if( ii .gt. 0 .and. jj .gt. 0 ) then

                  if( daxy(ii,jj) .gt. vmax .or.
     &                daxy(ii,jj) .le. vmin ) then

                     vval = vmax
                     if( daxy(ii,jj) .le. vmin ) vval = vmin

                  if( abs( ( daxy(ii,jj) - daxy(i,j) ) / vval )
     &                .gt. eps ) then

                     fac = ( daxy(ii,jj) - vval )
     &                   / ( daxy(ii,jj) - daxy(i,j) )

                  else

                     fac = 0.0

                  end if

                     x1 = vx(ii) + xdel * fac
                     y1 = vy(jj) + ydel * fac

                     icm = 1

                  end if

            end if

*-----------------------------------------------------------------------
*           ic = 7 : first down check
*-----------------------------------------------------------------------

            ic = 7

                  ii = i
                  jj = j - 1

            if( jj .gt. 0 ) then

                     vval = vmax
                     if( daxy(ii,jj) .le. vmin ) vval = vmin

                  if( abs( ( daxy(ii,jj) - daxy(i,j) ) / vval )
     &                .gt. eps ) then

                     fac = ( daxy(ii,jj) - vval )
     &                   / ( daxy(ii,jj) - daxy(i,j) )

                  else

                     fac = 0.0

                  end if

                     vxx = vx(ii) + xdel * fac
                     vyy = vy(jj)

            else

                     vxx = vx(i)
                     vyy = vy(j)

            end if

                  if( icm .eq. 0 ) then

                     x1 = x0
                     y1 = y0

                  end if

                     x2 = vxx
                     y2 = vyy

                     x12 = ( x2 - x1 )
                     y12 = ( y2 - y1 )
                     r12 = sqrt( x12**2 + y12**2 )

*-----------------------------------------------------------------------
*        write new region
*-----------------------------------------------------------------------

               imat(i,j) = -imat(i,j)

               ic  = 8
               kc  = 3
               icf = 0

*-----------------------------------------------------------------------
*        eight directions movement and final write information
*-----------------------------------------------------------------------

  300    continue

*-----------------------------------------------------------------------
*     write information
*-----------------------------------------------------------------------

      if( i .eq. i0 .and. j .eq. j0 .and.
     &  ( ic .eq. 5 .or. ic .eq. 6 .or. ic .eq. 7 ) ) then

               if( icf .eq. 0 ) then

                  x3 = x0
                  y3 = y0

                  icf = 1

                  goto 310

               else

                     ip = ip + 1
                     write(jhq) x1, y1
                     ip = ip + 1
                     write(jhq) x2, y2

               end if

*-----------------------------------------------------------------------

         if( ip .gt. 3 ) then

            ibmap = ibmap + 1
            nclmm = ip

            write(jhd) iclip, nfrmm, nclmm, nbmap, npath,
     &                 ihsb, iline, ilizt, rcol, bcol

               rewind(jhq)

            do jj = 1, ip

                   read(jhq) x1, y1
                  write(jhb) x1, y1

            end do

         end if

*-----------------------------------------------------------------------

               goto 100

      end if

*-----------------------------------------------------------------------
*        four directions movement
*-----------------------------------------------------------------------

  500    continue

               if( ic .eq. 1 ) goto 501
               if( ic .eq. 2 ) goto 502
               if( ic .eq. 3 ) goto 503
               if( ic .eq. 4 ) goto 504
               if( ic .eq. 5 ) goto 505
               if( ic .eq. 6 ) goto 506
               if( ic .eq. 7 ) goto 507
               if( ic .eq. 8 ) goto 508

*-----------------------------------------------------------------------
*        ic = 1
*-----------------------------------------------------------------------

  501    ic = 1

                        ii = i + 1
                        jj = j

            if( ii .le. nx ) then

               if( kk .le. abs(imat(ii,jj)) ) then

                  if( j + 1 .le. ny ) then

                           ic = 0

                     do jl = j + 1, ny

                        if( kk .le. abs(imat(ii,jl)) ) then

                           if( imat(ii,jl) .eq. kk .and.
     &                         ic .eq. 0 ) then

                              imat(ii,jl) = -kk

                           else if( kk .lt. abs(imat(ii,jl)) ) then

                              ic = 1

                           end if

                        else

                           goto 601

                        end if

                     end do

  601                continue

                  end if

                        i = ii
                        j = jj

                        imat(i,j) = -kk

                        ic = 7
                        goto 300

               else

                        vval = vmax
                        if( daxy(ii,jj) .le. vmin ) vval = vmin

                  if( abs( ( daxy(ii,jj) - daxy(i,j) ) / vval )
     &                .gt. eps ) then

                     fac = ( daxy(ii,jj) - vval )
     &                   / ( daxy(ii,jj) - daxy(i,j) )

                  else

                     fac = 0.0

                  end if

                        x3 = vx(ii) - xdel * fac
                        y3 = vy(jj)

                        ic = 2
                        goto 310

               end if

            else

                        x3 = vx(i)
                        y3 = vy(j)

                        ic = 3
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 2
*-----------------------------------------------------------------------

  502    ic = 2

                        ii = i + 1
                        jj = j + 1

            if( ii .le. nx .and. jj .le. ny ) then

               if( kk .gt. abs(imat(ii,jj)) .and.
     &             imat(ii,jj) .ne. 0 ) then

                        vval = vmax
                        if( daxy(ii,jj) .le. vmin ) vval = vmin

                  if( abs( ( daxy(ii,jj) - daxy(i,j) ) / vval )
     &                .gt. eps ) then

                     fac = ( daxy(ii,jj) - vval )
     &                   / ( daxy(ii,jj) - daxy(i,j) )

                  else

                     fac = 0.0

                  end if

                        x3 = vx(ii) - xdel * fac
                        y3 = vy(jj) - ydel * fac

                        ic = 3
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 3
*-----------------------------------------------------------------------

  503    ic = 3

                        ii = i
                        jj = j + 1

            if( jj .le. ny ) then

               if( kk .le. abs(imat(ii,jj)) ) then

                  if( i - 1 .ge. 1 ) then

                           ic = 0

                     do il = i - 1, 1, -1

                        if( kk .le. abs(imat(il,jj)) ) then

                           if( imat(il,jj) .eq. kk .and.
     &                         ic .eq. 0 ) then

                              imat(il,jj) = -kk

                           else if( kk .lt. abs(imat(il,jj)) ) then

                              ic = 1

                           end if

                        else

                           goto 603

                        end if

                     end do

  603                continue

                  end if

                        i = ii
                        j = jj

                        imat(i,j) = -kk

                        ic = 1
                        goto 300

               else

                        vval = vmax
                        if( daxy(ii,jj) .le. vmin ) vval = vmin

                  if( abs( ( daxy(ii,jj) - daxy(i,j) ) / vval )
     &                .gt. eps ) then

                     fac = ( daxy(ii,jj) - vval )
     &                   / ( daxy(ii,jj) - daxy(i,j) )

                  else

                     fac = 0.0

                  end if

                        x3 = vx(ii)
                        y3 = vy(jj) - ydel * fac

                        ic = 4
                        goto 310

               end if

            else

                        x3 = vx(i)
                        y3 = vy(j)

                        ic = 5
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 4
*-----------------------------------------------------------------------

  504    ic = 4

                        ii = i - 1
                        jj = j + 1

            if( ii .gt. 0 .and. jj .le. ny ) then

               if( kk .gt. abs(imat(ii,jj)) .and.
     &             imat(ii,jj) .ne. 0 ) then

                        vval = vmax
                        if( daxy(ii,jj) .le. vmin ) vval = vmin

                  if( abs( ( daxy(ii,jj) - daxy(i,j) ) / vval )
     &                .gt. eps ) then

                     fac = ( daxy(ii,jj) - vval )
     &                   / ( daxy(ii,jj) - daxy(i,j) )

                  else

                     fac = 0.0

                  end if

                        x3 = vx(ii) + xdel * fac
                        y3 = vy(jj) - ydel * fac

                        ic = 5
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 5
*-----------------------------------------------------------------------

  505    ic = 5

                        ii = i - 1
                        jj = j

            if( ii .gt. 0 ) then

               if( kk .le. abs(imat(ii,jj)) ) then

                  if( j - 1 .ge. 1 ) then

                           ic = 0

                     do jl = j - 1, 1, -1

                        if( kk .le. abs(imat(ii,jl)) ) then

                           if( imat(ii,jl) .eq. kk .and.
     &                         ic .eq. 0 ) then

                              imat(ii,jl) = -kk

                           else if( kk .lt. abs(imat(ii,jl)) ) then

                              ic = 1

                           end if

                        else

                           goto 605

                        end if

                     end do

  605                continue

                  end if

                        i = ii
                        j = jj

                        imat(i,j) = -kk

                        ic = 3
                        goto 300

               else

                        vval = vmax
                        if( daxy(ii,jj) .le. vmin ) vval = vmin

                  if( abs( ( daxy(ii,jj) - daxy(i,j) ) / vval )
     &                .gt. eps ) then

                     fac = ( daxy(ii,jj) - vval )
     &                   / ( daxy(ii,jj) - daxy(i,j) )

                  else

                     fac = 0.0

                  end if

                        x3 = vx(ii) + xdel * fac
                        y3 = vy(jj)

                        ic = 6
                        goto 310

               end if

            else

                        x3 = vx(i)
                        y3 = vy(j)

                        ic = 7
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 6
*-----------------------------------------------------------------------

  506    ic = 6

                        ii = i - 1
                        jj = j - 1

            if( ii .gt. 0 .and. jj .gt. 0 ) then

               if( kk .gt. abs(imat(ii,jj)) .and.
     &             imat(ii,jj) .ne. 0 ) then

                        vval = vmax
                        if( daxy(ii,jj) .le. vmin ) vval = vmin

                  if( abs( ( daxy(ii,jj) - daxy(i,j) ) / vval )
     &                .gt. eps ) then

                     fac = ( daxy(ii,jj) - vval )
     &                   / ( daxy(ii,jj) - daxy(i,j) )

                  else

                     fac = 0.0

                  end if

                        x3 = vx(ii) + xdel * fac
                        y3 = vy(jj) + ydel * fac

                        ic = 7
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 7
*-----------------------------------------------------------------------

  507    ic = 7

                        ii = i
                        jj = j - 1

            if( jj .gt. 0 ) then

               if( kk .le. abs(imat(ii,jj)) ) then

                  if( i + 1 .le. nx ) then

                           ic = 0

                     do il = i + 1, nx

                        if( kk .le. abs(imat(il,jj)) ) then

                           if( imat(il,jj) .eq. kk .and.
     &                         ic .eq. 0 ) then

                              imat(il,jj) = -kk

                           else if( kk .lt. abs(imat(il,jj)) ) then

                              ic = 1

                           end if

                        else

                           goto 607

                        end if

                     end do

  607                continue

                  end if

                        i = ii
                        j = jj

                        imat(i,j) = -kk

                        ic = 5
                        goto 300

               else

                        vval = vmax
                        if( daxy(ii,jj) .le. vmin ) vval = vmin

                  if( abs( ( daxy(ii,jj) - daxy(i,j) ) / vval )
     &                .gt. eps ) then

                     fac = ( daxy(ii,jj) - vval )
     &                   / ( daxy(ii,jj) - daxy(i,j) )

                  else

                     fac = 0.0

                  end if

                        x3 = vx(ii)
                        y3 = vy(jj) + ydel * fac

                        ic = 8
                        goto 310

               end if

            else

                        x3 = vx(i)
                        y3 = vy(j)

                        ic = 1
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 8
*-----------------------------------------------------------------------

  508    ic = 8

                     ii = i + 1
                     jj = j - 1

            if( ii .le. nx .and. jj .gt. 0 ) then

               if( kk .gt. abs(imat(ii,jj)) .and.
     &             imat(ii,jj) .ne. 0 ) then

                        vval = vmax
                        if( daxy(ii,jj) .le. vmin ) vval = vmin

                  if( abs( ( daxy(ii,jj) - daxy(i,j) ) / vval )
     &                .gt. eps ) then

                     fac = ( daxy(ii,jj) - vval )
     &                   / ( daxy(ii,jj) - daxy(i,j) )

                  else

                     fac = 0.0

                  end if

                        x3 = vx(ii) - xdel * fac
                        y3 = vy(jj) + ydel * fac

                        ic = 1
                        goto 310

               end if

            end if

                  goto 501

*-----------------------------------------------------------------------

  310       continue

                           x13 = ( x3 - x1 )
                           y13 = ( y3 - y1 )
                           r13 = sqrt( x13**2 + y13**2 )

                           if( r13 .lt. sddel ) goto 300

                           x23 = ( x3 - x2 )
                           y23 = ( y3 - y2 )
                           r23 = sqrt( x23**2 + y23**2 )

                           if( r23 .lt. sddel ) goto 300

                        rdcos = ( x13 * x23 + y13 * y23 ) / r13 / r23

                     if( rdcos .gt. 0.9999d0 ) then

                           x2 = x3
                           y2 = y3

                     else

                           ip = ip + 1
                           write(jhq) x1, y1

                           x1 = x2
                           y1 = y2
                           x2 = x3
                           y2 = y3

                     end if

            goto 300

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine cntdrw(nx,ny,dx,dy,dax,day,p,z,jof,
     &                  iwd2,ipd2,noned,
     &                  col1,col2,col3,rlin,itrace,ierr)
*                                                                      *
*         NX,NY   : DIMENSION OF MATRIX P( , )                         *
*         DX,DY   : MESH WIDTH                                         *
*         DAX,DAY : GRIT POINT COORDINATE                              *
*         P( , )  : DATA ON GRID                                       *
*         Z       : POTENTIAL                                          *
*                                                                      *
*         JOL     : OUTPUT FILE                                        *
*         IWD2    : WIDTH OF LINES                                     *
*         IPD2    : 1-> SPLINE 0-> NO                                  *
*                                                                      *
*         NONED   : NUMBER OF LINES FOR JOL                            *
*         COL     : COLOR OF LINE                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter(istx=0,isty=1,ists=2)

      dimension p(nx,ny)
      dimension dax(nx),day(ny)
      dimension itrace(nx,ny)
      dimension col(3)

      logical ods

      ods(i,j,ii,jj)=(p(i,j).lt.z .and. p(ii,jj).ge.z)
     &     .or. (p(i,j).ge.z .and. p(ii,jj).lt.z)

      ierr = 0
      col(1) = col1
      col(2) = col2
      col(3) = col3

      do 10 i=1,nx
      do 10 j=1,ny
          itrace(i,j)=1
          if(p(i,j).eq.z) p(i,j)=p(i,j)+abs(p(i,j)*1.e-7)
   10 continue

      do 120 ii=1,nx-1
        i=ii
        j=1
        if(itrace(i,j).eq.0) goto 120
        if(ods(i,j,i+1,j)) then
          istate=istx
        else
          goto 120
        endif

        call cntsub(p(1,1),nx,ny,dx,dy,dax,day,z,istate,30,i,j,
     &                  itrace,jof,iwd2,ipd2,noned,
     &                  rlin,col,ierr)
            if( ierr .ne. 0 ) return

  120 continue

      do 130 jj=1,ny-1
        i=1
        j=jj
        if(itrace(i,j).eq.0) goto 130
        if(ods(i,j,i,j+1)) then
          istate=isty
        else
          goto 130
        endif

        call cntsub(p(1,1),nx,ny,dx,dy,dax,day,z,istate,30,i,j,
     &                  itrace,jof,iwd2,ipd2,noned,
     &                  rlin,col,ierr)
            if( ierr .ne. 0 ) return

  130 continue

      do 140 ii=1,nx-1
        i=ii
        j=ny
        if(itrace(i,j).eq.0) goto 140
        if(ods(i,j,i+1,j)) then
          istate=istx
        else
          goto 140
        endif

        call cntsub(p(1,1),nx,ny,dx,dy,dax,day,z,istate,50,i,j,
     &                  itrace,jof,iwd2,ipd2,noned,
     &                  rlin,col,ierr)
            if( ierr .ne. 0 ) return

  140 continue

      do 150 jj=1,ny-1
        i=nx
        j=jj
        if(itrace(i,j).eq.0) goto 150
        if(ods(i,j,i,j+1)) then
          istate=isty
        else
          goto 150
        endif

        call cntsub(p(1,1),nx,ny,dx,dy,dax,day,z,istate,50,i,j,
     &                  itrace,jof,iwd2,ipd2,noned,
     &                  rlin,col,ierr)
            if( ierr .ne. 0 ) return

  150 continue

      do 160 ii=1,nx-1
      do 160 jj=1,ny-1
        i=ii
        j=jj
        if(itrace(i,j).eq.0) goto 160
        if(ods(i,j,i+1,j)) then
          istate=istx
        else if(ods(i,j,i,j+1)) then
          istate=isty
        else
          goto 160
        endif

        call cntsub(p(1,1),nx,ny,dx,dy,dax,day,z,istate,30,i,j,
     &                  itrace,jof,iwd2,ipd2,noned,
     &                  rlin,col,ierr)
            if( ierr .ne. 0 ) return

  160 continue

      return
      end
************************************************************************
*                                                                      *
      subroutine cntsub(p,nx,ny,dx,dy,dax,day,z,istate,ientry,i,j,
     &                  itrace,jof,iwd2,ipd2,noned,
     &                  rlin,col,ierr)
*                                                                      *
*         P( , )  : DATA ON GRID                                       *
*         NX,NY   : DIMENSION OF MATRIX P( , )                         *
*         DX,DY   : MESH WIDTH                                         *
*         DAX,DAY : GRIT POINT COORDINATE                              *
*         Z       : POTENTIAL                                          *
*         ISTATE  : = ISTX, ISTY, ISTS                                 *
*         IENTRY  : 30 (BOARDER) 50 (NORMAL)                           *
*                                                                      *
*         JOL     : OUTPUT FILE                                        *
*         IWD2    : WIDTH OF LINESZ                                    *
*         IPD2    : 1-> SPLINE 0-> NO                                  *
*                                                                      *
*         NONED   : NUMBER OF LINES FOR JOL                            *
*         rlin    : kind OF LINE                                       *
*         COL     : COLOR OF LINE                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter(istx=0,isty=1,ists=2)
      parameter ( eps = 1.0d-6 )

      dimension p(nx,ny)
      dimension dax(nx),day(ny)
      dimension itrace(nx,ny), ntrace(nx,ny)
      dimension col(3)
      real(8),allocatable:: dspl(:,:)
      logical ods

      ods(i,j,ii,jj)=(p(i,j).lt.z .and. p(ii,jj).ge.z)
     &            .or. (p(i,j).ge.z .and. p(ii,jj).lt.z)


      ic = 0
      ierr = 0


*-----------------------------------------------------------------------

         call pre_cntsub(p,nx,ny,z,istate,ientry,i,j,itrace,il)
               ipol  = ipd2
               if( ipol .lt. 0 .and. il .lt. 5 ) ipol = 0
               if( ipol .gt. 0 .and. il .lt. 3 ) ipol = 0
               ilas = il + ( il - 1 ) * abs(ipol)

         allocate(dspl(15,ilas))

*-----------------------------------------------------------------------

      if(ientry.eq.50) goto 50

 30   continue

        if(istate.eq.istx) then

            if( abs((p(i+1,j)-p(i,j))/z) .gt. eps ) then

              divxv=dax(i)+(z-p(i,j))/(p(i+1,j)-p(i,j))*dx

            else

              divxv=dax(i)

            end if

              call cntout(ic,divxv,day(j),jof,
     &                    iwd2,ipd2,noned,rlin,col,dspl,ierr)
              if( ierr .ne. 0 ) goto 99
              ic = ic + 1
        end if

        if(istate.eq.isty) then

            if( abs((p(i,j+1)-p(i,j))/z) .gt. eps ) then

              divyv=day(j)+(z-p(i,j))/(p(i,j+1)-p(i,j))*dy

            else

              divyv=day(j)

            end if

              call cntout(ic,dax(i),divyv,jof,
     &                    iwd2,ipd2,noned,rlin,col,dspl,ierr)
              if( ierr .ne. 0 ) goto 99
              ic = ic + 1
        end if

        if(itrace(i,j).eq.0) goto 40

        itrace(i,j)=0

        if(i.eq.nx .or. j.eq.ny) goto 40

        if(istate.eq.istx) then
          if(ods(i+1,j,i,j+1)) then
            istate=ists
          else
            istate=isty
          endif
        else if(istate.eq.isty) then

          if(ods(i+1,j,i,j+1)) then
            istate=ists
          else
            istate=istx
          endif

        else
          if(ods(i,j,i+1,j)) then
            istate=istx
          else
            istate=isty
          endif
        endif

   50 continue

        if(istate.eq.istx) then

            if( abs((p(i+1,j)-p(i,j))/z) .gt. eps ) then

              divxv=dax(i)+(z-p(i,j))/(p(i+1,j)-p(i,j))*dx

            else

              divxv=dax(i)

            end if

              call cntout(ic,divxv,day(j),jof,
     &                    iwd2,ipd2,noned,rlin,col,dspl,ierr)
              if( ierr .ne. 0 ) goto 99
              ic = ic + 1

        else if(istate.eq.isty) then

            if( abs((p(i,j+1)-p(i,j))/z) .gt. eps ) then

              divyv=day(j)+(z-p(i,j))/(p(i,j+1)-p(i,j))*dy

            else

              divyv=day(j)

            end if

              call cntout(ic,dax(i),divyv,jof,
     &                    iwd2,ipd2,noned,rlin,col,dspl,ierr)
              if( ierr .ne. 0 ) goto 99
              ic = ic + 1

        endif

        if(i.eq.1 .and. istate.eq.isty) goto 40
        if(j.eq.1 .and. istate.eq.istx) goto 40

        if(istate.eq.istx) then

          if(ods(i,j,i+1,j-1)) then
            istate=ists
            j=j-1
          else
            istate=isty
            i=i+1
            j=j-1
          endif

        else if(istate.eq.isty) then

          if(ods(i,j,i-1,j+1)) then
            istate=ists
            i=i-1
          else
            istate=istx
            i=i-1
            j=j+1
          endif

        else

          if(ods(i+1,j,i+1,j+1)) then
            istate=isty
            i=i+1
          else
            istate=istx
            j=j+1
          endif

        endif

      goto 30

   40 continue
              call cntout(-1,0.0d0,0.0d0,jof,
     &                    iwd2,ipd2,noned,rlin,col,dspl,ierr)
              if( ierr .ne. 0 ) goto 99

   99 deallocate(dspl)
      return
      end

************************************************************************
*                                                                      *
      subroutine pre_cntsub(p,nx,ny,z,istate,ientry,i,j,itrace,il)
*                                                                      *
*         P( , )  : DATA ON GRID                                       *
*         NX,NY   : DIMENSION OF MATRIX P( , )                         *
*         Z       : POTENTIAL                                          *
*         ISTATE  : = ISTX, ISTY, ISTS                                 *
*         IENTRY  : 30 (BOARDER) 50 (NORMAL)                           *
*         IL      : REQUEST SIZE OF 2ND RANK IN DSPL(15,*)             *
*                                                                      *
************************************************************************

      implicit none

*-----------------------------------------------------------------------

      integer,parameter:: istx=0,isty=1,ists=2
      integer,parameter:: eps = 1.0d-6

      integer,intent(in):: i,j,nx,ny
      integer,intent(in):: istate,ientry
      real(8),intent(in):: z
      real(8),intent(in):: p(nx,ny)
      integer,intent(in):: itrace(nx,ny)
      integer,intent(out):: il

      integer ii,jj,ni,nj,ic
      integer nstate,ntrace(nx,ny)

      logical ods

      ods(i,j,ii,jj)=(p(i,j).lt.z .and. p(ii,jj).ge.z)
     &            .or. (p(i,j).ge.z .and. p(ii,jj).lt.z)

      ntrace(:,:) = itrace(:,:)
      nstate = istate
      ni = i
      nj = j

      ic = 0

      if(ientry.eq.50) goto 50

 30   continue

        if(nstate.eq.istx .or. nstate.eq.isty) then

            if (ic.eq.0) il = 1
            if (ic.gt.0) il = il + 1
            ic = ic + 1

        end if

        if(ntrace(ni,nj).eq.0) goto 40

        ntrace(ni,nj)=0

        if(ni.eq.nx .or. nj.eq.ny) goto 40

        if(nstate.eq.istx) then
          if(ods(ni+1,nj,ni,nj+1)) then
            nstate=ists
          else
            nstate=isty
          endif
        else if(nstate.eq.isty) then
          if(ods(ni+1,nj,ni,nj+1)) then
            nstate=ists
          else
            nstate=istx
          endif
        else
          if(ods(ni,nj,ni+1,nj)) then
            nstate=istx
          else
            nstate=isty
          endif
        endif

   50 continue

        if(nstate.eq.istx .or. nstate.eq.isty) then

            if (ic.eq.0) il = 1
            if (ic.gt.0) il = il + 1
            ic = ic + 1

        endif

        if(ni.eq.1 .and. nstate.eq.isty) goto 40
        if(nj.eq.1 .and. nstate.eq.istx) goto 40

        if(nstate.eq.istx) then

          if(ods(ni,nj,ni+1,nj-1)) then
            nstate=ists
            nj=nj-1
          else
            nstate=isty
            ni=ni+1
            nj=nj-1
          endif

        else if(nstate.eq.isty) then

          if(ods(ni,nj,ni-1,nj+1)) then
            nstate=ists
            ni=ni-1
          else
            nstate=istx
            ni=ni-1
            nj=nj+1
          endif

        else

          if(ods(ni+1,nj,ni+1,nj+1)) then
            nstate=isty
            ni=ni+1
          else
            nstate=istx
            nj=nj+1
          endif

        endif

      goto 30

   40 continue

      !! do not increment il, becuase ic == -1

      return
      end

************************************************************************
*                                                                      *
      subroutine cntout(ic,x,y,jof,iwd2,ipd2,noned,rlin,rcol,dspl,ierr)
*                                                                      *
*         IC    :   0  -> START OF THE LINE                            *
*                  >0  -> SEQUENTIAL LINES                             *
*                  -1  -> END OF LINE                                  *
*                                                                      *
*         X, Y  :  COORDINATE OF THE POINTS                            *
*                                                                      *
*         JOL   :  OUTPUT FILE                                         *
*                                                                      *
*         IWD2  :  WIDTH OF LINES                                      *
*         IPD2  :  1-> SPLINE 0-> NO                                   *
*                                                                      *
*         NONED   : NUMBER OF LINES FOR JOL                            *
*         Rlin    : kind OF LINE                                       *
*         RCOL    : COLOR OF LINE                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

      dimension dspl(15,*)
      dimension rcol(3), rcob(3)

      save il

*-----------------------------------------------------------------------

      ierr = 0

      z   = 0.0
      iz  = 0
      w   = 1.0
      idc = 2

      rcob(1) = -1.0
      rcob(2) =  1.0
      rcob(3) =  1.0

*-----------------------------------------------------------------------

      if(ic.eq.0) then

          il=1
          dspl(1,il)=x
          dspl(2,il)=y

      else if(ic.gt.0) then

          il=il+1
          dspl(1,il)=x
          dspl(2,il)=y

      else

          ityp  = nint( rlin )
          iwi   = iwd2
          ipol  = ipd2
          imar  = 0
          icoma = 0
          facsiz = 0.0

*-----------------------------------------------------------------------
*           MINIMUM POINTS FOR SPLINE
*-----------------------------------------------------------------------

               if( ipol .lt. 0 .and. il .lt. 5 ) ipol = 0
               if( ipol .gt. 0 .and. il .lt. 3 ) ipol = 0

*-----------------------------------------------------------------------
*        SPLINE IPD2
*-----------------------------------------------------------------------
            if( ipol .ne. 0 ) then
*-----------------------------------------------------------------------

               syss = dble(ipol)

               call spln00(noned,noner,syss,il,ierr,
     &                     idc,iz,iz,ityp,iwi,0,
     &                     imar,rcol,rcob,icoma,facsiz,jof,jof,ityp,
     &                     dspl)

               if( ierr .ne. 0 ) goto 600

               goto 700

*-----------------------------------------------------------------------
            end if
*-----------------------------------------------------------------------

  600    continue

               noned = noned + 1

               ipol = 0

               write(jof) il,idc,iz,iz
               write(jof) ityp,iwi,ipol,imar,
     &                    rcol,rcob,icoma,facsiz

               do 100 i = 1, il

                  write(jof) dspl(1,i),dspl(2,i),z,z,z,z

  100          continue


  700    continue

      end if


      return
      end

************************************************************************
*                                                                      *
      subroutine writehd(jhf,idbg,ifon)
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character     today1*11, timestr1*5, today*100, timestr*100
      character     infnm(ichrl)*1, infn(ichrl)*1, avers*5
      common /dfil1/today, timestr, today1, timestr1, infnm, infn, avers
      common /dfil2/jtoday, jtime, jnm, mpage, iname, javers

*-----------------------------------------------------------------------

      common /con/  cm, dd

*-----------------------------------------------------------------------
*     HEADER OF EPS FILE
*-----------------------------------------------------------------------

      write(jhf,'(
     &''%!PS-Adobe-2.0''/
     &''%%Title: Input File = '',200a1
     &)')
     &  (infn(i),i=1,iname)

      write(jhf,'(
     &''%%CreationDate: '',A5,'' '',A11,/
     &''%%Creator: ANGEL ver.'',A5,/
     &''%%For: ANGEL users'',/
     &''%%DocumentFonts: (atend)'',/
     &''%%Pages: (atend)'',/
     &''%%BoundingBox: (atend)'',/
     &''%%EndComments''
     &)')
     &  timestr1, today1, avers


*-----------------------------------------------------------------------
*cKN 2023/08/23
*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine writehd1(jhf,idbg,ifon)
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character     today1*11, timestr1*5, today*100, timestr*100
      character     infnm(ichrl)*1, infn(ichrl)*1, avers*5
      common /dfil1/today, timestr, today1, timestr1, infnm, infn, avers
      common /dfil2/jtoday, jtime, jnm, mpage, iname, javers

*-----------------------------------------------------------------------

      common /con/  cm, dd

*cKN 2023/08/23
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*     WRITE BASIC OPERATERS AND CONSTANTS WITH COMMENTS
*-----------------------------------------------------------------------

      if( idbg .eq. 1 ) then

*-----------------------------------------------------------------------

      write(jhf,'(
     &''%%BeginProlog''
     &)')

      write(jhf,'(/''%'',71(''-''),/
     &             ''% Definition of Basic Operators''
     &           ,/''%'',71(''-'')/)')

      write(jhf,'(
     &''/N  {def} def''/
     &''/S  {exch} N''/
     &''/B  {bind N} N''/
     &)')

      write(jhf,'(
     &''/dv {transform round S round S itransform'',/
     &''     0.0967184 add S 0.0967184 add S} N''/
     &)')

      write(jhf,'(
     &''/mv {moveto} N''/
     &)')

      write(jhf,'(
     &''/L  {dv lineto} N'',/
     &''/M  {dv moveto} N'',/
     &''/TR {translate} N''/
     &)')

      write(jhf,'(''/x  {xal mul} N'')')
      write(jhf,'(''/y  {yal mul} N'')')
      write(jhf,'(''/xy {y S x S} N''/)')

      write(jhf,'(
     &''/l  {xy L} N'',/
     &''/m  {xy M} N''/
     &)')

      write(jhf,'(
     &''/np {newpath} N'',/
     &''/cp {closepath} N'',/
     &''/st {stroke} N'',/
     &''/fl {fill} N'',/
     &''/cl {clip} N'',/
     &''/gs {gsave} N'',/
     &''/gr {grestore} N''/
     &)')

      write(jhf,'(
     &''/sc {/ci3 S N /ci2 S N /ci1 S N ci1 0 lt {ci1 2 add setgray}''/
     &''    {ci1 1 sub ci2 ci3 sethsbcolor} ifelse} N''
     &/)')

      write(jhf,'(
     &''/sr {/ci3 S N /ci2 S N /ci1 S N ci1 0 lt {ci1 2 add setgray}''/
     &''    {ci1 1 sub ci2 ci3 setrgbcolor} ifelse} N''
     &/)')


*-----------------------------------------------------------------------
*     rotin0 : for axis, rotin: others
*-----------------------------------------------------------------------
*     rotin0  :(rx0 ry0 ang0) ; rx0 ry0 translate and ang1 rotate
*     rotout0 :               ; -ang0 rotate and -rx0 -ry0 translate

*     rotin   :(rx1 ry1 ang1) ; rx1 ry1 translate and ang1 rotate
*     rotout  :               ; -ang1 rotate and -rx1 -ry1 translate
*-----------------------------------------------------------------------

      write(jhf,'(
     &''/rotin0  {/ang0 S N /ry0 S N /rx0 S N '',
     &''rx0 ry0 TR ang0 rotate} N'',/
     &''/rotout0 { ang0 neg rotate} N''/
     &)')

      write(jhf,'(
     &''/rotin   {/ang1 S N /ry1 S N /rx1 S N '',
     &''rx1 ry1 TR ang1 rotate} N'',/
     &''/rotout  { ang1 neg rotate rx1 neg ry1 neg TR} N''/
     &)')

*-----------------------------------------------------------------------
*     scalin :(sx1 sy1 sclx scly) ; sx1 sy1 translate and
*                                   sclx scly scale and 0 0 moveto
*     scalout                     ; 0 0 transrate and
*                                   1/sclx 1/scly scale and
*                                   -sx1 -sy1 translate
*-----------------------------------------------------------------------

      write(jhf,'(
     &''/scalin  {/scly S N /sclx S N /sy1 S N '',
     &''/sx1 S N'',/
     &''           sx1 sy1 TR sclx scly scale 0 0 M} N'',/
     &''/scalout {1 sclx div 1 scly div scale sx1 neg '',
     &''sy1 neg TR} N''
     &)')

*-----------------------------------------------------------------------
*     LINE OPERATORS
*-----------------------------------------------------------------------

      write(jhf,'(/''%'',71(''-''),/
     &             ''% Line Operators''
     &           ,/''%'',71(''-'')/)')

      write(jhf,'(
     &''/lw  {setlinewidth} N'',/
     &''/sd  {setdash} N''/
     &)')

      write(jhf,'(
     &''/lc0 {0 setlinecap} N'',/
     &''/lc1 {1 setlinecap} N lc1'',/
     &''/lc2 {2 setlinecap} N''/
     &)')

      write(jhf,'(
     &''/lj0 {0 setlinejoin} N lj0'',/
     &''/lj1 {1 setlinejoin} N'',/
     &''/lj2 {2 setlinejoin} N''/
     &)')

      write(jhf,'(
     &''/sd0 {[] 0 sd} N /lm 0.001 N''
     &)')

*-----------------------------------------------------------------------
*     FONT OPERATORS
*-----------------------------------------------------------------------
*    /stsm :(str xp yp) ; Draw str1 at (xp,yp)
*    /stsw :(str xp yp ftn fts) ; Draw str1 at (xp,yp)
*                                       by scaling the chalacter
*    /stss :(str xp yp sx sy ftn fts) ; Draw str1 at (xp,yp)
*                                       by scaling the chalacter
*-----------------------------------------------------------------------

      write(jhf,'(/''%'',71(''-''),/
     &             ''% Font Operators''
     &           ,/''%'',71(''-'')/)')

      write(jhf,'(
     &''/stsm {mv show} N'',/
     &''/stsw {scalefont setfont mv show} N'',/
     &''/stss {scalefont setfont scalin show scalout} N''/
     &)')

*-----------------------------------------------------------------------
*     WRITE CONSTANTS
*-----------------------------------------------------------------------

*        CM : [ 28.346457 ] CM  to Point
*        DD : [ 0.24 ]      DPI to Point

*-----------------------------------------------------------------------

         write(jhf,'(/''%'',71(''-''),/
     &                ''% Constants''
     &              ,/''%'',71(''-'')/)')

         write(jhf,'(
     &   ''/cm {'',F10.6,'' mul} N'',/
     &   ''/dd {'',F5.2,'' mul} N''
     &   )') cm, dd

*-----------------------------------------------------------------------
*     CLIP PATH
*-----------------------------------------------------------------------

*        ax : axis path
*        ba : inside the screen

*-----------------------------------------------------------------------

         write(jhf,'(/''%'',71(''-''),/
     &                ''% Cliping Path Definition''
     &              ,/''%'',71(''-'')/)')

         write(jhf,'(''/ax { 0 0 M 0 yal L xal yal L '',
     &               ''xal 0 L cp} N'')')

         write(jhf,'(''/ba { xneg yneg M xneg ypst L xpst '',
     &               ''ypst L xpst yneg L cp} N'')')

         write(jhf,'(''/xpst 1000 cm N /ypst xpst N '',
     &               ''/xneg xpst neg N /yneg ypst neg N''/)')

*-----------------------------------------------------------------------

      write(jhf,'(
     &''%%EndProlog''
     &)')

*-----------------------------------------------------------------------
*     WITHOUT COMMENTS
*-----------------------------------------------------------------------

      else if( idbg .eq. 0 ) then

*-----------------------------------------------------------------------

      write(jhf,'(
     &''%%BeginProlog''
     &)')

*-----------------------------------------------------------------------

      write(jhf,'(
     &''/N {def} def /S {exch} N /B {bind N} N /dv {transform round '',
     &''S round S ''/
     &''itransform 0.0967184 add S 0.0967184 add S} N /mv {moveto} N'',
     &''/L {dv lineto} ''/
     &''N /M {dv moveto} N /TR {translate} N /x {xal mul} '',
     &''N /y {yal mul} N /xy {y S ''/
     &''x S} N /l {xy L} N /m {xy M} N /np {newpath} N /cp '',
     &''{closepath} N /st {stroke} ''
     &)')
      write(jhf,'(
     &''N /fl {fill} N /cl {clip} N /gs {gsave} N /gr {grestore} N '',
     &''/sc {/ci3 S N /ci2 ''/
     &''S N /ci1 S N ci1 0 lt {ci1 2 add setgray}{ci1 1 sub ci2 ci3 '',
     &''sethsbcolor} ''/
     &''ifelse} N /sr {/ci3 S N /ci2 S N /ci1 S N ci1 0 lt {ci1 2 '',
     &''add setgray}{ci1 1 ''/
     &''sub ci2 ci3 setrgbcolor} ifelse} N /rotin0 {/ang0 S N /ry0 '',
     &''S N /rx0 S N rx0 ''
     &)')
      write(jhf,'(
     &''ry0 TR ang0 rotate} N /rotout0 { ang0 neg rotate} N /rotin '',
     &''{/ang1 S N /ry1 S N ''/
     &''/rx1 S N rx1 ry1 TR ang1 rotate} N /rotout { ang1 neg '',
     &''rotate rx1 neg ry1 neg ''/
     &''TR} N /scalin {/scly S N /sclx S N /sy1 S N /sx1 S N sx1 '',
     &''sy1 TR sclx scly ''
     &)')
      write(jhf,'(
     &''scale 0 0 M} N /scalout {1 sclx div 1 scly div scale sx1 '',
     &''neg sy1 neg TR} N /lw ''/
     &''{setlinewidth} N /sd {setdash} N /lc0 {0 setlinecap} N /lc1 '',
     &''{1 setlinecap} N ''/
     &''lc1 /lc2 {2 setlinecap} N /lj0 {0 setlinejoin} N lj0 /lj1 '',
     &''{1 setlinejoin} N ''
     &)')
      write(jhf,'(
     &''/lj2 {2 setlinejoin} N /sd0 {[] 0 sd} N /lm 0.001 N /stsm '',
     &''{mv show} N /stsw ''/
     &''{scalefont setfont mv show} N /stss {scalefont setfont '',
     *''scalin show scalout} N''
     &)')
      write(jhf,'(
     &''/cm {28.346460 mul} N /dd {.24 mul} N /ax {0 0 M 0 yal L xal'',
     &'' yal L xal 0''/
     &''L cp} N /ba {xneg yneg M xneg ypst L xpst ypst L xpst yneg L'',
     &'' cp} N /xpst''/
     &''1000 cm N /ypst xpst N /xneg xpst neg N /yneg ypst neg N''
     &)')

*-----------------------------------------------------------------------

      write(jhf,'(
     &''%%EndProlog''
     &)')

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine inpexp(jsi,jsn,dsin,idsi,ierr)
*                                                                      *
*           expand input file on jhs                                   *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character dumr*10000
      character infn(ichrl)*1
      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------
*     INITIAL VALUE FOR INCLUDE FILE
*-----------------------------------------------------------------------

            ierr = 0

         do 30 i = 0, 9

            ill(i) = 0
            ilf(i) = 10000000

   30    continue

*-----------------------------------------------------------------------

  140       ill(jsn) = ill(jsn) + 1

               do ic = 1, icolm
                  dum(ic) = ' '
               end do

               read(jsi,'(10000a1)',iostat = ios ) (dum(ic),ic=1,icolm)
               if( ios .eq. -1 ) goto 251

*-----------------------------------------------------------------------

  240       continue

               if( ill(jsn) .gt. ilf(jsn) ) goto 251

               goto 250

  251          continue

                     if( jsn .le. 1 ) then

                        goto 141

                     else if( jsn .gt. 1 ) then

                        call closef(jsi,jsn)

                        goto 140

                     end if

  250       continue

*-----------------------------------------------------------------------
*        first column and blank line
*-----------------------------------------------------------------------

               k  = 1

            do 145 l = 1, icolm

               if(dum(l).ne.' '.and.dum(l).ne.tub) goto 146

  145       continue

               write(jhs,'()')

               goto 140

  146          k = l

*-----------------------------------------------------------------------
*        last column and blank line
*-----------------------------------------------------------------------

               kl = icolm

            do 165 l = icolm, 1, -1

               if(dum(l).ne.' '.and.dum(l).ne.tub ) goto 166

  165       continue

               write(jhs,'()')

               goto 140

  166          kl = l

*-----------------------------------------------------------------------

            if( ( dum(k  ) .eq. 'I' .or. dum(k  ) .eq. 'i' ) .and.
     &          ( dum(k+1) .eq. 'N' .or. dum(k+1) .eq. 'n' ) .and.
     &          ( dum(k+2) .eq. 'F' .or. dum(k+2) .eq. 'f' ) .and.
     &          ( dum(k+3) .eq. 'L' .or. dum(k+3) .eq. 'l' ) .and.
     &            dum(k+4) .eq. ':' ) then

                  k = k + 4

               call inclf(jsn,jsi,dsin,dum,k,icolm,ill,ilf,idsi,ierr)

                  if( ierr .ne. 0 ) goto 999

               goto 140

            end if

*-----------------------------------------------------------------------

            write(jhs,'(10000a1)') (dum(ic),ic=1,kl)

            goto 140

*-----------------------------------------------------------------------

  141 continue

               k  = 1
               kl = icolm

            do 265 l = icolm, 1, -1

               if(dum(l).ne.' '.and.dum(l).ne.tub ) goto 266

  265       continue

               return

  266          kl = l

            write(jhs,'(10000a1)') (dum(ic),ic=1,kl)

*-----------------------------------------------------------------------

  999 continue

      return
      end


************************************************************************
*                                                                      *
      subroutine ver1to2(jsi,infn,jfname,jsn)
*                                                                      *
*           CONVERT THE TEXT IN VERSION 1.70a TO VERSION 2.00          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character dumr*10000
      character infn(ichrl)*1
      character dsin(0:9)*200

      character yen*1
      character tub*1
      yen = char(92)
      tub = char(9)

*-----------------------------------------------------------------------

            do 20 i = 1, 200

                  dsin(jsn)(i:i) = ' '

   20       continue

*-----------------------------------------------------------------------

            do 58 i = 1, jfname

               dsin(jsn)(i:i) = infn(i)

   58       continue

*-----------------------------------------------------------------------

               i = jfname + 1
               dsin(jsn)(i:i+3) = '.an2'

            open( jhs, file = dsin(jsn), status = 'UNKNOWN' )

*-----------------------------------------------------------------------

               il = 0

  140          il = il + 1


               read(jsi,'(10000a1)', iostat = ios ) (dum(ic),ic=1,icolm)
               if( ios .eq. -1 ) goto 141

               if( ichar(dum(1)) .eq. 26 ) goto 140

*-----------------------------------------------------------------------
*     SKIP BLANK
*-----------------------------------------------------------------------

               k = 1

            do 155 l = 1, icolm

               if(dum(l).ne.' '.and.dum(l).ne.tub) goto 156

  155       continue

               goto 150

  156          k = l

*-----------------------------------------------------------------------
*     LAST COLUMN
*-----------------------------------------------------------------------

  150       continue

               kl = icolm

            do 165 l = icolm, 1, -1

               if(dum(l).ne.' '.and.dum(l).ne.tub.and.
     &            ichar(dum(l)) .ne. 26 ) goto 166

  165       continue

  166          kl = l

*-----------------------------------------------------------------------
*     TRANSFORM THE OLD FORMAT OF TEXT TO THAT OF NEW VERSION 2.00
*-----------------------------------------------------------------------

      if( ( (dum(k).eq.'X'.or.dum(k).eq.'x').and.dum(k+1).eq.':' ) .or.
     &    ( (dum(k).eq.'Y'.or.dum(k).eq.'y').and.dum(k+1).eq.':' ) .or.
     &    ( (dum(k).eq.'Y'.or.dum(k).eq.'y').and.
     &      (dum(k+1).eq.'R'.or.dum(k+1).eq.'r').and.
     &       dum(k+2).eq.':' ) .or.
     &      (dum(k).eq."'") .OR.
     &    ( (dum(k).eq.'H'.or.dum(k).eq.'h').and.dum(k+1).eq.':' ) .or.
     &    ( (dum(k).eq.'W'.or.dum(k).eq.'w').and.dum(k+1).eq.':' ) .or.
     &    ( (dum(k).eq.'A'.or.dum(k).eq.'a').and.
     &      (dum(k+1).eq.'W'.or.dum(k+1).eq.'w').and.
     &       dum(k+2).eq.':' ) ) then

*-----------------------------------------------------------------------

               j  = 0
               jr = 0

               it = 0

*-----------------------------------------------------------------------

  172          j  = j + 1
               if( j .gt. kl ) goto 170

               jr = jr + 1

*-----------------------------------------------------------------------

            if( dum(j) .eq. '$' ) then

               if( it .eq. 0 ) then

                  dumr(jr:jr+4)='{'//yen//'it '
                  jr = jr + 4
                  it = 1

               else

                  dumr(jr:jr)='}'
                  it = 0

               end if

            else if( dum(j) .eq. yen .and. dum(j+1) .eq. '[' ) then

                  dumr(jr:jr) = '['
                  j  = j  + 1

            else if( dum(j) .eq. yen .and. dum(j+1) .eq. ']' ) then

                  dumr(jr:jr) = ']'
                  j  = j  + 1

            else if( dum(j) .eq. '#' .and. dum(j+2) .eq. '#' ) then

               if( dum(j+1) .eq. 'a' ) then

                  dumr(jr:jr+5) = yen//'alpha'
                  jr = jr + 5

               else if( dum(j+1) .eq. 'b' ) then

                  dumr(jr:jr+4) = yen//'beta'
                  jr = jr + 4

               else if( dum(j+1) .eq. 'g' ) then

                  dumr(jr:jr+5) = yen//'gamma'
                  jr = jr + 5

               else if( dum(j+1) .eq. 'd' ) then

                  dumr(jr:jr+5) = yen//'delta'
                  jr = jr + 5

               else if( dum(j+1) .eq. 'e' ) then

                  dumr(jr:jr+7) = yen//'epsilon'
                  jr = jr + 7

               else if( dum(j+1) .eq. 'z' ) then

                  dumr(jr:jr+4) = yen//'zeta'
                  jr = jr + 4

               else if( dum(j+1) .eq. 'h' ) then

                  dumr(jr:jr+3) = yen//'eta'
                  jr = jr + 3

               else if( dum(j+1) .eq. 'q' ) then

                  dumr(jr:jr+5) = yen//'theta'
                  jr = jr + 5

               else if( dum(j+1) .eq. 'i' ) then

                  dumr(jr:jr+4) = yen//'iota'
                  jr = jr + 4

               else if( dum(j+1) .eq. 'k' ) then

                  dumr(jr:jr+5) = yen//'kappa'
                  jr = jr + 5

               else if( dum(j+1) .eq. 'l' ) then

                  dumr(jr:jr+6) = yen//'lambda'
                  jr = jr + 6

               else if( dum(j+1) .eq. 'm' ) then

                  dumr(jr:jr+2) = yen//'mu'
                  jr = jr + 2

               else if( dum(j+1) .eq. 'n' ) then

                  dumr(jr:jr+2) = yen//'nu'
                  jr = jr + 2

               else if( dum(j+1) .eq. 'x' ) then

                  dumr(jr:jr+2) = yen//'xi'
                  jr = jr + 2

               else if( dum(j+1) .eq. 'o' ) then

                  dumr(jr:jr) = 'o'

               else if( dum(j+1) .eq. 'p' ) then

                  dumr(jr:jr+2) = yen//'pi'
                  jr = jr + 2

               else if( dum(j+1) .eq. 'r' ) then

                  dumr(jr:jr+3) = yen//'rho'
                  jr = jr + 3

               else if( dum(j+1) .eq. 's' ) then

                  dumr(jr:jr+5) = yen//'sigma'
                  jr = jr + 5

               else if( dum(j+1) .eq. 't' ) then

                  dumr(jr:jr+3) = yen//'tau'
                  jr = jr + 3

               else if( dum(j+1) .eq. 'u' ) then

                  dumr(jr:jr+7) = yen//'upsilon'
                  jr = jr + 7

               else if( dum(j+1) .eq. 'f' ) then

                  dumr(jr:jr+3) = yen//'phi'
                  jr = jr + 3

               else if( dum(j+1) .eq. 'c' ) then

                  dumr(jr:jr+3) = yen//'chi'
                  jr = jr + 3

               else if( dum(j+1) .eq. 'y' ) then

                  dumr(jr:jr+3) = yen//'psi'
                  jr = jr + 3

               else if( dum(j+1) .eq. 'w' ) then

                  dumr(jr:jr+5) = yen//'omega'
                  jr = jr + 5

               else if( dum(j+1) .eq. 'v' ) then

                  dumr(jr:jr+5) = yen//'varpi'
                  jr = jr + 5

               else if( dum(j+1) .eq. 'V' ) then

                  dumr(jr:jr+8) = yen//'varsigma'
                  jr = jr + 8

               else if( dum(j+1) .eq. 'J' ) then

                  dumr(jr:jr+8) = yen//'vartheta'
                  jr = jr + 8

               else if( dum(j+1) .eq. 'j' ) then

                  dumr(jr:jr+6) = yen//'varphi'
                  jr = jr + 6

               else if( dum(j+1) .eq. 'G' ) then

                  dumr(jr:jr+5) = yen//'Gamma'
                  jr = jr + 5

               else if( dum(j+1) .eq. 'D' ) then

                  dumr(jr:jr+5) = yen//'Delta'
                  jr = jr + 5

               else if( dum(j+1) .eq. 'Q' ) then

                  dumr(jr:jr+5) = yen//'Theta'
                  jr = jr + 5

               else if( dum(j+1) .eq. 'L' ) then

                  dumr(jr:jr+6) = yen//'Lambda'
                  jr = jr + 6

               else if( dum(j+1) .eq. 'X' ) then

                  dumr(jr:jr+2) = yen//'Xi'
                  jr = jr + 2

               else if( dum(j+1) .eq. 'P' ) then

                  dumr(jr:jr+2) = yen//'Pi'
                  jr = jr + 2

               else if( dum(j+1) .eq. 'S' ) then

                  dumr(jr:jr+5) = yen//'Sigma'
                  jr = jr + 5

               else if( dum(j+1) .eq. 'U' ) then

                  dumr(jr:jr+7) = yen//'Upsilon'
                  jr = jr + 7

               else if( dum(j+1) .eq. 'F' ) then

                  dumr(jr:jr+3) = yen//'Phi'
                  jr = jr + 3

               else if( dum(j+1) .eq. 'Y' ) then

                  dumr(jr:jr+3) = yen//'Psi'
                  jr = jr + 3

               else if( dum(j+1) .eq. 'W' ) then

                  dumr(jr:jr+5) = yen//'Omega'
                  jr = jr + 5

               end if

                  j  = j  + 2

            else if( dum(j) .eq. '#' .and. dum(j+1) .eq. yen .and.
     &               dum(j+4) .eq. '#' ) then

               if( dum(j+2) .eq. '5' .and. dum(j+3) .eq. '2' ) then

                  dumr(jr:jr+3) = yen//'ast'
                  jr = jr + 3

               else if( dum(j+2) .eq. '4' .and. dum(j+3) .eq. '4' ) then

                  dumr(jr:jr+6) = yen//'exists'
                  jr = jr + 6

               end if

                  j  = j  + 4

            else if( dum(j) .eq. '#' .and. dum(j+1) .eq. yen .and.
     &               dum(j+5) .eq. '#' ) then

               if( dum(j+2) .eq. '2' .and.
     &             dum(j+3) .eq. '6' .and.
     &             dum(j+4) .eq. '1' ) then

                  dumr(jr:jr+2) = yen//'pm'
                  jr = jr + 2

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '6' .and.
     &                  dum(j+4) .eq. '0' ) then

                  dumr(jr:jr+4) = yen//'circ'
                  jr = jr + 4

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '0' .and.
     &                  dum(j+4) .eq. '5' ) then

                  dumr(jr:jr+5) = yen//'oplus'
                  jr = jr + 5

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '1' .and.
     &                  dum(j+4) .eq. '0' ) then

                  dumr(jr:jr+3) = yen//'cup'
                  jr = jr + 3

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '6' .and.
     &                  dum(j+4) .eq. '4' ) then

                  dumr(jr:jr+5) = yen//'times'
                  jr = jr + 5

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '6' .and.
     &                  dum(j+4) .eq. '7' ) then

                  dumr(jr:jr+6) = yen//'bullet'
                  jr = jr + 6

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '3' .and.
     &                  dum(j+4) .eq. '2' ) then

                  dumr(jr:jr+3) = yen//'vee'
                  jr = jr + 3

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '0' .and.
     &                  dum(j+4) .eq. '4' ) then

                  dumr(jr:jr+6) = yen//'otimes'
                  jr = jr + 6

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '7' .and.
     &                  dum(j+4) .eq. '0' ) then

                  dumr(jr:jr+3) = yen//'div'
                  jr = jr + 3

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '2' .and.
     &                  dum(j+4) .eq. '7' ) then

                  dumr(jr:jr+4) = yen//'cdot'
                  jr = jr + 4

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '3' .and.
     &                  dum(j+4) .eq. '1' ) then

                  dumr(jr:jr+5) = yen//'wedge'
                  jr = jr + 5

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '0' .and.
     &                  dum(j+4) .eq. '6' ) then

                  dumr(jr:jr+6) = yen//'oslash'
                  jr = jr + 6

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '0' .and.
     &                  dum(j+4) .eq. '7' ) then

                  dumr(jr:jr+3) = yen//'cap'
                  jr = jr + 3

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '4' .and.
     &                  dum(j+4) .eq. '0' ) then

                  dumr(jr:jr+7) = yen//'diamond'
                  jr = jr + 7

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '7' .and.
     &                  dum(j+4) .eq. '2' ) then

                  dumr(jr:jr+5) = yen//'equiv'
                  jr = jr + 5

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '7' .and.
     &                  dum(j+4) .eq. '3' ) then

                  dumr(jr:jr+6) = yen//'approx'
                  jr = jr + 6

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '7' .and.
     &                  dum(j+4) .eq. '1' ) then

                  dumr(jr:jr+3) = yen//'neq'
                  jr = jr + 3

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '6' .and.
     &                  dum(j+4) .eq. '5' ) then

                  dumr(jr:jr+6) = yen//'propto'
                  jr = jr + 6

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '4' .and.
     &                  dum(j+4) .eq. '3' ) then

                  dumr(jr:jr+2) = yen//'le'
                  jr = jr + 2

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '1' .and.
     &                  dum(j+4) .eq. '4' ) then

                  dumr(jr:jr+6) = yen//'subset'
                  jr = jr + 6

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '1' .and.
     &                  dum(j+4) .eq. '5' ) then

                  dumr(jr:jr+8) = yen//'subseteq'
                  jr = jr + 8

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '1' .and.
     &                  dum(j+4) .eq. '6' ) then

                  dumr(jr:jr+2) = yen//'in'
                  jr = jr + 2

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '1' .and.
     &                  dum(j+4) .eq. '7' ) then

                  dumr(jr:jr+5) = yen//'notin'
                  jr = jr + 5

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '6' .and.
     &                  dum(j+4) .eq. '3' ) then

                  dumr(jr:jr+2) = yen//'ge'
                  jr = jr + 2

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '1' .and.
     &                  dum(j+4) .eq. '1' ) then

                  dumr(jr:jr+6) = yen//'supset'
                  jr = jr + 6

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '1' .and.
     &                  dum(j+4) .eq. '2' ) then

                  dumr(jr:jr+8) = yen//'supseteq'
                  jr = jr + 8

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '5' .and.
     &                  dum(j+4) .eq. '3' ) then

                  dumr(jr:jr+14) = yen//'leftrightarrow'
                  jr = jr + 14

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '5' .and.
     &                  dum(j+4) .eq. '4' ) then

                  dumr(jr:jr+9) = yen//'leftarrow'
                  jr = jr + 9

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '5' .and.
     &                  dum(j+4) .eq. '6' ) then

                  dumr(jr:jr+2) = yen//'to'
                  jr = jr + 2

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '7' .and.
     &                  dum(j+4) .eq. '7' ) then

                  dumr(jr:jr+13) = yen//'hookleftarrow'
                  jr = jr + 13

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '3' .and.
     &                  dum(j+4) .eq. '3' ) then

                  dumr(jr:jr+14) = yen//'Leftrightarrow'
                  jr = jr + 14

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '3' .and.
     &                  dum(j+4) .eq. '4' ) then

                  dumr(jr:jr+9) = yen//'Leftarrow'
                  jr = jr + 9

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '3' .and.
     &                  dum(j+4) .eq. '6' ) then

                  dumr(jr:jr+10) = yen//'Rightarrow'
                  jr = jr + 10

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '7' .and.
     &                  dum(j+4) .eq. '5' ) then

                  dumr(jr:jr+3) = yen//'mid'
                  jr = jr + 3

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '5' .and.
     &                  dum(j+4) .eq. '5' ) then

                  dumr(jr:jr+7) = yen//'uparrow'
                  jr = jr + 7

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '5' .and.
     &                  dum(j+4) .eq. '7' ) then

                  dumr(jr:jr+9) = yen//'downarrow'
                  jr = jr + 9

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '3' .and.
     &                  dum(j+4) .eq. '5' ) then

                  dumr(jr:jr+7) = yen//'Uparrow'
                  jr = jr + 7

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '3' .and.
     &                  dum(j+4) .eq. '7' ) then

                  dumr(jr:jr+9) = yen//'Downarrow'
                  jr = jr + 9

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '0' .and.
     &                  dum(j+4) .eq. '0' ) then

                  dumr(jr:jr+5) = yen//'aleph'
                  jr = jr + 5

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '0' .and.
     &                  dum(j+4) .eq. '1' ) then

                  dumr(jr:jr+2) = yen//'Im'
                  jr = jr + 2

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '0' .and.
     &                  dum(j+4) .eq. '2' ) then

                  dumr(jr:jr+2) = yen//'Re'
                  jr = jr + 2

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '0' .and.
     &                  dum(j+4) .eq. '3' ) then

                  dumr(jr:jr+2) = yen//'wp'
                  jr = jr + 2

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '6' .and.
     &                  dum(j+4) .eq. '6' ) then

                  dumr(jr:jr+7) = yen//'partial'
                  jr = jr + 7

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '4' .and.
     &                  dum(j+4) .eq. '5' ) then

                  dumr(jr:jr+5) = yen//'infty'
                  jr = jr + 5

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '4' .and.
     &                  dum(j+4) .eq. '2' ) then

                  dumr(jr:jr+5) = yen//'prime'
                  jr = jr + 5

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '2' .and.
     &                  dum(j+4) .eq. '1' ) then

                  dumr(jr:jr+5) = yen//'nabla'
                  jr = jr + 5

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '2' .and.
     &                  dum(j+4) .eq. '6' ) then

                  dumr(jr:jr+4) = yen//'surd'
                  jr = jr + 4

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '2' .and.
     &                  dum(j+4) .eq. '0' ) then

                  dumr(jr:jr+5) = yen//'angle'
                  jr = jr + 5

               else if( dum(j+2) .eq. '3' .and.
     &                  dum(j+3) .eq. '3' .and.
     &                  dum(j+4) .eq. '0' ) then

                  dumr(jr:jr+3) = yen//'neg'
                  jr = jr + 3

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '4' .and.
     &                  dum(j+4) .eq. '7' ) then

                  dumr(jr:jr+8) = yen//'clubsuit'
                  jr = jr + 8

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '5' .and.
     &                  dum(j+4) .eq. '0' ) then

                  dumr(jr:jr+11) = yen//'diamondsuit'
                  jr = jr + 11

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '5' .and.
     &                  dum(j+4) .eq. '1' ) then

                  dumr(jr:jr+9) = yen//'heartsuit'
                  jr = jr + 9

               else if( dum(j+2) .eq. '2' .and.
     &                  dum(j+3) .eq. '5' .and.
     &                  dum(j+4) .eq. '2' ) then

                  dumr(jr:jr+9) = yen//'spadesuit'
                  jr = jr + 9

               end if

                  j  = j  + 5

            else if( dum(j) .eq. yen .and. dum(j+1) .eq. '3' .and.
     &             ( dum(j+2) .eq. '2' .or. dum(j+2) .eq. '3' ) ) then

               if( dum(j+2) .eq. '2' .and. dum(j+3) .eq. '1' ) then

                  dumr(jr:jr+2) = yen//'.a'

               else if( dum(j+2) .eq. '3' .and. dum(j+3) .eq. '1' ) then

                  dumr(jr:jr+2) = yen//'.A'

               else if( dum(j+2) .eq. '2' .and. dum(j+3) .eq. '2' ) then

                  dumr(jr:jr+2) = yen//'"a'

               else if( dum(j+2) .eq. '3' .and. dum(j+3) .eq. '2' ) then

                  dumr(jr:jr+2) = yen//'"A'

               else if( dum(j+2) .eq. '2' .and. dum(j+3) .eq. '3' ) then

                  dumr(jr:jr+2) = yen//'"o'

               else if( dum(j+2) .eq. '3' .and. dum(j+3) .eq. '3' ) then

                  dumr(jr:jr+2) = yen//'"O'

               else if( dum(j+2) .eq. '2' .and. dum(j+3) .eq. '4' ) then

                  dumr(jr:jr+2) = yen//'"u'

               else if( dum(j+2) .eq. '3' .and. dum(j+3) .eq. '4' ) then

                  dumr(jr:jr+2) = yen//'"U'

               end if

                  jr = jr + 2
                  j  = j  + 3

            else

               dumr(jr:jr) = dum(j)

            end if

               goto 172

*-----------------------------------------------------------------------

  170    continue

         do 171 j = 1, jr

               dum(j) = dumr(j:j)

  171    continue

               kl = jr

      end if

*-----------------------------------------------------------------------
*     REWITE THE SOURCE
*-----------------------------------------------------------------------

               write(jhs,'(10000a1)') (dum(ic),ic=1,kl)

               goto 140

*-----------------------------------------------------------------------

  141 continue

            close(jhs)

            call closef(jsi,jsn)
            call openf(jsi,jsn,dsin)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine strwhb(icn,chaf,fts,iftn,strw,strh,strb,ifon)
*                                                                      *
*           CALCULATE THE STRING WIDTH, HIGHT AND BOTTOM SIZE          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      common /stwhb/  stwhb(0:13,3,0:255)
      common /stwhbk/ stwhbk(3)

      character chaf(ichrl)*1

      logical spjpn
      logical numec

*-----------------------------------------------------------------------

      character yen*1

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------

            yen = char(92)

*-----------------------------------------------------------------------

         strw =  0.0
         strh = -100000.0
         strb =  100000.0

*-----------------------------------------------------------------------

         k = 0

  100    k = k + 1

*-----------------------------------------------------------------------

         if( spjpn(k,icn,icf,chaf,ifon,njpn) ) then

*              SKIP JAPANESE TWO BITE CHARACTERS

*-----------------------------------------------------------------------

               strw = strw + dble(njpn) * stwhbk(1) * fts / 1000.0

               strh = max( strh, stwhbk(2) * fts / 1000.0 )
               strb = min( strb, stwhbk(3) * fts / 1000.0 )

               k = icf

*-----------------------------------------------------------------------

         else

*-----------------------------------------------------------------------

               if     ( chaf(k)   .eq. yen .and.
     &                ( chaf(k+1) .eq. yen .or.
     &                  chaf(k+1) .eq. '(' .or.
     &                  chaf(k+1) .eq. ')' ) ) then

                  k     = k + 1
                  ichnm = ichar( chaf(k) )

               else if( chaf(k) .eq. yen .and.
     &                  numec(chaf(k+1)) .and.
     &                  numec(chaf(k+2)) .and.
     &                  numec(chaf(k+3)) ) then

                        read(chaf(k+1),'(I1)') ir2
                        read(chaf(k+2),'(I1)') ir1
                        read(chaf(k+3),'(I1)') ir0

                  k     = k + 3
                  ichnm = ir2 * 8**2 + ir1 * 8 + ir0

               else if( chaf(k) .eq. yen .and.
     &                  numec(chaf(k+1)) .and.
     &                  numec(chaf(k+2)) ) then

                        read(chaf(k+1),'(I1)') ir1
                        read(chaf(k+2),'(I1)') ir0

                  k     = k + 2
                  ichnm = ir1 * 8 + ir0

               else

                  ichnm = ichar( chaf(k) )

               end if

*-----------------------------------------------------------------------

                  ichmm = max(0,ichnm)
                  ichnm = min(255,ichmm)

                  strw = strw +     stwhb(iftn,1,ichnm) * fts / 1000.0
                  strh = max( strh, stwhb(iftn,2,ichnm) * fts / 1000.0 )
                  strb = min( strb, stwhb(iftn,3,ichnm) * fts / 1000.0 )

         end if

*-----------------------------------------------------------------------

            if( k .lt. icn ) goto 100

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      block data
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /stwhb/  stwhb(0:13,3,0:255)
      common /stwrl/  stwrl4(3,4)
      common /stwhbk/ stwhbk(3)

*-----------------------------------------------------------------------

      data  stwhbk/
     &    900.,  810.,  -90./

*-----------------------------------------------------------------------

      data  stwrl4/
     &    549.,   11.,  535.,
     &   1000.,  -60., 1050.,
     &    549.,   10.,  515.,
     &    500.,  480., 1090./

*-----------------------------------------------------------------------

      data  ((stwhb(0,i,j),i=1,3),j=0,31)/
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

      data  ((stwhb(0,i,j),i=1,3),j=32,63)/
     &    278.,    0.,    0.,
     &    278.,  718.,    0.,
     &    355.,  718.,  463.,
     &    556.,  688.,    0.,
     &    556.,  775., -115.,
     &    889.,  703.,  -19.,
     &    667.,  718.,  -15.,
     &    222.,  718.,  463.,
     &    333.,  733., -207.,
     &    333.,  733., -207.,
     &    389.,  718.,  431.,
     &    584.,  505.,    0.,
     &    278.,  106., -147.,
     &    333.,  322.,  232.,
     &    278.,  106.,    0.,
     &    278.,  737.,  -19.,
     &    556.,  703.,  -19.,
     &    556.,  703.,    0.,
     &    556.,  703.,    0.,
     &    556.,  703.,  -19.,
     &    556.,  703.,    0.,
     &    556.,  688.,  -19.,
     &    556.,  703.,  -19.,
     &    556.,  688.,    0.,
     &    556.,  703.,  -19.,
     &    556.,  703.,  -19.,
     &    278.,  516.,    0.,
     &    278.,  516., -147.,
     &    584.,  495.,   11.,
     &    584.,  390.,  115.,
     &    584.,  495.,   11.,
     &    556.,  727.,    0./

      data  ((stwhb(0,i,j),i=1,3),j=64,95)/
     &   1015.,  737.,  -19.,
     &    667.,  718.,    0.,
     &    667.,  718.,    0.,
     &    722.,  737.,  -19.,
     &    722.,  718.,    0.,
     &    667.,  718.,    0.,
     &    611.,  718.,    0.,
     &    778.,  737.,  -19.,
     &    722.,  718.,    0.,
     &    278.,  718.,    0.,
     &    500.,  718.,  -19.,
     &    667.,  718.,    0.,
     &    556.,  718.,    0.,
     &    833.,  718.,    0.,
     &    722.,  718.,    0.,
     &    778.,  737.,  -19.,
     &    667.,  718.,    0.,
     &    778.,  737.,  -56.,
     &    722.,  718.,    0.,
     &    667.,  737.,  -19.,
     &    611.,  718.,    0.,
     &    722.,  718.,  -19.,
     &    667.,  718.,    0.,
     &    944.,  718.,    0.,
     &    667.,  718.,    0.,
     &    667.,  718.,    0.,
     &    611.,  718.,    0.,
     &    278.,  722., -196.,
     &    278.,  737.,  -19.,
     &    278.,  722., -196.,
     &    469.,  688.,  264.,
     &    556.,  -75., -125./

      data  ((stwhb(0,i,j),i=1,3),j=96,127)/
     &    222.,  725.,  470.,
     &    556.,  538.,  -15.,
     &    556.,  718.,  -15.,
     &    500.,  538.,  -15.,
     &    556.,  718.,  -15.,
     &    556.,  538.,  -15.,
     &    278.,  728.,    0.,
     &    556.,  538., -220.,
     &    556.,  718.,    0.,
     &    222.,  718.,    0.,
     &    222.,  718., -220.,
     &    500.,  718.,    0.,
     &    222.,  718.,    0.,
     &    833.,  538.,    0.,
     &    556.,  538.,    0.,
     &    556.,  538.,  -14.,
     &    556.,  538., -207.,
     &    556.,  538., -207.,
     &    333.,  538.,    0.,
     &    500.,  538.,  -15.,
     &    278.,  669.,   -7.,
     &    556.,  523.,  -15.,
     &    500.,  523.,    0.,
     &    722.,  523.,    0.,
     &    500.,  523.,    0.,
     &    500.,  523., -214.,
     &    500.,  523.,    0.,
     &    334.,  722., -196.,
     &    260.,  737.,  -19.,
     &    334.,  722., -196.,
     &    584.,  326.,  180.,
     &    278.,    0.,    0./

      data  ((stwhb(0,i,j),i=1,3),j=128,159)/
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

      data  ((stwhb(0,i,j),i=1,3),j=160,191)/
     &    278.,    0.,    0.,
     &    333.,  523., -195.,
     &    556.,  623., -115.,
     &    556.,  718.,  -16.,
     &    167.,  703.,  -19.,
     &    556.,  688.,    0.,
     &    556.,  737., -207.,
     &    556.,  737., -191.,
     &    556.,  603.,   99.,
     &    191.,  718.,  463.,
     &    333.,  725.,  470.,
     &    556.,  446.,  108.,
     &    333.,  446.,  108.,
     &    333.,  446.,  108.,
     &    500.,  728.,    0.,
     &    500.,  728.,    0.,
     &    278.,    0.,    0.,
     &    556.,  313.,  240.,
     &    556.,  718., -159.,
     &    556.,  718., -159.,
     &    278.,  315.,  190.,
     &    278.,    0.,    0.,
     &    537.,  718., -173.,
     &    350.,  517.,  202.,
     &    220.,  106., -149.,
     &    333.,  106., -149.,
     &    333.,  718.,  463.,
     &    556.,  446.,  108.,
     &   1000.,  106.,    0.,
     &   1000.,  703.,  -19.,
     &    278.,    0.,    0.,
     &    611.,  525., -201./

      data  ((stwhb(0,i,j),i=1,3),j=192,223)/
     &    278.,    0.,    0.,
     &    333.,  734.,  593.,
     &    333.,  734.,  593.,
     &    333.,  734.,  593.,
     &    333.,  722.,  606.,
     &    333.,  684.,  627.,
     &    333.,  731.,  595.,
     &    333.,  706.,  604.,
     &    333.,  706.,  604.,
     &    278.,    0.,    0.,
     &    333.,  756.,  572.,
     &    333.,    0., -225.,
     &    278.,    0.,    0.,
     &    333.,  734.,  593.,
     &    333.,    0., -225.,
     &    333.,  734.,  593.,
     &   1000.,  313.,  240.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

      data  ((stwhb(0,i,j),i=1,3),j=224,255)/
     &    278.,    0.,    0.,
     &   1000.,  718.,    0.,
     &    278.,    0.,    0.,
     &    370.,  737.,  304.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    556.,  718.,    0.,
     &    778.,  737.,  -19.,
     &   1000.,  737.,  -19.,
     &    365.,  737.,  304.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    889.,  538.,  -15.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,  523.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    222.,  718.,    0.,
     &    611.,  545.,  -22.,
     &    944.,  538.,  -15.,
     &    611.,  728.,  -15.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(1,i,j),i=1,3),j=0,31)/
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

      data  ((stwhb(1,i,j),i=1,3),j=32,63)/
     &    278.,    0.,    0.,
     &    278.,  718.,    0.,
     &    355.,  718.,  463.,
     &    556.,  688.,    0.,
     &    556.,  775., -115.,
     &    889.,  703.,  -19.,
     &    667.,  718.,  -15.,
     &    222.,  718.,  463.,
     &    333.,  733., -207.,
     &    333.,  733., -207.,
     &    389.,  718.,  431.,
     &    584.,  505.,    0.,
     &    278.,  106., -147.,
     &    333.,  322.,  232.,
     &    278.,  106.,    0.,
     &    278.,  737.,  -19.,
     &    556.,  703.,  -19.,
     &    556.,  703.,    0.,
     &    556.,  703.,    0.,
     &    556.,  703.,  -19.,
     &    556.,  703.,    0.,
     &    556.,  688.,  -19.,
     &    556.,  703.,  -19.,
     &    556.,  688.,    0.,
     &    556.,  703.,  -19.,
     &    556.,  703.,  -19.,
     &    278.,  516.,    0.,
     &    278.,  516., -147.,
     &    584.,  495.,   11.,
     &    584.,  390.,  115.,
     &    584.,  495.,   11.,
     &    556.,  727.,    0./

      data  ((stwhb(1,i,j),i=1,3),j=64,95)/
     &   1015.,  737.,  -19.,
     &    667.,  718.,    0.,
     &    667.,  718.,    0.,
     &    722.,  737.,  -19.,
     &    722.,  718.,    0.,
     &    667.,  718.,    0.,
     &    611.,  718.,    0.,
     &    778.,  737.,  -19.,
     &    722.,  718.,    0.,
     &    278.,  718.,    0.,
     &    500.,  718.,  -19.,
     &    667.,  718.,    0.,
     &    556.,  718.,    0.,
     &    833.,  718.,    0.,
     &    722.,  718.,    0.,
     &    778.,  737.,  -19.,
     &    667.,  718.,    0.,
     &    778.,  737.,  -56.,
     &    722.,  718.,    0.,
     &    667.,  737.,  -19.,
     &    611.,  718.,    0.,
     &    722.,  718.,  -19.,
     &    667.,  718.,    0.,
     &    944.,  718.,    0.,
     &    667.,  718.,    0.,
     &    667.,  718.,    0.,
     &    611.,  718.,    0.,
     &    278.,  722., -196.,
     &    278.,  737.,  -19.,
     &    278.,  722., -196.,
     &    469.,  688.,  264.,
     &    556.,  -75., -125./

      data  ((stwhb(1,i,j),i=1,3),j=96,127)/
     &    222.,  725.,  470.,
     &    556.,  538.,  -15.,
     &    556.,  718.,  -15.,
     &    500.,  538.,  -15.,
     &    556.,  718.,  -15.,
     &    556.,  538.,  -15.,
     &    278.,  728.,    0.,
     &    556.,  538., -220.,
     &    556.,  718.,    0.,
     &    222.,  718.,    0.,
     &    222.,  718., -220.,
     &    500.,  718.,    0.,
     &    222.,  718.,    0.,
     &    833.,  538.,    0.,
     &    556.,  538.,    0.,
     &    556.,  538.,  -14.,
     &    556.,  538., -207.,
     &    556.,  538., -207.,
     &    333.,  538.,    0.,
     &    500.,  538.,  -15.,
     &    278.,  669.,   -7.,
     &    556.,  523.,  -15.,
     &    500.,  523.,    0.,
     &    722.,  523.,    0.,
     &    500.,  523.,    0.,
     &    500.,  523., -214.,
     &    500.,  523.,    0.,
     &    334.,  722., -196.,
     &    260.,  737.,  -19.,
     &    334.,  722., -196.,
     &    584.,  326.,  180.,
     &    278.,    0.,    0./

      data  ((stwhb(1,i,j),i=1,3),j=128,159)/
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

      data  ((stwhb(1,i,j),i=1,3),j=160,191)/
     &    278.,    0.,    0.,
     &    333.,  523., -195.,
     &    556.,  623., -115.,
     &    556.,  718.,  -16.,
     &    167.,  703.,  -19.,
     &    556.,  688.,    0.,
     &    556.,  737., -207.,
     &    556.,  737., -191.,
     &    556.,  603.,   99.,
     &    191.,  718.,  463.,
     &    333.,  725.,  470.,
     &    556.,  446.,  108.,
     &    333.,  446.,  108.,
     &    333.,  446.,  108.,
     &    500.,  728.,    0.,
     &    500.,  728.,    0.,
     &    278.,    0.,    0.,
     &    556.,  313.,  240.,
     &    556.,  718., -159.,
     &    556.,  718., -159.,
     &    278.,  315.,  190.,
     &    278.,    0.,    0.,
     &    537.,  718., -173.,
     &    350.,  517.,  202.,
     &    220.,  106., -149.,
     &    333.,  106., -149.,
     &    333.,  718.,  463.,
     &    556.,  446.,  108.,
     &   1000.,  106.,    0.,
     &   1000.,  703.,  -19.,
     &    278.,    0.,    0.,
     &    611.,  525., -201./

      data  ((stwhb(1,i,j),i=1,3),j=192,223)/
     &    278.,    0.,    0.,
     &    333.,  734.,  593.,
     &    333.,  734.,  593.,
     &    333.,  734.,  593.,
     &    333.,  722.,  606.,
     &    333.,  684.,  627.,
     &    333.,  731.,  595.,
     &    333.,  706.,  604.,
     &    333.,  706.,  604.,
     &    278.,    0.,    0.,
     &    333.,  756.,  572.,
     &    333.,    0., -225.,
     &    278.,    0.,    0.,
     &    333.,  734.,  593.,
     &    333.,    0., -225.,
     &    333.,  734.,  593.,
     &   1000.,  313.,  240.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

      data  ((stwhb(1,i,j),i=1,3),j=224,255)/
     &    278.,    0.,    0.,
     &   1000.,  718.,    0.,
     &    278.,    0.,    0.,
     &    370.,  737.,  304.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    556.,  718.,    0.,
     &    778.,  737.,  -19.,
     &   1000.,  737.,  -19.,
     &    365.,  737.,  304.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    889.,  538.,  -15.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,  523.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    222.,  718.,    0.,
     &    611.,  545.,  -22.,
     &    944.,  538.,  -15.,
     &    611.,  728.,  -15.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(2,i,j),i=1,3),j=0,31)/
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

      data  ((stwhb(2,i,j),i=1,3),j=32,63)/
     &    278.,    0.,    0.,
     &    333.,  718.,    0.,
     &    474.,  718.,  447.,
     &    556.,  698.,    0.,
     &    556.,  775., -115.,
     &    889.,  710.,  -19.,
     &    722.,  718.,  -19.,
     &    278.,  718.,  445.,
     &    333.,  734., -208.,
     &    333.,  734., -208.,
     &    389.,  718.,  387.,
     &    584.,  506.,    0.,
     &    278.,  146., -168.,
     &    333.,  345.,  215.,
     &    278.,  146.,    0.,
     &    278.,  737.,  -19.,
     &    556.,  710.,  -19.,
     &    556.,  710.,    0.,
     &    556.,  710.,    0.,
     &    556.,  710.,  -19.,
     &    556.,  710.,    0.,
     &    556.,  698.,  -19.,
     &    556.,  710.,  -19.,
     &    556.,  698.,    0.,
     &    556.,  710.,  -19.,
     &    556.,  710.,  -19.,
     &    333.,  512.,    0.,
     &    333.,  512., -168.,
     &    584.,  514.,   -8.,
     &    584.,  419.,   87.,
     &    584.,  514.,   -8.,
     &    611.,  727.,    0./

      data  ((stwhb(2,i,j),i=1,3),j=64,95)/
     &    975.,  737.,  -19.,
     &    722.,  718.,    0.,
     &    722.,  718.,    0.,
     &    722.,  737.,  -19.,
     &    722.,  718.,    0.,
     &    667.,  718.,    0.,
     &    611.,  718.,    0.,
     &    778.,  737.,  -19.,
     &    722.,  718.,    0.,
     &    278.,  718.,    0.,
     &    556.,  718.,  -18.,
     &    722.,  718.,    0.,
     &    611.,  718.,    0.,
     &    833.,  718.,    0.,
     &    722.,  718.,    0.,
     &    778.,  737.,  -19.,
     &    667.,  718.,    0.,
     &    778.,  737.,  -52.,
     &    722.,  718.,    0.,
     &    667.,  737.,  -19.,
     &    611.,  718.,    0.,
     &    722.,  718.,  -19.,
     &    667.,  718.,    0.,
     &    944.,  718.,    0.,
     &    667.,  718.,    0.,
     &    667.,  718.,    0.,
     &    611.,  718.,    0.,
     &    333.,  722., -196.,
     &    278.,  737.,  -19.,
     &    333.,  722., -196.,
     &    584.,  698.,  323.,
     &    556.,  -75., -125./

      data  ((stwhb(2,i,j),i=1,3),j=96,127)/
     &    278.,  727.,  454.,
     &    556.,  546.,  -14.,
     &    611.,  718.,  -14.,
     &    556.,  546.,  -14.,
     &    611.,  718.,  -14.,
     &    556.,  546.,  -14.,
     &    333.,  727.,    0.,
     &    611.,  546., -217.,
     &    611.,  718.,    0.,
     &    278.,  725.,    0.,
     &    278.,  725., -224.,
     &    556.,  718.,    0.,
     &    278.,  718.,    0.,
     &    889.,  546.,    0.,
     &    611.,  546.,    0.,
     &    611.,  546.,  -14.,
     &    611.,  546., -207.,
     &    611.,  546., -207.,
     &    389.,  546.,    0.,
     &    556.,  546.,  -14.,
     &    333.,  676.,   -6.,
     &    611.,  532.,  -14.,
     &    556.,  532.,    0.,
     &    778.,  532.,    0.,
     &    556.,  532.,    0.,
     &    556.,  532., -214.,
     &    500.,  532.,    0.,
     &    389.,  722., -196.,
     &    280.,  737.,  -19.,
     &    389.,  722., -196.,
     &    584.,  343.,  163.,
     &    278.,    0.,    0./

      data  ((stwhb(2,i,j),i=1,3),j=128,159)/
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

      data  ((stwhb(2,i,j),i=1,3),j=160,191)/
     &    278.,    0.,    0.,
     &    333.,  532., -186.,
     &    556.,  628., -118.,
     &    556.,  718.,  -16.,
     &    167.,  710.,  -19.,
     &    556.,  698.,    0.,
     &    556.,  737., -210.,
     &    556.,  727., -184.,
     &    556.,  636.,   76.,
     &    238.,  718.,  447.,
     &    500.,  727.,  454.,
     &    556.,  484.,   76.,
     &    333.,  484.,   76.,
     &    333.,  484.,   76.,
     &    611.,  727.,    0.,
     &    611.,  727.,    0.,
     &    278.,    0.,    0.,
     &    556.,  313.,  227.,
     &    556.,  718., -171.,
     &    556.,  718., -171.,
     &    278.,  334.,  172.,
     &    278.,    0.,    0.,
     &    556.,  700., -191.,
     &    350.,  524.,  194.,
     &    278.,  127., -146.,
     &    500.,  127., -146.,
     &    500.,  718.,  445.,
     &    556.,  484.,   76.,
     &   1000.,  146.,    0.,
     &   1000.,  710.,  -19.,
     &    278.,    0.,    0.,
     &    611.,  532., -195./

      data  ((stwhb(2,i,j),i=1,3),j=192,223)/
     &    278.,    0.,    0.,
     &    333.,  750.,  604.,
     &    333.,  750.,  604.,
     &    333.,  750.,  604.,
     &    333.,  737.,  610.,
     &    333.,  678.,  604.,
     &    333.,  750.,  604.,
     &    333.,  729.,  614.,
     &    333.,  729.,  614.,
     &    278.,    0.,    0.,
     &    333.,  776.,  568.,
     &    333.,    0., -228.,
     &    278.,    0.,    0.,
     &    333.,  750.,  604.,
     &    333.,    0., -228.,
     &    333.,  750.,  604.,
     &   1000.,  333.,  227.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

      data  ((stwhb(2,i,j),i=1,3),j=224,255)/
     &    278.,    0.,    0.,
     &   1000.,  718.,    0.,
     &    278.,    0.,    0.,
     &    370.,  737.,  276.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    611.,  718.,    0.,
     &    778.,  745.,  -27.,
     &   1000.,  737.,  -19.,
     &    365.,  737.,  276.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    889.,  546.,  -14.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,  532.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,  718.,    0.,
     &    611.,  560.,  -29.,
     &    944.,  546.,  -14.,
     &    611.,  731.,  -14.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(3,i,j),i=1,3),j=0,31)/
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

      data  ((stwhb(3,i,j),i=1,3),j=32,63)/
     &    278.,    0.,    0.,
     &    333.,  718.,    0.,
     &    474.,  718.,  447.,
     &    556.,  698.,    0.,
     &    556.,  775., -115.,
     &    889.,  710.,  -19.,
     &    722.,  718.,  -19.,
     &    278.,  718.,  445.,
     &    333.,  734., -208.,
     &    333.,  734., -208.,
     &    389.,  718.,  387.,
     &    584.,  506.,    0.,
     &    278.,  146., -168.,
     &    333.,  345.,  215.,
     &    278.,  146.,    0.,
     &    278.,  737.,  -19.,
     &    556.,  710.,  -19.,
     &    556.,  710.,    0.,
     &    556.,  710.,    0.,
     &    556.,  710.,  -19.,
     &    556.,  710.,    0.,
     &    556.,  698.,  -19.,
     &    556.,  710.,  -19.,
     &    556.,  698.,    0.,
     &    556.,  710.,  -19.,
     &    556.,  710.,  -19.,
     &    333.,  512.,    0.,
     &    333.,  512., -168.,
     &    584.,  514.,   -8.,
     &    584.,  419.,   87.,
     &    584.,  514.,   -8.,
     &    611.,  727.,    0./

      data  ((stwhb(3,i,j),i=1,3),j=64,95)/
     &    975.,  737.,  -19.,
     &    722.,  718.,    0.,
     &    722.,  718.,    0.,
     &    722.,  737.,  -19.,
     &    722.,  718.,    0.,
     &    667.,  718.,    0.,
     &    611.,  718.,    0.,
     &    778.,  737.,  -19.,
     &    722.,  718.,    0.,
     &    278.,  718.,    0.,
     &    556.,  718.,  -18.,
     &    722.,  718.,    0.,
     &    611.,  718.,    0.,
     &    833.,  718.,    0.,
     &    722.,  718.,    0.,
     &    778.,  737.,  -19.,
     &    667.,  718.,    0.,
     &    778.,  737.,  -52.,
     &    722.,  718.,    0.,
     &    667.,  737.,  -19.,
     &    611.,  718.,    0.,
     &    722.,  718.,  -19.,
     &    667.,  718.,    0.,
     &    944.,  718.,    0.,
     &    667.,  718.,    0.,
     &    667.,  718.,    0.,
     &    611.,  718.,    0.,
     &    333.,  722., -196.,
     &    278.,  737.,  -19.,
     &    333.,  722., -196.,
     &    584.,  698.,  323.,
     &    556.,  -75., -125./

      data  ((stwhb(3,i,j),i=1,3),j=96,127)/
     &    278.,  727.,  454.,
     &    556.,  546.,  -14.,
     &    611.,  718.,  -14.,
     &    556.,  546.,  -14.,
     &    611.,  718.,  -14.,
     &    556.,  546.,  -14.,
     &    333.,  727.,    0.,
     &    611.,  546., -217.,
     &    611.,  718.,    0.,
     &    278.,  725.,    0.,
     &    278.,  725., -224.,
     &    556.,  718.,    0.,
     &    278.,  718.,    0.,
     &    889.,  546.,    0.,
     &    611.,  546.,    0.,
     &    611.,  546.,  -14.,
     &    611.,  546., -207.,
     &    611.,  546., -207.,
     &    389.,  546.,    0.,
     &    556.,  546.,  -14.,
     &    333.,  676.,   -6.,
     &    611.,  532.,  -14.,
     &    556.,  532.,    0.,
     &    778.,  532.,    0.,
     &    556.,  532.,    0.,
     &    556.,  532., -214.,
     &    500.,  532.,    0.,
     &    389.,  722., -196.,
     &    280.,  737.,  -19.,
     &    389.,  722., -196.,
     &    584.,  343.,  163.,
     &    278.,    0.,    0./

      data  ((stwhb(3,i,j),i=1,3),j=128,159)/
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

      data  ((stwhb(3,i,j),i=1,3),j=160,191)/
     &    278.,    0.,    0.,
     &    333.,  532., -186.,
     &    556.,  628., -118.,
     &    556.,  718.,  -16.,
     &    167.,  710.,  -19.,
     &    556.,  698.,    0.,
     &    556.,  737., -210.,
     &    556.,  727., -184.,
     &    556.,  636.,   76.,
     &    238.,  718.,  447.,
     &    500.,  727.,  454.,
     &    556.,  484.,   76.,
     &    333.,  484.,   76.,
     &    333.,  484.,   76.,
     &    611.,  727.,    0.,
     &    611.,  727.,    0.,
     &    278.,    0.,    0.,
     &    556.,  313.,  227.,
     &    556.,  718., -171.,
     &    556.,  718., -171.,
     &    278.,  334.,  172.,
     &    278.,    0.,    0.,
     &    556.,  700., -191.,
     &    350.,  524.,  194.,
     &    278.,  127., -146.,
     &    500.,  127., -146.,
     &    500.,  718.,  445.,
     &    556.,  484.,   76.,
     &   1000.,  146.,    0.,
     &   1000.,  710.,  -19.,
     &    278.,    0.,    0.,
     &    611.,  532., -195./

      data  ((stwhb(3,i,j),i=1,3),j=192,223)/
     &    278.,    0.,    0.,
     &    333.,  750.,  604.,
     &    333.,  750.,  604.,
     &    333.,  750.,  604.,
     &    333.,  737.,  610.,
     &    333.,  678.,  604.,
     &    333.,  750.,  604.,
     &    333.,  729.,  614.,
     &    333.,  729.,  614.,
     &    278.,    0.,    0.,
     &    333.,  776.,  568.,
     &    333.,    0., -228.,
     &    278.,    0.,    0.,
     &    333.,  750.,  604.,
     &    333.,    0., -228.,
     &    333.,  750.,  604.,
     &   1000.,  333.,  227.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

      data  ((stwhb(3,i,j),i=1,3),j=224,255)/
     &    278.,    0.,    0.,
     &   1000.,  718.,    0.,
     &    278.,    0.,    0.,
     &    370.,  737.,  276.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    611.,  718.,    0.,
     &    778.,  745.,  -27.,
     &   1000.,  737.,  -19.,
     &    365.,  737.,  276.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    889.,  546.,  -14.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,  532.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,  718.,    0.,
     &    611.,  560.,  -29.,
     &    944.,  546.,  -14.,
     &    611.,  731.,  -14.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0.,
     &    278.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(4,i,j),i=1,3),j=0,31)/
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(4,i,j),i=1,3),j=32,63)/
     &    250.,    0.,    0.,
     &    333.,  672.,  -17.,
     &    713.,  705.,    0.,
     &    500.,  673.,  -16.,
     &    549.,  707.,    0.,
     &    833.,  655.,  -36.,
     &    778.,  662.,  -18.,
     &    439.,  500.,  -17.,
     &    333.,  673., -191.,
     &    333.,  673., -191.,
     &    500.,  551.,  134.,
     &    549.,  533.,    0.,
     &    250.,  105., -152.,
     &    549.,  288.,  233.,
     &    250.,   95.,  -17.,
     &    278.,  646.,  -18.,
     &    500.,  685.,  -17.,
     &    500.,  673.,    0.,
     &    500.,  686.,    0.,
     &    500.,  685.,  -17.,
     &    500.,  685.,    0.,
     &    500.,  685.,  -17.,
     &    500.,  685.,  -18.,
     &    500.,  673.,  -16.,
     &    500.,  685.,  -18.,
     &    500.,  685.,  -18.,
     &    278.,  460.,  -17.,
     &    278.,  460., -152.,
     &    549.,  522.,    0.,
     &    549.,  390.,  141.,
     &    549.,  522.,    0.,
     &    444.,  686.,  -17./

      data  ((stwhb(4,i,j),i=1,3),j=64,95)/
     &    549.,  475.,    0.,
     &    722.,  673.,    0.,
     &    667.,  673.,    0.,
     &    722.,  673.,    0.,
     &    612.,  688.,    0.,
     &    611.,  673.,    0.,
     &    763.,  673.,    0.,
     &    603.,  673.,    0.,
     &    722.,  673.,    0.,
     &    333.,  673.,    0.,
     &    631.,  689.,  -18.,
     &    722.,  673.,    0.,
     &    686.,  688.,    0.,
     &    889.,  673.,    0.,
     &    722.,  673.,   -8.,
     &    722.,  685.,  -17.,
     &    768.,  673.,    0.,
     &    741.,  685.,  -17.,
     &    556.,  673.,    0.,
     &    592.,  673.,    0.,
     &    611.,  673.,    0.,
     &    690.,  673.,    0.,
     &    439.,  500., -233.,
     &    768.,  688.,    0.,
     &    645.,  673.,    0.,
     &    795.,  684.,    0.,
     &    611.,  673.,    0.,
     &    333.,  674., -155.,
     &    863.,  478.,    0.,
     &    333.,  674., -155.,
     &    658.,  674.,    0.,
     &    500., -206., -252./

      data  ((stwhb(4,i,j),i=1,3),j=96,127)/
     &    500.,  917.,  881.,
     &    631.,  500.,  -18.,
     &    549.,  741., -223.,
     &    549.,  500., -231.,
     &    494.,  740.,  -19.,
     &    439.,  502.,  -19.,
     &    521.,  671., -224.,
     &    411.,  499., -225.,
     &    603.,  514., -202.,
     &    329.,  503.,  -17.,
     &    603.,  499., -224.,
     &    549.,  501.,    0.,
     &    549.,  740.,  -17.,
     &    576.,  500., -223.,
     &    521.,  507.,  -16.,
     &    549.,  499.,  -19.,
     &    549.,  487.,  -19.,
     &    521.,  690.,  -17.,
     &    549.,  499., -230.,
     &    603.,  500.,  -21.,
     &    439.,  500.,  -19.,
     &    576.,  507.,  -18.,
     &    713.,  583.,  -18.,
     &    686.,  500.,  -17.,
     &    493.,  766., -225.,
     &    686.,  500., -228.,
     &    494.,  756., -225.,
     &    480.,  673., -183.,
     &    200.,  673., -177.,
     &    480.,  673., -183.,
     &    549.,  307.,  203.,
     &    250.,    0.,    0./

      data  ((stwhb(4,i,j),i=1,3),j=128,159)/
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(4,i,j),i=1,3),j=160,191)/
     &    250.,    0.,    0.,
     &    620.,  685.,    0.,
     &    247.,  735.,  459.,
     &    549.,  639.,    0.,
     &    167.,  677.,  -12.,
     &    713.,  404.,  125.,
     &    500.,  687., -193.,
     &    753.,  533.,  -26.,
     &    753.,  550.,  -36.,
     &    753.,  532.,  -33.,
     &    753.,  548.,  -36.,
     &   1042.,  511.,  -15.,
     &    987.,  511.,  -15.,
     &    603.,  910.,    0.,
     &    987.,  511.,  -15.,
     &    603.,  888.,  -22.,
     &    400.,  685.,  385.,
     &    549.,  645.,    0.,
     &    411.,  737.,  459.,
     &    549.,  639.,    0.,
     &    549.,  524.,    8.,
     &    713.,  404.,  123.,
     &    494.,  746.,  -21.,
     &    460.,  473.,  113.,
     &    549.,  456.,   71.,
     &    549.,  549.,  -25.,
     &    549.,  443.,   82.,
     &    549.,  394.,  135.,
     &   1000.,   95.,  -17.,
     &    603., 1010., -120.,
     &   1000.,  276.,  220.,
     &    658.,  629.,  -16./

      data  ((stwhb(4,i,j),i=1,3),j=192,223)/
     &    823.,  658.,  -18.,
     &    686.,  740.,  -53.,
     &    795.,  734.,  -15.,
     &    987.,  573., -211.,
     &    768.,  673.,  -17.,
     &    768.,  675.,  -15.,
     &    823.,  719.,  -24.,
     &    768.,  509.,    0.,
     &    768.,  492.,  -17.,
     &    713.,  470.,    0.,
     &    713.,  470., -125.,
     &    713.,  540.,  -70.,
     &    713.,  470.,    0.,
     &    713.,  470., -125.,
     &    713.,  468.,    0.,
     &    713.,  555.,  -58.,
     &    768.,  673.,    0.,
     &    713.,  718.,  -19.,
     &    790.,  673.,  -17.,
     &    790.,  675.,  -15.,
     &    890.,  673.,  293.,
     &    823.,  751., -101.,
     &    549.,  917.,  -38.,
     &    250.,  310.,  210.,
     &    713.,  288.,    0.,
     &    603.,  454.,    0.,
     &    603.,  477.,    0.,
     &   1042.,  510.,  -20.,
     &    987.,  513.,  -15.,
     &    603.,  911.,    2.,
     &    987.,  508.,  -20.,
     &    603.,  890.,  -19./

      data  ((stwhb(4,i,j),i=1,3),j=224,255)/
     &    494.,  745.,    0.,
     &    329.,  746., -198.,
     &    790.,  670.,  -20.,
     &    790.,  675.,  -15.,
     &    786.,  673.,  293.,
     &    713.,  752., -108.,
     &    384.,  926., -293.,
     &    384.,  925.,  -85.,
     &    384.,  926., -293.,
     &    384.,  926.,  -80.,
     &    384.,  925.,  -79.,
     &    384.,  926.,  -80.,
     &    494.,  926.,  -75.,
     &    494.,  935.,  -85.,
     &    494.,  926.,  -70.,
     &    494.,  935.,  -80.,
     &    250.,    0.,    0.,
     &    329.,  746., -198.,
     &    274.,  916., -107.,
     &    686.,  922.,  -83.,
     &    686.,  975.,  -88.,
     &    686.,  921.,  -81.,
     &    384.,  926., -293.,
     &    384.,  925.,  -85.,
     &    384.,  926., -293.,
     &    384.,  926.,  -80.,
     &    384.,  925.,  -79.,
     &    384.,  926.,  -80.,
     &    494.,  926.,  -75.,
     &    494.,  935.,  -85.,
     &    494.,  926.,  -70.,
     &    250.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(5,i,j),i=1,3),j=0,31)/
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(5,i,j),i=1,3),j=32,63)/
     &    250.,    0.,    0.,
     &    333.,  676.,   -9.,
     &    408.,  676.,  431.,
     &    500.,  662.,    0.,
     &    500.,  727.,  -87.,
     &    833.,  676.,  -13.,
     &    778.,  676.,  -13.,
     &    333.,  676.,  433.,
     &    333.,  676., -177.,
     &    333.,  676., -177.,
     &    500.,  676.,  265.,
     &    564.,  506.,    0.,
     &    250.,  102., -141.,
     &    333.,  257.,  194.,
     &    250.,  100.,  -11.,
     &    278.,  676.,  -14.,
     &    500.,  676.,  -14.,
     &    500.,  676.,    0.,
     &    500.,  676.,    0.,
     &    500.,  676.,  -14.,
     &    500.,  676.,    0.,
     &    500.,  688.,  -14.,
     &    500.,  684.,  -14.,
     &    500.,  662.,   -8.,
     &    500.,  676.,  -14.,
     &    500.,  676.,  -22.,
     &    278.,  459.,  -11.,
     &    278.,  459., -141.,
     &    564.,  514.,   -8.,
     &    564.,  386.,  120.,
     &    564.,  514.,   -8.,
     &    444.,  676.,   -8./

      data  ((stwhb(5,i,j),i=1,3),j=64,95)/
     &    921.,  676.,  -14.,
     &    722.,  674.,    0.,
     &    667.,  662.,    0.,
     &    667.,  676.,  -14.,
     &    722.,  662.,    0.,
     &    611.,  662.,    0.,
     &    556.,  662.,    0.,
     &    722.,  676.,  -14.,
     &    722.,  662.,    0.,
     &    333.,  662.,    0.,
     &    389.,  662.,  -14.,
     &    722.,  662.,    0.,
     &    611.,  662.,    0.,
     &    889.,  662.,    0.,
     &    722.,  662.,  -11.,
     &    722.,  676.,  -14.,
     &    556.,  662.,    0.,
     &    722.,  676., -178.,
     &    667.,  662.,    0.,
     &    556.,  676.,  -14.,
     &    611.,  662.,    0.,
     &    722.,  662.,  -14.,
     &    722.,  662.,  -11.,
     &    944.,  662.,  -11.,
     &    722.,  662.,    0.,
     &    722.,  662.,    0.,
     &    611.,  662.,    0.,
     &    333.,  662., -156.,
     &    278.,  676.,  -14.,
     &    333.,  662., -156.,
     &    469.,  662.,  297.,
     &    500.,  -75., -125./

      data  ((stwhb(5,i,j),i=1,3),j=96,127)/
     &    333.,  676.,  433.,
     &    444.,  460.,  -10.,
     &    500.,  683.,  -10.,
     &    444.,  460.,  -10.,
     &    500.,  683.,  -10.,
     &    444.,  460.,  -10.,
     &    333.,  683.,    0.,
     &    500.,  460., -218.,
     &    500.,  683.,    0.,
     &    278.,  683.,    0.,
     &    278.,  683., -218.,
     &    500.,  683.,    0.,
     &    278.,  683.,    0.,
     &    778.,  460.,    0.,
     &    500.,  460.,    0.,
     &    500.,  460.,  -10.,
     &    500.,  460., -217.,
     &    500.,  460., -217.,
     &    333.,  460.,    0.,
     &    389.,  460.,  -10.,
     &    278.,  579.,  -10.,
     &    500.,  450.,  -10.,
     &    500.,  450.,  -14.,
     &    722.,  450.,  -14.,
     &    500.,  450.,    0.,
     &    500.,  450., -218.,
     &    444.,  450.,    0.,
     &    480.,  680., -181.,
     &    200.,  676.,  -14.,
     &    480.,  680., -181.,
     &    541.,  323.,  183.,
     &    250.,    0.,    0./

      data  ((stwhb(5,i,j),i=1,3),j=128,159)/
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(5,i,j),i=1,3),j=160,191)/
     &    250.,    0.,    0.,
     &    333.,  467., -218.,
     &    500.,  579., -138.,
     &    500.,  676.,   -8.,
     &    167.,  676.,  -14.,
     &    500.,  662.,    0.,
     &    500.,  676., -189.,
     &    500.,  676., -148.,
     &    500.,  602.,   58.,
     &    180.,  676.,  431.,
     &    444.,  676.,  433.,
     &    500.,  416.,   33.,
     &    333.,  416.,   33.,
     &    333.,  416.,   33.,
     &    556.,  683.,    0.,
     &    556.,  683.,    0.,
     &    250.,    0.,    0.,
     &    500.,  250.,  201.,
     &    500.,  676., -149.,
     &    500.,  676., -153.,
     &    250.,  310.,  199.,
     &    250.,    0.,    0.,
     &    453.,  662., -154.,
     &    350.,  466.,  196.,
     &    333.,  102., -141.,
     &    444.,  102., -141.,
     &    444.,  676.,  433.,
     &    500.,  416.,   33.,
     &   1000.,  100.,  -11.,
     &   1000.,  706.,  -19.,
     &    250.,    0.,    0.,
     &    444.,  466., -218./

      data  ((stwhb(5,i,j),i=1,3),j=192,223)/
     &    250.,    0.,    0.,
     &    333.,  678.,  507.,
     &    333.,  678.,  507.,
     &    333.,  674.,  507.,
     &    333.,  638.,  532.,
     &    333.,  601.,  547.,
     &    333.,  664.,  507.,
     &    333.,  623.,  523.,
     &    333.,  623.,  523.,
     &    250.,    0.,    0.,
     &    333.,  711.,  512.,
     &    333.,    0., -225.,
     &    250.,    0.,    0.,
     &    333.,  678.,  507.,
     &    333.,    0., -165.,
     &    333.,  674.,  507.,
     &   1000.,  250.,  201.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(5,i,j),i=1,3),j=224,255)/
     &    250.,    0.,    0.,
     &    889.,  662.,    0.,
     &    250.,    0.,    0.,
     &    276.,  676.,  394.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    611.,  662.,    0.,
     &    722.,  734.,  -80.,
     &    889.,  668.,   -6.,
     &    310.,  676.,  394.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    667.,  460.,  -10.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    278.,  460.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    278.,  683.,    0.,
     &    500.,  551., -112.,
     &    722.,  460.,  -10.,
     &    500.,  683.,   -9.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(6,i,j),i=1,3),j=0,31)/
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(6,i,j),i=1,3),j=32,63)/
     &    250.,    0.,    0.,
     &    333.,  667.,   -9.,
     &    420.,  666.,  421.,
     &    500.,  676.,    0.,
     &    500.,  731.,  -89.,
     &    833.,  676.,  -13.,
     &    778.,  666.,  -18.,
     &    333.,  666.,  436.,
     &    333.,  669., -181.,
     &    333.,  669., -180.,
     &    500.,  666.,  255.,
     &    675.,  506.,    0.,
     &    250.,  101., -129.,
     &    333.,  255.,  192.,
     &    250.,  100.,  -11.,
     &    278.,  666.,  -18.,
     &    500.,  676.,   -7.,
     &    500.,  676.,    0.,
     &    500.,  676.,    0.,
     &    500.,  676.,   -7.,
     &    500.,  676.,    0.,
     &    500.,  666.,   -7.,
     &    500.,  686.,   -7.,
     &    500.,  666.,   -8.,
     &    500.,  676.,   -7.,
     &    500.,  676.,  -17.,
     &    333.,  441.,  -11.,
     &    333.,  441., -129.,
     &    675.,  514.,   -8.,
     &    675.,  386.,  120.,
     &    675.,  514.,   -8.,
     &    500.,  664.,  -12./

      data  ((stwhb(6,i,j),i=1,3),j=64,95)/
     &    920.,  666.,  -18.,
     &    611.,  668.,    0.,
     &    611.,  653.,    0.,
     &    667.,  666.,  -18.,
     &    722.,  653.,    0.,
     &    611.,  653.,    0.,
     &    611.,  653.,    0.,
     &    722.,  666.,  -18.,
     &    722.,  653.,    0.,
     &    333.,  653.,    0.,
     &    444.,  653.,  -18.,
     &    667.,  653.,    0.,
     &    556.,  653.,    0.,
     &    833.,  653.,    0.,
     &    667.,  653.,  -15.,
     &    722.,  666.,  -18.,
     &    611.,  653.,    0.,
     &    722.,  666., -182.,
     &    611.,  653.,    0.,
     &    500.,  667.,  -18.,
     &    556.,  653.,    0.,
     &    722.,  653.,  -18.,
     &    611.,  653.,  -18.,
     &    833.,  653.,  -18.,
     &    611.,  653.,    0.,
     &    556.,  653.,    0.,
     &    556.,  653.,    0.,
     &    389.,  663., -153.,
     &    278.,  666.,  -18.,
     &    389.,  663., -153.,
     &    422.,  666.,  301.,
     &    500.,  -75., -125./

      data  ((stwhb(6,i,j),i=1,3),j=96,127)/
     &    333.,  666.,  436.,
     &    500.,  441.,  -11.,
     &    500.,  683.,  -11.,
     &    444.,  441.,  -11.,
     &    500.,  683.,  -13.,
     &    444.,  441.,  -11.,
     &    278.,  678., -207.,
     &    500.,  441., -206.,
     &    500.,  683.,   -9.,
     &    278.,  654.,    0.,
     &    278.,  654., -207.,
     &    444.,  683.,  -11.,
     &    278.,  683.,  -11.,
     &    722.,  441.,   -9.,
     &    500.,  441.,   -9.,
     &    500.,  441.,  -11.,
     &    500.,  441., -205.,
     &    500.,  441., -209.,
     &    389.,  441.,    0.,
     &    389.,  442.,  -13.,
     &    278.,  546.,  -11.,
     &    500.,  441.,  -11.,
     &    444.,  441.,  -18.,
     &    667.,  441.,  -18.,
     &    444.,  441.,  -11.,
     &    444.,  441., -206.,
     &    389.,  428.,  -81.,
     &    400.,  687., -177.,
     &    275.,  666.,  -18.,
     &    400.,  687., -177.,
     &    541.,  323.,  183.,
     &    250.,    0.,    0./

      data  ((stwhb(6,i,j),i=1,3),j=128,159)/
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(6,i,j),i=1,3),j=160,191)/
     &    250.,    0.,    0.,
     &    389.,  473., -205.,
     &    500.,  560., -143.,
     &    500.,  670.,   -6.,
     &    167.,  676.,  -10.,
     &    500.,  653.,    0.,
     &    500.,  682., -182.,
     &    500.,  666., -162.,
     &    500.,  697.,   53.,
     &    214.,  666.,  421.,
     &    556.,  666.,  436.,
     &    500.,  403.,   37.,
     &    333.,  403.,   37.,
     &    333.,  403.,   37.,
     &    500.,  681., -207.,
     &    500.,  682., -204.,
     &    250.,    0.,    0.,
     &    500.,  243.,  197.,
     &    500.,  666., -159.,
     &    500.,  666., -143.,
     &    250.,  310.,  199.,
     &    250.,    0.,    0.,
     &    523.,  653., -123.,
     &    350.,  461.,  191.,
     &    333.,  101., -129.,
     &    556.,  101., -129.,
     &    556.,  666.,  436.,
     &    500.,  403.,   37.,
     &    889.,  100.,  -11.,
     &   1000.,  706.,  -19.,
     &    250.,    0.,    0.,
     &    500.,  471., -205./

      data  ((stwhb(6,i,j),i=1,3),j=192,223)/
     &    250.,    0.,    0.,
     &    333.,  664.,  492.,
     &    333.,  664.,  494.,
     &    333.,  661.,  492.,
     &    333.,  624.,  517.,
     &    333.,  583.,  532.,
     &    333.,  650.,  492.,
     &    333.,  606.,  508.,
     &    333.,  606.,  508.,
     &    250.,    0.,    0.,
     &    333.,  691.,  492.,
     &    333.,    0., -227.,
     &    250.,    0.,    0.,
     &    333.,  664.,  494.,
     &    333.,   40., -169.,
     &    333.,  661.,  492.,
     &    889.,  243.,  197.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(6,i,j),i=1,3),j=224,255)/
     &    250.,    0.,    0.,
     &    889.,  653.,    0.,
     &    250.,    0.,    0.,
     &    276.,  676.,  406.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    556.,  653.,    0.,
     &    722.,  722., -105.,
     &    944.,  666.,   -8.,
     &    310.,  676.,  406.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    667.,  441.,  -11.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    278.,  441.,  -11.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    278.,  683.,  -11.,
     &    500.,  554., -135.,
     &    667.,  441.,  -12.,
     &    500.,  679., -207.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(7,i,j),i=1,3),j=0,31)/
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(7,i,j),i=1,3),j=32,63)/
     &    250.,    0.,    0.,
     &    333.,  691.,  -13.,
     &    555.,  691.,  404.,
     &    500.,  700.,    0.,
     &    500.,  750.,  -99.,
     &   1000.,  692.,  -14.,
     &    833.,  691.,  -16.,
     &    333.,  691.,  356.,
     &    333.,  694., -168.,
     &    333.,  694., -168.,
     &    500.,  691.,  255.,
     &    570.,  506.,    0.,
     &    250.,  155., -180.,
     &    333.,  287.,  171.,
     &    250.,  156.,  -13.,
     &    278.,  691.,  -19.,
     &    500.,  688.,  -13.,
     &    500.,  688.,    0.,
     &    500.,  688.,    0.,
     &    500.,  688.,  -14.,
     &    500.,  688.,    0.,
     &    500.,  676.,   -8.,
     &    500.,  688.,  -13.,
     &    500.,  676.,    0.,
     &    500.,  688.,  -13.,
     &    500.,  688.,  -13.,
     &    333.,  472.,  -13.,
     &    333.,  472., -180.,
     &    570.,  514.,   -8.,
     &    570.,  399.,  107.,
     &    570.,  514.,   -8.,
     &    500.,  689.,  -13./

      data  ((stwhb(7,i,j),i=1,3),j=64,95)/
     &    930.,  691.,  -19.,
     &    722.,  690.,    0.,
     &    667.,  676.,    0.,
     &    722.,  691.,  -19.,
     &    722.,  676.,    0.,
     &    667.,  676.,    0.,
     &    611.,  676.,    0.,
     &    778.,  691.,  -19.,
     &    778.,  676.,    0.,
     &    389.,  676.,    0.,
     &    500.,  676.,  -96.,
     &    778.,  676.,    0.,
     &    667.,  676.,    0.,
     &    944.,  676.,    0.,
     &    722.,  676.,  -18.,
     &    778.,  691.,  -19.,
     &    611.,  676.,    0.,
     &    778.,  691., -176.,
     &    722.,  676.,    0.,
     &    556.,  692.,  -19.,
     &    667.,  676.,    0.,
     &    722.,  676.,  -19.,
     &    722.,  676.,  -18.,
     &   1000.,  676.,  -15.,
     &    722.,  676.,    0.,
     &    722.,  676.,    0.,
     &    667.,  676.,    0.,
     &    333.,  678., -149.,
     &    278.,  691.,  -19.,
     &    333.,  678., -149.,
     &    581.,  676.,  311.,
     &    500.,  -75., -125./

      data  ((stwhb(7,i,j),i=1,3),j=96,127)/
     &    333.,  691.,  356.,
     &    500.,  473.,  -14.,
     &    556.,  676.,  -14.,
     &    444.,  473.,  -14.,
     &    556.,  676.,  -14.,
     &    444.,  473.,  -14.,
     &    333.,  691.,    0.,
     &    500.,  473., -206.,
     &    556.,  676.,    0.,
     &    278.,  691.,    0.,
     &    333.,  691., -203.,
     &    556.,  676.,    0.,
     &    278.,  676.,    0.,
     &    833.,  473.,    0.,
     &    556.,  473.,    0.,
     &    500.,  473.,  -14.,
     &    556.,  473., -205.,
     &    556.,  473., -205.,
     &    444.,  473.,    0.,
     &    389.,  473.,  -14.,
     &    333.,  630.,  -12.,
     &    556.,  461.,  -14.,
     &    500.,  461.,  -14.,
     &    722.,  461.,  -14.,
     &    500.,  461.,    0.,
     &    500.,  461., -205.,
     &    444.,  461.,    0.,
     &    394.,  698., -175.,
     &    220.,  691.,  -19.,
     &    394.,  698., -175.,
     &    520.,  333.,  173.,
     &    250.,    0.,    0./

      data  ((stwhb(7,i,j),i=1,3),j=128,159)/
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(7,i,j),i=1,3),j=160,191)/
     &    250.,    0.,    0.,
     &    333.,  501., -203.,
     &    500.,  588., -140.,
     &    500.,  684.,  -14.,
     &    167.,  688.,  -12.,
     &    500.,  676.,    0.,
     &    500.,  706., -155.,
     &    500.,  691., -132.,
     &    500.,  613.,   61.,
     &    278.,  691.,  404.,
     &    500.,  691.,  356.,
     &    500.,  415.,   36.,
     &    333.,  415.,   36.,
     &    333.,  415.,   36.,
     &    556.,  691.,    0.,
     &    556.,  691.,    0.,
     &    250.,    0.,    0.,
     &    500.,  271.,  181.,
     &    500.,  691., -134.,
     &    500.,  691., -132.,
     &    250.,  417.,  248.,
     &    250.,    0.,    0.,
     &    540.,  676., -186.,
     &    350.,  478.,  198.,
     &    333.,  155., -180.,
     &    500.,  155., -180.,
     &    500.,  691.,  356.,
     &    500.,  415.,   36.,
     &   1000.,  156.,  -13.,
     &   1000.,  706.,  -29.,
     &    250.,    0.,    0.,
     &    500.,  501., -201./

      data  ((stwhb(7,i,j),i=1,3),j=192,223)/
     &    250.,    0.,    0.,
     &    333.,  713.,  528.,
     &    333.,  713.,  528.,
     &    333.,  704.,  528.,
     &    333.,  674.,  547.,
     &    333.,  637.,  565.,
     &    333.,  691.,  528.,
     &    333.,  667.,  537.,
     &    333.,  667.,  537.,
     &    250.,    0.,    0.,
     &    333.,  740.,  527.,
     &    333.,    0., -218.,
     &    250.,    0.,    0.,
     &    333.,  713.,  528.,
     &    333.,   44., -173.,
     &    333.,  704.,  528.,
     &   1000.,  271.,  181.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(7,i,j),i=1,3),j=224,255)/
     &    250.,    0.,    0.,
     &   1000.,  676.,    0.,
     &    250.,    0.,    0.,
     &    300.,  688.,  397.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    667.,  676.,    0.,
     &    778.,  737.,  -74.,
     &   1000.,  684.,   -5.,
     &    330.,  688.,  397.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    722.,  473.,  -14.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    278.,  461.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    278.,  676.,    0.,
     &    500.,  549.,  -92.,
     &    722.,  473.,  -14.,
     &    556.,  691.,  -12.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(8,i,j),i=1,3),j=0,31)/
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(8,i,j),i=1,3),j=32,63)/
     &    250.,    0.,    0.,
     &    333.,  684.,  -13.,
     &    555.,  685.,  398.,
     &    500.,  700.,    0.,
     &    500.,  733., -100.,
     &    833.,  692.,  -10.,
     &    778.,  682.,  -19.,
     &    333.,  685.,  369.,
     &    333.,  685., -179.,
     &    333.,  685., -179.,
     &    500.,  685.,  249.,
     &    570.,  506.,    0.,
     &    250.,  134., -182.,
     &    333.,  282.,  166.,
     &    250.,  135.,  -13.,
     &    278.,  685.,  -18.,
     &    500.,  683.,  -14.,
     &    500.,  683.,    0.,
     &    500.,  683.,    0.,
     &    500.,  683.,  -13.,
     &    500.,  683.,    0.,
     &    500.,  669.,  -13.,
     &    500.,  679.,  -15.,
     &    500.,  669.,    0.,
     &    500.,  683.,  -13.,
     &    500.,  683.,  -10.,
     &    333.,  459.,  -13.,
     &    333.,  459., -183.,
     &    570.,  514.,   -8.,
     &    570.,  399.,  107.,
     &    570.,  514.,   -8.,
     &    500.,  684.,  -13./

      data  ((stwhb(8,i,j),i=1,3),j=64,95)/
     &    832.,  685.,  -18.,
     &    667.,  683.,    0.,
     &    667.,  669.,    0.,
     &    667.,  685.,  -18.,
     &    722.,  669.,    0.,
     &    667.,  669.,    0.,
     &    667.,  669.,    0.,
     &    722.,  685.,  -18.,
     &    778.,  669.,    0.,
     &    389.,  669.,    0.,
     &    500.,  669.,  -99.,
     &    667.,  669.,    0.,
     &    611.,  662.,    0.,
     &    889.,  669.,  -12.,
     &    722.,  669.,  -15.,
     &    722.,  685.,  -18.,
     &    611.,  669.,    0.,
     &    722.,  685., -208.,
     &    667.,  669.,    0.,
     &    556.,  685.,  -18.,
     &    611.,  669.,    0.,
     &    722.,  669.,  -18.,
     &    667.,  669.,  -18.,
     &    889.,  669.,  -18.,
     &    667.,  669.,    0.,
     &    611.,  669.,    0.,
     &    611.,  669.,    0.,
     &    333.,  674., -159.,
     &    278.,  685.,  -18.,
     &    333.,  674., -157.,
     &    570.,  669.,  304.,
     &    500.,  -75., -125./

      data  ((stwhb(8,i,j),i=1,3),j=96,127)/
     &    333.,  685.,  369.,
     &    500.,  462.,  -14.,
     &    500.,  699.,  -13.,
     &    444.,  462.,  -13.,
     &    500.,  699.,  -13.,
     &    444.,  462.,  -13.,
     &    333.,  698., -205.,
     &    500.,  462., -203.,
     &    556.,  699.,   -9.,
     &    278.,  684.,   -9.,
     &    278.,  684., -207.,
     &    500.,  699.,   -8.,
     &    278.,  699.,   -9.,
     &    778.,  462.,   -9.,
     &    556.,  462.,   -9.,
     &    500.,  462.,  -13.,
     &    500.,  462., -205.,
     &    500.,  462., -205.,
     &    389.,  462.,    0.,
     &    389.,  462.,  -13.,
     &    278.,  594.,   -9.,
     &    556.,  462.,   -9.,
     &    444.,  462.,  -13.,
     &    667.,  462.,  -13.,
     &    500.,  462.,  -13.,
     &    444.,  462., -205.,
     &    389.,  449.,  -78.,
     &    348.,  686., -187.,
     &    220.,  685.,  -18.,
     &    348.,  686., -187.,
     &    571.,  333.,  173.,
     &    250.,    0.,    0./

      data  ((stwhb(8,i,j),i=1,3),j=128,159)/
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(8,i,j),i=1,3),j=160,191)/
     &    250.,    0.,    0.,
     &    389.,  492., -205.,
     &    500.,  576., -143.,
     &    500.,  683.,  -12.,
     &    167.,  683.,  -14.,
     &    500.,  669.,    0.,
     &    500.,  707., -156.,
     &    500.,  685., -143.,
     &    500.,  586.,   34.,
     &    278.,  685.,  398.,
     &    500.,  685.,  369.,
     &    500.,  415.,   32.,
     &    333.,  415.,   32.,
     &    333.,  415.,   32.,
     &    556.,  703., -205.,
     &    556.,  704., -205.,
     &    250.,    0.,    0.,
     &    500.,  269.,  178.,
     &    500.,  685., -145.,
     &    500.,  685., -139.,
     &    250.,  405.,  257.,
     &    250.,    0.,    0.,
     &    500.,  669., -193.,
     &    350.,  525.,  175.,
     &    333.,  134., -182.,
     &    500.,  134., -182.,
     &    500.,  685.,  369.,
     &    500.,  415.,   32.,
     &   1000.,  135.,  -13.,
     &   1000.,  706.,  -29.,
     &    250.,    0.,    0.,
     &    500.,  492., -205./

      data  ((stwhb(8,i,j),i=1,3),j=192,223)/
     &    250.,    0.,    0.,
     &    333.,  697.,  516.,
     &    333.,  697.,  516.,
     &    333.,  690.,  516.,
     &    333.,  655.,  536.,
     &    333.,  623.,  553.,
     &    333.,  678.,  516.,
     &    333.,  655.,  525.,
     &    333.,  655.,  525.,
     &    250.,    0.,    0.,
     &    333.,  729.,  516.,
     &    333.,    5., -218.,
     &    250.,    0.,    0.,
     &    333.,  697.,  516.,
     &    333.,   44., -173.,
     &    333.,  690.,  516.,
     &   1000.,  269.,  178.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(8,i,j),i=1,3),j=224,255)/
     &    250.,    0.,    0.,
     &    944.,  669.,    0.,
     &    250.,    0.,    0.,
     &    266.,  685.,  399.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    611.,  669.,    0.,
     &    722.,  764., -125.,
     &    944.,  677.,   -8.,
     &    300.,  685.,  400.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    722.,  462.,  -13.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    278.,  462.,   -9.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    278.,  699.,   -9.,
     &    500.,  560., -119.,
     &    722.,  462.,  -13.,
     &    500.,  705., -200.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(9,i,j),i=1,3),j=0,31)/
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

      data  ((stwhb(9,i,j),i=1,3),j=32,63)/
     &    450.,    0.,    0.,
     &    450.,  572.,  -15.,
     &    450.,  562.,  328.,
     &    450.,  639.,  -32.,
     &    450.,  662., -126.,
     &    450.,  622.,  -15.,
     &    450.,  543.,  -15.,
     &    450.,  562.,  328.,
     &    450.,  622., -108.,
     &    450.,  622., -108.,
     &    450.,  607.,  257.,
     &    450.,  470.,   44.,
     &    450.,  122., -112.,
     &    450.,  285.,  231.,
     &    450.,  109.,  -15.,
     &    450.,  629.,  -80.,
     &    450.,  622.,  -15.,
     &    450.,  622.,    0.,
     &    450.,  622.,    0.,
     &    450.,  622.,  -15.,
     &    450.,  622.,    0.,
     &    450.,  607.,  -15.,
     &    450.,  622.,  -15.,
     &    450.,  607.,    0.,
     &    450.,  622.,  -15.,
     &    450.,  622.,  -15.,
     &    450.,  385.,  -15.,
     &    450.,  385., -112.,
     &    450.,  472.,   42.,
     &    450.,  376.,  138.,
     &    450.,  472.,   42.,
     &    450.,  572.,  -15./

      data  ((stwhb(9,i,j),i=1,3),j=64,95)/
     &    450.,  622.,  -15.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  580.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  580.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,  -13.,
     &    450.,  580.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  580., -138.,
     &    450.,  562.,    0.,
     &    450.,  580.,  -20.,
     &    450.,  562.,    0.,
     &    450.,  562.,  -18.,
     &    450.,  562.,  -13.,
     &    450.,  562.,  -13.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  622., -108.,
     &    450.,  629.,  -80.,
     &    450.,  622., -108.,
     &    450.,  622.,  354.,
     &    450.,  -75., -125./

      data  ((stwhb(9,i,j),i=1,3),j=96,127)/
     &    450.,  562.,  328.,
     &    450.,  441.,  -15.,
     &    450.,  629.,  -15.,
     &    450.,  441.,  -15.,
     &    450.,  629.,  -15.,
     &    450.,  441.,  -15.,
     &    450.,  629.,    0.,
     &    450.,  441., -157.,
     &    450.,  629.,    0.,
     &    450.,  657.,    0.,
     &    450.,  657., -157.,
     &    450.,  629.,    0.,
     &    450.,  629.,    0.,
     &    450.,  441.,    0.,
     &    450.,  441.,    0.,
     &    450.,  441.,  -15.,
     &    450.,  441., -157.,
     &    450.,  441., -157.,
     &    450.,  441.,    0.,
     &    450.,  441.,  -15.,
     &    450.,  561.,  -15.,
     &    450.,  426.,  -15.,
     &    450.,  426.,  -10.,
     &    450.,  426.,  -10.,
     &    450.,  426.,    0.,
     &    450.,  426., -157.,
     &    450.,  426.,    0.,
     &    450.,  622., -108.,
     &    450.,  750., -250.,
     &    450.,  622., -108.,
     &    450.,  320.,  197.,
     &    450.,    0.,    0./

      data  ((stwhb(9,i,j),i=1,3),j=128,159)/
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

      data  ((stwhb(9,i,j),i=1,3),j=160,191)/
     &    450.,    0.,    0.,
     &    450.,  430., -157.,
     &    450.,  614.,  -49.,
     &    450.,  611.,  -21.,
     &    450.,  665.,  -57.,
     &    450.,  562.,    0.,
     &    450.,  622., -143.,
     &    450.,  580.,  -78.,
     &    450.,  506.,   58.,
     &    450.,  562.,  328.,
     &    450.,  562.,  328.,
     &    450.,  446.,   70.,
     &    450.,  446.,   70.,
     &    450.,  446.,   70.,
     &    450.,  629.,    0.,
     &    450.,  629.,    0.,
     &    450.,    0.,    0.,
     &    450.,  285.,  231.,
     &    450.,  580.,  -78.,
     &    450.,  580.,  -78.,
     &    450.,  327.,  189.,
     &    450.,    0.,    0.,
     &    450.,  562.,  -78.,
     &    450.,  383.,  130.,
     &    450.,  100., -134.,
     &    450.,  100., -134.,
     &    450.,  562.,  328.,
     &    450.,  446.,   70.,
     &    450.,  111.,  -15.,
     &    450.,  622.,  -15.,
     &    450.,    0.,    0.,
     &    450.,  430., -157./

      data  ((stwhb(9,i,j),i=1,3),j=192,223)/
     &    450.,    0.,    0.,
     &    450.,  672.,  497.,
     &    450.,  672.,  497.,
     &    450.,  654.,  477.,
     &    450.,  606.,  489.,
     &    450.,  565.,  525.,
     &    450.,  609.,  501.,
     &    450.,  580.,  477.,
     &    450.,  595.,  492.,
     &    450.,    0.,    0.,
     &    450.,  627.,  463.,
     &    450.,   10., -151.,
     &    450.,    0.,    0.,
     &    450.,  672.,  497.,
     &    450.,    0., -151.,
     &    450.,  669.,  492.,
     &    450.,  285.,  231.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

      data  ((stwhb(9,i,j),i=1,3),j=224,255)/
     &    450.,    0.,    0.,
     &    450.,  562.,    0.,
     &    450.,    0.,    0.,
     &    450.,  580.,  249.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  562.,    0.,
     &    450.,  629.,  -80.,
     &    450.,  562.,    0.,
     &    450.,  580.,  249.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  441.,  -15.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  426.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  629.,    0.,
     &    450.,  506.,  -80.,
     &    450.,  441.,  -15.,
     &    450.,  629.,  -15.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(10,i,j),i=1,3),j=0,31)/
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

      data  ((stwhb(10,i,j),i=1,3),j=32,63)/
     &    450.,    0.,    0.,
     &    450.,  572.,  -15.,
     &    450.,  562.,  328.,
     &    450.,  639.,  -32.,
     &    450.,  662., -126.,
     &    450.,  622.,  -15.,
     &    450.,  543.,  -15.,
     &    450.,  562.,  328.,
     &    450.,  622., -108.,
     &    450.,  622., -108.,
     &    450.,  607.,  257.,
     &    450.,  470.,   44.,
     &    450.,  122., -112.,
     &    450.,  285.,  231.,
     &    450.,  109.,  -15.,
     &    450.,  629.,  -80.,
     &    450.,  622.,  -15.,
     &    450.,  622.,    0.,
     &    450.,  622.,    0.,
     &    450.,  622.,  -15.,
     &    450.,  622.,    0.,
     &    450.,  607.,  -15.,
     &    450.,  622.,  -15.,
     &    450.,  607.,    0.,
     &    450.,  622.,  -15.,
     &    450.,  622.,  -15.,
     &    450.,  385.,  -15.,
     &    450.,  385., -112.,
     &    450.,  472.,   42.,
     &    450.,  376.,  138.,
     &    450.,  472.,   42.,
     &    450.,  572.,  -15./

      data  ((stwhb(10,i,j),i=1,3),j=64,95)/
     &    450.,  622.,  -15.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  580.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  580.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,  -13.,
     &    450.,  580.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  580., -138.,
     &    450.,  562.,    0.,
     &    450.,  580.,  -20.,
     &    450.,  562.,    0.,
     &    450.,  562.,  -18.,
     &    450.,  562.,  -13.,
     &    450.,  562.,  -13.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  622., -108.,
     &    450.,  629.,  -80.,
     &    450.,  622., -108.,
     &    450.,  622.,  354.,
     &    450.,  -75., -125./

      data  ((stwhb(10,i,j),i=1,3),j=96,127)/
     &    450.,  562.,  328.,
     &    450.,  441.,  -15.,
     &    450.,  629.,  -15.,
     &    450.,  441.,  -15.,
     &    450.,  629.,  -15.,
     &    450.,  441.,  -15.,
     &    450.,  629.,    0.,
     &    450.,  441., -157.,
     &    450.,  629.,    0.,
     &    450.,  657.,    0.,
     &    450.,  657., -157.,
     &    450.,  629.,    0.,
     &    450.,  629.,    0.,
     &    450.,  441.,    0.,
     &    450.,  441.,    0.,
     &    450.,  441.,  -15.,
     &    450.,  441., -157.,
     &    450.,  441., -157.,
     &    450.,  441.,    0.,
     &    450.,  441.,  -15.,
     &    450.,  561.,  -15.,
     &    450.,  426.,  -15.,
     &    450.,  426.,  -10.,
     &    450.,  426.,  -10.,
     &    450.,  426.,    0.,
     &    450.,  426., -157.,
     &    450.,  426.,    0.,
     &    450.,  622., -108.,
     &    450.,  750., -250.,
     &    450.,  622., -108.,
     &    450.,  320.,  197.,
     &    450.,    0.,    0./

      data  ((stwhb(10,i,j),i=1,3),j=128,159)/
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

      data  ((stwhb(10,i,j),i=1,3),j=160,191)/
     &    450.,    0.,    0.,
     &    450.,  430., -157.,
     &    450.,  614.,  -49.,
     &    450.,  611.,  -21.,
     &    450.,  665.,  -57.,
     &    450.,  562.,    0.,
     &    450.,  622., -143.,
     &    450.,  580.,  -78.,
     &    450.,  506.,   58.,
     &    450.,  562.,  328.,
     &    450.,  562.,  328.,
     &    450.,  446.,   70.,
     &    450.,  446.,   70.,
     &    450.,  446.,   70.,
     &    450.,  629.,    0.,
     &    450.,  629.,    0.,
     &    450.,    0.,    0.,
     &    450.,  285.,  231.,
     &    450.,  580.,  -78.,
     &    450.,  580.,  -78.,
     &    450.,  327.,  189.,
     &    450.,    0.,    0.,
     &    450.,  562.,  -78.,
     &    450.,  383.,  130.,
     &    450.,  100., -134.,
     &    450.,  100., -134.,
     &    450.,  562.,  328.,
     &    450.,  446.,   70.,
     &    450.,  111.,  -15.,
     &    450.,  622.,  -15.,
     &    450.,    0.,    0.,
     &    450.,  430., -157./

      data  ((stwhb(10,i,j),i=1,3),j=192,223)/
     &    450.,    0.,    0.,
     &    450.,  672.,  497.,
     &    450.,  672.,  497.,
     &    450.,  654.,  477.,
     &    450.,  606.,  489.,
     &    450.,  565.,  525.,
     &    450.,  609.,  501.,
     &    450.,  580.,  477.,
     &    450.,  595.,  492.,
     &    450.,    0.,    0.,
     &    450.,  627.,  463.,
     &    450.,   10., -151.,
     &    450.,    0.,    0.,
     &    450.,  672.,  497.,
     &    450.,    0., -151.,
     &    450.,  669.,  492.,
     &    450.,  285.,  231.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

      data  ((stwhb(10,i,j),i=1,3),j=224,255)/
     &    450.,    0.,    0.,
     &    450.,  562.,    0.,
     &    450.,    0.,    0.,
     &    450.,  580.,  249.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  562.,    0.,
     &    450.,  629.,  -80.,
     &    450.,  562.,    0.,
     &    450.,  580.,  249.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  441.,  -15.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  426.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  629.,    0.,
     &    450.,  506.,  -80.,
     &    450.,  441.,  -15.,
     &    450.,  629.,  -15.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(11,i,j),i=1,3),j=0,31)/
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

      data  ((stwhb(11,i,j),i=1,3),j=32,63)/
     &    450.,    0.,    0.,
     &    450.,  572.,  -15.,
     &    450.,  562.,  277.,
     &    450.,  651.,  -45.,
     &    450.,  666., -126.,
     &    450.,  616.,  -15.,
     &    450.,  543.,  -15.,
     &    450.,  562.,  277.,
     &    450.,  616., -102.,
     &    450.,  616., -102.,
     &    450.,  601.,  219.,
     &    450.,  469.,   30.,
     &    450.,  174., -111.,
     &    450.,  313.,  203.,
     &    450.,  171.,  -15.,
     &    450.,  626.,  -77.,
     &    450.,  616.,  -15.,
     &    450.,  616.,    0.,
     &    450.,  616.,    0.,
     &    450.,  616.,  -15.,
     &    450.,  616.,    0.,
     &    450.,  601.,  -15.,
     &    450.,  616.,  -15.,
     &    450.,  601.,    0.,
     &    450.,  616.,  -15.,
     &    450.,  616.,  -15.,
     &    450.,  425.,  -15.,
     &    450.,  425., -110.,
     &    450.,  501.,   15.,
     &    450.,  389.,  109.,
     &    450.,  501.,   15.,
     &    450.,  580.,  -14./

      data  ((stwhb(11,i,j),i=1,3),j=64,95)/
     &    450.,  616.,  -15.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  580.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  580.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,  -12.,
     &    450.,  580.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  580., -138.,
     &    450.,  562.,    0.,
     &    450.,  580.,  -22.,
     &    450.,  562.,    0.,
     &    450.,  562.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  616., -102.,
     &    450.,  626.,  -77.,
     &    450.,  616., -102.,
     &    450.,  616.,  250.,
     &    450.,  -75., -125./

      data  ((stwhb(11,i,j),i=1,3),j=96,127)/
     &    450.,  562.,  277.,
     &    450.,  454.,  -15.,
     &    450.,  626.,  -15.,
     &    450.,  459.,  -15.,
     &    450.,  626.,  -15.,
     &    450.,  454.,  -15.,
     &    450.,  626.,    0.,
     &    450.,  454., -146.,
     &    450.,  626.,    0.,
     &    450.,  658.,    0.,
     &    450.,  658., -146.,
     &    450.,  626.,    0.,
     &    450.,  626.,    0.,
     &    450.,  454.,    0.,
     &    450.,  454.,    0.,
     &    450.,  454.,  -15.,
     &    450.,  454., -142.,
     &    450.,  454., -142.,
     &    450.,  454.,    0.,
     &    450.,  459.,  -17.,
     &    450.,  562.,  -15.,
     &    450.,  439.,  -15.,
     &    450.,  439.,    0.,
     &    450.,  439.,    0.,
     &    450.,  439.,    0.,
     &    450.,  439., -142.,
     &    450.,  439.,    0.,
     &    450.,  616., -102.,
     &    450.,  750., -250.,
     &    450.,  616., -102.,
     &    450.,  347.,  144.,
     &    450.,    0.,    0./

      data  ((stwhb(11,i,j),i=1,3),j=128,159)/
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

      data  ((stwhb(11,i,j),i=1,3),j=160,191)/
     &    450.,    0.,    0.,
     &    450.,  449., -146.,
     &    450.,  614.,  -49.,
     &    450.,  611.,  -28.,
     &    450.,  661.,  -60.,
     &    450.,  562.,    0.,
     &    450.,  616., -131.,
     &    450.,  580.,  -70.,
     &    450.,  517.,   49.,
     &    450.,  562.,  277.,
     &    450.,  562.,  277.,
     &    450.,  446.,   70.,
     &    450.,  446.,   70.,
     &    450.,  446.,   70.,
     &    450.,  626.,    0.,
     &    450.,  626.,    0.,
     &    450.,    0.,    0.,
     &    450.,  313.,  203.,
     &    450.,  580.,  -70.,
     &    450.,  580.,  -70.,
     &    450.,  342.,  156.,
     &    450.,    0.,    0.,
     &    450.,  580.,  -70.,
     &    450.,  430.,  132.,
     &    450.,  135., -150.,
     &    450.,  135., -150.,
     &    450.,  562.,  277.,
     &    450.,  446.,   70.,
     &    450.,  116.,  -15.,
     &    450.,  616.,  -15.,
     &    450.,    0.,    0.,
     &    450.,  449., -146./

      data  ((stwhb(11,i,j),i=1,3),j=192,223)/
     &    450.,    0.,    0.,
     &    450.,  661.,  508.,
     &    450.,  661.,  508.,
     &    450.,  657.,  483.,
     &    450.,  636.,  493.,
     &    450.,  585.,  505.,
     &    450.,  631.,  468.,
     &    450.,  625.,  485.,
     &    450.,  625.,  485.,
     &    450.,    0.,    0.,
     &    450.,  678.,  481.,
     &    450.,    0., -206.,
     &    450.,    0.,    0.,
     &    450.,  661.,  488.,
     &    450.,    0., -199.,
     &    450.,  667.,  493.,
     &    450.,  313.,  203.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

      data  ((stwhb(11,i,j),i=1,3),j=224,255)/
     &    450.,    0.,    0.,
     &    450.,  562.,    0.,
     &    450.,    0.,    0.,
     &    450.,  580.,  196.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  562.,    0.,
     &    450.,  584.,  -22.,
     &    450.,  562.,    0.,
     &    450.,  580.,  196.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  454.,  -15.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  439.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  626.,    0.,
     &    450.,  463.,  -24.,
     &    450.,  454.,  -15.,
     &    450.,  626.,  -15.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(12,i,j),i=1,3),j=0,31)/
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

      data  ((stwhb(12,i,j),i=1,3),j=32,63)/
     &    450.,    0.,    0.,
     &    450.,  572.,  -15.,
     &    450.,  562.,  277.,
     &    450.,  651.,  -45.,
     &    450.,  666., -126.,
     &    450.,  616.,  -15.,
     &    450.,  543.,  -15.,
     &    450.,  562.,  277.,
     &    450.,  616., -102.,
     &    450.,  616., -102.,
     &    450.,  601.,  219.,
     &    450.,  469.,   30.,
     &    450.,  174., -111.,
     &    450.,  313.,  203.,
     &    450.,  171.,  -15.,
     &    450.,  626.,  -77.,
     &    450.,  616.,  -15.,
     &    450.,  616.,    0.,
     &    450.,  616.,    0.,
     &    450.,  616.,  -15.,
     &    450.,  616.,    0.,
     &    450.,  601.,  -15.,
     &    450.,  616.,  -15.,
     &    450.,  601.,    0.,
     &    450.,  616.,  -15.,
     &    450.,  616.,  -15.,
     &    450.,  425.,  -15.,
     &    450.,  425., -110.,
     &    450.,  501.,   15.,
     &    450.,  389.,  109.,
     &    450.,  501.,   15.,
     &    450.,  580.,  -14./

      data  ((stwhb(12,i,j),i=1,3),j=64,95)/
     &    450.,  616.,  -15.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  580.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  580.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,  -12.,
     &    450.,  580.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  580., -138.,
     &    450.,  562.,    0.,
     &    450.,  580.,  -22.,
     &    450.,  562.,    0.,
     &    450.,  562.,  -18.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  562.,    0.,
     &    450.,  616., -102.,
     &    450.,  626.,  -77.,
     &    450.,  616., -102.,
     &    450.,  616.,  250.,
     &    450.,  -75., -125./

      data  ((stwhb(12,i,j),i=1,3),j=96,127)/
     &    450.,  562.,  277.,
     &    450.,  454.,  -15.,
     &    450.,  626.,  -15.,
     &    450.,  459.,  -15.,
     &    450.,  626.,  -15.,
     &    450.,  454.,  -15.,
     &    450.,  626.,    0.,
     &    450.,  454., -146.,
     &    450.,  626.,    0.,
     &    450.,  658.,    0.,
     &    450.,  658., -146.,
     &    450.,  626.,    0.,
     &    450.,  626.,    0.,
     &    450.,  454.,    0.,
     &    450.,  454.,    0.,
     &    450.,  454.,  -15.,
     &    450.,  454., -142.,
     &    450.,  454., -142.,
     &    450.,  454.,    0.,
     &    450.,  459.,  -17.,
     &    450.,  562.,  -15.,
     &    450.,  439.,  -15.,
     &    450.,  439.,    0.,
     &    450.,  439.,    0.,
     &    450.,  439.,    0.,
     &    450.,  439., -142.,
     &    450.,  439.,    0.,
     &    450.,  616., -102.,
     &    450.,  750., -250.,
     &    450.,  616., -102.,
     &    450.,  347.,  144.,
     &    450.,    0.,    0./

      data  ((stwhb(12,i,j),i=1,3),j=128,159)/
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

      data  ((stwhb(12,i,j),i=1,3),j=160,191)/
     &    450.,    0.,    0.,
     &    450.,  449., -146.,
     &    450.,  614.,  -49.,
     &    450.,  611.,  -28.,
     &    450.,  661.,  -60.,
     &    450.,  562.,    0.,
     &    450.,  616., -131.,
     &    450.,  580.,  -70.,
     &    450.,  517.,   49.,
     &    450.,  562.,  277.,
     &    450.,  562.,  277.,
     &    450.,  446.,   70.,
     &    450.,  446.,   70.,
     &    450.,  446.,   70.,
     &    450.,  626.,    0.,
     &    450.,  626.,    0.,
     &    450.,    0.,    0.,
     &    450.,  313.,  203.,
     &    450.,  580.,  -70.,
     &    450.,  580.,  -70.,
     &    450.,  342.,  156.,
     &    450.,    0.,    0.,
     &    450.,  580.,  -70.,
     &    450.,  430.,  132.,
     &    450.,  135., -150.,
     &    450.,  135., -150.,
     &    450.,  562.,  277.,
     &    450.,  446.,   70.,
     &    450.,  116.,  -15.,
     &    450.,  616.,  -15.,
     &    450.,    0.,    0.,
     &    450.,  449., -146./

      data  ((stwhb(12,i,j),i=1,3),j=192,223)/
     &    450.,    0.,    0.,
     &    450.,  661.,  508.,
     &    450.,  661.,  508.,
     &    450.,  657.,  483.,
     &    450.,  636.,  493.,
     &    450.,  585.,  505.,
     &    450.,  631.,  468.,
     &    450.,  625.,  485.,
     &    450.,  625.,  485.,
     &    450.,    0.,    0.,
     &    450.,  678.,  481.,
     &    450.,    0., -206.,
     &    450.,    0.,    0.,
     &    450.,  661.,  488.,
     &    450.,    0., -199.,
     &    450.,  667.,  493.,
     &    450.,  313.,  203.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

      data  ((stwhb(12,i,j),i=1,3),j=224,255)/
     &    450.,    0.,    0.,
     &    450.,  562.,    0.,
     &    450.,    0.,    0.,
     &    450.,  580.,  196.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  562.,    0.,
     &    450.,  584.,  -22.,
     &    450.,  562.,    0.,
     &    450.,  580.,  196.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  454.,  -15.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  439.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,  626.,    0.,
     &    450.,  463.,  -24.,
     &    450.,  454.,  -15.,
     &    450.,  626.,  -15.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0.,
     &    450.,    0.,    0./

*-----------------------------------------------------------------------

      data  ((stwhb(13,i,j),i=1,3),j=0,31)/
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(13,i,j),i=1,3),j=32,63)/
     &    250.,    0.,    0.,
     &    333.,  672.,  -17.,
     &    713.,  705.,    0.,
     &    500.,  673.,  -16.,
     &    549.,  707.,    0.,
     &    833.,  655.,  -36.,
     &    778.,  662.,  -18.,
     &    439.,  500.,  -17.,
     &    333.,  673., -191.,
     &    333.,  673., -191.,
     &    500.,  551.,  134.,
     &    549.,  533.,    0.,
     &    250.,  105., -152.,
     &    549.,  288.,  233.,
     &    250.,   95.,  -17.,
     &    278.,  646.,  -18.,
     &    500.,  685.,  -17.,
     &    500.,  673.,    0.,
     &    500.,  686.,    0.,
     &    500.,  685.,  -17.,
     &    500.,  685.,    0.,
     &    500.,  685.,  -17.,
     &    500.,  685.,  -18.,
     &    500.,  673.,  -16.,
     &    500.,  685.,  -18.,
     &    500.,  685.,  -18.,
     &    278.,  460.,  -17.,
     &    278.,  460., -152.,
     &    549.,  522.,    0.,
     &    549.,  390.,  141.,
     &    549.,  522.,    0.,
     &    444.,  686.,  -17./

      data  ((stwhb(13,i,j),i=1,3),j=64,95)/
     &    549.,  475.,    0.,
     &    722.,  673.,    0.,
     &    667.,  673.,    0.,
     &    722.,  673.,    0.,
     &    612.,  688.,    0.,
     &    611.,  673.,    0.,
     &    763.,  673.,    0.,
     &    603.,  673.,    0.,
     &    722.,  673.,    0.,
     &    333.,  673.,    0.,
     &    631.,  689.,  -18.,
     &    722.,  673.,    0.,
     &    686.,  688.,    0.,
     &    889.,  673.,    0.,
     &    722.,  673.,   -8.,
     &    722.,  685.,  -17.,
     &    768.,  673.,    0.,
     &    741.,  685.,  -17.,
     &    556.,  673.,    0.,
     &    592.,  673.,    0.,
     &    611.,  673.,    0.,
     &    690.,  673.,    0.,
     &    439.,  500., -233.,
     &    768.,  688.,    0.,
     &    645.,  673.,    0.,
     &    795.,  684.,    0.,
     &    611.,  673.,    0.,
     &    333.,  674., -155.,
     &    863.,  478.,    0.,
     &    333.,  674., -155.,
     &    658.,  674.,    0.,
     &    500., -206., -252./

      data  ((stwhb(13,i,j),i=1,3),j=96,127)/
     &    500.,  917.,  881.,
     &    631.,  500.,  -18.,
     &    549.,  741., -223.,
     &    549.,  500., -231.,
     &    494.,  740.,  -19.,
     &    439.,  502.,  -19.,
     &    521.,  671., -224.,
     &    411.,  499., -225.,
     &    603.,  514., -202.,
     &    329.,  503.,  -17.,
     &    603.,  499., -224.,
     &    549.,  501.,    0.,
     &    549.,  740.,  -17.,
     &    576.,  500., -223.,
     &    521.,  507.,  -16.,
     &    549.,  499.,  -19.,
     &    549.,  487.,  -19.,
     &    521.,  690.,  -17.,
     &    549.,  499., -230.,
     &    603.,  500.,  -21.,
     &    439.,  500.,  -19.,
     &    576.,  507.,  -18.,
     &    713.,  583.,  -18.,
     &    686.,  500.,  -17.,
     &    493.,  766., -225.,
     &    686.,  500., -228.,
     &    494.,  756., -225.,
     &    480.,  673., -183.,
     &    200.,  673., -177.,
     &    480.,  673., -183.,
     &    549.,  307.,  203.,
     &    250.,    0.,    0./

      data  ((stwhb(13,i,j),i=1,3),j=128,159)/
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0.,
     &    250.,    0.,    0./

      data  ((stwhb(13,i,j),i=1,3),j=160,191)/
     &    250.,    0.,    0.,
     &    620.,  685.,    0.,
     &    247.,  735.,  459.,
     &    549.,  639.,    0.,
     &    167.,  677.,  -12.,
     &    713.,  404.,  125.,
     &    500.,  687., -193.,
     &    753.,  533.,  -26.,
     &    753.,  550.,  -36.,
     &    753.,  532.,  -33.,
     &    753.,  548.,  -36.,
     &   1042.,  511.,  -15.,
     &    987.,  511.,  -15.,
     &    603.,  910.,    0.,
     &    987.,  511.,  -15.,
     &    603.,  888.,  -22.,
     &    400.,  685.,  385.,
     &    549.,  645.,    0.,
     &    411.,  737.,  459.,
     &    549.,  639.,    0.,
     &    549.,  524.,    8.,
     &    713.,  404.,  123.,
     &    494.,  746.,  -21.,
     &    460.,  473.,  113.,
     &    549.,  456.,   71.,
     &    549.,  549.,  -25.,
     &    549.,  443.,   82.,
     &    549.,  394.,  135.,
     &   1000.,   95.,  -17.,
     &    603., 1010., -120.,
     &   1000.,  276.,  220.,
     &    658.,  629.,  -16./

      data  ((stwhb(13,i,j),i=1,3),j=192,223)/
     &    823.,  658.,  -18.,
     &    686.,  740.,  -53.,
     &    795.,  734.,  -15.,
     &    987.,  573., -211.,
     &    768.,  673.,  -17.,
     &    768.,  675.,  -15.,
     &    823.,  719.,  -24.,
     &    768.,  509.,    0.,
     &    768.,  492.,  -17.,
     &    713.,  470.,    0.,
     &    713.,  470., -125.,
     &    713.,  540.,  -70.,
     &    713.,  470.,    0.,
     &    713.,  470., -125.,
     &    713.,  468.,    0.,
     &    713.,  555.,  -58.,
     &    768.,  673.,    0.,
     &    713.,  718.,  -19.,
     &    790.,  673.,  -17.,
     &    790.,  675.,  -15.,
     &    890.,  673.,  293.,
     &    823.,  751., -101.,
     &    549.,  917.,  -38.,
     &    250.,  310.,  210.,
     &    713.,  288.,    0.,
     &    603.,  454.,    0.,
     &    603.,  477.,    0.,
     &   1042.,  510.,  -20.,
     &    987.,  513.,  -15.,
     &    603.,  911.,    2.,
     &    987.,  508.,  -20.,
     &    603.,  890.,  -19./

      data  ((stwhb(13,i,j),i=1,3),j=224,255)/
     &    494.,  745.,    0.,
     &    329.,  746., -198.,
     &    790.,  670.,  -20.,
     &    790.,  675.,  -15.,
     &    786.,  673.,  293.,
     &    713.,  752., -108.,
     &    384.,  926., -293.,
     &    384.,  925.,  -85.,
     &    384.,  926., -293.,
     &    384.,  926.,  -80.,
     &    384.,  925.,  -79.,
     &    384.,  926.,  -80.,
     &    494.,  926.,  -75.,
     &    494.,  935.,  -85.,
     &    494.,  926.,  -70.,
     &    494.,  935.,  -80.,
     &    250.,    0.,    0.,
     &    329.,  746., -198.,
     &    274.,  916., -107.,
     &    686.,  922.,  -83.,
     &    686.,  975.,  -88.,
     &    686.,  921.,  -81.,
     &    384.,  926., -293.,
     &    384.,  925.,  -85.,
     &    384.,  926., -293.,
     &    384.,  926.,  -80.,
     &    384.,  925.,  -79.,
     &    384.,  926.,  -80.,
     &    494.,  926.,  -75.,
     &    494.,  935.,  -85.,
     &    494.,  926.,  -70.,
     &    250.,    0.,    0./

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine bbox(ipc,ic,xin0,yin0,ain)
*                                                                      *
*              PURPOSE : CALCULATE BOUNDING BOX                        *
*                                                                      *
*              IPC = 1 ;   pt                                          *
*              IPC = 2 ;   cm     * CM                                 *
*              IPC = 3 ;   x, y   * XAL(YAL) * CM                      *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb

      common /con/  cm, dd
      common /frm/  xal, yal

*-----------------------------------------------------------------------

            if( icb .eq. 1 ) return


            if( ipc .eq. 1 ) then

                  xin = xin0
                  yin = yin0

            else if( ipc .eq. 2 ) then

                  xin = xin0 * cm
                  yin = yin0 * cm

            else if( ipc .eq. 3 ) then

                  xin = xin0 * xal * cm
                  yin = yin0 * yal * cm

            end if

                  xin = xin * sxg0
                  yin = yin * syg0

                  dcs = cos( rog0 / 180.0 * pi )
                  dsn = sin( rog0 / 180.0 * pi )

                  xout = xog0 + ( xin * dcs - yin * dsn )
                  yout = yog0 + ( xin * dsn + yin * dcs )


            if( ic .eq. 0 ) then

               bx01 = min(bx01,xout)
               bx02 = max(bx02,xout)

               by01 = min(by01,yout)
               by02 = max(by02,yout)

            else if( ic .eq. 1 ) then

               xog0 = xout
               yog0 = yout

               rog0 = rog0 + ain

            end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine wrbox(jhf,idbg,ica,icbox,iccb,
     &                 b1x, b2x, b1y, b2y,
     &                 bdd, bcb, bcl, bdc, bcs)
*                                                                      *
*              PURPOSE : WRITE BOX ON jhf                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------
*     BOX OPERATORS
*-----------------------------------------------------------------------

*     /bxn ; single box      : ( x1 x2 y1 y2 dd cb cl bxn )
*     /bxo ; oval box        : ( x1 x2 y1 y2 dd cb cl dc bxo )
*     /bxd ; double box      : ( x1 x2 y1 y2 dd cb cl dc bxd )
*     /bxs ; shadow box      : ( x1 x2 y1 y2 dd cb cl dc cs bxs )
*     /bxp ; oval shadow box : ( x1 x2 y1 y2 dd cb cl dc cs bxp )

*      (x1,y1) and (x2,y2) : left down and right up coordinates
*      dd                  : width of line
*      cb                  : color of background
*      cl                  : color of line
*      dc                  : distance of double line or
*                            oval corner or shadow transform
*      cs                  : color of shadow

*    THE FOLLOWING BOX CONSTANTS ARE DEFINED IN BEGINING OF THE MAIN

*        BXWS =  0.25      : box small spacing
*        BXWL =  0.50      : box large spacing
*        BXSS =  0.15      : small obal box corner and shadow small size
*        BXSL =  0.25      : large obal box corner and shadow large size
*        BXDS =  0.10      : small double line distance
*        BXDL =  0.20      : large double line distance
*        BXLS =  0.028     : thin line width
*        BXLL =  0.056     : thick line width

*-----------------------------------------------------------------------

      common /fdfn/ ifd(0:13),ikan1,ikan2,ikan3,ikan4,
     &              idsm(9), idln(6), idxe, idye, ibd1, ibd2, ibd3,
     &              ibd4, ibd5, ibd6, ibd7, ibd8, ibd9, ibd10, ibd11

      dimension bcb(3), bcl(3), bcs(3)

*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

               if( icbox .ge. 10 .and. icbox .le. 13 ) then

                  if( ibd7 .eq. 0 ) then

                     if( ibd1 .eq. 0 ) then

                        call wbelm(jhf,1)

                        ibd1 = 1

                     end if

                     if( ibd4 .eq. 0 ) then

                        call wbelm(jhf,4)

                        ibd4 = 1

                     end if

                     if( ibd5 .eq. 0 ) then

                        call wbelm(jhf,5)

                        ibd5 = 1

                     end if

                        call wbelm(jhf,7)

                        ibd7 = 1

                     if( idbg .eq. 1 )  write(jhf,'()')

                  end if

                  if( iccb .eq. 0 ) then

                     if( ica .eq. 0 ) then

                        write(jhf,'(4g14.5,/g14.5,6f7.3,'' bxn'')')
     &                             b1x, b2x, b1y, b2y,
     &                             bdd, bcb, bcl

                     else

                        write(jhf,'(''xps '',g14.5,'' add xps '',g14.5,
     &                              '' add yps '',/g14.5,'' add yps '',
     &                              g14.5,'' add '',
     &                             /g14.5,6f7.3,'' bxn'')')
     &                             b1x, b2x, b1y, b2y,
     &                             bdd, bcb, bcl

                     end if

                  end if

               else if( icbox .ge. 20 .and. icbox .le. 23 ) then

                  if( ibd8 .eq. 0 ) then

                     if( ibd1 .eq. 0 ) then

                        call wbelm(jhf,1)

                        ibd1 = 1

                     end if

                     if( ibd2 .eq. 0 ) then

                        call wbelm(jhf,2)

                        ibd2 = 1

                     end if

                     if( ibd4 .eq. 0 ) then

                        call wbelm(jhf,4)

                        ibd4 = 1

                     end if

                     if( ibd6 .eq. 0 ) then

                        call wbelm(jhf,6)

                        ibd6 = 1

                     end if

                        call wbelm(jhf,8)

                        ibd8 = 1

                     if( idbg .eq. 1 )  write(jhf,'()')

                  end if

                  if( iccb .eq. 0 ) then

                     if( ica .eq. 0 ) then

                        write(jhf,'(4g14.5,
     &                             /g14.5,6f7.3,g14.5,'' bxo'')')
     &                             b1x, b2x, b1y, b2y,
     &                             bdd, bcb, bcl, bdc

                     else

                        write(jhf,'(''xps '',g14.5,'' add xps '',g14.5,
     &                              '' add yps '',/g14.5,'' add yps '',
     &                              g14.5,'' add '',
     &                             /g14.5,6f7.3,g14.5,'' bxo'')')
     &                             b1x, b2x, b1y, b2y,
     &                             bdd, bcb, bcl, bdc

                     end if

                  end if

               else if( icbox .ge. 30 .and. icbox .le. 33 ) then

                  if( ibd9 .eq. 0 ) then

                     if( ibd1 .eq. 0 ) then

                        call wbelm(jhf,1)

                        ibd1 = 1

                     end if

                     if( ibd2 .eq. 0 ) then

                        call wbelm(jhf,2)

                        ibd2 = 1

                     end if

                     if( ibd4 .eq. 0 ) then

                        call wbelm(jhf,4)

                        ibd4 = 1

                     end if

                     if( ibd5 .eq. 0 ) then

                        call wbelm(jhf,5)

                        ibd5 = 1

                     end if

                        call wbelm(jhf,9)

                        ibd9 = 1

                     if( idbg .eq. 1 )  write(jhf,'()')

                  end if

                  if( iccb .eq. 0 ) then

                     if( ica .eq. 0 ) then

                        write(jhf,'(4g14.5,
     &                             /g14.5,6f7.3,g14.5,'' bxd'')')
     &                             b1x, b2x, b1y, b2y,
     &                             bdd, bcb, bcl, bdc

                     else

                        write(jhf,'(''xps '',g14.5,'' add xps '',g14.5,
     &                              '' add yps '',/g14.5,'' add yps '',
     &                              g14.5,'' add '',
     &                             /g14.5,6f7.3,g14.5,'' bxd'')')
     &                             b1x, b2x, b1y, b2y,
     &                             bdd, bcb, bcl, bdc

                     end if

                  end if

               else if( icbox .ge. 40 .and. icbox .le. 43 ) then

                  if( ibd10 .eq. 0 ) then

                     if( ibd1 .eq. 0 ) then

                        call wbelm(jhf,1)

                        ibd1 = 1

                     end if

                     if( ibd2 .eq. 0 ) then

                        call wbelm(jhf,2)

                        ibd2 = 1

                     end if

                     if( ibd3 .eq. 0 ) then

                        call wbelm(jhf,3)

                        ibd3 = 1

                     end if

                     if( ibd4 .eq. 0 ) then

                        call wbelm(jhf,4)

                        ibd4 = 1

                     end if

                     if( ibd5 .eq. 0 ) then

                        call wbelm(jhf,5)

                        ibd5 = 1

                     end if

                        call wbelm(jhf,10)

                        ibd10 = 1

                     if( idbg .eq. 1 )  write(jhf,'()')

                  end if

                  if( iccb .eq. 0 ) then

                     if( ica .eq. 0 ) then

                        write(jhf,'(5g14.5,
     &                             /6f7.3,g14.5,3f7.3,'' bxs'')')
     &                             b1x, b2x, b1y, b2y,
     &                             bdd, bcb, bcl, bdc, bcs

                     else

                        write(jhf,'(''xps '',g14.5,'' add xps '',g14.5,
     &                              '' add yps '',/g14.5,'' add yps '',
     &                              g14.5,'' add '',g14.5,
     &                             /6f7.3,g14.5,3f6.3,'' bxs'')')
     &                             b1x, b2x, b1y, b2y,
     &                             bdd, bcb, bcl, bdc, bcs

                     end if

                  end if

               else if( icbox .ge. 50 .and. icbox .le. 53 ) then

                  if( ibd11 .eq. 0 ) then

                     if( ibd1 .eq. 0 ) then

                        call wbelm(jhf,1)

                        ibd1 = 1

                     end if

                     if( ibd2 .eq. 0 ) then

                        call wbelm(jhf,2)

                        ibd2 = 1

                     end if

                     if( ibd3 .eq. 0 ) then

                        call wbelm(jhf,3)

                        ibd3 = 1

                     end if

                     if( ibd4 .eq. 0 ) then

                        call wbelm(jhf,4)

                        ibd4 = 1

                     end if

                     if( ibd6 .eq. 0 ) then

                        call wbelm(jhf,6)

                        ibd6 = 1

                     end if

                        call wbelm(jhf,11)

                        ibd11 = 1

                     if( idbg .eq. 1 )  write(jhf,'()')

                  end if

                  if( iccb .eq. 0 ) then

                     if( ica .eq. 0 ) then

                        write(jhf,'(5g14.5,
     &                             /6f7.3,g14.5,3f7.3,'' bxp'')')
     &                             b1x, b2x, b1y, b2y,
     &                             bdd, bcb, bcl, bdc, bcs

                     else

                        write(jhf,'(''xps '',g14.5,'' add xps '',g14.5,
     &                              '' add yps '',/g14.5,'' add yps '',
     &                              g14.5,'' add '',g14.5,
     &                             /6f7.3,g14.5,3f7.3,'' bxp'')')
     &                             b1x, b2x, b1y, b2y,
     &                             bdd, bcb, bcl, bdc, bcs

                     end if

                  end if

               end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine wbelm(jhf,ibd)
*                                                                      *
*              PURPOSE : WRITE BOX OPERATORS                           *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      if(      ibd .eq.  1 ) then

         write(jhf,'(
     & '' /bxe {/bcl3 S N /bcl2 S N /bcl1 S N /bcb3 S N /bcb2 S N '',
     & ''/bcb1 S N /bdd S N ''/
     & ''/by2 S N /by1 S N /bx2 S N /bx1 S N} N''
     &   )')

      else if( ibd .eq.  2 ) then

         write(jhf,'(
     & ''/bxf {/bx3 bx1 bdc add N /bx4 bx2 bdc sub N /by3 by1 bdc '',
     & ''add N''/
     & ''/by4 by2 bdc sub N} N''
     &   )')

      else if( ibd .eq.  3 ) then

         write(jhf,'(
     & ''/bxg {/bx5 bx2 bdc add N /by5 by1 bdc sub N} N''
     &   )')

      else if( ibd .eq.  4 ) then

         write(jhf,'(
     & ''/bxh {/yp2 S N /yp1 S N /xp2 S N /xp1 S N} N''
     &   )')

      else if( ibd .eq.  5 ) then

         write(jhf,'(
     & ''/bxj {bxh xp1 yp1 M xp2 yp1 L xp2 yp2 L xp1 yp2 L cp} N''
     &   )')

      else if( ibd .eq.  6 ) then

         write(jhf,'(
     & ''/bxk {/yp4 S N /yp3 S N /xp4 S N /xp3 S N bxh xp4 yp1 M '',
     & ''xp4 yp3 dv''/
     & ''bdc 270 360 arc xp2 yp3 L xp2 yp4 L xp4 yp4 dv bdc 0 90 '',
     & ''arc xp4 yp2 L''/
     & ''xp3 yp2 L xp3 yp4 dv bdc 90 180 arc xp1 yp4 L xp1 yp3 L '',
     & ''xp3 yp3 dv''/
     & ''bdc 180 270 arc xp3 yp1 L cp} N''
     &   )')

      else if( ibd .eq.  7 ) then

         write(jhf,'(
     & ''/bxn {bxe bcb1 -2 ge {bcb1 bcb2 bcb3 sc bx1 bx2 by1 by2 '',
     & ''bxj fl} if ''/
     & ''bcl1 bcl2 bcl3 sc sd0 bdd lw bx1 bx2 by1 by2 bxj st} N''
     &   )')

      else if( ibd .eq.  8 ) then

         write(jhf,'(
     & ''/bxo {/bdc S N bxe bxf bcb1 -2 ge {bcb1 bcb2 bcb3 sc '',
     & ''bx1 bx2 by1 by2 ''/
     & ''bx3 bx4 by3 by4 bxk fl} if bcl1 bcl2 bcl3 sc sd0 bdd lw '',
     & ''bx1 bx2 by1 by2 ''/
     & ''bx3 bx4 by3 by4 bxk st} N''
     &   )')

      else if( ibd .eq.  9 ) then

         write(jhf,'(
     & ''/bxd {/bdc S N bxe bxf bcb1 -2 ge {bcb1 bcb2 bcb3 sc '',
     & ''bx1 bx2 by1 by2 ''/
     & ''bxj fl} if bcl1 bcl2 bcl3 sc sd0 bdd lw bx1 bx2 by1 by2 '',
     & ''bxj st bx3 bx4 ''/
     & ''by3 by4 bxj st} N''
     &   )')

      else if( ibd .eq. 10 ) then

         write(jhf,'(
     & ''/bxs {/bcs3 S N /bcs2 S N /bcs1 S N /bdc S N bxe bxf bxg '',
     & ''bcs1 bcs2 bcs3 sc ''/
     & ''bx3 bx5 by5 by4 bxj fl bcl1 bcl2 bcl3 sc sd0 bdd lw bx3 '',
     & ''bx5 by5 by4 bxj st ''/
     & ''bcb1 bcb2 bcb3 sc bx1 bx2 by1 by2 bxj fl bcl1 bcl2 bcl3 sc '',
     & ''sd0 bdd lw ''/
     & ''bx1 bx2 by1 by2 bxj st} N''
     &   )')

      else if( ibd .eq. 11 ) then

         write(jhf,'(
     & ''/bxp {/bcs3 S N /bcs2 S N /bcs1 S N /bdc S N bxe bxf bxg '',
     & ''/bx6 bx3 bdc add N ''/
     & ''/by6 by4 bdc sub N bcs1 bcs2 bcs3 sc bx3 bx5 by5 by4 bx6 '',
     & ''bx2 by1 by6 bxk fl ''/
     & ''bcl1 bcl2 bcl3 sc sd0 bdd lw bx3 bx5 by5 by4 bx6 bx2 by1 '',
     & ''by6 bxk st ''/
     & ''bcb1 bcb2 bcb3 sc bx1 bx2 by1 by2 bx3 bx4 by3 by4 bxk fl '',
     & ''bcl1 bcl2 bcl3 sc sd0 bdd lw bx1 bx2 by1 by2 bx3 bx4 by3 '',
     & ''by4 bxk st} N''
     &   )')

      end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine ftdfn(jhf,ifon,ifnn,idbg)
*                                                                      *
*              PURPOSE : FIRST DEFINITION OF FONTS                     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /fdfn/ ifd(0:13),ikan1,ikan2,ikan3,ikan4,
     &              idsm(9), idln(6), idxe, idye, ibd1, ibd2, ibd3,
     &              ibd4, ibd5, ibd6, ibd7, ibd8, ibd9, ibd10, ibd11

      character jfont(12)*100
      dimension jfin(12)
      character kfont(2)*9

      data jfont/
     &           'fnm2 /GothicBBB-Medium-RKSJ-H ',
     &           'fnm2 /GothicBBB-Medium-H ',
     &           'fnm2 /GothicBBB-Medium-EUC-H ',
     &           'fnm2 /Ryumin-Light-RKSJ-H ',
     &           'fnm2 /Ryumin-Light-H ',
     &           'fnm2 /Ryumin-Light-EUC-H ',
     &           'fnm3 /GothicBBB-Medium-RKSJ-H ',
     &           'fnm3 /GothicBBB-Medium-H ',
     &           'fnm3 /GothicBBB-Medium-EUC-H ',
     &           'fnm3 /Ryumin-Light-RKSJ-H ',
     &           'fnm3 /Ryumin-Light-H ',
     &           'fnm3 /Ryumin-Light-EUC-H '/

      data jfin/ 30,25,29,26,21,25,30,25,29,26,21,25/

      data kfont/
     &           'fdf1 fend',
     &           'fdf2 fend'/

*-----------------------------------------------------------------------
*     FONTS PARAMETERS
*-----------------------------------------------------------------------

*        IFON = 0 ; ENGLISH ONLY
*        IFON > 0 ; JAPANESE AND ENGLISH

*        IFON = 1 ; Shift JIS
*        IFON = 2 ; Normal JIS
*        IFON = 3 ; EUC JIS

*        Courier resize = 0.75
*        Kanji   resize = 0.9
*        Kanji   Italic = 0.2

*-----------------------------------------------------------------------

         if( idbg .eq. 1 )  write(jhf,'()')

         if( ifon .eq. 0 ) then

*-----------------------------------------------------------------------

            if(      ifnn .eq.  0 ) then

               write(jhf,'(
     &         ''/ft000 {/Helvetica findfont} N'')')

            else if( ifnn .eq.  1 ) then

               write(jhf,'(
     &         ''/ft001 {/Helvetica-Oblique findfont} N'')')

            else if( ifnn .eq.  2 ) then

               write(jhf,'(
     &         ''/ft002 {/Helvetica-Bold findfont} N'')')

            else if( ifnn .eq.  3 ) then

               write(jhf,'(
     &         ''/ft003 {/Helvetica-BoldOblique findfont} N'')')

            else if( ifnn .eq.  4 ) then

               write(jhf,'(
     &         ''/ft004 {/Symbol findfont} N'')')

            else if( ifnn .eq.  5 ) then

               write(jhf,'(
     &         ''/ft005 {/Times-Roman findfont} N'')')

            else if( ifnn .eq.  6 ) then

               write(jhf,'(
     &         ''/ft006 {/Times-Italic findfont} N'')')

            else if( ifnn .eq.  7 ) then

               write(jhf,'(
     &         ''/ft007 {/Times-Bold findfont} N'')')

            else if( ifnn .eq.  8 ) then

               write(jhf,'(
     &         ''/ft008 {/Times-BoldItalic findfont} N'')')

            else if( ifnn .eq.  9 ) then

               write(jhf,'(
     &         ''/ft009 {/Courier findfont '',
     &         ''[0.75 0 0 1 0 0] makefont} N'')')

            else if( ifnn .eq. 10 ) then

               write(jhf,'(
     &         ''/ft010 {/Courier-Oblique findfont '',
     &         ''[0.75 0 0 1 0 0] makefont} N'')')

            else if( ifnn .eq. 11 ) then

               write(jhf,'(
     &         ''/ft011 {/Courier-Bold findfont '',
     &         ''[0.75 0 0 1 0 0] makefont} N'')')

            else if( ifnn .eq. 12 ) then

               write(jhf,'(
     &         ''/ft012 {/Courier-BoldOblique findfont '',
     &         ''[0.75 0 0 1 0 0] makefont} N'')')

            else if( ifnn .eq. 13 ) then

               write(jhf,'(
     &         ''/ft013 {/Symbol findfont '',
     &         ''[1 0 0.2 1 0 0] makefont} N'')')

            end if

*-----------------------------------------------------------------------

         else if( ifon .ne. 0 ) then

*-----------------------------------------------------------------------

            if(      ifnn .eq.  0 ) then

               write(jhf,'(
     &         ''/ft000 {/Gothic-Helvetica findfont} N'')')

            else if( ifnn .eq.  1 ) then

               write(jhf,'(
     &         ''/ft001 {/Gothic-Helvetica-Oblique findfont} N'')')

            else if( ifnn .eq.  2 ) then

               write(jhf,'(
     &         ''/ft002 {/Gothic-Helvetica-Bold findfont} N'')')

            else if( ifnn .eq.  3 ) then

               write(jhf,'(
     &         ''/ft003 {/Gothic-Helvetica-BoldOblique findfont} N'')')

            else if( ifnn .eq.  4 ) then

               write(jhf,'(
     &         ''/ft004 {/Gothic-Symbol findfont} N'')')

            else if( ifnn .eq.  5 ) then

               write(jhf,'(
     &         ''/ft005 {/Ryumin-Times-Roman findfont} N'')')

            else if( ifnn .eq.  6 ) then

               write(jhf,'(
     &         ''/ft006 {/Ryumin-Times-Italic findfont} N'')')

            else if( ifnn .eq.  7 ) then

               write(jhf,'(
     &         ''/ft007 {/Gothic-Times-Bold findfont} N'')')

            else if( ifnn .eq.  8 ) then

               write(jhf,'(
     &         ''/ft008 {/Gothic-Times-BoldItalic findfont} N'')')

            else if( ifnn .eq.  9 ) then

               write(jhf,'(
     &         ''/ft009 {/Ryumin-Courier findfont} N'')')

            else if( ifnn .eq. 10 ) then

               write(jhf,'(
     &         ''/ft010 {/Ryumin-Courier-Oblique findfont} N'')')

            else if( ifnn .eq. 11 ) then

               write(jhf,'(
     &         ''/ft011 {/Gothic-Courier-Bold findfont} N'')')

            else if( ifnn .eq. 12 ) then

               write(jhf,'(
     &         ''/ft012 {/Gothic-Courier-BoldOblique findfont} N'')')

            else if( ifnn .eq. 13 ) then

               write(jhf,'(
     &         ''/ft013 {/Symbol findfont '',
     &         ''[1 0 0.2 1 0 0] makefont} N'')')

            end if

*-----------------------------------------------------------------------

            if( ikan1 .eq. 0 .and. ifnn .le. 8 ) then

               write(jhf,'(''/fnm2 {findfont} N'')')

               ikan1 = 1

            end if

            if( ikan2 .eq. 0 .and. ifnn .ge. 9 .and. ifnn .le. 12 ) then

               write(jhf,'(''/fnm3 {findfont [0.75 0 0 1 0 0] '',
     &                     ''makefont} N'')')

               ikan2 = 1

            end if

            if( ikan3 .eq.  0 .and.
     &        ( ifnn  .eq.  0 .or.
     &          ifnn  .eq.  2 .or.
     &          ifnn  .eq.  4 .or.
     &          ifnn  .eq.  5 .or.
     &          ifnn  .eq.  7 .or.
     &          ifnn  .eq.  9 .or.
     &          ifnn  .eq. 11 ) ) then

               write(jhf,'(''/fdf1 {findfont [0.9 0 0 0.9 0 0] '',
     &                     ''makefont] N} N'')')

               ikan3 = 1

            end if

            if( ikan4 .eq.  0 .and.
     &        ( ifnn  .eq.  1 .or.
     &          ifnn  .eq.  3 .or.
     &          ifnn  .eq.  6 .or.
     &          ifnn  .eq.  8 .or.
     &          ifnn  .eq. 10 .or.
     &          ifnn  .eq. 12 ) ) then

               write(jhf,'(''/fdf2 {findfont [0.9 0 0.2 0.9 0 0] '',
     &                     ''makefont] N} N'')')

               ikan4 = 1

            end if


*-----------------------------------------------------------------------

            ij = ifon

*-----------------------------------------------------------------------

            if(      ifnn .eq.  0 ) then

               write(jhf,'(
     &         ''ftnm /Gothic-Helvetica '',
     &         ''fnm1 /Helvetica''
     &         )')
               write(jhf,'(A)') jfont(ij)(1:jfin(ij))//kfont(1)

            else if( ifnn .eq.  1 ) then

               write(jhf,'(
     &         ''ftnm /Gothic-Helvetica-Oblique '',
     &         ''fnm1 /Helvetica-Oblique''
     &         )')
               write(jhf,'(A)') jfont(ij)(1:jfin(ij))//kfont(2)

            else if( ifnn .eq.  2 ) then

               write(jhf,'(
     &         ''ftnm /Gothic-Helvetica-Bold '',
     &         ''fnm1 /Helvetica-Bold''
     &         )')
               write(jhf,'(A)') jfont(ij)(1:jfin(ij))//kfont(1)

            else if( ifnn .eq.  3 ) then

               write(jhf,'(
     &         ''ftnm /Gothic-Helvetica-BoldOblique '',
     &         ''fnm1 /Helvetica-BoldOblique''
     &         )')
               write(jhf,'(A)') jfont(ij)(1:jfin(ij))//kfont(2)

            else if( ifnn .eq.  4 ) then

               write(jhf,'(
     &         ''ftnm /Gothic-Symbol '',
     &         ''fnm1 /Symbol''
     &         )')
               write(jhf,'(A)') jfont(ij)(1:jfin(ij))//kfont(1)

            else if( ifnn .eq.  5 ) then

               write(jhf,'(
     &         ''ftnm /Ryumin-Times-Roman '',
     &         ''fnm1 /Times-Roman''
     &         )')
               write(jhf,'(A)') jfont(ij+3)(1:jfin(ij+3))//kfont(1)

            else if( ifnn .eq.  6 ) then

               write(jhf,'(
     &         ''ftnm /Ryumin-Times-Italic '',
     &         ''fnm1 /Times-Italic''
     &         )')
               write(jhf,'(A)') jfont(ij+3)(1:jfin(ij+3))//kfont(2)

            else if( ifnn .eq.  7 ) then

               write(jhf,'(
     &         ''ftnm /Gothic-Times-Bold '',
     &         ''fnm1 /Times-Bold''
     &         )')
               write(jhf,'(A)') jfont(ij)(1:jfin(ij))//kfont(1)

            else if( ifnn .eq.  8 ) then

               write(jhf,'(
     &         ''ftnm /Gothic-Times-BoldItalic '',
     &         ''fnm1 /Times-BoldItalic''
     &         )')
               write(jhf,'(A)') jfont(ij)(1:jfin(ij))//kfont(2)

            else if( ifnn .eq.  9 ) then

               write(jhf,'(
     &         ''ftnm /Ryumin-Courier '',
     &         ''fnm1 /Courier''
     &         )')
               write(jhf,'(A)') jfont(ij+9)(1:jfin(ij+9))//kfont(1)

            else if( ifnn .eq. 10 ) then

               write(jhf,'(
     &         ''ftnm /Ryumin-Courier-Oblique '',
     &         ''fnm1 /Courier-Oblique''
     &         )')
               write(jhf,'(A)') jfont(ij+9)(1:jfin(ij+9))//kfont(2)

            else if( ifnn .eq. 11 ) then

               write(jhf,'(
     &         ''ftnm /Gothic-Courier-Bold '',
     &         ''fnm1 /Courier-Bold''
     &         )')
               write(jhf,'(A)') jfont(ij+6)(1:jfin(ij+6))//kfont(1)

            else if( ifnn .eq. 12 ) then

               write(jhf,'(
     &         ''ftnm /Gothic-Courier-BoldOblique '',
     &         ''fnm1 /Courier-BoldOblique''
     &         )')
               write(jhf,'(A)') jfont(ij+6)(1:jfin(ij+6))//kfont(2)

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

               ifd(ifnn) = 1

               if( idbg .eq. 1 )  write(jhf,'()')

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine lndfn(jhf,ityp,idbg)
*                                                                      *
*              PURPOSE : FIRST DEFINITION OF LINES                     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------
*     LINE OPERATORS
*-----------------------------------------------------------------------

*     /sd0 ; Solid Line
*     /sd1 ; Dot-Dashed Line
*     /sd2 ; Short-Dashed Line
*     /sd3 ; Long-Dashed Line
*     /sd4 ; Dotted Line
*     /sd5 ; Dot-Dot-Dashed Line
*     /sd6 ; Dot-Dot-Dot-Dashed Line

*     FOR LINES { /lu length N} IS NECESSARY
*     WHICH IS UNIT OF THE PATTERN

*-----------------------------------------------------------------------

         if(      ityp .eq. 1 ) then

            write(jhf,'(
     &      ''/sd1 {[ .375 lu dw sub .125 lu dw add lm .125 lu dw '',
     &      ''add ] 0 sd} N''
     &      )')

         else if( ityp .eq. 2 ) then

            write(jhf,'(
     &      ''/sd2 {[ .360 lu dw 2 div sub .280 lu dw add ] 0 sd} N''
     &      )')

         else if( ityp .eq. 3 ) then

            write(jhf,'(
     &      ''/sd3 {[ .440 lu dw 2 div sub .120 lu dw add ] 0 sd} N''
     &      )')

         else if( ityp .eq. 4 ) then

            write(jhf,'(
     &      ''/sd4 {[ lm .5 lu ] 0 sd} N''
     &      )')

         else if( ityp .eq. 5 ) then

            write(jhf,'(
     &      ''/sd5 {[ .3125 lu dw sub .125 lu dw 3 div 2 mul add '',
     &      ''lm .125 lu dw 3''/
     &      ''div 2 mul add lm .125 lu dw 3 div 2 mul add ] 0 sd} N''
     &      )')

         else if( ityp .eq. 6 ) then

            write(jhf,'(
     &      ''/sd6 {[ .250 lu dw sub .125 lu dw 2 div add lm .125 '',
     &      ''lu dw 2 div add lm''/
     &      ''.125 lu dw 2 div add lm .125 lu dw 2 div add ] 0 sd} N''
     &      )')

         end if

         if( idbg .eq. 1 )  write(jhf,'()')

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine smdfn(jhf,imar,idbg)
*                                                                      *
*              PURPOSE : FIRST DEFINITION OF SYMBOLS                   *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------
*     SYMBOL OPERATORS
*-----------------------------------------------------------------------

*     /s01 ; Dot
*     /s02 ; Plus
*     /s03 ; Open Circle
*     /s04 ; Solid Circle
*     /s05 ; Open Square
*     /s06 ; Solid Square
*     /s07 ; Open Triangle (up)
*     /s08 ; Solid Triangle (up)
*     /s09 ; Open Triangle (down)
*     /s10 ; Solid Triangle (down)
*     /s11 ; Open Diamond (slim)
*     /s12 ; Solid Diamond (slim)
*     /s13 ; Open Diamond (square)
*     /s14 ; Solid Diamond (square)
*     /s15 ; Cross
*     /s16 ; Asterisk
*     /s17 ; Open Circle with color inside
*     /s18 ; Open Square with color inside
*     /s19 ; Open Triangle (up) with color inside
*     /s20 ; Open Triangle (down) with color inside
*     /s21 ; Open Diamond (slim) with color inside
*     /s13 ; Open Diamond (square) with color inside

*     FOR SYMBOLS [ /ss size N /sw width N ] IS NECESSARY
*     WHICH ARE THE SIZE OF SYMBOLS AND WIDTH OF THE STROKE


      common /fdfn/ ifd(0:13),ikan1,ikan2,ikan3,ikan4,
     &              idsm(9), idln(6), idxe, idye, ibd1, ibd2, ibd3,
     &              ibd4, ibd5, ibd6, ibd7, ibd8, ibd9, ibd10, ibd11

*-----------------------------------------------------------------------

      if( imar .eq. 2 ) then

         if( idsm(1) .eq. 0 ) then

            write(jhf,'(
     &      ''/pls {/yps S y N /xps S x N xps yps sm 2 div add M'',
     &      '' xps yps sm 2 div''/
     &      ''sub L xps sm 2 div sub yps M xps sm 2 div add yps L} N''
     &      )')

            idsm(1) = 1

         end if

      else if( imar .eq.  1 .or. imar .eq.  3 .or. imar .eq.  4 .or.
     &         imar .eq. 17 ) then

         if( idsm(2) .eq. 0 ) then

            write(jhf,'(
     &      ''/cir {/yps S y N /xps S x N xps sm 2 div add yps M'',
     &      '' xps yps sm 2 div''/
     &      ''0 360 arc cp} N''
     &      )')

            idsm(2) = 1

         end if

      else if( imar .eq.  5 .or. imar .eq.  6 .or. imar .eq. 18 ) then

         if( idsm(3) .eq. 0 ) then

            write(jhf,'(
     &      ''/squ {/yps S y sm 2 div sub N /xps S x sm 2 div sub N'',
     &      '' xps yps M xps''/
     &      ''sm add yps L xps sm add yps sm add L xps yps sm add L'',
     &      '' cp} N''
     &      )')

            idsm(3) = 1

         end if

      else if( imar .eq.  7 .or. imar .eq.  8 .or. imar .eq. 19 ) then

         if( idsm(4) .eq. 0 ) then

            write(jhf,'(
     &      ''/tr1 {/yps S y sm 0.57735 mul add N /xps S x N xps'',
     &      '' yps M /yps yps sm''/
     &      ''0.866025 mul sub N xps sm 2 div sub yps L xps sm 2'',
     &      '' div add yps L cp} N''
     &      )')

            idsm(4) = 1

         end if

      else if( imar .eq.  9 .or. imar .eq. 10 .or. imar .eq. 20 ) then

         if( idsm(5) .eq. 0 ) then

            write(jhf,'(
     &      ''/tr2 {/yps S y sm 0.57735 mul sub N /xps S x N xps'',
     &      '' yps M /yps yps sm''/
     &      ''0.866025 mul add N xps sm 2 div add yps L xps sm 2'',
     &      '' div sub yps L cp} N''
     &      )')

            idsm(5) = 1

         end if

      else if( imar .eq. 11 .or. imar .eq. 12 .or. imar .eq. 21 ) then

         if( idsm(6) .eq. 0 ) then

            write(jhf,'(
     &      ''/di1 {/yps S y N /xps S x N xps yps sm 0.866025 mul'',
     &      '' add M xps sm 2 div''/
     &      ''sub yps L xps yps sm 0.866025 mul sub L xps sm 2 div'',
     &      '' add yps L cp} N''
     &      )')

            idsm(6) = 1

         end if

      else if( imar .eq. 13 .or. imar .eq. 14 .or. imar .eq. 22 ) then

         if( idsm(7) .eq. 0 ) then

            write(jhf,'(
     &      ''/di2 {/yps S y N /xps S x N xps yps sm add M xps sm'',
     &      '' sub yps L''/
     &      ''xps yps sm sub L xps sm add yps L cp} N''
     &      )')

            idsm(7) = 1

         end if

      else if( imar .eq. 15 ) then

         if( idsm(8) .eq. 0 ) then

            write(jhf,'(
     &      ''/crs {/yps S y sm 2 div add N /xps S x sm 2 div sub'',
     &      '' N xps yps M xps''/
     &      ''sm add yps sm sub L xps yps sm sub M xps sm add yps L} N''
     &      )')

            idsm(8) = 1

         end if

      else if( imar .eq. 16 ) then

         if( idsm(9) .eq. 0 ) then

            write(jhf,'(
     &      ''/ast {y sth sub S x stw sub S M (*) show} N /smc'',
     &      '' {/stw sm 0.25 mul N''/
     &      ''/sth 0.4705 sm mul N /Times-Roman findfont sm'',
     &      '' scalefont setfont} N''
     &      )')

            idsm(9) = 1

         end if

      end if


         if( idbg .eq. 1 )  write(jhf,'()')

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine inclp(inps,dpin,idpi,ill,jsn,dum,lum,k,ierr,
     &                 psbbx,xyps,ixps,iyps)
*                                                                      *
*              purpose : define include ps file                        *
*                                                                      *
*              inps: {file.name} x() y() ix() iy() s() sx() sy() a()   *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dpin(ipsm)*200
      dimension idpi(ipsm)
      dimension psbbx(ipsm,4)
      dimension xyps(ipsm,6)
      dimension ixps(ipsm), iyps(ipsm)

      character dum(ichrl)*1
      character lum(ichrl)*1
      character dup(ichrl)*1
      dimension ill(0:9)

      character c1*1, c2*1

      character m_err*200
      common /error/ m_err, l_err, k_err


      logical exex
      logical dnen2, deqn3

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

         c1 = '('
         c2 = ')'

*-----------------------------------------------------------------------

         ierr = 0

            inps = inps + 1

            if( inps .gt. ipsm ) then

               m_err = 'Number of Include PS File is too many. '//
     &                 ' Max Number is 100 in one Graph.'
               ErrCha = ''
               ErrID = 'L:19883/R:inclp/F:a-main1.f'
               l_err = ill(jsn)
               k_err = jsn

               ierr = 1

               return

            end if

*-----------------------------------------------------------------------

   10    k = k + 1

               if( k .gt. icolm ) goto 999

            if( lum(k) .eq. '{' ) then

               ini = k + 1

            else if( lum(k) .eq. '}' ) then

               inf = k - 1

               goto 11

            end if

               goto 10

*-----------------------------------------------------------------------

   11          continue

               if( ini .gt. inf ) goto 999


            do 20 i = 1, 200

                  dpin(inps)(i:i) = ' '

   20       continue

            do 30 i = ini, inf

                  j = i - ini + 1

                  dpin(inps)(j:j) = dum(i)

   30       continue

                  idpi(inps) = inf - ini + 1

*-----------------------------------------------------------------------
*     INPUT FILE DPIN EXIST ?
*-----------------------------------------------------------------------

            inquire( file = dpin(inps), exist = exex )

            if( exex .eqv. .false. ) then

               m_err = 'Input PS File Name Error in INPS: .'//
     &                 ' (File not exist)'
               ErrCha = ''
               ErrID = 'L:19947/R:inclp/F:a-main1.f'
               l_err = ill(jsn)
               k_err = jsn

               ierr  = 1

               return

            end if

*-----------------------------------------------------------------------
*     OPEN INCLUDE FILE
*-----------------------------------------------------------------------

            open(jps, file = dpin(inps), status = 'OLD' )

*-----------------------------------------------------------------------
*     READ THE INFORMATION OF THE BOUNDING BOX
*-----------------------------------------------------------------------

         ilp = 0
         ilb = 0

  100    ilp = ilp + 1

            read(jps,'(10000a1)',iostat=ios ) (dup(ic),ic=1,icolm)
            if( ios .eq. -1 ) goto 200

            if( dup( 1) .eq. '%' .and.
     &          dup( 2) .eq. '%' .and.
     &          dup( 3) .eq. 'B' .and.
     &          dup( 4) .eq. 'o' .and.
     &          dup( 5) .eq. 'u' .and.
     &          dup( 6) .eq. 'n' .and.
     &          dup( 7) .eq. 'd' .and.
     &          dup( 8) .eq. 'i' .and.
     &          dup( 9) .eq. 'n' .and.
     &          dup(10) .eq. 'g' .and.
     &          dup(11) .eq. 'B' .and.
     &          dup(12) .eq. 'o' .and.
     &          dup(13) .eq. 'x' .and.
     &          dup(14) .eq. ':' ) then

               inum =  0
               ic   = 14

  101          ic = ic + 1

                  if( ic .gt. icolm .and. inum .ne. 4 ) goto 100

                  if( ic .gt. icolm ) then

                     ilb = ilb + 1

                     goto 100

                  end if

                  if( dup(ic) .eq. ' ' ) goto 101

                  if( dnen2(dup(ic)) ) goto 100


               inum = inum + 1

               ici = ic

  102          ic = ic + 1

                  if( ic .gt. icolm .and. inum .ne. 4 ) goto 100

                  if( ic .gt. icolm ) then

                     call rnum(rrnm,dup,ici,ic-1,ierr)
                        if(ierr.ne.0) goto 200

                     psbbx(inps,inum) = rrnm

                     ilb = ilb + 1

                     goto 100

                  end if


                  if( dup(ic) .eq. ' ' ) then

                     call rnum(rrnm,dup,ici,ic-1,ierr)
                        if(ierr.ne.0) goto 200

                     psbbx(inps,inum) = rrnm

                     goto 101

                  end if

                  if( deqn3(dup(ic)) ) goto 102

                  goto 100

*-----------------------------------------------------------------------

            else

                  goto 100

            end if


  200      if( ilb .eq. 0 ) then

               m_err = 'Input PS File is not EPS File. '//
     &                 '%%BoundingBox: is missing.'
               ErrCha = ''
               ErrID = 'L:20061/R:inclp/F:a-main1.f'
               l_err = ill(jsn)
               k_err = jsn

               ierr  = 1

               return

            end if

*-----------------------------------------------------------------------
*     READ PARAMETERS
*-----------------------------------------------------------------------

      ic = inf + 1

  300 ic = ic + 1

         if( ic .gt. icolm ) goto 400

         if(lum(ic).eq.' '.or.lum(ic).eq.tub) goto 300

               if(lum(ic).ne.'x'.and.lum(ic).ne.'y'.and.
     &            lum(ic).ne.'s'.and.lum(ic).ne.'a'.and.
     &            lum(ic).ne.'i') then

                  m_err = 'INPS: {filename} X() Y() S() IX() IY()'//
     &                    ' A() SX() SY() : Unexpected Parameter.'
                  ErrCha = ''
                  ErrID = 'L:20090/R:inclp/F:a-main1.f'
                  goto 999

               end if


               if(lum(ic).eq.'x') then

                     ic=ic+1
                     call pnum(lum,ic,icolm,rnm,ierr)
                     if(ierr.ne.0) then
                        m_err = 'INPS: Number in X( ) is Wrong.'
                        ErrCha = ''
                        ErrID = 'L:20103/R:inclp/F:a-main1.f'
                        goto 999
                     end if
                     xyps(inps,1)=rnm
                     ic=ic+1

               else if(lum(ic).eq.'y') then

                     ic=ic+1
                     call pnum(lum,ic,icolm,rnm,ierr)
                     if(ierr.ne.0) then
                        m_err = 'INPS: Number in Y( ) is Wrong.'
                        ErrCha = ''
                        ErrID = 'L:20116/R:inclp/F:a-main1.f'
                        goto 999
                     end if
                     xyps(inps,2)=rnm
                     ic=ic+1

               else if(lum(ic).eq.'a') then

                     ic=ic+1
                     call pnum(lum,ic,icolm,rnm,ierr)
                     if(ierr.ne.0) then
                        m_err = 'INPS: Number in A( ) is Wrong.'
                        ErrCha = ''
                        ErrID = 'L:20129/R:inclp/F:a-main1.f'
                        goto 999
                     end if
                     xyps(inps,3)=rnm
                     ic=ic+1

               else if( lum(ic)//lum(ic+1) .eq. 's(' ) then

                     ic=ic+1
                     call pnum(lum,ic,icolm,rnm,ierr)
                     if(ierr.ne.0) then
                        m_err = 'INPS: Number in S( ) is Wrong.'
                        ErrCha = ''
                        ErrID = 'L:20142/R:inclp/F:a-main1.f'
                        goto 999
                     end if
                     xyps(inps,4)=rnm
                     ic=ic+1

               else if(lum(ic)//lum(ic+1).eq.'sx') then

                     ic=ic+2
                     call pnum(lum,ic,icolm,rnm,ierr)
                     if(ierr.ne.0) then
                        m_err = 'INPS: Number in SX( ) is Wrong.'
                        ErrCha = ''
                        ErrID = 'L:20155/R:inclp/F:a-main1.f'
                        goto 999
                     end if
                     xyps(inps,5)=rnm
                     ic=ic+1

               else if(lum(ic)//lum(ic+1).eq.'sy') then

                     ic=ic+2
                     call pnum(lum,ic,icolm,rnm,ierr)
                     if(ierr.ne.0) then
                        m_err = 'INPS: Number in SY( ) is Wrong.'
                        ErrCha = ''
                        ErrID = 'L:20168/R:inclp/F:a-main1.f'
                        goto 999
                     end if
                     xyps(inps,6)=rnm
                     ic=ic+1

               else if(lum(ic)//lum(ic+1).eq.'ix') then

                     ic=ic+2
                     call pnum(lum,ic,icolm,rnm,ierr)
                     if(ierr.ne.0) then
                        m_err = 'INPS: Number in IX( ) is Wrong.'
                        ErrCha = ''
                        ErrID = 'L:20181/R:inclp/F:a-main1.f'
                        goto 999
                     end if
                     if(nint(rnm).ne.1.and.nint(rnm).ne.2.and.
     &                  nint(rnm).ne.3) then
                        m_err = 'INPS: Number in IX() is 1(Left), '//
     &                          '2(Center), or 3(Right).'
                        ErrCha = ''
                        ErrID = 'L:20189/R:inclp/F:a-main1.f'
                        goto 999
                     end if
                     ixps(inps)=nint(rnm)
                     ic=ic+1

               else if(lum(ic)//lum(ic+1).eq.'iy') then

                     ic=ic+2
                     call pnum(lum,ic,icolm,rnm,ierr)
                     if(ierr.ne.0) then
                        m_err = 'INPS: Number in IY( ) is Wrong.'
                        ErrCha = ''
                        ErrID = 'L:20202/R:inclp/F:a-main1.f'
                        goto 999
                     end if
                     if(nint(rnm).ne.1.and.nint(rnm).ne.2.and.
     &                  nint(rnm).ne.3.and.nint(rnm).ne.0) then
                        m_err = 'INPS: Number in IY() is 1(Bottom), '//
     &                          '2(Center), 3(Top), or 0(BaseLine).'
                        ErrCha = ''
                        ErrID = 'L:20210/R:inclp/F:a-main1.f'
                        goto 999
                     end if
                     iyps(inps)=nint(rnm)
                     ic=ic+1

               end if

               goto 300

*-----------------------------------------------------------------------

  400       continue

            close(jps)

            return

*-----------------------------------------------------------------------

  999 continue

            l_err = ill(jsn)
            k_err = jsn
            ierr = 1

      return
      end


************************************************************************
*                                                                      *
      subroutine seqlin(dum,in,jsi,jsn,ill,ilf,ierr)
*                                                                      *
*                                                                      *
*      PURPOSE  : READ SEQUENTIAL LINE                                 *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character duml(ichrl)*1

      character m_err*200
      common /error/ m_err, l_err, k_err

      dimension ill(0:9), ilf(0:9)

      character yen*1
      character tub*1
      tub = char(9)
      yen = char(92)

*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------

                  do 407 i = 1, ichrl

                     if( i .le. icolm ) then

                        duml(i) = dum(i)

                     else

                        duml(i) = ' '

                     end if

  407             continue

                  do 101 i = icolm, in, -1

                     if( dum(i) .ne. ' ' .and. dum(i) .ne. tub )
     &               goto 111

  101             continue
  111             continue

*-----------------------------------------------------------------------

            ict = in - 1

         if( dum(i) .eq. yen ) then

  403          icf = i - 1

                  do 402 i = 1, icf - in + 1

                     duml( ict + i ) = dum( in + i - 1 )

  402             continue

                     ict = ict + icf - in + 1

*-----------------------------------------------------------------------

               read(jsi,'(10000a1)',iostat = ios ) (dum(ic),ic=1,icolm)
               if( ios .eq. -1 ) goto 406

*-----------------------------------------------------------------------

                     ill(jsn)  = ill(jsn) + 1

                     if( ill(jsn) .gt. ilf(jsn) ) then

                        m_err ='The Description of Secential Line '//
     &                         'is Wrong.'
                        ErrCha = ''
                        ErrID = 'L:20330/R:seqlin/F:a-main1.f'
                        l_err = ill(jsn)
                        k_err = jsn
                        goto 999

                     end if

                  in  = 1

                  do 404 i = icolm, 1, -1

                     if( dum(i) .ne. ' ' .and. dum(i) .ne. tub )
     &               goto 405

  404             continue
  405             continue

               if( dum(i) .eq. yen ) goto 403


                  icf = i

                  do 408 i = 1, icf - in + 1

                     if( ict + i .gt. ichrl ) then

                        m_err ='Length of This Sequential Line '//
     &                         'is Too Long.'
                        ErrCha = ''
                        ErrID = 'L:20359/R:seqlin/F:a-main1.f'
                        l_err = ill(jsn)
                        k_err = jsn
                        goto 999

                     end if

                     duml( ict + i ) = dum( in + i - 1 )

  408             continue

  406          continue

         end if

*-----------------------------------------------------------------------

                  do 409 i = 1, ichrl

                        dum(i) = duml(i)

  409             continue

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

  999 continue

      ierr = 1

      return
      end


************************************************************************
*                                                                      *
      subroutine wfrmg(jhf,idbg,ibfon,clms,imcm,wid,hgt,mcm,
     &                 clal,clmo,ifon,idat)
*                                                                      *
*      PURPOSE  :  WRITE FRAME MESSAGES                                *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character mcm(6,ichrl)*1
      dimension imcm(6)

      character mcm0(18)*1
      character mcm1(26)*1
      character mcm2(12)*1
      character mcm3(19)*1
      character mcm5(25)*1
      character mcm6(19)*1

      character ccmt(ichrl)*1

      data mcm0/ '{','@','i','t',' ','F','i','l','e',' ','=',' ','@',
     &           'f','i','l','e','}'/
      data mcm1/ '{','@','i','t',' ','D','a','t','e',' ','=',' ','@',
     &           't','i','m','e',' ',' ','@','t','o','d','a','y','}'/
      data mcm2/ 'F','i','l','e',' ','=',' ','@','f','i','l','e'/
      data mcm3/ 'D','a','t','e',' ','=',' ','@','t','i','m','e',' ',
     &           '@','t','o','d','a','y'/
      data mcm5/ '{','@','i','t',' ','D','a','t','e',' ','=',' ','@',
     &           't','o','d','a','y',' ','@','t','i','m','e','}'/
      data mcm6/ 'D','a','t','e',' ','=',' ','@','t','o','d','a','y',
     &           ' ','@','t','i','m','e'/

      character,allocatable:: chag(:,:)

*-----------------------------------------------------------------------

      character yen*1

      dimension clal(3), clms(3)

*-----------------------------------------------------------------------

      allocate(chag(inig,0:ichrl))

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------

            yen = char(92)

*-----------------------------------------------------------------------

         mcm0( 2) = yen
         mcm0(13) = yen

         mcm1( 2) = yen
         mcm1(13) = yen
         mcm1(20) = yen

         mcm2( 8) = yen

         mcm3( 8) = yen
         mcm3(14) = yen

         mcm5( 2) = yen
         mcm5(13) = yen
         mcm5(20) = yen

         mcm6( 8) = yen
         mcm6(15) = yen

*-----------------------------------------------------------------------
*     FRAME MESSAGE
*-----------------------------------------------------------------------

*        Default Font is IFMSG ( 10 pt )
*        for the massage font
*        Corner Margin is 3 pt

*-----------------------------------------------------------------------

               ifmsg = ibfon
               fsmsg = 10.0
               cmarg = 3.0

*-----------------------------------------------------------------------

            if( idbg .eq. 1 )
     &      write(jhf,'(/''%'',71(''-''),/
     &                   ''% Frame Message''
     &                 ,/''%'',71(''-'')/)')

*-----------------------------------------------------------------------

               if( clms(1) .lt. -r0max ) clms(1) = -2.0

*-----------------------------------------------------------------------

         do 501 imc = 1, 6

            if( imcm(imc) .eq. -2 .or. imcm(imc) .gt. 0 ) then

                     ixc = imc - imc / 4 * 3
                     iyc = abs( imc / 4 - 1 ) * 2 + 1

               if( imc .eq. 1 ) then

                     xpo = cmarg
                     ypo = hgt - cmarg

               else if( imc .eq. 2 ) then

                     xpo = wid / 2.0
                     ypo = hgt - cmarg

               else if( imc .eq. 3 ) then

                     xpo = wid - cmarg
                     ypo = hgt - cmarg

               else if( imc .eq. 4 ) then

                     xpo = cmarg
                     ypo = cmarg

               else if( imc .eq. 5 ) then

                     xpo = wid / 2.0
                     ypo = cmarg

               else if( imc .eq. 6 ) then

                     xpo = wid - cmarg
                     ypo = cmarg

               end if


               if( imcm(imc) .eq. -2 ) then

                  if( imc .eq. 1 ) then

                     if( ifmsg .ne. 9 ) then

                           istl = 18

                           do 300 i = 1, istl
                              ccmt(i) = mcm0(i)
  300                      continue

                     else

                           istl = 12

                           do 302 i = 1, istl
                              ccmt(i) = mcm2(i)
  302                      continue

                     end if

                  else if( imc .eq. 3 ) then

                     if( ifmsg .ne. 9 ) then

                        if( idat .eq. 0 ) then

                           istl = 26

                           do 301 i = 1, istl
                              ccmt(i) = mcm1(i)
  301                      continue

                        else

                           istl = 25

                           do 305 i = 1, istl
                              ccmt(i) = mcm5(i)
  305                      continue

                        end if

                     else

                        if( idat .eq. 0 ) then

                           istl = 19

                           do 303 i = 1, istl
                              ccmt(i) = mcm3(i)
  303                      continue

                        else

                           istl = 19

                           do 306 i = 1, istl
                              ccmt(i) = mcm6(i)
  306                      continue

                        end if

                     end if

                  end if

                     call wtext(jhf,chag,idbg,clal,clmo,0,xpo,ypo,
     &                          ccmt,istl,ixc,iyc,0,0.d0,ifmsg,ifon,
     &                          0,indt,fsmsg,clms)

               else if( imcm(imc) .gt. 0 ) then

                           do 304 i = 1, imcm(imc)
                              ccmt(i) = mcm(imc,i)
  304                      continue

                  call wtext(jhf,chag,idbg,clal,clmo,0,xpo,ypo,
     &                       ccmt,imcm(imc),ixc,iyc,0,0.d0,ifmsg,ifon,
     &                       0,indt,fsmsg,clms)

               end if

            end if

  501    continue

*-----------------------------------------------------------------------

      deallocate(chag)
      return
      end

************************************************************************
*                                                                      *
      subroutine wftat(jhf,idbg,noaf,dwaf,afac,clax,dwat,
     &                 itic,jxtic,jytic,nxlta,nxsta,nylta,nysta,
     &                 idecx,idecy,ixexp,iyexp,
cKN 2024/01/24
     &                 txl,txr,tyd,tyu,noxn,noyn,
     &                 ifon,itfon,clnm,ixlog,iylog,
     &                 ixnum,iynum,nxtch,nytch,xcta,ycta,
     &                 xlta,xsta,ylta,ysta,tmup,
     &                 xtt,ixtt,ytt,iytt,
     &                 fstl,fstc,ibfon,cltx,atxs,ix,iy,
     &                 ixtxt,xtxp,iytxt,ytxp,noxt,noyt)
*                                                                      *
*      purpose  : write axis frame, tics and axis text                 *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      common /con/  cm, dd
      common /frm/  xal, yal

      common /fdfn/ ifd(0:13),ikan1,ikan2,ikan3,ikan4,
     &              idsm(9), idln(6), idxe, idye, ibd1, ibd2, ibd3,
     &              ibd4, ibd5, ibd6, ibd7, ibd8, ibd9, ibd10, ibd11

      common /suf/  ftht, ftss, ftks, ftsv,
     &              yplu, ypld, ypku, ypkd, yplh, ypls

      common /wtval1/ strl0, strh0, strb0

      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb

      dimension xlta(numtic), xsta(numtic), ylta(numtic), ysta(numtic)
      dimension nxtch(numtic), nytch(numtic)
      character xcta(numtic,30)*1,ycta(numtic,30)*1

      character chaf(ichrl)*1
      character xtt(ichrl)*1,ytt(ichrl)*1
      character,allocatable:: chag(:,:)

      dimension clax(3), cltx(3), clnm(3)
      dimension clal(3)

*-----------------------------------------------------------------------
cKN 2024/01/24

      character dmms*30
      character rdum*100
      character d4*4
      character d8*8
      character d13*13
      dimension rnmm(30)
      character xctal(numtic,30)*1, yctal(numtic,30)*1

*-----------------------------------------------------------------------

      character xtn1(6)*1
      data xtn1/ ')',' ','x','t','n','1'/
      character xtn2(6)*1
      data xtn2/ ')',' ','x','t','n','2'/
      character xtn3(6)*1
      data xtn3/ ')',' ','x','t','n','3'/
      character xtn4(6)*1
      data xtn4/ ')',' ','x','t','n','4'/

      character ytn1(6)*1
      data ytn1/ ')',' ','y','t','n','1'/
      character ytn2(6)*1
      data ytn2/ ')',' ','y','t','n','2'/
      character ytn3(6)*1
      data ytn3/ ')',' ','y','t','n','3'/
      character ytn4(6)*1
      data ytn4/ ')',' ','y','t','n','4'/

*-----------------------------------------------------------------------

      allocate(chag(inig,0:ichrl))

*-----------------------------------------------------------------------

               xap = xal * cm
               yap = yal * cm

*-----------------------------------------------------------------------
*     WRITE AXIS FRAME, TICS AND AXIS TEXT
*-----------------------------------------------------------------------

*           Default Width of Axis Frame is DWAF = 5 dd

*=======================================================================

         if( noaf .eq. 1 ) then

               waxis = dwaf * afac

            if( idbg .eq. 1 )
     &      write(jhf,'(/''%'',71(''-''),/
     &                   ''% Axis Frame''
     &                 ,/''%'',71(''-'')/)')

               if( clax(1) .lt. -r0max ) clax(1) = -2.0

               write(jhf,'(3f7.3,'' sc sd0 '',g14.5,
     &                   '' dd lw ax st'')') clax, waxis

                  wwa = waxis * dd

                  call bbox(1,0,-wwa,-wwa,0.d0)
                  call bbox(1,0,xap+wwa,-wwa,0.d0)
                  call bbox(1,0,-wwa,yap+wwa,0.d0)
                  call bbox(1,0,xap+wwa,yap+wwa,0.d0)


         end if

*-----------------------------------------------------------------------
*        Axis Tics

*           WTICS : Default Width of Axis Tics is DWAT = 4 dd

*           TLL   : Default Length of Long Tics is 0.423 cm
*           TSL   : Default Length of Short Tics is 0.618 * TLL

*-----------------------------------------------------------------------

               wtics = dwat * afac
               tll   = 0.423 * sqrt( afac ) * cm
               tsl   = tll * 0.618

*-----------------------------------------------------------------------

      if( itic .ne. 0 .and.
     &  ( nxlta .gt. 0 .or. nxsta .gt. 0 .or.
     &    nylta .gt. 0 .or. nysta .gt. 0 ) ) then

*-----------------------------------------------------------------------

            if( idbg .eq. 1 )
     &      write(jhf,'(/''%'',71(''-''),/
     &                   ''% Axis Tics''
     &                 ,/''%'',71(''-'')/)')

               write(jhf,'(g14.5,'' dd lw sd0'')') wtics

               if( idbg .eq. 1 )  write(jhf,'()')

            if( nxlta .gt. 0 .or. nylta .gt. 0 ) then

                  write(jhf,'(''/tll '',g14.5,'' N'')') tll

            end if

            if( nxlta .gt. 0 .and. jxtic .ne. 0 ) then

               if( itic .eq. 1 ) then

                  if( jxtic .eq. 2 ) then

                     write(jhf,'(
     &               ''/xtl {/xp S x N xp 0 M xp tll L st '',
     &               ''xp yal tll sub M xp yal L st} N''
     &               )')

                  else if( jxtic .eq.  1 ) then

                     write(jhf,'(
     &               ''/xtl {/xp S x N xp 0 M xp tll L st } N'')')

                  else if( jxtic .eq. -1 ) then

                     write(jhf,'(
     &               ''/xtl {/xp S x N '',
     &               ''xp yal tll sub M xp yal L st} N''
     &               )')

                  end if

               else if( itic .eq. -1 ) then

                  if( jxtic .eq. 2 ) then

                     write(jhf,'(
     &               ''/xtl {/xp S x N xp 0 M xp tll neg L st '',
     &               ''xp yal M xp yal tll add L st} N''
     &               )')

                  else if( jxtic .eq.  1 ) then

                     write(jhf,'(
     &               ''/xtl {/xp S x N xp 0 M xp tll neg L st } N''
     &               )')

                  else if( jxtic .eq. -1 ) then

                     write(jhf,'(
     &               ''/xtl {/xp S x N '',
     &               ''xp yal M xp yal tll add L st} N''
     &               )')

                  end if

               end if

            end if

            if( nylta .gt. 0 .and. jytic .ne. 0 ) then

               if( itic .eq. 1 ) then

                  if( jytic .eq. 2 ) then

                     write(jhf,'(
     &               ''/ytl {/yp S y N 0 yp M tll yp L st '',
     &               ''xal tll sub yp M xal yp L st} N''
     &               )')

                  else if( jytic .eq.  1 ) then

                     write(jhf,'(
     &               ''/ytl {/yp S y N 0 yp M tll yp L st} N''
     &               )')

                  else if( jytic .eq. -1 ) then

                     write(jhf,'(
     &               ''/ytl {/yp S y N '',
     &               ''xal tll sub yp M xal yp L st} N''
     &               )')

                  end if

               else if( itic .eq. -1 ) then

                  if( jytic .eq. 2 ) then

                     write(jhf,'(
     &               ''/ytl {/yp S y N 0 yp M tll neg yp L st '',
     &               ''xal yp M xal tll add yp L st} N''
     &               )')

                  else if( jytic .eq.  1 ) then

                     write(jhf,'(
     &               ''/ytl {/yp S y N 0 yp M tll neg yp L st} N''
     &               )')

                  else if( jytic .eq. -1 ) then

                     write(jhf,'(
     &               ''/ytl {/yp S y N '',
     &               ''xal yp M xal tll add yp L st} N''
     &               )')

                  end if

               end if

            end if

               if( idbg .eq. 1 )  write(jhf,'()')

            if( nxsta .gt. 0 .or. nysta .gt. 0 ) then

                  write(jhf,'(''/tsl '',g14.5,'' N'')') tsl

            end if

            if( nxsta .gt. 0 .and. jxtic .ne. 0 ) then

               if( itic .eq. 1 ) then

                  if( jxtic .eq. 2 ) then

                     write(jhf,'(
     &               ''/xts {/xp S x N xp 0 M xp tsl L st '',
     &               ''xp yal tsl sub M xp yal L st} N''
     &               )')

                  else if( jxtic .eq.  1 ) then

                     write(jhf,'(
     &               ''/xts {/xp S x N xp 0 M xp tsl L st} N''
     &               )')

                  else if( jxtic .eq. -1 ) then

                     write(jhf,'(
     &               ''/xts {/xp S x N '',
     &               ''xp yal tsl sub M xp yal L st} N''
     &               )')

                  end if

               else if( itic .eq. -1 ) then

                  if( jxtic .eq. 2 ) then

                     write(jhf,'(
     &               ''/xts {/xp S x N xp 0 M xp tsl neg L st '',
     &               ''xp yal M xp yal tsl add L st} N''
     &               )')

                  else if( jxtic .eq.  1 ) then

                     write(jhf,'(
     &               ''/xts {/xp S x N xp 0 M xp tsl neg L st} N''
     &               )')

                  else if( jxtic .eq. -1 ) then

                     write(jhf,'(
     &               ''/xts {/xp S x N '',
     &               ''xp yal M xp yal tsl add L st} N''
     &               )')

                  end if

               end if

            end if

            if( nysta .gt. 0 .and. jytic .ne. 0 ) then

               if( itic .eq. 1 ) then

                  if( jytic .eq. 2 ) then

                     write(jhf,'(
     &               ''/yts {/yp S y N 0 yp M tsl yp L st '',
     &               ''xal tsl sub yp M xal yp L st} N''
     &               )')

                  else if( jytic .eq.  1 ) then

                     write(jhf,'(
     &               ''/yts {/yp S y N 0 yp M tsl yp L st} N''
     &               )')

                  else if( jytic .eq. -1 ) then

                     write(jhf,'(
     &               ''/yts {/yp S y N '',
     &               ''xal tsl sub yp M xal yp L st} N''
     &               )')

                  end if

               else if( itic .eq. -1 ) then

                  if( jytic .eq. 2 ) then

                     write(jhf,'(
     &               ''/yts {/yp S y N 0 yp M tsl neg yp L st '',
     &               ''xal yp M xal tsl add yp L st} N''
     &               )')

                  else if( jytic .eq.  1 ) then

                     write(jhf,'(
     &               ''/yts {/yp S y N 0 yp M tsl neg yp L st} N''
     &               )')

                  else if( jytic .eq. -1 ) then

                     write(jhf,'(
     &               ''/yts {/yp S y N '',
     &               ''xal yp M xal tsl add yp L st} N''
     &               )')

                  end if

               end if

            end if

               if( idbg .eq. 1 )  write(jhf,'()')

*-----------------------------------------------------------------------

         if( nxlta .gt. 0 .and. jxtic .ne. 0 ) then

               write(jhf,'(4(g14.5,'' xtl''))') (xlta(n),n=1,nxlta)

               if( idbg .eq. 1 )  write(jhf,'()')

         end if

         if( nxsta .gt. 0 .and. jxtic .ne. 0 ) then

               write(jhf,'(4(g14.5,'' xts''))') (xsta(n),n=1,nxsta)

               if( idbg .eq. 1 )  write(jhf,'()')

         end if

*-----------------------------------------------------------------------

         if( nylta .gt. 0 .and. jytic .ne. 0 ) then

               write(jhf,'(4(g14.5,'' ytl''))') (ylta(n),n=1,nylta)

               if( idbg .eq. 1 )  write(jhf,'()')

         end if

         if( nysta .gt. 0 .and. jytic .ne. 0 ) then

               write(jhf,'(4(g14.5,'' yts''))') (ysta(n),n=1,nysta)

               if( idbg .eq. 1 )  write(jhf,'()')

         end if

*-----------------------------------------------------------------------



*-----------------------------------------------------------------------

            if( itic .eq. -1 .and.
     &        ( jxtic .ne. 0 .or. jytic .ne. 0 .or.
     &          nxlta .gt. 0 .or. nxsta .gt. 0 .or.
     &          nylta .gt. 0 .or. nysta .gt. 0 ) ) then

               if( nxlta .gt. 0 .and. jxtic .ne. 0 ) then

                  call bbox(1,0,0.d0,-tll,0.d0)
                  call bbox(1,0,0.d0,yap+tll,0.d0)
                  call bbox(1,0,xap,-tll,0.d0)
                  call bbox(1,0,xap,yap+tll,0.d0)

               else if( nxsta .gt. 0 .and. jxtic .ne. 0 ) then

                  call bbox(1,0,0.d0,-tsl,0.d0)
                  call bbox(1,0,0.d0,yap+tsl,0.d0)
                  call bbox(1,0,xap,-tsl,0.d0)
                  call bbox(1,0,xap,yap+tsl,0.d0)

               end if

               if( nylta .gt. 0 .and. jytic .ne. 0 ) then

                  call bbox(1,0,-tll,0.d0,0.d0)
                  call bbox(1,0,-tll,yap,0.d0)
                  call bbox(1,0,xap+tll,0.d0,0.d0)
                  call bbox(1,0,xap+tll,yap,0.d0)

               else if( nysta .gt. 0 .and. jytic .ne. 0 ) then

                  call bbox(1,0,-tsl,0.d0,0.d0)
                  call bbox(1,0,-tsl,yap,0.d0)
                  call bbox(1,0,xap+tsl,0.d0,0.d0)
                  call bbox(1,0,xap+tsl,yap,0.d0)

               end if

            end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*     Write Tics Number
*-----------------------------------------------------------------------

*           Default Tic Font Size is 24 pt

*           TXYV = Margin of Tic number from Axis     : 0.55 * TFS
*           TXLV = Margin of X number (LOG) from Axis : 0.75 * TFS

*           TMUP = Position Up of Minus               : 0.08 * FNS

*           xtmn is Minimum X-Value of Y-Tic Number
*           xtma is Maximum X-Value of Y-Tic Number

*-----------------------------------------------------------------------
*           FONT FOR TIC NUMBER
*-----------------------------------------------------------------------

*              ITFON : TIC FONT
*              IFTMN : MINUS FONT (Symbol Fixed )

*-----------------------------------------------------------------------

            if( idbg .eq. 1 )
     &      write(jhf,'(/''%'',71(''-''),/
     &                   ''% Write Tics Number''
     &                 ,/''%'',71(''-'')/)')

*-----------------------------------------------------------------------
*     INITIAL VALUES OF MARGIN
*-----------------------------------------------------------------------

               txl  = 0.0
               txr  = xal * cm
               tyd  = 0.0
               tyu  = yal * cm
               xtmn = txl
               xtma = txr

*-----------------------------------------------------------------------

            if( itic .lt. 0 ) then

               if( nxlta .gt. 0 .and. jxtic .ne. 0 ) then

                     txl = txl - tll
                     txr = txr + tll
                     xtmn = txl
                     xtma = txr

               else if( nxsta .gt. 0 .and. jxtic .ne. 0 ) then

                     txl = txl - tsl
                     txr = txr + tsl
                     xtmn = txl
                     xtma = txr

               end if

               if( nylta .gt. 0 .and. jytic .ne. 0 ) then

                     tyd = tyd - tll
                     tyu = tyu + tll

               else if( nysta .gt. 0 .and. jytic .ne. 0 ) then

                     tyd = tyd - tsl
                     tyu = tyu + tsl

               end if

            end if

*-----------------------------------------------------------------------
*     CONSTANTS AND COLOR OF NUMBER
*-----------------------------------------------------------------------

         if( noxn .eq. 1 .or. noyn .eq. 1 ) then

            tfs = fstc * afac

            txyv = 0.55
            txlv = 0.75


               if( ifd(itfon) .eq. 0 )
     &         call ftdfn(jhf,ifon,itfon,idbg)

               write(jhf,'(''/fts '',g14.5,'' N /ftsu '',g14.5,
     &         '' N /ftn001 ft'',I3.3,'' N''
     &         )') tfs, tfs * ftss, itfon

               if( clnm(1) .lt. -r0max ) clnm(1) = -2.0

               write(jhf,'(3f7.3,'' sc'')') clnm

         end if

*-----------------------------------------------------------------------
*     DEFINE OF '10' FOR LOG AXIS
*-----------------------------------------------------------------------

            if( ( noxn .eq. 1 .and. ixlog .eq. 1 ) .or.
     &          ( noyn .eq. 1 .and. iylog .eq. 1 ) ) then

                  icn = 2
                  chaf(1) = '1'
                  chaf(2) = '0'

                  call strwhb(icn,chaf,tfs,itfon,str10w,str10h,str10b,
     &                        ifon)

                  ftsu = tfs * ftss
                  ypup = tfs * yplu

            end if

*-----------------------------------------------------------------------
*     WRITE X-NUMBER
*-----------------------------------------------------------------------
*     POSITION OF X-NUMBER DETERMINED WITH MARGIN FROM X-AXIS
*-----------------------------------------------------------------------
cKN 2024/01/24

            ieex = -100

      if( noxn .eq. 1 .and. nxlta .gt. 0 ) then

*-----------------------------------------------------------------------
*     x-linear and exp number
*-----------------------------------------------------------------------
cKN 2024/01/24

      if( ixlog .eq. 0 .and. ixexp .eq. 1 ) then

               wmax = 1.d-50

               do iww = 1, nxlta

                 do kw = 1, nxtch(iww)
                    dmms(kw:kw) = xcta(iww,kw)
                 end do

                    read(dmms(1:nxtch(iww)),*) rnmm(iww)

                 if( abs(rnmm(iww)) .gt. wmax ) then
                    wmax = abs(rnmm(iww))
                    iwmax = iww
                    if( rnmm(iww) .lt. 0.0 ) iwmax = -iwmax
                 end if

               end do

*-----------------------------------------------------------------------

         if( idecx .gt. 0 ) then

            if( wmax .lt. 1.e-3 ) then

                  write(rdum,'(1pe11.4e3)') wmax

                  d4=rdum(8:8)//rdum(9:9)//rdum(10:10)//rdum(11:11)
                  read(d4,'(i4)') ieex

                  idecxl = idecx + ieex

*-----------------------------------------------------------------------

            do iww = 1, nxlta

*-----------------------------------------------------------------------

               if( rnmm(iww) .ge. 0.d0 ) then

*-----------------------------------------------------------------------

                     do kw = 1, idecxl + 1

                        if( kw .eq. 1 ) then
                           xctal(iww,kw)   = xcta(iww,kw-ieex+1)
                           xctal(iww,kw+1) = '.'
                        else
                           xctal(iww,kw+1) = xcta(iww,kw-ieex+1)
                        end if

                     end do

                       if( idecxl .eq. 0 ) then
                           xctal(iww,3) = '0'
                       end if

                     nxtch(iww) = idecxl + 3

*-----------------------------------------------------------------------

               else if( rnmm(iww) .lt. 0.d0 ) then

                     do kw = 1, idecxl + 1

                        if( kw .eq. 1 ) then
                           xctal(iww,kw)   = xcta(iww,1)
                           xctal(iww,kw+1) = xcta(iww,kw-ieex+2)
                           xctal(iww,kw+2) = '.'
                        else
                           xctal(iww,kw+1) = xcta(iww,kw-ieex+2)
                        end if

                     end do

                       if( idecxl .eq. 0 ) then
                           xctal(iww,4) = '0'
                       end if

                     nxtch(iww) = idecxl + 4

*-----------------------------------------------------------------------

               end if

*-----------------------------------------------------------------------

            end do

*-----------------------------------------------------------------------

                  idecx = idecxl
                  if( idecxl .eq. 0 ) idecx = 1

               do iww = 1, nxlta
                  do kw = 1, nxtch(iww)
                     xcta(iww,kw) = xctal(iww,kw)
                  end do
               end do

cKN 2024/01/24 w
c               do iww = 1, nxlta
c                  write(*,*) (xcta(iww,j),j=1,nxtch(iww)),
c     &                        idecx, nxtch(iww), ieex
c               end do
cKN 2024/01/24 w

         end if

*-----------------------------------------------------------------------

         else if( idecx .lt. 0 ) then

         if( wmax .ge. 1.e+4 ) then

                  write(rdum,'(1pe11.4e3)') wmax

                  d4=rdum(8:8)//rdum(9:9)//rdum(10:10)//rdum(11:11)
                  read(d4,'(i4)') ieex

*-----------------------------------------------------------------------

            do iww = 1, nxlta

*-----------------------------------------------------------------------

               if( rnmm(iww) .ge. 0.d0 ) then

*-----------------------------------------------------------------------

                     icg = ieex - nxtch(iww)

                  if( icg .lt. 0 ) then

                        xctal(iww,1) = xcta(iww,1)
                        xctal(iww,2) = '.'

                     do ikk = 2, nxtch(iww)
                        xctal(iww,ikk+1) = xcta(iww,ikk)
                     end do

                        nxtch(iww) = nxtch(iww) + 1

                  else if( icg .ge. 0 ) then

                        xctal(iww,1) = '0'
                        xctal(iww,2) = '.'

                        ikkf = 2

                     do ikk = 1, icg
                        ikkf = ikkf + 1
                        xctal(iww,ikkf) = '0'
                     end do

                     do ikk = 1, nxtch(iww)
                        ikkf = ikkf + 1
                        xctal(iww,ikkf) = xcta(iww,ikk)
                     end do

                        nxtch(iww) = ikkf

                  end if

*-----------------------------------------------------------------------

               else if( rnmm(iww) .lt. 0.d0 ) then

*-----------------------------------------------------------------------

                     icg = ieex - nxtch(iww) + 1

                  if( icg .lt. 0 ) then

                        xctal(iww,1) = xcta(iww,1)
                        xctal(iww,2) = xcta(iww,2)
                        xctal(iww,3) = '.'

                     do ikk = 3, nxtch(iww)
                        xctal(iww,ikk+1) = xcta(iww,ikk)
                     end do

                        nxtch(iww) = nxtch(iww) + 1

                  else if( icg .ge. 0 ) then

                        xctal(iww,1) = xcta(iww,1)
                        xctal(iww,2) = '0'
                        xctal(iww,3) = '.'

                        ikkf = 3

                     do ikk = 1, icg
                        ikkf = ikkf + 1
                        xctal(iww,ikkf) = '0'
                     end do

                     do ikk = 2, nxtch(iww)
                        ikkf = ikkf + 1
                        xctal(iww,ikkf) = xcta(iww,ikk)
                     end do

                        nxtch(iww) = ikkf

                  end if

*-----------------------------------------------------------------------

               end if

*-----------------------------------------------------------------------

            end do

*-----------------------------------------------------------------------

                     idecm = 100

               do iww = 1, nxlta
                     nzeo = -1
                  do kw = nxtch(iww), 1, -1
                     nzeo = nzeo + 1
                     if( xctal(iww,kw) .eq. '.' ) then
                        nzeo = nzeo - 1
                        goto 100
                     end if
                     if( xctal(iww,kw) .ne. '0' ) goto 100
                  end do
  100                if( nzeo .lt. idecm ) idecm = nzeo
               end do

               do iww = 1, nxlta
                     nxtch(iww) = nxtch(iww) - idecm
                  do kw = 1, nxtch(iww)
                     xcta(iww,kw) = xctal(iww,kw)
                  end do
               end do

                  if( iwmax .gt. 0 ) then
                     idecx = nxtch(iwmax) - 2
                  else
                     idecx = nxtch(-iwmax) - 3
                  end if

cKN 2024/01/24 w
c               do iww = 1, nxlta
c                  write(*,*) (xcta(iww,j),j=1,nxtch(iww)),
c     &                        idecx, nxtch(iww), ieex
c               end do
cKN 2024/01/24 w
*-----------------------------------------------------------------------

         end if

         end if

      end if

cKN 2024/01/24
*-----------------------------------------------------------------------

               if( ixnum .eq. 1 ) then

                  if( ixlog .eq. 0 ) then

                     tyd = tyd - txyv * tfs

                  else

                     tyd = tyd - txlv * tfs

                  end if

               else if( ixnum .eq. -1 ) then

                     tyu = tyu + txyv * tfs

               else if( ixnum .eq. 2 ) then

                  if( ixlog .eq. 0 ) then

                     tyd = tyd - txyv * tfs

                  else

                     tyd = tyd - txlv * tfs

                  end if

                     tyu = tyu + txyv * tfs

               end if

*-----------------------------------------------------------------------
*        macro for tic numbers
*-----------------------------------------------------------------------

*           ixtn1 ; linear without minus
*           ixtn2 ; linear with munus
*           ixtn3 ; log without minus
*           ixtn4 ; log with minus

*-----------------------------------------------------------------------

            ixtn1 = 0
            ixtn2 = 0
            ixtn3 = 0
            ixtn4 = 0

*-----------------------------------------------------------------------

         do 710 n = 1, nxlta

*-----------------------------------------------------------------------
            if( ixlog .eq. 0 ) then
*-----------------------------------------------------------------------

               if( xcta(n,1) .ne. '-' ) then

                  if( ixtn1 .eq. 0 ) then

                        if( idbg .eq. 1 )  write(jhf,'()')

                     if( ixnum .eq. 1 ) then

                        yp1 = tyd - tfs * ftht

                        write(jhf,'(
     &                  ''/xtn1 {/str001 S N /xps S N str001 xps '',
     &                  g14.5,'' ftn001 fts stsw } N''
     &                  )') yp1

                     else if( ixnum .eq. -1 ) then

                        yp2 = tyu

                        write(jhf,'(
     &                  ''/xtn1 {/str001 S N /xps S N str001 xps '',
     &                  g14.5,'' ftn001 fts stsw } N''
     &                  )') yp2

                     else if( ixnum .eq. 2 ) then

                        yp1 = tyd - tfs * ftht
                        yp2 = tyu

                        write(jhf,'(
     &                  ''/xtn1 {/str001 S N /xps S N''
     &                  '' str001 xps '',g14.5,
     &                  '' ftn001 fts stsw''/
     &                  '' str001 xps '',g14.5,'' stsm } N''
     &                  )') yp1, yp2

                     end if

                        ixtn1 = 1

                        if( idbg .eq. 1 )  write(jhf,'()')

                  end if


                     icn = nxtch(n)

                     do 810 i = 1, nxtch(n)

                        chaf(i) = xcta(n,i)

  810                continue

                     call strwhb(icn,chaf,tfs,itfon,strwf,strhf,strbf,
     &                           ifon)

                     xp1 = xlta(n) * xal * cm - strwf / 2.0
                     xp2 = xp1 + strwf

                     write(jhf,'(g14.5,'' ('',40a1)')
     &               xp1, ( chaf(i), i = 1, icn ), xtn1

               else if( xcta(n,1) .eq. '-' ) then

                  if( ixtn2 .eq. 0 ) then

                        if( idbg .eq. 1 )  write(jhf,'()')

                        icn = 1
                        chaf(1) = '-'
                        call strwhb(icn,chaf,tfs,4,strm,strhf,strbf,
     &                              ifon)

                        rmup = tfs * tmup

                     if( ixnum .eq. 1 ) then

                        yp1 = tyd - tfs * ftht

                        write(jhf,'(
     &                  ''/xtn2 {/str001 S N /xps S N (-) xps '',g14.5,
     &                  '' /Symbol findfont fts stsw'',/
     &                  '' str001 xps '',g14.5,'' add '',g14.5,
     &                  '' ftn001 fts stsw } N''
     &                  )') yp1 + rmup, strm, yp1

                     else if( ixnum .eq. -1 ) then

                        yp2 = tyu

                        write(jhf,'(
     &                  ''/xtn2 {/str001 S N /xps S N (-) xps '',g14.5,
     &                  '' /Symbol findfont fts stsw'',/
     &                  '' str001 xps '',g14.5,'' add '',g14.5,
     &                  '' ftn001 fts stsw } N''
     &                  )') yp2 + rmup, strm, yp2

                     else if( ixnum .eq. 2 ) then

                        yp1 = tyd - tfs * ftht
                        yp2 = tyu

                        write(jhf,'(
     &                  ''/xtn2 {/str001 S N /xps S N'',/
     &                  ''(-) xps '',g14.5,
     &                  '' /Symbol findfont fts stsw'',/
     &                  ''(-) xps '',g14.5,'' stsm'',/
     &                  '' str001 xps '',g14.5,'' add '',g14.5,
     &                  '' ftn001 fts stsw''/
     &                  '' str001 xps '',g14.5,'' add '',g14.5,
     &                  '' stsm } N''
     &                  )') yp1 + rmup, yp2 + rmup,
     &                      strm, yp1, strm, yp2

                     end if

                        ixtn2 = 1

                        if( idbg .eq. 1 )  write(jhf,'()')

                  end if


                     icn = nxtch(n) - 1

                     do 811 i = 1, nxtch(n) - 1

                        chaf(i) = xcta(n,i+1)

  811                continue

                     call strwhb(icn,chaf,tfs,itfon,strwf,strhf,strbf,
     &                           ifon)

                     xp1 = xlta(n) * xal * cm - ( strwf + strm ) / 2.0
                     xp2 = xp1 + strwf + strm

                     write(jhf,'(g14.5,'' ('',40a1)')
     &               xp1, ( chaf(i), i = 1, icn ), xtn2

               end if


                     if( ixnum .eq. 1 ) then

                        yq1 = yp1 + strbf

                     else if( ixnum .eq. -1 ) then

                        yq2 = yp2 + strhf

                     else if( ixnum .eq. 2 ) then

                        yq1 = yp1 + strbf
                        yq2 = yp2 + strhf

                     end if

*-----------------------------------------------------------------------
            else if( ixlog .ne. 0 ) then
*-----------------------------------------------------------------------

               if( xcta(n,1) .ne. '-' ) then

                  if( ixtn3 .eq. 0 ) then

                        if( idbg .eq. 1 )  write(jhf,'()')

                     if( ixnum .eq. 1 ) then

                        yp1  = tyd - tfs * ftht
                        yp1u = yp1 + ypup

                        write(jhf,'(
     &                  ''/xtn3 {/str001 S N /xps S N''/
     &                  ''(10) xps '',g14.5,'' ftn001 fts stsw'',/
     &                  ''str001 xps '',g14.5,'' add '',g14.5,
     &                  '' ftn001 ftsu stsw } N''
     &                  )') yp1, str10w, yp1u

                     else if( ixnum .eq. -1 ) then

                        yp2  = tyu
                        yp2u = tyu + ypup

                        write(jhf,'(
     &                  ''/xtn3 {/str001 S N /xps S N''/
     &                  ''(10) xps '',g14.5,'' ftn001 fts stsw'',/
     &                  ''str001 xps '',g14.5,'' add '',g14.5,
     &                  '' ftn001 ftsu stsw } N''
     &                  )') yp2, str10w, yp2u

                     else if( ixnum .eq. 2 ) then

                        yp1  = tyd - tfs * ftht
                        yp1u = yp1 + ypup
                        yp2  = tyu
                        yp2u = tyu + ypup

                        write(jhf,'(
     &                  ''/xtn3 {/str001 S N /xps S N''/
     &                  ''(10) xps '',g14.5,'' ftn001 fts stsw'',/
     &                  ''(10) xps '',g14.5,'' stsm'',/
     &                  ''str001 xps '',g14.5,'' add '',g14.5,
     &                  '' ftn001 ftsu stsw'',/
     &                  ''str001 xps '',g14.5,'' add '',g14.5,
     &                  '' stsm } N''
     &                  )') yp1, yp2, str10w, yp1u, str10w, yp2u

                     end if

                        ixtn3 = 1

                        if( idbg .eq. 1 )  write(jhf,'()')

                  end if

                     icn = nxtch(n)

                     do 812 i = 1, nxtch(n)

                        chaf(i) = xcta(n,i)

  812                continue

                     call strwhb(icn,chaf,ftsu,itfon,strwf,strhf,strbf,
     &                           ifon)

                     xp1 = xlta(n) * xal * cm
     &                   - ( strwf + str10w ) / 2.0
                     xp2 = xp1 + strwf + str10w

                     write(jhf,'(g14.5,'' ('',40a1)')
     &               xp1, ( chaf(i), i = 1, icn ), xtn3


               else

                  if( ixtn4 .eq. 0 ) then

                        if( idbg .eq. 1 )  write(jhf,'()')

                        icn = 1
                        chaf(1) = '-'
                        call strwhb(icn,chaf,ftsu,4,strm,strhf,strbf,
     &                              ifon)

                        rmup = ftsu * tmup

                     if( ixnum .eq. 1 ) then

                        yp1  = tyd - tfs * ftht
                        yp1u = yp1 + ypup

                        write(jhf,'(
     &                  ''/xtn4 {/str001 S N /xps S N''/
     &                  ''  (10) xps '',g14.5,'' ftn001 fts stsw'',/
     &                  ''  (-)  xps '',g14.5,'' add '',g14.5,
     &                  '' /Symbol findfont ftsu stsw'',/
     &                  ''str001 xps '',g14.5,'' add '',g14.5,
     &                  '' ftn001 ftsu stsw } N''
     &                  )') yp1, str10w, yp1u + rmup, str10w+strm, yp1u

                     else if( ixnum .eq. -1 ) then

                        yp2  = tyu
                        yp2u = tyu + ypup

                        write(jhf,'(
     &                  ''/xtn4 {/str001 S N /xps S N''/
     &                  ''  (10) xps '',g14.5,'' ftn001 fts stsw'',/
     &                  ''  (-)  xps '',g14.5,'' add '',g14.5,
     &                  '' /Symbol findfont ftsu stsw'',/
     &                  ''str001 xps '',g14.5,'' add '',g14.5,
     &                  '' ftn001 ftsu stsw } N''
     &                  )') yp2, str10w, yp2u + rmup, str10w+strm, yp2u

                     else if( ixnum .eq. 2 ) then

                        yp1  = tyd - tfs * ftht
                        yp1u = yp1 + ypup
                        yp2  = tyu
                        yp2u = tyu + ypup

                        write(jhf,'(
     &                  ''/xtn4 {/str001 S N /xps S N''/
     &                  ''  (10) xps '',g14.5,'' ftn001 fts stsw'',/
     &                  ''  (10) xps '',g14.5,'' stsm'',/
     &                  ''  (-)  xps '',g14.5,'' add '',g14.5,
     &                  '' /Symbol findfont ftsu stsw'',/
     &                  ''  (-)  xps '',g14.5,'' add '',g14.5,
     &                  '' stsm'',/
     &                  ''str001 xps '',g14.5,'' add '',g14.5,
     &                  '' ftn001 ftsu stsw'',/
     &                  ''str001 xps '',g14.5,'' add '',g14.5,
     &                  '' stsm } N''
     &                  )') YP1, YP2, STR10W, YP1U + RMUP,
     &                                str10w, yp2u + rmup,
     &                      str10w+strm, yp1u, str10w+strm, yp2u

                     end if

                        ixtn4 = 1

                        if( idbg .eq. 1 )  write(jhf,'()')

                  end if

                     icn = nxtch(n) - 1

                     do 813 i = 1, nxtch(n) - 1

                        chaf(i) = xcta(n,i+1)

  813                continue

                     call strwhb(icn,chaf,ftsu,itfon,strwf,strhf,strbf,
     &                           ifon)

                     xp1 = xlta(n) * xal * cm
     &                   - ( strwf + str10w + strm) / 2.0
                     xp2 = xp1 + strwf + str10w + strm

                     write(jhf,'(g14.5,'' ('',40a1)')
     &               xp1, ( chaf(i), i = 1, icn ), xtn4

               end if


                     if( ixnum .eq. 1 ) then

                        yq1 = yp1 + str10b

                     else if( ixnum .eq. -1 ) then

                        yq2 = yp2u + strhf

                     else if( ixnum .eq. 2 ) then

                        yq1 = yp1 + str10b
                        yq2 = yp2u + strhf

                     end if

*-----------------------------------------------------------------------
            end if
*-----------------------------------------------------------------------

                     if( ixnum .eq. 1 ) then

                        call bbox(1,0,xp1,yq1,0.d0)
                        call bbox(1,0,xp2,yq1,0.d0)

                     else if( ixnum .eq. -1 ) then

                        call bbox(1,0,xp1,yq2,0.d0)
                        call bbox(1,0,xp2,yq2,0.d0)

                     else if( ixnum .eq. 2 ) then

                        call bbox(1,0,xp1,yq1,0.d0)
                        call bbox(1,0,xp2,yq1,0.d0)
                        call bbox(1,0,xp1,yq2,0.d0)
                        call bbox(1,0,xp2,yq2,0.d0)

                     end if

  710    continue

*-----------------------------------------------------------------------

            if( ixnum .eq. 1 .or. ixnum .eq. 2 ) then

                  tyd = tyd - tfs * ftht

            end if

            if( ixnum .eq. -1 .or. ixnum .eq. 2 ) then

               if( ixlog .eq. 0 ) then

                  tyu = tyu + tfs * ftht

               else

                  tyu = tyu + tfs * 0.8

               end if

            end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     WRITE Y-NUMBER
*-----------------------------------------------------------------------
*     POSITION OF Y-NUMBER DETERMINED WITH MARGIN FROM Y-AXIS
*-----------------------------------------------------------------------
cKN 2024/01/24

            ieey = -100

      if( noyn .eq. 1 .and. nylta .gt. 0 ) then

*-----------------------------------------------------------------------
*     y-linear and exp number
*-----------------------------------------------------------------------
cKN 2024/01/24

      if( iylog .eq. 0 .and. iyexp .eq. 1 ) then

               wmax = 1.d-50

               do iww = 1, nylta

                 do kw = 1, nytch(iww)
                    dmms(kw:kw) = ycta(iww,kw)
                 end do

                    read(dmms(1:nytch(iww)),*) rnmm(iww)

                 if( abs(rnmm(iww)) .gt. wmax ) then
                    wmax = abs(rnmm(iww))
                    iwmax = iww
                    if( rnmm(iww) .lt. 0.0 ) iwmax = -iwmax
                 end if

               end do

*-----------------------------------------------------------------------

         if( idecy .gt. 0 ) then

            if( wmax .lt. 1.e-3 ) then

                  write(rdum,'(1pe11.4e3)') wmax

                  d4=rdum(8:8)//rdum(9:9)//rdum(10:10)//rdum(11:11)
                  read(d4,'(i4)') ieey

                  idecyl = idecy + ieey

*-----------------------------------------------------------------------

            do iww = 1, nylta

*-----------------------------------------------------------------------

               if( rnmm(iww) .ge. 0.d0 ) then

*-----------------------------------------------------------------------

                     do kw = 1, idecyl + 1

                        if( kw .eq. 1 ) then
                           yctal(iww,kw)   = ycta(iww,kw-ieey+1)
                           yctal(iww,kw+1) = '.'
                        else
                           yctal(iww,kw+1) = ycta(iww,kw-ieey+1)
                        end if

                     end do

                       if( idecyl .eq. 0 ) then
                           yctal(iww,3) = '0'
                       end if

                     nytch(iww) = idecyl + 3

*-----------------------------------------------------------------------

               else if( rnmm(iww) .lt. 0.d0 ) then

                     do kw = 1, idecyl + 1

                        if( kw .eq. 1 ) then
                           yctal(iww,kw)   = ycta(iww,1)
                           yctal(iww,kw+1) = ycta(iww,kw-ieey+2)
                           yctal(iww,kw+2) = '.'
                        else
                           yctal(iww,kw+1) = ycta(iww,kw-ieey+2)
                        end if

                     end do

                       if( idecyl .eq. 0 ) then
                           yctal(iww,4) = '0'
                       end if

                     nytch(iww) = idecyl + 4

*-----------------------------------------------------------------------

               end if

*-----------------------------------------------------------------------

            end do

*-----------------------------------------------------------------------

                  idecy = idecyl
                  if( idecyl .eq. 0 ) idecy = 1

               do iww = 1, nylta
                  do kw = 1, nytch(iww)
                     ycta(iww,kw) = yctal(iww,kw)
                  end do
               end do

cKN 2024/01/24 w
c               do iww = 1, nylta
c                  write(*,*) (ycta(iww,j),j=1,nytch(iww)),
c     &                        idecy, nytch(iww), ieey
c               end do
cKN 2024/01/24 w

         end if

*-----------------------------------------------------------------------

         else if( idecy .lt. 0 ) then

         if( wmax .ge. 1.e+4 ) then

                  write(rdum,'(1pe11.4e3)') wmax

                  d4=rdum(8:8)//rdum(9:9)//rdum(10:10)//rdum(11:11)
                  read(d4,'(i4)') ieey

*-----------------------------------------------------------------------

            do iww = 1, nylta

*-----------------------------------------------------------------------

               if( rnmm(iww) .ge. 0.d0 ) then

*-----------------------------------------------------------------------

                     icg = ieey - nytch(iww)

                  if( icg .lt. 0 ) then

                        yctal(iww,1) = ycta(iww,1)
                        yctal(iww,2) = '.'

                     do ikk = 2, nytch(iww)
                        yctal(iww,ikk+1) = ycta(iww,ikk)
                     end do

                        nytch(iww) = nytch(iww) + 1

                  else if( icg .ge. 0 ) then

                        yctal(iww,1) = '0'
                        yctal(iww,2) = '.'

                        ikkf = 2

                     do ikk = 1, icg
                        ikkf = ikkf + 1
                        yctal(iww,ikkf) = '0'
                     end do

                     do ikk = 1, nytch(iww)
                        ikkf = ikkf + 1
                        yctal(iww,ikkf) = ycta(iww,ikk)
                     end do

                        nytch(iww) = ikkf

                  end if

*-----------------------------------------------------------------------

               else if( rnmm(iww) .lt. 0.d0 ) then

*-----------------------------------------------------------------------

                     icg = ieey - nytch(iww) + 1

                  if( icg .lt. 0 ) then

                        yctal(iww,1) = ycta(iww,1)
                        yctal(iww,2) = ycta(iww,2)
                        yctal(iww,3) = '.'

                     do ikk = 3, nytch(iww)
                        yctal(iww,ikk+1) = ycta(iww,ikk)
                     end do

                        nytch(iww) = nytch(iww) + 1

                  else if( icg .ge. 0 ) then

                        yctal(iww,1) = ycta(iww,1)
                        yctal(iww,2) = '0'
                        yctal(iww,3) = '.'

                        ikkf = 3

                     do ikk = 1, icg
                        ikkf = ikkf + 1
                        yctal(iww,ikkf) = '0'
                     end do

                     do ikk = 2, nytch(iww)
                        ikkf = ikkf + 1
                        yctal(iww,ikkf) = ycta(iww,ikk)
                     end do

                        nytch(iww) = ikkf

                  end if

*-----------------------------------------------------------------------

               end if

*-----------------------------------------------------------------------

            end do

*-----------------------------------------------------------------------

                     idecm = 100

               do iww = 1, nylta
                     nzeo = -1
                  do kw = nytch(iww), 1, -1
                     nzeo = nzeo + 1
                     if( yctal(iww,kw) .eq. '.' ) then
                        nzeo = nzeo - 1
                        goto 200
                     end if
                     if( yctal(iww,kw) .ne. '0' ) goto 200
                  end do
  200                if( nzeo .lt. idecm ) idecm = nzeo
               end do

               do iww = 1, nylta
                     nytch(iww) = nytch(iww) - idecm
                  do kw = 1, nytch(iww)
                     ycta(iww,kw) = yctal(iww,kw)
                  end do
               end do

                  if( iwmax .gt. 0 ) then
                     idecy = nytch(iwmax) - 2
                  else
                     idecy = nytch(-iwmax) - 3
                  end if

cKN 2024/01/24 w
c               do iww = 1, nylta
c                  write(*,*) (ycta(iww,j),j=1,nytch(iww)),
c     &                        idecy, nytch(iww), ieey
c               end do
cKN 2024/01/24 w
*-----------------------------------------------------------------------

         end if

         end if

      end if

cKN 2024/01/24
*-----------------------------------------------------------------------

               txl = txl - tfs * txyv
               txr = txr + tfs * txyv

*-----------------------------------------------------------------------
*        macro for tic numbers
*-----------------------------------------------------------------------

*           iytn1   ; linear without minus
*           iytn2   ; linear with munus
*           iytn3   ; log without minus
*           iytn4   ; log with minus

*-----------------------------------------------------------------------

            iytn1 = 0
            iytn2 = 0

*-----------------------------------------------------------------------
         if( iylog .eq. 0 ) then
*-----------------------------------------------------------------------

            do 720 n = 1, nylta

*-----------------------------------------------------------------------

               if( ycta(n,1) .ne. '-' ) then

                  if( iytn1 .eq. 0 ) then

                        if( idbg .eq. 1 )  write(jhf,'()')

                     if( iynum .eq. 1 ) then

                        write(jhf,'(
     &                  ''/ytn1 {/str001 S N /yps S N /xps S N ''
     &                  ''str001 xps yps ftn001 fts stsw } N''
     &                  )')

                     else if( iynum .eq. -1 ) then

                        txrm = txr

                        if( iytn2 .ne. 0 ) txrm = txrm + strm

                        write(jhf,'(
     &                  ''/ytn1 {/str001 S N /yps S N ''
     &                  ''str001 '',g14.5,'' yps ftn001 fts stsw } N''
     &                  )') txrm

                     else if( iynum .eq. 2 ) then

                        txrm = txr

                        if( iytn2 .ne. 0 ) txrm = txrm + strm

                        write(jhf,'(
     &                  ''/ytn1 {/str001 S N /yps S N /xps S N ''
     &                  ''str001 xps yps ftn001 fts stsw'',/
     &                  ''str001 '',g14.5,'' yps stsm } N''
     &                  )') txrm

                     end if

                        iytn1 = 1

                        if( idbg .eq. 1 )  write(jhf,'()')

                  end if

                     icn = nytch(n)

                     do 814 i = 1, nytch(n)

                        chaf(i) = ycta(n,i)

  814                continue

                     call strwhb(icn,chaf,tfs,itfon,strwf,strhf,strbf,
     &                           ifon)

                        yp1 = ylta(n) * yal * cm - tfs * ftht / 2.0
                        yp2 = yp1 + strhf

                     if( iynum .eq. 1 ) then

                        xp1 = txl - strwf
                        xp2 = txl

                        write(jhf,'(2g14.5,'' ('',40a1)')
     &                  xp1, yp1, ( chaf(i), i = 1, icn ), ytn1

                     else if( iynum .eq. -1 ) then

                        xp1 = txr
                        xp2 = txr + strwf

                        write(jhf,'(g14.5,'' ('',40a1)')
     &                  yp1, ( chaf(i), i = 1, icn ), ytn1

                     else if( iynum .eq.  2 ) then

                        xp1 = txl - strwf
                        xp2 = txr + strwf

                        write(jhf,'(2g14.5,'' ('',40a1)')
     &                  xp1, yp1, ( chaf(i), i = 1, icn ), ytn1

                     end if

               else

                  if( iytn2 .eq. 0 ) then

                        if( idbg .eq. 1 )  write(jhf,'()')

                        icn = 1
                        chaf(1) = '-'
                        call strwhb(icn,chaf,tfs,4,strm,strhf,strbf,
     &                              ifon)

                        rmup = tfs * tmup

                     if( iynum .eq. 1 ) then

                        write(jhf,'(
     &                  ''/ytn2 {/str001 S N /yps S N /xps S N''/
     &                  ''(-) xps yps '',g14.5,
     &                  '' add /Symbol findfont'',
     &                  '' fts stsw''/
     &                  '' str001 xps '',g14.5,
     &                  '' add yps ftn001 fts stsw } N''
     &                  )') rmup, strm

                     else if( iynum .eq. -1 ) then

                        write(jhf,'(
     &                  ''/ytn2 {/str001 S N /yps S N'',/
     &                  ''(-) '',g14.5,'' yps '',g14.5,
     &                  '' add /Symbol findfont fts stsw''/
     &                  ''str001 '',g14.5,'' yps ftn001 fts stsw } N''
     &                  )') txr, rmup, txr + strm

                     else if( iynum .eq. 2 ) then

                        write(jhf,'(
     &                  ''/ytn2 {/str001 S N /yps S N /xps S N''/
     &                  ''(-) xps yps '',g14.5,
     &                  '' add /Symbol findfont'',
     &                  '' fts stsw''/
     &                  ''(-) '',g14.5,'' yps '',g14.5,
     &                  '' add stsm''/
     &                  ''str001 xps '',g14.5,
     &                  '' add yps ftn001 fts stsw'',/
     &                  ''str001 '',g14.5,'' yps ftn001 fts stsw } N''
     &                  )') rmup, txr, rmup, strm, txr + strm

                     end if

                        iytn2 = 1

                        if( idbg .eq. 1 )  write(jhf,'()')

                  end if

                     icn = nytch(n) - 1

                     do 815 i = 1, nytch(n) - 1

                        chaf(i) = ycta(n,i+1)

  815                continue

                     call strwhb(icn,chaf,tfs,itfon,strwf,strhf,strbf,
     &                           ifon)

                        yp1 = ylta(n) * yal * cm - tfs * ftht / 2.0
                        yp2 = yp1 + strhf


                     if( iynum .eq. 1 ) then

                        xp1 = txl - strwf - strm
                        xp2 = txl

                        write(jhf,'(2g14.5,'' ('',40a1)')
     &                  xp1, yp1, ( chaf(i), i = 1, icn ), ytn2

                     else if( iynum .eq. -1 ) then

                        xp1 = txr
                        xp2 = txr + strwf + strm

                        write(jhf,'(g14.5,'' ('',40a1)')
     &                  yp1, ( chaf(i), i = 1, icn ), ytn2

                     else if( iynum .eq.  2 ) then

                        xp1 = txl - strwf - strm
                        xp2 = txr + strwf + strm

                        write(jhf,'(2g14.5,'' ('',40a1)')
     &                  xp1, yp1, ( chaf(i), i = 1, icn ), ytn2

                     end if


               end if

*-----------------------------------------------------------------------

                        yp1 = yp1 + strbf

                     if( iynum .eq. 1 ) then

                        xtmn = min(xtmn,xp1)

                        call bbox(1,0,xp1,yp1,0.d0)
                        call bbox(1,0,xp1,yp2,0.d0)

                     else if( iynum .eq. -1 ) then

                        xtma = max(xtma,xp2)

                        call bbox(1,0,xp2,yp1,0.d0)
                        call bbox(1,0,xp2,yp2,0.d0)

                     else if( iynum .eq.  2 ) then

                        xtmn = min(xtmn,xp1)
                        xtma = max(xtma,xp2)

                        call bbox(1,0,xp1,yp1,0.d0)
                        call bbox(1,0,xp1,yp2,0.d0)
                        call bbox(1,0,xp2,yp1,0.d0)
                        call bbox(1,0,xp2,yp2,0.d0)

                     end if

  720       continue

*-----------------------------------------------------------------------
         end if
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        YLOG CASE
*-----------------------------------------------------------------------

         if( iylog .eq. 1 ) then

                     ytlm = 0.0

            do 740 n = 1, nylta

               if( ycta(n,1) .ne. '-' ) then

                     icn = nytch(n)

                     do 816 i = 1, nytch(n)

                        chaf(i) = ycta(n,i)

  816                continue

                     call strwhb(icn,chaf,ftsu,itfon,strwf,strhf,strbf,
     &                           ifon)

                        ytlm = max(ytlm,strwf)

               else if( ycta(n,1) .eq. '-' ) then

                        icn = 1
                        chaf(1) = '-'
                        call strwhb(icn,chaf,ftsu,4,strm,strhf,strbf,
     &                              ifon)

                        rmup = ftsu * tmup

                     icn = nytch(n) - 1

                     do 817 i = 1, nytch(n) - 1

                        chaf(i) = ycta(n,i+1)

  817                continue

                     call strwhb(icn,chaf,ftsu,itfon,strwf,strhf,strbf,
     &                           ifon)

                        ytlm = max(ytlm,strwf+strm)

               end if

  740       continue

*-----------------------------------------------------------------------

                  iytn3 = 0
                  iytn4 = 0

            do 730 n = 1, nylta

               if( ycta(n,1) .ne. '-' ) then

                  if( iytn3 .eq. 0 ) then

                        if( idbg .eq. 1 )  write(jhf,'()')

                     if( iynum .eq. 1 ) then

                        write(jhf,'(
     &                  ''/ytn3 {/str001 S N /yps S N'',/
     &                  ''(10) '',g14.5,'' yps ftn001 fts stsw'',/
     &                  ''str001 '',g14.5,'' yps '',g14.5,
     &                  '' add ftn001 ftsu stsw } N''
     &                  )') txl - ytlm - str10w, txl - ytlm, ypup

                     else if( iynum .eq. -1 ) then

                        write(jhf,'(
     &                  ''/ytn3 {/str001 S N /yps S N'',/
     &                  ''(10) '',g14.5,'' yps ftn001 fts stsw'',/
     &                  ''str001 '',g14.5,'' yps '',g14.5,
     &                  '' add ftn001 ftsu stsw } N''
     &                  )') txr, txr + str10w, ypup

                     else if( iynum .eq. 2 ) then

                        write(jhf,'(
     &                  ''/ytn3 {/str001 S N /yps S N'',/
     &                  ''(10) '',g14.5,'' yps ftn001 fts stsw'',/
     &                  ''(10) '',g14.5,'' yps stsm'',/
     &                  ''str001 '',g14.5,'' yps '',g14.5,
     &                  '' add ftn001 ftsu stsw'',/
     &                  ''str001 '',g14.5,'' yps '',g14.5,
     &                  '' add stsm } N''
     &                  )') txl - ytlm - str10w, txr,
     &                      txl - ytlm, ypup, txr + str10w, ypup

                     end if

                        iytn3 = 1

                        if( idbg .eq. 1 )  write(jhf,'()')

                  end if

                     icn = nytch(n)

                     do 818 i = 1, nytch(n)

                        chaf(i) = ycta(n,i)

  818                continue

                     call strwhb(icn,chaf,ftsu,itfon,strwf,strhf,strbf,
     &                           ifon)

                        yp1 = ylta(n) * yal * cm - tfs * ftht / 2.0
                        yp2 = yp1 + ypup + strhf

                     if( iynum .eq. 1 ) then

                        xp1 = txl - ytlm - str10w
                        xp2 = txl

                        write(jhf,'(g14.5,'' ('',40a1)')
     &                  yp1, ( chaf(i), i = 1, icn ), ytn3

                     else if( iynum .eq. -1 ) then

                        xp1 = txr
                        xp2 = txr + str10w + strwf

                        write(jhf,'(g14.5,'' ('',40a1)')
     &                  yp1, ( chaf(i), i = 1, icn ), ytn3

                     else if( iynum .eq.  2 ) then

                        xp1 = txl - ytlm - str10w
                        xp2 = txr + str10w + strwf

                        write(jhf,'(g14.5,'' ('',40a1)')
     &                  yp1, ( chaf(i), i = 1, icn ), ytn3

                     end if

               else if( ycta(n,1) .eq. '-' ) then

                  if( iytn4 .eq. 0 ) then

                        if( idbg .eq. 1 )  write(jhf,'()')

                     if( iynum .eq. 1 ) then

                        write(jhf,'(
     &                  ''/ytn4 {/str001 S N /yps S N'',/
     &                  ''(10) '',g14.5,'' yps ftn001 fts stsw'',/
     &                  ''(-) '',g14.5,'' yps '',g14.5,
     &                  '' add /Symbol findfont ftsu stsw'',/
     &                  ''str001 '',g14.5,'' yps '',g14.5,
     &                  '' add ftn001 ftsu stsw } N''
     &                  )') txl - ytlm - str10w,
     &                      txl - ytlm, ypup + rmup,
     &                      txl - ytlm + strm, ypup

                     else if( iynum .eq. -1 ) then

                        write(jhf,'(
     &                  ''/ytn4 {/str001 S N /yps S N'',/
     &                  ''(10) '',g14.5,'' yps ftn001 fts stsw'',/
     &                  ''(-) '',g14.5,'' yps '',g14.5,
     &                  '' add /Symbol findfont ftsu stsw'',/
     &                  ''str001 '',g14.5,'' yps '',g14.5,
     &                  '' add ftn001 ftsu stsw } N''
     &                  )') txr,
     &                      txr + str10w, ypup + rmup,
     &                      txr + str10w + strm, ypup

                     else if( iynum .eq. 2 ) then

                        write(jhf,'(
     &                  ''/ytn4 {/str001 S N /yps S N'',/
     &                  ''(10) '',g14.5,'' yps ftn001 fts stsw'',/
     &                  ''(10) '',g14.5,'' yps stsm'',/
     &                  ''(-) '',g14.5,'' yps '',g14.5,
     &                  '' add /Symbol findfont ftsu stsw'',/
     &                  ''(-) '',g14.5,'' yps '',g14.5,
     &                  '' add stsm'',/
     &                  ''str001 '',g14.5,'' yps '',g14.5,
     &                  '' add ftn001 ftsu stsw'',/
     &                  ''str001 '',g14.5,'' yps '',g14.5,
     &                  '' add stsm } N''
     &                  )') txl - ytlm - str10w, txr,
     &                      txl - ytlm, ypup + rmup,
     &                      txr + str10w, ypup + rmup,
     &                      txl - ytlm + strm, ypup,
     &                      txr + str10w + strm, ypup

                     end if

                        iytn4 = 1

                        if( idbg .eq. 1 )  write(jhf,'()')

                  end if


                     icn = nytch(n) - 1

                     do 819 i = 1, nytch(n) - 1

                        chaf(i) = ycta(n,i+1)

  819                continue

                     call strwhb(icn,chaf,ftsu,itfon,strwf,strhf,strbf,
     &                           ifon)

                        yp1 = ylta(n) * yal * cm - tfs * ftht / 2.0
                        yp2 = yp1 + ypup + strhf

                     if( iynum .eq. 1 ) then

                        xp1 = txl - ytlm - str10w
                        xp2 = txl

                        write(jhf,'(g14.5,'' ('',40a1)')
     &                  yp1, ( chaf(i), i = 1, icn ), ytn4

                     else if( iynum .eq. -1 ) then

                        xp1 = txr
                        xp2 = txr + str10w + strwf + strm

                        write(jhf,'(g14.5,'' ('',40a1)')
     &                  yp1, ( chaf(i), i = 1, icn ), ytn4

                     else if( iynum .eq.  2 ) then

                        xp1 = txl - ytlm - str10w
                        xp2 = txr + str10w + strwf + strm

                        write(jhf,'(g14.5,'' ('',40a1)')
     &                  yp1, ( chaf(i), i = 1, icn ), ytn4

                     end if

               end if

*-----------------------------------------------------------------------

                        yp1 = yp1 + strbf

                     if( iynum .eq. 1 ) then

                        xtmn = min(xtmn,xp1)

                        call bbox(1,0,xp1,yp1,0.d0)
                        call bbox(1,0,xp1,yp2,0.d0)

                     else if( iynum .eq. -1 ) then

                        xtma = max(xtma,xp2)

                        call bbox(1,0,xp2,yp1,0.d0)
                        call bbox(1,0,xp2,yp2,0.d0)

                     else if( iynum .eq.  2 ) then

                        xtmn = min(xtmn,xp1)
                        xtma = max(xtma,xp2)

                        call bbox(1,0,xp1,yp1,0.d0)
                        call bbox(1,0,xp1,yp2,0.d0)
                        call bbox(1,0,xp2,yp1,0.d0)
                        call bbox(1,0,xp2,yp2,0.d0)

                     end if

  730       continue

*-----------------------------------------------------------------------
         end if
*-----------------------------------------------------------------------

                  if( iynum .eq. 1 ) then

                     txl = min(txl,xtmn)

                  else if( iynum .eq. -1 ) then

                     txr = max(txr,xtma)

                  else if( iynum .eq. 2 ) then

                     txl = min(txl,xtmn)
                     txr = max(txr,xtma)

                  end if

*-----------------------------------------------------------------------
      end if
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*     WRITE AXIS TEXT
*-----------------------------------------------------------------------

*           Default Tic Font Size is 24 pt

*           XMARG = Margin of Text from number or axis : 0.70 * TFS
*           YMARG = Margin of Text from number or axis : 0.55 * TFS

*-----------------------------------------------------------------------
cKN 2024/01/24

      if( noxt .eq. 1 .or. noyt .eq. 1 .or.
     &    ieex .gt. -50 .or. ieey .gt. -50 ) then

*-----------------------------------------------------------------------

            if( idbg .eq. 1 )
     &      write(jhf,'(/''%'',71(''-''),/
     &                   ''% Write Axis Text''
     &                 ,/''%'',71(''-'')/)')

*-----------------------------------------------------------------------

                  tfs   = fstl * afac * atxs

                  xmarg = 0.70
                  ymarg = 0.55

*-----------------------------------------------------------------------

                  if( ix .eq. 0 ) then
                     xtt(1) = 'X'
                     ixtt   =  1
                  end if

                  if( iy .eq. 0 ) then
                     ytt(1) = 'Y'
                     iytt   =  1
                  end if

*-----------------------------------------------------------------------
cKN 2024/01/24

         if( ieex .gt. -50 ) then

               ic377 = 0
               ic000 = 0

               do k = 1, ixtt-3
                  if( xtt(k)//xtt(k+1)//xtt(k+2) .eq. '377' )
     &               ic377 = ic377 + 1
               end do
            if( ic377 .eq. 1 ) then
               do k = 1, ixtt-3
                  if( xtt(k)//xtt(k+1)//xtt(k+2) .eq. '000' )
     &               ic000 = ic000 + 1
               end do
            end if

*-----------------------------------------------------------------------

               ii = ixtt

            if( ic377 .gt. ic000 ) then
               d8 = '\377\000'
               do jj = 1, 8
                  ii = ii + 1
                  xtt(ii) = d8(jj:jj)
               end do
            end if

               if( noxt .eq. 0 ) then
                  ii = 0
                  noxt = 1
               end if

               ixtt = ii

*-----------------------------------------------------------------------

               d13 = ' \times\,10^{'
            do jj = 1, 13
               ii = ii + 1
               xtt(ii) = d13(jj:jj)
            end do

               write(d4,'(i4)') ieex
            do jj = 1, 4
               if( d4(jj:jj) .ne. ' ' ) then
                  ii = ii + 1
                  xtt(ii) = d4(jj:jj)
               end if
            end do

               ii = ii + 1
               xtt(ii) = '}'

               ixtt = ii

         end if

*-----------------------------------------------------------------------

         if( ieey .gt. -50 ) then

               ic377 = 0
               ic000 = 0

               do k = 1, iytt-3
                  if( ytt(k)//ytt(k+1)//ytt(k+2) .eq. '377' )
     &               ic377 = ic377 + 1
               end do
            if( ic377 .eq. 1 ) then
               do k = 1, iytt-3
                  if( ytt(k)//ytt(k+1)//ytt(k+2) .eq. '000' )
     &               ic000 = ic000 + 1
               end do
            end if

*-----------------------------------------------------------------------

               ii = iytt

            if( ic377 .gt. ic000 ) then
               d8 = '\377\000'
               do jj = 1, 8
                  ii = ii + 1
                  ytt(ii) = d8(jj:jj)
               end do
            end if

               if( noyt .eq. 0 ) then
                  ii = 0
                  noyt = 1
               end if

               iytt = ii

*-----------------------------------------------------------------------

               d13 = ' \times\,10^{'
            do jj = 1, 13
               ii = ii + 1
               ytt(ii) = d13(jj:jj)
            end do

               write(d4,'(i4)') ieey
            do jj = 1, 4
               if( d4(jj:jj) .ne. ' ' ) then
                  ii = ii + 1
                  ytt(ii) = d4(jj:jj)
               end if
            end do

               ii = ii + 1
               ytt(ii) = '}'

               iytt = ii

         end if

*-----------------------------------------------------------------------
cKN 2024/01/24 w

c          write(*,'(50a1)') (xtt(ii),ii=1,ixtt)
c          write(*,'(50a1)') (ytt(ii),ii=1,iytt)

cKN 2024/01/24 w
*-----------------------------------------------------------------------
*           DEFAULT FONT FOR AXIS TEXT
*-----------------------------------------------------------------------

*              IFTTT : Default Font For Axis Text

               ifttt = ibfon

*-----------------------------------------------------------------------
*           COLOR OF AXIS TEXT
*-----------------------------------------------------------------------

               if( cltx(1) .lt. -r0max ) cltx(1) = -2.0

*-----------------------------------------------------------------------

         if( noxt .eq. 1 ) then

               if( ixtxt .eq. 1 ) then

                     xpo = xtxp * xal * cm
                     ypo = tyd - tfs * xmarg

                  call wtext(jhf,chag,idbg,clal,clmo,0,xpo,ypo,
     &                       xtt,ixtt,2,3,0,0.d0,ifttt,ifon,
     &                       0,indt,tfs,cltx)

                     if( idbg .eq. 1 )  write(jhf,'()')

                     tyd = ypo - strh0

               else if( ixtxt .eq. -1 ) then

                     xpo = xtxp * xal * cm
                     ypo = tyu + tfs * xmarg

                  call wtext(jhf,chag,idbg,clal,clmo,0,xpo,ypo,
     &                       xtt,ixtt,2,1,0,0.d0,ifttt,ifon,
     &                       0,indt,tfs,cltx)

                     if( idbg .eq. 1 )  write(jhf,'()')

                     tyu = ypo + strh0

               else if( ixtxt .eq. 2 ) then

                     xpo = xtxp * xal * cm
                     ypo = tyd - tfs * xmarg

                  call wtext(jhf,chag,idbg,clal,clmo,0,xpo,ypo,
     &                       xtt,ixtt,2,3,0,0.d0,ifttt,ifon,
     &                       0,indt,tfs,cltx)

                     if( idbg .eq. 1 )  write(jhf,'()')

                     tyd = ypo - strh0


                     xpo = xtxp * xal * cm
                     ypo = tyu + tfs * xmarg

                  call wtext(jhf,chag,idbg,clal,clmo,0,xpo,ypo,
     &                       xtt,ixtt,2,1,0,0.d0,ifttt,ifon,
     &                       0,indt,tfs,cltx)

                     if( idbg .eq. 1 )  write(jhf,'()')

                     tyu = ypo + strh0

               end if

         end if

*-----------------------------------------------------------------------

         if( noyt .eq. 1 ) then

               if( iytxt .eq. 1 ) then

                     xpo = txl - tfs * ymarg
                     ypo = ytxp * yal * cm

                  call wtext(jhf,chag,idbg,clal,clmo,0,xpo,ypo,
     &                       ytt,iytt,2,1,1,90.0d0,ifttt,ifon,
     &                       0,indt,tfs,cltx)

                     if( idbg .eq. 1 )  write(jhf,'()')

                     txl = xpo - strh0

               else if( iytxt .eq. -1 ) then

                     xpo = txr + tfs * ymarg
                     ypo = ytxp * yal * cm

                  call wtext(jhf,chag,idbg,clal,clmo,0,xpo,ypo,
     &                       ytt,iytt,2,3,1,90.0d0,ifttt,ifon,
     &                       0,indt,tfs,cltx)

                     if( idbg .eq. 1 )  write(jhf,'()')

                     txr = xpo + strh0

               else if( iytxt .eq. 2 ) then

                     xpo = txl - tfs * ymarg
                     ypo = ytxp * yal * cm

                  call wtext(jhf,chag,idbg,clal,clmo,0,xpo,ypo,
     &                       ytt,iytt,2,1,1,90.0d0,ifttt,ifon,
     &                       0,indt,tfs,cltx)

                     if( idbg .eq. 1 )  write(jhf,'()')

                     txl = xpo - strh0


                     xpo = txr + tfs * ymarg
                     ypo = ytxp * yal * cm

                  call wtext(jhf,chag,idbg,clal,clmo,0,xpo,ypo,
     &                       ytt,iytt,2,3,1,90.d0,ifttt,ifon,
     &                       0,indt,tfs,cltx)

                     if( idbg .eq. 1 )  write(jhf,'()')

                     txr = xpo + strh0

               end if

         end if

*-----------------------------------------------------------------------

      end if

*=======================================================================

      deallocate(chag)
      return
      end


************************************************************************
*                                                                      *
      subroutine wlgnd(jhf,idbg,txr,regx,regy,regs,regd,
     &                 ncom,ilbox,
     &                 xmin,xmax,ixlog,ymin,ymax,iylog,
     &                 ibfon,fscm,cllg,clal,clmo,
     &                 ifon,clgb,clgl,clgs,xfac,iycm,rycm,
     &                 rlptl,sybw,erwd)
cKN 2024/01/24
*                                                                      *
*                                                                      *
*      PURPOSE  :   WRITE LEGEND COMMENTS, BOX, LINES AND SYMBOLS      *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter ( ss0 = 0.254 )
      parameter ( regx0  = 0.6, regxl0 = 1.4 )
      parameter ( regxc0 = 0.4, regyd0 = 0.9 )

      common /con/  cm, dd
      common /frm/  xal, yal
      common /wtval1/ strl0, strh0, strb0
      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb
      common /box/  bxws, bxwl, bxss, bxsl, bxds, bxdl, bxls, bxll

      character ccmt(ichrl)*1
      dimension iycm(mc,8),rycm(mc,7)

      character,allocatable:: chag(:,:)

      dimension clal(3), cllg(3), clgb(3), clgl(3), clgs(3)
      dimension bcb(3), bcl(3), bcs(3)

*-----------------------------------------------------------------------

      allocate(chag(inig,0:ichrl))

*-----------------------------------------------------------------------
*        rewind comment file
*-----------------------------------------------------------------------

            rewind(jhy)

*-----------------------------------------------------------------------
*        WRITE LEGEND LINES
*-----------------------------------------------------------------------

                  if( regx .lt. -r0max ) then

                     regx = txr / cm / xal + regx0 * regs / xal

                  else

                     regx = zonep(xmin,xmax,regx,ixlog)

                  end if

                     regxl = regxl0 * regs / xal
                     regxc = regx + regxl + regxc0 * regs / xal
                     regyd = regyd0 * regs * regd / yal

                  if( regy .lt. -r0max ) then

                     if( ncom / 2 * 2 .eq. ncom ) then

                        regy = 0.5 + regyd / 2.0
     &                       + dble( ncom / 2 - 1 ) * regyd

                     else

                        regy = 0.5
     &                       + dble( ( ncom - 1 ) / 2 )  * regyd

                     end if

                  else

                        regy = zonep(ymin,ymax,regy,iylog)

                  end if


                  if( ilbox .ne. 0 ) then

                        rebl = regx * xal * cm
                        rebr = rebl + regxl * xal * cm
                        rebu = regy * yal * cm
                        rebd = ( regy - dble( ncom - 1 ) * regyd )
     &                       * yal * cm

                        yddl = 0.0

                     if( iycm(1,7) .ne. 0 ) then

                        yddl = max( ss0 * rycm(1,1) * cm,
     &                               regyd * 0.3 * yal * cm )

                     end if

                     if( iycm(1,4) .gt. 0 ) then

                        yddl = max( yddl,
     &                              ss0 * rycm(1,1) * cm / 2.0 )

                     end if

                     if( iycm(1,6) .ne. 0 ) then

                        yddl = max( yddl, regyd * 0.3 * yal * cm )

                     end if

                        rebu = rebu + yddl


                        yddl = 0.0

                     if( iycm(ncom,7) .ne. 0 ) then

                        yddl = max( ss0 * rycm(ncom,1) * cm,
     &                               regyd * 0.3 * yal * cm )

                     end if

                     if( iycm(ncom,4) .gt. 0 ) then

                        yddl = max( yddl,
     &                              ss0 * rycm(ncom,1) * cm / 2.0 )

                     end if

                     if( iycm(ncom,6) .ne. 0 ) then

                        yddl = max( yddl, regyd * 0.3 * yal * cm )

                     end if

                        rebd = rebd - yddl


                  end if

*-----------------------------------------------------------------------
*     WRITE LEGEND COMMENT
*-----------------------------------------------------------------------

                  if( idbg .eq. 1 )
     &            write(jhf,'(/''%'',71(''-''),/
     &                         ''% Write Legend Comment''
     &                       ,/''%'',71(''-'')/)')

*-----------------------------------------------------------------------

               do 455 i = 1, ncom

                     iflgc = ibfon

                     tfs = fscm * regs

                     if( cllg(1) .lt. -r0max ) cllg(1) = -2.0


                     xpo = regxc * xal * cm
                     ypo = ( regy - dble( i - 1 ) * regyd ) * yal * cm


                     ixc  = 1
                     iyc  = 2


                     read(jhy) ( ccmt(k), k = 1, iycm(i,1) )


                  if( ilbox .eq. 0 ) then

                     call wtext(jhf,chag,idbg,clal,clmo,0,xpo,ypo,
     &                       ccmt,iycm(i,1),ixc,iyc,0,0.d0,iflgc,ifon,
     &                       0,indt,tfs,cllg)

                  else

                     iccb = 1

                     call wrbox(jhf,idbg,0,ilbox,iccb,
     &                          b1x, b2x, b1y, b2y,
     &                          bdd, bcb, bcl, bdc, bcs)

                     ilmul = i

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(''/mtl'',I3.3,'' {'')') ilmul

                     call wtext(jhf,chag,idbg,clal,clmo,-2,xpo,ypo,
     &                       ccmt,iycm(i,1),ixc,iyc,0,0.d0,iflgc,ifon,
     &                       ilmul,indt,tfs,cllg)

                        rebr = max(rebr, xpo + strl0 )
                        rebu = max(rebu, ypo + strh0 / 2.0 )
                        rebd = min(rebd, ypo - strh0 / 2.0 )

                  end if

                     if( idbg .eq. 1 )  write(jhf,'()')


  455          continue

*-----------------------------------------------------------------------
*           WRITE LEGEND BOX
*-----------------------------------------------------------------------

               if( ilbox .ne. 0 ) then

                  if( idbg .eq. 1 )
     &            write(jhf,'(/''%'',71(''-''),/
     &                         ''% Legend Box''
     &                       ,/''%'',71(''-'')/)')

                        if( clgb(1) .lt. -r0max .or.
     &                      clmo .gt. -r0max .or.
     &                      clal(1) .gt. -r0max ) clgb(1) = -1.0

                        if( clal(1) .gt. -r0max ) then

                           clgl(1) = clal(1)
                           clgl(2) = clal(2)
                           clgl(3) = clal(3)
                           clgs(1) = clal(1)
                           clgs(2) = clal(2)
                           clgs(3) = clal(3)

                        end if

                        if( clgl(1) .lt. -r0max .or.
     &                      clmo .gt. -r0max ) clgl(1) = -2.0

                        if( clgs(1) .lt. -r0max .or.
     &                      clmo .gt. -r0max ) clgs(1) = -2.0

                        bcb(1) = clgb(1)
                        bcb(2) = clgb(2)
                        bcb(3) = clgb(3)
                        bcl(1) = clgl(1)
                        bcl(2) = clgl(2)
                        bcl(3) = clgl(3)
                        bcs(1) = clgs(1)
                        bcs(2) = clgs(2)
                        bcs(3) = clgs(3)


                  if( ilbox .eq. 10 .or. ilbox .eq. 11. or.
     &                ilbox .eq. 20 .or. ilbox .eq. 21. or.
     &                ilbox .eq. 30 .or. ilbox .eq. 31. or.
     &                ilbox .eq. 40 .or. ilbox .eq. 41. or.
     &                ilbox .eq. 50 .or. ilbox .eq. 51. ) then

                        bds = tfs * bxwl

                  else

                        bds = tfs * bxwl * 2.0

                  end if

                  if( ilbox .eq. 11 .or. ilbox .eq. 13. or.
     &                ilbox .eq. 21 .or. ilbox .eq. 23. or.
     &                ilbox .eq. 31 .or. ilbox .eq. 33. or.
     &                ilbox .eq. 41 .or. ilbox .eq. 43. or.
     &                ilbox .eq. 51 .or. ilbox .eq. 53. ) then

                        bdd = tfs * bxll

                  else

                        bdd = tfs * bxls

                  end if

                  if( ilbox .eq. 20 .or. ilbox .eq. 21. or.
     &                ilbox .eq. 40 .or. ilbox .eq. 41. or.
     &                ilbox .eq. 50 .or. ilbox .eq. 51. ) then

                        bdc = tfs * bxsl

                  else
     &            if( ilbox .eq. 22 .or. ilbox .eq. 23. or.
     &                ilbox .eq. 42 .or. ilbox .eq. 43. or.
     &                ilbox .eq. 52 .or. ilbox .eq. 53. ) then

                        bdc = tfs * bxsl * 2.0

                  end if

                  if( ilbox .eq. 30 .or. ilbox .eq. 31. ) then

                        bdc = tfs * bxdl

                  else
     &            if( ilbox .eq. 32 .or. ilbox .eq. 33. ) then

                        bdc = tfs * bxdl * 2.0

                  end if

                     b1x = rebl - bds
                     b2x = rebr + bds

                     b1y = rebd - bds
                     b2y = rebu + bds

*-----------------------------------------------------------------------

                     iccb = 0

                     call wrbox(jhf,idbg,0,ilbox,iccb,
     &                          b1x, b2x, b1y, b2y,
     &                          bdd, bcb, bcl, bdc, bcs)

*-----------------------------------------------------------------------

                  if( ilbox .ge. 40 ) then

                     b2x = b2x + bdc
                     b1y = b1y - bdc

                  end if

                     call bbox(1,0,b1x,b1y,0.d0)
                     call bbox(1,0,b1x,b2y,0.d0)
                     call bbox(1,0,b2x,b1y,0.d0)
                     call bbox(1,0,b2x,b2y,0.d0)

*-----------------------------------------------------------------------

                  do 457 i = 1, ncom

                     write(jhf,'(''mtl'',i3.3,'' mts'',i3.3)') i, i

  457             continue

               end if

*-----------------------------------------------------------------------
*           WRITE LEGEND LINE AND SYMBOL
*-----------------------------------------------------------------------

                  if( idbg .eq. 1 )
     &            write(jhf,'(/''%'',71(''-''),/
     &                         ''% Legend Line and Symbol''
     &                       ,/''%'',71(''-'')/)')

                  call lgline(ixlog,iylog,xmax,xmin,ymax,ymin,
     &                        regx,regy,regs,regd,idbg,xfac,
     &                        ncom,iycm,rycm,regxl,regxc,regyd,jhf,txr,
     &                        clal,clmo,rlptl,sybw,erwd)
cKN 2024/01/24


*-----------------------------------------------------------------------

      deallocate(chag)
      return
      end

************************************************************************
*                                                                      *
      subroutine commt0(jhf,idbg,ifon,
     &                  iw,wxys,iwx,iwy,iwn,iwf,iwb,ibox,cbox,
     &                  fscm,clal,clmo,xmax,xmin,ymax,ymin,
     &                  rgxm,rgym,
     &                  ixlog,iylog,ibfon)
*                                                                      *
*        PURPOSE    :  WRITE NORMAL COMMENTS                           *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      dimension wxys(mc,9),iwx(mc),iwy(mc),iwn(mc),iwf(mc),iwb(mc)
      dimension ibox(mc),cbox(mc,3,3)

      dimension jndt(mc)

      character ccmt(ichrl)*1

      common /frm/  xal, yal
      common /con/  cm, dd

      common /lbox/  bcb(3), bcl(3), bcs(3), bds, bdd, bdc

      common /box/  bxws, bxwl, bxss, bxsl, bxds, bxdl, bxls, bxll
      common /suf/  ftht, ftss, ftks, ftsv,
     &              yplu, ypld, ypku, ypkd, yplh, ypls

      common /wtval1/ strl0, strh0, strb0
      common /wtval7/ strhtp, strbtp, strhbt, strbbt, strhall
      common /wtval8/ xpsi(0:inig), ypsi(0:inig), strl(0:inig),
     &                xpsl(0:inig), xpsr(0:inig), rind(0:inig)

      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb

      dimension rxp(0:inig)

      character,allocatable:: chag(:,:)
      dimension clal(3)
      dimension rcol(3)

*-----------------------------------------------------------------------

      allocate(chag(inig,0:ichrl))

*-----------------------------------------------------------------------
*     rewind information file
*-----------------------------------------------------------------------

            rewind(jhm)

*-----------------------------------------------------------------------
*     WRITE COMMENTS
*-----------------------------------------------------------------------

            if( idbg .eq. 1 )
     &      write(jhf,'(/''%'',71(''-''),/
     &                   ''% Write Normal Comment''
     &                 ,/''%'',71(''-''))')

*-----------------------------------------------------------------------
*           BLSI  DEFAULT BASE LINE SKIP : JPN 1.5 FTS (NOT NOW)
*                                          ENG 1.2 FTS (ONLY)
*-----------------------------------------------------------------------

            if( ifon .eq. 0 ) then

               blsi = 1.2

            else

               blsi = 1.2

            end if

*-----------------------------------------------------------------------
*        DO LOOP FOR COMMENTS
*-----------------------------------------------------------------------

                  imull = 0

*-----------------------------------------------------------------------

         do 410 i = 1, iw

*-----------------------------------------------------------------------

                  ipin = 0
                  indt = 0

            if( iwn(i) .eq. 0 ) goto 410

            if( iwb(i) .eq. -3 .and. imull .eq. 0 ) then

                  iwb(i) = -4

            end if

*-----------------------------------------------------------------------
*        INITIALIZATION
*-----------------------------------------------------------------------

            if( wxys(i,7) .gt. -r0max ) then

                  bls = wxys(i,7)

            else

                  bls = blsi

            end if

            if( wxys(i,9) .gt. -r0max ) then

                  bbls = wxys(i,9)

            else

                  bbls = blsi

            end if

*-----------------------------------------------------------------------

                  tfs = fscm * wxys(i,3)

                  rcol(1) = wxys(i,4)
                  rcol(2) = wxys(i,5)
                  rcol(3) = wxys(i,6)

                  if( clal(1) .gt. -r0max ) then
                     rcol(1) = clal(1)
                     rcol(2) = clal(2)
                     rcol(3) = clal(3)
                  end if

                  if( clmo .gt. -r0max .and.
     &                rcol(1) .gt. 0.0 )    rcol(1) = -2.0
                  if( rcol(1) .lt. -r0max ) rcol(1) = -2.0

                  bls  = tfs * bls
                  bbls = tfs * bbls


                  if( wxys(i,1) .lt. -r0max ) wxys(i,1) = rgxm
                  if( wxys(i,2) .lt. -r0max ) wxys(i,2) = rgym

                  xpo = zonep(xmin,xmax,wxys(i,1),ixlog) * xal * cm
                  ypo = zonep(ymin,ymax,wxys(i,2),iylog) * yal * cm

*-----------------------------------------------------------------------

               if( ibox(i) .ne. 0 ) then

                     iccb = 1

                     call wrbox(jhf,idbg,0,ibox(i),iccb,
     &                          b1x, b2x, b1y, b2y,
     &                          bdd, bcb, bcl, bdc, bcs)

               end if

*-----------------------------------------------------------------------
*        ONE LINE COMENT
*-----------------------------------------------------------------------

            if( iwb(i) .eq. -1 .or. iwb(i) .eq. -4 ) then

                  ipin  = 0
                  imull = 0

               if( iwb(i) .eq. -4 .and. ibox(i) .ne. 0 ) then

                     ipin  = ibox(i)

               end if

                  if( idbg .eq. 1 )  write(jhf,'()')

*-----------------------------------------------------------------------
*        ONE LINE WITH BOX FROM WT:    IPIN > 4 = IBOX(I)
*-----------------------------------------------------------------------

               if( ipin .gt. 4 ) then

                        if( cbox(i,1,1) .lt. -r0max .or.
     &                      clmo .gt. -r0max .or.
     &                      clal(1) .gt. -r0max ) cbox(i,1,1) = -1.0

                        if( clal(1) .gt. -r0max ) then

                           cbox(i,2,1) = clal(1)
                           cbox(i,2,2) = clal(2)
                           cbox(i,2,3) = clal(3)
                           cbox(i,3,1) = clal(1)
                           cbox(i,3,2) = clal(2)
                           cbox(i,3,3) = clal(3)

                        end if

                        if( cbox(i,2,1) .lt. -r0max .or.
     &                      clmo .gt. -r0max ) cbox(i,2,1) = -2.0

                        if( cbox(i,3,1) .lt. -r0max .or.
     &                      clmo .gt. -r0max ) cbox(i,3,1) = -2.0

                        bcb(1) = cbox(i,1,1)
                        bcb(2) = cbox(i,1,2)
                        bcb(3) = cbox(i,1,3)
                        bcl(1) = cbox(i,2,1)
                        bcl(2) = cbox(i,2,2)
                        bcl(3) = cbox(i,2,3)
                        bcs(1) = cbox(i,3,1)
                        bcs(2) = cbox(i,3,2)
                        bcs(3) = cbox(i,3,3)


                  if( ibox(i) .eq. 10 .or. ibox(i) .eq. 11. or.
     &                ibox(i) .eq. 20 .or. ibox(i) .eq. 21. or.
     &                ibox(i) .eq. 30 .or. ibox(i) .eq. 31. or.
     &                ibox(i) .eq. 40 .or. ibox(i) .eq. 41. or.
     &                ibox(i) .eq. 50 .or. ibox(i) .eq. 51. ) then

                        bds = tfs * bxwl

                  else

                        bds = tfs * bxwl * 2.0

                  end if

                  if( ibox(i) .eq. 11 .or. ibox(i) .eq. 13. or.
     &                ibox(i) .eq. 21 .or. ibox(i) .eq. 23. or.
     &                ibox(i) .eq. 31 .or. ibox(i) .eq. 33. or.
     &                ibox(i) .eq. 41 .or. ibox(i) .eq. 43. or.
     &                ibox(i) .eq. 51 .or. ibox(i) .eq. 53. ) then

                        bdd = tfs * bxll

                  else

                        bdd = tfs * bxls

                  end if

                  if( ibox(i) .eq. 20 .or. ibox(i) .eq. 21. or.
     &                ibox(i) .eq. 40 .or. ibox(i) .eq. 41. or.
     &                ibox(i) .eq. 50 .or. ibox(i) .eq. 51. ) then

                        bdc = tfs * bxsl

                  else
     &            if( ibox(i) .eq. 22 .or. ibox(i) .eq. 23. or.
     &                ibox(i) .eq. 42 .or. ibox(i) .eq. 43. or.
     &                ibox(i) .eq. 52 .or. ibox(i) .eq. 53. ) then

                        bdc = tfs * bxsl * 2.0

                  end if

                  if( ibox(i) .eq. 30 .or. ibox(i) .eq. 31. ) then

                        bdc = tfs * bxdl

                  else
     &            if( ibox(i) .eq. 32 .or. ibox(i) .eq. 33. ) then

                        bdc = tfs * bxdl * 2.0

                  end if

               end if

*-----------------------------------------------------------------------
*        MULTI LINE COMMENT
*-----------------------------------------------------------------------

            else if( iwb(i) .eq. -2 .and. iwn(i) .lt. 0 ) then

                  imull = 1

                  kndt        = 0
                  jndt(imull) = 0

                  if( idbg .eq. 1 )  write(jhf,'()')

                  write(jhf,'(''/mtl001 {} B /mts001 {} B'')')

                  xpsi(imull) = 0.0
                  ypsi(imull) = 0.0

                  strhall = 0.0

                  strhtp = bbls / bls * tfs * ftht
                  strbtp = 0.0


            else if( iwb(i) .eq. -2 ) then

                  imull = 1

                  kndt        = 0
                  jndt(imull) = 0

                  ipin = 1

                  xpsi(imull) = 0.0
                  ypsi(imull) = 0.0

                  strhall = 0.0


                  if( idbg .eq. 1 )  write(jhf,'()')

                  write(jhf,'(''/mtl'',i3.3,'' {'')') imull


            else if( iwb(i) .eq. i - 1 .or. iwb(i) .eq. -3 ) then

               if( iwb(i) .eq. i - 1 ) then

                  ipin = 2

               else if( iwb(i) .eq. -3 ) then

                  ipin = 3

               end if

                  imull = imull + 1


                  if( idbg .eq. 1 )  write(jhf,'()')

                  write(jhf,'(''/mtl'',i3.3,'' {'')') imull


               if( iwn(i) .lt. 0 ) then

                  xpsi(imull) = 0.0
                  ypsi(imull) = ypsi(imull-1) - bbls

                  strhall = strhall + bbls

               else

                  xpsi(imull) = 0.0
                  ypsi(imull) = ypsi(imull-1) - bls

                  strhall = strhall + bls

               end if

            end if

*-----------------------------------------------------------------------

            if( iwn(i) .lt. 0 ) then

               if( ipin .eq. 2 .or. ipin .eq. 3 ) then

                  if( imull .eq. 2 ) then

                           jndt(1) = 0
                           jndt(2) = 0

                  else if( imull .gt. 2 ) then

                     if( jndt(imull-1) .eq. 0 ) then

                           jndt(imull) = 0

                     else if( jndt(imull-1) .gt. 0 ) then

                        if( jndt(imull-2) .eq. 0 ) then

                           jndt(imull-1) = 0
                           jndt(imull)   = 0

                        else if( jndt(imull-2) .gt. 0 ) then

                           jndt(imull-1) = -1
                           jndt(imull)   =  0
                           kndt          =  1

                        end if

                     end if

                  end if

               end if


               if( ipin .eq. 2 ) then

                  write(jhf,'(''} B /mts'',i3.3,'' {} B'')') imull

                  if( idbg .eq. 1 )  write(jhf,'()')

               else if( ipin .eq. 3 ) then

                  strhbt = strh0
                  strbbt = 0.0


                  write(jhf,'(''} B /mts'',i3.3,'' {} B'')') imull

                  if( idbg .eq. 1 )  write(jhf,'()')

                  goto 412

               end if

                  goto 410

            end if

*-----------------------------------------------------------------------

               ixc  = iwx(i)
               iyc  = iwy(i)

               read(jhm) (ccmt(k),k=1,iwn(i))

*-----------------------------------------------------------------------
*     FONT FOR NORMAL COMMENT
*-----------------------------------------------------------------------

            if( iwf(i) .lt. 0 ) then

               ifnoc = ibfon

            else

               ifnoc = iwf(i)

            end if

*-----------------------------------------------------------------------

            if( wxys(i,8) .lt. -r0max ) then

               call wtext(jhf,chag,idbg,clal,clmo,ipin,xpo,ypo,
     &                    ccmt,iwn(i),ixc,iyc,0,0.d0,ifnoc,ifon,
     &                    imull,indt,tfs,rcol)

            else

               call wtext(jhf,chag,idbg,clal,clmo,ipin,xpo,ypo,
     &                    ccmt,iwn(i),ixc,iyc,1,wxys(i,8),ifnoc,ifon,
     &                    imull,indt,tfs,rcol)

            end if


*-----------------------------------------------------------------------
*        SUMMARY OF INDENT FOR MULTI LINE
*-----------------------------------------------------------------------

            if( ipin .gt. 0 ) then

               if( indt .gt. 0 ) then

                  if( imull .eq. 1 ) then

                           jndt(imull) = 1

                  else if( imull .gt. 1 ) then

                     if( jndt(imull-1) .eq. 0 ) then

                           jndt(imull) = 1

                     else if( ipin .eq. 3 ) then

                           jndt(imull) = -1
                           kndt        =  1

                     else if( jndt(imull-1) .gt. 0 ) then

                           jndt(imull) = jndt(imull-1) + 1

                     end if

                  end if

               else if( indt .eq. 0 ) then

                  if( imull .eq. 1 ) then

                           jndt(imull) = 0

                  else if( imull .eq. 2 ) then

                           jndt(1) = 0
                           jndt(2) = 0

                  else if( imull .gt. 2 ) then

                     if( jndt(imull-1) .eq. 0 ) then

                           jndt(imull) = 0

                     else if( jndt(imull-1) .gt. 0 ) then

                        if( jndt(imull-2) .eq. 0 ) then

                           jndt(imull-1) = 0
                           jndt(imull)   = 0

                        else if( jndt(imull-2) .gt. 0 ) then

                           jndt(imull-1) = -1
                           jndt(imull)   =  0
                           kndt          =  1

                        end if

                     end if

                  end if

               end if

            end if

*-----------------------------------------------------------------------
*        SUMMARY OF MULTI LINE
*-----------------------------------------------------------------------

  412       if( ipin .eq. 3 ) then

*-----------------------------------------------------------------------
*              SUMMARY OF INDENT
*-----------------------------------------------------------------------

               if( kndt .ne. 0 ) then

                        iindt = 0
                        ifndt = 0

                        rind0 = 0.0

                  do 517 im = 1, imull

                     if( jndt(im) .eq. 0 ) goto 517

                        ii = i - imull + im

                     if( jndt(im) .eq. 1 ) then

                        if( iwx(ii) .eq. 1 ) then

                           rind0 = rind(im)

                        else if( iwx(ii) .eq. 2 ) then

                           rindr = strl(im) - rind(im)
                           rindl = rind(im)

                           rxp(im) = strl(im) / 2.0 - rind(im)

                        else if( iwx(ii) .eq. 3 ) then

                           rind0 = strl(im) - rind(im)

                        end if

                           iindt = im

                     else if( jndt(im) .gt.  1 .or.
     &                        jndt(im) .eq. -1 ) then

                        if( iwx(ii) .eq. 1 ) then

                           rind0 = max(rind0,rind(im))

                        else if( iwx(ii) .eq. 2 ) then

                           rindr = max(rindr,strl(im)-rind(im))
                           rindl = max(rindl,rind(im))

                           rxp(im) = strl(im) / 2.0 - rind(im)

                        else if( iwx(ii) .eq. 3 ) then

                           rind0 = max(rind0,strl(im)-rind(im))

                        end if

                     end if

                     if( jndt(im) .eq. -1 ) then

                           ifndt = im

                        if( iwx(ii) .eq. 2 ) then

                           rindc = ( rindl + rindr ) / 2.0 - rindl

                        end if

                      do 518 kn = iindt, ifndt

                        if( iwx(ii) .eq. 1 ) then

                           rxp(kn) = rind0 - rind(kn)

                        else if( iwx(ii) .eq. 2 ) then

                           rxp(kn) = rxp(kn) - rindc

                        else if( iwx(ii) .eq. 3 ) then

                           rxp(kn) = - rind0 + ( strl(kn) - rind(kn) )

                        end if

  518                 continue


                     end if

  517             continue

               end if

*-----------------------------------------------------------------------
*              SUMMURY OF LEFT AND RIGHT X POSITION FOR MULTI-LINE
*-----------------------------------------------------------------------

                           xpsl0 =  10000.0
                           xpsr0 = -10000.0

                  do 516 im = 1, imull

                        ii = i - imull + im

                     if( iwn(ii) .gt. 0 ) then

                        if( jndt(im) .eq. 0 ) then

                           xpsl0 = min(xpsl0,xpsl(im))
                           xpsr0 = max(xpsr0,xpsr(im))

                        else

                           xpsl0 = min(xpsl0,xpsl(im)+rxp(im))
                           xpsr0 = max(xpsr0,xpsr(im)+rxp(im))

                        end if

                     end if

  516             continue

*-----------------------------------------------------------------------
*           STARTING POSITION
*-----------------------------------------------------------------------

                  xpsin = xpo
                  ypsin = ypo

*-----------------------------------------------------------------------
*           ROTATION
*-----------------------------------------------------------------------

               if( wxys(i,8) .lt. -r0max ) then

                  iang = 0

               else

                  iang = 1

                  if( idbg .eq. 1 )  write(jhf,'()')

                  write(jhf,'(3g14.5,'' rotin'')')
     &                        xpsin, ypsin, wxys(i,8)

                  xpsain = xpsin
                  ypsain = ypsin
                  angain = wxys(i,8)

                  call bbox(1,1,xpsain,ypsain,angain)

                  xpsin = 0.0
                  ypsin = 0.0

               end if


                  xpss = xpsin
                  ypss = ypsin

*-----------------------------------------------------------------------
*           DETERMINATION OF THE INITIAL POSITION
*-----------------------------------------------------------------------

               if( ibox(i) .eq. 0 ) then

                  if( iwy(i) .eq. 1 ) then

                     bb1y = ypss

                     ypss = ypss - strbbt + strhall

                     bb2y = ypss + strhtp + strbtp

                  else if( iwy(i) .eq. 2 ) then

                     bb1y = ypss
     &                    - ( strhall - strbbt + strhtp + strbtp ) / 2.0

                     bb2y = ypss
     &                    + ( strhall - strbbt + strhtp + strbtp ) / 2.0

                     ypss = ypss
     &                    + ( strhall - strbbt + strhtp + strbtp ) / 2.0
     &                    - ( strhtp + strbtp )

                  else if( iwy(i) .eq. 3 ) then

                     bb2y = ypss

                     ypss = ypss - ( strhtp + strbtp )

                     bb1y = ypss - strhall + strbbt

                  end if

                     bb1x = xpsl0 + xpss
                     bb2x = xpsr0 + xpss


*-----------------------------------------------------------------------
*        BOX FOR WT: TEXT
*-----------------------------------------------------------------------

               else if( ibox(i) .ne. 0 ) then

                        if( cbox(i,1,1) .lt. -r0max .or.
     &                      clmo .gt. -r0max .or.
     &                      clal(1) .gt. -r0max ) cbox(i,1,1) = -1.0

                        if( clal(1) .gt. -r0max ) then

                           cbox(i,2,1) = clal(1)
                           cbox(i,2,2) = clal(2)
                           cbox(i,2,3) = clal(3)
                           cbox(i,3,1) = clal(1)
                           cbox(i,3,2) = clal(2)
                           cbox(i,3,3) = clal(3)

                        end if

                        if( cbox(i,2,1) .lt. -r0max .or.
     &                      clmo .gt. -r0max ) cbox(i,2,1) = -2.0

                        if( cbox(i,3,1) .lt. -r0max .or.
     &                      clmo .gt. -r0max ) cbox(i,3,1) = -2.0

                        bcb(1) = cbox(i,1,1)
                        bcb(2) = cbox(i,1,2)
                        bcb(3) = cbox(i,1,3)
                        bcl(1) = cbox(i,2,1)
                        bcl(2) = cbox(i,2,2)
                        bcl(3) = cbox(i,2,3)
                        bcs(1) = cbox(i,3,1)
                        bcs(2) = cbox(i,3,2)
                        bcs(3) = cbox(i,3,3)


                  if( ibox(i) .eq. 10 .or. ibox(i) .eq. 11. or.
     &                ibox(i) .eq. 20 .or. ibox(i) .eq. 21. or.
     &                ibox(i) .eq. 30 .or. ibox(i) .eq. 31. or.
     &                ibox(i) .eq. 40 .or. ibox(i) .eq. 41. or.
     &                ibox(i) .eq. 50 .or. ibox(i) .eq. 51. ) then

                        bds = tfs * bxwl

                  else

                        bds = tfs * bxwl * 2.0

                  end if

                  if( ibox(i) .eq. 11 .or. ibox(i) .eq. 13. or.
     &                ibox(i) .eq. 21 .or. ibox(i) .eq. 23. or.
     &                ibox(i) .eq. 31 .or. ibox(i) .eq. 33. or.
     &                ibox(i) .eq. 41 .or. ibox(i) .eq. 43. or.
     &                ibox(i) .eq. 51 .or. ibox(i) .eq. 53. ) then

                        bdd = tfs * bxll

                  else

                        bdd = tfs * bxls

                  end if

                  if( ibox(i) .eq. 20 .or. ibox(i) .eq. 21. or.
     &                ibox(i) .eq. 40 .or. ibox(i) .eq. 41. or.
     &                ibox(i) .eq. 50 .or. ibox(i) .eq. 51. ) then

                        bdc = tfs * bxsl

                  else
     &            if( ibox(i) .eq. 22 .or. ibox(i) .eq. 23. or.
     &                ibox(i) .eq. 42 .or. ibox(i) .eq. 43. or.
     &                ibox(i) .eq. 52 .or. ibox(i) .eq. 53. ) then

                        bdc = tfs * bxsl * 2.0

                  end if

                  if( ibox(i) .eq. 30 .or. ibox(i) .eq. 31. ) then

                        bdc = tfs * bxdl

                  else
     &            if( ibox(i) .eq. 32 .or. ibox(i) .eq. 33. ) then

                        bdc = tfs * bxdl * 2.0

                  end if


                  if( iwx(i) .eq. 1 ) then

                     xpss = xpss + bds

                  else if( iwx(i) .eq. 3 ) then

                     xpss = xpss - bds

                  end if

                     b1x = xpsl0 + xpss - bds
                     b2x = xpsr0 + xpss + bds

                  if( iwy(i) .eq. 1 ) then

                     b1y  = ypss
                     ypss = ypss + bds - strbbt + strhall
                     b2y  = ypss + strhtp + strbtp + bds

                  else if( iwy(i) .eq. 2 ) then

                     b1y  = ypss
     &                    - ( strhall - strbbt + strhtp + strbtp ) / 2.0
     &                    - bds

                     b2y  = ypss
     &                    + ( strhall - strbbt + strhtp + strbtp ) / 2.0
     &                    + bds

                     ypss = ypss
     &                    + ( strhall - strbbt + strhtp + strbtp ) / 2.0
     &                    - ( strhtp + strbtp )

                  else if( iwy(i) .eq. 3 ) then

                     b2y  = ypss
                     ypss = ypss - bds - ( strhtp + strbtp )

                     b1y  = ypss - strhall + strbbt - bds

                  end if

                     bb1x = b1x
                     bb2x = b2x

                     bb1y = b1y
                     bb2y = b2y

                  if( ibox(i) .ge. 40 ) then

                     bb2x = bb2x + bdc
                     bb1y = bb1y - bdc

                  end if

*-----------------------------------------------------------------------

                     iccb = 0

                     call wrbox(jhf,idbg,0,ibox(i),iccb,
     &                          b1x, b2x, b1y, b2y,
     &                          bdd, bcb, bcl, bdc, bcs)

*-----------------------------------------------------------------------

               end if

*-----------------------------------------------------------------------

                     call bbox(1,0,bb1x,bb1y,0.d0)
                     call bbox(1,0,bb1x,bb2y,0.d0)
                     call bbox(1,0,bb2x,bb1y,0.d0)
                     call bbox(1,0,bb2x,bb2y,0.d0)

*-----------------------------------------------------------------------
*           STARTING POINT
*-----------------------------------------------------------------------

                  iindt = 0

                  if( idbg .eq. 1 )  write(jhf,'()')

                  write(jhf,'(''/yps '',g14.5,'' N'')')
     &                           ypss

*-----------------------------------------------------------------------

               do 650 im = 1, imull

                  if( jndt(im) .eq. 0 .and. iindt .eq. 0 ) then

                        write(jhf,'(''/xps '',g14.5,'' N'')')
     &                                 xpss

                        iindt = 1

                  else if( jndt(im) .ne. 0 ) then

                        write(jhf,'(''/xps '',g14.5,'' N'')')
     &                                 xpss + rxp(im)

                        iindt = 0

                  end if

                     write(jhf,'(''mtl'',i3.3,'' mts'',i3.3)') im, im


  650          continue

*-----------------------------------------------------------------------

                  if( iang .eq. 1 ) then

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(''rotout'')')

                     call bbox(1,1,0.d0,0.d0,-angain)
                     call bbox(1,1,-xpsain,-ypsain,0.d0)

                  end if

                  imull = 0

*-----------------------------------------------------------------------
*        END OF SUMMARY OF MULTI LINE
*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------

  410    continue

*-----------------------------------------------------------------------

      deallocate(chag)
      return
      end


************************************************************************
*                                                                      *
      subroutine commt1(jhf,idbg,ifon,iv,
     &                  fscm,clal,clmo,xmax,xmin,ymax,ymin,
     &                  rgxm,rgym,
     &                  ixlog,iylog,ibfon)
*                                                                      *
*        PURPOSE    :  WRITE TABLE COMMENT                             *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character ccmt(ichrl)*1

      common /frm/  xal, yal
      common /con/  cm, dd

      common /lbox/  bcb(3), bcl(3), bcs(3), bds, bdd, bdc
      common /box/  bxws, bxwl, bxss, bxsl, bxds, bxdl, bxls, bxll
      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb

      common /wtval1/ strl0, strh0, strb0
      common /wtval3/ strla(0:inig), strha(0:inig), strba(0:inig)
      common /wtval7/ strhtp, strbtp, strhbt, strbbt, strhall

*-----------------------------------------------------------------------

      dimension vxys(13)
      dimension vspt(2,0:nr)
      dimension ittx(0:nr,0:nc), ittv(0:nr,0:nc), itth(0:nr,0:nc)
      dimension csl(0:nr,0:nc,3)
      dimension lvm(0:nr,0:nc)

*-----------------------------------------------------------------------

      dimension tsrl(nr,nc), tsrh(nr,nc), tsrb(nr,nc)
      dimension tsrll(nc), tsrhh(nr), tsrbb(nr)
      dimension txps(nr,nc), typs(0:nr), txpm(nc)

      character,allocatable:: chag(:,:)
      dimension clal(3), rcol(3), tcb(3), tcl(3)

*-----------------------------------------------------------------------

      allocate(chag(inig,0:ichrl))

*-----------------------------------------------------------------------
*     rewind file
*-----------------------------------------------------------------------

            rewind(jwt)

*-----------------------------------------------------------------------
*     WRITE TABLE COMMENT
*-----------------------------------------------------------------------


            if( idbg .eq. 1 )
     &      write(jhf,'(/''%'',71(''-''),/
     &                   ''% Write Table Comment''
     &                 ,/''%'',71(''-'')/)')

*-----------------------------------------------------------------------

         do 750 i = 1, iv

*-----------------------------------------------------------------------
*           read information from file
*-----------------------------------------------------------------------

                  read(jwt) ivx,ivy,ivr,ivc,ivf,ivs,ittl

               do ii = 1, 13
                  read(jwt) vxys(ii)
               end do

               do ii = 0, ivr
                  read(jwt) vspt(1,ii), vspt(2,ii)
               do jj = 0, ivc
                  read(jwt) ittx(ii,jj), ittv(ii,jj), itth(ii,jj)
               end do
               end do

               do ii = 1, ivr
               do jj = 1, ivc
                  read(jwt) csl(ii,jj,1), csl(ii,jj,2), csl(ii,jj,3)
                  read(jwt) lvm(ii,jj)
               end do
               end do

*-----------------------------------------------------------------------

                  tfs = fscm * vxys(4)

                  rcol(1) = vxys(5)
                  rcol(2) = vxys(6)
                  rcol(3) = vxys(7)

                  if( clal(1) .gt. -r0max ) then
                     rcol(1) = clal(1)
                     rcol(2) = clal(2)
                     rcol(3) = clal(3)
                  end if

                  if( clmo .gt. -r0max .and. rcol(1) .gt. 0.0 )
     &               rcol(1) = -2.0
                  if( rcol(1) .lt. -r0max ) rcol(1) = -2.0

                  if( vxys(1) .lt. -r0max ) vxys(1) = rgxm
                  if( vxys(2) .lt. -r0max ) vxys(2) = rgym

                  xpo = zonep(xmin,xmax,vxys(1),ixlog) * xal * cm
                  ypo = zonep(ymin,ymax,vxys(2),iylog) * yal * cm

                  xpoi = 0.0
                  ypoi = 0.0

*-----------------------------------------------------------------------
*        FONT FOR TABLE COMMENT
*-----------------------------------------------------------------------

            if( ivf .lt. 0 ) then

               iftab = ibfon

            else

               iftab = ivf

            end if

*-----------------------------------------------------------------------
*        WRITE EACH ELEMENTS
*-----------------------------------------------------------------------

                  ipin = -1

            do 752 l = 1, ivr
            do 753 m = 1, ivc

                  iclng = lvm(l,m)

                  imull = - ( l * 1000 + m )

*-----------------------------------------------------------------------

               if( iclng .gt. 0 ) then

                     read(jwt) (ccmt(k),k=1,iclng)

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(''/mtl'',i6.6,'' {'')') -imull

                  call wtext(jhf,chag,idbg,clal,clmo,ipin,xpoi,ypoi,
     &                       ccmt,iclng,1,1,0,0.d0,iftab,ifon,
     &                       imull,indt,tfs,rcol)


                     tsrl(l,m) = strla(0)
                     tsrh(l,m) = strha(0)
                     tsrb(l,m) = strba(0)

*-----------------------------------------------------------------------

               else if( iclng .le. 0 ) then

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(''/mtl'',i6.6,'' {} B /mts'',i6.6,
     &                           '' {} B'')')
     &                           -imull, -imull

                     tsrl(l,m) = 0.0
                     tsrh(l,m) = 0.0
                     tsrb(l,m) = 0.0

               end if

  753       continue
  752       continue

*-----------------------------------------------------------------------
*        SUMMARY OF THE TABLE
*-----------------------------------------------------------------------


               do 754 l = 1, ivr

                     tsrhh(l) = 0.0
                     tsrbb(l) = 0.0

                  do 755 m = 1, ivc

                     tsrhh(l) = max(tsrhh(l),tsrh(l,m))
                     tsrbb(l) = min(tsrbb(l),tsrb(l,m))

  755             continue

  754          continue


               do 756 m = 1, ivc

                     tsrll(m) = 0.0

                  do 757 l = 1, ivr

                     tsrll(m) = max(tsrll(m),tsrl(l,m))

  757             continue

  756          continue

*-----------------------------------------------------------------------
*           COLOR AND SPACE
*-----------------------------------------------------------------------

                  if( vxys(8) .lt. -r0max .or.
     &                clmo .gt. -r0max .or.
     &                clal(1) .gt. -r0max ) vxys(8) = -1.0

                  if( clal(1) .gt. -r0max ) then

                     vxys(11) = clal(1)
                     vxys(12) = clal(2)
                     vxys(13) = clal(3)

                  end if

                  if( vxys(11) .lt. -r0max .or.
     &                clmo .gt. -r0max ) vxys(11) = -2.0

                     tcb(1) = vxys( 8)
                     tcb(2) = vxys( 9)
                     tcb(3) = vxys(10)
                     tcl(1) = vxys(11)
                     tcl(2) = vxys(12)
                     tcl(3) = vxys(13)


                  if( ivs .eq. 0 ) then

                           tbl = tfs * bxwl
                           tbs = tfs * bxws
                           tbc = tfs * bxdl

                  else if( ivs .eq. 1 ) then

                           tbl = tfs * bxwl * 2.0
                           tbs = tfs * bxwl
                           tbc = tfs * bxdl * 2.0

                  end if


*-----------------------------------------------------------------------
*           TOTAL HIGHT
*-----------------------------------------------------------------------

                     strhall = 0.0

                     irdd = 0
                     rrdd = 0.0

                     if( abs(itth(0,1)) .gt. 1 ) irdd = irdd + 1

                  do 758 l = 1, ivr

                     strhall = strhall + tsrhh(l) - tsrbb(l)

                     if( abs(itth(l,1)) .gt. 1 ) irdd = irdd + 1

                     rrdd = rrdd + vspt(1,l) + vspt(2,l)

  758             continue

                     strhall = strhall
     &                       + dble(irdd) * tbc
     &                       + dble(ivr*2) * tbs
     &                       + tfs * rrdd

*-----------------------------------------------------------------------
*        TOTAL WIDTH
*-----------------------------------------------------------------------

                     strwall = 0.0

                     irdd = 0

                     if( abs(ittv(1,0)) .gt. 1 ) irdd = irdd + 1

                  do 759 m = 1, ivc

                     strwall = strwall + tsrll(m)

                     if( abs(ittv(1,m)) .gt. 1 ) irdd = irdd + 1

  759             continue

                     strwall = strwall
     &                       + dble(irdd) * tbc
     &                       + dble(ivc*2) * tbl

*-----------------------------------------------------------------------
*           STARTING POSITION
*-----------------------------------------------------------------------

                  xpsin = xpo
                  ypsin = ypo

*-----------------------------------------------------------------------
*           ROTATION
*-----------------------------------------------------------------------

               if( vxys(3) .lt. -r0max ) then

                  iang = 0

               else

                  iang = 1

                  if( idbg .eq. 1 )  write(jhf,'()')

                  write(jhf,'(3g14.5,'' rotin'')')
     &                        xpsin, ypsin, vxys(3)

                  xpsain = xpsin
                  ypsain = ypsin
                  angain = vxys(3)

                  call bbox(1,1,xpsain,ypsain,angain)

                  xpsin = 0.0
                  ypsin = 0.0

               end if

*-----------------------------------------------------------------------
*           DETERMINATION OF THE INITIAL LEFT TOP POSITION
*-----------------------------------------------------------------------

                  if( ivx .eq. 3 ) then

                     xpsin = xpsin - strwall

                  else if( ivx .eq. 2 ) then

                     xpsin = xpsin - strwall / 2.0

                  end if

                  if( ivy .eq. 1 ) then

                     ypsin = ypsin + strhall

                  else if( ivy .eq. 2 ) then

                     ypsin = ypsin + strhall / 2.0

                  end if

*-----------------------------------------------------------------------

                     bb1x = xpsin
                     bb2x = xpsin + strwall
                     bb1y = ypsin
                     bb2y = ypsin - strhall

                     call bbox(1,0,bb1x,bb1y,0.d0)
                     call bbox(1,0,bb1x,bb2y,0.d0)
                     call bbox(1,0,bb2x,bb1y,0.d0)
                     call bbox(1,0,bb2x,bb2y,0.d0)

*-----------------------------------------------------------------------
*           Y STRING STARTING POSITION
*-----------------------------------------------------------------------

                        typs(0) = ypsin

                  do 760 l = 1, ivr

                     if( abs(itth(l-1,1)) .gt. 1 ) then

                        irdd = 1

                     else

                        irdd = 0

                     end if

                        rrdd = vspt(1,l) + vspt(2,l)

                     if( l .eq. 1 ) then

                        ypsdn = tsrhh(l) + tbs
     &                        + dble(irdd) * tbc
     &                        + tfs * rrdd

                     else

                        ypsdn = tsrhh(l) + tbs * 2.0
     &                        + dble(irdd) * tbc
     &                        + tfs * rrdd
     &                        - tsrbb(l-1)

                     end if

                        typs(l)   = typs(l-1) - ypsdn

  760             continue

*-----------------------------------------------------------------------
*           X STRING STARTING POSITION
*-----------------------------------------------------------------------

               do 762 m = 1, ivc

                     if( abs(ittv(1,m-1)) .gt. 1 ) then

                        irdd = 1

                     else

                        irdd = 0

                     end if

                     if( m .eq. 1 ) then

                        txpm(m) = xpsin + tbl + dble(irdd) * tbc

                     else

                        txpm(m) = txpm(m-1)
     &                          + tsrll(m-1) + tbl * 2.0
     &                          + dble(irdd) * tbc

                     end if

                  do 761 l = 1, ivr

                     if( ittx(l,m) .eq. 1 ) then

                        txps(l,m) = txpm(m)

                     else if( ittx(l,m) .eq. 2 ) then

                        txps(l,m) = txpm(m) + tsrll(m) / 2.0
     &                            - tsrl(l,m) / 2.0

                     else if( ittx(l,m) .eq. 3 ) then

                        txps(l,m) = txpm(m) + tsrll(m)
     &                            - tsrl(l,m)

                     end if

  761             continue

  762          continue

*-----------------------------------------------------------------------
*           BACK GROUND OF TABLE
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(3f7.3,'' sc'')') tcb

                  if( abs(itth(0,1)) .gt. 1 ) then

                        write(jhf,'(2g14.5,'' M '',2g14.5,
     &                                     '' L'')')
     &                  xpsin,         ypsin-tbc,
     &                  xpsin+strwall, ypsin-tbc

                  else

                        write(jhf,'(2g14.5,'' M '',2g14.5,
     &                                     '' L'')')
     &                  xpsin,         ypsin,
     &                  xpsin+strwall, ypsin

                  end if

                     ilc = 1

               do 780 l = 1, ivr

                  if( abs(itth(l,1)) .gt. 1 .or. l .eq. ivr ) then

                        write(jhf,'(2g14.5,'' L '',2g14.5,
     &                                     '' L cp fl'')')
     &                  xpsin+strwall, typs(l)+tsrbb(l)-tbs,
     &                  xpsin,         typs(l)+tsrbb(l)-tbs

                  end if

                  if( abs(itth(l,1)) .gt. 1 .and. l .lt. ivr ) then

                           write(jhf,'(2g14.5,'' M '',2g14.5,
     &                                        '' L'')')
     &                     xpsin,         typs(l)+tsrbb(l)-tbs-tbc,
     &                     xpsin+strwall, typs(l)+tsrbb(l)-tbs-tbc

                  end if

  780          continue

*-----------------------------------------------------------------------
*           SEL COLOR OF TABLE
*-----------------------------------------------------------------------

               do 782 l = 1, ivr
               do 781 m = 1, ivc

                  if( csl(l,m,1) .lt. -r0max) goto 781

                     write(jhf,'(3f7.3,'' sc'')')
     &                  csl(l,m,1),csl(l,m,2),csl(l,m,3)

                  if( l .eq. 1 ) then

                     if( abs(itth(l-1,m)) .gt. 1 ) then

                        write(jhf,'(2g14.5,'' M '',2g14.5,
     &                                     '' L'')')
     &                  txpm(m)-tbl,          typs(l-1)-tbc,
     &                  txpm(m)+tsrll(m)+tbl, typs(l-1)-tbc

                     else

                        write(jhf,'(2g14.5,'' M '',2g14.5,
     &                                     '' L'')')
     &                  txpm(m)-tbl,          typs(l-1),
     &                  txpm(m)+tsrll(m)+tbl, typs(l-1)

                     end if

                  else

                     if( abs(itth(l-1,m)) .gt. 1 ) then

                        write(jhf,'(2g14.5,'' M '',2g14.5,
     &                                     '' L'')')
     &                  txpm(m)-tbl,
     &                  typs(l-1)+tsrbb(l-1)-tbs-tbc,
     &                  txpm(m)+tsrll(m)+tbl,
     &                  typs(l-1)+tsrbb(l-1)-tbs-tbc

                     else

                        write(jhf,'(2g14.5,'' M '',2g14.5,
     &                                     '' L'')')
     &                  txpm(m)-tbl,
     &                  typs(l-1)+tsrbb(l-1)-tbs,
     &                  txpm(m)+tsrll(m)+tbl,
     &                  typs(l-1)+tsrbb(l-1)-tbs

                     end if

                  end if

                        write(jhf,'(2g14.5,'' L '',2g14.5,
     &                                     '' L cp fl'')')
     &                  txpm(m)+tsrll(m)+tbl,
     &                  typs(l)+tsrbb(l)-tbs,
     &                  txpm(m)-tbl,
     &                  typs(l)+tsrbb(l)-tbs

  781          continue
  782          continue

*-----------------------------------------------------------------------
*           WRITE EACH STRING
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

               do 763 l = 1, ivr
               do 764 m = 1, ivc

                        write(jhf,'(''/xps '',g14.5,'' N '',
     &                              ''/yps '',g14.5,'' N ''
     &                              ''mtl'',2i3.3,'' mts'',2i3.3)')
     &           txps(l,m), typs(l) + tfs * vspt(1,l), l, m, l, m

  764          continue
  763          continue

*-----------------------------------------------------------------------
*           WRITE TABLE HORIZONTAL LINES
*-----------------------------------------------------------------------
*              CHECK OF FIRST AND LAST VERTICAL LINE WIDTH
*-----------------------------------------------------------------------

                     icfv = 0
                     iclv = 0

                     m = ivc

               do 770 l = 0, ivr

                  if( ittv(l,0) .gt. 0 .and. icfv .eq. 0 ) icfv = 1
                  if( ittv(l,0) .lt. 0 .and. icfv .le. 1 ) icfv = 2

                  if( ittv(l,m) .gt. 0 .and. iclv .eq. 0 ) iclv = 1
                  if( ittv(l,m) .lt. 0 .and. iclv .le. 1 ) iclv = 2

  770          continue

*-----------------------------------------------------------------------
*           LINE WIDTH AND SET COLOR, DASH AND LINE CAP
*-----------------------------------------------------------------------

            if( ittl .ne. 0 ) then

                     if( idbg .eq. 1 )  write(jhf,'()')

                     tds = tfs * bxls
                     tdl = tfs * bxls * 3.0

                     write(jhf,'(''sd0 lc0 '',3f7.3,'' sc'')') tcl

*-----------------------------------------------------------------------

               do 765 l = 0, ivr

                        ils  = 0
                        ilc  = 0
                        ild0 = 0

               do 766 m = 1, ivc

                        if( itth(l,m) .gt. 0 ) then

                           ild =  1

                        else if( itth(l,m) .lt. 0 ) then

                           ild = -1

                        else if( itth(l,m) .eq. 0 ) then

                           ild = 0

                        end if


                  if( itth(l,m) .ne. 0 .and. ils .eq. 0 ) then

                        if( ilc .eq. 0 ) then

                           if( l .eq. 0 ) then

                              ypsf = ypsin

                           else

                              ypsf = typs(l) + tsrbb(l) - tbs

                           end if

                           if( abs(itth(l,m)) .gt. 1 ) then

                              ypsd = ypsf - tbc

                           end if

                              ilc = 1

                        end if

                              ild0 = ild

                           if( ( l .eq. 0 .or. l .eq. ivr .or.
     &                           abs(itth(l,m)) .gt. 1 ) .and.
     &                           abs(ittv(l,m-1)) .gt. 1 ) then

                              xpss = txpm(m) - tbl - tbc

                           else

                              xpss = txpm(m) - tbl

                           end if

                           if( m .eq. 1 .and. icfv .ne. 0 ) then

                              if( icfv .eq. 1 ) then

                                 xpss = xpss - tds / 2.0

                              else if( icfv .eq. 2 ) then

                                 xpss = xpss - tdl / 2.0

                              end if

                           end if

                              ils = 1

                  end if


                  if( ils .eq. 1 ) then


                     if( m .eq. ivc ) then

                        if( ild .eq. ild0 ) then

                              mm = m

                        else if( ild .ne. ild0 ) then

                              mm = m - 1

                        end if

                           if( ( l .eq. 0 .or. l .eq. ivr .or.
     &                           abs(itth(l,m)) .gt. 1 ) .and.
     &                           abs(ittv(l,m)) .gt. 1 ) then

                              xpsf = txpm(mm) + tsrll(mm) + tbl + tbc

                           else

                              xpsf = txpm(mm) + tsrll(mm) + tbl

                           end if

                           if( iclv .ne. 0 ) then

                              if( iclv .eq. 1 ) then

                                 xpsf = xpsf + tds / 2.0

                              else if( iclv .eq. 2 ) then

                                 xpsf = xpsf + tdl / 2.0

                              end if

                           end if


                           if( ild0 .eq. 1 ) then

                              write(jhf,'(g14.5,'' lw'')') tds

                           else

                              write(jhf,'(g14.5,'' lw'')') tdl

                           end if

                              write(jhf,'(2g14.5,'' M '',
     &                                    2g14.5,'' L st'')')
     &                        xpss, ypsf, xpsf, ypsf


                           if( abs(itth(l,mm)) .gt. 1 ) then

                              write(jhf,'(2g14.5,'' M '',
     &                                    2g14.5,'' L st'')')
     &                        xpss, ypsd, xpsf, ypsd

                           end if


                        if( ild .ne. ild0 .and. ild .ne. 0 ) then

                              xpss = txpm(m) - tbl

                           if( ( l .eq. 0 .or. l .eq. ivr ) .and.
     &                           abs(ittv(l,m)) .gt. 1 ) then

                              xpss = txpm(m) + tsrll(m) + tbl + tbc

                           else

                              xpss = txpm(m) + tsrll(m) + tbl

                           end if

                           if( iclv .ne. 0 ) then

                              if( iclv .eq. 1 ) then

                                 xpsf = xpsf + tds / 2.0

                              else if( iclv .eq. 2 ) then

                                 xpsf = xpsf + tdl / 2.0

                              end if

                           end if


                           if( ild .eq. 1 ) then

                              write(jhf,'(g14.5,'' lw'')') tds

                           else

                              write(jhf,'(g14.5,'' lw'')') tdl

                           end if

                              write(jhf,'(2g14.5,'' M '',
     &                                    2g14.5,'' L st'')')
     &                        xpss, ypsf, xpsf, ypsf

                        end if


                     else if( ild .ne. ild0 ) then

                              xpsf = txpm(m-1) + tsrll(m-1) + tbl

                           if( ild0 .eq. 1 ) then

                              write(jhf,'(g14.5,'' lw'')') tds

                           else

                              write(jhf,'(g14.5,'' lw'')') tdl

                           end if

                              write(jhf,'(2g14.5,'' M '',
     &                                    2g14.5,'' L st'')')
     &                        xpss, ypsf, xpsf, ypsf

                        if( ild .ne. 0 ) then

                              xpss = txpm(m) - tbl

                        else if( ild .eq. 0 ) then

                              ils = 0

                        end if

                              ild0 = ild

                     end if


                  end if


  766          continue

  765          continue


*-----------------------------------------------------------------------
*           WRITE TABLE VERTICAL LINES
*-----------------------------------------------------------------------


               do 767 m = 0, ivc

                        ims  = 0
                        imc  = 0
                        imd0 = 0

               do 768 l = 1, ivr

                        if( ittv(l,m) .gt. 0 ) then

                           imd =  1

                        else if( ittv(l,m) .lt. 0 ) then

                           imd = -1

                        else if( ittv(l,m) .eq. 0 ) then

                           imd = 0

                        end if


                  if( imd .ne. 0 .and. ims .eq. 0 ) then

                        if( imc .eq. 0 ) then

                           if( m .eq. 0 ) then

                              xpsf = xpsin

                           else

                              xpsf = txpm(m) + tsrll(m) + tbl

                           end if

                           if( abs(ittv(l,m)) .gt. 1 ) then

                              xpsd = xpsf + tbc

                           end if

                              imc = 1

                        end if

                              imd0 = imd

                        if( l .gt. 1 ) then

                           if( abs(itth(l-1,m)) .gt. 1 ) then

                              ypss = typs(l-1) + tsrbb(l-1) - tbs - tbc

                           else

                              ypss = typs(l-1) + tsrbb(l-1) - tbs

                           end if

                        else

                           if( abs(itth(l-1,m)) .gt. 1 ) then

                              ypss = ypsin - tbc

                           else

                              ypss = ypsin

                           end if

                        end if

                              ims = 1


                  else if( ims .eq. 1 .and.
     &                   ( imd .ne. imd0 .or.
     &                   ( l .gt. 1 .and.
     &                     abs(itth(l-1,m)) .gt. 1 ) ) ) then

                              ypsf = typs(l-1) + tsrbb(l-1) - tbs

                           if( imd0 .eq. 1 ) then

                              write(jhf,'(g14.5,'' lw'')') tds

                           else

                              write(jhf,'(g14.5,'' lw'')') tdl

                           end if

                              write(jhf,'(2g14.5,'' M '',
     &                                    2g14.5,'' L st'')')
     &                        xpsf, ypss, xpsf, ypsf

                           if( abs(ittv(l,m)) .gt. 1 ) then

                              write(jhf,'(2g14.5,'' M '',
     &                                    2g14.5,'' L st'')')
     &                        xpsd, ypss, xpsd, ypsf

                           end if

                        if( imd .ne. 0 ) then

                           if( abs(itth(l-1,m)) .gt. 1 ) then

                              ypss = typs(l-1) + tsrbb(l-1) - tbs - tbc

                           else

                              ypss = typs(l-1) + tsrbb(l-1) - tbs

                           end if

                        else if( imd .eq. 0 ) then

                              ims = 0

                        end if

                              imd0 = imd

                  end if


                  if( l .eq. ivr .and. ims .eq. 1 ) then

                              ypsf = typs(l) + tsrbb(l) - tbs

                        if( imd0 .eq. 1 ) then

                              write(jhf,'(g14.5,'' lw'')') tds

                        else

                              write(jhf,'(g14.5,'' lw'')') tdl

                        end if

                              write(jhf,'(2g14.5,'' M '',
     &                                    2g14.5,'' L st'')')
     &                        xpsf, ypss, xpsf, ypsf

                        if( abs(ittv(l,m)) .gt. 1 ) then

                              write(jhf,'(2g14.5,'' M '',
     &                                    2g14.5,'' L st'')')
     &                        xpsd, ypss, xpsd, ypsf

                        end if

                  end if

  768          continue

  767          continue



                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(''lc1'')')

            end if

*-----------------------------------------------------------------------

                  if( iang .eq. 1 ) then

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(''rotout'')')

                     call bbox(1,1,0.d0,0.d0,-angain)
                     call bbox(1,1,-xpsain,-ypsain,0.d0)

                  end if

*-----------------------------------------------------------------------

  750   continue

*-----------------------------------------------------------------------

      deallocate(chag)
      return
      end


************************************************************************
*                                                                      *
      subroutine wparn(jhf,idbg,pacn,rstx,
     &                 xmin,xmax,ixlog,ymin,ymax,iylog,ifon,
     &                 ibfon,fscm,clal,clmo,
     &                 regx,regy,regs,ilbox,cllg,clgb,clgl,clgs)
*                                                                      *
*                                                                      *
*      PURPOSE  :   WRITE VALUES OF THE PARAMETTERS AND CONSTANTS      *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter ( regyd0 = 0.9 )

      common /con/  cm, dd
      common /frm/  xal, yal

      common /wtval1/ strl0, strh0, strb0

      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb

      common /box/  bxws, bxwl, bxss, bxsl, bxds, bxdl, bxls, bxll

      character ccmt(ichrl)*1
      character pacn*1

      common /rval1/ cval(mxcval), aval(mxcval)

      dimension bval(mxcval)

      dimension iaak(mxcval)

      character ycm(mxcval)*20

      character,allocatable:: chag(:,:)
      dimension clal(3), cllg(3), clgb(3), clgl(3), clgs(3)
      dimension bcb(3), bcl(3), bcs(3)

*-----------------------------------------------------------------------

      allocate(chag(inig,0:ichrl))

*-----------------------------------------------------------------------
*        INITIAL PROCEDURE
*-----------------------------------------------------------------------

            do i = 1, mxcval

               if( pacn .eq. 'A' ) then

                  bval(i) = aval(i)

               else

                  bval(i) = cval(i)

               end if

            end do


               ncom = 0

            do 100 i = 1, mxcval

               if( bval(i) .gt. -r0max ) then

                  ncom = ncom + 1

                  iaak(ncom) = i

               end if

  100       continue

            if( ncom .eq. 0 ) goto 999

*-----------------------------------------------------------------------
*        REAL NUMBER TO CHARACTERS
*-----------------------------------------------------------------------

            do 110 i = 1, ncom

               write( ycm(i), '(a1,i2.2,'' = '',g14.6)' )
     &         pacn, iaak(i), bval(iaak(i))

  110       continue

*-----------------------------------------------------------------------
*        INITIAL POSITION
*-----------------------------------------------------------------------

                  if( regx .lt. -r0max ) then

                        regx = rstx

                  else

                        regx = zonep(xmin,xmax,regx,ixlog)

                  end if

                        regyd = regyd0 * regs / yal

                  if( regy .lt. -r0max ) then

                        regy = 1.0 + dble( ncom ) * regyd

                  else

                        regy = zonep(ymin,ymax,regy,iylog)

                  end if


                  if( ilbox .ne. 0 ) then

                        rebl = regx * xal * cm
                        rebr = rebl
                        rebu = regy * yal * cm
                        rebd = ( regy - dble( ncom - 1 ) * regyd )
     &                       * yal * cm

                  end if

*-----------------------------------------------------------------------
*     WRITE VALUE OF THE PARAMETERS
*-----------------------------------------------------------------------

                  if( idbg .eq. 1 )
     &            write(jhf,'(/''%'',71(''-''),/
     &                         ''% Write Values of the Parameters''
     &                       ,/''%'',71(''-'')/)')

*-----------------------------------------------------------------------

               do 455 i = 1, ncom

                     iflgc = ibfon

                     tfs = fscm * regs

                     if( cllg(1) .lt. -r0max ) cllg(1) = -2.0


                     xpo = regx * xal * cm
                     ypo = ( regy - dble( i - 1 ) * regyd ) * yal * cm


                     ixc  = 1
                     iyc  = 2

                  do 456 k = 1, 20

                     ccmt(k) = ycm(i)(k:k)

  456             continue


                  if( ilbox .eq. 0 ) then

                     call wtext(jhf,chag,idbg,clal,clmo,0,xpo,ypo,
     &                          ccmt,20,ixc,iyc,0,0.d0,iflgc,ifon,
     &                          0,indt,tfs,cllg)

                  else

                     iccb = 1

                     call wrbox(jhf,idbg,0,ilbox,iccb,
     &                          b1x, b2x, b1y, b2y,
     &                          bdd, bcb, bcl, bdc, bcs)

                     ilmul = i

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(''/mtl'',i3.3,'' {'')') ilmul

                     call wtext(jhf,chag,idbg,clal,clmo,-2,xpo,ypo,
     &                          ccmt,20,ixc,iyc,0,0.d0,iflgc,ifon,
     &                          ilmul,indt,tfs,cllg)

                        rebr = max(rebr, xpo + strl0 )
                        rebu = max(rebu, ypo + strh0 / 2.0 )
                        rebd = min(rebd, ypo - strh0 / 2.0 )

                  end if

                     if( idbg .eq. 1 )  write(jhf,'()')


  455          continue

*-----------------------------------------------------------------------
*           WRITE PARAMETERS BOX
*-----------------------------------------------------------------------

               if( ilbox .ne. 0 ) then

                  if( idbg .eq. 1 )
     &            write(jhf,'(/''%'',71(''-''),/
     &                         ''% Parameters Box''
     &                       ,/''%'',71(''-'')/)')

                        if( clgb(1) .lt. -r0max .or.
     &                      clmo .gt. -r0max .or.
     &                      clal(1) .gt. -r0max ) clgb(1) = -1.0

                        if( clal(1) .gt. -r0max ) then

                           clgl(1) = clal(1)
                           clgl(2) = clal(2)
                           clgl(3) = clal(3)
                           clgs(1) = clal(1)
                           clgs(2) = clal(2)
                           clgs(3) = clal(3)

                        end if

                        if( clgl(1) .lt. -r0max .or.
     &                      clmo .gt. -r0max ) clgl(1) = -2.0

                        if( clgs(1) .lt. -r0max .or.
     &                      clmo .gt. -r0max ) clgs(1) = -2.0

                        bcb(1) = clgb(1)
                        bcb(2) = clgb(2)
                        bcb(3) = clgb(3)
                        bcl(1) = clgl(1)
                        bcl(2) = clgl(2)
                        bcl(3) = clgl(3)
                        bcs(1) = clgs(1)
                        bcs(2) = clgs(2)
                        bcs(3) = clgs(3)


                  if( ilbox .eq. 10 .or. ilbox .eq. 11. or.
     &                ilbox .eq. 20 .or. ilbox .eq. 21. or.
     &                ilbox .eq. 30 .or. ilbox .eq. 31. or.
     &                ilbox .eq. 40 .or. ilbox .eq. 41. or.
     &                ilbox .eq. 50 .or. ilbox .eq. 51. ) then

                        bds = tfs * bxwl

                  else

                        bds = tfs * bxwl * 2.0

                  end if

                  if( ilbox .eq. 11 .or. ilbox .eq. 13. or.
     &                ilbox .eq. 21 .or. ilbox .eq. 23. or.
     &                ilbox .eq. 31 .or. ilbox .eq. 33. or.
     &                ilbox .eq. 41 .or. ilbox .eq. 43. or.
     &                ilbox .eq. 51 .or. ilbox .eq. 53. ) then

                        bdd = tfs * bxll

                  else

                        bdd = tfs * bxls

                  end if

                  if( ilbox .eq. 20 .or. ilbox .eq. 21. or.
     &                ilbox .eq. 40 .or. ilbox .eq. 41. or.
     &                ilbox .eq. 50 .or. ilbox .eq. 51. ) then

                        bdc = tfs * bxsl

                  else
     &            if( ilbox .eq. 22 .or. ilbox .eq. 23. or.
     &                ilbox .eq. 42 .or. ilbox .eq. 43. or.
     &                ilbox .eq. 52 .or. ilbox .eq. 53. ) then

                        bdc = tfs * bxsl * 2.0

                  end if

                  if( ilbox .eq. 30 .or. ilbox .eq. 31. ) then

                        bdc = tfs * bxdl

                  else
     &            if( ilbox .eq. 32 .or. ilbox .eq. 33. ) then

                        bdc = tfs * bxdl * 2.0

                  end if

                     b1x = rebl - bds
                     b2x = rebr + bds

                     b1y = rebd - bds
                     b2y = rebu + bds

*-----------------------------------------------------------------------

                     iccb = 0

                     call wrbox(jhf,idbg,0,ilbox,iccb,
     &                          b1x, b2x, b1y, b2y,
     &                          bdd, bcb, bcl, bdc, bcs)

*-----------------------------------------------------------------------

                  if( ilbox .ge. 40 ) then

                     b2x = b2x + bdc
                     b1y = b1y - bdc

                  end if

                     call bbox(1,0,b1x,b1y,0.d0)
                     call bbox(1,0,b1x,b2y,0.d0)
                     call bbox(1,0,b2x,b1y,0.d0)
                     call bbox(1,0,b2x,b2y,0.d0)

*-----------------------------------------------------------------------

                  do 457 i = 1, ncom

                     write(jhf,'(''mtl'',i3.3,'' mts'',i3.3)') i, i

  457             continue

               end if

*-----------------------------------------------------------------------

  999 deallocate(chag)
      return
      end

************************************************************************
*                                                                      *
      subroutine winclp(jhf,idbg,inps,
     &                  dpin,idpi,
     &                  psbbx,xyps,ixps,iyps,
     &                  xmax,xmin,ymax,ymin,ixlog,iylog)
*                                                                      *
*        PURPOSE    :  WRITE INCLUDE PS FILE                           *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      common /con/  cm, dd
      common /frm/  xal, yal

      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb


      character dpin(ipsm)*200
      dimension idpi(ipsm)
      dimension psbbx(ipsm,4)
      dimension xyps(ipsm,6)
      dimension ixps(ipsm), iyps(ipsm)

      character dps(ichrl)*1

      data ipsd / 0 /

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------
*     INCLUDE PS FILE
*-----------------------------------------------------------------------

            if( idbg .eq. 1 )
     &      write(jhf,'(/''%'',71(''-''),/
     &                   ''% Include PS File''
     &                 ,/''%'',71(''-'')/)')


         do 800 i = 1, inps

*-----------------------------------------------------------------------

                     if( xyps(i,1) .lt. -r0max ) then

                        xyps(i,1) = xmin

                        ixps(i) = 1

                     end if

                     if( xyps(i,2) .lt. -r0max ) then

                        xyps(i,2) = ymax

                        iyps(i) = 1

                     end if


                  xpo = zonep(xmin,xmax,xyps(i,1),ixlog) * xal * cm
                  ypo = zonep(ymin,ymax,xyps(i,2),iylog) * yal * cm


                  if( xyps(i,3) .lt. -r0max ) then

                     iang = 0

                     write(jhf,'(2g14.5,'' TR'')') xpo, ypo

                     xpsain = xpo
                     ypsain = ypo
                     angain = 0.0

                     call bbox(1,1,xpsain,ypsain,angain)

                  else

                     iang = 1

                     write(jhf,'(2g14.5,'' TR'')') xpo, ypo

                     write(jhf,'(g14.5,'' rotate'')') xyps(i,3)

                     xpsain = xpo
                     ypsain = ypo
                     angain = xyps(i,3)

                     call bbox(1,1,xpsain,ypsain,angain)

                  end if

                  if( xyps(i,4) .lt. -r0max ) xyps(i,4) = 1.0

                  if( xyps(i,5) .lt. -r0max ) xyps(i,5) = xyps(i,4)
                  if( xyps(i,6) .lt. -r0max ) xyps(i,6) = xyps(i,4)

*-----------------------------------------------------------------------
*           SCALE
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' scale'')')
     &                  xyps(i,5), xyps(i,6)

                     sxg0 = sxg0 * xyps(i,5)
                     syg0 = syg0 * xyps(i,6)

*-----------------------------------------------------------------------
*           DETERMIN THE LEFT DOWN CONNER AND TRANSLATE
*-----------------------------------------------------------------------

                     xdis = psbbx(i,3) - psbbx(i,1)
                     ydis = psbbx(i,4) - psbbx(i,2)

                  if( ixps(i) .eq. 1 ) then

                     xpo = - psbbx(i,1)

                  else if( ixps(i) .eq. 2 ) then

                     xpo = - xdis / 2.0 - psbbx(i,1)

                  else if( ixps(i) .eq. 3 ) then

                     xpo = - xdis - psbbx(i,1)

                  end if

                  if( iyps(i) .eq. 1 ) then

                     ypo = - psbbx(i,2)

                  else if( iyps(i) .eq. 2 ) then

                     ypo = - ydis / 2.0 - psbbx(i,2)

                  else if( iyps(i) .eq. 3 ) then

                     ypo = - ydis - psbbx(i,2)

                  end if


                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' TR'')') xpo, ypo

                     xpsbin = xpo
                     ypsbin = ypo

                     call bbox(1,1,xpsbin,ypsbin,0.d0)


*-----------------------------------------------------------------------
*           BOUNDING BOX
*-----------------------------------------------------------------------

                     xpsm = psbbx(i,3)
                     ypsm = psbbx(i,4)

                     call bbox(1,0,0.d0 ,0.d0 ,0.d0)
                     call bbox(1,0,xpsm,0.d0 ,0.d0)
                     call bbox(1,0,0.d0 ,ypsm,0.d0)
                     call bbox(1,0,xpsm,ypsm,0.d0)

*-----------------------------------------------------------------------
*           STARTING DESCRIPTION OF INCLUDING PS FILE
*-----------------------------------------------------------------------

                  if( ipsd .eq. 0 ) then

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(''/workdict 12 dict N'')')

                     ipsd = 1

                  end if

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(
     &               ''workdict begin /Save_state save N '',
     &               ''gs count /pcount S N''/
     &               ''/qcount countdictstack N'',
     &               ''/showpage{}N /erasepage{}N /copypage{}N''/
     &               ''sd0 -2 1 1 sc lc1 lj0'')')

                     write(jhf,'(/
     &               ''%%BeginDocument: '',200A1/)')
     &               ( dpin(i)(k:k), k = 1, idpi(i) )

*-----------------------------------------------------------------------
*           OPEN INCLUDE PS FILE AND COPY IT ON jhf EXCEPT FOR
*           %% LINES, %! LINES
*-----------------------------------------------------------------------

               open(jps, file = dpin(i), status = 'OLD' )

*-----------------------------------------------------------------------

                  il = 0

  840             il = il + 1

                  read(jps,'(10000a1)',iostat=ios) (dps(ic),ic=1,ichrl)
                  if( ios .eq. -1 ) goto 841

                  if( ichar(dps(1)) .eq. 26 ) goto 840

*-----------------------------------------------------------------------
*              SKIP %% LINES, %! LINES
*-----------------------------------------------------------------------

                  if( dps(1) .eq. '%' .and. dps(2) .eq. '%' ) goto 840
                  if( dps(1) .eq. '%' .and. dps(2) .eq. '!' ) goto 840

*-----------------------------------------------------------------------
*              LAST COLUMN
*-----------------------------------------------------------------------

                  do 865 l = ichrl, 1, -1

                     if(dps(l).ne.' ' .and. dps(l).ne.tub.and.
     &                  ichar(dps(l)) .ne. 26 )  goto 866

  865             continue

  866                kl = l

                     write(jhf,'(1024A1)') (dps(ic),ic=1,kl)

                     goto 840

*-----------------------------------------------------------------------
*           END OF INCLUDE PS FILE
*-----------------------------------------------------------------------

  841             continue

                     write(jhf,'(/
     &               ''%%EndDocument'')')

*-----------------------------------------------------------------------
*           CLOSE INCLUDE PS FILE
*-----------------------------------------------------------------------

                     close(jps)

*-----------------------------------------------------------------------
*           ENDING DESCRIPTION OF INCLUDE PS FILE
*-----------------------------------------------------------------------

                     write(jhf,'(/
     &               ''count pcount sub {pop} repeat countdictstack '',
     &               ''qcount sub {end} repeat gr''/
     &               ''Save_state restore end'')')


*-----------------------------------------------------------------------
*           SCALE, ROTATE AND TRANSLATE BACK TO THE ORIGINAL
*-----------------------------------------------------------------------

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' TR'')') -xpsbin, -ypsbin

                     call bbox(1,1,-xpsbin,-ypsbin,0.d0)


                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' scale'')')
     &                  1./xyps(i,5), 1./xyps(i,6)

                     sxg0 = sxg0 / xyps(i,5)
                     syg0 = syg0 / xyps(i,6)


                  if( iang .eq. 0 ) then

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(2g14.5,'' TR'')') -xpsain, -ypsain

                     call bbox(1,1,-xpsain,-ypsain,0.d0)

                  else

                     if( idbg .eq. 1 )  write(jhf,'()')

                     write(jhf,'(g14.5,'' rotate'')') -angain

                     write(jhf,'(2g14.5,'' TR'')') -xpsain, -ypsain

                     call bbox(1,1,0.d0,0.d0,-angain)
                     call bbox(1,1,-xpsain,-ypsain,0.d0)

                  end if

*-----------------------------------------------------------------------

  800    continue

                     if( idbg .eq. 1 )  write(jhf,'()')


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine rmarin(ij,kl)
*                                                                      *
* This is the initialization routine for the random number generator   *
*     ranmar()                                                         *
* note: The seed variables can have values between: 0 <= IJ <= 31328   *
*                                                   0 <= KL <= 30081   *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      dimension U(97)
      integer I97, J97

      common /raset1/ U, C, CD, CM, I97, J97

         if( ij .lt. 0 .or. kl .lt. 0 ) stop 999

  100    if( ij .gt. 31328 ) then
            ij = ij / 2
            goto 100
         end if

  200    if( kl .gt. 30081 ) then
            kl = kl / 2
            goto 200
         end if

      i = mod(IJ/177, 177) + 2
      j = mod(IJ    , 177) + 2
      k = mod(KL/169, 178) + 1
      l = mod(KL,     169)
      do 2 ii = 1, 97
         s = 0.0
         t = 0.5
         do 3 jj = 1, 24
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) .ge. 32) then
               s = s + t
            endif
            t = 0.5 * t
3        continue
         U(ii) = s
2     continue
      C = 362436.0 / 16777216.0
      CD = 7654321.0 / 16777216.0
      CM = 16777213.0 /16777216.0
      I97 = 97
      J97 = 33
      return
      end

      function ranmar()

      implicit real*8 (a-h,o-z)

      dimension U(97)
      integer I97, J97
      common /raset1/ U, C, CD, CM, I97, J97

         uni = U(I97) - U(J97)
         if( uni .lt. 0.0 ) uni = uni + 1.0
         U(I97) = uni
         I97 = I97 - 1
         if(I97 .eq. 0) I97 = 97
         J97 = J97 - 1
         if(J97 .eq. 0) J97 = 97
         C = C - CD
         if( C .lt. 0.0 ) C = C + CM
         uni = uni - C
         if( uni .lt. 0.0 ) uni = uni + 1.0
         RANMAR = uni

      return
      END

************************************************************************
*                                                                      *
      subroutine inpdat(indat,infn,iname,iexp,
     &                  iangb,iver2,inpage,ifpage,
     &                  iofn,ioname,ierr)
*                                                                      *
*       read input information                                         *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character infn(ichrl)*1
      character iofn(ichrl)*1
      character indat(ichrl)*1

*-----------------------------------------------------------------------
*     default values
*-----------------------------------------------------------------------

            iname = 0
            iangb = 0
            iver2 = 0
            iexp  = 0
            inpage = -1
            ifpage = -1
            ioname = 0

*-----------------------------------------------------------------------

         do i = icolm, 1, -1

            if( indat(i) .ne. ' ' ) goto 50

         end do

            goto 999

   50    kmax = i + 1

*-----------------------------------------------------------------------

            k = 0

  100    continue

            k = k + 1

            if( k .gt. kmax ) goto 200

         if( indat(k) .eq. ' ' ) goto 100

         if( indat(k) .eq. '-' .and. indat(k+1) .eq. 'a' ) then

            iangb = 1
            k = k + 1
            goto 100

         else if( indat(k  ) .eq. '-' .and. indat(k+1) .eq. 'v' .and.
     &            indat(k+2) .eq. '1' .and. indat(k+3) .eq. '2' ) then

            iver2 = 1
            k = k + 3
            goto 100

         else if( indat(k  ) .eq. '-' .and. indat(k+1) .eq. 'e' .and.
     &            indat(k+2) .eq. 'x' .and. indat(k+3) .eq. 'p' ) then

            iexp = 1
            k = k + 3
            goto 100

         else if( indat(k  ) .eq. '-' .and. indat(k+1) .eq. 'f' ) then

               ini = -1
               ifi = -1

               k = k + 1
  101          k = k + 1

               if( k .gt. kmax ) goto 200

            if( ini .eq. -1 .and. indat(k) .eq. ' ' ) then

               goto 101

            else if( ini .eq. -1 .and. indat(k) .ne. ' ' ) then

               ini = k
               goto 101

            else if( ini .gt. 0 .and. indat(k) .ne. ' ' ) then

               goto 101

            else if( ini .gt. 0 .and. indat(k) .eq. ' ' ) then

               ifi = k - 1

               call rnum(rnm,indat,ini,ifi,ierr)

               if( ierr .ne. 0 ) goto 999

               inpage = nint( rnm )

               if( inpage .le. 0 ) goto 999

               goto 100

            end if

         else if( indat(k  ) .eq. '-' .and. indat(k+1) .eq. 't' ) then

               ini = -1
               ifi = -1

               k = k + 1
  102          k = k + 1

               if( k .gt. kmax ) goto 200

            if( ini .eq. -1 .and. indat(k) .eq. ' ' ) then

               goto 102

            else if( ini .eq. -1 .and. indat(k) .ne. ' ' ) then

               ini = k
               goto 102

            else if( ini .gt. 0 .and. indat(k) .ne. ' ' ) then

               goto 102

            else if( ini .gt. 0 .and. indat(k) .eq. ' ' ) then

               ifi = k - 1

               call rnum(rnm,indat,ini,ifi,ierr)

               if( ierr .ne. 0 ) goto 999

               ifpage = nint( rnm )

               if( ifpage .le. 0 ) goto 999

               goto 100

            end if

         else if( indat(k  ) .eq. '-' .and. indat(k+1) .eq. 'o' ) then

               ini = -1
               ifi = -1

               k = k + 1
  103          k = k + 1

               if( k .gt. kmax ) goto 200

            if( ini .eq. -1 .and. indat(k) .eq. ' ' ) then

               goto 103

            else if( ini .eq. -1 .and. indat(k) .ne. ' ' ) then

               ini = k
               goto 103

            else if( ini .gt. 0 .and. indat(k) .ne. ' ' ) then

               goto 103

            else if( ini .gt. 0 .and. indat(k) .eq. ' ' ) then

               ifi = k - 1

               ioname = ifi - ini + 1

               do i = 1, ioname

                  iofn(i) = indat(i+ini-1)

               end do

               goto 100

            end if

         else

               ini = k

               k = k - 1
  104          k = k + 1

               if( k .gt. kmax ) goto 200

            if( indat(k) .ne. ' ' ) then

               goto 104

            else if( indat(k) .eq. ' ' ) then

               ifi = k - 1

               iname = ifi - ini + 1

               do i = 1, iname

                  infn(i) = indat(i+ini-1)

               end do

               goto 100

            end if

         end if

*-----------------------------------------------------------------------

  200    continue

         ierr = 0
         return

*-----------------------------------------------------------------------

  999    continue

         ierr = 1

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine applyColorMap(cmap,ndis,icrev,rcol,rin0)
*     cmap (= colormap)
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)
      dimension rcol(3)
      character cmap*99

C     Hue[0,360], Saturation[0,100], Brightness[0,100]

      integer, parameter :: n_phits2 = 11
      real phits2(3,n_phits2)
      data phits2 /
     &  294, 71,  34,
     &  317, 100, 58,
     &  350, 100, 100,
     &  41,  100, 100,
     &  60,  100, 100,
     &  90,  100, 100,
     &  147, 91,  81,
     &  182, 100, 100,
     &  230, 100, 100,
     &  233, 52,  91,
     &  195, 13,  100/

! This software includes colormap data extracted from Matplotlib,
! which is licensed under the 3-Clause BSD License.

! ------------------ !
! Matplotlib License !
! ------------------ !

! License agreement for matplotlib versions 1.3.0 and later
! =========================================================

! 1. This LICENSE AGREEMENT is between the Matplotlib Development Team
! ("MDT"), and the Individual or Organization ("Licensee") accessing and
! otherwise using matplotlib software in source or binary form and its
! associated documentation.

! 2. Subject to the terms and conditions of this License Agreement, MDT
! hereby grants Licensee a nonexclusive, royalty-free, world-wide license
! to reproduce, analyze, test, perform and/or display publicly, prepare
! derivative works, distribute, and otherwise use matplotlib
! alone or in any derivative version, provided, however, that MDT's
! License Agreement and MDT's notice of copyright, i.e., "Copyright (c)
! 2012- Matplotlib Development Team; All Rights Reserved" are retained in
! matplotlib  alone or in any derivative version prepared by
! Licensee.

! 3. In the event Licensee prepares a derivative work that is based on or
! incorporates matplotlib or any part thereof, and wants to
! make the derivative work available to others as provided herein, then
! Licensee hereby agrees to include in any such work a brief summary of
! the changes made to matplotlib .

! 4. MDT is making matplotlib available to Licensee on an "AS
! IS" basis.  MDT MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR
! IMPLIED.  BY WAY OF EXAMPLE, BUT NOT LIMITATION, MDT MAKES NO AND
! DISCLAIMS ANY REPRESENTATION OR WARRANTY OF MERCHANTABILITY OR FITNESS
! FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF MATPLOTLIB
! WILL NOT INFRINGE ANY THIRD PARTY RIGHTS.

! 5. MDT SHALL NOT BE LIABLE TO LICENSEE OR ANY OTHER USERS OF MATPLOTLIB
!  FOR ANY INCIDENTAL, SPECIAL, OR CONSEQUENTIAL DAMAGES OR
! LOSS AS A RESULT OF MODIFYING, DISTRIBUTING, OR OTHERWISE USING
! MATPLOTLIB , OR ANY DERIVATIVE THEREOF, EVEN IF ADVISED OF
! THE POSSIBILITY THEREOF.

! 6. This License Agreement will automatically terminate upon a material
! breach of its terms and conditions.

! 7. Nothing in this License Agreement shall be deemed to create any
! relationship of agency, partnership, or joint venture between MDT and
! Licensee.  This License Agreement does not grant permission to use MDT
! trademarks or trade name in a trademark sense to endorse or promote
! products or services of Licensee, or any third party.

! 8. By copying, installing or otherwise using matplotlib ,
! Licensee agrees to be bound by the terms and conditions of this License
! Agreement.

! ----------------------------------------------------------------------
! For more information, see:
!   https://matplotlib.org/stable/project/license.html

! ----------- !
! ColorBrewer !
! ----------- !
! Apache-Style Software License for ColorBrewer software and ColorBrewer Color Schemes

! Copyright (c) 2002 Cynthia Brewer, Mark Harrower, and The Pennsylvania State University.

! Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
! You may obtain a copy of the License at

! http://www.apache.org/licenses/LICENSE-2.0

! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.

! https://github.com/altercation/solarized/blob/master/LICENSE
! Copyright (c) 2011 Ethan Schoonover

! --------- !
! Solarized !
! --------- !
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.

! ------ !
! Yorick !
! ------ !
! BSD-style license for gist/yorick colormaps.

! Copyright:

!   Copyright (c) 1996.  The Regents of the University of California.
!              All rights reserved.

! Permission to use, copy, modify, and distribute this software for any
! purpose without fee is hereby granted, provided that this entire
! notice is included in all copies of any software which is or includes
! a copy or modification of this software and in all copies of the
! supporting documentation for such software.

! This work was produced at the University of California, Lawrence
! Livermore National Laboratory under contract no. W-7405-ENG-48 between
! the U.S. Department of Energy and The Regents of the University of
! California for the operation of UC LLNL.


!                   DISCLAIMER

! This software was prepared as an account of work sponsored by an
! agency of the United States Government.  Neither the United States
! Government nor the University of California nor any of their
! employees, makes any warranty, express or implied, or assumes any
! liability or responsibility for the accuracy, completeness, or
! usefulness of any information, apparatus, product, or process
! disclosed, or represents that its use would not infringe
! privately-owned rights.  Reference herein to any specific commercial
! products, process, or service by trade name, trademark, manufacturer,
! or otherwise, does not necessarily constitute or imply its
! endorsement, recommendation, or favoring by the United States
! Government or the University of California.  The views and opinions of
! authors expressed herein do not necessarily state or reflect those of
! the United States Government or the University of California, and
! shall not be used for advertising or product endorsement purposes.


!                 AUTHOR

! David H. Munro wrote Yorick and Gist.  Berkeley Yacc (byacc) generated
! the Yorick parser.  The routines in Math are from LAPACK and FFTPACK;
! MathC contains C translations by David H. Munro.  The algorithms for
! Yorick's random number generator and several special functions in
! Yorick/include were taken from Numerical Recipes by Press, et. al.,
! although the Yorick implementations are unrelated to those in
! Numerical Recipes.  A small amount of code in Gist was adapted from
! the X11R4 release, copyright M.I.T. -- the complete copyright notice
! may be found in the (unused) file Gist/host.c.

      integer, parameter :: n_matp1 = 20
      real matplotlib1(3,n_matp1,39) ! (HSB,i,imatp)
      data matplotlib1 /
C      imatp1=1: viridis
     &  53, 85, 99,      61, 89, 88,      71, 82, 87,      86, 70, 84,
     & 104, 59, 81,     128, 56, 77,     146, 67, 73,     158, 76, 68,
     & 166, 80, 64,     174, 79, 58,     181, 75, 55,     190, 71, 55,
     & 198, 68, 55,     208, 64, 55,     219, 59, 54,     234, 53, 53,
     & 251, 57, 50,     265, 68, 46,     277, 80, 40,     288, 98, 32,
C      imatp1=2: plasma
     &  62, 86, 97,      53, 84, 97,      46, 85, 98,      39, 82, 99,
     &  33, 77, 98,      25, 72, 96,      17, 67, 94,       8, 60, 90,
     & 356, 57, 86,     343, 62, 81,     331, 67, 77,     319, 73, 71,
     & 307, 80, 65,     294, 89, 63,     286, 97, 65,     279, 99, 65,
     & 272, 99, 64,     264, 97, 61,     256, 96, 58,     242, 94, 52,
C      imatp1=3: inferno
     &  61, 35, 99,      58, 53, 94,      48, 72, 96,      42, 87, 98,
     &  37, 96, 98,      31, 94, 97,      23, 87, 94,      14, 78, 89,
     &   4, 70, 84,     351, 69, 76,     340, 71, 69,     328, 73, 60,
     & 316, 75, 52,     301, 77, 44,     288, 83, 43,     277, 90, 41,
     & 267, 88, 35,     255, 80, 23,     245, 81, 11,     244, 96,  1,
C      imatp1=4: magma
     &  61, 24, 99,      43, 34, 99,      33, 43, 99,      25, 51, 99,
     &  19, 57, 99,      13, 62, 98,       5, 62, 95,     353, 64, 90,
     & 343, 68, 84,     333, 69, 75,     323, 70, 67,     311, 70, 58,
     & 298, 71, 50,     286, 78, 50,     277, 84, 49,     268, 86, 44,
     & 259, 80, 35,     249, 74, 22,     243, 78, 10,     244, 96,  1,
C      imatp1=5: cividis
     &  53, 78, 99,      51, 72, 95,      51, 63, 89,      50, 55, 83,
     &  49, 48, 77,      49, 40, 72,      48, 32, 66,      47, 24, 61,
     &  47, 16, 56,      47,  7, 51,      69,  0, 46,     223,  6, 44,
     & 224, 13, 43,     224, 23, 42,     223, 34, 42,     222, 48, 42,
     & 220, 64, 43,     215, 89, 44,     213,100, 39,     213,100, 30,
C      imatp1=6: greys
     &   0,  0,  0,       0,  0,  5,       0,  0, 11,       0,  0, 19,
     &   0,  0, 26,       0,  0, 33,       0,  0, 38,       0,  0, 44,
     &   0,  0, 49,       0,  0, 56,       0,  0, 61,       0,  0, 68,
     &   0,  0, 74,       0,  0, 79,       0,  0, 83,       0,  0, 88,
     &   0,  0, 91,       0,  0, 95,       0,  0, 97,       0,  0,100,
C      imatp1=7: purples
     & 270,100, 49,     268, 87, 51,     266, 77, 54,     264, 66, 58,
     & 261, 57, 61,     257, 48, 64,     251, 40, 68,     243, 33, 72,
     & 243, 29, 74,     244, 24, 77,     244, 21, 80,     241, 17, 83,
     & 238, 14, 86,     238, 10, 89,     239,  7, 91,     242,  5, 93,
     & 248,  4, 95,     255,  2, 96,     258,  1, 97,     270,  0, 99,
C      imatp1=8: blues
     & 215, 92, 41,     213, 93, 49,     211, 94, 57,     209, 91, 63,
     & 208, 85, 67,     207, 80, 71,     205, 73, 74,     203, 67, 77,
     & 203, 60, 79,     202, 53, 82,     202, 45, 84,     201, 36, 86,
     & 200, 29, 88,     203, 23, 90,     207, 18, 93,     209, 14, 94,
     & 208, 11, 96,     208,  8, 97,     209,  5, 98,     210,  3,100,
C      imatp1=9: greens
     & 143,100, 26,     144,100, 33,     144,100, 39,     143, 92, 45,
     & 141, 82, 50,     139, 73, 55,     137, 67, 60,     136, 62, 66,
     & 131, 53, 70,     125, 44, 74,     120, 37, 78,     117, 32, 82,
     & 114, 27, 85,     112, 22, 88,     110, 18, 90,     108, 14, 92,
     & 107, 10, 94,     105,  7, 96,     104,  5, 97,     102,  2, 98,
C      imatp1=10: oranges
     &  17, 96, 49,      17, 97, 56,      18, 97, 62,      19, 98, 70,
     &  19, 99, 78,      20, 98, 86,      21, 95, 89,      23, 92, 94,
     &  23, 86, 96,      24, 79, 98,      25, 72, 99,      26, 64, 99,
     &  27, 56, 99,      28, 47, 99,      29, 38, 99,      30, 30, 99,
     &  30, 23, 99,      30, 16, 99,      30, 12, 99,      30,  7,100,
C      imatp1=11: reds
     & 352,100, 40,     355, 95, 50,     356, 92, 60,     357, 90, 68,
     & 358, 88, 74,     359, 86, 81,       1, 83, 86,       4, 81, 93,
     &   6, 77, 95,       9, 72, 97,      11, 67, 98,      12, 60, 98,
     &  14, 53, 98,      15, 45, 98,      16, 38, 98,      17, 29, 99,
     &  18, 22, 99,      19, 15, 99,      19, 10, 99,      20,  5,100,
C      imatp1=12: ylorbr
     &  19, 94, 40,      19, 95, 48,      19, 96, 56,      20, 97, 65,
     &  21, 98, 73,      22, 98, 81,      23, 94, 86,      25, 91, 91,
     &  27, 88, 95,      30, 85, 98,      33, 80, 99,      36, 74, 99,
     &  40, 67, 99,      41, 56, 99,      44, 45, 99,      46, 37, 99,
     &  50, 30, 99,      53, 23,100,      55, 16,100,      60, 10,100,
C      imatp1=13: ylorrd
     & 342,100, 50,     345,100, 59,     347,100, 69,     350, 96, 77,
     & 355, 91, 83,       0, 87, 90,       5, 85, 94,       9, 83, 98,
     &  15, 80, 98,      21, 77, 99,      26, 75, 99,      30, 72, 99,
     &  34, 69, 99,      38, 61, 99,      42, 55, 99,      44, 48, 99,
     &  47, 41, 99,      49, 34,100,      53, 27,100,      60, 19,100,
C      imatp1=14: orrd
     &   0,100, 49,       0,100, 58,       0,100, 66,       1, 95, 73,
     &   3, 89, 79,       5, 83, 85,       7, 77, 89,      10, 70, 93,
     &  13, 67, 95,      17, 65, 97,      20, 61, 98,      23, 53, 99,
     &  27, 47, 99,      30, 42, 99,      33, 38, 99,      34, 32, 99,
     &  35, 25, 99,      35, 18, 99,      35, 13, 99,      34,  7,100,
C      imatp1=15: purd
     & 341,100, 40,     337,100, 48,     334,100, 56,     334, 97, 65,
     & 336, 93, 73,     337, 90, 81,     333, 86, 85,     329, 82, 90,
     & 327, 72, 89,     324, 60, 88,     321, 49, 85,     313, 37, 81,
     & 301, 25, 79,     297, 20, 81,     291, 16, 84,     284, 11, 88,
     & 276,  8, 91,     266,  5, 94,     269,  3, 96,     276,  2, 97,
C      imatp1=16: rdpu
     & 281,100, 41,     290, 99, 43,     298, 99, 45,     306, 99, 53,
     & 312, 99, 61,     317, 96, 70,     320, 86, 77,     324, 77, 85,
     & 328, 69, 90,     333, 61, 94,     337, 53, 97,     341, 44, 97,
     & 346, 35, 98,     352, 29, 98,       2, 24, 98,       5, 20, 98,
     &   5, 15, 99,       6, 11, 99,      10,  7, 99,      20,  4,100,
C      imatp1=17: bupu
     & 301,100, 30,     302, 93, 38,     302, 89, 46,     299, 79, 51,
     & 292, 67, 57,     284, 56, 62,     278, 48, 65,     269, 40, 69,
     & 257, 33, 72,     238, 27, 75,     225, 28, 79,     216, 28, 82,
     & 209, 26, 85,     209, 22, 87,     209, 17, 89,     208, 13, 92,
     & 206, 10, 94,     203,  7, 96,     199,  4, 97,     189,  2, 99,
C      imatp1=18: gnbu
     & 212, 93, 50,     208, 94, 57,     205, 95, 64,     203, 90, 69,
     & 202, 83, 72,     199, 75, 75,     197, 69, 78,     194, 63, 82,
     & 189, 54, 80,     180, 42, 78,     168, 36, 81,     152, 29, 84,
     & 133, 23, 86,     124, 18, 89,     112, 16, 91,     108, 14, 93,
     & 108, 11, 94,     105,  8, 95,      98,  6, 97,      85,  4, 98,
C      imatp1=19: pubu
     & 202, 97, 34,     202, 97, 42,     202, 97, 51,     202, 97, 58,
     & 202, 97, 64,     202, 94, 69,     201, 83, 72,     200, 73, 74,
     & 201, 61, 77,     203, 49, 79,     206, 39, 82,     209, 30, 84,
     & 214, 23, 86,     220, 16, 87,     232, 10, 89,     243,  7, 91,
     & 255,  5, 93,     274,  3, 95,     299,  2, 96,     330,  3,100,
C      imatp1=20: ylgnbu
     & 224, 90, 34,     227, 82, 44,     230, 76, 53,     226, 76, 60,
     & 218, 78, 63,     211, 80, 66,     204, 82, 70,     197, 84, 74,
     & 193, 78, 75,     189, 70, 76,     183, 60, 76,     175, 46, 77,
     & 164, 36, 80,     142, 25, 85,     106, 21, 90,      86, 24, 93,
     &  75, 27, 95,      68, 26, 97,      65, 20, 98,      60, 14,100,
C      imatp1=21: pubugn
     & 166, 98, 27,     167, 98, 33,     168, 99, 39,     173, 98, 44,
     & 179, 98, 47,     185, 94, 56,     192, 82, 64,     200, 72, 74,
     & 201, 63, 77,     201, 54, 79,     203, 44, 82,     207, 33, 84,
     & 214, 23, 86,     220, 16, 87,     232, 10, 89,     248,  7, 91,
     & 267,  6, 93,     287,  5, 94,     305,  3, 96,     330,  3,100,
C      imatp1=22: bugn
     & 143,100, 26,     144,100, 33,     144,100, 39,     143, 92, 45,
     & 141, 82, 50,     140, 73, 55,     144, 67, 61,     148, 63, 67,
     & 152, 56, 71,     157, 50, 74,     161, 43, 77,     163, 35, 81,
     & 165, 28, 85,     166, 21, 88,     168, 15, 91,     174, 11, 93,
     & 184,  8, 95,     191,  6, 97,     191,  4, 98,     189,  2, 99,
C      imatp1=23: ylgn
     & 155,100, 27,     153,100, 32,     152,100, 38,     148, 91, 43,
     & 144, 81, 48,     139, 72, 53,     137, 66, 59,     136, 62, 66,
     & 131, 53, 70,     124, 43, 75,     115, 38, 79,     104, 36, 83,
     &  95, 35, 87,      87, 33, 90,      79, 32, 93,      73, 30, 95,
     &  68, 27, 97,      64, 23, 99,      62, 16, 99,      60, 10,100,
C      imatp1=24: binary
     &   0,  0,  0,       0,  0,  5,       0,  0, 10,       0,  0, 15,
     &   0,  0, 20,       0,  0, 26,       0,  0, 31,       0,  0, 36,
     &   0,  0, 41,       0,  0, 47,       0,  0, 52,       0,  0, 58,
     &   0,  0, 63,       0,  0, 68,       0,  0, 73,       0,  0, 79,
     &   0,  0, 84,       0,  0, 89,       0,  0, 94,       0,  0,100,
C      imatp1=25: gist_yarg
     &   0,  0,  0,       0,  0,  5,       0,  0, 10,       0,  0, 15,
     &   0,  0, 20,       0,  0, 26,       0,  0, 31,       0,  0, 36,
     &   0,  0, 41,       0,  0, 47,       0,  0, 52,       0,  0, 58,
     &   0,  0, 63,       0,  0, 68,       0,  0, 73,       0,  0, 79,
     &   0,  0, 84,       0,  0, 89,       0,  0, 94,       0,  0,100,
C      imatp1=26: gist_gray
     &   0,  0,100,       0,  0, 94,       0,  0, 89,       0,  0, 84,
     &   0,  0, 79,       0,  0, 73,       0,  0, 68,       0,  0, 63,
     &   0,  0, 58,       0,  0, 52,       0,  0, 47,       0,  0, 41,
     &   0,  0, 36,       0,  0, 31,       0,  0, 26,       0,  0, 20,
     &   0,  0, 15,       0,  0, 10,       0,  0,  5,       0,  0,  0,
C      imatp1=27: gray
     &   0,  0,100,       0,  0, 94,       0,  0, 89,       0,  0, 84,
     &   0,  0, 79,       0,  0, 73,       0,  0, 68,       0,  0, 63,
     &   0,  0, 58,       0,  0, 52,       0,  0, 47,       0,  0, 41,
     &   0,  0, 36,       0,  0, 31,       0,  0, 26,       0,  0, 20,
     &   0,  0, 15,       0,  0, 10,       0,  0,  5,       0,  0,  0,
C      imatp1=28: bone
     &   0,  0,100,     179,  2, 95,     179,  5, 91,     179,  8, 86,
     & 179, 12, 81,     181, 16, 77,     189, 17, 72,     198, 18, 67,
     & 206, 19, 63,     214, 21, 58,     222, 23, 54,     231, 25, 49,
     & 239, 27, 44,     240, 28, 38,     240, 28, 31,     240, 28, 25,
     & 240, 28, 19,     240, 28, 12,     240, 28,  6,       0,  0,  0,
C      imatp1=29: pink
     &   0,  0,100,      60,  3, 98,      60,  7, 96,      60, 12, 94,
     &  60, 17, 92,      58, 22, 90,      51, 23, 88,      43, 25, 86,
     &  36, 26, 84,      27, 28, 82,      19, 30, 80,      10, 32, 78,
     &   0, 34, 76,       0, 35, 70,       0, 35, 64,       0, 35, 57,
     &   0, 36, 50,       0, 36, 41,       0, 39, 30,       0,100, 11,
C      imatp1=30: spring
     &  60,100,100,      56, 94,100,      53, 89,100,      48, 84,100,
     &  44, 79,100,      38, 73,100,      32, 68,100,      24, 63,100,
     &  16, 58,100,       5, 52,100,     354, 52,100,     343, 58,100,
     & 335, 63,100,     327, 68,100,     321, 73,100,     315, 79,100,
     & 311, 84,100,     306, 89,100,     303, 94,100,     300,100,100,
C      imatp1=31: summer
     &  60, 60,100,      62, 58, 97,      65, 57, 94,      69, 56, 92,
     &  72, 55, 89,      76, 53, 86,      81, 52, 84,      86, 50, 81,
     &  92, 49, 79,      99, 47, 76,     106, 45, 73,     116, 43, 70,
     & 125, 46, 68,     135, 52, 65,     142, 58, 63,     149, 65, 60,
     & 154, 72, 57,     159, 81, 55,     164, 90, 52,     167,100, 50,
C      imatp1=32: autumn
     &  60,100,100,      56,100,100,      53,100,100,      50,100,100,
     &  47,100,100,      44,100,100,      41,100,100,      37,100,100,
     &  34,100,100,      31,100,100,      28,100,100,      25,100,100,
     &  22,100,100,      18,100,100,      15,100,100,      12,100,100,
     &   9,100,100,       6,100,100,       3,100,100,       0,100,100,
C      imatp1=33: winter
     & 150,100,100,     153,100, 94,     156,100, 89,     161,100, 84,
     & 165,100, 79,     171,100, 73,     177,100, 68,     184,100, 68,
     & 190,100, 70,     197,100, 73,     202,100, 76,     208,100, 79,
     & 212,100, 81,     217,100, 84,     221,100, 86,     226,100, 89,
     & 229,100, 92,     233,100, 94,     236,100, 97,     240,100,100,
C      imatp1=34: cool
     & 300,100,100,     296, 94,100,     293, 89,100,     288, 84,100,
     & 284, 79,100,     278, 73,100,     272, 68,100,     264, 63,100,
     & 256, 58,100,     245, 52,100,     234, 52,100,     223, 58,100,
     & 215, 63,100,     207, 68,100,     201, 73,100,     195, 79,100,
     & 191, 84,100,     186, 89,100,     183, 94,100,     180,100,100,
C      imatp1=35: wistia
     &  30,100, 98,      31,100, 99,      33,100, 99,      34,100, 99,
     &  36,100, 99,      37,100,100,      39,100,100,      40,100,100,
     &  42,100,100,      43,100,100,      45, 98,100,      47, 96,100,
     &  49, 94,100,      51, 92,100,      53, 90,100,      55, 83, 98,
     &  58, 74, 96,      62, 66, 96,      66, 59, 98,      72, 52,100,
C      imatp1=36: hot
     &   0,  0,100,      60, 20,100,      60, 40,100,      60, 61,100,
     &  60, 81,100,      58,100,100,      50,100,100,      41,100,100,
     &  33,100,100,      25,100,100,      17,100,100,       8,100,100,
     &   0,100,100,       0,100, 86,       0,100, 73,       0,100, 58,
     &   0,100, 45,       0,100, 30,       0,100, 17,       0,100,  4,
C      imatp1=37: afmhot
     &   0,  0,100,      60, 10,100,      60, 20,100,      60, 31,100,
     &  60, 41,100,      57, 52,100,      47, 62,100,      40, 73,100,
     &  35, 83,100,      31, 94,100,      28,100, 94,      24,100, 83,
     &  19,100, 73,      12,100, 62,       2,100, 52,       0,100, 41,
     &   0,100, 31,       0,100, 20,       0,100, 10,       0,  0,  0,
C      imatp1=38: gist_heat
     &   0,  0,100,      30, 20,100,      30, 40,100,      30, 62,100,
     &  30, 83,100,      28,100,100,      22,100,100,      16,100, 94,
     &  11,100, 87,       3,100, 78,       0,100, 71,       0,100, 62,
     &   0,100, 55,       0,100, 47,       0,100, 39,       0,100, 31,
     &   0,100, 23,       0,100, 15,       0,100,  7,       0,  0,  0,
C      imatp1=39: copper
     &  33, 50,100,      30, 52,100,      27, 55,100,      24, 58,100,
     &  23, 59, 97,      23, 59, 91,      23, 59, 84,      23, 59, 77,
     &  23, 59, 71,      23, 59, 64,      23, 59, 58,      23, 59, 51,
     &  23, 59, 45,      23, 59, 38,      23, 59, 32,      23, 59, 25,
     &  23, 59, 19,      23, 59, 12,      23, 59,  6,       0,  0,  0/

      integer, parameter :: n_matp2 = 32
      real matplotlib2(3,n_matp2,40) ! (HSB,i,imatp)
      data matplotlib2 /
C      imatp2=1: piyg
     & 108, 75, 39,     103, 75, 44,     100, 76, 50,      97, 77, 56,
     &  94, 73, 62,      92, 69, 67,      90, 66, 72,      89, 58, 77,
     &  88, 49, 82,      87, 42, 86,      86, 34, 89,      86, 26, 92,
     &  84, 17, 95,      84, 12, 96,      84,  7, 96,      84,  2, 96,
     & 328,  2, 97,     328,  5, 98,     328,  9, 98,     327, 12, 98,
     & 325, 17, 97,     324, 21, 95,     324, 26, 93,     325, 33, 91,
     & 327, 40, 88,     327, 49, 86,     326, 60, 83,     325, 73, 80,
     & 325, 86, 75,     325, 90, 69,     325, 94, 62,     325, 99, 55,
C      imatp2=2: prgn
     & 143,100, 26,     141, 89, 33,     139, 83, 39,     138, 78, 45,
     & 134, 66, 53,     130, 57, 59,     126, 49, 66,     123, 41, 72,
     & 119, 33, 78,     115, 28, 84,     113, 23, 87,     111, 18, 90,
     & 108, 13, 93,     107,  9, 94,     107,  5, 95,     107,  2, 96,
     & 296,  1, 95,     297,  4, 93,     297,  6, 92,     293,  9, 89,
     & 286, 13, 86,     283, 17, 83,     281, 21, 79,     281, 25, 75,
     & 281, 30, 70,     282, 36, 65,     286, 45, 60,     289, 56, 55,
     & 291, 69, 50,     291, 76, 43,     291, 85, 36,     291,100, 29,
C      imatp2=3: brbg
     & 167,100, 23,     171, 99, 28,     173, 99, 33,     174, 99, 39,
     & 175, 86, 45,     175, 75, 51,     175, 66, 57,     174, 56, 64,
     & 172, 46, 71,     171, 39, 78,     170, 32, 82,     170, 24, 86,
     & 171, 17, 90,     171, 11, 92,     171,  7, 93,     171,  2, 95,
     &  43,  3, 96,      43, 10, 96,      43, 16, 96,      43, 23, 95,
     &  42, 31, 92,      42, 38, 89,      41, 46, 86,      38, 56, 82,
     &  35, 67, 77,      34, 77, 73,      33, 81, 67,      33, 86, 60,
     &  32, 92, 53,      32, 93, 46,      32, 93, 39,      32, 94, 32,
C      imatp2=4: puor
     & 276,100, 29,     272, 87, 36,     270, 78, 44,     268, 72, 51,
     & 265, 58, 57,     261, 46, 61,     255, 35, 66,     253, 29, 70,
     & 252, 23, 76,     251, 19, 80,     248, 15, 84,     244, 12, 87,
     & 237,  8, 91,     233,  6, 93,     233,  3, 94,     233,  1, 96,
     &  35,  5, 97,      35, 14, 98,      34, 22, 99,      34, 32, 99,
     &  33, 43, 99,      33, 53, 99,      33, 63, 98,      32, 72, 94,
     &  32, 83, 90,      32, 91, 86,      30, 92, 80,      29, 94, 75,
     &  28, 96, 68,      27, 95, 62,      26, 94, 56,      25, 93, 49,
C      imatp2=5: rdgy
     &   0,  0, 10,       0,  0, 16,       0,  0, 22,       0,  0, 29,
     &   0,  0, 36,       0,  0, 44,       0,  0, 51,       0,  0, 57,
     &   0,  0, 64,       0,  0, 70,       0,  0, 76,       0,  0, 80,
     &   0,  0, 86,       0,  0, 90,       0,  0, 94,       0,  0, 97,
     &  22,  3, 99,      22, 10, 99,      22, 17, 99,      21, 24, 98,
     &  19, 33, 97,      18, 41, 96,      17, 48, 94,      13, 53, 90,
     &  10, 59, 86,       7, 64, 82,       1, 68, 78,     356, 76, 73,
     & 352, 86, 68,     349, 89, 58,     346, 94, 49,     341,100, 40,
C      imatp2=6: rdbu
     & 211, 94, 38,     211, 88, 47,     210, 84, 56,     210, 81, 65,
     & 208, 75, 70,     205, 71, 72,     203, 66, 75,     202, 57, 78,
     & 201, 46, 82,     200, 36, 86,     199, 29, 88,     200, 22, 90,
     & 200, 15, 93,     201, 10, 94,     201,  6, 95,     201,  2, 96,
     &  22,  3, 97,      22, 10, 98,      22, 17, 98,      21, 24, 98,
     &  19, 33, 97,      18, 41, 96,      17, 48, 94,      13, 53, 90,
     &  10, 59, 86,       7, 64, 82,       1, 68, 78,     356, 76, 73,
     & 352, 86, 68,     349, 89, 58,     346, 94, 49,     341,100, 40,
C      imatp2=7: rdylbu
     & 237, 67, 58,     229, 65, 62,     222, 63, 66,     215, 61, 69,
     & 211, 56, 73,     207, 50, 77,     204, 45, 81,     201, 39, 84,
     & 199, 33, 87,     196, 28, 90,     195, 22, 92,     194, 17, 94,
     & 193, 11, 96,     146,  6, 96,      78, 12, 97,      63, 20, 99,
     &  55, 28, 99,      49, 34, 99,      45, 39, 99,      41, 45, 99,
     &  35, 51, 99,      31, 57, 99,      27, 62, 98,      22, 66, 97,
     &  17, 69, 96,      13, 73, 94,       9, 75, 91,       6, 78, 87,
     &   2, 81, 83,     356, 84, 77,     351, 91, 70,     346,100, 64,
C      imatp2=8: rdylgn
     & 151,100, 40,     149, 93, 46,     147, 87, 52,     146, 83, 58,
     & 139, 70, 63,     131, 58, 68,     120, 48, 72,     109, 48, 76,
     &  98, 49, 80,      90, 50, 84,      84, 48, 86,      80, 45, 89,
     &  75, 42, 92,      71, 38, 94,      67, 33, 96,      63, 27, 98,
     &  55, 28, 99,      50, 34, 99,      46, 41, 99,      42, 47, 99,
     &  36, 52, 99,      32, 58, 99,      27, 62, 98,      22, 66, 97,
     &  17, 69, 96,      13, 73, 94,       9, 75, 91,       6, 78, 87,
     &   2, 81, 83,     356, 84, 77,     351, 91, 70,     346,100, 64,
C      imatp2=9: spectral
     & 250, 51, 63,     228, 52, 66,     214, 62, 70,     204, 71, 73,
     & 194, 64, 71,     182, 53, 68,     165, 48, 74,     153, 40, 78,
     & 139, 32, 82,     119, 24, 85,      99, 28, 88,      84, 32, 91,
     &  72, 36, 94,      68, 35, 96,      65, 31, 98,      62, 27, 99,
     &  55, 28, 99,      50, 34, 99,      46, 41, 99,      42, 47, 99,
     &  36, 52, 99,      32, 58, 99,      27, 62, 98,      22, 66, 97,
     &  17, 69, 96,      12, 71, 94,       7, 69, 90,       0, 66, 87,
     & 352, 72, 82,     346, 79, 75,     340, 88, 68,     335, 99, 61,
C      imatp2=10: coolwarm
     & 348, 97, 70,     357, 79, 75,       2, 72, 79,       5, 68, 83,
     &   8, 64, 87,      10, 61, 90,      12, 57, 92,      14, 53, 94,
     &  15, 48, 96,      16, 43, 96,      17, 38, 96,      18, 33, 96,
     &  19, 26, 95,      20, 20, 93,      21, 12, 91,      21,  4, 88,
     & 216,  4, 88,     216, 11, 91,     217, 17, 94,     218, 23, 96,
     & 218, 29, 98,     219, 34, 99,     220, 39, 99,     221, 43, 99,
     & 222, 47, 98,     223, 51, 97,     224, 54, 94,     225, 57, 92,
     & 227, 61, 88,     228, 64, 84,     230, 66, 80,     232, 69, 75,
C      imatp2=11: bwr
     &   0,100,100,       0, 93,100,       0, 87,100,       0, 81,100,
     &   0, 74,100,       0, 67,100,       0, 61,100,       0, 55,100,
     &   0, 48,100,       0, 41,100,       0, 35,100,       0, 29,100,
     &   0, 22,100,       0, 16,100,       0,  9,100,       0,  3,100,
     & 240,  3,100,     240,  9,100,     240, 16,100,     240, 22,100,
     & 240, 29,100,     240, 35,100,     240, 41,100,     240, 48,100,
     & 240, 55,100,     240, 61,100,     240, 67,100,     240, 74,100,
     & 240, 81,100,     240, 87,100,     240, 93,100,     240,100,100,
C      imatp2=12: seismic
     &   0,100, 50,       0,100, 56,       0,100, 62,       0,100, 68,
     &   0,100, 75,       0,100, 82,       0,100, 88,       0,100, 94,
     &   0, 96,100,       0, 83,100,       0, 71,100,       0, 58,100,
     &   0, 44,100,       0, 32,100,       0, 19,100,       0,  7,100,
     & 240,  7,100,     240, 19,100,     240, 32,100,     240, 44,100,
     & 240, 58,100,     240, 71,100,     240, 83,100,     240, 96,100,
     & 240,100, 92,     240,100, 83,     240,100, 75,     240,100, 66,
     & 240,100, 56,     240,100, 47,     240,100, 38,     240,100, 30,
C      imatp2=13: twilight
     & 299,  4, 88,     353,  5, 87,      20, 12, 84,      24, 23, 81,
     &  23, 34, 79,      19, 43, 77,      15, 50, 75,       8, 53, 72,
     &   0, 53, 68,     350, 60, 63,     340, 67, 57,     331, 72, 50,
     & 322, 75, 42,     314, 75, 34,     305, 73, 25,     294, 68, 21,
     & 286, 72, 24,     282, 78, 31,     277, 78, 42,     270, 73, 52,
     & 262, 64, 61,     252, 55, 66,     240, 47, 70,     228, 47, 72,
     & 220, 46, 74,     212, 41, 75,     207, 35, 76,     201, 26, 77,
     & 198, 17, 79,     200,  8, 83,     234,  2, 86,     296,  4, 88,
C      imatp2=14: twilight_shifted
     & 287, 62, 21,     300, 70, 22,     309, 74, 29,     318, 75, 38,
     & 326, 74, 46,     336, 69, 54,     345, 64, 60,     355, 56, 66,
     &   4, 53, 70,      12, 52, 74,      17, 47, 76,      21, 39, 78,
     &  24, 28, 80,      23, 18, 83,      12,  8, 86,     328,  4, 88,
     & 273,  3, 87,     209,  5, 84,     197, 12, 81,     199, 21, 78,
     & 204, 31, 76,     210, 38, 75,     216, 44, 74,     224, 47, 73,
     & 234, 47, 71,     246, 51, 68,     257, 60, 64,     266, 69, 57,
     & 274, 76, 48,     280, 79, 36,     284, 76, 27,     287, 64, 21,
C      imatp2=15: hsv
     & 354,100,100,     343,100,100,     332,100,100,     321,100,100,
     & 308,100,100,     297,100,100,     286,100,100,     275,100,100,
     & 262,100,100,     251,100,100,     240, 98,100,     229,100,100,
     & 216,100,100,     205,100,100,     194,100,100,     183,100,100,
     & 170,100,100,     159,100,100,     148,100,100,     137,100,100,
     & 125,100,100,     113,100,100,     102,100,100,      91,100,100,
     &  79,100,100,      68,100,100,      56,100, 99,      45,100,100,
     &  33,100,100,      22,100,100,      11,100,100,       0,100,100,
C      imatp2=16: pastel1
     &   0,  0, 94,       0,  0, 94,       0,  0, 94,       0,  0, 94,
     & 329, 13, 99,     329, 13, 99,     329, 13, 99,      40, 17, 89,
     &  40, 17, 89,      40, 17, 89,      40, 17, 89,      60, 19,100,
     &  60, 19,100,      60, 19,100,      34, 34, 99,      34, 34, 99,
     &  34, 34, 99,      34, 34, 99,     285, 10, 89,     285, 10, 89,
     & 285, 10, 89,     108, 16, 92,     108, 16, 92,     108, 16, 92,
     & 108, 16, 92,     207, 21, 89,     207, 21, 89,     207, 21, 89,
     &   4, 30, 98,       4, 30, 98,       4, 30, 98,       4, 30, 98,
C      imatp2=17: pastel2
     &   0,  0, 80,       0,  0, 80,       0,  0, 80,       0,  0, 80,
     &  35, 15, 94,      35, 15, 94,      35, 15, 94,      35, 15, 94,
     &  50, 31,100,      50, 31,100,      50, 31,100,      50, 31,100,
     &  80, 17, 96,      80, 17, 96,      80, 17, 96,      80, 17, 96,
     & 322, 17, 95,     322, 17, 95,     322, 17, 95,     322, 17, 95,
     & 219, 12, 90,     219, 12, 90,     219, 12, 90,     219, 12, 90,
     &  24, 32, 99,      24, 32, 99,      24, 32, 99,      24, 32, 99,
     & 153, 20, 88,     153, 20, 88,     153, 20, 88,     153, 20, 88,
C      imatp2=18: paired
     &  21, 77, 69,      21, 77, 69,      21, 77, 69,      60, 40,100,
     &  60, 40,100,      60, 40,100,     269, 60, 60,     269, 60, 60,
     & 279, 16, 83,     279, 16, 83,     279, 16, 83,      29,100,100,
     &  29,100,100,      33, 56, 99,      33, 56, 99,      33, 56, 99,
     & 359, 88, 89,     359, 88, 89,     359, 88, 89,       0, 39, 98,
     &   0, 39, 98,     116, 72, 62,     116, 72, 62,     116, 72, 62,
     &  91, 38, 87,      91, 38, 87,     204, 82, 70,     204, 82, 70,
     & 204, 82, 70,     200, 26, 89,     200, 26, 89,     200, 26, 89,
C      imatp2=19: accent
     &   0,  0, 40,       0,  0, 40,       0,  0, 40,       0,  0, 40,
     &  24, 87, 74,      24, 87, 74,      24, 87, 74,      24, 87, 74,
     & 328, 99, 94,     328, 99, 94,     328, 99, 94,     328, 99, 94,
     & 214, 68, 69,     214, 68, 69,     214, 68, 69,     214, 68, 69,
     &  60, 40,100,      60, 40,100,      60, 40,100,      60, 40,100,
     &  29, 47, 99,      29, 47, 99,      29, 47, 99,      29, 47, 99,
     & 265, 17, 83,     265, 17, 83,     265, 17, 83,     265, 17, 83,
     & 120, 36, 78,     120, 36, 78,     120, 36, 78,     120, 36, 78,
C      imatp2=20: dark2
     &   0,  0, 40,       0,  0, 40,       0,  0, 40,       0,  0, 40,
     &  38, 82, 65,      38, 82, 65,      38, 82, 65,      38, 82, 65,
     &  44, 99, 90,      44, 99, 90,      44, 99, 90,      44, 99, 90,
     &  88, 81, 65,      88, 81, 65,      88, 81, 65,      88, 81, 65,
     & 329, 82, 90,     329, 82, 90,     329, 82, 90,     329, 82, 90,
     & 244, 37, 70,     244, 37, 70,     244, 37, 70,     244, 37, 70,
     &  25, 99, 85,      25, 99, 85,      25, 99, 85,      25, 99, 85,
     & 162, 82, 61,     162, 82, 61,     162, 82, 61,     162, 82, 61,
C      imatp2=21: set1
     &   0,  0, 60,       0,  0, 60,       0,  0, 60,       0,  0, 60,
     & 328, 47, 96,     328, 47, 96,     328, 47, 96,      21, 75, 65,
     &  21, 75, 65,      21, 75, 65,      21, 75, 65,      60, 80,100,
     &  60, 80,100,      60, 80,100,      29,100,100,      29,100,100,
     &  29,100,100,      29,100,100,     292, 52, 63,     292, 52, 63,
     & 292, 52, 63,     118, 57, 68,     118, 57, 68,     118, 57, 68,
     & 118, 57, 68,     206, 70, 72,     206, 70, 72,     206, 70, 72,
     & 359, 88, 89,     359, 88, 89,     359, 88, 89,     359, 88, 89,
C      imatp2=22: set2
     &   0,  0, 70,       0,  0, 70,       0,  0, 70,       0,  0, 70,
     &  35, 35, 89,      35, 35, 89,      35, 35, 89,      35, 35, 89,
     &  49, 81,100,      49, 81,100,      49, 81,100,      49, 81,100,
     &  82, 61, 84,      82, 61, 84,      82, 61, 84,      82, 61, 84,
     & 323, 40, 90,     323, 40, 90,     323, 40, 90,     323, 40, 90,
     & 221, 30, 79,     221, 30, 79,     221, 30, 79,     221, 30, 79,
     &  16, 61, 98,      16, 61, 98,      16, 61, 98,      16, 61, 98,
     & 161, 47, 76,     161, 47, 76,     161, 47, 76,     161, 47, 76,
C      imatp2=23: set3
     &  52, 56,100,      52, 56,100,      52, 56,100,     108, 16, 92,
     & 108, 16, 92,     108, 16, 92,     299, 32, 74,     299, 32, 74,
     &   0,  0, 85,       0,  0, 85,       0,  0, 85,     329, 18, 98,
     & 329, 18, 98,      82, 52, 87,      82, 52, 87,      82, 52, 87,
     &  31, 61, 99,      31, 61, 99,      31, 61, 99,     204, 39, 82,
     & 204, 39, 82,       6, 54, 98,       6, 54, 98,       6, 54, 98,
     & 247, 14, 85,     247, 14, 85,      60, 29,100,      60, 29,100,
     &  60, 29,100,     169, 33, 82,     169, 33, 82,     169, 33, 82,
C      imatp2=24: tab10
     & 185, 88, 81,     185, 88, 81,     185, 88, 81,     185, 88, 81,
     &  60, 82, 74,      60, 82, 74,      60, 82, 74,       0,  0, 49,
     &   0,  0, 49,       0,  0, 49,     318, 47, 89,     318, 47, 89,
     & 318, 47, 89,      10, 46, 54,      10, 46, 54,      10, 46, 54,
     & 271, 45, 74,     271, 45, 74,     271, 45, 74,     359, 81, 83,
     & 359, 81, 83,     359, 81, 83,     120, 72, 62,     120, 72, 62,
     & 120, 72, 62,      28, 94,100,      28, 94,100,      28, 94,100,
     & 204, 82, 70,     204, 82, 70,     204, 82, 70,     204, 82, 70,
C      imatp2=25: tab20
     & 189, 31, 89,     189, 31, 89,     185, 88, 81,     185, 88, 81,
     &  60, 35, 85,      60, 82, 74,      60, 82, 74,       0,  0, 78,
     &   0,  0, 49,       0,  0, 49,     334, 26, 96,     318, 47, 89,
     & 318, 47, 89,       9, 24, 76,      10, 46, 54,      10, 46, 54,
     & 274, 17, 83,     274, 17, 83,     271, 45, 74,       1, 41,100,
     &   1, 41,100,     359, 81, 83,     110, 38, 87,     110, 38, 87,
     & 120, 72, 62,      29, 52,100,      29, 52,100,      28, 94,100,
     & 214, 24, 90,     214, 24, 90,     204, 82, 70,     204, 82, 70,
C      imatp2=26: ocean
     &   0,  0,100,     195,  6, 96,     195, 13, 93,     195, 20, 90,
     & 195, 29, 87,     195, 38, 83,     195, 47, 80,     195, 57, 77,
     & 195, 69, 74,     195, 81, 70,     195, 94, 67,     196,100, 64,
     & 199,100, 61,     201,100, 58,     204,100, 54,     207,100, 51,
     & 212,100, 48,     216,100, 45,     221,100, 41,     227,100, 38,
     & 235,100, 35,     236,100, 32,     226,100, 29,     214,100, 25,
     & 195,100, 22,     174,100, 21,     157,100, 25,     145,100, 30,
     & 135,100, 35,     129,100, 40,     124,100, 45,     120,100, 50,
C      imatp2=27: gist_earth
     &   0,  0, 99,       2,  7, 95,      11, 13, 92,      18, 21, 88,
     &  24, 28, 84,      30, 34, 81,      36, 41, 77,      43, 47, 75,
     &  50, 48, 73,      57, 48, 72,      64, 48, 70,      70, 49, 69,
     &  78, 49, 67,      85, 50, 66,      92, 50, 65,     102, 52, 63,
     & 112, 53, 61,     124, 55, 59,     134, 56, 58,     144, 57, 56,
     & 155, 59, 54,     166, 60, 52,     177, 62, 50,     185, 65, 49,
     & 194, 69, 48,     201, 73, 48,     209, 77, 47,     217, 81, 47,
     & 226, 86, 46,     234, 90, 46,     241, 98, 45,       0,  0,  0,
C      imatp2=28: terrain
     &   0,  0,100,      10,  2, 93,      10,  4, 87,      10,  7, 81,
     &  10, 11, 74,      10, 16, 67,      10, 21, 61,      10, 27, 55,
     &  14, 34, 51,      25, 35, 58,      34, 36, 64,      41, 37, 70,
     &  47, 38, 77,      51, 38, 83,      55, 39, 90,      58, 39, 96,
     &  68, 40, 98,      83, 41, 96,      98, 42, 93,     113, 43, 91,
     & 129, 53, 88,     137, 66, 85,     144, 80, 83,     148, 95, 80,
     & 164,100, 74,     185,100, 74,     200,100, 93,     206, 97, 94,
     & 210, 91, 85,     216, 84, 76,     225, 76, 68,     240, 66, 60,
C      imatp2=29: gist_stern
     &   0,  0,100,      60,  8, 96,      60, 18, 93,      60, 28, 90,
     &  60, 41, 87,      60, 53, 83,      60, 65, 80,      60, 79, 77,
     &  60, 96, 74,      60, 84, 70,      60, 64, 67,      60, 42, 64,
     &  60, 14, 61,     240, 11, 65,     240, 30, 79,     240, 44, 92,
     & 240, 50, 96,     240, 50, 90,     240, 50, 83,     240, 50, 77,
     & 240, 50, 70,     240, 50, 64,     240, 50, 58,     240, 50, 51,
     & 226, 64, 44,     278, 50, 38,     328, 65, 47,     344, 79, 62,
     & 352, 88, 80,     355, 93, 95,     356, 94, 57,       0,  0,  0,
C      imatp2=30: gnuplot
     &  60,100,100,      55,100, 98,      51,100, 96,      46,100, 95,
     &  42,100, 93,      38,100, 91,      35,100, 89,      31,100, 88,
     &  28,100, 86,      25,100, 84,      22,100, 82,      20,100, 80,
     &  17,100, 78,      15,100, 76,      13,100, 74,      11,100, 71,
     &   0, 84, 69,     338, 86, 67,     317, 88, 64,     297, 90, 64,
     & 283, 94, 79,     276, 96, 90,     272, 97, 96,     270, 98, 99,
     & 268, 98, 98,     267, 99, 93,     268, 99, 84,     269, 99, 72,
     & 272, 99, 55,     279, 99, 38,     294, 99, 19,       0,  0,  0,
C      imatp2=31: gnuplot2
     &   0,  0,100,      60, 39,100,      60, 78,100,      58, 97,100,
     &  53, 90,100,      48, 83,100,      42, 77,100,      35, 71,100,
     &  26, 64,100,      16, 57,100,       3, 51,100,     349, 54,100,
     & 337, 61,100,     328, 67,100,     317, 72, 93,     302, 76, 83,
     & 288, 85, 87,     278, 93, 93,     271,100,100,     265,100,100,
     & 259,100,100,     253,100,100,     247,100,100,     241,100,100,
     & 240,100, 89,     240,100, 76,     240,100, 64,     240,100, 51,
     & 240,100, 37,     240,100, 25,     240,100, 12,       0,  0,  0,
C      imatp2=32: cmrmap
     &   0,  0,100,      60, 10, 97,      60, 21, 94,      60, 32, 92,
     &  59, 46, 90,      54, 57, 90,      52, 68, 90,      50, 79, 90,
     &  47, 89, 90,      43, 92, 90,      39, 95, 90,      35, 98, 90,
     &  30, 98, 91,      24, 94, 93,      17, 90, 96,      10, 86, 98,
     &   3, 78, 94,     354, 72, 84,     342, 70, 74,     324, 68, 64,
     & 300, 65, 54,     281, 70, 60,     267, 75, 66,     257, 79, 73,
     & 252, 78, 69,     249, 76, 63,     246, 73, 57,     240, 70, 50,
     & 240, 70, 37,     240, 70, 25,     240, 70, 12,       0,  0,  0,
C      imatp2=33: cubehelix
     &   0,  0,100,     148,  4, 98,     163,  8, 96,     178, 11, 93,
     & 196, 16, 94,     210, 19, 95,     225, 20, 95,     244, 20, 94,
     & 265, 25, 91,     283, 27, 87,     304, 29, 82,     322, 35, 83,
     & 339, 38, 82,     354, 38, 79,      11, 44, 74,      25, 51, 67,
     &  41, 56, 57,      59, 57, 48,      78, 61, 47,      96, 60, 47,
     & 119, 54, 45,     138, 64, 42,     153, 71, 38,     168, 74, 33,
     & 187, 73, 30,     201, 70, 30,     215, 65, 27,     231, 56, 23,
     & 254, 57, 17,     272, 62, 11,     291, 64,  5,       0,  0,  0,
C      imatp2=34: brg
     & 120,100,100,     115,100, 93,     111,100, 87,     106,100, 81,
     &  99,100, 74,      91,100, 67,      82,100, 61,      71,100, 55,
     &  55,100, 51,      43,100, 58,      33,100, 64,      25,100, 70,
     &  17,100, 77,      11,100, 83,       6,100, 90,       2,100, 96,
     & 357,100, 96,     353,100, 90,     348,100, 83,     342,100, 77,
     & 335,100, 70,     326,100, 64,     316,100, 58,     304,100, 51,
     & 288,100, 55,     277,100, 61,     268,100, 67,     260,100, 74,
     & 253,100, 81,     248,100, 87,     244,100, 93,     240,100,100,
C      imatp2=35: gist_rainbow
     & 315,100,100,     304,100,100,     294,100,100,     284,100,100,
     & 272,100,100,     262,100,100,     252,100,100,     242,100,100,
     & 230,100,100,     220,100,100,     210,100,100,     199,100,100,
     & 188,100,100,     178,100,100,     168,100,100,     157,100,100,
     & 146,100,100,     136,100,100,     126,100,100,     116,100,100,
     & 104,100,100,      94,100,100,      84,100,100,      74,100,100,
     &  62,100,100,      52,100,100,      42,100,100,      32,100,100,
     &  20,100,100,      10,100,100,       0,100,100,     350,100,100,
C      imatp2=36: rainbow
     &   0, 99,100,       3, 95,100,       6, 90,100,      10, 85,100,
     &  14, 79,100,      18, 75,100,      23, 70,100,      27, 65,100,
     &  33, 59, 98,      43, 52, 91,      58, 43, 85,      76, 41, 89,
     &  95, 39, 93,     111, 36, 96,     128, 39, 98,     139, 46, 99,
     & 149, 53, 99,     156, 59, 98,     163, 64, 96,     169, 70, 93,
     & 176, 77, 89,     182, 83, 87,     187, 91, 89,     192, 98, 91,
     & 199, 94, 93,     207, 87, 95,     216, 81, 96,     227, 75, 97,
     & 241, 70, 98,     253, 80, 99,     262, 90, 99,     270,100,100,
C      imatp2=37: jet
     &   0,100, 50,       0,100, 64,       0,100, 78,       0,100, 92,
     &   8,100,100,      15,100,100,      22,100,100,      29,100,100,
     &  37,100,100,      44,100,100,      51,100,100,      60, 99, 97,
     &  70, 87,100,      79, 77,100,      91, 67,100,     108, 57,100,
     & 131, 57,100,     148, 67,100,     160, 77,100,     169, 87,100,
     & 182, 99, 95,     192,100,100,     200,100,100,     207,100,100,
     & 216,100,100,     223,100,100,     231,100,100,     238,100,100,
     & 240,100, 92,     240,100, 78,     240,100, 64,     240,100, 50,
C      imatp2=38: turbo
     &   0, 97, 47,       4, 99, 58,       7, 99, 67,      10, 98, 75,
     &  13, 97, 83,      15, 95, 88,      18, 93, 93,      22, 89, 96,
     &  27, 85, 98,      32, 81, 99,      37, 78, 98,      44, 76, 96,
     &  54, 75, 90,      65, 76, 89,      74, 78, 94,      83, 77, 97,
     &  93, 73, 99,     108, 64, 99,     129, 64, 98,     146, 75, 96,
     & 160, 85, 92,     169, 89, 87,     180, 87, 82,     191, 83, 90,
     & 203, 78, 97,     212, 74, 99,     219, 71, 97,     225, 69, 90,
     & 231, 65, 78,     239, 60, 63,     252, 63, 45,     284, 69, 23,
C      imatp2=39: nipy_spectral
     &   0,  0, 80,       0, 62, 80,       0,100, 81,       0,100, 85,
     &   0,100, 94,       7,100,100,      30,100,100,      41,100,100,
     &  49,100, 98,      57,100, 94,      67,100, 96,      78,100,100,
     & 109,100,100,     120,100, 94,     120,100, 86,     120,100, 78,
     & 120,100, 68,     120,100, 60,     150,100, 64,     170,100, 66,
     & 179,100, 66,     191,100, 78,     200,100, 86,     206,100, 86,
     & 224,100, 86,     240,100, 83,     240,100, 70,     260,100, 63,
     & 293,100, 59,     292,100, 55,     292,100, 33,       0,  0,  0,
C      imatp2=40: gist_ncar
     & 300,  2, 99,     298, 15, 97,     298, 28, 95,     298, 42, 93,
     & 280, 65, 96,     276, 83,100,     293, 94,100,     318,100,100,
     & 356, 98,100,      10,100,100,      20, 98,100,      35, 95,100,
     &  46, 96,100,      51,100,100,      55,100,100,      64, 93,100,
     &  77, 79,100,      86, 85,100,      90, 98, 97,      90,100, 87,
     & 100,100, 86,     115, 99, 98,     139,100, 98,     159,100, 98,
     & 172,100, 99,     186,100,100,     194,100,100,     220,100,100,
     & 233,100, 73,     164,100, 29,     174,100, 22,     240,100, 50/

      integer, parameter :: n_matp3 = 80
      real matplotlib3(3,n_matp3,4) ! (HSB,i,imatp)
      data matplotlib3 /
C      imatp3=1: flag
     &   0,  0,  0,     240,100, 77,     216, 52,100,      40, 18,100,
     &   3,100,100,       0,100,  1,     240,100, 72,     218, 57,100,
     &  43, 13,100,       0,100, 87,       0,100,  5,     240,100, 67,
     & 220, 63,100,      31, 35,100,       0,100, 92,       0,100,  9,
     & 240,100, 61,     210, 39,100,      33, 30,100,       0,100, 96,
     &   0,100, 14,     240,100, 84,     213, 45,100,      36, 24,100,
     &   0,100,100,       0,  0,  0,     240,100, 79,     215, 50,100,
     &  39, 19,100,       2,100,100,       0,  0,  0,     240,100, 73,
     & 217, 56,100,      42, 14,100,       0,100, 86,       0,100,  4,
     & 240,100, 68,     220, 61,100,      30, 36,100,       0,100, 91,
     &   0,100,  8,     240,100, 63,     209, 38,100,      33, 31,100,
     &   0,100, 95,       0,100, 13,     240,100, 85,     212, 43,100,
     &  35, 26,100,       0,100,100,       0,  0,  0,     240,100, 80,
     & 214, 49,100,      38, 20,100,       1,100,100,       0,  0,  0,
     & 240,100, 75,     217, 54,100,      41, 15,100,       0,100, 85,
     &   0,100,  3,     240,100, 69,     219, 60,100,      30, 38,100,
     &   0,100, 90,       0,100,  7,     240,100, 64,     209, 36,100,
     &  32, 32,100,       0,100, 94,       0,100, 12,     240,100, 86,
     & 211, 42,100,      35, 27,100,       0,100, 98,       0,  0,  0,
     & 240,100, 81,     214, 47,100,      38, 22,100,       0,100,100,
C      imatp3=2: prism
     & 100,100, 99,      66,100,100,      41,100,100,       7,100,100,
     & 342,100,100,     258,100,100,     232,100, 96,     148, 99, 67,
     &  93,100,100,      53,100,100,      23,100,100,       0,100,100,
     & 305,100, 75,     244,100,100,     202,100, 67,     106,100, 91,
     &  75,100,100,      39,100,100,       5,100,100,     339,100,100,
     & 267,100, 97,     230,100, 94,     142, 96, 71,      91,100,100,
     &  59,100,100,      21,100,100,       0,100,100,     298,100, 72,
     & 250,100,100,     221,100, 84,     105,100, 93,      73,100,100,
     &  47,100,100,      14,100,100,     335,100, 97,     264,100, 99,
     & 239,100,100,     171,100, 56,      89,100,100,      58,100,100,
     &  30,100,100,       0,100,100,     292,100, 76,     249,100,100,
     & 217,100, 80,     115, 93, 83,      71,100,100,      45,100,100,
     &  12,100,100,     349,100,100,     262,100,100,     237,100,100,
     & 163,100, 60,      98,100,100,      64,100,100,      28,100,100,
     &   0,100,100,     317,100, 83,     256,100,100,     213,100, 76,
     & 110,100, 86,      80,100,100,      52,100,100,       9,100,100,
     & 346,100,100,     274,100, 91,     243,100,100,     156,100, 63,
     &  96,100,100,      62,100,100,      37,100,100,       0,100,100,
     & 311,100, 79,     254,100,100,     227,100, 91,     108,100, 88,
     &  78,100,100,      50,100,100,      19,100,100,       0,100,100,
C      imatp3=3: tab20b
     & 307, 28, 87,     307, 28, 87,     307, 28, 87,     307, 28, 87,
     & 310, 47, 80,     310, 47, 80,     310, 47, 80,     310, 47, 80,
     & 312, 50, 64,     312, 50, 64,     312, 50, 64,     312, 50, 64,
     & 308, 47, 48,     308, 47, 48,     308, 47, 48,     308, 47, 48,
     & 355, 35, 90,     355, 35, 90,     355, 35, 90,     355, 35, 90,
     & 354, 54, 83,     354, 54, 83,     354, 54, 83,     354, 54, 83,
     & 359, 57, 67,     359, 57, 67,     359, 57, 67,     359, 57, 67,
     &   2, 56, 51,       2, 56, 51,       2, 56, 51,       2, 56, 51,
     &  39, 35, 90,      39, 35, 90,      39, 35, 90,      39, 35, 90,
     &  41, 64, 90,      41, 64, 90,      41, 64, 90,      41, 64, 90,
     &  45, 69, 74,      45, 69, 74,      45, 69, 74,      45, 69, 74,
     &  39, 65, 54,      39, 65, 54,      39, 65, 54,      39, 65, 54,
     &  72, 28, 85,      72, 28, 85,      72, 28, 85,      72, 28, 85,
     &  75, 48, 81,      75, 48, 81,      75, 48, 81,      75, 48, 81,
     &  76, 49, 63,      76, 49, 63,      76, 49, 63,      76, 49, 63,
     &  80, 52, 47,      80, 52, 47,      80, 52, 47,      80, 52, 47,
     & 238, 29, 87,     238, 29, 87,     238, 29, 87,     238, 29, 87,
     & 238, 48, 81,     238, 48, 81,     238, 48, 81,     238, 48, 81,
     & 238, 49, 63,     238, 49, 63,     238, 49, 63,     238, 49, 63,
     & 238, 52, 47,     238, 52, 47,     238, 52, 47,     238, 52, 47,
C      imatp3=4: tab20c
     &   0,  0, 85,       0,  0, 85,       0,  0, 85,       0,  0, 85,
     &   0,  0, 74,       0,  0, 74,       0,  0, 74,       0,  0, 74,
     &   0,  0, 58,       0,  0, 58,       0,  0, 58,       0,  0, 58,
     &   0,  0, 38,       0,  0, 38,       0,  0, 38,       0,  0, 38,
     & 240,  7, 92,     240,  7, 92,     240,  7, 92,     240,  7, 92,
     & 238, 14, 86,     238, 14, 86,     238, 14, 86,     238, 14, 86,
     & 245, 23, 78,     245, 23, 78,     245, 23, 78,     245, 23, 78,
     & 248, 39, 69,     248, 39, 69,     248, 39, 69,     248, 39, 69,
     & 109, 17, 91,     109, 17, 91,     109, 17, 91,     109, 17, 91,
     & 114, 28, 85,     114, 28, 85,     114, 28, 85,     114, 28, 85,
     & 121, 40, 76,     121, 40, 76,     121, 40, 76,     121, 40, 76,
     & 138, 69, 63,     138, 69, 63,     138, 69, 63,     138, 69, 63,
     &  30, 35, 99,      30, 35, 99,      30, 35, 99,      30, 35, 99,
     &  27, 57, 99,      27, 57, 99,      27, 57, 99,      27, 57, 99,
     &  25, 76, 99,      25, 76, 99,      25, 76, 99,      25, 76, 99,
     &  19, 94, 90,      19, 94, 90,      19, 94, 90,      19, 94, 90,
     & 209, 17, 93,     209, 17, 93,     209, 17, 93,     209, 17, 93,
     & 200, 29, 88,     200, 29, 88,     200, 29, 88,     200, 29, 88,
     & 202, 50, 83,     202, 50, 83,     202, 50, 83,     202, 50, 83,
     & 205, 74, 74,     205, 74, 74,     205, 74, 74,     205, 74, 74/


C     -----------------------------------------------------------------------

C     rin0 is from 1 to 3
C     X    is from 0 to 1
C     f: X -> Hue, Sat, Bri
C     rcol(1) is from 1 to 2 [Hue]
C     rcol(2) is from 0 to 1 [Saturation]
C     rcol(3) is from 0 to 1 [Brightness]

      if(icrev/=0) rin0=4d0-rin0
      X=(rin0-1d0)/2d0
      X=min(1d0,X)
      X=max(0d0,X)

      if(ndis .ne. 0) then
        if(X/=1d0) X=1d0*int(X*ndis)/(ndis-1)
      endif

      if(cmap=="phits" .or. cmap=="") then
        if(ndis/=0 .or. icrev/=0) then
          rin0=X*2d0+1d0
          rcol(1) = chue(rin0)
        endif
        return
      endif


C     ---------- change manually ----------
      n=0
      if(cmap=="phits2") n=n_phits2

      imatp1=0
      if(cmap=='viridis'  ) imatp1=1
      if(cmap=='plasma'   ) imatp1=2
      if(cmap=='inferno'  ) imatp1=3
      if(cmap=='magma'    ) imatp1=4
      if(cmap=='cividis'  ) imatp1=5
      if(cmap=='greys'    ) imatp1=6
      if(cmap=='purples'  ) imatp1=7
      if(cmap=='blues'    ) imatp1=8
      if(cmap=='greens'   ) imatp1=9
      if(cmap=='oranges'  ) imatp1=10
      if(cmap=='reds'     ) imatp1=11
      if(cmap=='ylorbr'   ) imatp1=12
      if(cmap=='ylorrd'   ) imatp1=13
      if(cmap=='orrd'     ) imatp1=14
      if(cmap=='purd'     ) imatp1=15
      if(cmap=='rdpu'     ) imatp1=16
      if(cmap=='bupu'     ) imatp1=17
      if(cmap=='gnbu'     ) imatp1=18
      if(cmap=='pubu'     ) imatp1=19
      if(cmap=='ylgnbu'   ) imatp1=20
      if(cmap=='pubugn'   ) imatp1=21
      if(cmap=='bugn'     ) imatp1=22
      if(cmap=='ylgn'     ) imatp1=23
      if(cmap=='binary'   ) imatp1=24
      if(cmap=='gist_yarg') imatp1=25
      if(cmap=='gist_gray') imatp1=26
      if(cmap=='gray'     ) imatp1=27
      if(cmap=='bone'     ) imatp1=28
      if(cmap=='pink'     ) imatp1=29
      if(cmap=='spring'   ) imatp1=30
      if(cmap=='summer'   ) imatp1=31
      if(cmap=='autumn'   ) imatp1=32
      if(cmap=='winter'   ) imatp1=33
      if(cmap=='cool'     ) imatp1=34
      if(cmap=='wistia'   ) imatp1=35
      if(cmap=='hot'      ) imatp1=36
      if(cmap=='afmhot'   ) imatp1=37
      if(cmap=='gist_heat') imatp1=38
      if(cmap=='copper'   ) imatp1=39
      if(imatp1/=0) n=n_matp1

      imatp2=0
      if(cmap=='piyg'            ) imatp2=1
      if(cmap=='prgn'            ) imatp2=2
      if(cmap=='brbg'            ) imatp2=3
      if(cmap=='puor'            ) imatp2=4
      if(cmap=='rdgy'            ) imatp2=5
      if(cmap=='rdbu'            ) imatp2=6
      if(cmap=='rdylbu'          ) imatp2=7
      if(cmap=='rdylgn'          ) imatp2=8
      if(cmap=='spectral'        ) imatp2=9
      if(cmap=='coolwarm'        ) imatp2=10
      if(cmap=='bwr'             ) imatp2=11
      if(cmap=='seismic'         ) imatp2=12
      if(cmap=='twilight'        ) imatp2=13
      if(cmap=='twilight_shifted') imatp2=14
      if(cmap=='hsv'             ) imatp2=15
      if(cmap=='pastel1'         ) imatp2=16
      if(cmap=='pastel2'         ) imatp2=17
      if(cmap=='paired'          ) imatp2=18
      if(cmap=='accent'          ) imatp2=19
      if(cmap=='dark2'           ) imatp2=20
      if(cmap=='set1'            ) imatp2=21
      if(cmap=='set2'            ) imatp2=22
      if(cmap=='set3'            ) imatp2=23
      if(cmap=='tab10'           ) imatp2=24
      if(cmap=='tab20'           ) imatp2=25
      if(cmap=='ocean'           ) imatp2=26
      if(cmap=='gist_earth'      ) imatp2=27
      if(cmap=='terrain'         ) imatp2=28
      if(cmap=='gist_stern'      ) imatp2=29
      if(cmap=='gnuplot'         ) imatp2=30
      if(cmap=='gnuplot2'        ) imatp2=31
      if(cmap=='cmrmap'          ) imatp2=32
      if(cmap=='cubehelix'       ) imatp2=33
      if(cmap=='brg'             ) imatp2=34
      if(cmap=='gist_rainbow'    ) imatp2=35
      if(cmap=='rainbow'         ) imatp2=36
      if(cmap=='jet'             ) imatp2=37
      if(cmap=='turbo'           ) imatp2=38
      if(cmap=='nipy_spectral'   ) imatp2=39
      if(cmap=='gist_ncar'       ) imatp2=40
      if(imatp2/=0) n=n_matp2

      imatp3=0
      if(cmap=='flag'  ) imatp3=1
      if(cmap=='prism' ) imatp3=2
      if(cmap=='tab20b') imatp3=3
      if(cmap=='tab20c') imatp3=4
      if(imatp3/=0) n=n_matp3


C     --------------
      if(n==0)return
C     ---------- change manually ----------
      i=X*(n-1)+1
      dI=X*(n-1d0)+1d0 - i*1d0

C     ---------- change manually ----------
      if(cmap=="phits2") then
        Ylow1=phits2(1,i); Yupp1=phits2(1,i+1)
        Ylow2=phits2(2,i); Yupp2=phits2(2,i+1)
        Ylow3=phits2(3,i); Yupp3=phits2(3,i+1)
      else if(imatp1/=0) then
        Ylow1=matplotlib1(1,i,imatp1); Yupp1=matplotlib1(1,i+1,imatp1)
        Ylow2=matplotlib1(2,i,imatp1); Yupp2=matplotlib1(2,i+1,imatp1)
        Ylow3=matplotlib1(3,i,imatp1); Yupp3=matplotlib1(3,i+1,imatp1)
      else if(imatp2/=0) then
        Ylow1=matplotlib2(1,i,imatp2); Yupp1=matplotlib2(1,i+1,imatp2)
        Ylow2=matplotlib2(2,i,imatp2); Yupp2=matplotlib2(2,i+1,imatp2)
        Ylow3=matplotlib2(3,i,imatp2); Yupp3=matplotlib2(3,i+1,imatp2)
      else if(imatp3/=0) then
        Ylow1=matplotlib3(1,i,imatp3); Yupp1=matplotlib3(1,i+1,imatp3)
        Ylow2=matplotlib3(2,i,imatp3); Yupp2=matplotlib3(2,i+1,imatp3)
        Ylow3=matplotlib3(3,i,imatp3); Yupp3=matplotlib3(3,i+1,imatp3)
      else
        return
      endif
C     ---------- change manually ----------

C     ----- Hue -----
      if    (Yupp1-Ylow1 > 180d0) then
        Ylow1=Ylow1+360d0
      elseif(Ylow1-Yupp1 > 180d0) then
        Yupp1=Yupp1+360d0
      endif
      rcol(1)=Ylow1+(Yupp1-Ylow1)*dI
      if(rcol(1) > 360d0) rcol(1)=rcol(1)-360d0
      rcol(1)=rcol(1)/360d0+1d0
C     ----- Saturaion -----
      rcol(2)=Ylow2+(Yupp2-Ylow2)*dI
      rcol(2)=rcol(2)/100d0
C     ----- Brightness -----
      rcol(3)=Ylow3+(Yupp3-Ylow3)*dI
      rcol(3)=rcol(3)/100d0

C     ----------------------------------------
C import matplotlib.pyplot as plt
C import matplotlib.colors as mcolors
C import colorsys
C import numpy as np

C turn=4
C cmap_list, n_divisions = {}, {}
C n_divisions[1]=turn*5
C cmap_list[1]=[
C             plt.cm.viridis,
C             plt.cm.plasma,
C             ]
C n_divisions[2]=turn*8
C cmap_list[2]=[
C             plt.cm.PiYG,
C             plt.cm.PRGn,
C             ]
C n_divisions[3]=turn*20
C cmap_list[3]=[
C             plt.cm.flag,
C             plt.cm.prism,
C             ]

C for i in [1,2,3]:
C     print("      integer, parameter :: n_matp{} = {}".format(i,n_divisions[i]))
C     print("      real matplotlib{}(3,n_matp{},{}) ! (HSB,i,imatp)".format(i,i,len(cmap_list[i])))
C     print("      data matplotlib{} /".format(i))
C     for idata, cmap in enumerate(cmap_list[i]):
C         values = np.linspace(1, 0, n_divisions[i])
C         colors_rgb = [cmap(v)[:3] for v in values]
C         # RGB -> HSV
C         colors_hsb = []
C         for r, g, b_ in colors_rgb:
C             h, s, b = colorsys.rgb_to_hsv(r, g, b_)
C             colors_hsb.append((360*h, 100*s, 100*b))
C         colors_hsb=np.array(colors_hsb).astype(int)
C         # Results
C         print("C      imatp{}={}: {}".format(i,idata+1,cmap.name.lower()))
C         for ihsb,hsb in enumerate(colors_hsb):
C             if ihsb%turn==0:
C                 print("     &",end="")
C             if idata==len(cmap_list[i])-1 and ihsb==len(colors_hsb)-1:
C                 sepa="/"
C                 end=""
C             else:
C                 sepa=","
C                 if (1+ihsb)%turn==0:
C                     end="\n"
C                 else:
C                     end="    "
C             print(" {:3},{:3},{:3}{}".format(hsb[0],hsb[1],hsb[2],sepa),end=end)
C     print("\n")

C for i in [1,2,3]:
C     l_name_max=max(len(cmap.name) for cmap in cmap_list[i])
C     print("      imatp{}=0".format(i))
C     for idata, cmap in enumerate(cmap_list[i]):
C         l_name=len(cmap.name)
C         print("      if(cmap=='{}'{}) imatp{}={}".format(cmap.name.lower() ," "*(l_name_max-l_name),i ,idata+1))
C     print("      if(imatp{}/=0) n=n_matp{}".format(i,i))
C     print()

C for i in [1,2,3]:
C     print("      else if(imatp{}/=0) then".format(i))
C     print("        Ylow1=matplotlib{0}(1,i,imatp{0}); Yupp1=matplotlib{0}(1,i+1,imatp{0})".format(i))
C     print("        Ylow2=matplotlib{0}(2,i,imatp{0}); Yupp2=matplotlib{0}(2,i+1,imatp{0})".format(i))
C     print("        Ylow3=matplotlib{0}(3,i,imatp{0}); Yupp3=matplotlib{0}(3,i+1,imatp{0})".format(i))
C     ----------------------------------------
      return
      end
