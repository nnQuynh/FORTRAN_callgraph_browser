************************************************************************
*                                                                      *
*     'mmbank.inc'                                                     *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
cKN 2015/03/19 nsos, nsosa, ibksos, iaksos

      common /bnkptb1/ibkblz,ibknam,ibknmd,ibknty,ibknkf,ibknfc,
     &                ibknct
      common /bnkptb2/ibkwt,ibku,ibkv,ibkw,ibke,ibkx,ibky,ibkz,ibkt,
     &                ibkwin,ibkwnz,ibkxfc
      common /bnkptb3/ibkspx,ibkspy,ibkspz,ibkzst,ibksos

      common /bnkptbc/ibkec,ibkxc,ibkyc,ibkzc,ibktc

      common /bnkpta1/iakblz,iaknam,iaknmd,iaknty,iaknkf,iaknfc,
     &                iaknct
      common /bnkpta2/iakwt,iaku,iakv,iakw,iake,iakx,iaky,iakz,iakt,
     &                iakwin,iakwnz,iakxfc
      common /bnkpta3/iakspx,iakspy,iakspz,iakzst,iaksos

*-----------------------------------------------------------------------

      dimension iblz(1),name(1),nmed(1),nty(1),nkf(1),nfcs(1),
     &          ncnt(3,1),nzst(1),nsos(1)
      dimension wt(1),u(1),v(1),w(1),e(1),x(1),y(1),z(1),t(1),
     &          wtin(1),wtnz(1),xfcs(1)
      dimension ec(1),xc(1),yc(1),zc(1),tc(1)
      dimension spx(1),spy(1),spz(1)

      dimension iblza(1),namea(1),nmeda(1),ntya(1),nkfa(1),nfcsa(1),
     &          ncnta(3,1),nzsta(1),nsosa(1)
      dimension wta(1),ua(1),va(1),wa(1),ea(1),xa(1),ya(1),za(1),ta(1),
     &          wtina(1),wtnza(1),xfcsa(1)
      dimension spxa(1),spya(1),spza(1)

*-----------------------------------------------------------------------

      equivalence ( das, iblz,name,nmed,nty,nkf,nfcs,ncnt,nzst,nsos )
      equivalence ( das, wt,u,v,w,e,x,y,z,t )
      equivalence ( das, wtin,wtnz,xfcs )
      equivalence ( das, ec,xc,yc,zc,tc )
      equivalence ( das, spx,spy,spz )

      equivalence ( das, iblza,namea,nmeda,ntya,nkfa,nfcsa,ncnta,
     &                   nzsta,nsosa)
      equivalence ( das, wta,ua,va,wa,ea,xa,ya,za,ta )
      equivalence ( das, wtina,wtnza,xfcsa )
      equivalence ( das, spxa,spya,spza )

*-----------------------------------------------------------------------
