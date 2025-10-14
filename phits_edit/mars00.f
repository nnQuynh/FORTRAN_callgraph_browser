      function aread ( ibr, iret )
c
c     entry to return hollerith data in four byte words
c     ibr = 0  return next 4 characters on card
c     ibr = 1  read a new card and return first 4 characters
c     ibr = 2  return 4 characters starting with the first non blank
c              character and stopping with 2 consecutive blanks
c     ibr = 3  return next 4 characters, but stop at a double blank
c     ibr = 4  same as ibr = 2, but stop at a single blank
c     ibr =-1  returns last character scanned
c
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     94/11/28
c     time the module was last permanently updated:     14:13:09
c     programmer name:                                  l.m.petrie
c     module name:                                      s7aread
c     current archiving level number:                   00011
c     current number of permanent updates:              00011
c     date of last access by librarian:                 94/11/28
c     dataset name:  x4s.scale4.master
c
c#######################################################################
      external qrdbfr
      logical lskip, lascn
      character aread*4, equals*1
c
      common /unit/ inpt, outpt, idum(19)
      common /qrdbuf/ unit, tens, dnum, ichr, lchr, irpt,
     *  wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
      common /yrdbuf/ char, card
      save /unit/, /qrdbuf/, /yrdbuf/
c
      integer inpt, outpt
      double precision dnum, unit, tens, ten
      character*1 card, char, blank, e, d, comma
      logical wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
c
      dimension card(252), char(28)
      dimension unit(10), tens(8)
c
      save infg
!$OMP THREADPRIVATE(infg)
c
c     data statements
c
      data iblnk/1/, equals/'='/
      save iblnk !FURUTA
!$OMP THREADPRIVATE(iblnk)
c
c     end specification statements
c
      aread  = char(17) !FURUTA
      if ( iret.eq.2) then
        if (ibr.eq.3 .and. wrdflg ) then
          return
        else
          call enfile
C
CTAT/MOD(97/Oct) TO AVOID SUDDEN DEATH ---------------------------(FROM)
          return
CTAT/MOD(97/Oct) TO AVOID SUDDEN DEATH ---------------------------( TO )
C
        end if
      end if
      lascn  = lscan .and. ibr.gt.1
      iret   = 0
      if (ibr.eq.-1) then
c
c     return single character
c
        aread  = char(17) !FURUTA
        if ( illchr ) then
          ichr   = ichr + 1
          illchr = .false.
        endif
        if (ichr.gt.0.and.ichr.le.lchr) aread  = card(ichr)
        return
      end if
c
      if (ibr.eq.2)            iblnk  = 1
      if (ibr.eq.4)            iblnk  = 0
      lskip  = ibr.eq.2 .or. ibr.eq.4
      if (ibr.ne.1)            go to 110
  100 continue
      ichr   = 0
      call y0read ( card, lchr, iret )
      if (iret.ne.0)           return
  110 continue
      if (.not.lskip)          go to 120
      ichr   = ichr + 1
      if (ichr.gt.lchr)        go to 100
      if (card(ichr).eq.char(17)) go to 110 !FURUTA
      ichr   = ichr - 1
  120 continue
      nchr   = 0
      aread  = char(17) !FURUTA
      irpt   = 0
      if (ibr.eq.3)            go to 130
      if (ichr.ge.lchr)        go to 100
      wrdflg = .false.
      infg   = 0
  130 continue
      nchr   = nchr + 1
      if (wrdflg.or.nchr.gt.4) go to 150
      ichr   = ichr + 1
      if (ichr.gt.lchr) then
        wrdflg = ibr.gt.1
        go to 150
      end if
      infg   = infg + 1
      if (card(ichr).ne.char(17)) infg   = 0 !FURUTA
      aread(nchr:nchr) = card(ichr)
      wrdflg = (infg.gt.iblnk .or. card(ichr).eq.equals) .and. ibr.gt.1
      go to 130
c
c     test for scan or return
c
  150 continue
      if (lascn) then
        jchr   = ichr
        jnfg   = infg
c
c     scan ahead for end
c
  160   continue
        if (ichr.ge.lchr) then
          wrdflg = ibr.gt.1
          ichr   = 0
          call y0read ( card, lchr, iret )
          if (iret.ne.0)           return
        end if
        if (card(ichr+1).eq.char(17)) then !FURUTA
          ichr   = ichr+1
          infg   = infg+1
          wrdflg = infg.gt.iblnk .or. wrdflg
        else
          if (ichr+3.gt.lchr)          return
          if (card(ichr+1)//card(ichr+2)//card(ichr+3).eq.'end') then
            iret   = 1
            return
          else
            iret   = 0
            if (.not.wrdflg) then
              ichr   = jchr
              infg   = jnfg
            end if
          end if
          return
        end if
        go to 160
      end if
      return
      end
      subroutine clear ( ld, n )
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     91/06/20
c     time the module was last permanently updated:     11:29:43
c     programmer name:                                  l.m.petrie
c     module name:                                      s7clear
c     current archiving level number:                   00004
c     current number of permanent updates:              00004
c     date of last access by librarian:                 93/02/24
c     dataset name:  x4s.scale4.master
c
c#######################################################################
      dimension ld(*)
      if (n.le.0) return
      do 100 i=1,n
         ld(i) = 0
  100    continue
      return
      end
      function cread ( ibr, iret )
c
c     cread returns an eight character word according to
c     the argument ibr ( see comments in aread )
c
c     type statements
c
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     94/11/28
c     time the module was last permanently updated:     15:28:54
c     programmer name:                                  l.m.petrie
c     module name:                                      s7cread
c     current archiving level number:                   00006
c     current number of permanent updates:              00006
c     date of last access by librarian:                 94/11/28
c     dataset name:  x4s.scale4.master
c
c#######################################################################
      character*8 cread
      character*4 aread,first,last
c
      ib2 = 0
      if (ibr.gt.1) ib2 = 3
      first = aread ( ibr, iret )
      last  = aread ( ib2, iret )
      cread = first // last
      return
      end
      function dread ( ibr, iret )
c
c     dread returns a double precision floating point number
c     see comments for aread
c
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     94/11/28
c     time the module was last permanently updated:     14:16:15
c     programmer name:                                  l.m.petrie
c     module name:                                      s7dread
c     current archiving level number:                   00025
c     current number of permanent updates:              00025
c     date of last access by librarian:                 94/11/28
c     dataset name:  x4s.scale4.master
c
c#######################################################################
c
      common /unit/ inpt, outpt, idum(19)
      common /qrdbuf/ unit, tens, dnum, ichr, lchr, irpt,
     *  wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
      common /yrdbuf/ char, card
      save /unit/, /qrdbuf/, /yrdbuf/
c
      integer inpt, outpt
      double precision dnum, unit, tens, ten
      character*1 card, char, blank, e, d, comma
      logical wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
c
      dimension card(252), char(28)
      dimension unit(10), tens(8)
c
      logical errflg, expflg, numsgn, expsgn, perflg, zflg
      character*4 end
      double precision dread, di, frac
      save frac, iper, num, expsgn, numsgn
!$OMP THREADPRIVATE(frac, iper, num, expsgn, numsgn)
c
c     data statements
c
      data nexp/3/, mant/17/, end/'end '/
c
c     end specification statements
c
C
CTAT/MOD(97/Oct) TO AVOID SUDDEN DEATH ---------------------------(FROM)
      if (iret.eq.2) then
          call enfile
          return
      endif
CTAT/MOD(97/Oct) TO AVOID SUDDEN DEATH ---------------------------( TO )
C
      infg   = 0
      nsig   = mant
      expflg = .false.
      errflg = .false.
      wrdflg = .false.
      zflg   = .false.
      if (ibr.ne.1)     go to 110
  100 continue
      if (infg.gt.0)    go to 190
      ichr   = 0
      iret   = 0
      irpt   = 0
      call y0read ( card, lchr, iret )
      if (iret.ne.0) return
  110 continue
      if (irpt.gt.0)    go to 270
      if (iret.ne.0)    go to 260
      dnum   = 0
      frac   = 0
      iper   = 0
      num    = 0
      numsgn = .false.
      expsgn = .false.
      perflg = .false.
      illchr = .false.
      lneg   = .false.
  120 continue
      ichr   = ichr + 1
      if (ichr.gt.lchr) go to 100
  130 continue
      do 140 i=1,28
         if (char(i).eq.card(ichr)) go to (150,150,150,150,150,150,150,
     *    150,150,150,160,160,160,170,170,160,180,190,220,220,220,200,
     *    200,210,240,250,190,230),i
  140    continue
      go to 160
c
c     numeric character
c
  150 continue
      infg   = infg + 1
      if (expflg) then
        num    = 10*num + i - 1
        if (infg.ge.nsig) then
          if (ichr.lt.lchr .and. card(ichr+1).eq.char(18)) !FURUTA
     &         ichr   = ichr+1                             !FURUTA
          go to 190
        end if
      else
        if (infg.le.nsig) then
          di     = i - 1
          if (perflg) then
            frac   = unit(2)*frac + di
            iper   = iper + 1
          else
            dnum   = unit(2)*dnum + di
          end if
        end if
      end if
      go to 120
c
c     illegal character
c
  160 continue
      illchr = .true.
      if ( .not.lallc) then
        if (.not.errflg) then
          erread = .true.
          errflg = .true.
          if (crdflg) write(outpt,10000) (card(ii),ii=1,lchr)
          crdflg = .false.
          write(outpt,10100) ichr, card(ichr)
          call errtra
        end if
        go to 120
      else
        ichr   = ichr - 1
        go to 190
      end if
c
c     exponent flag - d or e
c
  170 continue
      expflg = .true.
      nsig   = nexp
      infg   = 1
      if (ichr+1.le.lchr.and.card(ichr+1).eq.char(17)) ichr = ichr + 1 !FURUTA
      go to 120
c
c     blank
c
  180 continue
      if (infg.eq.0)    go to 120
c
c     comma - field terminator
c
  190 continue
      it     = iper/10
      iu     = mod(iper,10)
      frac   = frac/unit(iu+1)/tens(it+1)
      it     = num/10
      iu     = mod(num,10)
      dnum   = dnum + frac
      if (expsgn) then
        dnum = dnum/unit(iu+1)/tens(it+1)
      else
        dnum = dnum*unit(iu+1)*tens(it+1)
      end if
      if (numsgn) dnum   = -dnum
      go to 270
c
c     plus sign
c
  200 continue
      if (infg.le.0) then
        infg   = infg + 1
      else
        expflg = .true.
        nsig   = nexp
        infg   = 1
      end if
      go to 120
c
c     minus sign
c
  210 continue
      if (infg.le.0) then
        numsgn = .true.
        infg   = infg + 1
      else
        expflg = .true.
        expsgn = .true.
        nsig   = nexp
        infg   = 1
      end if
      go to 120
c
c     repeat flag - r,*, or $
c
  220 continue
      irpt   = dnum
      dnum   = 0
      frac   = 0
      numsgn = .false.
      expsgn = .false.
      expflg = .false.
      perflg = .false.
      iper   = 0
      if (.not.zflg) infg   = 0
      go to 120
c
c     p - alternating sign repeat
c
  230 continue
      lneg   = .true.
      go to 220
c
c     z - zero field
c
  240 continue
      zflg   = .true.
      infg   = max0(infg,1)
      go to 220
c
c     decimal point
c
  250 continue
      perflg = .true.
      go to 120
  260 continue
      dnum   = 0.0
  270 continue
      dread  = dnum
      irpt   = irpt - 1
      if (lneg) dnum   = -dnum
  280 continue
      if (iret.ne.0)  return
      if (irpt.gt.0)  return
      if (.not.lscan) return
      if (illchr)     return
      if (ichr.ge.lchr) then
        call y0read ( card, lchr, iret )
        if (iret.ne.0) return
        ichr   = 0
      end if
      ie     = ichr + 1
      if (card(ie).eq.char(17)) then !FURUTA
        ichr   = ichr + 1
        go to 280
      end if
      in     = ie + 1
      id     = in + 1
      if (id.gt.lchr) return
      if (card(ie)//card(in)//card(id).ne.end) return
      iret   = 1
      return
c
c     formats
c
10000 format ( '0*****error in input.',
     *         ' card image printed on next line *****' / 10x,80a1 )
10100 format ( ' on the above card, character number ', i2,
     *         ' (image=', a1, ') is not valid.' )
      end
      double precision function dset ( x )
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     91/06/20
c     time the module was last permanently updated:     11:29:47
c     programmer name:                                  nmg
c     module name:                                      s7dset
c     current archiving level number:                   00002
c     current number of permanent updates:              00002
c     date of last access by librarian:                 93/02/24
c     dataset name:  x4s.scale4.master
c
c#######################################################################
      double precision x
      dset   = x
      return
      end
      subroutine enfile
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     94/11/28
c     time the module was last permanently updated:     14:17:09
c     programmer name:                                  l.m.petrie
c     module name:                                      s7enfile
c     current archiving level number:                   00007
c     current number of permanent updates:              00007
c     date of last access by librarian:                 94/11/28
c     dataset name:  x4s.scale4.master
c
c#######################################################################
C
CTAT/ADD(97/Oct) TO AVOID SUDDEN DEATH ---------------------------(FROM)
      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)
CTAT/ADD(97/Oct) TO AVOID SUDDEN DEATH ---------------------------( TO )
C
      common /unit/ inpt, outpt, idum(19)
      integer inpt, outpt, idum
      save /unit/
c
c     print end of file message and stop
c
      write(outpt,10000) inpt
      call errtra
C
CTAT/MOD(97/Oct) TO AVOID SUDDEN DEATH ---------------------------(FROM)
      return
CTAT/MOD(97/Oct) TO AVOID SUDDEN DEATH ---------------------------( TO )
C
10000 format ( 33h0 ***** end of file read on unit ,i2,6h ***** )
      end
      subroutine errtra
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     none
c     time the module was last permanently updated:     none
c     programmer name:                                  l.m.petrie
c     module name:                                      ulerrtra
c     current archiving level number:                   00000
c     current number of permanent updates:              00000
c     date of last access by librarian:                 93/02/24
c     dataset name:  x4s.scale4.master
c
c#######################################################################
C
CTAT/MOD(97/Oct) TO AVOID SUDDEN DEATH ---------------------------(FROM)
      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)
      if(ierror.eq.0) ierror=1
CTAT/MOD(97/Oct) TO AVOID SUDDEN DEATH ---------------------------( TO )
      return
      end
      subroutine ffpack ( ibase,ichr,ipos )
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     91/06/20
c     time the module was last permanently updated:     11:29:53
c     programmer name:                                  l.m.petrie
c     module name:                                      s7ffpack
c     current archiving level number:                   00006
c     current number of permanent updates:              00006
c     date of last access by librarian:                 93/02/24
c     dataset name:  x4s.scale4.master
c
c#######################################################################
      character ibase*8,ichr*1
      ibase(ipos:ipos) = ichr
      return
      end
      subroutine ffread(in,k,v,nf,n5,n6,iprtrg)
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     93/02/24
c     time the module was last permanently updated:     07:25:53
c     programmer name:                                  l.m.petrie
c     module name:                                      s7ffread
c     current archiving level number:                   00015
c     current number of permanent updates:              00015
c     date of last access by librarian:                 93/02/24
c     dataset name:  x4s.scale4.master
c
c#######################################################################
c
c****** free- and fixed-field card image translator
c
c   *** ncpw is the characters per word
      parameter(ncpw=4)
      dimension ny(80),in(*),k(*),v(*)
      dimension numb(10)
      character*1 numb,nt,nbb,ndp,npl,npa,nmi,nss,ndl,ns,ne,nq,nm,nn
      character*1 nr,ni,nff,na,nl,nk,nsls,ncom,nsl,ndb,nh,k,ny,nyncc
      character*8 blank,buffer
c
c ***** double precision for ibm 360
c
      double precision v,dnum,vnum,sgn,dum,one,ten
      data one/1.0d0/,ten/10.0d0/
c
c ***** end of double precision for ibm 360
c
      data numb/'0','1','2','3','4','5','6','7','8','9'/
      data nt,nbb,ndp,npl,npa,nmi/'t',' ','.','+','&','-'/
      data nss,ndl,ns,ne/'s','$','*','e'/
      data nq,nm,nn/'q','m','n'/
      data nr,ni,nff,na,nl,nk/'r','i','f','a','l','k'/
      data nsls,ncom,nsl,ifree/'/',',','''',0/
      data ndb/'#'/,nh/'h'/,blank/'        '/
      save ifree !FURUTA

  100 do 110 i=1,37
         in(i)  = 0
         v(i)   = 0.0
  110    k(i)   = nbb
      iexp   = 0
      nscl   = 0
      ndpn   = 0
      iex    = 0
      vnum   = 0.0
      sgn    = one
      iblnk  = 0
      nf     = 1
      ncc    = 1
      ifc    = 0
      ids    = 0
      nman   = 0
  120 read(n5,10100) ny
      call y0trns ( ny, 80 )
      if (iprtrg.gt.0)                             write(n6,10000) ny
      if (ny(1).ne.nsl)                            go to 130
      write(n6,10100) nbb,(ny(i),i=2,80)
      go to 120
  130 nyncc  = ny(ncc)
      if (nyncc.eq.ncom)                           nyncc  = nbb
      if (ifree.eq.0)                              go to 140
      ifc    = ifc+1
      if (ifc.gt.12)                               ifc    = 1
      if (ifc.eq.3)                                go to 150
      if (ifc.eq.11.and.iex.ne.0.and.nyncc.eq.nbb) go to 320
      if (nyncc.eq.ne.and.ifc.gt.4)                go to 240
  140 if (nyncc.eq.nbb)                            go to 280
      if (nyncc.eq.npl.or.nyncc.eq.npa)            go to 230
      if (nyncc.eq.nmi)                            go to 250
      if (nyncc.eq.ne.and.iblnk.eq.1)              go to 240
      iblnk  = 1
      if (nyncc.eq.ndp)                            go to 200
  150 do 160 num=1,10
         if (numb(num).eq.nyncc)                   go to 190
  160    continue
  170 if (k(nf).ne.nbb.and.nyncc.eq.nbb)           go to 280
      k(nf)  = nyncc
      in(nf) = sgn*vnum
      vnum   = 0.0
      sgn    = one
      iblnk  = 0
      nman   = 0
      if (nyncc.eq.nbb)                            go to 320
      if (k(nf).ne.ndl.and.k(nf).ne.ns.and.k(nf).ne.ndb) go to 180
      ndbl   = 1
      if (k(nf).eq.ndb)                            ndbl   = 2
      ifree  = 0
      ids    = 1
      if (nf.eq.1)                                 go to 300
      if (k(nf).ne.k(nf-1))                        go to 300
      ids    = 0
      k(nf)  = nbb
      go to 320
  180 continue
      if (nyncc.eq.nh)                             go to 340
      if ( ifree.eq.0  .and. nyncc.ne.nsls .and. nyncc.ne.nr .and.
     *     nyncc.ne.ni .and. nyncc.ne.nff  .and. nyncc.ne.na .and.
     *     nyncc.ne.nl .and. nyncc.ne.nk   .and. nyncc.ne.nq .and.
     *     nyncc.ne.nm .and. nyncc.ne.nn )         go to 300
      if ( ifree.eq.1  .and. nyncc.eq.nss )        go to 300
      if (ncc.eq.72)                               write(n6,10200) ny
      go to 320
  190 continue
      if (iex.ne.0)                                go to 210
      dnum   = num-1
      vnum   = ten*vnum + dnum
      nscl   = nscl + ndpn
      nman   = nman + 1
      go to 220
  200 ndpn   = -1
      go to 220
  210 continue
      iexp   = iexp*10 + iex*(num-1)
  220 if (ifc.eq.12)                               go to 290
      go to 320
  230 if (nman.gt.0)                               go to 240
      sgn    = one
      go to 320
  240 iex    = 1
      go to 320
  250 if (nman.gt.0)                               go to 260
      sgn    = -one
      go to 320
  260 iex    = -1
      go to 320
  270 ifree  = 1
      ids    = 0
      go to 310
  280 if (ids.eq.1)                                go to 270
      if (ifree.eq.1 .and. ifc.eq.12 .and. k(nf).ne.nbb) go to 290
      if (iblnk.eq.0)                              go to 320
      if (ifree.eq.1 .and. ifc.ne.12)              go to 320
  290 vnum   = sgn*vnum
      dum    = ten**iabs(iexp+nscl)
      if (iexp+nscl.ge.0)                   v(nf)  = vnum*dum
      if (iexp+nscl.lt.0)                   v(nf)  = vnum/dum
      nman   = 0
  295 ndpn   = 0
      iex    = 0
      iexp   = 0
      nscl   = 0
      iblnk  = 0
      vnum   = 0.0
      sgn    = one
  300 nf     = nf + 1
      if (nyncc.eq.nss)                            go to 320
  310 ifc    = 0
      if (ifree.eq.1)                       ncc    = ((ncc+11)/12)*12
      if (k(nf-1).eq.nt)                           go to 330
  320 if (ncc.ge.73)                               go to 330
      ncc    = ncc+1
      if (ncc.le.72)                               go to 130
      go to 280
  330 nf     = nf-1
      if (nf.le.0)                                 go to 100
      return
  340 continue
      nhol   = min0(in(nf),72-ncc)
      if (nhol.ne.in(nf)) write(n6,10200) ny
      in(nf) = nhol
      if (nhol.le.0)                               go to 295
      nccmax = ncc + nhol
      kcpw   = ndbl*ncpw
      nwds   = (nhol+kcpw-1)/kcpw
      do 370 kh=1,nwds
         buffer = blank
         do 350 kc=1,kcpw
            ncc    = ncc+1
            if (ncc.gt.nccmax)                     go to 360
            nyncc  = ny(ncc)
            call ffpack ( buffer,nyncc,kc )
  350       continue
  360    continue
         read(buffer,10300) vnum
         k(nf)  = nh
         v(nf)  = vnum
         nyncc  = ny(ncc)
         if (ncc.gt.nccmax)                        go to 295
         nf     = nf+1
  370    continue
      nf     = nf-1
      go to 295
10000 format(' ',t53,80a1)
10100 format(80a1)
10200 format('0 ****** incomplete field at end of card,',
     * ' some data may be lost'/'         card image follows'/1x,80a1)
10300 format(a8)
      end
      subroutine fidas (d,ld,lld,ll2,j3,n5,n6)
      external fidcom
c
c   *** data array modification
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     93/09/20
c     time the module was last permanently updated:     07:31:12
c     programmer name:                                  l.m.petrie
c     module name:                                      s7fidas
c     current archiving level number:                   00020
c     current number of permanent updates:              00020
c     date of last access by librarian:                 93/09/20
c     dataset name:  x4s.scale4.master
c
c#######################################################################
      dimension in(37),k(37),v(37),ld(*),vmt(18),w(12),
     *          d(*),lld(*)
      common /fidasc/ iprtrg
      dimension prt(3)
      equivalence (bb,lbb)
      dimension vs(74)
      equivalence (v(1),vs(1))
      double precision v,vv,vu,del,delt,dif
      character*1 k,lbb,ldl,ls,lr1,li,lt,lss,lf,la,lap,lpl,lmi,lz,lv,lu
      character*1 bb,lq,lm,ln,le,lli,lki,lo,lc,lsls,ldb,lhl,lb
      character*4 prt,vmt
      dimension dv(2),vi(2)
      equivalence (dv(1),delt),(vi(1),vv)
      data ldl,ls,lr1,li,lt,lss/'$','*','r','i','t','s'/
      data lf,la,lap,lpl,lmi/'f','a','&','+','-'/
      data lz,lv,lu,bb/'z','v','u',' '/
      data lq,lm,ln/'q','m','n'/,le/'e'/,lli/'l'/
      data lki/'k'/
      data lo,lc,prt(1),prt(2),prt(3)/'o','c','off ','prt ','on  '/
      data lsls /'/'/
      data ldb/'#'/, lhl/'h'/
      kdbl   = 0
      ldbl   = 0
      j      = 0
      j3     = 0
      iii    = 0
      ivmt   = 0
      ncount = 0
  100 call ffread(in,k,v,nf,n5,n6,iprtrg)
      do 120 i=1,nf
         if (k(i).ne.lpl.and.k(i).ne.lap.and.k(i).ne.lmi) go to 120
c   *** exponentiation (&,+,-)
         l      = in(i)
         if (l.eq.0) go to 120
         e      = 10.0**l
         if (k(i).eq.lmi) go to 110
         v(i)   = v(i)*e
         go to 120
  110    v(i)   = v(i)/e
  120    continue
      i      = 1
  130 continue
         if (iii.eq.0)                 go to 140
         if (iii-2) 570,510,540
  140    iii    = 0
         if (k(i).eq.lbb)              go to 370
         if (k(i).eq.lhl)              go to 370
         if (k(i).eq.ldl)              go to 200
         if (k(i).eq.ls)               go to 190
         if (k(i).eq.ldb)              go to 180
         if (k(i).eq.lr1)              go to 400
         if (k(i).eq.li)               go to 560
         if (k(i).eq.lt)               go to 150
         if (k(i).eq.lss)              go to 340
         if (k(i).eq.lf)               go to 290
         if (k(i).eq.la)               go to 330
         if (k(i).eq.lz)               go to 390
         if (k(i).eq.le)               go to 320
         if (k(i).eq.lq)               go to 420
         if (k(i).eq.lm.or.k(i).eq.ln) go to 460
         if (k(i).eq.lu)               go to 160
         if (k(i).eq.lv)               go to 170
         if (k(i).eq.lli)              go to 500
         if (k(i).eq.lki)              go to 530
         if (k(i).eq.lo)               go to 350
         if (k(i).eq.lc)               go to 360
         if (k(i).eq.lsls)             go to 210
         go to 370
c   *** terminate (t)
  150    j2     = 0
         itest  = in(i)
         kdbl   = ldbl
         if (j.eq.0) go to 250
         go to 220
c   *** variable format control (u,v)
  160    read(n5,10900)vmt
  170    ivmt   = 1
         go to 190
c    *** begin new array (*,$,#)
  180    kdbl   = ldbl
         ldbl   = 1
         kkk    = 0
         go to 210
  190    kkk    = 0
         kdbl   = ldbl
         ldbl   = 0
         go to 210
  200    kkk    = 1
         kdbl   = ldbl
         ldbl   = 0
  210    if (ncount+j.eq.0) go to 270
         j2     = 1
  220    ij     = j
         ic     = ncount
         inc    = 1
         if (kdbl.eq.0) go to 230
         ij     = j/2
         ic     = ncount/2
         inc    = 2
  230    write(n6,10400) ll,lb,ij
         if (j.eq.ncount) go to 250
         imax   = j1 + j - 1
         if (j.gt.10000) imax   = j + 99
         if (kk.eq.1) then
           write(n6,'(1x,10i12)')      (ld(ii),ii=j1,imax,inc)
         else
           write(n6,'(1x,1p,10e12.5)') (d(ii),ii=j1,imax,inc)
         end if
         j3     = j3 + 1
         write(n6,10500) ic,ll,lb
  250    if (j2.ne.0) go to 270
         write(n6,10600)in(i),lt
  260    return
  270    if (k(i).eq.lsls) go to 590
         ll     = in(i)
         ll3    = ll + ll2 - 1
         kk     = kkk
         lb     = k(i)
         inc    = ldbl+1
         j      = inc*v(i)
         j1     = lld(ll3) + j
         ncount = lld(ll3+ 1) - j1
         j      = 0
         if (ivmt.eq.0) go to 600
c   *** variable format control (u,v)
  280    if (ncount.eq.0) go to 100
         j      = ncount
         ncc    = j1 + ncount - 1
         read(n5,vmt) (d(j2),j2=j1,ncc)
         ivmt   = 0
         go to 100
c   *** fill array (f)
  290    if (j.ge.ncount) go to 310
         ncc    = j+1
         do 300 ii=ncc,ncount,inc
            j2     = j1 + ii -1
            if (ldbl.eq.1) then
              d(j2)   = vs(2*i-1)
              d(j2+1) = vs(2*i)
            else
              d(j2)   = v(i)
            end if
            if (kk.ne.0)  ld(j2) = v(i)
  300       continue
         j      = ncount
         go to 600
  310    write(n6,10700) ll,lb
         go to 600
c   *** end array (e)
  320    if (j.le.ncount) j=ncount
         go to 600
c   *** address modification (a)
  330    j      = v(i)
         if (j.gt.ncount.or.j.le.0) write(n6,10800)j,ll,lb
         j      = inc*(j-1)
         go to 600
c   *** skip (s)
  340    j      = j + inc*in(i)
         go to 600
c   *** turn print trigger on/off
  350    iprtrg = -iprtrg
         write(n6,10100)prt(2),prt(iprtrg+2),in(i),lo
         go to 600
c   *** print count in current array
  360    write(n6,10200)ll,lb,j,in(i),lc
         go to 600
c   *** no modification
  370    continue
  380    j2     = j1 + j
         j      = j+inc
         if (ldbl.eq.1) then
           d(j2)   = vs(2*i-1)
           d(j2+1) = vs(2*i)
         else
           if (k(i).eq.lhl) then
             d(j2)   = vs(2*i-1)
             go to 600
           else
             d(j2)   = v(i)
           end if
         end if
         if (kk.ne.0) ld(j2) = v(i)
         go to 600
c   *** zero (z)
  390    in(i)  = v(i) + in(i)
         v(i)   = 0.0
c   *** repeat (r)
  400    l      = in(i)
         do 410 ii=1,max0(1,l)
            j2     = j1 + j
            if (ldbl.eq.1) then
              d(j2)   = vs(2*i-1)
              d(j2+1) = vs(2*i)
            else
              d(j2)   = v(i)
            end if
            if (kk.ne.0)  ld(j2) = v(i)
  410       j      = j+inc
         go to 600
c   *** sequence repeat (q)
  420    l      = v(i) + in(i)
         lseq   = 1
         if (v(i).eq.0.0.or.in(i).eq.0) go to 430
         l      = v(i)
         lseq   = in(i)
  430    do 450 lsq=1,lseq
            do 440 ii=1,l
               j2     = j1 + j
               j4     = j2 - l*inc
               d(j2)  = d(j4)
               if (ldbl.eq.1) d(j2+1) = d(j4+1)
  440          j      = j+inc
  450       continue
         go to 600
c   *** inverted sequence repeat (n)
  460    l      = v(i) + in(i)
         lseq   = 1
         if (v(i).eq.0.0.or.in(i).eq.0) go to 470
         l      = v(i)
         lseq   = in(i)
  470    do 490 lsq=1,lseq
            do 480 ii=1,l
               j2     = j1 + j
               j4     = j2 - 2*ii*inc + inc
c   *** reversed sign inverted sequence repeat (m)
               if (ldbl.eq.1)  then
                 vi(1) = d(j4)
                 vi(2) = d(j4+1)
               else
                 vv     = d(j4)
               end if
               if (k(i).eq.lm)     vv = -vv
               if (ldbl.eq.1) then
                 d(j2)   = vi(1)
                 d(j2+1) = vi(2)
               else
                 d(j2)  = vv
               end if
  480          j      = j+inc
  490       continue
         go to 600
c   *** logarithmic interpolation (l)
  500    l      = in(i) + 1
         vv     = v(i)
         iii    = 2
         go to 600
  510    vu     = v(i)/vv
         del    = dexp(dlog(vu)/dble(l))
         j2     = j1 + j
         if (ldbl.eq.1) then
           d(j2)   = vi(1)
           d(j2+1) = vi(2)
         else
           d(j2)   = vv
         end if
         j      = j+inc
         do 520 ii=2,l
            j2     = j1 + j
            vv     = del*vv
            if (ldbl.eq.1) then
              d(j2)   = vi(1)
              d(j2+1) = vi(2)
            else
              d(j2)  = vv
            end if
  520       j      = j+inc
         if (iii.ne.0) go to 140
         go to 600
c   *** fixed 10 log interp.  (k)
  530    l      = in(i) + 1
         vv     = v(i)
         iii    = 3
         go to 600
  540    dif    = (v(i)-vv)/9.0d0
         del    = dexp(2.302585092994045d0/dble(l))
         j2     = j1 + j
         if (ldbl.eq.1) then
           d(j2)   = vi(1)
           d(j2+1) = vi(2)
         else
           d(j2)   = vv
         end if
         vv     = vv-dif
         j      = j+inc
         do 550 ii=2,l
            j2     = j1 + j
            dif    = dif*del
            delt   = vv + dif
            if (ldbl.eq.1) then
              d(j2)   = dv(1)
              d(j2+1) = dv(2)
            else
              d(j2)   = delt
            end if
  550       j      = j+inc
         if (iii.ne.0) go to 140
         go to 600
c   *** interpolate (i)
  560    l      = in(i) + 1
         vv     = v(i)
         iii    = 1
         go to 600
  570    del    = (v(i)-vv)/float(l)
         do 580 ii=1,l
            j2     = j1 + j
            delt   = del*float(ii-1) + vv
            if (ldbl.eq.1) then
              d(j2)   = dv(1)
              d(j2+1) = dv(2)
            else
              d(j2)   = delt
            end if
            if (kk.ne.0)  ld(j2) = delt + sign(0.5d0,del)
  580       j      = j+inc
         if (iii.ne.0) go to 140
         go to 600
c   *** read alphanumeric
  590    ll     = in(i)
         ll3    = ll+ll2-1
         j2     = 1
         j      = 0
         ncount = 0
         mode   = 0
         if (v(i).eq.0) mode   = 1
         lstart = v(mode+i)
         lstop  = v(mode+i+1)
         lstar  = lld(ll3) + lstart - 1
         lsto   = lld(ll3) + lstop - 1
         read(n5,10900) (d(ii),ii=lstar,lsto)
         write(n6,10300)ll,lstart,lstop
         i      = mode+i+1
  600    continue
         i      = i+1
         if ( i.le.nf ) go to 130
      go to 100
c   *** format statements
10000 format(6(i2,a1,f9.0),t5,2a4,5(4x,2a4))
10100 format(' ',t53,2a4,i7,a1)
10200 format(t90,i5,a1,' array',i7,' entries read at',i4,a1)
10300 format('0',i5,'/ array read from',i5,' to',i5)
10400 format(1h0,i5,a1,' array',i7,' entries read')
10500 format('0****** error ',i7,' entries required in ',i3,a1,' array'/
     *  '0 data edit continues' )
10600 format(1h0,i5,a1)
10700 format('0****** fill option ignored in ',i2,a1,' array' )
10800 format('0****** warning  address',i7,' is beyond limits of',i3,a1,
     *  ' array')
10900 format(18a4)
      end
C
CTAT/MOD(97/Oct) TO COMPATIBLE ON USING WITH C PROGRAM        ----(FROM)
C     function fread ( ibr, iret )
      function ftread ( ibr, iret )
CTAT/MOD(97/Oct) TO COMPATIBLE ON USING WITH C PROGRAM        ----( TO )
C
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     94/11/28
c     time the module was last permanently updated:     15:30:26
c     programmer name:                                  l.m.petrie
c     module name:                                      s7fread
c     current archiving level number:                   00008
c     current number of permanent updates:              00008
c     date of last access by librarian:                 94/11/28
c     dataset name:  x4s.scale4.master
c
c#######################################################################
      double precision dread
c
c     fread returns a single precision floating point number
c
C
CTAT/MOD(97/Oct) TO COMPATIBLE ON USING WITH C PROGRAM        ----(FROM)
      ftread  = dread ( ibr, iret )
CTAT/MOD(97/Oct) TO COMPATIBLE ON USING WITH C PROGRAM        ----( TO )
C
      return
      end
      function iread ( ibr, iret )
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     94/11/28
c     time the module was last permanently updated:     15:31:16
c     programmer name:                                  l.m.petrie
c     module name:                                      s7iread
c     current archiving level number:                   00008
c     current number of permanent updates:              00008
c     date of last access by librarian:                 94/11/28
c     dataset name:  x4s.scale4.master
c
c#######################################################################
      double precision dread
c
c     iread returns a full word integer
c
      iread  = dread ( ibr, iret )
      return
      end
      logical function lread ( iz, ir )
c
c     lread returns the value true if the next character on the
c     card is a number, otherwise it returns false.
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     94/11/28
c     time the module was last permanently updated:     15:24:41
c     programmer name:                                  l.m.petrie
c     module name:                                      s7lread
c     current archiving level number:                   00010
c     current number of permanent updates:              00010
c     date of last access by librarian:                 94/11/28
c     dataset name:  x4s.scale4.master
c
c#######################################################################
      character*1 test
c
c
      common /unit/ inpt, outpt, idum(19)
      common /qrdbuf/ unit, tens, dnum, ichr, lchr, irpt,
     *  wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
      common /yrdbuf/ char, card
      save /unit/, /qrdbuf/, /yrdbuf/
c
      integer inpt, outpt
      double precision dnum, unit, tens, ten
      character*1 card, char, blank, e, d, comma
      logical wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
c
      dimension card(252), char(28)
      dimension unit(10), tens(8)
c
      if ( irpt.gt.0 ) then
c       if irpt > 0 then there is a repeated number to be returned
        lread = .true.
        return
      end if
C
CTAT/MOD(97/Oct) TO AVOID SUDDEN DEATH ---------------------------(FROM)
      if ( ir.eq.2 ) then
           call enfile
           return
      endif
CTAT/MOD(97/Oct) TO AVOID SUDDEN DEATH ---------------------------( TO )
C
      ii     = ichr
  100 continue
      if ( ii.ge.lchr ) then
        call y0read ( card, lchr, ir )
        ichr   = 0
        if ( ir.ne.0 ) then
          lread  = .false.
          return
        end if
        ii     = ichr
      end if
      ii    = ii + 1
      test  = card(ii)
      if ( test.eq.char(17) ) go to 100 !FURUTA
      lread = test.ge.char(01) .and. test.le.char(10)
      lread = lread .or.
     *        test.eq.char(23) .or.
     *        test.eq.char(24) .or.
     *        test.eq.char(26)
      return
      end

      subroutine scanof
      common /unit/ inpt, outpt, idum(19)
      common /qrdbuf/ unit, tens, dnum, ichr, lchr, irpt,
     *  wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
      common /yrdbuf/ char, card
      save /unit/, /qrdbuf/, /yrdbuf/
c
      integer inpt, outpt
      double precision dnum, unit, tens, ten
      character*1 card, char, blank, e, d, comma
      logical wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
c
      dimension card(252), char(28)
      dimension unit(10), tens(8)
c
      lscan = .false.
      return
      end


      subroutine scanon
c
c  scanon sets the reading mode to scan ahead for end
c  scanof turns off the scan ahead reading mode
c
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     94/11/28
c     time the module was last permanently updated:     15:25:49
c     programmer name:                                  l.m.petrie
c     module name:                                      s7scanon
c     current archiving level number:                   00006
c     current number of permanent updates:              00006
c     date of last access by librarian:                 94/11/28
c     dataset name:  x4s.scale4.master
c
c#######################################################################
c
      common /unit/ inpt, outpt, idum(19)
      common /qrdbuf/ unit, tens, dnum, ichr, lchr, irpt,
     *  wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
      common /yrdbuf/ char, card
      save /unit/, /qrdbuf/, /yrdbuf/
c
      integer inpt, outpt
      double precision dnum, unit, tens, ten
      character*1 card, char, blank, e, d, comma
      logical wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
c
      dimension card(252), char(28)
      dimension unit(10), tens(8)
c
      lscan=.true.
      return
c
c     entry scanof
c
      entry allowc
c
      lallc = .true.
      return
c
      entry resetc
c
      lallc = .false.
      return
c
      entry rstptr(i)
c
      ichr = i
      return
c
      entry getptr(i)
      i = ichr
      return
c
      entry setbin
c
      lbin = .true.
      return
c
      entry resetb
c
      lbin = .false.
      return
c
      entry ionums ( nin, nou, min, mou )
c
      min = inpt
      mou = outpt
      if (nin.gt.0) inpt = nin
      if (nou.gt.0) outpt = nou
      return
c
      entry rcrdln ( len,lnsv )
c
      lnsv = lchr
      lchr = len
      return
      end
      subroutine y0read ( cbuf, lbuf, iret )
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     94/11/28
c     time the module was last permanently updated:     15:26:47
c     programmer name:                                  l.m.petrie
c     module name:                                      s7y0read
c     current archiving level number:                   00013
c     current number of permanent updates:              00013
c     date of last access by librarian:                 94/11/28
c     dataset name:  x4s.scale4.master
c
c#######################################################################
c
c     specification statements
c
c
      common /unit/ inpt, outpt, idum(19)
      common /qrdbuf/ unit, tens, dnum, ichr, lchr, irpt,
     *  wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
      common /yrdbuf/ char, card
      save /unit/, /qrdbuf/, /yrdbuf/
c
      integer inpt, outpt
      double precision dnum, unit, tens, ten
      character*1 card, char, blank, e, d, comma
      logical wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
c
      dimension card(252), char(28)
      dimension unit(10), tens(8)
c
      character*1 cbuf
      dimension cbuf(lbuf)
      external qrdbfr
c
c     end of specification statements
c
      crdflg = .true.
  100 continue
      if (lbin) then
        read(inpt,iostat=ios) cbuf
      else
        read(inpt,10000,iostat=ios) cbuf
      end if
      if (ios.ne.0) then
        iret = 2
        if ( .not.lscan ) call enfile
        return
      end if
      if (cbuf(1).ne.'''') then
        call y0trns ( cbuf,lbuf )
        return
      end if
      write(outpt,10100) cbuf
      go to 100
c
c     format statements
c
10000 format(80a1)
10100 format('0',80a1)
      end
      subroutine y0trns ( card, lchr )
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     91/06/26
c     time the module was last permanently updated:     15:40:01
c     programmer name:                                  tza
c     module name:                                      s7y0trns
c     current archiving level number:                   00003
c     current number of permanent updates:              00003
c     date of last modification for scale-pc:           95/07/14
c
c#######################################################################
c
c     y0trns translates lowercase letters to uppercase letters          pc4.3
c     for an ascii character set on pc                                  pc4.3
c
      character*1 card(lchr)
c
      do 100 i=1,lchr
        ii     = ichar(card(i))
c       uppercase letters --> lowercase letters
        if ( ii.ge.65 .and. ii.le.90 ) ii     = ii + 32
        card(i) = char(ii)
  100   continue
      return
      end
      subroutine yread ( dd, ld, lim, it, lerr, iret)
c
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     95/05/09
c     time the module was last permanently updated:     14:04:48
c     programmer name:                                  l.m.petrie
c     module name:                                      s7yread
c     current archiving level number:                   00012
c     current number of permanent updates:              00012
c     date of last access by librarian:                 95/05/09
c     dataset name:  x4s.scale4.master
c
c#######################################################################
c
      common /unit/ inpt, outpt, idum(19)
      common /qrdbuf/ unit, tens, dnum, ichr, lchr, irpt,
     *  wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
      common /yrdbuf/ char, card
      save /unit/, /qrdbuf/, /yrdbuf/
c
      integer inpt, outpt
      double precision dnum, unit, tens, ten
      character*1 card, char, blank, e, d, comma
      logical wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
c
      dimension card(252), char(28)
      dimension unit(10), tens(8)
c
      double precision vnum, vstor, del, vl, dread
      character*1 c,bb,fl,ad,tt,sk,hi,qs,sm,gl,bk,cm
      character*4 aread
      logical lsave, lerr, lsall
      dimension dd(lim), ld(lim)
      dimension vn(2), vs(2)
      equivalence (vnum,vn(1)),(vlow,vn(2)),(vstor,vs(1)),(vslow,vs(2))
c     data statements
c
      data tt/'t'/, sk/'s'/, hi/'i'/, qs/'q'/, sm/'n'/, gl/'l'/, bk/'b'/
c
c     end of specification statements
c
c
      lerr   = .false.
      lsave  = lscan
      lsall  = lallc
      lscan  = .true.
      lallc  = .true.
      m1     = -1
      iret   = 0
      iz     = 0
      iard   = 0
      mult   = max0(1,it)
      llim   = mult*lim
      i      = 1
  100 continue
      if ( iret.ne.0 )  go to 270
      vnum   = dread(iz,iret)
      c      = aread(m1,iard)
  110 continue
      if (c.eq.tt)      go to 270
      if (c.eq.char(17) .or. c.eq.char(18)) then !FURUTA
c     no modification
        if (i.gt.llim) then
          lerr   = .true.
          write(outpt,10000)
        else
          dd(i)  = vnum
          if (it.eq.0) ld(i)  = vnum
          i      = i+1
          if (it.eq.2) then
            dd(i)  = vlow
            i      = i+1
          end if
        end if
      else
        in     = vnum
        if (c.eq.sk) then
c     skip
          i      = i+in*mult
          if (i.le.0 .or. i.gt.llim) then
            write(outpt,10100) in,i
            lerr   = .true.
            i      = max0(1,min0(i,llim-(mult-1)))
          end if
        else
          vnum   = dread(iz,iret)
          if (c.eq.char(11)) then !FURUTA
c     address modification
            nn     = vnum
            i      = nn*mult
            if (i.le.0 .or. i.gt.llim) then
              write(outpt,10200) nn
              lerr   = .true.
              i      = max0(1,min0(i,llim-(mult-1)))
            end if
          else if (c.eq.char(16)) then !FURUTA
c     fill array
  150       continue
            if (i.gt.llim) go to 100
            dd(i)  = vnum
            if (it.eq.0) ld(i)  = vnum
            i      = i+1
            if (it.eq.2) then
              dd(i)  = vlow
              i      = i+1
            end if
            go to 150
          else if (c.eq.qs .or. c.eq.sm .or. c.eq.bk) then
c     sequence repeat
            il     = vnum*mult
            in     = max0(in,1)
            ib     = 0
c     b is a reverse repeat of il entries with a backspace of in
            if (c.eq.bk) then
              ib     = in*mult
              in     = 1
            endif
            if (i-il-ib.le.0) then
              write(outpt,10300) c,il,i,ib
              lerr   = .true.
            else if (i-1+il*in.gt.llim) then
              write(outpt,10400) c,il,i,ib,in
              lerr   = .true.
            else
              do 180 lsq=1,in
                do 170 ii=1,il
c     q is a repeat of il entries
                  if (c.eq.qs) then
                    j      = i - il
c     n is a reverse repeat of il entries
                  else
                    j      = i - (2*ii-1) - ib
                  end if
                  dd(i)  = dd(j)
                  i      = i+1
                  if (it.eq.2) then
                    dd(i)  = dd(j+1)
                    i      = i+1
                  end if
  170             continue
  180           continue
            end if
          else if (c.eq.gl) then
c     logarithmic interpolation
            l      = in+1
            vl     = l
            vstor  = vnum
            vnum   = dread(iz,iret)
            c      = aread(m1,iard)
            if (c.ne.bb) then
              in     = vnum
              vnum   = dread(iz,iret)
            end if
            if (i+l*mult.gt.llim) then
              write(outpt,10500) gl,in
              lerr   = .true.
            else
              del    = exp(log(vnum/vstor)/vl)
              do 260 ii=1,l
                dd(i)  = vstor
                if (it.eq.0) ld(i)  = vstor
                i      = i+1
                if (it.eq.2) then
                  dd(i)  = vslow
                  i      = i+1
                end if
                vstor  = del*vstor
  260           continue
            end if
            go to 110
          else if (c.eq.hi) then
c     linear interpolation
            l      = in+1
            vl     = l
            vstor  = vnum
            vnum   = dread(iz,iret)
            c      = aread(m1,iard)
            if (c.ne.bb) then
              in     = vnum
              vnum   = dread(iz,iret)
            end if
            if (i+l*mult.gt.llim) then
              write(outpt,10500) hi,in
              lerr   = .true.
            else
              del    = (vnum-vstor)/vl
              if (it.eq.0) then
                 if (anint(del*vl).eq.(vnum-vstor)) del=1.00001*del
              end if
              if (abs(anint(del)-del).lt.0.001) del = anint(del)
              do 240 ii=1,l
                dd(i)  = vstor
                if (it.eq.0) ld(i)  = vstor
                i      = i+1
                if (it.eq.2) then
                  dd(i)  = vslow
                  i      = i+1
                end if
                vstor  = vstor + del
  240           continue
            end if
            go to 110
          end if
        end if
      end if
      if (c.ne.tt) go to 100
  270 continue
      lscan  = lsave
      lallc  = lsall
      return
10000 format('0***** error - attempt to store past the end of the ',
     *       'array in yread *****')
10100 format('0***** error - a skip of ',i5,' increments the array ',
     *       'index to ',i5,' which is outside the array in yread ',
     *       '*****')
10200 format('0***** error - an address of ',i5,' is outside the ',
     *       'array in yread *****')
10300 format('0***** error - a ',a,' sequence repeat of ',i5,
     *       ' items starting from an array index of ',i5,
     *       ' with a backspace of ',i5/14x,
     *       ' starts before the beginning of the array in yread *****')
10400 format('0***** error - a ',a,' sequence repeat of ',i5,
     *       ' items starting from an array index of ',i5,
     *       ' with a backspace of ',i5/14x,' repeated ',i5,' times',
     *       ' will store past the end of the array in yread *****')
10500 format('0***** error - interpolation type ',a,
     *       ' specifying ',i5,' entries',
     *       ' will store past the end of the array in yread *****')
      end


      block data qrdbfr
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     94/11/28
c     time the module was last permanently updated:     14:24:34
c     programmer name:                                  l.m.petrie
c     module name:                                      s7qrdbfr
c     current archiving level number:                   00010
c     current number of permanent updates:              00010
c     date of last access by librarian:                 94/11/28
c     dataset name:  x4s.scale4.master
c
c#######################################################################
c
      common /unit/ inpt, outpt, idum(19)
      common /qrdbuf/ unit, tens, dnum, ichr, lchr, irpt,
     *  wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
      common /yrdbuf/ char, card
      save /unit/, /qrdbuf/, /yrdbuf/
c
      integer inpt, outpt
      double precision dnum, unit, tens
      character*1 card, char !FURUTA, blank, e, d, comma
      logical wrdflg, crdflg, lscan, lallc, lneg, lbin, erread, illchr
c
      dimension card(252), char(28)
      dimension unit(10), tens(8)
c
      data lscan /.false./, lallc /.false./, lneg  /.false./
      data lbin  /.false./, erread/.false./
      data char/'0','1','2','3','4','5','6','7','8','9','a','b','c',
     *  'd','e','f',' ',',','r','*','$','&','+','-','z','.','o','p' /
      data unit/1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9/
      data tens/1d00,1d10,1d20,1d30,1d40,1d50,1d60,1d70/
      data irpt/0/, ichr/80/, lchr/80/
c
      data inpt/5/, outpt/6/
c
      end


      block data fidcom
c#######################################################################
c
c                audit trail information
c
c     date the module was last permanently updated:     91/07/08
c     time the module was last permanently updated:     10:25:15
c     programmer name:                                  l.m.petrie
c     module name:                                      s7fidcom
c     current archiving level number:                   00001
c     current number of permanent updates:              00001
c     date of last access by librarian:                 93/02/24
c     dataset name:  x4s.scale4.master
c
c#######################################################################
      common /fidasc/ iprtrg
      data iprtrg/-1/
      end

