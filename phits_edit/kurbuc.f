c [proton track]
c===========================
      subroutine ionts_setup
c===========================

      return
      end

c============================
      subroutine pts_flt(fpl)
c============================

      return
      end

c=====================================================
      subroutine ionts_trn(mark,markp,fpl,nbeta,itmak)
c=====================================================

      return
      end

c========================
      subroutine pts_reac
c========================

      return
      end

c=========================
      subroutine pts_input
c=========================

!      write(*,*)'[Track Structure] is not available'
!      stop

      return
      end



c [carbon ion track]
c==============================
      subroutine cts_flt(fpl)
c==============================

      return
      end


c==========================================
      subroutine cts_reac
c==========================================

      return
      end

c=========================
      subroutine cts_input
c=========================

!      write(*,*)'[Track Structure] is not available'
!      stop

      return
      end

c=========================
      subroutine KURBUCcheck(KBCflg)
c=========================
      include 'err.inc'
! T.Sato 2020/03/18
      common /tsmaxcom/tsmax
      write(ErrCha,'(''Warning :: tsmax > 1.0d-3, but this exe file '',
     &'' cannot run KURBUC. ITSART is used instead.'')')
      ErrID = 'L:76/R:KURBUCcheck/F:kurbuc.f' !W00_000_000
      call ErrWrite(ErrID,ErrCha)
      KBCflg = 1
      return

      end
