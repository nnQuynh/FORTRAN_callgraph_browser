!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas_mesh
      module moddas_mesh

! mesh data (sub.getmsh, getmshdc)
      double precision, allocatable :: gmsh(:)
      integer :: iaddress_gmsh(3) = 1
      integer :: icurrent_gmsh = 1

! tally : r, x, y, z
      double precision, allocatable :: das_itrrg(:)
      double precision, allocatable :: das_itxrg(:)
      double precision, allocatable :: das_ityrg(:)
      double precision, allocatable :: das_itzrg(:)

! tally : energy, angle, time
      double precision, allocatable :: das_iterg(:)
      double precision, allocatable :: das_iterg2(:)
      double precision, allocatable :: das_itarg(:)
      double precision, allocatable :: das_ittrg(:)

! weight window : x, y, z
      double precision, allocatable :: das_iwxrg(:)
      double precision, allocatable :: das_iwyrg(:)
      double precision, allocatable :: das_iwzrg(:)

! weight window bias : x, y, z
      double precision, allocatable :: das_iwbxrg(:)
      double precision, allocatable :: das_iwbyrg(:)
      double precision, allocatable :: das_iwbzrg(:)

! source : x, y, z
      double precision, allocatable :: das_istxx(:)
      double precision, allocatable :: das_istyy(:)
      double precision, allocatable :: das_istzz(:)

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_mesh_initialize()

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)

      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
      common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)
      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)
      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)

      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)

      common /wbxyz/ iwbxty(6), iwbxnm(6), iwbxrg(6),
     &               rwbxmi(6), rwbxma(6), rwbxdl(6),
     &               iwbyty(6), iwbynm(6), iwbyrg(6),
     &               rwbymi(6), rwbyma(6), rwbydl(6),
     &               iwbzty(6), iwbznm(6), iwbzrg(6),
     &               rwbzmi(6), rwbzma(6), rwbzdl(6)
      
      common /isormd/ istdd(isrc),istcc(isrc),
     &                isxtp(isrc),isinx(isrc),istxx(isrc),
     &                isytp(isrc),isiny(isrc),istyy(isrc),
     &                isztp(isrc),isinz(isrc),istzz(isrc),
     &                sxmin(isrc),sxmax(isrc),sxdel(isrc),
     &                symin(isrc),symax(isrc),sydel(isrc),
     &                szmin(isrc),szmax(isrc),szdel(isrc)

      itrrg(:) = 1
      itxrg(:) = 1
      ityrg(:) = 1
      itzrg(:) = 1

      iterg(:) = 1
      iterg2(:) = 1
      itarg(:) = 1
      ittrg(:) = 1

      iwxrg(:) = 1
      iwyrg(:) = 1
      iwzrg(:) = 1

      iwbxrg(:) = 1
      iwbyrg(:) = 1
      iwbzrg(:) = 1

      istxx(:) = 1
      istyy(:) = 1
      istzz(:) = 1

      allocate(das_itrrg(1), das_itxrg(1), das_ityrg(1), das_itzrg(1))

      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module
