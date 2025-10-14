C --------------------------------------------------------
C  set:%key%[value] parameter
C --------------------------------------------------------
      module CHARVARMOD
      implicit none

      include 'param.inc'
      include 'err.inc'

      type chapara
       character(len=100) :: key
       character(len=100) :: value
      end type
      type(chapara) :: cvariable(maxchapara)

      character(len=200) :: filnm,ifilnm
      integer :: irwt,irwtaddflag,irwtflag
      integer :: icpcount
      integer :: iendflag ,iinflflag

      integer,parameter :: inflmax = 9
      type inflnlog
      integer :: n1,n2
      end type
      type(inflnlog) :: infln(inflmax)
      character(len=200) :: nowfile,logfile(inflmax),logfile0

      integer :: iflog,inplog,infllog(inflmax)

      integer :: iflogr(maxchapara)
     1          ,inplogr(maxchapara)
     2          ,infllogr(maxchapara)
     3          ,cc
      character(len=200) :: filenames(maxchapara)

      contains
C --------------------------------------------------------

      subroutine set_Charpara(jsi,ifileflag,fname,ierr)
      implicit none
      integer,intent(in) :: jsi,ifileflag
      character(len=200),intent(in) :: fname
      integer,intent(out) :: ierr
      integer :: npe,me
      common /mpi00/ npe, me
      integer :: ios
      character(len=200) :: chin,chlw,chadjl
      integer :: i,j
      integer :: inum,inumtmp,iunit,inumlog,inumjsi

C - Initial --------------------
      ierr=0
      filnm  = "phits_rwt.inp"
      ifilnm = "phits.inp"
      logfile = " "

      icpcount = 0
      irwt     = 1
      irwtaddflag = 0
      iendflag = 0
      iinflflag = 0

      infln%n1 = 0
      infln%n2 = 0

      iflog   = 0
      inplog  = 0
      infllog = 0

      iflogr   = 1
      inplogr  = 1
      infllogr = 0
      filenames = ifilnm
      cc       = 1

      inumtmp = 901
      inumlog = 31
      inumjsi = 31

      cvariable%key   = ' '
      cvariable%value = ' '

C ------------------------------

C Header Parameter Read
      call set_CharparaFlag(jsi,ifileflag,fname,irwtflag)
      filenames(1) = ifilnm
      if(irwt==0) return
      rewind(jsi)
      inquire(jsi, number=inum)
      inumjsi = inum
      inumlog = inum

      if(me.eq.0)then ! Only for master

C rwt file write start
       open(901, file=filnm,status='unknown')

       if(irwt.gt.0) then
C Header Parameter Replace
        do
         read(jsi,'(a)', iostat = ios) chin
         if( ios == -1 ) exit

         iflog  = iflog  + 1
         inplog = inplog + 1
        
         do i = 1, len_trim(chin)
          if (chin(i:i) >= 'A' .and. chin(i:i) <= 'Z') then
           chin(i:i) = char(ichar(chin(i:i)) + 32)
          endif
         enddo

         call chcomp(chin,1,200,i)
         chadjl = adjustl(chin)

         if(index(chadjl(1:4),'file')/=0) then
          chin = "#"//trim(chin)
         endif

         if(index(chadjl(1:4),'$rwt')/=0) then
          chin = "$rwt=0 "
         endif

C First Section Search
c$$$         if(chadjl(1:1) == "[") then
         if(chadjl(1:1) == "[" .or.
     &        index(chadjl(1:4),'infl')/=0 .or.
     &        index(chadjl(1:4),'set:')/=0 ) then
          if( irwtflag == 0) then
           write(901,fmt='(a)') '$rwt=0'
           irwtaddflag = 1
          endif
c$$$          write(901,fmt='(a)') chin(1:len_trim(chin))
          backspace(inum)
          iflog  = iflog  - 1
          inplog = inplog - 1
          exit
         endif
         write(901,fmt='(a)') chin(1:len_trim(chin))
        enddo
       endif
 
       if(irwtaddflag == 1) then
         iflog = iflog  + 1
         inplog = inplog + 1
       endif   

C =================
C read main input  
C =================
       inquire(unit=inum,name=nowfile)
       logfile0 = nowfile
       
       do
        iinflflag = 0
        read(inum,'(a)', iostat = ios) chin

        inquire(unit=inum,name=nowfile) 
        if(nowfile/=logfile0) then    
         cc = cc + 1
         iflogr(cc)    = iflog + 1
         inplogr(cc)   = inplog
         if(inum == inumjsi) inplogr(cc)   = inplog + 1
         if(inum /= inumjsi) infllogr(cc)  = infllog(inum-901) + 1
         if(inum == inumjsi) filenames(cc) = ifilnm
         if(inum /= inumjsi) filenames(cc) = logfile(inum-901)
        endif    

        if( ios == -1 ) then
         if( inum == 902 ) then

          if(infllog(inum-901) == 0 ) then
           write(6,103) trim(filenames(cc)) 
  103 format(" ** Warning:"1x,a,1x"is empty")
            cc = cc - 1
            iflog = iflog + 1
            write(901,fmt='(a)') "       "
          endif  

          if( infln(inum-901)%n2 > infllog(inum-901) ) then
           write(6,101,advance='no')
           write(6,102) infln(inum-901)%n2
     1          ,infllog(inum-901)
     2          ,trim(filenames(cc))
 101       format(
     1     " ** Warning: the upper limit line specified in infl, ")
 102       format(
     1     i0,", exceeds, ",i0,", the number of lines contained in ",a)
          endif

          infllog(inum-901) = 0  
          close(inum)
          inumtmp = 901
          inum = inumjsi
          cycle  

         elseif( inum >= 903 ) then

          if(infllog(inum-901) == 0 ) then
            infllogr(cc) = 0
            cc = cc - 1
          endif  

          if( infln(inum-901)%n2 > infllog(inum-901) ) then
          write(6,101,advance='no')  
          write(6,102) trim(filenames(cc))
     1                ,infllog(inum-901)
     2                ,infln(inum-901)%n2    
          endif

          close(inum)   
          infllog(inum-901) = 0
          infllog(inum-902) = infllog(inum-902) + 1
          inum = inum - 1      
          inumtmp = inumtmp - 1      
          cycle

         else
          exit

         endif 
        endif 

C -------------------------------
        chlw = chin
        do i = 1, len_trim(chin)
         if (chin(i:i) >= 'A' .and. chin(i:i) <= 'Z') then
          chlw(i:i) = char(ichar(chin(i:i)) + 32)
         endif
        enddo
        call chcomp(chlw,1,200,i)
        chadjl = adjustl(chlw)

        if((index(chadjl(1:5),'[end]')/=0) .and.
     1     (index(chin,'off')==0)) then
         iendflag = 1   
        endif

        if(iendflag == 1 ) goto 999 

         if(index(chadjl(1:4),'infl')/=0) then
          call set_Replace_Call(chin,chadjl,chlw,ierr)
          if(ierr.ne.0)return
  
          call Search_infl(chin,inumtmp,ierr)
          if(ierr.ne.0)return
  
          if(iinflflag == 1 ) goto 991
          if(inumtmp==902) inplog = inplog + 1
          inum = inumtmp
          cycle 
         endif          
         call set_Function(chin,inum,ierr)
         if(ierr.ne.0)return

  991 continue
C -----------------------------
        iflog  = iflog  + 1
        if(inum == inumjsi) then
         inplog = inplog + 1
C            write(61,8) inum,iflog,inplog,ifilnm,chin  ! for check
        else
         infllog(inum-901) = infllog(inum-901) + 1  
C           inquire(unit=inum,name=nowfile) 
C           write(61,9) inum,iflog,inplog,infllog(inum-901),  ! for check
C      1   trim(nowfile),chin
        endif 
    8 format(i3,i3,i3,3x,1x,a10,1x,a20)    
    9 format(i3,i3,i3,i3,1x,a,1x,a20)

C -----------------------------
C infl skip line

        if(inum >= 902) then
         if(( infln(inum-901)%n1 /= 0) .and. 
     1      ( infln(inum-901)%n1 > infllog(inum-901) )) then
          chin = "$! "//trim(chin)
         endif  

         if(( infln(inum-901)%n2 /= 0) .and. 
     1      ( infln(inum-901)%n2 < infllog(inum-901) )) then
          chin = "$! "//trim(chin)
         endif
        endif    
C -----------------------------  
  999 continue
        write(901,fmt='(a)') chin(1:len_trim(chin)) 
        inumlog = inum  
        inquire(unit=inum,name=logfile0) 
       enddo
       close(901)

      endif ! Only for master

      close(jsi)
      open(unit=jsi,file=trim(filnm),status='unknown')

      return
      end subroutine set_Charpara

C -----------------------------------------------------

      subroutine set_CharparaFlag(jsi,ifileflag,fname,irwtflag)
      implicit none
      integer,intent(in) :: jsi,ifileflag
      character(len=200),intent(in) :: fname
      integer,intent(out) :: irwtflag
      character(len=200) :: chin,chlw
      character(len=200) :: chadjl
      integer :: i,ios
      integer :: jfileflag

      irwtflag = 0
      jfileflag = ifileflag
      if(jfileflag.eq.1)call Get_file0(fname)
      do
       read(jsi,'(a)', iostat = ios) chin
       if( ios == -1 ) exit

       chlw = chin
       do i = 1, len_trim(chin)
        if (chin(i:i) >= 'A' .and. chin(i:i) <= 'Z') then
         chlw(i:i) = char(ichar(chin(i:i)) + 32)
        endif
       enddo

       call chcomp(chlw,1,200,i)
       chadjl = adjustl(chlw)
       if(chadjl(1:1) == "[") exit

       if(jfileflag == 0) then
        if(index(chadjl(1:4),'file')/=0) then
         jfileflag = 1
         call Get_file(chadjl,jsi)
        endif
       endif

       if(index(chadjl(1:5),'$rwt=')/=0) then
        irwtflag = 1
        call set_rwtflag(chadjl)
       endif
      enddo

      return
      end subroutine set_CharparaFlag

C -----------------------------------------------------
      subroutine Get_file(chin,jsi)
      implicit none
      character(len=200),intent(in) :: chin
      integer,intent(in) :: jsi
      integer :: j1,j2
      integer :: j
      character(len=200) :: hfile

      j1=0
      j2=0
      j1 = index(chin(1:len_trim(chin)),"=")

      hfile = chin(j1+1:len_trim(chin))
      ifilnm = trim(hfile)

      j = index(hfile,".", back=.true.)
      if( j /= 0 ) then
       filnm = hfile(1:j-1)//"_rwt"//hfile(j:len_trim(hfile))
      else
       filnm = hfile//"_rwt"
      endif

      close(jsi)
      open(unit=jsi, file=trim(hfile),status='unknown')
      read(jsi,*)

      return
      end subroutine Get_file

C -----------------------------------------------------
      subroutine Get_file0(chin)
      implicit none
      character(len=200),intent(in) :: chin
      integer :: j
      character(len=200) :: hfile

      hfile = chin(1:len_trim(chin))
      ifilnm = trim(hfile)

      j = index(hfile,".", back=.true.)
      if( j /= 0 ) then
       filnm = hfile(1:j-1)//"_rwt"//hfile(j:len_trim(hfile))
      else
       filnm = hfile//"_rwt"
      endif

      return
      end subroutine Get_file0

C -----------------------------------------------------
      subroutine set_rwtflag(chin)
      implicit none
      character(len=200),intent(in) :: chin

      read(chin(6:len_trim(chin)),*) irwt

      return
      end subroutine set_rwtflag

C -----------------------------------------------------

      subroutine set_Function(chin,inumtmp,ierr)
      implicit none
      character(len=200),intent(inout) :: chin
      integer :: inumtmp
      integer,intent(out) :: ierr
      character(len=200) :: chlw,chadjl
      integer :: i

      ierr=0
      chlw = chin
      do i = 1, len_trim(chin)
       if (chin(i:i) >= 'A' .and. chin(i:i) <= 'Z') then
        chlw(i:i) = char(ichar(chin(i:i)) + 32)
       endif
      enddo
      call chcomp(chlw,1,200,i)
      chadjl = adjustl(chlw)

C ---------------------------------------------------------------

      if((chadjl(1:1) /= "#").and.
     1   (chadjl(1:1) /= "!").and.
     2   (chadjl(1:1) /= "$")) then
 
       if((index(chadjl(1:4),'set:')/=0)              
     1    .and.(index(chadjl(1:5),'set:c')==0 )) then
        call set_Search(chin,ierr)
        if(ierr.ne.0)return
        chin = "$"//trim(chin) ! T.Sato 2023/03/30 change from # to $
       endif

       call set_Replace_Call(chin,chadjl,chlw,ierr)
       if(ierr.ne.0)return

      endif 

      return
      end subroutine set_Function

C -----------------------------------------------------

      subroutine set_Search(chin,ierr)
      implicit none
      character(len=200),intent(in) :: chin
      integer,intent(out) :: ierr
      integer :: i,j
      integer :: j1,j2,j3,j4
      character(len=100) :: key
      character(len=100) :: value

      ierr=0
      j1=0
      j2=0
      j3=0
      j4=0
      j1 = index(chin(1:len_trim(chin)),"%") 
      j2 = j1 + index(chin(j1+1:len_trim(chin)),"%") 
      j3 = index(chin(1:len_trim(chin)),"[") 
      j4 = j3 + index(chin(j3+1:len_trim(chin)),"]", back=.true.) 
     
C -- Error ---
      if((j1==0).or.(j2==0).or.(j1==j2)) then
       call ErrorMessage(chin)
       write(ErrCha,
     1      '("error = Set Parameter is not enclosed in %")')
       ErrID = ""
       call ErrWrite(ErrID,ErrCha)
       call EndingTime
       ierr=1
       return
      endif

C -- Error ---
      if((j3==0).or.(j4==0).or.(j3==j4)) then
       call ErrorMessage(chin)
       write(ErrCha,
     1      '("error = Set Parameter is not enclosed in [ ]")')
       ErrID = ""
       call ErrWrite(ErrID,ErrCha)
       call EndingTime
       ierr=1
       return
      endif

      key   = chin(j1+1:j2-1)
      value = chin(j3+1:j4-1)

      do i=1,maxchapara
       if( trim(key) == cvariable(i)%key) then
        cvariable(i)%value = trim(value)
        exit
       endif

       if( trim(key) /= cvariable(i)%key) then
        if( i == maxchapara ) then
        icpcount = icpcount + 1
C Error --
        if( icpcount > 1000 ) then
         write(ErrCha,'("The maximum number of 
     1      chapara definitions is up to 1000")')
          ErrID = ""
          call ErrWrite(ErrID,ErrCha)
          ierr=1
          return
         endif
         cvariable(icpcount)%key = trim(key)
         cvariable(icpcount)%value = trim(value)
        endif
       endif
      enddo

      return
      end subroutine set_Search

C -----------------------------------------------------

      subroutine Search_infl(chin,inumtmp,ierr)
      implicit none
      character(len=200),intent(inout) :: chin
      integer,intent(inout)  :: inumtmp
      integer,intent(out) :: ierr
      
      character(len=200) :: chlw
      integer :: i
      logical :: I_EXIST,OPND
      integer :: j0,j1,j2,j3,j4,j5
      integer :: n1,n2
      character(len=200) :: ifile

      ierr=0
      chlw=chin
      call chcomp(chlw,1,200,i)

      n1 = 0
      n2 = 0

      j0=0
      j1=0
      j2=0
      j3=0
      j4=0
      j5=0
      j1 = index(chlw(1:len_trim(chlw)),"{")
      j2 = j1 + index(chlw(j1+1:len_trim(chlw)),"}")
      j3 = index(chlw(1:len_trim(chlw)),"[")
      j4 = j3 + index(chlw(j3+1:len_trim(chlw)),"]")

C T.Sato 2023/03/24 find comment remark
      j0 = index(chlw(1:len_trim(chlw)),"$")
      if(j0.eq.0) j0 = index(chlw(1:len_trim(chlw)),"#")
      if(j0.ne.0) then ! there is comment remark
       if(j3.gt.j0) j3=0  ! [ is behind comment
       if(j4.gt.j0) j4=0  ! ] is behind comment
      endif

C -- Error ---
      if( ( ((j3/=0).and.(j4/=0)).and.(j3==j4)).or.
     1      ((j3==0).and.(j4/=0))) then
       call ErrorMessage(chin)
       write(ErrCha,
     1      '("error = infl Specify number of lines in [ ]")')
       ErrID = ""
       call ErrWrite(ErrID,ErrCha)
       call EndingTime
       ierr=1
       return
      endif

C -- Error ---
      if((j3/=0).and.(j4/=0)) then
       j5 = index(chlw(j3+1:j4-1),"-")
       if(j5 == 0) then
        call ErrorMessage(chin)
        write(ErrCha,
     1      '("error = infl Specify number of lines in - ")')
        ErrID = ""
        call ErrWrite(ErrID,ErrCha)
        call EndingTime
        ierr=1
        return
       endif

       if(index(chlw(j3:j4),"-")-index(chlw(j3:j4),"[")/=1) then
        read(chlw(j3+1:j3+j5-1),*) n1
       endif

       if(index(chlw(j3:j4),"]")-index(chlw(j3:j4),"-")/=1) then
        read(chlw(j3+j5+1:j4-1),*) n2
       endif
      endif

      ifile = chlw(j1+1:j2-1)

      inquire(file=trim(ifile),exist=I_EXIST, opened=OPND)

      if( I_EXIST .eqv. .FALSE. ) then 
       write(6,*) 
     1'Warning in making rwt file: infl file is not found ->> '
     2,trim(ifile)
       iinflflag = 1
       return
      endif

C       if( I_EXIST .eqv. .FALSE. ) then 
C        call ErrorMessage(chin)
C        write(6,*) 
C      1'Error : infl file is not found ->> '
C      2,trim(ifile)
C        call EndingTime 
C        call close_rwtfile(1)
C        call parastop(999)
C       endif

      if( OPND .eqv. .FALSE. ) then 
       inumtmp = inumtmp + 1
C Error -
        if ( inumtmp >= (901+inflmax) ) then
         call ErrorMessage1(chin)
         write(ErrCha,'("Number of Include File is too deep. 
     ! Max Number is 8.")')
         ErrID = ""
         call ErrWrite(ErrID,ErrCha)
         call EndingTime
         ierr=1
         return
        endif 

       open(unit=inumtmp,file=ifile,status='old')
       logfile(inumtmp-901) = ifile
      else
       call ErrorMessage(chin)
       write(6,*)
     1 "Error = ",trim(ifile)," file is recursively called."
       call EndingTime
       ierr=1
       return
      endif

      infln(inumtmp-901)%n1 = n1
      infln(inumtmp-901)%n2 = n2

      return
      end subroutine Search_infl

C --------------------------------------------------------

      subroutine set_Replace_Call(chin,chadjl,chlw,ierr)
      implicit none
      character(len=200),intent(inout) :: chin,chadjl,chlw
      integer,intent(out) :: ierr
      integer :: iwarning

      ierr=0
      iwarning = 0
       do
        if((index(chadjl(1:3),'set')==0).and.(index(chlw,'%')/=0)) then   
         call set_Replace(chin,iwarning,ierr)
         if(ierr.ne.0)return
         chlw=chin
         if( iwarning == 1 ) exit
        elseif((index(chadjl(1:5),'set:c')/=0)
     1     .and.(index(chlw,'%')/=0)) then
         call set_Replace(chin,iwarning,ierr)
         if(ierr.ne.0)return
         chlw=chin
         if( iwarning == 1 ) exit         
        else
         exit 
        endif 
       enddo      

      return
      end subroutine set_Replace_Call

C --------------------------------------------------------

      subroutine set_Replace(chin,iwarning,ierr)
      implicit none
      character(len=200),intent(inout) :: chin
      integer,intent(out) :: ierr
      character(len=200) :: chtmp
      integer :: i
      integer :: j1,j2
      integer :: iwarning

      ierr=0
      chtmp = chin
      j1=0
      j2=0
      j1 = index(chin(1:len_trim(chin)),"%")
      j2 = index(chin(j1+1:len_trim(chin)),"%")

      if((j1==0).or.(j2==0)) then
       iwarning = 1
C -- Error ---
      elseif(icpcount==0) then
       call ErrorMessage(chin)
       j2 = j1 + index(chin(j1+1:len_trim(chin)),"%")
       write(6,101,advance='no') chin(j1+1:j2-1)
       write(ErrCha,
     1      '(": character variable is not found")')
       ErrID = ""
       call ErrWrite(ErrID,ErrCha)
       call EndingTime
       ierr=1
       return
      else
       j2 = j1 + index(chin(j1+1:len_trim(chin)),"%")
      endif

      if(iwarning /= 1) then
       if(icpcount /= 0) then
        do i=1,icpcount
         if(chin(j1+1:j2-1) == cvariable(i)%key ) then
          chin=" "
          chin=chtmp(1:j1-1)//trim(cvariable(i)%value)
     1  //chtmp(j2+1:len_trim(chtmp))
         exit
         endif

C -- Error ---
         if(i == icpcount) then 
          call ErrorMessage(chin)
          write(6,101,advance='no') chin(j1+1:j2-1)
          write(ErrCha,
     1      '(": character variable is not found")')
          ErrID = ""
          call ErrWrite(ErrID,ErrCha)
          call EndingTime
          ierr=1
          return
         endif
        enddo
       endif
      endif
  101 format("error = %",a,"%")

      return
      end subroutine set_Replace

C --------------------------------------------------------

      subroutine ErrLine_Adjust(ctmp,l_err)
      implicit none
      character(len=*),intent(out) :: ctmp
      integer,intent(inout) :: l_err
      integer :: i,ii

C       do i=1,cc
C        write(6,101) iflogr(i),inplogr(i),infllogr(i),filenames(i)
C       enddo
C   101 format(3(i3,1x),a40) 

      if(irwt.gt.0)then
       ii = 0
       do i=1,maxchapara
        if((i>=2).and.(iflogr(i)==1)) exit
        if( iflogr(i) <= l_err ) then        
         ii = ii + 1     
        endif     
       enddo    
      
       if(irwtaddflag == 0) then
        if( infllogr(ii) == 0 ) then
         l_err = l_err - iflogr(ii) + inplogr(ii)
        else  
         l_err = l_err - iflogr(ii) + infllogr(ii)
        endif
       else
        if( infllogr(ii) == 0 ) then
         l_err = l_err - iflogr(ii) + inplogr(ii) - 1
        else  
         l_err = l_err - iflogr(ii) + infllogr(ii)
        endif
       endif        
       ctmp = trim(filenames(ii))
      endif
    
      return
      end subroutine ErrLine_Adjust

C-----------------------------------------------------------------------
C     Ending Time
C-----------------------------------------------------------------------
      subroutine EndingTime
      implicit none
      integer :: iyer1,imon1,iday1,ihor1,imin1,isec1

            call date_a_time(iyer1,imon1,iday1,
     &                       ihor1,imin1,isec1)

            write(6,'(/79(''-''))')
            write(6,'('' job termination date : '',
     &            i4,''/'',i2.2,''/'',i2.2)') iyer1, imon1, iday1

            write(6,'(''                 time :   '',
     &              i2.2,'':'',i2.2,'':'',i2.2)') ihor1, imin1, isec1

            write(6,*)
            write(6,*) 'END'

      return
      end subroutine EndingTime
C-----------------------------------------------------------------------
C     Error Message
C-----------------------------------------------------------------------
      subroutine ErrorMessage(chin)
      implicit none
      character(len=200),intent(in) :: chin
      character(len=200) :: ctmp
   
      write(6,*)
      write(6,*) "***** Error Message from Input File *****"
      write(6,*)
      ctmp = ""
      call ErrLine_Adjust(ctmp,iflog)
      iflog = iflog + 1 
      write(6,100) trim(ctmp),iflog
   
  100 format(a," : ",i0)

      return
      end subroutine ErrorMessage
C-----------------------------------------------------------------------
C     Error Message for infl deep Error
C-----------------------------------------------------------------------
      subroutine ErrorMessage1(chin)
      implicit none
      character(len=200),intent(in) :: chin
      character(len=200) :: ctmp

      write(6,*)
      write(6,*) "***** Error Message from Input File *****"
      write(6,*)
      ctmp = ""
      call ErrLine_Adjust(ctmp,iflog)
      write(6,100) trim(ctmp)
      write(6,101) trim(chin)

  100 format(a)
  101 format("line : ",a)

      return
      end subroutine ErrorMessage1
*-----------------------------------------------------------------------
*  set char parameter rwt file delete flag
*-----------------------------------------------------------------------
      subroutine close_rwtfile(ierr)
      implicit none
      integer,intent(in) :: ierr
      integer :: npe,me
      common /mpi00/ npe, me
      integer :: num,icc

      call parabcsti(irwt)      ! Broadcast irwt to MPI processes
      if(irwt/=0) then
       if(me.eq.0)then
        if(npe.gt.1)then
         call paraiccp          ! Synchronize MPI processes
        endif
        inquire(file=trim(filnm),number=num)
        if(num.ne.-1) then ! file is still open
         if(irwt.eq.1) close(num,status="delete")
         if(irwt.ge.2) close(num)
        elseif(irwt.eq.1) then  ! make sure to delete *.rwt file
         open(901,file=trim(filnm),status='unknown')
         close(901,status="delete")
        endif
        if(irwt.eq.3.and.ierr.eq.0) then
         write(*,*)'RWT file was successfully generated!!'
         call parafin
         stop
        endif
       elseif(me.gt.0)then
        icc=0
        inquire(file=trim(filnm),number=num)
        if(num.ne.-1) close(num) ! file is still open
        call parasi(icc,1,0)    ! Synchronize MPI processes
        if(irwt.eq.3)then
         call parafin
         stop
        endif
       endif
      endif

      end subroutine close_rwtfile
*-----------------------------------------------------------------------

      end module CHARVARMOD

