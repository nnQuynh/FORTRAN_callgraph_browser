!
      character(len=256) :: ErrCha
      character(len=256) :: ErrID
      character(len=256) :: MsgID
!
	  common / Eflag / SOURCE_Eflag,tdcheck_CallFlag
	  integer :: SOURCE_Eflag,tdcheck_CallFlag
!
	  real(8) :: tbina_tmp
	  real(8) :: before_time,now_time
!
      integer :: nep_set
!
	  common / SPC / ZSPC_Eflag,ZenkakuSpaceLine,zspc_count
	  integer :: column,zspc_count,ZSPC_Eflag
	  character(len=200) :: ZenkakuSpaceLine(10000)

      common/errcom/ErrCha,ErrID

!  *************  How to output error message  ***********************
!  1. include 'err.inc' if it is not included in the subroutine
!  2. write error message like below:
!    write(ErrCha,'("Parameter a is wrong: a=",f10.2)') a
!    ErrID=''
!    call ErrWrite(ErrID,ErrCha)
!  *******************************************************************
