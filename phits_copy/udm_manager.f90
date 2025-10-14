module udm_Manager

!=======================================================================
! [udm_int]
use udm_int_sample,   caller_udm_int_sample   => caller
use udm_int_ee2mumu,  caller_udm_int_ee2mumu  => caller
use udm_int_kill,     caller_udm_int_kill     => caller
use udm_int_neutrino, caller_udm_int_neutrino => caller
! [udm_part]
use udm_part_sample,   caller_udm_part_sample  => caller
!=======================================================================

implicit none
contains




subroutine user_defined_interaction(action,index)
integer action,index
!=======================================================================
call caller_udm_int_sample(action,index)
call caller_udm_int_ee2mumu(action,index)
call caller_udm_int_kill(action,index)
call caller_udm_int_neutrino(action,index)
!=======================================================================
end subroutine user_defined_interaction




subroutine user_defined_particle(action)
use udm_Parameter, only: udm_part_num
integer action,index
do index=1,udm_part_num
!=======================================================================
call caller_udm_part_sample(action,index)
!=======================================================================
enddo
end subroutine user_defined_particle

end module udm_Manager
