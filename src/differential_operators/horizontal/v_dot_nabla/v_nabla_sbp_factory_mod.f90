module v_nabla_sbp_factory_mod

use v_nabla_sbp_mod,     only : v_nabla_sbp_operator_t
use sbp_factory_mod,     only : create_sbp_operator


implicit none

contains

function create_v_nabla_sbp_operator(sbp_operator_name) result(v_nabla_sbp_operator)
    character(len=*), intent(in) :: sbp_operator_name
    type(v_nabla_sbp_operator_t) :: v_nabla_sbp_operator

    v_nabla_sbp_operator%sbp_op = create_sbp_operator(sbp_operator_name)
end function

end module v_nabla_sbp_factory_mod
