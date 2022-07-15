module v_nabla_factory_mod

use abstract_v_nabla_mod, only : v_nabla_operator_t
use v_nabla_mod,          only : v_nabla_c2_operator_t, v_nabla_c4_operator_t,   &
                                 v_nabla_up1_operator_t, v_nabla_up3_operator_t, &
                                 v_nabla_up4_operator_t
use parcomm_mod,          only : parcomm_global

implicit none

contains

subroutine create_v_nabla_hor_operator(v_nabla_operator,halo_width,v_nabla_operator_name)
    character(len=*),          intent(in) :: v_nabla_operator_name
    class(v_nabla_operator_t), allocatable, &
                               intent(out) :: v_nabla_operator
    integer(kind=4),           intent(out) :: halo_width

    select case(v_nabla_operator_name)
    case("c2")
        v_nabla_operator = v_nabla_c2_operator_t()
        halo_width = 1
    case("c4")
        v_nabla_operator = v_nabla_c4_operator_t()
        halo_width = 2
    case("up1")
        v_nabla_operator = v_nabla_up1_operator_t()
        halo_width = 1
    case("up3")
        v_nabla_operator = v_nabla_up3_operator_t()
        halo_width = 2
    case("up4")
        v_nabla_operator = v_nabla_up4_operator_t()
        halo_width = 3
    case default
        call parcomm_global%abort("create_v_nabla_hor_operator, "//&
                                  "unknown horizontal advection operator:"//&
                                  v_nabla_operator_name)
    end select
end subroutine create_v_nabla_hor_operator

end module v_nabla_factory_mod
