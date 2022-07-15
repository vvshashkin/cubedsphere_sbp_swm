module adv_z_factory_mod

use abstract_adv_z_mod, only : adv_z_t
use adv_z_mod,          only : adv_z_c2_t, adv_z_c4_t, adv_z_up4_t
use parcomm_mod,        only : parcomm_global

implicit none

contains

subroutine create_adv_z_operator(adv_z_op, adv_z_op_name)
    class(adv_z_t), allocatable, intent(out) :: adv_z_op
    character(len=*), intent(in) :: adv_z_op_name

    select case(adv_z_op_name)
    case("adv_z_c2")
        adv_z_op = adv_z_c2_t()
    case("adv_z_c4")
        adv_z_op = adv_z_c4_t()
    case("adv_z_up4")
        adv_z_op = adv_z_up4_t()
    case default
        call parcomm_global%abort("create_adv_z_operator, unknown adv_z_operator_name:"//&
                                  adv_z_op_name)
    end select
end subroutine create_adv_z_operator

end module adv_z_factory_mod
