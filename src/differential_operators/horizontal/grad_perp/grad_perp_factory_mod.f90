module grad_perp_factory_mod

use domain_mod,              only : domain_t
use abstract_grad_perp_mod,  only : grad_perp_operator_t
use parcomm_mod,             only : parcomm_global

implicit none

contains

subroutine create_grad_perp_operator(grad_perp, grad_perp_operator_name, domain)

    class(grad_perp_operator_t), allocatable, intent(out) :: grad_perp
    character(len=*),                         intent(in)  :: grad_perp_operator_name
    type(domain_t),                           intent(in)  :: domain

    if(grad_perp_operator_name == "grad_perp_c_sbp42" .or. &
       grad_perp_operator_name == "grad_perp_c_sbp21") then
        call create_grad_perp_operator_c_sbp(grad_perp, grad_perp_operator_name, domain)
    else
        call parcomm_global%abort("unknown grad_perp operator name: "// grad_perp_operator_name)
    end if
end subroutine create_grad_perp_operator

subroutine create_grad_perp_operator_c_sbp(grad_perp, grad_perp_operator_name, domain)

    use grad_perp_c_sbp_mod,  only : grad_perp_c_sbp_t
    use sbp_factory_mod,      only : create_sbp_operator
    use exchange_factory_mod, only : create_xy_points_halo_exchange
    use halo_factory_mod,     only : create_halo_procedure

    class(grad_perp_operator_t), allocatable, intent(out) :: grad_perp
    character(len=*),                         intent(in)  :: grad_perp_operator_name
    type(domain_t),                           intent(in)  :: domain

    type(grad_perp_c_sbp_t), allocatable :: grad_perp_sbp
    integer(kind=4) :: halo_width

    allocate(grad_perp_sbp)

    select case(grad_perp_operator_name)
    case("grad_perp_c_sbp21")
        grad_perp_sbp%sbp_diff = create_sbp_operator("D21_staggered_i2c")
        halo_width = 2
    case("grad_perp_c_sbp42")
        grad_perp_sbp%sbp_diff = create_sbp_operator("D42_staggered_i2c")
        halo_width = 4
    case default
        call parcomm_global%abort("unknown sbp grad_perp operator: "//grad_perp_operator_name)
    end select
    grad_perp_sbp%exchange = create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                                          domain%topology, halo_width, 'full')

    call move_alloc(grad_perp_sbp, grad_perp)

end subroutine create_grad_perp_operator_c_sbp

end module grad_perp_factory_mod
