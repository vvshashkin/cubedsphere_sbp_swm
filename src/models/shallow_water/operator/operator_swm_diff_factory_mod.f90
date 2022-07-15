module operator_swm_diff_factory_mod

use domain_mod,     only : domain_t
use operator_mod,   only : operator_t
use config_swm_mod, only : config_swm_t

implicit none

contains

subroutine create_swm_diff_operator(operator, swm_config, domain)

    use operator_swm_diff_mod,  only : operator_swm_diff_t
    use hordiff_factory_mod,    only : create_hordiff_operator

    type(domain_t),                 intent(in)  :: domain
    type(config_swm_t),             intent(in)  :: swm_config
    type(operator_swm_diff_t),      allocatable, intent(out) :: operator

    !type(operator_swm_diff_t), allocatable :: swm_diff_op

    integer(kind=4) :: halo_width_xy, t

    !allocate(swm_diff_op)
    allocate(operator)

    call create_hordiff_operator(operator%hordiff_uv, swm_config%hordiff_uv_name, &
                                 swm_config%uv_diff_coeff, domain)
    call create_hordiff_operator(operator%hordiff, swm_config%hordiff_h_name, &
                                 swm_config%h_diff_coeff, domain)

    !call move_alloc(swm_diff_op, operator)

end subroutine create_swm_diff_operator

end module operator_swm_diff_factory_mod
