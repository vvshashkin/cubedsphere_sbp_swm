module operator_advection_factory_mod

use domain_mod,   only : domain_t
use operator_mod, only : operator_t

implicit none

contains

subroutine create_advection_operator(operator, massflux_operator_name, div_operator_name, domain)

    use operator_advection_mod, only : operator_advection_t
    use div_factory_mod,        only : create_div_operator
    use massflux_factory_mod,   only : create_massflux_operator
    use grid_field_factory_mod, only : create_grid_field

    type(domain_t),                 intent(in)  :: domain
    character(len=*),               intent(in)  :: massflux_operator_name
    character(len=*),               intent(in)  :: div_operator_name
    class(operator_t), allocatable, intent(out) :: operator

    type(operator_advection_t), allocatable :: advection_operator

    allocate(advection_operator)

    advection_operator%div_op  = create_div_operator(domain, div_operator_name)
    advection_operator%flux_op = create_massflux_operator(domain, massflux_operator_name)

    !WORKAROUND
    call create_grid_field(advection_operator%hu, 10, 0, domain%mesh_p)
    call create_grid_field(advection_operator%hv, 10, 0, domain%mesh_p)

    call move_alloc(advection_operator, operator)

end subroutine create_advection_operator

end module operator_advection_factory_mod
