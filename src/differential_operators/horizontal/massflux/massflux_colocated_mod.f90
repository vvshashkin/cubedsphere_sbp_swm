module massflux_colocated_mod

use abstract_massflux_mod, only : massflux_operator_t
use grid_field_mod,        only : grid_field_t
use domain_mod,            only : domain_t

implicit none

type, extends(massflux_operator_t), public :: massflux_colocated_t

contains

procedure :: calc_massflux => calc_colocated_massflux

end type massflux_colocated_t

contains

subroutine calc_colocated_massflux(this, fx, fy, f, u, v, domain)
    class(massflux_colocated_t), intent(inout) :: this
    type(domain_t),              intent(in)    :: domain
    type(grid_field_t),          intent(inout) :: f, u, v
    !output:
    type(grid_field_t),          intent(inout) :: fx, fy

    call fx%assign_prod(1.0_8, f, u, domain%mesh_p)
    call fy%assign_prod(1.0_8, f, v, domain%mesh_p)

end subroutine calc_colocated_massflux

end module massflux_colocated_mod
