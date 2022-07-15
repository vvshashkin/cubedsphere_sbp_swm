module abstract_scalar_advection3d_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

implicit none

type, abstract :: scalar_advection3d_t
    contains
    procedure(calc_adv3d), deferred :: calc_adv3d
end type scalar_advection3d_t

abstract interface
    subroutine calc_adv3d(this,f_tend,f,u,v,eta_dot,domain)
        import scalar_advection3d_t, grid_field_t, domain_t
        class(scalar_advection3d_t), intent(inout) :: this
        type(grid_field_t),          intent(inout) :: f, u, v, eta_dot
        type(domain_t),              intent(in)    :: domain
        !output
        type(grid_field_t),          intent(inout) :: f_tend
    end subroutine calc_adv3d
end interface

end module abstract_scalar_advection3d_mod
