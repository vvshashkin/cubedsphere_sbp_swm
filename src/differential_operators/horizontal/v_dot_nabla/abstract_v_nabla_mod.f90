module abstract_v_nabla_mod

use grid_field_mod,    only : grid_field_t
use mesh_mod,          only : mesh_t

implicit none

type, abstract, public :: v_nabla_operator_t
contains
    procedure(calc_v_nabla_i),         deferred :: calc_v_nabla
end type v_nabla_operator_t

abstract interface
    subroutine calc_v_nabla_i(this, f_tend, f, ut, vt, mesh)
        import v_nabla_operator_t, grid_field_t, mesh_t

        class(v_nabla_operator_t), intent(inout) :: this
        type(grid_field_t),        intent(in)    :: f      !advected field
        type(grid_field_t),        intent(in)    :: ut, vt !contravariant components
        type(mesh_t),              intent(in)    :: mesh
        type(grid_field_t),        intent(inout) :: f_tend !advective tendency
    end subroutine calc_v_nabla_i
end interface

contains

end module abstract_v_nabla_mod
