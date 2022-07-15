module abstract_adv_z_mod

use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t

implicit none

type, abstract :: adv_z_t
    contains
    procedure calc_z_adv
    procedure(calc_z_adv_tile), deferred :: calc_z_adv_tile
end type adv_z_t

abstract interface
    subroutine calc_z_adv_tile(this, f_tend, f, eta_dot, mesh, scale)
        import adv_z_t, tile_field_t, tile_mesh_t
        class(adv_z_t),     intent(in)    :: this
        type(tile_field_t), intent(in)    :: f, eta_dot
        type(tile_mesh_t),  intent(in)    :: mesh
        real(kind=8),       intent(in)    :: scale
        !output
        type(tile_field_t), intent(inout) :: f_tend
    end subroutine calc_z_adv_tile
end interface

contains

subroutine calc_z_adv(this, f_tend, f, eta_dot, mesh)
    class(adv_z_t),     intent(in)    :: this
    type(grid_field_t), intent(in)    :: f, eta_dot
    type(mesh_t),       intent(in)    :: mesh
    !output
    type(grid_field_t), intent(inout) :: f_tend

    integer(kind=4) :: t

    do t=mesh%ts, mesh%te
        call this%calc_z_adv_tile(f_tend%tile(t), f%tile(t), eta_dot%tile(t), &
                                  mesh%tile(t), mesh%vertical_scale)
    end do
end subroutine calc_z_adv

end module abstract_adv_z_mod
