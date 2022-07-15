module advection_w_3d_mod

use abstract_scalar_advection3d_mod, only : scalar_advection3d_t
use grid_field_mod,                  only : grid_field_t, tile_field_t
use domain_mod,                      only : domain_t
use mesh_mod,                        only : tile_mesh_t
use halo_mod,                        only : halo_t
use abstract_interpolators3d_mod,    only : interpolator_uv2w_t
use abstract_v_nabla_mod,            only : v_nabla_operator_t
use abstract_adv_z_mod,              only : adv_z_t

type, extends(scalar_advection3d_t) :: advection_w_C3d_t
    class(interpolator_uv2w_t), allocatable :: interp_uv2w_op
    class(halo_t), allocatable              :: halo_f
    integer(kind=4)                         :: halo_width
    class(v_nabla_operator_t), allocatable  :: v_nabla_op
    class(adv_z_t), allocatable             :: adv_z
    type(grid_field_t)                      :: uw, vw, f_tend_z
    contains
    procedure calc_adv3d => calc_adv3d_C
end type advection_w_C3d_t

contains

subroutine calc_adv3d_C(this,f_tend,f,u,v,eta_dot,domain)
    class(advection_w_C3d_t),    intent(inout) :: this
    type(grid_field_t),          intent(inout) :: f, u, v, eta_dot
    type(domain_t),              intent(in)    :: domain
    !output
    type(grid_field_t),          intent(inout) :: f_tend

    call this%interp_uv2w_op%interp_uv2w(this%uw,this%vw,u,v,domain)
    call this%halo_f%get_halo_scalar(f, domain, this%halo_width)

    call this%v_nabla_op%calc_v_nabla(f_tend, f, this%uw, this%vw, domain%mesh_w)
    call this%adv_z%calc_z_adv(this%f_tend_z,f,eta_dot,domain%mesh_w)

    call f_tend%update(1.0_8,this%f_tend_z,domain%mesh_w)
end subroutine calc_adv3d_C

end module advection_w_3d_mod
