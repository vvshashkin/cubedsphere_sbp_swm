module advection_p_3d_mod

use abstract_scalar_advection3d_mod, only : scalar_advection3d_t
use grid_field_mod,                  only : grid_field_t, tile_field_t
use domain_mod,                      only : domain_t
use mesh_mod,                        only : tile_mesh_t
use halo_mod,                        only : halo_t
use abstract_interpolators2d_mod,    only : interpolator2d_vec2vec_t
use abstract_vertical_operator_mod,  only : vertical_operator_t
use abstract_v_nabla_mod,            only : v_nabla_operator_t
use abstract_adv_z_mod,              only : adv_z_t

type, extends(scalar_advection3d_t) :: advection_p_C3d_t
    class(interpolator2d_vec2vec_t), allocatable :: interp_uv2p_op
    class(vertical_operator_t), allocatable      :: interp_w2p_op
    class(halo_t), allocatable                   :: halo_f
    integer(kind=4)                              :: halo_width
    class(v_nabla_operator_t), allocatable       :: v_nabla_op
    class(adv_z_t), allocatable                  :: adv_z
    type(grid_field_t)                           :: up, vp, eta_dot_p, f_tend_z
    contains
    procedure calc_adv3d => calc_adv3d_C
end type advection_p_C3d_t

contains

subroutine calc_adv3d_C(this,f_tend,f,u,v,eta_dot,domain)
    class(advection_p_C3d_t),    intent(inout) :: this
    type(grid_field_t),          intent(inout) :: f, u, v, eta_dot
    type(domain_t),              intent(in)    :: domain
    !output
    type(grid_field_t),          intent(inout) :: f_tend

    call this%interp_uv2p_op%interp2d_vec2vec(this%up,this%vp,u,v,domain)
    call this%interp_w2p_op%apply(this%eta_dot_p, eta_dot, domain)
    call this%halo_f%get_halo_scalar(f, domain, this%halo_width)

    call this%v_nabla_op%calc_v_nabla(f_tend, f, this%up, this%vp, domain%mesh_p)
    call this%adv_z%calc_z_adv(this%f_tend_z,f,this%eta_dot_p,domain%mesh_p)

    call f_tend%update(1.0_8,this%f_tend_z,domain%mesh_p)
end subroutine calc_adv3d_C

end module advection_p_3d_mod
