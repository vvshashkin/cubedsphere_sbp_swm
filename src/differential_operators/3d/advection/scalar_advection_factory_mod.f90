module scalar_advection_factory_mod

use abstract_scalar_advection3d_mod, only : scalar_advection3d_t
use domain_mod,                      only : domain_t
use parcomm_mod,                     only : parcomm_global

implicit none

contains

subroutine create_scalar_advection3d_operator(adv_op, scalar_advection_op_name, &
                                              hor_adv_op_name, z_adv_op_name, domain)
    class(scalar_advection3d_t), allocatable, intent(out) :: adv_op
    character(len=*), intent(in) :: scalar_advection_op_name, hor_adv_op_name, &
                                    z_adv_op_name
    type(domain_t),   intent(in) :: domain

    select case(scalar_advection_op_name)
    case("advection_p_staggered")
        call create_p_3d_advection(adv_op, hor_adv_op_name, z_adv_op_name, domain)
    case("advection_w_staggered")
        call create_w_3d_advection(adv_op, hor_adv_op_name, z_adv_op_name, domain)
    case default
        call parcomm_global%abort("create_scalar_advection3d_operator, unknown "//&
                                  "scalar_advection_op_name: "// scalar_advection_op_name)
    end select
end subroutine create_scalar_advection3d_operator

subroutine create_p_3d_advection(adv_op, hor_adv_op_name, z_adv_op_name, domain)

    use v_nabla_factory_mod, only : create_v_nabla_hor_operator
    use adv_z_factory_mod,   only : create_adv_z_operator
    use halo_factory_mod,    only : create_halo_procedure
    use advection_p_3d_mod,  only : advection_p_C3d_t

    use interpolator2d_factory_mod,    only : create_vec2vec_interpolator2d
    use vertical_operator_factory_mod, only : create_vertical_operator
    use grid_field_factory_mod,        only : create_grid_field

    class(scalar_advection3d_t), allocatable, intent(out) :: adv_op
    character(len=*), intent(in) :: hor_adv_op_name, z_adv_op_name
    type(domain_t),   intent(in) :: domain

    type(advection_p_C3d_t), allocatable :: adv_p3d
    integer(kind=4) :: halo_width

    allocate(adv_p3d)

    call create_vec2vec_interpolator2d(adv_p3d%interp_uv2p_op, "interp2d_uv2pvec_C_sbp42", domain)
    call create_vertical_operator(adv_p3d%interp_w2p_op, "vertical_interp_w2p_sbp42")

    call create_v_nabla_hor_operator(adv_p3d%v_nabla_op,halo_width,hor_adv_op_name)

    adv_p3d%halo_width = halo_width
    call create_halo_procedure(adv_p3d%halo_f,domain,max(halo_width,2),"ECS_O")

    call create_adv_z_operator(adv_p3d%adv_z, z_adv_op_name)
    call create_grid_field(adv_p3d%up, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%vp, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%eta_dot_p, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%f_tend_z, 0, 0, domain%mesh_p)

    call move_alloc(adv_p3d, adv_op)

end subroutine create_p_3d_advection

subroutine create_w_3d_advection(adv_op, hor_adv_op_name, z_adv_op_name, domain)

    use v_nabla_factory_mod, only : create_v_nabla_hor_operator
    use adv_z_factory_mod,   only : create_adv_z_operator
    use halo_factory_mod,    only : create_halo_procedure
    use advection_w_3d_mod,  only : advection_w_C3d_t

    use interpolator_uv2w_factory_mod, only : create_uv2w_interpolator
    use grid_field_factory_mod,        only : create_grid_field

    class(scalar_advection3d_t), allocatable, intent(out) :: adv_op
    character(len=*), intent(in) :: hor_adv_op_name, z_adv_op_name
    type(domain_t),   intent(in) :: domain

    type(advection_w_C3d_t), allocatable :: adv_w3d
    integer(kind=4) :: halo_width

    allocate(adv_w3d)

    call create_uv2w_interpolator(adv_w3d%interp_uv2w_op, "uv2w_staggered", &
                                  "interp2d_uv2pvec_C_sbp42", "vertical_interp_p2w_sbp42", &
                                  domain)

    call create_v_nabla_hor_operator(adv_w3d%v_nabla_op,halo_width,hor_adv_op_name)

    adv_w3d%halo_width = halo_width
    call create_halo_procedure(adv_w3d%halo_f,domain,max(halo_width,2),"ECS_Oz")

    call create_adv_z_operator(adv_w3d%adv_z, z_adv_op_name)
    call create_grid_field(adv_w3d%uw, 0, 0, domain%mesh_w)
    call create_grid_field(adv_w3d%vw, 0, 0, domain%mesh_w)
    call create_grid_field(adv_w3d%f_tend_z, 0, 0, domain%mesh_w)

    call move_alloc(adv_w3d, adv_op)

end subroutine create_w_3d_advection

end module scalar_advection_factory_mod
