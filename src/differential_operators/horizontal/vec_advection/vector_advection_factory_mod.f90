module vector_advection_factory_mod

use domain_mod,                     only : domain_t
use abstract_vector_advection_mod,  only : vector_advection_operator_t
use parcomm_mod,                    only : parcomm_global
use grid_field_mod,                 only : grid_field_t
use v_nabla_mod,                    only : v_nabla_up4_operator_t, &
                                           v_nabla_up3_operator_t, &
                                           v_nabla_up1_operator_t, &
                                           v_nabla_c2_operator_t,  &
                                           v_nabla_c4_operator_t

implicit none

private

public :: create_vector_advection_operator

contains

subroutine create_vector_advection_operator(vec_advection_op, vec_advection_op_name, domain)

    class(vector_advection_operator_t), allocatable, intent(out) :: vec_advection_op
    character(len=*),                                intent(in)  :: vec_advection_op_name
    type(domain_t),                                  intent(in)  :: domain

    select case(vec_advection_op_name)

    case("vector_advection_Ah21")
        call create_vector_advection_Ah_covariant(vec_advection_op, "d21", 1, domain)
    case("vector_advection_Ah42")
        call create_vector_advection_Ah_covariant(vec_advection_op, "d42", 3, domain)
    case("vector_advection_Ah63")
        call create_vector_advection_Ah_covariant(vec_advection_op, "d63", 5, domain)
    case("vector_advection_C_up4")
        call create_vector_advection_C(vec_advection_op, v_nabla_up4_operator_t(), &
                                       "interp2d_uv2pvec_C_sbp42", "interp2d_pvec2uv_C_sbp42", &
                                                                                     3,domain)
    case("vector_advection_C_up3")
        call create_vector_advection_C(vec_advection_op, v_nabla_up3_operator_t(), &
                                       "interp2d_uv2pvec_C_sbp42", "interp2d_pvec2uv_C_sbp42", &
                                                                                     3,domain)
    case("vector_advection_C_up1")
        call create_vector_advection_C(vec_advection_op, v_nabla_up1_operator_t(), &
                                       "interp2d_uv2pvec_C_sbp42", "interp2d_pvec2uv_C_sbp42", &
                                                                                     3,domain)
    case("vector_advection_C_c2")
        call create_vector_advection_C(vec_advection_op, v_nabla_c2_operator_t(), &
                                       "interp2d_uv2pvec_C_sbp42", "interp2d_pvec2uv_C_sbp42", &
                                                                                     3,domain)
    case("vector_advection_C_c4")
        call create_vector_advection_C(vec_advection_op, v_nabla_c4_operator_t(), &
                                       "interp2d_uv2pvec_C_sbp42", "interp2d_pvec2uv_C_sbp42", &
                                                                                     3,domain)
    case default
        call parcomm_global%abort("Unknown vector advection operator: "//vec_advection_op_name)
    end select

end subroutine create_vector_advection_operator

subroutine create_vector_advection_Ah_covariant(vec_advection_op,sbp_operator_name, &
                                                halo_width, domain)

    use vector_advection_Ah_mod,   only : vector_advection_Ah_t
    use exchange_factory_mod,      only : create_xy_points_halo_exchange
    use halo_factory_mod,          only : create_vector_halo_procedure
    use v_nabla_sbp_factory_mod,   only : create_v_nabla_sbp_operator

    class(vector_advection_operator_t), allocatable, intent(out) :: vec_advection_op
    character(len=*),                                intent(in)  :: sbp_operator_name
    integer(kind=4),                                 intent(in)  :: halo_width
    type(domain_t),                                  intent(in)  :: domain

    type(vector_advection_Ah_t), allocatable :: vec_advection_Ah_op

    allocate(vec_advection_Ah_op)

    vec_advection_Ah_op%v_nabla_op = create_v_nabla_sbp_operator(sbp_operator_name)
    call create_vector_halo_procedure(vec_advection_Ah_op%sync_edges_cov,domain,0, &
                                                    "ecs_Ah_vec_sync_covariant")
    call create_vector_halo_procedure(vec_advection_Ah_op%sync_edges_contra,domain,0, &
                                                    "ecs_Ah_vec_sync_contra")

    vec_advection_Ah_op%exch_uv_interior =  &
                    create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                                 domain%topology,  halo_width, 'full')

    call move_alloc(vec_advection_Ah_op, vec_advection_op)

end subroutine create_vector_advection_Ah_covariant

subroutine create_vector_advection_C(vec_advection_op, v_nabla_op, uv2p_interp_name, &
                                     p2uv_interp_name, halo_width, domain)

    use vector_advection_C_mod,       only : vector_advection_C_t
    use grid_field_factory_mod,       only : create_grid_field
    use interpolator2d_factory_mod,   only : create_vec2vec_interpolator2d
    use halo_factory_mod,             only : create_vector_halo_procedure
    use abstract_v_nabla_mod,         only : v_nabla_operator_t

    character(len=*),  intent(in)         :: uv2p_interp_name, p2uv_interp_name
    class(v_nabla_operator_t), intent(in) :: v_nabla_op
    integer(kind=4),   intent(in)         :: halo_width
    type(domain_t),    intent(in)         :: domain

    class(vector_advection_operator_t), allocatable, intent(out) :: vec_advection_op

    type(vector_advection_C_t), allocatable :: vec_advection_C_op

    allocate(vec_advection_C_op)

    vec_advection_C_op%v_nabla_op = v_nabla_op

    call create_grid_field(vec_advection_C_op%u_at_v, 0, 0, domain%mesh_v)
    call create_grid_field(vec_advection_C_op%v_at_u, 0, 0, domain%mesh_u)
    call create_grid_field(vec_advection_C_op%uh,halo_width+1, 0, domain%mesh_p)
    call create_grid_field(vec_advection_C_op%vh,halo_width+1, 0, domain%mesh_p)

    call create_vec2vec_interpolator2d(vec_advection_C_op%interp_v2h_op, uv2p_interp_name, domain)
    call create_vec2vec_interpolator2d(vec_advection_C_op%interp_h2v_op, p2uv_interp_name, domain)

    call create_vector_halo_procedure(vec_advection_C_op%halo_uv, domain,max(halo_width,2),"ecs_C_vec")
    call create_vector_halo_procedure(vec_advection_C_op%tendency_edge_sync, domain, &
                                                                              1,"C_vec_default")

    call move_alloc(vec_advection_C_op, vec_advection_op)

end subroutine create_vector_advection_C

end module vector_advection_factory_mod
