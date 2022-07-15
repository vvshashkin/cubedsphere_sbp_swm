module operator_swm_factory_mod

use domain_mod,      only : domain_t
use operator_mod,    only : operator_t
use config_swm_mod,  only : config_swm_t
use parcomm_mod,     only : parcomm_global
use test_fields_mod, only : scalar_field_generator_t, zero_scalar_field_generator, &
                            set_scalar_test_field

implicit none

contains

subroutine create_swm_operator(operator, grav, swm_config, domain, orography_generator)

    type(domain_t),                  intent(in)  :: domain
    type(config_swm_t),              intent(in)  :: swm_config
    real(kind=8),                    intent(in)  :: grav

    class(scalar_field_generator_t), intent(in), optional :: orography_generator

    class(operator_t), allocatable,  intent(out) :: operator

    class(scalar_field_generator_t), allocatable :: orography_generator_loc

    if(present(orography_generator)) then
        orography_generator_loc = orography_generator
    else
        orography_generator_loc = zero_scalar_field_generator
    end if

    if(swm_config%swm_op_type == "vector_invariant") then
        call create_vector_invariant_swm_operator(operator, grav, swm_config, domain, &
                                                  orography_generator_loc)
    else if(swm_config%swm_op_type == "advective") then
        call create_advective_swm_operator(operator, grav, swm_config, domain, &
                                           orography_generator_loc)
    else
        call parcomm_global%abort("unknown swm_op_type: "//swm_config%swm_op_type)
    end if
end subroutine create_swm_operator

subroutine create_vector_invariant_swm_operator(operator, grav, swm_config, domain, &
                                                orography_generator)

    use operator_swm_mod,       only : operator_swm_t
    use div_factory_mod,        only : create_div_operator
    use grad_factory_mod,       only : create_grad_operator
    use curl_factory_mod,       only : create_curl_operator
    use coriolis_factory_mod,   only : create_coriolis
    use KE_factory_mod,         only : create_KE_operator
    use massflux_factory_mod,   only : create_massflux_operator
    use co2contra_factory_mod,  only : create_co2contra_operator
    use quadrature_factory_mod, only : create_quadrature
    use hordiff_factory_mod,    only : create_hordiff_operator

    use grid_field_factory_mod, only : create_grid_field

    type(domain_t),                  intent(in)  :: domain
    type(config_swm_t),              intent(in)  :: swm_config
    real(kind=8),                    intent(in)  :: grav
    class(scalar_field_generator_t), intent(in)  :: orography_generator

    class(operator_t), allocatable,  intent(out) :: operator

    type(operator_swm_t), allocatable :: swm_op

    integer(kind=4) :: halo_width_xy, t

    !WORKAROUND
    halo_width_xy = 8

    allocate(swm_op)

    swm_op%div_op = create_div_operator(domain, swm_config%div_op_name)

    swm_op%grad_op =  create_grad_operator(domain, swm_config%grad_op_name)

    call create_coriolis(swm_op%coriolis_op, swm_config%coriolis_op_name, domain)

    call create_curl_operator(swm_op%curl_op, swm_config%curl_op_name, domain)

    call create_KE_operator(swm_op%KE_op, swm_config%KE_op_name, domain)

    swm_op%massflux_op = create_massflux_operator(domain, swm_config%massflux_op_name)

    swm_op%co2contra_op = create_co2contra_operator(domain, swm_config%co2contra_op_name)

    call create_quadrature(swm_op%quadrature_h, swm_config%quadrature_name, domain%mesh_p)
    call create_quadrature(swm_op%quadrature_u, swm_config%quadrature_name, domain%mesh_u)
    call create_quadrature(swm_op%quadrature_v, swm_config%quadrature_name, domain%mesh_v)
    call create_quadrature(swm_op%quadrature_w, swm_config%quadrature_name, domain%mesh_q)

    call create_grid_field(swm_op%KE,  halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%div, halo_width_xy, 0, domain%mesh_p)


    call create_grid_field(swm_op%curl, halo_width_xy, 0, domain%mesh_q)

    call create_grid_field(swm_op%hu, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%hv, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%hu_diag, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%hv_diag, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%cor_u, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%cor_v, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%ut, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%vt, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%grad_x, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%grad_y, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%h_surf,    halo_width_xy, 0, domain%mesh_p)
    call set_scalar_test_field(swm_op%h_surf, orography_generator, domain%mesh_p, 0)

    !WORKAROUND
    swm_op%grav = grav
    !WORKAROUND
    ! do t = domain%partition%ts, domain%partition%te
    !     swm_op%h_surf%tile(t)%p = 0.0_8
    ! end do

    call create_grid_field(swm_op%KE_diag_u,  0, 0, domain%mesh_u)
    call create_grid_field(swm_op%KE_diag_v,  0, 0, domain%mesh_v)
    call create_grid_field(swm_op%PE_diag,    0, 0, domain%mesh_p)

    call move_alloc(swm_op, operator)

end subroutine create_vector_invariant_swm_operator

subroutine create_advective_swm_operator(operator, grav, swm_config, domain, &
                                         orography_generator)

    use operator_adv_swm_mod,   only : operator_adv_swm_t
    use div_factory_mod,        only : create_div_operator
    use grad_factory_mod,       only : create_grad_operator
    use coriolis_factory_mod,   only : create_coriolis
    use massflux_factory_mod,   only : create_massflux_operator
    use co2contra_factory_mod,  only : create_co2contra_operator
    use quadrature_factory_mod, only : create_quadrature
    use hordiff_factory_mod,    only : create_hordiff_operator

    use grid_field_factory_mod, only : create_grid_field

    use vector_advection_factory_mod, only : create_vector_advection_operator

    type(domain_t),                  intent(in)  :: domain
    type(config_swm_t),              intent(in)  :: swm_config
    real(kind=8),                    intent(in)  :: grav
    class(operator_t), allocatable,  intent(out) :: operator
    class(scalar_field_generator_t), intent(in)  :: orography_generator

    type(operator_adv_swm_t), allocatable :: swm_op

    integer(kind=4) :: halo_width_xy, t

    !WORKAROUND
    halo_width_xy = 8

    allocate(swm_op)

    swm_op%v_components_type = swm_config%v_components_type

    swm_op%div_op = create_div_operator(domain, swm_config%div_op_name)

    swm_op%grad_op =  create_grad_operator(domain, swm_config%grad_op_name)

    call create_coriolis(swm_op%coriolis_op, swm_config%coriolis_op_name, domain)

    swm_op%massflux_op = create_massflux_operator(domain, swm_config%massflux_op_name)

    swm_op%co2contra_op = create_co2contra_operator(domain, swm_config%co2contra_op_name)

    call create_vector_advection_operator(swm_op%adv_uv, swm_config%vector_adv_op_name, domain)

    call create_quadrature(swm_op%quadrature_h, swm_config%quadrature_name, domain%mesh_p)
    call create_quadrature(swm_op%quadrature_u, swm_config%quadrature_name, domain%mesh_u)
    call create_quadrature(swm_op%quadrature_v, swm_config%quadrature_name, domain%mesh_v)

    call create_grid_field(swm_op%div, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%h_total, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%hu, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%hv, halo_width_xy, 0, domain%mesh_v)
    call create_grid_field(swm_op%cor_u, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%cor_v, halo_width_xy, 0, domain%mesh_v)
    call create_grid_field(swm_op%ut, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%vt, halo_width_xy, 0, domain%mesh_v)
    call create_grid_field(swm_op%grad_x, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%grad_y, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%h_surf,    halo_width_xy, 0, domain%mesh_p)
    call set_scalar_test_field(swm_op%h_surf, orography_generator, domain%mesh_p, 0)

    !WORKAROUND
    swm_op%grav = grav
    !WORKAROUND
    ! do t = domain%partition%ts, domain%partition%te
    !     swm_op%h_surf%tile(t)%p = 0.0_8
    ! end do

    call create_grid_field(swm_op%KE_diag_u,  0, 0, domain%mesh_u)
    call create_grid_field(swm_op%KE_diag_v,  0, 0, domain%mesh_v)
    call create_grid_field(swm_op%PE_diag,    0, 0, domain%mesh_p)

    call move_alloc(swm_op, operator)

end subroutine create_advective_swm_operator

end module operator_swm_factory_mod
