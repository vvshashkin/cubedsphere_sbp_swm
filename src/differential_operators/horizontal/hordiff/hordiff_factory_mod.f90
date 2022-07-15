module hordiff_factory_mod

use domain_mod,             only : domain_t
use grid_field_factory_mod, only : create_grid_field
use abstract_hordiff_mod,   only : hordiff_operator_t

implicit none

contains

subroutine create_hordiff_operator(hordiff_op, hordiff_op_name, hordiff_coeff, domain)

    class(hordiff_operator_t), allocatable, intent(out) :: hordiff_op
    character(len=*),                       intent(in)  :: hordiff_op_name
    real(kind=8),                           intent(in)  :: hordiff_coeff
    type(domain_t),                         intent(in)  :: domain

    select case(hordiff_op_name)
    case("hordiff_c_biharm_div")
        call create_Cgrid_hordiff_div_operator(hordiff_op, hordiff_coeff, domain)
    case("hordiff_c_biharm_curl")
        call create_Cgrid_hordiff_curl_operator(hordiff_op, hordiff_coeff, domain)
    case("hordiff_c_biharm_curl_div")
        call create_Cgrid_hordiff_curl_div_operator(hordiff_op, hordiff_coeff, domain)
!    case("hordiff_colocated")
!        hordiff_op = hordiff_colocated_t()
    case("hordiff_scalar_Ah_sbp_21_narrow")
        call create_laplace_based_hordiff_operator(hordiff_op, "laplace_ah_sbp21_narrow", hordiff_coeff, domain)
    case("hordiff_scalar_Ah_sbp_42_narrow")
        call create_laplace_based_hordiff_operator(hordiff_op, "laplace_ah_sbp42_narrow", hordiff_coeff, domain)
    case("hordiff_scalar_Ah_sbp_63_narrow")
        call create_laplace_based_hordiff_operator(hordiff_op, "laplace_ah_sbp63_narrow", hordiff_coeff, domain)
    case("hordiff_scalar_Ah")
        call create_scalar_hordiff_operator(hordiff_op, hordiff_coeff, domain, "Ah", isscalar=.true.)
    case("hordiff_scalar_C")
        call create_scalar_hordiff_operator(hordiff_op, hordiff_coeff, domain, "C", isscalar=.true.)
    case("hordiff_vec_Ah")
        call create_colocated_hordiff_operator(hordiff_op, hordiff_coeff, domain, "Ah")
    case("hordiff_vec_xyz_Ah")
        call create_colocated_vec_xyz_hordiff_operator(hordiff_op, hordiff_coeff, domain, "Ah")
    case("hordiff_vec_xyz_Ah_sbp_21_narrow")
        call create_Ah_lap_based_vec_xyz_hordiff_operator(hordiff_op, "laplace_ah_sbp21_narrow", hordiff_coeff, domain)
    case("hordiff_vec_xyz_Ah_sbp_42_narrow")
        call create_Ah_lap_based_vec_xyz_hordiff_operator(hordiff_op, "laplace_ah_sbp42_narrow", hordiff_coeff, domain)
    case("hordiff_vec_xyz_Ah_sbp_63_narrow")
        call create_Ah_lap_based_vec_xyz_hordiff_operator(hordiff_op, "laplace_ah_sbp63_narrow", hordiff_coeff, domain)
    case default
        call domain%parcomm%abort("Unknown hordiff_op_name "//hordiff_op_name)
    end select

end subroutine create_hordiff_operator

subroutine create_Cgrid_hordiff_curl_div_operator(hordiff_op, hordiff_coeff, domain)

    use hordiff_Cgrid_mod,     only : hordiff_c_curl_div_t

    class(hordiff_operator_t), allocatable, intent(out) :: hordiff_op
    real(kind=8),                           intent(in)  :: hordiff_coeff
    type(domain_t),                         intent(in)  :: domain

    type(hordiff_c_curl_div_t), allocatable :: hordiff_curl_div

    integer(kind=4) :: halo_width
    real(kind=8)    :: hx

    allocate(hordiff_curl_div)

    !WORKAROUND
    halo_width = 4

    call create_grid_field(hordiff_curl_div%u_tend,  halo_width, 0, domain%mesh_u)
    call create_grid_field(hordiff_curl_div%v_tend,  halo_width, 0, domain%mesh_v)

    call create_Cgrid_hordiff_curl_operator(hordiff_curl_div%curl_diff_op, hordiff_coeff, domain)
    call create_Cgrid_hordiff_div_operator(hordiff_curl_div%div_diff_op, hordiff_coeff, domain)

    call move_alloc(hordiff_curl_div, hordiff_op)

end subroutine create_Cgrid_hordiff_curl_div_operator

subroutine create_Cgrid_hordiff_div_operator(hordiff_op, hordiff_coeff, domain)

    use hordiff_Cgrid_mod,     only : hordiff_c_div_t
    use div_factory_mod,       only : create_div_operator
    use grad_factory_mod,      only : create_grad_operator
    use co2contra_factory_mod, only : create_co2contra_operator

    class(hordiff_operator_t), allocatable, intent(out) :: hordiff_op
    real(kind=8),                           intent(in)  :: hordiff_coeff
    type(domain_t),                         intent(in)  :: domain

    type(hordiff_c_div_t), allocatable :: hordiff_div

    integer(kind=4) :: halo_width
    real(kind=8)    :: hx

    allocate(hordiff_div)

    !WORKAROUND
    halo_width = 4

    call create_grid_field(hordiff_div%div, halo_width, 0, domain%mesh_p)
    call create_grid_field(hordiff_div%ut,  halo_width, 0, domain%mesh_u)
    call create_grid_field(hordiff_div%vt,  halo_width, 0, domain%mesh_v)

    hordiff_div%div_op = create_div_operator(domain, "divergence_c_sbp21")
    hordiff_div%grad_op = create_grad_operator(domain, "gradient_c_sbp21")
    hordiff_div%co2contra_op = create_co2contra_operator(domain, "co2contra_c_sbp21")

    hx = domain%mesh_p%tile(domain%mesh_p%ts)%hx

    hordiff_div%diff_coeff = hordiff_coeff*domain%mesh_p%scale*hx
    hordiff_div%diff_order = 2

    call move_alloc(hordiff_div, hordiff_op)

end subroutine create_Cgrid_hordiff_div_operator

subroutine create_Cgrid_hordiff_curl_operator(hordiff_op, hordiff_coeff, domain)

    use hordiff_Cgrid_mod,     only : hordiff_c_curl_t
    use curl_factory_mod,      only : create_curl_operator
    use grad_perp_factory_mod, only : create_grad_perp_operator
    use co2contra_factory_mod, only : create_co2contra_operator

    class(hordiff_operator_t), allocatable, intent(out) :: hordiff_op
    real(kind=8),                           intent(in)  :: hordiff_coeff
    type(domain_t),                         intent(in)  :: domain

    type(hordiff_c_curl_t), allocatable :: hordiff_curl

    integer(kind=4) :: halo_width
    real(kind=8)    :: hx

    allocate(hordiff_curl)

    !WORKAROUND
    halo_width = 4

    call create_grid_field(hordiff_curl%curl, halo_width, 0, domain%mesh_q)
    call create_grid_field(hordiff_curl%ut,   halo_width, 0, domain%mesh_u)
    call create_grid_field(hordiff_curl%vt,   halo_width, 0, domain%mesh_v)

    call create_curl_operator(hordiff_curl%curl_op, "curl_c_sbp21", domain)
    call create_grad_perp_operator(hordiff_curl%grad_perp_op, "grad_perp_c_sbp21", domain)
    hordiff_curl%co2contra_op = create_co2contra_operator(domain, "co2contra_c_sbp21_new")

    hx = domain%mesh_p%tile(domain%mesh_p%ts)%hx

    hordiff_curl%diff_coeff = hordiff_coeff*domain%mesh_p%scale*hx
    hordiff_curl%diff_order = 2

    call move_alloc(hordiff_curl, hordiff_op)

end subroutine create_Cgrid_hordiff_curl_operator

subroutine create_laplace_based_hordiff_operator(hordiff_op, laplace_op_name, hordiff_coeff, domain)

    use hordiff_scalar_mod,    only : hordiff_scalar_t
    use halo_factory_mod,      only : create_halo_procedure
    use laplace_factory_mod,   only : create_laplace_operator

    class(hordiff_operator_t), allocatable, intent(out) :: hordiff_op
    real(kind=8),                           intent(in)  :: hordiff_coeff
    type(domain_t),                         intent(in)  :: domain
    character(len=*),                       intent(in)  :: laplace_op_name

    type(hordiff_scalar_t), allocatable :: hordiff_scalar

    integer(kind=4) :: halo_width
    real(kind=8)    :: hx

    allocate(hordiff_scalar)

    !WORKAROUND
    halo_width = 9

    call create_laplace_operator(hordiff_scalar%laplace_op, laplace_op_name, domain)
    call create_grid_field(hordiff_scalar%f_tend_inter, halo_width, 0, domain%mesh_xy)

    hx = domain%mesh_xy%tile(domain%mesh_o%ts)%hx

    hordiff_scalar%diff_coeff = hordiff_coeff*domain%mesh_xy%scale*hx
    hordiff_scalar%diff_order = 2

    call move_alloc(hordiff_scalar, hordiff_op)

end subroutine create_laplace_based_hordiff_operator

subroutine create_scalar_hordiff_operator(hordiff_op, hordiff_coeff, domain, &
                                          staggering, isscalar)

    use hordiff_scalar_mod,    only : hordiff_scalar_t
    use div_factory_mod,       only : create_div_operator
    use grad_factory_mod,      only : create_grad_operator
    use co2contra_factory_mod, only : create_co2contra_operator
    use halo_factory_mod,      only : create_halo_procedure
    use laplace_factory_mod,   only : create_laplace_operator

    class(hordiff_operator_t), allocatable, intent(out) :: hordiff_op
    real(kind=8),                           intent(in)  :: hordiff_coeff
    type(domain_t),                         intent(in)  :: domain
    character(len=*),                       intent(in)  :: staggering
    logical,                                intent(in)  :: isscalar

    type(hordiff_scalar_t), allocatable :: hordiff_scalar

    integer(kind=4) :: halo_width
    real(kind=8)    :: hx

    allocate(hordiff_scalar)

    !WORKAROUND
    halo_width = 5

    select case(staggering)
    case ("Ah")
        !call create_laplace_operator(hordiff_scalar%laplace_op,"divgrad_laplace_ch_sbp21",domain)
        call create_laplace_operator(hordiff_scalar%laplace_op,"laplace_ch_halo4",domain)
        call create_grid_field(hordiff_scalar%f_tend_inter, 8, 0, domain%mesh_xy)
        ! call create_grid_field(hordiff_scalar%div, halo_width, 0, domain%mesh_xy)
        ! call create_grid_field(hordiff_scalar%gx,  halo_width, 0, domain%mesh_y)
        ! call create_grid_field(hordiff_scalar%gy,  halo_width, 0, domain%mesh_x)
        ! call create_grid_field(hordiff_scalar%gxt, halo_width, 0, domain%mesh_y)
        ! call create_grid_field(hordiff_scalar%gyt, halo_width, 0, domain%mesh_x)
        !
        ! hordiff_scalar%div_op = create_div_operator(domain, "divergence_ch_sbp21")
        ! hordiff_scalar%grad_op = create_grad_operator(domain, "gradient_ch_sbp21")
        ! hordiff_scalar%co2contra_op = create_co2contra_operator(domain, "co2contra_ch_sbp21")
        ! if(isscalar) then
        !     call create_halo_procedure(hordiff_scalar%edge_sync, domain, 1, "Ah_scalar_sync")
        ! end if
    case ("C", "A")
        call create_laplace_operator(hordiff_scalar%laplace_op,"divgrad_laplace_c_sbp21",domain)
        call create_grid_field(hordiff_scalar%f_tend_inter, halo_width, 0, domain%mesh_o)
        ! call create_grid_field(hordiff_scalar%div, halo_width, 0, domain%mesh_o)
        ! call create_grid_field(hordiff_scalar%gx,  halo_width, 0, domain%mesh_x)
        ! call create_grid_field(hordiff_scalar%gy,  halo_width, 0, domain%mesh_y)
        ! call create_grid_field(hordiff_scalar%gxt, halo_width, 0, domain%mesh_x)
        ! call create_grid_field(hordiff_scalar%gyt, halo_width, 0, domain%mesh_y)
        !
        ! hordiff_scalar%div_op = create_div_operator(domain, "divergence_c_sbp21")
        ! hordiff_scalar%grad_op = create_grad_operator(domain, "gradient_c_sbp21")
        ! hordiff_scalar%co2contra_op = create_co2contra_operator(domain, "co2contra_c_sbp21")
    case default
        call domain%parcomm%abort("This staggering: "//staggering//&
                                  "is not currently implemented in create_hordiff_scalar")
    end select

    hx = domain%mesh_xy%tile(domain%mesh_o%ts)%hx

    hordiff_scalar%diff_coeff = hordiff_coeff*domain%mesh_xy%scale*hx
    hordiff_scalar%diff_order = 2

    call move_alloc(hordiff_scalar, hordiff_op)

end subroutine create_scalar_hordiff_operator

subroutine create_colocated_hordiff_operator(hordiff_op, hordiff_coeff, domain, staggering)

    use hordiff_colocated_mod,    only : hordiff_colocated_t
    use halo_factory_mod,         only : create_vector_halo_procedure

    class(hordiff_operator_t), allocatable, intent(out) :: hordiff_op
    real(kind=8),                           intent(in)  :: hordiff_coeff
    character(len=*),                       intent(in)  :: staggering
    type(domain_t),                         intent(in)  :: domain

    type(hordiff_colocated_t), allocatable :: hordiff_uv

    allocate(hordiff_uv)

    select case(staggering)
    case ("Ah")
        call create_scalar_hordiff_operator(hordiff_uv%hordiff_1comp, hordiff_coeff, domain, &
                                           "Ah", isscalar=.false.) !WORKAROUND
        call create_vector_halo_procedure(hordiff_uv%edge_sync, domain, 1, "ecs_Ah_vec_sync_covariant")
    case default
        call domain%parcomm%abort("This staggering: "//staggering//&
                                  "is not currently implemented in create_colocated)hordiff_operator")
    end select

    call move_alloc(hordiff_uv, hordiff_op)

end subroutine create_colocated_hordiff_operator
subroutine create_colocated_vec_xyz_hordiff_operator(hordiff_op, hordiff_coeff, domain, staggering)

    use hordiff_colocated_mod,    only : hordiff_colocated_xyz_t
    use halo_factory_mod,         only : create_vector_halo_procedure
    use grid_field_factory_mod,   only : create_grid_field

    class(hordiff_operator_t), allocatable, intent(out) :: hordiff_op
    real(kind=8),                           intent(in)  :: hordiff_coeff
    character(len=*),                       intent(in)  :: staggering
    type(domain_t),                         intent(in)  :: domain

    type(hordiff_colocated_xyz_t), allocatable :: hordiff_uv

    integer(kind=4) :: halo_width

    !WORKAROUND
    halo_width = 5

    allocate(hordiff_uv)

    select case(staggering)
    case ("Ah")
        ! call create_laplace_based_hordiff_operator(hordiff_uv%hordiff_1comp, "laplace_ah_sbp21_narrow", hordiff_coeff, domain)
        call create_scalar_hordiff_operator(hordiff_uv%hordiff_1comp, hordiff_coeff, domain, &
                                           "Ah", isscalar=.false.) !WORKAROUND
        call create_vector_halo_procedure(hordiff_uv%edge_sync, domain, 1, "ecs_Ah_vec_sync_covariant")

        call create_grid_field(hordiff_uv%vx,      halo_width, 0, domain%mesh_p)
        call create_grid_field(hordiff_uv%vy,      halo_width, 0, domain%mesh_p)
        call create_grid_field(hordiff_uv%vz,      halo_width, 0, domain%mesh_p)
        call create_grid_field(hordiff_uv%vx_tend,          0, 0, domain%mesh_p)
        call create_grid_field(hordiff_uv%vy_tend,          0, 0, domain%mesh_p)
        call create_grid_field(hordiff_uv%vz_tend,          0, 0, domain%mesh_p)

    case default
        call domain%parcomm%abort("This staggering: "//staggering//&
                                  "is not currently implemented in create_colocated_vec_xyz_hordiff_operator")
    end select

    call move_alloc(hordiff_uv, hordiff_op)

end subroutine create_colocated_vec_xyz_hordiff_operator
subroutine create_Ah_lap_based_vec_xyz_hordiff_operator(hordiff_op, laplace_op_name, hordiff_coeff, domain)

    use hordiff_colocated_mod,    only : hordiff_colocated_xyz_t
    use halo_factory_mod,         only : create_vector_halo_procedure
    use grid_field_factory_mod,   only : create_grid_field

    class(hordiff_operator_t), allocatable, intent(out) :: hordiff_op
    real(kind=8),                           intent(in)  :: hordiff_coeff
    character(len=*),                       intent(in)  :: laplace_op_name
    type(domain_t),                         intent(in)  :: domain

    type(hordiff_colocated_xyz_t), allocatable :: hordiff_uv

    integer(kind=4) :: halo_width

    !WORKAROUND
    halo_width = 8

    allocate(hordiff_uv)

    call create_laplace_based_hordiff_operator(hordiff_uv%hordiff_1comp, laplace_op_name, hordiff_coeff, domain)

    call create_vector_halo_procedure(hordiff_uv%edge_sync, domain, 1, "ecs_Ah_vec_sync_covariant")

    call create_grid_field(hordiff_uv%vx,      halo_width, 0, domain%mesh_xy)
    call create_grid_field(hordiff_uv%vy,      halo_width, 0, domain%mesh_xy)
    call create_grid_field(hordiff_uv%vz,      halo_width, 0, domain%mesh_xy)
    call create_grid_field(hordiff_uv%vx_tend,          0, 0, domain%mesh_xy)
    call create_grid_field(hordiff_uv%vy_tend,          0, 0, domain%mesh_xy)
    call create_grid_field(hordiff_uv%vz_tend,          0, 0, domain%mesh_xy)

    call move_alloc(hordiff_uv, hordiff_op)

end subroutine create_Ah_lap_based_vec_xyz_hordiff_operator
end module hordiff_factory_mod
