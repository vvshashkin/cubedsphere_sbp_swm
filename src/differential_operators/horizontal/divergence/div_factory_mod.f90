module div_factory_mod

use domain_mod,        only : domain_t
use mesh_mod,          only : mesh_t
use abstract_div_mod,  only : div_operator_t
use parcomm_mod,       only : parcomm_global

implicit none

contains

function create_div_operator(domain, div_operator_name) result(div)
    type(domain_t),    intent(in)    :: domain
    character(len=*),  intent(in)    :: div_operator_name

    class(div_operator_t), allocatable :: div

    select case(div_operator_name)
    case("divergence_c2")
        div = create_div_c2_operator(domain)
    case("divergence_c_sbp21")
        div = create_div_c_sbp21_operator(domain)
    case("divergence_c_sbp42")
        div = create_div_c_sbp42_operator(domain)
    case("divergence_a2_ecs","divergence_a2_cons", "divergence_a2_fv")
        div = create_div_a2_operator(domain,div_operator_name)
    case("divergence_ah2")
        div = create_div_ah2_operator(domain)
    case("divergence_ch_sbp21", "divergence_ch_sbp42")
        div = create_div_ch_sbp_operator(domain,div_operator_name)
    case("divergence_ah42_sbp", "divergence_ah63_sbp", "divergence_ah43_sbp")
        div = create_div_ah_sbp_operator(domain,div_operator_name)
    case default
        call parcomm_global%abort("unknown divergence operator: "//div_operator_name)
    end select
end function create_div_operator

function create_div_c2_operator(domain) result(div)

    use div_c2_mod,             only : div_c2_t
    use halo_factory_mod,       only : create_vector_halo_procedure

    type(domain_t),   intent(in)  :: domain
    type(div_c2_t)                :: div

    integer(kind=4), parameter :: halo_width=1

    !div%exch_halo = create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
    !                                                     domain%topology, halo_width, 'full')
    call create_vector_halo_procedure(div%halo_procedure,domain,1,"C_vec_default")
end function create_div_c2_operator

function create_div_c_sbp21_operator(domain) result(div)

    use div_c2_mod,             only : div_c_sbp21_t
    use exchange_factory_mod,   only : create_symmetric_halo_vec_exchange_C

    type(domain_t),   intent(in)  :: domain
    type(div_c_sbp21_t)           :: div

    integer(kind=4), parameter :: halo_width=1

    div%exch_halo = create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
                                                         domain%topology, halo_width, 'full')
end function create_div_c_sbp21_operator

function create_div_c_sbp42_operator(domain) result(div)

    use div_c_sbp42_mod,        only : div_c_sbp42_t
    use exchange_factory_mod,   only : create_symmetric_halo_vec_exchange_C
    use sbp_factory_mod,        only : create_sbp_operator
    use grid_field_factory_mod, only : create_grid_field

    type(domain_t),   intent(in)  :: domain
    type(div_c_sbp42_t)           :: div

    integer(kind=4), parameter :: halo_width=3

    call create_grid_field(div%Gu,halo_width+1,0,domain%mesh_u)
    call create_grid_field(div%Gv,halo_width+1,0,domain%mesh_v)
    call create_grid_field(div%Dx,0,0,domain%mesh_p)
    call create_grid_field(div%Dy,0,0,domain%mesh_p)

    div%exch_halo = create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
                                                         domain%topology, halo_width, 'full')
    div%sbp_diff = create_sbp_operator("D42_staggered_i2c")
end function create_div_c_sbp42_operator

function create_div_a2_operator(domain, div_operator_name) result(div)

    use div_a2_mod,       only : div_a2_t
    use halo_factory_mod, only : create_vector_halo_procedure

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: div_operator_name
    type(div_a2_t) :: div

    integer(kind=4), parameter :: ecs_halo_width=2, default_halo_width=1

    if(div_operator_name == "divergence_a2_ecs") then
        call create_vector_halo_procedure(div%halo_procedure,domain,ecs_halo_width,"ecs_A_vec")
        div%subtype="default"
    else if(div_operator_name == "divergence_a2_cons") then
        call create_vector_halo_procedure(div%halo_procedure,domain,default_halo_width,"A_vec_default")
        div%subtype="cons"
    else if(div_operator_name == "divergence_a2_fv") then
        call create_vector_halo_procedure(div%halo_procedure,domain,default_halo_width,"A_vec_default")
        div%subtype="fv"
    else
        call parcomm_global%abort("div_factory_mod, create_div_a2_operator, "// &
                                  "unknown divergence_a2 operator subtype: "// div_operator_name)
    end if

end function create_div_a2_operator

function create_div_ah2_operator(domain) result(div)

    use div_ah2_mod,          only : div_ah2_t
    use exchange_factory_mod, only : create_xy_points_halo_exchange

    type(domain_t),   intent(in)  :: domain
    type(div_ah2_t)               :: div

    integer(kind=4), parameter :: halo_width=2

    div%exch_halo = create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                                 domain%topology,  halo_width, 'full')

end function create_div_ah2_operator

function create_div_ah_sbp_operator(domain, div_operator_name) result(div)

    use div_ah_sbp_mod,       only : div_ah_sbp_t
    use exchange_factory_mod, only : create_xy_points_halo_exchange
    use halo_factory_mod,     only : create_halo_procedure
    use sbp_factory_mod,      only : create_sbp_operator

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: div_operator_name
    type(div_ah_sbp_t)            :: div

    integer(kind=4)            :: halo_width_interior
    integer(kind=4), parameter :: halo_width_edges=1

    select case(div_operator_name)
    case ("divergence_ah42_sbp")
        halo_width_interior = 3
        div%sbp_op = create_sbp_operator("d42")
    case ("divergence_ah63_sbp")
        halo_width_interior = 5
        div%sbp_op = create_sbp_operator("d63")
    case ("divergence_ah43_sbp")
        halo_width_interior = 5
        div%sbp_op = create_sbp_operator("d43")
    case default
        call parcomm_global%abort("div_factory_mod, create_div_ah_sbp_operator"// &
                                  " - unknown SBP operator: "//div_operator_name)
    end select

    div%exch_uv_interior =  &
                    create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                                 domain%topology,  halo_width_interior, 'full')
    call create_halo_procedure(div%sync_edges,domain,1,"Ah_scalar_sync")

end function create_div_ah_sbp_operator

function create_div_ch_sbp_operator(domain, div_operator_name) result(div)

    use div_ch_sbp_mod,         only : div_ch_sbp_t
    use exchange_factory_mod,   only : create_symmetric_halo_vec_exchange_Ch
    use halo_factory_mod,       only : create_halo_procedure
    use sbp_factory_mod,        only : create_sbp_operator
    use grid_field_factory_mod, only : create_grid_field

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: div_operator_name
    type(div_ch_sbp_t)            :: div

    integer(kind=4)            :: halo_width_interior
    integer(kind=4), parameter :: halo_width_edges=1

    select case(div_operator_name)
    case ("divergence_ch_sbp21")
        halo_width_interior = 2
        div%sbp_op = create_sbp_operator("D21_staggered_c2i")
        call create_grid_field(div%Gu, halo_width_interior+1, 0, domain%mesh_y)
        call create_grid_field(div%Gv, halo_width_interior+1, 0, domain%mesh_x)
    case ("divergence_ch_sbp42")
        halo_width_interior = 3
        div%sbp_op = create_sbp_operator("D42_staggered_c2i")
        call create_grid_field(div%Gu, halo_width_interior+1, 0, domain%mesh_y)
        call create_grid_field(div%Gv, halo_width_interior+1, 0, domain%mesh_x)
    case default
        call parcomm_global%abort("div_factory_mod, create_div_ch_sbp_operator"// &
                                  " - unknown SBP operator: "//div_operator_name)
    end select

    div%exch_uv =  &
            create_symmetric_halo_vec_exchange_Ch(domain%partition, domain%parcomm, &
                                        domain%topology,  halo_width_interior, 'full')
    call create_halo_procedure(div%sync_edges,domain,1,"Ah_scalar_sync")

end function create_div_ch_sbp_operator

end module div_factory_mod
