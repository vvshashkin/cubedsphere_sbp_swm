module curl_factory_mod

use domain_mod,         only : domain_t
use abstract_curl_mod,  only : curl_operator_t
use parcomm_mod,        only : parcomm_global

implicit none

contains

subroutine create_curl_operator(curl, curl_operator_name, domain)

    class(curl_operator_t), allocatable, intent(out) :: curl
    character(len=*),                    intent(in)  :: curl_operator_name
    type(domain_t),                      intent(in)  :: domain

    if(curl_operator_name(1:8) == "curl_div") then
        call create_curl_operator_div_based(curl, curl_operator_name(6:), domain)
    else if(curl_operator_name == "curl_c_sbp42" .or. &
            curl_operator_name == "curl_c_sbp21") then
        call create_curl_operator_c_sbp(curl, curl_operator_name, domain)
    else
        call parcomm_global%abort("unknown curl operator name: "// curl_operator_name)
    end if
end subroutine create_curl_operator

subroutine create_curl_operator_div_based(curl, div_operator_name, domain)

    use curl_div_based_mod,     only : curl_div_based_t
    use div_factory_mod,        only : create_div_operator
    use grid_field_factory_mod, only : create_grid_field

    class(curl_operator_t), allocatable, intent(out) :: curl
    character(len=*),                    intent(in)  :: div_operator_name
    type(domain_t),                      intent(in)  :: domain

    type(curl_div_based_t), allocatable :: curl_div_based

    integer(kind=4) :: halo_width

    !This is temporary solution
    !We need a way to get information about required halo width from operator
    halo_width = 5

    allocate(curl_div_based)

    curl_div_based%div_op = create_div_operator(domain, div_operator_name)

    call create_grid_field(curl_div_based%ut, halo_width, 0, domain%mesh_u)
    call create_grid_field(curl_div_based%vt, halo_width, 0, domain%mesh_v)

    call move_alloc(curl_div_based, curl)

end subroutine create_curl_operator_div_based

subroutine create_curl_operator_c_sbp(curl_out, curl_operator_name, domain)

    use curl_c_sbp_mod,  only : curl_c_sbp_t
    use sbp_factory_mod, only : create_sbp_operator
    use exchange_factory_mod, only : create_symmetric_halo_vec_exchange_C
    use halo_factory_mod,     only : create_halo_procedure

    class(curl_operator_t), allocatable, intent(out) :: curl_out
    character(len=*),                    intent(in)  :: curl_operator_name
    type(domain_t),                      intent(in)  :: domain

    type(curl_c_sbp_t), allocatable :: curl
    integer(kind=4) :: halo_width

    allocate(curl)

    select case(curl_operator_name)
    case("curl_c_sbp21")
        curl%sbp_diff = create_sbp_operator("D21_staggered_c2i")
        halo_width = 2
    case("curl_c_sbp42")
        curl%sbp_diff = create_sbp_operator("D42_staggered_c2i")
        halo_width = 4
    case default
        call parcomm_global%abort("unknown sbp curl operator: "//curl_operator_name)
    end select
    curl%exchange = create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
                                                         domain%topology, halo_width, 'full')

    call create_halo_procedure(curl%sync_edges,domain,1,"Ah_scalar_sync")
    call move_alloc(curl, curl_out)

end subroutine create_curl_operator_c_sbp

end module curl_factory_mod
