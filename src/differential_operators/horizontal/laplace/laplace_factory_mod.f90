module laplace_factory_mod

use domain_mod,           only : domain_t
use abstract_laplace_mod, only : laplace_operator_t
use parcomm_mod,          only : parcomm_global

implicit none

contains

subroutine create_laplace_operator(laplace_operator,laplace_operator_name,domain)
    character(len=*), intent(in) :: laplace_operator_name
    type(domain_t),   intent(in) :: domain
    !output:
    class(laplace_operator_t), allocatable, intent(out) :: laplace_operator

    if(laplace_operator_name(1:15) == "divgrad_laplace") then
        call create_divgrad_laplace(laplace_operator,laplace_operator_name, domain)
    else if(laplace_operator_name == "laplace_ch_halo2" .or. &
            laplace_operator_name == "laplace_ch_halo4") then
        call create_laplace_ch_halo(laplace_operator,laplace_operator_name, domain)
    else if(laplace_operator_name == "laplace_ah_sbp21_narrow") then
        call create_laplace_ah_sbp21_narrow(laplace_operator,laplace_operator_name, domain)
    else if(laplace_operator_name == "laplace_ah_sbp42_narrow") then
        call create_laplace_ah_sbp42_narrow(laplace_operator,laplace_operator_name, domain)
    else if(laplace_operator_name == "laplace_ah_sbp63_narrow") then
        call create_laplace_ah_sbp63_narrow(laplace_operator, laplace_operator_name, domain)
    else
        call parcomm_global%abort("create laplace_operator, unknown laplace_operator_name:"//&
                                   laplace_operator_name)
    end if
end subroutine create_laplace_operator

subroutine create_divgrad_laplace(laplace_operator,laplace_operator_name, domain)

    use divgrad_laplace_mod,    only : divgrad_laplace_t
    use grad_factory_mod,       only : create_grad_operator
    use div_factory_mod,        only : create_div_operator
    use co2contra_factory_mod,  only : create_co2contra_operator
    use grid_field_factory_mod, only : create_grid_field

    character(len=*), intent(in) :: laplace_operator_name
    type(domain_t),   intent(in) :: domain
    !output:
    class(laplace_operator_t), allocatable, intent(out) :: laplace_operator
    !local
    type(divgrad_laplace_t), allocatable :: laplace
    character(len=:), allocatable :: grad_name, div_name, co2contra_name
    integer(kind=4)  :: halo_width

    allocate(laplace)
    select case(laplace_operator_name)
    case("divgrad_laplace_c_sbp21")
        div_name = "divergence_c_sbp21"
        co2contra_name = "co2contra_c_sbp21"
        grad_name = "gradient_c_sbp21"
    case("divgrad_laplace_c_sbp42")
        div_name = "divergence_c_sbp42"
        co2contra_name = "co2contra_c_sbp42"
        grad_name = "gradient_c_sbp42"
    case("divgrad_laplace_ch_sbp21")
        div_name = "divergence_ch_sbp21"
        co2contra_name = "co2contra_ch_sbp21"
        grad_name = "gradient_ch_sbp21"
    case("divgrad_laplace_ch_sbp42")
        div_name = "divergence_ch_sbp42"
        co2contra_name = "co2contra_ch_sbp42"
        grad_name = "gradient_ch_sbp42"
    case("divgrad_laplace_ah_sbp21")
        div_name = "divergence_ah2"
        co2contra_name = "co2contra_colocated"
        grad_name = "gradient_ah21_sbp_ecs"
    case("divgrad_laplace_ah_sbp42")
        div_name = "divergence_ah42_sbp"
        co2contra_name = "co2contra_colocated"
        grad_name = "gradient_ah42_sbp_ecs"
    case default
        call parcomm_global%abort("create_divgrad_laplace, incorrect laplace_operator_name:"//&
                                  laplace_operator_name)
    end select

    laplace%div_operator       = create_div_operator(domain, div_name)
    laplace%co2contra_operator = create_co2contra_operator(domain, co2contra_name)
    laplace%grad_operator      = create_grad_operator(domain,grad_name)

    !WORKAROUND
    halo_width = 6
    call create_grid_field(laplace%gx,  halo_width, 0, domain%mesh_u)
    call create_grid_field(laplace%gy,  halo_width, 0, domain%mesh_v)
    call create_grid_field(laplace%gxt, halo_width, 0, domain%mesh_u)
    call create_grid_field(laplace%gyt, halo_width, 0, domain%mesh_v)

    call move_alloc(laplace, laplace_operator)

end subroutine create_divgrad_laplace

subroutine create_laplace_ch_halo(laplace_operator,laplace_operator_name, domain)

    use laplace_ch_halo_mod, only : laplace_ch_halo_t
    use halo_factory_mod,    only : create_halo_procedure

    character(len=*), intent(in) :: laplace_operator_name
    type(domain_t),   intent(in) :: domain
    !output:
    class(laplace_operator_t), allocatable, intent(out) :: laplace_operator
    !local
    type(laplace_ch_halo_t), allocatable :: laplace

    allocate(laplace)
    select case(laplace_operator_name)
    case("laplace_ch_halo2")
        laplace%halo_width = 1
        laplace%order      = 2
    case("laplace_ch_halo4")
        laplace%halo_width = 3
        laplace%order      = 4
    case default
        call parcomm_global%abort("create_laplace_ch_halo, incorrect laplace_operator_name:"//&
                                  laplace_operator_name)
    end select


    call create_halo_procedure(laplace%halo_f,   domain,laplace%halo_width,"ECS_xy")
    call create_halo_procedure(laplace%edge_sync,domain,1,"Ah_scalar_sync")

    call move_alloc(laplace, laplace_operator)

    contains
    function find_position_in_str(char,str) result(pos)
        character(len=1), intent(in) :: char
        character(len=*), intent(in) :: str
        integer(kind=4)  :: pos
        integer(kind=4)  :: i
        pos = -1
        do i=1, len(str)
            if(str(i:i) == char) then
                pos = i
                return
            end if
        end do
        end
end subroutine create_laplace_ch_halo
subroutine create_laplace_ah_sbp21_narrow(laplace_operator, laplace_operator_name, domain)

    use laplace_Ah_sbp21_narrow_mod, only : laplace_Ah_sbp21_narrow_t
    use exchange_factory_mod,        only : create_xy_points_halo_exchange
    use grid_field_factory_mod,      only : create_grid_field
    use sbp_factory_mod,             only : create_sbp_operator
    use halo_factory_mod,            only : create_halo_procedure

    character(len=*), intent(in) :: laplace_operator_name
    type(domain_t),   intent(in) :: domain
    !output:
    class(laplace_operator_t), allocatable, intent(out) :: laplace_operator
    !local
    type(laplace_Ah_sbp21_narrow_t), allocatable :: laplace

    allocate(laplace)

    laplace%exchange = create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                               domain%topology,  1, 'full')

    call create_grid_field(laplace%d2f_x, 0, 0, domain%mesh_xy)
    call create_grid_field(laplace%d2f_y, 0, 0, domain%mesh_xy)

    call create_grid_field(laplace%d1f_x, 1, 0, domain%mesh_xy)
    call create_grid_field(laplace%d1f_y, 1, 0, domain%mesh_xy)

    call create_grid_field(laplace%d2f_xy, 0, 0, domain%mesh_xy)
    call create_grid_field(laplace%d2f_yx, 0, 0, domain%mesh_xy)

    laplace%sbp_d1_op = create_sbp_operator("d21")

    call create_halo_procedure(laplace%edge_sync,domain,1,"Ah_scalar_sync")

    call move_alloc(laplace, laplace_operator)

end subroutine create_laplace_ah_sbp21_narrow
subroutine create_laplace_ah_sbp42_narrow(laplace_operator, laplace_operator_name, domain)

    use laplace_Ah_sbp42_narrow_mod, only : laplace_Ah_sbp42_narrow_t
    use exchange_factory_mod,        only : create_xy_points_halo_exchange
    use grid_field_factory_mod,      only : create_grid_field
    use sbp_factory_mod,             only : create_sbp_operator
    use halo_factory_mod,            only : create_halo_procedure

    character(len=*), intent(in) :: laplace_operator_name
    type(domain_t),   intent(in) :: domain
    !output:
    class(laplace_operator_t), allocatable, intent(out) :: laplace_operator
    !local
    type(laplace_Ah_sbp42_narrow_t), allocatable :: laplace

    allocate(laplace)

    laplace%exchange = create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                               domain%topology,  6, 'full')

    call create_grid_field(laplace%d2f_x, 0, 0, domain%mesh_xy)
    call create_grid_field(laplace%d2f_y, 0, 0, domain%mesh_xy)

    call create_grid_field(laplace%d1f_x, 6, 0, domain%mesh_xy)
    call create_grid_field(laplace%d1f_y, 6, 0, domain%mesh_xy)

    call create_grid_field(laplace%d2f_xy, 0, 0, domain%mesh_xy)
    call create_grid_field(laplace%d2f_yx, 0, 0, domain%mesh_xy)

    laplace%sbp_d1_op = create_sbp_operator("d42")

    call create_halo_procedure(laplace%edge_sync,domain,1,"Ah_scalar_sync")

    call move_alloc(laplace, laplace_operator)

end subroutine create_laplace_ah_sbp42_narrow
subroutine create_laplace_ah_sbp63_narrow(laplace_operator, laplace_operator_name, domain)

    use laplace_Ah_sbp63_narrow_mod, only : laplace_Ah_sbp63_narrow_t
    use exchange_factory_mod,        only : create_xy_points_halo_exchange
    use grid_field_factory_mod,      only : create_grid_field
    use sbp_factory_mod,             only : create_sbp_operator
    use halo_factory_mod,            only : create_halo_procedure

    character(len=*), intent(in) :: laplace_operator_name
    type(domain_t),   intent(in) :: domain
    !output:
    class(laplace_operator_t), allocatable, intent(out) :: laplace_operator
    !local
    type(laplace_Ah_sbp63_narrow_t), allocatable :: laplace

    allocate(laplace)

    laplace%exchange = create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                               domain%topology,  9, 'full')

    call create_grid_field(laplace%d2f_x, 0, 0, domain%mesh_xy)
    call create_grid_field(laplace%d2f_y, 0, 0, domain%mesh_xy)

    call create_grid_field(laplace%d1f_x, 9, 0, domain%mesh_xy)
    call create_grid_field(laplace%d1f_y, 9, 0, domain%mesh_xy)

    call create_grid_field(laplace%d2f_xy, 0, 0, domain%mesh_xy)
    call create_grid_field(laplace%d2f_yx, 0, 0, domain%mesh_xy)

    laplace%sbp_d1_op = create_sbp_operator("d63")

    call create_halo_procedure(laplace%edge_sync,domain,1,"Ah_scalar_sync")

    call move_alloc(laplace, laplace_operator)

end subroutine create_laplace_ah_sbp63_narrow
end module laplace_factory_mod
