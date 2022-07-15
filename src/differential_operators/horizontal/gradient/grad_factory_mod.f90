module grad_factory_mod

use abstract_grad_mod, only : grad_operator_t
use domain_mod,        only : domain_t
use parcomm_mod,       only : parcomm_global

implicit none

contains

function create_grad_operator(domain, grad_operator_name) result(grad)
    type(domain_t),    intent(in)  :: domain
    character(len=*),  intent(in)  :: grad_operator_name

    class(grad_operator_t), allocatable :: grad

    if(grad_operator_name == 'gradient_c2_ecs') then
        grad = create_grad_c2_ecs_operator(domain)
    ! else if(grad_operator_name == 'gradient_c2_cons') then
    !     grad = create_grad_contra_c2_cons_operator(domain)
    else if(grad_operator_name == 'gradient_c_sbp21') then
        grad = create_grad_c_sbp21_operator(domain)
    else if(grad_operator_name == 'gradient_c_sbp42') then
        grad = create_grad_c_sbp42_operator(domain)
    else if(grad_operator_name == 'gradient_a2_ecs' .or. &
            grad_operator_name == 'gradient_a2_cons') then
        grad = create_grad_a2_operator(domain,grad_operator_name)
    else if(grad_operator_name == 'gradient_ah21_sbp_ecs' .or. &
            grad_operator_name == 'gradient_ah42_sbp_ecs' .or. &
            grad_operator_name == 'gradient_ah63_sbp_ecs' .or. &
            grad_operator_name == 'gradient_ah43_sbp_ecs') then
        grad = create_grad_ah_sbp_operator(domain, grad_operator_name)
    else if (grad_operator_name == 'gradient_ch_sbp21' .or. &
             grad_operator_name == 'gradient_ch_sbp42') then
        grad = create_grad_ch_sbp_operator(domain, grad_operator_name)
    else if(grad_operator_name == "gradient_ch_ecs_halo2" .or. &
            grad_operator_name == "gradient_ch_ecs_halo4") then
        grad = create_grad_ch_halo_operator(domain, grad_operator_name)
    else
        call parcomm_global%abort("unknown gradient operator: "//grad_operator_name)
    end if
end

function create_grad_c2_ecs_operator(domain) result(grad)

    use grad_c2_ecs_mod, only : grad_c2_ecs_t
    use halo_factory_mod,       only : create_halo_procedure

    type(domain_t),   intent(in)  :: domain
    type(grad_c2_ecs_t)           :: grad

    integer(kind=4), parameter :: ecs_halo_width=2

    grad = grad_c2_ecs_t()
    call create_halo_procedure(grad%halo_procedure,domain,ecs_halo_width,"ECS_O")

end function create_grad_c2_ecs_operator

function create_grad_c_sbp21_operator(domain) result(grad)

    use grad_c2_ecs_mod,      only : grad_c_sbp21_t
    use exchange_factory_mod, only : create_o_points_halo_exchange

    type(domain_t),   intent(in) :: domain
    type(grad_c_sbp21_t)         :: grad

    integer(kind=4), parameter :: halo_width=2

    grad%exch_halo = create_o_points_halo_exchange( &
                    domain%partition, domain%parcomm, domain%topology,  halo_width, 'full')
end function create_grad_c_sbp21_operator

function create_grad_c_sbp42_operator(domain) result(grad)

    use grad_c_sbp42_mod,        only : grad_c_sbp42_t
    use exchange_factory_mod,    only : create_o_points_halo_exchange
    use grid_field_factory_mod,  only : create_grid_field

    type(domain_t),   intent(in)      :: domain
    type(grad_c_sbp42_t)              :: grad

    integer(kind=4), parameter :: halo_width=3

    grad%exch_f = create_o_points_halo_exchange( &
                    domain%partition, domain%parcomm, domain%topology,  halo_width, 'cross')

end function create_grad_c_sbp42_operator

! function create_grad_contra_c2_cons_operator(domain) result(grad)
!
!     use grad_contra_c2_ecs_mod, only : grad_contra_c2_cons_t
!     use halo_factory_mod,   only : create_halo_procedure, create_vector_halo_procedure
!
!     type(domain_t),   intent(in)  :: domain
!     type(grad_contra_c2_cons_t)   :: grad
!
!     integer(kind=4), parameter :: halo_width=1
!
!     grad = grad_contra_c2_cons_t()
!     call create_halo_procedure(grad%halo_procedure,domain,halo_width,"A_default")
!     call create_vector_halo_procedure(grad%sync_procedure,domain,0,"C_vec_default")
!
! end function create_grad_contra_c2_cons_operator

function create_grad_a2_operator(domain, grad_operator_name) result(grad)

    use grad_a2_mod,        only : grad_a2_t
    use halo_factory_mod,   only : create_halo_procedure

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: grad_operator_name
    type(grad_a2_t)               :: grad

    integer(kind=4), parameter :: ecs_halo_width=2, default_halo_width=1

    grad = grad_a2_t()
    if(grad_operator_name=="gradient_a2_ecs") then
        call create_halo_procedure(grad%halo_procedure,domain,ecs_halo_width,"ECS_O")
    else if(grad_operator_name=="gradient_a2_cons") then
        call create_halo_procedure(grad%halo_procedure,domain,default_halo_width,"A_default")
    else
        call parcomm_global%abort("grad_factory_mod, create_grad_a2_operator "//&
                                  "unknown gradient_a2 subtype: "// grad_operator_name)
    end if

end function create_grad_a2_operator

function create_grad_ah_sbp_operator(domain, grad_operator_name) result(grad)

    use grad_ah_sbp_mod,           only : grad_ah_sbp_t
    use exchange_factory_mod,      only : create_xy_points_halo_exchange
    use halo_factory_mod,          only : create_vector_halo_procedure
    use sbp_factory_mod,           only : create_sbp_operator

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: grad_operator_name
    type(grad_ah_sbp_t)           :: grad

    integer(kind=4)               :: halo_width_interior
    integer(kind=4), parameter    :: halo_width_edges=1

    select case(grad_operator_name)
    case ("gradient_ah21_sbp_ecs")
        halo_width_interior = 1
        grad%sbp_op = create_sbp_operator("d21")
    case ("gradient_ah42_sbp_ecs")
        halo_width_interior = 3
        grad%sbp_op = create_sbp_operator("d42")
    case ("gradient_ah63_sbp_ecs")
        halo_width_interior = 5
        grad%sbp_op = create_sbp_operator("d63")
    case ("gradient_ah43_sbp_ecs")
        halo_width_interior = 5
        grad%sbp_op = create_sbp_operator("d43")
    case default
        call parcomm_global%abort("grad_factory_mod, create_grad_ah_sbp_operator"// &
                                  " - unknown SBP operator: "//grad_operator_name)
    end select

    grad%exch_scalar_interior =  &
              create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                         domain%topology,  halo_width_interior, 'full')

    call create_vector_halo_procedure(grad%sync_edges,domain,0,"ecs_Ah_vec_sync_covariant")

end function create_grad_ah_sbp_operator

function create_grad_ch_sbp_operator(domain, grad_operator_name) result(grad)

    use grad_ch_sbp_mod,      only : grad_ch_sbp_t
    use exchange_factory_mod, only : create_xy_points_halo_exchange
    use sbp_factory_mod,      only : create_sbp_operator

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: grad_operator_name
    type(grad_ch_sbp_t)           :: grad

    integer(kind=4)               :: halo_width_interior

    select case(grad_operator_name)
    case ("gradient_ch_sbp21")
        halo_width_interior = 1
        grad%sbp_op = create_sbp_operator("D21_staggered_i2c")
    case ("gradient_ch_sbp42")
        halo_width_interior = 3
        grad%sbp_op = create_sbp_operator("D42_staggered_i2c")
    case default
        call parcomm_global%abort("grad_factory_mod, create_grad_ch_sbp_operator"// &
                                  " - unknown SBP operator: "//grad_operator_name)
    end select

    grad%exch_scalar_interior =  &
              create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                         domain%topology,  halo_width_interior, 'full')

end function create_grad_ch_sbp_operator

function create_grad_ch_halo_operator(domain, grad_operator_name) result(grad)

    use grad_ch_halo_mod,     only : grad_ch_halo_t
    use halo_factory_mod,     only : create_halo_procedure

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: grad_operator_name
    type(grad_ch_halo_t)          :: grad

    integer(kind=4)               :: halo_width_interior

    select case (grad_operator_name)
    case("gradient_ch_ecs_halo2")
        grad%order             = 2
        grad%input_halo_width  = 0
        grad%output_halo_width = 0
        call create_halo_procedure(grad%halo,domain,grad%input_halo_width,"ECS_xy")
    case("gradient_ch_ecs_halo4")
        grad%order             = 4
        grad%input_halo_width  = 1
        grad%output_halo_width = 0
        call create_halo_procedure(grad%halo,domain,grad%input_halo_width,"ECS_xy")
    case default
        call parcomm_global%abort("create_grad_ch_halo, unknown operator: "//&
                                  grad_operator_name)
    end select

end function create_grad_ch_halo_operator

end module grad_factory_mod
