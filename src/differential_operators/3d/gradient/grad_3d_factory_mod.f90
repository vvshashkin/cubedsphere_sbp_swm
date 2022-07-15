module grad_3d_factory_mod

use abstract_grad_3d_mod, only : grad_3d_operator_t
use domain_mod,           only : domain_t
use parcomm_mod,          only : parcomm_global

implicit none

contains

subroutine create_grad_3d_operator(grad_3d, domain, grad_name, diff_eta_name)

    class(grad_3d_operator_t), allocatable, intent(out) :: grad_3d
    type(domain_t),                         intent(in)  :: domain
    character(len=*),                       intent(in)  :: grad_name
    character(len=*),             optional, intent(in)  :: diff_eta_name
    !if diff_eta_name present, then grad_name means name of horizontal grad operator

    if (present(diff_eta_name)) then
        call create_grad_3d_hor_vert(grad_3d, domain, grad_name, diff_eta_name)
    else
        select case(grad_name)
        case default
            call parcomm_global%abort("unknown 3d gradient operator: "//grad_name)
        end select
    end if

end subroutine create_grad_3d_operator

subroutine create_grad_3d_hor_vert(grad_3d, domain, horizontal_grad_name, diff_eta_name)

    use grad_3d_hor_vert_mod,          only : grad_3d_hor_vert_t
    use grad_factory_mod,              only : create_grad_operator
    use vertical_operator_factory_mod, only : create_vertical_operator

    class(grad_3d_operator_t), allocatable, intent(out) :: grad_3d
    type(domain_t),                         intent(in)  :: domain
    character(len=*),                       intent(in)  :: horizontal_grad_name, &
                                                           diff_eta_name
    type(grad_3d_hor_vert_t), allocatable :: grad_3d_hor_vert

    allocate(grad_3d_hor_vert)

    grad_3d_hor_vert%grad_xy = create_grad_operator(domain, horizontal_grad_name)
    call create_vertical_operator(grad_3d_hor_vert%diff_eta, diff_eta_name)

    call move_alloc(grad_3d_hor_vert, grad_3d)

end subroutine create_grad_3d_hor_vert

end module grad_3d_factory_mod
