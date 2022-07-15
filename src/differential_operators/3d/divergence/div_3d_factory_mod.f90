module div_3d_factory_mod

use abstract_div_3d_mod, only : div_3d_operator_t
use domain_mod,           only : domain_t
use parcomm_mod,          only : parcomm_global

implicit none

contains

subroutine create_div_3d_operator(div_3d, domain, div_name, diff_eta_name)

    class(div_3d_operator_t), allocatable, intent(out) :: div_3d
    type(domain_t),                         intent(in)  :: domain
    character(len=*),                       intent(in)  :: div_name
    character(len=*),             optional, intent(in)  :: diff_eta_name
    !if diff_eta_name present, then div_name means name of horizontal div operator

    if (present(diff_eta_name)) then
        call create_div_3d_hor_vert(div_3d, domain, div_name, diff_eta_name)
    else
        select case(div_name)
        case default
            call parcomm_global%abort("unknown 3d divient operator: "//div_name)
        end select
    end if

end subroutine create_div_3d_operator

subroutine create_div_3d_hor_vert(div_3d, domain, horizontal_div_name, diff_eta_name)

    use div_3d_hor_vert_mod,           only : div_3d_hor_vert_t
    use div_factory_mod,               only : create_div_operator
    use vertical_operator_factory_mod, only : create_vertical_operator
    use grid_field_factory_mod,        only : create_grid_field

    class(div_3d_operator_t), allocatable, intent(out) :: div_3d
    type(domain_t),                         intent(in)  :: domain
    character(len=*),                       intent(in)  :: horizontal_div_name, &
                                                           diff_eta_name
    type(div_3d_hor_vert_t), allocatable :: div_3d_hor_vert

    integer(kind=4) :: halo_width

    !WORKAROUND
    halo_width = 5

    allocate(div_3d_hor_vert)

    div_3d_hor_vert%div_uv_op = create_div_operator(domain, horizontal_div_name)
    call create_vertical_operator(div_3d_hor_vert%diff_eta_op, diff_eta_name)

    call create_grid_field(div_3d_hor_vert%Jw, halo_width, 0, domain%mesh_w)
    call create_grid_field(div_3d_hor_vert%diff_eta, 0, 0, domain%mesh_p)

    call move_alloc(div_3d_hor_vert, div_3d)

end subroutine create_div_3d_hor_vert

end module div_3d_factory_mod
