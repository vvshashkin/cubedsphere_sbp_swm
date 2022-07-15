module interpolator_uv2w_factory_mod

use abstract_interpolators3d_mod,  only : interpolator_uv2w_t
use interpolators_uv2w_mod,        only : uv2w_colocated_t, uv2w_hor_colocated_t, &
                                          uv2w_staggered_t
use vertical_operator_factory_mod, only : create_vertical_operator
use domain_mod,                    only : domain_t
use parcomm_mod,                   only : parcomm_global

implicit none

contains

subroutine create_uv2w_interpolator(uv2w_interpolator, uv2w_interpolator_name, &
                                    uv2w_hor_part_name, uv2w_vert_part_name, domain)
    class(interpolator_uv2w_t), allocatable, intent(out) :: uv2w_interpolator
    character(len=*), intent(in) :: uv2w_interpolator_name, uv2w_hor_part_name, &
                                    uv2w_vert_part_name
    type(domain_t),   intent(in) :: domain

    if(uv2w_interpolator_name == "uv2w_colocated") then
        uv2w_interpolator = uv2w_colocated_t()
    else if(uv2w_interpolator_name == "uv2w_hor_colocated") then
        call create_hor_colocated_uv2w(uv2w_interpolator, uv2w_interpolator_name, &
                                       uv2w_vert_part_name)
    else if(uv2w_interpolator_name == "uv2w_staggered") then
        call create_staggered_uv2w(uv2w_interpolator, uv2w_interpolator_name, &
                                   uv2w_hor_part_name, uv2w_vert_part_name, domain)
    else
        call parcomm_global%abort("create_uv2w_interpolator error, unknown uv2w_interpolator_name: "//&
                                   uv2w_interpolator_name)
    end if
end subroutine create_uv2w_interpolator

subroutine create_hor_colocated_uv2w(uv2w_interpolator, uv2w_interpolator_name, &
                                     uv2w_vert_part_name)

    class(interpolator_uv2w_t), allocatable, intent(out) :: uv2w_interpolator
    character(len=*), intent(in) :: uv2w_interpolator_name, uv2w_vert_part_name

    type(uv2w_hor_colocated_t), allocatable :: uv2w_hor_colocated

    allocate(uv2w_hor_colocated)

    call create_vertical_operator(uv2w_hor_colocated%p2w_oper,uv2w_vert_part_name)

    call move_alloc(uv2w_hor_colocated, uv2w_interpolator)
end subroutine create_hor_colocated_uv2w

subroutine create_staggered_uv2w(uv2w_interpolator, uv2w_interpolator_name, &
                                 uv2w_hor_part_name, uv2w_vert_part_name, domain)

    use interpolator2d_factory_mod,   only : create_vec2vec_interpolator2d
    use grid_field_factory_mod,       only : create_grid_field

    class(interpolator_uv2w_t), allocatable, intent(out) :: uv2w_interpolator
    character(len=*), intent(in) :: uv2w_interpolator_name, uv2w_hor_part_name, &
                                    uv2w_vert_part_name
    type(domain_t),   intent(in) :: domain

    type(uv2w_staggered_t), allocatable :: uv2w_staggered

    integer(kind=4) :: halo_width

    allocate(uv2w_staggered)

    !WORKAROUND
    halo_width  = 3
    call create_vec2vec_interpolator2d(uv2w_staggered%uv2p_oper, uv2w_hor_part_name, domain)

    call create_grid_field(uv2w_staggered%up,0,0,domain%mesh_p)
    call create_grid_field(uv2w_staggered%vp,0,0,domain%mesh_p)

    call create_vertical_operator(uv2w_staggered%p2w_oper,uv2w_vert_part_name)

    call move_alloc(uv2w_staggered, uv2w_interpolator)

    contains
    subroutine local_throw_error()
        call parcomm_global%abort("create_staggered_uv2w, unsupported uv2w_interpolator_name "//&
                                  uv2w_interpolator_name)
    end subroutine local_throw_error
end subroutine create_staggered_uv2w

end module interpolator_uv2w_factory_mod
