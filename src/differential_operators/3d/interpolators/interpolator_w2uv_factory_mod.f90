module interpolator_w2uv_factory_mod

use abstract_interpolators3d_mod,  only : interpolator_w2uv_t
use interpolators_w2uv_mod,        only : w2uv_colocated_t, w2uv_hor_colocated_t,&
                                          w2uv_staggered_t
use vertical_operator_factory_mod, only : create_vertical_operator
use domain_mod,                    only : domain_t
use parcomm_mod,                   only : parcomm_global


implicit none

contains

subroutine create_w2uv_interpolator(w2uv_interpolator, w2uv_interpolator_name, &
                                    w2uv_hor_part_name, w2uv_vert_part_name, domain)
    class(interpolator_w2uv_t), allocatable, intent(out) :: w2uv_interpolator
    character(len=*), intent(in) :: w2uv_interpolator_name, w2uv_hor_part_name, &
                                    w2uv_vert_part_name
    type(domain_t),   intent(in) :: domain

    if(w2uv_interpolator_name == "w2uv_colocated") then
        w2uv_interpolator = w2uv_colocated_t()
    else if(w2uv_interpolator_name == "w2uv_hor_colocated") then
        call create_hor_colocated_w2uv(w2uv_interpolator, w2uv_interpolator_name, &
                                       w2uv_vert_part_name)
    else if(w2uv_interpolator_name == "w2uv_staggered") then
        call create_staggered_w2uv(w2uv_interpolator, w2uv_interpolator_name, &
                                   w2uv_hor_part_name, w2uv_vert_part_name, domain)
    else
        call parcomm_global%abort("create_w2uv_interpolator error, unknown w2uv_interpolator_name: "//&
                                   w2uv_interpolator_name)
    end if
end subroutine create_w2uv_interpolator

subroutine create_hor_colocated_w2uv(w2uv_interpolator, w2uv_interpolator_name, &
                                     w2uv_vert_part_name)

    class(interpolator_w2uv_t), allocatable, intent(out) :: w2uv_interpolator
    character(len=*), intent(in) :: w2uv_interpolator_name, w2uv_vert_part_name

    type(w2uv_hor_colocated_t), allocatable :: w2uv_hor_colocated

    allocate(w2uv_hor_colocated)

    call create_vertical_operator(w2uv_hor_colocated%w2p_oper, w2uv_vert_part_name)

    call move_alloc(w2uv_hor_colocated, w2uv_interpolator)
end subroutine create_hor_colocated_w2uv

subroutine create_staggered_w2uv(w2uv_interpolator, w2uv_interpolator_name, &
                                 w2uv_hor_part_name, w2uv_vert_part_name, domain)

    use interpolator2d_factory_mod,   only : create_scalar2vec_interpolator2d
    use grid_field_factory_mod,       only : create_grid_field

    class(interpolator_w2uv_t), allocatable, intent(out) :: w2uv_interpolator
    character(len=*), intent(in) :: w2uv_interpolator_name, w2uv_hor_part_name, &
                                    w2uv_vert_part_name
    type(domain_t),   intent(in) :: domain

    type(w2uv_staggered_t), allocatable :: w2uv_staggered

    !WORKAROUND
    integer(kind=4) :: halo_width = 3

    allocate(w2uv_staggered)

    call create_scalar2vec_interpolator2d(w2uv_staggered%h2v_oper,  &
                                          w2uv_hor_part_name, domain)

    call create_vertical_operator(w2uv_staggered%w2p_oper,w2uv_vert_part_name)

    call create_grid_field(w2uv_staggered%wp,halo_width+1,0,domain%mesh_p)

    call move_alloc(w2uv_staggered, w2uv_interpolator)

end subroutine create_staggered_w2uv

end module interpolator_w2uv_factory_mod
