module interpolator_uv2q_sbp_mod

use abstract_interpolators2d_mod, only : interpolator2d_vec2vec_t
use grid_field_mod,               only : grid_field_t, tile_field_t
use domain_mod,                   only : domain_t
use exchange_abstract_mod,        only : exchange_t
use mesh_mod,                     only : tile_mesh_t
use sbp_operator_mod,             only : sbp_operator_t

implicit none

type, public, extends(interpolator2d_vec2vec_t) :: interpolator_uv2q_sbp_Ch_t
    class(exchange_t),     allocatable :: exchange
    class(sbp_operator_t), allocatable :: sbp_interp_uv2q
contains
    procedure, public :: interp2d_vec2vec => interp_uv2q_Ch_sbp
end type interpolator_uv2q_sbp_Ch_t

contains

subroutine interp_uv2q_Ch_sbp(this, u, v, u_source, v_source, domain)

    class(interpolator_uv2q_sbp_Ch_t), intent(inout) :: this
    type(grid_field_t),                intent(inout) :: u_source, v_source, u, v
    type(domain_t),                    intent(in)    :: domain

    integer(kind=4) :: t

    call this%exchange%do_vec(v_source, u_source, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call interp_uv2q_tile(u%tile(t), v%tile(t), u_source%tile(t), v_source%tile(t), &
                              this%sbp_interp_uv2q, domain%mesh_o%tile(t))
    end do

end subroutine interp_uv2q_Ch_sbp

subroutine interp_uv2q_tile(u_q, v_q, u, v, sbp_interp_uv2q, mesh_q)

    use tile_mod, only : tile_t

    type(tile_field_t),     intent(inout) :: u_q, v_q
    type(tile_field_t),     intent(inout) :: u, v
    class(sbp_operator_t),  intent(in)    :: sbp_interp_uv2q
    type(tile_mesh_t),      intent(in)    :: mesh_q

    type(tile_t) :: work

    work = tile_t(is = mesh_q%is, ie = mesh_q%ie, &
                  js = mesh_q%js, je = mesh_q%je, &
                  ks = mesh_q%ks, ke = mesh_q%ke)

    call sbp_interp_uv2q%apply(u_q, work, mesh_q%ny, 'y', u)
    call sbp_interp_uv2q%apply(v_q, work, mesh_q%nx, 'x', v)

end subroutine interp_uv2q_tile

end module interpolator_uv2q_sbp_mod
