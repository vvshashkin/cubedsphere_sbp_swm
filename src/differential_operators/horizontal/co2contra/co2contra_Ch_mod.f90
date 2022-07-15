module co2contra_ch_mod

use abstract_co2contra_mod,       only : co2contra_operator_t
use grid_field_mod,               only : grid_field_t, tile_field_t
use mesh_mod,                     only : tile_mesh_t
use domain_mod,                   only : domain_t
use exchange_abstract_mod,        only : exchange_t
use parcomm_mod,                  only : parcomm_global
use halo_mod,                     only : halo_vec_t
use sbp_operator_mod,             only : sbp_operator_t
use abstract_interpolators2d_mod, only : interpolator2d_vec2vec_t

implicit none

type, extends(co2contra_operator_t), public :: co2contra_ch_sbp_t
    character(len=:), allocatable :: operator_name
    class(interpolator2d_vec2vec_t), allocatable :: interp_q2uv_op
    class(interpolator2d_vec2vec_t), allocatable :: interp_uv2q_op
    type(grid_field_t)                           :: uq, vq
    contains
        procedure :: transform    => transform_co2contra_ch_sbp
!        procedure :: transform2co => transform_contra2co_ah_c_sbp
end type co2contra_ch_sbp_t

contains

subroutine transform_co2contra_ch_sbp(this, u_contra, v_contra, u_cov, v_cov, domain)
    class(co2contra_ch_sbp_t),    intent(inout)  :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u_cov, v_cov
    !output:
    type(grid_field_t),           intent(inout) :: u_contra, v_contra

    integer(kind=4) :: t

    call this%interp_uv2q_op%interp2d_vec2vec(this%uq, this%vq, u_cov, v_cov, domain)

    do t = domain%partition%ts, domain%partition%te
        call start_co2contra_transform_at_w_tile(this%uq%tile(t), this%vq%tile(t), &
                                                 domain%mesh_o%tile(t))
    end do

    call this%interp_q2uv_op%interp2d_vec2vec(u_contra, v_contra, this%uq, this%vq, domain)

    do t = domain%partition%ts, domain%partition%te
        call finalize_co2contra_transform_at_uv_tile(u_contra%tile(t), v_contra%tile(t), &
                   u_cov%tile(t), v_cov%tile(t), &
                domain%mesh_y%tile(t), domain%mesh_x%tile(t))
    end do

end subroutine transform_co2contra_ch_sbp
subroutine start_co2contra_transform_at_w_tile(uw, vw, mesh_q)

    type(tile_field_t),     intent(inout) :: uw, vw
    type(tile_mesh_t),      intent(in)    :: mesh_q

    real(kind=8)    :: u, v
    integer(kind=4) :: i, j, k

    do k = mesh_q%ks,mesh_q%ke
        do j = mesh_q%js, mesh_q%je
            do i = mesh_q%is, mesh_q%ie
                u = uw%p(i,j,k)*mesh_q%J(i,j,k)*mesh_q%Qi(2,i,j,k)
                v = vw%p(i,j,k)*mesh_q%J(i,j,k)*mesh_q%Qi(2,i,j,k)
                !Change components for (u,v) to (v,u) interpolation
                uw%p(i,j,k) = v
                vw%p(i,j,k) = u
            end do
        end do
    end do
end subroutine start_co2contra_transform_at_w_tile
subroutine finalize_co2contra_transform_at_uv_tile(u_contra, v_contra, u_cov, v_cov, mesh_u, mesh_v)

    type(tile_field_t),     intent(inout) :: u_contra, v_contra, u_cov, v_cov
    type(tile_mesh_t),      intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k

    do k = mesh_u%ks,mesh_u%ke
        do j = mesh_u%js, mesh_u%je
            do i = mesh_u%is, mesh_u%ie
                u_contra%p(i,j,k) = mesh_u%Qi(1,i,j,k)*u_cov%p(i,j,k)+u_contra%p(i,j,k)/mesh_u%J(i,j,k)
            end do
        end do
    end do
    do k = mesh_v%ks,mesh_v%ke
        do j = mesh_v%js, mesh_v%je
            do i = mesh_v%is, mesh_v%ie
                v_contra%p(i,j,k) = mesh_v%Qi(3,i,j,k)*v_cov%p(i,j,k)+v_contra%p(i,j,k)/mesh_v%J(i,j,k)
            end do
        end do
    end do

end subroutine finalize_co2contra_transform_at_uv_tile

end module co2contra_ch_mod
