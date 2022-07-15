module co2contra_3d_h_colocated_mod

use abstract_co2contra_3d_mod,      only : co2contra_3d_operator_t
use grid_field_mod,                 only : grid_field_t, tile_field_t
use mesh_mod,                       only : tile_mesh_t
use domain_mod,                     only : domain_t
use abstract_vertical_operator_mod, only : vertical_operator_t

implicit none

type, extends(co2contra_3d_operator_t), public :: co2contra_3d_h_colocated_t
    class(vertical_operator_t), allocatable :: interp_w2p, interp_p2w
    type(grid_field_t) :: wp
contains
    procedure :: transform    => transform_co2contra_3d_h_colocated
    ! procedure :: transform2co => transform_contra2co_3d_h_colocated
end type co2contra_3d_h_colocated_t

contains

subroutine transform_co2contra_3d_h_colocated(this, u_contra, v_contra, w_contra, &
                                                    u_cov, v_cov, w_cov, domain)
    class(co2contra_3d_h_colocated_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u_cov, v_cov, w_cov
    !output:
    type(grid_field_t),           intent(inout) :: u_contra, v_contra, w_contra

    integer(kind=4) :: t

    call this%interp_w2p%apply(this%wp, w_cov, domain)

    do t = domain%partition%ts, domain%partition%te
        call start_co2contra_transform_at_p_tile(u_contra%tile(t), v_contra%tile(t), this%wp%tile(t), &
                                                 u_cov%tile(t), v_cov%tile(t), domain%mesh_p%tile(t))
    end do

    call this%interp_p2w%apply(w_contra, this%wp, domain)

    do t = domain%partition%ts, domain%partition%te
        call finalize_co2contra_transform_at_w_tile(w_contra%tile(t), w_cov%tile(t), &
                                                    domain%mesh_w%tile(t))
    end do

end subroutine transform_co2contra_3d_h_colocated

subroutine start_co2contra_transform_at_p_tile(u_contra, v_contra, wp, u_cov, v_cov, mesh)

    type(tile_field_t),     intent(inout) :: u_cov, v_cov
    type(tile_field_t),     intent(inout) :: u_contra, v_contra, wp
    type(tile_mesh_t),      intent(in)    :: mesh

    real(kind=8)    :: uw2v, vw2u, uv2w
    integer(kind=4) :: i, j, k

    do k = mesh%ks,mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                u_contra%p(i,j,k) = u_cov%p(i,j,k)*mesh%Qi(1,i,j,k) + &
                                    v_cov%p(i,j,k)*mesh%Qi(2,i,j,k) + &
                                       wp%p(i,j,k)*mesh%Qi(4,i,j,k)

               v_contra%p(i,j,k) = u_cov%p(i,j,k)*mesh%Qi(2,i,j,k) + &
                                   v_cov%p(i,j,k)*mesh%Qi(3,i,j,k) + &
                                      wp%p(i,j,k)*mesh%Qi(5,i,j,k)
                !part of the w_contra at p points that are to be interpolated to w points
                wp%p(i,j,k) = (u_cov%p(i,j,k)*mesh%Qi(4,i,j,k)+v_cov%p(i,j,k)*mesh%Qi(5,i,j,k))*mesh%J(i,j,k)
            end do
        end do
    end do
end subroutine start_co2contra_transform_at_p_tile

subroutine finalize_co2contra_transform_at_w_tile(w_contra, w_cov, mesh_w)

    type(tile_field_t), intent(in)    :: w_cov
    type(tile_field_t), intent(inout) :: w_contra
    type(tile_mesh_t),  intent(in)    :: mesh_w

    integer(kind=4) :: i, j, k

    do k = mesh_w%ks,mesh_w%ke
        do j = mesh_w%js, mesh_w%je
            do i = mesh_w%is, mesh_w%ie
                w_contra%p(i,j,k) = mesh_w%Qi(6,i,j,k)*w_cov%p(i,j,k)+w_contra%p(i,j,k)/mesh_w%J(i,j,k)
            end do
        end do
    end do

end subroutine finalize_co2contra_transform_at_w_tile

end module co2contra_3d_h_colocated_mod
