module co2contra_3d_Cgrid_mod

use abstract_co2contra_3d_mod,      only : co2contra_3d_operator_t
use grid_field_mod,                 only : grid_field_t, tile_field_t
use mesh_mod,                       only : tile_mesh_t
use domain_mod,                     only : domain_t
use abstract_vertical_operator_mod, only : vertical_operator_t
use abstract_interpolators2d_mod,   only : interpolator2d_vec2vec_t
use abstract_interpolators2d_mod,   only : interpolator2d_vec2vec_t

implicit none

type, extends(co2contra_3d_operator_t), public :: co2contra_3d_Cgrid_t
    class(interpolator2d_vec2vec_t), allocatable :: interp_h2v
    class(interpolator2d_vec2vec_t), allocatable :: interp_v2h
    class(vertical_operator_t), allocatable      :: interp_w2p, interp_p2w
    type(grid_field_t)                           :: up, vp, wp
contains
    procedure :: transform    => transform_co2contra_3d_Cgrid
    ! procedure :: transform2co => transform_contra2co_3d_Cgrid
end type co2contra_3d_Cgrid_t

contains


subroutine transform_co2contra_3d_Cgrid(this, u_contra, v_contra, w_contra, &
                                                  u_cov, v_cov, w_cov, domain)
    class(co2contra_3d_Cgrid_t), intent(inout) :: this
    type(domain_t),              intent(in)    :: domain
    type(grid_field_t),          intent(inout) :: u_cov, v_cov, w_cov
    !output:
    type(grid_field_t),          intent(inout) :: u_contra, v_contra, w_contra

    integer(kind=4) :: t

    call this%interp_v2h%interp2d_vec2vec(this%up, this%vp, u_cov, v_cov, domain)
    call this%interp_w2p%apply(this%wp, w_cov, domain)

    do t = domain%partition%ts, domain%partition%te
        call start_co2contra_transform_at_p_tile(this%up%tile(t), this%vp%tile(t), &
                                                 this%wp%tile(t), domain%mesh_p%tile(t))
    end do

    !u_contra contains v part at u, v_contra contains u part at v
    call this%interp_h2v%interp2d_vec2vec(u_contra, v_contra, this%up, this%vp, domain)
    call this%interp_p2w%apply(w_contra, this%wp, domain)

    do t = domain%partition%ts, domain%partition%te
        call finalize_co2contra_transform_at_uvw_tile(             &
             u_contra%tile(t), v_contra%tile(t), w_contra%tile(t), &
             u_cov%tile(t),    v_cov%tile(t),    w_cov%tile(t),    &
             domain%mesh_u%tile(t), domain%mesh_v%tile(t), domain%mesh_w%tile(t))
    end do

end subroutine transform_co2contra_3d_Cgrid

subroutine start_co2contra_transform_at_p_tile(up, vp, wp, mesh_o)

    type(tile_field_t),     intent(inout) :: up, vp, wp
    type(tile_mesh_t),      intent(in)    :: mesh_o

    real(kind=8)    :: uw2v, vw2u, uv2w
    integer(kind=4) :: i, j, k

    do k = mesh_o%ks,mesh_o%ke
        do j = mesh_o%js, mesh_o%je
            do i = mesh_o%is, mesh_o%ie
                uw2v = (up%p(i,j,k)*mesh_o%Qi(2,i,j,k)+wp%p(i,j,k)*mesh_o%Qi(5,i,j,k))*mesh_o%J(i,j,k)
                vw2u = (vp%p(i,j,k)*mesh_o%Qi(2,i,j,k)+wp%p(i,j,k)*mesh_o%Qi(4,i,j,k))*mesh_o%J(i,j,k)
                uv2w = (up%p(i,j,k)*mesh_o%Qi(4,i,j,k)+vp%p(i,j,k)*mesh_o%Qi(5,i,j,k))*mesh_o%J(i,j,k)
                !Change components for (u,v) to (v,u) interpolation
                up%p(i,j,k) = vw2u
                vp%p(i,j,k) = uw2v
                wp%p(i,j,k) = uv2w
            end do
        end do
    end do
end subroutine start_co2contra_transform_at_p_tile
subroutine finalize_co2contra_transform_at_uvw_tile(u_contra, v_contra, w_contra, &
                                                    u_cov,    v_cov,    w_cov,    &
                                                    mesh_u, mesh_v, mesh_w)

    type(tile_field_t), intent(inout) :: u_contra, v_contra, w_contra
    type(tile_field_t), intent(in)    :: u_cov, v_cov, w_cov
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v, mesh_w

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
    do k = mesh_w%ks,mesh_w%ke
        do j = mesh_w%js, mesh_w%je
            do i = mesh_w%is, mesh_w%ie
                w_contra%p(i,j,k) = mesh_w%Qi(6,i,j,k)*w_cov%p(i,j,k)+w_contra%p(i,j,k)/mesh_w%J(i,j,k)
            end do
        end do
    end do
end subroutine finalize_co2contra_transform_at_uvw_tile
end module co2contra_3d_Cgrid_mod
