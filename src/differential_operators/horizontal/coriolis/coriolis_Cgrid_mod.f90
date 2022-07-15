module coriolis_Cgrid_mod

use grid_field_mod,                only : grid_field_t, tile_field_t
use domain_mod,                    only : domain_t
use abstract_coriolis_mod,         only : coriolis_operator_t
use mesh_mod,                      only : tile_mesh_t
use interpolator_w2h_mod,          only : interpolator_w2h_t
use abstract_interpolators2d_mod,  only : interpolator2d_vec2vec_t

implicit none

type, public, extends(coriolis_operator_t) :: coriolis_Cgrid_t
    type(grid_field_t) :: f !coriolis parameter
    type(grid_field_t) :: Ghu, Ghv, Ghu_p, Ghv_p, curl_p, cor_u_p, cor_v_p

    class(interpolator2d_vec2vec_t), allocatable :: interp_v2h_op
    type(interpolator_w2h_t)                     :: interp_w2h_op
    class(interpolator2d_vec2vec_t), allocatable :: interp_h2v_op

contains
    procedure, public :: calc_coriolis
    procedure, public :: calc_coriolis_vec_inv
end type coriolis_Cgrid_t

contains

subroutine calc_coriolis(this, cor_u, cor_v, ut, vt, domain)
    class(coriolis_Cgrid_t), intent(inout) :: this
    type(domain_t),          intent(in)    :: domain
    type(grid_field_t),      intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),      intent(inout) :: cor_u, cor_v

    integer(kind=4) :: t

    call domain%parcomm%abort("calc_coriolis is not implemented!")

end subroutine calc_coriolis

subroutine calc_coriolis_vec_inv(this, cor_u, cor_v, hu, hv, h, curl, domain)
    class(coriolis_Cgrid_t), intent(inout) :: this
    type(domain_t),          intent(in)    :: domain
    type(grid_field_t),      intent(inout) :: hu, hv! massflux contravariant components
    type(grid_field_t),      intent(inout) :: h, curl
    type(grid_field_t),      intent(inout) :: cor_u, cor_v

    integer(kind=4) :: t

    do t = domain%partition%ts, domain%partition%te
        call calc_Ghuv_tile(this%Ghu%tile(t), hu%tile(t), domain%mesh_u%tile(t))
        call calc_Ghuv_tile(this%Ghv%tile(t), hv%tile(t), domain%mesh_v%tile(t))
    end do

    call this%interp_v2h_op%interp2d_vec2vec(this%Ghu_p, this%Ghv_p, this%Ghu, this%Ghv, domain)
    call this%interp_w2h_op%interp_w2h(this%curl_p, curl, domain)

    do t = domain%partition%ts, domain%partition%te
        call calc_coriolis_vec_inv_tile_p(this%Ghu_p%tile(t), this%Ghv_p%tile(t), &
                                          h%tile(t), this%curl_p%tile(t), this%f%tile(t), &
                                          this%cor_u_p%tile(t), this%cor_v_p%tile(t), &
                                          domain%mesh_p%tile(t))
    end do

    call this%interp_h2v_op%interp2d_vec2vec(cor_u, cor_v, this%cor_u_p, this%cor_v_p, domain)

end subroutine calc_coriolis_vec_inv
subroutine calc_Ghuv_tile(Ghuv, uv, mesh)

    type(tile_field_t), intent(in)    :: uv
    type(tile_field_t), intent(inout) :: Ghuv
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    real(kind=8)    :: u_contra, v_contra

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                Ghuv%p(i,j,k) = mesh%J(i,j,k)*uv%p(i,j,k)
            end do
        end do
    end do

end subroutine calc_Ghuv_tile
subroutine calc_coriolis_vec_inv_tile_p(Ghu_p, Ghv_p, h, curl_p, f, cor_u_p, cor_v_p, mesh_p)

    type(tile_field_t), intent(in)    :: Ghu_p, Ghv_p, curl_p, f, h
    type(tile_field_t), intent(inout) :: cor_u_p, cor_v_p
    type(tile_mesh_t),  intent(in)    :: mesh_p

    integer(kind=4) :: i, j, k

!This implementation works only for unstaggered case, so mesh_u==mesh_v==mesh_p
    do k = mesh_p%ks, mesh_p%ke
        do j = mesh_p%js, mesh_p%je
            do i = mesh_p%is, mesh_p%ie
                cor_u_p%p(i,j,k) =  (f%p(i,j,1)+curl_p%p(i,j,k))*Ghv_p%p(i,j,k)/h%p(i,j,k)
                cor_v_p%p(i,j,k) = -(f%p(i,j,1)+curl_p%p(i,j,k))*Ghu_p%p(i,j,k)/h%p(i,j,k)
            end do
        end do
    end do


end subroutine calc_coriolis_vec_inv_tile_p

end module coriolis_Cgrid_mod
