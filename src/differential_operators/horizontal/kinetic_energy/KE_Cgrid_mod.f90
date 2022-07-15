module ke_Cgrid_mod

use abstract_KE_mod,              only : KE_operator_t
use grid_field_mod,               only : grid_field_t, tile_field_t
use domain_mod,                   only : domain_t
use mesh_mod,                     only : tile_mesh_t
use abstract_interpolators2d_mod, only : interpolator2d_vec2vec_t

implicit none

type, public, extends(KE_operator_t) :: ke_Cgrid_t
    type(grid_field_t)                           :: KE_u,  KE_v
    type(grid_field_t)                           :: KE_uh, KE_vh
    class(interpolator2d_vec2vec_t), allocatable :: interp_op
contains
    procedure, public :: calc_KE
end type ke_Cgrid_t

contains

subroutine calc_KE(this, KE, u, v, ut, vt, domain)
    class(KE_Cgrid_t),  intent(inout) :: this
    type(domain_t),     intent(in)    :: domain
    type(grid_field_t), intent(inout) :: u, v, ut, vt
    type(grid_field_t), intent(inout) :: KE

    integer(kind=4) :: t

    do t = domain%partition%ts, domain%partition%te
        call calc_KE_detG_uv_on_tile(this%KE_u%tile(t), u%tile(t), ut%tile(t), domain%mesh_u%tile(t))
        call calc_KE_detG_uv_on_tile(this%KE_v%tile(t), v%tile(t), vt%tile(t), domain%mesh_v%tile(t))
    end do

    call this%interp_op%interp2d_vec2vec(this%KE_uh, this%KE_vh, this%KE_u, this%KE_v, domain)

    do t = domain%partition%ts, domain%partition%te
        call calc_KE_on_tile(KE%tile(t), this%KE_uh%tile(t), this%KE_vh%tile(t), domain%mesh_p%tile(t))
    end do
end subroutine calc_KE

subroutine calc_KE_detG_uv_on_tile(KE_uv, uv, uvt, mesh)

    type(tile_field_t), intent(in)    :: uv, uvt
    type(tile_field_t), intent(inout) :: KE_uv
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    real(kind=8)    :: u_contra, v_contra

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                KE_uv%p(i,j,k) = mesh%J(i,j,k)*(uv%p(i,j,k)*uvt%p(i,j,k))/2.0_8
            end do
        end do
    end do

end subroutine calc_KE_detG_uv_on_tile
subroutine calc_KE_on_tile(KE, KE_uh, KE_vh, mesh)

    type(tile_field_t), intent(in)    :: KE_uh, KE_vh
    type(tile_field_t), intent(inout) :: KE
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    real(kind=8)    :: u_contra, v_contra

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                KE%p(i,j,k) = (KE_uh%p(i,j,k)+KE_vh%p(i,j,k))/mesh%J(i,j,k)
            end do
        end do
    end do

end subroutine calc_KE_on_tile
end module ke_Cgrid_mod
