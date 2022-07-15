module ke_colocated_mod

use abstract_KE_mod, only : KE_operator_t
use grid_field_mod,  only : grid_field_t, tile_field_t
use domain_mod,      only : domain_t
use mesh_mod,        only : tile_mesh_t

implicit none

type, public, extends(KE_operator_t) :: ke_colocated_t
contains
    procedure, public :: calc_KE
end type ke_colocated_t

contains

subroutine calc_KE(this, KE, u, v, ut, vt, domain)
    class(ke_colocated_t), intent(inout) :: this
    type(domain_t),          intent(in)    :: domain
    type(grid_field_t),      intent(inout) :: u, v, ut, vt
    type(grid_field_t),      intent(inout) :: KE

    integer(kind=4) :: t

    do t = domain%mesh_p%ts, domain%mesh_p%te
        call calc_KE_on_tile(KE%tile(t), u%tile(t), v%tile(t), &
                             ut%tile(t), vt%tile(t), domain%mesh_p%tile(t))
    end do

end subroutine calc_KE

subroutine calc_KE_on_tile(KE, u, v, ut, vt, mesh_p)

    type(tile_field_t), intent(in)    :: u, v, ut, vt
    type(tile_field_t), intent(inout) :: KE
    type(tile_mesh_t),  intent(in)    :: mesh_p

    integer(kind=4) :: i, j, k
    real(kind=8)    :: u_contra, v_contra

    do k = mesh_p%ks, mesh_p%ke
        do j = mesh_p%js, mesh_p%je
            do i = mesh_p%is, mesh_p%ie
                KE%p(i,j,k) = (u%p(i,j,k)*ut%p(i,j,k) + v%p(i,j,k)*vt%p(i,j,k))/2.0_8
            end do
        end do
    end do
end subroutine calc_KE_on_tile

end module ke_colocated_mod
