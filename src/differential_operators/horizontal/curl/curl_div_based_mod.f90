module curl_div_based_mod

use grid_field_mod,    only : grid_field_t, tile_field_t
use domain_mod,        only : domain_t
use abstract_div_mod,  only : div_operator_t
use abstract_curl_mod, only : curl_operator_t

implicit none

!This curl operator is based on the use of divergence operator
!Div in curvilinear coords div(u,v) = (d/dx(u*G) + d/dy(v*G))/G
!Where u,v - contravariant vector components
!Curl in curvilinear coords curl(U,V) = (d/dx(V)-d/dy(U))/G
!Here U,V - covariant vector components.
!Thus div(V/G,-U/G) = curl(U,V), or curl(\vec{u}) = div(\vec{u}^\perp)

type, public, extends(curl_operator_t) :: curl_div_based_t
    class(div_operator_t), allocatable :: div_op
    type(grid_field_t) :: ut, vt
contains
    procedure, public :: calc_curl
end type curl_div_based_t

contains

subroutine calc_curl(this, curl, u, v, domain)

    class(curl_div_based_t),  intent(inout) :: this
    type(domain_t),           intent(in)    :: domain
    type(grid_field_t),       intent(inout) :: u, v !covariant components
    type(grid_field_t),       intent(inout) :: curl

    integer(kind=4) :: t

    do t = domain%partition%ts, domain%partition%te
        call transform_vectors(u%tile(t), v%tile(t), this%ut%tile(t), this%vt%tile(t), &
                               domain%mesh_u%tile(t), domain%mesh_v%tile(t))
    end do

    call this%div_op%calc_div(curl, this%ut, this%vt, domain)

end subroutine calc_curl

subroutine transform_vectors(u, v, ut, vt, mesh_u, mesh_v)
!This procedure transforms velocity vector with contravariant coordinates (u, v)
!to the vector with covariant coordinates (V, -U) scaled by 1/G.
    use mesh_mod, only : tile_mesh_t

    type(tile_field_t), intent(in)    :: u, v
    type(tile_field_t), intent(inout) :: ut, vt !contravariant components of perp vec
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k

!This implementation works only for unstaggered case, so mesh_u==mesh_v==mesh_p
    do k = mesh_u%ks, mesh_u%ke
        do j = mesh_u%js, mesh_u%je
            do i = mesh_u%is, mesh_u%ie
                ut%p(i,j,k) =  v%p(i,j,k)/mesh_u%J(i,j,k)
                vt%p(i,j,k) = -u%p(i,j,k)/mesh_u%J(i,j,k)
            end do
        end do
    end do


end subroutine transform_vectors

end module curl_div_based_mod
