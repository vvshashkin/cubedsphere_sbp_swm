module curl_c_sbp_mod

use grid_field_mod,        only : grid_field_t, tile_field_t
use domain_mod,            only : domain_t
use abstract_curl_mod,     only : curl_operator_t
use sbp_operator_mod,      only : sbp_operator_t
use exchange_abstract_mod, only : exchange_t
use halo_mod,              only : halo_t

implicit none

type, public, extends(curl_operator_t) :: curl_c_sbp_t
    class(exchange_t),     allocatable :: exchange
    class(sbp_operator_t), allocatable :: sbp_diff
    class(halo_t),         allocatable :: sync_edges
contains
    procedure, public :: calc_curl
end type curl_c_sbp_t

contains

subroutine calc_curl(this, curl, u, v, domain)

    class(curl_c_sbp_t),   intent(inout) :: this
    type(domain_t),        intent(in)    :: domain
    type(grid_field_t),    intent(inout) :: u, v !covariant components
    type(grid_field_t),    intent(inout) :: curl

    integer(kind=4) :: t

    call this%exchange%do_vec(u,v,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_curl_c_sbp_tile(curl%tile(t), u%tile(t), v%tile(t),   &
                                  this%sbp_diff, domain%mesh_x%tile(t), &
                                  domain%mesh_y%tile(t), domain%mesh_xy%tile(t), &
                                  domain%mesh_xy%scale)
    end do

    call this%sync_edges%get_halo_scalar(curl,domain,1)

end subroutine calc_curl

subroutine calc_curl_c_sbp_tile(curl, u, v, sbp_diff, mesh_x, mesh_y, mesh_xy, scale)
    use tile_mod,       only : tile_t
    use mesh_mod,       only : tile_mesh_t
    use grid_field_mod, only : grid_field_t

    type(tile_field_t),    intent(inout)    :: u, v
    class(sbp_operator_t), intent(in)    :: sbp_diff
    type(tile_mesh_t),     intent(in)    :: mesh_x, mesh_y, mesh_xy
    real(kind=8),          intent(in)    :: scale

    type(tile_field_t),    intent(inout) :: curl

    integer(kind=4) :: i, j, k
    integer(kind=4) :: is, ie, js, je, ks, ke
    real(kind=8)    :: dvx(mesh_xy%is:mesh_xy%ie,mesh_xy%js:mesh_xy%je,1)
    real(kind=8)    :: duy(mesh_xy%is:mesh_xy%ie,mesh_xy%js:mesh_xy%je,1)
    type(tile_t)    :: dxdy_tile
    real(kind=8)    :: dv, du, dvu
    integer(kind=4) :: nw

    is = mesh_xy%is; ie = mesh_xy%ie
    js = mesh_xy%js; je = mesh_xy%je
    ks = mesh_xy%ks; ke = mesh_xy%ke

    dxdy_tile = tile_t(is=is, ie=ie, js=js, je=je, ks=1, ke=1)

    do k=ks, ke
        dxdy_tile%ks = k; dxdy_tile%ke = k

        call sbp_diff%apply(dvx, dxdy_tile, dxdy_tile, mesh_xy%nx, 'x', v)

        !Edges penalty procedure
        if(is == 1) then
            do j=js, je
                dv = 0.0_8
                do i=1, size(sbp_diff%proj_operator_l)
                    dv = dv + sbp_diff%proj_operator_l(i)*(v%p(i,j,k)-v%p(1-i,j,k))
                end do
                dvx(1,j,1) = dvx(1,j,1) + 0.5_8*dv / sbp_diff%Al_out(1)
            end do
        end if
        if(ie == mesh_xy%nx) then
            do j=js, je
                dv = 0.0_8
                nw = size(sbp_diff%proj_operator_r)
                do i = 0,nw-1
                    dv = dv + sbp_diff%proj_operator_r(nw-i)*(v%p(ie+i,j,k)-v%p(ie-i-1,j,k))
                    !if(j==1) print *, "P", i, sbp_diff%proj_operator_r(nw-i)
                end do
                dvx(ie,j,1) = dvx(ie,j,1) + 0.5_8*dv / sbp_diff%Ar_out(size(sbp_diff%Ar_out))
            end do
        end if

        call sbp_diff%apply(duy, dxdy_tile, dxdy_tile, mesh_xy%ny, 'y', u)

        if(js == 1) then
            do i=is, ie
                du = 0.0_8
                do j=1, size(sbp_diff%proj_operator_l)
                    du = du + sbp_diff%proj_operator_l(j)*(u%p(i,j,k)-u%p(i,1-j,k))
                end do
                duy(i,1,1) = duy(i,1,1) + 0.5_8*du / sbp_diff%Al_out(1)
            end do
        end if
        if(je == mesh_xy%ny) then
            do i=is, ie
                du = 0.0_8
                nw = size(sbp_diff%proj_operator_r)
                do j= 0,nw-1
                    du = du + sbp_diff%proj_operator_r(nw-j)*(u%p(i,je+j,k)-u%p(i,je-j-1,k))
                end do
                duy(i,je,1) = duy(i,je,1) + 0.5_8*du / sbp_diff%Ar_out(size(sbp_diff%Ar_out))
            end do
        end if

        do j=js, je
            do i=is, ie
                curl%p(i,j,k) = (dvx(i,j,1)-duy(i,j,1)) / (mesh_xy%J(i,j,k)*mesh_xy%hx*scale)
            end do
        end do
    end do

end subroutine calc_curl_c_sbp_tile


end module curl_c_sbp_mod
