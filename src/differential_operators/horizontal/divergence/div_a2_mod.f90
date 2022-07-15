module div_a2_mod

use domain_mod,         only : domain_t
use abstract_div_mod,   only : div_operator_t
use grid_field_mod,     only : grid_field_t, tile_field_t
use halo_mod,           only : halo_vec_t
use parcomm_mod,        only : parcomm_global

implicit none

type, public, extends(div_operator_t) :: div_a2_t
    class(halo_vec_t), allocatable :: halo_procedure
    character(:),      allocatable :: subtype
contains
    procedure, public :: calc_div => calc_div_a2
end type div_a2_t

contains

subroutine calc_div_a2(this, div, u, v, domain)
    class(div_a2_t),        intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: u, v
    !output
    type(grid_field_t),     intent(inout) :: div

    integer(kind=4), parameter :: halo_width = 1
    integer(kind=4) :: t

    call this%halo_procedure%get_halo_vector(u,v,domain,halo_width)
    select case(this%subtype)
    case("cons")
        do t = domain%partition%ts, domain%partition%te
            call calc_div_on_tile_cons(div%tile(t), u%tile(t), v%tile(t), &
                                       domain%mesh_o%tile(t), domain%partition%Nh,&
                                       domain%mesh_o%scale)
        end do
    case("fv")
        do t = domain%partition%ts, domain%partition%te
            call calc_div_on_tile_fv(div%tile(t), u%tile(t), v%tile(t),            &
                                     domain%mesh_o%tile(t), domain%mesh_x%tile(t), &
                                     domain%mesh_y%tile(t),&
                                     domain%mesh_o%scale)
        end do
    case("default")
        do t = domain%partition%ts, domain%partition%te
            call calc_div_on_tile(div%tile(t), u%tile(t), v%tile(t),  &
                                  domain%mesh_o%tile(t),&
                                  domain%mesh_o%scale)
        end do
    case default
        call parcomm_global%abort("div_a2_mod, calc_div_a2, unknown subtype: "// this%subtype)
    end select


end subroutine calc_div_a2

subroutine calc_div_on_tile(div, u, v, mesh, scale)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),     intent(inout) :: div
    type(tile_field_t),     intent(in)    :: u, v
    type(tile_mesh_t),      intent(in)    :: mesh
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: hx, mult_loc
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx
    do k = ks, ke
        do j = js, je
            do i = is, ie
                div%p(i,j,k) = (mesh%J(i+1,j,k)*u%p(i+1,j,k)-mesh%J(i-1,j,k)*u%p(i-1,j,k) +  &
                                mesh%J(i,j+1,k)*v%p(i,j+1,k)-mesh%J(i,j-1,k)*v%p(i,j-1,k))/  &
                                (2._8*mesh%J(i,j,k)*hx*scale)
            end do
        end do
    end do

end subroutine calc_div_on_tile

subroutine calc_div_on_tile_cons(div, u, v, mesh, nx, scale)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),     intent(inout) :: div
    type(tile_field_t),     intent(in)    :: u, v
    type(tile_mesh_t),      intent(in)    :: mesh
    integer(kind=4),        intent(in)    :: nx
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: hx, mult_loc
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k
    integer(kind=4) :: ip1, im1, jp1, jm1

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx
    do k = ks, ke
        do j = js, je
            jm1 = max(1,j-1)
            jp1 = min(nx,j+1)
            do i = is, ie
                im1 = max(1,i-1)
                ip1 = min(nx,i+1)
                div%p(i,j,k) = (mesh%J(ip1,j,k)*u%p(i+1,j,k)-mesh%J(im1,j,k)*u%p(i-1,j,k) +  &
                                mesh%J(i,jp1,k)*v%p(i,j+1,k)-mesh%J(i,jm1,k)*v%p(i,j-1,k))/  &
                                (2._8*mesh%J(i,j,k)*hx*scale)
            end do
        end do
    end do

end subroutine calc_div_on_tile_cons

subroutine calc_div_on_tile_fv(div, u, v, mesh_o, mesh_x, mesh_y, scale)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),     intent(inout) :: div
    type(tile_field_t),     intent(in)    :: u, v
    type(tile_mesh_t),      intent(in)    :: mesh_o, mesh_x, mesh_y
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: hx, mult_loc
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k

    is = mesh_o%is; ie = mesh_o%ie
    js = mesh_o%js; je = mesh_o%je
    ks = mesh_o%ks; ke = mesh_o%ke

    hx = mesh_o%hx
    do k = ks, ke
        do j = js, je
            do i = is, ie
                div%p(i,j,k) = (mesh_x%J(i+1,j,k)*(u%p(i+1,j,k)+u%p(i,j,k)) -     &
                                mesh_x%J(i  ,j,k)*(u%p(i  ,j,k)+u%p(i-1,j,k)) +   &
                                mesh_y%J(i,j+1,k)*(v%p(i,j+1,k)+v%p(i,j,k))-      &
                                mesh_y%J(i,j  ,k)*(v%p(i,j  ,k)+v%p(i,j-1,k))) /  &
                                (2._8*mesh_o%J(i,j,k)*hx*scale)
            end do
        end do
    end do

end subroutine calc_div_on_tile_fv

end module div_a2_mod
