module div_c2_mod

use abstract_div_mod,   only : div_operator_t
use domain_mod,         only : domain_t
use mesh_mod,           only : mesh_t, tile_mesh_t
use grid_field_mod,     only : grid_field_t, tile_field_t
use exchange_halo_mod,  only : exchange_t
use halo_mod,           only : halo_vec_t

implicit none

type, public, extends(div_operator_t) :: div_c2_t
    class(halo_vec_t), allocatable :: halo_procedure
contains
    procedure, public :: calc_div => calc_div_c2
end type div_c2_t

type, public, extends(div_operator_t) :: div_c_sbp21_t
    class(exchange_t), allocatable :: exch_halo
    class(halo_vec_t), allocatable :: halo_procedure
contains
    procedure, public :: calc_div => calc_div_c_sbp21
end type div_c_sbp21_t

contains

subroutine calc_div_c2(this, div, u, v, domain)

    class(div_c2_t),        intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: u, v
    !out put
    type(grid_field_t),     intent(inout) :: div

    integer(kind=4) :: i, j, k, t
    real(kind=8) hx

    call this%halo_procedure%get_halo_vector(u,v,domain,1)

    do t = domain%partition%ts, domain%partition%te
        call calc_div_on_tile(div%tile(t), u%tile(t), v%tile(t),            &
                              domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                              domain%mesh_o%tile(t),domain%mesh_o%scale)
    end do

end subroutine calc_div_c2

subroutine calc_div_on_tile(div, u, v, mesh_u, mesh_v, mesh_p, scale)

    type(tile_field_t),  intent(inout) :: div
    type(tile_field_t),  intent(in)    :: u, v
    type(tile_mesh_t),   intent(in)    :: mesh_u, mesh_v, mesh_p
    real(kind=8),        intent(in)    :: scale

    real(kind=8) :: hx
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k

    is = mesh_p%is; ie = mesh_p%ie
    js = mesh_p%js; je = mesh_p%je
    ks = mesh_p%ks; ke = mesh_p%ke

    hx = mesh_p%hx

    do k = ks, ke
        do j = js, je
            do i = is, ie
                div%p(i,j,k) = (mesh_u%J(i+1,j,k)*u%p(i+1,j,k)-mesh_u%J(i,j,k)*u%p(i,j,k)  +  &
                                mesh_v%J(i,j+1,k)*v%p(i,j+1,k)-mesh_v%J(i,j,k)*v%p(i,j,k))/  &
                                (mesh_p%J(i,j,k)*hx*scale)
            end do
        end do
    end do

end subroutine calc_div_on_tile

subroutine calc_div_c_sbp21(this, div, u, v, domain)

    class(div_c_sbp21_t),   intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: u, v
    !out put
    type(grid_field_t),     intent(inout) :: div

    integer(kind=4) :: i, j, k, t
    real(kind=8) hx

    call this%exch_halo%do_vec(u,v,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_div_on_tile_sbp21(div%tile(t), u%tile(t), v%tile(t),            &
                                    domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                                    domain%mesh_o%tile(t),domain%mesh_o%scale)
    end do

end subroutine calc_div_c_sbp21

subroutine calc_div_on_tile_sbp21(div, u, v, mesh_u, mesh_v, mesh_p, scale)

    type(tile_field_t),  intent(inout) :: div
    type(tile_field_t),  intent(in)    :: u, v
    type(tile_mesh_t),   intent(in)    :: mesh_u, mesh_v, mesh_p
    real(kind=8),        intent(in)    :: scale

    real(kind=8)    :: hx
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k
    real(kind=8)    :: dx(mesh_p%is:mesh_p%ie,mesh_p%js:mesh_p%je)
    real(kind=8)    :: dy(mesh_p%is:mesh_p%ie,mesh_p%js:mesh_p%je)

    is = mesh_p%is; ie = mesh_p%ie
    js = mesh_p%js; je = mesh_p%je
    ks = mesh_p%ks; ke = mesh_p%ke

    hx = mesh_p%hx

    do k = ks, ke

        do j=js,je
            if(is == 1) dx(1,j) = (mesh_u%J(2,j,k)*u%p(2,j,k)-&
                                   mesh_u%J(1,j,k)*(0.25_8*u%p(1,j,k)+0.75_8*u%p(0,j,k)))
            if(is<=2 .and. ie >= 2) &
                        dx(2,j) = (mesh_u%J(3,j,k)*u%p(3,j,k)-mesh_u%J(2,j,k)*u%p(2,j,k)- &
                                   0.25_8*mesh_u%J(1,j,k)*(u%p(1,j,k)-u%p(0,j,k)))
            do i=max(3,is),min(ie,mesh_p%nx-2)
                dx(i,j) = (mesh_u%J(i+1,j,k)*u%p(i+1,j,k)-mesh_u%J(i,j,k)*u%p(i,j,k))
            end do
            if(is <= mesh_p%nx-1 .and. ie >= mesh_p%nx-1) then
                i = mesh_p%nx-1
                dx(i,j) = (mesh_u%J(i+1,j,k)*u%p(i+1,j,k)-mesh_u%J(i,j,k)*u%p(i,j,k)+ &
                           0.25_8*mesh_u%J(i+2,j,k)*(u%p(i+2,j,k)-u%p(i+3,j,k)))
            end if
            if(ie == mesh_p%nx) &
                dx(ie,j) = (mesh_u%J(ie+1,j,k)*(0.75_8*u%p(ie+2,j,k)+0.25_8*u%p(ie+1,j,k))- &
                            mesh_u%J(ie,j,k)*u%p(ie,j,k))
        end do

        if(js==1) then
            j=1
            do i=is,ie
                dy(i,j) = (mesh_v%J(i,j+1,k)*v%p(i,j+1,k)- &
                           mesh_v%J(i,j,k)*(0.25_8*v%p(i,j,k)+0.75_8*v%p(i,j-1,k)))
            end do
        end if
        if(js<=2.and.je>=2) then
            j=2
            do i=is,ie
                dy(i,j) = (mesh_v%J(i,j+1,k)*v%p(i,j+1,k)-mesh_v%J(i,j,k)*v%p(i,j,k)- &
                                   0.25_8*mesh_v%J(i,j-1,k)*(v%p(i,j-1,k)-v%p(i,j-2,k)))
            end do
        end if
        do j=max(js,3),min(je,mesh_p%ny-2)
            do i=is,ie
                dy(i,j) = (mesh_v%J(i,j+1,k)*v%p(i,j+1,k)-mesh_v%J(i,j,k)*v%p(i,j,k))
            end do
        end do
        if(js<=mesh_p%ny-1.and.je>=mesh_p%ny-1) then
            j = mesh_p%ny-1
            do i=is,ie
                dy(i,j) = (mesh_v%J(i,j+1,k)*v%p(i,j+1,k)-mesh_v%J(i,j,k)*v%p(i,j,k)+ &
                            0.25_8*mesh_v%J(i,j+2,k)*(v%p(i,j+2,k)-v%p(i,j+3,k)))
            end do
        end if
        if(je==mesh_p%ny) then
            j = mesh_p%ny
            do i =is,ie
                dy(i,j) = mesh_v%J(i,j+1,k)*(0.75_8*v%p(i,j+2,k)+0.25_8*v%p(i,j+1,k))- &
                          mesh_v%J(i,j,k)*v%p(i,j,k)
            end do
        end if

        do j = js,je
            do i = is, ie
                div%p(i,j,k) =  (dx(i,j)+dy(i,j))/(mesh_p%J(i,j,k)*hx*scale)
            end do
        end do
    end do

end subroutine calc_div_on_tile_sbp21

end module div_c2_mod
