module laplace_ch_halo_mod

use abstract_laplace_mod,   only: laplace_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use mesh_mod,               only : tile_mesh_t
use domain_mod,             only : domain_t
use halo_mod,               only : halo_t

implicit none

type, extends(laplace_operator_t) :: laplace_ch_halo_t
    integer(kind=4) :: order, halo_width
    class(halo_t), allocatable :: halo_f, edge_sync
    contains
    procedure :: calc_laplace
end type laplace_ch_halo_t

contains

subroutine calc_laplace(this, f1, f, domain)
    class(laplace_ch_halo_t),  intent(inout) :: this
    type(domain_t),            intent(in)    :: domain
    type(grid_field_t),        intent(inout) :: f
    !output:
    type(grid_field_t),        intent(inout) :: f1

    integer(kind=4) :: t

    call this%halo_f%get_halo_scalar(f,domain,this%halo_width)

    do t=domain%mesh_xy%ts, domain%mesh_xy%te
        select case(this%order)
        case(2)
            call calc_laplace_tile_2(f1%tile(t),f%tile(t),domain%mesh_xy%tile(t),  &
                                     domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                                     domain%mesh_xy%scale)
        case(4)
            call calc_laplace_tile_4(f1%tile(t),f%tile(t),domain%mesh_xy%tile(t),  &
                                     domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                                     domain%mesh_xy%scale)
        end select
    end do

    call this%edge_sync%get_halo_scalar(f1,domain,0)

end subroutine calc_laplace

subroutine calc_laplace_tile_2(f1,f,mesh_xy,mesh_x,mesh_y,scale)
    !result
    type(tile_field_t), intent(inout) :: f1
    !input
    type(tile_field_t), intent(in)    :: f
    type(tile_mesh_t),  intent(in)    :: mesh_xy, mesh_x, mesh_y
    real(kind=8),       intent(in)    :: scale

    integer(kind=4) :: i, j, k
    real(kind=8)    :: gx(mesh_xy%is-1:mesh_xy%ie,mesh_xy%js-1:mesh_xy%je+1)
    real(kind=8)    :: gy(mesh_xy%is-1:mesh_xy%ie+1,mesh_xy%js-1:mesh_xy%je)
    real(kind=8)    :: gxt(mesh_xy%is-1:mesh_xy%ie,mesh_xy%js:mesh_xy%je)
    real(kind=8)    :: gyt(mesh_xy%is:mesh_xy%ie,mesh_xy%js-1:mesh_xy%je)
    real(kind=8)    :: gh1, gh0

    do k= mesh_xy%ks, mesh_xy%ke

        do j = mesh_xy%js-1, mesh_xy%je+1
            do i = mesh_xy%is-1, mesh_xy%ie
                gx(i,j) = (f%p(i+1,j,k)-f%p(i,j,k)) / (mesh_xy%hx*scale)
            end do
        end do
        do j = mesh_xy%js-1, mesh_xy%je
            do i = mesh_xy%is-1, mesh_xy%ie+1
                gy(i,j) = (f%p(i,j+1,k)-f%p(i,j,k)) / (mesh_xy%hy*scale)
            end do
        end do

        do j=mesh_xy%js, mesh_xy%je
            do i = mesh_xy%is-1, mesh_xy%ie
                gh0 = 0.5_8*(gy(i,j-1)+gy(i,j))*mesh_xy%Qi(2,i,j,k)*mesh_xy%J(i,j,k)
                gh1 = 0.5_8*(gy(i+1,j-1)+gy(i+1,j))*mesh_xy%Qi(2,i+1,j,k)*mesh_xy%J(i+1,j,k)
                gxt(i,j) = mesh_y%Qi(1,i,j,k)*gx(i,j) + 0.5_8*(gh0+gh1) / mesh_y%J(i,j,k)
            end do
        end do

        do j=mesh_xy%js-1, mesh_xy%je
            do i = mesh_xy%is, mesh_xy%ie
                gh0 = 0.5_8*(gx(i-1,j)+gx(i,j))*mesh_xy%Qi(2,i,j,k)*mesh_xy%J(i,j,k)
                gh1 = 0.5_8*(gx(i-1,j+1)+gx(i,j+1))*mesh_xy%Qi(2,i+1,j,k)*mesh_xy%J(i,j+1,k)
                gyt(i,j) = mesh_x%Qi(3,i,j,k)*gy(i,j) + 0.5_8*(gh0+gh1) / mesh_x%J(i,j,k)
            end do
        end do

        do j = mesh_xy%js, mesh_xy%je
            do i = mesh_xy%is, mesh_xy%ie
                f1%p(i,j,k) = (mesh_y%J(i,j,k)*gxt(i,j) - mesh_y%J(i-1,j,k)*gxt(i-1,j)  +&
                               mesh_x%J(i,j,k)*gyt(i,j) - mesh_x%J(i,j-1,k)*gyt(i,j-1)) /&
                               (mesh_xy%J(i,j,k)*mesh_xy%hx*scale)
            end do
        end do
    end do

end subroutine calc_laplace_tile_2

subroutine calc_laplace_tile_4(f1,f,mesh_xy,mesh_x,mesh_y,scale)
    !result
    type(tile_field_t), intent(inout) :: f1
    !input
    type(tile_field_t), intent(in)    :: f
    type(tile_mesh_t),  intent(in)    :: mesh_xy, mesh_x, mesh_y
    real(kind=8),       intent(in)    :: scale

    integer(kind=4) :: i, j, k
    real(kind=8)    :: gx(mesh_xy%is-2:mesh_xy%ie+1,mesh_xy%js-3:mesh_xy%je+3)
    real(kind=8)    :: gy(mesh_xy%is-3:mesh_xy%ie+3,mesh_xy%js-2:mesh_xy%je+1)
    real(kind=8)    :: gx_at_h(mesh_xy%is:mesh_xy%ie,mesh_xy%js-3:mesh_xy%je+3)
    real(kind=8)    :: gy_at_h(mesh_xy%is-3:mesh_xy%ie+3,mesh_xy%js:mesh_xy%je)
    real(kind=8)    :: gxt(mesh_xy%is-2:mesh_xy%ie+1,mesh_xy%js:mesh_xy%je)
    real(kind=8)    :: gyt(mesh_xy%is:mesh_xy%ie,mesh_xy%js-2:mesh_xy%je+1)
    real(kind=8)    :: div_x, div_y

    do k= mesh_xy%ks, mesh_xy%ke

        do j = mesh_xy%js-3, mesh_xy%je+3
            do i = mesh_xy%is-2, mesh_xy%ie+1
                gx(i,j) = (-f%p(i+2,j,k)+27.0_8*f%p(i+1,j,k)  - &
                              27.0_8*f%p(i,j,k)+f%p(i-1,j,k)) / &
                                            (24.0_8*mesh_xy%hx*scale)
            end do
        end do
        do j = mesh_xy%js-2, mesh_xy%je+1
            do i = mesh_xy%is-3, mesh_xy%ie+3
                gy(i,j) = (-f%p(i,j+2,k)+27.0_8*f%p(i,j+1,k)    - &
                                27.0_8*f%p(i,j,k)+f%p(i,j-1,k)) / &
                                            (24.0_8*mesh_xy%hy*scale)
            end do
        end do

        do j=mesh_xy%js-3, mesh_xy%je+3
            do i=mesh_xy%is, mesh_xy%ie
                gx_at_h(i,j) = (-gx(i-2,j)+9.0_8*gx(i-1,j)+                &
                                         9.0_8*gx(i,j)-gx(i+1,j))/16.0_8 * &
                                            mesh_xy%Qi(2,i,j,k)*mesh_xy%J(i,j,k)
            end do
        end do
        do j=mesh_xy%js, mesh_xy%je
            do i=mesh_xy%is-3, mesh_xy%ie+3
                gy_at_h(i,j) = (-gy(i,j-2)+9.0_8*gy(i,j-1)+                &
                                       9.0_8*gy(i,j)-gy(i,j+1))/16.0_8 *   &
                                       mesh_xy%Qi(2,i,j,k)*mesh_xy%J(i,j,k)
            end do
        end do

        do j=mesh_xy%js, mesh_xy%je
            do i = mesh_xy%is-2, mesh_xy%ie+1
                gxt(i,j) = mesh_y%Qi(1,i,j,k)*gx(i,j) + &
                           (-gy_at_h(i+2,j)+9.0_8*gy_at_h(i+1,j)+ &
                            9.0_8*gy_at_h(i,j)-gy_at_h(i-1,j)) / (16.0_8*mesh_y%J(i,j,k))
            end do
        end do

        do j=mesh_xy%js-2, mesh_xy%je+1
            do i = mesh_xy%is, mesh_xy%ie
                gyt(i,j) = mesh_x%Qi(3,i,j,k)*gy(i,j) + &
                           (-gx_at_h(i,j+2)+9.0_8*gx_at_h(i,j+1)+ &
                            9.0_8*gx_at_h(i,j)-gx_at_h(i,j-1))/ (16.0_8*mesh_x%J(i,j,k))
            end do
        end do

        do j = mesh_xy%js, mesh_xy%je
            do i = mesh_xy%is, mesh_xy%ie
                div_x = (-mesh_y%J(i+1,j,k)*gxt(i+1,j) + 27.0_8*mesh_y%J(i,j,k)*gxt(i,j) - &
                         27.0_8*mesh_y%J(i-1,j,k)*gxt(i-1,j)+mesh_y%J(i-2,j,k)*gxt(i-2,j)) / 24.0_8
                div_y = (-mesh_x%J(i,j+1,k)*gyt(i,j+1) + 27.0_8*mesh_x%J(i,j,k)*gyt(i,j) - &
                         27.0_8*mesh_x%J(i,j-1,k)*gyt(i,j-1)+mesh_x%J(i,j-2,k)*gyt(i,j-2)) / 24.0_8
                f1%p(i,j,k) = (div_x/mesh_xy%hx + div_y/mesh_xy%hy) / (mesh_xy%J(i,j,k)*scale)
            end do
        end do
    end do

end subroutine calc_laplace_tile_4

end module laplace_ch_halo_mod
