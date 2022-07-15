module grad_c_sbp42_mod

use abstract_grad_mod,      only : grad_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use mesh_mod,               only : mesh_t, tile_mesh_t
use domain_mod,             only : domain_t
use exchange_abstract_mod,  only : exchange_t

implicit none

type, public, extends(grad_operator_t) :: grad_c_sbp42_t
    class(exchange_t), allocatable  :: exch_f
contains
    procedure, public :: calc_grad => calc_grad_c_sbp42
end type grad_c_sbp42_t

!Differentiation constants
integer(kind=4), parameter :: n_edge = 4
real(kind=8), parameter :: dedge(5,4) = reshape( &
                           [-2.0_8,       3.0_8,     -1.0_8,       0.0_8,       0.0_8,       &
                            -1.0_8,       1.0_8,      0.0_8,       0.0_8,       0.0_8,       &
                             1._8/24._8, -9._8/8._8,  9._8/8._8,  -1._8/24._8,  0.0_8,       &
                            -1._8/71._8,  6._8/71._8,-83._8/71._8, 81._8/71._8,-3._8/71._8], &
                            [5,4])
integer(kind=4), parameter :: d_last_nonzero(4) = [3,2,4,5]
real(kind=8), parameter :: din(4) = [1.0_8/24.0_8, -9.0_8/8.0_8, 9.0_8/8.0_8, -1.0/24.0_8]
real(kind=8), parameter :: proj_op(3) = [15.0_8/8.0_8, -5.0_8/4.0_8, 3.0_8/8.0_8]
real(kind=8), parameter :: sbp_quad(4)= [7._8/18._8, 9._8/8._8, 1._8, 71._8/72._8]

contains

subroutine calc_grad_c_sbp42(this, gx, gy, f, domain)
    class(grad_c_sbp42_t), intent(inout) :: this
    type(grid_field_t),    intent(inout) :: gx
    type(grid_field_t),    intent(inout) :: gy
    type(grid_field_t),    intent(inout) :: f
    type(domain_t),        intent(in)    :: domain

    integer(kind=4) :: t

    call this%exch_f%do(f,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_covariant_grad_on_tile(gx%tile(t), gy%tile(t), f%tile(t),            &
                                         domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                                         domain%mesh_o%tile(t),domain%mesh_o%scale)
    end do

end subroutine calc_grad_c_sbp42

subroutine calc_covariant_grad_on_tile(gx, gy, f, mesh_x, mesh_y, mesh_o, scale)

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(in)    :: f
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_o
    real(kind=8),           intent(in)    :: scale

    integer(kind=4) :: i, i1, i2, j, j1, j2, k, ks, ke, is, ie, js, je, n
    real(kind=8)    :: dx1, dh, dhy(mesh_y%is:mesh_y%ie)

    dx1 = 1.0_8/(scale*mesh_o%hx)

    ks = mesh_o%ks; ke = mesh_o%ke

    do k = ks,ke
        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je
        do j=js,je

            do i = is,min(ie,n_edge)
                i1 = d_last_nonzero(i)
                gx%p(i,j,k) =sum(dedge(1:i1,i)*f%p(1:i1,j,k)) *dx1
            end do
            if(is == 1) then
                dh = (proj_op(1)*(f%p(1,j,k)-f%p( 0,j,k))+&
                       proj_op(2)*(f%p(2,j,k)-f%p(-1,j,k))+&
                       proj_op(3)*(f%p(3,j,k)-f%p(-2,j,k)))*dx1
                gx%p(1,j,k) = gx%p(1,j,k)+0.5_8*dh/sbp_quad(1)
            end if

            do i=max(is,n_edge+1),min(ie,mesh_o%nx+1-n_edge)
                gx%p(i,j,k) = sum(f%p(i-2:i+1,j,k)*din(1:4)) *dx1
            end do

            n = mesh_o%nx+1
            do i = max(is,n-n_edge+1),ie
                i2 = n-i+1
                i1 = d_last_nonzero(i2)
                gx%p(i,j,k) =-sum(dedge(1:i1,i2)*f%p(n-1:n-i1:-1,j,k)) *dx1
            end do
            if(ie == n) then
                dh = (proj_op(1)*(f%p(n  ,j,k)-f%p(n-1,j,k))+&
                       proj_op(2)*(f%p(n+1,j,k)-f%p(n-2,j,k))+&
                       proj_op(3)*(f%p(n+2,j,k)-f%p(n-3,j,k)))*dx1
                gx%p(n,j,k) = gx%p(n,j,k)+0.5_8*dh/sbp_quad(1)
            end if
        end do

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        do j = js,min(je,n_edge)
            j1 = d_last_nonzero(j)
            do i=is,ie
                gy%p(i,j,k) =sum(dedge(1:j1,j)*f%p(i,1:j1,k)) *dx1
            end do
        end do
        if(js==1) then
            do i=is,ie
                dh = (proj_op(1)*(f%p(i,1,k)-f%p(i, 0,k))+&
                      proj_op(2)*(f%p(i,2,k)-f%p(i,-1,k))+&
                      proj_op(3)*(f%p(i,3,k)-f%p(i,-2,k)))*dx1
                gy%p(i,1,k) = gy%p(i,1,k)+0.5_8*dh/sbp_quad(1)
            end do
        end if

        do j=max(js,n_edge+1),min(je,mesh_o%ny+1-n_edge)
            do i=is,ie
                gy%p(i,j,k) = sum(f%p(i,j-2:j+1,k)*din(1:4)) *dx1
            end do
        end do

        n = mesh_o%ny+1
        do j = max(js,n-n_edge+1),je
            j2 = n-j+1
            j1 = d_last_nonzero(j2)
            do i=is,ie
                gy%p(i,j,k) =-sum(dedge(1:j1,j2)*f%p(i,n-1:n-j1:-1,k)) *dx1
            end do
        end do
        if(je==n) then
            do i=is,ie
                dh = (proj_op(1)*(f%p(i,n  ,k)-f%p(i,n-1,k))+&
                      proj_op(2)*(f%p(i,n+1,k)-f%p(i,n-2,k))+&
                      proj_op(3)*(f%p(i,n+2,k)-f%p(i,n-3,k)))*dx1
                gy%p(i,n,k) = gy%p(i,n,k)+0.5_8*dh/sbp_quad(1)
            end do
        end if
    end do
end subroutine calc_covariant_grad_on_tile

end module grad_c_sbp42_mod
