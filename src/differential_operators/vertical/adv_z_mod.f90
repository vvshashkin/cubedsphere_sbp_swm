module adv_z_mod

use abstract_adv_z_mod, only : adv_z_t
use grid_field_mod,     only : grid_field_t, tile_field_t
use mesh_mod,           only : mesh_t, tile_mesh_t

implicit none

type, extends(adv_z_t) :: adv_z_c2_t
    contains
    procedure :: calc_z_adv_tile => calc_z_adv_c2_tile
end type adv_z_c2_t

type, extends(adv_z_t) :: adv_z_c4_t
    contains
    procedure :: calc_z_adv_tile => calc_z_adv_c4_tile
end type adv_z_c4_t

type, extends(adv_z_t) :: adv_z_up4_t
    contains
    procedure :: calc_z_adv_tile => calc_z_adv_up4_tile
end type adv_z_up4_t

contains

subroutine calc_z_adv_c2_tile(this, f_tend, f, eta_dot, mesh,scale)
    class(adv_z_c2_t),  intent(in)    :: this
    type(tile_field_t), intent(in)    :: f, eta_dot
    type(tile_mesh_t),  intent(in)    :: mesh
    real(kind=8),       intent(in)    :: scale
    !output
    type(tile_field_t), intent(inout) :: f_tend

    integer(kind=4) :: is, ie, js, je, ks, ke
    integer(kind=4) :: i, j, k
    integer(kind=4) :: km1, kp1

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks, ke
        !constant extension of field above/below the the upper/lower boundary
        km1 = max(k-1,ks)
        kp1 = min(k+1,ke)
        do j=js, je
            do i=is, ie
                f_tend%p(i,j,k) =-eta_dot%p(i,j,k)*(f%p(i,j,kp1)-f%p(i,j,km1)) / (2.0_8*mesh%hz*scale)
            end do
        end do
    end do
end subroutine calc_z_adv_c2_tile

subroutine calc_z_adv_c4_tile(this, f_tend, f, eta_dot, mesh,scale)
    class(adv_z_c4_t),  intent(in)    :: this
    type(tile_field_t), intent(in)    :: f, eta_dot
    type(tile_mesh_t),  intent(in)    :: mesh
    real(kind=8),       intent(in)    :: scale
    !output
    type(tile_field_t), intent(inout) :: f_tend

    integer(kind=4) :: is, ie, js, je, ks, ke
    integer(kind=4) :: i, j, k
    integer(kind=4) :: km1, kp1, km2, kp2
    real(kind=8)    :: df

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks, ke
        !constant extension of field above/below the the upper/lower boundary
        km1 = max(k-1,ks); km2 = max(k-2,ks)
        kp1 = min(k+1,ke); kp2 = min(k+2,ke)
        do j=js, je
            do i=is, ie
                df = (-f%p(i,j,kp2)+8._8*f%p(i,j,kp1)-8._8*f%p(i,j,km1)+f%p(i,j,km2)) &
                                                         / (12.0_8*mesh%hz*scale)
                f_tend%p(i,j,k) =-eta_dot%p(i,j,k)*df
            end do
        end do
    end do
end subroutine calc_z_adv_c4_tile

subroutine calc_z_adv_up4_tile(this, f_tend, f, eta_dot, mesh,scale)
    class(adv_z_up4_t), intent(in)    :: this
    type(tile_field_t), intent(in)    :: f, eta_dot
    type(tile_mesh_t),  intent(in)    :: mesh
    real(kind=8),       intent(in)    :: scale
    !output
    type(tile_field_t), intent(inout) :: f_tend

    integer(kind=4) :: is, ie, js, je, ks, ke
    integer(kind=4) :: i, j, k
    integer(kind=4) :: km1, kp1, km2, kp2, km3, kp3
    real(kind=8)    :: df, left, right

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks, ke
        !constant extension of field above/below the the upper/lower boundary
        km1 = max(k-1,ks); km2 = max(k-2,ks); km3 = max(k-3,ks)
        kp1 = min(k+1,ke); kp2 = min(k+2,ke); kp3 = min(k+3,ke)
        do j=js, je
            do i=is, ie
                left  = 0.5_8+sign(0.5_8,eta_dot%p(i,j,k))
                right = 1.0_8-left
                df =  (left*( 3._8*f%p(i,j,kp1)+10._8*f%p(i,j,k)- &
                             18._8*f%p(i,j,km1)+6._8*f%p(i,j,km2)-f%p(i,j,km3))-&
                      right*( 3._8*f%p(i,j,km1)+10._8*f%p(i,j,k)- &
                             18._8*f%p(i,j,kp1)+6._8*f%p(i,j,kp2)-f%p(i,j,kp3))) &
                                                         / (12.0_8*mesh%hz*scale)
                f_tend%p(i,j,k) =-eta_dot%p(i,j,k)*df
            end do
        end do
    end do
end subroutine calc_z_adv_up4_tile

end module adv_z_mod
