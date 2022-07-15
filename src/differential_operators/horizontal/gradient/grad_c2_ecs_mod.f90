module grad_c2_ecs_mod

use abstract_grad_mod,      only : grad_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use mesh_mod,               only : mesh_t, tile_mesh_t
use halo_mod,               only : halo_t, halo_vec_t
use domain_mod,             only : domain_t
use exchange_abstract_mod,  only : exchange_t

implicit none

type, public, extends(grad_operator_t) :: grad_c2_ecs_t
    class(halo_t), allocatable :: halo_procedure
contains
    procedure, public :: calc_grad => calc_grad_c2_ecs
end type grad_c2_ecs_t

! type, public, extends(grad_operator_t) :: grad_contra_c2_cons_t
!     class(halo_t),     allocatable :: halo_procedure
!     class(halo_vec_t), allocatable :: sync_procedure
! contains
!     procedure, public :: calc_grad => calc_grad_contra_c2_cons
! end type grad_contra_c2_cons_t

type, public, extends(grad_operator_t) :: grad_c_sbp21_t
    class(exchange_t), allocatable     :: exch_halo
contains
    procedure, public :: calc_grad => calc_grad_c_sbp21
end type grad_c_sbp21_t

contains

subroutine calc_grad_c2_ecs(this, gx, gy, f, domain)
    class(grad_c2_ecs_t),  intent(inout) :: this
    type(grid_field_t),    intent(inout) :: gx
    type(grid_field_t),    intent(inout) :: gy
    type(grid_field_t),    intent(inout) :: f
    type(domain_t),        intent(in)    :: domain

    integer(kind=4) :: t
    integer(kind=4), parameter :: halo_width=1

    call this%halo_procedure%get_halo_scalar(f,domain,halo_width)

    do t = domain%partition%ts, domain%partition%te
        call calc_grad_on_tile(gx%tile(t), gy%tile(t), f%tile(t),            &
                               domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                               domain%mesh_o%tile(t),domain%mesh_o%scale)
    end do

end subroutine calc_grad_c2_ecs

! subroutine calc_grad_contra_c2_cons(this, gx, gy, f, domain)
!     class(grad_contra_c2_cons_t), intent(inout) :: this
!     type(grid_field_t),           intent(inout) :: gx
!     type(grid_field_t),           intent(inout) :: gy
!     type(grid_field_t),           intent(inout) :: f
!     type(domain_t),               intent(in)    :: domain
!
!     integer(kind=4) :: t
!     integer(kind=4), parameter :: halo_width=1
!
!     call this%halo_procedure%get_halo_scalar(f,domain,halo_width)
!
!     do t = domain%partition%ts, domain%partition%te
!         call calc_grad_on_tile_cons(gx%tile(t), gy%tile(t), f%tile(t),            &
!                                domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
!                                domain%mesh_o%tile(t),domain%mesh_o%scale)
!     end do
!
!     call this%sync_procedure%get_halo_vector(gx,gy,domain,0)
!
! end subroutine calc_grad_contra_c2_cons

subroutine calc_grad_c_sbp21(this, gx, gy, f, domain)
    class(grad_c_sbp21_t), intent(inout) :: this
    type(grid_field_t),    intent(inout) :: gx
    type(grid_field_t),    intent(inout) :: gy
    type(grid_field_t),    intent(inout) :: f
    type(domain_t),        intent(in)    :: domain

    integer(kind=4) :: t
    integer(kind=4), parameter :: halo_width=1

    call this%exch_halo%do(f,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_grad_on_tile_sbp21(gx%tile(t), gy%tile(t), f%tile(t),            &
                                     domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                                     domain%mesh_o%tile(t),domain%mesh_o%scale)
    end do

end subroutine calc_grad_c_sbp21

subroutine calc_grad_on_tile(gx, gy, f, mesh_x, mesh_y, mesh_o, scale)

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(in)    :: f
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_o
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: hx, mult_loc
    integer(kind=4) :: ks, ke
    integer(kind=4) :: jsu, jeu, isu, ieu
    integer(kind=4) :: jsv, jev, isv, iev
    integer(kind=4) :: jsx, jex, isx, iex
    integer(kind=4) :: jsy, jey, isy, iey
    integer(kind=4) :: i, j, k
    real(kind=8), allocatable :: fdx(:,:), fdy(:,:)
    real(kind=8)    :: fdx_at_y, fdy_at_x

    ks = mesh_o%ks; ke = mesh_o%ke

    isu = mesh_x%is; ieu = mesh_x%ie
    jsu = mesh_x%js; jeu = mesh_x%je
    isv = mesh_y%is; iev = mesh_y%ie
    jsv = mesh_y%js; jev = mesh_y%je

    isx = min(isu,isv);   iex = max(ieu,iev+1)
    jsx = min(jsu,jsv-1); jex = max(jeu,jev)
    allocate(fdx(isx:iex,jsx:jex))
    isy = min(isu-1,isv); iey = max(ieu,iev)
    jsy = min(jsu,jsv);   jey = max(jeu+1,jev)
    allocate(fdy(isy:iey,jsy:jey))

    hx = mesh_o%hx

    do k = ks, ke

        do j= jsx, jex
            do i= isx, iex
                gx%p(i,j,k) = (f%p(i,j,k)-f%p(i-1,j,k))/(hx*scale)
            end do
        end do

        do j= jsy, jey
            do i= isy, iey
                gy%p(i,j,k) = (f%p(i,j,k)-f%p(i,j-1,k))/(hx*scale)
            end do
        end do

    end do

end subroutine calc_grad_on_tile

! subroutine calc_grad_on_tile_cons(gx, gy, f, mesh_x, mesh_y, mesh_o, scale)
!
!     type(tile_field_t),     intent(inout) :: gx, gy
!     type(tile_field_t),     intent(in)    :: f
!     type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_o
!     real(kind=8),           intent(in)    :: scale
!
!     real(kind=8)    :: hx, mult_loc
!     integer(kind=4) :: ks, ke
!     integer(kind=4) :: jsu, jeu, isu, ieu
!     integer(kind=4) :: jsv, jev, isv, iev
!     integer(kind=4) :: jsx, jex, isx, iex
!     integer(kind=4) :: jsy, jey, isy, iey
!     integer(kind=4) :: nx, ny
!     integer(kind=4) :: i, j, k, im1, i0, jm1, j0
!     real(kind=8), allocatable :: fdx(:,:), fdy(:,:)
!     real(kind=8), allocatable :: fdx_at_o(:,:), fdy_at_o(:,:)
!     real(kind=8)    :: fdx_at_y, fdy_at_x
!
!     nx = mesh_o%nx
!     ny = mesh_o%ny
!
!     ks = mesh_o%ks; ke = mesh_o%ke
!
!     isu = mesh_x%is; ieu = mesh_x%ie
!     jsu = mesh_x%js; jeu = mesh_x%je
!     isv = mesh_y%is; iev = mesh_y%ie
!     jsv = mesh_y%js; jev = mesh_y%je
!
!     isx = min(isu,isv);          iex = min(max(ieu,iev+1),nx+1)
!     jsx = max(min(jsu,jsv-1),1); jex = min(max(jeu,jev),ny)
!     allocate(fdx(isx:iex,jsx:jex))
!     isy = max(min(isu-1,isv),1); iey = min(max(ieu,iev),nx)
!     jsy = min(jsu,jsv);          jey = min(max(jeu+1,jev),ny+1)
!     allocate(fdy(isy:iey,jsy:jey))
!
!     allocate(fdx_at_o(isv:iev,jsx:jex))
!     allocate(fdy_at_o(isy:iey,jsu:jeu))
!
!     hx = mesh_o%hx
!
!     do k = ks, ke
!
!         do j= jsx, jex
!             do i= isx,iex
!                 fdx(i,j) = (f%p(i,j,k)-f%p(i-1,j,k))/(hx*scale)
!             end do
!         end do
!
!         do j= jsy,jey
!             do i= isy, iey
!                 fdy(i,j) = (f%p(i,j,k)-f%p(i,j-1,k))/(hx*scale)
!             end do
!         end do
!
!         do j=jsx,jex
!             do i=isv,iev
!                 fdx_at_o(i,j) = 0.5*(fdx(i,j)+fdx(i+1,j))*mesh_o%G(i,j)*mesh_o%Qi(2,i,j)
!             end do
!         end do
!
!         do j=jsu,jeu
!             do i=isy,iey
!                 fdy_at_o(i,j) = 0.5*(fdy(i,j)+fdy(i,j+1))*mesh_o%G(i,j)*mesh_o%Qi(2,i,j)
!             end do
!         end do
!
!         !transform to contravariant components
!         do j=jsu,jeu
!             do i=isu,ieu
!                 i0  = min(i,nx)
!                 im1 = max(i-1,1)
!                 fdy_at_x = 0.5_8*(fdy_at_o(i0,j)+fdy_at_o(im1,j)) / mesh_x%G(i,j)
!                 gx%p(i,j,k) = mesh_x%Qi(1,i,j)*fdx(i,j) + fdy_at_x
!             end do
!         end do
!         do j=jsv,jev
!             j0  = min(j,ny)
!             jm1 = max(j-1,1)
!             do i=isv,iev
!                 fdx_at_y = 0.5_8*(fdx_at_o(i,j0)+fdx_at_o(i,jm1)) / mesh_y%G(i,j)
!                 gy%p(i,j,k) = mesh_y%Qi(3,i,j)*fdy(i,j) + fdx_at_y
!             end do
!         end do
!     end do
!
! end subroutine calc_grad_on_tile_cons

subroutine calc_grad_on_tile_sbp21(gx, gy, f, mesh_x, mesh_y, mesh_o, scale)

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(in)    :: f
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_o
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: hx, mult_loc
    integer(kind=4) :: ks, ke
    integer(kind=4) :: is, ie, js, je
    integer(kind=4) :: nx, ny
    integer(kind=4) :: i, j, k

    nx = mesh_o%nx
    ny = mesh_o%ny

    ks = mesh_o%ks; ke = mesh_o%ke

    hx = mesh_o%hx

    do k = ks, ke

        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je
        do j= js, je
            if(is == 1) &
                gx%p(1,j,k) = (0.5_8*f%p(2,j,k)+0.5_8*f%p(1,j,k)- &
                               1.5_8*f%p(0,j,k)+0.5_8*f%p(-1,j,k))/(hx*scale)
            do i=max(is,2), min(ie,nx)
                gx%p(i,j,k) = (f%p(i,j,k)-f%p(i-1,j,k))/(hx*scale)
            end do
            if(ie == nx+1) &
                gx%p(ie,j,k) = (-0.5_8*f%p(ie-2,j,k)-0.5_8*f%p(ie-1,j,k)+&
                                 1.5_8*f%p(ie,j,k)  -0.5_8*f%p(ie+1,j,k))/(hx*scale)
        end do

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        if(js==1) then
            do i= is, ie
                gy%p(i,1,k) = (0.5_8*f%p(i,2,k)+0.5_8*f%p(i,1,k)-&
                               1.5_8*f%p(i,0,k)+0.5_8*f%p(i,-1,k))/(hx*scale)
            end do
        end if
        do j= max(js,2), min(je,ny)
            do i= is, ie
                gy%p(i,j,k) = (f%p(i,j,k)-f%p(i,j-1,k))/(hx*scale)
            end do
        end do
        if(je == ny+1) then
            do i= is, ie
                gy%p(i,je,k) = (-0.5_8*f%p(i,je-2,k)-0.5_8*f%p(i,je-1,k)+&
                                 1.5_8*f%p(i,je,k)  -0.5_8*f%p(i,je+1,k))/(hx*scale)
            end do
        end if

    end do

end subroutine calc_grad_on_tile_sbp21

end module grad_c2_ecs_mod
