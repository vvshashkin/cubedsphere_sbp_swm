!Module to interpolate values to virtual points beyond face edge
!Edge interpolations are treated as one following case:
!    __________
!   / source  /
!  / /^b     /
! / +-->a   /
!/_________/     ^z
!|         |     |  /^y
!|  ^b     |     | /
!|  | targ |     |/
!|  +-->a  |     +------->x
!|_________|

!Edge and corners numeration:
!4__________3
!| edge2    |
!|e        e|
!|d        d|
!|g        g|
!|3        4|
!1----------2
!   edge1

module ecs_halo_xy_mod

use halo_mod,          only : halo_t
use exchange_halo_mod, only : exchange_t
use parcomm_mod,       only : parcomm_global

implicit none

type, extends(halo_t) :: ecs_halo_xy_t

    integer(kind=4)                       :: ts, te
    class(exchange_t),        allocatable :: exch_halo
    type(ecs_tile_halo_xy_t), allocatable :: tile(:)
    logical :: is_z_interfaces

    contains

    procedure :: get_halo_scalar => get_ecs_halo

end type

type ecs_tile_halo_xy_t
  integer(kind=4) n          !number of real grid points along cubed-sphere face edge
  integer(kind=4) halo_width !number of rows in halo-zone
  integer(kind=4) wsx, wex     !starting and ending index of interpolation weights x-edges
  integer(kind=4) wsy, wey     !starting and ending index of interpolation weights y-edges
  logical lhalo(4), lcorn(4) !flags to perform haloprocedures at specific edge and corner
  real(kind=8)   , allocatable :: wx(:,:,:) !interpolation weights for edges along x-axis
  integer(kind=4), allocatable :: indx(:,:) !interpolation stencil base point index for edges along x-axis
                                           !   s-----s-t----s------s
                                           !   base--^
  real(kind=8)   , allocatable :: wy(:,:,:) !interpolation weights, edges along y-axis
  integer(kind=4), allocatable :: indy(:,:) !interpolation stencil base point index, edgws along x-axis

  contains
  procedure, public :: interp  => ext_ecs_tile_halo

end type ecs_tile_halo_xy_t

contains


subroutine get_ecs_halo(this,f,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(ecs_halo_xy_t),     intent(inout) :: this
    class(grid_field_t),      intent(inout) :: f
    type(domain_t),           intent(in)    :: domain
    integer(kind=4),          intent(in)    :: halo_width

    integer(kind=4) t

    call this%exch_halo%do(f, domain%parcomm)

    do t=this%ts,this%te
        if(this%is_z_interfaces) then
            call this%tile(t)%interp(f%tile(t),domain%partition%tiles_xyz%tile(t),halo_width)
        else
            call this%tile(t)%interp(f%tile(t),domain%partition%tiles_xy%tile(t),halo_width)
        end if
    end do
end subroutine get_ecs_halo

subroutine ext_ecs_tile_halo(this, f, tile, halo_width)
!interpolate source face values at halo zones to target face virtual points
use grid_field_mod, only : tile_field_t
use tile_mod,       only : tile_t

class(ecs_tile_halo_xy_t), intent(in)    :: this
type(tile_field_t),        intent(inout) :: f
type(tile_t),              intent(in)    :: tile
integer(kind=4),           intent(in)    :: halo_width

!local
integer(kind=4) :: i, j, k, n, i1, i2
real(kind=8)    :: zbufx(f%is:f%ie)
real(kind=8)    :: zbufy(f%js:f%je)

!set corner diagonals & remove discontinuity between adjacent edges
n = this%n
if(this%lcorn(1)) then
    do k=tile%ks, tile%ke
        f%p(1,1,k) = (f%p(1,1,k)+f%p(1,0,k)+f%p(0,1,k)) / 3.0_8
        do j=1, halo_width
            f%p(1-j,1-j,k) = 0.5_8*(f%p(-j,1,k)+f%p(1,-j,k))
            f%p(-j,1,k) = f%p(1-j,1-j,k)
            f%p(1,-j,k) = f%p(1-j,1-j,k)
        end do
    end do
end if
if(this%lcorn(2)) then
    do k=tile%ks, tile%ke
        f%p(n,1,k) = (f%p(n,1,k)+f%p(n,0,k)+f%p(n+1,1,k)) / 3.0_8
        do j=1, halo_width
            f%p(n+j,1-j,k) = 0.5_8*(f%p(n+j+1,1,k)+f%p(n,-j,k))
            f%p(n+j+1,1,k) = f%p(n+j,1-j,k)
            f%p(n,-j,k)    = f%p(n+j,1-j,k)
        end do
    end do
end if
if(this%lcorn(3)) then
    do k=tile%ks, tile%ke
        f%p(n,n,k) = (f%p(n,n,k)+f%p(n,n+1,k)+f%p(n+1,n,k)) / 3.0_8
        do j=1, halo_width
            f%p(n+j,n+j,k) = 0.5_8*(f%p(n+j+1,n,k)+f%p(n,n+j+1,k))
            f%p(n+j+1,n,k) = f%p(n+j,n+j,k)
            f%p(n,n+j+1,k) = f%p(n+j,n+j,k)
        end do
    end do
end if
if(this%lcorn(4)) then
    do k=tile%ks, tile%ke
        f%p(1,n,k) = (f%p(1,n,k)+f%p(0,n,k)+f%p(1,n+1,k)) / 3.0_8
        do j=1, halo_width
            f%p(1-j,n+j,k) = 0.5_8*(f%p(1,n+j+1,k)+f%p(-j,n,k))
            f%p(1,n+j+1,k) = f%p(1-j,n+j,k)
            f%p(-j,n,k)    = f%p(1-j,n+j,k)
        end do
    end do
end if
!Edges
if(this%lhalo(1)) then
    do k=tile%ks,tile%ke
        do i = max(tile%is-halo_width,2),min(tile%ie+halo_width,n-1)
            f%p(i,1,k) = 0.5_8*(f%p(i,1,k)+f%p(i,0,k))
        end do
        do j=1, halo_width
            zbufx(f%is:f%ie) = f%p(f%is:f%ie,-j,k)
            i1 = max(tile%is-halo_width,2-j)
            i2 = min(tile%ie+halo_width,n+j-1)
            call ecs_ext_halo_1e(zbufx,f%is,f%ie, i1, i2,this%wx(-1:2,this%wsx:this%wex,j),&
                                 this%wsx,this%wex,this%indx(this%wsx:this%wex,j))
            f%p(i1:i2,1-j,k) = zbufx(i1:i2)
        end do
    end do
end if
if(this%lhalo(2)) then
    do k=tile%ks,tile%ke
        do i = max(tile%is-halo_width,2),min(tile%ie+halo_width,n-1)
            f%p(i,n,k) = 0.5_8*(f%p(i,n,k)+f%p(i,n+1,k))
        end do
        do j=1, halo_width
            zbufx(f%is:f%ie) = f%p(f%is:f%ie,n+j+1,k)
            i1 = max(tile%is-halo_width,2-j)
            i2 = min(tile%ie+halo_width,n+j-1)
            call ecs_ext_halo_1e(zbufx,f%is,f%ie, i1, i2,this%wx(-1:2,this%wsx:this%wex,j),&
                                 this%wsx,this%wex,this%indx(this%wsx:this%wex,j))
            f%p(i1:i2,n+j,k) = zbufx(i1:i2)
        end do
    end do
end if
if(this%lhalo(3)) then
    do k=tile%ks,tile%ke
        do j = max(tile%js-halo_width,2),min(tile%je+halo_width,n-1)
            f%p(1,j,k) = 0.5_8*(f%p(1,j,k)+f%p(0,j,k))
        end do
        do i=1, halo_width
            do j=f%js,f%je
                zbufy(j) = f%p(-i,j,k)
            end do
            i1 = max(tile%js-halo_width,2-i)
            i2 = min(tile%je+halo_width,n+i-1)
            call ecs_ext_halo_1e(zbufy,f%js,f%je, i1, i2,this%wy(-1:2,this%wsy:this%wey,i),&
                                 this%wsy,this%wey,this%indy(this%wsy:this%wey,i))
            do j=i1,i2
                f%p(1-i,j,k) = zbufy(j)
            end do
        end do
    end do
end if
if(this%lhalo(4)) then
    do k=tile%ks,tile%ke
        do j = max(tile%js-halo_width,2),min(tile%je+halo_width,n-1)
            f%p(n,j,k) = 0.5_8*(f%p(n,j,k)+f%p(n+1,j,k))
        end do
        do i=1, halo_width
            do j=f%js,f%je
                zbufy(j) = f%p(n+i+1,j,k)
            end do
            i1 = max(tile%js-halo_width,2-i)
            i2 = min(tile%je+halo_width,n+i-1)
            call ecs_ext_halo_1e(zbufy,f%js,f%je, i1, i2,this%wy(-1:2,this%wsy:this%wey,i),&
                                 this%wsy,this%wey,this%indy(this%wsy:this%wey,i))
            do j=i1,i2
                f%p(n+i,j,k) = zbufy(j)
            end do
        end do
    end do
end if
end subroutine ext_ecs_tile_halo

subroutine ecs_ext_halo_1e(zf,i1v,i2v,i1,i2,w,ws,we,ind)
!input
integer(kind=4), intent(in) :: i1v,i2v,i1,i2
!in-output:
real(kind=8), intent(inout) :: zf(i1v:i2v)!input: source face values, output: interpolated target face values
!input
integer(kind=4), intent(in) :: ws,we
real(kind=8),    intent(in) :: w(-1:2,ws:we)
integer(kind=4), intent(in) :: ind(ws:we)
!locals
real(kind=8) zh(i1:i2)!buffer for interpolated values
integer i, ii

do i = i1,i2
    ii = ind(i)
    zh(i) = sum(w(:,i)*zf(ii-1:ii+2))
end do
zf(i1:i2) = zh(i1:i2)

end subroutine ecs_ext_halo_1e

end module ecs_halo_xy_mod
