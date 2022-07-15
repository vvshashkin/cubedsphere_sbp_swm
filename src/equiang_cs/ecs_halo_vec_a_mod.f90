!Object to interpolate vector values to virtual points beyond face edge
!non-staggered (A) grid. The algorithm is as follows:
!1)transform adjacent face's vector components to this face components
!2)make halo-interpolation for u,v as for scalars (non-staggered grid!)

module ecs_halo_vec_a_mod

use halo_mod,          only : halo_t, halo_vec_t
use exchange_halo_mod, only : exchange_t
use ecs_halo_mod,      only : ecs_tile_halo_t

implicit none

type, extends(halo_vec_t) :: ecs_halo_vec_t

    integer(kind=4)                        :: ts, te
    class(exchange_t),         allocatable :: exch_halo
    type(ecs_tile_halo_vec_t), allocatable :: tile(:)

    contains

    procedure :: get_halo_vector => get_ecs_a_vector_halo

end type

type ecs_tile_halo_vec_t
  integer(kind=4) n          !number of real grid points along cubed-sphere face edge
  integer(kind=4) panel_ind
  integer(kind=4) ish, ieh, jsh, jeh
  integer(kind=4) halo_width !number of rows in halo-zone (maximum)
  integer(kind=4) corner_halo_width !number of rows needed to compute corner areas
  logical lhalo(4)   !flags to perform haloprocedures at specific edge and corner!

  !Transformation matrices source panel vect ->target panel vect
  real(kind=8)   , allocatable :: TM1(:,:,:) !lower edge
  real(kind=8)   , allocatable :: TM2(:,:,:) !upper edge
  real(kind=8)   , allocatable :: TM3(:,:,:) !left edge
  real(kind=8)   , allocatable :: TM4(:,:,:) !right edge

  class(ecs_tile_halo_t), allocatable :: scalar_halo !scalar halo to interpolate individual
                                                     !vector components after transformation
  contains
  procedure, public :: interpv  => ext_ecs_a_vector_tile_halo

end type ecs_tile_halo_vec_t

contains

subroutine get_ecs_a_vector_halo(this,u,v,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(ecs_halo_vec_t),    intent(inout) :: this
    class(grid_field_t),      intent(inout) :: u,v
    type(domain_t),           intent(in)    :: domain
    integer(kind=4),          intent(in)    :: halo_width

    integer(kind=4) t

    call this%exch_halo%do_vec(u,v, domain%parcomm)

    do t=this%ts,this%te
        call this%tile(t)%interpv(u%tile(t),v%tile(t),domain%partition%tiles_o%tile(t),halo_width)
    end do
end subroutine get_ecs_a_vector_halo

subroutine ext_ecs_a_vector_tile_halo(this, u, v, tile, halo_width)
!interpolate source face vectors at halo zones to target face virtual points
use tile_mod,       only : tile_t
use grid_field_mod, only: tile_field_t

class(ecs_tile_halo_vec_t), intent(in)    :: this
type(tile_field_t),         intent(inout) :: u, v
type(tile_t),               intent(in)    :: tile
integer(kind=4),            intent(in)    :: halo_width
!locals
integer(kind=4) i, j, k, ks, ke
real(kind=8)    zu, zv

ks = u%ks
ke = u%ke

!Apply transform matrix
if(this%lhalo(1)) then
    do k = ks, ke
        do j=1,halo_width!max(halo_width,this%halo_width)
            do i=this%ish, this%ieh
                zu = u%p(i,1-j,k)
                zv = v%p(i,1-j,k)
                u%p(i,1-j,k) = this%TM1(1,i,j)*zu + this%TM1(2,i,j)*zv
                v%p(i,1-j,k) = this%TM1(3,i,j)*zu + this%TM1(4,i,j)*zv
            end do
        end do
        do j=halo_width+1,this%corner_halo_width
            do i=this%ish, 2
                zu = u%p(i,1-j,k)
                zv = v%p(i,1-j,k)
                u%p(i,1-j,k) = this%TM1(1,i,j)*zu + this%TM1(2,i,j)*zv
                v%p(i,1-j,k) = this%TM1(3,i,j)*zu + this%TM1(4,i,j)*zv
            end do
            do i=this%n-1,this%ieh
                zu = u%p(i,1-j,k)
                zv = v%p(i,1-j,k)
                u%p(i,1-j,k) = this%TM1(1,i,j)*zu + this%TM1(2,i,j)*zv
                v%p(i,1-j,k) = this%TM1(3,i,j)*zu + this%TM1(4,i,j)*zv
            end do
        end do
    end do
end if
if(this%lhalo(2)) then
    do k = ks, ke
        do j=1,halo_width
            do i=this%ish, this%ieh
                zu = u%p(i,this%n+j,k)
                zv = v%p(i,this%n+j,k)
                u%p(i,this%n+j,k) = this%TM2(1,i,j)*zu + this%TM2(2,i,j)*zv
                v%p(i,this%n+j,k) = this%TM2(3,i,j)*zu + this%TM2(4,i,j)*zv
            end do
        end do
        do j=halo_width+1,this%corner_halo_width
            do i=this%ish, 2
                zu = u%p(i,this%n+j,k)
                zv = v%p(i,this%n+j,k)
                u%p(i,this%n+j,k) = this%TM2(1,i,j)*zu + this%TM2(2,i,j)*zv
                v%p(i,this%n+j,k) = this%TM2(3,i,j)*zu + this%TM2(4,i,j)*zv
            end do
            do i=this%n-1,this%ieh
                zu = u%p(i,this%n+j,k)
                zv = v%p(i,this%n+j,k)
                u%p(i,this%n+j,k) = this%TM2(1,i,j)*zu + this%TM2(2,i,j)*zv
                v%p(i,this%n+j,k) = this%TM2(3,i,j)*zu + this%TM2(4,i,j)*zv
            end do
        end do
    end do
end if
if(this%lhalo(3)) then
    do k = ks, ke
        do j=1,halo_width
            do i=this%jsh, this%jeh
                zu = u%p(1-j,i,k)
                zv = v%p(1-j,i,k)
                u%p(1-j,i,k) = this%TM3(1,i,j)*zu + this%TM3(2,i,j)*zv
                v%p(1-j,i,k) = this%TM3(3,i,j)*zu + this%TM3(4,i,j)*zv
            end do
        end do
        do j=halo_width+1,this%corner_halo_width
            do i=this%jsh, 2
                zu = u%p(1-j,i,k)
                zv = v%p(1-j,i,k)
                u%p(1-j,i,k) = this%TM3(1,i,j)*zu + this%TM3(2,i,j)*zv
                v%p(1-j,i,k) = this%TM3(3,i,j)*zu + this%TM3(4,i,j)*zv
            end do
            do i=this%n-1,this%jeh
                zu = u%p(1-j,i,k)
                zv = v%p(1-j,i,k)
                u%p(1-j,i,k) = this%TM3(1,i,j)*zu + this%TM3(2,i,j)*zv
                v%p(1-j,i,k) = this%TM3(3,i,j)*zu + this%TM3(4,i,j)*zv
            end do
        end do
    end do
end if
if(this%lhalo(4)) then
    do k = ks, ke
        do j=1,halo_width
            do i=this%jsh, this%jeh
                zu = u%p(this%n+j,i,k)
                zv = v%p(this%n+j,i,k)
                u%p(this%n+j,i,k) = this%TM4(1,i,j)*zu + this%TM4(2,i,j)*zv
                v%p(this%n+j,i,k) = this%TM4(3,i,j)*zu + this%TM4(4,i,j)*zv
            end do
        end do
        do j=halo_width+1,this%corner_halo_width
            do i=this%jsh, 2
                zu = u%p(this%n+j,i,k)
                zv = v%p(this%n+j,i,k)
                u%p(this%n+j,i,k) = this%TM4(1,i,j)*zu + this%TM4(2,i,j)*zv
                v%p(this%n+j,i,k) = this%TM4(3,i,j)*zu + this%TM4(4,i,j)*zv
            end do
            do i=this%n-1,this%jeh
                zu = u%p(this%n+j,i,k)
                zv = v%p(this%n+j,i,k)
                u%p(this%n+j,i,k) = this%TM4(1,i,j)*zu + this%TM4(2,i,j)*zv
                v%p(this%n+j,i,k) = this%TM4(3,i,j)*zu + this%TM4(4,i,j)*zv
            end do
        end do
    end do
end if

call this%scalar_halo%interp(u,tile,halo_width)
call this%scalar_halo%interp(v,tile,halo_width)

end subroutine ext_ecs_a_vector_tile_halo

end module ecs_halo_vec_a_mod
