module coriolis_colocated_mod

use grid_field_mod,        only : grid_field_t, tile_field_t
use domain_mod,            only : domain_t
use abstract_coriolis_mod, only : coriolis_operator_t

implicit none

type, public, extends(coriolis_operator_t) :: coriolis_colocated_t
    type(grid_field_t) :: f !coriolis parameter
contains
    procedure, public :: calc_coriolis
    procedure, public :: calc_coriolis_contra
    procedure, public :: calc_coriolis_vec_inv
end type coriolis_colocated_t

contains

subroutine calc_coriolis(this, cor_u, cor_v, ut, vt, domain)
    class(coriolis_colocated_t), intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),     intent(inout) :: cor_u, cor_v

    integer(kind=4) :: t

    do t = domain%partition%ts, domain%partition%te
        call calc_coriolis_on_tile(ut%tile(t), vt%tile(t), this%f%tile(t), &
                                   cor_u%tile(t), cor_v%tile(t), &
                                   domain%mesh_u%tile(t), domain%mesh_v%tile(t))
    end do


end subroutine calc_coriolis

subroutine calc_coriolis_on_tile(ut, vt, f, cor_u, cor_v, mesh_u, mesh_v)

    !coriolis(\vec{u}) = f*\vec{u}^\perp

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t), intent(in)    :: ut, vt, f
    type(tile_field_t), intent(inout) :: cor_u, cor_v
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k
    real(kind=8)    :: u_covariant, v_covariant

!This implementation works only for colocatedgered case, so mesh_u==mesh_v==mesh_p
    do k = mesh_u%ks, mesh_u%ke
        do j = mesh_u%js, mesh_u%je
            do i = mesh_u%is, mesh_u%ie
                cor_u%p(i,j,k) =  f%p(i,j,1)*vt%p(i,j,k)*mesh_u%J(i,j,k)
                cor_v%p(i,j,k) = -f%p(i,j,1)*ut%p(i,j,k)*mesh_u%J(i,j,k)
            end do
        end do
    end do


end subroutine calc_coriolis_on_tile
subroutine calc_coriolis_vec_inv(this, cor_u, cor_v, hu, hv, h, curl, domain)
    class(coriolis_colocated_t), intent(inout) :: this
    type(domain_t),              intent(in)    :: domain
    type(grid_field_t),          intent(inout) :: hu, hv! massflux contravariant components
    type(grid_field_t),          intent(inout) :: h, curl
    type(grid_field_t),          intent(inout) :: cor_u, cor_v

    integer(kind=4) :: t

    do t = domain%partition%ts, domain%partition%te
        call calc_coriolis_vec_inv_tile(hu%tile(t), hv%tile(t), h%tile(t), curl%tile(t), this%f%tile(t), &
                                        cor_u%tile(t), cor_v%tile(t), &
                                        domain%mesh_u%tile(t), domain%mesh_v%tile(t))
    end do


end subroutine calc_coriolis_vec_inv

subroutine calc_coriolis_vec_inv_tile(hu, hv, h, curl, f, cor_u, cor_v, mesh_u, mesh_v)

    !coriolis(\vec{u}) = f*\vec{u}^\perp

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t), intent(in)    :: hu, hv, curl, f, h
    type(tile_field_t), intent(inout) :: cor_u, cor_v
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k
    real(kind=8)    :: u_covariant, v_covariant

!This implementation works only for unstaggered case, so mesh_u==mesh_v==mesh_p
    do k = mesh_u%ks, mesh_u%ke
        do j = mesh_u%js, mesh_u%je
            do i = mesh_u%is, mesh_u%ie
                cor_u%p(i,j,k) =  (f%p(i,j,1)+curl%p(i,j,k))*hv%p(i,j,k)/h%p(i,j,k)*mesh_u%J(i,j,k)
                cor_v%p(i,j,k) = -(f%p(i,j,1)+curl%p(i,j,k))*hu%p(i,j,k)/h%p(i,j,k)*mesh_u%J(i,j,k)
            end do
        end do
    end do


end subroutine calc_coriolis_vec_inv_tile

subroutine calc_coriolis_contra(this, cor_ut, cor_vt, ut, vt, domain)
    class(coriolis_colocated_t), intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: ut, vt!covariant components
    type(grid_field_t),     intent(inout) :: cor_ut, cor_vt !contravariant

    integer(kind=4) :: t

    do t = domain%partition%ts, domain%partition%te
        call calc_coriolis_contra_on_tile(ut%tile(t), vt%tile(t), this%f%tile(t), &
                                          cor_ut%tile(t), cor_vt%tile(t), &
                                          domain%mesh_u%tile(t), domain%mesh_v%tile(t))
    end do


end subroutine calc_coriolis_contra

subroutine calc_coriolis_contra_on_tile(ut, vt, f, cor_ut, cor_vt, mesh_u, mesh_v)

    !coriolis(\vec{u}) = f*\vec{u}^\perp

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t), intent(in)    :: ut, vt, f
    type(tile_field_t), intent(inout) :: cor_ut, cor_vt
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k
    real(kind=8)    :: u_covariant, v_covariant

!This implementation works only for colocatedgered case, so mesh_u==mesh_v==mesh_p
    do k = mesh_u%ks, mesh_u%ke
        do j = mesh_u%js, mesh_u%je
            do i = mesh_u%is, mesh_u%ie
                u_covariant = mesh_u%Q(1,i,j,k)*ut%p(i,j,k)+mesh_u%Q(2,i,j,k)*vt%p(i,j,k)
                v_covariant = mesh_v%Q(2,i,j,k)*ut%p(i,j,k)+mesh_v%Q(3,i,j,k)*vt%p(i,j,k)
                cor_ut%p(i,j,k) =  f%p(i,j,1)*v_covariant/mesh_u%J(i,j,k)
                cor_vt%p(i,j,k) = -f%p(i,j,1)*u_covariant/mesh_v%J(i,j,k)
            end do
        end do
    end do


end subroutine calc_coriolis_contra_on_tile

end module coriolis_colocated_mod
