module coriolis_Cgrid_noncons_mod

use grid_field_mod,                only : grid_field_t, tile_field_t
use domain_mod,                    only : domain_t
use abstract_coriolis_mod,         only : coriolis_operator_t
use mesh_mod,                      only : tile_mesh_t
use abstract_interpolators2d_mod,  only : interpolator2d_vec2vec_t
use abstract_co2contra_mod,        only : co2contra_operator_t

implicit none

type, public, extends(coriolis_operator_t) :: coriolis_Cgrid_noncons_t
    type(grid_field_t) :: f, u, v, uh, vh !coriolis parameter
    class(co2contra_operator_t), allocatable     :: co2contra
    class(interpolator2d_vec2vec_t), allocatable :: interp_v2h_op
    class(interpolator2d_vec2vec_t), allocatable :: interp_h2v_op

contains
    procedure, public :: calc_coriolis
    procedure, public :: calc_coriolis_vec_inv
    procedure, public :: calc_coriolis_contra
end type coriolis_Cgrid_noncons_t

contains

subroutine calc_coriolis_contra(this, cor_ut, cor_vt, ut, vt, domain)
    class(coriolis_Cgrid_noncons_t), intent(inout) :: this
    type(domain_t),                  intent(in)    :: domain
    type(grid_field_t),              intent(inout) :: ut, vt!covariant components
    type(grid_field_t),              intent(inout) :: cor_ut, cor_vt !contravariant

    integer(kind=4) :: t

    call this%co2contra%transform2co(this%u,this%v, ut, vt, domain)
    call this%interp_v2h_op%interp2d_vec2vec(this%uh, this%vh, this%u, this%v, domain)

    do t = domain%partition%ts, domain%partition%te
        call calc_coriolis_contra_on_tile(this%uh%tile(t), this%vh%tile(t), this%f%tile(t), &
                                          domain%mesh_p%tile(t))
    end do

    call this%interp_h2v_op%interp2d_vec2vec(cor_ut, cor_vt, this%uh, this%vh, domain)

    do t = domain%partition%ts, domain%partition%te
        call divide_by_G(cor_ut%tile(t), domain%mesh_u%tile(t))
        call divide_by_G(cor_vt%tile(t), domain%mesh_v%tile(t))
    end do

end subroutine calc_coriolis_contra

subroutine calc_coriolis_contra_on_tile(u, v, f, mesh)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),   intent(inout) :: u, v
    type(tile_field_t),   intent(in)    :: f
    type(tile_mesh_t),    intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    real(kind=8)    :: ut_tend, vt_tend

    do k=mesh%ks, mesh%ke
        do j=mesh%js, mesh%je
            do i=mesh%is, mesh%ie
                ut_tend = v%p(i,j,k)*f%p(i,j,1)
                vt_tend =-u%p(i,j,k)*f%p(i,j,1)

                u%p(i,j,k) = ut_tend
                v%p(i,j,k) = vt_tend
            end do
        end do
    end do

end subroutine calc_coriolis_contra_on_tile

subroutine divide_by_G(uv, mesh)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),   intent(inout) :: uv
    type(tile_mesh_t),    intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    real(kind=8)    :: ut_tend, vt_tend

    do k=mesh%ks, mesh%ke
        do j=mesh%js, mesh%je
            do i=mesh%is, mesh%ie
                uv%p(i,j,k) = uv%p(i,j,k) / mesh%J(i,j,k)
            end do
        end do
    end do

end subroutine divide_by_G

subroutine calc_coriolis(this, cor_u, cor_v, ut, vt, domain)
    class(coriolis_Cgrid_noncons_t), intent(inout) :: this
    type(domain_t),                  intent(in)    :: domain
    type(grid_field_t),              intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),              intent(inout) :: cor_u, cor_v

    integer(kind=4) :: t

    call domain%parcomm%abort("calc_coriolis is not implemented!")

end subroutine calc_coriolis

subroutine calc_coriolis_vec_inv(this, cor_u, cor_v, hu, hv, h, curl, domain)
    class(coriolis_Cgrid_noncons_t), intent(inout) :: this
    type(domain_t),                  intent(in)    :: domain
    type(grid_field_t),              intent(inout) :: hu, hv! massflux contravariant components
    type(grid_field_t),              intent(inout) :: h, curl
    type(grid_field_t),              intent(inout) :: cor_u, cor_v

    call domain%parcomm%abort("calc_coriolis_vec_inv is not implemented!")

end subroutine calc_coriolis_vec_inv


end module coriolis_Cgrid_noncons_mod
