module hot_bubble_testcase_mod

use parcomm_mod,        only : parcomm_global
use domain_mod,         only : domain_t
use stvec_mod,          only : stvec_t
use stvec_nh_mod,       only : stvec_nh_t
use test_fields_3d_mod, only : scalar_field3d_t
use const_mod,          only : Cp, Earth_grav, pi

implicit none

type, extends(scalar_field3d_t) :: hot_bubble_theta_t
    real(kind=8) :: amp = 10.0_8, theta0 = 300._8
    real(kind=8) :: lam0 = pi, phi0 = 0.0_8, h0 = 3000.0_8
    real(kind=8) :: dx = 2000._8, dy = 2000._8, dz = 2000._8
    real(kind=8) :: scale
    contains
    procedure :: get_scalar_field_tile => get_hot_bubble_theta_tile
end type hot_bubble_theta_t

type, extends(scalar_field3d_t) :: hot_bubble_PExner_t
    real(kind=8) :: theta0 = 300._8
    contains
    procedure :: get_scalar_field_tile => get_hot_bubble_PExner_tile
end type hot_bubble_PExner_t

contains

subroutine get_hot_bubble_initial_conditions(stvec, domain)
    class(stvec_t),   intent(inout) :: stvec
    type(domain_t),   intent(in)    :: domain

    select type(stvec)
    type is (stvec_nh_t)
        call get_hot_bubble_initial_conditions_nh_stvec(stvec, domain)
    class default
        call parcomm_global%abort("hot_bubble testcases does not support this type of stvec")
    end select
end subroutine get_hot_bubble_initial_conditions

subroutine get_hot_bubble_initial_conditions_nh_stvec(stvec, domain)
    class(stvec_nh_t), intent(inout) :: stvec
    type(domain_t),    intent(in)    :: domain

    type(hot_bubble_theta_t)  :: theta_generator
    type(hot_bubble_PExner_t) :: P0_generator

    call stvec%u%assign(0.0_8,domain%mesh_u)
    call stvec%v%assign(0.0_8,domain%mesh_v)
    call stvec%eta_dot%assign(0.0_8,domain%mesh_w)

    theta_generator%scale = domain%mesh_p%scale
    call theta_generator%get_scalar_field(stvec%theta,domain%mesh_w,0)
    call P0_generator%get_scalar_field(stvec%P,domain%mesh_p,0)

end subroutine get_hot_bubble_initial_conditions_nh_stvec

subroutine get_hot_bubble_theta_tile(this,f,mesh,halo_width)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t
    use const_mod,      only : pi
    use sph_coords_mod, only : cart2sph

    class(hot_bubble_theta_t),  intent(in)    :: this
    type(tile_field_t),         intent(inout) :: f
    type(tile_mesh_t),          intent(in)    :: mesh
    integer(kind=4),            intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: d, lam, phi, p, H0

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    H0 = (Cp*this%theta0) / Earth_grav

    do k=ks,ke
        do j=js,je
            do i=is,ie
                d = acos(-mesh%rx(i,j,k))
                d = sqrt((d*this%scale/this%dx)**2+&
                         ((mesh%h(i,j,k)-this%h0)/this%dz)**2)
                d = min(1.0_8,d)
                f%p(i,j,k) = this%theta0+this%amp*0.5_8*(1.0_8+cos(pi*d))
                ! if(d<=1.0_8) then
                !     f%p(i,j,k) = this%theta0+this%amp!*0.5_8*(1.0_8+cos(pi*d))!/p
                ! else
                !     f%p(i,j,k) = this%theta0
                ! end if
            end do
        end do
    end do

end subroutine get_hot_bubble_theta_tile

subroutine get_hot_bubble_PExner_tile(this,f,mesh,halo_width)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    class(hot_bubble_PExner_t), intent(in)    :: this
    type(tile_field_t),         intent(inout) :: f
    type(tile_mesh_t),          intent(in)    :: mesh
    integer(kind=4),            intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: H0

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    H0 = (Cp*this%theta0) / Earth_grav

    do k=ks,ke
        do j=js,je
            do i=is,ie
                f%p(i,j,k) = 1.0_8 - mesh%h(i,j,k) / H0
            end do
        end do
    end do

end subroutine get_hot_bubble_PExner_tile

end module hot_bubble_testcase_mod
