module advection3d_testcases_mod

use parcomm_mod,        only : parcomm_global
use domain_mod,         only : domain_t
use stvec_mod,          only : stvec_t
use stvec_nh_mod,       only : stvec_nh_t
use parcomm_mod,        only : parcomm_global
use test_fields_3d_mod, only : scalar_field3d_t
use const_mod,          only : pi

implicit none

type, extends(scalar_field3d_t) :: cosine_bell_t
    real(kind=8) :: h0 = 5e3
    real(kind=8) :: r = 7.0_8*pi/64.0, rz = 2.5e3_8
contains
    procedure :: get_scalar_field_tile => get_cosine_bell_tile
end type cosine_bell_t

contains

subroutine get_advection3d_initial_conditions(stvec, domain, testcase_name)
    class(stvec_t),   intent(inout) :: stvec
    type(domain_t),   intent(in)    :: domain
    character(len=*), intent(in)    :: testcase_name

    select case(testcase_name)
    case("advection3d_solid_rotation")
        call get_advection3d_solid_rotation_initial_conditions(stvec, domain)
    case("Straka_buble")
        call get_Straka_buble_initial_conditions(stvec,domain)
    case default
        call parcomm_global%abort("unknown advection3d testcase name: "//testcase_name)
    end select

end subroutine get_advection3d_initial_conditions

subroutine get_advection3d_solid_rotation_initial_conditions(stvec, domain)
    class(stvec_t),   intent(inout) :: stvec
    type(domain_t),   intent(in)    :: domain

    type(cosine_bell_t) :: f_generator


    select type(stvec)
    type is (stvec_nh_t)
        call stvec%u%assign(0.0_8,domain%mesh_u)
        call stvec%v%assign(0.0_8,domain%mesh_v)
        call stvec%eta_dot%assign(0.0_8,domain%mesh_w)
        f_generator = cosine_bell_t(r = 7.0_8*pi/64.0_8,rz = 1.5e3_8)
        call f_generator%get_scalar_field(stvec%theta,domain%mesh_w,0)
        call f_generator%get_scalar_field(stvec%P,domain%mesh_p,0)
    class default
    call parcomm_global%abort("unsupported stvec type in get_advection3d_solid_rotation_initial_conditions")
    end select
end subroutine get_advection3d_solid_rotation_initial_conditions

subroutine get_Straka_buble_initial_conditions(stvec, domain)
    use Straka_testcase_mod, only : Straka_theta_t

    class(stvec_t),   intent(inout) :: stvec
    type(domain_t),   intent(in)    :: domain

    type(Straka_theta_t) :: f_generator

    select type(stvec)
    type is (stvec_nh_t)
        call stvec%u%assign(0.0_8,domain%mesh_u)
        call stvec%v%assign(0.0_8,domain%mesh_v)
        call stvec%eta_dot%assign(0.0_8,domain%mesh_w)
        f_generator%scale = domain%mesh_p%scale
        f_generator%h0 = 2500.0_8
        call f_generator%get_scalar_field(stvec%theta,domain%mesh_w,0)
        call f_generator%get_scalar_field(stvec%P,domain%mesh_p,0)
    class default
    call parcomm_global%abort("unsupported stvec type in get_advection3d_solid_rotation_initial_conditions")
    end select
end subroutine get_Straka_buble_initial_conditions

subroutine get_cosine_bell_tile(this,f,mesh,halo_width)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t
    use const_mod,      only : pi

    class(cosine_bell_t), intent(in)    :: this
    type(tile_field_t),   intent(inout) :: f
    type(tile_mesh_t),    intent(in)    :: mesh
    integer(kind=4),      intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: d

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                d = sqrt(acos(mesh%rx(i,j,k))**2 / this%r**2 + &
                         (mesh%h(i,j,k)-this%h0)**2/this%rz**2)
                d = min(d,1.0_8)
                f%p(i,j,k) = 0.5_8*(1.0_8+cos(pi*d))
            end do
        end do
    end do

end subroutine get_cosine_bell_tile

end module advection3d_testcases_mod
