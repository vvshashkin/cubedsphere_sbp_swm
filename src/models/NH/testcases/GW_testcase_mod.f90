module NH_GW_testcase_mod

use parcomm_mod,        only : parcomm_global
use domain_mod,         only : domain_t
use stvec_mod,          only : stvec_t
use stvec_nh_mod,       only : stvec_nh_t
use test_fields_3d_mod, only : scalar_field3d_t

use solid_rotation_fields_factory_mod, only : create_solid_rotation_field_generators

implicit none

type, extends(scalar_field3d_t) :: GW_theta_t
    real(kind=8) :: amp = 0.01_8
    real(kind=8) :: a = 5e3, Lz = 10e3
    real(kind=8) :: hor_scale
    contains
    procedure :: get_scalar_field_tile => get_GWtheta_tile
end type GW_theta_t

contains

subroutine get_GW_initial_conditions(stvec, domain, testcase_name)
    class(stvec_t),   intent(inout) :: stvec
    type(domain_t),   intent(in)    :: domain
    character(len=*), intent(in)    :: testcase_name

    logical :: is_nonlin

    select case(testcase_name)
    case ("GW")
        is_nonlin = .true.
    case("GW_linear")
        is_nonlin = .false.
    case default
        call parcomm_global%abort("unknown type of NH GW test: "// testcase_name)
    end select

    select type(stvec)
    type is (stvec_nh_t)
        call get_GW_initial_conditions_nh_stvec(stvec, domain, is_nonlin)
    class default
        call parcomm_global%abort("NH GW testcases does not support this type of stvec")
    end select

end subroutine get_GW_initial_conditions

subroutine get_GW_initial_conditions_nh_stvec(stvec, domain, is_nonlin)

    use grid_field_mod,                only: grid_field_t
    use grid_field_factory_mod,        only: create_grid_field
    use test_fields_3d_mod,            only: scalar_field3d_t, vector_field3d_t
    use const_mod,                     only: Earth_grav

    class(stvec_nh_t),   intent(inout) :: stvec
    type(domain_t),      intent(in)    :: domain
    logical,             intent(in)    :: is_nonlin

    real(kind=8), parameter :: Nb = 0.01_8, T0 = 300.0_8, pref = 1e5_8, U0 = 20.0_8

    type(GW_theta_t) :: theta_generator
    type(grid_field_t) :: P0, theta0

    class(scalar_field3d_t), allocatable :: theta0_generator
    class(scalar_field3d_t), allocatable :: P0_generator
    class(vector_field3d_t), allocatable :: wind_gen

    call stvec%u%assign(0.0_8,domain%mesh_u)
    call stvec%v%assign(0.0_8,domain%mesh_v)
    call stvec%eta_dot%assign(0.0_8,domain%mesh_w)
    call stvec%P%assign(0.0_8,domain%mesh_p)

    theta_generator = GW_theta_t(amp=1._8,a=5e3,Lz=10e3,hor_scale=domain%mesh_w%scale)
    call theta_generator%get_scalar_field(stvec%theta,domain%mesh_w,0)

    if(is_nonlin) then
        call create_grid_field(P0,0,0,domain%mesh_p)
        call create_grid_field(theta0,0,0,domain%mesh_w)
        call create_solid_rotation_field_generators("Nb_const",u0=U0,omega=0.0_8,       &
                                                    sphere_rad = domain%mesh_p%scale, Nb = Nb, &
                                                    grav = Earth_grav,alpha=0.0_8,      &
                                                    theta_gen  = theta0_generator, &
                                                    PExner_gen = P0_generator,  &
                                                    wind_gen   = wind_gen)

        call wind_gen%get_vector_field(stvec%u,stvec%v,stvec%eta_dot, &
                                       domain%mesh_u,domain%mesh_v,domain%mesh_w,0,"contravariant")

        call P0_generator%get_scalar_field(P0,domain%mesh_p,0)
        call theta0_generator%get_scalar_field(theta0,domain%mesh_w,0)

        call stvec%P%update(1.0_8,P0,domain%mesh_p)
        call stvec%theta%update(1.0_8,theta0,domain%mesh_w)
    end if

end subroutine get_GW_initial_conditions_nh_stvec

subroutine get_GWtheta_tile(this,f,mesh,halo_width)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t
    use const_mod,      only : pi

    class(GW_theta_t),    intent(in)    :: this
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
                d = acos(mesh%rx(i,j,k))*this%hor_scale
                f%p(i,j,k) = this%amp*this%a**2 / (this%a**2+d**2)*sin(pi*mesh%h(i,j,k)/this%Lz)
            end do
        end do
    end do

end subroutine get_GWtheta_tile

end module NH_GW_testcase_mod
