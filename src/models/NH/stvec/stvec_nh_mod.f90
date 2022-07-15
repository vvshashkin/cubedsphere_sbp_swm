module stvec_nh_mod

use stvec_mod,      only : stvec_t
use grid_field_mod, only : grid_field_t
use parcomm_mod,    only : parcomm_global
use domain_mod,     only : domain_t

implicit none

type, extends(stvec_t) :: stvec_nh_t
    type(grid_field_t) :: u, v, eta_dot, theta, P
    real(kind=8)       :: model_time
contains
    procedure, public :: create_similar
    procedure, public :: assign_s1v1s2v2
    procedure, public :: update_s1v1s2v2
    procedure, public :: update_s1v1
end type stvec_nh_t

contains

subroutine create_similar(this, destination)
    class(stvec_nh_t),           intent(in)    :: this
    class(stvec_t), allocatable, intent(inout) :: destination

    allocate(stvec_nh_t :: destination)
    select type (destination)
    class is (stvec_nh_t)
        destination%u = this%u%create_similar()
        destination%v = this%v%create_similar()
        destination%eta_dot = this%eta_dot%create_similar()
        destination%theta = this%theta%create_similar()
        destination%P = this%P%create_similar()
        destination%model_time = 0.0_8
    class default
        call parcomm_global%abort("stvec_nh_t%create_similar: types error")
    end select
end subroutine create_similar

subroutine assign_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

    class(stvec_nh_t), intent(inout) :: this
    real(kind=8),       intent(in)       :: scalar1, scalar2
    class(stvec_t),     intent(in)       :: v1, v2
    class(domain_t),    intent(in)       :: domain

    integer(kind=4) :: N

    select type (v1)
    class is (stvec_nh_t)
        select type(v2)
        class is (stvec_nh_t)
            call this%u%assign(scalar1, v1%u, scalar2, v2%u, domain%mesh_u)
            call this%v%assign(scalar1, v1%v, scalar2, v2%v, domain%mesh_v)
            call this%eta_dot%assign(scalar1, v1%eta_dot, scalar2, v2%eta_dot, domain%mesh_w)
            call this%theta%assign(scalar1, v1%theta, scalar2, v2%theta, domain%mesh_w)
            call this%P%assign(scalar1, v1%P, scalar2, v2%P, domain%mesh_p)
            this%model_time = scalar1*v1%model_time+scalar2*v2%model_time
    class default
        call parcomm_global%abort("stvec_nh_t%assign_s1v1s2v2: types error")
    end select
    class default
        call parcomm_global%abort("stvec_nh_t%assign_s1v1s2v2: types error")
    end select

end subroutine assign_s1v1s2v2
subroutine update_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

    class(stvec_nh_t), intent(inout) :: this
    real(kind=8),       intent(in)       :: scalar1, scalar2
    class(stvec_t),     intent(in)       :: v1, v2
    class(domain_t),    intent(in)       :: domain

    integer(kind=4) :: N

    select type (v1)
    class is (stvec_nh_t)
        select type(v2)
        class is (stvec_nh_t)
            call this%u%update(scalar1, v1%u, scalar2, v2%u, domain%mesh_u)
            call this%v%update(scalar1, v1%v, scalar2, v2%v, domain%mesh_v)
            call this%eta_dot%update(scalar1, v1%eta_dot, scalar2, v2%eta_dot, domain%mesh_w)
            call this%theta%update(scalar1, v1%theta, scalar2, v2%theta, domain%mesh_w)
            call this%P%update(scalar1, v1%P, scalar2, v2%P, domain%mesh_p)
            this%model_time = this%model_time+scalar1*v1%model_time+scalar2*v2%model_time
    class default
        call parcomm_global%abort("stvec_nh_t%update_s1v1s2v2: types error")
    end select
    class default
        call parcomm_global%abort("stvec_nh_t%update_s1v1s2v2: types error")
    end select

end subroutine update_s1v1s2v2
subroutine update_s1v1(this, scalar1, v1, domain)

    class(stvec_nh_t), intent(inout) :: this
    real(kind=8),       intent(in)       :: scalar1
    class(stvec_t),     intent(in)       :: v1
    class(domain_t),    intent(in)       :: domain

    integer(kind=4) :: N

    select type (v1)
    class is (stvec_nh_t)
        call this%u%update(scalar1, v1%u, domain%mesh_u)
        call this%v%update(scalar1, v1%v, domain%mesh_v)
        call this%eta_dot%update(scalar1, v1%eta_dot, domain%mesh_w)
        call this%theta%update(scalar1, v1%theta, domain%mesh_w)
        call this%P%update(scalar1, v1%P, domain%mesh_p)
        this%model_time = this%model_time+scalar1*v1%model_time
    class default
        call parcomm_global%abort("stvec_nh_t%update_s1v1: types error")
    end select

end subroutine update_s1v1


end module stvec_nh_mod
