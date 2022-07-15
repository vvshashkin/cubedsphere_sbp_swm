module stvec_swm_mod

use stvec_mod,      only : stvec_t
use grid_field_mod, only : grid_field_t
use parcomm_mod,    only : parcomm_global
use domain_mod,     only : domain_t

implicit none

type, extends(stvec_t) :: stvec_swm_t
    type(grid_field_t) :: h, u, v
contains
    procedure, public :: create_similar
    procedure, public :: assign_s1v1s2v2
    procedure, public :: update_s1v1s2v2
    procedure, public :: update_s1v1
end type stvec_swm_t

contains



subroutine create_similar(this, destination)
    class(stvec_swm_t),          intent(in)    :: this
    class(stvec_t), allocatable, intent(inout) :: destination

    allocate(stvec_swm_t :: destination)
    select type (destination)
    class is (stvec_swm_t)
        destination%h = this%h%create_similar()
        destination%u = this%u%create_similar()
        destination%v = this%v%create_similar()
    class default
        call parcomm_global%abort("stvec_swm_t%create_similar: types error")
    end select
end subroutine create_similar

subroutine assign_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

    class(stvec_swm_t), intent(inout) :: this
    real(kind=8),       intent(in)       :: scalar1, scalar2
    class(stvec_t),     intent(in)       :: v1, v2
    class(domain_t),    intent(in)       :: domain

    integer(kind=4) :: N

    select type (v1)
    class is (stvec_swm_t)
        select type(v2)
        class is (stvec_swm_t)
            call this%h%assign(scalar1, v1%h, scalar2, v2%h, domain%mesh_p)
            call this%u%assign(scalar1, v1%u, scalar2, v2%u, domain%mesh_u)
            call this%v%assign(scalar1, v1%v, scalar2, v2%v, domain%mesh_v)

    class default
        call parcomm_global%abort("stvec_swm_t%assign_s1v1s2v2: types error")
    end select
    class default
        call parcomm_global%abort("stvec_swm_t%assign_s1v1s2v2: types error")
    end select

end subroutine assign_s1v1s2v2
subroutine update_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

    class(stvec_swm_t), intent(inout) :: this
    real(kind=8),       intent(in)       :: scalar1, scalar2
    class(stvec_t),     intent(in)       :: v1, v2
    class(domain_t),    intent(in)       :: domain

    integer(kind=4) :: N

    select type (v1)
    class is (stvec_swm_t)
        select type(v2)
        class is (stvec_swm_t)
            call this%h%update(scalar1, v1%h, scalar2, v2%h, domain%mesh_p)
            call this%u%update(scalar1, v1%u, scalar2, v2%u, domain%mesh_u)
            call this%v%update(scalar1, v1%v, scalar2, v2%v, domain%mesh_v)

    class default
        call parcomm_global%abort("stvec_swm_t%update_s1v1s2v2: types error")
    end select
    class default
        call parcomm_global%abort("stvec_swm_t%update_s1v1s2v2: types error")
    end select

end subroutine update_s1v1s2v2
subroutine update_s1v1(this, scalar1, v1, domain)

    class(stvec_swm_t), intent(inout) :: this
    real(kind=8),       intent(in)       :: scalar1
    class(stvec_t),     intent(in)       :: v1
    class(domain_t),    intent(in)       :: domain

    integer(kind=4) :: N

    select type (v1)
    class is (stvec_swm_t)
        call this%h%update(scalar1, v1%h, domain%mesh_p)
        call this%u%update(scalar1, v1%u, domain%mesh_u)
        call this%v%update(scalar1, v1%v, domain%mesh_v)
    class default
        call parcomm_global%abort("stvec_swm_t%update_s1v1: types error")
    end select

end subroutine update_s1v1


end module stvec_swm_mod
