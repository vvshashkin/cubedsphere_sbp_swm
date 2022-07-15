module stvec_mod

use container_abstract_mod, only : state_abstract_t
use parcomm_mod,            only : parcomm_global
use domain_mod,             only : domain_t

implicit none

private
type, abstract, public, extends(state_abstract_t) :: stvec_t

contains
    procedure, public :: copy_to => copy_stvec_to
    procedure, public :: create_similar => create_similar_stvec

    procedure, public :: update_s1v1 => update_stvec_s1v1 !v = v + s1*v1
    procedure, public :: update_s1v1s2v2 => update_stvec_s1v1s2v2 !v = v + s1*v1+s2*v2
    generic :: update => update_s1v1, update_s1v1s2v2

    procedure, public :: assign_s1         => assign_stvec_s1 !v = s1
    procedure, public :: assign_s1v1       => assign_stvec_s1v1 !v = s1*v1
    procedure, public :: assign_s1v1s2v2   => assign_stvec_s1v1s2v2 !v = s1*v1+s2*v2
    generic :: assign => assign_s1v1, assign_s1, assign_s1v1s2v2

    procedure, public :: algebraic_norm2 => compute_stvec_algebraic_norm2
    procedure, public :: algebraic_dot   => compute_stvec_algebraic_dot

end type stvec_t

contains

subroutine copy_stvec_to(this, destination)
    class(stvec_t),              intent(in)    :: this
    class(stvec_t), allocatable, intent(inout) :: destination

    call parcomm_global%abort("copy_stvec_to not implemented for specific stvec class")

end subroutine copy_stvec_to

subroutine create_similar_stvec(this,destination)
    class(stvec_t),              intent(in)    :: this
    class(stvec_t), allocatable, intent(inout) :: destination

    call parcomm_global%abort("create_similar_stvec not implemented for specific stvec class")

end subroutine create_similar_stvec

subroutine update_stvec_s1v1(this, scalar1, v1, domain)

    class(stvec_t),      intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    class(stvec_t),      intent(in)    :: v1
    class(domain_t),     intent(in)    :: domain

    call parcomm_global%abort("update_stvec_s1v1 not implemented for specific stvec class")

end subroutine update_stvec_s1v1

subroutine update_stvec_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

    class(stvec_t),      intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    class(stvec_t),      intent(in)    :: v1
    real(kind=8),        intent(in)    :: scalar2
    class(stvec_t),      intent(in)    :: v2
    class(domain_t),     intent(in)    :: domain

    call parcomm_global%abort("update_stvec_s1v1s2v2 not implemented for specific stvec class")

end subroutine update_stvec_s1v1s2v2

subroutine assign_stvec_s1(this, scalar1, domain)

    class(stvec_t),      intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    class(domain_t),     intent(in)    :: domain

    call parcomm_global%abort("assign_stvec_s1 not implemented for specific stvec class")

end subroutine assign_stvec_s1

subroutine assign_stvec_s1v1(this, scalar1, v1, domain)

    class(stvec_t),      intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    class(stvec_t),      intent(in)    :: v1
    class(domain_t),     intent(in)    :: domain

    call parcomm_global%abort("assign_stvec_s1v1 not implemented for specific stvec class")

end subroutine assign_stvec_s1v1

subroutine assign_stvec_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

    class(stvec_t),      intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    class(stvec_t),      intent(in)    :: v1
    real(kind=8),        intent(in)    :: scalar2
    class(stvec_t),      intent(in)    :: v2
    class(domain_t),     intent(in)    :: domain

    call parcomm_global%abort("assign_stvec_s1v1s2v2 not implemented for specific stvec class")

end subroutine assign_stvec_s1v1s2v2

function compute_stvec_algebraic_norm2(this, domain) result(norm2)
    use parcomm_mod, only : parcomm_t
    use mpi

    class(stvec_t),      intent(in)  :: this
    class(domain_t),     intent(in)  :: domain
    real(kind=8)                     :: norm2

    norm2 = 0.0_8
    call parcomm_global%abort("compute_stvec_algebraic_norm2 not implemented for specific stvec class")

end function compute_stvec_algebraic_norm2

function compute_stvec_algebraic_dot(this, other, domain) result(dot_product)

    class(stvec_t),      intent(in)  :: this
    class(stvec_t),      intent(in)  :: other
    class(domain_t),     intent(in)  :: domain
    real(kind=8)                     :: dot_product

    dot_product = 0.0_8
    call parcomm_global%abort("compute_stvec_algebraic_dot not implemented for specific stvec class")

end function compute_stvec_algebraic_dot

end module stvec_mod
