module stvec_abstract_mod

use container_abstract_mod, only : state_abstract_t
use parcomm_mod,            only : parcomm_global

implicit none

type, abstract, public, extends(state_abstract_t) :: stvec_abstract_t

contains

    procedure(add),   deferred :: add
    procedure(copy),  deferred :: copy
    procedure                  :: dot
    procedure                  :: norm
    procedure                  :: smult

end type stvec_abstract_t

abstract interface
    subroutine add(this, other, alpha, beta)
        !calculates linear combination alpha*this+beta*other
        import stvec_abstract_t
        class(stvec_abstract_t), intent(inout) :: this
        class(stvec_abstract_t), intent(in)    :: other
        real(kind=8),            intent(in)    :: alpha, beta
    end subroutine add

    subroutine copy(this, source_stvec)
        !copies information from source_stvec to this
        import stvec_abstract_t
        class(stvec_abstract_t), intent(inout) :: this
        class(stvec_abstract_t), intent(in)    :: source_stvec
    end subroutine copy
end interface

contains

real(kind=8) function norm(this) result(l2)
    class(stvec_abstract_t), intent(in) :: this
    !calculates norm of state vector
    l2 = 0._8
    call parcomm_global%abort("norm function not implemented for specific stvec class")
end function norm

subroutine smult(this, alpha)
    !multiplicates state vector by scalar
    class(stvec_abstract_t), intent(inout) :: this
    real(kind=8),            intent(in)    :: alpha

    call parcomm_global%abort("scalar multiplication not implemented for specific stvec class")
end subroutine smult

real(kind=8) function dot(this,other) result(dot_prod)
    !calculates dot(this,other)
    class(stvec_abstract_t), intent(in) :: this
    class(stvec_abstract_t), intent(in) :: other

    dot_prod = 0._8
    call parcomm_global%abort("dot product function is not implemented for specific stvec class")
end function dot

end module stvec_abstract_mod
