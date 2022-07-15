module stvec_iomega_mod

use stvec_mod, only: stvec_t
use parcomm_mod,            only : parcomm_global
use domain_mod,             only : domain_t

implicit none

type, extends(stvec_t) :: stvec_iomega_t
    integer(kind=4) N
    complex(kind=8), allocatable :: f(:)

    contains
    procedure, public :: copy_to => copy_stvec_to
    procedure, public :: create_similar => create_similar_stvec
    procedure, public :: update_s1v1 => update_stvec_s1v1 !v = v + s1*v1
    procedure, public :: update_s1v1s2v2 => update_stvec_s1v1s2v2 !v = v + s1*v1+s2*v

    procedure, public :: assign_s1         => assign_stvec_s1 !v = s1
    procedure, public :: assign_s1v1       => assign_stvec_s1v1 !v = s1*v1
    procedure, public :: assign_s1v1s2v2   => assign_stvec_s1v1s2v2 !v = s1*v1+s2*v2

    procedure, public :: algebraic_norm2 => compute_stvec_algebraic_norm2
    procedure, public :: algebraic_dot   => compute_stvec_algebraic_dot

end type stvec_iomega_t

contains

subroutine init_stvec_iomega(new_stvec, N, f)
    type(stvec_iomega_t),      intent(out) :: new_stvec
    integer(kind=4),           intent(in)  :: N
    complex(kind=8), optional, intent(in)  :: f(1:N)

    if(allocated(new_stvec%f) .and. new_stvec%N /= N) then
        deallocate(new_stvec%f)
    end if
    if (.not. allocated(new_stvec%f)) then
        allocate(new_stvec%f(1:N))
    end if

    new_stvec%N = N
    if(present(f)) then
        new_stvec%f(1:N) = f(1:N)
    end if

end subroutine init_stvec_iomega

subroutine copy_stvec_to(this,destination)
    class(stvec_iomega_t),       intent(in)    :: this
    class(stvec_t), allocatable, intent(inout) :: destination

    allocate(stvec_iomega_t :: destination)
    select type (destination)
    class is (stvec_iomega_t)
        call init_stvec_iomega(destination,this%N,this%f)
    class default
        call parcomm_global%abort("stvec_iomega_t%copy_to: impossible error")
    end select

end subroutine copy_stvec_to

subroutine create_similar_stvec(this,destination)
    class(stvec_iomega_t),       intent(in)    :: this
    class(stvec_t), allocatable, intent(inout) :: destination

    allocate(stvec_iomega_t :: destination)
    select type (destination)
    class is (stvec_iomega_t)
        call init_stvec_iomega(destination,this%N)
    class default
        call parcomm_global%abort("stvec_iomega_t%copy_to: types error")
    end select

end subroutine create_similar_stvec

subroutine update_stvec_s1v1(this, scalar1, v1, domain)

    class(stvec_iomega_t), intent(inout) :: this
    real(kind=8),          intent(in)    :: scalar1
    class(stvec_t),        intent(in)    :: v1
    class(domain_t),       intent(in)    :: domain

    integer(kind=4) :: N

    select type (v1)
    class is (stvec_iomega_t)
        N = this%N
        this%f(1:N) = this%f(1:N)+scalar1*v1%f(1:N)
    class default
        call parcomm_global%abort("stvec_iomega_t%update_s1v1: types error")
    end select
end subroutine update_stvec_s1v1

subroutine update_stvec_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

    class(stvec_iomega_t), intent(inout) :: this
    real(kind=8),          intent(in)    :: scalar1
    class(stvec_t),        intent(in)    :: v1
    real(kind=8),          intent(in)    :: scalar2
    class(stvec_t),        intent(in)    :: v2
    class(domain_t),       intent(in)    :: domain

    integer(kind=4) :: N

    select type (v1)
    class is (stvec_iomega_t)
    select type(v2)
    class is (stvec_iomega_t)

    N = this%N
    this%f(1:N) = this%f(1:N)+scalar1*v1%f(1:N)+scalar2*v2%f(1:N)
    class default
        call parcomm_global%abort("stvec_iomega_t%update_s1v1s2v2: types error")
    end select
    class default
        call parcomm_global%abort("stvec_iomega_t%update_s1v1s2v2: types error")
    end select


end subroutine update_stvec_s1v1s2v2

subroutine assign_stvec_s1(this, scalar1, domain)

    class(stvec_iomega_t), intent(inout) :: this
    real(kind=8),          intent(in)    :: scalar1
    class(domain_t),       intent(in)    :: domain

    integer(kind=4) :: N

    N = this%N
    this%f(1:N) = scalar1
end subroutine assign_stvec_s1


subroutine assign_stvec_s1v1(this, scalar1, v1, domain)

    class(stvec_iomega_t), intent(inout) :: this
    real(kind=8),          intent(in)    :: scalar1
    class(stvec_t),        intent(in)    :: v1
    class(domain_t),       intent(in)    :: domain

    integer(kind=4) :: N

    select type (v1)
    class is (stvec_iomega_t)
        N = this%N
        this%f(1:N) = scalar1*v1%f(1:N)
    class default
        call parcomm_global%abort("stvec_iomega_t%assign_s1v1: types error")
    end select
end subroutine assign_stvec_s1v1

subroutine assign_stvec_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

    class(stvec_iomega_t), intent(inout) :: this
    real(kind=8),          intent(in)    :: scalar1
    class(stvec_t),        intent(in)    :: v1
    real(kind=8),          intent(in)    :: scalar2
    class(stvec_t),        intent(in)    :: v2
    class(domain_t),       intent(in)    :: domain

    integer(kind=4) :: N

    select type (v1)
    class is (stvec_iomega_t)
    select type(v2)
    class is (stvec_iomega_t)

    N = this%N
    this%f(1:N) = scalar1*v1%f(1:N)+scalar2*v2%f(1:N)
    class default
        call parcomm_global%abort("stvec_iomega_t%assign_s1v1s2v2: types error")
    end select
    class default
        call parcomm_global%abort("stvec_iomega_t%assign_s1v1s2v2: types error")
    end select

end subroutine assign_stvec_s1v1s2v2

function compute_stvec_algebraic_norm2(this, domain) result(norm2)
    use parcomm_mod, only : parcomm_t
    use mpi

    class(stvec_iomega_t), intent(in)  :: this
    class(domain_t),       intent(in)  :: domain
    real(kind=8)                       :: norm2

    norm2 = sqrt(sum(real(this%f(1:this%n))**2+imag(this%f(1:this%n))**2))

end function compute_stvec_algebraic_norm2

function compute_stvec_algebraic_dot(this, other, domain) result(dot_product)

    class(stvec_iomega_t), intent(in)  :: this
    class(stvec_t),        intent(in)  :: other
    class(domain_t),       intent(in)  :: domain
    real(kind=8)                       :: dot_product

    integer(kind=4) :: N

    dot_product = 0.0_8
    select type (other)
    class is (stvec_iomega_t)
        N = this%N
        dot_product = sum(real(this%f(1:N))*real(other%f(1:N))+ &
                          imag(this%f(1:N))*imag(other%f(1:N)))
    class default
        call parcomm_global%abort("stvec_iomega_t%algebraic_dot: types error")
    end select
end function compute_stvec_algebraic_dot

end module stvec_iomega_mod
