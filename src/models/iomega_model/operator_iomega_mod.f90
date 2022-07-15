module operator_iomega_mod

use operator_mod,     only: operator_t
use stvec_mod,        only: stvec_t
use stvec_iomega_mod, only: stvec_iomega_t
use domain_mod,       only: domain_t
use parcomm_mod,      only : parcomm_global


implicit none

private
public :: init_iomega_operator

type, public, extends(operator_t) :: operator_iomega_t
    integer(kind=4) N
    complex(kind=8), allocatable :: omega(:) !eigen values, imag == oscillation frequency,
                                             !              real == amplification/decay
    contains
    procedure, public :: apply => apply_iomega
    procedure, public :: solve => solve_iomega
end type operator_iomega_t

contains

subroutine init_iomega_operator(operator, omega)
    type(operator_iomega_t), allocatable :: operator
    complex(kind=8)                      :: omega(:)

    allocate(operator)
    operator%N = size(omega)
    operator%omega = omega
end

subroutine apply_iomega(this,vout,vin,domain)
    class(operator_iomega_t),  intent(inout) :: this
    class(stvec_t),            intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),            intent(inout) :: vin
    type(domain_t),            intent(in)    :: domain

    integer i

    select type (vout)
    class is (stvec_iomega_t)
    select type (vin)
    class is (stvec_iomega_t)
        do i=1,this%N
            vout%f(i) = this%omega(i)*vin%f(i)
        end do
    class default
        call parcomm_global%abort("iomega operator failure: vin of wrong type")
    end select
    class default
        call parcomm_global%abort("iomega operator failure: vout of wrong type")
    end select

end subroutine apply_iomega

subroutine solve_iomega(this, vout, rhs, dt, domain)
    class(operator_iomega_t), intent(inout) :: this
    class(stvec_t),           intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),           intent(inout) :: rhs
    real(kind=8),             intent(in)    :: dt
    type(domain_t),           intent(in)    :: domain

    integer(kind=4) i

    select type (vout)
    class is (stvec_iomega_t)
    select type (rhs)
    class is (stvec_iomega_t)
        do i=1,this%N
            vout%f(i) = rhs%f(i) / (1.0_8 - dt*this%omega(i))
        end do
    class default
        call parcomm_global%abort("iomega solve-operator failure: rhs of wrong type")
    end select
    class default
        call parcomm_global%abort("iomega solve-operator failure: vout of wrong type")
    end select
end subroutine solve_iomega

end module operator_iomega_mod
