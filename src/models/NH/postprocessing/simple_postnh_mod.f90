module simple_postnh_mod

use domain_mod,                  only : domain_t
use stvec_mod,                   only : stvec_t
use abstract_postprocessing_mod, only : postprocessing_t
use outputer_abstract_mod,       only : outputer_t, outputer_vector_t
use stvec_nh_mod,                only : stvec_nh_t
use parcomm_mod,                 only : parcomm_global

implicit none

type, extends(postprocessing_t) :: simple_postnh_t
    integer(kind=4) :: Nlon, Nlat
    class(outputer_t), allocatable :: outputer_theta
    class(outputer_t), allocatable :: outputer_P
    class(outputer_vector_t), allocatable :: outputer_uv
contains
    procedure :: write_fields
end type simple_postnh_t

contains
subroutine write_fields(this, irec, stvec, domain)
    class(simple_postnh_t),  intent(inout) :: this
    integer(kind=4),         intent(in)    :: irec
    class(stvec_t),          intent(inout) :: stvec
    type(domain_t),          intent(in)    :: domain

    select type(stvec)
    class is (stvec_nh_t)
        call this%outputer_P%write(stvec%P, domain, 'P.dat', irec)
        call this%outputer_theta%write(stvec%theta, domain, 'theta.dat', irec)
        call this%outputer_theta%write(stvec%eta_dot, domain, 'w.dat', irec)
        call this%outputer_uv%write(stvec%u, stvec%v, domain, 'u.dat', 'v.dat', irec)
    class default
        call parcomm_global%abort("unsupported state vector type in simple_postnh_t%write")
    end select
end subroutine write_fields

end module simple_postnh_mod
