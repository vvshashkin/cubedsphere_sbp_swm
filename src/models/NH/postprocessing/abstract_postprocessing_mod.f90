module abstract_postprocessing_mod

use domain_mod, only : domain_t
use stvec_mod,  only : stvec_t

implicit none

type, abstract :: postprocessing_t
contains
    procedure(write_fields), deferred :: write_fields
end type postprocessing_t

abstract interface
    subroutine write_fields(this, irec, stvec, domain)
        import domain_t, stvec_t, postprocessing_t
        class(postprocessing_t), intent(inout) :: this
        integer(kind=4),         intent(in)    :: irec
        class(stvec_t),          intent(inout) :: stvec
        type(domain_t),          intent(in)    :: domain
    end subroutine write_fields
end interface

end module abstract_postprocessing_mod
