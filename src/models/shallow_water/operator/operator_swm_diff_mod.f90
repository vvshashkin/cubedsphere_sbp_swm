module operator_swm_diff_mod

use stvec_mod,      only : stvec_t
use domain_mod,     only : domain_t
use operator_mod,   only : operator_t

use abstract_hordiff_mod,    only : hordiff_operator_t

use stvec_swm_mod, only : stvec_swm_t
use parcomm_mod,   only : parcomm_global

implicit none

type, public, extends(operator_t) :: operator_swm_diff_t

    class(hordiff_operator_t),   allocatable :: hordiff_uv
    class(hordiff_operator_t),   allocatable :: hordiff

contains
    procedure, public :: apply
end type operator_swm_diff_t

contains

subroutine apply(this, vout, vin, domain)
    class(operator_swm_diff_t), intent(inout) :: this
    class(stvec_t),             intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),             intent(inout) :: vin
    type(domain_t),             intent(in)    :: domain

    select type (vout)
    class is (stvec_swm_t)
        select type (vin)
        class is (stvec_swm_t)

            call this%hordiff_uv%calc_diff_vec(vout%u, vout%v, vin%u, vin%v, domain)
            call this%hordiff%calc_diff(vout%h, vin%h, domain%mesh_p, domain)

        class default
            call parcomm_global%abort("swm operator failure: vin of wrong type")
        end select
    class default
        call parcomm_global%abort("swm operator failure: vout of wrong type")
    end select

end subroutine apply

end module operator_swm_diff_mod
