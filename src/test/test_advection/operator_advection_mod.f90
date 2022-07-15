module operator_advection_mod

use stvec_mod,             only : stvec_t
use domain_mod,            only : domain_t
use operator_mod,          only : operator_t
use grid_field_mod,        only : grid_field_t
use abstract_div_mod,      only : div_operator_t
use abstract_massflux_mod, only : massflux_operator_t
use stvec_advection_mod,   only : stvec_advection_t
use parcomm_mod,           only : parcomm_global

implicit none

type, public, extends(operator_t) :: operator_advection_t
    class(div_operator_t),      allocatable :: div_op
    class(massflux_operator_t), allocatable :: flux_op
    type(grid_field_t) :: hu, hv
contains
    procedure, public :: apply
end type operator_advection_t

contains

subroutine apply(this, vout, vin, domain)
    class(operator_advection_t), intent(inout) :: this
    class(stvec_t),              intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),              intent(inout) :: vin
    type(domain_t),              intent(in)    :: domain

    select type (vout)
    class is (stvec_advection_t)
        select type (vin)
        class is (stvec_advection_t)

            !call this%hu%assign_prod(-1.0_8, vin%h, vin%u, domain%mesh_p)
            !call this%hv%assign_prod(-1.0_8, vin%h, vin%v, domain%mesh_p)
            call this%flux_op%calc_massflux(this%hu, this%hv, vin%h, vin%u, vin%v, domain)

            call this%div_op%calc_div(vout%h, this%hu, this%hv, domain)

            call vout%h%assign(-1.0_8,vout%h,domain%mesh_p)
            call vout%u%assign(0.0_8, domain%mesh_u)
            call vout%v%assign(0.0_8, domain%mesh_v)

        class default
            call parcomm_global%abort("Advection operator failure: vin of wrong type")
        end select
    class default
        call parcomm_global%abort("Advection operator failure: vout of wrong type")
    end select
end subroutine apply

end module operator_advection_mod
