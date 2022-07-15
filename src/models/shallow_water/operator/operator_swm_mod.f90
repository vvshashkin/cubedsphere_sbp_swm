module operator_swm_mod

use stvec_mod,      only : stvec_t
use domain_mod,     only : domain_t
use operator_mod,   only : operator_t
use grid_field_mod, only : grid_field_t

use abstract_div_mod,        only : div_operator_t
use abstract_grad_mod,       only : grad_operator_t
use abstract_coriolis_mod,   only : coriolis_operator_t
use abstract_curl_mod,       only : curl_operator_t
use abstract_KE_mod,         only : KE_operator_t
use abstract_massflux_mod,   only : massflux_operator_t
use abstract_co2contra_mod,  only : co2contra_operator_t
use abstract_quadrature_mod, only : quadrature_t

use stvec_swm_mod, only : stvec_swm_t
use parcomm_mod,   only : parcomm_global

use key_value_mod, only : key_value_r8_t

implicit none

type, public, extends(operator_t) :: operator_swm_t

    class(div_operator_t),       allocatable :: div_op
    class(grad_operator_t),      allocatable :: grad_op
    class(coriolis_operator_t),  allocatable :: coriolis_op
    class(curl_operator_t),      allocatable :: curl_op
    class(KE_operator_t),        allocatable :: KE_op
    class(massflux_operator_t),  allocatable :: massflux_op
    class(co2contra_operator_t), allocatable :: co2contra_op
    class(quadrature_t),         allocatable :: quadrature_h, quadrature_u, &
                                                quadrature_v, quadrature_w

    !work fields for operator
    type(grid_field_t) :: h_surf !orography
    type(grid_field_t) :: div, grad_x, grad_y, curl
    type(grid_field_t) :: cor_u, cor_v
    type(grid_field_t) :: KE !kinetic energy
    type(grid_field_t) :: ut, vt !contravariant components
    type(grid_field_t) :: hu, hv !mass fluxes in continuty eq.

    !work fields for diagnostics
    type(grid_field_t) :: KE_diag_u, KE_diag_v !kinetic energy
    type(grid_field_t) :: PE_diag, hu_diag, hv_diag !kinetic energy

!gravity acceleration. Should be moved to mesh?
    real(kind=8) :: grav

contains
    procedure, public :: apply
    procedure, public :: get_diagnostics
    procedure, public :: get_diagnostics_tend
    procedure, public :: calc_energy
    procedure, public :: calc_enstrophy
end type operator_swm_t

contains

subroutine apply(this, vout, vin, domain)
    class(operator_swm_t), intent(inout) :: this
    class(stvec_t),        intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),        intent(inout) :: vin
    type(domain_t),        intent(in)    :: domain

    !uncomment for enegry tendency diagnostics
    !type(key_value_r8_t) :: te_tend

    real(kind=8) :: ke_u, ke_v, pe, te

    select type (vout)
    class is (stvec_swm_t)
        select type (vin)
        class is (stvec_swm_t)

            call this%co2contra_op%transform(this%ut, this%vt, vin%u, vin%v, domain)
            call this%massflux_op%calc_massflux(this%hu, this%hv, &
                         vin%h, this%ut, this%vt, domain)

            !momentum eq part
            call this%KE_op%calc_KE(this%KE, vin%u, vin%v, this%ut, this%vt, domain)

            !store grav*grad(h+hs+KE) in KE array
            call this%KE%update(this%grav, vin%h, this%grav, this%h_surf, domain%mesh_p)

            call this%grad_op%calc_grad(this%grad_x, this%grad_y, this%KE, domain)

            call this%curl_op%calc_curl(this%curl, vin%u, vin%v, domain)

            call this%coriolis_op%calc_coriolis_vec_inv(this%cor_u, this%cor_v, &
                                            this%hu, this%hv, vin%h, this%curl, domain)


            call vout%u%assign(-1.0_8, this%grad_x, 1.0_8, this%cor_u, domain%mesh_u)
            call vout%v%assign(-1.0_8, this%grad_y, 1.0_8, this%cor_v, domain%mesh_v)

            !continuty eq part
            call this%div_op%calc_div(this%div, this%hu, this%hv, domain)

            call vout%h%assign(-1.0_8, this%div, domain%mesh_p)

        class default
            call parcomm_global%abort("swm operator failure: vin of wrong type")
        end select
    class default
        call parcomm_global%abort("swm operator failure: vout of wrong type")
    end select

    !uncomment for enegry tendency diagnostics
    !te_tend = this%get_diagnostics_tend(vin, vout, domain)
    !call te_tend%print()
end subroutine apply

function get_diagnostics(this, v, domain) result(diagnostics)

    use key_value_mod, only : key_value_r8_t

    class(operator_swm_t), intent(inout) :: this
    class(stvec_t),        intent(inout) :: v
    type(domain_t),        intent(in)    :: domain

    type(key_value_r8_t)  :: diagnostics

    integer(kind=4), parameter :: ndiag = 7
    real(kind=8) :: te, ke, pe, enstrophy

    allocate(diagnostics%keys(ndiag))
    allocate(diagnostics%values(ndiag))

    select type(v)
    type is (stvec_swm_t)
        diagnostics%keys(1)%str = "hmin"
        diagnostics%values(1) = v%h%minimum(domain%mesh_p, domain%parcomm)
        diagnostics%keys(2)%str = "hmax"
        diagnostics%values(2) = v%h%maximum(domain%mesh_p, domain%parcomm)
        diagnostics%keys(3)%str = "mass"

        diagnostics%values(3) = this%quadrature_h%mass(v%h, domain%mesh_p, domain%parcomm)
        call this%calc_energy(te,ke,pe,v, domain)
        diagnostics%keys(4)%str = "TE"
        diagnostics%values(4) = te
        diagnostics%keys(5)%str = "KE"
        diagnostics%values(5) = ke
        diagnostics%keys(6)%str = "PE"
        diagnostics%values(6) = pe

        call this%calc_enstrophy(enstrophy,v, domain)
        diagnostics%keys(7)%str = "Enstrophy"
        diagnostics%values(7) = enstrophy
    class default
        call parcomm_global%abort("wrong type of v in operator_swm_t get_diagnostics")
    end select

end function get_diagnostics

function get_diagnostics_tend(this, v, vtend, domain) result(diagnostics)

    use key_value_mod, only : key_value_r8_t

    class(operator_swm_t), intent(inout) :: this
    class(stvec_t),        intent(inout) :: v, vtend
    type(domain_t),        intent(in)    :: domain

    type(key_value_r8_t)  :: diagnostics

    real(kind=8) :: ke_u_tend, ke_v_tend, pe_tend, te, ke, pe

    allocate(diagnostics%keys(2))
    allocate(diagnostics%values(2))

    select type(v)
    type is (stvec_swm_t)
    select type(vtend)
    type is (stvec_swm_t)

        call this%co2contra_op%transform(this%ut, this%vt, v%u, v%v, domain)

        call this%massflux_op%calc_massflux(this%hu_diag, this%hv_diag, &
                                            vtend%h, this%ut, this%vt, domain)

        call this%hu_diag%assign_prod(0.5_8, this%hu_diag, v%u, domain%mesh_u)
        call this%hv_diag%assign_prod(0.5_8, this%hv_diag, v%v, domain%mesh_v)

        call this%PE_diag%assign_prod(this%grav, v%h, vtend%h, domain%mesh_p)

        call this%KE_diag_u%assign_prod(1.0_8, this%hu, vtend%u, domain%mesh_u)
        call this%KE_diag_v%assign_prod(1.0_8, this%hv, vtend%v, domain%mesh_v)

        ke_u_tend = this%quadrature_u%mass(this%KE_diag_u, domain%mesh_u, domain%parcomm) +&
                    this%quadrature_u%mass(this%hu_diag, domain%mesh_u, domain%parcomm)
        ke_v_tend = this%quadrature_v%mass(this%KE_diag_v, domain%mesh_v, domain%parcomm) +&
                    this%quadrature_v%mass(this%hv_diag, domain%mesh_v, domain%parcomm)
        pe_tend   = this%quadrature_h%mass(this%PE_diag, domain%mesh_p, domain%parcomm)

        diagnostics%keys(1)%str = "total energy tendency"
        diagnostics%values(1)   = ke_u_tend+ke_v_tend+pe_tend

        call this%calc_energy(te,ke,pe,v, domain)
        diagnostics%keys(2)%str = "total energy tendency / te"
        diagnostics%values(2)   = (ke_u_tend+ke_v_tend+pe_tend) / te
    class default
        call parcomm_global%abort("wrong type of vtend in operator_swm_t get_diagnostics")
    end select
    class default
        call parcomm_global%abort("wrong type of v in operator_swm_t get_diagnostics")
    end select

end function get_diagnostics_tend

subroutine calc_energy(this, te, ke, pe, vin, domain)


    class(operator_swm_t), intent(inout) :: this
    class(stvec_swm_t),    intent(inout) :: vin
    type(domain_t),        intent(in)    :: domain

    real(kind=8), intent(out) :: te, ke, pe

    call this%co2contra_op%transform(this%ut, this%vt, vin%u, vin%v, domain)
    call this%massflux_op%calc_massflux(this%hu, this%hv, &
                                             vin%h, this%ut, this%vt, domain)
    call this%massflux_op%calc_massflux(this%hu, this%hv, &
                                             vin%h, this%ut, this%vt, domain)

    call this%KE_diag_u%assign_prod(0.5_8,this%hu,vin%u,domain%mesh_u)
    call this%KE_diag_v%assign_prod(0.5_8,this%hv,vin%v,domain%mesh_v)
    call this%PE_diag%assign(1.0_8,vin%h,1.0_8,this%h_surf,domain%mesh_p)
    call this%PE_diag%assign_prod(0.5_8*this%grav,this%PE_diag,this%PE_diag,domain%mesh_p)

    ke = this%quadrature_u%mass(this%KE_diag_u,domain%mesh_u,domain%parcomm)+&
         this%quadrature_v%mass(this%KE_diag_v,domain%mesh_v,domain%parcomm)
    pe   = this%quadrature_h%mass(this%PE_diag,domain%mesh_p,domain%parcomm)

    te = ke+pe
end subroutine calc_energy

subroutine calc_enstrophy(this, enstrophy, vin, domain)

    use coriolis_factory_mod, only : calc_coriolis_parameter

    !ATTENTION: will not work correctly for C-grid

    class(operator_swm_t), intent(inout) :: this
    class(stvec_swm_t),    intent(inout) :: vin
    type(domain_t),        intent(in)    :: domain

    real(kind=8), intent(out) :: enstrophy

    enstrophy = 0.0_8

    call this%curl_op%calc_curl(this%curl, vin%u, vin%v, domain)
    call calc_coriolis_parameter(this%PE_diag,domain%mesh_q)
    call this%curl%update(1.0_8,this%PE_diag,domain%mesh_q)
    call this%curl%assign_prod(0.5_8,this%curl, this%curl, domain%mesh_q)
    call this%curl%assign_ratio(1.0_8,this%curl, vin%h, domain%mesh_q)

    enstrophy = this%quadrature_w%mass(this%curl, domain%mesh_q, domain%parcomm)

end subroutine calc_enstrophy

end module operator_swm_mod
