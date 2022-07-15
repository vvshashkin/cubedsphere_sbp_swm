module massflux_factory_mod

use domain_mod,             only : domain_t
use abstract_massflux_mod,  only : massflux_operator_t
use massflux_colocated_mod, only : massflux_colocated_t
use massflux_Cgrid_mod,     only : massflux_chalo_t, massflux_c_sbp21_t, massflux_c_sbp42_t, &
                                   massflux_c_up4_t
use parcomm_mod,            only : parcomm_global
use halo_factory_mod,       only : create_halo_procedure, create_vector_halo_procedure
use sbp_factory_mod,        only : create_sbp_operator


implicit none

contains

function create_massflux_operator(domain, massflux_operator_name) result(massflux)
    type(domain_t),    intent(in)    :: domain
    character(len=*),  intent(in)    :: massflux_operator_name

    class(massflux_operator_t), allocatable :: massflux

    if(massflux_operator_name == 'massflux_colocated') then
        massflux = massflux_colocated_t()
    elseif(massflux_operator_name == 'massflux_c2') then
        massflux = massflux_chalo_t(order=2)
        select type(massflux)
        type is (massflux_chalo_t)
            call create_halo_procedure(massflux%halo,domain,2,"ECS_O")
            call create_vector_halo_procedure(massflux%halo_flux,domain,2,"C_vec_default")
        end select
    elseif(massflux_operator_name == 'massflux_c4') then
        massflux = massflux_chalo_t(order=4)
        select type(massflux)
        type is (massflux_chalo_t)
            call create_halo_procedure(massflux%halo,domain,2,"ECS_O")
            call create_vector_halo_procedure(massflux%halo_flux,domain,2,"C_vec_default")
        end select
    elseif(massflux_operator_name == 'massflux_c_sbp21') then
        massflux = massflux_c_sbp21_t()
    elseif(massflux_operator_name == 'massflux_c_sbp42') then
        massflux = create_c_sbp42_massflux()
    elseif(massflux_operator_name == 'massflux_c_up4') then
        massflux = massflux_c_up4_t()
        select type(massflux)
        type is (massflux_c_up4_t)
            call create_halo_procedure(massflux%halo,domain,3,"ECS_O")
            call create_vector_halo_procedure(massflux%halo_flux,domain,2,"C_vec_default")
        end select
    else
        call parcomm_global%abort("unknown massflux operator: "//massflux_operator_name)
    end if

end function create_massflux_operator

function create_c_sbp42_massflux() result(massflux)
    type(massflux_c_sbp42_t) massflux

    massflux%sbp_interp_h2v = create_sbp_operator("W42_stagered_interp_c2i")

end function create_c_sbp42_massflux

end module massflux_factory_mod
