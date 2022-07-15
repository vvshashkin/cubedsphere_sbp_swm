module nh_model_factory_mod

use nh_model_mod,                  only : nh_model_t
use config_nh_model_mod,           only : config_nh_model_t
use domain_mod,                    only : domain_t
use domain_factory_mod,            only : create_domain
use stvec_nh_factory_mod,          only : create_stvec_nh
use timescheme_factory_mod,        only : create_timescheme
use nh_testcases_mod,              only : get_initial_conditions
use postprocessing_nh_factory_mod, only : create_nh_postprocessing
use nh_operator_factory_mod,       only : create_nh_operator

implicit none

contains

subroutine create_nh_model(nh_model, config)

    !WORKAROUND
    integer(kind=4), parameter :: halo_width_hor = 8, halo_width_ver = 0

    class(nh_model_t),       intent(out)   :: nh_model
    type(config_nh_model_t), intent(in)    :: config

    call create_domain(nh_model%domain, config%config_domain)

    call create_stvec_nh(nh_model%stvec, nh_model%domain, halo_width_hor, halo_width_ver)
    call get_initial_conditions(nh_model%stvec, nh_model%domain, config%testcase_name)

    call create_timescheme(nh_model%timescheme, nh_model%stvec, config%timescheme_name)

    call create_nh_operator(nh_model%operator, config%config_operator, nh_model%domain)

    nh_model%simulation_time = config%simulation_time
    nh_model%dt              = config%dt
    nh_model%tau_write       = config%tau_write
    nh_model%tau_diagnostics = config%tau_diagnostics

    call create_nh_postprocessing(nh_model%postproc, config%config_postprocessing, &
                                  nh_model%domain)

end subroutine create_nh_model

end module nh_model_factory_mod
