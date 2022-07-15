module postprocessing_nh_factory_mod

use abstract_postprocessing_mod, only : postprocessing_t
use simple_postnh_mod,           only : simple_postnh_t
use config_postnh_mod,           only : config_postnh_t
use domain_mod,                  only : domain_t
use outputer_factory_mod,        only : create_latlon_outputer,      &
                                        create_latlon_vec_outputer,  &
                                        create_master_paneled_outputer
use parcomm_mod,                 only : parcomm_global

implicit none

contains

subroutine create_nh_postprocessing(post,config,domain)
    class(postprocessing_t), allocatable, intent(out) :: post
    type(config_postnh_t), intent(in) :: config
    type(domain_t),        intent(in) :: domain

    type(simple_postnh_t), allocatable :: post_loc

    allocate(post_loc)
    post_loc%Nlon = config%Nlon
    post_loc%Nlat = config%Nlat

    select case(config%outputer_name)
    case("latlon")
        if(domain%horizontal_staggering == "C") then
            call create_latlon_outputer(post_loc%outputer_theta,config%Nlat,config%Nlon, &
                                        "A",domain,is_z_interfaces = .true.)
            call create_latlon_outputer(post_loc%outputer_P,config%Nlat,config%Nlon, &
                                        "A",domain,is_z_interfaces = .false.)
            call create_latlon_vec_outputer(post_loc%outputer_uv,  config%Nlat, config%Nlon, "C", &
                                            "contravariant", domain)
        else if(domain%horizontal_staggering == "Ah") then
            call create_latlon_outputer(post_loc%outputer_theta,config%Nlat,config%Nlon, &
                                        "Ah", domain,is_z_interfaces=.true.)
            call create_latlon_outputer(post_loc%outputer_P,config%Nlat,config%Nlon, &
                                        "Ah", domain,is_z_interfaces=.false.)
            call create_latlon_vec_outputer(post_loc%outputer_uv,  config%Nlat, config%Nlon, "Ah", &
                                            "contravariant", domain)
        end if
    case("paneled")
        call create_master_paneled_outputer(post_loc%outputer_theta, "w", domain)
        call create_master_paneled_outputer(post_loc%outputer_P, "p", domain)
    case default
        call parcomm_global%abort("create_nh_postprocessing, unknown outputer: "//&
                                   config%outputer_name)
    end select

    call move_alloc(post_loc,post)
end subroutine create_nh_postprocessing

end module postprocessing_nh_factory_mod
