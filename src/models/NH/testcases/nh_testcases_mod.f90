module nh_testcases_mod

use parcomm_mod,  only : parcomm_global
use domain_mod,   only : domain_t
use stvec_mod,    only : stvec_t

use NH_GW_testcase_mod,             only : get_GW_initial_conditions
use NH_solid_rotation_testcase_mod, only : get_nh_solid_rotation_initial_conditions
use advection3d_testcases_mod,      only : get_advection3d_initial_conditions
use Straka_testcase_mod,            only : get_Straka_initial_conditions
use hot_bubble_testcase_mod,        only : get_hot_bubble_initial_conditions

implicit none

contains

subroutine get_initial_conditions(stvec, domain, testcase_name)
    class(stvec_t),   intent(inout) :: stvec
    type(domain_t),   intent(in)    :: domain
    character(len=*), intent(in)    :: testcase_name

    select case(testcase_name)
    case("GW", "GW_linear")
        call get_GW_initial_conditions(stvec,domain,testcase_name)
    case("solid_rotation")
        call get_nh_solid_rotation_initial_conditions(stvec,domain)
    case("Straka","straka")
        call get_Straka_initial_conditions(stvec,domain)
    case("hot_bubble")
        call get_hot_bubble_initial_conditions(stvec,domain)
    case("advection3d_solid_rotation","Straka_buble")
        call get_advection3d_initial_conditions(stvec,domain,testcase_name)
    case default
        call parcomm_global%abort("unknown NH testcases name: "//testcase_name)
    end select
end subroutine get_initial_conditions

end module nh_testcases_mod
