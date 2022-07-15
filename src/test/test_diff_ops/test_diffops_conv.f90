program test_diffops

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment, &
                                parcomm_global
use test_diffops_mod, only: test_div, test_grad, test_conv, test_curl

implicit none

integer(kind=4), parameter :: Ns(3) = [32,64,128]

call init_global_parallel_enviroment()

call test_conv(operator_name="gradient_c_sbp21",staggering="C",Ns=Ns)
call test_conv(operator_name="gradient_c_sbp42",staggering="C",Ns=Ns)
call test_conv(operator_name="gradient_ch_sbp21",staggering="Ch",Ns=Ns)
call test_conv(operator_name="gradient_ch_sbp42",staggering="Ch",Ns=Ns)
call test_conv(operator_name="gradient_ah_c21_sbp_ecs",staggering="Ah_C",Ns=Ns)
call test_conv(operator_name="gradient_ah21_sbp_ecs",staggering="Ah",Ns=Ns)
call test_conv(operator_name="gradient_ah42_sbp_ecs",staggering="Ah",Ns=Ns)
call test_conv(operator_name="gradient_ah63_sbp_ecs",staggering="Ah",Ns=Ns)
call test_conv(operator_name="gradient_ah43_sbp_ecs",staggering="Ah",Ns=Ns)
call test_conv(operator_name="divergence_c_sbp21",staggering="C",Ns=Ns)
call test_conv(operator_name="divergence_c_sbp42",staggering="C",Ns=Ns)
call test_conv(operator_name="divergence_ch_sbp21",staggering="Ch",Ns=Ns)
call test_conv(operator_name="divergence_ch_sbp42",staggering="Ch",Ns=Ns)
call test_conv(operator_name="divergence_ah_c_sbp21",staggering="Ah_C",Ns=Ns)
call test_conv(operator_name="divergence_ah42_sbp",staggering="Ah",Ns=Ns)
call test_conv(operator_name="divergence_ah63_sbp",staggering="Ah",Ns=Ns)
call test_conv(operator_name="divergence_ah43_sbp",staggering="Ah",Ns=Ns)
call test_conv(operator_name="curl_divergence_ah42_sbp",staggering="Ah",Ns=Ns)
call test_conv(operator_name="curl_c_sbp21",staggering="C",Ns=Ns)
call test_conv(operator_name="curl_c_sbp42",staggering="C",Ns=Ns)
call test_conv(operator_name="curl_divergence_ah43_sbp",staggering="Ah",Ns=Ns)

call deinit_global_parallel_enviroment()

end program test_diffops
