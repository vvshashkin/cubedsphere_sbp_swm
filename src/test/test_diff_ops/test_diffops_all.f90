program test_diffops

use parcomm_mod,      only : init_global_parallel_enviroment, &
                             deinit_global_parallel_enviroment, &
                             parcomm_global
use test_diffops_mod, only: test_div, test_grad, test_conv, test_curl,      &
                            test_coriolis, test_KE, test_coriolis_vec_inv,  &
                            test_co2contra, test_grad_perp, test_contra2co, &
                            test_vec_advection

use key_value_mod,    only : key_value_r8_t

implicit none

real(kind=8) :: err
type(key_value_r8_t)  :: errs

call init_global_parallel_enviroment()

errs = test_div(N=32,div_oper_name="divergence_a2_ecs",staggering="A")
if(parcomm_global%myid == 0) then
    print *, "divergence_a2_ecs"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_div(N=32,div_oper_name="divergence_a2_cons",staggering="A")
if(parcomm_global%myid == 0) then
    print *, "divergence_a2_cons"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_div(N=32,div_oper_name="divergence_a2_fv",staggering="A")
if(parcomm_global%myid == 0) then
    print *, "divergence_a2_fv"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_div(N=32,div_oper_name="divergence_c2",staggering="C")
if(parcomm_global%myid == 0) then
    print *, "divergence_c2"
    print "(A,5E15.7)", "Err: ", errs%values
end if

errs = test_div(N=32,div_oper_name="divergence_ch_sbp21",staggering="Ch")
if(parcomm_global%myid == 0) then
    print *, "divergence_ch_sbp21"
    print "(A,5E15.7)", "Err: ", errs%values
end if

errs = test_div(N=32,div_oper_name="divergence_ch_sbp42",staggering="Ch")
if(parcomm_global%myid == 0) then
    print *, "divergence_ch_sbp42"
    print "(A,5E15.7)", "Err: ", errs%values
end if

errs = test_div(N=32,div_oper_name="divergence_c_sbp21",staggering="C")
if(parcomm_global%myid == 0) then
    print *, "divergence_c_sbp21"
    print "(A,5E15.7)", "Err: ", errs%values
end if

errs = test_div(N=32,div_oper_name="divergence_c_sbp42",staggering="C")
if(parcomm_global%myid == 0) then
    print *, "divergence_c_sbp42"
    print "(A,5E15.7)", "Err: ", errs%values
end if

errs = test_div(N=32,div_oper_name="divergence_ah2",staggering="Ah")
if(parcomm_global%myid == 0) then
    print *, "divergence_ah2"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_div(N=32,div_oper_name="divergence_ah42_sbp",staggering="Ah")
if(parcomm_global%myid == 0) then
    print *, "divergence_ah42_sbp"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_div(N=32,div_oper_name="divergence_ah43_sbp",staggering="Ah")
if(parcomm_global%myid == 0) then
    print *, "divergence_ah43_sbp"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_grad(N=32,grad_oper_name="gradient_a2_ecs",staggering="A")
if(parcomm_global%myid == 0) then
    print *, "gradient_a2_ecs"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_grad(N=32,grad_oper_name="gradient_a2_cons",staggering="A")
if(parcomm_global%myid == 0) then
    print *, "gradient_a2_cons"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_grad(N=32,grad_oper_name="gradient_ch_sbp21",staggering="Ch")
if(parcomm_global%myid == 0) then
    print *, "gradient_ch_sbp21"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_grad(N=32,grad_oper_name="gradient_ch_sbp42",staggering="Ch")
if(parcomm_global%myid == 0) then
    print *, "gradient_ch_sbp42"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_grad(N=32,grad_oper_name="gradient_ah21_sbp_ecs",staggering="Ah")
if(parcomm_global%myid == 0) then
    print *, "gradient_ah21_sbp_ecs"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_grad(N=32,grad_oper_name="gradient_ah42_sbp_ecs",staggering="Ah")
if(parcomm_global%myid == 0) then
    print *, "gradient_ah42_sbp_ecs"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_grad(N=32,grad_oper_name="gradient_ah43_sbp_ecs",staggering="Ah")
if(parcomm_global%myid == 0) then
    print *, "gradient_ah43_sbp_ecs"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_grad(N=32,grad_oper_name="gradient_c2_ecs",staggering="C")
if(parcomm_global%myid == 0) then
    print *, "gradient_c2_ecs"
    print "(A,4E15.7)", "Err: ", errs%values
end if

! errs = test_grad(N=32,grad_oper_name="gradient_c2_cons",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_c2_cons"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

errs = test_grad(N=32,grad_oper_name="gradient_c_sbp21",staggering="C")
if(parcomm_global%myid == 0) then
    print *, "gradient_c_sbp21"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_grad(N=32,grad_oper_name="gradient_c_sbp42",staggering="C")
if(parcomm_global%myid == 0) then
    print *, "gradient_c_sbp42"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_curl(N=32,curl_oper_name="curl_divergence_ah42_sbp",staggering="Ah")
if(parcomm_global%myid == 0) then
    print *, "curl_divergence_ah42_sbp"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_curl(N=32,curl_oper_name="curl_divergence_ah43_sbp",staggering="Ah")
if(parcomm_global%myid == 0) then
    print *, "curl_divergence_ah43_sbp"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_curl(N=32,curl_oper_name="curl_c_sbp42",staggering="C")
if(parcomm_global%myid == 0) then
    print *, "curl_c_sbp42"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_curl(N=32,curl_oper_name="curl_c_sbp21",staggering="C")
if(parcomm_global%myid == 0) then
    print *, "curl_c_sbp21"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_coriolis(N=32, coriolis_op_name="coriolis_colocated", staggering="A")
if (parcomm_global%myid==0) then
    print *, "coriolis_colocated"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_coriolis_vec_inv(N=32, coriolis_op_name="coriolis_Cgrid_sbp42", staggering="C")
if (parcomm_global%myid==0) then
    print *, "coriolis_Cgrid_sbp42"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_coriolis_vec_inv(N=32, coriolis_op_name="coriolis_Cgrid_sbp21", staggering="C")
if (parcomm_global%myid==0) then
    print *, "coriolis_Cgrid_sbp21"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_KE(N=32, KE_oper_name="KE_colocated", staggering="A")
if (parcomm_global%myid==0) then
    print *, "KE_colocated"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_KE(N=32, KE_oper_name="KE_Cgrid_sbp42", staggering="C")
if (parcomm_global%myid==0) then
    print *, "KE_Cgrid_sbp42"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_KE(N=32, KE_oper_name="KE_Cgrid_sbp21", staggering="C")
if (parcomm_global%myid==0) then
    print *, "KE_Cgrid_sbp21"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_co2contra(N=32,co2contra_oper_name="co2contra_c_sbp21", staggering="C")
if (parcomm_global%myid==0) then
    print *, "co2contra_c_sbp21"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_co2contra(N=32,co2contra_oper_name="co2contra_c_sbp21_new", staggering="C")
if (parcomm_global%myid==0) then
    print *, "co2contra_c_sbp21_new"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_contra2co(N=32,contra2co_oper_name="co2contra_c_sbp21_new", staggering="C")
if (parcomm_global%myid==0) then
    print *, "contra2co_c_sbp21_new"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_co2contra(N=32,co2contra_oper_name="co2contra_c_sbp42", staggering="C")
if (parcomm_global%myid==0) then
    print *, "co2contra_c_sbp42"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_co2contra(N=32,co2contra_oper_name="co2contra_c_sbp42_new", staggering="C")
if (parcomm_global%myid==0) then
    print *, "co2contra_c_sbp42_new"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_co2contra(N=32,co2contra_oper_name="co2contra_ch_sbp21", staggering="Ch")
if (parcomm_global%myid==0) then
    print *, "co2contra_ch_sbp21"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_co2contra(N=32,co2contra_oper_name="co2contra_ch_sbp42", staggering="Ch")
if (parcomm_global%myid==0) then
    print *, "co2contra_ch_sbp42"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_co2contra(N=32,co2contra_oper_name="co2contra_colocated", staggering="Ah")
if (parcomm_global%myid==0) then
    print *, "co2contra_colocated Ah grid"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_contra2co(N=32,contra2co_oper_name="co2contra_c_sbp42_new", staggering="C")
if (parcomm_global%myid==0) then
    print *, "contra2co_c_sbp42_new"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_contra2co(N=32,contra2co_oper_name="co2contra_colocated", staggering="Ah")
if (parcomm_global%myid==0) then
    print *, "contra2co_colocated"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_C_up4", staggering="C")
if (parcomm_global%myid==0) then
    print *, "vector_advection_C_up4"
    print "(A,4E15.7)", "Err: ", errs%values
end if
errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_C_up3", staggering="C")
if (parcomm_global%myid==0) then
    print *, "vector_advection_C_up3"
    print "(A,4E15.7)", "Err: ", errs%values
end if
errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_C_up1", staggering="C")
if (parcomm_global%myid==0) then
    print *, "vector_advection_C_up1"
    print "(A,4E15.7)", "Err: ", errs%values
end if
errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_C_c2", staggering="C")
if (parcomm_global%myid==0) then
    print *, "vector_advection_C_c2"
    print "(A,4E15.7)", "Err: ", errs%values
end if
errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_C_c4", staggering="C")
if (parcomm_global%myid==0) then
    print *, "vector_advection_C_c4"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_Ah21", staggering="Ah")
if (parcomm_global%myid==0) then
    print *, "vector_advection_Ah21"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_Ah42", staggering="Ah")
if (parcomm_global%myid==0) then
    print *, "vector_advection_Ah42"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_Ah63", staggering="Ah")
if (parcomm_global%myid==0) then
    print *, "vector_advection_Ah63"
    print "(A,4E15.7)", "Err: ", errs%values
end if

! errs = test_grad_perp(N=32, grad_perp_oper_name="grad_perp_c_sbp42", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "grad_perp_c_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

!
! errs = test_KE(N=32, KE_oper_name="KE_colocated", staggering="A")
! if (parcomm_global%myid==0) then
!     print *, "KE_colocated"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_KE(N=32, KE_oper_name="KE_Cgrid", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "KE_Cgrid"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

call deinit_global_parallel_enviroment()

end program test_diffops
