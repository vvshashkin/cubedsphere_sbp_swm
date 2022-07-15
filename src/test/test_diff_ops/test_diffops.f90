!Use for research purposes
!Auto test is test_diffops_all.f90
program test_diffops

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment, parcomm_global
use test_diffops_mod, only: test_div, test_grad, test_conv, test_curl, test_grad_perp, &
                            test_coriolis, test_coriolis_vec_inv, test_curl_grad,      &
                            test_co2contra, test_compatibility, test_KE, &
                            test_vec_advection, test_grad_3d, test_div_3d,  &
                            test_co2contra_3d, test_laplace
use key_value_mod,    only : key_value_r8_t

implicit none

real(kind=8) :: err
type(key_value_r8_t)  :: errs
integer(kind=4), parameter :: Ns(3) = [32,64,128]

call init_global_parallel_enviroment()

! errs = test_laplace(N=32,laplace_oper_name="divgrad_laplace_c_sbp21",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "divgrad_laplace_c_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_laplace(N=32,laplace_oper_name="divgrad_laplace_c_sbp21",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "divgrad_laplace_c_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
! errs = test_laplace(N=64,laplace_oper_name="divgrad_laplace_c_sbp21",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "divgrad_laplace_c_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_laplace(N=32,laplace_oper_name="divgrad_laplace_c_sbp42",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "divgrad_laplace_c_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
! errs = test_laplace(N=64,laplace_oper_name="divgrad_laplace_c_sbp42",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "divgrad_laplace_c_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
! errs = test_laplace(N=32,laplace_oper_name="divgrad_laplace_ch_sbp42",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "divgrad_laplace_ch_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
! errs = test_laplace(N=32,laplace_oper_name="divgrad_laplace_ah_sbp21",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "divgrad_laplace_ah_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
! errs = test_laplace(N=64,laplace_oper_name="divgrad_laplace_ah_sbp21",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "divgrad_laplace_ah_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
! errs = test_laplace(N=32,laplace_oper_name="divgrad_laplace_ah_sbp42",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "divgrad_laplace_ah_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
! errs = test_laplace(N=64,laplace_oper_name="divgrad_laplace_ah_sbp42",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "divgrad_laplace_ah_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
! errs = test_laplace(N=32,laplace_oper_name="laplace_ch_halo2",staggering="Ch")
! if(parcomm_global%myid == 0) then
!     print *, "laplace_ch_halo2"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
! errs = test_laplace(N=32,laplace_oper_name="laplace_ch_halo4",staggering="Ch")
! if(parcomm_global%myid == 0) then
!     print *, "laplace_ch_halo4"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

errs = test_laplace(N=32,laplace_oper_name="laplace_ah_sbp63_narrow",staggering="Ah")
if(parcomm_global%myid == 0) then
    print *, "laplace_ah_sbp63_narrow"
    print "(A,4E15.7)", "Err: ", errs%values
end if
errs = test_laplace(N=64,laplace_oper_name="laplace_ah_sbp63_narrow",staggering="Ah")
if(parcomm_global%myid == 0) then
    print *, "laplace_ah_sbp63_narrow"
    print "(A,4E15.7)", "Err: ", errs%values
end if
errs = test_laplace(N=128,laplace_oper_name="laplace_ah_sbp63_narrow",staggering="Ah")
if(parcomm_global%myid == 0) then
    print *, "laplace_ah_sbp63_narrow"
    print "(A,4E15.7)", "Err: ", errs%values
end if

! errs = test_laplace(N=32,laplace_oper_name="laplace_ah_sbp21_narrow",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "laplace_Ah_sbp21_narrow"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_laplace(N=64,laplace_oper_name="laplace_ah_sbp21_narrow",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "laplace_Ah_sbp21_narrow"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_laplace(N=128,laplace_oper_name="laplace_ah_sbp21_narrow",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "laplace_Ah_sbp21_narrow"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

! errs = test_div(N=32,div_oper_name="divergence_a2_ecs",staggering="A")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_a2_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

! errs = test_div(N=32,div_oper_name="divergence_a2_cons",staggering="A")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_a2_cons"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_a2_fv",staggering="A")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_a2_fv"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_c2",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_c2"
!     print "(A,5E15.7)", "Err: ", errs%values
! end if

! errs = test_div(N=32,div_oper_name="divergence_c_sbp21",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_c_sbp21"
!     print "(A,5E15.7)", "Err: ", errs%values
! end if
! errs = test_div(N=32,div_oper_name="divergence_ch_sbp21",staggering="Ch")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_ch_sbp21"
!     print "(A,5E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_ch_sbp42",staggering="Ch")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_ch_sbp21"
!     print "(A,5E15.7)", "Err: ", errs%values
! end if
!
!
!
! errs = test_div(N=32,div_oper_name="divergence_c_sbp42",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_c_sbp42"
!     print "(A,5E15.7)", "Err: ", errs%values
! end if

! errs = test_div(N=32,div_oper_name="divergence_ah2",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_ah2"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_ah42_sbp",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_ah42_sbp"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_ah43_sbp",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_ah43_sbp"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_a2_ecs",staggering="A")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_a2_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_a2_cons",staggering="A")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_a2_cons"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_ah21_sbp_ecs",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_ah21_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_ah42_sbp_ecs",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_ah42_sbp_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_ah43_sbp_ecs",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_ah43_sbp_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

! errs = test_grad(N=32,grad_oper_name="gradient_c2_ecs",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_c2_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

! errs = test_grad(N=32,grad_oper_name="gradient_c2_cons",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_c2_cons"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_c_sbp21",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_c_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_c_sbp42",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_c_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_ch_sbp21",staggering="Ch")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_ch_sbp21_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_ch_sbp42",staggering="Ch")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_ch_sbp42_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_ch_ecs_halo2",staggering="Ch")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_ch_ecs_halo2"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_ch_ecs_halo4",staggering="Ch")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_ch_ecs_halo4"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

! errs = test_grad(N=32,grad_oper_name="gradient_ah_c42_sbp_ecs",staggering="Ah_C")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_ah_c21_sbp_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

! errs = test_curl(N=32,curl_oper_name="curl_divergence_ah2",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "curl_divergence_ah42_sbp"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_curl(N=32,curl_oper_name="curl_divergence_ah42_sbp",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "curl_divergence_ah42_sbp"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_curl(N=32,curl_oper_name="curl_c_sbp21",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "curl_c_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_curl(N=32,curl_oper_name="curl_c_sbp42",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "curl_c_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_coriolis(N=32, coriolis_op_name="coriolis_colocated", staggering="A")
! if (parcomm_global%myid==0) then
!     print *, "coriolis_colocated"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_coriolis_vec_inv(N=32, coriolis_op_name="coriolis_Cgrid_sbp21", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "coriolis_Cgrid_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_coriolis_vec_inv(N=32, coriolis_op_name="coriolis_Cgrid_sbp42", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "coriolis_Cgrid_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!

! errs = test_KE(N=32, KE_oper_name="KE_Cgrid_sbp42", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "KE_Cgrid_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_KE(N=32, KE_oper_name="KE_Cgrid_sbp21", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "KE_Cgrid_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

! errs = test_coriolis(N=32, coriolis_op_name="coriolis_Cgrid_noncons_sbp21", &
!                      staggering="C", result_components = "contravariant")
! if (parcomm_global%myid==0) then
!     print *, "coriolis_Cgrid_noncons_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_coriolis(N=32, coriolis_op_name="coriolis_Cgrid_noncons_sbp42", &
!                      staggering="C", result_components = "contravariant")
! if (parcomm_global%myid==0) then
!     print *, "coriolis_Cgrid_noncons_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

!
! errs = test_co2contra(N=32,co2contra_oper_name="co2contra_colocated",staggering="A")
! if(parcomm_global%myid == 0) then
!     print *, "co2contra_colocated, A-grid"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_co2contra(N=32,co2contra_oper_name="co2contra_colocated",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "co2contra_colocated, Ah-grid"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

! errs = test_co2contra(N=32,co2contra_oper_name="co2contra_c_sbp21_new",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "co2contra c sbp21, C-grid"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_co2contra(N=32,co2contra_oper_name="co2contra_c_sbp42_new",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "co2contra c sbp42, C-grid"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_co2contra(N=32,co2contra_oper_name="co2contra_ch_sbp21", staggering="Ch")
! if (parcomm_global%myid==0) then
!     print *, "co2contra_ch_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_co2contra(N=32,co2contra_oper_name="co2contra_ch_sbp42", staggering="Ch")
! if (parcomm_global%myid==0) then
!     print *, "co2contra_ch_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

! call test_conv(operator_name="gradient_c_sbp21",staggering="C",Ns=Ns)
! call test_conv(operator_name="divergence_c_sbp21",staggering="C",Ns=Ns)
! call test_conv(operator_name="gradient_c_sbp42",staggering="C",Ns=Ns)
! call test_conv(operator_name="divergence_c_sbp42",staggering="C",Ns=Ns)
! call test_conv(operator_name="divergence_c2",staggering="C",Ns=Ns)
! call test_conv(operator_name="divergence_ah42_sbp",staggering="Ah",Ns=Ns)
! call test_conv(operator_name="divergence_ah43_sbp",staggering="Ah",Ns=Ns)
! call test_conv(operator_name="curl_divergence_ah42_sbp",staggering="Ah",Ns=Ns)
! call test_conv(operator_name="curl_divergence_ah43_sbp",staggering="Ah",Ns=Ns)

! errs = test_curl_grad(N=32,grad_oper_name="gradient_ah42_sbp_ecs",&
!                            curl_oper_name="curl_divergence_ah42_sbp",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "curl of grad _ah42_sbp_ecs test"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
! errs = test_curl_grad(N=32,grad_oper_name="gradient_c_sbp21",&
!                            curl_oper_name="curl_c_sbp21",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "curl of grad _c_sbp42_ecs test"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_curl_grad(N=32,grad_oper_name="gradient_c_sbp42",&
!                            curl_oper_name="curl_c_sbp42",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "curl of grad _c_sbp42_ecs test"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

!
! errs = test_curl_grad(N=32,grad_oper_name="gradient_ah2_ecs",&
!                            div_oper_name="divergence_ah2",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "curl of grad _ah2_ecs test"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

! print *, "Compatibility of Ah21:"
! call test_compatibility(div_operator_name  ="divergence_ah2",  &
!                         grad_operator_name = "gradient_ah21_sbp_ecs",  &
!                         co2contra_operator_name = "co2contra_colocated", &
!                         quadrature_name = "SBP_Ah21_quadrature", staggering="Ah")
!
! print *, "Compatibility of Ah42:"
! call test_compatibility(div_operator_name  ="divergence_ah42_sbp",  &
!                         grad_operator_name = "gradient_ah42_sbp_ecs",  &
!                         co2contra_operator_name = "co2contra_colocated", &
!                         quadrature_name = "SBP_Ah42_quadrature", staggering="Ah")
!
! print *, "Compatibility of C21"
! call test_compatibility(div_operator_name  ="divergence_c_sbp21",  &
!                         grad_operator_name = "gradient_c_sbp21",  &
!                         co2contra_operator_name = "co2contra_c_sbp21", &
!                         quadrature_name = "SBP_C21_quadrature", staggering="C")
!
! print *, "Compatibility of C42"
! call test_compatibility(div_operator_name  ="divergence_c_sbp42",  &
!                         grad_operator_name = "gradient_c_sbp42",  &
!                         co2contra_operator_name = "co2contra_c_sbp42", &
!                         quadrature_name = "SBP_C42_quadrature", staggering="C")

! errs = test_grad_perp(N=32, grad_perp_oper_name="grad_perp_c_sbp21", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "grad_perp_c_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad_perp(N=32, grad_perp_oper_name="grad_perp_c_sbp42", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "grad_perp_c_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

! errs =  test_div_3d(Nh = 32, Nz = 20, &
!                     hor_div_name = "divergence_c_sbp42", &
!                     diff_eta_name = "eta_diff_w2p_sbp42", &
!                     horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "div_3d_c_sbp42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_grad_3d(Nh = 32, Nz = 8, &
!                     hor_grad_name = "gradient_c_sbp42", &
!                     diff_eta_name = "eta_diff_p2w_sbp42", &
!                     horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "grad_3d_c_sbp42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
! errs =  test_co2contra_3d(Nh = 32, Nz = 8, &
!                     co2contra_3d_oper_name = "co2contra_3d_colocated", &
!                     horizontal_staggering = "Ah", vertical_staggering = "None")
! if (parcomm_global%myid==0) then
!     print *, "co2contra_3d_colocated"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
! errs =  test_co2contra_3d(Nh = 32, Nz = 8, &
!                     co2contra_3d_oper_name = "co2contra_3d_Cgrid_h_sbp42_v_sbp42", &
!                     horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "co2contra_3d_C_CharneyPhilips"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
! errs =  test_co2contra_3d(Nh = 32, Nz = 8, &
!                     co2contra_3d_oper_name = "co2contra_3d_h_colocated_v_sbp42", &
!                     horizontal_staggering = "Ah", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "co2contra_3d_h_colocated_v_sbp42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if

!
! errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_C_up4", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "vector_advection_C_up4"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
! errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_C_up3", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "vector_advection_C_up3"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
! errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_C_up1", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "vector_advection_C_up1"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
! errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_C_c2", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "vector_advection_C_c2"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
! errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_C_c4", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "vector_advection_C_c4"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_Ah21", staggering="Ah")
! if (parcomm_global%myid==0) then
!     print *, "vector_advection_Ah21"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_Ah42", staggering="Ah")
! if (parcomm_global%myid==0) then
!     print *, "vector_advection_Ah42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs = test_vec_advection(N=32, vecadv_oper_name="vector_advection_Ah63", staggering="Ah")
! if (parcomm_global%myid==0) then
!     print *, "vector_advection_Ah63"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if

call deinit_global_parallel_enviroment()

end program test_diffops
