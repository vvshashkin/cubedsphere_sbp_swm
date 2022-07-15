program test_diffops_3d

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment, parcomm_global
use test_diffops_3d_mod, only : test_w2uv_interp, test_uv2w_interp, test_scalar_advection_3d
use key_value_mod,       only : key_value_r8_t

implicit none

type(key_value_r8_t)  :: errs

call init_global_parallel_enviroment()

! errs =  test_w2uv_interp(Nh = 32, Nz = 8, &
!                          w2uv_interpolator_name = "w2uv_colocated", &
!                          w2uv_hor_part_name     = "None",           &
!                          w2uv_vert_part_name    = "None",           &
!                          horizontal_staggering  = "Ah", vertical_staggering = "None")
! if (parcomm_global%myid==0) then
!     print *, "w2uv_colocated"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_w2uv_interp(Nh = 32, Nz = 8, &
!                          w2uv_interpolator_name = "w2uv_hor_colocated",         &
!                          w2uv_hor_part_name     = "None",                       &
!                          w2uv_vert_part_name    = "vertical_interp_w2p_sbp21",  &
!                          horizontal_staggering = "Ah", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "w2uv_hor_colocated_sbp21"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_w2uv_interp(Nh = 32, Nz = 8, &
!                          w2uv_interpolator_name = "w2uv_hor_colocated",         &
!                          w2uv_hor_part_name     = "None",                       &
!                          w2uv_vert_part_name    = "vertical_interp_w2p_sbp42",  &
!                          horizontal_staggering = "Ah", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "w2uv_hor_colocated_sbp42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_w2uv_interp(Nh = 32, Nz = 8, &
!                          w2uv_interpolator_name = "w2uv_staggered",             &
!                          w2uv_hor_part_name     = "interp2d_p2uv_C_sbp42",     &
!                          w2uv_vert_part_name    = "vertical_interp_w2p_sbp21",  &
!                          horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "w2uv_staggered_C_sbp42_v_sbp21"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_w2uv_interp(Nh = 32, Nz = 8, &
!                          w2uv_interpolator_name = "w2uv_staggered",             &
!                          w2uv_hor_part_name     = "interp2d_p2uv_C_sbp42",      &
!                          w2uv_vert_part_name    = "vertical_interp_w2p_sbp42",  &
!                          horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "w2uv_staggered_C_sbp42_v_sbp42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if

errs =  test_uv2w_interp(Nh = 32, Nz = 8, &
                         uv2w_interpolator_name = "uv2w_colocated", &
                         uv2w_hor_part_name     = "",               &
                         uv2w_vert_part_name    = "",               &
                         horizontal_staggering = "Ah", vertical_staggering = "None")
if (parcomm_global%myid==0) then
    print *, "uv2w_colocated"
    print "(A,4E25.16)", "Err: ", errs%values
end if

errs =  test_uv2w_interp(Nh = 32, Nz = 8, &
                         uv2w_interpolator_name = "uv2w_hor_colocated",  &
                         uv2w_hor_part_name     = "",                          &
                         uv2w_vert_part_name    = "vertical_interp_p2w_sbp21", &
                         horizontal_staggering = "Ah", vertical_staggering = "CharneyPhilips")
if (parcomm_global%myid==0) then
    print *, "uv2w_hor_colocated_sbp21"
    print "(A,4E25.16)", "Err: ", errs%values
end if

errs =  test_uv2w_interp(Nh = 32, Nz = 8, &
                         uv2w_interpolator_name = "uv2w_hor_colocated",  &
                         uv2w_hor_part_name     = "",                          &
                         uv2w_vert_part_name    = "vertical_interp_p2w_sbp42", &
                         horizontal_staggering = "Ah", vertical_staggering = "CharneyPhilips")
if (parcomm_global%myid==0) then
    print *, "uv2w_hor_colocated_sbp42"
    print "(A,4E25.16)", "Err: ", errs%values
end if

errs =  test_uv2w_interp(Nh = 32, Nz = 8, &
                         uv2w_interpolator_name = "uv2w_staggered",                  &
                         uv2w_hor_part_name     = "interp2d_uv2pvec_C_sbp42",           &
                         uv2w_vert_part_name    = "vertical_interp_p2w_sbp21",       &
                         horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
if (parcomm_global%myid==0) then
    print *, "uv2w_staggered_C_sbp42_v_sbp21"
    print "(A,4E25.16)", "Err: ", errs%values
end if

errs =  test_uv2w_interp(Nh = 32, Nz = 8, &
                         uv2w_interpolator_name = "uv2w_staggered",                  &
                         uv2w_hor_part_name     = "interp2d_uv2pvec_C_sbp42",           &
                         uv2w_vert_part_name    = "vertical_interp_p2w_sbp42",       &
                         horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
if (parcomm_global%myid==0) then
    print *, "uv2w_staggered_C_sbp42_v_sbp42"
    print "(A,4E25.16)", "Err: ", errs%values
end if

errs =  test_scalar_advection_3d(Nh = 32, Nz = 10, &
                                 advection_oper_name      = "advection_w_staggered", &
                                 hor_advection_oper_name  = "up4",                   &
                                 vert_advection_oper_name = "adv_z_c2",              &
                                 points_type              = "w",                     &
                                 horizontal_staggering    = "C",                     &
                                 vertical_staggering      = "CharneyPhilips")
if (parcomm_global%myid==0) then
    print *, "advection_w_staggered, up4, adv_z_c2"
    print "(A,4E25.16)", "Err: ", errs%values
end if

errs =  test_scalar_advection_3d(Nh = 32, Nz = 10, &
                                 advection_oper_name      = "advection_p_staggered", &
                                 hor_advection_oper_name  = "up4",                   &
                                 vert_advection_oper_name = "adv_z_c2",              &
                                 points_type              = "p",                     &
                                 horizontal_staggering    = "C",                     &
                                 vertical_staggering      = "CharneyPhilips")
if (parcomm_global%myid==0) then
    print *, "advection_p_staggered, up4, adv_z_c2"
    print "(A,4E25.16)", "Err: ", errs%values
end if


call deinit_global_parallel_enviroment()

end program test_diffops_3d
