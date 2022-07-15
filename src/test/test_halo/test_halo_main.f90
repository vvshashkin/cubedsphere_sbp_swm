program test_halo_main

    use test_ecs_halo_mod,   only : test_ecs_halo
    use test_ecs_halo_c_mod, only : test_ecs_cvec_halo
    use test_halo_mod,       only : test_halo
    use parcomm_mod,         only : init_global_parallel_enviroment, &
                                    deinit_global_parallel_enviroment

    use test_generic_halo_mod, only : test_scalar_halo

    call init_global_parallel_enviroment()

    call test_halo()
    ! call test_ecs_halo
    call test_scalar_halo(nh=32,halo_width=2,staggering="A",halo_procedure_name="ECS_O", &
                          check_corners=.true.)
    call test_scalar_halo(nh=32,halo_width=3,staggering="Ah",halo_procedure_name="ECS_xy", &
                          check_corners=.true.)
    call test_ecs_cvec_halo("ecs_C_vec", "contravariant")
    call test_ecs_cvec_halo("ecs_C_vec_covariant", "covariant")

    call deinit_global_parallel_enviroment()

end
