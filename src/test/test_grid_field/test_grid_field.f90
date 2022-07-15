program test_grid_field_main

    use test_grid_field_mod, only : test_grid_field
    use parcomm_mod,         only : init_global_parallel_enviroment, &
                                    deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    call test_grid_field()

    call deinit_global_parallel_enviroment()

end
