program main

    use test_namelist_mod, only : test_namelist
    use parcomm_mod,       only : init_global_parallel_enviroment, &
                                  deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    call test_namelist

    call deinit_global_parallel_enviroment()

end
