program test_ts
    use test_rk4_mod, only: test_rk4

    use parcomm_mod,  only : init_global_parallel_enviroment, &
                             deinit_global_parallel_enviroment

implicit none

integer(kind=4) ierr

    call init_global_parallel_enviroment()

    call test_rk4()

    call deinit_global_parallel_enviroment()

end program test_ts
