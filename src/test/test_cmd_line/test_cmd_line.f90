program test_cmd_line

    use cmd_args_test_mod, only : test_cmd_args
    use parcomm_mod,       only : init_global_parallel_enviroment, &
                                  deinit_global_parallel_enviroment

    implicit none

    integer(kind=4) ierr

    call init_global_parallel_enviroment()

    call test_cmd_args()

    call deinit_global_parallel_enviroment()

end program test_cmd_line
