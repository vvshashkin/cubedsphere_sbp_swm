program paneled_output_test

    use test_paneled_output_mod, only : test_master_paneled_output!, &
    !                                    test_mpi_paneled_output
    use parcomm_mod,             only : init_global_parallel_enviroment, &
                                        deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    call test_master_paneled_output()
    ! call test_mpi_paneled_output()

    call deinit_global_parallel_enviroment()

end program paneled_output_test
