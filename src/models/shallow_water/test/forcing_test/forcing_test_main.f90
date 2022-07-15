program forcing_test

use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment, parcomm_global

use forcing_test_mod, only : run_forcing_test

call init_global_parallel_enviroment()


call run_forcing_test()


call deinit_global_parallel_enviroment()

end program forcing_test
